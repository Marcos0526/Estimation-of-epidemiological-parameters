"""Run and record the Mexico and synthetic SIRD ensemble size sweeps."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from hashlib import sha256
from importlib.metadata import version
from pathlib import Path
import platform
import subprocess
from tempfile import TemporaryDirectory
from time import perf_counter
import traceback

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

from codes import gaussian_rate_walk, run_sird_enkf
from enkf_sweep_data import (
    CHUNK_SIZE, DATA_DIR, FILTER_SETTINGS, LAST_DAY, OBSERVED, PHASE_FILE, RATES,
    RECOVERY_DAYS, SOURCE_FILE, START_DATE, STATES, SYNTHETIC_SETTINGS, mexico_data,
    synthetic_data,
)
from enkf_sweep_report import write_report


EXPERIMENT = "SIRD EnKF ensemble size"
REFERENCE = 160
SIZES = tuple(range(10, 201, 10))
SEEDS = tuple(range(2026, 2036))


def error_metrics(series):
    """Evaluate aligned observations, forecasts, and analyses from day 1."""
    evaluated = series.iloc[1:]
    errors = {}
    for stage in ("forecast", "analysis"):
        for state in OBSERVED:
            errors[f"{stage}_{state}"] = (
                evaluated[f"{stage}_{state}"] - evaluated[f"observed_{state}"]
            )
    errors["forecast_daily_deaths"] = (
        evaluated["forecast_daily_deaths"] - evaluated["observed_daily_deaths"]
    )
    errors = pd.DataFrame(errors)
    groups = [("Overall", errors)]
    groups.extend(errors.groupby(evaluated["period"], sort=False))
    rows = []
    for period, group in groups:
        for metric, values in group.abs().mean().items():
            rows.append((period, f"{metric}_mae", values, len(group)))
        for metric, values in np.sqrt(group.pow(2).mean()).items():
            rows.append((period, f"{metric}_rmse", values, len(group)))

    for rate in RATES:
        truth_column = f"true_{rate}"
        if truth_column in evaluated:
            truth = evaluated[truth_column]
            error = np.log(evaluated[f"analysis_{rate}"]) - np.log(truth)
            covered = truth.between(evaluated[f"{rate}_p05"], evaluated[f"{rate}_p95"])
            rows.append(("Overall", f"{rate}_log_rmse", np.sqrt(np.mean(error**2)), len(truth)))
            rows.append(("Overall", f"{rate}_p05_p95_coverage", covered.mean(), len(truth)))
    return pd.DataFrame(rows, columns=["period", "metric", "value", "days"])


def evaluate(data, size, seed):
    """Run the filter and return daily and aggregate metrics for one configuration."""
    settings = FILTER_SETTINGS.copy()
    rate_step = partial(gaussian_rate_walk, log_step_std=settings.pop("log_step_std"))
    mean, std, paths, forecast = run_sird_enkf(
        data[OBSERVED].to_numpy(dtype=float), ensemble_size=size, seed=seed,
        rate_step=rate_step, days=data["day"].to_numpy(), return_forecast=True, **settings,
    )
    population = settings["population"]
    if not all(np.isfinite(array).all() for array in (mean, std, paths, forecast[1:])):
        raise FloatingPointError("The filter returned non-finite results.")
    if (paths <= 0).any() or (mean[:, :4] < 0).any():
        raise FloatingPointError("The filter returned non-positive rates or negative counts.")
    if not np.allclose(mean[:, :4].sum(axis=1), population, rtol=0, atol=1e-5):
        raise FloatingPointError("The filtered mean does not conserve the population.")

    columns = [*STATES, *RATES]
    series = data.rename(columns={state: f"observed_{state}" for state in OBSERVED}).copy()
    series[[f"analysis_{column}" for column in columns]] = mean
    series[[f"analysis_sd_{column}" for column in columns]] = std
    series[[f"forecast_{state}" for state in OBSERVED]] = forecast
    for percentile, values in zip((5, 50, 95), np.percentile(paths, [5, 50, 95], axis=1)):
        series[[f"{rate}_p{percentile:02d}" for rate in RATES]] = values
    series["observed_daily_deaths"] = series["observed_D"].diff()
    series["forecast_daily_deaths"] = series["forecast_D"] - series["analysis_D"].shift(1)
    metrics = error_metrics(series)
    if not np.isfinite(metrics["value"]).all():
        raise FloatingPointError("The calculated metrics are not finite.")
    return series, metrics


def run_job(data, size, seed):
    start = perf_counter()
    try:
        series, metrics = evaluate(data, size, seed)
    except Exception:
        return perf_counter() - start, None, None, traceback.format_exc()
    return perf_counter() - start, series, metrics, None


def completed_jobs(jobs, datasets, workers):
    if workers == 1:
        for scenario, size, seed in jobs:
            yield (scenario, size, seed), run_job(datasets[scenario], size, seed)
        return
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(run_job, datasets[scenario], size, seed): (scenario, size, seed)
            for scenario, size, seed in jobs
        }
        for future in as_completed(futures):
            yield futures.pop(future), future.result()


def file_hash(path):
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def code_revision():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=DATA_DIR.parent, text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def prepare_inputs(mlflow, scenarios, folder):
    configuration = {
        "filter": FILTER_SETTINGS,
        "code_revision": code_revision(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "dependencies": {
            name: version(name) for name in ("numpy", "pandas", "scipy", "matplotlib", "mlflow")
        },
        "datasets": {},
        "code_sha256": {},
    }
    for name in ("codes.py", "enkf_sweep.py", "enkf_sweep_data.py", "enkf_sweep_report.py"):
        source = Path(__file__).with_name(name)
        configuration["code_sha256"][name] = file_hash(source)
        mlflow.log_artifact(str(source), "code")

    datasets = {}
    for scenario in scenarios:
        print(f"Preparing {scenario} observations", flush=True)
        if scenario == "mexico":
            details = {
                "start_date": START_DATE,
                "last_day": LAST_DAY,
                "recovery_days": RECOVERY_DAYS,
                "chunk_size": CHUNK_SIZE,
                "source_sha256": file_hash(SOURCE_FILE),
                "phase_calendar_sha256": file_hash(PHASE_FILE),
            }
            data = mexico_data()
        else:
            details = SYNTHETIC_SETTINGS.copy()
            data = synthetic_data()
        file = folder / f"{scenario}.csv"
        data.to_csv(file, index=False)
        configuration["datasets"][scenario] = {**details, "daily_inputs_sha256": file_hash(file)}
        mlflow.log_artifact(str(file), "inputs")
        datasets[scenario] = data
    return datasets, configuration


def unique_values(values):
    return tuple(dict.fromkeys(values))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", action="store_true", help="Run sizes 10, 160, and 200 with seed 2026.")
    parser.add_argument("--sizes", nargs="+", type=int, choices=SIZES)
    parser.add_argument("--seeds", nargs="+", type=int)
    parser.add_argument("--scenario", choices=("both", "mexico", "synthetic"), default="both")
    parser.add_argument("--workers", type=int, default=1, help="Number of independent filter processes.")
    parser.add_argument("--output-dir", type=Path, default=DATA_DIR / "enkf_sweep")
    args = parser.parse_args()
    if args.workers < 1 or (args.seeds and min(args.seeds) < 0):
        parser.error("Workers must be positive and seeds must be non-negative.")
    if args.pilot and (args.sizes or args.seeds):
        parser.error("Use either --pilot or custom sizes and seeds.")

    sizes = (10, REFERENCE, 200) if args.pilot else unique_values(args.sizes or SIZES)
    seeds = (SEEDS[0],) if args.pilot else unique_values(args.seeds or SEEDS)
    scenarios = ("mexico", "synthetic") if args.scenario == "both" else (args.scenario,)
    jobs = [(scenario, size, seed) for scenario in scenarios for size in sizes for seed in seeds]
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    import mlflow

    mlflow.set_tracking_uri(f"sqlite:///{output / 'mlflow.db'}")
    if mlflow.get_experiment_by_name(EXPERIMENT) is None:
        mlflow.create_experiment(EXPERIMENT, artifact_location=(output / "artifacts").as_uri())
    mlflow.set_experiment(EXPERIMENT)
    label = "pilot" if args.pilot else "sweep"
    rows, trajectories = [], []
    with mlflow.start_run(run_name=label, tags={"kind": "sweep"}) as parent:
        print(f"Sweep ID: {parent.info.run_id} | SQLite: {output / 'mlflow.db'}", flush=True)
        with TemporaryDirectory(prefix="enkf-sweep-") as temporary:
            folder = Path(temporary)
            datasets, configuration = prepare_inputs(mlflow, scenarios, folder)
            mlflow.log_dict({
                **configuration, "sizes": sizes, "seeds": seeds, "scenarios": scenarios,
                "workers": args.workers, "reference_size": REFERENCE,
            }, "configuration.json")
            mlflow.log_params({"reference_size": REFERENCE, "scheduled_runs": len(jobs)})

            for number, (job, result) in enumerate(completed_jobs(jobs, datasets, args.workers), 1):
                scenario, size, seed = job
                duration, series, metrics, error = result
                status = "FAILED" if error else "FINISHED"
                params = {"scenario": scenario, "ensemble_size": size, "seed": seed}
                values = {"runtime_seconds": duration}
                with mlflow.start_run(
                    run_name=f"{scenario}-{size}-seed-{seed}", nested=True, tags={"kind": "filter"},
                ) as run:
                    mlflow.log_params({
                        **params, **FILTER_SETTINGS,
                        "inputs_sha256": configuration["datasets"][scenario]["daily_inputs_sha256"],
                    })
                    if error:
                        mlflow.log_text(error, "error.txt")
                        mlflow.set_tag("error_message", error.strip().splitlines()[-1])
                    else:
                        series.to_csv(folder / "daily.csv", index=False)
                        metrics.to_csv(folder / "metrics.csv", index=False)
                        mlflow.log_artifact(str(folder / "daily.csv"))
                        mlflow.log_artifact(str(folder / "metrics.csv"))
                        values.update(metrics.loc[metrics["period"] == "Overall"].set_index("metric")["value"])
                        trajectories.append(series[["day", *[f"analysis_{r}" for r in RATES]]].assign(**params))
                    mlflow.log_metrics(values)
                    mlflow.end_run(status=status)
                rows.append({
                    **params, **values, "status": status, "run_id": run.info.run_id,
                    "error": error.strip().splitlines()[-1] if error else "",
                })
                print(f"[{number}/{len(jobs)}] {scenario} M={size} seed={seed}: {status} ({duration:.2f}s)", flush=True)

        report = output / "reports" / parent.info.run_id
        write_report(report, pd.DataFrame(rows), trajectories, REFERENCE)
        mlflow.log_artifacts(str(report), "report")
        failures = sum(row["status"] == "FAILED" for row in rows)
        mlflow.log_metrics({"completed_runs": len(rows) - failures, "failed_runs": failures})
        print(f"Report: {report} | {len(rows) - failures} finished, {failures} failed", flush=True)
    return int(len(rows) == failures)


if __name__ == "__main__":
    raise SystemExit(main())
