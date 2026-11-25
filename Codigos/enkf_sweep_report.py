"""Compare successful replicates and report numerical failures by ensemble size."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from enkf_sweep_data import RATES


def metric_plot(summary, scenario, folder):
    metrics = [
        ("forecast_I_rmse", "Active cases RMSE (people)"),
        ("forecast_R_rmse", "Cumulative recoveries RMSE (people)"),
        ("forecast_D_rmse", "Cumulative deaths RMSE (people)"),
        ("forecast_daily_deaths_rmse", "Daily deaths RMSE (people)"),
    ]
    if scenario == "synthetic":
        metrics += [(f"{rate}_log_rmse", f"log({rate}) RMSE") for rate in RATES]
    metrics += [("runtime_seconds", "Runtime (s)")]
    fig, axes = plt.subplots((len(metrics) + 1) // 2, 2, figsize=(10, 2.7 * ((len(metrics) + 1) // 2)),
                             layout="constrained", squeeze=False)
    for ax, (name, title) in zip(axes.flat, metrics):
        values = summary.loc[summary["metric"] == name].sort_values("ensemble_size")
        ax.plot(values["ensemble_size"], values["median"], marker="o", markersize=3)
        ax.fill_between(values["ensemble_size"], values["q25"], values["q75"], alpha=0.22)
        ax.set(title=title, xlabel="Ensemble members")
        ax.grid(axis="y", alpha=0.2)
    for ax in list(axes.flat)[len(metrics):]:
        ax.set_visible(False)
    fig.suptitle(f"{scenario.title()} | median and interquartile range across seeds")
    fig.savefig(folder / f"{scenario}_metrics.png", dpi=160)
    plt.close(fig)


def rate_spread_plot(spread, scenario, folder):
    fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True, layout="constrained")
    for ax, rate in zip(axes, RATES):
        values = spread.loc[spread["rate"] == rate].copy()
        values["log_iqr"] = np.log10(values["q75"]) - np.log10(values["q25"])
        grid = values.pivot(index="ensemble_size", columns="day", values="log_iqr")
        mesh = ax.pcolormesh(grid.columns, np.arange(len(grid)), grid.to_numpy(), shading="nearest")
        ax.set_yticks(np.arange(len(grid))[::2], labels=grid.index[::2])
        ax.set(title=rate, ylabel="Members")
        fig.colorbar(mesh, ax=ax, label="IQR in log10(rate)")
    axes[-1].set_xlabel("Day")
    fig.suptitle(f"{scenario.title()} | variation of the filtered rate across seeds")
    fig.savefig(folder / f"{scenario}_rate_spread.png", dpi=160)
    plt.close(fig)


def write_report(folder, runs, trajectories, reference):
    folder.mkdir(parents=True, exist_ok=True)
    keys = ["scenario", "ensemble_size"]
    runs = runs.sort_values([*keys, "seed"])
    runs.to_csv(folder / "runs.csv", index=False)
    finished = runs["status"].eq("FINISHED")
    counts = runs.assign(finished=finished, failed=~finished).groupby(keys)[["finished", "failed"]].sum()
    counts.to_csv(folder / "counts.csv")

    metric_columns = runs.select_dtypes("number").columns.difference(["ensemble_size", "seed"])
    rows = []
    for scenario, size in counts.index:
        group = runs.loc[finished & runs["scenario"].eq(scenario) & runs["ensemble_size"].eq(size)]
        for metric in metric_columns:
            values = group[metric].dropna()
            # Synthetic truth metrics are undefined for the Mexico scenario.
            if scenario == "mexico" and metric.startswith(tuple(f"{rate}_" for rate in RATES)):
                continue
            q25, median, q75 = values.quantile([0.25, 0.5, 0.75])
            rows.append((scenario, size, metric, q25, median, q75, len(values)))
    summary = pd.DataFrame(rows, columns=[*keys, "metric", "q25", "median", "q75", "count"])
    summary = summary.merge(counts.reset_index(), on=keys, how="left")
    baseline = summary.loc[summary["ensemble_size"].eq(reference)].set_index(["scenario", "metric"])["median"]
    summary["median_difference_from_reference"] = (
        summary["median"] - summary.set_index(["scenario", "metric"]).index.map(baseline)
    )
    summary.to_csv(folder / "summary.csv", index=False)

    for scenario, group in summary.groupby("scenario"):
        if group["count"].gt(0).any():
            metric_plot(group, scenario, folder)
    if not trajectories:
        return

    rate_data = pd.concat(trajectories, ignore_index=True).rename(
        columns={f"analysis_{rate}": rate for rate in RATES}
    )
    long_rates = rate_data.melt(id_vars=[*keys, "day"], value_vars=RATES,
                                var_name="rate", value_name="estimate")
    grouped = long_rates.groupby([*keys, "day", "rate"])["estimate"]
    spread = grouped.quantile([0.25, 0.5, 0.75]).unstack().set_axis(["q25", "median", "q75"], axis=1)
    spread["count"] = grouped.count()
    spread = spread.reset_index()
    spread.to_csv(folder / "rate_spread.csv", index=False)
    for scenario, group in spread.loc[spread["count"] >= 2].groupby("scenario"):
        rate_spread_plot(group, scenario, folder)
