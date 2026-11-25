"""Build the fixed datasets and configure the ensemble size experiment."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from codes import SIRD, get_Y_obs


DATA_DIR = Path(__file__).resolve().parents[1] / "Datos"
SOURCE_FILE = DATA_DIR / "raw_data" / "COVID19MEXICO2020_2021_2022_2023.csv"
PHASE_FILE = DATA_DIR / "metadata" / "covid_waves_parameters.csv"
START_DATE = "2020-02-28"
LAST_DAY = 1166
CHUNK_SIZE = 250_000
RECOVERY_DAYS = 14
OBS_STD = (10.0, 4.0, 1.5)
FILTER_SETTINGS = {
    "population": 126_014_024,
    "state_std": (1.0, 0.0, 0.0),
    "initial_rates": (0.3508, 0.0714, 0.0015),
    "initial_log_std": (0.35, 0.40, 0.50),
    "log_step_std": (0.12, 0.12, 0.15),
    "observation_cov": np.diag(np.square(OBS_STD)).tolist(),
}
OBSERVED = ["I", "R", "D"]
STATES = ["S", *OBSERVED]
RATES = ["beta", "gamma", "delta"]
SYNTHETIC_SETTINGS = {
    "last_day": 240,
    "initial_I_R_D": (10_000, 50_000, 1_000),
    "base_rates": (0.09, 0.0714, 0.0015),
    "log_amplitudes": (0.35, 0.20, 0.30),
    "periods_days": (120, 80, 60),
    "phase_radians": (0.0, 1.0, 2.0),
    "observation_noise_seed": 2025,
    "solver_rtol": 1e-10,
    "solver_atol": 1e-9,
}


def mexico_data(source=SOURCE_FILE, phase_file=PHASE_FILE):
    """Reconstruct the observations with the source notebook's rules."""
    dates = pd.date_range(START_DATE, periods=LAST_DAY + 1, name="date")
    start, end = dates[0], dates[-1]
    phases = pd.read_csv(phase_file, parse_dates=["Start Date", "End Date"])
    periods = pd.Series("Between waves", index=dates, name="period")
    for name, first, last in phases[
        ["Epidemiological Phase", "Start Date", "End Date"]
    ].itertuples(index=False, name=None):
        periods.loc[first:last] = name
    periods.loc[dates > phases["End Date"].max()] = "After sixth wave"
    if not set(phases["Epidemiological Phase"]) <= set(periods):
        raise ValueError("At least one source phase has no dates in the calendar.")

    observed = pd.DataFrame(0, index=dates, columns=[*OBSERVED, "C"], dtype="int64")
    record_counts = pd.Series(0, index=dates, dtype="int64")
    with pd.read_csv(
        source, usecols=["FECHA_INGRESO", "FECHA_DEF", "CLASIFICACION_FINAL"],
        chunksize=CHUNK_SIZE, na_values={"FECHA_DEF": ["9999-99-99"]},
    ) as reader:
        for chunk in reader:
            admission = pd.to_datetime(chunk["FECHA_INGRESO"], format="%Y-%m-%d")
            records = chunk.loc[admission.between(start, end)].copy()
            if records.empty:
                continue
            counts = admission.loc[records.index].value_counts()
            record_counts.loc[counts.index] += counts
            records["FECHA_DEF"] = pd.to_datetime(records["FECHA_DEF"], format="%Y-%m-%d")
            _, daily = get_Y_obs(records, start, end, gamma_time=RECOVERY_DAYS)
            observed += daily[observed.columns]

    if not record_counts.gt(0).all():
        raise ValueError("The source has dates without admission records.")
    if not (observed["C"] == observed[OBSERVED].sum(axis=1)).all():
        raise ValueError("Cumulative cases do not equal I + R + D.")
    if not observed[["C", "R", "D"]].diff().iloc[1:].ge(0).all().all():
        raise ValueError("Observed cumulative counts decreased.")
    return observed.assign(
        day=np.arange(len(dates)), period=periods, admission_records=record_counts
    ).reset_index()


def synthetic_data():
    """Generate one continuous SIRD trajectory and one fixed observation series."""
    config = SYNTHETIC_SETTINGS
    days = np.arange(config["last_day"] + 1)
    population = FILTER_SETTINGS["population"]
    initial_counts = config["initial_I_R_D"]
    base = np.array(config["base_rates"])
    amplitude = np.array(config["log_amplitudes"])
    frequency = 2 * np.pi / np.array(config["periods_days"])
    phase = np.array(config["phase_radians"])

    def rates(time):
        return base * np.exp(amplitude * np.sin(np.asarray(time)[..., None] * frequency + phase))

    def derivative(time, state):
        return SIRD(time, state, *rates(time), population)

    solution = solve_ivp(
        derivative, (days[0], days[-1]), (population - sum(initial_counts), *initial_counts),
        t_eval=days, rtol=config["solver_rtol"], atol=config["solver_atol"],
    )
    truth = solution.y.T
    if not solution.success or not np.isfinite(truth).all():
        raise RuntimeError(f"Synthetic integration failed: {solution.message}")
    if (truth < 0).any() or not np.allclose(truth.sum(axis=1), population, rtol=0, atol=1e-5):
        raise RuntimeError("Synthetic compartments violate the population constraints.")

    rng = np.random.default_rng(config["observation_noise_seed"])
    observed = truth[:, 1:] + rng.normal(scale=OBS_STD, size=(len(days), len(OBSERVED)))
    data = pd.DataFrame(observed, columns=OBSERVED).assign(day=days, period="Synthetic")
    data[[f"true_{state}" for state in STATES]] = truth
    data[[f"true_{rate}" for rate in RATES]] = rates(days)
    return data
