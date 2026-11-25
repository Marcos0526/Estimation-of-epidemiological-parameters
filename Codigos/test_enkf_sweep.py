"""Regression checks for forecast timing, reproducibility, and sweep metrics."""

from functools import partial
import unittest

import numpy as np
import pandas as pd

from codes import gaussian_rate_walk, run_sird_enkf
from enkf_sweep import error_metrics, evaluate
from enkf_sweep_data import FILTER_SETTINGS, OBSERVED, synthetic_data


class EnkfSweepTests(unittest.TestCase):
    def test_forecast_is_opt_in_and_precedes_analysis(self):
        observations = synthetic_data()[OBSERVED].to_numpy()[:5]
        settings = FILTER_SETTINGS.copy()
        rate_step = partial(gaussian_rate_walk, log_step_std=settings.pop("log_step_std"))
        settings.update(ensemble_size=10, seed=2026, rate_step=rate_step)
        original = run_sird_enkf(observations, **settings)
        with_forecast = run_sird_enkf(observations, **settings, return_forecast=True)
        self.assertEqual(len(original), 3)
        for before, after in zip(original, with_forecast[:3]):
            np.testing.assert_array_equal(before, after)
        self.assertTrue(np.isnan(with_forecast[3][0]).all())

        changed = observations.copy()
        changed[2, 0] += 300
        updated = run_sird_enkf(changed, **settings, return_forecast=True)
        np.testing.assert_array_equal(with_forecast[3][:3], updated[3][:3])
        self.assertFalse(np.array_equal(with_forecast[0][2], updated[0][2]))

    def test_same_seed_reproduces_the_daily_series_and_metrics(self):
        data = synthetic_data().iloc[:20]
        first = evaluate(data, 10, 2026)
        second = evaluate(data, 10, 2026)
        for left, right in zip(first, second):
            pd.testing.assert_frame_equal(left, right)

    def test_metrics_skip_initial_day_and_use_calendar_groups(self):
        data = synthetic_data().iloc[:3].copy()
        data["period"] = ["First", "First", "Second"]
        series, _ = evaluate(data, 10, 2026)
        series["forecast_I"] = series["observed_I"] + [1000, 3, 4]
        series["forecast_daily_deaths"] = series["observed_daily_deaths"] + [0, 1, -3]
        series["analysis_beta"] = 2 * series["true_beta"]
        series["beta_p05"] = series["true_beta"] * [0.8, 0.8, 1.1]
        series["beta_p95"] = series["true_beta"] * 1.2
        metrics = error_metrics(series).set_index(["period", "metric"])
        self.assertAlmostEqual(metrics.loc[("Overall", "forecast_I_rmse"), "value"], np.sqrt(12.5))
        self.assertAlmostEqual(metrics.loc[("Second", "forecast_I_mae"), "value"], 4)
        self.assertAlmostEqual(metrics.loc[("Overall", "forecast_daily_deaths_mae"), "value"], 2)
        self.assertAlmostEqual(metrics.loc[("Overall", "beta_log_rmse"), "value"], np.log(2))
        self.assertAlmostEqual(metrics.loc[("Overall", "beta_p05_p95_coverage"), "value"], 0.5)


if __name__ == "__main__":
    unittest.main()
