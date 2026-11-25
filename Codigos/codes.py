import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates  
from scipy.integrate import solve_ivp

################ 
# To CSV

def dates_to_datetime(df, columns): # Transform date columns to datetime format
    df['FECHA_DEF'] = pd.to_datetime(df['FECHA_DEF'], errors='coerce')
    for col in columns:
        df[col] = pd.to_datetime(df[col], format='%Y-%m-%d')
    return df

def covid_positive(df): # Filter the DataFrame to include only positive COVID-19 cases
    return df[(df['CLASIFICACION_FINAL'] == 1) | (df['CLASIFICACION_FINAL'] == 2) | (df['CLASIFICACION_FINAL'] == 3)]

def filter_by_date(df, first_date, last_date, column_name='FECHA_INGRESO',positive=True, saved=False): # Filter the DataFrame by date range
    if positive:
        df = covid_positive(df)
    df = dates_to_datetime(df, [column_name])
    first_date = pd.to_datetime(first_date)
    last_date = pd.to_datetime(last_date)

    filtro = df[column_name].between(
        first_date,
        last_date,
        inclusive="both"
    )

    out = df.loc[filtro].copy()

    if saved:
        filename = (
            f'filtered_data_{first_date:%Y-%m-%d}_'
            f'{last_date:%Y-%m-%d}.csv'
        )
        out.to_csv(filename, index=False)

    return out

def formatear_fechas(ax, intervalo_meses=2): # Format dates in plots
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=intervalo_meses))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))   # e.g., "Mar 2020"
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.tick_params(axis="x", rotation=45)
    ax.grid(alpha=0.3)

def plot_log(block, date_1, date_f, intervalo_meses=1): # Plot cumulative positive cases on a logarithmic scale

    block = dates_to_datetime(block, ['FECHA_INGRESO'])

    pos = covid_positive(block)
    idx = pd.date_range(date_1, date_f, freq='D')

    new = pos.groupby('FECHA_INGRESO').size().reindex(idx, fill_value=0)

    C = new.cumsum()      # Cumulative cases (= I + R + D)

    print('Casos positivos acumulados: ', C.iloc[-1])

    acumulado = pd.DataFrame({'Fecha': idx, 'Casos positivos acumulados': C.values})


    fig, ax = plt.subplots(figsize=(12, 4))
    ax.set_yscale('log')
    ax.plot(C)
    ax.set_title('Casos acumulados')
    ax.set_xlabel('Fecha')
    ax.set_ylabel('Número de casos')
    formatear_fechas(ax, intervalo_meses)
    plt.tight_layout()
    plt.show()

    return C,idx,acumulado


def plot_cumulated(block, date_1, date_f, intervalo_meses=1):# Plot cumulative positive cases and deaths

        block = dates_to_datetime(block, ['FECHA_INGRESO', 'FECHA_SINTOMAS', 'FECHA_DEF'])

        pos = covid_positive(block)
        idx = pd.date_range(date_1, date_f, freq='D')

        new    = pos.groupby('FECHA_INGRESO').size().reindex(idx, fill_value=0)
        deaths = pos['FECHA_DEF'].dropna().value_counts().reindex(idx, fill_value=0)

        C = new.cumsum()      # Cumulative cases (= I + R + D)
        D = deaths.cumsum()   # Cumulative deaths

        print('Casos positivos acumulados: ', C.iloc[-1])
        print('Defunciones acumuladas: ', D.iloc[-1])

        fig, ax = plt.subplots(1, 2, figsize=(12, 4))
        ax[0].plot(C)
        ax[0].set_title('Casos acumulados')
        ax[1].plot(D)
        ax[1].set_title('Defunciones acumuladas')
        for a in ax:
            formatear_fechas(a, intervalo_meses)
        plt.tight_layout()
        plt.show()

        return C, D,idx

#################
# SIRD simulation

def SIRD(t, y, beta, gamma, delta, N): # Define the SIRD model differential equations
    S, I, R, D = y
    infecciones = beta * S * I / N
    dS = -infecciones
    dI = infecciones - gamma * I - delta * I
    dR = gamma * I
    dD = delta * I

    return dS, dI, dR, dD


def SIRD_simulation(beta, gamma, delta, N, y0, t_max, dt=1.0, verbose=True):


    """
    Returns
    --------
    t  : np.ndarray, shape (n,)    -> days
    y  : np.ndarray, shape (n, 4)  -> columns [S, I, R, D]
    df : pd.DataFrame              -> same data with named columns
    """
    t_eval = np.arange(0, t_max + dt, dt)
    sol = solve_ivp(
        SIRD,
        t_span=(0, t_max),
        y0=y0,
        args=(beta, gamma, delta, N),
        t_eval=t_eval,
        method='RK45',
        rtol=1e-8,
        atol=1e-8
    )

    if not sol.success:
        raise RuntimeError(f"solve_ivp falló: {sol.message}")

    t, y = sol.t, sol.y.T
    df = pd.DataFrame(y, columns=["S", "I", "R", "D"])
    df.insert(0, "dia", t)

    if verbose:
        R0 = beta / (gamma + delta)
        error_N = np.abs(y.sum(axis=1) - N).max()
        print(f"R0 = {R0:.3f}")
        print(f"Pico de infectados: {df['I'].max():,.0f} (día {df.loc[df['I'].idxmax(), 'dia']:.0f})")
        print(f"Muertes finales: {df['D'].iloc[-1]:,.0f}")
        print(f"Error máximo de conservación |S+I+R+D - N|: {error_N:.2e}")

    return t, y, df

###############
# SIRD observations

def get_Y_obs(df, date_1, date_f, N=None, gamma_intervals=None,
              gamma_time=14, date_col='FECHA_INGRESO'):
    """
    Y      : np.ndarray (T, 3) with columns [I, R, D]
    df_obs : pd.DataFrame with date index and columns [I, R, D, C] (+ S if N given)

    gamma_intervals : list of (start, end, gamma) tuples, e.g.
                      [('2020-02-28', '2020-03-22', 0.0714), ...]
                      If None, falls back to a constant gamma_time (in days).
    gamma_time      : constant recovery period in days, used only when
                      gamma_intervals is None.
    """
    df = dates_to_datetime(df.copy(), [date_col, 'FECHA_DEF'])
    pos = covid_positive(df)
    idx = pd.date_range(date_1, date_f, freq='D')

    # Cumulative cases C = I + R + D
    new = pos.groupby(date_col).size().reindex(idx, fill_value=0)
    C = new.cumsum()

    # Cumulative deaths
    D = (pos['FECHA_DEF'].dropna()
             .value_counts()
             .reindex(idx, fill_value=0)
             .cumsum())

    # Survivors: cases without a recorded death
    surv = pos[pos['FECHA_DEF'].isna()].copy()

    if gamma_intervals is not None:
        # Build a Series: date -> T_rec (in days)
        t_series = pd.Series(np.nan, index=idx, dtype=float)
        for start, end, gamma in gamma_intervals:
            mask = (idx >= pd.Timestamp(start)) & (idx <= pd.Timestamp(end))
            t_series[mask] = round(1.0 / gamma)

        if t_series.isna().any():
            missing = idx[t_series.isna()]
            raise ValueError(
                f"No gamma defined for dates: "
                f"{missing[0].date()} ... {missing[-1].date()}"
            )

        # Assign T_rec to each survivor based on their registration date
        surv_dates = surv[date_col].clip(
            upper=pd.Timestamp(date_f)
        )
        t_map = t_series.to_dict()
        t_rec_per_case = surv_dates.map(
            lambda d: t_map.get(d, np.nan)
        )

        # Recovery date = registration date + T_rec for that interval
        rec_dates = surv_dates + pd.to_timedelta(t_rec_per_case, unit='D')

    else:
        # Constant recovery period
        surv_dates = surv[date_col].clip(upper=pd.Timestamp(date_f))
        rec_dates = surv_dates + pd.Timedelta(days=gamma_time)

    R = (rec_dates.value_counts()
             .reindex(idx, fill_value=0)
             .cumsum())

    I = C - R - D

    df_obs = pd.DataFrame({'I': I, 'R': R, 'D': D, 'C': C}, index=idx)
    if N is not None:
        df_obs['S'] = N - C

    Y = df_obs[['I', 'R', 'D']].to_numpy(dtype=float)
    return Y, df_obs


# Filtro de Kalman por ensamble para [S, I, R, D, β, γ, δ].
SIRD_OBSERVATION_MATRIX = np.eye(7)[1:4]


def _rate_triplet(values, name, *, positive=False):
    """Valida tres valores en el orden β, γ, δ."""
    values = np.asarray(values, dtype=float)
    if values.shape != (3,) or not np.isfinite(values).all():
        raise ValueError(f"{name} debe contener tres valores finitos.")
    if (values < 0).any() or (positive and (values == 0).any()):
        raise ValueError(f"Revise el signo de {name}.")
    return values


def _balance_compartments(ensemble, population):
    """Reconstruye S por balance después de imponer I, R, D ≥ 0."""
    observed = ensemble[:, 1:4]
    np.maximum(observed, 0, out=observed)
    totals = observed.sum(axis=1)
    excess = totals > population
    observed[excess] *= (population / totals[excess])[:, None]
    ensemble[:, 0] = np.maximum(population - observed.sum(axis=1), 0)


def initialize_sird_ensemble(first_observation, population, ensemble_size,
                             state_std, initial_rates, initial_log_std, rng):
    """Genera un ensamble (M, 7) centrado antes de las restricciones en el primer dato."""
    first = np.asarray(first_observation, dtype=float)
    state_std = np.asarray(state_std, dtype=float)
    if not np.isfinite(population) or population <= 0:
        raise ValueError("La población debe ser finita y positiva.")
    if not isinstance(ensemble_size, (int, np.integer)) or ensemble_size < 2:
        raise ValueError("Se necesitan al menos dos miembros.")
    if first.shape != (3,) or not np.isfinite(first).all() or (first < 0).any() or first.sum() > population:
        raise ValueError("La primera observación [I, R, D] es incompatible con la población.")
    if state_std.shape != (3,) or not np.isfinite(state_std).all() or (state_std < 0).any():
        raise ValueError("Se requieren tres dispersiones no negativas para I, R y D.")
    rates = _rate_triplet(initial_rates, "las medianas iniciales", positive=True)
    log_std = _rate_triplet(initial_log_std, "la dispersión inicial de las log-tasas")

    ensemble = np.empty((ensemble_size, 7))
    ensemble[:, 1:4] = rng.normal(first, state_std, size=(ensemble_size, 3))
    _balance_compartments(ensemble, population)
    # En escala logarítmica la dispersión es relativa y las tasas no se acumulan
    # en un límite artificial por truncamiento. La salida conserva sus unidades físicas.
    ensemble[:, 4:] = rng.lognormal(np.log(rates), log_std, size=(ensemble_size, 3))
    return ensemble


def gaussian_rate_walk(rates, start_day, end_day, rng, *, log_step_std):
    """Aplica incrementos gaussianos independientes a las log-tasas por raíz de día."""
    sigma = _rate_triplet(log_step_std, "la dispersión del camino aleatorio")
    dt = end_day - start_day
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("El paso de tiempo debe ser finito y positivo.")
    increments = rng.normal(0, sigma * np.sqrt(dt), size=rates.shape)
    return rates * np.exp(increments)


def predict_sird_ensemble(ensemble, population, start_day, end_day, rate_step, rng):
    """Avanza las tasas con la ley suministrada e integra el SIRD de cada miembro."""
    forecast = np.empty_like(ensemble, dtype=float)
    rates = np.asarray(rate_step(ensemble[:, 4:].copy(), start_day, end_day, rng),
                       dtype=float)
    if rates.shape != ensemble[:, 4:].shape or not np.isfinite(rates).all() or (rates <= 0).any():
        raise ValueError("El avance debe devolver una matriz (M, 3) de tasas finitas y positivas.")
    forecast[:, 4:] = rates
    for member, (state, parameters) in enumerate(zip(ensemble[:, :4], rates)):
        solution = solve_ivp(SIRD, (start_day, end_day), state,
                             args=(*parameters, population), t_eval=[end_day],
                             rtol=1e-7, atol=1e-8)
        if not solution.success:
            raise RuntimeError(f"Falló la integración del miembro {member} en el día {end_day}.")
        forecast[member, :4] = solution.y[:, -1]
    return forecast


def analyze_sird_ensemble(forecast, observation, population, observation_cov, rng):
    """Actualiza estados y log-tasas con un EnKF de observaciones perturbadas."""
    covariance = np.asarray(observation_cov, dtype=float)
    if (covariance.shape != (3, 3) or not np.isfinite(covariance).all()
            or not np.allclose(covariance, covariance.T, rtol=1e-12, atol=0)):
        raise ValueError("La covarianza de observación debe ser simétrica de tamaño 3 × 3.")
    try:
        noise_factor = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as error:
        raise ValueError("La covarianza de observación debe ser definida positiva.") from error

    analysis = forecast.copy()
    analysis[:, 4:] = np.log(forecast[:, 4:])
    predicted_obs = analysis @ SIRD_OBSERVATION_MATRIX.T
    anomalies = analysis - analysis.mean(axis=0)
    obs_anomalies = predicted_obs - predicted_obs.mean(axis=0)
    divisor = len(forecast) - 1
    cross_cov = anomalies.T @ obs_anomalies / divisor
    innovation_cov = obs_anomalies.T @ obs_anomalies / divisor + covariance
    gain = np.linalg.solve(innovation_cov, cross_cov.T).T

    noise = rng.normal(size=predicted_obs.shape) @ noise_factor.T
    analysis += (observation + noise - predicted_obs) @ gain.T
    analysis[:, 4:] = np.exp(analysis[:, 4:])
    if not np.isfinite(analysis).all() or (analysis[:, 4:] <= 0).any():
        raise FloatingPointError("La actualización excedió el rango numérico de las tasas o estados.")
    _balance_compartments(analysis, population)
    return analysis


def run_sird_enkf(observations, population, *, ensemble_size, state_std, initial_rates,
                  initial_log_std, rate_step, observation_cov, seed, days=None,
                  return_forecast=False):
    """Filtra Y (T, 3). Devuelve medias, desviaciones (T, 7) y tasas (T, M, 3).

    Los resúmenes usan tasas físicas. La primera fila corresponde al ensamble
    inicial condicionado al primer dato. days mide tiempo transcurrido en días.

    With return_forecast=True, append the (T, 3) forecast mean for I, R, D.
    The first row is NaN because initialization already uses the first observation.
    """
    observations = np.asarray(observations, dtype=float)
    if observations.ndim != 2 or observations.shape[1] != 3 or len(observations) == 0:
        raise ValueError("Se espera una serie no vacía de observaciones (T, 3).")
    if not np.isfinite(observations).all() or (observations < 0).any() or (observations.sum(axis=1) > population).any():
        raise ValueError("Las observaciones deben ser finitas y compatibles con la población.")
    if days is None:
        days = np.arange(len(observations))
    days = np.asarray(days, dtype=float)
    if days.shape != (len(observations),) or not np.isfinite(days).all() or (np.diff(days) <= 0).any():
        raise ValueError("Se requiere un día por observación, en orden estrictamente creciente.")

    rng = np.random.default_rng(seed)
    ensemble = initialize_sird_ensemble(observations[0], population, ensemble_size,
                                        state_std, initial_rates, initial_log_std, rng)
    mean = np.empty((len(days), 7))
    std = np.empty_like(mean)
    parameters = np.empty((len(days), ensemble_size, 3))
    forecast_mean = np.full((len(days), 3), np.nan) if return_forecast else None
    for k, day in enumerate(days):
        # El primer dato ya intervino en la inicialización.
        if k > 0:
            forecast = predict_sird_ensemble(ensemble, population, days[k - 1], day,
                                             rate_step, rng)
            if forecast_mean is not None:
                forecast_mean[k] = forecast[:, 1:4].mean(axis=0)
            ensemble = analyze_sird_ensemble(forecast, observations[k], population,
                                             observation_cov, rng)
        mean[k] = ensemble.mean(axis=0)
        std[k] = ensemble.std(axis=0, ddof=1)
        parameters[k] = ensemble[:, 4:]
    if forecast_mean is not None:
        return mean, std, parameters, forecast_mean
    return mean, std, parameters
