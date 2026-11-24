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
