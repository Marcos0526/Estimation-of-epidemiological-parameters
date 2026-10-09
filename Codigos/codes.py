import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates  
from scipy.integrate import solve_ivp

def dates_to_datetime(df, columns): #tranform date columns to datetime format
    df['FECHA_DEF'] = pd.to_datetime(df['FECHA_DEF'], errors='coerce')
    for col in columns:
        df[col] = pd.to_datetime(df[col], format='%Y-%m-%d')
    return df

def covid_positive(df): # filter dataframe to include only positive COVID cases
    return df[(df['CLASIFICACION_FINAL'] == 1) | (df['CLASIFICACION_FINAL'] == 2) | (df['CLASIFICACION_FINAL'] == 3)]

def filter_by_date(df, first_date, last_date, column_name='FECHA_INGRESO', saved=False): # filter dataframe by date range
    filtro = (df[column_name] >= first_date) & (df[column_name] <= last_date)
    out = df.loc[filtro]
    if saved:
        out.to_csv(f'filtered_data_{first_date}_{last_date}.csv', index=False)
    return out

def formatear_fechas(ax, intervalo_meses=2): #date formatting for plots
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=intervalo_meses))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))   # ej. "mar 2020"
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.tick_params(axis="x", rotation=45)
    ax.grid(alpha=0.3)

def plot_log(block, date_1, date_f, intervalo_meses=1): #plot cumulative positive cases in log scale

    block = dates_to_datetime(block, ['FECHA_INGRESO'])

    pos = covid_positive(block)
    idx = pd.date_range(date_1, date_f, freq='D')

    new = pos.groupby('FECHA_INGRESO').size().reindex(idx, fill_value=0)

    C = new.cumsum()      # casos acumulados (= I+R+D)

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


def plot_cumulated(block, date_1, date_f, intervalo_meses=1):# plot cumulative positive cases and deaths

        block = dates_to_datetime(block, ['FECHA_INGRESO', 'FECHA_SINTOMAS', 'FECHA_DEF'])

        pos = covid_positive(block)
        idx = pd.date_range(date_1, date_f, freq='D')

        new    = pos.groupby('FECHA_INGRESO').size().reindex(idx, fill_value=0)
        deaths = pos['FECHA_DEF'].dropna().value_counts().reindex(idx, fill_value=0)

        C = new.cumsum()      # casos acumulados (= I+R+D)
        D = deaths.cumsum()   # defunciones acumuladas

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

def SIRD(t, y, beta, gamma, delta, N): #Define the SIRD model differential equations
    S, I, R, D = y
    infecciones = beta * S * I / N
    dS = -infecciones
    dI = infecciones - gamma * I - delta * I
    dR = gamma * I
    dD = delta * I

    return dS, dI, dR, dD


def SIRD_simulation(beta, gamma, delta, N, y0, t_max, dt=1.0, verbose=True):
    """
    Simula el modelo SIRD con parámetros fijos.

    Devuelve
    --------
    t  : np.ndarray, shape (n,)    -> días
    y  : np.ndarray, shape (n, 4)  -> columnas [S, I, R, D]
    df : pd.DataFrame              -> mismos datos, con columnas con nombre
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