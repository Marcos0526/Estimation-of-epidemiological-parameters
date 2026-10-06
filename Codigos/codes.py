import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates  

def covid_positive(df):
    return df[(df['CLASIFICACION_FINAL'] == 1) | (df['CLASIFICACION_FINAL'] == 2) | (df['CLASIFICACION_FINAL'] == 3)]

def filter_by_date(df, first_date, last_date, column_name='FECHA_INGRESO', saved=False):
    filtro = (df[column_name] >= first_date) & (df[column_name] <= last_date)
    out = df.loc[filtro]
    if saved:
        out.to_csv(f'filtered_data_{first_date}_{last_date}.csv', index=False)
    return out

def formatear_fechas(ax, intervalo_meses=2):
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=intervalo_meses))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))   # ej. "mar 2020"
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.tick_params(axis="x", rotation=45)
    ax.grid(alpha=0.3)

def plot_log(block, date_1, date_f, intervalo_meses=1):

    block['FECHA_INGRESO'] = pd.to_datetime(block['FECHA_INGRESO'], format='%Y-%m-%d')
    block['FECHA_SINTOMAS'] = pd.to_datetime(block['FECHA_SINTOMAS'], format='%Y-%m-%d')
    block['FECHA_DEF'] = pd.to_datetime(block['FECHA_DEF'], errors='coerce')
    block['FECHA_DEF'] = pd.to_datetime(block['FECHA_DEF'], format='%Y-%m-%d')


    pos = covid_positive(block)
    idx = pd.date_range(date_1, date_f, freq='D')

    new    = pos.groupby('FECHA_INGRESO').size().reindex(idx, fill_value=0)
    deaths = pos['FECHA_DEF'].dropna().value_counts().reindex(idx, fill_value=0)

    C = new.cumsum()      # casos acumulados (= I+R+D)
    D = deaths.cumsum()   # defunciones acumuladas

    print('Casos positivos acumulados: ', C.iloc[-1])
    print('Defunciones acumuladas: ', D.iloc[-1])

    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    ax[0].semilogy(C.where(C > 0))      # evita log(0) al inicio
    ax[0].set_title('Casos acumulados (log)')
    ax[1].plot(D)
    ax[1].set_title('Defunciones acumuladas')
    for a in ax:
        formatear_fechas(a, intervalo_meses)
    plt.tight_layout()
    plt.show()

    return C, D,idx



