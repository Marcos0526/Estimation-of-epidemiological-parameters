# Estimation of epidemiological parameters

Proyecto para explorar datos epidemiológicos de COVID-19 en México y estudiar
la evolución de una epidemia mediante modelos compartimentales SIRD y un filtro
de Kalman. El código de análisis está organizado principalmente en notebooks de
Jupyter y en las funciones reutilizables de `Codigos/codes.py`.

## Estructura del proyecto

```text
.
├── Codigos/
│   ├── codes.py
│   ├── Exploracion.ipynb
│   ├── Lectura.ipynb
│   ├── first_simulation.ipynb
│   ├── get_filtered_by_date.ipynb
│   ├── plots.ipynb
│   └── random.ipynb
├── Datos/
│   ├── 240708 Catalogos.xlsx
│   ├── 240708 Descriptores_.xlsx
│   └── Actualizaciones en la presentación de información
│       referente a casos de COVID.pdf
├── Seminar/
│   ├── example_kalman_filter.ipynb
│   ├── example_kalman_filter copy(español).ipynb
│   ├── kalman_1_datos.png
│   └── kalman_2_resultados.png
├── .gitattributes
├── .gitignore
├── LICENSE
└── README.md
```

### Directorios y archivos principales

- **`Codigos/`**: código de preparación, exploración, visualización y
  modelado. `codes.py` concentra las funciones reutilizables.
- **`Codigos/Lectura.ipynb`**: lectura de los archivos de casos, comparación de
  columnas, selección de variables y concatenación de datos.
- **`Codigos/Exploracion.ipynb`**: revisión de tipos, valores nulos,
  duplicados, estadísticas descriptivas, casos positivos y defunciones.
- **`Codigos/get_filtered_by_date.ipynb`**: filtrado de casos por un intervalo
  de fechas.
- **`Codigos/plots.ipynb`**: generación de series acumuladas de casos y
  defunciones.
- **`Codigos/first_simulation.ipynb`**: primera simulación del modelo SIRD y
  construcción de observaciones para el filtro de Kalman.
- **`Codigos/random.ipynb`**: notebook auxiliar para pruebas o experimentos.
- **`Datos/`**: catálogos, descriptores y documentación de los datos de
  vigilancia epidemiológica.
- **`Seminar/`**: ejemplos didácticos del filtro de Kalman, en inglés y
  español, junto con imágenes de resultados.

Los notebooks hacen referencia a archivos CSV de casos de COVID-19 y a datos
filtrados que no forman parte del árbol versionado actual. El patrón
`Datos/*.csv` está excluido mediante `.gitignore`; por ello, esos archivos
deben generarse o incorporarse localmente antes de ejecutar los notebooks que
los utilizan.

## Funciones de `Codigos/codes.py`

### Preparación y filtrado de datos

#### `dates_to_datetime(df, columns)`

Convierte `FECHA_DEF` y las columnas indicadas en `columns` a fechas de
`pandas`. Para `FECHA_DEF` usa `errors='coerce'`, por lo que valores inválidos
se convierten en `NaT`; las demás columnas se interpretan con el formato
`%Y-%m-%d`. Devuelve el mismo `DataFrame` con las columnas convertidas.

#### `covid_positive(df)`

Filtra el `DataFrame` para conservar únicamente los registros cuya
`CLASIFICACION_FINAL` sea `1`, `2` o `3`, categorías consideradas casos
positivos. Devuelve un nuevo subconjunto del `DataFrame`.

#### `filter_by_date(df, first_date, last_date, column_name='FECHA_INGRESO', positive=True, saved=False)`

Filtra los registros entre `first_date` y `last_date`, con ambos extremos
incluidos. De forma predeterminada conserva solo los casos positivos y utiliza
`FECHA_INGRESO` como fecha de referencia. Si `saved=True`, guarda el resultado
como `filtered_data_<fecha-inicial>_<fecha-final>.csv` en el directorio de
trabajo actual. Devuelve el `DataFrame` filtrado.

### Visualización y series acumuladas

#### `formatear_fechas(ax, intervalo_meses=2)`

Configura un eje de Matplotlib para mostrar fechas por mes, con el intervalo
indicado entre marcas principales. También agrega marcas secundarias, gira las
etiquetas 45 grados y activa una cuadrícula.

#### `plot_log(block, date_1, date_f, intervalo_meses=1)`

Calcula y grafica en escala logarítmica los casos positivos acumulados entre
`date_1` y `date_f`. Completa con cero los días sin registros. Imprime el total
acumulado y devuelve:

```python
C, idx, acumulado
```

`C` es la serie acumulada, `idx` el índice diario y `acumulado` un
`DataFrame` con las columnas `Fecha` y `Casos positivos acumulados`.

#### `plot_cumulated(block, date_1, date_f, intervalo_meses=1)`

Calcula y grafica en paneles separados los casos positivos acumulados y las
defunciones acumuladas. Devuelve:

```python
C, D, idx
```

donde `C` representa casos, `D` defunciones e `idx` el rango diario.

### Modelo SIRD

#### `SIRD(t, y, beta, gamma, delta, N)`

Define las ecuaciones diferenciales del modelo compartimental SIRD. El estado
`y` contiene `[S, I, R, D]`: susceptibles, infectados, recuperados y
defunciones. `beta` controla la transmisión, `gamma` la recuperación,
`delta` la mortalidad y `N` el tamaño de la población. Devuelve las derivadas
`(dS, dI, dR, dD)`.

#### `SIRD_simulation(beta, gamma, delta, N, y0, t_max, dt=1.0, verbose=True)`

Integra numéricamente el modelo SIRD con `scipy.integrate.solve_ivp` y el
método `RK45`. Genera tiempos desde `0` hasta `t_max` con paso `dt`, valida
que el integrador haya terminado correctamente y devuelve:

```python
t, y, df
```

`t` contiene los tiempos, `y` es una matriz con las columnas `[S, I, R, D]`
y `df` es el mismo resultado en un `DataFrame` con la columna `dia`. Si
`verbose=True`, imprime `R0`, el pico de infectados, las muertes finales y el
error máximo de conservación de la población.

### Filtro de Kalman

#### `SIRD_lineal(beta, gamma, delta, s=1.0, dt=1.0)`

Construye la matriz continua `A` del modelo SIRD linealizado para el estado
`[I, R, D]`, donde `s` representa la fracción susceptible `S/N`. Calcula la
matriz de transición discreta `F` mediante la exponencial matricial
`expm(A * dt)`. Devuelve `(A, F)`.

#### `kalman_filter(z, F, H, Q, R, x0, P0)`

Implementa un filtro de Kalman discreto básico. En cada observación realiza
una predicción y una actualización, almacenando el estado y la covarianza
filtrados. Devuelve:

```python
x_filt, P_filt
```

`x_filt` tiene una fila por observación y `P_filt` contiene la matriz de
covarianza correspondiente a cada paso.

#### `get_Y_obs(df, date_1, date_f, N=None, gamma_time=14, date_col='FECHA_INGRESO')`

Construye las observaciones diarias `[I, R, D]` a partir de los registros
positivos:

- `C`: casos acumulados según `date_col`.
- `D`: defunciones acumuladas según `FECHA_DEF`.
- `R`: recuperados estimados desplazando `gamma_time` días los sobrevivientes.
- `I`: infectados activos calculados como `C - R - D`.

Devuelve `Y` como un arreglo NumPy de dimensiones `(T, 3)` y `df_obs`, un
`DataFrame` indexado por fecha con las columnas `I`, `R`, `D`, `C` y,
opcionalmente, `S` cuando se proporciona `N`.

## Dependencias

El módulo utiliza:

- `pandas`
- `numpy`
- `matplotlib`
- `scipy`


El flujo habitual es cargar los datos, convertir las fechas, filtrar los casos
positivos, generar las series observadas y finalmente ejecutar la simulación
SIRD o el filtro de Kalman.