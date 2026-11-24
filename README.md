# Estimation of epidemiological parameters

Proyecto para explorar datos epidemiológicos de COVID-19 en México y estudiar
la evolución de una epidemia. La tesis estudia el modelo compartimental SIRD
mediante el filtro de Kalman por ensamble (EnKF). El código de análisis está
organizado en notebooks de Jupyter y en las funciones reutilizables de
`Codigos/codes.py`.

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
├── kalman base/
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
  construcción de observaciones SIRD.
- **`Codigos/random.ipynb`**: notebook auxiliar para pruebas o experimentos.
- **`Datos/`**: catálogos, descriptores y documentación de los datos de
  vigilancia epidemiológica.
- **`kalman base/`** es el archivo local excluido de Git. Contiene los cuadernos del
  filtro clásico, su implementación en `classical_kalman.py` y los ejemplos
  didácticos con sus figuras en `Seminar/`.

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

### Observaciones SIRD

#### `get_Y_obs(df, date_1, date_f, N=None, gamma_intervals=None, gamma_time=14, date_col='FECHA_INGRESO')`

Construye las observaciones diarias `[I, R, D]` a partir de los registros
positivos:

- `C`: casos acumulados según `date_col`.
- `D`: defunciones acumuladas según `FECHA_DEF`.
- `R`: recuperados estimados al desplazar las fechas de ingreso de los
  sobrevivientes. El retraso es `gamma_time` días o el inverso redondeado de
  la tasa de recuperación correspondiente si se proporciona `gamma_intervals`.
- `I`: infectados activos calculados como `C - R - D`.

Devuelve `Y` como un arreglo NumPy de dimensiones `(T, 3)` y `df_obs`, un
`DataFrame` indexado por fecha con las columnas `I`, `R`, `D`, `C` y,
opcionalmente, `S` cuando se proporciona `N`.

### Estimación de estados y tasas con el EnKF

La función principal, `run_sird_enkf`, recibe las observaciones de
`get_Y_obs` en el orden `[I, R, D]`. Cada uno de los `M` miembros del
ensamble representa una trayectoria de `[S, I, R, D, beta, gamma, delta]`.
Los cuatro primeros valores se expresan en personas y las tres tasas en
inversos de día.

| Función | Qué hace y qué devuelve |
| --- | --- |
| `initialize_sird_ensemble` | Usa la primera observación como centro de los sorteos de `I, R, D`, con desviaciones `state_std`. Calcula `S` por balance y sortea tasas positivas con medianas `initial_rates` y desviaciones logarítmicas `initial_log_std`. Devuelve el ensamble inicial de forma `(M, 7)`. |
| `gaussian_rate_walk` | Avanza las tasas de cada miembro con incrementos gaussianos independientes en escala logarítmica. Multiplica `log_step_std` por la raíz de la duración del paso en días. Recibe y devuelve una matriz `(M, 3)` de tasas en unidades físicas. |
| `predict_sird_ensemble` | Llama a `rate_step` para obtener las tasas del intervalo. Con ellas fijas, integra `SIRD` mediante `solve_ivp` desde los compartimientos previos de cada miembro. Devuelve los cuatro compartimientos pronosticados y las tres tasas de cada miembro. |
| `analyze_sird_ensemble` | Usa `SIRD_OBSERVATION_MATRIX` para extraer `I, R, D` del pronóstico. Calcula la ganancia de Kalman con las covarianzas del ensamble y la covarianza de observación `observation_cov` de forma `(3, 3)`. Para cada miembro, compara su pronóstico con la observación perturbada, actualiza compartimientos y logaritmos de las tasas, y devuelve el ensamble corregido con tasas positivas. |
| `run_sird_enkf` | Inicializa el ensamble con el primer dato y alterna predicción y análisis desde la segunda fecha. `ensemble_size` fija el número de miembros, `seed` controla los sorteos y `days` permite indicar los tiempos cuando los pasos difieren de un día. |

Para `T` fechas, `run_sird_enkf` devuelve `mean, std, parameters`. Los dos
primeros arreglos tienen forma `(T, 7)` y contienen la media y la desviación
estándar muestral entre miembros. El tercero tiene forma `(T, M, 3)`. Por
ejemplo, `parameters[:, 0, 1]` contiene la trayectoria de recuperación
$\gamma(t)$ del primer miembro. Para cambiar la evolución de las tasas, se
pasa otra función con la firma `rate_step(rates, start_day, end_day, rng)`.

Dos funciones internas mantienen las restricciones del modelo.
`_rate_triplet` valida la longitud, la finitud y el signo de las ternas de
configuración. `_balance_compartments` corrige los recuentos negativos de
`I, R, D`, los reduce proporcionalmente si su suma supera `N` y recalcula
`S = N - I - R - D` en el propio ensamble.

### Funciones de los cuadernos de estimación

- En la estimación de 1167 fechas, `mark_periods(ax)` marca el comienzo de las olas y sombrea los intervalos fuera de ellas. `plot_band(ax, column)` dibuja la media del componente elegido con una banda de una desviación estándar entre miembros. [Cuaderno de estimación](Codigos/first_estimation.ipynb)
- En el recorrido de 30 fechas con 50 miembros, `analysis_step(forecast, day)` calcula la corrección diaria y devuelve el ensamble actualizado. Durante los primeros diez días, `record(day, operation, values, components, *, shared=False)` guarda los valores de cada operación para construir el DataFrame `operation_log`. `member_trace(member)` recupera las filas del miembro elegido junto con las magnitudes compartidas, como la observación y la ganancia. Por ejemplo, `member_trace(1).query("day == 1")` muestra cómo se actualizó el primer miembro el 29 de febrero de 2020. [Cuaderno paso a paso](Codigos/enkf_30_day_walkthrough.ipynb)

## Dependencias

El módulo utiliza:

- `pandas`
- `numpy`
- `matplotlib`
- `scipy`


El flujo habitual es cargar los datos, convertir las fechas, filtrar los casos
positivos, generar las series observadas y finalmente ejecutar la simulación
SIRD.
