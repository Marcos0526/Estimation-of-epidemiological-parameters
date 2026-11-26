# Estimation of epidemiological parameters

Proyecto para explorar datos epidemiológicos de COVID-19 en México y estudiar
la evolución de una epidemia. La tesis estudia el modelo compartimental SIRD
mediante el filtro de Kalman por ensamble (EnKF). El código de análisis está
organizado en notebooks de Jupyter y en las funciones reutilizables de
`Codigos/codes.py`.

Consulta el [contexto de investigación](markdowns/resume.md), la [nota de
implementación](markdowns/Codigos/first_estmation.md) y el [recorrido de la
primera estimación](markdowns/Codigos/recorrido_primera_estimacion.md). Las
instrucciones locales para agentes están en `markdowns/AGENTS.md`.

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
├── kalman_base/
├── markdowns/
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
- **`kalman_base/`** es el archivo local excluido de Git. Contiene los cuadernos del
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
Con `return_forecast=True`, se añade un cuarto arreglo de forma `(T, 3)`
con las medias pronosticadas de `I, R, D` antes de asimilar cada observación.
La primera fila contiene `NaN` porque la inicialización ya usa ese dato.

Dos funciones internas mantienen las restricciones del modelo.
`_rate_triplet` valida la longitud, la finitud y el signo de las ternas de
configuración. `_balance_compartments` corrige los recuentos negativos de
`I, R, D`, los reduce proporcionalmente si su suma supera `N` y recalcula
`S = N - I - R - D` en el propio ensamble.

### Funciones de los cuadernos de estimación

- En la estimación de 1167 fechas, `mark_periods(ax)` marca el comienzo de las olas y sombrea los intervalos fuera de ellas. `plot_band(ax, column)` dibuja la media del componente elegido con una banda de una desviación estándar entre miembros. [Cuaderno de estimación](Codigos/first_estimation.ipynb)
- En el recorrido de 30 fechas con 50 miembros, `analysis_step(forecast, day)` calcula la corrección diaria y devuelve el ensamble actualizado. Durante los primeros diez días, `record(day, operation, values, components, *, shared=False)` guarda los valores de cada operación para construir el DataFrame `operation_log`. `member_trace(member)` recupera las filas del miembro elegido junto con las magnitudes compartidas, como la observación y la ganancia. Por ejemplo, `member_trace(1).query("day == 1")` muestra cómo se actualizó el primer miembro el 29 de febrero de 2020. [Cuaderno paso a paso](Codigos/enkf_30_day_walkthrough.ipynb)

## Barrido del tamaño del ensamble

### ¿Qué cambia al usar más miembros?

Cada miembro del ensamble representa una trayectoria posible del modelo SIRD.
El filtro usa las diferencias entre miembros para estimar cómo se relacionan
los compartimientos con las tasas de contagio, recuperación y fallecimiento.
Esas relaciones determinan cuánto se corrige cada tasa al incorporar una
observación. Aumentar el número de miembros cambia la estimación de esas
relaciones y exige integrar más trayectorias en cada paso diario.

Para medir ese efecto, probamos **20 tamaños**, de 10 a 200 en pasos de 10, con
**10 semillas por tamaño** (2026 a 2035) y **dos escenarios**. Cada combinación
de escenario, tamaño y semilla es una corrida. En total se programaron
**400 corridas**. Se tomó como referencia un ensamble de 160 miembros porque
es el tamaño de la primera estimación. Las semillas permiten repetir los
sorteos del filtro para una configuración dada y medir cuánto varía el
resultado entre corridas del mismo tamaño. [Primera estimación](Codigos/first_estimation.ipynb) · [Código del barrido](Codigos/enkf_sweep.py)

| Escenario | ¿Qué recibe el filtro? | ¿Qué se puede evaluar? |
| --- | --- | --- |
| México | Series diarias del 28 de febrero de 2020 al 9 de mayo de 2023. Los casos activos se reconstruyen como $I=C-R-D$, donde $C$ son los casos acumulados. Se imputan recuperaciones a los 14 días para casos sin defunción registrada y se acumulan las defunciones por fecha registrada. | El error respecto a las series reconstruidas, en todo el periodo y por fase. La interpretación de la tasa de recuperación depende de la imputación. |
| Simulación | Una trayectoria SIRD de 240 días con tasas conocidas que varían suavemente. Se añade ruido a $I$, $R$ y $D$ con una semilla fija. | La cercanía de las tasas estimadas a las tasas conocidas en ese escenario. |

En cada escenario, todas las corridas reciben la misma serie de observaciones.
Se conservan la inicialización, la covarianza de observación y la regla de
evolución de las tasas. En México el filtro continúa entre olas y después de
la sexta ola. Las fases sirven para agrupar métricas.
[Preparación de los dos escenarios](Codigos/enkf_sweep_data.py)

### Cómo leer el error

En el día $t$, el filtro avanza los miembros del ensamble desde el análisis
del día anterior. La media de esos miembros constituye el **pronóstico a un
día**. Después incorpora las observaciones del día $t$ y obtiene el **análisis**.
Por ejemplo, en `daily.csv` las columnas `forecast_I`, `observed_I` y
`analysis_I` muestran, respectivamente, el pronóstico, el dato reconstruido
y la estimación corregida de casos activos. El día 0 se usa para inicializar
el filtro y queda fuera del cálculo de errores.

En la corrida de México con 160 miembros y semilla 2026, el 22 de marzo de
2023 se pronosticaron **642 309** casos activos frente a **35 345** reconstruidos.
Tras asimilar el dato, el análisis quedó en **35 345** al redondear a personas.
Ese día ilustra por qué conviene examinar por separado el ajuste después de
la corrección y el error del pronóstico. El ajuste de los compartimientos por
sí solo es insuficiente para determinar si las tasas estimadas se pueden
distinguir entre sí.

| Indicador | Cómo interpretarlo |
| --- | --- |
| MAE y RMSE | Ambos miden el error en personas. El MAE promedia los errores absolutos. El RMSE eleva al cuadrado cada diferencia antes de promediar y aplica una raíz al final, por lo que responde más a los errores grandes. |
| Error de defunciones diarias | El flujo pronosticado en el día $t$ se obtiene restando el acumulado del análisis de $t-1$ al acumulado pronosticado de $t$. Se compara con el incremento del acumulado observado entre esas dos fechas. |
| RMSE de las tasas en escala logarítmica | En la simulación, compara el logaritmo de la media estimada de cada tasa con el logaritmo de su valor conocido. |
| Cobertura de las tasas | En la simulación, indica la fracción de días en que la tasa conocida quedó entre los percentiles 5 y 95 de los miembros. |

Cada indicador se calcula primero para una corrida y después se resume por
tamaño con la mediana y los cuartiles de las corridas exitosas. El rango
intercuartílico describe la variación **entre semillas**. Los percentiles 5 y
95 por fecha describen la variación **entre miembros de una corrida**.

### Qué observamos en el barrido

- **Estabilidad numérica.** Terminaron 395 de las 400 corridas. En México
  fallaron tres de diez con 10 miembros y dos de diez con 20 miembros.
  El resto de los tamaños completó las diez réplicas en ambos escenarios.
  El conteo de fallos se conserva junto con las métricas de las corridas exitosas.
- **Datos de México.** El menor RMSE mediano de pronóstico de casos activos
  se obtuvo con 30 miembros (15 990 personas). Con 160 fue de 18 372. Para
  ese tamaño, del 22 al 24 de marzo de 2023 se concentró aproximadamente el
  66 % del error cuadrático total de casos activos, como mediana entre
  semillas. El menor RMSE mediano de defunciones diarias se obtuvo con
  140 miembros (32,3 defunciones).
- **Simulación.** El RMSE mediano de pronóstico de casos activos bajó de
  1 934 a 593 personas entre 10 y 200 miembros. Los mínimos del RMSE de las
  tasas conocidas se dieron con tamaños distintos para $\beta$, $\gamma$ y
  $\delta$. La comparación depende de la variable de interés.

La corrida de México con 160 miembros y semilla 2026 coincidió con la
primera estimación del cuaderno. Las seis corridas del piloto reprodujeron
sus series y métricas en el barrido completo. Los resultados descritos
corresponden al barrido local `5bda922378a240eb90eeed7cbe943e4e`.

### Cómo usar MLflow con este estudio

El script registra los experimentos automáticamente al ejecutarse. MLflow
guarda sus parámetros, métricas y estados en `Datos/enkf_sweep/mlflow.db`,
y sus archivos en `Datos/enkf_sweep/artifacts/`. El servidor local abre
la interfaz para consultar esos resultados, sin volver a ejecutar el filtro.

#### Activar el entorno y abrir la interfaz

Desde la raíz del repositorio, activa el entorno que usaste para el barrido.
Si todavía no tienes uno, créalo con un Python 3 reciente e instala las
dependencias. En macOS y Linux puedes usar estos comandos la primera vez.

```sh
python3 -m venv "$HOME/.venvs/enkf-sweep"
source "$HOME/.venvs/enkf-sweep/bin/activate"
python -m pip install -r Codigos/requirements-enkf-sweep.txt
```

Con el entorno activo, comprueba que existe la base del estudio e inicia
MLflow desde la raíz del repositorio.

```sh
ls -lh Datos/enkf_sweep/mlflow.db
mlflow server --host 127.0.0.1 --port 8081 \
  --backend-store-uri "sqlite:///$PWD/Datos/enkf_sweep/mlflow.db"
```

Deja la terminal abierta y selecciona el experimento **SIRD EnKF ensemble
size** en el navegador. Puedes detener el servidor con `Ctrl+C`. Si el
puerto 8081 está ocupado, cambia el número tanto en el comando como en la
dirección del navegador. Si la base todavía no existe, ejecuta primero el
piloto descrito en la sección de reproducción. [Abrir MLflow](http://127.0.0.1:8081)

#### Qué incluye el registro

| Dónde mirar | Qué se puede consultar |
| --- | --- |
| Corrida principal `sweep` | `configuration.json` identifica el código, las versiones, las semillas y las huellas SHA-256. `inputs/` contiene los datos diarios agregados y las tasas verdaderas de la simulación. `report/` contiene los resúmenes y las gráficas del barrido. |
| Corrida individual | **Parameters** muestra el escenario, tamaño, semilla y configuración del filtro. **Metrics** muestra valores globales de error, cobertura en la simulación y tiempo de cálculo. El estado indica si terminó o falló. |
| **Artifacts** de una corrida | `daily.csv` muestra las observaciones, el pronóstico y el análisis por día. `metrics.csv` contiene el error global (`Overall`) y por fase, como `Sixth Wave`. En una corrida fallida, `error.txt` contiene la traza. |

#### Cómo medir y comparar resultados

En la tabla del experimento, pega este filtro para ver las corridas
exitosas de México con semilla 2026 dentro del barrido completo.

```text
tags.mlflow.parentRunId = '5bda922378a240eb90eeed7cbe943e4e' AND params.scenario = 'mexico' AND params.seed = '2026' AND attributes.status = 'FINISHED'
```

Selecciona, por ejemplo, las filas de 30, 160 y 200 miembros y usa **Compare**
para revisar sus valores lado a lado. Para estudiar el efecto del tamaño
entre semillas, quita la condición sobre `params.seed` y consulta
`report/summary.csv` en la corrida principal. Ese archivo contiene la
mediana y los cuartiles por tamaño. [Guía de búsqueda de corridas de MLflow](https://mlflow.org/docs/latest/ml/search/search-runs/)

| Métrica en MLflow | Qué mide |
| --- | --- |
| `forecast_I_rmse` | Error de pronóstico a un día de casos activos reconstruidos, en personas. |
| `analysis_I_rmse` | Error de casos activos después de incorporar la observación, en personas. |
| `forecast_daily_deaths_rmse` | Error de pronóstico del número diario de defunciones, en personas. |
| `runtime_seconds` | Segundos de cálculo del filtro y de sus métricas para esa corrida. |
| `beta_log_rmse`, `gamma_log_rmse`, `delta_log_rmse` | Error de cada tasa en escala logarítmica frente a su valor conocido en la simulación. |
| `beta_p05_p95_coverage` y sus equivalentes para `gamma` y `delta` | Fracción de días de la simulación en que la tasa conocida quedó dentro de la banda del ensamble. Por ejemplo, 0,80 equivale al 80 % de los días. |

Para verificar el ejemplo del 22 de marzo de 2023, abre la corrida de
160 miembros y semilla 2026, descarga `daily.csv` desde **Artifacts** y busca
la fecha en la columna `date`. Compara `forecast_I`, `observed_I` y
`analysis_I` en esa fila. Las métricas de la pantalla resumen la corrida
completa, y los CSV permiten examinar días y fases concretos.
[Abrir la corrida del ejemplo](http://127.0.0.1:8081/#/experiments/1/runs/6462e2bf13c1497abd7530c9fd35d98c)

Para consultar los resultados desde el sistema de archivos, abre
`Datos/enkf_sweep/reports/<id-del-barrido>/`.

| Archivo | Para qué sirve |
| --- | --- |
| `runs.csv` | Localizar una corrida por escenario, tamaño y semilla, y consultar su estado y métricas globales. |
| `summary.csv` | Comparar medianas, cuartiles y diferencias respecto a 160 miembros. El campo `runtime_seconds` mide el cálculo por corrida y está condicionado por la carga de CPU. |
| `counts.csv` | Comprobar cuántas semillas terminaron o fallaron en cada escenario y tamaño. |
| `rate_spread.csv` y figuras `*_rate_spread.png` | Examinar cómo cambia con el tiempo la dispersión de las tasas estimadas entre semillas. |
| Figuras `*_metrics.png` | Comparar error y tiempo de cómputo entre tamaños. La línea representa la mediana y la banda, el rango intercuartílico. |

La base, los artefactos y los informes se guardan en `Datos/`, que está
excluido de Git. Para consultar las corridas desde otra copia del repositorio
se necesita acceso al directorio local de resultados o ejecutar un nuevo
barrido.

### Cómo repetir el experimento

Desde la raíz del repositorio, activa un entorno virtual e instala las
dependencias. El piloto evalúa 10, 160 y 200 miembros con la semilla 2026.
La opción `--scenario synthetic` permite probar el flujo sin el archivo
nacional.

```sh
python -m pip install -r Codigos/requirements-enkf-sweep.txt
python Codigos/enkf_sweep.py --pilot --scenario synthetic
```

Para incluir México, coloca los archivos
`Datos/raw_data/COVID19MEXICO2020_2021_2022_2023.csv` y
`Datos/metadata/covid_waves_parameters.csv` en sus rutas correspondientes.
Los siguientes comandos ejecutan, respectivamente, el piloto de seis
corridas, una sola configuración y el barrido completo.

```sh
python Codigos/enkf_sweep.py --pilot
python Codigos/enkf_sweep.py --scenario mexico --sizes 160 --seeds 2026
python Codigos/enkf_sweep.py --workers 4
```

Las pruebas de regresión se ejecutan con
`python -m unittest discover -s Codigos -p test_enkf_sweep.py`.

#####
El flujo habitual es cargar los datos, convertir las fechas, filtrar los casos
positivos, generar las series observadas y finalmente ejecutar la simulación
SIRD.
