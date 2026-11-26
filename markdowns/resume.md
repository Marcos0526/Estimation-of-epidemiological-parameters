# Contexto de investigación de la tesis SIRD–EnKF

Estado de referencia del repositorio al 25 de noviembre de 2026.

## Cómo retomar el proyecto

Al iniciar una sesión, lee [AGENTS.md](AGENTS.md) para conocer las instrucciones de trabajo y este archivo para reconstruir el contexto científico. Revisa el estado de Git y abre la fuente pertinente a la solicitud. El [README](../README.md) documenta el código y la reproducción del barrido. El [recorrido de la primera estimación](Codigos/recorrido_primera_estimacion.md) explica la actualización numérica de un miembro y la comparación por fases. [La nota de implementación](Codigos/first_estmation.md) detalla el algoritmo y las salidas. Consulta las celdas fuente de los cuadernos, los archivos de Python o los `.tex` para resolver discrepancias. Usa los resultados guardados para comprobar cifras y procedencia.

La estructura actual sitúa el código versionado en `Codigos/`. Los materiales del Kalman clásico están en `kalman_base/`. Los datos, borradores, materiales del Kalman clásico y animaciones se conservan localmente según `.gitignore`. En otra copia, comprueba qué archivos están disponibles antes de citar rutas o intentar reproducir resultados.

## Pregunta, alcance y aportación

La tesis de maestría estudia cómo el **filtro de Kalman por ensamble (EnKF)** permite estimar conjuntamente los compartimientos SIRD y la evolución de las tasas efectivas de contagio $\beta(t)$, recuperación $\gamma(t)$ y fallecimiento por COVID-19 $\delta(t)$ a partir de los registros de México. La pregunta que organiza el análisis es **en qué medida los datos permiten distinguir los cambios de una tasa de los de las demás**.

El análisis es agregado, nacional y retrospectivo. La corrida nacional cubre del **28 de febrero de 2020 al 9 de mayo de 2023**, ambas fechas incluidas. El SIRD supone mezcla homogénea y población cerrada y constante. Quedan fuera del modelo explícito la migración, los nacimientos, las muertes por otras causas y la heterogeneidad por región, edad o variante. El alcance excluye la atribución causal de cambios a medidas sanitarias o variantes y el pronóstico fuera del periodo de estudio.

El trabajo construye series a partir de registros de vigilancia, estima secuencialmente los estados y las tasas, y evalúa qué información sostiene esas estimaciones. El ajuste de los compartimientos, el error de pronóstico a un día y la recuperación de tasas conocidas responden a preguntas distintas.

## Modelo y cantidades observadas

El estado compartimental es $(S,I,R,D)$, en personas, con $S+I+R+D=N$. En la corrida nacional $N=126\,014\,024$. Las tres tasas se expresan en $\mathrm{día}^{-1}$. Si $f=\beta SI/N$, $g=\gamma I$ y $m=\delta I$ son los flujos instantáneos en personas por día, entonces

$$
(\dot S,\dot I,\dot R,\dot D)=(-f,\ f-g-m,\ g,\ m).
$$

- $I(t)$ es un compartimiento puntual. Las infecciones nuevas en un intervalo son su flujo de entrada $\int f\,dt$. El registro por ingreso es una serie de eventos registrados cuya fecha puede diferir de la fecha de infección y de la de síntomas.
- $D(t)$ es un acumulado. Las defunciones diarias corresponden al incremento de $D$ durante el intervalo. $\gamma$ y $\delta$ son tasas competitivas de salida de $I$. Cuando son constantes, el tiempo medio en $I$ es $1/(\gamma+\delta)$. La tasa $\delta$ se expresa por persona infectada y por día, mientras que la letalidad es una proporción.
- En el SIRD, $R_t=[\beta(t)/(\gamma(t)+\delta(t))]S(t)/N$ bajo condiciones congeladas en $t$. La columna histórica `R0` de la primera estimación guarda el cociente de las **tasas medias del ensamble** con susceptibilidad completa, $\bar\beta/(\bar\gamma+\bar\delta)$. El factor $S/N$ y el orden de las operaciones importan al comparar cocientes, promedios entre miembros y medianas en el tiempo.
- Para $I>0$, $\dot I/I=\beta S/N-\gamma-\delta$. La serie de activos por sí sola informa sobre esa combinación. Los datos de defunción aportan información sobre $\delta I$. En el canal de recuperación, la información sobre $\gamma I$ está condicionada por la imputación.

## Datos y reconstrucciones

Los insumos locales son los CSV anuales de `Datos/raw_data/`, el consolidado `COVID19MEXICO2020_2021_2022_2023.csv`, los catálogos y descriptores en `Datos/metadata/` y el calendario de fases `Datos/metadata/covid_waves_parameters.csv`. Los CSV en `Datos/waves/` contienen subconjuntos de registros por fecha de ingreso.

**Primera estimación nacional y escenario México del barrido.** Se lee el consolidado en bloques de 250 000 registros, con `FECHA_INGRESO`, `FECHA_DEF` y `CLASIFICACION_FINAL`. Se seleccionan ingresos en el intervalo y positivos con clasificación 1, 2 o 3. `9999-99-99` representa fecha de defunción ausente. Sobre un calendario diario común:

1. Se acumulan ingresos positivos por fecha de ingreso para formar $C$.
2. Se acumulan las defunciones de los casos seleccionados por su fecha registrada para formar $D$.
3. Se imputa recuperación a los 14 días del ingreso para los casos sin defunción registrada y se acumula para formar $R$.
4. Se reconstruyen $I=C-R-D$ y, cuando se solicita, $S=N-C$. Cada fila de entrada del filtro representa $(I,R,D)$. La primera fecha aporta $(2,0,0)$.
5. Los conteos de cada bloque se suman sobre el mismo calendario. Un desenlace posterior a la fase del ingreso conserva su fecha dentro del horizonte. Los intervalos entre fases se mantienen en la estimación.

La hipótesis de recuperación a 14 días produce ceros iniciales por construcción y aporta información indirecta a $\gamma$. Además, $I$ comparte información con $R$ y $D$ por la identidad anterior, aunque la primera configuración use covarianza diagonal. La fecha de defunción ausente es la condición del algoritmo para imputar.

**Reconstrucción exploratoria por fases.** `observations.ipynb` lee los CSV de fases, llega hasta el 31 de marzo de 2023 e imputa $\operatorname{round}(1/\gamma_{\rm medio})$ días según la fase del ingreso. Su calendario rellena los metadatos en abril de 2022 con los valores previos; sus registros fuente se limitan a los CSV por fase. El código aplica `drop_duplicates()` después de seleccionar tres variables y agregar la fase; dos personas con los mismos valores pueden quedar representadas por una sola fila. El producto `Datos/sird_observations.csv` alimenta `observations_plots.ipynb`. La corrida nacional y el barrido reconstruyen desde el consolidado con la regla fija de 14 días.

### Calendario de la corrida nacional

| Fase | Inicio | Fin |
| --- | --- | --- |
| Preconfinamiento | 2020-02-28 | 2020-03-22 |
| Primera ola | 2020-03-23 | 2020-09-30 |
| Segunda ola | 2020-10-01 | 2021-04-30 |
| Tercera ola | 2021-05-01 | 2021-11-30 |
| Cuarta ola | 2021-12-01 | 2022-03-31 |
| Entre olas | 2022-04-01 | 2022-04-30 |
| Quinta ola | 2022-05-01 | 2022-10-31 |
| Sexta ola | 2022-11-01 | 2023-03-31 |
| Después de la sexta ola | 2023-04-01 | 2023-05-09 |

Las fases agrupan resúmenes y métricas. Los 160 miembros continúan durante las **1 167 fechas**, incluidas las 30 de abril de 2022 y las 39 posteriores a la sexta ola. Las etiquetas de variantes y rangos del CSV son metadatos descriptivos del proyecto. Su procedencia y la comparabilidad con cada tasa requieren revisión por separado.

## La primera estimación y el EnKF implementado

La corrida de referencia de [first_estimation.ipynb](../Codigos/first_estimation.ipynb) usa **160 miembros y semilla 2026**. Un generador de NumPy controla la inicialización, la evolución de las tasas y la perturbación de observaciones.

| Elemento | Configuración de referencia |
| --- | --- |
| Compartimientos iniciales | Sorteos normales centrados en el primer $(I,R,D)$, con desviaciones $(1,0,0)$ personas. Se obtiene $S$ por balance. |
| Tasas iniciales | Lognormales independientes con medianas $(0.3508,0.0714,0.0015)$ por día, tomadas de `first_simulation.ipynb`, y desviaciones logarítmicas $(0.35,0.40,0.50)$. |
| Evolución de las tasas | $\log p^f_k=\log p^a_{k-1}+\sigma_p\sqrt{\Delta t}\,Z$, con normales estándar independientes y $\sigma=(0.12,0.12,0.15)$ por raíz de día. |
| Covarianza de observación | $V=\operatorname{diag}(100,16,2.25)$ personas² para $I,R,D$, con desviaciones $(10,4,1.5)$ personas. |

[La implementación](../Codigos/codes.py) sigue este ciclo:

1. `initialize_sird_ensemble` usa el dato del día 0 para construir $M$ vectores $(S,I,R,D,\beta,\gamma,\delta)$. La primera corrección ocurre en el día 1.
2. `predict_sird_ensemble` llama a `rate_step` y mantiene las tasas resultantes fijas durante cada paso. Integra el SIRD con `solve_ivp` (RK45) para cada miembro desde su análisis anterior. La interfaz de `rate_step(rates, start_day, end_day, rng)` recibe y devuelve tasas físicas estrictamente positivas de forma $(M,3)$.
3. `analyze_sird_ensemble` trabaja con $u=(S,I,R,D,\log\beta,\log\gamma,\log\delta)$. La matriz $H$ selecciona $h=Hu=(I,R,D)$. Con covarianzas muestrales calculadas sobre el pronóstico, resuelve $K(P_{hh}+V)=P_{uh}$ y aplica a cada miembro $u^a=u^f+K(Y_k+\varepsilon-h)$, con $\varepsilon\sim\mathcal N(0,V)$. Las tasas se corrigen mediante sus covarianzas con los canales observables.
4. Se recuperan las tasas mediante exponenciación y se corrigen los compartimientos con `_balance_compartments`. Se llevan a cero los valores negativos de $I,R,D$, se reescalan proporcionalmente si suman más que $N$ y se reconstruye $S=N-I-R-D$. La integración genera flujos no negativos. En la secuencia de análisis, la corrección puede revisar a la baja $R$ o $D$ respecto del día anterior.

`run_sird_enkf` devuelve medias y desviaciones estándar muestrales de forma $(T,7)$ y las trayectorias de tasas de forma $(T,M,3)$. Con `return_forecast=True` agrega la media pronosticada de $(I,R,D)$, de forma $(T,3)$ y con `NaN` inicial. La configuración experimental usa observaciones reconstruidas completas; la rutina valida una matriz finita y no negativa. El análisis del día $k$ usa los datos hasta $k$, que es filtrado. En suavizado se emplearían además datos de fechas posteriores para estimar ese estado.

**Lectura de resultados.** Las bandas de una desviación estándar describen dispersión entre miembros bajo los supuestos elegidos. El cuaderno por fases compara el pronóstico y el análisis del mismo filtro y exporta las tasas diarias. El código escribe `Datos/first_estimation_daily_parameters.csv`. La copia local consultada se encuentra en `Datos/simulated or estimated/first_estimation_daily_parameters.csv` y contiene 1 167 filas.

Al cierre de la sexta ola (día 1127), el cociente de las tasas medias fue 0.8262. En el día 1166 las tasas medias fueron $(0.087842,0.224029,0.00032134)\ \mathrm{día}^{-1}$, con cociente 0.3915. Las medianas por fase del recorrido se calculan tomando primero el resumen del ensamble por fecha y después la mediana temporal de cada columna. Por ejemplo, el cociente mediano del preconfinamiento fue 10.63, acompañado de una recuperación mediana de $0.0191\ \mathrm{día}^{-1}$, en un tramo donde la imputación impone ceros iniciales en $R$.

La comparación del recorrido con De León y colaboradores corresponde a un ajuste SEIARD con transmisión sintomática y asintomática separadas, otra ventana y otra relación con los registros. El recorrido documenta qué coeficientes del ajuste de la primera ola representan amplitudes o límites de curvas y qué rangos de la tabla siguen con procedencia incierta. Consulta esa sección para la comparación académica. Evita tratar las variantes del calendario como explicación causal de las trayectorias o los rangos de referencia como tasas verdaderas.

## Mapa de cuadernos y archivos

| Fuente en `Codigos/` | Qué hace y qué utiliza |
| --- | --- |
| [Lectura.ipynb](../Codigos/Lectura.ipynb) | Lee los cuatro CSV anuales, selecciona columnas, concatena y exporta el consolidado. |
| [Exploracion.ipynb](../Codigos/Exploracion.ipynb) | Revisa tipos, fechas, nulos, duplicados, clasificaciones, positivos y defunciones del consolidado. |
| [get_filtered_by_date.ipynb](../Codigos/get_filtered_by_date.ipynb) | Genera los CSV de preconfinamiento y seis olas por fecha de ingreso y el CSV con metadatos de fases. |
| [plots.ipynb](../Codigos/plots.ipynb) | Grafica casos positivos y defunciones acumulados del preconfinamiento, incluida una vista logarítmica de casos. |
| [observations.ipynb](../Codigos/observations.ipynb) y [observations_plots.ipynb](../Codigos/observations_plots.ipynb) | Generan y visualizan la reconstrucción exploratoria con $\gamma$ por fase. |
| [first_simulation.ipynb](../Codigos/first_simulation.ipynb) | Simulación determinista hasta el día 23 con $(0.3508,0.0714,0.0015)$ y construcción de observaciones del preconfinamiento con retraso de 14 días. |
| [pre-lockdown_sird_simulation.ipynb](../Codigos/pre-lockdown_sird_simulation.ipynb) | Parte del primer estado reconstruido, compara escenarios de tasas fijas hasta el día 100 y varía una tasa a la vez alrededor de las otras dos en sus puntos medios. Exporta tablas y GIF. El preconfinamiento histórico termina en el día 23. |
| [complete_sird_simulation.ipynb](../Codigos/complete_sird_simulation.ipynb) | Proyecta tres combinaciones fijas (mínimos, puntos medios y máximos simultáneos) hasta el día 1166. Exporta trayectorias, resumen y una animación. Es un ejercicio de escenarios deterministas. |
| [first_estimation.ipynb](../Codigos/first_estimation.ipynb) | Corre el EnKF nacional continuo, grafica compartimientos, tasas, dispersión y cociente reproductivo y muestra los cierres por periodo. |
| [first_estimation_by_phase.ipynb](../Codigos/first_estimation_by_phase.ipynb) | Repite la corrida continua conservando pronóstico y análisis completos. Muestra ambos y su corrección por fase y exporta las estimaciones diarias. |
| [enkf_30_day_walkthrough.ipynb](../Codigos/enkf_30_day_walkthrough.ipynb) | Ejercicio didáctico separado de 30 fechas (28 de febrero a 28 de marzo de 2020) y 50 miembros. Combina los registros de dos fases y guarda las operaciones de las primeras diez fechas, incluida la inicialización. `member_trace` muestra los pasos de un miembro y las magnitudes compartidas. |
| [random.ipynb](../Codigos/random.ipynb) | Inspección auxiliar de CSV exploratorios mediante rutas absolutas locales. |

El [recorrido en Markdown](Codigos/recorrido_primera_estimacion.md) sigue al miembro 1 de la corrida **nacional de 160 miembros** durante los tres primeros avances y muestra pronóstico, innovación, ganancia y corrección. El cuaderno de 50 miembros tiene su propio ensamble y su propia ganancia. `kalman_base/` conserva ejemplos didácticos del filtro clásico; `animations/` contiene utilidades y resultados de visualización.

## Barrido del tamaño del ensamble

[enkf_sweep.py](../Codigos/enkf_sweep.py) ejecuta el experimento, [enkf_sweep_data.py](../Codigos/enkf_sweep_data.py) prepara los dos escenarios y su configuración, y [enkf_sweep_report.py](../Codigos/enkf_sweep_report.py) agrega y grafica resultados. Se probaron tamaños de 10 a 200 en pasos de 10, semillas de 2026 a 2035, y dos escenarios, para 400 corridas planificadas. La referencia es 160 miembros. Dentro de cada escenario se mantienen las observaciones y la configuración del filtro.

- **México.** Mismas observaciones reconstruidas y calendario que la primera estimación. Las métricas se calculan por fase y globalmente.
- **Simulación.** Una trayectoria SIRD de 240 días transcurridos (241 fechas, días 0 a 240) con tasas positivas sinusoidales en escala logarítmica. Se observa $I,R,D$ con ruido gaussiano fijo entre corridas (semilla 2025). Las tasas verdaderas se guardan para evaluar recuperación de parámetros en ese escenario.
- **Métricas.** Se excluye el día 0 porque inicializa el filtro. Se comparan MAE y RMSE del pronóstico a un día y del análisis frente a las observaciones. El pronóstico de defunciones diarias se calcula como $D^f_k-D^a_{k-1}$ y se compara con $\Delta D$ observado. En sintético se calcula el RMSE entre el logaritmo de la tasa media estimada y el logaritmo de la tasa verdadera, y la fracción de días en que la tasa verdadera queda entre los percentiles 5 y 95 de los miembros. Las medianas y los cuartiles de métricas por tamaño se calculan entre semillas exitosas.

**Resultados del barrido local `5bda922378a240eb90eeed7cbe943e4e`.** Terminaron 395 de 400 corridas. Los cinco fallos ocurrieron en México (tres con 10 miembros y dos con 20). En México, el menor RMSE mediano de pronóstico de activos se obtuvo con 30 miembros (15 990 personas, frente a 18 372 con 160); el menor RMSE mediano de defunciones diarias se obtuvo con 140 (32.3 personas). Para 160 miembros, del 22 al 24 de marzo de 2023 se concentró cerca del 66 % del error cuadrático de activos, como mediana entre semillas según el README. En la corrida de referencia, el 22 de marzo se pronosticaron 642 309 activos frente a 35 345 reconstruidos y el análisis fue 35 345 al redondear. En sintético el RMSE mediano de pronóstico de activos pasó de 1 934 a 593 personas entre 10 y 200 miembros. Los mínimos de RMSE de tasas en esa simulación correspondieron a 10, 20 y 140 miembros para $\beta,\gamma,\delta$, respectivamente. El criterio, el escenario y la variación entre semillas condicionan cualquier comparación de tamaños.

Para reproducir o verificar, consulta los comandos del README. El piloto sintético se ejecuta con `python Codigos/enkf_sweep.py --pilot --scenario synthetic` y usa 10, 160 y 200 miembros con semilla 2026. Las dependencias están en `Codigos/requirements-enkf-sweep.txt`. Las pruebas se ejecutan con `python -m unittest discover -s Codigos -p test_enkf_sweep.py`. MLflow guarda configuraciones, huellas, series, informes y fallos bajo `Datos/enkf_sweep/`. Consulta `mlflow.db`, los artefactos de la corrida y `reports/<id>/` para recuperar evidencia existente. El directorio del informe contiene `counts.csv`, `summary.csv` y los resúmenes de dispersión.

## Estado de la tesis y límites de interpretación

- En `ayuda/formato-cap1.tex` se encuentran la pregunta, los objetivos y la voz del borrador. Los capítulos 2 y 3 se trabajan en `ayuda/formato-cap2.tex` y `ayuda/formato-cap3.tex`, con bibliografía en `ayuda/tesis-capitulos.bib`. El capítulo 2 pasa de epidemias al SIRD y al problema de observación. El capítulo 3 desarrolla Markov, filtrado y suavizado bayesianos, filtro de Kalman y EnKF, y enlaza con el capítulo 4, «Datos». La propiedad de Markov se refiere a la dependencia condicional del futuro respecto del estado presente; la homogeneidad temporal se refiere a la constancia del núcleo de transición.
- El borrador local del capítulo 3 incluye una sección adicional de **formulación del SIRD en espacio de estados**. Ilustra un estado de nueve componentes que añade dos flujos diarios, observación de conteos, canales ausentes y una posible proyección de compartimientos. Es una formulación teórica distinta de la implementación nacional de siete componentes que asimila niveles $(I,R,D)$. Contrasta el índice fuente con la meta de seis secciones indicada en los archivos de instrucciones antes de editar o compilar.
- Para redacción o revisión consulta las [instrucciones locales de tesis](ayuda/AGENTS.md), el capítulo 1 y las revisiones de `markdowns/ayuda/`. Conserva la voz académica en español, `\cite{}`, `\label` y el objetivo bibliográfico de seis referencias en el capítulo 2 y nueve en el 3. El código y sus comentarios nuevos se escriben en inglés. `Tesis_Mestria/` es el destino original o publicado; cualquier edición allí requiere autorización explícita.
- La orden documentada de compilación es `latexmk -pdf -outdir=.compilacion -interaction=nonstopmode -halt-on-error capitulos-2-y-3.tex` desde `ayuda/`. En el árbol consultado falta ese archivo maestro, y `Tesis_Mestria/` carece de los `.tex` originales mencionados en las instrucciones. Si la tarea requiere compilar, verifica primero la disponibilidad del maestro y la plantilla. Comprueba índice, referencias, bibliografía y el objetivo aproximado de 10 a 15 páginas de cuerpo.
- **Identificabilidad.** La estructural depende del modelo y de lo que podría observarse bajo condiciones ideales. La práctica depende además de la calidad de los registros, su cobertura y de los supuestos de reconstrucción, error y evolución de las tasas. Un buen ajuste tras la corrección y una banda estrecha entre miembros son compatibles con distintas trayectorias de tasas. El barrido aporta sensibilidad al tamaño del ensamble y resultados condicionados a un sintético con verdad conocida. La capacidad de separar las tres tasas en los registros reales sigue siendo la pregunta de investigación.
