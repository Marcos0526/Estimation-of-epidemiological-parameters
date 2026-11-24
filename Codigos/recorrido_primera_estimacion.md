# Del registro de un caso a la estimación diaria de β, γ y δ

¿Cómo puede un ingreso hospitalario o ambulatorio cambiar una estimación de contagio? El registro aporta una fecha y una clasificación. Al agrupar los registros por fecha se construye una serie diaria. El SIRD proyecta cuánto debería valer esa serie al día siguiente para cada una de 160 combinaciones de estados y tasas. El EnKF usa los errores de esas proyecciones y sus covarianzas para corregir tanto los compartimientos como las tasas. Seguiremos a **un mismo miembro** desde el viernes 28 de febrero hasta el lunes 2 de marzo de 2020. Luego situaremos esos tres avances en las 1167 fechas de la corrida nacional. [Cuaderno de la primera estimación](first_estimation.ipynb) · [Implementación](codes.py)

## 1. ¿Qué información entra al filtro?

El cuaderno lee el consolidado de México por bloques de 250 000 registros, selecciona ingresos en el periodo de estudio y conserva los positivos con clasificación final 1, 2 o 3. Cada bloque aporta conteos al mismo calendario. Si una persona tiene fecha de defunción, se acumula en $D$ en esa fecha. Si carece de ella, se imputa su recuperación 14 días después del ingreso. La suma de los ingresos positivos hasta el día da $C$. Así se obtiene $Y=(I,R,D)=(C-R-D,R,D)$. La fecha de ingreso representa el registro del caso y puede diferir de la fecha de infección.

| Fecha | Nuevos ingresos positivos | Acumulado C | Recuperaciones acumuladas R | Defunciones acumuladas D | Activos I |
|---|---:|---:|---:|---:|---:|
| 2020-02-28 | 2 | 2 | 0 | 0 | 2 |
| 2020-02-29 | 2 | 4 | 0 | 0 | 4 |
| 2020-03-01 | 2 | 6 | 0 | 0 | 6 |
| 2020-03-02 | 5 | 11 | 0 | 0 | 11 |

El 2 de marzo hay cinco ingresos nuevos y once casos activos reconstruidos. El filtro asimila los once activos de ese día. Los ceros de recuperación de estas primeras fechas se explican por el plazo de imputación y por el inicio del calendario.

## 2. ¿Cómo se construye la creencia inicial?

Cada miembro contiene $x=(S,I,R,D,β,γ,δ)$. Los cuatro primeros componentes se miden en personas y los tres últimos en día⁻¹. Fijamos $N=126 014 024$ y la semilla 2026. El 28 de febrero se sortean los tres recuentos alrededor de $(2,0,0)$ con desviaciones $(1,0,0)$ personas. Se llevan a cero los recuentos negativos y se calcula $S=N-I-R-D$. Las tasas iniciales se sortean independientemente como $p=\exp(\log m+sZ)$, con $Z$ normal estándar. Las medianas $m$ son $(0.3508,0.0714,0.0015)$ y las dispersiones logarítmicas $s$ son $(0.35,0.40,0.50)$. Las medianas proceden del ejemplo previo del proyecto. [Primera simulación](first_simulation.ipynb)

Para el **miembro 1** se obtiene `(126014022.793122, 1.206878, 0.000000, 0.000000, 0.213153, 0.076420, 0.001755)`. La media de las tasas de los 160 miembros es `(0.373143, 0.079829, 0.001657)`. El primer dato ya intervino en el sorteo inicial, por lo que la primera corrección ocurre el 29 de febrero.

## 3. ¿Qué operaciones se repiten cada día?

**Primero avanzan las tasas.** Para cada miembro y tasa se sortea un incremento independiente $e_p=\sigma_p\sqrt{\Delta t}\,Z_p$. Se obtiene $p^f=p^a_{\rm anterior}\exp(e_p)$. El superíndice $f$ denota pronóstico y $a$ análisis corregido. Para un día, las desviaciones de los incrementos son $(0.12,0.12,0.15)$. Esta caminata tiene deriva cero en el logaritmo. La multiplicación por $\exp(e_p)$ mantiene positivas las tasas y hace que el tamaño del cambio dependa de su valor anterior.

**Después avanzan los compartimientos.** Con las nuevas tasas constantes durante ese día, el SIRD calcula tres flujos instantáneos en personas por día:

$$
f=\beta SI/N,\qquad g=\gamma I,\qquad m=\delta I,\qquad
(\dot S,\dot I,\dot R,\dot D)=(-f,\ f-g-m,\ g,\ m).
$$

Un contagio pasa de $S$ a $I$. Las recuperaciones y los fallecimientos compiten por la salida de $I$. La suma de las cuatro derivadas es cero. El integrador RK45 evalúa estos flujos a lo largo del día y controla el error con pasos internos adaptativos. Si $F,G,M$ son las integrales diarias de los tres flujos, el avance es $(-F,F-G-M,G,M)$. En las cuentas de abajo recuperamos esas integrales de las diferencias entre los compartimientos pronosticados e iniciales.

**Por último se asimila el dato.** Para cada miembro se usa $u=(S,I,R,D,\log\beta,\log\gamma,\log\delta)$. La matriz $H$ toma sus componentes 2, 3 y 4 (contando desde 1), así que $h=Hu=(I,R,D)$. Se centran las 160 filas respecto de sus medias para formar $A_m=u_m-\bar u$ y $B_m=h_m-\bar h$. La ganancia usa a todos los miembros:

$$
P_{uh}=A^\top B/159,\quad W=B^\top B/159+V,\quad
V=\operatorname{diag}(100,16,2.25),\quad KW=P_{uh}.
$$

$V$ es la covarianza de observación en personas². El código resuelve el sistema para $K$. Para el miembro que seguimos se sortea $\epsilon=\operatorname{diag}(10,4,1.5)z$ con tres normales estándar independientes. La innovación es $v=Y+\epsilon-h$. El incremento de análisis es $Kv$ y el estado corregido en escala de trabajo es $u^a=u^f+Kv$. Se exponencian las tres últimas componentes, se llevan a cero los recuentos negativos de $I,R,D$ y, si exceden juntos a $N$, se reescalan proporcionalmente. Finalmente se reconstruye $S=N-I-R-D$. La perturbación expresa el error supuesto para la medición. El registro público conserva su valor original.

## 4. Tres avances consecutivos con los datos mexicanos

Para leer las tablas, sigue al miembro 1 desde su estado corregido de ayer hasta su nuevo estado corregido de hoy. Los vectores y matrices se redondean para lectura. El generador mantiene la precisión completa y el orden de sorteos del cuaderno. Las observaciones son conteos de los registros, mientras que los decimales en los estados y tasas son resultados de la simulación.

### 28/02 → 29/02/2020

| Operación del miembro 1 | Sustitución o resultado |
|---|---|
| Punto de partida $(S,I,R,D)$ | `(126014022.793122, 1.206878, 0.000000, 0.000000)` |
| Tasas anteriores × exponencial del incremento, componente a componente | `(0.213153, 0.076420, 0.001755) × exp(-0.013917, 0.011288, 0.190488) = (0.210207, 0.077287, 0.002123)` |
| Flujos iniciales $(βSI/N,γI,δI)$ | `(0.253695, 0.093277, 0.002563)` personas/día |
| Flujos integrados $(F,G,M)$ | `(0.271034, 0.099652, 0.002738)` personas |
| Pronóstico $(S-F,I+F-G-M,R+G,D+M)$ | `(126014022.522089, 1.375522, 0.099652, 0.002738)` |
| Pronóstico observable del miembro $h$ y media observable de los 160 miembros | `(1.375522, 0.099652, 0.002738)` y `(2.757384, 0.185658, 0.004022)` |
| Observación real + ruido − pronóstico = innovación | `(4.000000, 0.000000, 0.000000) + (-1.048411, 3.086765, -0.899415) − (1.375522, 0.099652, 0.002738) = (1.576068, 2.987113, -0.902153)` |

Las columnas de las matrices siguen el orden $I,R,D$. Las filas de $P_{uh}$ y $K$ siguen el orden de $u$. Cada fila de la derecha resulta de resolver el sistema con la fila correspondiente de la izquierda.

```text
W = cov(h,h) + V
(102.288404, 0.122049, 0.003368)
(0.122049, 16.013803, 0.000186)
(0.003368, 0.000186, 2.250012)
componente       P_uh                                  K
S           (-2.413821, -0.136039, -0.003566)      (-0.023588, -0.008315, -0.001549)
I           (2.288404, 0.122049, 0.003368)         (0.022363, 0.007451, 0.001463)
R           (0.122049, 0.013803, 0.000186)         (0.001192, 0.000853, 0.000081)
D           (0.003368, 0.000186, 0.000012)         (0.000033, 0.000011, 0.000005)
log β       (0.195496, 0.004410, 0.000213)         (0.001911, 0.000261, 0.000092)
log γ       (-0.063485, 0.026965, -0.000032)       (-0.000623, 0.001689, -0.000013)
log δ       (0.075323, 0.004563, 0.001187)         (0.000736, 0.000279, 0.000527)
```

Por ejemplo, la fila de $I$ da la corrección **+0.056183 personas** al multiplicar `(0.022363, 0.007451, 0.001463)` por la innovación `(1.576068, 2.987113, -0.902153)` y sumar los tres productos. El resultado para las siete componentes aparece a continuación.

| Operación | Valor del miembro 1 |
|---|---|
| Incremento $Kv$ en los cuatro recuentos y las tres log-tasas | `(-0.060618, 0.056183, 0.004354, 0.000081, 0.003708, 0.004075, 0.001519)` |
| Suma $u^f+Kv$, antes de aplicar las restricciones | `(126014022.461471, 1.431705, 0.104005, 0.002819, -1.555952, -2.556149, -6.153184)` |
| Resultado en unidades físicas, después de exponenciar las tasas y equilibrar los recuentos | `(126014022.461471, 1.431705, 0.104005, 0.002819, 0.210988, 0.077603, 0.002127)` |
| Media de los 160 miembros corregidos, en el mismo orden | `(126014021.055900, 2.754780, 0.185302, 0.004017, 0.369445, 0.079841, 0.001651)` |


¿Por qué el miembro 1 termina con 1.431705 activos cuando el registro indica 4? El pronóstico fue $1.206878+0.271034-0.099652-0.002738=1.375522$. Para este miembro, la innovación en activos fue $4-1.048411-1.375522≈1.576068$. Al ponderar las tres innovaciones, la corrección fue $0.022363(1.576068)+0.007451(2.987113)+0.001463(-0.902153)≈0.05618$. La varianza pronosticada de activos entre miembros fue 2.288404 personas², frente a las 100 personas² asignadas al error de observación. Con estos supuestos se obtiene una corrección pequeña. Al día siguiente, el miembro vuelve a integrarse desde $1.375522+0.056183=1.431705$ activos.

### 29/02 → 01/03/2020

| Operación del miembro 1 | Sustitución o resultado |
|---|---|
| Punto de partida $(S,I,R,D)$ | `(126014022.461471, 1.431705, 0.104005, 0.002819)` |
| Tasas anteriores × exponencial del incremento, componente a componente | `(0.210988, 0.077603, 0.002127) × exp(-0.244124, -0.100859, -0.092635) = (0.165286, 0.070158, 0.001939)` |
| Flujos iniciales $(βSI/N,γI,δI)$ | `(0.236641, 0.100445, 0.002775)` personas/día |
| Flujos integrados $(F,G,M)$ | `(0.248018, 0.105274, 0.002909)` personas |
| Pronóstico $(S-F,I+F-G-M,R+G,D+M)$ | `(126014022.213453, 1.571540, 0.209280, 0.005728)` |
| Pronóstico observable del miembro $h$ y media observable de los 160 miembros | `(1.571540, 0.209280, 0.005728)` y `(3.834462, 0.436607, 0.009567)` |
| Observación real + ruido − pronóstico = innovación | `(6.000000, 0.000000, 0.000000) + (17.583458, -5.085011, -0.283825) − (1.571540, 0.209280, 0.005728) = (22.011918, -5.294291, -0.289553)` |

Las columnas de las matrices siguen el orden $I,R,D$. Las filas de $P_{uh}$ y $K$ siguen el orden de $u$. Cada fila de la derecha resulta de resolver el sistema con la fila correspondiente de la izquierda.

```text
W = cov(h,h) + V
(106.762878, 0.496388, 0.015147)
(0.496388, 16.078427, 0.001162)
(0.015147, 0.001162, 2.250068)
componente       P_uh                                  K
S           (-7.274413, -0.575976, -0.016377)      (-0.067978, -0.033724, -0.006803)
I           (6.762878, 0.496388, 0.015147)         (0.063209, 0.028921, 0.006291)
R           (0.496388, 0.078427, 0.001162)         (0.004627, 0.004735, 0.000483)
D           (0.015147, 0.001162, 0.000068)         (0.000142, 0.000068, 0.000029)
log β       (0.554249, 0.029911, 0.000966)         (0.005183, 0.001700, 0.000393)
log γ       (-0.143646, 0.061370, -0.000211)       (-0.001363, 0.003859, -0.000086)
log δ       (0.119407, 0.004786, 0.002731)         (0.001117, 0.000263, 0.001206)
```

Por ejemplo, la fila de $I$ da la corrección **+1.236424 personas** al multiplicar `(0.063209, 0.028921, 0.006291)` por la innovación `(22.011918, -5.294291, -0.289553)` y sumar los tres productos. El resultado para las siete componentes aparece a continuación.

| Operación | Valor del miembro 1 |
|---|---|
| Incremento $Kv$ en los cuatro recuentos y las tres log-tasas | `(-1.315822, 1.236424, 0.076650, 0.002748, 0.104982, -0.050417, 0.022846)` |
| Suma $u^f+Kv$, antes de aplicar las restricciones | `(126014020.897631, 2.807964, 0.285929, 0.008476, -1.695095, -2.707424, -6.222973)` |
| Resultado en unidades físicas, después de exponenciar las tasas y equilibrar los recuentos | `(126014020.897631, 2.807964, 0.285929, 0.008476, 0.183582, 0.066708, 0.001983)` |
| Media de los 160 miembros corregidos, en el mismo orden | `(126014019.622481, 3.927037, 0.440708, 0.009774, 0.369635, 0.079449, 0.001668)` |


### 01/03 → 02/03/2020

| Operación del miembro 1 | Sustitución o resultado |
|---|---|
| Punto de partida $(S,I,R,D)$ | `(126014020.897631, 2.807964, 0.285929, 0.008476)` |
| Tasas anteriores × exponencial del incremento, componente a componente | `(0.183582, 0.066708, 0.001983) × exp(0.024445, -0.004978, 0.012881) = (0.188125, 0.066377, 0.002009)` |
| Flujos iniciales $(βSI/N,γI,δI)$ | `(0.528248, 0.186385, 0.005641)` personas/día |
| Flujos integrados $(F,G,M)$ | `(0.561175, 0.198002, 0.005993)` personas |
| Pronóstico $(S-F,I+F-G-M,R+G,D+M)$ | `(126014020.336456, 3.165143, 0.483931, 0.014469)` |
| Pronóstico observable del miembro $h$ y media observable de los 160 miembros | `(3.165143, 0.483931, 0.014469)` y `(5.679589, 0.804987, 0.017904)` |
| Observación real + ruido − pronóstico = innovación | `(11.000000, 0.000000, 0.000000) + (-14.702779, 7.152277, 1.107975) − (3.165143, 0.483931, 0.014469) = (-6.867922, 6.668346, 1.093506)` |

Las columnas de las matrices siguen el orden $I,R,D$. Las filas de $P_{uh}$ y $K$ siguen el orden de $u$. Cada fila de la derecha resulta de resolver el sistema con la fila correspondiente de la izquierda.

```text
W = cov(h,h) + V
(123.588047, 1.841624, 0.057395)
(1.841624, 16.284766, 0.004692)
(0.057395, 0.004692, 2.250259)
componente       P_uh                                  K
S           (-25.487066, -2.131082, -0.062346)     (-0.204611, -0.107718, -0.022263)
I           (23.588047, 1.841624, 0.057395)        (0.189485, 0.091654, 0.020482)
R           (1.841624, 0.284766, 0.004692)         (0.014665, 0.015828, 0.001678)
D           (0.057395, 0.004692, 0.000259)         (0.000461, 0.000236, 0.000103)
log β       (1.242281, 0.084002, 0.002459)         (0.009991, 0.004028, 0.000830)
log γ       (-0.322190, 0.108580, -0.000680)       (-0.002711, 0.006974, -0.000248)
log δ       (0.157294, 0.000896, 0.005160)         (0.001273, -0.000090, 0.002261)
```

Por ejemplo, la fila de $I$ da la corrección **-0.667789 personas** al multiplicar `(0.189485, 0.091654, 0.020482)` por la innovación `(-6.867922, 6.668346, 1.093506)` y sumar los tres productos. El resultado para las siete componentes aparece a continuación.

| Operación | Valor del miembro 1 |
|---|---|
| Incremento $Kv$ en los cuatro recuentos y las tres log-tasas | `(0.662604, -0.667789, 0.006664, -0.001479, -0.040852, 0.064853, -0.006868)` |
| Suma $u^f+Kv$, antes de aplicar las restricciones | `(126014020.999060, 2.497354, 0.490595, 0.012990, -1.711501, -2.647550, -6.216960)` |
| Resultado en unidades físicas, después de exponenciar las tasas y equilibrar los recuentos | `(126014020.999060, 2.497354, 0.490595, 0.012990, 0.180595, 0.070825, 0.001995)` |
| Media de los 160 miembros corregidos, en el mismo orden | `(126014016.403309, 6.699530, 0.876784, 0.020378, 0.392664, 0.079766, 0.001697)` |

La última línea de cada tabla es la estimación que aparece en el cuaderno. Una tasa se corrige cuando las variaciones de su logaritmo covarían con los pronósticos observables entre los miembros. Incluso con cero defunciones observadas, el ensamble puede mantener una tasa de mortalidad positiva. En esas primeras fechas la información directa sobre δ es escasa.

## 5. Mapa de funciones de la primera estimación

| Función o llamada | Papel en la cadena de cálculo |
|---|---|
| `get_Y_obs` | Llama a `dates_to_datetime` para normalizar fechas y a `covid_positive` para seleccionar clasificaciones 1–3. Cuenta ingresos y desenlaces, acumula por fecha y devuelve la matriz diaria y un DataFrame. |
| `run_sird_enkf` | Crea el generador con la semilla, inicializa, alterna pronóstico y análisis y devuelve medias y desviaciones muestrales de las siete componentes, además de los parámetros de cada miembro. |
| `initialize_sird_ensemble` | Usa `rng.normal` para los recuentos y `rng.lognormal` para las tasas. |
| `_rate_triplet` y `_balance_compartments` | La primera valida forma, finitud y signo de las ternas. La segunda impone recuentos factibles y reconstruye susceptibles durante inicialización y análisis. |
| `gaussian_rate_walk` | Sortea los incrementos de log-tasa con `rng.normal` y aplica el multiplicador exponencial. `functools.partial` fija las dispersiones y entrega la función intercambiable `rate_step`. |
| `predict_sird_ensemble` | Llama a `rate_step` y a `scipy.integrate.solve_ivp` para cada miembro. El integrador llama a `SIRD`, que devuelve las cuatro derivadas. |
| `analyze_sird_ensemble` | Proyecta con `SIRD_OBSERVATION_MATRIX`, centra con `mean`, forma las covarianzas por productos matriciales, resuelve la ganancia con `numpy.linalg.solve` y usa `numpy.linalg.cholesky` y `rng.normal` para perturbar la observación. Aplica `log`, `exp` y el balance. |
| Preparación y resumen en el cuaderno | `pandas.read_csv`, `to_datetime`, `between` y `date_range` delimitan bloques y calendario. `groupby`, `value_counts`, `reindex` y `cumsum` construyen conteos. `mean` y `std(ddof=1)` resumen el ensamble. |
| `plot_band` y `mark_periods` | La primera dibuja media y dispersión mediante `plot` y `fill_between`. La segunda marca comienzos de fase y periodos entre olas o posteriores con `axvline` y `axvspan`. Los puntos de observación usan `scatter`. |

## 6. ¿Qué estimó la corrida completa y cómo se compara con la tabla?

El CSV contiene una media de ensamble por fecha. Para representar **todos los días de cada fase**, calculamos la mediana temporal de cada columna. Las medianas son marginales, por lo que la mediana de los cocientes puede diferir del cociente de las medianas. Cada celda presenta **mediana EnKF / valor suministrado**. Las flechas conservan el sentido temporal indicado en la tabla de referencia. Los rangos suministrados sirven como referencias descriptivas y su procedencia debe evaluarse por separado. [CSV diario](../Datos/first_estimation_daily_parameters.csv) · [Metadatos de las fases](../Datos/metadata/covid_waves_parameters.csv)

| Fase y fechas de la corrida | β (día⁻¹) | γ (día⁻¹) | δ (día⁻¹) | Cociente reproductivo |
|---|---|---|---|---|
| Preconfinamiento · 2020-02-28 a 2020-03-22 | 0.2740 / 0.5534 – 1.1400 | 0.0191 / 0.0708 – 0.2000 | 0.001577 / 0.0561 – 0.1200 | 10.630 / 1.762 – 2.860 |
| Primera · 2020-03-23 a 2020-09-30 | 0.0951 / 0.5534 → 0.1069 | 0.0609 / 0.0708 → 0.1383 | 0.009066 / 0.0561 → 0.0211 | 1.159 / 1.100 - 1.300 |
| Segunda · 2020-10-01 a 2021-04-30 | 0.0881 / 0.3500 - 0.5000 | 0.0845 / 0.1200 - 0.1500 | 0.006078 / 0.0150 - 0.0250 | 0.959 / 1.200 - 1.500 |
| Tercera · 2021-05-01 a 2021-11-30 | 0.0855 / 0.6000 - 0.8500 | 0.0868 / 0.1400 - 0.1800 | 0.003507 / 0.0080 - 0.0150 | 0.942 / 1.400 - 1.800 |
| Cuarta · 2021-12-01 a 2022-03-31 | 0.0753 / 1.1000 - 1.6000 | 0.1030 / 0.2000 - 0.3300 | 0.001555 / 0.0010 - 0.0030 | 0.683 / 2.000 - 3.500 |
| Quinta · 2022-05-01 a 2022-10-31 | 0.0806 / 0.9000 - 1.3000 | 0.0906 / 0.2200 - 0.3000 | 0.000395 / 0.0005 - 0.0015 | 0.844 / 1.300 - 1.600 |
| Sexta · 2022-11-01 a 2023-03-31 | 0.1070 / 0.8000 - 1.2000 | 0.0849 / 0.2500 - 0.3300 | 0.000474 / 0.0003 - 0.0010 | 1.152 / 1.100 - 1.350 |

Abril de 2022 aporta 30 fechas entre la cuarta y la quinta ola. Del 1 de abril al 9 de mayo de 2023 hay otras 39 fechas después del cierre de la sexta. El filtro continúa durante ambos intervalos. En el día 1166 termina con β = 0.087842, γ = 0.224029 y δ = 0.00032134 día⁻¹. El cociente de las tasas medias es 0.3915. La columna `R0` usa susceptibilidad completa y congela las tasas en el valor del día. Para el número efectivo del SIRD se incorpora la fracción susceptible mediante $R_t=[β/(γ+δ)]S/N$. La última columna de la tabla aportada mezcla las etiquetas $R_0$ y $R(t)$, de modo que una coincidencia numérica requiere también verificar la definición.

### Qué se pudo verificar en artículos primarios

De León y colaboradores ajustaron un **SEIARD** para México con datos del 12 de marzo al 12 de agosto de 2020. En su sección 3.2 y tabla 1, la transmisión sintomática es $β_I(t)=0.5534 e^{-t/31.9108}+0.1069$. Por tanto, al inicio del reloj del artículo vale **0.6603** y a largo plazo tiende a **0.1069**. Para la mortalidad, la misma distinción da **0.0561 + 0.0211 = 0.0772** al inicio y **0.0211** como límite. La recuperación sube aproximadamente de **0.0708** a **0.1383**. El artículo modela por separado las transmisiones sintomática y asintomática. Su β sintomática y la β efectiva única del SIRD describen poblaciones distintas. Esta lectura corrige la interpretación de las dos flechas de transmisión y mortalidad de la primera ola. [De León et al., sección 3.2 y tabla 1](https://www.aimspress.com/fileOther/PDF/BDIA/bigdia-05-002.pdf)

Al alinear el calendario, el **12 de agosto de 2020** (153 días después del inicio del ajuste publicado), las funciones del SEIARD dan aproximadamente $(β_I,γ,δ)=(0.1115,0.1383,0.0211)$ día⁻¹. Ese mismo día el SIRD–EnKF obtuvo $(0.0838,0.0893,0.0078)$ día⁻¹. Las tres estimaciones quedan por debajo de sus contrapartes del SEIARD en esa fecha. El SEIARD ajustó sus curvas al intervalo completo, mientras que el EnKF actualizó las tasas secuencialmente sobre las series reconstruidas. [Funciones del SEIARD, sección 3.2](https://www.aimspress.com/fileOther/PDF/BDIA/bigdia-05-002.pdf) · [Estimaciones diarias del EnKF](../Datos/first_estimation_daily_parameters.csv)

El estudio *Transmission dynamics and forecasts of the COVID-19 pandemic in Mexico* estimó un número reproductivo efectivo con series de incidencia por inicio de síntomas y describió fluctuaciones alrededor de uno desde finales de julio de 2020. La mediana de nuestro cociente durante toda la primera ola fue 1.159. La ventana, el reloj de los casos y la definición del cociente son diferentes. Otro estudio de cinco olas mexicanas examinó casos confirmados y letalidad intrahospitalaria. Esta última expresa una proporción entre pacientes hospitalizados, mientras que δ en el SIRD tiene unidades de riesgo por infectado y por día. La documentación consultada no permite atribuir a estos artículos cada intervalo de las olas segunda a sexta ni los extremos del preconfinamiento de la tabla suministrada. [Estudio de transmisión en México](https://doi.org/10.1371/journal.pone.0254826) · [Estudio epidemiológico de cinco olas](https://doi.org/10.1186/s12879-022-07800-w)

### Conclusión

La concordancia depende del parámetro y de la fase. En la primera ola, el cociente mediano de 1.159 cae dentro de 1.100–1.300, aunque las medianas de las tres tasas quedan por debajo de ambos valores de sus respectivas flechas. En la sexta ola, la mortalidad mediana de 0.000474 día⁻¹ sí cae dentro de 0.0003–0.0010, pero la transmisión mediana de 0.107 día⁻¹ queda muy por debajo de 0.800–1.200. En el preconfinamiento el cociente mediano alcanza 10.63, acompañado de una recuperación mediana de apenas 0.0191 día⁻¹. Durante las primeras dos semanas la regla de recuperación impuesta por construcción produce ceros y condiciona la información con que el filtro separa las salidas de infectados.

En estos datos agregados, $I$ informa especialmente sobre el crecimiento neto $βS/N-γ-δ$. Las defunciones contribuyen información sobre $δI$, y la serie imputada de recuperaciones condiciona la estimación de $γI$. Además, el supuesto de errores de observación independientes es aproximado porque $I=C-R-D$ comparte información con las otras dos series. Por ello podemos describir el ajuste retrospectivo y localizar discrepancias con la literatura, pero el ajuste por sí solo deja abierta la identificabilidad de las tres tasas y la atribución de sus cambios a variantes específicas.
