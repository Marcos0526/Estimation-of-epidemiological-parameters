# Primera estimación de las tasas del SIRD

El [cuaderno de estimación](first_estimation.ipynb) aplica un EnKF continuo desde el 28 de febrero de 2020 (día 0) hasta el 9 de mayo de 2023 (día 1166). Conserva el ensamble al pasar del periodo previo al confinamiento a cada una de las seis olas. Para reproducirlo, abre el cuaderno con el intérprete de `uni` y ejecuta las celdas en orden. La primera celda de código reúne el intervalo, el archivo y los hiperparámetros. Las funciones del filtro están al final de [codes.py](codes.py).

## De los registros a las tres observaciones

Los CSV de `Datos/raw_data/` contienen registros individuales. En la muestra del archivo consolidado, la fecha de defunción ausente aparece como `9999-99-99`. Los CSV de `Datos/waves/` contienen registros individuales seleccionados por fecha de ingreso. El calendario de `Datos/metadata/covid_waves_parameters.csv` deja abril de 2022 entre dos olas y sitúa el cierre de la sexta el 31 de marzo de 2023 (día 1127). Para cubrir las 1167 fechas, el cuaderno lee tres columnas del consolidado en bloques de 250 000 registros y filtra las fechas de ingreso. Así incluye también los registros de abril de 2022 y los 39 días posteriores al cierre de la sexta ola.

La función existente `get_Y_obs` selecciona los casos positivos (clasificación final 1, 2 o 3) y construye una fila por día, incluidos los días sin eventos. Cuenta los ingresos positivos acumulados $C$, acumula las defunciones por su fecha registrada para obtener $D$ y asigna la recuperación a los 14 días del ingreso a quienes carecen de fecha de defunción. Al acumular esas recuperaciones obtiene $R$. Finalmente calcula los casos activos $I=C-R-D$. El filtro recibe estas tres series en el orden $[I,R,D]$. Cada bloque de lectura se agrega sobre el mismo calendario completo y luego se suman sus conteos diarios de activos y acumulados. Una defunción o recuperación de una persona ingresada al final de una ola se contabiliza en la fecha de su desenlace, incluso si pertenece a la ola siguiente. La primera tabla del cuaderno indica el inicio, el fin y los días de simulación de cada periodo.

## Construcción del ensamble inicial, paso a paso

1. Se fijan la población $N=126\,014\,024$, el tamaño del ensamble $M=160$ y la semilla 2026. `run_sird_enkf` crea un solo generador aleatorio y lo comparte entre inicialización, predicción y análisis.
2. `initialize_sird_ensemble` reserva una matriz de 160 filas y siete columnas. Cada fila representa una realización completa de $[S,I,R,D,\beta,\gamma,\delta]$. El primer dato del intervalo es $[2,0,0]$.
3. Para cada miembro se sortean $I,R,D$ de forma independiente alrededor de ese primer dato. `STATE_STD = [1, 0, 0]` da una desviación de una persona a los activos y fija ambos acumulados iniciales en cero. Las realizaciones son continuas, como los compartimientos de la ecuación diferencial.
4. `_balance_compartments` lleva a cero los recuentos negativos. Si $I+R+D>N$, reduce esos tres valores en la misma proporción hasta que sumen $N$. A continuación obtiene $S=N-I-R-D$. Con ello cada miembro comienza con la misma población. Estos ajustes de frontera pueden desplazar la media respecto al centro del sorteo.
5. Para cada miembro y cada tasa se sortea un valor normal en escala logarítmica, centrado en el logaritmo de la mediana inicial. Las medianas son $0.3508$, $0.0714$ y $0.0015$ por día, tomadas de [first_simulation.ipynb](first_simulation.ipynb). Las desviaciones iniciales en esa escala son $0.35$, $0.40$ y $0.50$. Al exponenciar, cada tasa es positiva. `_rate_triplet` comprueba la forma, finitud y signo de las ternas de configuración.
6. Se registran la media y la desviación estándar muestral de las siete columnas y se guardan las tres tasas de cada miembro. Como el primer dato ya intervino en la construcción del ensamble, el primer análisis de Kalman corresponde al segundo día.

## Qué ocurre entre dos fechas

**Avance de las tasas.** `gaussian_rate_walk` recibe una matriz de tasas físicas de tamaño $(M,3)$, los días inicial y final del paso y el generador aleatorio. Para cada miembro y tasa $p$, genera una innovación independiente $z\sim\mathcal N(0,1)$ y aplica

$$
\log p_k^f=\log p_{k-1}^a+\sigma_p\sqrt{t_k-t_{k-1}}\,z .
$$

Los superíndices $f$ y $a$ indican pronóstico y análisis. Las tres dispersiones del paso son $0.12$, $0.12$ y $0.15$ por raíz de día. Controlan cuánto pueden cambiar las log-tasas entre observaciones. Las dispersiones de inicialización controlan, en cambio, cuánto difieren los miembros en la primera fecha.

**Predicción del SIRD.** `predict_sird_ensemble` llama a la ley de avance y mantiene constantes las nuevas tasas durante ese paso. Integra cada miembro con `solve_ivp` y la función existente `SIRD`, partiendo de sus propios compartimientos analizados. El pronóstico contiene los cuatro compartimientos al final del paso y las tres tasas que los produjeron. Así puede surgir covarianza entre una tasa y una cantidad observable.

**Análisis.** `analyze_sird_ensemble` trabaja con $u=[S,I,R,D,\log\beta,\log\gamma,\log\delta]$. La matriz $H$ de tamaño $3\times7$ selecciona las columnas de $I,R,D$. Para cada miembro $m$ se obtiene $h_m=Hu_m^f$. A partir de las anomalías respecto de las medias del ensamble se calculan las covarianzas muestrales $P_{uh}$ y $P_{hh}$, con divisor $M-1$. La actualización es

$$
K=P_{uh}(P_{hh}+R_{\mathrm{obs}})^{-1},\qquad
u_m^a=u_m^f+K\bigl(Y_k+\varepsilon_m-h_m\bigr),\qquad
\varepsilon_m\sim\mathcal N(0,R_{\mathrm{obs}}).
$$

El código resuelve un sistema lineal para obtener la ganancia. Genera perturbaciones independientes entre miembros con la factorización de Cholesky de $R_{\mathrm{obs}}$. El cuaderno supone errores de observación con desviaciones de 10, 4 y 1.5 personas para $I,R,D$ y usa sus **varianzas** en la diagonal. Se trata de una hipótesis de trabajo. La reconstrucción $I=C-R-D$ puede correlacionar los errores de las tres series. La función admite también una covarianza simétrica definida positiva con términos fuera de la diagonal.

Tras el análisis se exponencian las log-tasas y se vuelve a aplicar el balance de compartimientos. La actualización puede revisar a la baja un acumulado estimado respecto al análisis del día anterior, aunque el SIRD produzca flujos diarios no negativos en cada integración.

### Cómo cambiar la ley de avance

El argumento `rate_step` acepta cualquier función con la firma `rate_step(rates, start_day, end_day, rng)`. Recibe y devuelve tasas **en unidades físicas**, finitas, estrictamente positivas y con forma $(M,3)$. Los tiempos se expresan en días desde un origen común. El cuaderno configura la ley gaussiana mediante `partial(gaussian_rate_walk, log_step_std=LOG_STEP_STD)`. Por ejemplo, para mantener las tasas de cada miembro durante la predicción basta con pasar esta función a `run_sird_enkf`.

```python
def constant_rate_step(rates, start_day, end_day, rng):
    return rates.copy()
```

El análisis aún podrá actualizar esas tasas con las observaciones. Una ley de reversión a la media puede usar el tiempo transcurrido y el generador. Una ley por tramos puede seleccionar sus valores a partir de los días del paso.

## Lectura de las salidas

`run_sird_enkf` recorre las fechas y devuelve tres arreglos, con $T$ igual al número de observaciones.

| Salida | Forma | Contenido |
| --- | --- | --- |
| `mean` | $(T,7)$ | Media por fecha de cada componente en el orden del estado aumentado y en sus unidades físicas. |
| `std` | $(T,7)$ | Desviación estándar entre miembros, con divisor $M-1$. |
| `parameter_paths` | $(T,M,3)$ | Trayectoria filtrada de las tres tasas de cada miembro. Por ejemplo, `parameter_paths[:, 0, 2]` contiene la serie de $\delta$ del primer miembro. |

En el cuaderno, `plot_band` dibuja la media y la banda de una desviación estándar. La primera figura superpone los datos a los compartimientos estimados y la segunda muestra las tasas. Las bandas describen la dispersión obtenida bajo los supuestos del ensamble.

La última celda grafica para todas las fechas el cociente $R_0(t)=\bar\beta(t)/(\bar\gamma(t)+\bar\delta(t))$, calculado con las tasas medias del análisis. Su interpretación supone susceptibilidad completa y tasas congeladas en el valor de esa fecha. Para el número reproductivo efectivo se incluye además el factor $S(t)/N$. La tabla de esa celda muestra las tasas y el cociente al cierre de cada periodo, incluidos el fin de la sexta ola y el día 1166.

Al cierre de la sexta ola (día 1127), la corrida obtiene $R_0(t)=0.8262$. Al cierre del horizonte solicitado (día 1166), las medias diarias son $\bar\beta=0.087842$, $\bar\gamma=0.224029$ y $\bar\delta=0.00032134$, con $R_0(t)=0.3915$. El plazo de recuperación imputado y las dispersiones elegidas influyen en estas estimaciones.
