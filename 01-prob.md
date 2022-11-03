
```r
library(tidyverse)
#> ── Attaching packages ─────────────────── tidyverse 1.3.2 ──
#> ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
#> ✔ tibble  3.1.8     ✔ dplyr   1.0.9
#> ✔ tidyr   1.2.0     ✔ stringr 1.4.1
#> ✔ readr   2.1.2     ✔ forcats 0.5.2
#> ── Conflicts ────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
library(kableExtra)
#> 
#> Attaching package: 'kableExtra'
#> 
#> The following object is masked from 'package:dplyr':
#> 
#>     group_rows
```

# Probabilidad {#probabilidad}

La estadística Bayesiana nace de Thomas Bayes (1702-1761), clérigo y científico-matemático vocacional dedicado a la probabilidad inversa en seguros.

Hasta el siglo XX y especialmente con el desarrollo computacional de la era moderna, no se pudo desarrollar por los problemas técnicos que implica su implementación (involucra integraciones irresolubles de forma manual), por lo que se impuso la versión frecuentista basada en el muestreo repetido. Ya en los últimos 20-30 años podríamos decir que despegó incorporándose en múltiples ámbitos de conocimiento y especialmente en los últimos 10 por los avances en hardware y sobre todo software, dados los requisitos de computación intensiva y simulación.

Realmente, el modo de razonar "en bayesiano" se corresponde con el modo que usamos habitualmente para razonar, con lo que podríamos decir que a priori todos somos bayesianos. Nuestras percepciones las solemos formular en términos de probabilidades e incertidumbres; así hablamos de fenómenos extraordinarios (con poca probabilidad), o de reacciones de sentido común (acontecen con alta probabilidad... en un entorno sensato), y actualizamos nuestras percepciones con la información que recibimos u observamos. 


::: {.exm #ex01a}
Por ejemplo, el pronóstico de lluvias se da con la "probabilidad de lluvia", y la incertidumbre te la intentan transmitir con los habituales avisos amarillos, naranjas... Cuando tú te asomas por la ventana, actualizas la información recibida según tus percepciones de las nubes, viento, temperatura... y tu experiencia (generalmente los trabajadores del campo suelen tener mucha más puntería para afinar la predicción). 
:::

::: {.exm #ex01b}
Cuando vas por la calle, sueles reconocer a los turistas extranjeros por su aspecto, especialmente los que vienen de países del norte de Europa, por supuesto los asiáticos,... Si te enseñan una foto de una persona, y te piden adivinar su procedencia, tu predicción la darás en términos de probabilidades: probabilidad de que sea inglesa, alemana, rusa, ....
:::

::: {.exm #ex01c}
Cuando prueban un nuevo fármaco para una afección sueles estar interesado en la probabilidad de que nos cure. Ante una enfermedad complicada solemos interesarnos por el tiempo de supervivencia o la probabilidad de superarla bajo nuestras circunstancias, y no tanto en valores genéricos basados en todos los enfermos que la han sufrido.
:::

## Fundamentos

1. Está construida sobre la base de que la información y la incertidumbre que tenemos sobre todo lo que conocemos y desconocemos es expresable en términos de  probabilidad.

1. No se acude a la idea del muestreo repetido para interpretar las propiedades de los estimadores y en consecuencia las estimaciones obtenidas.

1. La información observada actualiza el conocimiento sobre lo desconocido.

1. Desaparecen los estimadores y se infiere en términos de distribuciones de probabilidad sobre los parámetros de interés.


## Conceptos
Es importante concretar los conceptos básicos sobre los que trabajemos, para asegurarnos de que todos entendemos lo mismo.

Un **parámetro** es una característica poblacional de interés, desconocida, sobre la que pretendemos inferir a partir de datos observados.

**Ejemplo**

1. La probabilidad de que una planta se infecte por un hongo. 
1. El valor medio/máximo de cierto índice de contaminación medioambiental en una ubicación. 


Las **variables** identifican las características que hemos de observar en la población de interés para obtener información sobre el parámetro de interés.

**Ejemplo**

1. En cada planta hay que observar si tiene/no tiene el hongo; se trata de una variable dicotómica.
1. En la ubicación dada hay que registrar el índice de interés, o las características con las que se calcula; se trata de una variable numérica (continua). 

Los **datos** disponibles provienen de observaciones de la característica de interés sobre un sector/muestra de la población de interés.

**Ejemplo**

1. En una finca revisamos al azar un conjunto de plantas y registramos si tiene o no el hongo.  Obtendremos un conjunto de valores 0/1.
1. En la ubicación dada elegiremos un periodo de cadencia para ir registrando el índice de interés a lo largo del tiempo.


La **verosimilitud** da información sobre qué valores del parámetro son más/menos verosímiles a la vista de los datos. Es la información "probabilística" que proporciona sobre el parámetro de interés la distribución asumida para modelizar la variable que observamos, utilizando los datos observados.


En Bayesiano jugaremos siempre en términos de probabilidad. Solemos tener cierto conocimiento a priori sobre los parámetros que nos interesan, que expresaremos en términos de probabilidad. Modelizaremos con probabilidad los datos que vayamos a observar. Combinaremos la probabilidad previa con la información que dan los datos y actualizaremos la probabilidad posterior sobre el parámetro de interés. Con dicha probabilidad posterior extraeremos conclusiones y hablaremos de cuántas garantías les otorgamos.

## Probabilidad 

- ¿Qué entiendes por probabilidad?

$$\text{Probabilidad(A)}=\frac{\text{Casos favorables aconteciendo A}}{\text{Casos posibles/totales}}$$

- ¿Qué entiendes por probabilidad condicionada?


$$\text{Probabilidad (A|B)}=\frac{\text{Casos favorables aconteciendo A y B}}{\text{Casos posibles/totales aconteciendo B}}$$


El teorema de Bayes se formula a partir de la probabilidad condicional. Dados dos eventos, A y B, 
$$Pr(A|B)=\frac{Pr(A,B)}{Pr(B)}=\frac{Pr(B|A)\cdot Pr(A)}{Pr(B)}$$

que nos permite igualmente calcular las probabilidades conjuntas:
$$Pr(A,B)=Pr(A|B)\cdot Pr(B)=Pr(B|A)\cdot Pr(A).$$


::: {.exm #buques}
**Ejemplo**
Tenemos una tabla con el número de buques pesqueros que pescaron frente a las costas de Nueva Zelanda (básicamente pesqueros de arrastre de calamar), y de los que accidentalmente pescaron leones marinos (LM), protegidos y en peligro de extinción, durante las temporadas de pesca de 1987/88 a 1995/96 (\@ref{link2010}). Los buques están clasificados según su bandera.

|  | Japón | NZ | Rusia | Total|
| --- | --- | --- | --- | --- |
Pesca accidental LM | 1 | 6 | 23 | 30 |
Total | 19 | 96 | 123 | 238 |
Pr(país) |  |  |  |  |
"Pr(PA|país)" |  |  |  |  |
"Pr(país|PA)" |  |  |  |  |
Odds ratio |  |  |  |  |


```r
total = c(19,96,123)
paises = c("Japón", "NZ", "Rusia")
nacc = c(1,6,23)

datos = data.frame(paises, nacc, total)
datos
#>   paises nacc total
#> 1  Japón    1    19
#> 2     NZ    6    96
#> 3  Rusia   23   123
```


1. ¿Qué porcentaje de barcos pescan con cada una de las tres banderas? Si salimos un día cualquiera a navegar, ¿con qué probabilidad nos encontraremos un buque pesquero ruso? ¿Y neocelandés? ¿Y japonés?


```r
# probabilidades marginales país
datos$pais= round(datos$total/sum(datos$total),4)
datos
#>   paises nacc total   pais
#> 1  Japón    1    19 0.0798
#> 2     NZ    6    96 0.4034
#> 3  Rusia   23   123 0.5168
```

1. ¿Qué proporción de los barcos con bandera japonesa han realizado alguna captura accidental de leones marinos? ¿Y para el resto de banderas?

```r
# probabilidades condicionadas por país
datos$acc_pais= round(datos$nacc/datos$total,4)
datos
#>   paises nacc total   pais acc_pais
#> 1  Japón    1    19 0.0798   0.0526
#> 2     NZ    6    96 0.4034   0.0625
#> 3  Rusia   23   123 0.5168   0.1870
```

1. ¿Con qué probabilidad podríamos decir que se da una captura accidental de león marino en la costa neocelandesa?

```r
# probabilidades marginales pesca accidental
datos$acc= rep(round(sum(datos$nacc)/sum(datos$total),4),3)
datos
#>   paises nacc total   pais acc_pais    acc
#> 1  Japón    1    19 0.0798   0.0526 0.1261
#> 2     NZ    6    96 0.4034   0.0625 0.1261
#> 3  Rusia   23   123 0.5168   0.1870 0.1261
```

1. Si nos reportan una captura accidental de león marino, ¿con qué probabilidad te atreverías a vaticinar que el pesquero era ruso? ¿Y japonés? ¿Y neocelandés?

```r
# probabilidades condicionadas de pesca accidental por país
datos$pais_acc= round(datos$nacc/sum(datos$nacc),4)
datos
#>   paises nacc total   pais acc_pais    acc pais_acc
#> 1  Japón    1    19 0.0798   0.0526 0.1261   0.0333
#> 2     NZ    6    96 0.4034   0.0625 0.1261   0.2000
#> 3  Rusia   23   123 0.5168   0.1870 0.1261   0.7667
```
1. Si nos reportan una captura accidental de león marino, ¿qué ventaja tendremos si apostamos a que fue ruso? 

```r
# odds de pesca accidental por país
datos$odds.acc= round(datos$pais_acc/(1-datos$pais_acc),4)
datos
#>   paises nacc total   pais acc_pais    acc pais_acc
#> 1  Japón    1    19 0.0798   0.0526 0.1261   0.0333
#> 2     NZ    6    96 0.4034   0.0625 0.1261   0.2000
#> 3  Rusia   23   123 0.5168   0.1870 0.1261   0.7667
#>   odds.acc
#> 1   0.0344
#> 2   0.2500
#> 3   3.2863
```
1. Si salimos un día cualquiera a navegar, ¿con qué probabilidad nos encontraremos un buque pesquero ruso que ha realizado alguna captura accidental de león marino? ¿Y neocelandés? ¿Y japonés?

```r
# probabilidades conjuntas
datos$acc.pais = round(datos$nacc/sum(datos$total),4)
datos
#>   paises nacc total   pais acc_pais    acc pais_acc
#> 1  Japón    1    19 0.0798   0.0526 0.1261   0.0333
#> 2     NZ    6    96 0.4034   0.0625 0.1261   0.2000
#> 3  Rusia   23   123 0.5168   0.1870 0.1261   0.7667
#>   odds.acc acc.pais
#> 1   0.0344   0.0042
#> 2   0.2500   0.0252
#> 3   3.2863   0.0966
```


```r
# Toda la sintaxis reunida
total = c(19,96,123)
paises = c("Japón", "NZ", "Rusia")
nacc = c(1,6,23)

datos = data.frame(paises, nacc, total)
# probabilidades marginales país
datos$pais= round(datos$total/sum(datos$total),4)
# probabilidades condicionadas por país
datos$acc_pais= round(datos$nacc/datos$total,4)
# probabilidades marginales pesca accidental
datos$acc= rep(round(sum(datos$nacc)/sum(datos$total),4),3)
# probabilidades condicionadas de pesca accidental por país
datos$pais_acc= round(datos$nacc/sum(datos$nacc),4)
# odds de pesca accidental por país
datos$odds.acc= round(datos$pais_acc/(1-datos$pais_acc),4)
# probabilidades conjuntas
datos$acc.pais = round(datos$nacc/sum(datos$total),4)

kbl(datos) %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> paises </th>
   <th style="text-align:right;"> nacc </th>
   <th style="text-align:right;"> total </th>
   <th style="text-align:right;"> pais </th>
   <th style="text-align:right;"> acc_pais </th>
   <th style="text-align:right;"> acc </th>
   <th style="text-align:right;"> pais_acc </th>
   <th style="text-align:right;"> odds.acc </th>
   <th style="text-align:right;"> acc.pais </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Japón </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 0.0798 </td>
   <td style="text-align:right;"> 0.0526 </td>
   <td style="text-align:right;"> 0.1261 </td>
   <td style="text-align:right;"> 0.0333 </td>
   <td style="text-align:right;"> 0.0344 </td>
   <td style="text-align:right;"> 0.0042 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NZ </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 96 </td>
   <td style="text-align:right;"> 0.4034 </td>
   <td style="text-align:right;"> 0.0625 </td>
   <td style="text-align:right;"> 0.1261 </td>
   <td style="text-align:right;"> 0.2000 </td>
   <td style="text-align:right;"> 0.2500 </td>
   <td style="text-align:right;"> 0.0252 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rusia </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 123 </td>
   <td style="text-align:right;"> 0.5168 </td>
   <td style="text-align:right;"> 0.1870 </td>
   <td style="text-align:right;"> 0.1261 </td>
   <td style="text-align:right;"> 0.7667 </td>
   <td style="text-align:right;"> 3.2863 </td>
   <td style="text-align:right;"> 0.0966 </td>
  </tr>
</tbody>
</table>
:::


Repetimos los cálculos con el teorema de Bayes.

Probabilidad condicionada:
$$Pr(pais|acc)=\frac{Pr(pais,acc)}{Pr(acc)}=\frac{Pr(acc|pais) \cdot Pr(pais)}{Pr(acc)}$$
Probabilidad conjunta:
$$Pr(pais,acc)=Pr(acc|pais) \cdot Pr(pais)$$

```r
pais_acc=datos$acc_pais*datos$pais/datos$acc; pais_acc
#> [1] 0.03328692 0.19994052 0.76638858
acc.pais=datos$acc_pais*datos$pais; acc.pais
#> [1] 0.00419748 0.02521250 0.09664160
```

### Probabilidad y Variables aleatorias

Toda la información e incertidumbre que tenemos sobre una característica observable la podemos formular en términos de una distribución de probabilidad. Cuando hablamos de variables aleatorias, hablamos de características observables cuyas observaciones varían (son distintas) en la población en que se observan, y el modo en que varían o se distribuyen, se puede aproximar o describir matemáticamente con una distribución de probabilidad.

Una distribución de probabilidad viene descrita por una función matemática que permite calcular las probabilidades asociadas a la variable estadística que representa.

Básicamente podemos diferenciar entre variables discretas (categóricas, cualitativas) y continuas (numéricas, cuantitativas). Las primeras se identifican porque solo pueden tomar ciertos valores y las segundas pueden tomar cualquier valor entre dos valores cualesquiera.

**Ejemplo**

Discreta dicotómica. Se registra si un ave capturada está anillada.
Discreta entera. Se contabilizan el número de nidos en cada una de las parcelaciones de un terreno.
Continua.




### La aproximación frecuentista

