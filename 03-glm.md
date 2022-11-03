# Modelos lineales generalizados {#glm}

Los modelos lineales generalizados (Generalized Linear Models or GLM), son una clase de modelos introducidos por Nelder y Wedderburn (1972) y McCullagh y Nelder (1989), con el objetivo de extender la regresión lineal al caso en el que la variable dependiente no se distribuye necesariamente según una normal, pero su distribución todavía pertenece a la familia exponencial. Trabajamos a continuación con dos de los GLM más comunes en epidemiología y ciencias sociales: la regresión logística y la de Poisson, mostrando cómo usar R-INLA.

Un modelo lineal general (MLG) está basado en asumir una relación lineal entre cierta transformación del valor esperado de la respuesta $E(y_i)$ y los predictores disponibles, sean covariables, efectos fijos, efectos aleatorios, o incluso alguna función de estos.



Si $y$ representa una respuesta observada, $x_1,x_2,..$ una serie de covariables o efectos fijos, y $z_1,z_2,...$ efectos aleatorios, el modelo lineal general asume cierto modelo probabilístico para $y$ dentro de la familia Gausiana, dependiente de los predictores y ciertos parámetros desconocidos $\theta$, que se especifica en un primer nivel de la modelización.

$$E(y_i|x,z,\theta)=\gamma_i=g(\eta_i)$$ 

$\theta$ engloba todos los parámetros de la distribución de los datos, esto es, los involucrados en el predictor lineal, pero también otros como por ejemplo la varianza $\sigma^2$ del modelo normal.


El predictor lineal $\eta$ está relacionado linealmente con los predictores según:
$$\eta_i=\mu + \beta_1 x_{1i} +   \beta_2 x_{2i} +...+ z_{1i} + z_{2i}+...$$
relación que se suele representar en forma matricial como 
$$\eta=X\beta + Z u$$
Los modelos lineales que conocemos (regresión, anova, ancova, glm,modelos mixtos, aditivos,...) se engloban dentro del modelo lineal general.

Para obtener la distribución marginal de los valores ajustados y predichos necesitamos incorporar  a la función `inla` el argumento `control.compute=list(return.marginals.predictor=TRUE)`
El argumento `control.predictor=list(compute=TRUE)` en la función `inla` permite obtener las distribuciones predictivas para el predictor lineal, que en estos casos es distinto a los valores ajustados, `fit$summary.fitted`. Además:

-   `fit$summary.linear.predictor` resume  la inferencia posterior sobre los
    predictores lineales (distintos a los fitted cuando hay una función *link*)
-   `fit$marginals.linear.predictor` da las distribuciones posteriores marginales para los predictores lineales



### Modelos jerárquicos bayesianos
Un modelo bayesiano se modeliza a través de un modelo jerárquico en el que en el nivel 1 se define la distribución asumida sobre la variable respuesta y que determina la verosimilitud. Esta variable depende de unos parámetros que definen los efectos fijos y aleatorios, y para los que hay que proporcionar la información previa disponible a través de una distribución a priori en el segundo nivel del modelo jerárquico. La distribución a priori para los efectos fijos generalmente será común a todos ellos, mientras que la distribución a priori para los efectos aleatorios estará vinculada a otros hiperparámetros para los que también será preciso especificar una distribución a priori en un tercer nivel de la modelización.

El argumento `control.predictor=list(compute=TRUE)` en la función `inla` permite obtener las distribuciones predictivas para el predictor lineal. Así para el modelo anterior tendremos

\begin{eqnarray*}
(I)& y | X, Z, \theta &\sim f(y|x,z,\theta) \\
&& E(y_i|x,z,\theta)=\gamma_i=g(\eta_i) \\
&& \eta=X\beta + Z u \\
(II)&  \beta &\sim N(0,\sigma_{\beta}), \text{ con un valor dado  para } \sigma_{\beta} \\
&& u_i|\mu,\sigma \sim_{iid}  N(\mu,{\sigma}) \\
(III)& \mu  &\sim N(0,\sigma_{\mu}), \text{ con un valor dado  para } \sigma_{\mu} \\
&& \sigma \sim  GaI(\alpha_{\sigma},\beta_{\sigma}), \text{ con un valor dado  para }\alpha_{\sigma} y \beta_{\sigma}.
\end{eqnarray*}


## Regresión logística
La regresión logística es el modelo estándar para respuestas binarias (éxitos/fracasos). Tiene dos variaciones, en función de si la respuesta representa observaciones individuales (0/1) o conteos (de éxitos) en grupos de sujetos.
Si las observaciones son individualizadas, entonces
$$y_1,...,y_n \sim Ber(\pi_i), \ i=1,...,n$$
En el caso de que sean conteos en grupos, 
$$y_1,...,y_n \sim Bin(n_i,\pi_i), \ i=1,...,n$$
siendo $n_i$ el tamaño de cada uno de los $n$ grupos disponibles, y $\pi_i$ la probabilidad de éxito (output de interés).

La relación entre el predictor lineal $\eta_i$ construido con los predictores disponibles $x_i=(x_{i1},...x_{iM})$ y la probabilidad $\pi_i$ se especifica a través de la función *logit*:
$$\eta_i=logit(\pi_i)=log\left(\frac{\pi_i}{1-\pi_i}\right)=x_i\beta=\beta_0+\sum_{m=1}^M \beta_mx_{im}$$
de forma que
$$\pi_i=logit^{-1}(x_i\beta)=\frac{exp(x_i\beta)}{1-exp(x_i\beta)}$$
Una vez especificado el modelo, si no hay información previa disponible sobre los efectos (fijos) ${\beta_o,\beta_1,...\beta_M}$, se asumen distribuciones a priori independientes y  normales con media cero y varianza muy grande. 


### Interpretación de los coeficientes en la regresión logit

Puesto que $x_i\beta=\beta_0+\sum_{m=1}^M \beta_mx_{im}$, 
- cuando todas las $x$ son categóricas y su valor coincide con el nivel de referencias
- o cuando las $x$ son continuas y su valor es cero, 
el valor del predictor lineal $\eta_i=\beta_0=logit(\pi_i)$. En consecuencia, el logit inverso de $\beta_0$ se interpreta como la probabilidad de éxito $\pi_i$ cuando los predictores están en su nivel de referencia o son cero.
$$logit^{-1}(\beta_0)=Pr(y=1|x=0).$$ 

En cuanto a la interpretación de cualquier coeficiente de regresión, como $\beta_1$, echamos mano del concepto de odds y odds ratio. Los odds ratio, OR, representan el cociente de las posibilidades a favor de un evento $E$ bajo condiciones A y de las posibilidades del mismo evento bajo condiciones B. Nos sirve para evaluar cuánto afecta a dicho evento el hecho de variar las condiciones de B a A.
$$OR(A,B)=\frac{Pr(E|A)/(1-Pr(E|A))}{Pr(E|B)/(1-Pr(E|B))}.$$

En el modelo logístico, nos interesa saber el efecto que tiene sobre la respuesta (sobre las probabilidad de éxito en este caso) el incremento de una unidad en la variable $x$, y para ello consideramos los odds bajo $x$ y los odds bajo $x+1$, y en particular el logaritmo de los odds, log-odds:
$$log.odds(x+1)=\frac{P(y=1|x+1)}{P(y=0|x+1)}=\beta_0+\beta_1 (x+1)$$
$$log.odds(x)=\frac{P(y=1|x)}{P(y=0|x)}=\beta_0+\beta_1 x$$
de modo que 
$$log.ods (x+1)/x = \beta_1 \rightarrow exp(\beta_1)=\frac{odds(x+1)}{odds(x)}=OR(x+1,x).$$
Así $\beta_1$ representa el cambio en los odds a favor de un éxito cuando se incrementa en una unidad el predictor $x$ al que acompaña en el predictor lineal. Esta interpretación es muy común en Epidemiología.


## Ejemplo modelo logit. Mortalidad por infarto en Sheffield.

Utilizamos los datos *stroke*, disponibles en [datasets in SSTM-RINLA](https://sites.google.com/a/r-inla.org/stbook/datasets). Queremos evaluar la presencia de cierta asociación entre los niveles de NOx y el infarto en Sheffield, UK. Se dispone de la concentración anual de NOx medida en $\mu g/m^3$ y categorizada en quintiles, promediada durante el periodo 1994-1999, en la variable *NOx.class* y su análoga $NOx$, y el número de muertes por infarto  $y$ en cada distrito identificado por el índice de desventajas y privación *Townsend* (categorizado en quintiles). Se dispone igualmente del tamaño de la población para cada registro, en la variable $pop$. La respuesta $y$ se puede considerar entonces como conteos (de muertes) sobre la población de cada distrito,
$$y_i|\pi_i \sim Bin(n_i, \pi_i)$$
y el predictor lineal es función del nivel de NOx y del distrito:
$$\eta_i=logit(\pi_i)=\beta_0 + \sum_{k=2}^5 \beta_{1k} I(NOx_i=k)+ \sum_{h=2}^5 \beta_{2h} I(Townsend_i=h)+logit(\tilde{p_i})$$
siendo $n_i$ la población (número total de habitantes) del distrito en el que se ubica el registro $i$ (disponible en la variable $pop$). El término $\tilde{p_i}$ representa el resgo ajustado por sexo y edad de la mortalidad por infarto, calculada utilizando estandarización indirecta con ratios de referencia internos basados en 18 estratos (9 para edad y 2 para género), y que se usa como un riesgo base en el modelo (Maheswaran et al.2006). En el ejemplo se calcula como el ratio de la mortalidad dividido por la población de cada registro.

```r
my.dir="~/Dropbox/ESTADISTICA/BAYESIAN/VARIOS/"
Stroke <- read.csv(paste0(my.dir,"Stroke.csv"),sep=",",dec=".",header=TRUE)
#riesgo base: ajuste por tamaño de la población
Stroke$Adjusted.prob <- Stroke$stroke_exp/Stroke$pop
# logit del riesgo base
Stroke$logit.adjusted.prob <- log(Stroke$Adjusted.prob/(1-Stroke$Adjusted.prob))                          
```

Ajustamos ya el modelo, en el que las variables NOx y Townsend actúan como factores (efectos fijos) y el riesgo base se introduce como offset, para estandarizar los riesgos en función del tamaño de la población y poder equiparar así todos los distritos:

```r
formula.inla <- y ~ 1 + factor(NOx) + factor(Townsend) + offset(logit.adjusted.prob)
model.logistic <- inla(formula.inla, family="binomial", Ntrials=pop, data=Stroke)
round(model.logistic$summary.fixed[,1:5],3)
#>                     mean    sd 0.025quant 0.5quant
#> (Intercept)       -0.181 0.057     -0.293   -0.180
#> factor(NOx)2       0.132 0.059      0.016    0.132
#> factor(NOx)3       0.105 0.061     -0.014    0.105
#> factor(NOx)4       0.261 0.059      0.144    0.260
#> factor(NOx)5       0.425 0.062      0.302    0.425
#> factor(Townsend)2  0.077 0.061     -0.043    0.077
#> factor(Townsend)3  0.137 0.060      0.020    0.137
#> factor(Townsend)4 -0.132 0.063     -0.255   -0.132
#> factor(Townsend)5 -0.118 0.067     -0.250   -0.118
#>                   0.975quant
#> (Intercept)           -0.071
#> factor(NOx)2           0.248
#> factor(NOx)3           0.225
#> factor(NOx)4           0.377
#> factor(NOx)5           0.547
#> factor(Townsend)2      0.198
#> factor(Townsend)3      0.255
#> factor(Townsend)4     -0.009
#> factor(Townsend)5      0.014
```

Para obtener la probabilidad promedio de muerte por infarto en el distrito Townsend=1 y para el nivel NOx=1, que son los niveles base, nos apoyamos en la interceptación $\beta_0$. Sobre sus simulaciones será preciso deshacer el logit (con la función logit-inversa):

```r
prob.stroke <- inla.tmarginal(function(x) exp(x)/(1+exp(x)), model.logistic$marginals.fixed[[1]])
inla.zmarginal(prob.stroke)
#> Mean            0.455034 
#> Stdev           0.0139168 
#> Quantile  0.025 0.427486 
#> Quantile  0.25  0.445615 
#> Quantile  0.5   0.455082 
#> Quantile  0.75  0.464485 
#> Quantile  0.975 0.482137
ggplot(data.frame(prob.stroke),aes(x=x,y=y))+geom_line()+
         labs(x=expression(pi),y= expression(tilde(p)(paste(pi,"|",y,",",NOx[1],",",TS[1]))))
```

![](03-glm_files/figure-latex/unnamed-chunk-4-1.pdf)<!-- --> 

El efecto $\beta_{12}$ representa el efecto en los log.odds de la mortalidad por infarto de estar en el nivel $NOx=2$ frente al de estar en el nivel $NOx=1$. Si queremos evaluar el odds-ratio, simplemente calculamos la distribución posterior de $exp(\beta_{12})$ con `inla.tmarginal`. Si sólo estamos interesados en su media, bastaría utilizar `inla.emarginal':

```r
odds.nox21 <- inla.tmarginal(function(x) exp(x), model.logistic$marginals.fixed$"factor(NOx)2")
e<-inla.emarginal(exp, model.logistic$marginals.fixed$"factor(NOx)2")
ggplot(data.frame(odds.nox21),aes(x=x,y=y))+geom_line()+
        labs(x=expression(OR(NOx[21])),y= expression(tilde(p)(paste(OR(NOx[21]),"|",y))))+
        geom_vline(xintercept=e,color="pink")+ geom_text(x=e,y=1,label=paste("mean=",round(e,3)),color="red")
```

![](03-glm_files/figure-latex/unnamed-chunk-5-1.pdf)<!-- --> 
La probabilidad de muerte por infarto se incrementa en un 14.3% cuando la exposición de NOx cambia del primer al segundo nivel. Los odds-ratio para el resto de los niveles resultan realmente significativos: 11,3% el odds-ratio NOx3/NOx1, 30% NOx4/NOx1 y 53.2% NOx5/NOx1.

```r
inla.emarginal(exp, model.logistic$marginals.fixed$"factor(NOx)3")
#> [1] 1.113128
inla.emarginal(exp, model.logistic$marginals.fixed$"factor(NOx)4")
#> [1] 1.299909
inla.emarginal(exp, model.logistic$marginals.fixed$"factor(NOx)5")
#> [1] 1.532265
```


## Regresión de Poisson

La regresión de Poisson es útil cuando la variable respuesta representa conteos y estos toman valores discretos entre 0 y $+\infty$, sin una cota superior de referencia. El parámetro de interés es el número promedio de eventos $\lambda_i=E(y_i)$ y el link natural es el logaritmo, de modo que el predictor lineal está ligado con las covariables y factores según:
$$\eta_i=log(\lambda_i)=x_i \beta, \ \ \mbox{ y } \ \ \lambda_i=exp(x_i \beta)$$ 

Un modelo de Poisson puede especificarse según
$$y_i \sim Po(\lambda_i), i=1,...,n$$
$$\eta_i=log(\lambda_i)=\beta_0+\sum_{m=1}^M \beta_mx_{im}.$$

Para completar el modelo se especifican distribuciones a priori para $\beta$, típicamente como normales con media cero y una varianza grande cuando no hay información disponible de estudios previos u opinión de expertos.

Los coeficientes se interpretan a través de la función exponencial:

* $exp(\beta_0)=\lambda_i$ cuando todas las $x=0$ si son continuas, o para el primer nivel de las categorías posibles si son categóricas.
* $exp(\beta_m)$ es el cambio que se produce en la respuesta promedio $y$ cuando $x_m$ se incrementa en una unidad.

La mayoría de las veces que se utiliza la regresión de Poisson, el interés recáe en las ratios o riesgos relativos, más que en el número promedio de casos $\lambda_i$. Para cambiar la escala en términos de riesgo, ha de utilizarse un offset como factor de corrección en la especificación del modelo. Este offset representa el denominador del riesgo y entra en la regresión en una escala logarítmica, asumiendo que tiene un coeficiente de regresión fijado a 1:
$$\eta_i=log(\lambda_i)=\beta_0+\sum_{m=1}^M \beta_mx_{im}+log(Offset_i)$$
donde el riesgo relativo de que se produzca un evento se obtiene según
$$log\left(\frac{\lambda_i}{Offset_i}\right)=\beta_0+\sum_{m=1}^M \beta_mx_{im}$$
y los coeficientes entonces se interpretan en una escala de riesgo. En este caso al exponenciar la interceptación obtenemos el riesgo base, mientras que $exp(\beta_m)$ representa el cambio en el riesgo relativo debido a un cambio de unidad en el predictor correspondiente.

### Ejemplo modelo Poisson: incidentes en barcos
Utilizamos los datos *ships.csv* en [datasets for SSTM-RINLA](https://sites.google.com/a/r-inla.org/stbook/datasets) para estimar el riesgo mensual de incidentes en barcos. Los factores potenciales del riesgo son el periodo de construcción (*built*), el periodo de operación (*oper*) y el tipo de barco (*type*).
El modelo se escribe en INLA a continuación, utilizando como offset el *log(months)*, que son los meses que ha navegado y ponderan en consecuencia el riesgo de incidentes. El modelo con el offset será entonces 
$$y_i \sim Poisson (E_i \rho_i),$$
donde $\eta_i=log(\rho_i)$ es el predictor lineal y el promedio del número de incidentes $\lambda_i=E_i \rho_i$. El offset no se incluye en esta formulación en el predictor lineal.

```r
my.dir="~/Dropbox/ESTADISTICA/BAYESIAN/VARIOS/"
ShipsIncidents <- read.csv(file=paste0(my.dir,"Ships.csv"),sep=",") 

formula.inla <- y ~ 1 + built + oper + type
model.poisson <- inla(formula.inla,family="poisson", data=ShipsIncidents, offset=log(months))

round(model.poisson$summary.fixed[,1:5],3)
#>               mean    sd 0.025quant 0.5quant 0.975quant
#> (Intercept) -6.416 0.217     -6.852   -6.413     -5.998
#> built65-69   0.696 0.150      0.406    0.695      0.994
#> built70-74   0.819 0.170      0.487    0.818      1.153
#> built75-79   0.453 0.233     -0.012    0.455      0.903
#> oper75-79    0.384 0.118      0.153    0.384      0.617
#> typeB       -0.543 0.178     -0.882   -0.546     -0.185
#> typeC       -0.688 0.329     -1.366   -0.678     -0.075
#> typeD       -0.075 0.291     -0.664   -0.069      0.476
#> typeE        0.326 0.236     -0.141    0.327      0.785
```


```r
names(model.poisson$marginals.fixed)
#> [1] "(Intercept)" "built65-69"  "built70-74"  "built75-79" 
#> [5] "oper75-79"   "typeB"       "typeC"       "typeD"      
#> [9] "typeE"
# ratio medio de incidentes por mes en las categorías base
inla.emarginal(exp,model.poisson$marginals.fixed[[1]])
#> [1] 0.001674164
# riesgo relativo de barcos tipo E
inla.emarginal(exp,model.poisson$marginals.fixed$typeE)
#> [1] 1.424414
```
Así, la media de $exp(\beta_0)$, 0.0018, representa el ratio medio de incidentes por mes entre los barcos que fueron construidos entre el 60 y el 64, han operado entre el 60 y el 74 y son de tipo A (las categorías de referencia). El ratio en 1000 meses sería del 1.8.Para los barcos de tipo E el incremento en el ratio mensual de incidentes, comparado con los de tipo A es del 42,4%.



Otros datos modelizables con una regresión de Poisson son los que provienen del libro de Andrews and Herzberg y están descritos en  ([randomservices.org/random](http://www.randomservices.org/random/data/HorseKicks.html)) consistentes en el número de soldados muertos por coces de caballo en diversos cuerpos de caballería del ejército prusiano, entre 1875 y 1894.

Este modelo se implementa en INLA, a partir de datos simulados, con el siguiente código

```r
my.dir="~/Dropbox/ESTADISTICA/BAYESIAN/VARIOS/"
horse<-read.csv(paste0(my.dir,"HorseKicks.txt"),sep="", dec=".",header=TRUE)
horse$sum<-apply(horse[,2:ncol(horse)],1,sum)

fit=inla(sum~1,data=horse,family="poisson",control.predictor=list(compute=TRUE))
summary(fit)
#> 
#> Call:
#>    c("inla.core(formula = formula, family = family, 
#>    contrasts = contrasts, ", " data = data, quantiles = 
#>    quantiles, E = E, offset = offset, ", " scale = 
#>    scale, weights = weights, Ntrials = Ntrials, strata = 
#>    strata, ", " lp.scale = lp.scale, link.covariates = 
#>    link.covariates, verbose = verbose, ", " lincomb = 
#>    lincomb, selection = selection, control.compute = 
#>    control.compute, ", " control.predictor = 
#>    control.predictor, control.family = control.family, 
#>    ", " control.inla = control.inla, control.fixed = 
#>    control.fixed, ", " control.mode = control.mode, 
#>    control.expert = control.expert, ", " control.hazard 
#>    = control.hazard, control.lincomb = control.lincomb, 
#>    ", " control.update = control.update, 
#>    control.lp.scale = control.lp.scale, ", " 
#>    control.pardiso = control.pardiso, only.hyperparam = 
#>    only.hyperparam, ", " inla.call = inla.call, inla.arg 
#>    = inla.arg, num.threads = num.threads, ", " 
#>    blas.num.threads = blas.num.threads, keep = keep, 
#>    working.directory = working.directory, ", " silent = 
#>    silent, inla.mode = inla.mode, safe = FALSE, debug = 
#>    debug, ", " .parent.frame = .parent.frame)") 
#> Time used:
#>     Pre = 2.08, Running = 0.135, Post = 0.00666, Total = 2.22 
#> Fixed effects:
#>              mean    sd 0.025quant 0.5quant 0.975quant mode
#> (Intercept) 2.282 0.071      2.139    2.283       2.42   NA
#>             kld
#> (Intercept)   0
#> 
#> Marginal log-Likelihood:  -61.30 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
```
