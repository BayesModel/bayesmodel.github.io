# Regresión lineal {#inlabasics}



## Introducción

Una vez presentados los fundamentos de INLA vamos a utilizarlo para trabajar progresivamente desde los modelos más sencillos a los más sofisticados. Empezamos aquí con el modelo de regresión lineal simple, para continuar generalizando con el de regresión lineal múltiple.

Vamos a empezar a trabajar con INLA modelizando y ajustando un problema
sencillo de regresión lineal simple. La base de datos Davis (en la librería `carData`) contiene 200 registros de 5 variables relacionadas con un estudio sobre habituación de hombres y mujeres a la realización de ejercicio físico de forma regular. Las variables que se registraron son sexo, peso y altura (reales y reportados). Vamos a indagar la relación entre el peso real y el reportado a través de un análisis de regresión lineal simple.


```r
library("carData")
summary(Davis)
#>  sex         weight          height          repwt       
#>  F:112   Min.   : 39.0   Min.   : 57.0   Min.   : 41.00  
#>  M: 88   1st Qu.: 55.0   1st Qu.:164.0   1st Qu.: 55.00  
#>          Median : 63.0   Median :169.5   Median : 63.00  
#>          Mean   : 65.8   Mean   :170.0   Mean   : 65.62  
#>          3rd Qu.: 74.0   3rd Qu.:177.2   3rd Qu.: 73.50  
#>          Max.   :166.0   Max.   :197.0   Max.   :124.00  
#>                                          NA's   :17      
#>      repht      
#>  Min.   :148.0  
#>  1st Qu.:160.5  
#>  Median :168.0  
#>  Mean   :168.5  
#>  3rd Qu.:175.0  
#>  Max.   :200.0  
#>  NA's   :17
plot(weight ~ repwt, data=Davis)
```

![](01-regresion_files/figure-latex/unnamed-chunk-2-1.pdf)<!-- --> 

```r
# Excluimos el valor outlier y los NA
davis=Davis %>%
  filter(weight<160) %>%
  slice(which(!is.na(repwt)))
```


## Variables y relaciones

Entendemos como variable respuesta $y=weight$, de tipo numérico
(continua), y como variable explicativa o covariable, `repwt`, que pretendemos relacionar con un modelo de regresión lineal. 

La especificación de respuesta $y$ y predictores $x_1,x_2,...$, en INLA se registra en una fórmula del tipo:


```r
formula= y ~ 1 + x_1 + x_2 +...
```

en la que podemos prescindir del $1$ que identifica la interceptación, pues el ajuste, por defecto, siempre se resolverá con su estimación, salvo que en su lugar se escriba un '-1'.

En nuestro problema tendríamos pues,


```r
formula = weight ~ repwt
```

A continuación es procedente elegir el modelo sobre la respuesta, o lo
que es lo mismo, la verosimilitud.

## Verosimilitud

En principio es razonable asumir normalidad en la respuesta, además de independencia entre todas las observaciones. Así, el modelo propuesto para la respuesta es:

$$(y_i|\theta,\beta,\sigma^2) \sim N(\theta+\beta x_i,\sigma^2), i=1,...,n$$ 

donde $y=weight$ y $x=repwt$, e $i=1,...,n$ es un subíndice que identifica a cada uno de los registros disponibles en el banco de datos. $\theta$ es la interceptación de la recta de regresión y $\beta$ el coeficiente que explica la relación lineal entre $s$ e $y$. Veamos cómo especificar este modelo con INLA.


En INLA el vector $(\theta,\beta)$ identifica los **efectos latentes**, en cuya inferencia posterior estamos interesados. El modelo, o lo que es lo mismo, la verosimilitud, depende también de un parámetro de dispersión $\sigma^2$ sobre el que también querremos inferir. 


La función `name(inla.models())` proporciona un listado de todos los
tipos de modelos posibles para los datos (*likelihood*), parámetros
(*prior*), hiperparámetros (*latent*), ... Es más, el listado completo
de todas las distribuciones disponibles para cada uno de los tipos de
modelos lo obtenemos con el comando `inla.list.models()`

En particular, si ejecutamos `names(inla.models()$likelihood)` obtenemos
todas las distribuciones disponibles para modelizar los datos. La
distribución `gaussian` identifica la distribución normal que buscamos para resolver nuestro problema. Para obtener información sobre cómo está parametrizada y cuáles son los
hiperparámetros por defecto, basta consultar la documentación *Gaussian* con el comando:


```r
# documentación (parametrización y valores por defecto)
inla.doc("gaussian")
```

Para ajustar un modelo sencillo en INLA hay que echar mano de la función `inla`, en la que introducimos en primer lugar la `formula`, con la relación entre las variables, a continuación el modelo, en el argumento `family`, y después el banco de datos sobre el que trabajamos. Si no especificamos el argumento `family`, la función `inla`interpreta por defecto la opción `gaussian`, esto es, normalidad para los datos, de modo que podríamos excluir dicha especificación cuando modelizamos datos normales.


```r
# Asumiendo datos normales
fit=inla(formula,family="gaussian", data)
# equivalente a
fit=inla(formula, data)
```

Adelantamos pues un paso más en nuestro problema, añadiendo la verosimilitud normal y la base de datos. Haciendo acopio de lo que hemos visto antes, bastaría con ejecutar:

```r
# ajuste del modelo
formula = weight ~ 1+ repwt
fit=inla(formula,family="gaussian",data=davis)
```


## Parámetros e hiperparámetros

Por defecto INLA asume unas distribuciones difusas sobre los efectos latentes $\theta,\beta$ y la varianza $\sigma^2$. Veamos cuáles son esas distribuciones, y cómo modificarlas si tenemos información previa.

Pensando en el significado de los parámetros (reconocidos como *hyperpar* o hiperparámetros) que aparecen en la verosimilitud, es lógico asumir, ante ausencia de información:

- para los efectos latentes $\theta,\beta$ distribuciones normales independientes difusas y centradas en el cero 
- para la varianza $\sigma^2$ es habitual asumir una gamma inversa, también difusa (con media y varianza grandes) si no tenemos información previa.

Vamos en primer lugar con el parámetro $\sigma^2$ en la verosimilitud. 

En INLA, en lugar de asignar distribuciones a priori sobre las
varianzas, se hace sobre el logaritmo de las precisiones, para facilitar
el cálculo del máximo de la log-posterior (obtenida de la
log-verosimilitud y la log-prior). Así, asumir una gamma inversa difusa
$GaI(\alpha,\beta)$ para la varianza es equivalente a una Gamma difusa
$Ga(\alpha,\beta)$ para la precisión $\tau=1/\sigma^2$, y una
log-gamma difusa $Log-Gamma(\alpha,\beta)$ para la log-precisión
$log(\tau)$.

Por defecto, como ya verificamos en la documentación de la verosimilitud gausiana, (con `inla.doc("gaussian")`), la distribución a priori por defecto sobre la log-precisión ($\theta_1$ en la ayuda) es la log-gamma difusa $LogGa(1,5\cdot 10^{-5})$, lo que da un valor esperado para la precisión $\tau$ de 20000 y una varianza de $4\cdot 10^8$. 

Para modificar la distribución a priori de un parámetro podemos utilizar cualquiera de las distribuciones que ofrece INLA en su listado `names(inla.models()$prior)` (siempre buscando coherencia con la información sobre el parámetro en cuestión). 

Para definir una prior para los parámetros o hiperparámetros en INLA hay
que definir los siguientes argumentos:

-   prior, el nombre de la distribución a priori (alguna de las opciones de `names(inla.models()$prior)`)
-   param, los valores de los parámetros de la prior
-   initial, el valor inicial para el hiperparámetro
-   fixed, variable booleana para decir si el hiperparámetro es fijo o aleatorio.

La modificación la haremos con el argumento `control.family=list(hyper=list(...))` en la función `inla`, al que le proporcionaremos una lista con el nombre de los parámetros (el *short name* con el que los identifica INLA en la documentación), que apunta a una lista con la distribución (*prior*) y los parámetros (*param*) a utilizar.

En nuestro problema, si tenemos información previa sobre la precisión del modelo,
podemos modificarla con la sintaxis:


```r
prec.info = list(prior="loggamma", param =c(1,0.001))
fit2=inla(formula,family="gaussian",data=davis,
      control.family = list(hyper = list(prec = prec.info)))
```

Si nos conformamos con la previa por defecto de INLA, el modelo que estamos asumiendo en nuestro problema de regresión lineal simple será:

\begin{eqnarray*}
(y_i|\theta,\beta,\sigma^2) &\sim & N(\theta+\beta x_i,\sigma^2), i=1,...,n \\
\tau=1/\sigma^2 & \sim & Ga(1,0.00005)
\end{eqnarray*}


## Efectos fijos

Una variable explicativa entra en el modelo como **efecto fijo** cuando
se piensa que afecta a todas las observaciones del mismo modo (de un
modo lineal), y que su efecto es de interés primario en el estudio. 

En un modelo de regresión lineal todos los efectos latentes,  interceptación y coeficientes de covariables, son efectos fijos. Por defecto en INLA y ante ausencia de información, las priors para los efectos fijos,  esto es, ($\theta,\beta$ en nuestro caso), son gausianas centradas en cero y con varianzas grandes. Así la interceptación $\beta$ tiene una prior gausiana con media y precisión igual a cero (una distribución plana objetiva, que no integra 1), y los coeficientes $\beta$ también tienen una prior gausiana con media cero y precisión igual a 0.001. Estos valores por defecto se pueden consultar con el comando `inla.set.control.fixed.default()`, que da la siguiente información:

- *mean=0* y *mean.intercept=0* son las medias de la distribución normal para los coeficientes $\beta$ y la interceptación $\theta$ respectivamente
- *prec=0.001* y *prec.intercept=0* son las precisiones respectivas de las normales para $\beta$ y $\theta$.

Con todo, los parámetros de las priors sobre los efectos fijos se pueden modificar a través del argumento `control.fixed=list(...)` en la función `inla`, utilizando siempre los nombres que atribuye INLA a los diferentes parámetros e hiperparámetros (*short name*).


```r
formula = weight ~ 1+ repwt
fit=inla(formula,family="gaussian",data=davis,
         control.fixed=list(prec=0.0001,prec.intercept=0.01))
```

Volvemos sobre nuestro ejemplo, y asumiendo las a priori por defecto de INLA tendremos:
\begin{eqnarray*}
(y_i|\theta,\beta,\sigma^2) & \sim & N(\theta+\beta x_i,\sigma^2), i=1,...,n \\
\theta & \sim & N(0,\infty) \\
\beta & \sim & N(0,1000) \\
\tau=1/\sigma^2 & \sim & Ga(1,0.00005)
\end{eqnarray*}


## Resultados

Para controlar qué resultados ha de mostrar INLA tras realizar un ajuste, hemos de especificar los ajustes del argumento `control.results` de la función `inla`.


Para mostrar una descriptiva de los resultados del ajuste obtenido con `fit=inla(...)`,
utilizamos la sintaxis siguiente:

-   `summary(fit)` proporciona una descriptiva del ajuste
-  `names(fit$marginals.fixed)` lista los nombres de todos los efectos fijos
-   `fit$summary.fixed` resume  la inferencia posterior sobre los efectos fijos
- `names(fit$marginals.hyperpar)` lista los nombres de todos los hiperparámetros
-   `fit$summary.hyperpar` da un resumen de la inferencia posterior de
    los parámetros e hiperparámetros
-   `fit$summary.fitted.values` resume  la inferencia posterior sobre los valores ajustados
-   `fit$mlik` da la estimación de la log-verosimilitud marginal, útil para evaluar y comparar modelos.
  
  
Veamos los resultados para nuestro problema de regresión.


```r
# ajuste del modelo 
formula = weight ~ 1+ repwt
fit=inla(formula,family="gaussian",data=davis)
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
#>     Pre = 2.35, Running = 0.512, Post = 0.018, Total = 2.88 
#> Fixed effects:
#>              mean    sd 0.025quant 0.5quant 0.975quant mode
#> (Intercept) 2.734 0.814      1.135    2.734      4.333   NA
#> repwt       0.958 0.012      0.935    0.958      0.982   NA
#>             kld
#> (Intercept)   0
#> repwt         0
#> 
#> Model hyperparameters:
#>                                          mean    sd
#> Precision for the Gaussian observations 0.199 0.021
#>                                         0.025quant 0.5quant
#> Precision for the Gaussian observations      0.161    0.198
#>                                         0.975quant mode
#> Precision for the Gaussian observations      0.242   NA
#> 
#> Marginal log-Likelihood:  -426.73 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
fit$mlik
#>                                            [,1]
#> log marginal-likelihood (integration) -426.7325
#> log marginal-likelihood (Gaussian)    -426.7339
```

Obtenemos pues la salida con un descriptivo de la distribución posterior para los efectos fijos interceptación, $\theta$, y el coeficiente del regresor `repwt`, $\beta$, con la media, desviación típica y cuantiles con los que podemos evaluar la región creíble al 95%.

También muestra a continuación una tabla con los descriptivos de la distribución posterior para la precisión $\tau=1/\sigma^2$ de los datos.

Si queremos obtener los nombres con los que INLA reconoce los efectos fijos (interceptación y coeficiente del regresor) e hiperparámetros (precisión de los datos), llamamos a 


```r
names(fit$marginals.fixed)
#> [1] "(Intercept)" "repwt"
names(fit$marginals.hyperpar)
#> [1] "Precision for the Gaussian observations"
```

Y para describir la marginal posterior sobre cada uno de los datos ajustados:


```r
head(fit$summary.fitted.values)
#>                          mean        sd 0.025quant 0.5quant
#> fitted.Predictor.001 76.52862 0.2162777   76.10404 76.52862
#> fitted.Predictor.002 51.61089 0.2441583   51.13157 51.61089
#> fitted.Predictor.003 54.48601 0.2190148   54.05606 54.48601
#> fitted.Predictor.004 69.82000 0.1750418   69.47637 69.82000
#> fitted.Predictor.005 59.27789 0.1856080   58.91351 59.27789
#> fitted.Predictor.006 75.57025 0.2087749   75.16039 75.57025
#>                      0.975quant mode
#> fitted.Predictor.001   76.95321   NA
#> fitted.Predictor.002   52.09021   NA
#> fitted.Predictor.003   54.91597   NA
#> fitted.Predictor.004   70.16364   NA
#> fitted.Predictor.005   59.64226   NA
#> fitted.Predictor.006   75.98011   NA
```

Si queremos hacer un **análisis de sensibilidad** sobre las distribuciones a priori, reajustamos el modelo con otras priors y comparamos los resultados.


```r
# ajuste del modelo 
formula = weight ~ 1+ repwt
fit2=inla(formula,family="gaussian",data=davis,
         control.fixed = list(mean.intercept = 100, 
                                 prec.intercept = 10^(-2),
                                 mean = 0, 
                                 prec = 1), 
            control.family = list(hyper = list(
              prec = list(prior="loggamma", param =c(1,0.001)))))
fit2$summary.fixed
#>                  mean         sd 0.025quant  0.5quant
#> (Intercept) 3.3859342 0.81577662  1.7941226 3.3824734
#> repwt       0.9488563 0.01215735  0.9248445 0.9489068
#>             0.975quant mode          kld
#> (Intercept)  4.9973585   NA 2.204325e-09
#> repwt        0.9725817   NA 2.141045e-09
fit2$summary.hyperpar
#>                                              mean
#> Precision for the Gaussian observations 0.1983981
#>                                                 sd
#> Precision for the Gaussian observations 0.02084352
#>                                         0.025quant
#> Precision for the Gaussian observations  0.1602208
#>                                          0.5quant
#> Precision for the Gaussian observations 0.1976839
#>                                         0.975quant mode
#> Precision for the Gaussian observations  0.2410992   NA
```

## Distribuciones posteriores

Para obtener la distribución marginal de los valores ajustados y predichos necesitamos incorporar  a la función `inla` el argumento `control.compute=list(return.marginals.predictor=TRUE)`.
 
Tras conseguir un ajuste con `inla`, podemos acceder a todas las distribuciones marginales posteriores y predictivas a través de:

-   `fit$marginals.fixed` da las distribuciones posteriores marginales de los efectos fijos
-   `fit$marginals.fixed$x` da la distribución del efecto fijo `x` (`x` el nombre del efecto), y también se puede seleccionar con su ordinal en el conjunto de efectos fijos
`fit$marginals.fixed[[i]]`
-     `fit.marginals.hyperpar` da las distribuciones posteriores marginales de los parámetros e hiperparámetros
-   `fit$marginals.fitted.values` da las distribuciones posteriores marginales para los valores ajustados 

Con estas distribuciones, reconocidas como `marginal`, podemos hacer cálculos y gráficos de interés a través de estas funciones que operan sobre las distribuciones y que podemos consultar con `?inla.marginal`:

-   `inla.dmarginal(x, marginal, ...)` da la densidad en x
-   `inla.pmarginal(q, marginal, ...)` da las probabilidades o función de distribución
-   `inla.qmarginal(p, marginal,...)` da los cuantiles
-   `inla.rmarginal(n, marginal)` permite obtener $n$ simulaciones
-   `inla.hpdmarginal(p, marginal,...)` da la región HPD 
-   `inla.smarginal(marginal, ...)` da un suavizado con splines de la distribución marginal
-   `inla.emarginal(fun, marginal, ...)` calcula el valor esperado de una función
-   `inla.mmarginal(marginal,...)` calcula la moda posterior
-   `inla.tmarginal(fun, marginal,...)` transforma la distribución marginal de una función del parámetro
-   `inla.zmarginal(marginal,...)` calcula descriptivos de la marginal.

Veamos algunos ejemplos sobre nuestro problema. Vamos a mostrar a continuación, en un único gráfico, las distribuciones posteriores de los efectos fijos y el parámetro de dispersión de los datos $\sigma*, con líneas verticales que marquen el valor esperado posterior (en azul) y el HPD95% en rojo.


```r
# library(gridExtra)
# library(ggplot2)

g=list() # lista en que almacenamos los gráficos con d.posteriores
# Efectos fijos
names.fixed=names(fit$marginals.fixed)
n.fixed=length(names.fixed)
names=c(expression(theta),expression(beta))
for (i in 1:n.fixed){
 g[[i]] = ggplot(as.data.frame(fit$marginals.fixed[[i]])) + 
  geom_line(aes(x = x, y = y)) +
  labs(x=names[i],y="")+
   geom_vline(xintercept=inla.hpdmarginal(0.95,fit$marginals.fixed[[i]]),
             linetype="dashed",color="red")+
  geom_vline(xintercept=inla.emarginal(function(x) x,fit$marginals.fixed[[i]]),
             linetype="dashed",color="blue")
}


# Parámetros
g[[3]]= ggplot(as.data.frame(fit$marginals.hyperpar[[1]])) + 
  geom_line(aes(x = x, y = y)) +
  labs(x=expression(tau),y="")+
  geom_vline(xintercept=inla.hpdmarginal(0.95,fit$marginals.hyperpar[[1]]),
             linetype="dashed",color="red")+
  geom_vline(xintercept=inla.emarginal(function(x) x,fit$marginals.hyperpar[[1]]),
             linetype="dashed",color="blue")

# Transformamos la posterior en tau para obtener la posterior de sigma
sigma.post=inla.tmarginal(function(tau) tau^(-1/2),
  fit$marginals.hyperpar[[1]])
# y la pintamos
g[[4]]= ggplot(as.data.frame(sigma.post)) + 
  geom_line(aes(x = x, y = y)) +
  labs(x=expression(sigma),y="")+
  geom_vline(xintercept=inla.hpdmarginal(0.95,sigma.post),
             linetype="dashed",color="red")+
  geom_vline(xintercept=inla.emarginal(function(x) x,sigma.post),
             linetype="dashed",color="blue")

library(gridExtra)
grid.arrange(g[[1]],g[[2]],g[[3]],g[[4]],ncol=2)
```

![](01-regresion_files/figure-latex/unnamed-chunk-14-1.pdf)<!-- --> 

## Simulación de la posterior

Cuando queremos inferir sobre funciones de los efectos latentes e hiperparámetros que no proporciona por defecto `INLA`, podemos recurrir a simular de las distribuciones posteriores de los efectos involucrados, y con ellas evaluar la función que nos interesa para conseguir una muestra de su distribución posterior.
Para ello es preciso que al ajustar el modelo con `inla` hayamos incluido el argumento `control.compute=list(config=TRUE)`.

Para obtener simulaciones utilizamos las funciones:

- `inla.posterior.sample(n, fit,selection)`, para simular los efectos latentes, donde $n$ es el número de simulaciones pretendido, `fit` es el ajuste obtenido con `inla` y `selection` es una lista con el nombre de las componentes a simular.  
- `inla.hyperpar.sample(n,fit,improve.marginals=TRUE)` para simular de parámetros e hiperparámetros.

Para describir las distribuciones de los nuevos parámetros que queremos evaluar con dichas simulaciones, utilizaremos:

- `inla.posterior.sample()` para funciones sobre los efectos fijos
- `inla.hyperpar.sample()` para funciones sobre los hiperparámetros


Imaginemos que queremos obtener la distribución posterior de $(\theta+\beta)$. Hemos de simular pues, de las distribuciones posteriores de $\theta$ y de $\beta$, para luego sumar las simulaciones y obtener una muestra posterior de $\theta+\beta$.


```r
fit <- inla(formula, data = davis,
  control.compute = list(config = TRUE))
sims = inla.posterior.sample(100, fit, selection = list("(Intercept)"=1,repwt = 1))
```

Esto nos devuelve una lista de la dimensión del número de simulaciones, y en cada uno de los elementos tenemos: los valores simulados de los hiperparámetros (`hyper`), de los efectos latentes (`latent`)


```r
length(sims)
#> [1] 100
names(sims[[1]])
#> [1] "hyperpar" "latent"   "logdens"
sims[[1]]$hyperpar
#> Precision for the Gaussian observations 
#>                               0.2514044
sims[[1]]$latent
#>                sample:1
#> (Intercept):1 2.9201273
#> repwt:1       0.9577521
```

y la log-densidad de la posterior en esos valores (`logdens`)


```r
 sims[[1]]$logdens
#> $hyperpar
#> [1] 0.1424985
#> 
#> $latent
#> [1] 1008.068
#> 
#> $joint
#> [1] 1008.21
```

Ahora con la función `inla.posterior.sample.eval`, dado que nuestra función depende de efectos fijos (latentes), $\theta+\beta$, evaluamos dicha suma, aproximamos los descriptivos de la posterior, y la graficamos con un histograma:


```r
theta_beta=inla.posterior.sample.eval(function(...) {(Intercept)+repwt},sims)
theta_beta=as.vector(theta_beta)
summary(theta_beta)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.9381  3.2327  3.8525  3.8418  4.4880  5.5453
ggplot(data.frame(tb=theta_beta),aes(x=tb))+
  geom_histogram(stat="density")
#> Warning: Ignoring unknown parameters: binwidth, bins, pad
```

![](01-regresion_files/figure-latex/unnamed-chunk-18-1.pdf)<!-- --> 


## Regresión lineal múltiple con INLA

Vemos a continuación cómo se trabaja la regresión múltiple con INLA, simplemente añadiendo más predictores. 

El modelo de regresión asume una distribución normal para los datos, datos los efectos latentes (todos los relacionados con el valor esperado o predictor lineal) y el resto de parámetros o hiperparámetros del modelo (en nuestro caso la varianza de los datos).

Si tenemos $n$ observaciones
$$ (y_i|\mu_i,\sigma^2) \sim N(\mu_i,\sigma^2), \ \ i=1,...n$$
donde $\mu_i$ representa la media y $\sigma^2$ la varianza de los datos.
$$E(y_i|\mu_i,\sigma^2)=\mu_i, \ Var(y_i|\mu_i,\sigma^2)=\sigma^2$$
Si tenemos varios regresores $x_1,x_2,...,x_J$, la media $\mu_i$ se relaciona linealmente con un predictor lineal $\eta_i$ que se construye a partir de una combinación lineal de los predictores:
$$\eta_i=\theta+ \sum_{j=1}^J \beta_j x_{ij}$$
En el modelo de regresión la media de los datos $\mu$ coincide con el predictor lineal $\eta$, $\mu_i=\eta_i, i=1,...,n$. 

Nos queda a continuación especificar las distribuciones a priori sobre el vector de efectos latentes, en nuestro caso efectos fijos, $(\theta,\beta_1,...,\beta_J)$, y sobre el parámetro de dispersión o hiperparámetro $\sigma^2$. 

Cuando no tenemos información previa especificamos distribuciones difusas (vagas) sobre los parámetros. En INLA por defecto tendremos:

\begin{eqnarray*}
\theta &\sim& N(0,\sigma_{\theta}^2)  \\
\beta_j &\sim& N(0,\sigma_{\beta}^2) \\
log(\tau=1/\sigma^2) &\sim& Log-Ga(1,0.00005)
\end{eqnarray*}
con $\sigma_{theta}^2=\infty$ y $\sigma_{\beta}=1000$.

Ejemplificamos el análisis de regresión lineal múltiple sobre la base de datos `usair`, en la librería `brinla`. Esta BD contiene datos recopilados para investigar los factores determinantes de la polución, utilizando el nivel de SO2 como variable dependiente y las restantes como variables explicativas potenciales. Las relaciones entre las variables que incluye se muestra en la Figura \@ref(fig:regmul01) a continuación.


```r
data(usair, package = "brinla")
library(GGally)
pairs.chart <- ggpairs(usair, lower = list(continuous = "cor"), 
                       upper = list(continuous = "points", combo = "dot"))
ggplot2::theme(axis.text = element_text(size = 6)) 
#> List of 1
#>  $ axis.text:List of 11
#>   ..$ family       : NULL
#>   ..$ face         : NULL
#>   ..$ colour       : NULL
#>   ..$ size         : num 6
#>   ..$ hjust        : NULL
#>   ..$ vjust        : NULL
#>   ..$ angle        : NULL
#>   ..$ lineheight   : NULL
#>   ..$ margin       : NULL
#>   ..$ debug        : NULL
#>   ..$ inherit.blank: logi FALSE
#>   ..- attr(*, "class")= chr [1:2] "element_text" "element"
#>  - attr(*, "class")= chr [1:2] "theme" "gg"
#>  - attr(*, "complete")= logi FALSE
#>  - attr(*, "validate")= logi TRUE
pairs.chart
```

![(\#fig:regmul01)Relaciones entre variables en la base de datos usair(brinla).](01-regresion_files/figure-latex/regmul01-1.pdf) 

Apreciamos ya en el gráfico una correlación positiva muy alta entre las variables `pop` y `manuf`, y relevante para `negtemp` y `days` y también para `precip` y `days`, lo que posteriormente condicionará la selección de variables. 

Cuando trabajamos con más de una variable predictora en regresión lineal (realmente en cualquier modelo) surge un problema adicional, que es la **selección de variables**, o selección del mejor modelo de predicción. Esto se resuelve en INLA utilizando diversos criterios de selección entre los que destacamos:
  
- la verosimilitud marginal (valor de la log-verosimilitud): `mlik`; al cambiarle el signo tendremos $-loglikelihood$
- el criterio de información de la deviance (DIC): `dic`
- el criterio de información bayesiana ampliado (WAIC): `waic`.

El procedimiento a utilizar para la selección del modelo (de variables) es el siguiente:

1. ajustar todos los modelos resultantes de todas las combinaciones posibles de predictoras, 
1. calcular los índices de selección para cada uno de ellos 
1. proceder con la selección en base a dichos valores. 

Siempre se prefieren modelos con los valores más bajos para estos criterios y se descartan los que proporcionan valores más altos. Cuando no es el mismo modelo el que proporciona el valor mínimo en estos criterios, habremos de optar por alguno de ellos.

Por defecto, al ajustar un modelo con `inla`, nos devuelve la log-verosimilitud (accesible con `fit$mlik` si el ajuste se guardó en el objeto `fit`). Para obtener las otras medidas de selección, hemos de incluir como argumento de la función `inla`, la opción `control.compute = list(dic = TRUE, waic = TRUE))`. Los valores por defecto de `control.compute` los podemos consultar ejecutando el comando `inla.set.control.compute.default()`. 

Vayamos pues con el ajuste del modelo de regresión con todas las variables. Recordemos que si no especificamos el argumento `family`, interpreta por defecto la opción `gaussian`, esto es, normalidad para los datos. Así ajustamos el modelo y obtenemos las siguientes inferencias sobre los efectos fijos.


```r
formula <-  SO2 ~ negtemp + manuf + wind + precip + days
fit= inla(formula, data = usair, control.compute = list(dic = TRUE, waic = TRUE))
#summary(fit)
round(fit$summary.fixed,3)
#>                mean     sd 0.025quant 0.5quant 0.975quant
#> (Intercept) 135.490 49.850     37.178  135.497    233.760
#> negtemp       1.769  0.634      0.518    1.769      3.019
#> manuf         0.026  0.005      0.017    0.026      0.035
#> wind         -3.723  1.934     -7.536   -3.723      0.092
#> precip        0.625  0.387     -0.139    0.625      1.388
#> days         -0.057  0.174     -0.400   -0.057      0.287
#>             mode kld
#> (Intercept)   NA   0
#> negtemp       NA   0
#> manuf         NA   0
#> wind          NA   0
#> precip        NA   0
#> days          NA   0
```

La inferencia posterior sobre la desviación típica de los datos, $\sigma$, la resolvemos con sus descriptivos.


```r
sigma.post= inla.tmarginal(function(tau) tau^(-1/2),
                           fit$marginals.hyperpar[[1]])
inla.zmarginal(sigma.post)
#> Mean            15.6706 
#> Stdev           1.85534 
#> Quantile  0.025 12.5293 
#> Quantile  0.25  14.347 
#> Quantile  0.5   15.4919 
#> Quantile  0.75  16.7975 
#> Quantile  0.975 19.817
```

Para seleccionar las variables relevantes seguimos el procedimiento descrito anteriormente. Añadimos también el ajuste del modelo de regresión frecuentista, para el que calculamos como criterio de bondad de ajuste el AIC.


```r
# variables en la bd
vars=names(usair)[-1] # variables predictoras (excluye v.dpte)
nvars=length(vars) # nº variables predictoras
# Selección de variables. truco para concatenar todos los modelos posibles en una fórmula
listcombo <- unlist(sapply(1:nvars,function(x) combn(nvars, x, simplify=FALSE)),recursive=FALSE)
predterms <- lapply(listcombo, function(x) paste(vars[x],collapse="+"))
nmodels <- length(listcombo)
coefm <- matrix(NA,length(listcombo),4,dimnames=list(predterms,c("AIC","DIC","WAIC","MLIK")))
for(i in 1:nmodels){
  formula <- as.formula(paste("SO2 ~ ",predterms[[i]]))
  # modelo frecuentista
  lmi <- lm(formula, data=usair)
  # modelo bayesiano
  result <- inla(formula, family="gaussian", data=usair, control.compute=list(dic=TRUE, waic=TRUE))
  coefm[i,1] <- AIC(lmi)
  coefm[i,2] <- result$dic$dic
  coefm[i,3] <- result$waic$waic
  coefm[i,4] <- -result$mlik[1]
}
```

Ya solo resta comparar los resultados, respecto de cada uno de los criterios, para seleccionar con qué variables nos quedamos, y por lo tanto con qué modelo de predicción. Basta con encontrar el modelo que proporciona el mínimo valor en cada uno de los criterios.


```r
gana.aic = predterms[which.min(coefm[,1])]
gana.dic = predterms[which.min(coefm[,2])]
gana.waic = predterms[which.min(coefm[,3])]
gana.mlik = predterms[which.min(coefm[,4])]
gana.aic;gana.dic
#> [[1]]
#> [1] "negtemp+manuf+pop+wind+precip"
#> [[1]]
#> [1] "negtemp+manuf+pop+wind+precip"
gana.waic;gana.mlik
#> [[1]]
#> [1] "negtemp+manuf+pop+wind+precip"
#> [[1]]
#> [1] "manuf"
```

Concluimos pues que, tanto el criterio AIC en el modelo de regresión frecuentista, como los criterios DIC y WAIC en el modelo de regresión bayesiano, proporcionan el mejor ajuste. Este mejor modelo incluye como predictores las variables 'negtemp+manuf+pop+wind+precip'. Reajustamos el modelo con estas variables para derivar las inferencias y predicciones.

Las inferencias sobre los efectos fijos y la precisión de los datos se muestran a continuación.


```r
formula=SO2 ~ negtemp+manuf+pop+wind+precip
fit=inla(formula, family="gaussian", data=usair, control.compute=list(dic=TRUE, waic=TRUE))
fit$summary.fixed
#>                     mean          sd   0.025quant
#> (Intercept) 100.01549484 30.15083880 40.555055784
#> negtemp       1.12026612  0.41444673  0.302988558
#> manuf         0.06489651  0.01549632  0.034341608
#> pop          -0.03934787  0.01488967 -0.068708294
#> wind         -3.07258456  1.75727508 -6.536740141
#> precip        0.41922295  0.21555728 -0.005830357
#>                 0.5quant    0.975quant mode          kld
#> (Intercept) 100.01906558 159.455494342   NA 1.152499e-08
#> negtemp       1.12029301   1.937389819   NA 1.155082e-08
#> manuf         0.06489625   0.095452868   NA 1.156335e-08
#> pop          -0.03934751  -0.009989461   NA 1.156286e-08
#> wind         -3.07283892   0.393027080   NA 1.150245e-08
#> precip        0.41922935   0.844239641   NA 1.156060e-08
fit$summary.hyperpar
#>                                                mean
#> Precision for the Gaussian observations 0.005073499
#>                                                  sd
#> Precision for the Gaussian observations 0.001196977
#>                                          0.025quant
#> Precision for the Gaussian observations 0.003032416
#>                                            0.5quant
#> Precision for the Gaussian observations 0.004968958
#>                                          0.975quant mode
#> Precision for the Gaussian observations 0.007512641   NA
```

Y en la Figura \@ref(fig:regmul02) se muestran las distribuciones posteriores de todos los efectos latentes (efectos fijos) en el modelo: interceptación y coeficientes de los regresores.


```r
nfixed=length(fit$names.fixed)
g=list()
for(i in 1:nfixed){
  g[[i]]=ggplot(as.data.frame(fit$marginals.fixed[[i]])) + 
    geom_line(aes(x = x, y = y)) +
    labs(x=fit$names.fixed[i],y="")+
    geom_vline(xintercept=inla.hpdmarginal(0.95,fit$marginals.fixed[[i]]),
               linetype="dashed",color="red")+
    geom_vline(xintercept=inla.emarginal(function(x) x,fit$marginals.fixed[[i]]),
               linetype="dashed",color="blue")
}
grid.arrange(g[[1]],g[[2]],g[[3]],g[[4]],g[[5]],g[[6]],ncol=2)
```

![(\#fig:regmul02)Distribuciones posteriores de los efectos latentes.](01-regresion_files/figure-latex/regmul02-1.pdf) 

Por último, si queremos predecir nuevas observaciones con la **distribución predictiva a posteriori**, y ciertos valores de las variables explicativas, utilizamos el argumento `control.predictor`en la función `inla`, especificando en una lista los valores a predecir. Veamos cómo hacerlo con `inla`, ajustando un modelo sobre un `data.frame` combinado, en el que añadimos tantas filas como predicciones queremos conseguir, con los valores deseados para los predictores (y valores faltantes en la respuesta), y especificamos los valores a predecir en un vector indicador, a través del argumento `control.predictor(list(link=vector.indicador))`. Las predicciones se obtienen después, con el resumen de los datos ajustados `fitted.values`. Recordemos que para poder mostrar las distribuciones predictivas, hemos de añadir en `inla` el argumento `control.compute=list(return.marginals.predictor=TRUE)`. 


```r
# Predicción con INLA
## valores de los predictores en los que predecir: 3 escenarios
formula=SO2 ~ negtemp+manuf+pop+wind+precip
new.data <- data.frame(negtemp = c(-50, -60, -40), 
                       manuf = c(150, 100, 400), 
                       pop = c(200, 100, 300), 
                       wind = c(6, 7, 8), 
                       precip = c(10, 30, 20),
                       days=c(NA, NA,NA))

## añadimos los tres escenarios de predicción a la bd original, 
## dejando como faltantes los valores a predecir de la v.dpte
usair.combinado <- rbind(usair, data.frame(SO2=c(NA,NA,NA),new.data))
## creamos un vector con NA's para observaciones y 1's para predicciones
usair.indicador <- c(rep(NA, nrow(usair)), rep(1, nrow(new.data)))
## reajustamos el modelo añadiendo la opción de predicción de datos
fit.pred <- inla(formula, data = usair.combinado, 
                 control.compute=list(return.marginals.predictor=TRUE),
                 control.predictor = list(link = usair.indicador))
## y describimos los valores ajustados para los tres escenarios añadidos
fit.pred$summary.fitted.values[(nrow(usair)+1):nrow(usair.combinado),]
#>                         mean       sd 0.025quant 0.5quant
#> fitted.Predictor.42 31.62397 8.209857   15.43620 31.62473
#> fitted.Predictor.43 26.42298 5.477699   15.62271 26.42334
#> fitted.Predictor.44 53.16307 7.096415   39.17095 53.16362
#>                     0.975quant mode
#> fitted.Predictor.42   47.80712   NA
#> fitted.Predictor.43   37.22103   NA
#> fitted.Predictor.44   67.15178   NA
```

Así graficamos en la Figura \@ref(fig:regmul03) la distribución predictiva del nivel de SO2 para una combinación dada de valores de las variables predictivas, negtem=-50, manuf=150, pop=200, wind=6 y precip=10.


```r
pred=fit.pred$marginals.fitted.values[[42]]
ggplot(as.data.frame(pred)) + 
  geom_line(aes(x = x, y = y)) +
  labs(x=expression(eta),y="")+
  geom_vline(xintercept=inla.hpdmarginal(0.95,pred),
             linetype="dashed",color="red")+
  geom_vline(xintercept=inla.emarginal(function(x) x,pred),
             linetype="dashed",color="blue")
```

![(\#fig:regmul03)Distribución predictiva a posteriori de SO2 para una configuración dada de los predictores.](01-regresion_files/figure-latex/regmul03-1.pdf) 

## Conclusiones

En este tema hemos trabajado con el ajuste con INLA de modelos lineales de regresión, esto es, con efectos fijos. En posteriores temas trabajaremos modelos más sofisticados en los que incluiremos los efectos aleatorios, generalizaremos el modelo lineal y entenderemos el planteamiento de modelos a través de modelos jerárquicos bayesianos.
