# Modelo de ANOVA {#anova}


## Introducción


El modelo de ANOVA se plantea para comparar poblaciones normales, especialmente cuando son más de dos las poblaciones a comparar. Las poblaciones a comparar se identifican a través de una variable clasificadora (de tipo categórico) que actúa como predictora para estimar respuestas medias supuestamente distintas a comparar. 

## El modelo de ANOVA

Consideremos una variable respuesta $Y$ que se distribuye normal, y que viene afectada por una variable de clasificación $A$ con $a$ niveles de respuesta distintos (uno por cada una de las poblaciones a comparar). Supongamos que tenemos $n_i$ observaciones de la respuesta para cada uno de los niveles de respuesta de la variable clasificadora, $i=1,...,a$. El modelo de ANOVA se plantea asumiendo que en cada nivel o subpoblación, esperamos una valor distinto para la respuesta,
$$(y_{ij}|\mu_i,\sigma^2) \sim N(\mu_i,\sigma^2)$$
de modo que 
$$E(y_{ij}|\mu_i,\sigma^2)=\mu_i; \ Var(y_{ij}|\mu_i,\sigma^2)=\sigma^2, \ \ i=1,...,a; \ j=1,...,n_i$$
La formulación habitual de este modelo se suele dar en términos de un efecto global y común a todas las observaciones, $\theta$, y un efecto diferencial respecto del primer nivel del factor de clasificación $A$, $\alpha_i$, con los que se construye la media (identificada generalmente por $\mu$) o predictor lineal (identificada generalmente por $\eta$) y que en el modelo lineal coinciden:
$$\mu_{ij}=\eta_{ij}=\theta + \alpha_i$$
donde $\alpha_i=E(y_{ij}|\mu,\sigma^2)-E(y_{1j}|\mu,\sigma^2), i\geq 1$, esto es, $\alpha=1=0$.

En la modelización bayesiana es preciso añadir distribuciones a priori para cada uno de los parámetros del modelo: los efectos fijos $\theta,\alpha_i$, y la varianza  $\sigma^2$ de los datos. Ante ausencia de información, se asumirán distribuciones difusas:

\begin{eqnarray*}
(Y_{ij}|\mu_i,\sigma^2) & \sim & N(\mu_i,\sigma^2) \\
&& \mu_i = \theta + \alpha_i, i=1,...,a \\
\theta & \sim & N(0,\sigma_{\theta}^2), \ \sigma_{\theta}^2=0 \\
\alpha_i & \sim & N(0,1000), i\geq 1\\
\tau=1/\sigma^2 & \sim & Ga(1,0.00005)
\end{eqnarray*}

## Anova de una vía

Vamos a ilustrar el análisis ANOVA en INLA a través de la base de datos `coagulation`, en la librería [`faraway`](https://cran.r-project.org/web/packages/faraway/faraway.pdf), referidos a un estudio de tiempos de coagulación de la sangre en 24 animales a los que aleatoriamente se les asignó una de entre tres dietas distintas (variable clasificadora `diet`). Posteriormente, para estudiar el efecto de dichas dietas en la coagulación, se tomaron muestras de los tiempos de coagulación (en la variable `coag`, que es la respuesta).


```r
data(coagulation, package="faraway")
str(coagulation)
#> 'data.frame':	24 obs. of  2 variables:
#>  $ coag: num  62 60 63 59 63 67 71 64 65 66 ...
#>  $ diet: Factor w/ 4 levels "A","B","C","D": 1 1 1 1 2 2 2 2 2 2 ...
ggplot(coagulation,aes(x=diet,y=coag))+
  geom_boxplot()
```

![](02-anova_files/figure-latex/unnamed-chunk-2-1.pdf)<!-- --> 

Estamos planteando un modelo de Anova como el propuesto en la sección anterior, donde $\alpha_i$ identifica el efecto diferencial sobre la respuesta con la dieta A, para el resto de las dietas B y C. Los parámetros del modelo son, como en regresión, los efectos fijos $(\theta,\alpha_i)$ y la varianza $\sigma^2$, para los que asumimos las priors difusas que por defecto propone INLA. Ajustamos el modelo y obtenemos las inferencias a posteriori


```r
formula=coag ~ diet
fit=inla(formula,family="gaussian",data=coagulation,
         control.compute=list(return.marginals.predictor=TRUE))
fijos=round(fit$summary.fixed,3);fijos
#>               mean    sd 0.025quant 0.5quant 0.975quant
#> (Intercept) 61.016 1.172     58.700   61.016     63.337
#> dietB        4.979 1.513      1.983    4.980      7.970
#> dietC        6.977 1.513      3.981    6.978      9.968
#> dietD       -0.016 1.435     -2.859   -0.016      2.822
#>             mode kld
#> (Intercept)   NA   0
#> dietB         NA   0
#> dietC         NA   0
#> dietD         NA   0
tau=round(fit$summary.hyperpar,3);tau
#>                                          mean    sd
#> Precision for the Gaussian observations 0.197 0.059
#>                                         0.025quant 0.5quant
#> Precision for the Gaussian observations      0.099    0.191
#>                                         0.975quant mode
#> Precision for the Gaussian observations      0.328   NA
medias=round(fit$summary.linear.predictor,4)
```
Atendiendo a los descriptivos de la distribución posterior para los efectos fijos, concluimos:

- El tiempo esperado de coagulación para los animales que han seguido la dieta A es de 61.016(58.7,63.337).
- Los animales que han seguido la dieta B tienen un tiempo de coagulación esperado  superior en 4.979 unidades a los de la dieta A, y dicha diferencia es *significativamente distinta de cero* en el contexto bayesiano, dado que su HPD no incluye al cero, (1.983,7.97). 
- Una conclusión similar se deriva para la dieta C, que da un tiempo de coagulación esperado  superior en 6.977 unidades a los de la dieta A, y un HPD (3.981,9.968).
- Las diferencias en los tiempos de coagulación de seguir una dieta D frente a la dieta A no son relevantes. De hecho, la diferencia entre ellos es de -0.016 y el intervalo HPD contiene al cero, (-2.859,2.822).


Pintamos a continuación en la Figura \@ref(fig:anova01) la distribución posterior de las medias $\mu_i$ (o predictores lineales) para cada una de las dietas.


```r
dietas=levels(coagulation$diet)
pred=NULL
for(i in 1:length(dietas)){
index=which(coagulation$diet==dietas[i])[1]
# distrib. posterior
post=fit$marginals.fitted.values[[index]]
# media
e=fit$summary.fitted.values[index,1]
hpd.low=fit$summary.fitted.values[index,3]
hpd.up=fit$summary.fitted.values[index,5]
pred=rbind(pred,data.frame(dieta=dietas[i],
                           post,e=e,
                           hpd.low=hpd.low,hpd.up=hpd.up))
}

ggplot(pred, aes(x = x, y =y)) + 
  geom_line(aes(color=dieta))+
  labs(x=expression(paste("Tiempo medio de coagulación:",mu)),
       y="D.Posterior")
```

![(\#fig:anova01)Distribución posterior del tiempo medio de coagulación para las 4 dietas.](02-anova_files/figure-latex/anova01-1.pdf) 


```r
ggplot(pred, aes(x = x, y =y)) + 
  geom_line(aes(color=dieta))+
  geom_vline(aes(xintercept=e,color=dieta),linetype="dashed")+
  geom_vline(aes(xintercept=hpd.low,color=dieta),linetype="dotted")+
  geom_vline(aes(xintercept=hpd.up,color=dieta),linetype="dotted")+
  facet_wrap(vars(dieta))+
  labs(x=expression(paste("Tiempo medio de coagulación:",mu)),
       y="D.Posterior",title="D.Posterior, medias y HPD95%")
```

![(\#fig:anova02)Distribuciones posteriores, medias y HPD](02-anova_files/figure-latex/anova02-1.pdf) 

Como ya hacíamos en regresión, podemos inferir sobre la desviación típica de los datos, $\sigma$, transformando la distribución para la precisión $\tau$.


```r
sigma.post=inla.tmarginal(function(tau) tau^(-1/2),
  fit$marginals.hyperpar[[1]])
# y la pintamos
ggplot(as.data.frame(sigma.post)) + 
  geom_line(aes(x = x, y = y)) +
  labs(x=expression(sigma),y="D.Posterior")+
  geom_vline(xintercept=inla.hpdmarginal(0.95,sigma.post),
             linetype="dotted",color="blue")+
  geom_vline(xintercept=inla.emarginal(function(x) x,sigma.post),
             linetype="dashed",color="blue")
```

![](02-anova_files/figure-latex/unnamed-chunk-4-1.pdf)<!-- --> 

```r

# Valor esperado
sigma.e=round(inla.emarginal(function(tau) tau^(-1/2),
  fit$marginals.hyperpar[[1]]),4)
# HPD95%
sigma.hpd=round(inla.hpdmarginal(0.95,sigma.post),3)
paste("E(sigma.post)=",sigma.e,"HPD95%=(",sigma.hpd[1],",",sigma.hpd[2],")")
#> [1] "E(sigma.post)= 2.3347 HPD95%=( 1.681 , 3.069 )"
```

## Anova de varias vías

Generalmente, y en especial cuando trabajamos con experimentación, son varios los factores que controlamos para investigar el efecto que producen en una respuesta continua. Hablamos de modelos de Anova de varias vías.

Utilizamos la base de datos `butterfat` en la librería [`faraway`](https://cran.r-project.org/web/packages/faraway/faraway.pdf) para ilustrar el ajuste con INLA de un modelo de Anova de varias vías. Esta base de datos contiene 100 registros del contenido en grasa láctea, `Butterfat`, para muestras aleatorias de 20 vacas (10 de ellas de 2 años y 10 maduras, con más de 4 años -en la variable `Age`) de cada una de cinco razas (en la variable `Breed`).

El objetivo es investigar las diferencias en materia grasa entre razas y edad, con el fin último de identificar cuáles producen más materia grasa y cuáles menos. Veamos los datos en la Figura \@ref(fig:butterfat1).


```r
data(butterfat,package="faraway")
str(butterfat)
#> 'data.frame':	100 obs. of  3 variables:
#>  $ Butterfat: num  3.74 4.01 3.77 3.78 4.1 4.06 4.27 3.94 4.11 4.25 ...
#>  $ Breed    : Factor w/ 5 levels "Ayrshire","Canadian",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ Age      : Factor w/ 2 levels "2year","Mature": 2 1 2 1 2 1 2 1 2 1 ...
ggplot(butterfat,aes(x=Breed,y=Butterfat))+
  geom_boxplot(aes(color=Age))+
  theme(axis.text.x = element_text(angle = 45))
```

![(\#fig:butterfat1)Base de datos butterfat, en la librería Faraway](02-anova_files/figure-latex/butterfat1-1.pdf) 

A la vista del gráfico, apreciamos que por lo general, en la mayoría de las razas, las vacas más jóvenes tienen menor contenido en materia grasa que las más viejas. Sin embargo, tal afirmación no parece tan clara en las razas *Guernsey* y *Holstein-Fresian*, de modo que para modelizar nuestros datos vamos a considerar a priori, la posibilidad de interacciones entre los factores de clasificación `Breed` y `Age`.


Cuando nos enfrentamos a varios factores de clasificación, cabe la posibilidad de que interaccionen entre ellos, esto es, que en algunos niveles de un factor actúen de forma diferente a los otros cuando se combinan con los niveles de algún otro factor. El orden de una interacción viene dado por el número de factores de clasificación que involucra, de modo que hablamos de interacciones de orden 2 si consideramos la interacción entre dos factores, de orden 3 si consideramos la interacción entre tres factores, etc. Generalmente trabajamos con interacciones de orden bajo, dada la complejidad de las conclusiones en interacciones de orden alto. Por otro lado, siempre es importante tener en cuenta de cuántos datos disponemos para conocer a priori la posibilidad de estimar con fiabilidad los distintos efectos de interacción (una interacción de dos factores con $n_1$ y $n_2$ niveles de clasificación respectivamente, revierte en la estimación de $(n_1-1)\times(n_2-1)$ efectos de interacción).

Así, en nuestro problema si estamos planteando la posibilidad de que haya interacciones entre los dos factores de clasificación, estamos asumiendo un modelo de tipo siguiente, asumiendo normalidad en la respuesta:

$$(y_{ijk}|\eta_{ij},\sigma^2) \sim N(\eta_{ij},\sigma^2)$$
con 
$$\eta_{ij}=\theta+ \alpha_i + \beta_j + \alpha\beta_{ij}$$
donde $\alpha_i$ es el efecto diferencial (respecto del primer nivel) que aporta el nivel $i$ de la variable `Breed`, $\beta_j$ el efecto asociado a la variable `Age`, y $\alpha\beta$ la correspondiente interacción entre ellas. En *R* una interacción de orden 2 entre dos variables $f_1$ y $f_2$ se especifica con $f_1:f_2$; los efectos principales y la interacción también se pueden especificar con $f_1+f_2+f_1:f_2=f_1*f_2=(f_1+f_2)$^2.

Veamos cómo ajustar con INLA este modelo.


```r
formula=Butterfat ~ Breed * Age
fit=inla(formula,data=butterfat,
         control.compute=list(dic = TRUE, waic = TRUE))
fit$summary.fixed
#>                                        mean        sd
#> (Intercept)                      3.96605050 0.1314176
#> BreedCanadian                    0.52193553 0.1858549
#> BreedGuernsey                    0.93293190 0.1858549
#> BreedHolstein-Fresian           -0.30304829 0.1858549
#> BreedJersey                      1.16693161 0.1858549
#> AgeMature                        0.18793906 0.1858476
#> BreedCanadian:AgeMature         -0.28692013 0.2628334
#> BreedGuernsey:AgeMature         -0.08591997 0.2628334
#> BreedHolstein-Fresian:AgeMature -0.17493825 0.2628334
#> BreedJersey:AgeMature            0.13107657 0.2628334
#>                                 0.025quant    0.5quant
#> (Intercept)                      3.7077565  3.96604997
#> BreedCanadian                    0.1566437  0.52193621
#> BreedGuernsey                    0.5676400  0.93293262
#> BreedHolstein-Fresian           -0.6683396 -0.30304778
#> BreedJersey                      0.8016397  1.16693233
#> AgeMature                       -0.1773380  0.18793970
#> BreedCanadian:AgeMature         -0.8035053 -0.28692097
#> BreedGuernsey:AgeMature         -0.6025052 -0.08592082
#> BreedHolstein-Fresian:AgeMature -0.6915240 -0.17493890
#> BreedJersey:AgeMature           -0.3855087  0.13107576
#>                                 0.975quant mode
#> (Intercept)                     4.22434753   NA
#> BreedCanadian                   0.88722350   NA
#> BreedGuernsey                   1.29821977   NA
#> BreedHolstein-Fresian           0.06224017   NA
#> BreedJersey                     1.53221946   NA
#> AgeMature                       0.55321251   NA
#> BreedCanadian:AgeMature         0.22966986   NA
#> BreedGuernsey:AgeMature         0.43067002   NA
#> BreedHolstein-Fresian:AgeMature 0.34165120   NA
#> BreedJersey:AgeMature           0.64766646   NA
#>                                          kld
#> (Intercept)                     2.741551e-09
#> BreedCanadian                   2.741827e-09
#> BreedGuernsey                   2.741827e-09
#> BreedHolstein-Fresian           2.741827e-09
#> BreedJersey                     2.741832e-09
#> AgeMature                       2.740990e-09
#> BreedCanadian:AgeMature         2.741409e-09
#> BreedGuernsey:AgeMature         2.741408e-09
#> BreedHolstein-Fresian:AgeMature 2.741408e-09
#> BreedJersey:AgeMature           2.741408e-09
fit$dic$dic
#> [1] 120.4909
fit$waic$waic
#> [1] 121.7376
```

Observamos en la inferencia posterior para los efectos fijos, que todas las HPD asociadas a los efectos de interacción contienen al cero, lo que descarta la relevancia de la interacción a la hora de predecir la respuesta. Reajustamos pues el modelo eliminando la interacción, y comprobamos que efectivamente al eliminarla conseguimos reducir los valores del DIC y WAIC que usamos habitualmente para la selección de variables.


```r
formula=Butterfat ~ Breed + Age
fit=inla(formula,data=butterfat,
         control.predictor=list(compute=TRUE),
         control.compute=list(return.marginals.predictor=TRUE,
                              dic = TRUE, waic = TRUE))
fit$summary.fixed
#>                             mean         sd  0.025quant
#> (Intercept)            4.0077184 0.10124856  3.80873708
#> BreedCanadian          0.3784787 0.13071131  0.12159372
#> BreedGuernsey          0.8899744 0.13071131  0.63308923
#> BreedHolstein-Fresian -0.3905147 0.13071131 -0.64739952
#> BreedJersey            1.2324714 0.13071131  0.97558622
#> AgeMature              0.1045993 0.08267006 -0.05787063
#>                         0.5quant 0.975quant mode
#> (Intercept)            4.0077182  4.2067007   NA
#> BreedCanadian          0.3784790  0.6353625   NA
#> BreedGuernsey          0.8899746  1.1468580   NA
#> BreedHolstein-Fresian -0.3905145 -0.1336307   NA
#> BreedJersey            1.2324717  1.4893550   NA
#> AgeMature              0.1045993  0.2670692   NA
#>                                kld
#> (Intercept)           2.524753e-09
#> BreedCanadian         2.524905e-09
#> BreedGuernsey         2.524905e-09
#> BreedHolstein-Fresian 2.524905e-09
#> BreedJersey           2.524907e-09
#> AgeMature             2.525137e-09
fit$waic$waic
#> [1] 116.3439
fit$dic$dic
#> [1] 115.3808
```

Observamos ya a partir del modelo ajustado, que el efecto de la edad no es relevante (su HPD incluye al cero), pero sin embargo sí que hay diferencias debido a las razas.

Reajustamos de nuevo el modelo, excluyendo la variable `Age`, y verificamos la reducción (ligera) del DIC/WAIC, lo cual justifica usar este modelo para la predicción.


```r
formula=Butterfat ~ Breed 
fit=inla(formula,data=butterfat,
         control.predictor=list(compute=TRUE),
         control.compute=list(return.marginals.predictor=TRUE,
                              dic = TRUE, waic = TRUE))
fit$summary.fixed
#>                             mean         sd 0.025quant
#> (Intercept)            4.0600181 0.09271799  3.8778057
#> BreedCanadian          0.3784786 0.13112332  0.1207894
#> BreedGuernsey          0.8899742 0.13112332  0.6322849
#> BreedHolstein-Fresian -0.3905148 0.13112332 -0.6482038
#> BreedJersey            1.2324713 0.13112332  0.9747819
#>                         0.5quant 0.975quant mode
#> (Intercept)            4.0600180  4.2422316   NA
#> BreedCanadian          0.3784788  0.6361666   NA
#> BreedGuernsey          0.8899745  1.1476620   NA
#> BreedHolstein-Fresian -0.3905146 -0.1328266   NA
#> BreedJersey            1.2324715  1.4901590   NA
#>                                kld
#> (Intercept)           2.473535e-09
#> BreedCanadian         2.473520e-09
#> BreedGuernsey         2.473520e-09
#> BreedHolstein-Fresian 2.473520e-09
#> BreedJersey           2.473527e-09
fit$waic$waic
#> [1] 115.9279
fit$dic$dic
#> [1] 115.0074
```

Procederíamos igual que en el modelo de Anova de una vía para la representación de las distribuciones posteriores sobre las medias o predictores lineales en cada una de las razas. Igualmente representaremos la distribución posterior del parámetro de dispersión de los datos $\sigma$.




## Análisis de ANCOVA

En ocasiones tenemos una variable respuesta de tipo numérico, y como posibles predictores variables de tipo numérico y también variables clasificadoras o factores. Surge entonces la posibilidad de que los predictores numéricos afecten a la respuesta de modo distinto en diferentes niveles de clasificación de los factores; hablamos entonces de interacción entre covariables y factores. Veamos un ejemplo para comprender cómo funcionan estos modelos y cómo se ajustan con INLA.

Consideramos los datos de Galton sobre la regresión de las alturas de los hijos sobre la de los padres (Fte: [Galton's Height Data](http://www.randomservices.org/random/)). Tenemos la estatura del padre, de la madre y del hijo/a, identificado/a por su sexo.
Vamos a formular un modelo de regresión de la estatura de los hijos en función de la de sus padres y su género.


```r
my.dir="~/Dropbox/ESTADISTICA/BAYESIAN/VARIOS/"
datos<-read.csv(file=paste0(my.dir,"Galton.txt"),header=TRUE,dec=".", sep="")
str(datos)
#> 'data.frame':	898 obs. of  6 variables:
#>  $ Family: chr  "1" "1" "1" "1" ...
#>  $ Father: num  78.5 78.5 78.5 78.5 75.5 75.5 75.5 75.5 75 75 ...
#>  $ Mother: num  67 67 67 67 66.5 66.5 66.5 66.5 64 64 ...
#>  $ Gender: chr  "M" "F" "F" "F" ...
#>  $ Height: num  73.2 69.2 69 69 73.5 72.5 65.5 65.5 71 68 ...
#>  $ Kids  : int  4 4 4 4 4 4 4 4 2 2 ...
```

Asumimos pues como respuesta la variable `y=Height`, como regresores las variables $x_1=$`Father` y $x_2=$`Mother` con las estaturas del padre y la madre respectivamente, y con factor de clasificación la variable $G=$`Gender`, con niveles M/F. En principio cabrían posibles interacciones entre los regresores y los factores de clasificación, por lo que planteamos el modelo:

$$(y_{ij}|\eta_{ij},\sigma^2) \sim N(\eta_{ij},\sigma^2)$$
con el predictor lineal 
$$\eta_{ij}=\mu_{ij}=\beta_0+(\beta_1 + \alpha_{1M}) x_{1j} + (\beta_2+ \alpha_{2M}) x_{2j} + \alpha_M;\ \  j =1,...,n_i; i=M,F$$
donde $\alpha_M$ es el efecto diferencial global de los hombres frente a las mujeres al predecir la estatura, y $\alpha_{1M},\alpha_{2M}$ los efectos diferenciales que afectan a los regresores.

Asumimos una distribución vaga sobre $\beta_0, \beta_1$ y $\tau=1/\sigma^2$ y ajustamos el modelo Gausiano en INLA:

```r
formula = Height ~ 1+(Father+Mother)*Gender
fit = inla(formula,family = "gaussian",data=datos)
round(fit$summary.fixed,3)
#>                  mean    sd 0.025quant 0.5quant 0.975quant
#> (Intercept)    16.652 3.887      9.028   16.652     24.276
#> Father          0.400 0.039      0.324    0.400      0.477
#> Mother          0.307 0.045      0.218    0.307      0.396
#> GenderM         2.707 5.427     -7.939    2.707     13.353
#> Father:GenderM  0.012 0.058     -0.103    0.012      0.126
#> Mother:GenderM  0.027 0.062     -0.096    0.027      0.149
#>                mode kld
#> (Intercept)      NA   0
#> Father           NA   0
#> Mother           NA   0
#> GenderM          NA   0
#> Father:GenderM   NA   0
#> Mother:GenderM   NA   0
```

Observamos que ninguna de las interacciones tienen un efecto a considerar (su HPD posterior incluye al cero), de modo que las descartamos y reajustamos el modelo sin ellas. 

$$\eta_{ij}=\mu_{ij}=\beta_0+ \beta_1  x_{1j} + \beta_2 x_{2j} + \alpha_M;\ \  j =1,...,n_i; i=M,F.$$


```r
formula = Height ~ Father+Mother+Gender
fit = inla(formula,family = "gaussian",data=datos,
        control.predictor=list(compute=TRUE),
         control.compute=list(return.marginals.predictor=TRUE,
                              dic = TRUE, waic = TRUE))
round(fit$summary.fixed,3)
#>               mean    sd 0.025quant 0.5quant 0.975quant
#> (Intercept) 15.345 2.747      9.957   15.345     20.733
#> Father       0.406 0.029      0.349    0.406      0.463
#> Mother       0.321 0.031      0.260    0.321      0.383
#> GenderM      5.226 0.144      4.943    5.226      5.508
#>             mode kld
#> (Intercept)   NA   0
#> Father        NA   0
#> Mother        NA   0
#> GenderM       NA   0
```

Ahora todos los efectos fijos son relevantes para predecir la estatura de los hijos. Utilizamos este modelo para derivar las inferencias.


Representamos a continuación en la Figura \ref(fig:galton1) las distribuciones posteriores de los efectos fijos:


```r
fixed=names(fit$marginals.fixed)
g=list()
for(i in 1:length(fixed)){
  g[[i]]=ggplot(as.data.frame(fit$marginals.fixed[[i]]),aes(x=x,y=y))+
    geom_line()+
    geom_vline(xintercept=fit$summary.fixed$mean[i],linetype="dashed")+
    geom_vline(xintercept=fit$summary.fixed[i,3],linetype="dotted")+
    geom_vline(xintercept=fit$summary.fixed[i,5],linetype="dotted")+
    labs(x=fixed[i],y="D.posterior")
}
grid.arrange(g[[1]],g[[2]],g[[3]],g[[4]],ncol=2)
```

![(\#fig:galton1)Distribución posterior de los efectos fijos](02-anova_files/figure-latex/galton1-1.pdf) 

En media observamos que la estatura de los hombres es 5.23 unidades superior a la de las mujeres.

Con estas distribuciones podemos calcular cualquier probabilidad, como por ejemplo,la probabilidad de que la estatura de un hombre, sean como sean sus ancestros, supere en 5 unidades a la de una mujer:


```r
1-inla.pmarginal(5,fit$marginals.fixed$"GenderM")
#> [1] 0.9414368
```


Podemos acceder a las distribuciones posteriores de la estatura esperada de un sujeto y posicionar las estaturas de sus padres, que se muestran en la Figura \@ref(fig:galton2)


```r
# la predicción del predictor lineal para cada sujeto es:
pred<-fit$marginals.linear.predictor
# que en este caso coincide con los valores ajustados
fitted<-fit$marginals.fitted.values
ggplot(as.data.frame(pred$Predictor.1),aes(x=x,y=y))+
  geom_line()+
  labs(x="Estatura media del sujeto 1",y="D.posterior")+
  geom_vline(xintercept=fit$summary.fitted.values$mean[1],linetype="dashed")+
  geom_vline(xintercept=datos$Father[1],linetype="dashed",color="blue")+
  geom_vline(xintercept=datos$Mother[1],linetype="dashed",color="pink")+
  annotate("text",x=datos$Mother[1]+1,y=1,label="Madre")+
  annotate("text",x=datos$Father[1]-1,y=1,label="Padre")
```

![(\#fig:galton2)Predicción de la estatura del primer sujeto en la muestra](02-anova_files/figure-latex/galton2-1.pdf) 

Podemos ir más allá, prediciendo la estatura de un sujeto, sea hombre o mujer, cuando su padre mide 1.75m (68.9 pulgadas) y su madre 1.70m (66.9 pulgadas). Expresamos los resultados en centímetros.


```r
new.data=data.frame(Father=c(68.9,68.9),
                    Mother=c(66.9,66.9),
                    Gender=c("M","F"),
                    Height=c(NA,NA))
datos.combinado <- rbind(datos, data.frame(Family=c(NA,NA),new.data,Kids=c(NA,NA)))

## creamos un vector con NA's para observaciones y 1's para predicciones
datos.indicador <- c(rep(NA, nrow(datos)), rep(1, nrow(new.data)))
## reajustamos el modelo añadiendo la opción de predicción de datos
fit.pred <- inla(formula, data = datos.combinado, 
                 control.compute=list(return.marginals.predictor=TRUE),
                 control.predictor = list(link = datos.indicador))
## y describimos los valores ajustados para los escenarios añadidos
round(fit.pred$summary.fitted.values[(nrow(datos)+1):nrow(datos.combinado),]*2.54,1)
#>                       mean  sd 0.025quant 0.5quant
#> fitted.Predictor.899 177.9 0.3      177.3    177.9
#> fitted.Predictor.900 164.7 0.3      164.0    164.7
#>                      0.975quant mode
#> fitted.Predictor.899      178.6   NA
#> fitted.Predictor.900      165.3   NA
```

También graficar las distribuciones predictivas y probabilidades (Figura \@ref(fig:galton3):


```r
# Distribuciones predictivas
pred.M=as.data.frame(fit.pred$marginals.fitted.values[[(nrow(datos)+1)]])*2.54
pred.F=as.data.frame(fit.pred$marginals.fitted.values[[(nrow(datos)+2)]])*2.54
d.pred=rbind(pred.M,pred.F)
# atributo Gender
d.pred$Gender=rep(c("M","F"),c(nrow(pred.M),nrow(pred.F)))
# objetivo de estatura
d.pred$obj=rep(c(178,165),c(nrow(pred.M),nrow(pred.F)))

ggplot(d.pred,aes(x=x,y=y))+
  geom_line()+
  geom_vline(aes(xintercept=obj),linetype="dashed")+
  facet_wrap(vars(Gender),scales="free")+
  labs(x="Estatura",y="D.posterior")
```

![(\#fig:galton3)Distribución predictiva de la estatura de un sujeto cuyo padre mide 1,75cm y madre 1,07cm.](02-anova_files/figure-latex/galton3-1.pdf) 

Y calcular probabilidades, como la probabilidad de que dicho sujeto supere el 1.65m si es mujer, o el 1.78m si es hombre.


```r
# cálculo de probabilidades
p165F=round(1-inla.pmarginal(165,pred.F),2)
cat(paste("Pr(estatura>165|mujer,padre=175,madre=170)=",p165F))
#> Pr(estatura>165|mujer,padre=175,madre=170)= 0.16
cat("\n")
p178M=round(1-inla.pmarginal(178,pred.M),2)
cat(paste("Pr(estatura>178|hombre,padre=175,madre=170)=",p178M))
#> Pr(estatura>178|hombre,padre=175,madre=170)= 0.42
```

Podríamos también, modificar las especificaciones a priori sobre los parámetros $\beta_0$ y $\beta_1$ mediante el comando `control.fixed`. Por ejemplo, queremos asumir a priori $\beta_0\sim N(0,10^4)$ y $\beta_1\sim N(0,100)$ y ver cómo afecta a las inferencias.

```r
fit<-inla(formula,family="gaussian",data=datos,
                   control.fixed=list(mean=0,prec=0.01,
                   mean.intercept=0, prec.intercept=0.0001))
round(fit$summary.fixed,3)
#>               mean    sd 0.025quant 0.5quant 0.975quant
#> (Intercept) 15.335 2.745      9.949   15.335     20.720
#> Father       0.406 0.029      0.349    0.406      0.463
#> Mother       0.322 0.031      0.260    0.322      0.383
#> GenderM      5.225 0.144      4.943    5.225      5.507
#>             mode kld
#> (Intercept)   NA   0
#> Father        NA   0
#> Mother        NA   0
#> GenderM       NA   0
```

Si queremos especificar medias a priori diferentes para los coeficientes de los distintos regresores, hemos de especificarlos con listas. 


```r
fit = inla(formula,family = "gaussian",data=datos,
                   control.fixed=list(mean=list(Father=0.2,Mother=0.1)))
round(fit$summary.fixed,3)
#>               mean    sd 0.025quant 0.5quant 0.975quant
#> (Intercept) 15.345 2.747      9.957   15.345     20.733
#> Father       0.406 0.029      0.349    0.406      0.463
#> Mother       0.321 0.031      0.260    0.321      0.383
#> GenderM      5.226 0.144      4.943    5.226      5.508
#>             mode kld
#> (Intercept)   NA   0
#> Father        NA   0
#> Mother        NA   0
#> GenderM       NA   0
```

Si queremos modificar la especificación de la prior en $\sigma^2$, o lo que es equivalente, en la precisión $\tau$, con la distribución $log(\tau) \sim N(0,1)$ en lugar de $\tau \sim Ga(1,10^{-5})$, vemos cómo afecta a la inferencia posterior sobre la precisión.


```r
fit_n = inla(formula,family="gaussian", data=datos,
                   control.family=list(hyper=list(
                     prec=list(prior="gaussian",param=c(0,1)))))
# con el modelo log-gamma para precisión
round(fit$summary.hyperpar,3)
#>                                          mean   sd
#> Precision for the Gaussian observations 0.216 0.01
#>                                         0.025quant 0.5quant
#> Precision for the Gaussian observations      0.197    0.216
#>                                         0.975quant mode
#> Precision for the Gaussian observations      0.236   NA
# con el modelo normal para precisión
round(fit_n$summary.hyperpar,3)
#>                                          mean   sd
#> Precision for the Gaussian observations 0.216 0.01
#>                                         0.025quant 0.5quant
#> Precision for the Gaussian observations      0.197    0.216
#>                                         0.975quant mode
#> Precision for the Gaussian observations      0.236   NA
```

Cuando tenemos información previa disponible sobre la variación de los datos, será generalmente más intuitivo expresarla en términos de la desviación estándar $\sigma$. Bastará con conseguir la equivalencia en la escala de $log(\tau)$ para incluirla en el modelo. Por ejemplo, si sabemos que la desviación típica está entre 2 y 14, $\sigma \sim Unif(2,14)$, podemos calcular una prior para la log-precisión del siguiente modo:

1. simular una muestra de $\sigma \sim Unif(2,14)$
1. transformar a precisiones
1. calcular los parámetros de la Gamma para la precisión, a partir de su media y varianza

Hacemos los cálculos y graficamos la prior en la Figura \@ref(fig:galton4).


```r
# parámetros para sigma~Un(a1,b1)
a1<-2
b1<-14
# simulamos sigma de una distribución Unif(a1,b1)
sigma<-runif(n=10000,min=a1,max=b1)
# obtenemos la precisión
tau<-1/sigma^2
# Calculamos los parámetros alpha,beta de una distrib. Gamma para la precisión
# mean=alpha/beta; var=alpha/beta^2
beta= mean(tau)/var(tau)
alpha<-mean(tau)*beta
# dibujamos los valores de la precisión
tau.seq=sort(tau)
  # seq(min(tau),max(tau),length=1000)
prior=data.frame(tau=tau.seq,dprior=dgamma(tau.seq,alpha,beta))
ggplot(prior, aes(x=tau))+
  geom_histogram(aes(y=..density..),color="grey",fill="white")+
  geom_line(aes(y=dprior))
#> `stat_bin()` using `bins = 30`. Pick better value with
#> `binwidth`.
```

![(\#fig:galton4)Distribución a prior para tau con sigma ~ Uniforme(2,14)](02-anova_files/figure-latex/galton4-1.pdf) 

Utilicemos pues esos parámetros para especificar la prior sobre $\tau$ en INLA:

```r
fit = inla(formula,family="gaussian",data=datos,
                   control.family=list(hyper=list(
                     prec=list(prior="loggamma",param=c(alpha,beta)))))
round(fit$summary.hyperpar,3)
#>                                          mean   sd
#> Precision for the Gaussian observations 0.214 0.01
#>                                         0.025quant 0.5quant
#> Precision for the Gaussian observations      0.195    0.214
#>                                         0.975quant mode
#> Precision for the Gaussian observations      0.234   NA
```



## Efectos aleatorios

Desde una perspectiva frecuentista un modelo básico de Anova podría ser un modelo de efectos fijos, pero también de efectos aleatorios. Así por ejemplo el tratamiento dado en un ensayo clínico importa para comparar y diferenciar el efecto que provoca sobre un paciente; tratamiento en entonces es un efecto fijo, de interés primario. En otro ejemplo, se han aplicado varios fertilizantes a cultivos en fincas distintas; el interés primario será comparar los fertilizantes, pero no las fincas, por lo que fertilizante será un efecto fijo; sin embargo, el factor finca sólo tiene interés por cuanto aporta variabilidad en la respuesta, y no para comparar las fincas, por lo que se considerará como un efecto aleatorio.

Una variable predictiva, numérica o categórica, entra en el modelo como **efecto fijo** cuando se piensa que afecta a todas las observaciones del mismo modo (de un modo lineal), y que su efecto es de interés primario en el estudio. 
En un contexto bayesiano un efecto fijo tendrá un coeficiente asociado al que se le asigna a menudo una distribución a priori vaga (mínimo informativa), como una gausiana con media cero y varianza (conocida) grande. En cualquier caso, la distribución a priori que se asume para los efectos fijos es siempre una distribución conocida.

Un **efecto aleatorio** identifica a variables de tipo categórico que no son de interés primario en la investigación, pero que se considera que añaden incertidumbre y por lo tanto variabilidad a la respuesta. La modelización habitual de los efectos aleatorios es una prior gausiana con media cero y una precisión desconocida, para la que será preciso asignar una distribución a priori. La distribución a priori de los efectos aleatorios tiene parámetros desconocidos, llamados **hiperparámetros**, a los que hay que asignar asimismo una distribución a priori.

Puesto que no salimos del modelo lineal, seguiremos asumiendo una respuesta normal, *gaussian*, con media igual a un predictor lineal $\eta=\theta+ Z u$, donde $Z$ es la correspondiente matriz de diseño para los efectos aleatorios  $z_1, z_2,...$. Se asume además una varianza desconocida $\sigma^2$.

En INLA la fórmula de predicción de una respuesta $y$ a partir de un conjunto de efectos aleatorios z1,z2,... se especifica como:


```r
formula = y ~ 1  + f(z1, model="") + f(z2,model="") 
```

donde la función $f()$ especifica la relación entre el predictor lineal de la respuesta y los efectos aleatorios $z$. El tipo de relación asumida se incluye en el argumento `model` o modelo latente, que tiene como posibilidades `names(inla.models()$latent)`, si bien en el modelo lineal la opción habitual es `model="iid"`, que asume efectos aleatorios independientes e idénticamente distribuidos. La función $f()$ tiene muchos argumentos, que se pueden consultar con el comando `?f`. 


```r
names(inla.models()$latent)
#>  [1] "linear"       "iid"          "mec"         
#>  [4] "meb"          "rgeneric"     "cgeneric"    
#>  [7] "rw1"          "rw2"          "crw2"        
#> [10] "seasonal"     "besag"        "besag2"      
#> [13] "bym"          "bym2"         "besagproper" 
#> [16] "besagproper2" "fgn"          "fgn2"        
#> [19] "ar1"          "ar1c"         "ar"          
#> [22] "ou"           "intslope"     "generic"     
#> [25] "generic0"     "generic1"     "generic2"    
#> [28] "generic3"     "spde"         "spde2"       
#> [31] "spde3"        "iid1d"        "iid2d"       
#> [34] "iid3d"        "iid4d"        "iid5d"       
#> [37] "iidkd"        "2diid"        "z"           
#> [40] "rw2d"         "rw2diid"      "slm"         
#> [43] "matern2d"     "dmatern"      "copy"        
#> [46] "clinear"      "sigm"         "revsigm"     
#> [49] "log1exp"      "logdist"
```

    
Veamos cómo ajustar un modelo de efectos aleatorios a partir de un ejemplo sencillo. Comenzamos con la base de datos `broccoli` en la librería [`faraway`](https://cran.r-project.org/web/packages/faraway/faraway.pdf). Varios cultivadores suministran brócoli a una planta de procesamiento de alimentos. La planta da instrucciones a los cultivadores para que empaquen el brócoli en cajas de tamaño estándar. Debe haber 18 racimos de brócoli por caja y cada racimo debe pesar entre 1,33 y 1,5 libras. Debido a que los productores utilizan diferentes variedades, métodos de cultivo, etc., hay cierta variación en el peso de los racimos. El responsable de la planta seleccionó 3 cultivadores al azar y luego 4 cajas al azar suministradas por estos cultivadores. Se seleccionaron 3 racimos de brocoli de cada caja.

La variable de interés es el peso del racimo de brócoli, en la variable `wt`. Sin embargo, dado cómo se ha seleccionado la muestra, el objetivo no es ni la comparación entre cultivadores (`grower`), ni entre cajas (`box`) ni entre racimos (`cluster`). Sin embargo, de manera lógica intuimos que habrá variabilidad también entre cajas (efecto aleatorio `box`) y también entre cultivadores (efecto aleatorio `grower`), lo que nos conduce a un modelo en el que todos los predictores, `box`  y `grower`  intervienen como efectos aleatorios; la variable `cluster` la aprovechamos a modo de repeticiones de medidas en una misma caja de un mismo cultivador.

La base de datos cuenta con 36 registros (3 observaciones en cada combinación `grower-box`.

$$(y_{ijk}|\eta_{ij},\sigma^2 ) \sim N(\eta_{ijk},\sigma^2)$$

con $$\eta_{ijk} = \theta + \alpha_i^G + \beta_j^B; \ \  i=2,3; j=2,3,4; k=1,2,3$$
donde $\alpha^G$ representa el efecto aleatorio asociado al cultivador y $\beta^B$ a la caja.

Así el vector de efectos latentes está compuesto por el efecto fijo de interceptación $theta$ y los efectos aleatorios $u=(\alpha_2^G,\alpha_3^G,\beta_2^B,\beta_3^B,\beta_4^B)$. 

El siguiente paso es especificar una distribución a priori sobre los parámetros. INLA por defecto asigna una prior difusa sobre la interceptación $\theta$ y también sobre la precisión de los datos $\tau=1/\sigma^2$. Dado que los $\alpha_i^G$ representan el efecto diferencial asociado al cultivador, es razonable asumir independencia entre todos estos parámetros y una distribución idéntica, centrada en el cero (ante ausencia de información) y con una varianza desconocida. Del mismo modo, se asume que los $\beta_j^B$ son a priori independientes e idénticamente distribuidos (iid) con una normal centrada en el cero (ante ausencia de información) y con varianza desconocida.

\begin{eqnarray*}
\theta &\sim & N(0,\sigma_{\theta}^2), \ \sigma_{\theta}^2=\infty \\
log(\tau) &\sim & Log-Ga(1,5\cdot 10^{-5})\\
\alpha_i^G & \sim_{iid} & N(0,\sigma_{G}^2), i=2,3 \\
\beta_j^B & \sim_{iid} & N(0,\sigma_{B}^2), j=2,3,4
\end{eqnarray*}

Surgen pues, dos nuevos parámetros en las a priori, o hiperparámetros, $\sigma_{G}^2$ y $\sigma_{B}^2$, a los que también habrá que asignar una distribución a priori. Dado que se trata de varianzas, por defecto INLA asume gammas inversas difusas, o lo que es lo mismo, log-gammas difusas para las precisiones

\begin{eqnarray*}
\tau_{G}=1/\sigma_{G}^2 &\sim & Ga(1,5\cdot 10^{-5}) \\
\tau_{B}=1/\sigma_{B}^2 &\sim & Ga(1,5\cdot 10^{-5})
\end{eqnarray*}

Surgen pues, tres niveles de especificación del modelo: datos, parámetros e hiperparámetros, que generan un modelo jerárquico, y sobre el que hablaremos más adelante.


```r
data(broccoli, package="faraway")
str(broccoli)
#> 'data.frame':	36 obs. of  4 variables:
#>  $ wt     : num  352 369 383 339 367 328 376 359 388 365 ...
#>  $ grower : Factor w/ 3 levels "1","2","3": 1 1 1 2 2 2 3 3 3 1 ...
#>  $ box    : Factor w/ 4 levels "1","2","3","4": 1 1 1 1 1 1 1 1 1 2 ...
#>  $ cluster: Factor w/ 3 levels "1","2","3": 1 2 3 1 2 3 1 2 3 1 ...
formula = wt ~ f(grower,model="iid")+ f(box,model="iid")
fit = inla(formula, family="gaussian",data=broccoli,
           control.compute = list(dic=TRUE,waic=TRUE))  
```

Cuando queremos mostrar los resultados a posteriori sobre los efectos aleatorios a partir de un ajuste `fit` con `inla`, tenemos las siguientes opciones:

-   `fit$summary.random` resume  la inferencia posterior sobre los efectos
    aleatorios
-   `names(fit$marginals.random)` lista los nombre de todos los efectos aleatorios
-   `fit$marginals.random` da las distribuciones posteriores marginales de los efectos aleatorios


```r
fit$summary.random
#> $grower
#>   ID          mean         sd  0.025quant      0.5quant
#> 1  1  1.985980e-06 0.01315983 -0.02996227  8.981825e-07
#> 2  2 -1.390190e-05 0.01315986 -0.03004579 -6.287420e-06
#> 3  3  1.191598e-05 0.01315985 -0.02991009  5.389257e-06
#>   0.975quant mode         kld
#> 1 0.02998314   NA 0.006097810
#> 2 0.02989966   NA 0.006097838
#> 3 0.03003535   NA 0.006097830
#> 
#> $box
#>   ID        mean       sd 0.025quant    0.5quant 0.975quant
#> 1  1  0.29220708 1.504071  -2.476249  0.23643442   3.563340
#> 2  2 -0.16749264 1.487971  -3.271176 -0.13580544   2.656885
#> 3  3 -0.07400084 1.481622  -3.076984 -0.06005145   2.807207
#> 4  4 -0.05063106 1.480802  -3.031491 -0.04109151   2.847049
#>   mode          kld
#> 1   NA 0.0003927971
#> 2   NA 0.0003020735
#> 3   NA 0.0002670739
#> 4   NA 0.0002626010
```

Más que la inferencia sobre los efectos aleatorios, es importante la que se hace sobre las varianzas asociadas:


```r
fit$summary.hyperpar
#>                                                 mean
#> Precision for the Gaussian observations 3.098359e-03
#> Precision for grower                    1.089611e+04
#> Precision for box                       2.613325e-01
#>                                                   sd
#> Precision for the Gaussian observations 8.762370e-04
#> Precision for grower                    2.139107e+04
#> Precision for box                       7.181807e-01
#>                                           0.025quant
#> Precision for the Gaussian observations 1.965925e-03
#> Precision for grower                    3.341883e+02
#> Precision for box                       1.381779e-02
#>                                             0.5quant
#> Precision for the Gaussian observations 2.973353e-03
#> Precision for grower                    4.852999e+03
#> Precision for box                       9.855093e-02
#>                                           0.975quant mode
#> Precision for the Gaussian observations 4.956842e-03   NA
#> Precision for grower                    5.957611e+04   NA
#> Precision for box                       1.549269e+00   NA
```


Vemos que tanto la precisión asociada al efecto aleatorio caja (`box`), como al efecto cultivador, `grower`, son muy grandes, lo que implica varianzas muy pequeñas que posiblemente nos permita prescindir de dichos efectos aleatorios para ajustar un mejor modelo.Cuando transformamos a escala de desviaciones estándar, tenemos la distribución posterior para los tres tipo de error (en la Figura \@ref(fig:brocoli2)).


```r
nombres=c("sigma","grower","box")
sigma.post=as.data.frame(inla.tmarginal(function(tau) tau^(-1/2),
  fit$marginals.hyperpar[[1]]))
sigma.grower.post =as.data.frame(inla.tmarginal(function(tau) tau^(-1/2),
  fit$marginals.hyperpar[[2]]))
sigma.box.post = as.data.frame(inla.tmarginal(function(tau) tau^(-1/2),
  fit$marginals.hyperpar[[3]]))

sigma=rbind(sigma.post,sigma.grower.post,sigma.box.post)
sigma$efecto=rep(c("sigma","grower","box"),
                 c(nrow(sigma.post),nrow(sigma.grower.post),nrow(sigma.box.post))) 

ggplot(sigma,aes(x=x,y=y)) + 
  geom_line(aes(color=efecto)) +
  labs(x=expression(sigma),y="D.Posterior")+
  facet_wrap(vars(efecto),scales = "free")+
  theme(axis.text.x = element_text(angle = 45))
```

![(\#fig:brocoli2)Distribución posterior de la desviación típica para las tres fuentes de error: datos, caja y cultivador](02-anova_files/figure-latex/brocoli2-1.pdf) 

No obstante, antes de tomar una decisión sobre la exclusión de los efectos aleatorios, vamos a hacer una aproximación del porcentaje de varianza explicada por cada una de estas fuentes de variación. Utilizando simulaciones de las distribuciones posteriores de $\sigma^2, \sigma_{G}^2$ y $\sigma_{B}^2$ vamos a calcular la contribución a la varianza del efecto cultivador, $ \sigma_{G}^2/(\sigma^2 + \sigma_{G}^2)$ y la contribución a la varianza del efecto caja, $ \sigma_{B}^2/(\sigma^2 + \sigma_{B}^2)$, y calcular con ellas un porcentaje promedio.


```r
n=1000
tau=as.data.frame(inla.hyperpar.sample(n,fit,improve.marginals=TRUE))
sigma2=apply(tau,2,function(x) 1/x)
colnames(sigma2)=c("sigma2d","sigma2G","sigma2B")
# contribución a la varianza de grower
cG= sigma2[,2]/apply(sigma2,1,sum)
# contribución a la varianza de grower
cB= sigma2[,3]/apply(sigma2,1,sum)
cat(paste("Contribución media de grower a la varianza:",round(mean(cG)*100,6),"por 100","\n"))
#> Contribución media de grower a la varianza: 0.000147 por 100
cat(paste("Contribución media de box a la varianza:",round(mean(cB)*100,6),"por 100"))
#> Contribución media de box a la varianza: 4.415927 por 100
```


Ante estos resultados, y dados los valores del DIC (305.5829015) y del WAIC (306.7749495), se justifica la opción de prescindir de los efectos `grower` y `box` como efectos aleatorios y ajustar el modelo con un único efecto fijo global. 


```r
formula = wt ~ 1
fit = inla(formula, family="gaussian",data=broccoli,
           control.compute = list(dic=TRUE,waic=TRUE))  
fit$summary.fixed
#>                 mean      sd 0.025quant 0.5quant 0.975quant
#> (Intercept) 358.1667 2.79132   352.6611 358.1667   363.6722
#>             mode          kld
#> (Intercept)   NA 1.246858e-08
fit$summary.hyperpar
#>                                                mean
#> Precision for the Gaussian observations 0.003757592
#>                                                   sd
#> Precision for the Gaussian observations 0.0008836951
#>                                          0.025quant
#> Precision for the Gaussian observations 0.002258693
#>                                           0.5quant
#> Precision for the Gaussian observations 0.00367462
#>                                          0.975quant mode
#> Precision for the Gaussian observations 0.005689411   NA
```

Vemos que la variación en los indicadores DIC (308.1333981) y WAIC (307.8685218) es despreciable para este nuevo modelo.

Inferimos a continuación con las posteriores para la media global y la varianza de los datos.


```r
theta.post = as.data.frame(fit$marginals.fixed[[1]])
sigma.post=as.data.frame(inla.tmarginal(function(tau) tau^(-1/2),
  fit$marginals.hyperpar[[1]]))

posterior=rbind(theta.post,sigma.post)
posterior$efecto=rep(c("theta","sigma"),
                     c(nrow(theta.post),nrow(sigma.post)))
ggplot(posterior,aes(x=x,y=y)) + 
  geom_line(aes(color=efecto)) +
  labs(x="",y="D.Posterior")+
  facet_wrap(vars(efecto),scales = "free")+
  theme(legend.position = "none")
```

![](02-anova_files/figure-latex/unnamed-chunk-25-1.pdf)<!-- --> 

## Modelos mixtos

En ocasiones cuando ajustamos un modelo lineal tendremos algunos factores de clasificación que operan como efectos fijos y otros que operan como efectos aleatorios. Estaremos ante **modelos lineales mixtos**. Siendo estrictos, realmente el modelo con solo efectos aleatorios ya es un modelo mixto, puesto que incluye como efecto fijo una interceptación global.

En un modelo lineal mixto seguimos asumiendo una respuesta normal, *gaussian*, con media igual a un predictor lineal $\eta=X\beta + Z u$, donde $X$ es una matriz de diseño con los efectos fijos $x_1,x_2,...$, y $Z$ la correspondiente para los efectos aleatorios  $z_1, z_2,...$. Se asume además una varianza desconocida que puede ser distinta para distintos niveles de los predictores, y que en general se suele expresar a través de una matriz de covarianzas $\Sigma$, $(y|\eta,\Sigma) \sim N(\eta,\Sigma)$.

En INLA la fórmula de predicción de una respuesta $y$ a partir de un conjunto de efectos fijos x1,x2,..., y un conjunto de efectos aleatorios z1,z2,... se especifica como:


```r
formula = y ~ 1 + x1 + x2  + f(z1, model="") + f(z2,model="") 
```

De nuevo mencionar que la opción más habitual para los efectos aleatorios en un modelo lineal es `model="iid"`. 


Veamos cómo resolver las inferencias a través de un ejemplo disponible en [R-bloggers](https://www.r-bloggers.com/2019/09/bayesian-linear-mixed-models-random-intercepts-slopes-and-missing-data/), proporcionado por [Patrick Curran](https://curran.web.unc.edu/) y descargable desde [Github](https://github.com/MultiLevelAnalysis/Datasets-third-edition-Multilevel-book/blob/master/chapter%205/Curran/CurranLong.sav). Se refieren estos datos, a un estudio con 405 niños en los dos primeros años de la escuela infantil, medidos a lo largo de cuatro instantes equidistantes (no medidos todos en todos los sujetos) para registrar su progreso en lectura y en comportamiento antisocial. Nos centramos aquí exclusivamente en intentar predecir los progresos en lectura (variable `read`) a partir de su estimulación cognitiva en casa (`homecog`), teniendo en cuenta la existencia de medidas repetidas para cada sujeto (identificado como `id`) en los 4 instantes de medición (`occasion`).

Cargamos los datos y los inspeccionamos en la Figura \@ref(fig:curran1).


```r
# librería para leer archivos .sav
library(haven)
my.dir="~/Dropbox/ESTADISTICA/BAYESIAN/VARIOS/"
curran_dat = read_sav(paste0(my.dir,"CurranLong.sav")) %>%
  select(id, occasion, read, homecog) %>%
  filter(complete.cases(.))
curran_dat$id=as.factor(curran_dat$id)
# 'occasion' la convertimos en factor y la preservamos como numérica en 'time'
curran_dat$time=curran_dat$occasion
curran_dat$occasion=as.factor(curran_dat$occasion)
# Relaciones
g1=ggplot(curran_dat, aes(x=occasion,y=read))+
  geom_boxplot()
g2=ggplot(curran_dat, aes(x=homecog,y=read))+
  geom_point()+
  facet_wrap(vars(occasion))
grid.arrange(g1,g2,ncol=2)
```

![(\#fig:curran1)Descripción de la BD CurranLong sobre desarrollo de las habilidades lectoras en niños.](02-anova_files/figure-latex/curran1-1.pdf) 

Como base vamos a asumir normalidad en la respuesta de un sujeto $i$ en un instante $j$, y plantear un modelo lineal para obtener nuestras conclusiones.
$$( y_{ij}|\eta_{ij},\sigma^2 ) \sim N(\eta_{ij},\sigma^2);  \ i=1,...,450; j=1,2,3,4$$

A continuación, planteamos distintas alternativas de modelización que dan lugar a diversos predictores lineales. En los primeros modelos prescindimos de momento del efecto de la estimulación cognitiva.

### M0: interceptaciones distintas por sujetos

Solo estamos interesados en cuantificar el nivel general de habilidades lectoras, pero teniendo en cuenta posibles diferencias entre los niños. No nos interesan sin embargo, dichas diferencias. Pensamos pues en utilizar la variable `id` como efecto aleatorio y modelizar el predictor lineal con:
$$\eta_{ij}=\eta_i=\theta + \alpha_i^{id}$$
donde $\alpha_i^{id} \sim N(0,\sigma_{\alpha}^2)$ y a priori $\tau_{\alpha} \sim Ga(0.001,0.001)$.

El modelo que ajustamos genera la inferencia posterior para un efecto fijo de interacción, $\theta$, y dos varianzas (precisiones) que explican la variabilidad existente: debida a los datos $\sigma^2$, y debida a la variabilidad entre los sujetos $\sigma_{\alpha}^2$. En la Figura \@ref(fig:curranm0) se muestran las distribuciones posteriores obtenidas sobre efectos fijos y varianzas.


```r
prec.prior=list(prec=list(param=c(0.001,0.001)))
formula0= read ~ f(id,model="iid",hyper = prec.prior) 
fit0=inla(formula0,family="gaussian",data=curran_dat)
fit=fit0
fit$summary.fixed
#>                 mean         sd 0.025quant 0.5quant
#> (Intercept) 4.114053 0.05109306   4.013248 4.114248
#>             0.975quant mode          kld
#> (Intercept)   4.213756   NA 1.440557e-09
fit$summary.hyperpar
#>                                              mean
#> Precision for the Gaussian observations 0.4182785
#> Precision for id                        3.7326816
#>                                                 sd
#> Precision for the Gaussian observations 0.01924282
#> Precision for id                        1.07716432
#>                                         0.025quant
#> Precision for the Gaussian observations  0.3810323
#> Precision for id                         2.1871261
#>                                          0.5quant
#> Precision for the Gaussian observations 0.4180809
#> Precision for id                        3.5412659
#>                                         0.975quant mode
#> Precision for the Gaussian observations  0.4568346   NA
#> Precision for id                         6.3715903   NA

nfixed=length(names(fit$marginals.fixed))
nhyp=length(names(fit$marginals.hyperpar))
res=NULL
for(i in 1:nfixed){
res=rbind(res,data.frame(fit$marginals.fixed[[i]],
                         id=names(fit$marginals.fixed)[i],
                          tipo="fixed"))
}
for(j in 1:nhyp){
  res=rbind(res,data.frame(fit$marginals.hyperpar[[j]],
                           id=names(fit$marginals.hyperpar)[j],
                            tipo="prec"))
}
ggplot(res,aes(x=x,y=y))+
  geom_line(aes(color=id))+
  facet_wrap(vars(tipo),scales="free")+
  theme(legend.position="top",
        legend.title=element_blank(),
        legend.text = element_text(size=5))
```

![(\#fig:curranm0)Distribuciones posteriores para el modelo con interceptaciones aleatorias por sujeto.](02-anova_files/figure-latex/curranm0-1.pdf) 


La matriz de diseño $Z$ asociada a los efectos aleatorios se obtiene con la función `model.matrix`, y la podemos utilizar también para ajustar el modelo con `inla`, sustituyendo la modelización `model="iid"` por `model="z"` en la especificación de los efectos aleatorios. Para ello habremos de definir un nuevo índice para todos los registros de la base de datos, y aplicar sobre ellos la matriz de efectos aleatorios:


```r
Z <- as(model.matrix(~ 0 + id, data = curran_dat), "Matrix")
# índice nuevo
ID=1:nrow(curran_dat)
formula00= read ~ f(ID,model="z",Z=Z,hyper = prec.prior) 
fit00=inla(formula00,family="gaussian",data=curran_dat)
#> Warning in inla.model.properties.generic(inla.trim.family(model), mm[names(mm) == : Model 'z' in section 'latent' is marked as 'experimental'; changes may appear at any time.
#>   Use this model with extra care!!! Further warnings are disabled.
fit00$summary.hyperpar
#>                                              mean
#> Precision for the Gaussian observations 0.4180811
#> Precision for ID                        3.7326069
#>                                                 sd
#> Precision for the Gaussian observations 0.01955321
#> Precision for ID                        1.07788425
#>                                         0.025quant
#> Precision for the Gaussian observations  0.3801053
#> Precision for ID                         2.1934543
#>                                          0.5quant
#> Precision for the Gaussian observations 0.4179313
#> Precision for ID                        3.5390516
#>                                         0.975quant mode
#> Precision for the Gaussian observations  0.4571309   NA
#> Precision for ID                         6.3797445   NA
```
La diferencia entre ajustar los efectos aleatorios con "iid" o con "z" es que con la primera opción tendremos tantos efectos aleatorios como están definidos en el modelo. Sin embargo, la modelización con "z" producirá tantos efectos aleatorios como datos, con una única precisión/varianza asociada. Generalmente ambos modelos producen unas medias posterioris de los efectos aleatorios bastante parecidas, como se muestra en la Figura \@ref(fig:curranm0b). 


```r
# medias posteriori de los efectos aleatorias con model="iid"
random0=fit0$summary.random$id$mean
# medias posteriori de los efectos aleatorios con model="z"
random00=fit00$summary.random$ID$mean
plot(random00,random00,type="l")
points(random0,random0)
```

![(\#fig:curranm0b)Concordancia entre las medias posteriores de los efectos aleatorios para el modelo iid y el especificado con la matriz de diseño para los efectos aleatorios.](02-anova_files/figure-latex/curranm0b-1.pdf) 


### M1: interceptaciones distintas por sujetos e instantes

Dado que tenemos varias mediciones de un mismo sujeto, es razonable asumir que, por defecto, las habilidades cognitivas propias de un sujeto, $\alpha_i^{id}$, influyen en sus habilidades lectoras. 

$$\eta_{ij}=\theta + \alpha_i^{id}$$

Además, es de esperar que el tiempo que transcurre afecte de modo similar a la evolución de todos los sujetos (en la Figura \@ref(fig:curran1) se apreciaba cierto crecimiento). Hablamos pues de un modelo en el que predecimos las habilidades lectoras con un efecto fijo común $\theta$ afectado de cierta variación extra en función del sujeto $i$ y del instante de medición $j$.

$$\eta_{ij}=\theta + \alpha_i^{id} + \beta_j^{oc}$$
Si no nos interesa evaluar las habilidades cognitivas propias de cada sujeto, pero queremos reconocer de algún modo la variabilidad entre sujetos, $\sigma_{\alpha}^2$, estamos pensando en unos efectos aleatorios para el sujeto, `id`, $\alpha_i^{id} \sim N(0,\sigma_{\alpha}^2)$.

Si no nos interesa evaluar el efecto del tiempo sobre las habilidades lectoras, pero queremos reconocer la variabilidad entre los distintos periodos de tiempo, $\sigma_{\beta}^2$, estamos pensando en unos efectos aleatorios asociados al tiempo, `occasion`, $\beta_i^{oc} \sim N(0,\sigma_{\beta}^2)$.

Si asumimos una prior para las precisiones $\tau_{\alpha}$ y $\tau_{\beta}$ de media 1 y varianza 1000, podemos usar una $Ga(0.001,0.001)$.

con 
\begin{eqnarray*}
\text{Nivel II} && \\
\theta &\sim & N(0,100) \\
\alpha_i^{id} &\sim& N(0,\sigma_{\alpha}^2) \\
\beta_j^{oc} &\sim& N(0,\sigma_{\beta}) \\
\text{Nivel III} && \\
\tau_{\alpha}=1/\sigma_{\alpha}^2 &\sim& Ga(0.001,0.001) \\
\tau_{\beta}=1/\sigma_{\beta}^2 &\sim&Ga(0.001,0.001)
\end{eqnarray*}

En la Figura \@ref(fig:curran2) se muestran las distribuciones posteriores obtenidas sobre efectos fijos y varianzas.


```r
prec.prior=list(prec=list(param=c(0.001,0.001)))
formula1= read ~ f(id,model="iid",hyper = prec.prior) + f(occasion,model="iid",hyper = prec.prior)
fit1=inla(formula1,family="gaussian",data=curran_dat)
fit=fit1
fit$summary.fixed
#>                 mean     sd 0.025quant 0.5quant 0.975quant
#> (Intercept) 4.353429 0.8501   2.582487 4.353381   6.124708
#>             mode          kld
#> (Intercept)   NA 0.0004859045
fit$summary.hyperpar
#>                                              mean        sd
#> Precision for the Gaussian observations 2.4704987 0.1102818
#> Precision for id                        1.2807519 0.1012373
#> Precision for occasion                  0.4165207 0.3047343
#>                                         0.025quant
#> Precision for the Gaussian observations 2.26386738
#> Precision for id                        1.09977732
#> Precision for occasion                  0.04845845
#>                                          0.5quant
#> Precision for the Gaussian observations 2.4666179
#> Precision for id                        1.2741962
#> Precision for occasion                  0.3354912
#>                                         0.975quant mode
#> Precision for the Gaussian observations   2.698610   NA
#> Precision for id                          1.498455   NA
#> Precision for occasion                    1.190375   NA

nfixed=length(names(fit$marginals.fixed))
nhyp=length(names(fit$marginals.hyperpar))
res=NULL
for(i in 1:nfixed){
res=rbind(res,data.frame(fit$marginals.fixed[[i]],
                         id=names(fit$marginals.fixed)[i],
                          tipo="fixed"))
}
for(j in 1:nhyp){
  res=rbind(res,data.frame(fit$marginals.hyperpar[[j]],
                           id=names(fit$marginals.hyperpar)[j],
                            tipo="prec"))
}
ggplot(res,aes(x=x,y=y))+
  geom_line(aes(color=id))+
  facet_wrap(vars(tipo),scales="free")+
  theme(legend.position="top",
        legend.title=element_blank(),
        legend.text = element_text(size=5))
```

![(\#fig:curran2)Distribuciones posteriores para el modelo con interceptaciones aleatorias por sujeto y tiempo.](02-anova_files/figure-latex/curran2-1.pdf) 


### M2: interceptaciones distintas por sujeto y pendiente común

Por otro lado, y en base a la Figura \@ref(fig:curran3) podríamos considerar el efecto del tiempo que transcurre desde el inicio del estudio (`t=time`), como una covariable numérica que afecta de modo lineal a las habilidades lectoras.



```r
ggplot(curran_dat, aes(x=occasion,y=read))+
  geom_line(aes(group=id),color="grey",size=0.4)
```

![(\#fig:curran3)Relación entre el tiempo y las habilidades lectoras para cada sujeto (líneas).](02-anova_files/figure-latex/curran3-1.pdf) 

Podríamos seguir planteando un efecto aleatorio del sujeto sobre sus resultados lectores, y un efecto fijo asociado al tiempo transcurrido hasta el instante $t_j=j$, esto es, una interceptación aleatoria y una pendiente fija.

$$\eta_{ij}=\theta + \alpha_i^{id} + \beta \cdot t_{ij} $$
con 
\begin{eqnarray*}
\text{Nivel II} && \\
\theta &\sim & N(0,100) \\
\beta &\sim& N(0,100) \\
\alpha_i^{id} &\sim& N(0,\sigma_{\alpha}^2) \\
\text{Nivel III} && \\
\tau_{\alpha}=1/\sigma_{\alpha}^2 &\sim& Ga(0.001,0.001) 
\end{eqnarray*}

En la Figura \@ref(fig:curran4) se muestran las distribuciones posteriores obtenidas sobre efectos fijos y varianzas.


```r
prec.prior=list(prec=list(param=c(0.001,0.001)))
formula2= read ~ time + f(id,model="iid",hyper = prec.prior) 
fit2=inla(formula2,family="gaussian",data=curran_dat)
fit=fit2
fit$summary.fixed
#>                 mean         sd 0.025quant 0.5quant
#> (Intercept) 2.703744 0.05265822   2.600420 2.703750
#> time        1.101341 0.01758964   1.066827 1.101345
#>             0.975quant mode          kld
#> (Intercept)   2.807034   NA 1.225182e-11
#> time          1.135834   NA 5.561958e-12
fit$summary.hyperpar
#>                                             mean        sd
#> Precision for the Gaussian observations 2.174387 0.1013879
#> Precision for id                        1.285625 0.1091213
#>                                         0.025quant 0.5quant
#> Precision for the Gaussian observations   1.980400 2.172463
#> Precision for id                          1.083291 1.281336
#>                                         0.975quant mode
#> Precision for the Gaussian observations   2.379789   NA
#> Precision for id                          1.512927   NA

nfixed=length(names(fit$marginals.fixed))
nhyp=length(names(fit$marginals.hyperpar))
res=NULL
for(i in 1:nfixed){
res=rbind(res,data.frame(fit$marginals.fixed[[i]],
                         id=names(fit$marginals.fixed)[i],
                          tipo="fixed"))
}
for(j in 1:nhyp){
  res=rbind(res,data.frame(fit$marginals.hyperpar[[j]],
                         id=names(fit$marginals.hyperpar)[j],
                          tipo="prec"))
}
ggplot(res,aes(x=x,y=y))+
  geom_line(aes(color=id))+
  facet_wrap(vars(tipo),scales="free")+
  theme(legend.position="top",
        legend.title=element_blank(),
        legend.text = element_text(size=5))
```

![(\#fig:curran4)Distribuciones posteriores para el modelo con interceptaciones aleatorias por sujeto y efecto fijo del tiempo.](02-anova_files/figure-latex/curran4-1.pdf) 

## Efectos anidados

Hablamos de efectos anidados cuando cada miembro de un grupo está contenido completamente dentro de una única unidad de otro grupo. Que un factor A esté anidado en otro B, implica que cada nivel de B contiene niveles distintos de A, esto es, cada nivel de A está vinculado solo a algún nivel de B.

La base de datos `eggs` en la librería [`faraway`](https://cran.r-project.org/web/packages/faraway/faraway.pdf) nos resulta útil para describir este tipo de modelos con efectos anidados. Estos datos son los resultantes de un experimento para testar la consistencia en los tests de laboratorio que realizaban laboratorios distintos, técnicos distintos. Para ello se dividió en varias muestras un tarro de polvo de huevo seco homogeneizado (con idéntico contenido graso). Se enviaron 4 muestras a cada uno de los 6 laboratorios. De esas 4 muestras, 2 se etiquetaron como G y 2 como H (aun siendo idénticas).  Se dieron instrucciones a los laboratorios de dar dos muestras a dos técnicos distintos. Los técnicos recibieron instrucciones de dividir sus muestras en dos partes y medir el contenido graso de cada una. Así, cada laboratorio reportó 8 mediciones del contenido graso (`Fat`), cada técnico 4 mediciones, con 2 réplicas en cada una de las dos muestras.


```r
data(eggs,package="faraway")
```

Tenemos así en este ejemplo, a los técnicos (`Technician`) anidados en los laboratorios (`Lab`). En la Figura \@ref(fig:eggs1) se muestra claramente la variación entre laboratorios, entre técnicos, y debida al efecto irreal de tener dos muestras distintas G y H (`Sample`).


```r
ggplot(eggs,aes(x=Lab,y=Fat))+
  geom_boxplot(aes(color=Technician))+
  facet_wrap(vars(Sample))
```

![(\#fig:eggs1)Descripción de eggs: variación entre laboratorios y técnicos.](02-anova_files/figure-latex/eggs1-1.pdf) 

El modelo que planteamos para estimar la respuesta $y_{ijk}$, contenido graso de la muestra $k$ ($k=1,2$) del laboratorio $i$ ($i=1,...,6$), por el técnico $j$ ($j=1,2$) está basado como siempre, en el modelo normal, $(y_{ijk}\mu_{ijk},\sigma^2) \sim N(\mu_{ijk},\sigma^2)$, con una media o predictor lineal representado por:

$$ \mu_{ijk}= \theta + \alpha_i^{lab} + \beta_{j:i}^{tec} + \gamma_{k:(j:i)}^{sam}$$
y asumiendo las distribuciones a priori
\begin{eqnarray*}
\theta &\sim& N(0,1000) \\
\tau=1/\sigma^2 &\sim& Ga(0.001,0.001) \\
\alpha_i^{lab}&\sim& N(0,\sigma_{lab}^2); \  i = 1,...,4 \\
\beta_{j:i}^{tec}&\sim& N(0,\sigma_{tec}^2); \  j:i = 1,...,12 \\
\gamma_{k:(j:i)}^{sam}&\sim& N(0,\sigma_{sam}^2); \ k:(j:i)=1,...,24
\end{eqnarray*}

Para especificar en INLA los efectos anidados hemos de recurrir a la matriz del modelo, `model.matrix()`, para crear las matrices de los efectos aleatorios anidados.

A continuación hemos de crear los correspondientes índices, de longitud similar a la del número de registros en la base de datos, para aplicarles las correspondientes matrices de efectos anidados, y ya proceder con el ajuste como habitualmente hacemos.


```r
# matrices de efectos aleatorios anidados
Zlt <- as(model.matrix( ~ 0 + Lab:Technician, data = eggs), "Matrix")
Zlts <- as(model.matrix( ~ 0 + Lab:Technician:Sample, data = eggs), "Matrix")

# índices para aplicar los efectos aleatorios
eggs$IDt = eggs$IDts = 1:nrow(eggs)

# Ajuste
prec.prior=list(prec=list(param=c(0.001,0.001)))
formula = Fat ~ 1 + f(Lab,model="iid",hyper=prec.prior) +
  f(IDt,model="z",Z=Zlt,hyper=prec.prior)+
  f(IDts,model="z",Z=Zlts,hyper=prec.prior)

fit <- inla(formula,data = eggs, 
            control.predictor = list(compute = TRUE), 
            control.family=list(hyper=prec.prior),
            control.fixed=list(prec.intercept=0.001))
round(fit$summary.fixed,4)
#>               mean     sd 0.025quant 0.5quant 0.975quant
#> (Intercept) 0.3875 0.0554     0.2756   0.3875     0.4994
#>             mode kld
#> (Intercept)   NA   0
round(fit$summary.hyperpar,4)
#>                                             mean       sd
#> Precision for the Gaussian observations 142.0216  39.6765
#> Precision for Lab                       349.8811 651.0716
#> Precision for IDt                       206.0598 210.5203
#> Precision for IDts                      366.9448 335.2806
#>                                         0.025quant 0.5quant
#> Precision for the Gaussian observations    77.8346 137.4878
#> Precision for Lab                          22.4909 170.1948
#> Precision for IDt                          26.6541 144.1514
#> Precision for IDts                         59.8953 270.7312
#>                                         0.975quant mode
#> Precision for the Gaussian observations   232.9088   NA
#> Precision for Lab                        1808.4478   NA
#> Precision for IDt                         762.2898   NA
#> Precision for IDts                       1256.5560   NA
```

Alternativamente podríamos crear una variable índice a partir de las matrices de efectos aleatorios, para utilizarlas con `model="iid"` para describir los efectos aleatorios:


```r
eggs$labtech <- as.factor(apply(Zlt, 1, function(x){names(x)[x == 1]}))
eggs$labtechsamp <- as.factor(apply(Zlts, 1, function(x){names(x)[x == 1]}))

formula=Fat ~ 1 + f(Lab, model = "iid", hyper = prec.prior) +
    f(labtech, model = "iid", hyper = prec.prior) +
  f(labtechsamp, model = "iid", hyper = prec.prior)
fit=inla(formula, data = eggs, 
            control.predictor = list(compute = TRUE), 
            control.family=list(hyper=prec.prior),
         control.fixed=list(prec.intercept=0.001))
round(fit$summary.fixed,4)
#>               mean     sd 0.025quant 0.5quant 0.975quant
#> (Intercept) 0.3875 0.0554     0.2756   0.3875     0.4994
#>             mode kld
#> (Intercept)   NA   0
round(fit$summary.hyperpar,4)
#>                                             mean       sd
#> Precision for the Gaussian observations 141.9188  39.6195
#> Precision for Lab                       349.8671 651.1093
#> Precision for labtech                   206.0269 210.4740
#> Precision for labtechsamp               366.9583 335.2964
#>                                         0.025quant 0.5quant
#> Precision for the Gaussian observations    77.8090 137.3964
#> Precision for Lab                          22.4874 170.1767
#> Precision for labtech                      26.6550 144.1327
#> Precision for labtechsamp                  59.8972 270.7399
#>                                         0.975quant mode
#> Precision for the Gaussian observations   232.6619   NA
#> Precision for Lab                        1808.4558   NA
#> Precision for labtech                     762.1349   NA
#> Precision for labtechsamp                1256.6104   NA
```

El hecho de que las estimaciones de la precisión sean tan altas podría indicar que estos parámetros están identificados de modo pobre en el modelo y pueden requerir del uso de priors menos vagas.


## Datos longitudinales

Belenky et al. (2003) describen un estudio de los tiempos de reacción en pacientes a los que se ha privado de sueño durante 10 días; cada día se ha ido registrando la respuesta para cada uno de los 18 sujetos en el estudio. Los datos están disponibles como `sleepstudy` en la librería `lme4` y tienen como variables el tiempo medio de reacción en microsegundos (`Reaction`), el número de días con privación de sueño (`Days`) y un id para cada sujeto (`Subject`). Los tiempos de reacción se transforman a segundos para tener mayor estabilidad. Aun así, en la Figura \@ref(fig:sleep1) se aprecia que el número de días de falta de sueño afecta de modo distinto a cada sujeto.


```r
data(sleepstudy,package="lme4")
sleepstudy$Reaction <- sleepstudy$Reaction / 1000
ggplot(sleepstudy,aes(x=Days,y=Reaction))+
  geom_point(size=0.5)+
  geom_smooth(method="lm",color="blue",size=0.5)+
  facet_wrap(vars(Subject),ncol=6)+
  theme(axis.text.x = element_text(size=5),
        axis.text.y = element_text(size=5))
#> `geom_smooth()` using formula 'y ~ x'
```

![(\#fig:sleep1)Tiempos de reacción en función del número de días con falta de sueño (sleepstudy) para los 18 sujetos en el estudio](02-anova_files/figure-latex/sleep1-1.pdf) 

Un modelo razonable para estos datos es un modelo lineal que relacione los tiempos de reacción con los días, pero que tenga interceptaciones y pendientes diferentes para cada sujeto. El efecto sujeto entraría en el modelo como un efecto aleatorio y para relacionar todos los datos del mismo sujeto sin perder la asunción de independencia entre las observaciones de sujetos distintos. Si llamamos $y=Reaction$, estaríamos planteando el siguiente modelo:

$$ y_{ij}|\mu_{ij},\sigma^2 \sim N(\mu,\sigma^2), i=1,...,18; j=1, ...,10$$

con 

$$\mu_{ij}=\theta + \alpha_i + \beta \cdot x_{ij} + \gamma_{ij}$$

donde $x=Days$, ($theta,\beta$) se tratarían como efectos fijos con a prioris difusas ante falta de información, y $(\alpha_i,\gamma_{ij})$ como efectos aleatorios, con normales centradas en cero y una varianza desconocida a la que habría que asignar así mismo, una distribución a priori.

con 
\begin{eqnarray*}
\text{Nivel II} && \\
\theta &\sim & N(0,1000) \\
\beta &\sim& N(0,1000) \\
\alpha_i &\sim& N(0,\sigma_{\alpha}^2) \\
\gamma_{ij} &\sim & N(0,\sigma_{\gamma})^2 \\
\tau=1/\sigma^2 &\sim& Ga(0.001,0.001)\\
\text{Nivel III} && \\
\tau_{\alpha}=1/\sigma_{\alpha}^2 &\sim& Ga(0.001,0.001) \\
\tau_{\gamma}=1/\sigma_{\alpha}^2 &\sim& Ga(0.001,0.001) \end{eqnarray*}

En INLA modelizamos interceptaciones y pendientes distintas y vinculadas a un efecto aleatorio, incluyendo la interacción entre la covariable y el efecto aleatorio, como un efecto aleatorio en sí mismo, y asumiendo un modelo `iid` como habitualmente. A la hora de especificar tal interacción, es preciso ubicar en primer lugar el efecto aleatorio y detrás la covariable. Así la primera variable define el número de grupos y la segunda el valor de la covariable.


```r
prec.prior=list(prec=list(param=c(0.001,0.001)))
formula = Reaction ~ 1 + f(Subject, Days, model = "iid",hyper=prec.prior)
fit <- inla(formula,data = sleepstudy, 
            control.predictor = list(compute = TRUE), 
            control.family=list(hyper=prec.prior),
            control.fixed=list(prec=0.001,prec.intercept=0.001))
round(fit$summary.fixed,4)
#>               mean     sd 0.025quant 0.5quant 0.975quant
#> (Intercept) 0.2532 0.0041     0.2452   0.2532     0.2612
#>             mode kld
#> (Intercept)   NA   0
round(fit$summary.hyperpar,4)
#>                                             mean        sd
#> Precision for the Gaussian observations 1170.449  130.7965
#> Precision for Subject                   3732.841 1271.1589
#>                                         0.025quant 0.5quant
#> Precision for the Gaussian observations   930.7816 1164.475
#> Precision for Subject                    1757.6093 3563.494
#>                                         0.975quant mode
#> Precision for the Gaussian observations   1445.477   NA
#> Precision for Subject                     6704.976   NA
```

En la Figura \@ref(fig:sleep2) mostramos los datos y también los valores ajustados para las rectas, en términos de las interceptaciones y pendientes medias de las correspondientes distribuciones posteriores, además de la banda de estimación que construimos con los correspondientes percentiles de las posterioris.


```r
sleepstudy.pred = sleepstudy %>%
  mutate(fitted=fit$summary.fitted.values$mean,
         hpd.low=fit$summary.fitted.values$"0.025quant",
         hpd.up=fit$summary.fitted.values$"0.975quant") 

ggplot(sleepstudy.pred,aes(x=Days,y=Reaction))+
  geom_point(size=0.5)+
  geom_line(aes(y=fitted),color="blue")+
  geom_line(aes(y= hpd.low),color="skyblue")+
  geom_line(aes(y=hpd.up),color="skyblue")+
  facet_wrap(vars(Subject),ncol=6)+
  theme(axis.text.x = element_text(size=5),
        axis.text.y = element_text(size=5))
```

![](02-anova_files/figure-latex/unnamed-chunk-31-1.pdf)<!-- --> 



## Conclusiones

Hasta aquí desarrollamos los modelos lineales basados en Anova, esto es, en la integración de factores de clasificación como variables que van a explicar diferencias en la respuesta, como los efectos fijos, o variabilidad extra en los datos, como los efectos aleatorios, a veces incluso con otros predictores de tipo numérico, e incluso interaccionando con ellos.
