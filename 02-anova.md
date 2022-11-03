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
#> (Intercept) 61.016 1.173     58.698   61.016     63.339
#> dietB        4.979 1.514      1.980    4.980      7.973
#> dietC        6.977 1.514      3.978    6.978      9.971
#> dietD       -0.016 1.437     -2.861   -0.016      2.824
#>             mode kld
#> (Intercept)   NA   0
#> dietB         NA   0
#> dietC         NA   0
#> dietD         NA   0
tau=round(fit$summary.hyperpar,3);tau
#>                                          mean   sd
#> Precision for the Gaussian observations 0.197 0.06
#>                                         0.025quant 0.5quant
#> Precision for the Gaussian observations      0.099     0.19
#>                                         0.975quant mode
#> Precision for the Gaussian observations      0.323   NA
medias=round(fit$summary.linear.predictor,4)
```
Atendiendo a los descriptivos de la distribución posterior para los efectos fijos, concluimos:

- El tiempo esperado de coagulación para los animales que han seguido la dieta A es de 61.016(58.698,63.339).
- Los animales que han seguido la dieta B tienen un tiempo de coagulación esperado  superior en 4.979 unidades a los de la dieta A, y dicha diferencia es *significativamente distinta de cero* en el contexto bayesiano, dado que su HPD no incluye al cero, (1.98,7.973). 
- Una conclusión similar se deriva para la dieta C, que da un tiempo de coagulación esperado  superior en 6.977 unidades a los de la dieta A, y un HPD (3.978,9.971).
- Las diferencias en los tiempos de coagulación de seguir una dieta D frente a la dieta A no son relevantes. De hecho, la diferencia entre ellos es de -0.016 y el intervalo HPD contiene al cero, (-2.861,2.824).


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
#> [1] "E(sigma.post)= 2.3368 HPD95%=( 1.687 , 3.064 )"
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
#> (Intercept)                      3.96605050 0.1314191
#> BreedCanadian                    0.52193553 0.1858571
#> BreedGuernsey                    0.93293190 0.1858571
#> BreedHolstein-Fresian           -0.30304829 0.1858571
#> BreedJersey                      1.16693161 0.1858571
#> AgeMature                        0.18793906 0.1858497
#> BreedCanadian:AgeMature         -0.28692013 0.2628364
#> BreedGuernsey:AgeMature         -0.08591997 0.2628364
#> BreedHolstein-Fresian:AgeMature -0.17493825 0.2628364
#> BreedJersey:AgeMature            0.13107657 0.2628364
#>                                 0.025quant    0.5quant
#> (Intercept)                      3.7077533  3.96604997
#> BreedCanadian                    0.1566392  0.52193621
#> BreedGuernsey                    0.5676355  0.93293262
#> BreedHolstein-Fresian           -0.6683441 -0.30304778
#> BreedJersey                      0.8016352  1.16693233
#> AgeMature                       -0.1773425  0.18793970
#> BreedCanadian:AgeMature         -0.8035116 -0.28692097
#> BreedGuernsey:AgeMature         -0.6025115 -0.08592082
#> BreedHolstein-Fresian:AgeMature -0.6915303 -0.17493890
#> BreedJersey:AgeMature           -0.3855150  0.13107577
#>                                 0.975quant mode
#> (Intercept)                     4.22435069   NA
#> BreedCanadian                   0.88722796   NA
#> BreedGuernsey                   1.29822422   NA
#> BreedHolstein-Fresian           0.06224463   NA
#> BreedJersey                     1.53222392   NA
#> AgeMature                       0.55321697   NA
#> BreedCanadian:AgeMature         0.22967617   NA
#> BreedGuernsey:AgeMature         0.43067633   NA
#> BreedHolstein-Fresian:AgeMature 0.34165750   NA
#> BreedJersey:AgeMature           0.64767277   NA
#>                                          kld
#> (Intercept)                     2.743621e-09
#> BreedCanadian                   2.743955e-09
#> BreedGuernsey                   2.743953e-09
#> BreedHolstein-Fresian           2.743954e-09
#> BreedJersey                     2.743957e-09
#> AgeMature                       2.743117e-09
#> BreedCanadian:AgeMature         2.743536e-09
#> BreedGuernsey:AgeMature         2.743536e-09
#> BreedHolstein-Fresian:AgeMature 2.743535e-09
#> BreedJersey:AgeMature           2.743536e-09
fit$dic$dic
#> [1] 120.486
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
#> (Intercept)            4.0077184 0.10124855  3.80873711
#> BreedCanadian          0.3784787 0.13071130  0.12159375
#> BreedGuernsey          0.8899744 0.13071130  0.63308926
#> BreedHolstein-Fresian -0.3905147 0.13071130 -0.64739949
#> BreedJersey            1.2324714 0.13071130  0.97558625
#> AgeMature              0.1045993 0.08267006 -0.05787061
#>                         0.5quant 0.975quant mode
#> (Intercept)            4.0077182  4.2067007   NA
#> BreedCanadian          0.3784790  0.6353625   NA
#> BreedGuernsey          0.8899746  1.1468580   NA
#> BreedHolstein-Fresian -0.3905145 -0.1336307   NA
#> BreedJersey            1.2324717  1.4893550   NA
#> AgeMature              0.1045993  0.2670691   NA
#>                                kld
#> (Intercept)           2.524841e-09
#> BreedCanadian         2.524864e-09
#> BreedGuernsey         2.524860e-09
#> BreedHolstein-Fresian 2.524862e-09
#> BreedJersey           2.524862e-09
#> AgeMature             2.525095e-09
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
#> (Intercept)            4.0600181 0.09271686  3.8778080
#> BreedCanadian          0.3784786 0.13112173  0.1207927
#> BreedGuernsey          0.8899742 0.13112173  0.6322881
#> BreedHolstein-Fresian -0.3905148 0.13112173 -0.6482005
#> BreedJersey            1.2324713 0.13112173  0.9747851
#>                         0.5quant 0.975quant mode
#> (Intercept)            4.0600180  4.2422293   NA
#> BreedCanadian          0.3784788  0.6361633   NA
#> BreedGuernsey          0.8899745  1.1476588   NA
#> BreedHolstein-Fresian -0.3905146 -0.1328299   NA
#> BreedJersey            1.2324715  1.4901558   NA
#>                                kld
#> (Intercept)           2.471631e-09
#> BreedCanadian         2.471474e-09
#> BreedGuernsey         2.471475e-09
#> BreedHolstein-Fresian 2.471474e-09
#> BreedJersey           2.471469e-09
fit$waic$waic
#> [1] 115.9279
fit$dic$dic
#> [1] 115.0094
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
#> [1] 0.941422
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
#> (Intercept) 15.335 2.745      9.950   15.335     20.721
#> Father       0.406 0.029      0.349    0.406      0.463
#> Mother       0.322 0.031      0.260    0.322      0.383
#> GenderM      5.225 0.144      4.942    5.225      5.507
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
#> (Intercept) 15.345 2.747      9.957   15.345     20.734
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
#> Precision for the Gaussian observations      0.237   NA
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
#>     Pre = 2.56, Running = 0.238, Post = 0.0176, Total = 2.81 
#> Fixed effects:
#>                mean    sd 0.025quant 0.5quant 0.975quant
#> (Intercept) 358.167 2.734    352.778  358.167    363.555
#>             mode kld
#> (Intercept)   NA   0
#> 
#> Random effects:
#>   Name	  Model
#>     grower IID model
#>    box IID model
#> 
#> Model hyperparameters:
#>                                             mean       sd
#> Precision for the Gaussian observations 4.00e-03 1.00e-03
#> Precision for grower                    1.55e+04 1.56e+04
#> Precision for box                       2.11e+04 2.23e+04
#>                                         0.025quant 0.5quant
#> Precision for the Gaussian observations      0.002 4.00e-03
#> Precision for grower                      1765.697 1.09e+04
#> Precision for box                         2122.090 1.44e+04
#>                                         0.975quant mode
#> Precision for the Gaussian observations   6.00e-03   NA
#> Precision for grower                      5.70e+04   NA
#> Precision for box                         8.04e+04   NA
#> 
#> Deviance Information Criterion (DIC) ...............: 307.75
#> Deviance Information Criterion (DIC, saturated) ....: 26089.22
#> Effective number of parameters .....................: 1.84
#> 
#> Watanabe-Akaike information criterion (WAIC) ...: 307.66
#> Effective number of parameters .................: 1.67
#> 
#> Marginal log-Likelihood:  -166.64 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
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
#> 1  1  1.492374e-06 0.01130853 -0.02439995  8.736887e-07
#> 2  2 -1.044657e-05 0.01130854 -0.02442924 -6.115819e-06
#> 3  3  8.954184e-06 0.01130853 -0.02438167  5.242123e-06
#>   0.975quant mode          kld
#> 1 0.02440727   NA 0.0001310971
#> 2 0.02437801   NA 0.0001310977
#> 3 0.02442558   NA 0.0001310975
#> 
#> $box
#>   ID          mean         sd  0.025quant      0.5quant
#> 1  1  1.634571e-05 0.01047978 -0.02247619  9.734234e-06
#> 2  2 -9.371546e-06 0.01047977 -0.02254043 -5.580916e-06
#> 3  3 -4.140913e-06 0.01047977 -0.02252735 -2.465972e-06
#> 4  4 -2.833264e-06 0.01047977 -0.02252408 -1.687243e-06
#>   0.975quant mode          kld
#> 1 0.02255789   NA 0.0001322533
#> 2 0.02249359   NA 0.0001322526
#> 3 0.02250665   NA 0.0001322523
#> 4 0.02250992   NA 0.0001322523
```
Más que la inferencia sobre los efectos aleatorios, es importante la que se hace sobre las varianzas asociadas:

```r
fit$summary.hyperpar
#>                                                 mean
#> Precision for the Gaussian observations 3.831722e-03
#> Precision for grower                    1.547070e+04
#> Precision for box                       2.106179e+04
#>                                                   sd
#> Precision for the Gaussian observations 8.988398e-04
#> Precision for grower                    1.562618e+04
#> Precision for box                       2.234736e+04
#>                                           0.025quant
#> Precision for the Gaussian observations 2.271079e-03
#> Precision for grower                    1.765697e+03
#> Precision for box                       2.122090e+03
#>                                             0.5quant
#> Precision for the Gaussian observations 3.767080e-03
#> Precision for grower                    1.087754e+04
#> Precision for box                       1.440273e+04
#>                                           0.975quant mode
#> Precision for the Gaussian observations 5.782886e-03   NA
#> Precision for grower                    5.702135e+04   NA
#> Precision for box                       8.039870e+04   NA
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
#> Contribución media de grower a la varianza: 5.1e-05 por 100
cat(paste("Contribución media de box a la varianza:",round(mean(cB)*100,6),"por 100"))
#> Contribución media de box a la varianza: 4.5e-05 por 100
```


Ante estos resultados, y dados los valores del DIC (307.7519675) y del WAIC (307.6586291), se justifica la opción de prescindir de los efectos `grower` y `box` como efectos aleatorios y ajustar el modelo con un único efecto fijo global. 


```r
formula = wt ~ 1
fit = inla(formula, family="gaussian",data=broccoli,
           control.compute = list(dic=TRUE,waic=TRUE))  
fit$summary.fixed
#>                 mean       sd 0.025quant 0.5quant
#> (Intercept) 358.1667 2.793582   352.6576 358.1667
#>             0.975quant mode          kld
#> (Intercept)   363.6758   NA 1.090747e-08
fit$summary.hyperpar
#>                                                mean
#> Precision for the Gaussian observations 0.003790063
#>                                                   sd
#> Precision for the Gaussian observations 0.0008736111
#>                                          0.025quant
#> Precision for the Gaussian observations 0.002286081
#>                                            0.5quant
#> Precision for the Gaussian observations 0.003715627
#>                                          0.975quant mode
#> Precision for the Gaussian observations 0.005654613   NA
```
Vemos que la variación en los indicadores DIC (308.1460381) y WAIC (307.8947917) es despreciable para este nuevo modelo.

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

![](02-anova_files/figure-latex/unnamed-chunk-24-1.pdf)<!-- --> 

## Modelos mixtos

En ocasiones cuando ajustamos un modelo lineal tendremos algunos factores de clasificación que operan como efectos fijos y otros que operan como efectos aleatorios. Estaremos ante **modelos lineales mixtos**. Siendo estrictos, realmente el modelo con solo efectos aleatorios ya es un modelo mixto, puesto que incluye como efecto fijo una interceptación global.

En un modelo lineal mixto seguimos asumiendo una respuesta normal, *gaussian*, con media igual a un predictor lineal $\eta=X\beta + Z u$, donde $X$ es una matriz de diseño con los efectos fijos $x_1,x_2,...$, y $Z$ la correspondiente para los efectos aleatorios  $z_1, z_2,...$. Se asume además una varianza desconocida que puede ser distinta para distintos niveles de los predictores, y que en general se suele expresar a través de una matriz de covarianzas $\Sigma$, $(y|\eta,\Sigma) \sim N(\eta,\Sigma)$.

En INLA la fórmula de predicción de una respuesta $y$ a partir de un conjunto de efectos fijos x1,x2,..., y un conjunto de efectos aleatorios z1,z2,... se especifica como:


```r
formula = y ~ 1 + x1 + x2  + f(z1, model="") + f(z2,model="") 
```

De nuevo mencionar que la opción más habitual para los efectos aleatorios en un modelo lineal es `model="iid"`. 

Veamos cómo resolver las inferencias con la base de datos `penicillin` en la librería (`faraway`](https://cran.r-project.org/web/packages/faraway/faraway.pdf). Se recogen datos de la producción de penicilina (`yield`) para cuatro procesos de fabricación distintos (`treat`) y con varias mezclas de la materia prima, que son bastante variables. El objetivo es investigar las diferencias entre los procesos de fabricación, pero teniendo en cuenta la variabilidad extra que podrían introducir las mezclas. Es decir, estamos pensando en un modelo lineal para predecir `yield`, en el que `treat` interviene como efecto fijo y `blend` como efecto aleatorio.


```r
data(penicillin,package="faraway")
ggplot(penicillin,aes(x=treat,y=yield))+
  geom_point(aes(color=blend))
```

![](02-anova_files/figure-latex/unnamed-chunk-26-1.pdf)<!-- --> 
Modelizamos pues con INLA:

```r
penicillin$treat = relevel(penicillin$treat,"D")
prec.prior=list(prec=list(list(param = c(0.001, 0.001))))
  formula = yield ~ -1+ treat + f(blend,model="iid",hyper=prec.prior)
fit=inla(formula,family="gaussian",data=penicillin,
         control.predictor=list(compute=TRUE),
         control.compute=list(dic=TRUE,waic=TRUE))
fit$summary.fixed
#>            mean       sd 0.025quant 0.5quant 0.975quant
#> treatD 85.51299 2.383184   80.74987 85.53106   90.17464
#> treatA 83.52432 2.383016   78.76269 83.54196   88.18683
#> treatB 84.51865 2.383099   79.75628 84.53651   89.18073
#> treatC 88.49600 2.383442   83.73063 88.51470   93.15637
#>        mode          kld
#> treatD   NA 2.457502e-08
#> treatA   NA 2.447486e-08
#> treatB   NA 2.452448e-08
#> treatC   NA 2.472949e-08
fit$summary.random
#> $blend
#>       ID          mean         sd  0.025quant      0.5quant
#> 1 Blend1  1.594697e-04 0.01274530 -0.02672789  8.430042e-05
#> 2 Blend2 -6.276773e-05 0.01274325 -0.02736032 -3.324530e-05
#> 3 Blend3 -1.338306e-05 0.01274287 -0.02721817 -7.129237e-06
#> 4 Blend4  6.069374e-05 0.01274320 -0.02700668  3.204153e-05
#> 5 Blend5 -8.746062e-05 0.01274361 -0.02743175 -4.630681e-05
#>   0.975quant mode         kld
#> 1 0.02764182   NA 0.009996808
#> 2 0.02700140   NA 0.009959472
#> 3 0.02714209   NA 0.009952884
#> 4 0.02735487   NA 0.009958947
#> 5 0.02693139   NA 0.009965956
fit$summary.hyperpar
#>                                                mean
#> Precision for the Gaussian observations 0.062132647
#> Precision for blend                     0.004159303
#>                                                 sd
#> Precision for the Gaussian observations 0.02359525
#> Precision for blend                     0.02345196
#>                                           0.025quant
#> Precision for the Gaussian observations 0.0263688200
#> Precision for blend                     0.0000938719
#>                                             0.5quant
#> Precision for the Gaussian observations 0.0587500392
#> Precision for blend                     0.0008409634
#>                                         0.975quant mode
#> Precision for the Gaussian observations  0.1178322   NA
#> Precision for blend                      0.0272079   NA
fit$dic$dic; fit$waic$waic
#> [1] 128.3473
#> [1] 131.1274
```
Representamos en la figura \@ref(fig:penicillin2) las distribuciones posteriores de los efectos fijos y de los efectos aleatorios.

```r
fixed=NULL
for(i in 1:4){
  fixed=rbind(fixed,as.data.frame(fit$marginals.fixed[[i]]))
}
fixed
#>            x            y
#> 1   74.33659 2.030611e-05
#> 2   75.78705 1.357323e-04
#> 3   77.62177 1.164368e-03
#> 4   79.76863 9.779470e-03
#> 5   80.74987 2.220118e-02
#> 6   81.56739 4.037222e-02
#> 7   82.48144 7.096317e-02
#> 8   83.08280 9.613228e-02
#> 9   83.55309 1.170101e-01
#> 10  83.95194 1.341712e-01
#> 11  84.30699 1.479577e-01
#> 12  84.63378 1.585861e-01
#> 13  84.78980 1.627630e-01
#> 14  84.94221 1.661987e-01
#> 15  85.09186 1.689024e-01
#> 16  85.23944 1.708801e-01
#> 17  85.29804 1.714690e-01
#> 18  85.35646 1.719426e-01
#> 19  85.38562 1.721363e-01
#> 20  85.41475 1.723011e-01
#> 21  85.47294 1.725445e-01
#> 22  85.53106 1.726726e-01
#> 23  85.58916 1.726855e-01
#> 24  85.64728 1.725830e-01
#> 25  85.67635 1.724884e-01
#> 26  85.70544 1.723650e-01
#> 27  85.76370 1.720312e-01
#> 28  85.82210 1.715812e-01
#> 29  85.96892 1.699456e-01
#> 30  86.11751 1.675737e-01
#> 31  86.26853 1.644555e-01
#> 32  86.42280 1.605787e-01
#> 33  86.74485 1.504805e-01
#> 34  87.09320 1.370977e-01
#> 35  87.48262 1.201661e-01
#> 36  87.93938 9.927727e-02
#> 37  88.51993 7.375585e-02
#> 38  89.39597 4.228720e-02
#> 39  90.17464 2.336039e-02
#> 40  91.10584 1.031342e-02
#> 41  93.14554 1.219365e-03
#> 42  94.91260 1.383652e-04
#> 43  96.40410 1.813512e-05
#> 44  72.35137 2.027775e-05
#> 45  73.80241 1.357513e-04
#> 46  75.63644 1.164949e-03
#> 47  77.78206 9.785538e-03
#> 48  78.76269 2.221485e-02
#> 49  79.57972 4.039563e-02
#> 50  80.49327 7.099902e-02
#> 51  81.09435 9.617439e-02
#> 52  81.56445 1.170543e-01
#> 53  81.96316 1.342144e-01
#> 54  82.31811 1.479975e-01
#> 55  82.64481 1.586208e-01
#> 56  82.80080 1.627947e-01
#> 57  82.95319 1.662271e-01
#> 58  83.10282 1.689272e-01
#> 59  83.25037 1.709011e-01
#> 60  83.30897 1.714885e-01
#> 61  83.36738 1.719605e-01
#> 62  83.39654 1.721533e-01
#> 63  83.42567 1.723174e-01
#> 64  83.48386 1.725592e-01
#> 65  83.54196 1.726857e-01
#> 66  83.60006 1.726969e-01
#> 67  83.65818 1.725928e-01
#> 68  83.68725 1.724974e-01
#> 69  83.71634 1.723731e-01
#> 70  83.77459 1.720376e-01
#> 71  83.83299 1.715860e-01
#> 72  83.97981 1.699462e-01
#> 73  84.12839 1.675703e-01
#> 74  84.27942 1.644481e-01
#> 75  84.43370 1.605674e-01
#> 76  84.75578 1.504621e-01
#> 77  85.10418 1.370733e-01
#> 78  85.49369 1.201373e-01
#> 79  85.95057 9.924673e-02
#> 80  86.53132 7.372720e-02
#> 81  87.40775 4.226638e-02
#> 82  88.18683 2.334728e-02
#> 83  89.11857 1.030715e-02
#> 84  91.15952 1.218668e-03
#> 85  92.92745 1.383207e-04
#> 86  94.41865 1.815604e-05
#> 87  73.34399 2.029191e-05
#> 88  74.79474 1.357418e-04
#> 89  76.62911 1.164659e-03
#> 90  78.77534 9.782508e-03
#> 91  79.75628 2.220803e-02
#> 92  80.57356 4.038394e-02
#> 93  81.48736 7.098113e-02
#> 94  82.08857 9.615338e-02
#> 95  82.55878 1.170322e-01
#> 96  82.95755 1.341928e-01
#> 97  83.31255 1.479777e-01
#> 98  83.63930 1.586035e-01
#> 99  83.79530 1.627789e-01
#> 100 83.94770 1.662130e-01
#> 101 84.09734 1.689149e-01
#> 102 84.24491 1.708907e-01
#> 103 84.30351 1.714788e-01
#> 104 84.36192 1.719516e-01
#> 105 84.39108 1.721449e-01
#> 106 84.42021 1.723093e-01
#> 107 84.47840 1.725519e-01
#> 108 84.53651 1.726792e-01
#> 109 84.59461 1.726913e-01
#> 110 84.65273 1.725880e-01
#> 111 84.68180 1.724930e-01
#> 112 84.71089 1.723691e-01
#> 113 84.76914 1.720345e-01
#> 114 84.82754 1.715836e-01
#> 115 84.97437 1.699459e-01
#> 116 85.12295 1.675720e-01
#> 117 85.27398 1.644519e-01
#> 118 85.42825 1.605731e-01
#> 119 85.75032 1.504713e-01
#> 120 86.09869 1.370856e-01
#> 121 86.48815 1.201518e-01
#> 122 86.94498 9.926203e-02
#> 123 87.52562 7.374155e-02
#> 124 88.40186 4.227681e-02
#> 125 89.18073 2.335384e-02
#> 126 90.11220 1.031029e-02
#> 127 92.15253 1.219017e-03
#> 128 93.92002 1.383430e-04
#> 129 95.41137 1.814556e-05
#> 130 77.31444 2.034869e-05
#> 131 78.76401 1.357042e-04
#> 132 80.59977 1.163500e-03
#> 133 82.74847 9.770391e-03
#> 134 83.73063 2.218069e-02
#> 135 84.54889 4.033710e-02
#> 136 85.46368 7.090934e-02
#> 137 86.06546 9.606897e-02
#> 138 86.53605 1.169436e-01
#> 139 86.93510 1.341062e-01
#> 140 87.29031 1.478976e-01
#> 141 87.61722 1.585336e-01
#> 142 87.77329 1.627150e-01
#> 143 87.92574 1.661558e-01
#> 144 88.07543 1.688648e-01
#> 145 88.22304 1.708482e-01
#> 146 88.28165 1.714394e-01
#> 147 88.34008 1.719153e-01
#> 148 88.36924 1.721102e-01
#> 149 88.39838 1.722762e-01
#> 150 88.45658 1.725220e-01
#> 151 88.51470 1.726526e-01
#> 152 88.57281 1.726679e-01
#> 153 88.63093 1.725679e-01
#> 154 88.66001 1.724746e-01
#> 155 88.68910 1.723523e-01
#> 156 88.74736 1.720210e-01
#> 157 88.80576 1.715735e-01
#> 158 88.95259 1.699441e-01
#> 159 89.10117 1.675784e-01
#> 160 89.25219 1.644662e-01
#> 161 89.40645 1.605952e-01
#> 162 89.72845 1.505078e-01
#> 163 90.07673 1.371341e-01
#> 164 90.46603 1.202089e-01
#> 165 90.92261 9.932285e-02
#> 166 91.50285 7.379869e-02
#> 167 92.37832 4.231838e-02
#> 168 93.15637 2.338006e-02
#> 169 94.08675 1.032285e-02
#> 170 96.12458 1.220413e-03
#> 171 97.89032 1.384329e-04
#> 172 99.38228 1.810403e-05
```



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

**M0**. Solo estamos interesados en cuantificar el nivel general de habilidades lectoras, pero teniendo en cuenta posibles diferencias entre los niños. No nos interesan sin embargo, dichas diferencias. Pensamos pues en utilizar la variable `id` como efecto aleatorio y modelizar
$$\eta_{ij}=\eta_i=\theta + \alpha_i^{id}$$
con $\alpha_i^{id} \sim N(0,\sigma_{\alpha}^2)$ y a priori $\tau_{\alpha} \sim Ga(0.001,0.001)$.


```r
prec.prior=list(prec=list(param=c(0.001,0.001)))
formula0= read ~ f(id,model="iid",hyper = prec.prior) 
fit0=inla(formula0,family="gaussian",data=curran_dat)
fit=fit0
fit$summary.fixed
#>                 mean         sd 0.025quant 0.5quant
#> (Intercept) 4.114053 0.05109346   4.013248 4.114248
#>             0.975quant mode          kld
#> (Intercept)   4.213757   NA 1.440195e-09
fit$summary.hyperpar
#>                                              mean
#> Precision for the Gaussian observations 0.4181767
#> Precision for id                        3.7355177
#>                                                 sd
#> Precision for the Gaussian observations 0.01924267
#> Precision for id                        1.07498923
#>                                         0.025quant
#> Precision for the Gaussian observations  0.3809534
#> Precision for id                         2.1867716
#>                                          0.5quant
#> Precision for the Gaussian observations 0.4179702
#> Precision for id                        3.5464628
#>                                         0.975quant mode
#> Precision for the Gaussian observations  0.4567565   NA
#> Precision for id                         6.3647887   NA

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

![](02-anova_files/figure-latex/curranm0-1.pdf)<!-- --> 


La matriz de diseño $Z$ asociada a los efectos aleatorios se obtiene con la función `model.matrix`, y la podemos utilizar para ajustar el modelo con `inla`. Para ello habremos de definir un nuevo índice para todos los registros de la base de datos, y aplicar sobre ellos la matriz de efectos aleatorios:


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
#> Precision for the Gaussian observations 0.4182615
#> Precision for ID                        3.7319300
#>                                                 sd
#> Precision for the Gaussian observations 0.01928087
#> Precision for ID                        1.07622200
#>                                         0.025quant
#> Precision for the Gaussian observations  0.3809454
#> Precision for ID                         2.1865180
#>                                          0.5quant
#> Precision for the Gaussian observations 0.4180623
#> Precision for ID                        3.5410230
#>                                         0.975quant mode
#> Precision for the Gaussian observations  0.4568972   NA
#> Precision for ID                         6.3675389   NA
```

La diferencia entre ajustar los efectos aleatorios con "iid" o con "z" es que con la primera opción tendremos tantos efectos aleatorios como están definidos en el modelo. Sin embargo, la modelización con "z" producirá tantos efectos aleatorios como datos, con una única precisión/varianza asociada. Generalmente ambos modelos producen unas medias posterioris de los efectos aleatorios bastante parecidas. 


```r
# medias posteriori de los efectos aleatorias con model="iid"
random0=fit0$summary.random$id$mean
# medias posteriori de los efectos aleatorios con model="z"
random00=fit00$summary.random$ID$mean
plot(random00,random00,type="l")
points(random0,random0)
```

![](02-anova_files/figure-latex/curranm0b-1.pdf)<!-- --> 


**M1**. Dado que tenemos varias mediciones de un mismo sujeto, es razonable asumir que, por defecto, las habilidades cognitivas propias de un sujeto, $\alpha_i^{id}$, influyen en sus habilidades lectoras. 

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


```r
prec.prior=list(prec=list(param=c(0.001,0.001)))
formula1= read ~ f(id,model="iid",hyper = prec.prior) + f(occasion,model="iid",hyper = prec.prior)
fit1=inla(formula1,family="gaussian",data=curran_dat)
fit=fit1
fit$summary.fixed
#>                 mean       sd 0.025quant 0.5quant
#> (Intercept) 4.353434 0.889756    2.47961 4.353368
#>             0.975quant mode          kld
#> (Intercept)   6.227631   NA 0.0002595765
fit$summary.hyperpar
#>                                              mean        sd
#> Precision for the Gaussian observations 2.4595751 0.1147768
#> Precision for id                        1.2681469 0.1053370
#> Precision for occasion                  0.4871471 0.3919959
#>                                         0.025quant
#> Precision for the Gaussian observations 2.23998977
#> Precision for id                        1.07269487
#> Precision for occasion                  0.05849255
#>                                          0.5quant
#> Precision for the Gaussian observations 2.4573888
#> Precision for id                        1.2640450
#> Precision for occasion                  0.3835039
#>                                         0.975quant mode
#> Precision for the Gaussian observations   2.692123   NA
#> Precision for id                          1.487437   NA
#> Precision for occasion                    1.509043   NA

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

![](02-anova_files/figure-latex/curran2-1.pdf)<!-- --> 


- **M2**. Por otro lado, y en base a la Figura \@ref(fig:curran3) podríamos considerar el efecto del tiempo que transcurre desde el inicio del estudio (`t=time`), como una covariable numérica que afecta de modo lineal a las habilidades lectoras.



```r
ggplot(curran_dat, aes(x=occasion,y=read))+
  geom_line(aes(group=id),color="grey",size=0.4)
```

![(\#fig:curran3)Relación entre el tiempo y las habilidades lectoras para cada sujeto (líneas).](02-anova_files/figure-latex/curran3-1.pdf) 

Podríamos seguir planteando un efecto aleatorio del sujeto sobre sus resultados lectores, y un efecto fijo asociado al tiempo transcurrido hasta el instante $t_j=j$, esto es, una interceptación aleatoria y una pendiente fija.

$$\eta_{ij}=\theta + \alpha_i^{id} + \beta \cdot t_j $$
con 
\begin{eqnarray*}
\text{Nivel II} && \\
\theta &\sim & N(0,100) \\
\beta &\sim& N(0,100) \\
\alpha_i^{id} &\sim& N(0,\sigma_{\alpha}^2) \\
\text{Nivel III} && \\
\tau_{\alpha}=1/\sigma_{\alpha}^2 &\sim& Ga(0.001,0.001) 
\end{eqnarray*}


```r
prec.prior=list(prec=list(param=c(0.001,0.001)))
formula2= read ~ time + f(id,model="iid",hyper = prec.prior) 
fit2=inla(formula2,family="gaussian",data=curran_dat)
fit=fit2
fit$summary.fixed
#>                 mean         sd 0.025quant 0.5quant
#> (Intercept) 2.703744 0.05265826   2.600420 2.703750
#> time        1.101341 0.01758962   1.066827 1.101345
#>             0.975quant mode          kld
#> (Intercept)   2.807034   NA 1.185142e-11
#> time          1.135834   NA 5.561966e-12
fit$summary.hyperpar
#>                                             mean        sd
#> Precision for the Gaussian observations 2.174407 0.1013782
#> Precision for id                        1.285604 0.1091225
#>                                         0.025quant 0.5quant
#> Precision for the Gaussian observations   1.980419 2.172490
#> Precision for id                          1.083308 1.281299
#>                                         0.975quant mode
#> Precision for the Gaussian observations   2.379770   NA
#> Precision for id                          1.512947   NA

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

![](02-anova_files/figure-latex/curran4-1.pdf)<!-- --> 

- **M3**. Siendo más estrictos al mirar la Figura \@ref(fig:curran3), podríamos pensar que, además de tener diferentes orígenes en las habilidades lectoras para cada sujeto (interceptaciones distintas), dado el efecto lineal ascendente del tiempo sobre las habilidades lectoras de los sujetos, no todos los sujetos evolucionan al mismo ritmo, esto es, no todas las rectas son paralelas. Tendríamos entonces efectos aleatorios asociados a los sujetos, y vinculados a las interceptaciones y a las pendientes de las rectas.

$$\eta_{ij}=\theta + \alpha_i^{id} + \beta \cdot t_j + \beta_{ij}$$

con 
\begin{eqnarray*}
\text{Nivel II} && \\
\theta &\sim & N(0,100) \\
\beta &\sim& N(0,100) \\
\alpha_i^{id} &\sim& N(0,\sigma_{\alpha}^2) \\
\beta_{ij} &\sim & N(0,\sigma_{j})^2 \\
\text{Nivel III} && \\
\tau_{\alpha}=1/\sigma_{\alpha}^2 &\sim& Ga(0.001,0.001) \\
\tau_{j}=1/\sigma_{j}^2 &\sim&Ga(0.001,0.001)
\end{eqnarray*}

Además, dado que no disponemos de todas las mediciones temporales para todos los sujetos, tendremos que la variable tiempo (`time`) está anidada en la variable `id`. Este efecto de anidación que nos va a reportar la variabilidad que hay entre las pendientes distintas de los sujetos para cada instante de tiempo.


```r
prec.prior=list(prec=list(param=c(0.001,0.001)))
Z <- as(model.matrix( ~ 0 + id:time, data = curran_dat), "Matrix")
ID=1:nrow(curran_dat)
formula3= read ~ time + f(id,model="iid",hyper=prec.prior)+
                        f(ID,model="z",Z=Z,hyper = prec.prior) 
fit3=inla(formula3,family="gaussian",data=curran_dat)
fit=fit3
fit$summary.fixed
#>                 mean         sd 0.025quant 0.5quant
#> (Intercept) 2.695707 0.04697528   2.603547 2.695705
#> time        1.120191 0.02267335   1.075918 1.120114
#>             0.975quant mode          kld
#> (Intercept)   2.787883   NA 3.994798e-11
#> time          1.164903   NA 1.193638e-09
fit$summary.hyperpar
#>                                              mean        sd
#> Precision for the Gaussian observations  3.016271 0.1690637
#> Precision for id                         1.565002 0.1432518
#> Precision for ID                        11.393950 1.7120438
#>                                         0.025quant
#> Precision for the Gaussian observations   2.695117
#> Precision for id                          1.301779
#> Precision for ID                          8.461430
#>                                          0.5quant
#> Precision for the Gaussian observations  3.012269
#> Precision for id                         1.558507
#> Precision for ID                        11.241642
#>                                         0.975quant mode
#> Precision for the Gaussian observations   3.361027   NA
#> Precision for id                          1.865692   NA
#> Precision for ID                         15.189679   NA

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

![](02-anova_files/figure-latex/curran5-1.pdf)<!-- --> 
