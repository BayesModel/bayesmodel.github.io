# Modelos lineales generalizados {#glm}


Los modelos lineales generalizados (Generalized Linear Models or GLM), son una clase de modelos introducidos por Nelder y Wedderburn (1972) y McCullagh y Nelder (1989), con el objetivo de extender la regresión lineal al caso en el que la variable dependiente no se distribuya necesariamente según una normal, pero su distribución todavía pertenezca a la familia exponencial (Binomial, Poisson, Gamma, Gausiana inversa básicamente). Trabajamos a continuación con dos de los GLM más comunes en epidemiología y ciencias sociales: la regresión logística y la de Poisson, mostrando cómo usar R-INLA.

Un modelo lineal generalizado está basado en asumir, además de una distribución de los datos dentro de la familia exponencial, una relación lineal entre cierta transformación del valor esperado de la respuesta $E(y_i)$ y los predictores disponibles, sean covariables, efectos fijos, efectos aleatorios, o incluso alguna función de estos.


Si $y$ representa una respuesta observada, $x_1,x_2,..$ una serie de covariables o efectos fijos, y $z_1,z_2,...$ efectos aleatorios, el valor esperado de la respuesta lo denotamos como $\mu$, que 

$$E(y_i|x,z,\theta)=\mu_i$$ 
Pues bien, la relación entre esta media $\mu$ y un predictor lineal $\eta$ que construimos a partir de los predictores disponibles, viene dado por una función link $g$ tal que:

$$g(\mu_i)=\eta_i=\mu + \beta_1 x_{1i} +   \beta_2 x_{2i} +...+ z_{1i} + z_{2i}+...$$

Todos los parámetros involucrados en el predictor lineal $\eta$ son los efectos latentes del modelo (fijos o aleatorios). Estos, junto con el resto de parámetros definidos en este primer nivel de la modelización (nivel de datos), han de modelizarse a continuación, en un segundo nivel del modelo, con sus correspondientes distribuciones a priori. 


El predictor lineal $\eta$ está relacionado linealmente con los predictores según:
$$\eta_i=\mu + \beta_1 x_{1i} +   \beta_2 x_{2i} +...+ z_{1i} + z_{2i}+...$$
relación que se suele representar en forma matricial como 
$$\eta=X\beta + Z u$$

Los modelos lineales que hemos visto antes (regresión, anova, ancova, modelos mixtos) se engloban dentro del modelo lineal generalizado.


El argumento `control.predictor=list(compute=TRUE)` en la función `inla` permite obtener las distribuciones predictivas para el predictor lineal, que en estos casos es distinto a los valores ajustados, `fit$summary.fitted`. Además para obtener la distribución marginal de los valores ajustados y predichos necesitamos incorporar a la función `inla` el argumento `control.compute=list(return.marginals.predictor=TRUE)`, y ya con todo ello podemos:

-   `fit$summary.linear.predictor` resumir  la inferencia posterior sobre los
    predictores lineales (distintos a los fitted cuando hay una función *link*)
-   `fit$marginals.linear.predictor` graficar y describir las distribuciones posteriores marginales para los predictores lineales


##  Modelos jerárquicos bayesianos

A lo largo del curso ya hemos ido comentando en ocasiones, algo sobre la especificación de un modelo en varios niveles. Presentamos ya de lleno estos modelos lineales generalizados como modelos multi-nivel o modelos jerárquicos, denominados así porque se va especificando por niveles (o jerarquías) la información disponible sobre todo aquello que es desconocido, distribución de los datos y parámetros. 

Un modelo bayesiano se modeliza a través de un modelo jerárquico o multinivel en el que en el nivel I se define la distribución asumida sobre la variable respuesta y que determina la verosimilitud. Esta variable depende de unos parámetros que definen los efectos fijos y aleatorios, y para los que hay que proporcionar la información previa disponible a través de una distribución a priori en el segundo nivel del modelo jerárquico. La distribución a priori para los efectos fijos generalmente será común a todos ellos, mientras que la distribución a priori para los efectos aleatorios estará vinculada a otros hiperparámetros para los que también será preciso especificar una distribución a priori en un tercer nivel de la modelización, y así sucesivamente.


\begin{eqnarray*}
Nivel I &&\\
( y | X, Z, \theta) &\sim & f(y|x,z,\theta) \text{f en fam.exponencial}\\
&& E(y|x,z,\theta)=\mu; Var(y|x,z,\theta)=\Sigma \\
&& g(\mu)=\eta=X\beta + Z u \\
Nivel II &&\\
\beta &\sim & N(0,\sigma_{\beta}), \text{ con un valor dado  para } \sigma_{\beta} \\
u|\sigma_u^2 &\sim_{iid}&  N(0,{\sigma_u^2}) \\
\Sigma|s &\sim& F_{\Sigma|s} \\
Nivel III &&\\
\sigma_u^2 &\sim&  F_{\sigma} \\
s &\sim&  F_{s}
\end{eqnarray*}


## Regresión logística
La regresión logística es el modelo estándar para respuestas binarias (éxitos/fracasos). Tiene dos variaciones, en función de si la respuesta representa observaciones individuales (0/1) o conteos (de éxitos) en grupos de sujetos.

Si las observaciones son individualizadas, entonces

$$y_i|\pi_i \sim Ber(\pi_i), \ i=1,...,n$$
En el caso de que sean conteos en grupos, 

$$y_i|\pi_i\sim Bin(n_i,\pi_i), \ i=1,...,n$$
siendo $n_i$ el tamaño de cada uno de los $n$ grupos disponibles, y $\pi_i$ la probabilidad de éxito (output de interés).

La relación entre el predictor lineal $\eta$ construido con los predictores disponibles $x=(x_{1},...x_{M})$ y la probabilidad $\pi$ se especifica a través de la función *logit*:
$$logit(\pi)=log\left(\frac{\pi}{1-\pi}\right)=\eta=X\beta=\beta_0+\sum_{j=1}^M \beta_j x_{j}$$
de forma que
$$\pi=logit^{-1}(X\beta)=\frac{exp(X\beta)}{1-exp(X\beta)}$$
Una vez especificado el modelo, si no hay información previa disponible sobre los efectos (fijos) ${\beta_o,\beta_1,...\beta_M}$, se asumen distribuciones a priori independientes y  normales con media cero y varianza muy grande. 


### Interpretación de los coeficientes en la regresión logit {-}

Puesto que $X\beta=\beta_0+\sum_{j=1}^M \beta_j x_{j}$, la interceptación del predictor lineal $\beta_0$ se interpreta como los predictores toman el valor cero si son numéricos, o están en el nivel de referencia (para la estimación) si son categóricos, $\eta(X=0)=\beta_0=logit(\pi)$. En consecuencia, el logit inverso de $\beta_0$ se interpreta como la probabilidad de éxito $\pi_i$ cuando los predictores están en su nivel de referencia o son cero.
$$logit^{-1}(\beta_0)=Pr(y=1|X=0).$$ 

En cuanto a la interpretación de cualquier otro coeficiente de regresión en el predictor lineal, como $\beta_1$, echamos mano del concepto de *odds* y *odds ratio*. 

Los odds ratio, OR, comparan, a través de un cociente, las posibilidades a favor de un evento $E$ bajo condiciones A y de las posibilidades del mismo evento bajo condiciones B. Nos sirve para evaluar cuánto afecta a dicho evento el hecho de variar las condiciones de B a A.
$$OR(A,B)=\frac{Pr(E|A)/(1-Pr(E|A))}{Pr(E|B)/(1-Pr(E|B))}.$$

En el modelo logístico, nos interesa saber el efecto que tiene sobre la respuesta (realmente sobre las probabilidad de éxito) el incremento de una unidad en la variable $x$, y para ello consideramos los odds bajo $x$ y los odds bajo $x+1$, y en particular el logaritmo de los odds, log-odds:
$$log.odds(x+1)=log \left( \frac{P(y=1|x+1)}{P(y=0|x+1)} \right)=log \left(\frac{\pi}{1-\pi} | x+1\right)=\beta_0+\beta_1 (x+1)$$
$$log.odds(x)=log \left(\frac{P(y=1|x)}{P(y=0|x)} \right)=log \left(\frac{\pi}{1-\pi} | x\right)=\beta_0+\beta_1 x$$
de modo que 
$$log(OR(x+1,x))=log \left(\frac{ods (x+1)}{ods(x)} \right) = log.odds(x+1)-log.ods(x)= \beta_1$$
y tenemos entonces que la exponencial del coeficiente $\beta_1$ nos da el odds ratio asociado a dicha covariable.
$$exp(\beta_1)=\frac{odds(x+1)}{odds(x)}=OR(x+1,x).$$
Es decir, el coeficiente $\beta_1$ representa el cambio en los odds a favor de un éxito cuando se incrementa en una unidad el predictor $x$ al que acompaña en el predictor linea. Esta interpretación es muy común en Epidemiología.


### Intención de voto feb2022

Tenemos acceso a los datos completos obtenidos en la encuesta encargada por El País y la Cadena Ser a la empresa "40dB", en febrero de 2022, sobre la intención de voto nacional ([fuente](https://elpais.com/espana/2022-02-07/consulte-todos-los-datos-internos-de-la-encuesta-de-el-pais-cuestionarios-cruces-y-respuestas-individuales.html)).

Queremos predecir la probabilidad de votar al partido que gobierna mayoritariamente en la actualidad, PSOE, registrado en la variable `psoe`. Vamos a utilizar como predictores dos factores que nos dicen si el sujeto tiene simpatía por ese partido, `psoe_sim`, y si votó PSOE en las últimas elecciones `psoe_past`; también utilizaremos la comunidad autónoma `ccaa` como un efecto aleatorio, para contabilizar posible variación extra.


```r
library(readxl)
barometro_feb22 <- read_excel("~/Desktop/barometro/03_Datos_barómetro_febrero.xlsx")
names(barometro_feb22)
#>  [1] "id"                  "sexo"               
#>  [3] "edad"                "edad_r"             
#>  [5] "hab"                 "prov"               
#>  [7] "ccaa"                "edu"                
#>  [9] "cs"                  "p1"                 
#> [11] "p2"                  "p3"                 
#> [13] "p4_1"                "p4_2"               
#> [15] "p4_3"                "p4_4"               
#> [17] "p5"                  "p6"                 
#> [19] "p7"                  "hab_r"              
#> [21] "clase_social_r"      "situacion_laboral_r"
#> [23] "educacion_r"         "ponde"
datos=barometro_feb22 %>%
  select(id,p2,p3,p5,ccaa) %>%
  mutate(psoe=1*(p2=="PSOE (Partido Socialista Obrero Español)"),
         psoe_simp=1*(p3=="PSOE (Partido Socialista Obrero Español)"),
         psoe_past=1*(p5=="PSOE (Partido Socialista Obrero Español)")) 
#summary(datos)
```


Para ajustar el modelo utilizamos el argumento `family=binomial` y la opción `control.predictor = list(link = 1)` para establecer la función link apropiada para tener los valores ajustados en la escala correcta.


```r
prec.prior=list(prec=list(param=c(0.001,0.001)))
formula = psoe ~ psoe_simp + psoe_past+ f(ccaa,model="iid",hyper=prec.prior)
fit=inla(formula,family="binomial",data=datos,control.predictor = list(link = 1))
fit$summary.fixed
#>                  mean        sd 0.025quant  0.5quant
#> (Intercept) -3.590576 0.1897251  -3.998945 -3.578281
#> psoe_simp    3.085802 0.1865732   2.723932  3.084366
#> psoe_past    2.352353 0.1856829   1.989732  2.351807
#>             0.975quant mode          kld
#> (Intercept)  -3.251389   NA 5.812191e-07
#> psoe_simp     3.455821   NA 3.434611e-07
#> psoe_past     2.718074   NA 9.426592e-07
fit$summary.hyperpar
#>                        mean       sd 0.025quant 0.5quant
#> Precision for ccaa 70.76511 250.8882    1.97876 10.70226
#>                    0.975quant mode
#> Precision for ccaa   588.4032   NA
```

Las distribuciones posteriores de los exponenciales de los efectos aleatorios se muestran en la Figura \@ref(fig:eleccion2), identificadas en verde las de efectos positivos en la media posterior del predictor lineal (log-odds $>1$), y en rojo las de efectos negativos, e interpretables como más y menos favorables a votar por el PSOE.


```r
random = as.data.frame(fit$summary.random)
random$pro=1*exp(random$ccaa.mean)>1
ggplot(random,aes(x=exp(ccaa.mean),y=ccaa.ID)) +
  geom_point(aes(color=pro))+
  geom_errorbarh(aes(xmin=exp(ccaa.0.025quant),xmax=exp(ccaa.0.975quant),color=pro))+
  geom_vline(xintercept=1,linetype="dotted")+
  labs(x="Medias y HPD posteriores para el log-odds",y="Comunidad autónoma")+
  theme(legend.position="none")
```

![(\#fig:eleccion2)Medias y HPD posterioris para el log-odds de los efectos aleatorios](03-glm_files/figure-latex/eleccion2-1.pdf) 

La distribución posterior de la probabilidad de voto para el PSOE para la comunidad con más variabilidad en el efecto aleatorio (Castilla-La Mancha), viene representada en la Figura \@ref(fig:eleccion3) para las cuatro combinaciones posibles de valores para los predictores de simpatía y voto en el pasado.


```r
datos_pred = datos %>%
  mutate(f.post=round(fit$summary.fitted.values$mean,3),
        f.hpd.low=round(fit$summary.fitted.values$"0.025quant",3),
         f.hpd.up=round(fit$summary.fitted.values$"0.975quant",3)) %>%
  distinct(f.post,.keep_all = TRUE) %>%
  filter(ccaa=="Castilla - La Mancha") %>%
  mutate(simpast=str_c(psoe_simp,psoe_past))
ggplot(datos_pred,aes(x=f.post,y=simpast))+
  geom_point(aes(color=simpast))+
  geom_errorbarh(aes(xmin=f.hpd.low,xmax=f.hpd.up,color=simpast),height=0.2)+
  labs(x="Medias y HPD posteriores para la probabilidad de voto PSOE",title="Castilla - La Mancha")+
  scale_y_discrete(name="", 
                   labels=c("No simpatía/No votó PSOE","No simpatía/Votó PSOE",
                            "Simpatía/No votó PSOE","Simpatía/Votó PSOE"))+
  theme(legend.position="none")
```

![(\#fig:eleccion3,)Medias y HPD posterioris para la probabilidad de voto PSOE](03-glm_files/figure-latex/eleccion3,-1.pdf) 



### Mortalidad por infarto en Sheffield

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

![](03-glm_files/figure-latex/unnamed-chunk-6-1.pdf)<!-- --> 

El efecto $\beta_{12}$ representa el efecto en los log.odds de la mortalidad por infarto de estar en el nivel $NOx=2$ frente al de estar en el nivel $NOx=1$. Si queremos evaluar el odds-ratio, simplemente calculamos la distribución posterior de $exp(\beta_{12})$ con `inla.tmarginal`. Si sólo estamos interesados en su media, bastaría utilizar `inla.emarginal':


```r
odds.nox21 <- inla.tmarginal(function(x) exp(x), model.logistic$marginals.fixed$"factor(NOx)2")
e<-inla.emarginal(exp, model.logistic$marginals.fixed$"factor(NOx)2")
ggplot(data.frame(odds.nox21),aes(x=x,y=y))+geom_line()+
        labs(x=expression(OR(NOx[21])),y= expression(tilde(p)(paste(OR(NOx[21]),"|",y))))+
        geom_vline(xintercept=e,color="pink")+ geom_text(x=e,y=1,label=paste("mean=",round(e,3)),color="red")
```

![](03-glm_files/figure-latex/unnamed-chunk-7-1.pdf)<!-- --> 

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

La regresión de Poisson es útil cuando la variable respuesta representa conteos y estos toman valores discretos entre 0 y $+\infty$, sin una cota superior de referencia. El parámetro de interés es el número promedio de eventos $\lambda=E(y)$ y el link natural es el logaritmo, de modo que el predictor lineal está ligado con las covariables y factores según:
$$\eta=log(\lambda)=x \beta, \ \ \mbox{ y } \ \ \lambda_i=exp(X \beta)$$  

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

### Incidentes en barcos

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
#>     Pre = 2.13, Running = 0.129, Post = 0.00655, Total = 2.26 
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




