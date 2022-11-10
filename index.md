--- 
title: "Una primera aproximación a la Estadística Bayesiana"
author: "Asunción M. Mayoral y Javier Morales. IUI CIO-UMH"
date: "`Noviembre 2022`"
site: bookdown::bookdown_site
documentclass: book
bibliography: [book.bib, packages.bib]
# url: your book url like https://bookdown.org/yihui/bookdown
# cover-image: path to the social sharing image like images/cover.jpg
description: 
  Se trata de un curso de 10 horas ofrecido a alumnado de doctorado con formación BIO.
biblio-style: "apacite"
link-citations: TRUE
csl: chicago-fullnote-bibliography.csl
---



# Contexto

El curso "Una primera aproximación a la Estadística Bayesiana", con una duración de 10 horas, está ofertado a alumnado de doctorado con una formación científica en el ámbito BIO. 

El curso se desarrollará íntegramente online a través de videoconferencia síncrona, durante un total de 10 horas, seccionadas en 4 sesiones de dos horas y media cada una de ellas, los días 8, 10, 15 y 17 de noviembre de 2022, de 16 a 18:30h.

Se presentarán los contenidos a través de ejemplos prácticos programados en R, para que el estudiantado pueda ir generando resultados y comentando las interpretaciones derivadas del análisis.

La evaluación es continua basada en la interacción con el profesorado durante las sesiones de trabajo.

## Objetivos de aprendizaje 

- Conocer los conceptos básicos en el planteamiento bayesiano de la Estadística.
- Identificar la relevancia de la información previa y de expertos, y la proporcionada por los datos.
- Conocer los procedimientos básicos para conjugar la información disponible.
- Aplicar los conocimientos básicos en problemas sencillos.
- Descubrir las dificultades computacionales en la inferencia bayesiana.


## Contenidos propuestos

1. De probabilidad va la historia: la relevancia de las probabilidades condicionadas y el teorema de Bayes.
1. La jerga base: incertidumbre, a priori, a posteriori y verosímil.
1. Manos en la masa 1: ¿con qué probabilidad ocurrió?
1. Manos en la masa 2: ¿con qué abundancia ocurrió?
1. Curioseando para saber más: manuales y software.

## Contenidos definitivos
1. SESIÓN 1: Probabilidades, Bayes y las proporciones. 
1. SESIÓN 2: INLA y la regresión.
1. SESIÓN 3: INLA y el ANOVA.
1. SESIÓN 4: INLA y los GLM.

## INLA
INLA es una librería de R que aproxima la inferencia Bayesiana para
modelos gausianos latentes (LGM). Sus siglas provienen de Integrated
Nested Laplace Approximation (INLA), que es un método para aproximar las
inferencias bayesianas a través de la aproximación de Laplace.

Aunque la metodología INLA se ha desarrollado sobre modelos que se
pueden expresar como campos aleatorios markovianos gausianos (\*Gaussian
Markov random fields, GMRF), es viable para una gran familia de modelos
habituales en la práctica estadística.

Disponemos de referencias múltiples y documentación de esta librería en
la web [r-inla.org](https://www.r-inla.org), y en particular en el
manual de referencia de Gómez-Rubio (2021) titulado [Bayesian inference
with
INLA](https://becarioprecario.bitbucket.io/inla-gitbook/index.html),
también publicado por [Chapman & Hall-CRC
Press](https://www.routledge.com/Bayesian-inference-with-INLA/Gomez-Rubio/p/book/9781138039872).

### Instalación

Para instalar la librería INLA hemos de ejecutar, desde R, el comando


```r
install.packages("INLA", repos=c(getOption("repos"), 
                INLA="https://inla.r-inla-download.org/R/stable"), 
                dep=TRUE)
# y a continuación la cargamos con:
library(INLA)
```


Para instalar actualizaciones, basta con ejecutar


```r
options(repos = c(getOption("repos"), 
                  INLA="https://inla.r-inla-download.org/R/testing"))
update.packages("INLA", dep=TRUE)
```

Las descargas y documentación completa sobre INLA está disponible en
[R-INLA home](http://www.r-inla.org/home).

Ya desde R, para pedir ayuda sobre funciones en INLA, basta usar el
comando `inla.doc()`, especificando dentro y entrecomillada, la
función/objeto sobre el que se solicita ayuda. Por ejemplo,
`inla.doc("ar1")` o `inla.doc("loggamma")`.


También utilizaremos una librería accesoria de INLA, `brinla`,  desarrollada por [Faraway, Yue y Wan, 2022](https://rdrr.io/github/julianfaraway/brinla/).


```r
install.packages("remotes")
remotes::install_github("julianfaraway/brinla")
```

### Fundamentos

INLA está basado en la resolución de integrales vía la aproximación de
Laplace, que aproxima el integrando a través de una expansión de Taylor
de segundo grado que permite calcular la integral analíticamente.
\begin{eqnarray*}
I_n&=&\int_x exp[nf(x)]dx \\
&\approx& \int_x exp[n(f(x_0)+1/2 (x-x_0)^2 f''(x_0))] dx \\
&=& exp[nf(x_0)] \cdot \sqrt{\frac{2\pi}{-n f''(x_0)}}
\end{eqnarray*}

Evita así los largos tiempos de simulación de las cadenas de Markov
Monte Carlo. Cuando las distribuciones a integrar son Gausianas, Laplace
da órdenes buenos de aproximación. Y este es el principio que usa para
modelizar la mayoría de los modelos habituales, que se integran dentro
de la amplia clase de los modelos gausianos latentes, en los que se
aplica INLA.

INLA se ha aplicado en mapeo estadístico, modelos de cohorte
multidimensionales, modelos de asociación espacial, genética, análisis
medioambientales, salud y epidemiología, dinámicas de infecciones,
estudios agronómicos, meta-análisis, impacto del cambio climático y
muchos más ámbitos (Rue et al, 2017).



