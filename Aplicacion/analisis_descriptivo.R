#vamos a crear el ambiente de simulacion para el modelo 
#slash con n = 50. Los valores iniciales se proponen completando los datos faltantes
#por medio del algoritmo MDA propuesto por schafer(1997), luego
#utilizamos los pesos iguales a uno, suponiendo normalidad de los
#datos. Luego, corremos una iteración para extraer la matriz de 
#coeficiente y la matriz de varcov para utilizarlos como datos
#iniciales

setwd("ruta")
#-----------------------------------------------------------------------------
#Datos para aplicar el algoritmo MDA (Comuna Robledo, Medellín)
#-----------------------------------------------------------------------------
load(r"ruta")

#-----------------------------------------------------------------------------
#librería a usar
#------------------------------------------------------------------------------
library(readxl) 
library(beepr)
library(ggplot2)
library(tidyverse)
library(aplpack)
library(ggplotify)
#necesario para generar los pesos del modelo Slash
#install.packages("gamlss")
#install.packages("gamlss.dist")
#install.packages("gamlss.tr")
#library(gamlss)
#library(gamlss.tr)
#library(gamlss.dist)

#necesario para generar los datos iniciales
#library(norm)
#library(norm2)
#necesario para aplicar el proceso de imputación
library(mnorm)
#para medio vectorizar  a una matriz
library(fBasics)
library(mvtnorm)
#------------------------------------------------------------------------------


#-----------------------------------------------------------------------------
#matriz de datos (respuestas y covariables)
#---------------------------------------------------------------------------
datos <- base_de_datos #datos completos
datos <- as.matrix(datos)

#respuestas y covariables
respuestas <- datos[, 1:3] #matriz de variables respuesta
covariables <- cbind(rep(1, times = nrow(respuestas)), datos[,9], 
                    datos[, 6], datos[,8]) #matriz de covariables
colnames(covariables) <- c("X_0", "edad", "sexo", "tiempo_lechem")

datos_obs <- data.frame(respuestas, covariables) #conjunto de datos observados
#incluyendo las covariables

datos_1 <- datos_obs[-(which(is.na(datos_obs$per_bran) == T)), ]
datos_completos <- datos_1[-(which(is.na(datos_1$peso_act) == T)), ]
#-------------------------------------------------------------------

#-------------------------Analisis_Descriptivo-----------------------
#--------------------------------------------------------------------
# Crear el bagplot
per_bran <- respuestas[, 1]
peso <- respuestas[, 2]
talla <- respuestas[, 3]

#-----------------------------construcción_bagplots----------------------------

#Bagplot entre peso vs perímetro braquial
#par(mfrow = c(1, 1), pty = "s")
bagplot(per_bran, peso, show.whiskers = F, axes = F, xlab = "", ylab = "", show.looppoints = F, 
        show.bagpoints = F, show.baghull = TRUE, pch = 1, col.baghull = "gray", 
        col.loophull = "lightgray", cex = 1.5)
box()
at.Y1 <- seq(from = 6, to = 18, by = 2)
at.Y2 <- seq(from = 4, to = 16, by = 2)
axis(1, cex.axis = 1.2, at = at.Y1, line = -0.8, tick = F)
axis(2, cex.axis = 1.2, at = at.Y2, line = -0.8, tick = F)
title(xlab = "Perímetro braquial", line = 1.4, cex.lab = 1.4)
title(ylab = "Peso", line = 1.4, cex.lab = 1.4)
#bagplot_perbravsweight <- recordPlot()

#Bagplot entre talla vs perímetro braquial
#par(mfrow = c(1, 1), pty = "s")
bagplot(per_bran, talla, show.whiskers = F, axes = F, xlab = "", ylab = "", show.looppoints = F, 
        show.bagpoints = F, show.baghull = TRUE, pch = 1, col.baghull = "gray", 
        col.loophull = "lightgray", cex = 1.5)
box()
at.Y1 <- seq(from = 6, to = 18, by = 2)
at.Y2 <- seq(from = 50, to = 130, by = 10)
axis(1, cex.axis = 1.2, at = at.Y1, line = -0.8, tick = F)
axis(2, cex.axis = 1.2, at = at.Y2, line = -0.8, tick = F)
title(xlab = "Perímetro braquial", line = 1.4, cex.lab = 1.4)
title(ylab = "Talla", line = 1.4, cex.lab = 1.4)
#bagplot_perbravslength <- recordPlot()


#Bagplot entre talla vs peso
#par(mfrow = c(1, 1), pty = "s")
bagplot(peso, talla, show.whiskers = F, axes = F, xlab = "", ylab = "", show.looppoints = F, 
        show.bagpoints = F, show.baghull = TRUE, pch = 1, col.baghull = "gray", 
        col.loophull = "lightgray", cex = 1.5)
box()
at.Y1 <- seq(from = 4, to = 16, by = 2)
at.Y2 <- seq(from = 50, to = 130, by = 10)
axis(1, cex.axis = 1.2, at = at.Y1, line = -0.8, tick = F)
axis(2, cex.axis = 1.2, at = at.Y2, line = -0.8, tick = F)
title(xlab = "Peso", line = 1.4, cex.lab = 1.4)
title(ylab = "Talla", line = 1.4, cex.lab = 1.4)
#bagplot_lengthvsweight <- recordPlot()


#------------------------------------------------------------------------------
#-------------------------preparacion_boxplot----------------------------------
outcome.variable <- as.data.frame(respuestas)
matriz_modelo <- as.data.frame(covariables)
matriz_modelo$sexo <- as.factor(matriz_modelo$sexo)
matriz_modelo$age_cut <- as.factor(cut(matriz_modelo$edad,
                                       include.lowest = T,
                                       breaks = seq(0,4,1),
                                       labels=c("0 - 1","1 - 2","2 - 3",
                                                "3 - 4")))
eliminar <- which(is.na(outcome.variable$peso_act) == T)

outcome.variable <- outcome.variable[-eliminar, ]
matriz_modelo <- matriz_modelo[-eliminar, ]



#---------------------------creacion_boxplot----------------------------------

#box_plots peso vs edad y rango de edad
data.boxplot <- data.frame(outcome.variable$peso_act,matriz_modelo$sexo,
                           matriz_modelo$age_cut)
str(data.boxplot)
names(data.boxplot) <- c("weight", "gender", "age_cut")

#pdf("boxplot_logsn.pdf", height = 8.2, width = 16.4)
#par(mfrow = c(1, 1))
cols <- c("white", 'gray')
boxplot(weight~gender+age_cut,data=data.boxplot, at = c(1:2, 4:5, 7:8, 10:11),
        col = cols, xlab = "", ylab = "", xaxt = "n", yaxt="n")
legend("topleft", fill = cols, legend = c("Niña","Niño"), 
       horiz = F, bty="n",cex=2)
at_1.box <- levels(matriz_modelo$age_cut)
at_2.box <- seq(from = 4, to = 16, by = 2)
axis(side = 1, at = c(1.5,4.5,7.5,10.5),
     labels=at_1.box, cex.axis = 2.2, line = 0.1, tick = F)
axis(side = 2, at = at_2.box, cex.axis = 2.2, line = -0.8, tick = F)
title(xlab = "Edad", line = 2.6, cex.lab = 2.5)
title(ylab = "Peso", line = 2.4, cex.lab = 2.5)
#boxplot_peso_age <- recordPlot()
#------------------------------------------------------------------------------

#boxplot talla por edad y rango de edad
data.boxplot1 <- data.frame(outcome.variable$talla_act,matriz_modelo$sexo,
                           matriz_modelo$age_cut)
str(data.boxplot1)
names(data.boxplot1) <- c("length", "gender", "age_cut")

#pdf("boxplot_logsn.pdf", height = 8.2, width = 16.4)
#par(mfrow = c(1, 1))
cols <- c("white", 'gray')
boxplot(length~gender+age_cut,data=data.boxplot1, at = c(1:2, 4:5, 7:8, 10:11),
        col = cols, xlab = "", ylab = "", xaxt = "n", yaxt="n")
legend("topleft", fill = cols, legend = c("Niña","Niño"), 
       horiz = F, bty="n",cex=2)
at_1.box <- levels(matriz_modelo$age_cut)
at_2.box <- seq(from = 50, to = 130, by = 10)
axis(side = 1, at = c(1.5,4.5,7.5,10.5),
     labels=at_1.box, cex.axis = 2.2, line = 0.1, tick = F)
axis(side = 2, at = at_2.box, cex.axis = 2.2, line = -0.8, tick = F)
title(xlab = "Edad", line = 2.6, cex.lab = 2.5)
title(ylab = "Talla", line = 2.4, cex.lab = 2.5)
#boxplot_peso_age <- recordPlot()

#------------------------------------------------------------------------------
#boxplot perimetro braquial por edad y rango de edad
data.boxplot2 <- data.frame(outcome.variable$per_bran,matriz_modelo$sexo,
                            matriz_modelo$age_cut)
data.boxplot2 <- data.boxplot2[which(is.na(data.boxplot2$outcome.variable.per_bran) 
                                     == F), ] 
str(data.boxplot2)
names(data.boxplot2) <- c("armcircun", "gender", "age_cut")

#pdf("boxplot_logsn.pdf", height = 8.2, width = 16.4)
#par(mfrow = c(1, 1))
cols <- c("white", 'gray')
boxplot(armcircun~gender+age_cut,data=data.boxplot2, at = c(1:2, 4:5, 7:8, 10:11),
        col = cols, xlab = "", ylab = "", xaxt = "n", yaxt="n")
legend("bottomright", fill = cols, legend = c("Niña","Niño"), 
       horiz = F, bty="n",cex=2)
at_1.box <- levels(matriz_modelo$age_cut)
at_2.box <- seq(from = 6, to = 18, by = 2)
axis(side = 1, at = c(1.5,4.5,7.5,10.5),
     labels=at_1.box, cex.axis = 2.2, line = 0.1, tick = F)
axis(side = 2, at = at_2.box, cex.axis = 2.2, line = -0.8, tick = F)
title(xlab = "Edad", line = 2.6, cex.lab = 2.5)
title(ylab = "Perímetro braquial", line = 2.4, cex.lab = 2.5)

#-------------------------------------------------------------------------
#----------graficos variables respuesta vs tiempo leche materna------------
#--------------------------------------------------------------------------
# grafico perimetro braquial vs tiempo leche materna

ggplot(datos_completos, aes(x = datos_completos$tiempo_lechem, 
                            y = datos_completos$per_bran)) +
  geom_point(color = "black",  # Color de los puntos
             size = 2           # Tamaño de los puntos
             # Transparencia de los puntos
  ) +
  theme_minimal() +  # Utilizar un tema minimalista
  labs(
    #title = "",
    x = "Tiempo leche materna",
    y = "Perímetro braquial"
  ) +
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    axis.line = element_line(color = "black", size = 1),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", size = 0.5)#,
    #plot.background = element_rect(color = "black", size = 1, fill = NA)
  ) 

#---------------------------------------------------------------------------
# grafico peso vs tiempo leche materna

ggplot(datos_completos, aes(x = datos_completos$tiempo_lechem, 
                            y = datos_completos$peso_act)) +
  geom_point(
    color = "black",  # Color de los puntos
    size = 2           # Tamaño de los puntos
    # Transparencia de los puntos
  ) +
  theme_minimal() +  # Utilizar un tema minimalista
  labs(
    #title = "",
    x = "Tiempo leche materna",
    y = "Peso"
  ) +
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    axis.line = element_line(color = "black", size = 1),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", size = 0.5)#,
    #plot.background = element_rect(color = "black", size = 1, fill = NA)
  )

#---------------------------------------------------------------------------

# grafico peso vs tiempo leche materna


ggplot(datos_completos, aes(x = datos_completos$tiempo_lechem, 
                            y = datos_completos$talla_act)) +
  geom_point(
    color = "black",  # Color de los puntos
    size = 2           # Tamaño de los puntos
    # Transparencia de los puntos
  ) +
  theme_minimal() +  # Utilizar un tema minimalista
  labs(
    #title = "",
    x = "Tiempo leche materna",
    y = "Talla"
  ) +
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    axis.line = element_line(color = "black", size = 1),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", size = 0.5)#,
    #plot.background = element_rect(color = "black", size = 1, fill = NA)
  )

