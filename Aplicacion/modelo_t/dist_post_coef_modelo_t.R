
#distribución posterior de los parámetros del modelo slash
setwd("ruta")
load("parametros_reales_modelo_t.rdata")
load("aplicacion_modelo_t.Rdata")
library(fBasics)
library(ggplot2)

#----------------------preparacion------------------------------------------
#extraemos la cadena de cada uno de los coeficientes
simulaciones <- aplicacion_t
#obtenidas
Nsim = length(simulaciones)
matriz_coeficientes <- list()
matriz_varianzas <- list()
nu <- c()
for (i in 1:Nsim) {
  matriz_coeficientes[[i]] <- simulaciones[[i]][[1]] #matriz de betas 
  matriz_varianzas[[i]] <- simulaciones[[i]][[2]] #matriz de varianzas y covarianzas
  nu[i] <- simulaciones[[i]][[3]]
}

A <- matrix(0, ncol = nrow(as.matrix(matriz_coeficientes[[1]]))*
              ncol(as.matrix(matriz_coeficientes[[1]])), 
            nrow = length(matriz_coeficientes))
colnames(A) <- c("B_01", "B_11", "B_21", "B_31",
                 "B_02", "B_12", "B_22", "B_32",
                 "B_03", "B_13", "B_23", "B_33")

for (i in 1:length(matriz_coeficientes)) {
  coef <-  vec(matriz_coeficientes[[i]])
  A[i, ] <- coef
}


A <- as.data.frame(A)
nu <- as.data.frame(nu)
#----------------grafico de la distirbucion posterior de los betas------------
#---------------------------------------------------------------------------
#b_01
ggplot(A, aes(x = B_01)) +
  geom_density(fill = NA, color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1)) +
  labs(x = expression(paste(beta[0][1])), y = "Densidad")

#B_11

ggplot(A, aes(x = B_11)) +
  geom_density(fill = NA, color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1)) +
  labs(x = expression(paste(beta[1][1])), y = "Densidad")
#B_21

ggplot(A, aes(x = B_21)) +
  geom_density(fill = NA, color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1)) +
  labs(x = expression(paste(beta[2][1])), y = "Densidad")
#B_31

ggplot(A, aes(x = B_31)) +
  geom_density(fill = NA, color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1)) +
  labs(x = expression(paste(beta[3][1])), y = "Densidad")
#B_02

ggplot(A, aes(x = B_02)) +
  geom_density(fill = NA, color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1)) +
  labs(x = expression(paste(beta[0][2])), y = "Densidad")
#B_12

ggplot(A, aes(x = B_12)) +
  geom_density(fill = NA, color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1)) +
  labs(x = expression(paste(beta[1][2])), y = "Densidad")
#B_22

ggplot(A, aes(x = B_22)) +
  geom_density(fill = NA, color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1)) +
  labs(x = expression(paste(beta[2][2])), y = "Densidad")
#B_32

ggplot(A, aes(x = B_32)) +
  geom_density(fill = NA, color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1)) +
  labs(x = expression(paste(beta[3][2])), y = "Densidad")
#B_03

ggplot(A, aes(x = B_03)) +
  geom_density(fill = NA, color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1)) +
  labs(x = expression(paste(beta[0][3])), y = "Densidad")
#B_13

ggplot(A, aes(x = B_13)) +
  geom_density(fill = NA, color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1)) +
  labs(x = expression(paste(beta[1][3])), y = "Densidad")
#B_23

ggplot(A, aes(x = B_23)) +
  geom_density(fill = NA, color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1)) +
  labs(x = expression(paste(beta[2][3])), y = "Densidad")
#B_33
ggplot(A, aes(x = B_33)) +
  geom_density(fill = NA, color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1)) +
  labs(x = expression(paste(beta[3][3])), y = "Densidad")

#nu

ggplot(nu, aes(x = nu)) +
  geom_density(fill = NA, color = "black", alpha = 0.7) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        axis.line = element_line(color = "black", size = 1)) +
  labs(x = expression(paste(nu)), y = "Densidad")
