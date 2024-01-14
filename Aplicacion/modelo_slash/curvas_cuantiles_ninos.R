setwd("ruta")
#-----------------------------------------------------------------------------
#Datos para aplicar el algoritmo MDA (Comuna Robledo, Medellín)
#-----------------------------------------------------------------------------
load(r"ruta")#ruta base de datos real
load("parametros_reales_slash.Rdata")

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

coef.slash <- estim[[1]]
matdisp.slash <- estim[[2]]
nu.slash <- estim[[3]]
log_datos <- log(respuestas)


imputacion_total <- function(log_datos, x_df, B_1, M_1, nu_1){
  n <- nrow(log_datos)
  p <- ncol(log_datos) 
  #aquí mantengo las filas y corro por columnas.
  for (i in 1:n){ #corredor de individuos
    #este me mantiene la secuencia de las filas en orden
    #para que no se me actualicen cada vez que j comienza a recorrer las filas
    for (j in 1:p){#corren las columnas
      if(is.na(log_datos[i,j]) == T){
        w_i <- rbeta(1, nu_1, 1)
        imp <- mvrnorm(1, t(B_1[, 1:p, drop = F])%*%x_df[i,], 
                       M_1/w_i)
        #aquí es suponiendo que conozco mv y betas
        log_datos[i,j] <- imp[j]
      }
    }
  }
  return(log_datos)
}

set.seed(123)
datos_imputados <- exp(imputacion_total(log_datos = log_datos, x_df = covariables,
                                        B_1 = coef.slash, M_1 = matdisp.slash, nu_1 = nu.slash))




data.children <- data.frame(datos_imputados, covariables)

colnames(data.children)[colnames(data.children) == "peso_act"] <- "weight"
colnames(data.children)[colnames(data.children) == "talla_act"] <- "length"
colnames(data.children)[colnames(data.children) == "sexo"] <- "gender"
colnames(data.children)[colnames(data.children) == "edad"] <- "age"
colnames(data.children)[colnames(data.children) == "per_bran"] <- "ac"
colnames(data.children)[colnames(data.children) == "tiempo_lechem"] <- "tl"



data.children$gender <- as.factor(data.children$gender)
levels(data.children$gender)[1] <- "female"
levels(data.children$gender)[2] <- "male"
attach(data.children)
data.female <- data.children[data.children$gender == "female", ]
data.male <- data.children[data.children$gender == "male", ]
weight.female <- data.female$weight
weight.male <- data.male$weight
length.female <- data.female$length
length.male <- data.male$length
age.female <- data.female$age
age.male <- data.male$age
armcirc.female <- data.female$ac
armcirc.male <- data.male$ac


# Quantile function of the univariate slash distribution
qslash <- function(p, mu, sigma2, nu, m){
  # p: probability
  # mu: location parameter
  # sigma2: dispersion parameter
  # nu: tail parameter
  # m: independent draws for Monte Carlo integration
  qtl.func <- function(w) {
    set.seed(123)
    u <- rbeta(m, nu, 1)
    integrand <- Vectorize(function(u) {
      stats::pnorm(q = w, mean = mu, sd = sqrt(sigma2/u))
    })
    return(mean(integrand(u)) - p)
  }
  return(uniroot(qtl.func, c(-1e+16,1e+16))$root)
}
#------------------------calculo_cuantiles_slash------------------------------
q.stdslash_0.005 <- qslash(p = 0.005, mu = 0, sigma2 = 1, 
                           nu = nu.slash, m = 3000)
q.stdslash_0.05 <- qslash(p = 0.05, mu = 0, sigma2 = 1, 
                          nu = nu.slash, m = 3000)
q.stdslash_0.25 <- qslash(p = 0.25, mu = 0, sigma2 = 1, 
                          nu = nu.slash, m = 3000)
q.stdslash_0.50 <- 0
q.stdslash_0.75 <- -q.stdslash_0.25
q.stdslash_0.95 <- -q.stdslash_0.05
q.stdslash_0.995 <- -q.stdslash_0.005
sqrt.sigma11 <- sqrt(matdisp.slash[1,1])
sqrt.sigma22 <- sqrt(matdisp.slash[2,2])
sqrt.sigma33 <- sqrt(matdisp.slash[3,3])
mean.tlm <- round(mean(data.children$tl), 2)
mean.age <- mean(data.children$age)


#---------------------------------quantiles male-------------------------------

# Weight (male)
quant.weight.male.t1wm <- list()
quant.weight.male.t2wm <- list()
t1wm <- list(seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05),
             seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05))
t2wm <- list(seq(3.05, 4, 0.05), seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),
             seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),seq(3.05, 4, 0.05))

quant.weight.male.t1wm[[1]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                t1wm[[1]] + coef.slash[3,2] + coef.slash[4,2] * mean.tlm + sqrt.sigma22* q.stdslash_0.005)
quant.weight.male.t2wm[[1]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                     t2wm[[1]] + coef.slash[3,2] +coef.slash[4,2] * mean.tlm + sqrt.sigma22*  q.stdslash_0.005)

quant.weight.male.t1wm[[2]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                       t1wm[[1]] + coef.slash[3,2]+ coef.slash[4,2] * mean.tlm + sqrt.sigma22* 
                                       q.stdslash_0.05)
quant.weight.male.t2wm[[2]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                       t2wm[[1]] + coef.slash[3,2] +coef.slash[4,2] * mean.tlm + sqrt.sigma22* 
                                       q.stdslash_0.05)

quant.weight.male.t1wm[[3]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                     t1wm[[1]] + coef.slash[3,2] +coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.25)
quant.weight.male.t2wm[[3]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                     t2wm[[1]] + coef.slash[3,2]+coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.25)

quant.weight.male.t1wm[[4]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                     t1wm[[1]] +coef.slash[3,2] +coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.50)
quant.weight.male.t2wm[[4]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                     t2wm[[1]] + coef.slash[3,2]+coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.50)
quant.weight.male.t1wm[[5]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                     t1wm[[1]] + coef.slash[3,2]+coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.75)
quant.weight.male.t2wm[[5]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                     t2wm[[1]] + coef.slash[3,2]+coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.75)

quant.weight.male.t1wm[[6]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                      t1wm[[1]] + coef.slash[3,2]+coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.95)
quant.weight.male.t2wm[[6]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                     t2wm[[1]] + coef.slash[3,2]+coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.95)

quant.weight.male.t1wm[[7]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                     t1wm[[1]] + coef.slash[3,2]+coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.995)
quant.weight.male.t2wm[[7]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                     t2wm[[1]] + coef.slash[3,2] + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.995)
#par(mfrow = c(1, 1), pty = "s")
plot(age.male, weight.male, xaxt = "n", yaxt = "n", xlab = "", 
     ylab = "", xlim = c(0, 4), ylim = c(0, 28), lty = 1, main = "", 
     cex = 1.4, col = "white")
at1wm.weight.male <- seq(from = 0, to = 4, by = 2)
at2wm.weight.male <- seq(from = 0, to = 28, by = 14)
axis(side = 1, at = at1wm.weight.male, cex.axis = 2.7, line = 0.3, 
     tick = F)
axis(side = 2, at = at2wm.weight.male, cex.axis = 2.7, line = -0.8, 
     tick = F)
title(xlab = "Edad", line = 3, cex.lab = 2.5)
title(ylab = "Peso (niño)", line = 2.4, cex.lab = 2.5)
for (i in 1:7) {
  lines(t1wm[[i]], quant.weight.male.t1wm[[i]], lwd = 3)
  lines(t2wm[[i]], quant.weight.male.t2wm[[i]], lwd = 3)
}
est.weight.male0.005 <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                mean.age + coef.slash[3,2] + coef.slash[4,2] * mean.tlm + sqrt.sigma22*q.stdslash_0.005)
text(2.9, est.weight.male0.005 + 1.9, expression(paste(0.5, th)), 
     srt = 14, cex = 1.1)
est.weight.male0.05 <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                               mean.age + coef.slash[3,2] +coef.slash[4,2] * mean.tlm + sqrt.sigma22* q.stdslash_0.05)
text(2.9, est.weight.male0.05 + 2.2, expression(paste(5, th)), srt = 14, 
     cex = 1.2)
est.weight.male0.25 <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                               mean.age + coef.slash[3,2] + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * q.stdslash_0.25)
text(2.9, est.weight.male0.25 + 2.5,expression(paste(25, th)), 
     srt = 14, cex = 1.2)
est.weight.male0.50 <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                               mean.age + coef.slash[3,2] +coef.slash[4,2] * mean.tlm + sqrt.sigma22* q.stdslash_0.50)
text(2.9, est.weight.male0.50 + 2.8, expression(paste(50, th)), 
     srt = 14, cex = 1.2)
est.weight.male0.75 <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                               mean.age + coef.slash[3,2] +coef.slash[4,2] * mean.tlm + sqrt.sigma22 * q.stdslash_0.75)
text(2.9, est.weight.male0.75 + 3.0, expression(paste(75, th)), 
     srt = 14, cex = 1.2)
est.weight.male0.95 <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                               mean.age + coef.slash[3,2] +coef.slash[4,2] * mean.tlm + sqrt.sigma22 * q.stdslash_0.95)
text(2.9, est.weight.male0.95 + 3.6, expression(paste(95, th)), 
     srt = 14, cex = 1.2)
est.weight.male0.995 <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                mean.age + coef.slash[3,2]+coef.slash[4,2] * mean.tlm + sqrt.sigma22* q.stdslash_0.995)
text(2.9, est.weight.male0.995 + 4.0, expression(paste(99.5, th)), 
     srt = 14, cex = 1.0)



# Length (male)
quant.length.male.t1lm <- list()
quant.length.male.t2lm <- list()
t1lm <- list(seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05),
             seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05))
t2lm <- list(seq(3.05, 4, 0.05), seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),
             seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),seq(3.05, 4, 0.05))

quant.length.male.t1lm[[1]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                     t1lm[[1]] + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.005)
quant.length.male.t2lm[[1]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                     t2lm[[1]] + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.005)

quant.length.male.t1lm[[2]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                     t1lm[[1]] + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.05)
quant.length.male.t2lm[[2]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                     t2lm[[1]] + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.05)

quant.length.male.t1lm[[3]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                     t1lm[[1]] + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.25)
quant.length.male.t2lm[[3]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                     t2lm[[1]] + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.25)

quant.length.male.t1lm[[4]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                     t1lm[[1]] + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.50)
quant.length.male.t2lm[[4]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                     t2lm[[1]] + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.50)

quant.length.male.t1lm[[5]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                     t1lm[[1]] + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.75)
quant.length.male.t2lm[[5]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                     t2lm[[1]] + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.75)

quant.length.male.t1lm[[6]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                     t1lm[[1]] + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.95)
quant.length.male.t2lm[[6]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                     t2lm[[1]] + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.95)

quant.length.male.t1lm[[7]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                     t1lm[[1]] + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.995)
quant.length.male.t2lm[[7]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                     t2lm[[1]] + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.995)

#par(mfrow = c(1, 1), pty = "s")
plot(age.male, length.male, xaxt = "n", yaxt = "n", xlab = "", 
     ylab = "", xlim = c(0, 4), ylim = c(44, 150), lty = 1, main = "", 
     cex = 1.4, col = "white")
at1lm.length.male <- seq(from = 0, to = 4, by = 2)
at2lm.length.male <- seq(from = 44, to = 150, by = 53)
axis(side = 1, at = at1lm.length.male, cex.axis = 2.7, line = 0.1, 
     tick = F)
axis(side = 2, at = at2lm.length.male, cex.axis = 2.7, line = -0.8, 
     tick = F)
title(xlab = "Edad", line = 3, cex.lab = 2.5)
title(ylab = "Talla (niño)", line = 2.4, cex.lab = 2.5)
for (i in 1:7) {
  lines(t1lm[[i]], quant.length.male.t1lm[[i]], lwd = 3)
  lines(t2lm[[i]], quant.length.male.t2lm[[i]], lwd = 3)
}
est.length.male0.005 <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                mean.age + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                q.stdslash_0.005)
text(2.9, est.length.male0.005 + 11.0, expression(paste(0.5, th)), 
     srt = 14, cex = 1.2)
est.length.male0.05 <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                               mean.age + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.05)
text(2.9, est.length.male0.05 + 12, expression(paste(5, th)), srt = 14, 
     cex = 1.2)
est.length.male0.25 <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                               mean.age + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.25)
text(2.9, est.length.male0.25 + 13.1, expression(paste(25, th)), 
     srt = 14, cex = 1.2)
est.length.male0.50 <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                               mean.age + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.50)
text(2.9, est.length.male0.50 + 13.1, expression(paste(50, th)), 
     srt = 14, cex = 1.2)
est.length.male0.75 <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                               mean.age + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.75)
text(2.9, est.length.male0.75 + 13.4, expression(paste(75, th)), srt = 14, 
     cex = 1.2)
est.length.male0.95 <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                               mean.age + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.95)
text(2.9, est.length.male0.95 + 15.2, expression(paste(95, th)), 
     srt = 14, cex = 1.2)
est.length.male0.995 <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                mean.age + coef.slash[3,3]+coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.995)
text(2.9, est.length.male0.995 + 16.2, expression(paste(99.5, th)), 
     srt = 14, cex = 0.9)



#arm circunference(male)
quant.armcirunf.male.t1acm <- list()
quant.armcirunf.male.t2acm <- list()
t1acm <- list(seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05),
              seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05))
t2acm <- list(seq(3.05, 4, 0.05), seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),
              seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),seq(3.05, 4, 0.05))

quant.armcirunf.male.t1acm[[1]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                         t1acm[[1]] + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                           q.stdslash_0.005)
quant.armcirunf.male.t2acm[[1]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                         t2acm[[1]] + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                           q.stdslash_0.005)

quant.armcirunf.male.t1acm[[2]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                         t1acm[[1]] + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                           q.stdslash_0.05)
quant.armcirunf.male.t2acm[[2]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                         t2acm[[1]] + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                           q.stdslash_0.05)

quant.armcirunf.male.t1acm[[3]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                         t1acm[[1]] + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                           q.stdslash_0.25)
quant.armcirunf.male.t2acm[[3]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                         t2acm[[1]] + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                           q.stdslash_0.25)

quant.armcirunf.male.t1acm[[4]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                         t1acm[[1]] + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                           q.stdslash_0.50)
quant.armcirunf.male.t2acm[[4]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                         t2acm[[1]] + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                           q.stdslash_0.50)

quant.armcirunf.male.t1acm[[5]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                         t1acm[[1]] + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                           q.stdslash_0.75)
quant.armcirunf.male.t2acm[[5]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                         t2acm[[1]] + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                           q.stdslash_0.75)

quant.armcirunf.male.t1acm[[6]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                         t1acm[[1]] + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                           q.stdslash_0.95)
quant.armcirunf.male.t2acm[[6]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                         t2acm[[1]] + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                           q.stdslash_0.95)

quant.armcirunf.male.t1acm[[7]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                         t1acm[[1]] + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                           q.stdslash_0.995)
quant.armcirunf.male.t2acm[[7]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                         t2acm[[1]] + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                           q.stdslash_0.995)

#par(mfrow = c(1, 1), pty = "s")
plot(age.male, armcirc.male, xaxt = "n", yaxt = "n", xlab = "", 
     ylab = "", xlim = c(0, 4), ylim = c(2, 26), lty = 1, main = "", 
     cex = 1.4, col = "white")
at1acm.ac.male <- seq(from = 0, to = 4, by = 2)
at2acm.ac.male <- seq(from = 2, to = 26, by = 12)
axis(side = 1, at = at1acm.ac.male, cex.axis = 2.7, line = 0.1, 
     tick = F)
axis(side = 2, at = at2acm.ac.male, cex.axis = 2.7, line = -0.8, 
     tick = F)
title(xlab = "Edad", line = 3, cex.lab = 2.5)
title(ylab = "Perímetro braquial (niño)", line = 2.4, cex.lab = 2.5)
for (i in 1:7) {
  lines(t1acm[[i]], quant.armcirunf.male.t1acm[[i]], lwd = 3)
  lines(t2acm[[i]], quant.armcirunf.male.t2acm[[i]], lwd = 3)
}
est.armcirunf.male0.005 <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                   mean.age + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                   q.stdslash_0.005)
text(2.9, est.armcirunf.male0.005 + 0.6, expression(paste(0.5, th)), 
     srt = 14, cex = 1.0)
est.armcirunf.male0.05 <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                  mean.age + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.05)
text(2.9, est.armcirunf.male0.05 + 0.6, expression(paste(5, th)), srt = 14, 
     cex = 1.2)
est.armcirunf.male0.25 <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                  mean.age + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.25)
text(2.9, est.armcirunf.male0.25 + 0.6, expression(paste(25, th)), 
     srt = 14, cex = 1.2)
est.armcirunf.male0.50 <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                  mean.age + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.50)
text(2.9, est.armcirunf.male0.50 + 0.6, expression(paste(50, th)), 
     srt = 14, cex = 1.2)
est.armcirunf.male0.75 <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                  mean.age + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.75)
text(2.9, est.armcirunf.male0.75 + 0.6, expression(paste(75, th)), srt = 14, 
     cex = 1.2)
est.armcirunf.male0.95 <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                  mean.age + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.95)
text(2.9, est.armcirunf.male0.95 + 0.6, expression(paste(95, th)), 
     srt = 14, cex = 1.2)
est.armcirunf.male0.995 <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                   mean.age + coef.slash[3,1]+coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.995)
text(2.9, est.armcirunf.male0.995 + 0.6, expression(paste(99.5, th)), 
     srt = 14, cex = 0.9)




