setwd("ruta")
#-----------------------------------------------------------------------------
#Datos para aplicar el algoritmo MDA (Comuna Robledo, Medellín)
#-----------------------------------------------------------------------------
load(r"ruta")
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

#---------------------------------quantiles female-----------------------------

# Weight (female)
quant.weight.female.t1wf <- list()
quant.weight.female.t2wf <- list()
t1wf <- list(seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05),
             seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05))
t2wf <- list(seq(3.05, 4, 0.05), seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),
             seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),seq(3.05, 4, 0.05))

quant.weight.female.t1wf[[1]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
      t1wf[[1]] + coef.slash[4,2] * mean.tlm + sqrt.sigma22* 
        q.stdslash_0.005)
quant.weight.female.t2wf[[1]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                       t2wf[[1]] + coef.slash[4,2] * mean.tlm + sqrt.sigma22* 
                                       q.stdslash_0.005)

quant.weight.female.t1wf[[2]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                       t1wf[[1]] + coef.slash[4,2] * mean.tlm + sqrt.sigma22* 
                                       q.stdslash_0.05)
quant.weight.female.t2wf[[2]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                       t2wf[[1]] + coef.slash[4,2] * mean.tlm + sqrt.sigma22* 
                                       q.stdslash_0.05)

quant.weight.female.t1wf[[3]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                       t1wf[[1]] + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.25)
quant.weight.female.t2wf[[3]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                       t2wf[[1]] + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.25)

quant.weight.female.t1wf[[4]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                       t1wf[[1]] + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.50)
quant.weight.female.t2wf[[4]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                       t2wf[[1]] + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.50)
quant.weight.female.t1wf[[5]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                       t1wf[[1]] + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.75)
quant.weight.female.t2wf[[5]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                       t2wf[[1]] + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.75)

quant.weight.female.t1wf[[6]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                       t1wf[[1]] + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.95)
quant.weight.female.t2wf[[6]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                       t2wf[[1]] + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.95)

quant.weight.female.t1wf[[7]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                       t1wf[[1]] + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.995)
quant.weight.female.t2wf[[7]] <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                       t2wf[[1]] + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * 
                                       q.stdslash_0.995)
#par(mfrow = c(1, 1), pty = "s")
plot(age.female, weight.female, xaxt = "n", yaxt = "n", xlab = "", 
     ylab = "", xlim = c(0, 4), ylim = c(0, 28), lty = 1, main = "", 
     cex = 1.4, col = "white")
at1wf.weight.female <- seq(from = 0, to = 4, by = 2)
at2wf.weight.female <- seq(from = 0, to = 28, by = 14)
axis(side = 1, at = at1wf.weight.female, cex.axis = 2.7, line = 0.3, 
     tick = F)
axis(side = 2, at = at2wf.weight.female, cex.axis = 2.7, line = -0.8, 
     tick = F)
title(xlab = "Edad", line = 3, cex.lab = 2.5)
title(ylab = "Peso (niña)", line = 2.2, cex.lab = 2.5)
for (i in 1:7) {
  lines(t1wf[[i]], quant.weight.female.t1wf[[i]], lwd = 3)
  lines(t2wf[[i]], quant.weight.female.t2wf[[i]], lwd = 3)
}
est.weight.female0.005 <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                mean.age + coef.slash[4,2] * mean.tlm + sqrt.sigma22*q.stdslash_0.005)
text(2.9, est.weight.female0.005 + 1.9, expression(paste(0.5, th)), 
     srt = 14, cex = 1.1)
est.weight.female0.05 <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                               mean.age + coef.slash[4,2] * mean.tlm + sqrt.sigma22* q.stdslash_0.05)
text(2.9, est.weight.female0.05 + 2.2, expression(paste(5, th)), srt = 14, 
     cex = 1.1)
est.weight.female0.25 <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                               mean.age + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * q.stdslash_0.25)
text(2.9, est.weight.female0.25 + 2.5,expression(paste(25, th)), 
     srt = 14, cex = 1.1)
est.weight.female0.50 <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                 mean.age + coef.slash[4,2] * mean.tlm + sqrt.sigma22* q.stdslash_0.50)
text(2.9, est.weight.female0.50 + 3.0, expression(paste(50, th)), 
     srt = 14, cex = 1.1)
est.weight.female0.75 <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                   mean.age + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * q.stdslash_0.75)
text(2.9, est.weight.female0.75 + 3.0, expression(paste(75, th)), 
     srt = 14, cex = 1.1)
est.weight.female0.95 <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                               mean.age + coef.slash[4,2] * mean.tlm + sqrt.sigma22 * q.stdslash_0.95)
text(2.9, est.weight.female0.95 + 3.6, expression(paste(95, th)), 
     srt = 14, cex = 1.1)
est.weight.female0.995 <- exp(coef.slash[1,2] + coef.slash[2,2]* 
                                mean.age + coef.slash[4,2] * mean.tlm + sqrt.sigma22* q.stdslash_0.995)
text(2.9, est.weight.female0.995 + 4.0, expression(paste(99.5, th)), 
     srt = 14, cex = 1.1)



# Length (female)
quant.length.female.t1lf <- list()
quant.length.female.t2lf <- list()
t1lf <- list(seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05),
             seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05))
t2lf <- list(seq(3.05, 4, 0.05), seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),
             seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),seq(3.05, 4, 0.05))

quant.length.female.t1lf[[1]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                       t1lf[[1]] + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.005)
quant.length.female.t2lf[[1]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                       t2lf[[1]] + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.005)

quant.length.female.t1lf[[2]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                       t1lf[[1]] + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.05)
quant.length.female.t2lf[[2]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                       t2lf[[1]] + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.05)

quant.length.female.t1lf[[3]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                       t1lf[[1]] + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.25)
quant.length.female.t2lf[[3]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                       t2lf[[1]] + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.25)

quant.length.female.t1lf[[4]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                       t1lf[[1]] + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.50)
quant.length.female.t2lf[[4]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                       t2lf[[1]] + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.50)

quant.length.female.t1lf[[5]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                       t1lf[[1]] + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.75)
quant.length.female.t2lf[[5]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                       t2lf[[1]] + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.75)

quant.length.female.t1lf[[6]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                       t1lf[[1]] + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.95)
quant.length.female.t2lf[[6]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                       t2lf[[1]] + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.95)

quant.length.female.t1lf[[7]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                       t1lf[[1]] + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.995)
quant.length.female.t2lf[[7]] <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                       t2lf[[1]] + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                       q.stdslash_0.995)

#par(mfrow = c(1, 1), pty = "s")
plot(age.female, length.female, xaxt = "n", yaxt = "n", xlab = "", 
     ylab = "", xlim = c(0, 4), ylim = c(44, 150), lty = 1, main = "", 
     cex = 1.4, col = "white")
at1lf.length.female <- seq(from = 0, to = 4, by = 2)
at2lf.length.female <- seq(from = 44, to = 150, by = 53)
axis(side = 1, at = at1lf.length.female, cex.axis = 2.7, line = 0.1, 
     tick = F)
axis(side = 2, at = at2lf.length.female, cex.axis = 2.7, line = -0.8, 
     tick = F)
title(xlab = "Edad", line = 3, cex.lab = 2.5)
title(ylab = "Talla (niña)", line = 2.4, cex.lab = 2.5)
for (i in 1:7) {
  lines(t1lf[[i]], quant.length.female.t1lf[[i]], lwd = 3)
  lines(t2lf[[i]], quant.length.female.t2lf[[i]], lwd = 3)
}
est.length.female0.005 <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                mean.age + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * 
                                q.stdslash_0.005)
text(2.9, est.length.female0.005 + 11.2, expression(paste(0.5, th)), 
     srt = 14, cex = 1.2)
est.length.female0.05 <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                               mean.age + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.05)
text(2.9, est.length.female0.05 + 12, expression(paste(5, th)), srt = 14, 
     cex = 1.3)
est.length.female0.25 <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                               mean.age + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.25)
text(2.9, est.length.female0.25 + 13.1, expression(paste(25, th)), 
     srt = 14, cex = 1.3)
est.length.female0.50 <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                               mean.age + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.50)
text(2.9, est.length.female0.50 + 13.1, expression(paste(50, th)), 
     srt = 14, cex = 1.3)
est.length.female0.75 <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                               mean.age + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.75)
text(2.9, est.length.female0.75 + 13.4, expression(paste(75, th)), srt = 14, 
     cex = 1.3)
est.length.female0.95 <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                               mean.age + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.95)
text(2.9, est.length.female0.95 + 15.2, expression(paste(95, th)), 
     srt = 14, cex = 1.3)
est.length.female0.995 <- exp(coef.slash[1,3] + coef.slash[2,3] * 
                                mean.age + coef.slash[4,3] * mean.tlm + sqrt.sigma33 * q.stdslash_0.995)
text(2.9, est.length.female0.995 + 16.2, expression(paste(99.5, th)), 
     srt = 14, cex = 1.3)



#arm circunference(female)
quant.armcirunf.female.t1acf <- list()
quant.armcirunf.female.t2acf <- list()
t1acf <- list(seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05),
              seq(0, 2.75, 0.05), seq(0, 2.75, 0.05), seq(0, 2.75, 0.05))
t2acf <- list(seq(3.05, 4, 0.05), seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),
              seq(3.05, 4, 0.05),seq(3.05, 4, 0.05),seq(3.05, 4, 0.05))

quant.armcirunf.female.t1acf[[1]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                           t1acf[[1]] + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                       q.stdslash_0.005)
quant.armcirunf.female.t2acf[[1]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                           t2acf[[1]] + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                       q.stdslash_0.005)

quant.armcirunf.female.t1acf[[2]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                           t1acf[[1]] + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                       q.stdslash_0.05)
quant.armcirunf.female.t2acf[[2]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                           t2acf[[1]] + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                       q.stdslash_0.05)

quant.armcirunf.female.t1acf[[3]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                           t1acf[[1]] + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                       q.stdslash_0.25)
quant.armcirunf.female.t2acf[[3]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                           t2acf[[1]] + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                       q.stdslash_0.25)

quant.armcirunf.female.t1acf[[4]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                           t1acf[[1]] + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                       q.stdslash_0.50)
quant.armcirunf.female.t2acf[[4]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                           t2acf[[1]] + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                       q.stdslash_0.50)

quant.armcirunf.female.t1acf[[5]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                           t1acf[[1]] + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                       q.stdslash_0.75)
quant.armcirunf.female.t2acf[[5]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                           t2acf[[1]] + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                       q.stdslash_0.75)

quant.armcirunf.female.t1acf[[6]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                           t1acf[[1]] + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                       q.stdslash_0.95)
quant.armcirunf.female.t2acf[[6]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                           t2acf[[1]] + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                       q.stdslash_0.95)

quant.armcirunf.female.t1acf[[7]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                           t1acf[[1]] + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                       q.stdslash_0.995)
quant.armcirunf.female.t2acf[[7]] <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                           t2acf[[1]] + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                       q.stdslash_0.995)

#par(mfrow = c(1, 1), pty = "s")
plot(age.female, armcirc.female, xaxt = "n", yaxt = "n", xlab = "", 
     ylab = "", xlim = c(0, 4), ylim = c(2, 26), lty = 1, main = "", 
     cex = 1.4, col = "white")
at1acf.ac.female <- seq(from = 0, to = 4, by = 2)
at2acf.ac.female <- seq(from = 2, to = 26, by = 12)
axis(side = 1, at = at1acf.ac.female, cex.axis = 2.7, line = 0.1, 
     tick = F)
axis(side = 2, at = at2acf.ac.female, cex.axis = 2.7, line = -0.8, 
     tick = F)
title(xlab = "Edad", line = 3, cex.lab = 2.5)
title(ylab = "Perímetro braquial (niña)", line = 2.4, cex.lab = 2.5)
for (i in 1:7) {
  lines(t1acf[[i]], quant.armcirunf.female.t1acf[[i]], lwd = 3)
  lines(t2acf[[i]], quant.armcirunf.female.t2acf[[i]], lwd = 3)
}
est.armcirunf.female0.005 <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                   mean.age + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * 
                                q.stdslash_0.005)
text(2.9, est.armcirunf.female0.005 + 0.6, expression(paste(0.5, th)), 
     srt = 14, cex = 1.2)
est.armcirunf.female0.05 <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                  mean.age + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.05)
text(2.9, est.armcirunf.female0.05 + 0.6, expression(paste(5, th)), srt = 14, 
     cex = 1.3)
est.armcirunf.female0.25 <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                  mean.age + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.25)
text(2.9, est.armcirunf.female0.25 + 0.6, expression(paste(25, th)), 
     srt = 14, cex = 1.3)
est.armcirunf.female0.50 <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                  mean.age + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.50)
text(2.9, est.armcirunf.female0.50 + 0.6, expression(paste(50, th)), 
     srt = 14, cex = 1.3)
est.armcirunf.female0.75 <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                  mean.age + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.75)
text(2.9, est.armcirunf.female0.75 + 0.6, expression(paste(75, th)), srt = 14, 
     cex = 1.3)
est.armcirunf.female0.95 <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                  mean.age + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.95)
text(2.9, est.armcirunf.female0.95 + 0.6, expression(paste(95, th)), 
     srt = 14, cex = 1.3)
est.armcirunf.female0.995 <- exp(coef.slash[1,1] + coef.slash[2,1] * 
                                   mean.age + coef.slash[4,1] * mean.tlm + sqrt.sigma11 * q.stdslash_0.995)
text(2.9, est.armcirunf.female0.995 + 0.7, expression(paste(99.5, th)), 
     srt = 14, cex = 1.2)




