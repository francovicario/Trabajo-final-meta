library(meta)
library(metafor)
library(dmetar)

#Heterogeneidad

dat <- dat.raudenbush1985
dat

dat$seTE <- sqrt(dat$vi)
m.gen <- metagen(TE = yi,
                 seTE = seTE,
                 data = dat,
                 studlab = paste(study, author, sep=": "),
                 sm = "SMD",
                 method.tau = "REML",    
                 prediction = TRUE)     
summary(m.gen)

forest(m.gen, leftcols = c("studlab"), rightcols = c("effect","ci"))
funnel(m.gen)
rma_obj <- rma(yi = dat$yi, vi = dat$vi, method = "REML")
regtest(rma_obj, model = "rma", predictor = "sei")

outliers <- find.outliers(m.gen)
print(outliers) 

m.gen.inf <- InfluenceAnalysis(m.gen, random = TRUE)
plot(m.gen.inf, "baujat")     
plot(m.gen.inf, "influence")  
plot(m.gen.inf, "ES")         
plot(m.gen.inf, "I2") 

m.rma <- rma(yi = dat$yi, sei = dat$seTE, method = "REML", test = "knha")
set.seed(123)
res.gosh <- gosh(m.rma)
plot(res.gosh, alpha = 0.02)

res.gosh.diagnos <- gosh.diagnostics(res.gosh, km.params = list(centers = 2), db.params = list(eps=0.08, MinPts=50))
res.gosh.diagnos
plot(res.gosh.diagnos)

#subgrupos

# A) Subgrupo por 'setting'
#Primero con tau.common = FALSE (estimaciones separadas de tau^2 por subgrupo)
m.setting_separateTau <- update(m.gen, subgroup = setting, tau.common = FALSE)
summary(m.setting_separateTau)

#Luego con tau.common = TRUE (misma estimación tau^2 para todos los subgrupos)
m.setting_commonTau <- update(m.gen, subgroup = setting, tau.common = TRUE)
summary(m.setting_commonTau)

# Forest plot por subgrupos (ajustá xlim/cex si es necesario)
forest(m.setting_separateTau, xlab = "Diferencia de Media Estandarizada")
       print.subgroup.name = TRUE)

# B) Subgrupo por 'tester' (aware vs blind)
m.tester_separateTau <- update(m.gen, subgroup = tester, tau.common = FALSE)
summary(m.tester_separateTau)

m.tester_commonTau <- update(m.gen, subgroup = tester, tau.common = TRUE)
summary(m.tester_commonTau)

forest(m.tester_separateTau, xlab ="Diferencia de Media Estandarizada")

#metaregresion
mr.weeks <- metareg(m.gen, ~ weeks)
summary(mr.weeks)

bubble(mr.weeks,
       xlab = "Semanas",
       ylab = "Tamnho de Efecto)",
       studlab = TRUE) 

bubble(mr.weeks,
       xlab = "Weeks",
       ylab = "Effect size (SMD)",
       pred = TRUE)

#Sesgo

funnel(m.gen, main = "Funnel plot – Raudenbush 1985")

funnel(m.gen, main = "Funnel plot – Raudenbush 1985", contour = c(0.90, 0.95, 0.99),col.contour = c("darkgray", "gray", "lightgray"))

radial(m.gen, main = "Radial plot – Raudenbush 1985")

### Radial con Egger-
radial(m.gen, main = "Radial plot con línea de Egger")

### Crear variables de Egger
yi_over_se  <- m.gen$TE / m.gen$seTE      # yi / SEi
inv_se      <- 1 / m.gen$seTE             # 1 / SEi

### Ajustar modelo de Egger (regresión lineal)
egger_model <- lm(yi_over_se ~ inv_se)

### Agregar la línea de Egger al radial plot
abline(egger_model, lwd = 2)


metabias(m.gen, method = "rank")

metabias(m.gen, method = "linreg")

metabias(m.gen, method = "mm")

tf <- trimfill(m.gen)
summary(tf)
funnel(tf, main = "Trim and fill – Raudenbush 1985")

#datos faltantes
#al no haber datos faltantes seguimos el procedimiento del libro y modificamos los datos agregandole datos flatantes

set.seed(123)
dat2 <- dat

# columnas básicas
dat2$n1i <- as.integer(dat2$n1i)
dat2$n2i <- as.integer(dat2$n2i)
dat2$TE   <- dat2$yi
dat2$seTE <- sqrt(dat2$vi)

# Elegimos estudios más precisos
ord <- order(dat2$seTE)
sel7 <- ord[1:7]
sel_missing <- sel7[c(3,5)]

# Definir pérdidas
dat2$Nem <- 0
dat2$Ncm <- 0
dat2$Nem[sel_missing] <- c(8,5)
dat2$Ncm[sel_missing] <- dat2$Nem[sel_missing]

# Observados por brazo
dat2$Neo <- dat2$n1i - dat2$Nem
dat2$Nco <- dat2$n2i - dat2$Ncm

# Proporciones faltantes
dat2$Pe <- dat2$Nem / dat2$n1i
dat2$Pc <- dat2$Ncm / dat2$n2i

# Metaanálisis basado en los datos observados (pre-ajuste)
m.gen.miss <- metagen(
  TE = dat2$TE,
  seTE = dat2$seTE,
  data = dat2,
  studlab = paste(dat2$study, dat2$author, sep=": "),
  sm = "SMD",
  method.tau = "REML",
  prediction = TRUE
)

summary(m.gen.miss)

# semiss: función del libro Schwarzer et al. (2015), cap. 6

semiss <- function(TE, seTE, mu.e, mu.c, nu.e, nu.c, rho){
  
  # Delta para cada estudio
  delta.e <- rnorm(length(TE), mean = mu.e, sd = nu.e)
  delta.c <- rnorm(length(TE), mean = mu.c, sd = nu.c)
  
  # Correlación entre brazos
  if (rho != 0) {
    delta.c <- rho * delta.e + sqrt(1 - rho^2) * delta.c
  }
  
  # Ajustar TE
  TE.star <- TE + delta.e - delta.c
  
  # Ajustar varianza
  seTE.star <- sqrt(seTE^2 + nu.e^2 + nu.c^2 - 2*rho*nu.e*nu.c)
  
  return(list(TE = TE.star, seTE = seTE.star))
}

mu <- 8
nu <- sqrt(5)

# Convenience:
TE  <- dat2$yi
seT <- dat2$seTE

#fixed equal
S1 <- semiss(TE, seT, mu.e = mu, mu.c = mu, nu.e = 0, nu.c = 0, rho = 0)

m.S1 <- metagen(
  TE = S1$TE,
  seTE = S1$seTE,
  studlab = paste(dat2$study, dat2$author),
  sm = "SMD",
  method.tau = "REML"
)
summary(m.S1)

#fixed opposite
S2 <- semiss(TE, seT, mu.e = mu, mu.c = -mu, nu.e = 0, nu.c = 0, rho = 0)

m.S2 <- metagen(
  TE = S2$TE,
  seTE = S2$seTE,
  studlab = paste(dat2$study, dat2$author),
  sm = "SMD",
  method.tau = "REML"
)
summary(m.S2)

#random equal
S3 <- semiss(TE, seT, mu.e = mu, mu.c = mu, nu.e = nu, nu.c = nu, rho = 1)

m.S3 <- metagen(
  TE = S3$TE,
  seTE = S3$seTE,
  studlab = paste(dat2$study, dat2$author),
  sm = "SMD",
  method.tau = "REML"
)
summary(m.S3)

#random equal
S4 <- semiss(TE, seT, mu.e = mu, mu.c = -mu, nu.e = nu, nu.c = nu, rho = 0)

m.S4 <- metagen(
  TE = S4$TE,
  seTE = S4$seTE,
  studlab = paste(dat2$study, dat2$author),
  sm = "SMD",
  method.tau = "REML"
)
summary(m.S4)

#graficos
forest(m.S1, main="Missing Data — Fixed Equal")
forest(m.S2, main="Missing Data — Fixed Opposite")
forest(m.S3, main="Missing Data — Random Equal")
forest(m.S4, main="Missing Data — Random Uncorrelated")

funnel(m.S1,
       main = "Funnel plot — Missing Data: Fixed Equal",
       xlab = "Effect size (SMD)")

funnel(m.S2,
       main = "Funnel plot — Missing Data: Fixed Opposite",
       xlab = "Effect size (SMD)")

funnel(m.S3,
       main = "Funnel plot — Missing Data: Random Equal",
       xlab = "Effect size (SMD)")

funnel(m.S4,
       main = "Funnel plot — Missing Data: Random Uncorrelated",
       xlab = "Effect size (SMD)")

#potencia efecto
library(metapower)

n_exp  <- mean((dat$n1i + dat$n2i)/2)

my_power1<- mpower(effect_size = .1, study_size =n_exp, k = 19,i2 = .50, es_type = "d")
my_power1
plot_mpower(my_power1)

my_power2<- mpower(effect_size = .2, study_size =n_exp, k = 19, i2 = .50, es_type = "d")
my_power2
plot_mpower(my_power2)

my_power3<- mpower(effect_size = .1, study_size =n_exp, k = 19, i2 = .75, es_type = "d")
my_power3
plot_mpower(my_power3)

#potencia homogeneidad

my_homogen_power1<-  homogen_power(effect_size = .10, study_size = n_exp, , k = 19, i2 = .50, es_type = "d")
my_homogen_power1

my_homogen_power <- homogen_power(effect_size = .20, study_size = n_exp, , k = 19, i2 = .50, es_type = "d")
my_homogen_power

my_homogen_power3<- homogen_power(effect_size = .10, study_size = n_exp, , k = 19, i2 = .60, es_type = "d")
my_homogen_power3

my_homogen_power4<- homogen_power(effect_size = .10, study_size = n_exp, , k = 19, i2 = .25, es_type = "d")
my_homogen_power4



