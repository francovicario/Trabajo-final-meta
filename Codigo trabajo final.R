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
       print.subgroup.name = TRUE

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

radial(m.gen, main = "Radial plot – Raudenbush 1985")

### Radial con Egger
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





