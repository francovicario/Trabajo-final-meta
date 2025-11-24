dat <- dat.raudenbush1985
dat

res <- rma(yi, vi, data = dat)
res

res_re <- rma(yi, vi, data=dat, method="REML")
summary(res_re) 

predict(res_re, digits=3) 

m.gen <- metagen(TE = yi, seTE = sqrt(vi), data=dat, sm="SMD", method.tau="REML", prediction=TRUE)
summary(m.gen)

m.gen.inf <- InfluenceAnalysis(m.gen, random = TRUE)
plot(m.gen.inf, "baujat")
plot(m.gen.inf, "influence")
plot(m.gen.inf, "ES") 
plot(m.gen.inf, "I2")

m.rma <- rma(yi=dat$yi, vi=dat$vi, method="REML", test="knha")
res.gosh <- gosh(m.rma)         
plot(res.gosh, alpha=0.03)
res.gosh.diag <- gosh.diagnostics(res.gosh, km.params=list(centers=2))
plot(res.gosh.diag)