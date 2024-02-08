#multi
model = '
  Eta1 =~ simuvar1 + simuvar2 + simuvar3 
  Eta2 =~ simuvar4 + simuvar5 + simuvar6 
  Eta3 =~ simuvar7 + simuvar8 + simuvar9'
fits_random <- datagen(model = model, schwellen = 6, ID=1000, times=1, items=3, latvar = 3) #multi
simu=fits_random[["data"]][["data1"]]
fit_ord <- lavaan::cfa(model, data = simu, ordered = TRUE, estimator = "WLS",std.lv=TRUE )
object=fit_ord
summary(fit_ord)

#univ
fits_random <- datagen(schwellen = 4, ID=1000, times=1, items=5) #univ
simu=fits_random[["data"]][["data1"]]
fit_ord <- lavaan::cfa(fits_random[["model"]][[1]], data = simu, ordered = TRUE, estimator = "WLS",std.lv=TRUE )
object=fit_ord



################
model=fits_random[["model"]][["model1"]]
fluc_ord <- strucchange::gefp(fit_ord, fit=NULL, scores = estfun.WLS, order.by = simu[,"rand1"]) #ordinal/categorical
fluc_num <- strucchange::gefp(fit_ord, fit=NULL, scores = estfun.WLS, order.by = simu[,"rand2"]) #numeric
fluc_cat <- strucchange::gefp(fit_ord, fit=NULL, scores = estfun.WLS, order.by = simu[,"rand3"]) #nominal
p_maxBB    <- strucchange::sctest(fluc_num, functional = strucchange::maxBB)$p.value #dm
p_meanL2BB <- strucchange::sctest(fluc_num, functional = strucchange::meanL2BB)$p.value #cvm 
p_sumlm    <- strucchange::sctest(fluc_num, functional=  strucchange::supLM(from = 0.15, to = NULL))$p.value #suplm
p_ordL2BB  <- strucchange::sctest(fluc_ord, functional = strucchange::ordL2BB(fluc_ord))$p.value[1]  #LM #istead of fluc1 ... factor(simu[,"rand1"])
p_ordwmax  <- strucchange::sctest(fluc_ord, functional = strucchange::ordwmax(fluc_ord))$p.value[1]  #wdmo 
p_catl2bb  <- strucchange::sctest(fluc_cat, functional = strucchange::catL2BB(fluc_cat))$p.value  #LMuo
  