## mirt:
library(mirt)
ls7 <- mirt::expand.table(LSAT7)
mirtmod <- mirt(ls7[,1:5], 1, itemtype = "Rasch", SE = TRUE)


## reshape data and fit with glmer():
library(reshape2)
library(lme4)
ls7$person <- 1:nrow(ls7)
ls7long <- melt(ls7, id = "person")
lme4mod <- glmer(value ~ -1 + variable + (1 | person), family = binomial,
                 data = ls7long, nAGQ = 5L) #random effect: latent var; fixed effects: parameters (here: thresholds...)

## fit with lavaan WLS estimator
library(lavaan)
wlsmod <- cfa("f1 =~ 1*Item.1 + 1*Item.2 + 1*Item.3 + 1*Item.4 + 1*Item.5", estimator = "WLS", data=ls7, ordered = TRUE)

## score calculation:
library(merDeriv)
mirtsc <- estfun.AllModelClass(mirtmod)
lme4sc <- estfun.glmerMod(lme4mod, ranpar = "var")
