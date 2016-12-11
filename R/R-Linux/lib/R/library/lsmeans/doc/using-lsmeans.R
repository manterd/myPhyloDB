### R code from vignette source 'using-lsmeans.rnw'

###################################################
### code chunk number 1: using-lsmeans.rnw:69-71
###################################################
options(show.signif.stars=FALSE, prompt="R> ", continue="   ", 
    useFancyQuotes=FALSE, width=100, digits=6)


###################################################
### code chunk number 2: using-lsmeans.rnw:155-158
###################################################
library("lsmeans")
oranges.lm1 <- lm(sales1 ~ price1 + price2 + day + store, data = oranges)
anova(oranges.lm1)


###################################################
### code chunk number 3: using-lsmeans.rnw:161-162
###################################################
( oranges.rg1 <- ref.grid(oranges.lm1) )


###################################################
### code chunk number 4: using-lsmeans.rnw:165-166
###################################################
summary(oranges.rg1)


###################################################
### code chunk number 5: using-lsmeans.rnw:171-172
###################################################
lsmeans(oranges.rg1, "day")   ## or lsmeans(oranges.lm1, "day")


###################################################
### code chunk number 6: using-lsmeans.rnw:175-176
###################################################
with(oranges, tapply(sales1, day, mean))


###################################################
### code chunk number 7: using-lsmeans.rnw:187-189
###################################################
lsmeans(oranges.lm1, "day", at = list(price1 = 50, 
    price2 = c(40,60), day = c("2","3","4")) )


###################################################
### code chunk number 8: using-lsmeans.rnw:198-201
###################################################
org.lsm <- lsmeans(oranges.lm1, "day", by = "price2", 
    at = list(price1 = 50, price2 = c(40,60), day = c("2","3","4")) )
org.lsm


###################################################
### code chunk number 9: using-lsmeans.rnw:204-207 (eval = FALSE)
###################################################
## lsmeans(oranges.lm1, ~ day | price, at = ... )         # Ex 1
## lsmeans(oranges.lm1, c("day","price2"), at = ... )     # Ex 2
## lsmeans(oranges.lm1, ~ day * price, at = ... )         # Ex 3


###################################################
### code chunk number 10: using-lsmeans.rnw:217-218
###################################################
str(org.lsm)


###################################################
### code chunk number 11: using-lsmeans.rnw:224-226
###################################################
( org.sum <- summary(org.lsm, infer = c(TRUE,TRUE), 
                    level = .90, adjust = "bon", by = "day") )


###################################################
### code chunk number 12: using-lsmeans.rnw:231-232
###################################################
class(org.sum)


###################################################
### code chunk number 13: using-lsmeans.rnw:236-237
###################################################
transform(org.sum, lsrubles = lsmean * 34.2)


###################################################
### code chunk number 14: using-lsmeans.rnw:245-247
###################################################
org.lsm2 <- update(org.lsm, by.vars = NULL, level = .99)
org.lsm2


###################################################
### code chunk number 15: org-plot
###################################################
plot(org.lsm, by = "price2")


###################################################
### code chunk number 16: using-lsmeans.rnw:270-271
###################################################
contrast(org.lsm, method = "eff")


###################################################
### code chunk number 17: using-lsmeans.rnw:276-278
###################################################
days.lsm <- lsmeans(oranges.rg1, "day")
( days_contr.lsm <- contrast(days.lsm, "trt.vs.ctrl", ref = c(5,6)) )


###################################################
### code chunk number 18: using-lsmeans.rnw:283-284 (eval = FALSE)
###################################################
## confint(contrast(days.lsm, "trt.vs.ctrlk"))


###################################################
### code chunk number 19: using-lsmeans.rnw:292-293
###################################################
pairs(org.lsm)


###################################################
### code chunk number 20: using-lsmeans.rnw:296-297
###################################################
cld(days.lsm, alpha = .10)


###################################################
### code chunk number 21: days-cmp
###################################################
plot(days.lsm, comparisons = TRUE, alpha = .10)


###################################################
### code chunk number 22: using-lsmeans.rnw:324-325
###################################################
rbind(pairs(org.lsm)[1:3], pairs(org.lsm, by = "day")[1])


###################################################
### code chunk number 23: using-lsmeans.rnw:330-331
###################################################
rbind(pairs(lsmeans(org.lsm, "day")), pairs(lsmeans(org.lsm, "price2")))


###################################################
### code chunk number 24: using-lsmeans.rnw:337-340
###################################################
oranges.mlm <- lm(cbind(sales1,sales2) ~ price1 + price2 + day + store, 
                 data = oranges)
ref.grid(oranges.mlm)


###################################################
### code chunk number 25: using-lsmeans.rnw:343-345
###################################################
org.mlsm <- lsmeans(oranges.mlm, ~ day | variety, mult.name = "variety")
cld(org.mlsm, sort = FALSE)


###################################################
### code chunk number 26: using-lsmeans.rnw:353-354
###################################################
org.vardiff <- update(pairs(org.mlsm, by = "day"), by = NULL)


###################################################
### code chunk number 27: using-lsmeans.rnw:357-358
###################################################
cld(org.vardiff)


###################################################
### code chunk number 28: using-lsmeans.rnw:364-365
###################################################
contrast(org.mlsm, interaction = c("poly", "pairwise"))


###################################################
### code chunk number 29: using-lsmeans.rnw:373-375
###################################################
# Ensure we see the same results each time
set.seed(123454321)


###################################################
### code chunk number 30: using-lsmeans.rnw:377-380
###################################################
library("multcomp")
days.glht <- as.glht(days_contr.lsm)
summary(days.glht, test = adjusted("Westfall"))


###################################################
### code chunk number 31: using-lsmeans.rnw:383-385 (eval = FALSE)
###################################################
## days.glht1 <- glht(oranges.lm1, 
##                    lsm("day", contr = "trt.vs.ctrl", ref = c(5,6)))


###################################################
### code chunk number 32: using-lsmeans.rnw:389-391 (eval = FALSE)
###################################################
## summary(days_contr.lsm, adjust = "mvt")
## summary(days.glht)


###################################################
### code chunk number 33: using-lsmeans.rnw:397-398 (eval = FALSE)
###################################################
## summary(as.glht(pairs(org.lsm)))


###################################################
### code chunk number 34: using-lsmeans.rnw:401-402 (eval = FALSE)
###################################################
## summary(as.glht(pairs(org.lsm), by = NULL))


###################################################
### code chunk number 35: using-lsmeans.rnw:405-406 (eval = FALSE)
###################################################
## summary(as.glht(pairs(org.lsm, by = NULL)))


###################################################
### code chunk number 36: using-lsmeans.rnw:418-423
###################################################
data("Oats", package = "nlme")
library("lme4")
Oats.lmer <- lmer(log(yield) ~ Variety*factor(nitro) + (1|Block/Variety), 
                 data = Oats)
anova(Oats.lmer)


###################################################
### code chunk number 37: oatcontr (eval = FALSE)
###################################################
## contrast(lsmeans(Oats.lmer, "nitro"), "poly")


###################################################
### code chunk number 38: using-lsmeans.rnw:430-431
###################################################
cat("NOTE: Results may be misleading due to involvement in interactions")


###################################################
### code chunk number 39: using-lsmeans.rnw:433-434
###################################################
contrast(lsmeans(Oats.lmer, "nitro"), "poly")


###################################################
### code chunk number 40: using-lsmeans.rnw:438-440
###################################################
Oats.lmer2 <- lmer(log(yield) ~ Variety + poly(nitro,2) 
                                + (1|Block/Variety),  data = Oats)


###################################################
### code chunk number 41: using-lsmeans.rnw:444-445
###################################################
Oats.lsm2 <- lsmeans(Oats.lmer2, ~ nitro | Variety, cov.reduce = FALSE)


###################################################
### code chunk number 42: using-lsmeans.rnw:457-462
###################################################
library("xtable")
xtbl <- xtable(Oats.lsm2, caption = "Example using \\texttt{xtable}",
    label = "xtable:example")
print(xtbl, table.placement = "t")    
cat("See Table~\\ref{xtable:example}.\n")


###################################################
### code chunk number 43: oatslmer
###################################################
lsmip(Oats.lmer, Variety ~ nitro, ylab = "Observed log(yield)")


###################################################
### code chunk number 44: oatslmer2
###################################################
lsmip(Oats.lsm2, Variety ~ nitro, ylab = "Predicted log(yield)")


###################################################
### code chunk number 45: using-lsmeans.rnw:498-499
###################################################
str(Oats.lsm2)


###################################################
### code chunk number 46: using-lsmeans.rnw:502-503
###################################################
summary(Oats.lsm2, type = "response")


###################################################
### code chunk number 47: using-lsmeans.rnw:514-518
###################################################
Oats.log1 <- lmer(log(yield + 5) ~ Variety + factor(nitro) 
                  + (1|Block/Variety), data = Oats)
( Oats.rg1 <- update(ref.grid(Oats.log1), 
                    tran = make.tran("genlog", 5)) )


###################################################
### code chunk number 48: using-lsmeans.rnw:523-524
###################################################
round(predict(Oats.rg1, type = "response"), 1)


###################################################
### code chunk number 49: using-lsmeans.rnw:529-531
###################################################
my.tran <- make.tran("boxcox", c(.567, 10))
my.tran$linkfun(10:15)


###################################################
### code chunk number 50: using-lsmeans.rnw:539-543
###################################################
Oats.bc <- with(my.tran, lmer(linkfun(yield) ~ Variety + factor(nitro)
                              + (1|Block/Variety), data = Oats))
( rg.bc <- ref.grid(Oats.bc) )
round(predict(rg.bc, type = "response"), 1)


###################################################
### code chunk number 51: using-lsmeans.rnw:548-549
###################################################
rg.bc.regrid <- regrid(rg.bc)


###################################################
### code chunk number 52: using-lsmeans.rnw:552-553
###################################################
round(rg.bc.regrid@bhat, 1)


###################################################
### code chunk number 53: using-lsmeans.rnw:558-560
###################################################
summary(lsmeans(rg.bc, "Variety"), type = "response")
lsmeans(rg.bc.regrid, "Variety")


###################################################
### code chunk number 54: using-lsmeans.rnw:579-583
###################################################
rg.log <- regrid(rg.bc, "log")
lsm.log <- lsmeans(rg.log, "Variety")
summary(lsm.log, type = "response")
summary(pairs(lsm.log), type = "response")


###################################################
### code chunk number 55: using-lsmeans.rnw:598-600
###################################################
Oats.Vlsm = lsmeans(Oats.lmer2, "Variety")
test(Oats.Vlsm, null = log(100), type = "response")


###################################################
### code chunk number 56: using-lsmeans.rnw:612-613
###################################################
test(Oats.Vlsm, null = log(100), delta = 0.20, type = "r")


###################################################
### code chunk number 57: using-lsmeans.rnw:620-621
###################################################
test(contrast(Oats.Vlsm, "trt.vs.ctrlk"), side = ">")


###################################################
### code chunk number 58: using-lsmeans.rnw:625-626
###################################################
test(contrast(Oats.Vlsm, "trt.vs.ctrlk"), side = "nonsup", delta = .25)


###################################################
### code chunk number 59: chick-plot
###################################################
require("lattice")
xyplot(weight ~ Time | Diet, groups = ~ Chick, data = ChickWeight, 
    type = "o", layout=c(4, 1))


###################################################
### code chunk number 60: using-lsmeans.rnw:648-650
###################################################
Chick.lmer <- lmer(sqrt(weight) ~ Diet * Time + (0 + Time | Chick), 
    data = ChickWeight)


###################################################
### code chunk number 61: using-lsmeans.rnw:653-654
###################################################
Chick.lst <- lstrends (Chick.lmer, ~ Diet, var = "Time")


###################################################
### code chunk number 62: using-lsmeans.rnw:657-658
###################################################
cld (Chick.lst)


###################################################
### code chunk number 63: using-lsmeans.rnw:663-665
###################################################
lstrends(Chick.lmer, ~ Diet | Time, var = "Time", 
    transform = "response", at = list(Time = c(5, 15)))


###################################################
### code chunk number 64: using-lsmeans.rnw:677-680
###################################################
lsm.options(ref.grid = list(level = .90),
            lsmeans = list(),
            contrast = list(infer = c(TRUE, TRUE)))


###################################################
### code chunk number 65: using-lsmeans.rnw:685-686
###################################################
get.lsm.option("estble.tol")


###################################################
### code chunk number 66: using-lsmeans.rnw:703-704
###################################################
lsmeans(Oats.lmer2, pairwise ~ Variety)


###################################################
### code chunk number 67: using-lsmeans.rnw:708-709
###################################################
lsm.options(ref.grid = NULL, contrast = NULL)


###################################################
### code chunk number 68: using-lsmeans.rnw:719-722
###################################################
nutr.lm <- lm(gain ~ (age + group + race)^2, data = nutrition)
library("car")
Anova(nutr.lm)


###################################################
### code chunk number 69: nutr-intplot
###################################################
lsmip(nutr.lm, race ~ age | group)
lsmeans(nutr.lm, ~ group*race)


###################################################
### code chunk number 70: using-lsmeans.rnw:739-741
###################################################
nutr.lsm <- lsmeans(nutr.lm, ~ group * race, weights = "proportional",
    at = list(age = c("2","3"), race = c("Black","White")))


###################################################
### code chunk number 71: using-lsmeans.rnw:744-747
###################################################
nutr.lsm    
summary(pairs(nutr.lsm, by = "race"), by = NULL)
summary(pairs(nutr.lsm, by = "group"), by = NULL)


###################################################
### code chunk number 72: using-lsmeans.rnw:760-764
###################################################
lsmeans(nutr.lm, "race", weights = "equal")
lsmeans(nutr.lm, "race", weights = "prop")
lsmeans(nutr.lm, "race", weights = "outer")
lsmeans(nutr.lm, "race", weights = "cells")


###################################################
### code chunk number 73: using-lsmeans.rnw:773-775
###################################################
temp = lsmeans(nutr.lm, c("group","race"), weights = "prop")
lsmeans(temp, "race", weights = "prop")


###################################################
### code chunk number 74: using-lsmeans.rnw:780-781
###################################################
with(nutrition, tapply(gain, race, mean))


###################################################
### code chunk number 75: using-lsmeans.rnw:789-793
###################################################
library("mediation")
levels(framing$educ) = c("NA","Ref","< HS", "HS", "> HS","Coll +")
framing.glm = glm(cong_mesg ~ age + income + educ + emo + gender * factor(treat),
                  family = binomial, data = framing)


###################################################
### code chunk number 76: framinga
###################################################
lsmip(framing.glm, treat ~ educ | gender, type = "response")


###################################################
### code chunk number 77: framingb
###################################################
lsmip(framing.glm, treat ~ educ | gender, type = "response",
      cov.reduce = emo ~ treat*gender + age + educ + income)


###################################################
### code chunk number 78: using-lsmeans.rnw:821-823
###################################################
ref.grid(framing.glm, 
    cov.reduce = emo ~ treat*gender + age + educ + income)@grid


###################################################
### code chunk number 79: using-lsmeans.rnw:854-856 (eval = FALSE)
###################################################
## rg <- ref.grid(my.model, at = list(x1 = c(5,10,15)),
##                cov.reduce = list(x2 ~ x1,  x3 ~ x1 + x2))


###################################################
### code chunk number 80: housing-plot
###################################################
library("ordinal")
data(housing, package = "MASS")
housing.clm <- clm(Sat ~ (Infl + Type + Cont)^2,
                   data = housing, weights = Freq, link = "probit")
lsmip(housing.clm, Cont ~ Infl | Type, layout = c(4,1))


###################################################
### code chunk number 81: using-lsmeans.rnw:908-909
###################################################
test(pairs(lsmeans(housing.clm, ~ Infl | Type)), joint = TRUE)


###################################################
### code chunk number 82: using-lsmeans.rnw:912-913
###################################################
test(pairs(lsmeans(housing.clm, ~ Cont | Type)), joint = TRUE)


###################################################
### code chunk number 83: using-lsmeans.rnw:918-919
###################################################
ref.grid(housing.clm, mode = "cum.prob")


###################################################
### code chunk number 84: using-lsmeans.rnw:922-924
###################################################
lsmeans(housing.clm, ~ Infl, at = list(cut = "Medium|High"), 
        mode = "cum.prob")


###################################################
### code chunk number 85: using-lsmeans.rnw:927-929
###################################################
summary(lsmeans(housing.clm, ~ Infl, at = list(cut = "Medium|High"), 
                mode = "linear.predictor"), type = "response")


###################################################
### code chunk number 86: using-lsmeans.rnw:937-945
###################################################
require("nlme")
options(contrasts = c("contr.treatment", "contr.poly"))
Chick.nlme = nlme(weight ~ SSlogis(Time, asym, xmid, scal), 
    data = ChickWeight,
    fixed = list(asym + xmid ~ Diet, scal ~ 1),
    random = asym ~ 1 | Chick, 
    start = c(200, 100, 200, 100,   10, 0, 0, 0,   7))
Chick.nlme


###################################################
### code chunk number 87: using-lsmeans.rnw:948-950
###################################################
cld(lsmeans(Chick.nlme, ~ Diet, param = "asym"))    
cld(lsmeans(Chick.nlme, ~ Diet, param = "xmid"))    


###################################################
### code chunk number 88: using-lsmeans.rnw:962-967
###################################################
library("MCMCpack")
counts <- c(18, 17, 15,   20, 10, 20,   25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
posterior <- MCMCpoisson(counts ~ outcome + treatment, mcmc = 1000)


###################################################
### code chunk number 89: using-lsmeans.rnw:970-971
###################################################
( post.lsm <- lsmeans(posterior, "treatment") )


###################################################
### code chunk number 90: using-lsmeans.rnw:974-976
###################################################
library("coda")
summary(as.mcmc(post.lsm))


