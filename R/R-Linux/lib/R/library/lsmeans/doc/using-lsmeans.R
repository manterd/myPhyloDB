### R code from vignette source 'using-lsmeans.rnw'

###################################################
### code chunk number 1: using-lsmeans.rnw:49-51
###################################################
options(show.signif.stars=FALSE, prompt="R> ", continue="   ", 
    useFancyQuotes=FALSE, width=100, digits=6)


###################################################
### code chunk number 2: using-lsmeans.rnw:135-138
###################################################
library("lsmeans")
oranges.lm1 <- lm(sales1 ~ price1 + price2 + day + store, data = oranges)
anova(oranges.lm1)


###################################################
### code chunk number 3: using-lsmeans.rnw:141-142
###################################################
( oranges.rg1 <- ref.grid(oranges.lm1) )


###################################################
### code chunk number 4: using-lsmeans.rnw:145-146
###################################################
summary(oranges.rg1)


###################################################
### code chunk number 5: using-lsmeans.rnw:151-152
###################################################
lsmeans(oranges.rg1, "day")   ## or lsmeans(oranges.lm1, "day")


###################################################
### code chunk number 6: using-lsmeans.rnw:155-156
###################################################
with(oranges, tapply(sales1, day, mean))


###################################################
### code chunk number 7: using-lsmeans.rnw:165-167
###################################################
lsmeans(oranges.lm1, "day", at = list(price1 = 50, 
    price2 = c(40,60), day = c("2","3","4")) )


###################################################
### code chunk number 8: using-lsmeans.rnw:176-179
###################################################
org.lsm <- lsmeans(oranges.lm1, "day", by = "price2", 
    at = list(price1 = 50, price2 = c(40,60), day = c("2","3","4")) )
org.lsm


###################################################
### code chunk number 9: using-lsmeans.rnw:182-185 (eval = FALSE)
###################################################
## lsmeans(oranges.lm1, ~ day | price, at = ... )         # Ex 1
## lsmeans(oranges.lm1, c("day","price2"), at = ... )     # Ex 2
## lsmeans(oranges.lm1, ~ day * price, at = ... )         # Ex 3


###################################################
### code chunk number 10: using-lsmeans.rnw:195-196
###################################################
str(org.lsm)


###################################################
### code chunk number 11: using-lsmeans.rnw:202-204
###################################################
( org.sum <- summary(org.lsm, infer = c(TRUE,TRUE), 
                    level = .90, adjust = "bon", by = "day") )


###################################################
### code chunk number 12: using-lsmeans.rnw:209-210
###################################################
class(org.sum)


###################################################
### code chunk number 13: using-lsmeans.rnw:214-215
###################################################
transform(org.sum, lsrubles = lsmean * 34.2)


###################################################
### code chunk number 14: using-lsmeans.rnw:223-225
###################################################
org.lsm2 <- update(org.lsm, by.vars = NULL, level = .99)
org.lsm2


###################################################
### code chunk number 15: org-plot
###################################################
plot(org.lsm, by = "price2")


###################################################
### code chunk number 16: using-lsmeans.rnw:247-248
###################################################
contrast(org.lsm, method = "eff")


###################################################
### code chunk number 17: using-lsmeans.rnw:253-255
###################################################
days.lsm <- lsmeans(oranges.rg1, "day")
( days_contr.lsm <- contrast(days.lsm, "trt.vs.ctrl", ref = c(5,6)) )


###################################################
### code chunk number 18: using-lsmeans.rnw:260-261 (eval = FALSE)
###################################################
## confint(contrast(days.lsm, "trt.vs.ctrlk"))


###################################################
### code chunk number 19: using-lsmeans.rnw:269-270
###################################################
pairs(org.lsm)


###################################################
### code chunk number 20: using-lsmeans.rnw:273-274
###################################################
cld(days.lsm, alpha = .10)


###################################################
### code chunk number 21: days-cmp
###################################################
plot(days.lsm, comparisons = TRUE, alpha = .10)


###################################################
### code chunk number 22: using-lsmeans.rnw:297-298
###################################################
rbind(pairs(org.lsm)[1:3], pairs(org.lsm, by = "day")[1])


###################################################
### code chunk number 23: using-lsmeans.rnw:303-304
###################################################
rbind(pairs(lsmeans(org.lsm, "day")), pairs(lsmeans(org.lsm, "price2")))


###################################################
### code chunk number 24: using-lsmeans.rnw:310-313
###################################################
oranges.mlm <- lm(cbind(sales1,sales2) ~ price1 + price2 + day + store, 
                 data = oranges)
ref.grid(oranges.mlm)


###################################################
### code chunk number 25: using-lsmeans.rnw:316-318
###################################################
org.mlsm <- lsmeans(oranges.mlm, ~ day | variety, mult.name = "variety")
cld(org.mlsm, sort = FALSE)


###################################################
### code chunk number 26: using-lsmeans.rnw:326-327
###################################################
org.vardiff <- update(pairs(org.mlsm, by = "day"), by = NULL)


###################################################
### code chunk number 27: using-lsmeans.rnw:330-331
###################################################
cld(org.vardiff)


###################################################
### code chunk number 28: using-lsmeans.rnw:336-337
###################################################
contrast(org.mlsm, interaction = c("poly", "pairwise"))


###################################################
### code chunk number 29: using-lsmeans.rnw:344-346
###################################################
# Ensure we see the same results each time
set.seed(123454321)


###################################################
### code chunk number 30: using-lsmeans.rnw:348-351
###################################################
library("multcomp")
days.glht <- as.glht(days_contr.lsm)
summary(days.glht, test = adjusted("Westfall"))


###################################################
### code chunk number 31: using-lsmeans.rnw:354-356 (eval = FALSE)
###################################################
## days.glht1 <- glht(oranges.lm1, 
##                    lsm("day", contr = "trt.vs.ctrl", ref = c(5,6)))


###################################################
### code chunk number 32: using-lsmeans.rnw:360-362 (eval = FALSE)
###################################################
## summary(days_contr.lsm, adjust = "mvt")
## summary(days.glht)


###################################################
### code chunk number 33: using-lsmeans.rnw:368-369 (eval = FALSE)
###################################################
## summary(as.glht(pairs(org.lsm)))


###################################################
### code chunk number 34: using-lsmeans.rnw:372-373 (eval = FALSE)
###################################################
## summary(as.glht(pairs(org.lsm), by = NULL))


###################################################
### code chunk number 35: using-lsmeans.rnw:376-377 (eval = FALSE)
###################################################
## summary(as.glht(pairs(org.lsm, by = NULL)))


###################################################
### code chunk number 36: using-lsmeans.rnw:389-394
###################################################
data("Oats", package = "nlme")
library("lme4")
Oats.lmer <- lmer(log(yield) ~ Variety*factor(nitro) + (1|Block/Variety), 
                 data = Oats)
anova(Oats.lmer)


###################################################
### code chunk number 37: oatcontr (eval = FALSE)
###################################################
## ### We disable pbkrtest temporarily due to a bug in version 0.4.5.
## lsm.options(disable.pbkrtest = TRUE) 
## contrast(lsmeans(Oats.lmer, "nitro"), "poly")


###################################################
### code chunk number 38: using-lsmeans.rnw:403-404
###################################################
cat("NOTE: Results may be misleading due to involvement in interactions")


###################################################
### code chunk number 39: using-lsmeans.rnw:406-407
###################################################
### We disable pbkrtest temporarily due to a bug in version 0.4.5.
lsm.options(disable.pbkrtest = TRUE) 
contrast(lsmeans(Oats.lmer, "nitro"), "poly")


###################################################
### code chunk number 40: using-lsmeans.rnw:410-412
###################################################
Oats.lmer2 <- lmer(log(yield) ~ Variety + poly(nitro,2) 
                                + (1|Block/Variety),  data = Oats)


###################################################
### code chunk number 41: using-lsmeans.rnw:415-417
###################################################
Oats.lsm2 <- lsmeans(Oats.lmer2, ~ nitro | Variety, cov.reduce = FALSE)
Oats.lsm2


###################################################
### code chunk number 42: oatslmer
###################################################
lsmip(Oats.lmer, Variety ~ nitro, ylab = "Observed log(yield)")


###################################################
### code chunk number 43: oatslmer2
###################################################
lsmip(Oats.lsm2, Variety ~ nitro, ylab = "Predicted log(yield)")


###################################################
### code chunk number 44: using-lsmeans.rnw:458-459
###################################################
str(Oats.lsm2)


###################################################
### code chunk number 45: using-lsmeans.rnw:462-463
###################################################
summary(Oats.lsm2, type = "response")


###################################################
### code chunk number 46: using-lsmeans.rnw:479-481
###################################################
Oats.Vlsm = lsmeans(Oats.lmer2, "Variety")
test(Oats.Vlsm, null = log(100), type = "response")


###################################################
### code chunk number 47: using-lsmeans.rnw:493-494
###################################################
test(Oats.Vlsm, null = log(100), delta = 0.20, type = "r")


###################################################
### code chunk number 48: using-lsmeans.rnw:500-501
###################################################
test(contrast(Oats.Vlsm, "trt.vs.ctrlk"), side = ">")


###################################################
### code chunk number 49: using-lsmeans.rnw:504-505
###################################################
test(contrast(Oats.Vlsm, "trt.vs.ctrlk"), side = "nonsup", delta = .25)


###################################################
### code chunk number 50: chick-plot
###################################################
require("lattice")
xyplot(weight~Time | Diet, groups = ~ Chick, data = ChickWeight, 
    type = "o", layout=c(4, 1))


###################################################
### code chunk number 51: using-lsmeans.rnw:526-528
###################################################
Chick.lmer <- lmer(weight ~ Diet * Time + (0 + Time | Chick), 
    data = ChickWeight)


###################################################
### code chunk number 52: using-lsmeans.rnw:531-532
###################################################
( Chick.lst <- lstrends (Chick.lmer, ~ Diet, var = "Time") )


###################################################
### code chunk number 53: using-lsmeans.rnw:535-536
###################################################
cld (Chick.lst)


###################################################
### code chunk number 54: using-lsmeans.rnw:547-550
###################################################
lsm.options(ref.grid = list(level = .90),
            lsmeans = list(),
            contrast = list(infer = c(TRUE, TRUE)))


###################################################
### code chunk number 55: using-lsmeans.rnw:555-556
###################################################
get.lsm.option("estble.tol")


###################################################
### code chunk number 56: using-lsmeans.rnw:567-568
###################################################
lsmeans(Oats.lmer2, pairwise ~ Variety)


###################################################
### code chunk number 57: using-lsmeans.rnw:572-573
###################################################
lsm.options(ref.grid = NULL, contrast = NULL)


###################################################
### code chunk number 58: using-lsmeans.rnw:583-586
###################################################
nutr.lm <- lm(gain ~ (age + group + race)^2, data = nutrition)
library("car")
Anova(nutr.lm)


###################################################
### code chunk number 59: nutr-intplot
###################################################
lsmip(nutr.lm, race ~ age | group)
lsmeans(nutr.lm, ~ group*race)


###################################################
### code chunk number 60: using-lsmeans.rnw:603-605
###################################################
nutr.lsm <- lsmeans(nutr.lm, ~ group * race, weights = "proportional",
    at = list(age = c("2","3"), race = c("Black","White")))


###################################################
### code chunk number 61: using-lsmeans.rnw:608-611
###################################################
nutr.lsm    
summary(pairs(nutr.lsm, by = "race"), by = NULL)
summary(pairs(nutr.lsm, by = "group"), by = NULL)


###################################################
### code chunk number 62: using-lsmeans.rnw:624-628
###################################################
lsmeans(nutr.lm, "race", weights = "equal")
lsmeans(nutr.lm, "race", weights = "prop")
lsmeans(nutr.lm, "race", weights = "outer")
lsmeans(nutr.lm, "race", weights = "cells")


###################################################
### code chunk number 63: using-lsmeans.rnw:637-639
###################################################
temp = lsmeans(nutr.lm, c("group","race"), weights = "prop")
lsmeans(temp, "race", weights = "prop")


###################################################
### code chunk number 64: using-lsmeans.rnw:644-645
###################################################
with(nutrition, tapply(gain, race, mean))


###################################################
### code chunk number 65: using-lsmeans.rnw:652-656
###################################################
library("mediation")
levels(framing$educ) = c("NA","Ref","< HS", "HS", "> HS","Coll +")
framing.glm = glm(cong_mesg ~ age + income + educ + emo + gender * factor(treat),
                  family = binomial, data = framing)


###################################################
### code chunk number 66: framinga
###################################################
lsmip(framing.glm, treat ~ educ | gender, type = "response")


###################################################
### code chunk number 67: framingb
###################################################
lsmip(framing.glm, treat ~ educ | gender, type = "response",
      cov.reduce = emo ~ treat*gender + age + educ + income)


###################################################
### code chunk number 68: using-lsmeans.rnw:683-685
###################################################
ref.grid(framing.glm, 
    cov.reduce = emo ~ treat*gender + age + educ + income)@grid


###################################################
### code chunk number 69: using-lsmeans.rnw:716-718 (eval = FALSE)
###################################################
## rg <- ref.grid(my.model, at = list(x1 = c(5,10,15)),
##                cov.reduce = list(x2 ~ x1,  x3 ~ x1 + x2))


###################################################
### code chunk number 70: housing-plot
###################################################
library("ordinal")
data(housing, package = "MASS")
housing.clm <- clm(Sat ~ (Infl + Type + Cont)^2,
                   data = housing, weights = Freq, link = "probit")
lsmip(housing.clm, Cont ~ Infl | Type, layout = c(4,1))


###################################################
### code chunk number 71: using-lsmeans.rnw:763-764
###################################################
test(pairs(lsmeans(housing.clm, ~ Infl | Type)), joint = TRUE)


###################################################
### code chunk number 72: using-lsmeans.rnw:767-768
###################################################
test(pairs(lsmeans(housing.clm, ~ Cont | Type)), joint = TRUE)


###################################################
### code chunk number 73: using-lsmeans.rnw:773-774
###################################################
ref.grid(housing.clm, mode = "cum.prob")


###################################################
### code chunk number 74: using-lsmeans.rnw:777-779
###################################################
lsmeans(housing.clm, ~ Infl, at = list(cut = "Medium|High"), 
        mode = "cum.prob")


###################################################
### code chunk number 75: using-lsmeans.rnw:782-784
###################################################
summary(lsmeans(housing.clm, ~ Infl, at = list(cut = "Medium|High"), 
                mode = "linear.predictor"), type = "response")


###################################################
### code chunk number 76: using-lsmeans.rnw:792-800
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
### code chunk number 77: using-lsmeans.rnw:803-805
###################################################
cld(lsmeans(Chick.nlme, ~ Diet, param = "asym"))    
cld(lsmeans(Chick.nlme, ~ Diet, param = "xmid"))    


###################################################
### code chunk number 78: using-lsmeans.rnw:817-822
###################################################
library("MCMCpack")
counts <- c(18, 17, 15,   20, 10, 20,   25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
posterior <- MCMCpoisson(counts ~ outcome + treatment, mcmc = 1000)


###################################################
### code chunk number 79: using-lsmeans.rnw:825-826
###################################################
( post.lsm <- lsmeans(posterior, "treatment") )


###################################################
### code chunk number 80: using-lsmeans.rnw:829-831
###################################################
library("coda")
summary(as.mcmc(post.lsm))


