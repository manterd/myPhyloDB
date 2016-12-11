### R code from vignette source 'extending.rnw'

###################################################
### code chunk number 1: extending.rnw:42-44
###################################################
options(show.signif.stars=FALSE, prompt="R> ", continue="   ")
set.seed(271828)


###################################################
### code chunk number 2: extending.rnw:55-56 (eval = FALSE)
###################################################
## help("extending-lsmeans", package="lsmeans")


###################################################
### code chunk number 3: extending.rnw:64-68
###################################################
fake = expand.grid(rep = 1:5, A = c("a1","a2"), B = c("b1","b2","b3"))
fake$y = c(11.46,12.93,11.87,11.01,11.92,17.80,13.41,13.96,14.27,15.82,
           23.14,23.75,-2.09,28.43,23.01,24.11,25.51,24.11,23.95,30.37,
           17.75,18.28,17.82,18.52,16.33,20.58,20.55,20.77,21.21,20.10)


###################################################
### code chunk number 4: extending.rnw:74-79
###################################################
library(MASS)
fake.rlm = rlm(y ~ A * B, data = fake)

library(lsmeans)
lsmeans(fake.rlm, ~B | A)


###################################################
### code chunk number 5: extending.rnw:85-86
###################################################
fake.lts = ltsreg(y ~ A * B, data = fake)


###################################################
### code chunk number 6: extending.rnw:91-92
###################################################
lsmeans:::recover.data.lm


###################################################
### code chunk number 7: extending.rnw:95-96
###################################################
recover.data.lqs = lsmeans:::recover.data.lm


###################################################
### code chunk number 8: extending.rnw:99-101
###################################################
rec.fake = recover.data(fake.lts)
head(rec.fake)


###################################################
### code chunk number 9: extending.rnw:115-116
###################################################
args(lsmeans:::lsm.basis.lm)


###################################################
### code chunk number 10: extending.rnw:123-124
###################################################
MASS:::predict.lqs


###################################################
### code chunk number 11: extending.rnw:128-139
###################################################
lsm.basis.lqs = function(object, trms, xlev, grid, ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    bhat = coef(object)
    Xmat = model.matrix(trms, data=object$model)
    V = rev(object$scale)[1]^2 * solve(t(Xmat) %*% Xmat)
    nbasis = matrix(NA)
    dfargs = list(df = nrow(Xmat) - ncol(Xmat))
    dffun = function(k, dfargs) dfargs$df
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs)
}


###################################################
### code chunk number 12: extending.rnw:143-144
###################################################
lsmeans(fake.lts, ~ B | A)


###################################################
### code chunk number 13: extending.rnw:155-156 (eval = FALSE)
###################################################
## nbasis = estimability::nonest.basis(Xmat)


###################################################
### code chunk number 14: extending.rnw:177-180
###################################################
form = ~ data$x + data[[5]]
base::all.vars(form)
lsmeans::.all.vars(form)


###################################################
### code chunk number 15: extending.rnw:187-188
###################################################
.get.offset(terms(~ speed + offset(.03*breaks)), head(warpbreaks))


###################################################
### code chunk number 16: extending.rnw:205-229 (eval = FALSE)
###################################################
## recover.data.rsm = function(object, data, mode = c("asis", "coded", "decoded"), ...) {
##     mode = match.arg(mode)
##     cod = codings(object)
##     fcall = object$call
##     if(is.null(data))
##         data = lsmeans::recover.data(fcall, delete.response(terms(object)), object$na.action, ...)
##     if (!is.null(cod) && (mode == "decoded")) {
##         pred = cpred = attr(data, "predictors")
##         trms = attr(data, "terms")
##         data = decode.data(as.coded.data(data, formulas = cod))
##         for (form in cod) {
##             vn = all.vars(form)
##             if (!is.na(idx <- grep(vn[1], pred))) {
##                 pred[idx] = vn[2]
##                 cpred = setdiff(cpred, vn[1])
##             }
##         }
##         attr(data, "predictors") = pred
##         new.trms = update(trms, reformulate(c("1", cpred)))   # excludes coded variables
##         attr(new.trms, "orig") = trms       # save orig terms as an attribute
##         attr(data, "terms") = new.trms
##     }
##     data
## }


###################################################
### code chunk number 17: extending.rnw:241-265 (eval = FALSE)
###################################################
## lsm.basis.rsm = function(object, trms, xlev, grid, 
##                          mode = c("asis", "coded", "decoded"), ...) {
##     mode = match.arg(mode)
##     cod = codings(object)
##     if(!is.null(cod) && mode == "decoded") {
##         grid = coded.data(grid, formulas = cod)
##         trms = attr(trms, "orig")   # get back the original terms we saved
##     }
##     
##     m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
##     X = model.matrix(trms, m, contrasts.arg = object$contrasts)
##     bhat = as.numeric(object$coefficients) 
##     V = lsmeans::.my.vcov(object, ...)
##     
##     if (sum(is.na(bhat)) > 0)
##         nbasis = estimability::nonest.basis(object$qr)
##     else
##         nbasis = estimability::all.estble
##     dfargs = list(df = object$df.residual)
##     dffun = function(k, dfargs) dfargs$df
## 
##     list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
##          dffun = dffun, dfargs = dfargs, misc = list())
## }


###################################################
### code chunk number 18: extending.rnw:274-280 (eval = FALSE)
###################################################
## if (requireNamespace("lsmeans", quietly = TRUE)) {
##     importFrom("lsmeans", "recover.data", "lsm.basis")
##     importFrom("estimability", "all.estble", "nonest.basis")
##     S3method(recover.data, rsm)
##     S3method(lsm.basis, rsm)
## }


###################################################
### code chunk number 19: extending.rnw:290-292
###################################################
library("rsm")
example("rsm")   ### (output is not shown) ###


###################################################
### code chunk number 20: extending.rnw:295-297
###################################################
lsmeans(CR.rs2, ~ x1 * x2, mode = "coded", 
        at = list(x1 = c(-1, 0, 1), x2 = c(-2, 2)))


###################################################
### code chunk number 21: extending.rnw:300-301
###################################################
codings(CR.rs1)


###################################################
### code chunk number 22: extending.rnw:304-306
###################################################
lsmeans(CR.rs2, ~ Time * Temp, mode = "decoded", 
        at = list(Time = c(80, 85, 90), Temp = c(165, 185)))


