## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(fig.path='figdir/fig', debug=TRUE, echo=TRUE)

## ----results='asis'------------------------------------------------------
library(xtable)
options(xtable.floating = FALSE)

## ----results='asis'------------------------------------------------------
data(tli)
xtable(tli[1:10, ])

## ----results='asis'------------------------------------------------------
design.matrix <- model.matrix(~ sex*grade, data = tli[1:10, ])
xtable(design.matrix, digits = 0)

## ----results='asis'------------------------------------------------------
fm1 <- aov(tlimth ~ sex + ethnicty + grade + disadvg, data = tli)
xtable(fm1)

## ----results='asis'------------------------------------------------------
fm2 <- lm(tlimth ~ sex*ethnicty, data = tli)
xtable(fm2)

## ----results='asis'------------------------------------------------------
xtable(anova(fm2))

## ----results='asis'------------------------------------------------------
fm2b <- lm(tlimth ~ ethnicty, data = tli)
xtable(anova(fm2b, fm2))

## ----results='asis'------------------------------------------------------
fm3 <- glm(disadvg ~ ethnicty*grade, data = tli, family = binomial)
xtable(fm3)

## ----results='asis'------------------------------------------------------
xtable(anova(fm3))

## ----results='asis'------------------------------------------------------
pr1 <- prcomp(USArrests)
xtable(pr1)

## ----results='asis'------------------------------------------------------
xtable(summary(pr1))

## ----include=FALSE-------------------------------------------------------
# pr2 <- princomp(USArrests)
# xtable(pr2)

## ----results='asis'------------------------------------------------------
temp.ts <- ts(cumsum(1 + round(rnorm(100), 0)),
              start = c(1954, 7), frequency = 12)
temp.table <- xtable(temp.ts, digits = 0)
temp.table

## ----include=FALSE-------------------------------------------------------
# ## Demonstrate saving to file
# for(i in c("latex", "html")) {
#   outFileName <- paste("xtable.", ifelse(i=="latex", "tex", i), sep = "")
#   print(xtable(lm.D9), type = i, file = outFileName, append = TRUE,
#         latex.environments = NULL)
#   print(xtable(lm.D9), type = i, file = outFileName, append = TRUE,
#         latex.environments = "")
#   print(xtable(lm.D9), type = i, file = outFileName, append = TRUE,
#         latex.environments = "center")
#   print(xtable(anova(glm.D93, test = "Chisq")),
#         type = i, file = outFileName,
#         append = TRUE)
#   print(xtable(anova(glm.D93)), hline.after = c(1),
#         size = "small", type = i,
#         file = outFileName, append = TRUE)
#   # print(xtable(pr2), type = i, file = outFileName, append = TRUE)
# }

## ----results='asis'------------------------------------------------------
dat <- mtcars[1:3, 1:6]
x <- xtable(dat)
x

## ----results='asis'------------------------------------------------------
align(x) <- xalign(x)
digits(x) <- xdigits(x)
display(x) <- xdisplay(x)
x

## ----results='asis'------------------------------------------------------
xtable(dat, auto = TRUE)

## ----results='asis'------------------------------------------------------
x <- xtable(dat)
autoformat(x)

## ----results='asis'------------------------------------------------------
insane <- data.frame(Name = c("Ampersand","Greater than","Less than",
                            "Underscore","Per cent","Dollar",
                            "Backslash","Hash","Caret","Tilde",
                            "Left brace","Right brace"),
                     Character = I(c("&",">","<","_","%","$",
                                     "\\","#","^","~","{","}")))
colnames(insane)[2] <- paste(insane[, 2], collapse = "")
xtable(insane)

## ----results='asis'------------------------------------------------------
wanttex <- xtable(data.frame(Column =
                             paste("Value_is $10^{-",1:3,"}$", sep = "")))
print(wanttex, sanitize.text.function =
      function(str) gsub("_", "\\_", str, fixed = TRUE))

## ----results='asis'------------------------------------------------------
mat <- round(matrix(c(0.9, 0.89, 200, 0.045, 2.0), c(1, 5)), 4)
rownames(mat) <- "$y_{t-1}$"
colnames(mat) <- c("$R^2$", "$\\bar{x}$", "F-stat", "S.E.E", "DW")
mat <- xtable(mat)
print(mat, sanitize.text.function = function(x) {x})

## ----results='asis'------------------------------------------------------
money <- matrix(c("$1,000","$900","$100"), ncol = 3,
                dimnames = list("$\\alpha$",
                                c("Income (US$)","Expenses (US$)",
                                  "Profit (US$)")))
print(xtable(money), sanitize.rownames.function = function(x) {x})

## ----results='asis'------------------------------------------------------
print(xtable(anova(fm3), caption = "\\tt latex.environments = \"\""),
      floating = TRUE, latex.environments = "")
print(xtable(anova(fm3), caption = "\\tt latex.environments = \"center\""),
      floating = TRUE, latex.environments = "center")

## ----results='asis'------------------------------------------------------
tli.table <- xtable(tli[1:10, ])
align(tli.table) <- rep("r", 6)
tli.table

## ----results='asis'------------------------------------------------------
align(tli.table) <- "|rrl|l|lr|"
tli.table

## ----results='asis'------------------------------------------------------
align(tli.table) <- "|rr|lp{3cm}l|r|"
tli.table

## ----results='asis'------------------------------------------------------
display(tli.table)[c(2,6)] <- "f"
digits(tli.table) <- 3
tli.table

## ----results='asis'------------------------------------------------------
digits(tli.table) <- 1:(ncol(tli)+1)
tli.table

## ----results='asis'------------------------------------------------------
digits(tli.table) <- matrix(0:4, nrow = 10, ncol = ncol(tli)+1)
tli.table

## ----results='asis'------------------------------------------------------
tli.table <- xtable(tli[1:10, ])
print(tli.table, include.rownames = FALSE)

## ----results='asis'------------------------------------------------------
align(tli.table) <- "|r|r|lp{3cm}l|r|"
print(tli.table, include.rownames = FALSE)

## ------------------------------------------------------------------------
align(tli.table) <- "|rr|lp{3cm}l|r|"

## ----results='asis'------------------------------------------------------
print(tli.table, include.colnames = FALSE)

## ----results='asis'------------------------------------------------------
print(tli.table, include.colnames = FALSE,
      hline.after = c(0,nrow(tli.table)))

## ----results='asis'------------------------------------------------------
print(tli.table, include.colnames = FALSE, include.rownames = FALSE)

## ----results='asis'------------------------------------------------------
print(tli.table, rotate.rownames = TRUE, rotate.colnames = TRUE)

## ----results='asis'------------------------------------------------------
print(xtable(anova(fm3)), hline.after = c(1))

## ----results='asis'------------------------------------------------------
tli.table <- xtable(tli[1:10, ])
print(tli.table, include.rownames = FALSE, booktabs = TRUE)

## ----results='asis'------------------------------------------------------
bktbs <- xtable(matrix(1:10, ncol = 2))
hlines <- c(-1, 0, 1, nrow(bktbs))
print(bktbs, booktabs = TRUE, hline.after = hlines)

## ----results='asis'------------------------------------------------------
print(xtable(anova(fm3)), size = "large")

## ----results='asis'------------------------------------------------------
print(xtable(anova(fm3)), size = "\\setlength{\\tabcolsep}{12pt}")

## ----results='asis'------------------------------------------------------
x <- matrix(rnorm(1000), ncol = 10)
x.big <- xtable(x, caption = "A \\code{longtable} spanning several pages")
print(x.big, hline.after=c(-1, 0), tabular.environment = "longtable")

## ------------------------------------------------------------------------
Grade3 <- c("A","B","B","A","B","C","C","D","A","B",
            "C","C","C","D","B","B","D","C","C","D")
Grade6 <- c("A","A","A","B","B","B","B","B","C","C",
            "A","C","C","C","D","D","D","D","D","D")
Cohort <- table(Grade3, Grade6)
Cohort

## ----results='asis'------------------------------------------------------
xtable(Cohort)

## ----results='asis'------------------------------------------------------
addtorow <- list()
addtorow$pos <- list(0, 0)
addtorow$command <- c("& \\multicolumn{4}{c}{Grade 6} \\\\\n",
                      "Grade 3 & A & B & C & D \\\\\n")
print(xtable(Cohort), add.to.row = addtorow, include.colnames = FALSE)

## ----results='asis'------------------------------------------------------
x <- x[1:30, ]
x.side <- xtable(x, caption = "A sideways table")
print(x.side, floating = TRUE, floating.environment = "sidewaystable")

## ----results='asis'------------------------------------------------------
x <- x[1:20, ]
x.rescale <- xtable(x)
print(x.rescale, scalebox = 0.7)

## ----results='asis'------------------------------------------------------
df <- data.frame(name = c("A","B"), right = c(1.4, 34.6),
                 left = c(1.4, 34.6), text = c("txt1","txt2"))
print(xtable(df, align = c("l", "|c", "|R{3cm}", "|L{3cm}", "| p{3cm}|")),
      floating = FALSE, include.rownames = FALSE)

## ----results='asis'------------------------------------------------------
df.width <- data.frame(One = c("item 1", "A"), Two = c("item 2", "B"),
                       Three = c("item 3", "C"), Four = c("item 4", "D"))
x.width <- xtable(df.width)
align(x.width) <- "|l|X|l|l|l|"
print(x.width, tabular.environment = "tabularx", width = "\\textwidth")

## ------------------------------------------------------------------------
x.out <- print(tli.table, print.results = FALSE)

## ------------------------------------------------------------------------
x.ltx <- toLatex(tli.table)
class(x.ltx)
x.ltx

## ----results='asis'------------------------------------------------------
toLatex(sessionInfo())

