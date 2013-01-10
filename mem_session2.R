library(lme4)
library(reshape)
library(plyr)
library(ggplot2)
library(arm)
std.err.perc <- function(p, n) sqrt((p*(1-p))/n)

#########
# More coefficient fun
# centering, transformations, contrasts etc.

ldt <- read.csv("data/russianLDTclean.csv")
ldt <- ldt[ldt$PrimeType == "semi-regular", ]
ldt$Condition <- relevel(ldt$Condition, "unmatched prime")
ldt$Target.Freq <- relevel(ldt$Target.Freq, "low frequency")

ldt.lm1 <- lm(RT ~ Condition, data = ldt)
(ldt.lmer1 <- lmer(RT ~ Condition + (1|Subject) + (1|Item), data = ldt, REML = TRUE))
(ldt.lmer2 <- lmer(RT ~ Condition + ILR + (1|Subject) + (1|Item), data = ldt, REML = TRUE))
(ldt.lmer3 <- lmer(RT ~ Condition + ILR + Condition:ILR + (1|Subject) + (1|Item), data = ldt, REML = TRUE))

ldt$ILR.sumc <- ldt$ILR
contrasts(ldt$ILR)
contrasts(ldt$ILR.sumc) <- contr.sum(4)
contrasts(ldt$ILR.sumc)

(ldt.lmer3b <- lmer(RT ~ Condition + ILR.sumc + Condition:ILR.sumc + (1|Subject) + (1|Item), data = ldt, REML = TRUE))

(ldt.lmer4 <- lmer(RT ~ Condition * ILR + Target.Freq  + (1|Subject) + (1|Item), data = ldt, REML = TRUE))
(ldt.lmer5 <- lmer(RT ~ Condition * ILR * Target.Freq  + (1|Subject) + (1|Item), data = ldt[ldt$PrimeType == "semi-regular", ], REML = TRUE))

std.err <- function(x) sd(x)/sqrt(length(x))
ldt.means <- ddply(ldt, c("Condition", "ILR", "Target.Freq"), summarize, meanRT = mean(RT), SE <- std.err(RT))
colnames(ldt.means)[5] <- "SE"
ggplot(ldt.means, aes(ILR, meanRT)) + geom_pointrange(aes(ymax = meanRT + 1.96*SE, ymin = meanRT - 1.96*SE, linetype = Condition), position = position_dodge(width = .5)) + facet_wrap(~ Target.Freq)

(ldt.lmer6 <- lmer(RT ~ Condition * ILR * logFreq + (1|Subject) + (1|Item), data = ldt, REML = TRUE))
(ldt.lmer7 <- lmer(RT ~ Condition * ILR * scale(logFreq) + (1|Subject) + (1|Item), data = ldt, REML = TRUE))

library(multcomp)


######
# hypothesis testing





#######
# model building


##################
### UNDER CONSTRUCTION
##################
ch <- read.csv("data/LASRS Data for WS 2013-01-09.csv")

ch$FW.all <- ch$FW
ch$FW.all[is.na(ch$FW.all)] <- 0

ch.melt <- as.data.frame(melt(ch, measure.vars = c("Target", "FW", "FW.all")))

ch.means <- ddply(ch.melt, c("xDistalRate", "xCondition", "variable"), summarize, mean = mean(value, na.rm = TRUE), n = length(value))
  
ch.means$std.err.perc <- std.err.perc(ch.means$mean, ch.means$n)
  
ggplot(ch.means, aes(xDistalRate, mean)) + geom_pointrange(aes(ymin = mean - 1.96*std.err.perc, ymax = mean + 1.96*std.err.perc, color = factor(xCondition)), position = position_dodge(width = .3)) + scale_color_brewer(palette = "Set1") + facet_wrap(~ variable)
  
ch.subjmeans <- ddply(ch.melt, c("Subject", "xDistalRate", "xCondition", "variable"), summarize, mean = mean(value, na.rm = TRUE), n = length(value))
ggplot(ch.subjmeans, aes(xDistalRate, mean)) + geom_point() + geom_smooth() + geom_smooth(method = "lm")+ facet_wrap(~ variable)
  
ch.glmer1 <- glmer(FW.all ~ xDistalRate + (1+xDistalRate|Subject) + (1+xDistalRate|xBasename), data = ch, family = "binomial")
  
ch$fCondition <- as.factor(ch$xCondition)
ch$fSubject <- as.factor(ch$Subject)
ch.glmer1b <- glmer(FW.all ~ xDistalRate + (1+xDistalRate|fSubject:fCondition) + (1+xDistalRate|xBasename), data = ch, family = "binomial")
  
ch.glmer1c <- glmer(FW.all ~ xDistalRate + (1+xDistalRate|fSubject/fCondition) + (1+xDistalRate|xBasename), data = ch, family = "binomial")
  
ch.glmer2 <- glmer(FW.all ~ factor(xDistalRate) + (1+xDistalRate|Subject) + (1+xDistalRate|xBasename), data = ch, family = "binomial")
  
ch.glmer3 <- glmer(FW.all ~ xDistalRate * factor(xCondition) + (1+xDistalRate|Subject) + (1+xDistalRate|xBasename), data = ch, family = "binomial", verbose = TRUE)

##############################
pc1 <- read.csv("pc1.data.processed")
  
ggplot(pc1, aes(factor(response))) + geom_bar() + facet_grid(ptype ~ pred)
 
ggplot(pc1[pc1$ptype != "filler", ], aes(factor(response))) + geom_bar() + facet_grid(ptype ~ pred)
  
pc1.clmm <- clmm(factor(response) ~ ptype * pred + (1|subj) + (1|item), link = "probit", threshold = "equidistant", data = pc1)
  
pc1.clmm2 <- clmm(factor(response) ~ ptype * pred + (1|subj) + (1|item), link = "probit", threshold = "flexible", data = pc1)
  
pc1.clmm3 <- clmm(factor(response) ~ ptype * pred + (1|subj) + (1|item), link = "logit", threshold = "flexible", data = pc1)










##############################  
ldt.lm <- lm(RT ~ Condition * ILR, data = ldt)
summary(ldt.lm)
anova(ldt.lm)
  
ldt.bysubj <- cast(ldt, Subject + Condition + ILR ~ ., value = "RT", c(mean, sd))
  
bysubj.aov <- aov(mean ~ Condition * ILR + Error(Subject/Condition), data = ldt.bysubj)
  
bysubj.lmer <- lmer(mean ~ Condition * ILR + (1 + Condition|Subject), data = ldt.bysubj)
#ldt.lmer <- lmer(RT ~ Condition * ILR + (1+Condition|Subject) + (1+ILR|Item), data = ldt)
#print(ldt.lmer, cor = FALSE)
  
summary(bysubj.aov)
anova(bysubj.lmer)
  
  x <- 1:10
  err <- rnorm(10, sd = 2)
  effect <- 2
  y <- effect*x + err
  plot(x, y, ylim = c(0, max(y)))
  xy.lm <- lm(y ~ x)
  summary(xy.lm)
  abline(xy.lm)
  resid(xy.lm)
  segments(resid(xy.lm), 1.75*x)
  

