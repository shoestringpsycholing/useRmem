library(lme4)
library(reshape)
library(plyr)
library(ggplot2)
library(arm)
std.err <- function(x) sd(x)/sqrt(length(x))
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

library(multcomp)

my.contrasts <- rbind(
                  "ILR2 vs. grand mean for unmatched" = c(0, 0, 1/4, 1/4, 1/4, 0, 0, 0),
                  "ILR2+ vs. grand mean for unmatched" = c(0, 0, -3/4, 1/4, 1/4, 0, 0, 0),
                  "ILR3 vs. grand mean for unmatched" = c(0, 0, 1/4, -3/4, 1/4, 0, 0, 0),
                  "grand mean for unmatched" = c(1, 0, 1/4, 1/4, 1/4, 0, 0, 0))

summary(glht(ldt.lmer3, my.contrasts))

#########################
# Soup to nuts with logistic data

# 0. Read data in
imp <- read.csv("data/implicit_learning_mixed_effects_all_3Nov12.csv")
head(imp)
summary(imp)

# some factor coding for reading ease
imp$group <- NA
imp$group[imp$GROUP == 1] <- "incidental"
imp$group[imp$GROUP == 2] <- "intentional"
imp$group[imp$GROUP == 3] <- "control"
imp$group <- as.factor(imp$group)
summary(imp$group)
contrasts(imp$group)

imp$rule <- NA
imp$rule[imp$RULE == 1] <- "noun-det"
imp$rule[imp$RULE == 2] <- "adj-det-noun"
imp$rule[imp$RULE == 3] <- "adj-det-adj-noun"
imp$rule <- as.factor(imp$rule)
summary(imp$rule)
contrasts(imp$rule)
# re-order
imp$rule <- reorder(imp$rule, imp$RULE)
summary(imp$rule)
contrasts(imp$rule)
## # alternative, using relevel()
## imp$rule <- relevel(imp$rule, "adj-det-noun")
## imp$rule <- relevel(imp$rule, "noun-det")

imp$gramm <- NA
imp$gramm[imp$GRAMM == 0] <- "ungrammatical"
imp$gramm[imp$GRAMM == 1] <- "grammatical"
imp$gramm <- as.factor(imp$gramm)
summary(imp$gramm)
contrasts(imp$gramm)

# since we want to analyze ACC, explicitly subset, just for clarity and consistency
imp.all <- imp
imp <- imp[!is.na(imp$ACC), ]
nrow(imp.all)
nrow(imp)
summary(imp)

# 1. Plot the data!
imp.means <- ddply(imp, c("group", "rule", "gramm"), summarize, mean.ACC = mean(ACC), N = length(ACC))
imp.means$se <- std.err.perc(imp.means$mean.ACC, imp.means$N)

# proportions, so barplot meaningful
meansplot <- ggplot(imp.means, aes(rule, mean.ACC)) + geom_bar(aes(fill = group), stat = "identity", position = "dodge") + geom_errorbar(aes(ymin = mean.ACC - 1.96*se, ymax = mean.ACC + 1.96*se, group = group), width = 0.3, position = position_dodge(width = 0.9)) + facet_wrap(~ gramm)

# 2. Think about structure & hypotheses

# hmmmmm.....

# 3. Fit models with simple random effects, increasing fixed effect complexity
(imp.fit1 <- glmer(ACC ~ group + (1|SUBJECT) + (1|ITEM), data = imp, family = "binomial"))
(imp.fit2 <- glmer(ACC ~ group + rule + (1|SUBJECT) + (1|ITEM), data = imp, family = "binomial"))
(imp.fit3 <- glmer(ACC ~ group * rule + (1|SUBJECT) + (1|ITEM), data = imp, family = "binomial"))

imp.means.fit3 <- ddply(imp, c("group", "rule"), summarize, mean.ACC = mean(ACC), N = length(ACC))
imp.means.fit3$se <- std.err.perc(imp.means.fit3$mean.ACC, imp.means.fit3$N)
meansplot.collapsegramm <- ggplot(imp.means.fit3, aes(rule, mean.ACC)) + geom_bar(aes(fill = group), stat = "identity", position = "dodge") + geom_errorbar(aes(ymin = mean.ACC - 1.96*se, ymax = mean.ACC + 1.96*se, group = group), width = 0.3, position = position_dodge(width = 0.9))
print(meansplot.collapsegramm)

# contrasts!
fit3.contrasts <- rbind(
                    "control v. incidental in noun-det" = c(0, 1, 0, 0, 0, 0, 0, 0, 0),
                    "control v. incidental in adj-det-noun" = c(0, 1, 0, 0, 0, 1, 0, 0, 0),
                    "control v. incidental in adj-det-adj-noun" = c(0, 1, 0, 0, 0, 0, 0, 1, 0),
                    "control v. intentional in noun-det" = c(0, 0, 1, 0, 0, 0, 0, 0, 0),
                    "control v. intentional in adj-det-noun" = c(0, 0, 1, 0, 0, 0, 1, 0, 0),
                    "control v. intentional in adj-det-adj-noun" = c(0, 0, 1, 0, 0, 0, 0, 0, 1),
                    "incidental v. intentional in noun-det" = c(0, 1, -1, 0, 0, 0, 0, 0, 0),
                    "incidental v. intentional in adj-det-noun" = c(0, 1, -1, 0, 0, 1, -1, 0, 0),
                    "incidental v. intentional in adj-det-adj-noun" = c(0, 1, -1, 0, 0, 0, 0, 1, -1))
summary(glht(imp.fit3, fit3.contrasts))

pval.z <- function(z) 2*(1-pnorm(abs(z)))
pval.z(-2.893)
                    
# relevel to check
imp$group.relevel <- relevel(imp$group, "intentional")
imp$rule.relevel <- relevel(imp$rule, "adj-det-noun")

(imp.fit3.relevel <- glmer(ACC ~ group.relevel * rule.relevel + (1|SUBJECT) + (1|ITEM), data = imp, family = "binomial"))

print(meansplot)
(imp.fit4 <- glmer(ACC ~ group * rule * gramm + (1|SUBJECT) + (1|ITEM), data = imp, family = "binomial"))

imp.fit5 <- glmer(ACC ~ group * rule + group * gramm + rule * gramm + (1|SUBJECT) + (1|ITEM), data = imp, family = "binomial")
anova(imp.fit4, imp.fit5)

## # code to bootstrap model comparison, stolen from Faraway (2006) "Extending the Linear Model with R"
## n.sims <- 1000 # this could take a while!
## lrstat <- numeric(n.sims)
## for(i in 1:n.sims) {
##   y <- unlist(simulate(imp.fit4))
##   bsmaller <- glmer(y ~ group * rule + group * gramm + rule * gramm + (1|SUBJECT) + (1|ITEM), data = imp, family = "binomial")
##   balt <- glmer(y ~ group * rule * gramm + (1|SUBJECT) + (1|ITEM), data = imp, family = "binomial")
##   lrstat[i] <- as.numeric(2*(logLik(balt)-logLik(bsmaller)))
## }

## mean(lrstat > 8.3816)

# 4. Add random effects
#(imp.fit6 <- glmer(ACC ~ group * rule * gramm + (1 + rule * gramm|SUBJECT) + (1 + group|ITEM), data = imp, family = "binomial"))

# code to save and recall model later, for time-saving
#save(imp.fit6, file = "data/imp.fit6.RData")
load("data/imp.fit6.RData")

(imp.fit7 <- glmer(ACC ~ group * rule * gramm + (1 + rule|SUBJECT) + (1|ITEM), data = imp, family = "binomial"))
anova(imp.fit4, imp.fit7)
print(imp.fit4, cor = FALSE)
print(imp.fit7, cor = FALSE)

(imp.fit8 <- glmer(ACC ~ group * rule * gramm + (1 + rule + gramm|SUBJECT) + (1|ITEM), data = imp, family = "binomial"))
anova(imp.fit7, imp.fit8)
print(imp.fit7, cor = FALSE)
print(imp.fit8, cor = FALSE)

#(imp.fit9 <- glmer(ACC ~ group * rule * gramm + (1 + rule + gramm|SUBJECT) + (1+group|ITEM), data = imp, family = "binomial"))
save(imp.fit9, file = "data/imp.fit9.RData")
load("data/imp.fit9.RData")
anova(imp.fit8, imp.fit9)
print(imp.fit8, cor = FALSE)

anova(imp.fit4, imp.fit9)
print(imp.fit4, cor = FALSE)
print(imp.fit9, cor = FALSE)
print(meansplot)

anova(imp.fit6, imp.fit9)
#imp.fit11 <- glmer(ACC ~ group * rule + group * gramm + rule * gramm + (1 + rule * gramm|SUBJECT) + (1 + group|ITEM), data = imp, family = "binomial")
#save(imp.fit11, file = "data/imp.fit11.RData")
load("data/imp.fit11.RData") # false convergence error message!
anova(imp.fit6, imp.fit11)
print(imp.fit11, cor = FALSE)

#imp.fit10 <- glmer(ACC ~ group * rule + rule * gramm + group * gramm + (1 + rule + gramm|SUBJECT) + (1+group|ITEM), data = imp, family = "binomial")
#save(imp.fit10, file = "data/imp.fit10.RData")
load("data/imp.fit10.RData")

# 5. Center and transform variables as necessary
imp$cGramm <- imp$gramm
contrasts(imp$cGramm)
contrasts(imp$cGramm) <- c(-1, 1) # why 1's?  Gives interpretation for interactions of "diff from grand mean"
contrasts(imp$cGramm)
imp.fit4b <- glmer(ACC ~ group * rule * cGramm + (1 |SUBJECT) + (1 |ITEM), data = imp, family = "binomial")
print(imp.fit4, cor = FALSE)
print(imp.fit4b, cor = FALSE)

# 6. Try to interpret, try additional models as needed

# ???

# 7. Try "influence.ME" package
# use simple ranef model, for speed reasons
# if wanting to use for inferential purposes, use with full "final" model -- may take a long time!
library(influence.ME)
#influence.fit4.bysubj <- influence(imp.fit4, "SUBJECT")
#save(influence.fit4.bysubj, file = "data/influence.fit4.bysubj.RData")
load("data/influence.fit4.bysubj.RData")


cooks.fit4.bysubj <- as.data.frame(cooks.distance(influence.fit4.bysubj))
cooks.fit4.bysubj$Subject <- as.factor(rownames(cooks.fit4.bysubj))
cooks.fit4.bysubj$Subject <- reorder(cooks.fit4.bysubj$Subject, cooks.fit4.bysubj$V1)

ggplot(cooks.fit4.bysubj, aes(V1, Subject)) +
  geom_text(aes(label = Subject)) +
  geom_vline(xintercept = 4/nrow(cooks.fit4.bysubj), color = "red", linetype = 2)

imp.fit4.rem2 <- glmer(ACC ~ group * rule * gramm + (1|SUBJECT) + (1|ITEM), data = imp[!imp$SUBJECT %in% c(47, 59), ], family = "binomial")
print(imp.fit4.rem2, cor = FALSE)

#influence.fit4.rem2.bysubj <- influence(imp.fit4.rem2, "SUBJECT")
#save(influence.fit4.rem2.bysubj, file = "data/influence.fit4.rem2.bysubj.RData")
load("data/influence.fit4.rem2.bysubj.RData")


cooks.fit4.rem2.bysubj <- as.data.frame(cooks.distance(influence.fit4.rem2.bysubj))
cooks.fit4.rem2.bysubj$Subject <- as.factor(rownames(cooks.fit4.rem2.bysubj))
cooks.fit4.rem2.bysubj$Subject <- reorder(cooks.fit4.rem2.bysubj$Subject, cooks.fit4.rem2.bysubj$V1)

ggplot(cooks.fit4.rem2.bysubj, aes(V1, Subject)) +
  geom_text(aes(label = Subject)) +
  geom_vline(xintercept = 4/nrow(cooks.fit4.rem2.bysubj), color = "red", linetype = 2)


# 8. Final reporting!


##############################
# Ordinal data

pc1 <- read.csv("data/pc1.data.processed")
  
ggplot(pc1, aes(factor(response))) + geom_bar() + facet_grid(ptype ~ pred)
 
ggplot(pc1[pc1$ptype != "filler", ], aes(factor(response))) + geom_bar() + facet_grid(ptype ~ pred)

library(ordinal)
pc1.clmm <- clmm(factor(response) ~ ptype * pred + (1|subj) + (1|item), link = "probit", threshold = "equidistant", data = pc1)
  
pc1.clmm2 <- clmm(factor(response) ~ ptype * pred + (1|subj) + (1|item), link = "probit", threshold = "flexible", data = pc1)
  
pc1.clmm3 <- clmm(factor(response) ~ ptype * pred + (1|subj) + (1|item), link = "logit", threshold = "flexible", data = pc1)

