# R script: Causual inference 1: Matching
# Author: Yanhan Si
# Nov 11, 2020
# 79: -------------------------------------------------------------------------
# install packages
# install.packages("tableone")
# install.packages("Matching")
# library: --------------------------------------------------------------------
library(tableone)
library(Matching)
library(tidyverse)
library(readr)
# data: -----------------------------------------------------------------------
path = "./data"
file_name = sprintf('%s/student-mat.csv', path)
math = read_delim(file_name, ";", escape_double = FALSE, trim_ws = TRUE)
mydata = math %>% 
  select(-G1, -G2) %>% 
  mutate(
         schoolsup = if_else(schoolsup == "yes", 1, 0, missing = 0),
         famsup = if_else(famsup == "yes", 1, 0, missing = 0),
         paid = if_else(paid == "yes", 1, 0, missing = 0),
         activities = if_else(activities == "yes", 1, 0, missing = 0),
         nursery = if_else(nursery == "yes", 1, 0, missing = 0),
         higher = if_else(higher == "yes", 1, 0, missing = 0),
         internet = if_else(internet == "yes", 1, 0, missing = 0),
         romantic = if_else(romantic == "yes", 1, 0, missing = 0),
         school = if_else(school == "GP", 1, 0, missing = 0),
         sex = if_else(sex == "F", 1, 0, missing = 0),
         address = if_else(address == "U", 1, 0, missing = 0),
         famsize = if_else(famsize == "LE3", 1, 0, missing = 0),
         Pstatus = if_else(Pstatus == "T", 1, 0, missing = 0),
         Mjob = if_else(Mjob == "at_home", 1, 0, missing = 0),
         Fjob = if_else(Fjob == "at_home", 1, 0, missing = 0),
         #reason_home = if_else(reason == "home", 1, 0, missing = 0),
         #reason_course = if_else(reason == "course", 1, 0, missing = 0),
         #reason_reputation = if_else(reason == "reputation", 1, 0, missing = 0),
         guardian = if_else(guardian == "other", 0, 1, missing = 0),
         treatment = if_else(absences <3, 1, 0),
         grade = G3) %>%
  select(-reason, -absences, -G3) %>% 
  filter(grade > 0)

nrow(mydata)

t.test(formula = age ~ treatment ,data= mydata)

#hist(mydata$grade)
#summary(mydata$grade)
#table(mydata$treatment)
#write_csv(mydata, "cleaned_math.csv")

summary(mydata)
table(mydata$treatment)
nrow(mydata)
names(mydata)[1:28]

xvars<-names(mydata)[c(1:28)]

mydata %>% filter(treatment == 1 ) %>% select(grade) %>% summary()
mydata %>% filter(treatment == 0 ) %>% select(grade) %>% summary()

table1<- CreateTableOne(vars=xvars,strata="treatment", data=mydata, test=TRUE)
print(table1,smd=TRUE)


greedymatch<-Match(Y=mydata$grade, Tr=mydata$treatment,M=1,X=mydata[xvars], replace=FALSE)
matched<-mydata[unlist(greedymatch[c("index.treated","index.control")]), ]

#get table 1 for matched data with standardized differences
matchedtab1<-CreateTableOne(vars=xvars, strata ="treatment", 
                            data=matched, test = TRUE)
print(matchedtab1, smd = TRUE)

#outcome analysis
y_trt1<-matched$grade[matched$treatment==1]
y_con1<-matched$grade[matched$treatment==0]

#pairwise difference
diffy1<-y_trt1-y_con1

#paired t-test
t.test(diffy1)

##########################
#propensity score matching
#########################

small_data = mydata %>% select(-grade)

psmodel<-glm(treatment ~ ., family=binomial(link ="logit"),data=small_data)

#show coefficients etc
summary(psmodel)
#create propensity score
pscore<-psmodel$fitted.values


#do greedy matching on logit(PS) using Match with a caliper

logit <- function(p) {log(p)-log(1-p)}
psmatch<-Match(Y=mydata$grade,Tr=mydata$treatment,M=1,X=logit(pscore),replace=FALSE,caliper=.2)
matched<-mydata[unlist(psmatch[c("index.treated","index.control")]), ]


#get standardized differences
matchedtab2<-CreateTableOne(vars=xvars, strata ="treatment", 
                            data=matched, test = TRUE)
print(matchedtab2, smd = TRUE)

#outcome analysis
y_trt2<-matched$grade[matched$treatment==1]
y_con2<-matched$grade[matched$treatment==0]

#pairwise difference
diffy2<-y_trt2-y_con2

#paired t-test
t.test(diffy2)

# -----------------------------------------------------------------------------------


#create weights
weight<-ifelse(mydata$treatment==1,1/(pscore),1/(1-pscore))
#hist(weight)
library(survey)
#design <- svydesign(ids=~1, weights=~boosted, data=lalonde)
weighteddata<-svydesign(ids = ~ 1, weights = ~ weight, data = mydata)

glm1 <- svyglm(grade ~ treatment, design=weighteddata)
summary(glm1)
#apply weights to data


#weighted table 1
weightedtable <-svyCreateTableOne(vars = xvars, strata = "treatment", 
                                  data = weighteddata, test = FALSE)
## Show table with SMD
print(weightedtable, smd = TRUE)

#get causal risk difference
lm.obj<-lm(grade~treatment,weights=weight, data = mydata)

summary(lm.obj)
betaiptw<-coef(lm.obj)
SE<-sqrt(diag(vcovHC(glm.obj, type="HC0")))

causalrd<-(betaiptw[2])
lcl<-(betaiptw[2]-1.96*SE[2])
ucl<-(betaiptw[2]+1.96*SE[2])
c(lcl,causalrd,ucl)

#get causal relative risk. Weighted GLM
glm.obj<-glm(died~treatment,weights=weight,family=quasibinomial(link=log))
#summary(glm.obj)
betaiptw<-coef(glm.obj)
#to properly account for weighting, use asymptotic (sandwich) variance
SE<-sqrt(diag(vcovHC(glm.obj, type="HC0")))

#get point estimate and CI for relative risk (need to exponentiate)
causalrr<-exp(betaiptw[2])
lcl<-exp(betaiptw[2]-1.96*SE[2])
ucl<-exp(betaiptw[2]+1.96*SE[2])
c(lcl,causalrr,ucl)



