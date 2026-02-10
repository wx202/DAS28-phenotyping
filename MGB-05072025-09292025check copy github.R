rm(list=ls())
options(warn = -1)
library(tidyr)
library(dplyr)
library(readxl)
library(readr)
library(pROC)
library(parsnip)
library(glmnet)
library(lubridate)
library('SuperLearner')
library(brulee)
library('tidyverse')
# setwd("/Users/xuanwang/Library/CloudStorage/Box-Box/RACAT/MGB")
# list = read.csv('/Users/xuanwang/Library/CloudStorage/Box-Box/RACAT/All_Features_n2875_0505_2025.csv', stringsAsFactor = F)

#### Set constants
days_per_month = 30
radius_features_months = 3 # radius around target month to extract features
date.screen <- as.Date("2008-12-31")
days.per.year <- 365.25

#### feature screening 
{
load(paste0("data/derived-data-2025-05-07.RData"))
 
# d.proximal, one feature file for each of lab, med, diag, proc, nlp
# PatientICN target_month features

# d.vara.cln, file with outcome
# PatientICN Date.Of.Observation Das28
  
# demo
# PatientICN target_month age_target sex   race     hu

## screen features begin here
# fitting a multivariable logistic regression model for moderate-high disease activity at target months 
# against each individual candidate feature while additionally adjusting by age, sex, race, and healthcare utilization, 
# Codified and narrative features with p-values ≤0.1 and prevalence >5% were then selected to be features 
# used for training the disease activity algorithm. 
# We used target months before 2009 for feature screening and reserved the data from 2009 onward for training and validation. 
screenFeatures = function(d.proximal, threshold=0.1, threshold.pre=0.05) {

  feature.screen <- d.proximal %>%
    inner_join(targets %>% select(PatientICN,target_month,screen_month), by=c("PatientICN","target_month")) %>%
    subset(target_month <= screen_month) %>% # pre-09 data for screening
    select(-screen_month) %>%
    replace(is.na(.), 0) 
  
  analytic.screen <- d.vara.cln %>% subset(Date.Of.Observation <= date.screen) %>%
    select(PatientICN,target_date=Date.Of.Observation,Das28) %>%
    left_join(demo, by=c("PatientICN","target_date")) %>%
    select(PatientICN,target_month,sex,race,age_target,hu,Das28) %>%
    inner_join(feature.screen, by=c("PatientICN","target_month"))
    
  # aa=analytic.screen %>% left_join(demo,by=c('PatientICN',"target_date")) 
  # merged=aa %>% inner_join(feature.screen, by=c("PatientICN","target_date")) 
  feat.nme = names(analytic.screen)[!(names(analytic.screen) %in% c("PatientICN","target_month","Das28",
                                                  "age_target","sex","race",'hu' ))]

  p.screen = sapply(feat.nme, function(fn) {
    f.screen = paste0("Das28~sex+race+age_target+hu+","`",fn,"`")
    m.screen = lm(as.formula(f.screen),data=analytic.screen)
    # summary(m.screen)$coefficient[2,4]
    summ.screen <- summary(m.screen)$coefficient
    if(nrow(summ.screen)==6+1) { # linear model cannot be fit likely due to rarity of feature
      rv <- 1
    } else rv <- summ.screen[7+1,"Pr(>|t|)"]
    rv
  })
  p.screen.adjusted = p.adjust(p.screen, method="BH", n=length(feat.nme)) # TODO: check what # features to correct

  if (threshold==1){threshold=quantile(p.screen.adjusted, 0.1)}
  # d.proximal.screened = d.proximal %>%replace(is.na(.), 0) %>%subset(target_date > date.screen) %>%
  #   select(c("PatientICN","target_date",feat.nme[p.screen.adjusted < threshold]))
  d.proximal.screened =d.proximal %>%
    inner_join(targets %>% select(PatientICN,target_month,screen_month), by=c("PatientICN","target_month")) %>%
    subset(target_month> screen_month) %>% # pre-09 data for screening
    select(c("PatientICN","target_month",feat.nme[p.screen< threshold])) %>%
    replace(is.na(.), 0)
  # d.proximal.screened= d.proximal.screened %>%mutate(across(c(3:ncol(.)), ~if_else(.>0,1,0)))
  d.proximal.screened <- d.proximal.screened %>%
    select(c(1:2,(which(apply(d.proximal.screened[,-c(1:2)]>0,2,mean)>threshold.pre)+2)))
  
  # fea.screen=p.screen[names(p.screen)%in%names(d.proximal.screened)]
  # save(fea.screen,file='data/pnlp.featurescreen.rda')
}

}

#### missing data imputation
{
  load(paste0("data/data-screened-2025-05-07.RData"))
  # # dat=dat%>%select(-c(target_month)) 
  # dat$target_month=NULL
  # list.screen=names(dat)[12:(ncol(dat)-1)]
  # list.screen=list[which(list$Code%in%list.screen),]
  # write_csv(list.screen,'list.screen.csv')
  # dat=dat[,c(1:11,which(names(dat)%in%list$Code),ncol(dat))]
  
## dat, after date.screen
# "PatientICN"     "target_date"    "Das28"          "age_target"     "sex"            "race"          
#  "hu"             "crp_median"     "esr_median"     "ccp_binary"     "rf_binary"
  
## impute missing crp and esr using codified features and nlp feature
## for binary ccp and rf, add a missing indicator
  dat=dat%>%mutate(across(c(8:9), ~log(.+1)))
  for (v in c( "crp_median", "esr_median")){
    ind.tr=which(!is.na(dat[,v]))
    ind.imp=which(is.na(dat[,v]))
    cvfit = cv.glmnet(x=as.matrix(dat[ind.tr,12:(ncol(dat)-1)]),
                      y=as.numeric((dat[ind.tr,v])[[1]]), nfolds=10)
    pred=predict(cvfit, newx = as.matrix(dat[ind.imp,12:(ncol(dat)-1)]),s = "lambda.min")
    dat[ind.imp,v]=pred
  }
  tmp=dat[,1:11]
  tmp$ccp_binary_ind=as.numeric(is.na(tmp$ccp_binary))
  tmp$ccp_binary[is.na(tmp$ccp_binary)]=0
  tmp$rf_binary_ind=as.numeric(is.na(tmp$rf_binary))
  tmp$rf_binary[is.na(tmp$rf_binary)]=0
  dat=cbind(tmp,dat[,12:ncol(dat)])
  save(dat,file=paste0('data/dat-analysis-',Sys.Date(),'.rda'))
}

#### Analysis
{
load(paste0("data/dat-analysis-2025-05-07.rda"))
# dat
# "PatientICN"     "target_date"    "Das28"          "age_target"     "sex"            "race"          
# "hu"             "crp_median"     "esr_median"     "ccp_binary"     "rf_binary"      "ccp_binary_ind"
# "rf_binary_ind"

## define codified and nlp features to be binary
dat=dat %>%mutate(across(c(14:(ncol(.)-1)), ~if_else(.>0,1,0)))
d.cln = dat %>% select(-c(PatientICN,target_date,Das28))
d.cln$Das28_LDA=as.factor(d.cln$Das28_LDA)
names(d.cln)[1:20]=gsub(":","",names(d.cln)[1:20])
names(d.cln)[1:20]=gsub("\\.","",names(d.cln)[1:20])
## index for lab, codified features, nlp, try different sets of features
ind.lab=5:10
ind.co=11:28
ind.nlp=29:(ncol(d.cln)-1)
dat.co=d.cln[,-c(ind.lab,ind.nlp)]
dat.nlp=d.cln[,-c(ind.lab,ind.co)]
dat.co.lab= d.cln[,-c(ind.nlp)]
dat.co.nlp=d.cln[,-c(ind.lab)]
dat.all=d.cln
set=0
if (set==1){d.cln=dat.co} 
if (set==2){d.cln=dat.co.lab}
if (set==3){d.cln=dat.co.nlp} 
if (set==4){d.cln=dat.nlp} 
if (set==0){d.cln=dat.all}

## use common features with the other site when do transportation
if (common==1){
load("~/Dropbox (Harvard University)/Xuan/Projects/RACAT/covar.common-2025-05-22.rda")
d.cln=d.cln[,c(covar.common,"Das28_LDA")] }

## write model
covar=names(d.cln)[-ncol(d.cln)]
ind.all=1:length(covar)
# write.csv(covar,file=paste0('out/covariates-',Sys.Date(),'.csv'))
names(d.cln)[-ncol(d.cln)]=paste0('z',1:(ncol(d.cln)-1) )
outcome.name='Das28_LDA'
f <- as.formula(paste0(outcome.name,"~",paste0(names(d.cln)[ind.all],collapse="+")))

## one time training and validation: training using 80% of data
train_prop = .8 
if (common==1){train_prop =1}
i=1
set.seed(2024+i)
train_index = sample(1:nrow(d.cln), ceiling(train_prop*nrow(d.cln)), replace=F)
valid_index = setdiff(1:nrow(d.cln),train_index)
d.train = d.cln %>% slice(train_index)
d.valid = d.cln %>% slice(valid_index)

## superlearner with mean, random forest, XGBoost, and neural 
SL.brulee <- function(Y, X, newX, family, ...) {
  f <- as.formula(paste0("~",paste0(colnames(X),collapse="+")))
  fit.brulee <- brulee_mlp(x = model.matrix(f, data=X), y = as.factor(Y) %>% relevel(ref="0"), penalty = 0.03, hidden_units = c(64,32,10), dropout=0)
  
  pred <- predict(fit.brulee, new_data=model.matrix(f, data=newX), type="prob") %>% pull(".pred_1")
  fit <- list(object=fit.brulee)
  
  out <- list(pred=pred, fit=fit)
  class(out$fit) <- c("SL.brulee")
  return(out)
}

predict.SL.brulee <- function(object, newdata, ...) {
  f <- as.formula(paste0("~",paste0(colnames(newdata),collapse="+")))
  
  pred <- predict(object$object, new_data=model.matrix(f, data=newdata), type="prob") %>% pull(".pred_1")
  return(pred)
}

#
  time.start=Sys.time()
  m.sl.lib <- c("SL.mean","SL.ranger", "SL.xgboost","SL.brulee")
  # m.sl.lib <- c("SL.ranger")
  m.sl <- SuperLearner(Y=as.numeric(d.train$Das28_LDA)-1, X=d.train %>% select(-c(Das28_LDA)), 
                       family = binomial(), SL.library = m.sl.lib)
  pred <- predict(m.sl, d.valid %>% select(-c(Das28_LDA)), onlySL = TRUE)$pred
  AUC.sl=auc(response=d.valid$Das28_LDA,predictor=pred)[[1]]
  library(pROC)
  roc_obj <- roc(response = d.valid$Das28_LDA, predictor = pred)
  auc_value <- auc(roc_obj)
  print(auc_value)
  ci_auc <- ci.auc(roc_obj)
  print(ci_auc)
  time.end=Sys.time()
  time=time.end-time.start
  print(time);print(AUC.sl) 
  res=list('time'=time,'fit.sl'=m.sl,'covar'=covar,'AUC.sl'=AUC.sl,'pred'=pred ) 
  save(res,file=paste0("out/sl",set,'-rf-',Sys.Date(),".rda"))
  if (common==1){
    res=list('fit.sl'=m.sl,'covar'=covar ) 
    save(res,file=paste0("out/MGBcommon-sl",set,'-',Sys.Date(),".rda")) }
# }
  
## permutation importance
permute_importance <- function(model, X, Y, metric = function(y, yhat) mean((y - yhat)^2)) {
  baseline <- metric(Y, predict(model, newdata = X)$pred)
  sapply(names(X), function(var) {
    X_perm <- X
    X_perm[[var]] <- sample(X_perm[[var]])
    permuted_perf <- metric(Y, predict(model, newdata = X_perm)$pred)
    permuted_perf - baseline
  })
}
important=permute_importance(res$fit.sl, X=d.train %>% select(-c(Das28_LDA)), Y=as.numeric(d.train$Das28_LDA)-1)
importance=data.frame('Code'=covar,'importance'=important)
# list.screen$Code[239:248]=gsub(":","",list.screen$Code[239:248])
# list.screen$Code[239:248]=gsub("\\.","",list.screen$Code[239:248])
# list.screen.importance=left_join(list.screen,importance,by='Code')
# write_csv(list.screen.importance,'list.screen.importance.csv')

## sharp importance 
library(iml)
predict_sl <- function(model, newdata) {
  predict(model, newdata = newdata)$pred
}
predictor <- Predictor$new(
  model = m.sl,
  data = d.train %>% select(-c(Das28_LDA)),
  y =as.numeric(d.train$Das28_LDA)-1,
  predict.function = predict_sl
)
# shapley <- Shapley$new(predictor, x.interest = (d.train %>% select(-c(Das28_LDA)))[1, ])  # explanation for first obs
# print(shapley$results)
shap_imp <- FeatureImp$new(predictor, loss = "mae")
plot(shap_imp)
importance=data.frame('Code'=covar[match(shap_imp$results$feature,names(d.cln)[-ncol(d.cln)])],
                      'importance'=shap_imp$results$importance,'importance.05'=shap_imp$results$importance.05,'importance.95'=shap_imp$results$importance.95)
# list.screen$Code[239:248]=gsub(":","",list.screen$Code[239:248])
# list.screen$Code[239:248]=gsub("\\.","",list.screen$Code[239:248])
# list.screen.importance=left_join(importance,list.screen,by='Code')
# write_csv(list.screen.importance,file=paste0("out/list.screen.importance-",Sys.Date(),".csv"))

## plot top features
top_feats <- list.screen.importance %>%
  arrange(desc(importance)) %>%
  slice_head(n = 30)
# top_feats$Description[is.na(top_feats$Description)]=top_feats$Code[is.na(top_feats$Description)]
top_feats$Code=paste0(top_feats$Code,':',top_feats$Description)
pdf(file=paste0("out/Imp",Sys.Date(),"_common.pdf"))
ggplot(top_feats, aes(x = reorder(Code, importance), y = importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    x = "Feature",
    y = "Importance (increase in MAE)",
    title = "Top 30 Features by Permutation Importance"
  ) +
  theme_minimal()
dev.off()
}
 

#### transportability investigation
{
load("/Users/xuanwang/Library/CloudStorage/Box-Box/RACAT/VA/VAcommon.sl0-2025-05-23.rda")
fit=res$fit.sl
pred=predict(fit,d.cln %>% select(-c(Das28_LDA)),onlySL=TRUE)$pred
AUC.sl=auc(response=d.cln$Das28_LDA,predictor=pred)[[1]]
print(AUC.sl)
roc_obj <- roc(response = d.cln$Das28_LDA, predictor = pred)
auc_value <- auc(roc_obj)
print(auc_value)
ci_auc <- ci.auc(roc_obj)
print(ci_auc)

}




