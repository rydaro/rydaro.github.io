# ER OUTCOMES POISSON 180 DAYS

fulltib180<-fulltib2 %>%
  filter(enrolltime >=180)

ae <-fulltib180[fulltib180$Brnd_Nm %in% c('PROVENGE','ZYTIGA','XTANDI'),]

ae <-ae %>%
  select(Patid,provnum,agecat,racecat,educat,housecat,Division,
         Product,met,Aso,diabetes,hypertension,CHF,osteoporosis,arrythmia,uro,er60,ercount180,enrolltime) %>% na.omit

ae$timetreat <-ae$enrolltime +1

ae$Division<-as.numeric(ae$Division)
ae$Product<-as.numeric(ae$Product)


aepois <-glm(ercount180~provnum,data=ae,family=poisson(link="log"))

summary(aepois)


for(i in 2:length(aepois$coefficients)){
  out <-data.frame(t(exp(summary(aepois)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(aepois)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(aepois$coefficients[i])
  print(out)
}


aepois<-glm(ercount180~provnum +as.factor(agecat)+as.factor(racecat)+as.factor(educat)+
      as.factor(housecat)+as.factor(Division)+
      as.factor(Product)+as.factor(met)+as.factor(Aso)+as.factor(diabetes)+
      as.factor(hypertension)+
      as.factor(CHF)+as.factor(osteoporosis)+as.factor(arrythmia)+
      as.factor(uro),data=ae,family=poisson(link=log))

summary(aepois)


for(i in 2:length(aepois$coefficients)){
  out <-data.frame(t(exp(summary(aepois)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(aepois)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(aepois$coefficients[i])
  print(out)
}

matched <- matchit(provnum~as.factor(agecat)+as.factor(racecat)+as.factor(educat)+as.factor(housecat)+as.factor(Division)+
                     as.factor(Product)+as.factor(met)+as.factor(Aso)+as.factor(diabetes)+as.factor(hypertension)+
                     as.factor(CHF)+as.factor(osteoporosis)+as.factor(arrythmia)+
                     as.factor(uro),data=ae,distance = 'logit',
                   #method='optimal',caliper=.15)
                   method = 'full',discard='both',max.controls=5,caliper=.15)
#mehtod='nearest',caliper=.15)
summary<-summary(matched)
summary$nn

dta_m <- match.data(matched)
dim(dta_m)
colnames(dta_m)
dta_m$Division<-as.numeric(dta_m$Division)
dta_m$Product<-as.numeric(dta_m$Product)

matchedmodreg <-glm(ercount180~provnum,data=dta_m,family=poisson(link=log))
summary(matchedmodreg)
for(i in 2:length(matchedmodreg$coefficients)){
  out <-data.frame(t(exp(summary(matchedmodreg)$coefficients[i,1] + qnorm(c(0.5,0.025,0.975)) * summary(matchedmodreg)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(matchedmodreg$coefficients[i])
  print(out)
}


stratae <-glm(provnum~as.factor(agecat)+as.factor(racecat)+as.factor(educat)+as.factor(housecat)+as.factor(Division)+
                as.factor(Product)+as.factor(met)+as.factor(Aso)+as.factor(diabetes)+as.factor(hypertension)+
                as.factor(CHF)+as.factor(osteoporosis)+as.factor(arrythmia)+
                as.factor(uro),data=ae,family=binomial(link='logit'))

ae$pr_score <-predict(stratae, type = "response")


summary(ae$pr_score)

ae$pr_score_trim <-ifelse(ae$pr_score<.01,.01,ae$pr_score)

ae <- within(ae, quintile <- as.integer(cut(pr_score, quantile(pr_score, probs=0:5/5), include.lowest=TRUE)))

summary(ae$quintile)

round(quantile(ae$pr_score,probs=0:5/5),3)

ae$IPTW <-ae$provnum/ae$pr_score_trim + (1-ae$provnum)/(1-ae$pr_score_trim)

prep <-glm(ercount180~provnum+pr_score,data=ae,family=poisson(link=log))
summary(prep)

for(i in 2:length(prep$coefficients)){
  out <-data.frame(t(exp(summary(prep)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(prep)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(prep$coefficients[i])
  print(out)
}

###########

aeiptw<-glm(ercount180~provnum +as.factor(agecat)+as.factor(racecat)+as.factor(educat)+
              as.factor(housecat)+as.factor(Division)+
              as.factor(Product)+as.factor(met)+as.factor(Aso)+as.factor(diabetes)+
              as.factor(hypertension)+
              as.factor(CHF)+as.factor(osteoporosis)+as.factor(arrythmia)+
              as.factor(uro),data=ae,family=poisson(link=log),weights=as.vector(IPTW))

summary(aeiptw)



for(i in 2:length(aeiptw$coefficients)){
  out <-data.frame(t(exp(summary(aeiptw)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(aeiptw)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(aeiptw$coefficients[i])
  print(out)
}

#############################################################
#################################################################
################################################################
################################################################

