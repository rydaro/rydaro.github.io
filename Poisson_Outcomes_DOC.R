#poisson Outcome Docetaxel
#check to make sure you are using right matched data set

fulltib180<-fulltib2 %>%
  filter(enrolltime >=180)

doc <-fulltib180[fulltib180$Brnd_Nm %in% c('PROVENGE','DOCETAXEL'),]

doc <-doc %>%
  select(Patid,provnum,agecat,racecat,educat,housecat,Division,
         Product,met,Aso,diabetes,hypertension,CHF,osteoporosis,arrythmia,uro,er60,ercount180,enrolltime) %>% na.omit

doc$timetreat <-doc$enrolltime +1

doc$Division<-as.numeric(doc$Division)
doc$Product<-as.numeric(doc$Product)

docpois <-glm(ercount180~provnum,data=doc,family=poisson(link="log"))



summary(docpois)
pchisq(docpois$deviance, df=docpois$df.residual, lower.tail=FALSE)

for(i in 2:length(docpois$coefficients)){
  out <-data.frame(t(exp(summary(docpois)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(docpois)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(docpois$coefficients[i])
  print(out)
}


docpois<-glm(ercount180~provnum +as.factor(agecat)+as.factor(racecat)+as.factor(educat)+
               as.factor(housecat)+as.factor(Division)+
               as.factor(Product)+as.factor(met)+as.factor(Aso)+as.factor(diabetes)+
               as.factor(hypertension)+
               as.factor(CHF)+as.factor(osteoporosis)+as.factor(arrythmia)+
               as.factor(uro),data=doc,family=poisson(link=log))

summary(docpois)


for(i in 2:length(docpois$coefficients)){
  out <-data.frame(t(exp(summary(docpois)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(docpois)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(docpois$coefficients[i])
  print(out)
}


matched <- matchit(provnum~as.factor(agecat)+as.factor(racecat)+as.factor(educat)+as.factor(housecat)+as.factor(Division)+
                     as.factor(Product)+as.factor(met)+as.factor(Aso)+as.factor(diabetes)+as.factor(hypertension)+
                     as.factor(CHF)+as.factor(osteoporosis)+as.factor(arrythmia),data=doc,
                   distance = 'logit',
                   #method='optimal',caliper=.15)
                   method = 'full',discard='both',max.controls=4,caliper=.15)
#mehtod='nearest',caliper=.15)
summary<-summary(matched)
summary$nn
summary$reduction


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


stratdoc <-glm(provnum~as.factor(agecat)+as.factor(racecat)+as.factor(educat)+as.factor(housecat)+as.factor(Division)+
                 as.factor(Product)+as.factor(met)+as.factor(Aso)+as.factor(diabetes)+as.factor(hypertension)+
                 as.factor(CHF)+as.factor(osteoporosis)+as.factor(arrythmia)+
                 as.factor(uro),data=doc,family=binomial(link='logit'))

doc$pr_score <-predict(stratdoc, type = "response")

summary(doc$pr_score)
doc$pr_score_trim <-ifelse(doc$pr_score<.01,.01,doc$pr_score)
doc$pr_score_trim <-ifelse(doc$pr_score>.99,.99,doc$pr_score_trim)

doc$IPTW <-doc$provnum/doc$pr_score_trim + (1-doc$provnum)/(1-doc$pr_score_trim)
summary(doc$IPTW)

summary(doc$pr_score)



prep <-glm(ercount180~provnum+pr_score,data=doc,family=poisson(link=log))
summary(prep)

for(i in 2:length(prep$coefficients)){
  out <-data.frame(t(exp(summary(prep)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(prep)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(prep$coefficients[i])
  print(out)
}

###########

dociptw<-glm(ercount180~provnum +as.factor(agecat)+as.factor(racecat)+as.factor(educat)+
              as.factor(housecat)+as.factor(Division)+
              as.factor(Product)+as.factor(met)+as.factor(Aso)+as.factor(diabetes)+
              as.factor(hypertension)+
              as.factor(CHF)+as.factor(osteoporosis)+as.factor(arrythmia)+
              as.factor(uro),data=doc,family=poisson(link=log),weights=as.vector(IPTW))

summary(dociptw)



for(i in 2:length(dociptw$coefficients)){
  out <-data.frame(t(exp(summary(dociptw)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(dociptw)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(dociptw$coefficients[i])
  print(out)
}

#############################################################
#################################################################
################################################################
################################################################