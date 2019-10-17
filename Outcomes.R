#outcomes
#ERMerge and firstline_regression (at least beginning) must be ran prior
library(survival)
#set.seed(17891)
library(dplyr)
library(mgcv)
library(MatchIt)
library(ggplot2)
#firstline regression 
firstline <-read.csv('/Users/ryanross/Downloads/firstline5575patients.csv')
pats <-unique(firstline$Patid)


lasts <-read.csv("~/OptumProstate/enroltime.csv")
#lasts$lastdate <-as.Date(lasts$lastdate)

lasts$lastdate<-as.Date(lasts$lastdate)
lasts$Patid <-as.character(lasts$Patid)

lasts<-lasts[,1:3]


#firstline <-left_join(firstline,lasts,by="Patid")

prolast <-left_join(pro,lasts,by="Patid")

enrolltime <- prolast %>%
  group_by(Patid) %>%
  summarize(enrolltime = max(lastdate,as.Date(Date)) - min(as.Date(Date))) %>%
  filter(Patid %in% pats)

enrolltime$enrolltime<-as.numeric(enrolltime$enrolltime)
enrolltime$Patid <-as.numeric(enrolltime$Patid)

firstline <-left_join(firstline,enrolltime,by="Patid")
#firstline$enrolltime<-as.numeric(firstline$enrolltime)
firstline$enrolltime<-firstline$enrolltime +1

#summary(firstline$enrolltime[as.Date(firstline$Date)>="2014-01-01"])

summary(firstline$enrolltime[firstline$Brnd_Nm=="PROVENGE"])

summary(firstline$enrolltime[firstline$Brnd_Nm=="DOCETAXEL"])

summary(firstline$enrolltime[firstline$Brnd_Nm=="ZYTIGA"|firstline$Brnd_Nm=="XTANDI"])

enroll_drug <- pro %>%
  group_by(Patid) %>%
  summarize(enroll_drug = max(as.Date(Date)) - min(as.Date(Date))) %>%
  filter(Patid %in% pats)

enroll_drug$Patid<-as.numeric(enroll_drug$Patid)
enroll_drug$enroll_drug<-as.numeric(enroll_drug$enroll_drug)

firstline <-left_join(firstline,enroll_drug,by="Patid")

firstline$enroll_drug<-firstline$enroll_drug +1

summary(firstline$enroll_drug[firstline$Brnd_Nm=="PROVENGE"])

summary(firstline$enroll_drug[firstline$Brnd_Nm=="DOCETAXEL"])

summary(firstline$enroll_drug[firstline$Brnd_Nm=="ZYTIGA"|firstline$Brnd_Nm=="XTANDI"])


tibble <-firstline %>%
  group_by(Patid) %>%  
  arrange(Year,Date)%>%
  select(Patid,Year,Date,Brnd_Nm,Yrdob,Aso,Bus,Cdhp,Prov_Type,Gdr_Cd,
         Product,RACE_CD,Division,D_EDUCATION_LEVEL_CODE,
         D_HOUSEHOLD_INCOME_RANGE_CODE,class,diabetes,hypertension,
         arrythmia,CHF,osteoporosis,met,enrolltime,enroll_drug) %>%
  mutate(age = Year - Yrdob) %>%
  filter(row_number()==1)

tibble2 <- firstline %>%
  group_by(Patid) %>%
  summarise(keep=any(Brnd_Nm %in% c('PROVENGE','DOCETAXEL','ZYTIGA','XTANDI')))



fulltib <-inner_join(tibble2,tibble,by="Patid")
fulltib <-fulltib[fulltib$keep==1,]

fulltib$housecat <-ifelse(fulltib$D_HOUSEHOLD_INCOME_RANGE_CODE==1 |
                            fulltib$D_HOUSEHOLD_INCOME_RANGE_CODE==2,1,
                          ifelse(fulltib$D_HOUSEHOLD_INCOME_RANGE_CODE==3 |
                                   fulltib$D_HOUSEHOLD_INCOME_RANGE_CODE==4 |
                                   fulltib$D_HOUSEHOLD_INCOME_RANGE_CODE==5,2,
                                 ifelse(fulltib$D_HOUSEHOLD_INCOME_RANGE_CODE==6,3,4)))

fulltib$housecat <-factor(fulltib$housecat,levels=c(1,2,3,4),labels =c("<50k","50-99k",">99k","Unknown"))

fulltib$provider<-ifelse(fulltib$class=="Facility" | fulltib$class=="Other individuals" | 
                           fulltib$class=="Radiation oncologist",'Other',
                         as.character(fulltib$class))
fulltib$provnum <-as.numeric(ifelse(fulltib$Brnd_Nm=='PROVENGE',1,0))


fulltib$agecat <-ifelse(fulltib$age<55,0,ifelse(fulltib$age>=55&fulltib$age<65,1,ifelse(fulltib$age>=65&fulltib$age<75,2,3)))
fulltib$agecat <- factor(fulltib$agecat,levels = c(0,1,2,3),labels=c("<55","55-64","65-74",">75"))

fulltib$racecat <-ifelse(fulltib$RACE_CD==6,0,ifelse(fulltib$RACE_CD==2,1,ifelse(fulltib$RACE_CD==3,2,
                                                                                 ifelse(fulltib$RACE_CD==4,3,4))))

fulltib$racecat <-factor(fulltib$racecat,levels=c(0,1,2,3,4),labels=c("White","Asian","Black","Hispanic","Unknown"))

fulltib$educat <-ifelse(fulltib$D_EDUCATION_LEVEL_CODE %in% c('A','B'),0,ifelse(fulltib$D_EDUCATION_LEVEL_CODE %in% c('C','D'),1,2))
fulltib$educat <-factor(fulltib$educat, levels=c(0,1,2), labels=c('No College',"Some College or More","Uknown"))


fulltib$uro <-ifelse(is.na(fulltib$class),0,ifelse(fulltib$class=='Urologist',1,0))


fulltib1 <-full_join(fulltib,ercounts60,by='Patid')
fulltib1$ercount60[is.na(fulltib1$ercount60)] <-0
fulltib1$er60 <-ifelse(fulltib1$ercount60>=1,1,0)
fulltib2 <-full_join(fulltib1,ercounts180,by='Patid')
fulltib2$ercount180[is.na(fulltib2$ercount180)]<-0


fulltib1 <-full_join(fulltib,ercounts60,by='Patid')
fulltib1$ercount60[is.na(fulltib1$ercount60)] <-0
fulltib1$er60 <-ifelse(fulltib1$ercount60>=1,1,0)

fulltib2 <-full_join(fulltib1,ercounts180,by='Patid')
fulltib2$ercount180[is.na(fulltib2$ercount180)]<-0

fulltib2$pae <-ifelse(fulltib2$Brnd_Nm=="PROVENGE",0,ifelse(fulltib2$Brnd_Nm %in% c("ZYTIGA","XTANDI"),1,99))
fulltib2$pdoc <-ifelse(fulltib2$Brnd_Nm=="PROVENGE",0,ifelse(fulltib2$Brnd_Nm=="DOCETAXEL",1,99))

fulltib2$Brand <-factor(fulltib2$pdoc,levels=c(0,1,99),labels = c("Sipuleucel-T","Docetaxel","Abiraterone or Enzalutamide"))



firstline_single <-fulltib2 %>%
  select(Patid,Brand,pae,pdoc,agecat,racecat,educat,housecat,Division,
         Product,met,Aso,diabetes,hypertension,CHF,osteoporosis,arrythmia,uro,enrolltime,
         er60,ercount180,enroll_drug) %>%
  na.omit

#cols <- c("agecat", "racecat", "educat", "housecat","met","Aso","diabetes","hypertension","CHF","osteoporosis","arrythmia","uro")

#firstline_single[cols] <- lapply(firstline_single[cols], factor)

firstline_single

save(firstline_single,file="firstline.Rdata")


fulltib60<-firstline_single %>%
  filter(enroll_drug >=60)

######################################################################

ae <-fulltib60[fulltib60$Brnd_Nm %in% c('PROVENGE','ZYTIGA','XTANDI'),]

ae <-ae %>%
  select(Patid,provnum,agecat,racecat,educat,housecat,Division,
         Product,met,Aso,diabetes,hypertension,CHF,osteoporosis,arrythmia,uro,er60,ercount180,enrolltime)


ae$Division<-as.numeric(ae$Division)
ae$Product<-as.numeric(ae$Product)



matched <- matchit(provnum~as.factor(agecat)+as.factor(racecat)+as.factor(educat)+as.factor(housecat)+as.factor(Division)+
                     as.factor(Product)+as.factor(met)+as.factor(Aso)+as.factor(diabetes)+as.factor(hypertension)+
                     as.factor(CHF)+as.factor(osteoporosis)+as.factor(arrythmia)+
                     as.factor(uro),data=ae,distance = 'logit',
                   #method='optimal',caliper=.15)
                   method = 'full',discard='both',caliper=.2,max.controls=20)
#mehtod='nearest',caliper=.15)
summary<-summary(matched)
summary$nn
summary$reduction

s.out <- summary(matched, standardize=TRUE)
#plot(s.out)


dta_m <- match.data(matched)
dim(dta_m)
colnames(dta_m)
dta_m$Division<-as.numeric(dta_m$Division)
dta_m$Product<-as.numeric(dta_m$Product)




#Bias
#pre-match
control <-ae[ae$provnum==0,]
ncon<-nrow(control)
treat <-ae[ae$provnum==1,]
ntre<-nrow(treat)
var <-NULL
varval <-NULL
bias<-NULL
label<-NULL
counter <-0
for(j in 2:16){
  for(k in 1:nrow(unique(ae[,j]))){
    counter <-counter +1
    tab <-unique(ae[,j])
    var[counter] <-names(tab)
    varval[counter]<-as.numeric(tab[k,])
    bias[counter]<-(mean(treat[,j]==as.numeric(tab[k,])) - 
                      mean(control[,j]==as.numeric(tab[k,]))) /sd(treat[,j]==as.numeric(tab[k,])) 
    label[counter]<-paste0(var[counter],varval[counter])
  }
}

#post match 
control1 <-dta_m[dta_m$provnum==0,]
ncon1<-nrow(control)
treat1 <-dta_m[dta_m$provnum==1,]
ntre1<-nrow(treat)
var1 <-NULL
varval1 <-NULL
bias1<-NULL
label1<-NULL
counter1 <-0
for(j in 2:17){
  for(k in 1:length(unique(dta_m[,j]))){
    counter1 <-counter1 +1
    tab <-unique(dta_m[,j])
    var1[counter1] <-colnames(dta_m)[j]
    varval1[counter1]<-tab[k]
    bias1[counter1]<-(mean(treat1[,j]==tab[k]) - 
                        weighted.mean(control1[,j]==tab[k],control1$weights)) / sd(treat[,j]==as.numeric(tab[k]))
    label1[counter1]<-paste0(var1[counter1],varval1[counter1])
  }
}


forest <-data.frame(label,bias,match=rep("unmatched",length(label)))
forest1 <-data.frame(label=label1,bias=bias1,match=rep("matched",length(label1)))

forestall <-rbind(forest,forest1)

fp <- ggplot(data =forestall,aes(x=label, y=bias,color=match,shape=match)) +
  geom_point(size=2.5) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Label") + ylab("Bias") +
  theme_minimal()+ theme(axis.ticks.y=element_blank(),
                         panel.grid.minor=element_blank(),
                         legend.title=element_blank()) +
  ggtitle("Standardized Bias of Each Covariate")
print(fp)

#Now look at outcomes
#unmatched first
summary(ae$er60[ae$provnum==0])
summary(ae$er60[ae$provnum==1])

t.test(ae$er60[ae$provnum==0],ae$er60[ae$provnum==1])

sip <-ae[ae$provnum==1,]
nosip <-ae[ae$provnum==0,]

res <- prop.test(x = c(sum(nosip$er60), sum(sip$er60)), 
                 n = c(nrow(nosip), nrow(sip)))

res

modfullgroup <-glm(er60~provnum,data=ae,family=binomial(link=logit))
summary(modfullgroup)
for(i in 2:length(modfullgroup$coefficients)){
  out <-data.frame(t(exp(summary(modfullgroup)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(modfullgroup)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(modfullgroup$coefficients[i])
  print(out)
}

#AUC
# prs_df <- data.frame(pr_score = predict(modfullgroup, type = "response"),
#                      er = modfullgroup$model$er)
# head(prs_df)
# 
# roc_full_resolution <- roc(prs_df$er, prs_df$pr_score)
# plot(roc_full_resolution, print.auc=TRUE)



modfullgroup2 <-glm(er60~provnum +as.factor(agecat)+as.factor(racecat)+as.factor(educat)+
                     as.factor(housecat)+as.factor(Division)+
                     as.factor(Product)+as.factor(met)+as.factor(Aso)+as.factor(diabetes)+
                     as.factor(hypertension)+
                     as.factor(CHF)+as.factor(osteoporosis)+as.factor(arrythmia)+
                     as.factor(uro),data=ae,family=binomial(link=logit))
summary(modfullgroup2)
for(i in 2:length(modfullgroup2$coefficients)){
  out <-data.frame(t(exp(summary(modfullgroup2)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(modfullgroup2)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(modfullgroup2$coefficients[i])
  print(out)
}
anova(modfullgroup2)

#AUC
# prs_df <- data.frame(pr_score = predict(modfullgroup2, type = "response"),
#                      er = modfullgroup2$model$er)
# head(prs_df)
# 
# roc_full_resolution <- roc(prs_df$er, prs_df$pr_score)
# plot(roc_full_resolution, print.auc=TRUE)
# 
# 
# #Confusion
# confusionMatrix(as.factor(round(prs_df$pr_score)),as.factor(prs_df$er))

#matched
summary(dta_m$er60[dta_m$provnum==0])
summary(dta_m$er60[dta_m$provnum==1])

t.test(dta_m$er60[dta_m$provnum==0],dta_m$er60[dta_m$provnum==1])

sip <-dta_m[dta_m$provnum==1,]
nosip <-dta_m[dta_m$provnum==0,]

res <- prop.test(x = c(sum(nosip$er60), sum(sip$er60)), 
                 n = c(nrow(nosip), nrow(sip)))

res

matchedmodreg <-glm(er60~provnum,data=dta_m,family=binomial(link='logit'))
summary(matchedmodreg)
for(i in 2:length(matchedmodreg$coefficients)){
  out <-data.frame(t(exp(summary(matchedmodreg)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(matchedmodreg)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(matchedmodreg$coefficients[i])
  print(out)
}

#matchedmod <-clogit(er60~provnum,data=dta_m,method="approximate",subset = subclass )
#summary(matchedmod)


######
# prs_df <- data.frame(pr_score = predict(matchedmod, type = "expected"),
#                      er = dta_m$er)
# head(prs_df)
# 
# roc_full_resolution <- roc(prs_df$er, prs_df$pr_score)
# plot(roc_full_resolution, print.auc=TRUE)
# 
# 
# #Confusion
# confusionMatrix(as.factor(round(prs_df$pr_score)),as.factor(prs_df$er))

######    STRATIFICATION    #####################################

#A&E vs. SIP
#define strata cuts

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

quint1<-t.test(ae$er60[ae$provnum==1&ae$quintile==1],ae$er60[ae$provnum==0&ae$quintile==1])
quint1<-data.frame(quint1[4])
quint2<-t.test(ae$er60[ae$provnum==1&ae$quintile==2],ae$er60[ae$provnum==0&ae$quintile==2])
quint2<-data.frame(quint2[4])
quint3<-t.test(ae$er60[ae$provnum==1&ae$quintile==3],ae$er60[ae$provnum==0&ae$quintile==3])
quint3<-data.frame(quint3[4])
quint4<-t.test(ae$er60[ae$provnum==1&ae$quintile==4],ae$er60[ae$provnum==0&ae$quintile==4])
quint4<-data.frame(quint4[4])
quint5<-t.test(ae$er60[ae$provnum==1&ae$quintile==5],ae$er60[ae$provnum==0&ae$quintile==5])
quint5<-data.frame(quint5[4])

quintile_df <-data.frame(rbind(c(mean(quint1[,1]),max(quint1[,1]) - mean(quint1[,1]),1),
                               c(mean(quint2[,1]),max(quint2[,1]) - mean(quint2[,1]),2),
                               c(mean(quint3[,1]),max(quint3[,1]) - mean(quint3[,1]),3),
                               c(mean(quint4[,1]),max(quint4[,1]) - mean(quint4[,1]),4),
                               c(mean(quint5[,1]),max(quint5[,1]) - mean(quint5[,1]),5)))
colnames(quintile_df)<-c('mean','sd','quintile')

p<- ggplot(quintile_df, aes(x=quintile, y=mean)) + 
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2)

print(p)

#save to talk about, but results insignificant
# strataeglm <-glm(er60~provnum+quintile+provnum*quintile,data=ae,family = binomial(link='logit'))
# 
# summary(strataeglm)
# 
# for(i in 2:length(strataeglm$coefficients)){
#   out <-data.frame(t(exp(summary(strataeglm)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(strataeglm)$coefficients[i,2])))
#   colnames(out) <-c('OR','Lower','Upper')
#   rownames(out)<-names(strataeglm$coefficients[i])
#   print(out)
# }

prep <-glm(er60~provnum+pr_score,data=ae,family = binomial(link="logit"))
summary(prep)

for(i in 2:length(prep$coefficients)){
  out <-data.frame(t(exp(summary(prep)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(prep)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(prep$coefficients[i])
  print(out)
}

splineae <-gam(er60~provnum+s(pr_score,k=10,m=3),data=ae,family=binomial(link='logit'))
summary(splineae)
plot.gam(splineae,shade= TRUE)
gam.check(splineae)

exp(.1040 +qnorm(c(0.5,0.025,0.975)) *(.2567) )


summary(ae$IPTW)

summary(ae$pr_score)


ggplot(ae,aes(x = IPTW, color=as.factor(provnum),fill=as.factor(provnum))) +
  geom_density(alpha=.47) +
  xlab("IPTW of getting Sip-T")


aeiptw<-glm(er60~provnum +as.factor(agecat)+as.factor(racecat)+as.factor(educat)+
              as.factor(housecat)+as.factor(Division)+
              as.factor(Product)+as.factor(met)+as.factor(Aso)+as.factor(diabetes)+
              as.factor(hypertension)+
              as.factor(CHF)+as.factor(osteoporosis)+as.factor(arrythmia)+
              as.factor(uro),data=ae,family=binomial(link=logit),weights=as.vector(IPTW))

summary(aeiptw)



for(i in 2:length(aeiptw$coefficients)){
  out <-data.frame(t(exp(summary(aeiptw)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(aeiptw)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(aeiptw$coefficients[i])
  print(out)
}

####################################################################
###################################################################
###################################################################
###################################################################
###################################################################








