##Outcomes Docetaxel Group

#run first chunk of Outcomes for dataset creations
#Docetaxel
#MOVE TO NEW SCRIPT to avoid confusion

doc <-fulltib60[fulltib60$Brnd_Nm %in% c('PROVENGE','DOCETAXEL'),]

doc <-doc %>%
  select(Patid,provnum,agecat,racecat,educat,housecat,Division,
         Product,met,Aso,diabetes,hypertension,CHF,osteoporosis,arrythmia,uro,er60,ercount180,enrolltime) %>% na.omit

doc$timetreat <-doc$enrolltime

doc$Division<-as.numeric(doc$Division)
doc$Product<-as.numeric(doc$Product)


matched <- matchit(provnum~as.factor(agecat)+as.factor(racecat)+as.factor(educat)+as.factor(housecat)+as.factor(Division)+
                     as.factor(Product)+as.factor(met)+as.factor(Aso)+as.factor(diabetes)+as.factor(hypertension)+
                     as.factor(CHF)+as.factor(osteoporosis)+as.factor(arrythmia),
                   data=doc,
                   distance = 'logit',
                   #method='optimal',caliper=.15)
                   method = 'full',discard='both',max.controls=10,caliper=.15)
#mehtod='nearest',caliper=.15)
summary<-summary(matched)
summary$nn
summary$reduction

#full model
mod <-glm(provnum~as.factor(agecat)+as.factor(racecat)+as.factor(educat)+as.factor(housecat)+as.factor(Division)+
            as.factor(Product)+as.factor(met)+as.factor(Aso)+as.factor(diabetes)+as.factor(hypertension)+
            as.factor(CHF)+as.factor(osteoporosis)+as.factor(arrythmia),
          data=doc,family=binomial(link=logit))
summary(mod)

for(i in 2:length(mod$coefficients)){
  out <-data.frame(t(exp(summary(mod)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(mod)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(mod$coefficients[i])
  print(out)
}

anova(mod,test="Chisq")


prs_df <- data.frame(pr_score = predict(mod, type = "response"),
                     provenge = mod$model$provnum)
head(prs_df)


ggplot(prs_df,aes(x = pr_score, color=as.factor(provenge),fill=as.factor(provenge))) +
  geom_density(alpha=.47) +
  xlab("Probability of getting Sip-T")



s.out <- summary(matched, standardize=TRUE)
#plot(s.out)


dta_m <- match.data(matched)
dim(dta_m)
colnames(dta_m)
dta_m$Division<-as.numeric(dta_m$Division)
dta_m$Product<-as.numeric(dta_m$Product)




#Bias
#pre-match
control <-doc[doc$provnum==0,]
ncon<-nrow(control)
treat <-doc[doc$provnum==1,]
ntre<-nrow(treat)
var <-NULL
varval <-NULL
bias<-NULL
label<-NULL
counter <-0
for(j in 2:16){
  for(k in 1:nrow(unique(doc[,j]))){
    counter <-counter +1
    tab <-unique(doc[,j])
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
  ggtitle("Standardized Bias of Each Covariate Full")
print(fp)


#Now look at outcomes
#unmatched first
summary(doc$er60[doc$provnum==0])
summary(doc$er60[doc$provnum==1])

t.test(doc$er60[doc$provnum==0],doc$er60[doc$provnum==1])

sip <-doc[doc$provnum==1,]
nosip <-doc[doc$provnum==0,]

res <- prop.test(x = c(sum(nosip$er60), sum(sip$er60)), 
                 n = c(nrow(nosip), nrow(sip)))

res

modfullgroup <-glm(er60~provnum,data=doc,family=binomial(link=logit))
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
                      as.factor(CHF)+as.factor(osteoporosis)+as.factor(arrythmia),
                    data=doc,family=binomial(link=logit))
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






#docetaxel
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


doc <- within(doc, quintile <- as.integer(cut(pr_score, quantile(pr_score, probs=0:5/5), include.lowest=TRUE)))

summary(doc$quintile)

round(quantile(doc$pr_score,probs=0:5/5),3)


quint1<-t.test(doc$er60[doc$provnum==1&doc$quintile==1],doc$er60[doc$provnum==0&doc$quintile==1])
quint1<-data.frame(quint1[4])
quint2<-t.test(doc$er60[doc$provnum==1&doc$quintile==2],doc$er60[doc$provnum==0&doc$quintile==2])
quint2<-data.frame(quint2[4])
quint3<-t.test(doc$er60[doc$provnum==1&doc$quintile==3],doc$er60[doc$provnum==0&doc$quintile==3])
quint3<-data.frame(quint3[4])
quint4<-t.test(doc$er60[doc$provnum==1&doc$quintile==4],doc$er60[doc$provnum==0&doc$quintile==4])
quint4<-data.frame(quint4[4])
quint5<-t.test(doc$er60[doc$provnum==1&doc$quintile==5],doc$er60[doc$provnum==0&doc$quintile==5])
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

# stratdocglm <-glm(er60~provnum+quintile+provnum*quintile,data=doc,family = binomial(link='logit'))
# 
# summary(stratdocglm)
# 
# for(i in 2:length(stratdocglm$coefficients)){
#   out <-data.frame(t(exp(summary(stratdocglm)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(stratdocglm)$coefficients[i,2])))
#   colnames(out) <-c('OR','Lower','Upper')
#   rownames(out)<-names(stratdocglm$coefficients[i])
#   print(out)
# }

prep <-glm(er60~provnum+pr_score,data=doc,family = binomial(link="logit"))

for(i in 2:length(prep$coefficients)){
  out <-data.frame(t(exp(summary(prep)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(prep)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(prep$coefficients[i])
  print(out)
}

splinedoc <-gam(er60~provnum+s(pr_score,k=10,m=3),data=doc,family=binomial(link='logit'))
summary(splinedoc)
plot.gam(splinedoc,all.terms = TRUE, shade= TRUE)
gam.check(splinedoc)





ggplot(doc,aes(x = IPTW, color=as.factor(provnum),fill=as.factor(provnum))) +
  geom_density(alpha=.47) +
  xlab("IPTW of getting Sip-T")


dociptw<-glm(er60~provnum +as.factor(agecat)+as.factor(racecat)+as.factor(educat)+
               as.factor(housecat)+as.factor(Division)+
               as.factor(Product)+as.factor(met)+as.factor(Aso)+as.factor(diabetes)+
               as.factor(hypertension)+
               as.factor(CHF)+as.factor(osteoporosis)+as.factor(arrythmia),data=doc,family=binomial(link=logit),weights=as.vector(IPTW))

summary(dociptw)



for(i in 2:length(dociptw$coefficients)){
  out <-data.frame(t(exp(summary(dociptw)$coefficients[i,1] +     qnorm(c(0.5,0.025,0.975)) * summary(dociptw)$coefficients[i,2])))
  colnames(out) <-c('OR','Lower','Upper')
  rownames(out)<-names(dociptw$coefficients[i])
  print(out)
}
