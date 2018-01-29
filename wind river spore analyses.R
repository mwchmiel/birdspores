library(MASS)
library(ggplot2)
library(broom)
library(dplyr)

birds<-read.csv("bandingdata.csv")
head(birds)
birds<-droplevels(subset(birds, Alpha != "DESTROYED"))
birds<-droplevels(subset(birds, totgreen != "NA"))

birds$Tarsal.L<-as.numeric(as.character(birds$Tarsal.L))

birds$wing.chord<-as.numeric(as.character(birds$wing.chord))
birds$mass<-as.numeric(as.character(birds$mass))

summary(lm(wing.chord~Tarsal.L, data=birds)) #body size metrics correlated, shouldn't
#use both in a model together

#####################################################
##############Starting with total spores by
##############behavior and tarsal length
##############and total by taxon
#####################################################


birdsmod<-glm(totgreen~behavior+Tarsal.L, family = poisson, data=birds)
summary(birdsmod)
anova(birdsmod, test="Chisq")
pchisq(deviance(birdsmod), df.residual(birdsmod), lower=F) #!

plot(rstudent(birdsmod)~birdsmod$fitted.values)

############################
####Tarsal length influence on spore count, drop NA values for Tarsal length
############################

plot(totgreen~Tarsal.L, data=subset(birds, Tarsal.L != "NA"))
totgreenline<-lm(totgreen~Tarsal.L, data=subset(birds, Tarsal.L != "NA"))
abline(totgreenline)

plot(leggreen~Tarsal.L, data=subset(birds, Tarsal.L != "NA"))
leggreenline<-lm(leggreen~Tarsal.L, data=subset(birds, Tarsal.L != "NA"))
abline(leggreenline)

plot(rectgreen~Tarsal.L, data=subset(birds, Tarsal.L != "NA"))
rectgreenline<-lm(rectgreen~Tarsal.L, data=subset(birds, Tarsal.L != "NA"))
abline(rectgreenline)

##########################
#Tarsal length influence on spore count, drop spore zeroes
##########################

plot(totgreen~Tarsal.L, data=subset(birds, totgreen > 0 ))
totgreenline<-lm(totgreen~Tarsal.L, data=subset(birds, totgreen > 0))
abline(totgreenline)

plot(leggreen~Tarsal.L, data=subset(birds, leggreen > 0))
leggreenline<-lm(leggreen~Tarsal.L, data=subset(birds, leggreen > 0))
abline(leggreenline)

plot(rectgreen~Tarsal.L, data=subset(birds, rectgreen > 0))
rectgreenline<-lm(rectgreen~Tarsal.L, data=subset(birds, rectgreen > 0))
abline(rectgreenline)

#########################
#both drop tarsal NA and zero values
#########################

plot(totgreen~Tarsal.L, data=subset(birds, Tarsal.L != "NA" & totgreen > 0))
totgreenline<-lm(totgreen~Tarsal.L, data=subset(birds, Tarsal.L != "NA" & totgreen > 0))
abline(totgreenline)
summary(totgreenline)


plot(leggreen~Tarsal.L, data=subset(birds, Tarsal.L != "NA" & leggreen > 0))
leggreenline<-lm(leggreen~Tarsal.L, data=subset(birds, Tarsal.L != "NA" & leggreen > 0))
abline(leggreenline)

plot(rectgreen~Tarsal.L, data=subset(birds, Tarsal.L != "NA" & rectgreen > 0))
rectgreenline<-lm(rectgreen~Tarsal.L, data=subset(birds, Tarsal.L != "NA" & rectgreen > 0))
abline(rectgreenline)


##########################
#GLM Models
##########################


library(MASS)
# negative binomial
sum(residuals(birdsmod, type = "pearson")^2)

#quasipoisson model
birdsquasi<-glm(totgreen~behavior+Tarsal.L, family = quasipoisson, data=birds)
summary(birdsquasi)
anova(birdsquasi, test="Chisq")
pchisq(deviance(birdsquasi), df.residual(birdsquasi), lower=F) 
plot(rstudent(birdsquasi)~birdsquasi$fitted.values)


birdsnegbi<-glm.nb(totgreen~behavior+Tarsal.L, data=birds)
summary(birdsnegbi)
anova(birdsnegbi, test="Chisq")
pchisq(deviance(birdsnegbi), df.residual(birdsnegbi), lower=F) 
birdsnegbi0<-glm.nb(totgreen~1, data=subset(birds, Tarsal.L != "NA"))
1-(deviance(birdsnegbi)/deviance(birdsnegbi0)) #0.010 (goodness of fit ~R2)
m1<-update(birdsnegbi, .~. - Tarsal.L)
m1
anova(m1, birdsnegbi)
mod1overall<-anova(birdsnegbi, birdsnegbi0) #significantly better than null
plot(rstudent(birdsnegbi)~birdsnegbi$fitted.values)

model1factors<-drop1(birdsnegbi, test="LRT") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
mod1overall
drop1(birdsnegbi, test="LRT")


model1factors
birdsnegbi_tidy<-tidy(birdsnegbi)
model1factors<-tidy(model1factors)
mod1overalltidy<-tidy(mod1overall)

write.csv(rbind(birdsnegbi_tidy,model1tidy), "GLM summaries/test1.csv")
write.csv(model1factors, "GLM summaries/mod1factors.csv")
write.csv(mod1overalltidy, "GLM summaries/mod1overall.csv")

###################################################################
#totbytaxonandleg<-glm.nb(totgreen~Taxon.group+Tarsal.L, data=birds)
#summary(totbytaxonandleg)

# OVERALL SPORE NUMBER NOT SIGNIFICANT BY TAXONOMIC GROUP WHEN INCLUDING 
#TARSAL LENGTH, EXAMINE ONLY TAXONOMIC GROUP

totbytaxon<-glm.nb(totgreen~Taxon.group, data=birds)
summary(totbytaxon)
anova(totbytaxon, test="Chisq")
pchisq(deviance(totbytaxon), df.residual(totbytaxon), lower=F) #
totbytaxon0<-glm.nb(totgreen~1, data=birds)
1-(deviance(totbytaxon)/deviance(totbytaxon0)) 
mod2overall<-anova(totbytaxon, totbytaxon0, test="Chisq") #significantly better than null

mod2factors<-drop1(totbytaxon, test="LRT") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
plot(rstudent(totbytaxon)~totbytaxon$fitted.values)


totbytaxon_tidy<-tidy(totbytaxon)
totbytaxon_tidy
mod2factorstidy<-tidy(mod2factors)
mod2overalltidy<-tidy(mod2overall)
write.csv(totbytaxon_tidy, "GLM summaries/totbytaxon.csv")
write.csv(mod2factorstidy, "GLM summaries/mod2factors.csv")
write.csv(mod2overalltidy, "GLM summaries/mod2overall.csv")

##################################################################
###############leg spores############################
#################################################################

legsporebytarsalbehavior<-glm.nb(leggreen~behavior+Tarsal.L, data=birds)
summary(legsporebytarsalbehavior)
anova(legsporebytarsalbehavior, test="Chisq")
pchisq(deviance(legsporebytarsalbehavior), df.residual(legsporebytarsalbehavior), lower=F) #
legsporebytarsalbehavior0<-glm.nb(leggreen~1, data=subset(birds, Tarsal.L != "NA"))
1-(deviance(legsporebytarsalbehavior)/deviance(legsporebytarsalbehavior0))
anova(legsporebytarsalbehavior, legsporebytarsalbehavior0, test="Chisq") #significantly better than null


drop1(legsporebytarsalbehavior, test="LRT") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
#should drop behavior
##next model

legsporebytarsal<-glm.nb(leggreen~Tarsal.L, data=birds)
summary(legsporebytarsal)
anova(legsporebytarsal, test="Chisq")
pchisq(deviance(legsporebytarsal), df.residual(legsporebytarsal), lower=F) 
legsporebytarsal0<-glm.nb(leggreen~1, data=subset(birds, Tarsal.L != "NA"))
1-(deviance(legsporebytarsal)/deviance(legsporebytarsal0)) 
legmod1overall<-anova(legsporebytarsal0, legsporebytarsal, test="LR") #significantly better than null

legmod1parameters<-drop1(legsporebytarsal, test="LRT") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
plot(rstudent(legsporebytarsal)~legsporebytarsal$fitted.values)
legmod1overall
legsporebytarsal_tidy<-tidy(legsporebytarsal)
legmod1overalltidy<-tidy(legmod1overall)
legmod1parameterstidy<-tidy(legmod1parameters)
write.csv(legsporebytarsal_tidy, "GLM summaries/leggreenbytarsalfull.csv")
write.csv(legmod1parameterstidy, "GLM summaries/legmod1parameters.csv")
write.csv(legmod1overalltidy, "GLM summaries/legmod1overall.csv")



legsporebytaxonleg<-glm.nb(leggreen~Taxon.group+Tarsal.L, data=birds)
summary(legsporebytaxonleg)
anova(legsporebytaxonleg, test="Chisq")
pchisq(deviance(legsporebytaxonleg), df.residual(legsporebytaxonleg), lower=F) #
legsporebytaxonleg0<-glm.nb(leggreen~1, data=subset(birds, Tarsal.L != "NA"))
1-(deviance(legsporebytaxonleg)/deviance(legsporebytaxonleg0)) 
anova(legsporebytaxonleg, legsporebytaxonleg0, test="Chisq") #significantly better than null

drop1(legsporebytaxonleg, test="F") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
###neither significant, next model

legsporebytaxon<-glm.nb(leggreen~Taxon.group, data=birds)
summary(legsporebytaxon)
anova(legsporebytaxon, test="Chisq")
pchisq(deviance(legsporebytaxon), df.residual(legsporebytaxon), lower=F) 
legsporebytaxon0<-glm.nb(leggreen~1, data=birds)
1-(deviance(legsporebytaxon)/deviance(legsporebytaxon0)) 
legmod2overall<-anova(legsporebytaxon, legsporebytaxon0, test="Chisq") #significantly better than null
legmod2overall
legmod2parameters<-drop1(legsporebytaxon, test="LRT") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
plot(rstudent(legsporebytarsal)~legsporebytarsal$fitted.values)

legsporebytaxon_tidy<-tidy(legsporebytaxon)
legmod2overalltidy<-tidy(legmod2overall)
legmod2parameterstidy<-tidy(legmod2parameters)
write.csv(legsporebytaxon_tidy, "GLM summaries/leggreenbytaxonfull.csv")
write.csv(legmod2parameterstidy, "GLM summaries/legmod2parameters.csv")
write.csv(legmod2overalltidy, "GLM summaries/legmod2overall.csv")

####################################################################
###########################rect spores##############################
####################################################################

rectsporebytarsalbehavior<-glm.nb(rectgreen~behavior+Tarsal.L, data=birds)
summary(rectsporebytarsalbehavior)
anova(rectsporebytarsalbehavior, test="Chisq")
pchisq(deviance(rectsporebytarsalbehavior), df.residual(rectsporebytarsalbehavior), lower=F) 
rectsporebytarsalbehavior0<-glm.nb(rectgreen~1, data=subset(birds, Tarsal.L != "NA"))
1-(deviance(rectsporebytarsalbehavior)/deviance(rectsporebytarsalbehavior0)) 
anova(rectsporebytarsalbehavior, rectsporebytarsalbehavior0, test="Chisq") #significantly better than null

drop1(rectsporebytarsalbehavior, test="F") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
#should drop tarsal.L
##next model

rectsporebybehavior<-glm.nb(rectgreen~behavior, data=birds)
summary(rectsporebybehavior)
anova(rectsporebybehavior, test="Chisq")
pchisq(deviance(rectsporebybehavior), df.residual(rectsporebybehavior), lower=F) 
rectsporebybehavior0<-glm.nb(rectgreen~1, data=birds)
1-(deviance(rectsporebybehavior)/deviance(rectsporebybehavior0)) 
rectmod1overall<-anova(rectsporebybehavior, rectsporebybehavior0, test="Chisq") #significantly better than null
rectmod1parameters
rectmod1parameters<-drop1(rectsporebybehavior, test="LRT") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
plot(rstudent(rectsporebybehavior)~rectsporebybehavior$fitted.values)

rectsporebybehavior_tidy<-tidy(rectsporebybehavior)
rectmod1parameterstidy<-tidy(rectmod1parameters)
rectmod1overalltidy<-tidy(rectmod1overall)

write.csv(rectsporebybehavior_tidy, "GLM summaries/rectgreenbybehaviorfull.csv")
write.csv(rectmod1parameterstidy, "GLM summaries/rectmod1parameters.csv")
write.csv(rectmod1overalltidy, "GLM summaries/rectmod1overall.csv")



rectsporebytaxonrect<-glm.nb(rectgreen~Taxon.group+Tarsal.L, data=birds)
summary(rectsporebytaxonrect)
anova(rectsporebytaxonrect, test="Chisq")
pchisq(deviance(rectsporebytaxonrect), df.residual(rectsporebytaxonrect), lower=F) 
rectsporebytaxonrect0<-glm.nb(rectgreen~1, data=subset(birds, Tarsal.L != "NA"))
1-(deviance(rectsporebytaxonrect)/deviance(rectsporebytaxonrect0)) 
anova(rectsporebytaxonrect, rectsporebytaxonrect0, test="Chisq") #significantly better than null

drop1(rectsporebytaxonrect, test="F") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
###drop tarsal length next model

rectsporebytaxon<-glm.nb(rectgreen~Taxon.group, data=birds)
summary(rectsporebytaxon)
anova(rectsporebytaxon, test="Chisq")
pchisq(deviance(rectsporebytaxon), df.residual(rectsporebytaxon), lower=F) 
rectsporebytaxon0<-glm.nb(rectgreen~1, data=birds)
1-(deviance(rectsporebytaxon)/deviance(rectsporebytaxon0))
rectmod2overall<-anova(rectsporebytaxon, rectsporebytaxon0, test="Chisq") #significantly better than null

rectmod2overall
rectmod2parameters<-drop1(rectsporebytaxon, test="LRT") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
plot(rstudent(rectsporebytarsal)~rectsporebytarsal$fitted.values)

rectsporebytaxon_tidy<-tidy(rectsporebytaxon)
rectmod2parameterstidy<-tidy(rectmod2parameters)
rectmod2overalltidy<-tidy(rectmod2overall)

write.csv(rectsporebytaxon_tidy, "GLM summaries/rectgreenbytaxonfull.csv")
write.csv(rectmod2parameterstidy, "GLM summaries/rectmod2parameters.csv")
write.csv(rectmod2overalltidy, "GLM summaries/rectmod2overall.csv")


###################################################################



library(plyr)


sumbyspp<- ddply(birds, "Alpha", summarise,
                N    = length(Alpha),
                mean = mean(totgreen),
                sd   = sd(totgreen),
                se   = sd / sqrt(N)
                
)
sumbyspp$CI<-sumbyspp$se*1.96
head(sumbyspp)
sumbyspp

#make a  subset that only includes spp w at least 3 captures
smalist<-droplevels(sumbyspp$Alpha[c(which(sumbyspp$N < 3))])
smalist
smallset<-droplevels(subset(birds, Alpha != "BHGR" & Alpha != "FLIN" & Alpha != "GCKI" & Alpha != "HAWO" & Alpha != "HEWA" & Alpha != "PUFI" & Alpha != "RECR" & Alpha != "RSFL" & Alpha != "RUGR" & Alpha != "STJA" & Alpha!= "RUHU" & Alpha != "VASW" & Alpha != "VATH" & Alpha != "WCSP"))
length(levels(smallset$Alpha))
which(smallset$Tarsal.L != "NA")
length(birds$Tarsal.L)


sumsmall<- ddply(smallset, "Alpha", summarise,
                 N = length(Alpha),
                 mean = mean(totgreen),
                 sd = sd(totgreen),
                 se = sd / sqrt(N)
                 )
sumsmall

##################################################
###############Alpha with Tarsal.L################
##################################################
## all models suggest dropping tarsal.l
#################################################

totalalpha<-glm.nb(totgreen~Alpha+Tarsal.L, data=smallset)
summary(totalalpha)
anova(totalalpha, test="Chisq")
pchisq(deviance(totalalpha), df.residual(totalalpha), lower=F) 
totalalpha0<-glm.nb(totgreen~1, data=subset(smallset, Tarsal.L != "NA"))
1-(deviance(totalalpha)/deviance(totalalpha0)) 
anova(totalalpha, totalalpha0, test="Chisq") #significantly better than null

drop1(totalalpha, test="F") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
##drop tarsal length
plot(rstudent(totalalpha)~totalalpha$fitted.values)

legalpha<-glm.nb(leggreen~Alpha+Tarsal.L, data=smallset)
summary(legalpha)
anova(legalpha, test="Chisq")
pchisq(deviance(legalpha), df.residual(legalpha), lower=F) 
legalpha0<-glm.nb(leggreen~1, data=subset(smallset, Tarsal.L != "NA"))
1-(deviance(legalpha)/deviance(legalpha0))
anova(legalpha, legalpha0, test="Chisq") #significantly better than null

drop1(legalpha, test="F") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
#drop tarsal length
plot(rstudent(legalpha)~legalpha$fitted.values)

rectalpha<-glm.nb(rectgreen~Alpha+Tarsal.L, data=smallset)
summary(rectalpha)
anova(rectalpha, test="Chisq")
pchisq(deviance(rectalpha), df.residual(rectalpha), lower=F) 
rectalpha0<-glm.nb(rectgreen~1, data=subset(smallset, Tarsal.L != "NA"))
1-(deviance(rectalpha)/deviance(rectalpha0))
anova(rectalpha, rectalpha0, test="Chisq") #significantly better than null

drop1(rectalpha, test="F") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
#drop tarsal length
plot(rstudent(rectalpha)~rectalpha$fitted.values)

######################################################
###################Alpha##############################
######################################################

totalalpha1<-glm.nb(totgreen~Alpha, data=smallset)
summary(totalalpha1)
anova(totalalpha1, test="Chisq")
pchisq(deviance(totalalpha1), df.residual(totalalpha1), lower=F) 
totalalpha10<-glm.nb(totgreen~1, data=smallset)
1-(deviance(totalalpha1)/deviance(totalalpha10)) 
alphamod1overall<-anova(totalalpha1, totalalpha10, test="Chisq") #significantly better than null

alphamod1parameters<-drop1(totalalpha1, test="LRT") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
plot(rstudent(totalalpha1)~totalalpha1$fitted.values)

summary(totalalpha1)
totalalpha1tidy<-tidy(totalalpha1)
alphamod1parameterstidy<-tidy(alphamod1parameters)
alphamod1overalltidy<-tidy(alphamod1overall)

write.csv(totalalpha1tidy, "GLM summaries/totalalpha.csv")
write.csv(alphamod1parameterstidy, "GLM summaries/totalalphaparameters.csv")
write.csv(alphamod1overalltidy, "GLM summaries/totalalphaoverall.csv")

legalpha1<-glm.nb(leggreen~Alpha, data=smallset)
summary(legalpha1)
anova(legalpha1, test="Chisq")
pchisq(deviance(legalpha1), df.residual(legalpha1), lower=F) 
legalpha10<-glm.nb(leggreen~1, data=smallset)
1-(deviance(legalpha1)/deviance(legalpha10)) 
alpha1legmodoverall<-anova(legalpha1, legalpha10, test="Chisq") #significantly better than null

alpha1legmodparameters<-drop1(legalpha1, test="LRT") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
plot(rstudent(legalpha1)~legalpha1$fitted.values)


legalpha1tidy<-tidy(legalpha1)
alpha1legmodparameterstidy<-tidy(alpha1legmodparameters)
alpha1legmodoveralltidy<-tidy(alpha1legmodoverall)

write.csv(legalpha1tidy, "GLM summaries/legalpha.csv")
write.csv(alpha1legmodparameterstidy, "GLM summaries/legalphaparameters.csv")
write.csv(alpha1legmodoveralltidy, "GLM summaries/legalphaoverall.csv")



rectalpha1<-glm.nb(rectgreen~Alpha, data=smallset)
summary(rectalpha1)
anova(rectalpha1, test="Chisq")
pchisq(deviance(rectalpha1), df.residual(rectalpha1), lower=F) 
rectalpha10<-glm.nb(rectgreen~1, data=smallset)
1-(deviance(rectalpha1)/deviance(rectalpha10)) 
rectalphaoverall<-anova(rectalpha1, rectalpha10, test="Chisq") #significantly better than null

rectalphaparameters<-drop1(rectalpha1, test="LRT") #TWO WAY ANOVA DROPS EACH PREDICTOR SEPARATELY
plot(rstudent(rectalpha1)~rectalpha1$fitted.values)

rectalpha1tidy<-tidy(rectalpha1)
rectalphaoveralltidy<-tidy(rectalphaoverall)
rectalphaparameterstidy<-tidy(rectalphaparameters)

write.csv(rectalpha1tidy, "GLM summaries/rectalpha.csv")
write.csv(rectalphaparameterstidy, "GLM summaries/rectalphaparameters.csv")
write.csv(rectalphaoveralltidy, "GLM summaries/rectalphaoverall.csv")

for(i in levels(birds$Alpha)){
  print(c(i, mean(birds$totgreen[which(birds$Alpha == i)])), sep = ",")
}

nonzero<-droplevels(subset(birds, totgreen != 0))
min(nonzero$totgreen)

for(i in levels(nonzero$Alpha)){
  print(c(i, mean(nonzero$totgreen[which(nonzero$Alpha == i)])), sep = ",")
}

speciesmod<-glm.nb(totgreen~Alpha, data = nonzero)
summary(speciesmod)
anova(speciesmod)
pchisq(deviance(speciesmod), df.residual(speciesmod), lower= F) 
plot(rstudent(speciesmod))


library(multcomp)
library(lsmeans)

###################################
#posthoc pairwise tests
###################################

#ALPHA
alphatotalpairwise<-lsmeans(totalalpha1, pairwise~Alpha, adjust = "tukey")
alphatotalpairwise$contrasts

alphatotalpairs<-cld(alphatotalpairwise, alpha= 0.05, Letters= letters, adjust = "tukey")
alphatotalpairs<-arrange(alphatotalpairs, Alpha)

letterstotalalpha<-alphatotalpairs$.group

alphalegpairwise<-lsmeans(legalpha1, pairwise~Alpha, adjust = "tukey")
alphalegpairwise$contrasts

alphalegpairs<-cld(alphalegpairwise, alpha= 0.05, Letters= letters, adjust = "tukey")
alphalegpairs<-arrange(alphalegpairs, Alpha)

lettersalphaleg<-alphalegpairs$.group

alpharectpairwise<-lsmeans(rectalpha1, pairwise~Alpha, adjust = "tukey")
alpharectpairwise$contrasts

alpharectpairs<-cld(alpharectpairwise, alpha= 0.05, Letters= letters, adjust = "tukey")
alpharectpairs<-arrange(alpharectpairs, Alpha)

lettersalpharect<-alpharectpairs$.group

#total

totalbehaviorpairwise<-lsmeans(birdsnegbi, pairwise~behavior, adjust = "tukey")
totalbehaviorpairwise$contrasts

totalbehaviorpairs<-cld(totalbehaviorpairwise, alpha= 0.05, Letters= letters, adjust = "tukey")
totalbehaviorpairs<-arrange(totalbehaviorpairs, behavior)

letterstotalbehavior<-totalbehaviorpairs$.group

totaltaxonpairwise<-lsmeans(totbytaxon, pairwise~Taxon.group, adjust = "tukey")
totaltaxonpairwise$contrasts

totaltaxonpairs<-cld(totaltaxonpairwise, alpha= 0.05, Letters= letters, adjust = "tukey")
totaltaxonpairs<-arrange(totaltaxonpairs, Taxon.group)

letterstotaltaxon<-totaltaxonpairs$.group

#leg

legbehaviorpairwise<-lsmeans(legsporebytarsalbehavior, pairwise~behavior, adjust = "tukey")
legbehaviorpairwise$contrasts

legbehaviorpairs<-cld(legbehaviorpairwise, alpha= 0.05, Letters= letters, adjust = "tukey")
legbehaviorpairs<-arrange(legbehaviorpairs, behavior)

letterslegbehavior<-legbehaviorpairs$.group

legtaxonpairwise<-lsmeans(legsporebytaxon, pairwise~Taxon.group, adjust = "tukey")
legtaxonpairwise$contrasts

legtaxonpairs<-cld(legtaxonpairwise, alpha= 0.05, Letters= letters, adjust = "tukey")
legtaxonpairs<-arrange(legtaxonpairs, Taxon.group)

letterslegtaxon<-legtaxonpairs$.group

#rect

recttaxonpairwise<-lsmeans(rectsporebytaxon, pairwise~Taxon.group, adjust = "tukey")
recttaxonpairwise$contrasts

recttaxonpairs<-cld(recttaxonpairwise, alpha= 0.05, Letters= letters, adjust = "tukey")
recttaxonpairs<-arrange(recttaxonpairs, Taxon.group)

lettersrecttaxon<-recttaxonpairs$.group
lettersrecttaxon

rectbehaviorpairwise<-lsmeans(rectsporebybehavior, pairwise~behavior, adjust = "tukey")
rectbehaviorpairwise$contrasts

rectbehaviorpairs<-cld(rectbehaviorpairwise, alpha= 0.05, Letters= letters, adjust = "tukey")
rectbehaviorpairs<-arrange(rectbehaviorpairs, behavior)

lettersrectbehavior<-rectbehaviorpairs$.group
lettersrectbehavior

##############
###figures for paper
##############


totalpha<-ggplot(data=subset(smallset, totgreen > 0), aes(Alpha, totgreen)) +
  geom_boxplot(fill="gray") +
  theme_bw()+
  theme(#axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "inch"),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.major.y=element_blank())+
  xlab("")+
  ylab("Total Number of Spores")+
  annotate("text", y = 80, x = c(1,6, 10, 17), label = c("bc","ab", "a", "c"), fontface = "italic", size= 3)
  #annotate("text", y=100, x=1:20, label= letterstotalalpha, fontface = "italic")
#ggsave(path ="../Desktop/WIRI paper 1/WIRI Manuscript 1 - Ecology Letters/plots/", filename = "totALPHAboxfew.png", width=10, height = 5)
totalpha

legalpha<-ggplot(data=subset(smallset, leggreen > 0), aes(Alpha, leggreen)) +
  geom_boxplot(fill="gray") +
  theme_bw()+
  theme(#axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "inch"),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.major.y=element_blank())+
  xlab("")+
  ylab("Number of Leg Spores")
#annotate("text", y=40, x=1:20, label= lettersalphaleg, fontface = "italic")
#ggsave(path ="../Desktop/WIRI paper 1/WIRI Manuscript 1 - Ecology Letters/plots/", filename = "legALPHAboxfew.png", width=10, height = 5)
legalpha


rectalpha<-ggplot(data=subset(smallset, rectgreen > 0), aes(Alpha, rectgreen)) +
  geom_boxplot(fill="gray") +
  theme_bw()+
  theme(#axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "inch"),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.major.y=element_blank())+
  xlab("")+
  ylab("Number of Rectrix Spores")+
  annotate("text", y = 50, x = c(10, 16), label = c("a", "b"), fontface = "italic", size =3)
  #annotate("text", y=50, x=1:20, label= lettersalpharect, fontface = "italic")
#ggsave(path ="../Desktop/WIRI paper 1/WIRI Manuscript 1 - Ecology Letters/plots/", filename = "rectALPHAboxfew.png", width=10, height = 5)
rectalpha

totbehavior<-ggplot(data=subset(birds, behavior != "flower" & totgreen >0), aes(behavior, totgreen)) +
  geom_boxplot(fill="gray") +
  theme_bw()+
  theme(#axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.major.y=element_blank())+
  xlab("")+
  ylab("Total Number of Spores")+
  annotate("text", x=1:4, y=70, label= letterstotalbehavior, fontface = "italic", size =3)
#ggsave(path = "../Desktop/WIRI paper 1/WIRI Manuscript 1 - Ecology Letters/plots/", filename = "totbehaviorboxfew.png", width=10, height = 5)
totbehavior


legbehavior<-ggplot(data=subset(birds, behavior != "flower" & leggreen > 0), aes(behavior, leggreen)) +
  geom_boxplot(fill="gray") +
  theme_bw()+
  theme(#axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.major.y=element_blank())+
  xlab("")+
  ylab("Number of Leg Spores")
#annotate("text", x=1:4, y=10, label= letterslegbehavior)
#ggsave(path ="../Desktop/WIRI paper 1/WIRI Manuscript 1 - Ecology Letters/plots/", filename = "legbehaviorboxfew.png", width=10, height = 5)
legbehavior


rectbehavior<-ggplot(data=subset(birds, behavior != "flower" & rectgreen > 0), aes(behavior, rectgreen)) +
  geom_boxplot(fill="gray") +
  theme_bw()+
  theme(#axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.major.y=element_blank())+
  xlab("")+
  ylab("Number of Rectrix Spores")+
  annotate("text", x=1:4, y=80, label= lettersrectbehavior, fontface = "italic", size =3)
#ggsave(path ="../Desktop/WIRI paper 1/WIRI Manuscript 1 - Ecology Letters/plots/", filename = "rectbehaviorboxfew.png", width=10, height = 5)
rectbehavior

tottaxon<-ggplot(data=subset(birds, totgreen >0), aes(Taxon.group, totgreen)) +
  geom_boxplot(fill="gray") +
  theme_bw()+
  theme(#axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.major.y=element_blank())+
  xlab("")+
  ylab("Total Number of Spores")+
  annotate("text", x=c(2,4,5,7,8,9,11,12), y=85, label= c("ab","bc","bc","bc","bc","a","c","ab"), fontface = "italic", size =3)
#ggsave(path ="../Desktop/WIRI paper 1/WIRI Manuscript 1 - Ecology Letters/plots/", filename = "tottaxonboxfew.png", width=10, height = 5)
tottaxon

totaltaxonpairs

legtaxon<-ggplot(data=subset(birds, leggreen > 0), aes(Taxon.group, leggreen)) +
  geom_boxplot(fill="gray") +
  theme_bw()+
  theme(#axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.major.y=element_blank())+
  xlab("")+
  ylab("Number of Leg Spores")
  #annotate("text", x=1:16, y=80, label= letterslegtaxon, fontface = "italic", size =3)
#ggsave(path ="../Desktop/WIRI paper 1/WIRI Manuscript 1 - Ecology Letters/plots/", filename = "legtaxonboxfew.png", width=10, height = 5)
legtaxon

recttaxon<-ggplot(data=subset(birds, rectgreen > 0), aes(Taxon.group, rectgreen)) +
  geom_boxplot(fill="gray") +
  theme_bw()+
  theme(#axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.major.y=element_blank())+
  xlab("")+
  ylab("Number of Rectrix Spores")+
  annotate("text", y=50, x=c(2,10,11), label= c("a","b","a"), fontface = "italic", size =3)
#ggsave(path ="../Desktop/WIRI paper 1/WIRI Manuscript 1 - Ecology Letters/plots/", filename = "recttaxonboxfew.png", width=10, height = 5)
recttaxon

recttaxonpairs

library(cowplot)

totalspores<-plot_grid(totalpha, tottaxon, totbehavior, labels = c("A", "B", "C"), ncol =2, nrow =2, align="h")
legspores<-plot_grid(legalpha, legtaxon, legbehavior, labels = c("A", "B", "C"), ncol =2, nrow =2, align="h")
rectspores<-plot_grid(rectalpha, recttaxon, rectbehavior, labels = c("A", "B", "C"), ncol =2, nrow =2, align="h")

save_plot("plots/totalmulti.pdf", totalspores, ncol=2, nrow=2)
save_plot("plots/legmulti.pdf", legspores, ncol=2, nrow=2)
save_plot("plots/rectmulti.pdf", rectspores, ncol=2, nrow=2)

#pdf generation of full size plots with no margin
pdf("plots/totalmultifullsize.pdf")
plot(totalspores)
dev.off()

pdf("plots/legmultifullsize.pdf")
plot(legspores)
dev.off()

pdf("plots/rectmultifullsize.pdf")
plot(rectspores)
dev.off()

#pdfs of plots with 1 inch margins
pdf("plots/totalmultimargins.pdf", paper="a4r")
par(omi=c(1,1,1,1))
plot(totalspores)
dev.off()

pdf("plots/legmultimargins.pdf", paper="a4r")
par(omi=c(1,1,1,1))
plot(legspores)
dev.off()

pdf("plots/rectmultimargins.pdf", paper="a4r")
par(omi=c(1,1,1,1))
plot(rectspores)
dev.off()

#png of plots for presentation
png("plots/totalmultifullsize.png", width=1150, height=1150, res=150)
plot(totalspores)
dev.off()

png("plots/legmultifullsize.png", width=1150, height=1150, res=150)
plot(legspores)
dev.off()

png("plots/rectmultifullsize.png", width=1150, height=1150, res=150)
plot(rectspores)
dev.off()


tarsaltaxon<-ggplot(data=subset(birds, Tarsal.L != "NA"), aes(Taxon.group, Tarsal.L)) +
  geom_boxplot() +
  theme_bw()+
  theme(#axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.major.y=element_blank())+
  xlab("Taxon")+
  ylab("Leg length")

legtaxon
plot_grid(legtaxon, tarsaltaxon, nrow=1, ncol=2 )

######## Total spores by tarsal length

detach("package:cowplot", unload=TRUE)

testing<-glm.nb(totgreen~Tarsal.L, data=subset(birds, Tarsal.L != "NA" & totgreen > 0))
summary(testing)

figure1legall<-ggplot(data=subset(birds, Tarsal.L != "NA"), aes(Tarsal.L, totgreen)) +
  geom_point() +
  scale_x_continuous(breaks = seq(min(0), max(50), by = 5)) +
  scale_y_continuous(breaks = seq(min(0), max(100), by = 10)) +
  geom_smooth(method="glm.nb", color="black")+
  theme_classic()+
  xlab("Tarsal Length (mm)")+
  ylab("Total Spore Count")

figure1legall

#fullsizenomargins
pdf("plots/Figure1full.pdf")
plot(figure1legall)
dev.off()

#with one inch margins
pdf("plots/Figure1margins.pdf", paper="a4r")
par(omi=c(1,1,1,1))
plot(figure1legall)
dev.off()

#png file for presentation
png("plots/Figure1full.png", width=1150, height=1150, res=150)
plot(figure1legall)
dev.off()

