library(metaplus)

# speed execution this uses just 99 bootstrap replications
# for final results this should be changed to 999
nosims <- 99

data(mag)
mag1 <- metaplus(mag$yi,mag$sei,slab=mag$study)

summary(mag1)

plot(mag1)

# can modify plot to obtain Odds Ratios
plot(mag1,atransf=exp, at=log(c(.01,.1,1,10,100)),xlab="Odds Ratio")

mag2 <- metaplus(mag$yi,mag$sei,slab=mag$study,random="t-dist")

summary(mag2)

summary(test.outliers(mag2,nsim=nosims))

mag3 <- metaplus(mag$yi,mag$sei,slab=mag$study,random="mixture")
 
summary(mag3)

summary(test.outliers(mag3,nsim=nosims))

data(cdp)

cdp1 <- metaplus(cdp$yi,cdp$sei,plotci=TRUE,slab=cdp$study)
summary(cdp1)

plot(cdp1)

cdp2 <- metaplus(cdp$yi,cdp$sei,plotci=TRUE,slab=cdp$study,random="t-dist")
summary(cdp2)

summary(test.outliers(cdp2,nsim=nosims))

cdp3 <- metaplus(cdp$yi,cdp$sei,plotci=TRUE,slab=cdp$study,random="mixture")
summary(cdp3)

cdp3.outlier.probs <- outlier.probs(cdp3)
plot(cdp3.outlier.probs)

summary(test.outliers(cdp3,nsim=nosims))

plot(cdp1,extrameta=list(cdp2,cdp3))

data(marinho)

marinho1 <- metaplus(marinho$meaneffect,marinho$seeffect,plotci=TRUE,slab=marinho$study)
summary(marinho1)

marinho2 <- metaplus(marinho$meaneffect,marinho$seeffect,plotci=TRUE,slab=marinho$study,random="t-dist")
summary(marinho2)

marinho3 <- metaplus(marinho$meaneffect,marinho$seeffect,plotci=TRUE,slab=marinho$study,random="mixture")
summary(marinho3)

plot(marinho1,extrameta=list(marinho2,marinho3))

marinho3.outlier.probs <- outlier.probs(marinho3)
plot(marinho3.outlier.probs)

data(exercise)

exercise1 <- metaplus(exercise$smd,sqrt(exercise$varsmd),mods=exercise[,c("duration"),drop=FALSE],plotci=TRUE,slab=exercise$study)
summary(exercise1)

exercise2 <- metaplus(exercise$smd,sqrt(exercise$varsmd),mods=exercise[,c("duration"),drop=FALSE],plotci=TRUE,random="t-dist")
summary(exercise2)

summary(test.outliers(exercise2,nsim=nosims))

exercise3 <- metaplus(exercise$smd,sqrt(exercise$varsmd),mods=exercise[,c("duration"),drop=FALSE],plotci=TRUE,random="mixture")
summary(exercise3)

summary(test.outliers(exercise3,nsim=nosims))

exercise3.outlier.probs <- outlier.probs(exercise3)
plot(exercise3.outlier.probs)

exercise$duration4 <- exercise$duration-4
exercise$duration8 <- exercise$duration-8
exercise$duration12 <- exercise$duration-12

exercise.wk4 <- metaplus(exercise$smd,sqrt(exercise$varsmd),mods=exercise[,c("duration4"),drop=FALSE],plotci=TRUE,
                         label="Random Mixture (Week 4)",slab=exercise$study,random="mixture")
exercise.wk8 <- metaplus(exercise$smd,sqrt(exercise$varsmd),mods=exercise[,c("duration8"),drop=FALSE],plotci=TRUE,
                         label="Random Mixture (Week 8)",slab=exercise$study,random="mixture")
exercise.wk12 <- metaplus(exercise$smd,sqrt(exercise$varsmd),mods=exercise[,c("duration12"),drop=FALSE],plotci=TRUE,
                          label="Random Mixture (Week 12)",slab=exercise$study,random="mixture")

exercise.nodurn <- metaplus(exercise$smd,sqrt(exercise$varsmd),plotci=TRUE,
                          label="Random Mixture (No Duration)",slab=exercise$study,random="mixture")

plot(exercise.nodurn,extrameta=list(exercise.wk4,exercise.wk8,exercise.wk12))
