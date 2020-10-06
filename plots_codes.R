# Nucleosides project#

# Figures #

# Fig.1 #
#Residuals #
pseu <- logCobsPred[,1:178]; p <- as.data.frame(colMeans(pseu))
X1mA <- logCobsPred[,179:356]; p2 <- as.data.frame(colMeans(X1mA))
X3mC <- logCobsPred[,357:534]; p3 <- as.data.frame(colMeans(X3mC))
X1 <- logCobsPred[,535:712]; p4 <- as.data.frame(colMeans(X1))
A <- logCobsPred[,713:890]; p5 <- as.data.frame(colMeans(A))
N4ac <- logCobsPred[,891:1068]; p6 <- as.data.frame(colMeans(N4ac))
N2N2 <- logCobsPred[,1069:1246]; p7 <- as.data.frame(colMeans(N2N2))
MTA <- logCobsPred[,1247:1424]; p8 <- as.data.frame(colMeans(MTA))
X7mG <- logCobsPred[,1425:1602]; p9 <- as.data.frame(colMeans(X7mG))
I <- logCobsPred[,1603:1780]; p10 <- as.data.frame(colMeans(I))
G <- logCobsPred[,1781:1958]; p11 <- as.data.frame(colMeans(G))
X5mU <- logCobsPred[,1959:2136]; p12 <- as.data.frame(colMeans(X5mU))
X6mA <- logCobsPred[,2137:2314]; p13 <- as.data.frame(colMeans(X6mA))

xx = cbind(p,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13)

residuals = as.data.frame(logCobs) - xx

tiff("residuals.tiff", width = 4000, height = 4000, units = 'px', res = 300)

par(mfrow=c(4,4))
plot(tage, residuals$pseu, xlab = "Age", ylab="Residuals",ylim=c(-4,4),cex.main = 1.8, cex.sub = 1.5, cex.lab = 1.5,cex.sub = 1.5,cex=2,las=1, pch=16, main="Pseu",col=c('blue', 'red')); abline(h=0, col=17)
plot(tage, residuals$X1mA, xlab = "Age", ylab="Residuals",ylim=c(-4,4),cex.main = 1.8, cex.sub = 1.5, cex.lab = 1.5,cex.sub = 1.5,cex=2, las=1, pch=16,main="1mA",col=c('blue', 'red')); abline(h=0, col=17)
plot(tage, residuals$X3mC, xlab = "Age", ylab="Residuals",ylim=c(-4,4),cex.main = 1.8, cex.sub = 1.5, cex.lab = 1.5,cex.sub = 1.5,cex=2, las=1, pch=16,main="3mC",col=c('blue', 'red'));abline(h=0, col=17)
plot(tage, residuals$X, xlab = "Age", ylab="Residuals",ylim=c(-4,4),cex.main = 1.8, cex.sub = 1.5, cex.lab = 1.5,cex.sub = 1.5,cex=2, las=1, pch=16,main="X",col=c('blue', 'red'));abline(h=0, col=17)
plot(tage, residuals$A, xlab = "Age", ylab="Residuals",ylim=c(-4,4),cex.main = 1.8, cex.sub = 1.5, cex.lab = 1.5,cex.sub = 1.5,cex=2,las=1, pch=16,main="A",col=c('blue', 'red'));abline(h=0, col=17)
plot(tage, residuals$N4ac, xlab = "Age", ylab="Residuals",ylim=c(-4,4),cex.main = 1.8, cex.sub = 1.5, cex.lab = 1.5,cex.sub = 1.5,cex=2, las=1, pch=16,main="N4ac",col=c('blue', 'red'));abline(h=0, col=17)
plot(tage, residuals$N2N2, xlab = "Age", ylab="Residuals",ylim=c(-4,4),cex.main = 1.8, cex.sub = 1.5, cex.lab = 1.5,cex.sub = 1.5,cex=2, las=1, pch=16,main="N2N2",col=c('blue', 'red'));abline(h=0, col=17)
plot(tage, residuals$MTA, xlab = "Age", ylab="Residuals", ylim=c(-4,4),cex.main = 1.8, cex.sub = 1.5, cex.lab = 1.5,cex.sub = 1.5,cex=2, las=1, pch=16,main="MTA",col=c('blue', 'red'));abline(h=0, col=17)
plot(tage, residuals$X7mG, xlab = "Age", ylab="Residuals", ylim=c(-4,4),cex.main = 1.8, cex.sub = 1.5, cex.lab = 1.5,cex.sub = 1.5,cex=2,las=1, pch=16,main="7mG",col=c('blue', 'red'));abline(h=0, col=17)
plot(tage, residuals$I, xlab = "Age", ylab="Residuals", ylim=c(-4,4),cex.main = 1.8, cex.sub = 1.5, cex.lab = 1.5,cex.sub = 1.5,cex=2,las=1, pch=16,main="I",col=c('blue', 'red'));abline(h=0,col=17)
plot(tage, residuals$G, xlab = "Age", ylab="Residuals", ylim=c(-4,4),cex.main = 1.8, cex.sub = 1.5, cex.lab = 1.5,cex.sub = 1.5,cex=2,las=1, pch=16,main="G",col=c('blue', 'red'));abline(h=0, col=17)
plot(tage, residuals$X5mU, xlab = "Age", ylab="Residuals",ylim=c(-4,4),cex.main = 1.8, cex.sub = 1.5, cex.lab = 1.5,cex.sub = 1.5,cex=2,las=1, pch=16,main="5mU",col=c('blue', 'red'));abline(h=0, col=17)
plot(tage, residuals$X6mA, xlab = "Age", ylab="Residuals",ylim=c(-4,4), cex.main = 1.8, cex.sub = 1.5, cex.lab = 1.5,cex.sub = 1.5,cex=2,las=1, pch=16,main="6mA",col=c('blue', 'red'));abline(h=0, col=17)


# Fig.2 #
# Visual predictive check (VPC) -actual data versus data predicted by the model

## ACTUAL DATA ###
# Actual
t <- logCobs 
actual <- as.numeric(t[,1])# 
med <- median(actual)#
hist(actual, main="Pseu", xlab="",cex.main = 2, font.main= 2, col.main= "black");
lines( c(med,med), c(0,70), col = "black", lwd = 2) #median
a.5 <- quantile(actual, probs = c(0.05))
lines( c(a.5,a.5), c(0,70), col = "red", lwd = 2) #5th
a.95 <- quantile(actual, probs = c(0.95))
lines( c(a.95,a.95), c(0,70), col = "blue", lwd = 2) #5th

#Predicted 
pseu <- logCobsPred[,1:178];
pseu5 <- as.data.frame(apply(pseu, c(1), function(x) quantile(x, probs = c(0.05)))); dim(pseu5)
colnames(pseu5) <- "rr"
a.5 <- quantile(pseu5$rr, probs = c(0.05))
lines( c(a.5,a.5), c(0,70), col = "red", lwd = 2, lty=2) #5th
a.95 <- quantile(pseu5$rr, probs = c(0.95))
lines( c(a.95,a.95), c(0,70), col = "red", lwd = 2, lty=2) #5th

pseu95 <- as.data.frame(apply(pseu, c(1), function(x) quantile(x, probs = c(0.95)))); dim(pseu95)
colnames(pseu95) <- "rr"
a.5 <- quantile(pseu95$rr, probs = c(0.05))
lines( c(a.5,a.5), c(0,70), col = "blue", lwd = 2, lty=2) #5th
a.95 <- quantile(pseu95$rr, probs = c(0.95))
lines( c(a.95,a.95), c(0,70), col = "blue", lwd = 2, lty=2) #5th

pseu50 <- as.data.frame(apply(pseu, c(1), function(x) quantile(x, probs = c(0.5)))); dim(pseu50)
colnames(pseu50) <- "rr"
a.5 <- quantile(pseu50$rr, probs = c(0.05))
lines( c(a.5,a.5), c(0,70), col = "black", lwd = 2, lty=2) #5th
a.95 <- quantile(pseu50$rr, probs = c(0.95))
lines( c(a.95,a.95), c(0,70), col = "black", lwd = 2, lty=2) #


#Fig.3 #
# Prior effect of covariates #

par(mfrow=c(2,2))

#########
# Theta #
#########
par(mar=c(5.1,5.1,4.1,2.1))
colnames(theta1) = c('Pseu',	'1mA'	,'3mC',	'X','A'	,'N4Ac',	'N2N2','MTA','7mG','I','G',	'5mU','6mA');
boxplot(lapply(as.data.frame((theta1)), quantile, probs = c(0.95, 0.5, 0.05)), horizontal = FALSE, outline = FALSE, medcol = "red", whisklty = 1,staplelwd = 3, outpch = 8,
        outcex = .7, cex.lab=1.5, lwd =0.5, cex.axis = 1.4,cex.main = 1.6, main ="Typical nucleoside/creatinine excretion rates",mgp=c(1,1,0), las=2)
title(ylab = expression(paste(" (",theta[j],")"),c(5, 8, 9, 5) + 0.1 ,cex.lab=1))
abline(h = 0, col = "darkgray", lwd=1) 


###################
# Casecontrol-frc #
###################
colnames(FRC1) = c('Pseu',	'1mA'	,'3mC',	'X','A'	,'N4Ac',	'N2N2','MTA','7mG','I','G',	'5mU','6mA');
boxplot(lapply(exp(as.data.frame(FRC1)), quantile, probs = c(0.95, 0.5, 0.05)),horizontal = FALSE, outline = FALSE, ylim=c(0.5,2.5),medcol = "red", whisklty = 1,staplelwd = 3, outpch = 8,
        outcex = .7, cex.lab=1.5, lwd =0.5, cex.axis = 1.4,cex.main = 1.6, main ="Case/control effect",mgp=c(1,1,0), las=2)
title(ylab = expression(paste("exp (",  beta[case/control], ")"), c(5, 8, 9, 5) + 0.1 ,cex.lab=3))
abline(h = 1, col = "darkgray", lwd=1) 


#######
# FRA # age effect
#######
colnames(FRA1) = c('Pseu',	'1mA'	,'3mC',	'X','A'	,'N4Ac',	'N2N2','MTA','7mG','I','G',	'5mU','6mA');
boxplot(lapply(exp(as.data.frame(FRA1)), quantile, probs = c(0.95, 0.5, 0.05)),horizontal = FALSE, outline = FALSE, ylim=c(0.5,2.5),medcol = "red", whisklty = 1,staplelwd = 3, outpch = 8,
        outcex = .7, cex.lab=1.5, lwd =0.5, cex.axis = 1.4,cex.main = 1.6, main ="Age effect",mgp=c(1,1,0), las=2)
title(ylab = expression(paste("exp (",  beta[age], ")"), c(5, 8, 9, 5) + 0.1 ,cex.lab=3))
abline(h = 1.22, col = "darkgray", lwd=1) 
abline(h = 1.0, col = "black", lwd=1)

#######
# FRG # sex effect
########
colnames(FRG1) = c('Pseu',	'1mA'	,'3mC',	'X','A'	,'N4Ac',	'N2N2','MTA','7mG','I','G',	'5mU','6mA');
boxplot(lapply(exp(as.data.frame(FRG1)), quantile, probs = c(0.95, 0.5, 0.05)),horizontal = FALSE, outline = FALSE, ylim=c(0.5,2.5),medcol = "red", whisklty = 1,staplelwd = 3, outpch = 8,
        outcex = .7, cex.lab=1.5, lwd =0.5, cex.axis = 1.4,cex.main = 1.6, main ="Sex effect",mgp=c(1,1,0), las=2)
title(ylab = expression(paste("exp (",  beta[sex], ")"),c(5, 8, 9, 5) + 0.1 ,cex.lab=3))
abline(h = 1.35, col = "darkgray", lwd=1)
abline(h = 1.0, col = "black", lwd=1)


dev.off()

#Fig.4#

###Correlation matrix ##

pseu <- omega[,1:13];p1 <- as.matrix(rowMeans(pseu)); dim(p1)
X1mA <- omega[,14:26];p2 <- as.matrix(rowMeans(X1mA)); dim(X1mA)
X3mC <- omega[,27:39];p3 <- as.matrix(rowMeans(X3mC)); dim(X3mC)
X1 <- omega[,40:52];p4 <- as.matrix(rowMeans(X1));dim(X1)
A <- omega[,53:65];p5 <- as.matrix(rowMeans(A));dim(A)
N4ac <- omega[,66:78];p6 <- as.matrix(rowMeans(N4ac ));dim(N4ac)
N2N2 <- omega[,79:91];p7 <- as.matrix(rowMeans(N2N2));dim(N2N2)
MTA <- omega[,92:104];p8 <- as.matrix(rowMeans(MTA));dim(MTA)
X7mG <- omega[,105:117];p9 <- as.matrix(rowMeans(X7mG));dim(X7mG)
I <- omega[,118:130];p10 <- as.matrix(rowMeans(I));dim(I)
G <- omega[,131:143];p11 <- as.matrix(rowMeans(G));dim(G)
X5mU <- omega[,144:156];p12 <- as.matrix(rowMeans(X5mU));dim(X5mU)
X6mA <- omega[,157:169];p13 <- as.matrix(rowMeans(X6mA));dim(X6mA)

cor.mat <- cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13); dim(cor.mat)
colnames(cor.mat) <-  c('Pseu',	'1mA'	,'3mC',	'X','A'	,'N4Ac',	'N2N2','MTA','7mG','I','G',	'5mU','6mA')
cor.plot <- rcorr(cor.mat);

library(RColorBrewer);library("corrplot")

setwd("D:/Bayesian/Bayesian stat_nucleosides/codes/modelowanie stezen_22092016/prior_onAgeSexCaCo/ploty")
tiff("corr.matrix.tiff",width = 4000, height = 4000, units = 'px', res = 300)

cor.mat1 <-cor(cor.mat)
corrplot(cor.mat1,type="upper", order="hclust", cl.lim=c(0,1),  
         col=brewer.pal(n=7, name="PuOr"),cl.length=20,addCoef.col="grey",pch=.5,pch.cex = 0.5,
         insig="pch", tl.col = "black",tl.pos="td",tl.srt = 22); 
dev.off()

#Fig.5 #
## Individual probability of disease #

tiff("Ind_prob.tiff", width = 4000, height = 4000, units = "px", res = 300)

par(mfrow=c(3,4))


plot (theta,dbeta(theta,a1,b1), type ="l", lty=2, ylab="",ylim=c(0,2.5), col=2, lwd=2, main="Cancer",
      xlab="Individual 1",axes=T); #PRIOR
par(new=T)
hist(ppp[,1], prob=TRUE, main = "", xlab="Individual 1",ylim=c(0,2.5), axes=F);
par(new=T);
abline(v = 0.5, col= "darkblue", lwd=3); 
par(new=T)
abline(v=mean(ppp[,1]), col=33, lwd=2, lty=2)

plot (theta,dbeta(theta,a1,b1), type ="l", lty=2,ylab="",ylim=c(0,2.5), col=2, lwd=2, main="Cancer",
      xlab="Individual 5") #PRIOR
par(new=T)
hist(ppp[,5], prob=TRUE, main = "", xlab="Individual 5",ylim=c(0,2.5))
par(new=T);
abline(v = 0.5, col= "darkblue", lwd=3)
par(new=T)
abline(v=mean(ppp[,5]), col=33, lwd=2, lty=2)

plot (theta,dbeta(theta,a1,b1), type ="l", lty=2,ylab="",ylim=c(0,2.5), col=2, lwd=2, main="Cancer",
      xlab="Individual 10") #PRIOR
par(new=T)
hist(ppp[,10], prob=TRUE, main = "", xlab="Individual 10",ylim=c(0,2.5))
par(new=T);
abline(v = 0.5, col= "darkblue", lwd=3)
par(new=T)
abline(v=mean(ppp[,10]), col=33, lwd=2, lty=2)

plot (theta,dbeta(theta,a1,b1), type ="l", lty=2,ylab="",ylim=c(0,2.5), col=2, lwd=2, main="Healthy",
      xlab="Individual 13") #PRIOR
par(new=T)
hist(ppp[,13], prob=TRUE, main = "", xlab="Individual 13",ylim=c(0,2.5))
par(new=T);
abline(v = 0.5, col= "darkblue", lwd=3)
par(new=T)
abline(v=mean(ppp[,13]), col=33, lwd=2, lty=2)


plot (theta,dbeta(theta,a1,b1), type ="l", lty=2,ylab="",ylim=c(0,2.5), col=2, lwd=2, main="Cancer",
      xlab="Individual 30") #PRIOR
par(new=T)
hist(ppp[,30], prob=TRUE, main = "", xlab="Individual 30",ylim=c(0,2.5))
par(new=T);
abline(v = 0.5, col= "darkblue", lwd=3)
par(new=T)
abline(v=mean(ppp[,30]), col=33, lwd=2, lty=2)

plot (theta,dbeta(theta,a1,b1), type ="l", lty=2,ylab="",ylim=c(0,2.5), col=2, lwd=2, main="Healthy",
      xlab="Individual 54") #PRIOR
par(new=T)
hist(ppp[,54], prob=TRUE, main = "", xlab="Individual 54",ylim=c(0,2.5))
par(new=T);
abline(v = 0.5, col= "darkblue", lwd=3)
par(new=T)
abline(v=mean(ppp[,54]), col=33, lwd=2, lty=2)


plot (theta,dbeta(theta,a1,b1), type ="l", lty=2,ylab="",ylim=c(0,2.5), col=2, lwd=2, main="Healthy",
      xlab="Individual 59") #PRIOR
par(new=T)
hist(ppp[,59], prob=TRUE, main = "", xlab="Individual 59",ylim=c(0,2.5))
par(new=T);
abline(v = 0.5, col= "darkblue", lwd=3)
par(new=T)
abline(v=mean(ppp[,59]), col=33, lwd=2, lty=2)


plot (theta,dbeta(theta,a1,b1), type ="l", lty=2,ylab="",ylim=c(0,2.5), col=2, lwd=2, main="Cancer",
      xlab="Individual 63") #PRIOR
par(new=T)
hist(ppp[,63], prob=TRUE, main = "", xlab="Individual 63",ylim=c(0,2.5))
par(new=T);
abline(v = 0.5, col= "darkblue", lwd=3)
par(new=T)
abline(v=mean(ppp[,63]), col=33, lwd=2, lty=2)


plot (theta,dbeta(theta,a1,b1), type ="l", lty=2,ylab="",ylim=c(0,2.5), col=2, lwd=2, main="Cancer",
      xlab="Individual 68") #PRIOR
par(new=T)
hist(ppp[,68], prob=TRUE, main = "", xlab="Individual 68",ylim=c(0,2.5))
par(new=T);
abline(v = 0.5, col= "darkblue", lwd=3)
par(new=T)
abline(v=mean(ppp[,68]), col=33, lwd=2, lty=2)


plot (theta,dbeta(theta,a1,b1), type ="l", lty=2,ylab="",ylim=c(0,2.5), col=2, lwd=2, main="Healthy",
      xlab="Individual 38") #PRIOR
par(new=T)
hist(ppp[,38], prob=TRUE, main = "", xlab="Individual 38",ylim=c(0,2.5))
par(new=T);
abline(v = 0.5, col= "darkblue", lwd=3)
par(new=T)
abline(v=mean(ppp[,38]), col=33, lwd=2, lty=2)


plot (theta,dbeta(theta,a1,b1), type ="l", lty=2,ylab="",ylim=c(0,2.5), col=2, lwd=2, main="Healthy",
      xlab="Individual 29") #PRIOR
par(new=T)
hist(ppp[,29], prob=TRUE, main = "", xlab="Individual 29",ylim=c(0,2.5))
par(new=T);
abline(v = 0.5, col= "darkblue", lwd=3)
par(new=T)
abline(v=mean(ppp[,29]), col=33, lwd=2, lty=2)


plot (theta,dbeta(theta,a1,b1), type ="l", lty=2,ylab="",ylim=c(0,2.5), col=2, lwd=2, main="Cancer",
      xlab="Individual 33") #PRIOR
par(new=T)
hist(ppp[,33], prob=TRUE, main = "", xlab="Individual 33",ylim=c(0,2.5))
par(new=T);
abline(v = 0.5, col= "darkblue", lwd=3)
par(new=T)
abline(v=mean(ppp[,33]), col=33, lwd=2, lty=2)
dev.off()
