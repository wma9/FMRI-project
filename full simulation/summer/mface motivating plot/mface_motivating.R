# to get rid of the AsIs attributes
library(refund)
data(DTI2)
unAsIs <- function(X) {
  if("AsIs" %in% class(X)) {
    class(X) <- class(X)[-match("AsIs", class(X))]
  }
  X
}
# 100 subjects, visits 1~8
# get rid of the "AsIs" format

uncca <- unAsIs(DTI2$cca)
unrcst <- unAsIs(DTI2$rcst)
# plot
par(mfrow=c(1,1), mar=c(2.8,2.5,1,1), mgp = c(1.4,0.5,0))
matplot(1:dim(uncca)[2], t(uncca), col="darkgrey", type="l", main="CCA",
        xaxt='n', xlab="Profile", ylab="Value")
lines(uncca[which(DTI2$visit==4&DTI2$id==20001),], col="red", lwd=1.5)# first subject, 1~2 visits
lines(uncca[which(DTI2$visit==5&DTI2$id==20001),], col="red", lwd=1.5)
lines(uncca[which(DTI2$visit==2&DTI2$id==20003),], col="blue", lwd=1.5)# last subject, 1~2 visits
lines(uncca[which(DTI2$visit==3&DTI2$id==20003),], col="blue", lwd=1.5)
axis(1, at=c(1,19,37,56,74,93),labels=c(0,0.2,0.4,0.6,0.8,1), cex.lab=0.8, cex.axis=0.8)

matplot(1:dim(unrcst)[2], t(unrcst), col="darkgrey", type="l", main="RCST",
        xaxt='n', xlab="Profile", ylab="Value")
lines(unrcst[which(DTI2$visit==4&DTI2$id==20001),], col="red", lwd=1.5)
lines(unrcst[which(DTI2$visit==5&DTI2$id==20001),], col="red", lwd=1.5)
lines(unrcst[which(DTI2$visit==2&DTI2$id==20003),], col="blue", lwd=1.5)
lines(unrcst[which(DTI2$visit==3&DTI2$id==20003),], col="blue", lwd=1.5)
axis(1, at=c(1,11,22,33,44,55),labels=c(0,0.2,0.4,0.6,0.8,1), cex.lab=0.8, cex.axis=0.8)

# fpca
results.rcst <- fpca.face(unrcst, pve = 0.95, center = TRUE) # pve need to be chosen!
results.cca <- fpca.face(uncca, pve = 0.95, center = TRUE) # pve need to be chosen!
matplot(1:dim(unrcst)[2], (results.rcst$efunctions*sqrt(dim(unrcst)[2]))[,1:3], type="l", 
        col=c(1,2,3), lty=c(1,2,3),lwd=1.5,
        xaxt='n', xlab="Profile", ylab="",main="Univariate fPCA on RCST",
        ylim=c(min((results.rcst$efunctions*sqrt(dim(unrcst)[2]))[,1:3], 
                 (results.cca$efunctions*sqrt(dim(uncca)[2]))[,1:3]),
        max((results.rcst$efunctions*sqrt(dim(unrcst)[2]))[,1:3], 
            (results.cca$efunctions*sqrt(dim(uncca)[2]))[,1:3])))
legend("topright", c("1th PC", "2nd PC", "3rd PC"), col=c(1,2,3), lty=c(1,2,3), lwd=c(1.5,1.5,1.5), bty='n', cex=0.8)
axis(1, at=c(1,11,22,33,44,55),labels=c(0,0.2,0.4,0.6,0.8,1), cex.lab=0.8, cex.axis=0.8)

  
matplot(1:dim(uncca)[2], (results.cca$efunctions*sqrt(dim(uncca)[2]))[,1:3], type="l", 
        col=c(1,2,3),lty=c(1,2,3),lwd=1.5,
        xaxt='n', xlab="Profile", ylab="", main="Univariate fPCA on CCA",
        ylim=c(min((results.rcst$efunctions*sqrt(dim(unrcst)[2]))[,1:3], 
                   (results.cca$efunctions*sqrt(dim(uncca)[2]))[,1:3]),
               max((results.rcst$efunctions*sqrt(dim(unrcst)[2]))[,1:3], 
                   (results.cca$efunctions*sqrt(dim(uncca)[2]))[,1:3])))
legend("topright", c("1th PC", "2nd PC", "3rd PC"), col=c(1,2,3), lty=c(1,2,3),  lwd=c(1.5,1.5,1.5),bty='n', cex=0.8)
axis(1, at=c(1,19,37,56,74,93),labels=c(0,0.2,0.4,0.6,0.8,1), cex.lab=0.8, cex.axis=0.8)

# score correlation
scores.cca <- results.cca$scores[, 1:3]/sqrt(dim(uncca)[2])
scores.rcst <- results.rcst$scores[, 1:3]/sqrt(dim(unrcst)[2])
COR <- cor(scores.cca,scores.rcst)
image(x=seq(dim(scores.cca)[2]), y=seq(dim(scores.rcst)[2]), z=COR, xlab="CCA", ylab="RCST")
text(expand.grid(x=seq(dim(scores.cca)[2]), y=seq(dim(scores.rcst)[2])), labels=round(c(COR),2))
