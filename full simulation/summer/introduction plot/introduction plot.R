load("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/03.13.2018/real data application/realdata.Rdata")
# 1th voxel 
voxel=1
fmridata <- (center_whole_fmri)[(1+943*(voxel-1)):(943*voxel),]
sub1 <- dplyr::filter(fmridata, ID==10)[1:3,]
sub2 <- dplyr::filter(fmridata, ID==20)[1:3,]

par(mfrow=c(1,3), mar=c(2.8,2.5,1,1), mgp = c(1.4,0.5,0))
for(i in 1:3){
  plot(1:23,sub1[i,4:26], type="l", xlab="Time (Seconds)", ylab="fMRI time series", xaxt='n', ylim = c(-11,11))
  axis(1, at=c(5,10,15,20),labels=c(10,20,30,40))
}

par(mfrow=c(1,3), mar=c(3.4,2.6,0.3,1), mgp = c(1.4,0.5,0))
for(i in 1:3){
  plot(1:23,sub2[i,4:26], type="l", xlab="Time (Seconds)", ylab="fMRI time series", xaxt='n', ylim = c(-5,5), lty=2)
  axis(1, at=c(5,10,15,20),labels=c(10,20,30,40))
  if(i==2){
    mtext("(c)",side=1, cex=0.8, outer=F, line= 2.5)
  }
}

par(mfrow=c(1,1), mar=c(2.2,2.4,1.1,1), mgp = c(1.3,0.4,0))
plot(2:4, sub1$rating, type="o", pch = 16, cex=0.8, axes = F, xaxt='n',
     xlim=c(1.9,4.15), ylim=c(100,600),
     xlab='', ylab="Pain Rating", cex.lab=0.7)
lines(2:4, sub2$rating, type="o", pch = 16, cex=0.8, lty = 2, xaxt='n')
axis(1, at=c(2,3,4), labels=c("Subject2; Visit1", "Subject2; Visit2","Subject2; Visit3"), tck = -0.08,cex.axis=0.7)
axis(2, at=c(100,200,300,400,500), labels=c(100,200,300,400,500), tck = -0.08, cex.axis=0.7)
axis(3, at=c(2,3,4), labels=c("Subject1; Visit1", "Subject1; Visit2","Subject1; Visit3"), cex.lab = 0.8, tck = -0.08,cex.axis=0.7)
text(3,550,"(a)", cex=0.8)
mtext("(b)",side=1, cex=0.8, outer=F, line= 1.1, adj=c(0.5))

