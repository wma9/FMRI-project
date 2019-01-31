for( i in c(0.0025,0.0049,0.01,0.0144,0.0225,0.04,0.0529,0.0625,0.0729,0.0841,0.0961,0.1225)){
  for(k in c(1,0)){
    cat(paste0("smooth=",k," b.var=",i," 2020 seeds: \n"))
    p=0
    for(j in 1:4){
      cat(paste0("seed",j,"\n"))
      #path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/power.bonf/rest/20050/bonf_test_power_",k,"_",i,"_seed",j,"_grp200-rep50.RData")
      path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/power.bonf.homo/result/2020/test_homo_power_",k,"_",i,"_seed",j,"_grp20-rep20.RData")
      load(path)
      p=p+power2.sim[[1]]$Power
      print(power2.sim[[1]]$Power)
      print(power2.sim[[1]]$realTau)
    }
    cat(paste0("smooth=",k," b.var=",i," 2020 mean is: ", p/4))
    cat("--------------------------\n")
  }
}

index=1
plist=list()
for(k in c(1,0)){
  cat(paste0("smooth=",k))
  pp=c()
  for( i in c(0.0025,0.0049,0.01,0.0144,0.0225,0.04,0.0529,0.0625,0.0729,0.0841)){
    p=0
    for(j in 1:4){
      #path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/power.bonf/rest/20050/bonf_test_power_",k,"_",i,"_seed",j,"_grp200-rep50.RData")
      path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/power.bonf.homo/result/20050/bonf_homo_test_power_",k,"_",i,"_seed",j,"_grp200-rep50.RData")
      load(path)
      p=p+power2.sim[[1]]$Power
    }
    p=p/4
    pp=c(pp,p)
    cat("--------------------------\n")
}
  plist[[index]]=pp
  index=index+1
}


index=1
plist=list()
for(k in c(1,0)){
  cat(paste0("smooth=",k))
  pp=c()
  for( i in c(0.009)){
    p=0
    for(j in 1:4){
      #path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/power.bonf/rest/20050/bonf_test_power_",k,"_",i,"_seed",j,"_grp200-rep50.RData")
      #path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/power.bonf.homo/result/20050/bonf_homo_test_power_",k,"_",i,"_seed",j,"_grp20-rep20.RData")
      #path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/addpoint/bonf_homo_addpoint/bonf_homo_test_power_",k,"_",i,"_seed",j,"_grp200-rep50.RData")
      path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/addpoint/hetero_power_addpoint/heter_test_power_",k,"_",i,"_seed",j,"_grp200-rep50.RData")
      #path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/addpoint/bonf_heter_addpoint/bonf_heter_test_power_",k,"_",i,"_seed",j,"_grp50-rep50.RData")
      #path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/addpoint/homo_power_addpoint/homo_test_power_",k,"_",i,"_seed",j,"_grp50-rep50.RData")
      #path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/addpoint/bonf_homo_addpoint/bonf_homo_test_power_",k,"_",i,"_seed",j,"_grp50-rep50.RData")
      load(path)
      p=p+power2.sim[[1]]$Power
      print(power2.sim[[1]]$Power)
    }
    p=p/4
    pp=c(pp,p)
    cat("--------------------------\n")
  }
  plist[[index]]=pp
  index=index+1
}

