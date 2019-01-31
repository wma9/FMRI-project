

  for( i in c(0.08)){
    for(sub in c(20,50,200)){
      for(vis in c(20,50)){
        for(k in c(1,0)){
        #path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/power.bonf/rest/20050/bonf_test_power_",k,"_",i,"_seed",j,"_grp200-rep50.RData")
        #path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/prediction/homo/homo_pred_",k,"_",i,"_seed1_grp",sub,"-rep",vis,".RData")
        #path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/prediction/hetero/test/heter_pred_",k,"_",i,"_seed1_grp",sub,"-rep",vis,".RData")
        #path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/prediction/homo/test/homo_pred_",k,"_",i,"_seed1_grp",sub,"-rep",vis,".RData")
        #path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/prediction/addtest/heter/heter_pred_",k,"_",i,"_seed1_grp",sub,"-rep",vis,".RData")
        #path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/prediction/addtest/heter/moreSimu/heter_pred_",k,"_",i,"_seed1_grp",sub,"-rep",vis,".RData")
        #path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/prediction/dummy_addtest/heter/heter_pred_",k,"_",i,"_seed1_grp",sub,"-rep",vis,".RData")
        path <- paste0("/Users/mawanying/Downloads/fwdsubmissionoffmriprojectfolder (1)/Simulation/summer/prediction/dummy_addtest/homo/homo_pred_",k,"_",i,"_seed1_grp",sub,"-rep",vis,".RData")
        load(path)
        cat(paste0(i,"-",sub,"_",vis,"_smooth=",k,"\n"))
        print(est.sim$est.mean)
        cat("--------------------------\n")
      }
    }
  }
  }
  