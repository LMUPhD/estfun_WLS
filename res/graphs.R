
library(ggpubr)
library(ggplot2)
library(reshape2)
setwd("C:/Users/classe/Desktop/Diss/Paper3/estfun_WLS/res")

#catL2BB  --> LMuo             
#ordmax   --> WDMo             
#ordL2BB  --> maxLMo               
#supLM    --> maxLM            
#meanL2BB --> CvM             
#maxBB    --> DM



################################################################################
######################## Multivariate Plot #####################################
################################################################################

cond = c("mirt","WLS")
stat = c("maxBB","meanL2BB","supLM","ordL2BB","ordmax","catL2BB")
for(c in cond){
  for(t in stat){
    assign(paste0("tab_",c,"_",t),matrix(NA,nrow=4,ncol=3))
  }
}


anzahl = c("500","1000","2000")
schwellen = c("1","2","4","6")
for(n in 1:3){
  for(s in 1:4){

    load(paste0("multivariate/only betas/240210_",schwellen[s],"_",anzahl[n],"_multi_betas.RData"))
    
    
    for(c in 1:2){
      for(t in 1:6){
        mat = get(paste0("tab_",cond[c],"_",stat[t]))
        mat[s,n] <- resu[t,c]
        if(n==1 & s==1){
          colnames(mat) <- c("n=500","n=1000","n=2000")
          rownames(mat) <- c("k=1","k=2","k=4","k=6")
        }
        assign(paste0("tab_",cond[c],"_",stat[t]),mat)
      }
    }
  }
}




############################### Plot ########################################
#jpeg(filename="multidim_plot.jpg", width=11000, height=10500, res=1200)
jpeg(filename="multidim_plot.jpg", width=5500, height=5250, res=600)

g = c("tab_mirt_maxBB","tab_WLS_maxBB","tab_mirt_meanL2BB","tab_WLS_meanL2BB","tab_mirt_supLM","tab_WLS_supLM","tab_mirt_ordL2BB","tab_WLS_ordL2BB","tab_mirt_ordmax","tab_WLS_ordmax","tab_mirt_catL2BB","tab_WLS_catL2BB")
for(i in 1:length(g)){
  tabb <- 1-get(g[i])
  tab = melt(tabb, id.vars = rownames(tab), variable.name = colnames(tab))
  
  
  ##### correct multivariate
  if(g[i] %in% c("tab_mirt_meanL2BB","tab_WLS_meanL2BB")){
    tab[which(tab$Var1 %in% c("k=6","k=4","k=2")),"value"] = NA
  }
  if(g[i] %in% c("tab_mirt_supLM","tab_WLS_supLM")){
    tab[which(tab$Var1 %in% c("k=6","k=4")),"value"] = NA
  }
  tab = tab[order(tab$Var2),]
  
  pl <- ggplot(tab, aes(x = Var2, y = value, group = Var1)) +  geom_point(aes(color = Var1),size=2) +
    geom_line(aes(color = Var1, linetype = Var1), linewidth = 0.5) + ylim(c(0,1)) +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.title=element_blank()) 
  assign(paste0("plot",i),pl)
}
pl_arr <- ggarrange(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10,plot11,plot12, ncol = 2, nrow = 6,common.legend=T)
annotate_figure(pl_arr, 
                left = text_grob(c("LMuo                    WDMo                maxLMo                 maxLM                    CvM                    DM "), rot = 90, face = "bold"),
                bottom = text_grob(c("      MML                                                                                     WLS"),face="bold"),
                fig.lab = "Multidimensional Model", fig.lab.face = "bold", fig.lab.size = 20)


dev.off()













################################################################################
######################## Univariate Plot #######################################
################################################################################


cond = c("mirt","WLS")
stat = c("maxBB","meanL2BB","supLM","ordL2BB","ordmax","catL2BB")
for(c in cond){
  for(t in stat){
    assign(paste0("tab_",c,"_",t),matrix(NA,nrow=4,ncol=3))
  }
}


anzahl = c("500","1000","2000")
schwellen = c("1","2","4","6")
for(n in 1:3){
  for(s in 1:4){
    
    load(paste0("univariate/only betas/240210_",anzahl[n],"_",schwellen[s],"_univ_betas.RData"))
    
    
    for(c in 1:2){
      for(t in 1:6){
        mat = get(paste0("tab_",cond[c],"_",stat[t]))
        mat[s,n] <- resu[t,c]
        if(n==1 & s==1){
          colnames(mat) <- c("n=500","n=1000","n=2000")
          rownames(mat) <- c("k=1","k=2","k=4","k=6")
        }
        assign(paste0("tab_",cond[c],"_",stat[t]),mat)
      }
    }
  }
}


############################### Plot ########################################

#jpeg(filename="unidim_plot.jpg", width=11000, height=10500, res=1200)
jpeg(filename="unidim_plot.jpg", width=5500, height=5250, res=600)

g = c("tab_mirt_maxBB","tab_WLS_maxBB","tab_mirt_meanL2BB","tab_WLS_meanL2BB","tab_mirt_supLM","tab_WLS_supLM","tab_mirt_ordL2BB","tab_WLS_ordL2BB","tab_mirt_ordmax","tab_WLS_ordmax","tab_mirt_catL2BB","tab_WLS_catL2BB")
for(i in 1:length(g)){
  tabb <- 1-get(g[i])
  tab = melt(tabb, id.vars = rownames(tab), variable.name = colnames(tab))
  
  
  ##### correct univariate
  if(g[i] %in% c("tab_mirt_meanL2BB","tab_WLS_meanL2BB")){
    tab[which(tab$Var1 %in% c("k=6")),"value"] = NA
  }
  tab = tab[order(tab$Var2),]
  
  
  
  pl <- ggplot(tab, aes(x = Var2, y = value, group = Var1)) +  geom_point(aes(color = Var1),size=2) +
    geom_line(aes(color = Var1, linetype = Var1), linewidth = 0.5) + ylim(c(0,1)) +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.title=element_blank()) 
  assign(paste0("plot",i),pl)
}
pl_arr <- ggarrange(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10,plot11,plot12, ncol = 2, nrow = 6,common.legend=T)
annotate_figure(pl_arr, 
                left = text_grob(c("LMuo                    WDMo                maxLMo                 maxLM                    CvM                    DM "), rot = 90, face = "bold"),
                bottom = text_grob(c("      MML                                                                                      WLS"),face="bold"),
                fig.lab = "Unidimensional Model", fig.lab.face = "bold", fig.lab.size = 20)


dev.off()



################################################################################
######################## Computation time #####################################
################################################################################

cond = c("mirt","WLS")
for(c in cond){
  assign(paste0("tab_",c),matrix(NA,nrow=4,ncol=3))
}


anzahl = c("500","1000","2000")
schwellen = c("1","2","4","6")
for(n in 1:3){
  for(s in 1:4){
    
    load(paste0("univariate/all parameters/231124_",anzahl[n],"_",schwellen[s],"_result_univ_simpleWLS.RData"))
    
    
    for(c in 1:2){
      
      mat = get(paste0("tab_",cond[c]))
      mat[s,n] <- resu_t[4,c]
      if(n==1 & s==1){
        colnames(mat) <- c("n=500","n=1000","n=2000")
        rownames(mat) <- c("k=1","k=2","k=4","k=6")
      }
      assign(paste0("tab_",cond[c]),mat)
      
    }
  }
}


wb=loadWorkbook("computation_time.xlsx")
sh=addWorksheet(wb,"computation_time univ")
writeData(wb, sh, tab_mirt, startCol = 1, startRow = 1, rowNames = T)
writeData(wb, sh, tab_WLS, startCol = 1, startRow = 10, rowNames = T)


saveWorkbook(wb, "computation_time.xlsx", overwrite = TRUE)
















