
library(ggpubr)
library(ggplot2)
library(reshape2)
setwd("C:/Users/classe/Desktop/Diss/Paper3/estfun_WLS/res")

################################################################################
######################## Multivariate Plot #####################################
################################################################################

cond = c("mirt","WLS")
stat = c("maxBB","ordmax","ordL2BB","catL2BB")
for(c in cond){
  for(t in stat){
    assign(paste0("tab_",c,"_",t),matrix(NA,nrow=4,ncol=3))
  }
}


anzahl = c("500","1000","2000")
schwellen = c("1","2","4","6")
for(n in 1:3){
  for(s in 1:4){

    load(paste0("simpleWLS/multivariate/231125_",schwellen[s],"_",anzahl[n],"_result_multi_simpleWLS.RData"))
    
    
    for(c in 1:2){
      for(t in 1:4){
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


################################################################################
############################### Sandbox ########################################
################################################################################

g = c("tab_mirt_maxBB","tab_WLS_maxBB","tab_mirt_ordmax","tab_WLS_ordmax","tab_mirt_ordL2BB","tab_WLS_ordL2BB","tab_mirt_catL2BB","tab_WLS_catL2BB")
for(i in 1:length(g)){
  tab = melt(get(g[i]), id.vars = rownames(tab), variable.name = colnames(tab))
  
  pl <- ggplot(tab, aes(x = Var2, y = value, group = Var1)) +  geom_point(aes(color = Var1),size=2) +
    geom_line(aes(color = Var1, linetype = Var1), size = 0.5) + ylim(c(0,0.1)) +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.title=element_blank()) 
    #labs(linetype="number of\nthresholds",color="number of\nthresholds") #x="sample size",y="alpha error", 
  assign(paste0("plot",i),pl)
}
pl_arr <- ggarrange(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8, ncol = 2, nrow = 4,common.legend=T)
annotate_figure(pl_arr, 
                left = text_grob(c("catL2BB                     ordL2BB                       ordmax                        maxBB "), rot = 90, face = "bold"),
                bottom = text_grob(c("mirt                                                                   WLS"),face="bold"))

