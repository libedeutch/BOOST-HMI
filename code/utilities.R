perf_metric <- function(gamma, pval, cutoff=0.05,mistbin=FALSE, bf_cutoff){
  # pval for BOOST-HM needs only positive bayes factor 
  # need table(re,gamma) true is columns 
  if(mistbin){
    re <- ifelse(pval>=bf_cutoff, TRUE,FALSE)
  }
  else{
    re <- ifelse(pval<=cutoff, TRUE,FALSE)
  }
  confuse_matrix = table(re,gamma)
  if(nrow(confuse_matrix)==1 & rownames(confuse_matrix)[1]=='TRUE'){
    confuse_matrix = rbind( c(0,0),confuse_matrix )
    row.names(confuse_matrix) = c("FALSE","TRUE")
  }
  if(nrow(confuse_matrix)==1& rownames(confuse_matrix)[1]=='FALSE'){
    confuse_matrix = rbind( confuse_matrix,c(0,0) )
    row.names(confuse_matrix) = c("FALSE","TRUE")
  }
  TN = confuse_matrix[1]
  FP = confuse_matrix[2]
  FN = confuse_matrix[3]
  TP = confuse_matrix[4]
  # colnames(result_sub) = c("Sensitivity","Specificity" ,"F1 score" ,"FDR", "AUC" ,"MCC")
  Sensitivity = ifelse(TP==0,0,TP/(TP+FN))
  Specificity = ifelse(TN==0,0,TN/(TN+FP))
  F1_score = ifelse(TP==0,0,2*TP/(2*TP+FP+FN))
  FDR = ifelse(FP==0,0,FP/(FP+TP))
  if(mistbin) {AUC = auc(gamma,pval,direction = "<")} #,direction = "<"
  else{AUC = auc(gamma,pval,direction = ">");  }
  MCC = ifelse((TP+FP)==0|(TP+FN)==0|(TN+FP)==0|(TN+FN)==0,0,
               (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  return(list(Sensitivity =Sensitivity,Specificity = Specificity,F1_score = F1_score,FDR = FDR, AUC = AUC,MCC = MCC ))
}

morani <- function(loc, x){
  weight<- as.matrix(dist(cbind(loc[,1], loc[,2])))
  weight <- 1/weight
  diag(weight) <- 0
  xbar = mean(x)
  n = nrow(loc)
  #weight = matrix(1, nrow = n, ncol = n)
  s0 = sum(weight)
  denom = sum((x - xbar)^2)
  num = 0
  for(i in 1:n){
    for(j in 1:n){
      if(j !=i){
        num = num+weight[i,j]*(x[i] -xbar)*(x[j]-xbar)
      }
    }
  }
  return (n*num/s0/denom)
} # it works the same as MoranI in ape package 

gearyc <- function(loc, x){ # convert gearyc to (-1,1)
  weight <- as.matrix(dist(cbind(loc[,1], loc[,2])))
  weight <- 1/weight
  diag(weight) <- 0
  xbar = mean(x)
  n = nrow(loc)
  s0 = sum(weight)
  denom = sum((x - xbar)^2)
  num = 0
  for(i in 1:n){
    for(j in 1:n){
      if(j !=i){
        num = num+weight[i,j]*(x[i] -x[j])^2
      }
    }
  }
  tmp = (n-1)*num/2/s0/denom
  return (1-tmp)
}


bf<- function(rre){
  iter = nrow(rre[[1]]$mu)
  burn = iter/2; 
  bf_neg <- rep(NA, length(rre))
  bf_pos <- rep(NA, length(rre))
  for (i in 1:length(rre)){
    theta = rre[[i]]$theta[burn:iter,2]
    #bf.pos <- (sum(theta>=0/length(theta))/((sum(theta<0))/length(theta)+0.001)
    bf_neg[i] <- (sum(theta<=0)/length(theta))/((sum(theta>0))/length(theta)+0.0001)
    bf_pos[i] = (sum(theta>=0)/length(theta))/((sum(theta<0))/length(theta)+0.001)
  }
  return(list(bf_neg = bf_neg, bf_pos = bf_pos))
}
bf2<- function(rre){
  iter = nrow(theta)
  burn = iter/2; 
  bf_neg <- rep(NA, length(rre))
  bf_pos <- rep(NA, length(rre))
  for (i in 1:length(rre)){
    theta = theta[burn:iter,2]
    #bf.pos <- (sum(theta>=0/length(theta))/((sum(theta<0))/length(theta)+0.001)
    bf_neg[i] <- (sum(theta<=0)/length(theta))/((sum(theta>0))/length(theta)+0.0001)
    bf_pos[i] = (sum(theta>=0)/length(theta))/((sum(theta<0))/length(theta)+0.001)
  }
  return(list(bf_neg = bf_neg, bf_pos = bf_pos))
}

bf2<- function(rre,rep = 3){
  iter = nrow(rre[[1]]$theta)
  burn = iter/2; 
  bf_neg <- rep(NA, length(rre)/rep)
  bf_pos <- rep(NA, length(rre)/rep)
  for (i in 1:(length(rre)/3)){
    theta = c(rre[[i]]$theta[burn:iter,2],rre[[i+249]]$theta[burn:iter,2],rre[[i+249*2]]$theta[burn:iter,2])
    #bf.pos <- (sum(theta>=0/length(theta))/((sum(theta<0))/length(theta)+0.001)
    bf_neg[i] <- (sum(theta<=0)/length(theta))/((sum(theta>0))/length(theta)+0.0001)
    bf_pos[i] = (sum(theta>=0)/length(theta))/((sum(theta<0))/length(theta)+0.001)
  }
  return(list(bf_neg = bf_neg, bf_pos = bf_pos))
}

bf3<- function(thetaT){
  iter = nrow(thetaT)
  burn = iter/2; 
  bf_neg <- rep(NA, ncol(thetaT))
  bf_pos <- rep(NA, ncol(thetaT))
  for (i in 1:(ncol(thetaT))){
    theta = c(thetaT[burn:iter,i])
    #bf.pos <- (sum(theta>=0/length(theta))/((sum(theta<0))/length(theta)+0.001)
    bf_neg[i] <- (sum(theta<=0)/length(theta))/((sum(theta>0))/length(theta)+0.0001)
    bf_pos[i] = (sum(theta>=0)/length(theta))/((sum(theta<0))/length(theta)+0.001)
  }
  return(list(bf_neg = bf_neg, bf_pos = bf_pos))
}
thetaF = cbind(theta2,theta3,theta4,theta5)
bf4<- function(theta1,theta2,theta3,theta4,theta5){ # four chains thetaF <- cbind(thta2,3,4,5)
  iter = nrow(theta2)
  burn = iter/2; 
  bf_neg <- rep(NA, ncol(theta2))
  bf_pos <- rep(NA, ncol(theta2))
  for (i in 1:(ncol(theta2))){
    theta = c(theta1[burn:iter,i],theta2[burn:iter,i],theta3[burn:iter,i],theta4[burn:iter,i],theta5[burn:iter,i])
    #bf.pos <- (sum(theta>=0/length(theta))/((sum(theta<0))/length(theta)+0.001)
    bf_neg[i] <- (sum(theta<=0)/length(theta))/((sum(theta>0))/length(theta)+0.0001)
    bf_pos[i] = (sum(theta>=0)/length(theta))/((sum(theta<0))/length(theta)+0.001)
  }
  return(list(bf_neg = bf_neg, bf_pos = bf_pos))
}


bf5<- function(theta1,theta2,theta3,theta4){ # four chains thetaF <- cbind(thta2,3,4,5)
  iter = nrow(theta2)
  burn = iter/2; 
  bf_neg <- rep(NA, ncol(theta2))
  bf_pos <- rep(NA, ncol(theta2))
  for (i in 1:(ncol(theta2))){
    theta = c(theta1[burn:iter,i],theta2[burn:iter,i],theta3[burn:iter,i],theta4[burn:iter,i])
    #bf.pos <- (sum(theta>=0/length(theta))/((sum(theta<0))/length(theta)+0.001)
    bf_neg[i] <- (sum(theta<=0)/length(theta))/((sum(theta>0))/length(theta)+0.0001)
    bf_pos[i] = (sum(theta>=0)/length(theta))/((sum(theta<0))/length(theta)+0.001)
  }
  return(list(bf_neg = bf_neg, bf_pos = bf_pos))
}

###################### BinSpect plots for seqfish ##################

idx_bskmeans = idx_spk[!idx_spk%in%idx_hm]# spk only 
n = length(idx_bskmean)
x = rep(loc[,1], n )
y = rep(loc[,2], n)
col = c()
for(i in idx_bskmean){
  #col = c(col,colMeans(rbind(rre2[[i]]$z_colmeans,rre2[[i+249]]$z_colmeans,rre2[[i+249*2]]$z_colmeans )))  
  col = c(col,colMeans(rbind(rre1[[i]]$z_colmeans,rre2[[i]]$z_colmeans,rre3[[i]]$z_colmeans,rre4[[i]]$z_colmeans, rre5[[i]]$z_colmeans)))  
}
index = rep(idx_bskmean, each = nrow(loc))
label_bskmean = rep(label[idx_bskmean],each = nrow(loc))
# ggdata = tibble(x=x,y=y,col = col,index = index,label_spko = label_spko )%>% 
#   mutate (index = factor(index,levels = idx_spko,
#                          labels = paste0(gene_symbol[idx_spko]," BF ", round(gg$bf_neg[idx_spko],3))))

ggdata = tibble(x=x,y=y,col = col,index = index,label_bskmean = label_bskmean )%>% 
  mutate (index = factor(index,levels = idx_bskmean,
                         labels = paste0(gene_symbol[idx_bskmean]," (", 
                                         #sprintf("%.3f", round(gg$bf_neg[idx_spko],3)),")")))
                                         formatC(round(gg$bf_neg[idx_bskmean],3), format = "e",digits = 3),")")))
ggplot(ggdata, aes(x, y, col = col,group = index)) +
  geom_point(aes(col = col),size = 2.5) +
  #facet_wrap(~label_spko)+
  facet_wrap(~index)+
  scale_color_gradient2(low = "black",high = "red",midpoint=0.5, mid = "gray",name = "Relative expression",limits = c(0,1))+
  theme(
    panel.border = element_rect(colour = "gray", fill=NA, linewidth = 1),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(face = "bold",size = 10),
    legend.position = c(0.55,0.08),
    legend.direction = "horizontal",
    text = element_text(size=12),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )+ xlab("") +ylab("")+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))

ggsave("/Users/gilly/Library/CloudStorage/OneDrive-TheUniversityofTexasatDallas/fish/seqfish/paper_figure/bskmeans.png",
       dpi = 300, width = 20,height = 18)
ggsave("/Users/gilly/Library/CloudStorage/OneDrive-TheUniversityofTexasatDallas/fish/seqfish/paper_figure/bskmeans.pdf",
       dpi = 300, width = 20,height = 18)








n = length(idx_bsrank)
x = rep(loc[,1], n )
y = rep(loc[,2], n)
col = c()
for(i in idx_bsrank){
  #col = c(col,colMeans(rbind(rre2[[i]]$z_colmeans,rre2[[i+249]]$z_colmeans,rre2[[i+249*2]]$z_colmeans )))  
  col = c(col,colMeans(rbind(rre1[[i]]$z_colmeans,rre2[[i]]$z_colmeans,rre3[[i]]$z_colmeans,rre4[[i]]$z_colmeans, rre5[[i]]$z_colmeans)))  
}
index = rep(idx_bsrank, each = nrow(loc))
label_bsrank = rep(label[idx_bsrank],each = nrow(loc))
# ggdata = tibble(x=x,y=y,col = col,index = index,label_spko = label_spko )%>% 
#   mutate (index = factor(index,levels = idx_spko,
#                          labels = paste0(gene_symbol[idx_spko]," BF ", round(gg$bf_neg[idx_spko],3))))

ggdata = tibble(x=x,y=y,col = col,index = index,label_bsrank = label_bsrank )%>% 
  mutate (index = factor(index,levels = idx_bsrank,
                         labels = paste0(gene_symbol[idx_bsrank]," (", 
                                         #sprintf("%.3f", round(gg$bf_neg[idx_spko],3)),")")))
                                         formatC(round(gg$bf_neg[idx_bsrank],3), format = "e",digits = 3),")")))
ggplot(ggdata, aes(x, y, col = col,group = index)) +
  geom_point(aes(col = col),size = 2.5) +
  #facet_wrap(~label_spko)+
  facet_wrap(~index)+
  scale_color_gradient2(low = "black",high = "red",midpoint=0.5, mid = "gray",name = "Relative expression",limits = c(0,1))+
  theme(
    panel.border = element_rect(colour = "gray", fill=NA, linewidth = 1),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(face = "bold",size = 10),
    legend.position = c(0.4,0.08),
    legend.direction = "horizontal",
    text = element_text(size=12),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )+ xlab("") +ylab("")+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))

ggsave("/Users/gilly/Library/CloudStorage/OneDrive-TheUniversityofTexasatDallas/fish/seqfish/paper_figure/bsrank.png",
       dpi = 300, width = 20,height = 20)
ggsave("/Users/gilly/Library/CloudStorage/OneDrive-TheUniversityofTexasatDallas/fish/seqfish/paper_figure/bsrank.pdf",
       dpi = 300, width = 20,height = 20)






##################### Venn plot ####################

library(ggplot2)
library(gridExtra)
library(VennDiagram)
x = list(BOOST_HMI = gene_symbol[idx_hm], 
         BinSpect_kmeans = gene_symbol[idx_bskmean], BinSpect_rank = gene_symbol[idx_bsrank],SPARK = gene_symbol[idx_spk])
names(x) <- c("BOOST-HMI", "BinSpect-kmeans", "BinSpect-rank","SPARK")

ggVennDiagram(x,edge_size = 0,set_size = 4,cat.just=list(c(0.6,1) , c(0,0) , c(0,0) , c(0.2,0.5)))+ # , "#D3E2B7","#F7C97E"
  scale_fill_gradient(low =  "gray87", high = "#ECA8A9")+
  theme(legend.position = "None")

ggsave("/Users/gilly/Library/CloudStorage/OneDrive-TheUniversityofTexasatDallas/fish/seqfish/paper_figure/venn4.png",
       dpi = 300, width = 12,height = 6)
ggsave("/Users/gilly/Library/CloudStorage/OneDrive-TheUniversityofTexasatDallas/fish/seqfish/paper_figure/venn4.pdf",
       dpi = 300, width = 12,height = 6)