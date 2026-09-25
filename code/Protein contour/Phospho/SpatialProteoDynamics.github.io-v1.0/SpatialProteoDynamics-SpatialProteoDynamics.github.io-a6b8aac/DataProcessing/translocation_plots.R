
library(metap)
library(dplyr)
library(ggplot2)
library(ggrepel)
setwd("C:/Users/svg993/Desktop/testR/translocation")
full_table<-read.table("merged_table_HeLa_EGF_prot.txt", head=T, sep="\t")
pval_FR1<-read.table("EGF_prot_limma_CONDX2minEGFvsCFR1.csv", sep=",", head=T)
pval_FR2<-read.table("EGF_prot_limma_CONDX2minEGFvsCFR2.csv", sep=",", head=T)
pval_FR3<-read.table("EGF_prot_limma_CONDX2minEGFvsCFR3.csv", sep=",", head=T)
pval_FR4<-read.table("EGF_prot_limma_CONDX2minEGFvsCFR4.csv", sep=",", head=T)
pval_FR5<-read.table("EGF_prot_limma_CONDX2minEGFvsCFR5.csv", sep=",", head=T)
pval_FR6<-read.table("EGF_prot_limma_CONDX2minEGFvsCFR6.csv", sep=",", head=T)

groups<-read.table("Groups.txt", head=T, sep="\t")
prot_2min<-dplyr::select(full_table,contains("2min"))
prot_ctrl<-dplyr::select(full_table,contains("CTRL"))

prot_2min<-2^prot_2min
prot_ctrl<-2^prot_ctrl


FR<-c(rep("FR1",4), rep("FR2", 4),  rep("FR3", 4),  rep("FR4", 4),  rep("FR5", 4),  rep("FR6", 4))
prot_2min_mean<-t(apply(prot_2min, 1, function(x) tapply(x,FR,function(x) {mean(x)})))
prot_CTRL_mean<-t(apply(prot_ctrl, 1, function(x) tapply(x,FR,function(x) {mean(x)})))

prot_2min_mean_raw<-prot_2min_mean
prot_CTRL_mean_raw<-prot_CTRL_mean

prot_2min_scaled<-t(apply(prot_2min_mean_raw, 1, function(x) {x/sum(x)}))
prot_Ctrl_scaled<-t(apply(prot_CTRL_mean_raw, 1, function(x) {x/sum(x)}))

MS_2min<-abs(prot_2min_scaled-prot_Ctrl_scaled)
MS_2min_mean<-apply(MS_2min, 1, function(x) mean(x))

MS_max1<-(apply(MS_2min, 1, function(x) {match(max(x),x)}))

maxN <- function(x, N=2){
  len <- length(x)
  if(N>len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x,partial=len-N+1)[len-N+1]
}

MS_max2<-(apply(MS_2min, 1, function(x) {match(maxN(x),x)}))


isEmpty <- function(x) {
  return(identical(x, numeric(0)))
}

temp_table<-cbind(prot_2min_scaled, prot_Ctrl_scaled, MS_2min_mean, MS_max1, MS_max2)
temp_table<-as.data.frame(temp_table)
row.names(temp_table)<-full_table$Gene.names

pval<-c()

for(i in 1:nrow(temp_table)){
  
  if(temp_table$MS_max1[i]=="1"){
    x=pval_FR1[pval_FR1$X==row.names(temp_table)[i],3]
    if(isEmpty(x)){
      pval[i]=1}
    else{pval[i]<-pval_FR1[pval_FR1$X==row.names(temp_table)[i],3]}
  }
  
  else if(temp_table$MS_max1[i]=="2"){
    x=pval_FR2[pval_FR2$X==row.names(temp_table)[i],3]
    if(isEmpty(x)){
      pval[i]=1}
    else{pval[i]<-pval_FR2[pval_FR2$X==row.names(temp_table)[i],3]}
  } 
  else if(temp_table$MS_max1[i]=="3"){
    x=pval_FR3[pval_FR3$X==row.names(temp_table)[i],3]
    if(isEmpty(x)){
      pval[i]=1}    
    else{pval[i]<-pval_FR3[pval_FR3$X==row.names(temp_table)[i],3]}
  }
    
  else if(temp_table$MS_max1[i]=="4"){
    x=pval_FR4[pval_FR4$X==row.names(temp_table)[i],3]
    if(isEmpty(x)){
      pval[i]=1}
    else{pval[i]<-pval_FR4[pval_FR4$X==row.names(temp_table)[i],3]}
  }
    
  else if(temp_table$MS_max1[i]=="5"){
    x=pval_FR5[pval_FR5$X==row.names(temp_table)[i],3]
    if(isEmpty(x)){
      pval[i]=1}
    else{pval[i]<-pval_FR5[pval_FR5$X==row.names(temp_table)[i],3]}
  }
    
  else if(temp_table$MS_max1[i]=="6"){
    x=pval_FR6[pval_FR6$X==row.names(temp_table)[i],3]
    if(isEmpty(x)){
      pval[i]=1}
    else{pval[i]<-pval_FR6[pval_FR6$X==row.names(temp_table)[i],3]}
  }
}

pval_2<-c()

for(i in 1:nrow(temp_table)){
  
  if(temp_table$MS_max2[i]=="1"){
    x=pval_FR1[pval_FR1$X==row.names(temp_table)[i],3]
    if(isEmpty(x)){
      pval_2[i]=1}
    else{pval_2[i]<-pval_FR1[pval_FR1$X==row.names(temp_table)[i],3]}
  }
  
  else if(temp_table$MS_max2[i]=="2"){
    x=pval_FR2[pval_FR2$X==row.names(temp_table)[i],3]
    if(isEmpty(x)){
      pval_2[i]=1}
    else{pval_2[i]<-pval_FR2[pval_FR2$X==row.names(temp_table)[i],3]}
  } 
  else if(temp_table$MS_max2[i]=="3"){
    x=pval_FR3[pval_FR3$X==row.names(temp_table)[i],3]
    if(isEmpty(x)){
      pval_2[i]=1}    
    else{pval_2[i]<-pval_FR3[pval_FR3$X==row.names(temp_table)[i],3]}
  }
  
  else if(temp_table$MS_max2[i]=="4"){
    x=pval_FR4[pval_FR4$X==row.names(temp_table)[i],3]
    if(isEmpty(x)){
      pval_2[i]=1}
    else{pval_2[i]<-pval_FR4[pval_FR4$X==row.names(temp_table)[i],3]}
  }
  
  else if(temp_table$MS_max2[i]=="5"){
    x=pval_FR5[pval_FR5$X==row.names(temp_table)[i],3]
    if(isEmpty(x)){
      pval_2[i]=1}
    else{pval_2[i]<-pval_FR5[pval_FR5$X==row.names(temp_table)[i],3]}
  }
  
  else if(temp_table$MS_max2[i]=="6"){
    x=pval_FR6[pval_FR6$X==row.names(temp_table)[i],3]
    if(isEmpty(x)){
      pval_2[i]=1}
    else{pval_2[i]<-pval_FR6[pval_FR6$X==row.names(temp_table)[i],3]}
  }
}

pval_combi<-c()
for(i in 1:length(pval)){
	pval_combi<-c(pval_combi, sumlog(c(pval[i],pval_2[i]))$p)
}

color<-c()

for(i in 1:nrow(temp_table)){
  if((temp_table$MS_max1[i]==1 || temp_table$MS_max1[i]==2) && (temp_table$MS_max2[i]==3 || temp_table$MS_max2[i]==4)){
    color[i]="limegreen"
  }else if((temp_table$MS_max1[i]==3 || temp_table$MS_max1[i]==4) && (temp_table$MS_max2[i]==1 || temp_table$MS_max2[i]==2)){
    color[i]="limegreen"
  }else if((temp_table$MS_max1[i]==1 || temp_table$MS_max1[i]==2) && (temp_table$MS_max2[i]==5 || temp_table$MS_max2[i]==6)){
    color[i]="dodgerblue"
  }else if((temp_table$MS_max1[i]==5 || temp_table$MS_max1[i]==6) && (temp_table$MS_max2[i]==1 || temp_table$MS_max2[i]==2)){
    color[i]="dodgerblue"
  }else if((temp_table$MS_max1[i]==3 || temp_table$MS_max1[i]==4) && (temp_table$MS_max2[i]==5 || temp_table$MS_max2[i]==6)){
    color[i]="coral"
  }else if((temp_table$MS_max1[i]==5 || temp_table$MS_max1[i]==6) && (temp_table$MS_max2[i]==3 || temp_table$MS_max2[i]==4)){
    color[i]="coral"
  }else{color[i]="gray"}
}
temp_table<-cbind(temp_table, color, pval_combi)
pval_combi_FDR<-p.adjust(temp_table$pval_combi, method="BH")
temp_table<-cbind(temp_table, pval_combi_FDR)
temp_table<-temp_table[,13:18]
temp_table<-cbind(temp_table, row.names(temp_table))
colnames(temp_table)[6]<-"pval_combi_FDR"
colnames(temp_table)[7]<-"Gene"
plot_2min<-ggplot(temp_table, aes(x=MS_2min_mean, y=-log10(pval_combi_FDR)))+
  geom_point(color=color)+
  geom_text_repel(aes(label=Gene), data=temp_table[(temp_table$pval_combi_FDR<0.05) & (temp_table$MS_2min_mean>0.1),])+
  scale_x_continuous(limits=c(0,0.3))+
  scale_y_continuous(limits=c(0,6))+
  theme_bw()

