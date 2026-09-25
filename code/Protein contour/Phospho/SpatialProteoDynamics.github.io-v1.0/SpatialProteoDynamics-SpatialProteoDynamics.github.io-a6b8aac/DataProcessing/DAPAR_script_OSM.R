library(DAPAR)
library(Prostar)
library(dplyr)
library(MSnbase)
library(imp4p)
library(limma)
setwd("D:/00-PROJECTS/01-AMV/02-Subcellular fractionation/04b-OSMOTIC STRESS v14/01-PROTEOME directDIA")
exprsFile<-c("Report_Proteins-directDIA_v14_SUBCELL_OSMSTRSS_PROTEOME.txt")
exprsData<-read.table(exprsFile, head=T, sep="\t", quote="")
exprsData<-exprsData[which(exprsData$PG.Genes!=""),]
exprsData$PG.Genes<-gsub("^;", "", exprsData$PG.Genes)
exprsData$PG.Genes<-gsub(";.*", "", exprsData$PG.Genes)
exprsData<-exprsData[!duplicated(exprsData$PG.Genes),]
colnames_exprsData<-colnames(exprsData)
colnames_exprsData<-gsub(".raw.PG.Quantity", "", colnames_exprsData)
colnames_exprsData<-gsub(".*OSMSTR_", "OSMSTR", colnames_exprsData)
colnames(exprsData)<-colnames_exprsData
metadata<-read.table("Groups.txt", head=T, sep="\t")
indExpData<-c(1:20)
indiceID<-21
fractions<-c("FR1", "FR2", "FR3", "FR4", "FR5", "FR6")
prot_tables<-list()
group_tables<-list()
for(i in 1:6){prot_tables[[i]]<-dplyr::select(exprsData,contains(fractions[i]));
			prot_tables[[i]]<-cbind(prot_tables[[i]],exprsData$PG.Genes)}
for(i in 1:6){group_tables[[i]]<-dplyr::filter(metadata,grepl(fractions[i], Sample.name));
				group_tables[[i]][,2]<-as.character(group_tables[[i]][,2]);
				group_tables[[i]][,2]<-as.factor(group_tables[[i]][,2])}


id_table<-as.character(exprsData$PG.Genes)
merge_table<-as.data.frame(id_table)
# merge_table_average<-as.data.frame(id_table)
# merge_table_sd<-as.data.frame(id_table)
row.names(merge_table)<-merge_table[,1]
# row.names(merge_table_average)<-merge_table[,1]
# row.names(merge_table_sd)<-merge_table[,1]

#FR1
for(j in 1:6){
prot<-createMSnset(prot_tables[[j]], group_tables[[j]],indExpData, indiceID=21, indexForOriginOfValue = NULL, pep_prot_data = "protein", logData=T)
prot_mv<-mvFilter(prot, type="atLeastOneCond", th=3)
m<-group_tables[[j]][2]								
mv<-t(apply(exprs(prot_mv), 1, function(x) tapply(x,m,function(x) {length(x[!is.na(x)])})))
mv<-(mv[,c(paste("1hSORB3hREL_FR",j,sep=""), paste("1hSORB24hREL_FR",j,sep=""), paste("1hSORB30minREL_FR",j,sep=""),  paste("1hSORBITOL_FR",j,sep=""), paste("CTRL_FR",j,sep=""))])
mv_keep<-data.frame(matrix(vector(),nrow(exprs(prot_mv)),4, dimnames=list(c(), c( "1hSORB3hvsC","24hvsC", "1hSORB30minvsC", "1hSORBITOLvsC"))))
for(i in 1:nrow(exprs(prot_mv))){
	if(mv[i,5]>=3 || mv[i,1]>=3){mv_keep[i,1]<-"KEEP"}
}
for(i in 1:nrow(exprs(prot_mv))){
	if (mv[i,5]>=3 || mv[i,2]>=3){mv_keep[i,2]<-"KEEP"}
}
for(i in 1:nrow(exprs(prot_mv))){
	if (mv[i,5]>=3 || mv[i,3]>=3){mv_keep[i,3]<-"KEEP"}
}
for(i in 1:nrow(exprs(prot_mv))){
	if (mv[i,5]>=3 || mv[i,4]>=3){mv_keep[i,4]<-"KEEP"}
}

# prot_norm<-wrapper.normalizeD(prot_mv, method="LOESS", type="overall")
# prot_imp_1<-impute.slsa(exprs(prot_norm), group_tables[[j]]$Condition)
# qData <- (prot_imp_1)
# values <- getQuantile4Imp(qData)$shiftedImpVal
# prot_imp_2<-impute.detQuant(qData, values)


prot_norm<-wrapper.normalizeD(prot_mv, method="LOESS", type="overall")
prot_imp_KNN<-wrapper.impute.KNN(prot_norm, 15)
MEC<-findMECBlock(prot_norm)
prot_MEC<-reIntroduceMEC(prot_imp_KNN, MEC)
qData <- exprs(prot_MEC)
values <- getQuantile4Imp(qData)$shiftedImpVal
prot_imp_2<-impute.detQuant(qData, values)


for(k in 1:4){
	input_limma<-cbind(prot_imp_2[,(4*k-3):(k*4)],prot_imp_2[,17:20])
	input_limma<-input_limma[which(mv_keep[,k]=="KEEP"),]
	conditions<-rbind(pData(prot)[(4*k-3):(k*4),],pData(prot)[17:20,])
	prot_limma<-limmaCompleteTest(input_limma, conditions, comp.type="OnevsOne")
	BH_pvalue<-p.adjust(as.numeric(unlist(prot_limma$P_Value)), method = "BH", n = length(as.numeric(unlist(prot_limma$P_Value))))
	write.table(cbind(prot_limma$logFC, prot_limma$P_Value, BH_pvalue), paste("OSM_STRSS_limma_COND",colnames(mv_keep)[k],"FR",j,".csv", sep=""),col.names=NA, sep=",") 
	rnk<-prot_limma$logFC
	names(rnk)<-NULL
	write.table(rnk, paste("OSM_STRSS_RNK_COND",colnames(mv_keep)[k],"FR",j,".rnk", sep=""),col.names=F, sep="\t", quote=FALSE) 

} 

ave_prot_imp<-apply((prot_imp_2), 1, function(x) tapply(x,m, mean))
raw_prot_imp<-2^(t(ave_prot_imp))
prot_imp_sd<-apply(2^(prot_imp_2), 1, function(x) tapply(x,m, sd))

merge_table<-merge(merge_table, prot_imp_2, by=0, all=T)
row.names(merge_table)<-merge_table[,1]
merge_table<-dplyr::select(merge_table,contains("OSMSTR"))

write.table((raw_prot_imp), paste("OSM_STRSS_imputed_average_raw_FR",j,".csv", sep=""), col.names=NA, sep=",")
write.table(t(ave_prot_imp), paste("OSM_STRSS_imputed_average_log2_FR",j,".csv", sep=""), col.names=NA, sep=",")
write.table(t(prot_imp_sd), paste("OSM_STRSS_imputed_sd_FR",j,".csv", sep=""), col.names=NA, sep=",")

prot_results<-cbind(mv_keep, prot_imp_2)
write.table(prot_results, paste("OSM_STRSS_imputed_FR",j,".csv", sep=""), col.names=NA, sep=",")
}

merge_table_na<-merge_table[rowSums(is.na(merge_table))<120,]
#merge_table_norm<-normalizeMedianValues(merge_table_na)
write.table(merge_table_na, "merged_table_OSMSTRSS_prot.txt", sep="\t")
