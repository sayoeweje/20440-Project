## run from within the monomer folder
## there can only be one *_scoresEnergetics and one Centroid/*_scoresCentroid in this folder
## args[1] is the folder for the multimer
## NOTE: results may differ slightly from publication because of missing data in the mutagenesis experiments used in teh paper

args <- commandArgs(trailingOnly = TRUE)

#Load monomer
data1 = read.table(dir("./","_scoresEnergetics$"),header=T,row.names=1)
print(data1[,1:2])
data1Z = as.data.frame(apply(data1,2,function(x){return(scale(as.numeric(x)))}))
rownames(data1Z) = rownames(data1)
data1Z[is.na(data1Z)] = 0
data2 = read.table(dir("./Centroid","_scoresCentroid$", full.names=TRUE),header=T,row.names=1)
data2Z = as.data.frame(apply(data2,2,function(x){return(scale(as.numeric(x)))}))
rownames(data2Z) = rownames(data2)
data2Z[is.na(data2Z)] = 0
colnames(data2Z) = paste(colnames(data2Z),"CENTROIDSC",sep="")
data2Z = data2Z[rownames(data1Z),]
data2 = data2[rownames(data1),]

#Add RSA
data1 = cbind(data2[,"RSA"],data1)
colnames(data1)[1]="RSA"
dataZ = cbind(data1Z,data2Z)

#Identify max NodeEdgeBetweennessSTRIDE_sidechain
acidsMonomer = c()
for (i in 1:length(rownames(dataZ))) {
    node = rownames(dataZ)[i]
    node = unlist(strsplit(node,split=""))
    node = paste(node[1:(length(node)-1)],collapse="")
    acidsMonomer = c(acidsMonomer,node)
}
allAcids = unique(c(acidsMonomer))

#NodeEdgeBetweenness is the max of the monomer/multimer
NodeEdgeBetweennessSTRIDE_AVERAGE_MONOMER = (dataZ[,"NodeEdgeBetweennessSTRIDE"]+dataZ[,"NodeEdgeBetweennessSTRIDE_sidechain"])/2
NodeEdgeBetweennessSTRIDE_sidechain_ALL = c(NodeEdgeBetweennessSTRIDE_AVERAGE_MONOMER)
NodeEdgeBetweennessSTRIDE_sidechain_MAX = tapply(NodeEdgeBetweennessSTRIDE_sidechain_ALL,as.factor(c(acidsMonomer)),max,na.rm=TRUE)
NodeEdgeBetweennessSTRIDE_sidechain_MAX = NodeEdgeBetweennessSTRIDE_sidechain_MAX[allAcids]

#RSA
RSA_ALL = c(data1[,"RSA"])
RSA_MIN = tapply(RSA_ALL,as.factor(c(acidsMonomer)),min,na.rm=TRUE)
RSA_MIN = RSA_MIN[allAcids]

y = NodeEdgeBetweennessSTRIDE_sidechain_MAX

out=cbind(names(y),y)
write.table(out,file="FinalSum",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
