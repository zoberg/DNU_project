#_______________________________________________________________________________#
dataset="alon"
source("../Functions/functions.R")
source("local_functions.R")
data=loadCSV(name=dataset,path="../../Datasets/")
load("../../Datasets/seed.RData")
set.seed(saved.seed)
exprs=data$exprs
labels=data$labels
#_______________________________________________________________________________#

prop=.9; n=20
bcr=c();
size=round(length(labels)*prop)
print(dim(exprs))

gamma = 1
best_models=NULL
sampling=matrix(0, nrow=n, ncol=size)

for(i in 1:n)
	sampling[i,]=sample(c(1:length(labels)),size,replace=FALSE)

tstart=Sys.time()
for(i in 1:n){
	norm=normTrain(x=exprs[sampling[i,],], robust=FALSE)
	trainData=norm$values
	m=norm$ms
	sd=norm$sds
	trainLabels=labels[sampling[i,]]

	test=c(1:length(labels))[-sampling[i,]]

	testData=normTest(x=exprs[test,], ms=m, sds=sd)
	testLabels=labels[test]
	
	model = SparseLogReg(trainData, trainLabels, gamma, tol = 10^-3)

	best_models=c(best_models, list(model))
	if(length(test)==1){res=table(predict=linear(model,t(as.matrix(testData))),true=testLabels)}
	else{res=table(predict=linear(model,testData),true=testLabels)}
	res=formatConfusion(res)
		
	if(sum(res[,1])==0){bcri=res[2,2]/(res[1,2]+res[2,2])}
	else if(sum(res[,2])==0){bcri=res[1,1]/(res[1,1]+res[2,1])}
	else{bcri=((res[1,1]/(res[1,1]+res[2,1]))+(res[2,2]/(res[1,2]+res[2,2])))/2}
	bcr=c(bcr,bcri)
	print(bcri)
}
protocol_time=Sys.time()-tstart
print(protocol_time)
print("Kuncheva index")
KI_stat=KI(best_models, length(exprs[1,]))
print(KI_stat)
print("Jaccard index")
JIv=JI(best_models, length(exprs[1,]))
print(JIv)
print("average BCR")
print(mean(bcr))

cat(c("Dimensions:", length(labels), length(exprs[1,]),"\n"), file="readme")
cat(c("Dataset", dataset,"\n"), file="readme", append=TRUE)
cat(c("Resamplings", n,"\n"), file="readme", append=TRUE)
cat(c("Protocol time", protocol_time, "\n"), file="readme", append=TRUE)
cat(c("Average BCR", mean(bcr), "\n"), file="readme", append=TRUE)
cat(c("Kuncheva index", KI_stat$KI, "\n"), file="readme", append=TRUE)
cat(c("Jaccard index", JIv, "\n"), file="readme", append=TRUE)
cat(c("Minimal set size", KI_stat$minsetsize, "\n"), file="readme", append=TRUE)
cat(c("Average set size", KI_stat$averagesetsize, "\n"), file="readme", append=TRUE)
unlink("readme")
