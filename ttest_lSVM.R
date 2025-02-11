dataset = "alon"
source("/Functions/functions.R")
library(genefilter)
library(e1071)

data=loadCSV(name=dataset,path="../datasets/")
load("../datasets/seed.RData")

set.seed(saved.seed)
exprs=data$exprs
labels=data$labels

prop=.9                  # ????????? ????????? ? ???????? ?????? ??????
n=200                    # ????? ????????????
ngenes = dim(exprs)[2]                # ????? ?????, ??????????????? ? ???? ?????? alon
size=round(length(labels)*prop)    #????? ????????? ???????? ? ?????? ???????????
print(dim(exprs))	

#?????? ???????? ???????? ?????, ?? ??????? ?? ????? ??????? ??????????????
sigSizes = c(1000, 896, 768, 640, 512, 448, 384, 320, 256, 224, 192, 160, 128, 112, 96, 80, 64, 56, 48, 40, 32, 28, 24, 20, 16, 8, 4, 2)

#?????? ???? list, ?????? ??????? ????? ??????? ????? ??????? ??????? ? ??????????? ?????? ??? 200 ???????????? ??????? ? ????????????? ???????? ????????
sigArray.list <- vector("list", n)

sampling=matrix(0, nrow=n, ncol=size)  #??????? ? ?????????????

for(i in 1:n)
	sampling[i,]=sample(c(1:length(labels)),size,replace=FALSE)

bcr = matrix(0, nrow = n, ncol = length(sigSizes)) #??????? ? BCR

tstart=Sys.time()
for(i in 1:n){
	print(i)
	sigArray_i = matrix(0, nrow = length(sigSizes), ncol = ngenes)

	test=c(1:length(labels))[-sampling[i,]]
	testLabels=labels[test]

	#???????????? ??????
	norm=normTrain(x=exprs[sampling[i,],], robust=TRUE)#try TRUE?
	trainData=norm$values

	m=norm$ms
	sd=norm$sds
	trainLabels=labels[sampling[i,]]

	testData=normTest(x=exprs[test,], ms=m, sds=sd)

	#?????????? ???????? t-?????
	w_index = rankTTest(trainData,trainLabels)[,1]

	#?????????? ?????????????? ?????????? BCR ??? ??????? ?? ???????????? ??? ???????? ?????? ?????
	for (k in 1: length(sigSizes)){
		cur_sigArray = w_index[1:sigSizes[k]]

		#?????????? ?????? ? ??????? SVM ? ???????? ?????
                model=svm(trainData[,cur_sigArray],trainLabels,cost=1,type="C-classification",kernel="linear",scale=FALSE)

		#?????????? BCR		
		if(length(test)==1){res=table(predict=predict(model,testData[,cur_sigArray]),true=testLabels)} else{res=table(predict=predict(model,testData[,cur_sigArray]),true=testLabels)}
		res=formatConfusion(res)

		if(sum(res[,1])==0){bcri=res[2,2]/(res[1,2]+res[2,2])} else if(sum(res[,2])==0){bcri=res[1,1]/(res[1,1]+res[2,1])} else{bcri=((res[1,1]/(res[1,1]+res[2,1]))+(res[2,2]/(res[1,2]+res[2,2])))/2}

		bcr[i,k]=bcri
		sigArray_i[k, cur_sigArray] = 1
		}
	#?????????? ????? ???????? ? ?????? ???????????? ???? ????????
	sigArray.list[[i]] = sigArray_i
}
#?????????? ??????? ???????? ? ???????? BCR
KI = array(0, dim = length(sigSizes))
avg_bcr = array(0, dim = length(sigSizes))
for (i in 1:length(sigSizes)){
	featsLists = matrix(0, nrow=n, ncol=sigSizes[i])
	for (j in 1:n) {temp = which(sigArray.list[[j]][i,]!=0);featsLists[j,] = temp}
	KI[i] = stability(featsLists,ngenes,type = "Kuncheva")	
	avg_bcr[i] = mean(bcr[,i])	
}

#???????? ?????????? ??????????? ?? ????
save(KI, file="KI.RData")
save(avg_bcr, file="avg_bcr.RData")
save(bcr, file = "bcr.RData")
protocol_time = Sys.time() - tstart
print(protocol_time)
print(KI)
print(avg_bcr)

cat(c("Dimensions:", length(labels), length(exprs[1,]),"\n"), file="readme")
cat(c("Dataset", dataset,"\n"), file="readme", append=TRUE)
cat(c("Resamplings", n,"\n"), file="readme", append=TRUE)
cat(c("Protocol time", protocol_time, "\n"), file="readme", append=TRUE)