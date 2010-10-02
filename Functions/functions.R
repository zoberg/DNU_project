# name:	String. Default is 'alon'. Must be one of 'alon', 'golub', 'immuno', 'prostate', 'prostateWelsh', 'prostateChandran', 'prostateSingh' or 'lymphomaNoNA'.
# path: String. Path to the folder containing the files alon.csv, golub.csv, immuno.csv, prostate.csv, prostateWelsh.csv, prostateChandran.csv, prostateSingh.csv and lymphomaNoNA.csv.
#
# return:	A list containing to entries: 
#			- $labels contains a column matrix nx1 with the labels of the n samples.
#			- $exprs contains a nxd matrix with the expression data of n samples and d probesets. 
#
loadCSV<-function(name="alon",path=PATH$DATADIR){
	if(name!="alon" & name!="golub" & name!="immuno" & name!="prostate" & name!="prostateWelsh" & name!="prostateChandran" & name!="prostateSingh" & name!="lymphomaNoNA" & name!="DLBCL"){
		cat("You have to specify one of the following for 'name': 'alon', 'golub', 'immuno', 'prostate', 'DLBCL', 'prostateWelsh', 'prostateChandran', 'prostateSingh' or 'lymphomaNoNA'.\n")
		return(NULL)
	}
	f=paste(path,name,".csv",sep="")
	d=read.table(file=f,sep=",",quote="",dec=".",header=TRUE,row.names=1)
	data=list(labels=as.matrix(d[,1]),exprs=as.matrix(d[,-1]))
	return(data)
}
#-------------------------------------------------------------------------------------
# Args:
#	x: A nxm matrix
#	robust: Boolean.
#
# Return:	A list with three elements:
#				$ms: A 1xm matrix. If robust is TRUE, $ms elements are the median of each column of 'x'. If robust is FALSE, the mean is used instead of the median.
#				$sds: A 1xm matrix. If robust is TRUE, $sds elements are the IQRs of each column of 'x'. If robust is FALSE, the s.d. is used instead of the IQRs.
#				$values: A nxm matrix. $values(i,j)=(x(i,j)-$ms(j))/$sds(j)
normTrain<-function(x,robust=TRUE){
	m=matrix(nr=1,nc=dim(x)[2])
	sd=matrix(nr=1,nc=dim(x)[2])
	for(i in 1:dim(x)[2]){
		if(robust){
			mediani=median(x[,i])
			# iqr = inter-quantile range = quantile(3/4)-quantile(1/4)
			# divide by 1.35 (normality assumption)
			iqri=IQR(x[,i])/1.35
			x[,i]=(x[,i]-mediani)/iqri
			m[i]=mediani
			sd[i]=iqri
		}
		else{
			meani=mean(x[,i])
			sdi=sd(x[,i])
			if(sdi!=0){x[,i]=(x[,i]-meani)/sdi}
			m[i]=meani
			sd[i]=sdi
		}
	}
	params=list(ms=m,sds=sd,values=x)
	return(params)
}
#-------------------------------------------------------------------------------------
# Args:
#	x: A nxm matrix
#	ms: A m-length vector or a 1xm matrix.
#	sds: A m-length vector or a 1xm matrix.
#
# Return: A nxm matrix with element(i,j)=(x(i,j)-ms(j))/sds(j)
normTest<-function(x,ms,sds){
	if(!is.null(dim(x))){
		for(i in 1:dim(x)[2]){
			x[,i]=(x[,i]-ms[i])/sds[i]
		}
	}
	else{
		for(i in 1:length(x)){
			x[i]=(x[i]-ms[i])/sds[i]
		}
	}
	return(x)
}
#-------------------------------------------------------------------------------------
# Args: 
#	res: the result of table() function. Must be a 1x2, 2x2 or 2x1 matrix with dimnames in the set ("-1","1").
#
# Return: a 2x2 matrix with (("-1","1"),("-1","1")) dimnames and the corresponding values from res at the right position in the matrix. This function standardize 2x2 confusion matrices. 
formatConfusion<-function(res){
	#print(res)
	if(dim(res)[1]==2 & dim(res)[2]==2){
		res=res[c("-1","1"),c("-1","1")]
	}
	else if(dim(res)[2]==2){
		res=res[c(rownames(res),rownames(res)),c("-1","1")]
		tmp=matrix(nc=2,nr=2,data=0)
		colnames(tmp)=c("-1","1")
		rownames(tmp)=c("-1","1")
		if(rownames(res)[1]=="-1"){
			tmp[1,1]=res[1,1]
			tmp[1,2]=res[1,2]
			tmp[2,1]=0
			tmp[2,2]=0
		}
		else{
			tmp[1,1]=0
			tmp[1,2]=0
			tmp[2,1]=res[1,1]
			tmp[2,2]=res[1,2]
		}
		res=tmp
	}
	else if(dim(res)[1]==2){
		res=res[c("-1","1"),c(colnames(res),colnames(res))]
		tmp=matrix(nc=2,nr=2,data=0)
		colnames(tmp)=c("-1","1")
		rownames(tmp)=c("-1","1")
		if(colnames(res)[1]=="-1"){
			tmp[1,1]=res[1,1]
			tmp[2,1]=res[2,1]
			tmp[1,2]=0
			tmp[2,2]=0
		}
		else{
			tmp[1,1]=0
			tmp[2,1]=0
			tmp[1,2]=res[1,1]
			tmp[2,2]=res[2,1]
		}
		res=tmp
	}
	else{
		tmp=matrix(nc=2,nr=2,data=0)
		if(colnames(res)=="-1" && rownames(res)=="-1"){
			tmp[1,1]=res
			tmp[2,1]=0
			tmp[1,2]=0
			tmp[2,2]=0
			res=tmp
		}
		else if(colnames(res)=="1" && rownames(res)=="1"){
			tmp[1,1]=0
			tmp[2,1]=0
			tmp[1,2]=0
			tmp[2,2]=res
			res=tmp
		}
		else if(colnames(res)=="-1" && rownames(res)=="1"){
			tmp[1,1]=0
			tmp[2,1]=res
			tmp[1,2]=0
			tmp[2,2]=0
			res=tmp
		}
		else{
			tmp[1,1]=0
			tmp[2,1]=0
			tmp[1,2]=res
			tmp[2,2]=0
			res=tmp
		}
	}
	dimnames(res)=list(predict=c("-1","1"),true=c("-1","1"))
	return(res)	
}
#-------------------------------------------------------------------------------------
#_______________________________________________________________________________#
JI <- function(models_list, input_dim)
{
	partitions = length(models_list)
	best_sets=matrix(0,nrow=partitions,ncol=input_dim)
	for (i in 1:partitions)
		best_sets[i,] = unlist(models_list[i])[1:input_dim]
	sum = 0
	for (i in 1:(partitions-1)){
		for(j in (i+1):partitions){
			tempi = which(best_sets[i,]!=0)
			tempj = which(best_sets[j,]!=0)
			intersect = length(tempi[match(tempj, tempi, nomatch = 0)])
			sum = sum + intersect / (length(tempi) + length(tempj) -intersect)
		}
	}
	return(2*sum/(partitions*(partitions - 1)))
}
#_______________________________________________________________________________#
KI <- function(models_list, input_dim)
{
	partitions = length(models_list)
	min_set=input_dim
	best_sets=matrix(0,nrow=partitions,ncol=input_dim)
	for (i in 1:partitions){
		best_sets[i,] = unlist(models_list[i])[1:input_dim]
		if (length(which(best_sets[i,]!=0))<min_set) min_set=length(which(best_sets[i,]!=0))}

	sum=0
	tempi=array(0, dim=min_set)
	tempj=array(0, dim=min_set)

	for (i in 1:(partitions-1)){
		for (k in 1:min_set){
		    temp = which(abs(best_sets[i,])==sort(abs(best_sets[i,]),decreasing=TRUE)[k])
		    if (length(temp)>1){ 
				count = 1
				for (m in 1:length(temp))
					if (tempi[k-m]==temp[1]) count=count+1
				tempi[k]=temp[count]
				}
			    else tempi[k] = temp}
	
		for (j in i:partitions){
			for (k in 1:min_set){
			    temp = which(abs(best_sets[j,])==sort(abs(best_sets[j,]),decreasing=TRUE)[k])	    
			    if (length(temp)>1){ 
				count = 1
				for (m in 1:length(temp))
					if (tempj[k-m]==temp[1]) count=count+1
				tempj[k]=temp[count]
				}
			    else tempj[k] = temp}
			sum = sum + (length(tempi[match(tempi,tempj, nomatch=0)])-min_set*min_set/input_dim)/(min_set-min_set*min_set/input_dim)
		}
	}
	return(list(KI=sum*2/(partitions*(partitions-1)), minsetsize=min_set, averagesetsize=(length(which(best_sets!=0))/partitions)))
}
