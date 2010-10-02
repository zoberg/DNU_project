dyn.load("logreg.so")
#_______________________________________________________________________________#
standardization <- function(X)
{
	m = dim(X)[1]
	n = dim(X)[2]
	mu = array(0, dim = n)
	sigma = mu
	Xstd = matrix(0, ncol = n, nrow = m)
	for(i in 1:n){
		mu[i] = mean(X[,i])
		sigma[i] = sd(X[,i])
		Xstd[, i] = (X[, i] - mu[i])/sigma[i]
	}
	#X = diag(y)%*%(X - diag(array(1, dim = m))%*%t(mu))%*%ginv(diag(sigma))
out = list(values = Xstd, mu = mu, sigma = sigma)
return(out)
}
#_______________________________________________________________________________#
linear <- function(model, data)
{
	alpha = as.vector(model$weights)
	intercept = as.numeric(model$intercept)
	out=array(0, dim=length(data[,1]))
	out=sign(intercept + data%*%alpha)
	return(out)
}
#_______________________________________________________________________________#
linear2 <- function(model, data)
{
	alpha = as.vector(model$weights)
	intercept = as.numeric(model$intercept)
	out=array(0, dim=length(data[,1]))
	out=sign(data%*%alpha)
	return(out)
}
#_______________________________________________________________________________#
SparseLogReg <- function(Xb, y, gamma, tol, normalize = FALSE)
{
	if (missing(tol)) tol=0.000001
	if (missing(gamma)) gamma=1
	if (normalize==TRUE){
		Xs = normTrain(Xb, robust = FALSE)
		X = Xs$values
		sigma = Xs$sds
		mu = Xs$ms}
	else{X = Xb}
	n = length(X[1,])
	m = length(y)
	#if (gamma>=maxregval(X,y)*m) stop("The regularization parameter is too big, no solution. Maximal value of the parameter is: ", maxregval(X,y)*m, call. = FALSE)
	time=Sys.time()
	out<-.C("my_SLR", 
		as.numeric(tol), 
		as.integer(n),
		as.numeric(X), 
		as.integer(c(-13, y)), 
		as.integer(c(-13, seq(1:m))), 
		as.integer(m), 
		as.numeric(gamma), 
		as.numeric(0), 
		as.numeric(array(0, dim=n+1)), 
		as.numeric(0))
	time=Sys.time()-time
	model = list(weights = as.vector(unlist(out[9])[2:(n+1)]), intercept = - as.numeric(out[10]), gamma = gamma/m, gap = 0, l1n = 0, time = time)

	model$l1n = sum(abs(model$weights))*model$gamma
	if (normalize==TRUE){
	temp = solve(diag(as.vector(sigma)))
	model$intercept = as.numeric(- model$intercept - t(model$weights)%*%temp%*%as.vector(mu))
	model$weights = temp%*%model$weights
	model$gap = l_avg(Xb, y, model) + model$l1n - G(gap(Xb,y,model)$theta)}
return(model)
}
#_______________________________________________________________________________#
