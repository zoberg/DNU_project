library(lars)
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
#_______________________________________________________________________________#
cholinsert <- function(R = NULL, xnew,  xold, eps = .Machine$double.eps, Gram = FALSE) 
{
    xtx <- if (Gram) 
        xnew
    else sum(xnew^2)
    norm.xnew <- sqrt(xtx)
    if (is.null(R)) {
        R <- matrix(norm.xnew, 1, 1)
        #attr(R, "rank") <- 1
        return(R)
    }
    Xtx <- if (Gram) 
        xold
    else drop(t(xnew) %*% xold)
    r <- backsolvet(R, Xtx)
    rpp <- norm.xnew^2 - sum(r^2)
    #rank <- attr(R, "rank")
    if (rpp <= eps) 
        rpp <- eps
    else {
        rpp <- sqrt(rpp)
        #rank <- rank + 1
    }
    R <- cbind(rbind(R, 0), c(r, rpp))
    #attr(R, "rank") <- rank
    R
}
#_______________________________________________________________________________#
choldelete <- function (R, k) 
{
    p <- dim(R)[1]
    if (p == 1) 
        return(NULL)
    R <- delcol(R, rep(1, p), k)[[1]][-p, , drop = FALSE]
    #attr(R, "rank") <- p - 1
    R
}
#_______________________________________________________________________________#
lars2 <- function(X, y, stop, maxk)
{
	#library(lars)
	#load("xy.Rdata")
	#stop = par
	#maxk = 400

	Xt = t(X)
	
	n = dim(X)[1]
	p = dim(X)[2]
	nvars = min(n-1,p)
	
	if (missing(maxk)) maxk = 8 * nvars

	beta = array(0, dim=p)
	beta_new = array(0, dim=p)

	mu = y*0
	I = seq(1:p)
	A = c()
	R = c()

	lassocond = 0
	stopcond = 0
	vars = 0
	k = 0

	while((vars < nvars) & (!stopcond) & (k < maxk)){
		k = k + 1
		c = Xt%*%(y - mu)
		C = max(abs(c[I]))
		j = which(abs(c[I])==C)[1]

		if(!lassocond){
			R = cholinsert(R, X[,I[j]], X[,A])
			A = c(A, I[j])
			I = I[-j]
			vars = vars + 1}

		s = sign(c[A])

		GA1 <- backsolve(R, backsolvet(R, s))
		AA = 1/sqrt(sum(GA1%*%s))
		w = AA*GA1
		if (length(w)==1) u = X[, A]*w
		else u = X[, A]%*%w

		if (vars == nvars)
			gamma = C/AA
		else{
			a = Xt%*%u
	        temp = !is.na(c((C - c[I])/(AA - a[I]), (C + c[I])/(AA + a[I])))
			gamma=min(c(temp[temp > 0], C/AA))}

		lassocond = 0
		temp = - beta[A]/w
		gamma_tilde = min(c(temp[temp>0], gamma))
		j = which(temp == gamma_tilde)

		if (gamma_tilde < gamma){
			gamma = gamma_tilde
			lassocond = 1}

		mu = mu + gamma*u
		beta_new[A] = beta[A] + gamma*w

		if (stop > 0){
			t2 = sum(abs(beta_new))
			if (t2 >= stop){
				t1 = sum(abs(beta))
				ss = (stop - t1)/(t2 - t1)
				beta_new = beta + ss*(beta_new - beta)
				stopcond = 1}}

		if (lassocond == 1){
			R = choldelete(R, j)
			I = c(I, A[j])
			A = A[-j]
			vars = vars - 1}

		beta[A] = beta_new[A]
		beta[I] = 0

		if ((stop < 0) & (vars >= -stop))
			stopcond = 1
	}
return(beta)
}
#_______________________________________________________________________________#
irls_lars2 = function(X, y, stop, tol, normalize = FALSE, pobj = 0)
{
	if (missing(stop)) stop=1
	if (missing(tol)) tol=0.0000001
	if (normalize==TRUE){
		x = normTrain(X, robust = FALSE)
		sd = x$sds
		ms = x$ms
		x = x$values}
	else x = X
	m = dim(x)[1]
	n = dim(x)[2]
	w = array(0, dim=n)
	v = log(sum(y>0)/sum(y<0))
	model=list(weights=as.vector(w), intercept = as.numeric(v), gamma = as.numeric(stop), time = Sys.time(), pobj = 0, gap = 0)

	A = diag(y)%*%x

	lc = X%*%w + v				#z2
	z3 = 1 + exp(-lc*y)			#z3

	beta = 0.9
	alpha = 0.3
	prog = 1
	pre_nll = sum(log(z3))/m
	gap = 1

	#max_iter = 10000
	#for(iter in 1:max_iter){
	while(prog > tol){
		#if ((iter %% 3)==0) print(post_nll)
		p = 1/(1 + exp(-lc)) 	#z0
		L = p*(1 - p)			#z1
		logp = (1 - 1/z3)*y 	#z5	plog
		#logp = plog(X,y,model)
		z = lc + (logp/L)		#z4

		A1 = diag(as.vector(sqrt(L)))

		atmp = A1%*%x
		btmp = A1%*%(z - v)

		preobj = sum(log(z3))
		vstep = sum(logp)/sum(L)

		w_ = lars2(atmp, btmp, stop, 400)
		step = w_ - w

		t = 1
		eps = sum((t(logp)%*%x)%*%step)
		for (biter in 1:500){
			w_ = w + t*step
			postobj = sum(log(1+exp(-(x%*%w_ + v)*y)))
			if (postobj < (preobj - alpha*t*eps))
				break
			else
				t = beta*t
		}

		w = w + t*step
		model$weights = as.vector(w)

		t = 1
		eps = sum(logp)*vstep
		for (biter in 1:500){
			v_ = v + t*vstep
			postobj = sum(log(1 + exp(-(x%*%w + v_)*y)))
			if (postobj < (preobj - alpha*t*eps))
				break
			else
				t = beta*t
		}
		v = v + t*vstep
		model$intercept = as.numeric(v)

		lc = x%*%w + v
		z3 = 1 + exp(-lc*y)
		post_nll = sum(log(z3))/m
		if(missing(pobj)) prog = pre_nll - post_nll
		else prog = post_nll - pobj  
		pre_nll = post_nll
	}
	if (pobj!=0) model$gap = pre_nll - pobj
	if (missing(pobj)) model$pobj=pre_nll
	model$time=Sys.time()-model$time
	if (normalize == TRUE){
		temp = solve(diag(as.vector(sd)))
		model$intercept = as.numeric(model$intercept - t(model$weights)%*%temp%*%as.vector(ms))
		model$weights = as.vector(temp%*%model$weights)}
return(model)
}
#_______________________________________________________________________________#
