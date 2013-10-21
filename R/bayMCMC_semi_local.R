bayMCMC_semi_local <-
function(data_x, data_y, data_xnew, Xvar, Xvarpred, warm=1000, M=1000, mutprob=0.44, errorprob=0.44, epsilonprob = 0.44, mutsizp=1.0, errorsizp=1.0, epsilonsizp = 1.0, 
	prior_alpha=1.0, prior_beta=0.05, err_int = c(-5,5), err_ngrid=10001, num_batch=20, step=10, alpha=0.95, ...)
{	
	data_y <- as.vector(data_y)	
	 if (is.vector(data_xnew)) 
        		data_xnew <- as.matrix(t(data_xnew))
	    testfordim <- sum(dim(data_x) == dim(data_xnew)) == 2
	    twodatasets <- TRUE
	    if (testfordim) 
	        twodatasets <- sum(data_x == data_xnew) != prod(dim(data_x))
	
	SPECURVES1 = data_x
	Specresp1 = data_y

	# prior density IG(1.0,0.05)

	logpriorh2 = function(h2)
	{
		logp = prior_alpha * log(prior_beta) - lgamma(prior_alpha)
		logp = logp - 1.0*(prior_alpha + 1.0)*log(h2) - prior_beta/h2	
		return(logp)
	}
	
	########################
	# negative log posterior
	########################

	cost = function(xp)
	{
		n = length(Specresp1)		
		Wmat = funopare.kernel(Specresp1, SPECURVES1, SPECURVES1, bandwidth = exp(xp[1]), ...)$NWweit
		if(is.null(Wmat))
		{	
			result = -100000
		}
		else
		{
			tildacurves = (diag(n) - Wmat) %*% Xvar
			tilday = (diag(n) - Wmat) %*% as.matrix(Specresp1)
			betahat = solve(t(tildacurves)%*%tildacurves)%*%t(tildacurves)%*%tilday
		
			resi = Specresp1 - Xvar%*%betahat
			NW = funopare.kernel(resi, SPECURVES1, SPECURVES1, bandwidth = exp(xp[1]), ...)$Estimated.values
			yhat = NW + Xvar%*%betahat	
			resid = as.vector(Specresp1 - yhat)
			epsilon = scale(resid)
			std = sd(resid)
			cont = (2.0*pi)^(-0.5)
		
			b = exp(xp[2])
			badj = exp(xp[3])/(1.0+exp(xp[3]))
			logf = vector(,length(resid))
			for(i in 1:length(resid))
			{
				temp = epsilon[i] - epsilon[-i]
				res = sum(cont*exp(-0.5*((temp/(b*(1+badj*abs(epsilon[-i]))))^2))/(b*(1+badj*abs(epsilon[-i]))))
				logf[i] = log(res/length(temp)/std)
			}
			sumlogf = sum(logf)
			#log Jacobi and log prior
			priorJacobi = vector(,2)
			for(i in 1:2)
			{
				priorJacobi[i] = xp[i] + logpriorh2((exp(xp[i]))^2)
			}
			result = sumlogf + sum(priorJacobi) + log(badj) + log(1.0 - badj)
		}
		return(-result)
	}
	
# sampling bandwidth h used in the functional NW estimator
	gibbs_mean = function(xp, k, mutsizp)
	{
		fx = xp[4]
		dv = rnorm(1)*mutsizp
		
		xp[1] = xp[1] + dv
		fy = cost(xp)
		rvalue = fx - fy
		
		if(is.nan(rvalue)) accept = 0
		else
		{
			if(fx > fy) accept = 1
			else
			{
				un = runif(1)
				if(un < exp(rvalue)) accept = 1
				else accept = 0
			}
		}
		accept_mean=0
		mutsizc = mutsizp/(mutprob * (1.0 - mutprob))
		if(accept == 1)
		{
			accept_mean = accept_mean + 1
			xp[4] = fy
			mutsizp = mutsizp + mutsizc*(1.0 -mutprob)/k
		}	
		else
		{
			xp[1] = xp[1] - dv
			mutsizp = mutsizp - mutsizc*mutprob/k
		}		
		return(list(xpnew = xp, mutsizpnew = mutsizp, mutsizcnew = mutsizc, acceptnw = accept_mean)) 
	}	
	# sampling bandwidth b used in the kernel-form error density
	gibbs_erro = function(xp, k, errorsizp)
	{
		fx = xp[4]
		dv = rnorm(1)*errorsizp
		
		xp[2] = xp[2] + dv
		fy = cost(xp)
		rvalue = fx - fy
		if(is.nan(rvalue)) accept = 0
		else
		{
			if(fx > fy) accept = 1
			else
			{
				un = runif(1)
				if(un < exp(rvalue)) accept = 1
				else accept = 0
			}
		}
		accept_erro=0
		errorsizc = errorsizp/(errorprob * (1.0 - errorprob))
		if(accept == 1)
		{
			accept_erro = accept_erro + 1
			xp[4] = fy
			errorsizp = errorsizp + errorsizc*(1-errorprob)/k
		}
		else	
		{
			xp[2] = xp[2] - dv
			errorsizp = errorsizp - errorsizc*errorprob/k
		}	
		return(list(xperronew = xp, errorsizpnew = errorsizp, errorsizcnew = errorsizc, accepterro = accept_erro))
	}
	# sampling localized bandwidth adjustment tau_{\varepsilon}
	gibbs_epsilon = function(xp, k, epsilonsizp)
	{
		fx = xp[4]
		dv = rnorm(1)*epsilonsizp
		
		xp[3] = xp[3] + dv
		fy = cost(xp)
		rvalue = fx - fy
		if(is.nan(rvalue)) accept = 0
		else
		{
			if(fx > fy) accept = 1
			else
			{
				un = runif(1)
				if(un < exp(rvalue)) accept = 1
				else accept = 0
			}
		}
		accept_epsilon = 0
		epsilonsizc = epsilonsizp/(epsilonprob * (1.0 - epsilonprob))
		if(accept == 1)
		{
			accept_epsilon = accept_epsilon + 1
			xp[4] = fy
			epsilonsizp = epsilonsizp + epsilonsizc*(1-epsilonprob)/k
		}
		else
		{
			xp[3] = xp[3] - dv
			epsilonsizp = epsilonsizp - epsilonsizc*epsilonprob/k
		}
		return(list(xpepsilonnew = xp, epsilonsizpnew = epsilonsizp, epsilonsizcnew = epsilonsizc, acceptepsilon = accept_epsilon))
	}
	## warm-up stage		
	# initial values		
	ini_val = runif(3,min=1,max=3)
	xp = c(ini_val, cost(xp = ini_val))
	acceptnw = accepterro = acceptepsilon = vector(,warm)
	xpwarm = matrix(,warm,4)
	# burn-in period
	for(k in 1:warm)
	{
			dum = gibbs_mean(xp, k, mutsizp)
			xp = dum$xpnew
			acceptnw[k] = dum$acceptnw
			mutsizp = dum$mutsizpnew
			
			dum2 = gibbs_erro(xp, k, errorsizp)
			xp = dum2$xperronew
			accepterro[k] = dum2$accepterro
			errorsizp = dum2$errorsizpnew
			
			dum3 = gibbs_epsilon(xp, k, epsilonsizp)
			xp = dum3$xpepsilonnew
			acceptepsilon[k] = dum3$acceptepsilon
			epsilonsizp = dum3$epsilonsizpnew
			xpwarm[k,] = xp		
	}
	# MCMC recording
	acceptnwMCMC = accepterroMCMC = acceptepsilonMCMC = vector(,M)
	xpM = xpMsquare = matrix(,M,4)
	cpost = matrix(, M/step, 4)
	for(k in 1:M)
	{
			dumMCMC = gibbs_mean(xp, k+warm, mutsizp)
			xp = dumMCMC$xpnew
			acceptnwMCMC[k] = dumMCMC$acceptnw
			mutsizp = dumMCMC$mutsizpnew
			
			dum2MCMC = gibbs_erro(xp, k+warm, errorsizp)
			xp = dum2MCMC$xperronew
			accepterroMCMC[k] = dum2MCMC$accepterro
			errorsizp = dum2MCMC$errorsizpnew
			
			dum3MCMC = gibbs_epsilon(xp, k+warm, epsilonsizp)
			xp = dum3MCMC$xpepsilonnew
			acceptepsilonMCMC[k] = dum3MCMC$acceptepsilon
			epsilonsizp = dum3MCMC$epsilonsizpnew
			xpM[k,] = xp	
			xpMsquare[k,] = exp(xp)^2
			index = ceiling(k/step)
			cpost[index,] = exp(xp)^2	
	}

	# ergodic average
	xpfinalres = colMeans(xpM)
	# obtaining the bandwidth of regression and residuals,
	kernelestfinal = funopare.kernel(Specresp1, SPECURVES1, SPECURVES1, bandwidth = exp(xpfinalres[1]), ...)
	Wmatkernelestfinal = kernelestfinal$NWweit							
	n = dim(Wmatkernelestfinal)[1]
		
	tildacurves = (diag(n) - Wmatkernelestfinal) %*% Xvar
	tilday = (diag(n) - Wmatkernelestfinal) %*% as.matrix(Specresp1)
	betahat = solve(t(tildacurves)%*%tildacurves)%*%t(tildacurves)%*%tilday
		
	resi = Specresp1 - Xvar%*%betahat
	NW = funopare.kernel(Response = resi, CURVES = SPECURVES1, PRED = data_xnew, bandwidth = exp(xpfinalres[1]), ...)
	regressionfunctionestimate = (Xvar%*%betahat + NW$Estimated.values)	
	residfinal = Specresp1 - regressionfunctionestimate

	SIF <- function(BAND_MATRIX, NUM_ITERATIONS, NUM_BATCH) 
	{
	        	size_batch = NUM_ITERATIONS/NUM_BATCH
	        h_mean <- colMeans(BAND_MATRIX)
	        h_mean <-  matrix(rep(h_mean,times= dim(BAND_MATRIX)[1]), nrow = dim(BAND_MATRIX)[1], 
					ncol = dim(BAND_MATRIX)[2], byrow = TRUE)
	        sigma_square_tilde <- (1/(NUM_ITERATIONS - 1)) * colSums((BAND_MATRIX - h_mean)^2) 
	        sum_par_mean = rep(0, times = dim(BAND_MATRIX)[2])        
        		for (i in 1: NUM_BATCH)
	        {                
        	        	sum_par_mean = sum_par_mean + (colMeans(BAND_MATRIX[((i-1) * size_batch + 1):(size_batch * i),]) - h_mean[1,] )^2            
	        }
        	sigma_square_hat = size_batch/(NUM_BATCH - 1) * sum_par_mean
	        var_h = sqrt(sum_par_mean/(NUM_BATCH^2 - NUM_BATCH))
	        return(list(batch_se = round(sqrt(sigma_square_hat),4), 
			total_se = round(sqrt(sigma_square_tilde),4), 
			SIF = round(sigma_square_hat/sigma_square_tilde,4), VAR_H = round(var_h,4)))
	}
	sif_value = SIF(exp(xpM[,1:(ncol(xpM)-1)]), M, num_batch)

	logpriors_admkr = function(h2)
	{
		logf = 0
		dm = length(h2)
		for(i in 1:(dm-1))
		{
			logf = logf + logpriorh2(h2[i])
		}
		logf = logf + logpriorh2(h2[dm])
		return(logf)
	}
	
	loglikelihood_admkr = function(h, resid)
	{
	    dm = length(h)
		b = h[2]
		badj = h[3]
		    
		epsilon = scale(resid)
		std = sd(resid)
		cont = (2.0*pi)^(-0.5)
		
		logf = vector(,length(resid))
		for(i in 1:length(resid))
		{
			temp = epsilon[i] - epsilon[-i]
			res = sum(cont*exp(-0.5*((temp/(b*(1+badj*abs(epsilon[-i]))))^2))/(b*(1+badj*abs(epsilon[-i]))))
			logf[i] = log(res/length(temp)/std)
		}
		sumlogf = sum(logf)			
  		return(sumlogf)
	}


	logdensity_admkr = function(tau2, cpost)
	{
	    dm = ncol(cpost)
	    len = nrow(cpost)
	    band = vector(, dm)
	    for(j in 1:dm)
	    {
        	temp = tem2 = 0
	        for(i in 1:len)
	        {
        	    temp = temp + cpost[i,j]
	            tem2 = tem2 + cpost[i,j]^2
	        }
	        sigma = sqrt(tem2/len - (temp/len)^2)
	        temp = exp(1.0/(dm + 4) * log(4/(dm + 2)))
	        band[j] = temp * sigma * exp(-1/(dm + 4) * log(len))
	    }
	    hprod = prod(band)
	    cont = exp(-0.5 * (dm) * log(2.0 * pi))
	    xsum = 0
	    for(i in 1:len)
	    {
        	temp = 0
	        for(j in 1:dm)
	        {
	            tem2 = (tau2[j] - cpost[i,j])/band[j]
         	   temp = temp + tem2^2
	        }
	        xsum = xsum + cont * exp(-0.5 * temp)/hprod
	    }
	    hatf = log(xsum/len)
	    return(hatf)
	}

        mlikeres = loglikelihood_admkr(exp(xpfinalres[1:3]), residfinal) + logpriors_admkr(exp(xpfinalres[1:3])^2) - 
			logdensity_admkr(colMeans(xpMsquare[,1:3]), cpost[,1:3])

	# kernel error density estimation
	admkr.denadj = function(ban, badj, eps, res.data)
	{
		res = as.numeric(res.data)
		data.num = length(res)
		std.res = sd(res)
		epsilon = (res-mean(res))/std.res
		eps.std = (eps-mean(res))/std.res
		tem = tem2 = vector(,length(epsilon))
		for(i in 1:(length(epsilon)))
		{
			tem[i]=(eps.std-epsilon[i])/(ban*(1+badj*abs(res[i])))
			tem2[i]=dnorm(tem[i])/(ban*(1+badj*abs(res[i])))
		}		
		hatf=(sum(tem2)/data.num)/std.res
		return(hatf)		
	}
	admkr.cdfadj = function(ban, badj, eps, res.data)
	{
		res = as.numeric(res.data)
		data.num = length(res)
		std.res = sd(res)
		epsilon = (res-mean(res))/std.res
		eps.std = (eps-mean(res))/std.res
		tem = tem2 = vector(,length(epsilon))
		for(i in 1:(length(epsilon)))
		{
			tem[i]=(eps.std-epsilon[i])/(ban*(1+badj*abs(res[i])))
			tem2[i]=pnorm(tem[i])
		}
		hatf = sum(tem2)/data.num
		return(hatf)
	}			
	# approximate ISE	
	y = seq(err_int[1], err_int[2], by = diff(err_int)/(err_ngrid-1))
	fore.den.mkr = fore.cdf.mkr = vector(,length(y))
	for(i in 1:(length(y)))
	{
		eps = y[i]
		fore.den.mkr[i] = admkr.denadj(exp(xpfinalres[2]), exp(xpfinalres[3]), eps, residfinal)
		fore.cdf.mkr[i] = admkr.cdfadj(exp(xpfinalres[2]), exp(xpfinalres[3]), eps, residfinal)
	}
	if (twodatasets) 
	{
		NWpred = NW$Predicted.values + Xvarpred%*%betahat					
		PI = NWpred + c(y[which.min(abs(fore.cdf.mkr - (1-alpha)/2))], y[which.min(abs(fore.cdf.mkr - (1+alpha)/2))])
		return(list(xpfinalres = c(exp(xpfinalres[1:2]), exp(xpfinalres[3])/(1.0+exp(xpfinalres[3]))), mhat = regressionfunctionestimate, betahat = betahat,
			sif_value = sif_value, mlikeres = mlikeres,
			acceptnwMCMC = mean(acceptnwMCMC), accepterroMCMC = mean(accepterroMCMC), 
			acceptepsilonMCMC = mean(acceptepsilonMCMC),
			fore.den.mkr = fore.den.mkr, fore.cdf.mkr = fore.cdf.mkr, pointforecast=NWpred, PI = PI))
	}
	else
	{
		return(list(xpfinalres = c(exp(xpfinalres[1:2]), exp(xpfinalres[3])/(1.0+exp(xpfinalres[3]))), mhat = regressionfunctionestimate, betahat = betahat,
			sif_value = sif_value, mlikeres = mlikeres,
			acceptnwMCMC = mean(acceptnwMCMC), accepterroMCMC = mean(accepterroMCMC), 
			acceptepsilonMCMC = mean(acceptepsilonMCMC),
			fore.den.mkr = fore.den.mkr, fore.cdf.mkr = fore.cdf.mkr))
	}
}
