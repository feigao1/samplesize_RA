#### Sample Size Calculation based on Z ####
samplesize_CF <- function(para, R1, R0 = 1, alpha=0.05, beta=0.8){
	gamma0 = rec_gamma(para)
	gamma1 = trial_gamma(para,para$lambda0*R1)
	gamma00 = gamma0$gamma00; gamma01 = gamma0$gamma01; gamma1 = gamma1$gamma1
	vZ = varZ(para, R1, R0)
	Zsum = qnorm(1-alpha/2) + qnorm(beta,0,sqrt(vZ))
	return(ceiling((gamma00+gamma1)/(((log(R1)-log(R0))/Zsum)^2-gamma01)))
}

power_CF <- function(para, R1, R0 = 1, alpha=0.05, N){
	gamma0 = rec_gamma(para,N)
	gamma1 = trial_gamma(para,para$lambda0*R1,N)
	gamma00 = gamma0$gamma00; gamma01 = gamma0$gamma01; gamma1 = gamma1$gamma1
	qbeta = abs((log(R1)-log(R0))/(sqrt((gamma00+gamma1)/N+gamma01)))-qnorm(1-alpha/2)
	vZ = varZ(para, R1, R0)
	beta = pnorm(qbeta,0,sqrt(vZ))
	return(beta)
}

#### Data Generation ####
# N: total screening sample size
# para: lambda0, p, q, r, OmegaT, sigmaOmegaT, betaT, sigmabetaT, Time, tau
# R: risk ratio, default=1
generatedata <- function(N, para, R=1){
	p = para$p; q = para$q; r = para$r; lambda0 = para$lambda0; lambda1 = lambda0*R; tau = para$tau
	OmegaT = para$OmegaT; sigmaOmegaT = para$sigmaOmegaT; betaT = para$betaT; sigmabetaT = para$sigmabetaT; Time = para$Time
	Npos = rbinom(1,N,p); Nneg = N - Npos
	Npos_test = rbinom(1,Npos,q); Nneg_trial = rbinom(1,Nneg,r)
	PR = betaT + lambda0*(1-p)/p*(OmegaT-betaT*Time)
	NR = rbinom(1,Npos_test,PR)
	wbetaT = rnorm(1,betaT,sigmabetaT); wOmegaT = rnorm(1,OmegaT,sigmaOmegaT)
	Nevent = rpois(1,lambda1*tau*Nneg_trial)
	return(list(Npos = Npos, Nneg = Nneg, Npos_test = Npos_test, Nneg_trial = Nneg_trial, NR = NR, Nevent = Nevent, PR = PR, wbetaT = wbetaT, wOmegaT = wOmegaT))
}

### lambda0 estimation ###
rec_est <- function(N,Npos,Npos_test,NR,para){
	pest = Npos/N; qest = Npos_test/Npos; Nneg = N - Npos
	betaT = para$betaT; OmegaT = para$OmegaT; Time = para$Time
	lambdaest = (NR/qest - betaT*Npos) / (Nneg*(OmegaT-betaT*Time))
	para_est = list(lambda0=lambdaest, p=pest, q=qest, OmegaT=OmegaT, sigmaOmegaT = para$sigmaOmegaT, betaT=betaT, sigmabetaT = para$sigmabetaT, Time=Time)
	varlog = rec_gamma(para_est,N)$varlog
	Est=lambdaest; SE = sqrt(varlog)*lambdaest;
	CI = c(lambdaest-qnorm(0.975)*SE, lambdaest+qnorm(0.975)*SE)
	CI_log = c(lambdaest * exp(-qnorm(0.975)*sqrt(varlog)), lambdaest * exp(qnorm(0.975)*sqrt(varlog)))
	return(list(Est=Est,SE = SE,CI = CI, CI_log = CI_log))
}

### lambda1 estimation ###
trial_est <- function(N,Nneg,Nneg_trial,Nevent,tau){
	pest = 1-Nneg/N; rest = Nneg_trial/Nneg
	lambdaest = Nevent /Nneg_trial/tau
	para_est = list(p=pest,r=rest,tau=tau)
	varlog = trial_gamma(para_est,lambdaest, N)$varlog
	Est=lambdaest; SE = sqrt(varlog)*lambdaest;
	CI = c(lambdaest - qnorm(0.975)*SE, lambdaest + qnorm(0.975)*SE)
	CI_log = c(lambdaest * exp(-qnorm(0.975)*sqrt(varlog)), lambdaest * exp(qnorm(0.975)*sqrt(varlog)))
	return(list(Est=Est,SE = SE,CI = CI, CI_log = CI_log))
}

#### Rejection Indicator based on Z ####
Rejection <- function(lambda0_est, lambda1_est, R0 = 1, Zcut = qnorm(0.975)){
	if (lambda0_est$Est<=0) return(F)
	if (lambda1_est$Est==0) {
		Zstat = (log(lambda0_est$Est))/(lambda0_est$SE/lambda0_est$Est)
		return(abs(Zstat)>Zcut)
	}
	Zstat = efficacy_test(lambda0_est,lambda1_est, R0)
	return(abs(Zstat)>Zcut)
}

### Z calculation ###
efficacy_test<-function(lambda0_est,lambda1_est, R0 = 1){
	lambda0 = lambda0_est$Est; lambda0_SE = lambda0_est$SE; V0 = (lambda0_SE/lambda0)^2
	lambda1 = lambda1_est$Est; lambda1_SE = lambda1_est$SE; V1 = (lambda1_SE/lambda1)^2
	wR = lambda1/lambda0; V = V0+V1
	Z = (log(wR)-log(R0))/sqrt(V)
	return(Z)
}

#### Expected Numbers - Recency Assay ####
rec_num <- function(N,para){
	lambda0 = para$lambda0; p = para$p; q = para$q
	OmegaT = para$OmegaT; betaT = para$betaT; Time = para$Time
	PR = betaT + lambda0 *(1-p)/p*(OmegaT-betaT*Time)
	ENR = N*p*PR*q
	return(ENR)
}
#### Expected Numbers - Active Arm Trial ####
inc_num <- function(N,lambda1,para){
	p = para$p; r = para$r; tau = para$tau
	Event = N * tau * (1-p) * r * lambda1
	return(Event)	
}

#### Expected Numbers - Summary ####
summarize_expectnum <-function(N, para, R){
	p = para$p; q = para$q; lambda0 = para$lambda0; MDRI = para$OmegaT; FRR = para$betaT; Time = para$Time;
	r = para$r; tau = para$tau
	return(c(N,N*p*q,rec_num(N,para),N*(1-p)*r,inc_num(N,lambda0*R,para)))	
}

rec_gamma<-function(para, N=NULL){
	lambda0 = para$lambda0; p = para$p; q = para$q; OmegaT = para$OmegaT; betaT = para$betaT; Time = para$Time
	sigma2Omega = para$sigmaOmegaT^2;sigma2beta = para$sigmabetaT^2
	PR = betaT + lambda0 *(1-p)/p*(OmegaT-betaT*Time)
	gamma00 = (PR*(1-PR)/(PR-betaT)^2/q + 1/(1-p) + (1-p*q)*sigma2beta/(PR-betaT)^2/q)/p
	gamma01 = sigma2Omega/(OmegaT-betaT*Time)^2 + sigma2beta*((OmegaT-PR*Time)/(PR-betaT)/(OmegaT-betaT*Time))^2
	return(list(gamma00 = gamma00, gamma01 = gamma01, varlog = gamma00/N+gamma01))
}

trial_gamma<-function(para,lambda1, N = NULL){
	p = para$p; r = para$r; tau = para$tau
	gamma1 = 1 / (lambda1*(1-p)*r*tau)
	return(list(gamma1 = gamma1, varlog = gamma1/N))
}

#### NR - N+test*wbetaT, N+test, N+, wOmega -wbeta*T, Nevent, N-enroll, NR
varW <-function(p,q,PR,OmegaT,betaT,r,lambda1,tau){
	variance = matrix(0,7,7)
	variance[1,1] = p*q*PR*(1-PR) + p*q*(1-p*q)*(PR-betaT)^2
	variance[2,2] = p*q*(1-p*q)
	variance[3,3] = p*(1-p)
	variance[5,5] = (1-p)*r*lambda1*tau*(1+lambda1*tau*(1-r)+lambda1*tau*r*p)
	variance[6,6] = (1-p)*r*(1-r+p*r)
	variance[7,7] = p*q*PR*(1-PR*q*p)
	variance[1,2] = variance[2,1] = p*q*(1-p*q)*(PR-betaT)
	variance[1,3] = variance[3,1] = p*(1-p)*q*(PR-betaT)
	variance[1,5] = variance[5,1] = -p*(1-p)*(PR-betaT)*q*r*lambda1*tau
	variance[1,6] = variance[6,1] = -p*(1-p)*(PR-betaT)*q*r
	variance[1,7] = variance[7,1] = p*q*PR*(1-PR) + p*q*(1-p*q)*(PR-betaT)*PR
	variance[2,3] = variance[3,2] = p*(1-p)*q
	variance[2,5] = variance[5,2] = -p*(1-p)*q*r*lambda1*tau
	variance[2,6] = variance[6,2] = -p*(1-p)*q*r
	variance[2,7] = variance[7,2] = p*q*(1-p*q)*PR
	variance[3,5] = variance[5,3] = -p*(1-p)*r*lambda1*tau
	variance[3,6] = variance[6,3] = -p*(1-p)*r
	variance[3,7] = variance[7,3] = p*(1-p)*q*PR
	variance[5,6] = variance[6,5] = (1-p)*r*(1-r+p*r)*lambda1*tau
	variance[5,7] = variance[7,5] = -p*(1-p)*PR*q*r*lambda1*tau
	variance[6,7] = variance[7,6] = -p*(1-p)*PR*q*r
	return(variance)
}
#### NR - N+test*wbetaT, N+, Nevent, N-enroll, NR
varW5 <-function(p,PR,OmegaT,betaT,r,lambda1,tau){
	variance = matrix(0,5,5)
	variance[1,1] = p*PR*(1-PR) + p*(1-p)*(PR-betaT)^2
	variance[2,2] = p*(1-p)
	variance[3,3] = (1-p)*r*lambda1*tau*(1+lambda1*tau*(1-r)+lambda1*tau*r*p)
	variance[4,4] = (1-p)*r*(1-r+p*r)
	variance[5,5] = p*PR*(1-PR*p)
	variance[1,2] = variance[2,1] = p*(1-p)*(PR-betaT)
	variance[1,3] = variance[3,1] = -p*(1-p)*(PR-betaT)*r*lambda1*tau
	variance[1,4] = variance[4,1] = -p*(1-p)*(PR-betaT)*r
	variance[1,5] = variance[5,1] = p*PR*(1-PR) + p*(1-p)*(PR-betaT)*PR
	variance[2,3] = variance[3,2] = -p*(1-p)*r*lambda1*tau
	variance[2,4] = variance[4,2] = -p*(1-p)*r
	variance[2,5] = variance[5,2] = p*(1-p)*PR
	variance[3,4] = variance[4,3] = (1-p)*r*(1-r+p*r)*lambda1*tau
	variance[3,5] = variance[5,3] = -p*(1-p)*PR*r*lambda1*tau
	variance[4,5] = variance[5,4] = -p*(1-p)*PR*r
	return(variance)
}

#### NR - N+test*wbetaT, N+, Nevent, N-enroll, NR, wOmegaT-wbetaT
varW6 <-function(p,PR,OmegaT,betaT,r,lambda1,tau,sigmaOmegaT, sigmabetaT,Time,N){
	variance = matrix(0,6,6)
	variance[1,1] = p*PR*(1-PR) + p*(1-p)*(PR-betaT)^2
	variance[2,2] = p*(1-p)
	variance[3,3] = (1-p)*r*lambda1*tau*(1+lambda1*tau*(1-r)+lambda1*tau*r*p)
	variance[4,4] = (1-p)*r*(1-r+p*r)
	variance[5,5] = p*PR*(1-PR*p)
	variance[1,2] = variance[2,1] = p*(1-p)*(PR-betaT)
	variance[1,3] = variance[3,1] = -p*(1-p)*(PR-betaT)*r*lambda1*tau
	variance[1,4] = variance[4,1] = -p*(1-p)*(PR-betaT)*r
	variance[1,5] = variance[5,1] = p*PR*(1-PR) + p*(1-p)*(PR-betaT)*PR
	variance[2,3] = variance[3,2] = -p*(1-p)*r*lambda1*tau
	variance[2,4] = variance[4,2] = -p*(1-p)*r
	variance[2,5] = variance[5,2] = p*(1-p)*PR
	variance[3,4] = variance[4,3] = (1-p)*r*(1-r+p*r)*lambda1*tau
	variance[3,5] = variance[5,3] = -p*(1-p)*PR*r*lambda1*tau
	variance[4,5] = variance[5,4] = -p*(1-p)*PR*r
	variance[1,6] = variance[6,1] = p*sigmabetaT^2*Time
	variance = variance * N
	variance[6,6] = sigmaOmegaT^2 + sigmabetaT^2*Time^2
	return(variance)
}

var_AB <- function(p,q,PR,OmegaT, betaT,r,lambda1,tau){
	a = matrix(0,7,2)
	a[1,1] = -1/(p*q*(PR-betaT))
	a[2,1] = 1/(p*q)
	a[3,1] = -1/(p*(1-p))
	a[4,1] = 1/OmegaT
	a[5,1] = 1/((1-p)*r*lambda1*tau)
	a[6,1] = -1/((1-p)*r)
	a[1,2] = -2*PR*(1-PR)/(p*q)^2/(PR-betaT)^3
	a[2,2] = PR^2/(p*q*(PR-betaT))^2
	a[3,2] = -(p^2+(1-p)^2)/(p*(1-p))^2
	a[4,2] = 0
	a[5,2] = -1/((1-p)*r*lambda1*tau)^2
	a[6,2] = 0
	a[7,2] = (1-2*PR)/(p*q*(PR-betaT))^2
	cov_mat = varW(p,q,PR,OmegaT, betaT, r,lambda1,tau)
	return(t(a)%*%cov_mat%*%a)
}


var_AB_sigmabeta <- function(p,q,PR,OmegaT, betaT,r,lambda1,tau,sigmabetaT){
	a = matrix(0,7,2)
	a[1,1] = -1/(p*q*(PR-betaT))
	a[2,1] = 1/(p*q)
	a[3,1] = -1/(p*(1-p))
	a[4,1] = 1/OmegaT
	a[5,1] = 1/((1-p)*r*lambda1*tau)
	a[6,1] = -1/((1-p)*r)
	a[1,2] = -2*PR*(1-PR)/(p*q)^2/(PR-betaT)^3-2*sigmabetaT^2*(1-p)/(p*q)^2/(PR-betaT)^3
	a[2,2] = (PR^2+sigmabetaT^2*(1-2*p))/(p*q*(PR-betaT))^2
	a[3,2] = -(p^2+(1-p)^2)/(p*(1-p))^2
	a[4,2] = 0
	a[5,2] = -1/((1-p)*r*lambda1*tau)^2
	a[6,2] = 0
	a[7,2] = (1-2*PR)/(p*q*(PR-betaT))^2
	cov_mat = varW(p,q,PR,OmegaT, betaT, r,lambda1,tau)
	return(t(a)%*%cov_mat%*%a)
}

var_Z_sigmaOmega <- function(p,PR,OmegaT, betaT,r,lambda1,tau,sigmaOmegaT,Time,N){
	V = varW5(p,PR,OmegaT, betaT, r,lambda1,tau)*N
	a = matrix(0,5,2)
	a[1,1] = -1/(p*(PR-betaT))
	a[2,1] = -1/(1-p)
	a[3,1] = 1/((1-p)*r*lambda1*tau)
	a[4,1] = -1/((1-p)*r)
	
	a[1,2] = -2*PR*(1-PR)/p^2/(PR-betaT)^3
	a[2,2] = (PR^2)/(p*(PR-betaT))^2 -(p^2+(1-p)^2)/(p*(1-p))^2
	a[3,2] = -1/((1-p)*r*lambda1*tau)^2
	a[5,2] = (1-2*PR)/(p*(PR-betaT))^2
	
	a[,1] = a[,1]/N; a[,2] = a[,2]/N^2
	mat = t(a)%*%V%*%a
	lastterm = c(1/(OmegaT-betaT*Time),-2*sigmaOmegaT^2/(OmegaT-betaT*Time)^3)
	mat = mat + lastterm%*%t(lastterm)*sigmaOmegaT^2
	tA = log(R1)-log(R0)
	tB = (PR*(1-PR)/(p*q*(PR-betaT)^2) + 1/p + 1/(1-p) + 1/((1-p)*r*lambda1*tau))/N + sigmaOmegaT^2/(OmegaT - betaT*Time)^2 
	c = rep(0,2); c[1] = 1/sqrt(tB); c[2] = -tA/(2*tB^(3/2))
	return(t(c)%*%mat%*%c)
}

var_Z_both <- function(p,PR,OmegaT, betaT,r,lambda1,tau,sigmaOmegaT,sigmabetaT,Time,N){
	V = varW6(p,PR,OmegaT, betaT, r,lambda1,tau,sigmaOmegaT,sigmabetaT,Time,N)
	a = matrix(0,6,2)
	a[1,1] = -1/(p*(PR-betaT)) / N
	a[2,1] = -1/(1-p) / N
	a[3,1] = 1/((1-p)*r*lambda1*tau) / N
	a[4,1] = -1/((1-p)*r) / N
	a[6,1] = 1/(OmegaT-betaT*Time)
	
	a[1,2] = -2*(PR*(1-PR)+sigmabetaT^2*(1-p))/p^2/(PR-betaT)^3/N^2 - 2*sigmabetaT^2*(OmegaT-PR*Time)/p/(PR-betaT)^3/(OmegaT-betaT*Time)/N
	a[2,2] = ((PR^2+sigmabetaT^2*(1-2*p))/(p*(PR-betaT))^2 -(p^2+(1-p)^2)/(p*(1-p))^2)/N^2 + sigmabetaT^2*(OmegaT-PR*Time)/p/(PR-betaT)^2/(OmegaT-betaT*Time)/N
	a[3,2] = -1/((1-p)*r*lambda1*tau)^2/N^2
	a[5,2] = (1-2*PR)/(p*(PR-betaT))^2/N^2
	a[6,2] = -2*sigmaOmegaT^2/(OmegaT-betaT*Time)^2 - 2*sigmabetaT^2*Time*(OmegaT-PR*Time)/(PR-betaT)/(OmegaT-betaT*Time)^3
	
	mat = t(a)%*%V%*%a
	tA = log(R1)-log(R0)
	tB = ((PR*(1-PR)+sigmabetaT^2*(1-p))/(p*q*(PR-betaT)^2) + 1/p + 1/(1-p) + 1/((1-p)*r*lambda1*tau))/N + sigmaOmegaT^2/(OmegaT - betaT*Time)^2 + sigmabetaT^2*((OmegaT-PR*Time)/(PR-betaT)/(OmegaT-betaT*Time))^2
	c = rep(0,2); c[1] = 1/sqrt(tB); c[2] = -tA/(2*tB^(3/2))
	return(t(c)%*%mat%*%c)
}


varZ <-function(para, R1, R0 = 1){
	p=para$p; q= para$q; OmegaT = para$OmegaT;betaT=para$betaT; r=para$r;tau =para$tau; Time = para$Time
	lambda0 = para$lambda0; lambda1 = lambda0*R1
	PR = betaT + lambda0 *(1-p)/p*(OmegaT-betaT*Time)
	varAB = var_AB(p,q,PR,OmegaT,betaT,r,lambda1,tau)
	tA = log(R1)-log(R0)
	tB = PR*(1-PR)/(p*q*(PR-betaT)^2) + 1/p + 1/(1-p) + 1/((1-p)*r*lambda1*tau)
	c = rep(0,2); c[1] = 1/sqrt(tB); c[2] = -tA/(2*tB^(3/2))
	return(t(c)%*%varAB%*%c)
}

varZ_sigmabeta <-function(para, R1, R0 = 1, N = N){
	p=para$p; q= para$q; OmegaT = para$OmegaT;betaT=para$betaT; r=para$r;tau =para$tau; Time = para$Time
	sigmabetaT = para$sigmabetaT
	lambda0 = para$lambda0; lambda1 = lambda0*R1
	PR = betaT + lambda0 *(1-p)/p*(OmegaT-betaT*Time)
	varAB = var_AB_sigmabeta(p,q,PR,OmegaT,betaT,r,lambda1,tau,sigmabetaT)
	tA = log(R1)-log(R0)
	tB = (PR*(1-PR) + sigmabetaT^2*(1-p))/(p*q*(PR-betaT)^2) + 1/p + 1/(1-p) + 1/((1-p)*r*lambda1*tau)
	c = rep(0,2); c[1] = 1/sqrt(tB); c[2] = -tA/(2*tB^(3/2))
	return(t(c)%*%varAB%*%c)
}

varZ_sigmaOmega <-function(para, R1, R0 = 1, N = N){
	p=para$p; q= para$q; OmegaT = para$OmegaT;betaT=para$betaT; r=para$r;tau =para$tau; Time = para$Time
	sigmaOmegaT = para$sigmaOmegaT
	lambda0 = para$lambda0; lambda1 = lambda0*R1
	PR = betaT + lambda0 *(1-p)/p*(OmegaT-betaT*Time)
	return(var_Z_sigmaOmega(p,PR,OmegaT, betaT,r,lambda1,tau,sigmaOmegaT,Time,N))
}
varZ_both <-function(para, R1, R0 = 1, N = N){
	p=para$p; q= para$q; OmegaT = para$OmegaT;betaT=para$betaT; r=para$r;tau =para$tau; Time = para$Time
	sigmaOmegaT = para$sigmaOmegaT; sigmabetaT = para$sigmabetaT
	lambda0 = para$lambda0; lambda1 = lambda0*R1
	PR = betaT + lambda0 *(1-p)/p*(OmegaT-betaT*Time)
	return(var_Z_both(p,PR,OmegaT, betaT,r,lambda1,tau,sigmaOmegaT,sigmabetaT,Time,N))
}
