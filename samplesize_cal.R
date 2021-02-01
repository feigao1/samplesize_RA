source('recencyassay_code.R')
source('setting_code.R')
tau_set = c(1,2); R0 = 0.5; R1 = 0.15
#### MSM ####
for (tau in tau_set){
	para = SettingMSM_para(q=1, r=0.85, tau=tau, N=NULL, alpha=0.05, beta=0.9)
	N  = samplesize_CF(para, R1 = R1, R0 = R0)
	N_crit = summarize_expectnum(N,para,R1) # N, N+, Nrec, N-enroll, Nevent
	print(c(tau,N_crit))
}

##### Simulations ######
nrep = 10000;  set.seed(123)
R0 = 0.5; R1 = 0.15; tau_set = c(1,2)
for (tau in tau_set){
	para = SettingMSM_para(q=1, r=0.85, tau=tau, N=NULL, alpha=0.05, beta=0.9)
	N  = samplesize_CF(para, R1 = R1, R0 = R0)
	Rej0 = Rej0_I = Rej1 = Rej1_I = lambda0 = lambda1_null = lambda1_alt = NULL
	for (rep in 1:nrep){
		data = generatedata(N,para,R=R0)
		para_est = para; para_est$OmegaT = data$wOmegaT;  para_est$betaT = data$wbetaT
		lambda0_est = rec_est(N,data$Npos,data$Npos_test,data$NR, para_est)
		lambda1_est = trial_est(N,data$Nneg,data$Nneg_trial,data$Nevent,para$tau)
		Rej0 = c(Rej0,Rejection(lambda0_est,lambda1_est,R0=R0))
						
		data = generatedata(N,para,R=R1)
		para_est = para; para_est$OmegaT = data$wOmegaT;  para_est$betaT = data$wbetaT
		lambda0_est = rec_est(N,data$Npos,data$Npos_test,data$NR, para_est)
		lambda1_est = trial_est(N,data$Nneg,data$Nneg_trial,data$Nevent,para$tau)
		Rej1 = c(Rej1,Rejection(lambda0_est,lambda1_est,R0=R0))
	}
	print(c(tau,N,mean(Rej0),mean(Rej1)))
}

#### Women Trial ####
for (tau in tau_set){
	para = SettingWomen_para(q=1, r=0.85, tau=tau, N=NULL, alpha=0.05, beta=0.9)
	N  = samplesize_CF(para, R1 = R1, R0 = R0)
	N_crit = summarize_expectnum(N,para,R1) # N, N+, Nrec, N-enroll, Nevent
	print(c(tau,N_crit))
}

##### Simulations ######
nrep = 10000;  set.seed(123)
R0 = 0.5; R1 = 0.15; tau_set = c(1,2)
for (tau in tau_set){
	para = SettingWomen_para(q=1, r=0.85, tau=tau, N=NULL, alpha=0.05, beta=0.9)
	N  = samplesize_CF(para, R1 = R1, R0 = R0)
	Rej0 = Rej0_I = Rej1 = Rej1_I = lambda0 = lambda1_null = lambda1_alt = NULL
	for (rep in 1:nrep){
		data = generatedata(N,para,R=R0)
		para_est = para; para_est$OmegaT = data$wOmegaT;  para_est$betaT = data$wbetaT
		lambda0_est = rec_est(N,data$Npos,data$Npos_test,data$NR, para_est)
		lambda1_est = trial_est(N,data$Nneg,data$Nneg_trial,data$Nevent,para$tau)
		Rej0 = c(Rej0,Rejection(lambda0_est,lambda1_est,R0=R0))
						
		data = generatedata(N,para,R=R1)
		para_est = para; para_est$OmegaT = data$wOmegaT;  para_est$betaT = data$wbetaT
		lambda0_est = rec_est(N,data$Npos,data$Npos_test,data$NR, para_est)
		lambda1_est = trial_est(N,data$Nneg,data$Nneg_trial,data$Nevent,para$tau)
		Rej1 = c(Rej1,Rejection(lambda0_est,lambda1_est,R0=R0))
	}
	print(c(tau,N,mean(Rej0),mean(Rej1)))
}