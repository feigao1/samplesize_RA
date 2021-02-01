#### Returning Values:
# lambda0, p, q, r, OmegaT, sigmaOmegaT, betaT, sigmabetaT, Time, tau, N, alpha, beta
Setting1_para <- function(q, r, tau, N=NULL, alpha=0.05, beta=0.8){  ##### MOZ #####
	MDRI = 118/365.25; RSE_MDRI = 0.07; FRR = 0.015; RSE_FRR = 0.25; Time = 2
	lambda0 = 0.0089; p = 0.124;
	sigmaMDRI = MDRI * RSE_MDRI; sigmaFRR = FRR * RSE_FRR
	return(list(lambda0=lambda0, p=p, q=q, r=r, OmegaT=MDRI, sigmaOmegaT=sigmaMDRI, betaT=FRR, sigmabetaT=sigmaFRR, Time=Time, tau=tau, N=N, alpha=alpha, beta=beta))
}
Setting2_para <- function(q, r, tau, N=NULL, alpha=0.05, beta=0.8){ ##### SA AGYW 14-17 #####
	MDRI = 118/365.25; RSE_MDRI = 0.07; FRR = 0.015; RSE_FRR = 0.25; Time = 2
	lambda0 = 0.047; p = 0.276;
	sigmaMDRI = MDRI * RSE_MDRI; sigmaFRR = FRR * RSE_FRR
	return(list(lambda0=lambda0, p=p, q=q, r=r, OmegaT=MDRI, sigmaOmegaT=sigmaMDRI, betaT=FRR, sigmabetaT=sigmaFRR, Time=Time, tau=tau, N=N, alpha=alpha, beta=beta))
}
Setting3_para <- function(q, r, tau, N=NULL, alpha=0.05, beta=0.8){ ##### SA AGYW >18 #####
	MDRI = 118/365.25; RSE_MDRI = 0.07; FRR = 0.015; RSE_FRR = 0.25; Time = 2
	lambda0 = 0.069; p = 0.52;
	sigmaMDRI = MDRI * RSE_MDRI; sigmaFRR = FRR * RSE_FRR
	return(list(lambda0=lambda0, p=p, q=q, r=r, OmegaT=MDRI, sigmaOmegaT=sigmaMDRI, betaT=FRR, sigmabetaT=sigmaFRR, Time=Time, tau=tau, N=N, alpha=alpha, beta=beta))
}
Setting4_para <- function(q, r, tau, N=NULL, alpha=0.05, beta=0.8){ ##### USA MSM #####
	MDRI = 142/365.25; RSE_MDRI = 0.1; FRR = 0.01; RSE_FRR = 0.25; Time = 2
	lambda0 = 0.0342; p = 0.145; 
	sigmaMDRI = MDRI * RSE_MDRI; sigmaFRR = FRR * RSE_FRR
	return(list(lambda0=lambda0, p=p, q=q, r=r, OmegaT=MDRI, sigmaOmegaT=sigmaMDRI, betaT=FRR, sigmabetaT=sigmaFRR, Time=Time, tau=tau, N=N, alpha=alpha, beta=beta))
}