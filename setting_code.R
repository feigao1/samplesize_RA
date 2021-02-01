#### Returning Values:
# lambda0, p, q, r, OmegaT, sigmaOmegaT, betaT, sigmabetaT, Time, tau, N, alpha, beta
#### Hypothetical MSM Trial
### Mimicing HPTN083
SettingMSM_para <- function(q, r, tau, N=NULL, alpha=0.05, beta=0.9){
	popnames <- c('US-BMSM', 'US-other','Brazil','Peru','Buenos Aires','CapeTown','Bangkok','ChiangMai','India','Hanoi','TGW')
	Inc <- c(5.9,1.3, 5, 3.5, 6.4, 4.7, 5.2, 8.2, 0.87, 4, 8)/100
	Subtype <- c('B','B','B','B','B','C','A/E','A/E','C','A/E',NA)
	nsub <- c(843,854,796,829,335,152,413,140,0,199,NA)
	prevalence <- c(0.15,0.15,0.15,0.15,0.15,0.25,0.15,0.15,0.25,0.15,NA)
	RA_subtype <- c('A','B','C','D','A/E')
	MDRI <- c(159,142,118,182,NA)/365.25
	RSE_MDRI <- c(0.13,0.10,0.07,0.19,NA)
	FRR <- c(0.003,0.015,0.010,0.039,NA)
	
	pop <- data.frame(popnames, Subtype, Inc, prevalence,nsub)
	RA <- data.frame(RA_subtype, MDRI,RSE_MDRI,FRR)
	colnames(RA)[1] = 'Subtype'
	pop_inc = sum(nsub*Inc,na.rm=T)/sum(nsub,na.rm=T)
	pop_p = sum(nsub*prevalence,na.rm=T)/sum(nsub,na.rm=T)
	
	pop = merge(pop,RA,by='Subtype')
	pop$nsub_NA = pop$nsub; pop$nsub_NA[is.na(pop$MDRI)] = NA 
	pop$MDRI_sigma2 = (pop$MDRI * pop$RSE_MDRI)^2
	pop_MDRI = sum(pop$nsub_NA*pop$MDRI,na.rm=T)/sum(pop$nsub_NA,na.rm=T)
	pop_MDRI_sigma = sqrt(sum(pop$nsub_NA*pop$MDRI_sigma2,na.rm=T)/sum(pop$nsub_NA,na.rm=T))
	pop_FRR = sum(pop$nsub_NA*pop$FRR,na.rm=T)/sum(pop$nsub_NA,na.rm=T)
	
	RSE_FRR = 0.25; sigmaFRR = pop_FRR * RSE_FRR; Time = 2
	return(list(lambda0=pop_inc, p=pop_p, q=q, r=r, OmegaT=pop_MDRI, sigmaOmegaT=pop_MDRI_sigma, betaT=pop_FRR, sigmabetaT=sigmaFRR, Time=Time, tau=tau, N=N, alpha=alpha, beta=beta))
}
#### Hypothetical Women Trial
### Mimicing HPTN084
SettingWomen_para <- function(q, r, tau, N=NULL, alpha=0.05, beta=0.9){
	MDRI = 118/365.25; RSE_MDRI = 0.07; FRR = 0.015; RSE_FRR = 0.25; Time = 2
	lambda0 = 0.035; p = 0.25;
	sigmaMDRI = MDRI * RSE_MDRI; sigmaFRR = FRR * RSE_FRR
	return(list(lambda0=lambda0, p=p, q=q, r=r, OmegaT=MDRI, sigmaOmegaT=sigmaMDRI, betaT=FRR, sigmabetaT=sigmaFRR, Time=Time, tau=tau, N=N, alpha=alpha, beta=beta))
}