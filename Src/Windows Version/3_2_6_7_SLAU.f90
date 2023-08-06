SUBROUTINE SLAU_1D()
    
    USE BasicData,ONLY:gamma,id,one,two,half,p2,zero
    USE BasicData_1D,ONLY:rho_1D,u_1D,p_1D,rhs
	
	IMPLICIT NONE

    REAL(p2):: delta_plus,delta_minus,delta

    INTEGER i
	REAL(p2):: rhoL,uL,pL,hL,cL,EL
	REAL(p2):: rhoR,uR,pR,hR,cR,ER
    REAL(p2):: cFace,delP,delRho,maL,maR,alphaL,alphaR,betaL,betaR,vtFace,Mcap,Xi,Vnabs,VnabsL,VnabsR,fnG,pBar,mass

    REAL(p2),DIMENSION(:,:),ALLOCATABLE :: xFlux

    REAL(p2),DIMENSION(:),ALLOCATABLE :: flux_plus,flux_minus
	
	ALLOCATE(xFlux(3,id))
    ALLOCATE(flux_plus(3),flux_minus(3))


    	DO i=1,id
!Step 1:	Define the left and right states of cell interface with MUSCL.
		!Left state
		CALL Reconstruction_L(i,rhoL,rho_1D)

        CALL Reconstruction_L(i,uL,u_1D)

        CALL Reconstruction_L(i,pL,p_1D)
        
        !		Right state
        CALL Reconstruction_R(i,rhoR,rho_1D)
        
        CALL Reconstruction_R(i,uR,u_1D)

        CALL Reconstruction_R(i,pR,p_1D)
        
 !Step 2:	Calculate the numerical flux
        EL=pL/(gamma-one)+half*rhoL*uL*uL
        ER=pR/(gamma-one)+half*rhoR*uR*uR
	
		cL=SQRT(gamma*pL/rhoL)
		cR=SQRT(gamma*pR/rhoR)
        cFace=half*(cL+cR)

        hL=(gamma/(gamma-one))*pL/rhoL+(uL*uL)/two
		hR=(gamma/(gamma-one))*pR/rhoR+(uR*uR)/two
        
        delP=pR-pL
        delRho=rhoR-rhoL
        
        maL=uL/cFace
        maR=uR/cFace
            
        alphaL=MAX(zero, one-FLOOR(ABS(maL)))
        alphaR=MAX(zero, one-FLOOR(ABS(maR)))
        
        betaL = (one-alphaL)*half*(one+SIGN(one,maL)) + (alphaL)*0.25*(two-maL)*((maL+one)**2)
        betaR = (one-alphaR)*half*(one-SIGN(one,maR)) + (alphaR)*0.25*(two+maR)*((maR-one)**2)
        
        vtface = SQRT(half*((uL*uL)+ (uR*uR)))
        Mcap   = MIN(one, vtface/cFace)
        Xi     = (one - Mcap)**2
            
        Vnabs = (rhoL *ABS(uL) + rhoR*ABS(uR))/(rhoL + rhoR)
        
        fnG = -one*MAX(MIN(maL,zero),-one)*MIN(MAX(maR,zero),one)
        
        pbar = half*((pL+pR) + (betaL-betaR)*(pL-pR) + (one-xi)*(betaL+betaR-one)*(pL+pR))
        
        VnabsL = (one - fnG)*Vnabs + fnG*ABS(uL)
        !VnabsL = (one - fnG)*Vnabs - fnG*abs(uL)
        VnabsR = (one - fnG)*Vnabs + fnG*ABS(uR)
        
        mass = half*((rhoL*(uL+VnabsL) + rhoR*(uR-VnabsR)) - (Xi*delp/cFace))

        flux_plus(1)=half*(mass + ABS(mass))
        flux_plus(2)=flux_plus(1)*uL
        flux_plus(3)=flux_plus(1)*hL

        flux_minus(1)=half*(mass - ABS(mass))
        flux_minus(2)=flux_minus(1)*uR
        flux_minus(3)=flux_minus(1)*hR

        xFlux(:,i)=flux_plus+flux_minus 
        xFlux(2,i)=xFlux(2,i)+pbar

        END DO

        DO i=1,id-1
		    rhs(:,i)=rhs(:,i)-(xFlux(:,i+1)-xFlux(:,i))
	    END DO

        DEALLOCATE(xFlux)

END SUBROUTINE SLAU_1D