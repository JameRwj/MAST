    SUBROUTINE SLAU(variableL,variableR,nx,ny,flux)
    USE BasicData,ONLY:p2,zero,half,one,two,gamma

    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: variableL,variableR,flux
    REAL(p2),DIMENSION(4) :: FluxL,FluxR,flux_plus,flux_minus
    REAL(p2) :: rhoL,uL,vL,pL,cL,EL,hL,qnL
    REAL(p2) :: rhoR,uR,vR,pR,cR,ER,hR,qnR
    REAL(p2):: cFace,delP,delRho,maL,maR,alphaL,alphaR,betaL,betaR,vtFace,Mcap,Xi,Vnabs,VnabsL,VnabsR,fnG,pBar,mass

    REAL(p2) nx,ny

    rhoL = variableL(1)
    uL   = variableL(2)
    vL   = variableL(3)
    pL   = variableL(4)

    rhoR = variableR(1)
    uR   = variableR(2)
    vR   = variableR(3)
    pR   = variableR(4)

    cL=sqrt(gamma*pL/rhoL)
    cR=sqrt(gamma*pR/rhoR)
    cFace=half*(cL+cR)
    hL=(gamma/(gamma-one))*pL/rhoL+(uL*uL+vL*vL)/two
    hR=(gamma/(gamma-one))*pR/rhoR+(uR*uR+vR*vR)/two
    delP=pR-pL
    delRho=rhoR-rhoL

    qnL=uL*nx+vL*ny
    qnR=uR*nx+vR*ny

    maL=qnL/cFace
    maR=qnR/cFace

    alphaL=max(zero, one-floor(abs(maL)))
    alphaR=max(zero, one-floor(abs(maR)))

    betaL = (one-alphaL)*half*(one+sign(one,maL)) + (alphaL)*0.25*(two-maL)*((maL+one)**2)
    betaR = (one-alphaR)*half*(one-sign(one,maR)) + (alphaR)*0.25*(two+maR)*((maR-one)**2)

    vtface = sqrt(half*((uL*uL) + (vL*vL) + (uR*uR) + (vR*vR)))
    Mcap   = min(one, vtface/cFace)
    Xi     = (one - Mcap)**2

    Vnabs = (rhoL *abs(qnL) + rhoR*abs(qnR))/(rhoL + rhoR)

    fnG = -one*max(min(maL,zero),-one)*min(max(maR,zero),one)

    pbar = half*((pL+pR) + (betaL-betaR)*(pL-pR) + (one-xi)*(betaL+betaR-one)*(pL+pR))

    VnabsL = (one - fnG)*Vnabs + fnG*abs(qnL)
    VnabsR = (one - fnG)*Vnabs + fnG*abs(qnR)

    mass = half*((rhoL*(qnL+VnabsL) + rhoR*(qnR-VnabsR)) - (Xi*delp/cFace))

    flux_plus(1)=half*(mass + abs(mass))
    flux_plus(2)=flux_plus(1)*uL
    flux_plus(3)=flux_plus(1)*vL
    flux_plus(4)=flux_plus(1)*hL

    flux_minus(1)=half*(mass - abs(mass))
    flux_minus(2)=flux_minus(1)*uR
    flux_minus(3)=flux_minus(1)*vR
    flux_minus(4)=flux_minus(1)*hR


    Flux(:)=flux_plus+flux_minus
    Flux(2)=Flux(2)+pbar*nx
    Flux(3)=Flux(3)+pbar*ny
    
    END SUBROUTINE SLAU
