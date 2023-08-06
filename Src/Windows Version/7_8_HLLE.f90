    SUBROUTINE HLLE(variableL,variableR,nx,ny,flux)
    USE BasicData,ONLY:p2,zero,half,one,two,gamma

    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: variableL,variableR,flux
    REAL(p2),DIMENSION(4) :: FluxL,FluxR,qL,qR
    REAL(p2) :: rhoL,uL,vL,pL,cL,EL,hL,qnL
    REAL(p2) :: rhoR,uR,vR,pR,cR,ER,hR,qnR
    real(p2) :: rhoFace,cFace,uFace,vFace,hFace,qnFace
    real(p2) :: sL,sR,sStar

    REAL(p2) nx,ny


    rhoL = variableL(1)
    uL   = variableL(2)
    vL   = variableL(3)
    pL   = variableL(4)

    rhoR = variableR(1)
    uR   = variableR(2)
    vR   = variableR(3)
    pR   = variableR(4)

    cL=SQRT(gamma*pL/rhoL)
    cR=SQRT(gamma*pR/rhoR)

    EL=pL/(gamma-one)+half*rhoL*(uL*uL+vL*vL)
    ER=pR/(gamma-one)+half*rhoR*(uR*uR+vR*vR)

    hL=(EL+pL)/rhoL
    hR=(ER+pR)/rhoR

    qL(1)=rhoL
    qL(2)=rhoL*uL
    qL(3)=rhoL*vL
    qL(4)=EL

    qR(1)=rhoR
    qR(2)=rhoR*uR
    qR(3)=rhoR*vR
    qR(4)=ER

    rhoFace=SQRT(rhoL*rhoR)
    uFace=(SQRT(rhoL)*uL+SQRT(rhoR)*uR)/(SQRT(rhoL)+SQRT(rhoR))
    vFace=(SQRT(rhoL)*vL+SQRT(rhoR)*vR)/(SQRT(rhoL)+SQRT(rhoR))
    hFace=(SQRT(rhoL)*hL+SQRT(rhoR)*hR)/(SQRT(rhoL)+SQRT(rhoR))
    cFace=SQRT( (gamma-one)*(hFace-half*(uFace*uFace+vFace*vFace)) )

    qnL=uL*nx+vL*ny
    qnR=uR*nx+vR*ny
    qnFace=(SQRT(rhoL)*qnL+SQRT(rhoR)*qnR)/(SQRT(rhoL)+SQRT(rhoR))

    sL=min(zero,qnL-cL,qnFace-cFace)
    sR=max(zero,qnR+cR,qnFace+cFace)

    FluxL(1)=rhoL*qnL
    FluxL(2)=rhoL*qnL*uL+nx*pL
    FluxL(3)=rhoL*qnL*vL+ny*pL
    FluxL(4)=rhoL*qnL*hL

    FluxR(1)=rhoR*qnR
    FluxR(2)=rhoR*qnR*uR+nx*pR
    FluxR(3)=rhoR*qnR*vR+ny*pR
    FluxR(4)=rhoR*qnR*hR

    Flux(:)=sR*FluxL/(sR-sL)-sL*FluxR/(sR-sL)+sL*sR*(qR-qL)/(sR-sL)

    END SUBROUTINE HLLE
