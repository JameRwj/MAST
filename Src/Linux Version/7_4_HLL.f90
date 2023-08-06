    SUBROUTINE HLL(variableL,variableR,nx,ny,flux)
    USE BasicData,ONLY:p2,zero,half,one,two,gamma
    
    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: variableL,variableR,flux
    REAL(p2),DIMENSION(4) :: FluxL,FluxR,qL,qR
    REAL(p2) :: rhoL,uL,vL,pL,cL,EL,hL,qnL
    REAL(p2) :: rhoR,uR,vR,pR,cR,ER,hR,qnR

    REAL(p2) sL,sR
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

    eL=pL/(gamma-one)/rhoL+half*(uL**2+vL**2)
    eR=pR/(gamma-one)/rhoR+half*(uR**2+vR**2)

    HL=gamma/(gamma-one)*pL/rhoL+half*(uL**2+vL**2)
    HR=gamma/(gamma-one)*pR/rhoR+half*(uR**2+vR**2)

    qnL=uL*nx+vL*ny
    qnR=uR*nx+vR*ny

    FluxL(1)=rhoL*qnL
    FluxL(2)=rhoL*qnL*uL+nx*pL
    FluxL(3)=rhoL*qnL*vL+ny*pL
    FluxL(4)=rhoL*qnL*hL

    FluxR(1)=rhoR*qnR
    FluxR(2)=rhoR*qnR*uR+nx*pR
    FluxR(3)=rhoR*qnR*vR+ny*pR
    FluxR(4)=rhoR*qnR*hR

    sL=min(qnR-cR,qnL-cL)
    sR=max(qnL+cL,qnR+cR)

    qL(1)=rhoL
    qL(2)=rhoL*uL
    qL(3)=rhoL*vL
    qL(4)=rhoL*eL

    qR(1)=rhoR
    qR(2)=rhoR*uR
    qR(3)=rhoR*vR
    qR(4)=rhoR*eR

    IF(sL .GE. zero)THEN
        Flux(:)=FluxL
    ENDIF

    IF( (sL .LE. zero) .and. (sR .GE. zero))THEN
        Flux(:)=sR/(sR-sL)*FluxL(:)-sL/(sR-sL)*FluxR(:)+sR*sL/(sR-sL)*(qR(:)-qL(:))
    ENDIF

    IF( sR .LE. zero)THEN
        Flux(:)=FluxR
    END IF

    END SUBROUTINE HLL
