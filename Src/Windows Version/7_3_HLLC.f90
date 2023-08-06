    SUBROUTINE HLLC(variableL,variableR,nx,ny,flux)
    USE BasicData,ONLY:p2,zero,half,one,two,gamma
    
    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: variableL,variableR,flux
    REAL(p2),DIMENSION(4) :: FluxL,FluxR,FluxStarL,FluxStarR
    REAL(p2) :: rhoL,uL,vL,pL,cL,EL,hL,qnL
    REAL(p2) :: rhoR,uR,vR,pR,cR,ER,hR,qnR
    REAL(p2) cFace,uFace,vFace,hFace,qnFace,rhoFace
    REAL(p2) rhoStarL,rhoStarR,eStarL,eStarR

    REAL(p2) sL,sR,sStar,pStar
    REAL(p2) nx,ny
    REAL(p2) alfaL,alfaR


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

    rhoFace=SQRT(rhoL*rhoR)
    uFace=(SQRT(rhoL)*uL+SQRT(rhoR)*uR)/(SQRT(rhoL)+SQRT(rhoR))
    vFace=(SQRT(rhoL)*vL+SQRT(rhoR)*vR)/(SQRT(rhoL)+SQRT(rhoR))
    hFace=(SQRT(rhoL)*hL+SQRT(rhoR)*hR)/(SQRT(rhoL)+SQRT(rhoR))
    cFace=SQRT( (gamma-one)*(hFace-half*(uFace**2+vFace**2)) )

    qnL=uL*nx+vL*ny
    qnR=uR*nx+vR*ny
    qnFace=(SQRT(rhoL)*qnL+SQRT(rhoR)*qnR)/(SQRT(rhoL)+SQRT(rhoR))

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
    !sL=min(qnFace-cFace,qnL-cL)
    !sR=max(qnFace+cFace,qnR+cR)

    alfaL=rhoL*(sL-qnL)
    alfaR=rhoR*(sR-qnR)

    sStar=(alfaR*qnR-alfaL*qnL+pL-pR)/(alfaR-alfaL)

    rhoStarL=alfaL/(sL-sStar)
    rhoStarR=alfaR/(sR-sStar)

    eStarL=eL+(sStar-qnL)*(sStar+pL/alfaL)
    eStarR=eR+(sStar-qnR)*(sStar+pR/alfaR)

    pStar=(alfaR*pL-alfaL*pR-alfaL*alfaR*(qnL-qnR))/(alfaR-alfaL)

    FluxStarL(1)=rhoStarL*sStar
    FluxStarL(2)=rhoStarL*sStar*(uL+nx*(sStar-qnL))+nx*pStar
    FluxStarL(3)=rhoStarL*sStar*(vL+ny*(sStar-qnL))+ny*pStar
    FluxStarL(4)=sStar*(rhoStarL*eStarL+pStar)

    FluxStarR(1)=rhoStarR*sStar
    FluxStarR(2)=rhoStarR*sStar*(uR+nx*(sStar-qnR))+nx*pStar
    FluxStarR(3)=rhoStarR*sStar*(vR+ny*(sStar-qnR))+ny*pStar
    FluxStarR(4)=sStar*(rhoStarR*eStarR+pStar)

    IF(sL .GE. zero)THEN
        Flux(:)=FluxL
    ENDIF

    IF( (sL .LE. zero) .and. (sStar .GE. zero))THEN
        Flux(:)=FluxStarL
    ENDIF

    IF( (sStar .LE. zero) .and. (sR .GE. zero))THEN
        Flux(:)=FluxStarR
    ENDIF

    IF( sR .LE. zero)THEN
        Flux(:)=FluxR
    END IF

    END SUBROUTINE HLLC
