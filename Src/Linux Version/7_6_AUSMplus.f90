    SUBROUTINE AUSMplus(variableL,variableR,nx,ny,flux)
    USE BasicData,ONLY:p2,zero,half,one,two,gamma

    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: variableL,variableR,flux
    REAL(p2):: rhoL,uL,vL,pL,hL,cL,maL,qnL
    REAL(p2):: rhoR,uR,vR,pR,hR,cR,maR,qnR
    REAL(p2):: cCriticalL,cCriticalR,cFace,maFace
    REAL(p2):: maPlus,maMinus
    REAL(p2):: pPlus,pMinus

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

    hL=(gamma/(gamma-one))*pL/rhoL+uL*uL/two+vL*vL/two
    hR=(gamma/(gamma-one))*pR/rhoR+uR*uR/two+vR*vR/two

    qnL=uL*nx+vL*ny
    qnR=uR*nx+vR*ny

    !------------------------------------------------------------------------------------------------------------
    cCriticalL=SQRT(two*(gamma-one)/(gamma+one)*hL)
    cCriticalR=SQRT(two*(gamma-one)/(gamma+one)*hR)

    cL=(cCriticalL*cCriticalL)/MAX(cCriticalL,ABS(qnL))
    cR=(cCriticalR*cCriticalR)/MAX(cCriticalR,ABS(qnR))

    !cFace=MIN(cL,cR)
    
    cFace=half*(cL+cR)

    maL=qnL/cFace
    maR=qnR/cFace

    IF(ABS(maL) .LE. one)THEN
        maPlus=0.25_p2*(maL+one)**2+one/8.0_p2*(maL**2-one)**2
    ELSE
        maPlus=half*(maL+ABS(maL))
    END IF

    IF(ABS(maR) .LE. one)THEN
        maMinus=-0.25_p2*(maR-one)**2-one/8.0_p2*(maR**2-one)**2
    ELSE
        maMinus=half*(maR-ABS(maR))
    END IF

    IF(ABS(maL) .LE. one)THEN
        pPlus=0.25_p2*(maL+one)**2*(two-maL)+3.0_p2/16.0_p2*maL*(maL**2-one)**2
    ELSE
        pPlus=half*(one+SIGN(one,maL))
    END IF

    IF(ABS(maR) .LE. one)THEN
        pMinus=0.25_p2*(maR-one)**2*(two+maR)-3.0_p2/16.0_p2*maR*(maR**2-one)**2
    ELSE
        pMinus=half*(one-SIGN(one,maR))
    END IF

    maFace=maPlus+maMinus

    IF(maFace .GT. zero)THEN
        Flux(1)=maFace*rhoL*cFace
        Flux(2)=maFace*rhoL*uL*cFace+(pPlus*pL+pMinus*pR)*nx
        Flux(3)=maFace*rhoL*vL*cFace+(pPlus*pL+pMinus*pR)*ny
        Flux(4)=maFace*rhoL*hL*cFace
    ELSE
        Flux(1)=maFace*rhoR*cFace
        Flux(2)=maFace*rhoR*uR*cFace+(pPlus*pL+pMinus*pR)*nx
        Flux(3)=maFace*rhoR*vR*cFace+(pPlus*pL+pMinus*pR)*ny
        Flux(4)=maFace*rhoR*hR*cFace
    END IF

    END SUBROUTINE AUSMplus
