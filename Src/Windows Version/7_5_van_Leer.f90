
    SUBROUTINE van_Leer(variableL,variableR,nx,ny,flux)
    USE BasicData,ONLY:p2,zero,half,one,two,gamma

    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: variableL,variableR,flux
    REAL(p2),DIMENSION(4) :: FluxL,FluxR,qL,qR
    REAL(p2) :: rhoL,uL,vL,pL,cL,EL,hL,qnL,maL
    REAL(p2) :: rhoR,uR,vR,pR,cR,ER,hR,qnR,maR

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

    maL=qnL/cL
    maR=qnR/cR

    IF(maL .GE. one)THEN
        fluxL(1)=rhoL*qnL
        fluxL(2)=rhoL*qnL*uL+pL*nx
        fluxL(3)=rhoL*qnL*vL+pL*ny
        fluxL(4)=(pL*gamma/(gamma-one)+rhoL*(uL*uL+vL*vL)*half)*qnL
    ELSE IF(maL .LE. -one)THEN
        fluxL(1)=zero
        fluxL(2)=zero
        fluxL(3)=zero
        fluxL(4)=zero
    ELSE
        fluxL(1)=0.25*rhoL*cL*(MaL+one)*(MaL+one)
        fluxL(2)=fluxL(1)*(nx*(-qnL+two*cL)/gamma+uL)
        fluxL(3)=fluxL(1)*(ny*(-qnL+two*cL)/gamma+vL)
        fluxL(4)=fluxL(1)*half*(((gamma-one)*qnL+two*cL)**2/(gamma*gamma-one)&
            +(uL*uL+vL*vL-qnL*qnL))
    END IF

    IF(maR .GE. one)THEN
        fluxR(1)=zero
        fluxR(2)=zero
        fluxR(3)=zero
        fluxR(4)=zero
    ELSE IF(maR .LE. -one)THEN
        fluxR(1)=rhoR*qnR
        fluxR(2)=rhoR*qnR*uR+pR*nx
        fluxR(3)=rhoR*qnR*vR+pR*ny
        fluxR(4)=(pR*gamma/(gamma-one)+rhoR*(uR*uR+vR*vR)*half)*qnR
    ELSE
        fluxR(1)=-0.25*rhoR*cR*(MaR-one)*(MaR-one)
        fluxR(2)=fluxR(1)*(nx*(-qnR-two*cR)/gamma+uR)
        fluxR(3)=fluxR(1)*(ny*(-qnR-two*cR)/gamma+vR)
        fluxR(4)=fluxR(1)*half*(((gamma-one)*qnR-two*cR)**2/(gamma*gamma-one)+(uR*uR+vR*vR-qnR*qnR))
    END IF

    Flux(:)=fluxL(:)+fluxR(:)

    END SUBROUTINE van_Leer