    SUBROUTINE HLLC_1D()

    USE BasicData,ONLY:gamma,id,one,two,half,p2,zero
    USE BasicData_1D,ONLY:rho_1D,u_1D,p_1D,rhs

    IMPLICIT NONE

    INTEGER i

    REAL(p2):: delta_plus,delta_minus,delta

    REAL(p2):: rhoL,uL,pL,cL,eL,hL
    REAL(p2):: rhoR,uR,pR,cR,eR,hR
    REAL(p2):: RT,rhoFace,uFace,cFace,hFace
    REAL(p2):: sL,sR,sStar
    REAL(p2):: rhoAvg,cAvg
    REAL(p2):: pStar,pPvrs,qqL,qqR
    REAL(p2):: fluxL(3),fluxR(3),fluxStarL(3),fluxStarR(3)
    REAL(p2):: qL(3),qR(3),qStarL(3),qStarR(3)
    REAL(p2):: xFlux(3,id)

    DO i=1,id
        !Left state
        CALL Reconstruction_L(i,rhoL,rho_1D)

        CALL Reconstruction_L(i,uL,u_1D)

        CALL Reconstruction_L(i,pL,p_1D)
        
        !		Right state
        CALL Reconstruction_R(i,rhoR,rho_1D)
        
        CALL Reconstruction_R(i,uR,u_1D)

        CALL Reconstruction_R(i,pR,p_1D)
        
        cL=SQRT(gamma*pL/rhoL)
        cR=SQRT(gamma*pR/rhoR)

        eL=pL/(gamma-one)+half*rhoL*uL*uL
        eR=pR/(gamma-one)+half*rhoR*uR*uR

        hL=(EL+pL)/rhoL
        hR=(ER+pR)/rhoR

        ! Roe's averaged variables

        hFace=(SQRT(rhoL)*hL+SQRT(rhoR)*hR)/(SQRT(rhoL)+SQRT(rhoR))
        uFace=(SQRT(rhoL)*uL+SQRT(rhoR)*uR)/(SQRT(rhoL)+SQRT(rhoR))
        cFace=SQRT((gamma-one)*(hFace-half*uFace*uFace))

        !calculate the conservative variables

        qL(1)=rhoL
        qL(2)=rhoL*uL
        qL(3)=eL

        qR(1)=rhoR
        qR(2)=rhoR*uR
        qR(3)=eR


        !sL=min(uL-cL,uFace-cFace)
        !sR=max(uR+cR,uFace+cFace)

        sL=min(uL-cL,uR-cR)
        sR=max(uL+cL,uR+cR)

        sStar=( pR-pL+rhoL*uL*(sL-uL)-rhoR*uR*(sR-uR) )/( rhoL*(sL-uL)-rhoR*(sR-uR) )

        !compute the flux

        fluxL(1)=rhoL*uL
        fluxL(2)=rhoL*uL*uL+pL
        fluxL(3)=rhoL*uL*hL

        fluxR(1)=rhoR*uR
        fluxR(2)=rhoR*uR*uR+pR
        fluxR(3)=rhoR*uR*hR

        !calculate the star conservative states

        qStarL(1)=rhoL*(sL-uL)/(sL-sStar)
        qStarL(2)=rhoL*(sL-uL)/(sL-sStar)*sStar
        qStarL(3)=rhoL*(sL-uL)/(sL-sStar)*(eL/rhoL+(sStar-uL)*(sStar+pL/rhoL/(sL-uL)))

        qStarR(1)=rhoR*(sR-uR)/(sR-sStar)
        qStarR(2)=rhoR*(sR-uR)/(sR-sStar)*sStar
        qStarR(3)=rhoR*(sR-uR)/(sR-sStar)*(eR/rhoR+(sStar-uR)*(sStar+pR/rhoR/(sR-uR)))

        fluxStarL(:)=fluxL(:)+sL*(qStarL(:)-qL(:))
        fluxStarR(:)=fluxR(:)+sR*(qStarR(:)-qR(:))

        IF(sL>=0)THEN
            xFlux(:,i)=fluxL(:)
        ELSEIF(sL<=zero .and. sStar>=zero)THEN
            xFlux(:,i)=fluxStarL(:)
        ELSEIF(sStar<=zero .and. sR>=zero)THEN
            xFlux(:,i)=fluxStarR(:)
        ELSEIF(sR<=zero)THEN
            xFlux(:,i)=fluxR(:)
        ENDIF

    END DO

    DO i=1,id-1
        rhs(:,i)=rhs(:,i)-(xFlux(:,i+1)-xFlux(:,i))
    END DO

    end subroutine HLLC_1D