    SUBROUTINE AUSMplus_1D()

    USE BasicData,ONLY:gamma,id,one,two,half,p2,zero
    USE BasicData_1D,ONLY:rho_1D,u_1D,p_1D,rhs

    IMPLICIT NONE

    REAL(p2):: delta_plus,delta_minus,delta

    INTEGER i
    REAL(p2):: rhoL,uL,pL,hL,cL,maL,maR
    REAL(p2):: rhoR,uR,pR,hR,cR
    REAL(p2):: cCriticalL,cCriticalR,cFace,maFace
    REAL(p2):: maPlus,maMinus
    REAL(p2):: pPlus,pMinus

    REAL(p2):: sL,sR,del,fp

    REAL(p2),DIMENSION(:,:),ALLOCATABLE :: xFlux

    ALLOCATE(xFlux(3,id))

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

        hL=(gamma/(gamma-one))*pL/rhoL+uL*uL/two
        hR=(gamma/(gamma-one))*pR/rhoR+uR*uR/two

        !------------------------------------------------------------------------------------------------------------
        cCriticalL=SQRT(two*(gamma-one)/(gamma+one)*hL)
        cCriticalR=SQRT(two*(gamma-one)/(gamma+one)*hR)

        cL=(cCriticalL*cCriticalL)/MAX(cCriticalL,ABS(uL))
        cR=(cCriticalR*cCriticalR)/MAX(cCriticalR,ABS(uR))

        cFace=MIN(cL,cR)

        maL=uL/cFace
        maR=uR/cFace

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
            xFlux(1,i)=maFace*rhoL*cFace
            xFlux(2,i)=maFace*rhoL*uL*cFace+(pPlus*pL+pMinus*pR)
            xFlux(3,i)=maFace*rhoL*hL*cFace
        ELSE
            xFlux(1,i)=maFace*rhoR*cFace
            xFlux(2,i)=maFace*rhoR*uR*cFace+(pPlus*pL+pMinus*pR)
            xFlux(3,i)=maFace*rhoR*hR*cFace
        END IF

    END DO

    DO i=1,id-1
        rhs(:,i)=rhs(:,i)-(xFlux(:,i+1)-xFlux(:,i))
    END DO

    DEALLOCATE(xFlux)

    END SUBROUTINE AUSMplus_1D
