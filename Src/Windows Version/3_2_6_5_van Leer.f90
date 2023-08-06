    SUBROUTINE van_Leer_1D()

    USE BasicData,ONLY:gamma,id,one,two,half,p2,zero
    USE BasicData_1D,ONLY:rho_1D,u_1D,p_1D,rhs

    IMPLICIT NONE

    REAL(p2):: delta_plus,delta_minus,delta

    INTEGER i
    REAL(p2):: rhoL,uL,pL,hL,cL,maL,maPlus
    REAL(p2):: rhoR,uR,pR,hR,cR,maR,maMinus
    REAL(p2):: maN

    REAL(p2),DIMENSION(:),ALLOCATABLE :: fluxL,fluxR
    REAL(p2),DIMENSION(:,:),ALLOCATABLE :: xFlux

    ALLOCATE(fluxL(3),fluxR(3))
    ALLOCATE(xFlux(3,id))


    DO i=1,id
        !Step 1:	Define the left and right states of cell interface with MUSCL.
        !		!Left state
        CALL Reconstruction_L(i,rhoL,rho_1D)

        CALL Reconstruction_L(i,uL,u_1D)

        CALL Reconstruction_L(i,pL,p_1D)
        
        !		Right state
        CALL Reconstruction_R(i,rhoR,rho_1D)
        
        CALL Reconstruction_R(i,uR,u_1D)

        CALL Reconstruction_R(i,pR,p_1D)

        !Step 2:	Calculate the numerical flux
        cL=SQRT(gamma*pL/rhoL)
        maL=uL/cL
        IF(maL .GE. one)THEN
            maPlus=maL
        ELSE IF(maL .LE. -one)THEN
            maPlus=zero
        ELSE
            maPlus=0.25_p2*(maL+one)**2
        END IF

        cR=SQRT(gamma*pR/rhoR)
        maR=uR/cR
        IF(maR .GE. one)THEN
            maMinus=zero
        ELSE IF(maR .LE. -one)THEN
            maMinus=maR
        ELSE
            maMinus=-0.25_p2*(maR-one)**2
        END IF

        maN=maPlus+maMinus

        IF(maL .GE. one)THEN
            fluxL(1)=rhoL*uL
            fluxL(2)=rhoL*uL*uL+pL
            fluxL(3)=(pL*gamma/(gamma-one)+half*rhoL*(uL*uL))*uL
        ELSE IF(maL .LE. -one)THEN
            fluxL(1)=zero
            fluxL(2)=zero
            fluxL(3)=zero
        ELSE
            fluxL(1)=0.25_p2*rhoL*cL*(MaL+one)*(MaL+one)
            fluxL(2)=fluxL(1)*((-uL+two*cL)/gamma+uL)
            fluxL(3)=fluxL(1)*half*(((gamma-one)*uL+two*cL)**2/(gamma*gamma-one))
        END IF

        IF(maR .GE. one)THEN
            fluxR(1)=zero
            fluxR(2)=zero
            fluxR(3)=zero
        ELSE IF(maR .LE. -one)THEN
            fluxR(1)=rhoR*uR
            fluxR(2)=rhoR*uR*uR+pR
            fluxR(3)=(pR*gamma/(gamma-one)+half*rhoR*(uR*uR))*uR
        ELSE
            fluxR(1)=-0.25_p2*rhoR*cR*(MaR-one)*(MaR-one)
            fluxR(2)=fluxR(1)*((-uR-two*cR)/gamma+uR)
            fluxR(3)=fluxR(1)*half*(((gamma-one)*uR-two*cR)**2/(gamma*gamma-one))
        END IF

        xFlux(:,i)=fluxL(:)+fluxR(:)

    END DO

    DO i=1,id-1
        rhs(:,i)=rhs(:,i)-(xFlux(:,i+1)-xFlux(:,i))
    END DO

    DEALLOCATE(fluxL,fluxR)
    DEALLOCATE(xFlux)

    END SUBROUTINE van_Leer_1D