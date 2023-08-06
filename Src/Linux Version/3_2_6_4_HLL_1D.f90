    SUBROUTINE HLL_1D()

    USE BasicData,ONLY:gamma,id,one,two,half,p2,zero
    USE BasicData_1D,ONLY:rho_1D,u_1D,p_1D,rhs

    IMPLICIT NONE

    REAL(p2)::delta_plus,delta_minus,delta

    INTEGER i
    real(p2)::rhoL,uL,pL,hL,cL,EL
    real(p2)::rhoR,uR,pR,hR,cR,ER
    real(p2)::cFace,uFace,hFace
    real(p2)::sL,sR

    REAL(p2),DIMENSION(:,:),ALLOCATABLE :: xFlux

    REAL(p2),DIMENSION(:),ALLOCATABLE :: qL,qR,FL,FR

    ALLOCATE(xFlux(3,id))
    ALLOCATE(qL(3),qR(3),FL(3),FR(3))


    DO i=1,id
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

        hL=(gamma/(gamma-one))*pL/rhoL+(uL*uL)/two
        hR=(gamma/(gamma-one))*pR/rhoR+(uR*uR)/two

        cL=SQRT(gamma*pL/rhoL)
        cR=SQRT(gamma*pR/rhoR)

        qL(1)=rhoL
        qL(2)=rhoL*uL
        qL(3)=EL

        qR(1)=rhoR
        qR(2)=rhoR*uR
        qR(3)=ER

        FL(1)=rhoL*uL
        FL(2)=rhoL*uL*uL+pL
        FL(3)=uL*(EL+pL)

        FR(1)=rhoR*uR
        FR(2)=rhoR*uR*uR+pR
        FR(3)=uR*(ER+pR)

        sR=max(zero,uL+cL,uR+cR)
        sL=min(zero,uR-cR,uL-cL)

        xFlux(:,i)=sR/(sR-sL)*FL-sL/(sR-sL)*FR+sR*sL/(sR-sL)*(qR-qL)

    END DO

    DO i=1,id-1
        rhs(:,i)=rhs(:,i)-(xFlux(:,i+1)-xFlux(:,i))
    END DO

    DEALLOCATE(xFlux,qL,qR,FL,FR)

    END SUBROUTINE HLL_1D