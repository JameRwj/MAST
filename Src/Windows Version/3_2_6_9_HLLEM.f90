    SUBROUTINE HLLEM_1D()

    USE BasicData,ONLY:gamma,id,one,two,half,p2,zero
    USE BasicData_1D,ONLY:rho_1D,u_1D,p_1D,rhs

    IMPLICIT NONE

    REAL(p2)::delta_plus,delta_minus,delta

    INTEGER i
    real(p2)::rhoL,uL,pL,hL,cL,EL
    real(p2)::rhoR,uR,pR,hR,cR,ER
    real(p2)::cFace,uFace,hFace
    real(p2)::sL,sR
    real(p2)::del,alfa2
    real(p2)::fp

    REAL(p2),DIMENSION(:,:),ALLOCATABLE :: xFlux

    REAL(p2),DIMENSION(:),ALLOCATABLE :: qL,qR,FL,FR,R2

    ALLOCATE(xFlux(3,id))
    ALLOCATE(R2(3),qL(3),qR(3),FL(3),FR(3))


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
        !Step 2:	Calculate the numerical flux
        EL=pL/(gamma-one)+half*rhoL*uL*uL
        ER=pR/(gamma-one)+half*rhoR*uR*uR

        hL=(gamma/(gamma-one))*pL/rhoL+(uL*uL)/two
        hR=(gamma/(gamma-one))*pR/rhoR+(uR*uR)/two

        cL=SQRT(gamma*pL/rhoL)
        cR=SQRT(gamma*pR/rhoR)

        hFace=(SQRT(rhoL)*hL+SQRT(rhoR)*hR)/(SQRT(rhoL)+SQRT(rhoR))
        uFace=(SQRT(rhoL)*uL+SQRT(rhoR)*uR)/(SQRT(rhoL)+SQRT(rhoR))
        cFace=SQRT((gamma-one)*(hFace-half*uFace*uFace))

        qL(1)=rhoL
        qL(2)=rhoL*uL
        qL(3)=EL

        qR(1)=rhoR
        qR(2)=rhoR*uR
        qR(3)=ER

        R2(1)=one
        R2(2)=uFace
        R2(3)=half*uFace*uFace

        FL(1)=rhoL*uL
        FL(2)=rhoL*uL*uL+pL
        FL(3)=uL*(EL+pL)

        FR(1)=rhoR*uR
        FR(2)=rhoR*uR*uR+pR
        FR(3)=uR*(ER+pR)

        !        sR=max(uL+cL,uR+cR)
        !        sL=min(uR-cR,uL-cL)
        !
        !sR=max(uFace+cFace,uR+cR)
        !sL=min(uFace-cFace,uL-cL)

        sR=uFace+cFace
        sL=uFace-cFace

        del=cFace/(ABS(uFace)+cFace)

        alfa2=rhoR-rhoL-(pR-pL)/cFace/cFace




        if(sL>=zero)then
            xFlux(:,i)=FL
        elseif(sL<zero .and. sR>zero)then
            xFlux(:,i)=sR/(sR-sL)*FL-sL/(sR-sL)*FR+sR*sL/(sR-sL)*(qR-qL-del*alfa2*R2)
        elseif(sR<=zero)then
            xFlux(:,i)=FR
        endif


        !        xFlux(:,i)=sR/(sR-sL)*FL-sL/(sR-sL)*FR+sR*sL/(sR-sL)*(qR-qL-del*alfa2*R2)

    END DO

    DO i=1,id-1
        rhs(:,i)=rhs(:,i)-(xFlux(:,i+1)-xFlux(:,i))
    END DO

    DEALLOCATE(xFlux)

    END SUBROUTINE HLLEM_1D