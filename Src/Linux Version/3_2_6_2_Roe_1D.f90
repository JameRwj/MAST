    SUBROUTINE Roe_1D()

    USE BasicData,ONLY:gamma,id,one,two,half,p2,zero
    USE BasicData_1D,ONLY:rho_1D,u_1D,p_1D,rhs

    IMPLICIT NONE

    REAL(p2):: delta_plus,delta_minus,delta

    INTEGER i
    REAL(p2):: rhoL,uL,pL,hL,cL,EL
    REAL(p2):: rhoR,uR,pR,hR,cR,ER
    REAL(p2):: rhoFace,cFace,uFace,hFace
    REAL(p2):: sL,sR
    REAL(p2):: alfa1,alfa2,alfa3
    REAL(p2):: fp
    REAL(p2):: lamda1,lamda2,lamda3,ee

    REAL(p2),DIMENSION(:,:),ALLOCATABLE :: xFlux
    REAL(p2),DIMENSION(:),ALLOCATABLE :: qL,qR,FL,FR,R1,R2,R3

    ALLOCATE(xFlux(3,id))
    ALLOCATE(R1(3),R2(3),R3(3),qL(3),qR(3),FL(3),FR(3))

    DO i=1,id

        !		Left state
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

        rhoFace=sqrt(rhoL*rhoR)
        hFace=(SQRT(rhoL)*hL+SQRT(rhoR)*hR)/(SQRT(rhoL)+SQRT(rhoR))
        uFace=(SQRT(rhoL)*uL+SQRT(rhoR)*uR)/(SQRT(rhoL)+SQRT(rhoR))
        cFace=SQRT((gamma-one)*(hFace-half*uFace**2))


        alfa1=(pR-pL-rhoFace*cFace*(uR-uL))/two/cFace/cFace
        alfa2= rhoR-rhoL-(pR-pL)/cFace/cFace
        alfa3=(pR-pL+rhoFace*cFace*(uR-uL))/two/cFace/cFace

        qL(1)=rhoL
        qL(2)=rhoL*uL
        qL(3)=EL

        qR(1)=rhoR
        qR(2)=rhoR*uR
        qR(3)=ER

        R1(1) = one
        R1(2) = uFace - cFace
        R1(3) = hFace - uFace*cFace

        R2(1)=  one
        R2(2)=  uFace
        R2(3)=  half*uFace*uFace

        R3(1) = one
        R3(2) = uFace + cFace
        R3(3) = hFace + uFace*cFace

        FL(1)=rhoL*uL
        FL(2)=rhoL*uL*uL+pL
        FL(3)=uL*(EL+pL)

        FR(1)=rhoR*uR
        FR(2)=rhoR*uR*uR+pR
        FR(3)=uR*(ER+pR)

        lamda1=abs(uFace-cFace)
        lamda2=abs(uFace)
        lamda3=abs(uFace+cFace)

        xFlux(:,i)=half*(FL(:)+FR(:)-lamda1*alfa1*R1(:)-lamda2*alfa2*R2(:)-lamda3*alfa3*R3(:))

    END DO

    DO i=1,id-1
        rhs(:,i)=rhs(:,i)-(xFlux(:,i+1)-xFlux(:,i))
    END DO

    DEALLOCATE(xFlux)

    END SUBROUTINE Roe_1D