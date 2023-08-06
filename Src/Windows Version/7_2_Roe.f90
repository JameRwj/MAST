    SUBROUTINE Roe(variableL,variableR,nx,ny,flux)
    USE BasicData,ONLY:p2,zero,half,one,two,gamma

    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: variableL,variableR,flux
    REAL(p2),DIMENSION(4) :: R1,R2,R3,R4,qL,qR,FluxL,FluxR,diss
    REAL(p2) :: rhoL,uL,vL,pL,cL,EL,hL,qnL,tnL
    REAL(p2) :: rhoR,uR,vR,pR,cR,ER,hR,qnR,tnR
    REAL(p2) :: rhoFace,cFace,uFace,vFace,hFace,qnFace,tnFace
    REAL(p2) :: lamda1,lamda2,lamda3,lamda4,ee,ee2
    REAL(p2) :: alfa1,alfa2,alfa3,alfa4

    REAL(p2):: nx,ny,tx,ty

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

    eL=pL/(gamma-one)+half*rhoL*(uL**2+vL**2)
    eR=pR/(gamma-one)+half*rhoR*(uR**2+vR**2)

    hL=(EL+pL)/rhoL
    hR=(ER+pR)/rhoR

    rhoFace=SQRT(rhoL*rhoR)
    uFace=(SQRT(rhoL)*uL+SQRT(rhoR)*uR)/(SQRT(rhoL)+SQRT(rhoR))
    vFace=(SQRT(rhoL)*vL+SQRT(rhoR)*vR)/(SQRT(rhoL)+SQRT(rhoR))
    hFace=(SQRT(rhoL)*hL+SQRT(rhoR)*hR)/(SQRT(rhoL)+SQRT(rhoR))
    cFace=SQRT( (gamma-one)*(hFace-half*(uFace**2+vFace**2)) )

    tx=-ny
    ty=nx

    qnL=uL*nx+vL*ny
    qnR=uR*nx+vR*ny
    qnFace=(SQRT(rhoL)*qnL+SQRT(rhoR)*qnR)/(SQRT(rhoL)+SQRT(rhoR))

    tnL=uL*tx+vL*ty
    tnR=uR*tx+vR*ty
    tnFace=(SQRT(rhoL)*tnL+SQRT(rhoR)*tnR)/(SQRT(rhoL)+SQRT(rhoR))

    lamda1=abs(qnFace-cFace)
    lamda2=abs(qnFace)
    lamda3=abs(qnFace)
    lamda4=abs(qnFace+cFace)

    !ee=0.2
    !if( lamda1<ee )  lamda1=half*(lamda1*lamda1/ee+ee)
    !if( lamda4<ee )  lamda4=half*(lamda4*lamda4/ee+ee)

    !ee2=0.9*cFace
    !if( lamda2<ee2 )  lamda2=half*(lamda2*lamda2/ee2+ee2)

    !ee2=2.0*cFace
    !if( lamda3<ee2 )  lamda3=half*(lamda3*lamda3/ee2+ee2)

    alfa1=(pR-pL-rhoFace*cFace*(qnR-qnL))/(2.0*cFace*cFace)
    alfa2=rhoR-rhoL-(pR-pL)/cFace/cFace
    alfa3=rhoFace*(tnR-tnL)/cFace
    alfa4=(pR-pL+rhoFace*cFace*(qnR-qnL))/(2.0*cFace*cFace)

    R1(1)=one
    R1(2)=uFace-cFace*nx
    R1(3)=vFace-cFace*ny
    R1(4)=hFace-qnFace*cFace

    R2(1)=one
    R2(2)=uFace
    R2(3)=vFace
    R2(4)=half*(uFace**2+vFace**2)

    R3(1)=zero
    R3(2)=cFace*tx
    R3(3)=cFace*ty
    R3(4)=cFace*tnFace

    R4(1)=one
    R4(2)=uFace+cFace*nx
    R4(3)=vFace+cFace*ny
    R4(4)=hFace+qnFace*cFace

    diss(:)=lamda1*alfa1*R1(:)+lamda2*alfa2*R2(:)+lamda3*alfa3*R3(:)+lamda4*alfa4*R4(:)

    FluxL(1)=rhoL*qnL
    FluxL(2)=rhoL*qnL*uL+nx*pL
    FluxL(3)=rhoL*qnL*vL+ny*pL
    FluxL(4)=rhoL*qnL*hL

    FluxR(1)=rhoR*qnR
    FluxR(2)=rhoR*qnR*uR+nx*pR
    FluxR(3)=rhoR*qnR*vR+ny*pR
    FluxR(4)=rhoR*qnR*hR

    flux(:)=half*(FluxL(:)+FluxR(:)-diss(:))
    END SUBROUTINE Roe
