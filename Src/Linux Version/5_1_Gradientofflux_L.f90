    SUBROUTINE Gradientofflux_L(i,j,Face,nx,ny,L,faiL_1,faiL_2,faiL_3)
    USE BasicData,ONLY:p2,rho,u,v,p,zero,two,Reconstruction_Method

    IMPLICIT NONE

    INTEGER i,j,Face,m,n
    REAL(p2),DIMENSION(4,4) :: L,faiL_1,faiL_2,faiL_3
    REAL(p2),DIMENSION(4) :: flux1_plus ,flux2_plus ,flux3_plus ,flux4_plus
    REAL(p2),DIMENSION(4) :: flux1_minus,flux2_minus,flux3_minus,flux4_minus
    REAL(p2),DIMENSION(4) :: UL,UR

    REAL(p2) :: faiL1_1,faiL2_1,faiL3_1,faiL4_1,&
                faiR1_1,faiR2_1,faiR3_1,faiR4_1
    REAL(p2) :: faiL1_2,faiL2_2,faiL3_2,faiL4_2,&
                faiR1_2,faiR2_2,faiR3_2,faiR4_2
    REAL(p2) :: faiL1_3,faiL2_3,faiL3_3,faiL4_3,&
                faiR1_3,faiR2_3,faiR3_3,faiR4_3
    
    REAL(p2) :: UL1,UL2,UL3,UL4,UR1,UR2,UR3,UR4
    REAL(p2) :: delta,deltaL,deltaR
    REAL(p2) :: nx,ny

    delta = 1.0E-7_p2
    
    faiL_1 = zero
    faiL_2 = zero
    faiL_3 = zero

    !===========重构===========
    CALL Reconstruct_L(i,j,rho,Face,faiL1_1,faiL1_2,faiL1_3,UL1)
    CALL Reconstruct_L(i,j,u,  Face,faiL2_1,faiL2_2,faiL2_3,UL2)
    CALL Reconstruct_L(i,j,v,  Face,faiL3_1,faiL3_2,faiL3_3,UL3)
    CALL Reconstruct_L(i,j,p,  Face,faiL4_1,faiL4_2,faiL4_3,UL4)
                                         
    CALL Reconstruct_R(i,j,rho,Face,faiR1_1,faiR1_2,faiR1_3,UR1)
    CALL Reconstruct_R(i,j,u,  Face,faiR2_1,faiR2_2,faiR2_3,UR2)
    CALL Reconstruct_R(i,j,v,  Face,faiR3_1,faiR3_2,faiR3_3,UR3)
    CALL Reconstruct_R(i,j,p,  Face,faiR4_1,faiR4_2,faiR4_3,UR4)

    deltaL = delta
    deltaR = zero

    !===========左侧加减扰动计算通量===========
    !第一项加扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(1) = UL(1)+deltaL
    UR(1) = UR(1)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux1_plus)

    !第一项减扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(1) = UL(1)-deltaL
    UR(1) = UR(1)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux1_minus)

    !第二项加扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(2) = UL(2)+deltaL
    UR(2) = UR(2)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux2_plus)

    !第二项减扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(2) = UL(2)-deltaL
    UR(2) = UR(2)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux2_minus)

    !第三项加扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(3) = UL(3)+deltaL
    UR(3) = UR(3)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux3_plus)

    !第三项减扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(3) = UL(3)-deltaL
    UR(3) = UR(3)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux3_minus)

    !第四项加扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(4) = UL(4)+deltaL
    UR(4) = UR(4)+deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux4_plus)

    !第四项减扰动
    UL = [UL1,UL2,UL3,UL4]
    UR = [UR1,UR2,UR3,UR4]

    UL(4) = UL(4)-deltaL
    UR(4) = UR(4)-deltaR

    CALL Calculateflux(UL,UR,nx,ny,flux4_minus)

    !=======求每一个通量梯度=======

    L(:,1) = (flux1_plus-flux1_minus)/two/delta
    L(:,2) = (flux2_plus-flux2_minus)/two/delta
    L(:,3) = (flux3_plus-flux3_minus)/two/delta
    L(:,4) = (flux4_plus-flux4_minus)/two/delta

    IF ((Face == 3) .OR. (Face == 4)) THEN
        DO m = 1,4
            DO n = 1,4
                L(m,n) = -L(m,n)
            END DO
        END DO
    END IF

    SELECT CASE(Reconstruction_Method)
    CASE (1)   
        faiL_1(1,1) = faiL1_1
        faiL_1(2,2) = faiL2_1
        faiL_1(3,3) = faiL3_1
        faiL_1(4,4) = faiL4_1
        
        faiL_2 = zero
        faiL_3 = zero
    CASE (2)
        faiL_1(1,1) = faiL1_1
        faiL_1(2,2) = faiL2_1
        faiL_1(3,3) = faiL3_1
        faiL_1(4,4) = faiL4_1
        
        faiL_2(1,1) = faiL1_2
        faiL_2(2,2) = faiL2_2
        faiL_2(3,3) = faiL3_2
        faiL_2(4,4) = faiL4_2
        
        faiL_3(1,1) = faiL1_3
        faiL_3(2,2) = faiL2_3
        faiL_3(3,3) = faiL3_3
        faiL_3(4,4) = faiL4_3
    END SELECT
    
    END SUBROUTINE Gradientofflux_L