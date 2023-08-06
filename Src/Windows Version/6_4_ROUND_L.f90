    SUBROUTINE ROUND_L(i,j,U,Face,alpha1,alpha2,alpha3,UL)
    USE BasicData,ONLY: p2,id,jd,ghostLayers,zero,one,two,half

    IMPLICIT NONE

    INTEGER i,j,Face

    REAL(p2) :: u1,u2,u3
    REAL(p2) :: alpha1,alpha2,alpha3
    REAL(p2) :: faiL,UL
    REAL(p2) :: faiL_bar,UL_bar
    REAL(p2) :: omega_0,omega_1
    REAL(p2) :: gamma_0,gamma_1,lambda_1
    REAL(p2) :: Temp1,Temp2,Temp3,Temp4
    REAL(p2) :: delta

    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)

    delta = 1.0E-16

    SELECT CASE(Face)
    CASE(1)
        u1 = U(i+1,j)
        u2 = U(i-0,j)
        u3 = U(i-1,j)
    CASE(2)
        u1 = U(i,j+1)
        u2 = U(i,j-0)
        u3 = U(i,j-1)
    CASE(3)
        u1 = U(i-0,j)
        u2 = U(i-1,j)
        u3 = U(i-2,j)
    CASE(4)
        u1 = U(i,j-0)
        u2 = U(i,j-1)
        u3 = U(i,j-2)
    END SELECT

    !faiL_bar = (u2-u3+delta)/(u1-u3+delta)
    
    faiL_bar = (u2-u3)/(u1-u3)

    !gamma_0 = 7.8_p2
    !gamma_1 = 12.0_p2
    !
    !omega_0 = one/(one+gamma_0*(faiL_bar-one)**2)**6
    !omega_1 = one/(one+gamma_1*(faiL_bar-one)**2)**4
    !
    !alpha1 = one/3.0_p2*(one-omega_0)*(1-omega_1) + 3.0_p2/two*(one-omega_0)*omega_1 + one/two*omega_0
    !alpha2 = 5.0_p2/6.0_p2*(one-omega_0)*(1-omega_1) + one/two*omega_0
    !alpha3 = -(one/3.0_p2*(one-omega_0)*(1-omega_1) + 3.0_p2/two*(one-omega_0)*omega_1 + one/two*omega_0&
    !         -5.0_p2/6.0_p2*(one-omega_0)*(1-omega_1) + one/two*omega_0 - one)
    
    gamma_0  = 1100.0_p2
    gamma_1  = 800.0_p2
    lambda_1 = 0.15_p2
    
    omega_0 = one/(one+gamma_0*(faiL_bar-one)**4)**2
    omega_1 = one/(one+gamma_1*(faiL_bar-one)**4)**2
    
    Temp1 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiL_bar)*omega_0 + two*faiL_bar*(one-omega_0)
    Temp2 = two*faiL_bar
    Temp3 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiL_bar)*omega_1 + (lambda_1*faiL_bar-lambda_1+one)*(one-omega_1)
    Temp4 = lambda_1*faiL_bar-lambda_1+one
    
    IF ((faiL_bar > zero) .AND. (faiL_bar <= half)) THEN
        IF (Temp1 <= Temp2) THEN
            alpha1 = one/3.0_p2*omega_0
            alpha2 = two-7.0_p2/6.0_p2*omega_0
            alpha3 = 5.0_p2/6.0_p2*omega_0-one
        ELSE
            alpha1 = zero
            alpha2 = two
            alpha3 = -one
        END IF
    ELSEIF ((faiL_bar > half) .AND. (faiL_bar <= one)) THEN
        IF (Temp3 <= Temp4) THEN
            alpha1 = one/3.0_p2*omega_1 + (one-omega_1)*(one-lambda_1)
            alpha2 = 5.0_p2/6.0_p2*omega_1 + lambda_1*(one-omega_1)
            alpha3 = -one/6.0_p2*omega_1
        ELSE
            alpha1 = one-lambda_1
            alpha2 = lambda_1
            alpha3 = zero
        END IF
    ELSE
        alpha1 = zero
        alpha2 = one
        alpha3 = zero
    END IF

    IF (u1-u3==0.0) THEN
        alpha1 = zero
        alpha2 = one
        alpha3 = zero
    END IF
    
    UL = alpha1*u1+alpha2*u2+alpha3*u3
    
    END SUBROUTINE ROUND_L