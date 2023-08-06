    SUBROUTINE ROUND_L_1D(i,variableL,variable)

    USE BasicData ,ONLY: p2,id,ghostLayers,half,zero,one,two

    IMPLICIT NONE

    INTEGER i

    REAL(p2) :: variableL,variableL1
    REAL(p2) :: delta
    REAL(p2) :: faiL_bar,variableL_bar
    REAL(p2) :: omega0,omega1
    REAL(p2) :: gamma0,gamma1,lambda1
    REAL(p2) :: Temp1,Temp2,Temp3,Temp4
    REAL(p2) :: alpha1,alpha2,alpha3
    
    REAL(p2) :: variable(1-ghostLayers:id+ghostLayers)

    delta = 1.0E-16
    
    gamma0  = 1100.0_p2
    gamma1  = 800.0_p2
    lambda1 = 0.15_p2
    
    !faiL_bar = (variable(i-1)-variable(i-2)+delta)/(variable(i)-variable(i-2)+delta)
    
    faiL_bar = (variable(i-1)-variable(i-2))/(variable(i)-variable(i-2))
    
    omega0 = one/(one+gamma0*(faiL_bar-one)**4)**2
    omega1 = one/(one+gamma1*(faiL_bar-one)**4)**2
    
    Temp1 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiL_bar)*omega0 + two*faiL_bar*(one-omega0)
    Temp2 = two*faiL_bar
    Temp3 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiL_bar)*omega1 + (lambda1*faiL_bar-lambda1+one)*(one-omega1)
    Temp4 = lambda1*faiL_bar-lambda1+one
    
    IF ((faiL_bar > zero) .AND. (faiL_bar <= half)) THEN
        variableL_bar = MIN(Temp1,Temp2)
    ELSEIF ((faiL_bar > half) .AND. (faiL_bar <= one)) THEN
        variableL_bar = MIN(Temp3,Temp4)
    ELSE
        variableL_bar = faiL_bar
    END IF
    
    variableL = variableL_bar*(variable(i)-variable(i-2))+variable(i-2)
    
    IF (variable(i)-variable(i-2) == 0) THEN
        variableL = variable(i-1)
    END IF
        
    !IF ((faiL_bar > zero) .AND. (faiL_bar <= half)) THEN
    !    IF (Temp1 <= Temp2) THEN
    !        alpha1 = one/3.0_p2*omega0
    !        alpha2 = two-7.0_p2/6.0_p2*omega0
    !        alpha3 = 5.0_p2/6.0_p2*omega0-one
    !    ELSE
    !        alpha1 = zero
    !        alpha2 = two
    !        alpha3 = -one
    !    END IF
    !ELSEIF ((faiL_bar > half) .AND. (faiL_bar <= one)) THEN
    !    IF (Temp3 <= Temp4) THEN
    !        alpha1 = one/3.0_p2*omega1 + (one-omega1)*(one-lambda1)
    !        alpha2 = 5.0_p2/6.0_p2*omega1 + lambda1*(one-omega1)
    !        alpha3 = -one/6.0_p2*omega1
    !    ELSE
    !        alpha1 = one-lambda1
    !        alpha2 = lambda1
    !        alpha3 = zero
    !    END IF
    !ELSE
    !    alpha1 = zero
    !    alpha2 = one
    !    alpha3 = zero
    !END IF
    !
    !variableL = alpha1*variable(i) + alpha2*variable(i-1) + alpha3*variable(i-2)
    
    END SUBROUTINE ROUND_L_1D