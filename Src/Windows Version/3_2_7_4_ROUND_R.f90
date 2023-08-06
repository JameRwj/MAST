    SUBROUTINE ROUND_R_1D(i,variableR,variable)

    USE BasicData ,ONLY: p2,id,ghostLayers,half,zero,one,two

    IMPLICIT NONE

    INTEGER i

    REAL(p2) :: variableR
    REAL(p2) :: delta
    REAL(p2) :: faiR_bar,variableR_bar
    REAL(p2) :: omega0,omega1
    REAL(p2) :: gamma0,gamma1,lambda1
    REAL(p2) :: Temp1,Temp2,Temp3,Temp4
    REAL(p2) :: variable(1-ghostLayers:id+ghostLayers)
    REAL(p2) :: alpha1,alpha2,alpha3

    delta = 1.0E-16
    
    gamma0  = 1100.0_p2
    gamma1  = 800.0_p2
    lambda1 = 0.15_p2
    
    !faiR_bar = (variable(i)-variable(i+1)+delta)/(variable(i-1)-variable(i+1)+delta)
    
    faiR_bar = (variable(i)-variable(i+1))/(variable(i-1)-variable(i+1))
    
    omega0 = one/(one+gamma0*(faiR_bar-one)**4)**2
    omega1 = one/(one+gamma1*(faiR_bar-one)**4)**2
    
    Temp1 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiR_bar)*omega0 + two*faiR_bar*(one-omega0)
    Temp2 = two*faiR_bar
    Temp3 = (one/3.0_p2 + 5.0_p2/6.0_p2*faiR_bar)*omega1 + (lambda1*faiR_bar-lambda1+one)*(one-omega1)
    Temp4 = lambda1*faiR_bar-lambda1+one
    
    IF ((faiR_bar > 0.0) .AND. (faiR_bar <= half)) THEN
        variableR_bar = MIN(Temp1,Temp2)
    ELSEIF ((faiR_bar > half) .AND. (faiR_bar <= one)) THEN
        variableR_bar = MIN(Temp3,Temp4)
    ELSE
        variableR_bar = faiR_bar
    END IF
    
    variableR = variableR_bar*(variable(i-1)-variable(i+1))+variable(i+1)
    
    IF (variable(i-1)-variable(i+1)==0)THEN
        variableR = variable(i)
    END IF
    
    !
    !IF ((faiR_bar > zero) .AND. (faiR_bar <= half)) THEN
    !    IF (Temp1 <= Temp2) THEN
    !        alpha1 = one/3.0_p2*omega0
    !        alpha2 = two-7.0_p2/6.0_p2*omega0
    !        alpha3 = 5.0_p2/6.0_p2*omega0-one
    !    ELSE
    !        alpha1 = zero
    !        alpha2 = two
    !        alpha3 = -one
    !    END IF
    !ELSEIF ((faiR_bar > half) .AND. (faiR_bar <= one)) THEN
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
    !variableR = alpha1*variable(i-1)+alpha2*variable(i)+alpha3*variable(i+1)
    END SUBROUTINE ROUND_R_1D