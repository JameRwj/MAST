    SUBROUTINE Initialization_1D()

    USE BasicData,ONLY:id,p2,epsilon,Mainf,one,two,half,gamma
    USE BasicData_1D,ONLY:rho_1D,u_1D,p_1D

    IMPLICIT NONE

    INTEGER i,ii


    REAL(p2)::rhoL,rhoR,uL,uR,pL,pR

    REAL(p2)::deltaU,deltaP

    !×ó±ßµÄ×´Ì¬³õÖµ
    rhoL = one
    uL   = one
    pL   = one/gamma/mainf/mainf

    !ÓÒ±ßµÄ×´Ì¬³õÖµ
    rhoR = one/( two/( (gamma+one)*mainf*mainf )+(gamma-one)/(gamma+one) )
    uR   = two/( (gamma+one)*mainf*mainf )+(gamma-one)/(gamma+one)
    pR   = ( two*gamma*mainf*mainf/(gamma+one)-(gamma-one)/(gamma+one) )/gamma/mainf/mainf

    ii = FLOOR(0.5*id)
    DO i = 1,ii-1
        rho_1D(i) = rhoL
        u_1D(i)   = uL
        p_1D(i)   = pL
    END DO

    DO i = ii+1,id-1
        rho_1D(i) = rhoR
        u_1D(i)   = uR
        p_1D(i)   = pR
    END DO

    deltaU = one-(one-epsilon)/SQRT(one+epsilon*(mainf*mainf-1.0)/(one+(gamma-one)*mainf*mainf/two) )/SQRT(one+epsilon*(mainf*mainf-one)/(one-two*gamma*mainf*mainf/(gamma-one)) )
    deltaP = epsilon/SQRT(one+(one-epsilon)*(gamma+one)/(gamma-one)*(mainf*mainf-one)/mainf/mainf)

    rho_1D(ii)= (one-epsilon)*rhoL+epsilon*rhoR
    u_1D(ii)  = (one-deltaU)*uL+deltaU*uR
    p_1D(ii)  = (one-deltaP)*pL+deltaP*pR

    END SUBROUTINE Initialization_1D