    SUBROUTINE Limiter_function(delta_plus,delta_minus,delta,fai)

    USE BasicData,ONLY:p2,half,zero,one,Limiter,two

    IMPLICIT NONE

    REAL(p2):: delta_plus,delta_minus,delta
    REAL(p2):: r,psi,fai
    
    r=delta_plus/delta_minus
    
    IF(ABS(r) .GT. 1.0E+16)THEN
        r=SIGN(one,r)*1.0E+16
    END IF
    IF(delta_minus .EQ. zero)THEN
        r=SIGN(one,delta_plus)*1.0E+16
    END IF
    IF((delta_plus .EQ. zero))THEN
        r=zero
    END IF
    
    SELECT CASE(Limiter)
    CASE(0)		!no limiter
        delta = zero
        psi   = zero
    CASE(1)		!superbee
        psi   = MAX(MIN(two*r,one),MIN(r,two))
        psi   = MAX(zero,psi)
        delta = half*psi*delta_minus
    CASE(2)		!van Leer
        psi   = (r+ABS(r))/(one+ABS(r))
        delta = half*psi*delta_minus
    CASE(3)		!van Albada
        psi   = (r*r+r)/(one+r*r)
        delta = half*psi*delta_minus
    CASE(4)		!minmod
        psi   = MAX(MIN(r,one),zero)
        delta = half*psi*delta_minus
    CASE(5)     !From the paper of Deng Xi in JCP
        IF (r >= zero) THEN
            psi   = (two*r+two*two*r**2)/(one+two*r+3.0_p2*r**2)
        ELSE
            psi = zero
        END IF
        delta = half*psi*delta_minus
    END SELECT
    fai = half*psi

    END SUBROUTINE Limiter_function