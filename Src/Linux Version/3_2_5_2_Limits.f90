    SUBROUTINE Limits()

    USE BasicData,ONLY:p2,id,gamma,one,half
    USE BasicData_1D,ONLY:rho_1D,u_1D,p_1D,qq

    IMPLICIT NONE

    INTEGER i

    DO i=1,id-1
        rho_1D(i) = qq(1,i)
        u_1D(i)   = qq(2,i)/rho_1D(i)
        p_1D(i)   = (gamma-one)*(qq(3,i)-half*rho_1D(i)*u_1D(i)*u_1D(i))
    END DO

    DO i=1,id-1
        IF(ISNAN(rho_1D(i)))THEN
            PAUSE
        END IF
        IF(ISNAN(u_1D(i)))THEN
            PAUSE
        END IF
        IF(ISNAN(p_1D(i)))THEN
            PAUSE
        END IF
    END DO

    END SUBROUTINE Limits