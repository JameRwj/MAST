    SUBROUTINE CalConserv()

    USE BasicData,ONLY:id,p2,gamma,one,half
    USE BasicData_1D,ONLY:rho_1D,u_1D,p_1D,qq,factor

    IMPLICIT NONE

    INTEGER i,ii

    DO i=1,id-1
        qq(1,i)=rho_1D(i)
        qq(2,i)=rho_1D(i)*u_1D(i)
        qq(3,i)=p_1D(i)/(gamma-one)+half*rho_1D(i)*u_1D(i)*u_1D(i)
    END DO

    ii=FLOOR((id)*0.5)

    IF (factor==1) THEN
        qq(2,ii+1)=qq(2,ii-1)
    END IF
    !IF (factor == 2) THEN
    !    qq(2,ii+2)=qq(2,ii-1)
    !END IF
    

    END SUBROUTINE CalConserv