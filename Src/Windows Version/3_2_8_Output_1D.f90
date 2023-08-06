    SUBROUTINE output_1D()

    USE BasicData,ONLY:gamma,id,one,half
    USE BasicData_1D,ONLY:x_1D,rho_1D,u_1D,p_1D

    IMPLICIT NONE

    INTEGER::i

    OPEN(101,FILE="Results\Result_1D.plt")
    WRITE(101,"(1X,A)") "VARIABLES= ""X"", ""Mass Density"", ""Velocity"", ""Pressure"",""Total Energy"",""Mach Number"",""Mass flux"""
    WRITE(101,"(1X,A)") "ZONE DATAPACKING=POINT"
    DO i=1,id-1
        WRITE(101,"(7E20.12)") half*(x_1D(i)+x_1D(i+1)),&
            rho_1D(i),&
            u_1D(i),  &
            p_1D(i),  &
            p_1D(i)/(gamma-one)+half*rho_1D(i)*u_1D(i)*u_1D(i),&
            u_1D(i)/SQRT(gamma*p_1D(i)/rho_1D(i)),&
            rho_1D(i)*u_1D(i)
    END DO
    CLOSE(101)

    END SUBROUTINE output_1D