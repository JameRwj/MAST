SUBROUTINE CalTimeStep()
	
	USE BasicData,ONLY:p2,id,gamma
    USE BasicData_1D,ONLY:p_1D,rho_1D,u_1D,x_1D,deltaT,dt,cfl_1D,deltaX
	
	IMPLICIT NONE
	
	INTEGER ::i
    REAL(p2)::maxspeed

    DO i=1,id-1
        maxspeed = ABS(u_1D(i))+SQRT(gamma*p_1D(i)/rho_1D(i))
        deltaT(i)= cfl_1D*deltaX(i)/maxspeed
    END DO

    dt = deltaT(1)
    DO i = 1,id-1
        IF (deltaT(i) <= dt) THEN
            dt=deltaT(i)
        END IF
    END DO

       	
END SUBROUTINE CalTimeStep