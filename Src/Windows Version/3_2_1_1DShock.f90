    SUBROUTINE OneDShock()

    USE BasicData,ONLY:p2,Limiter,Riemann_Solver,Initialization_Method,id,x,ghostLayers,zero,epsilon
    USE BasicData_1D,ONLY:x_1D,rho_1D,u_1D,p_1D,deltaT,factor,qq,dqq,rhs,qq_0,residual,resid,ncount,maxCount,DeltaX,nst,dt,stageCoefficient

    IMPLICIT NONE

    INTEGER :: i,n

    WRITE(*,"(A80)")"========================Performing the 1D computation...========================"
    
    ALLOCATE(x_1D(id),rho_1D(1-ghostLayers:id-1+ghostLayers),u_1D(1-ghostLayers:id-1+ghostLayers),p_1D(1-ghostLayers:id-1+ghostLayers),deltaT(id-1),DeltaX(id-1))
    ALLOCATE(qq(3,1:id-1),dqq(3,1:id-1),rhs(3,1:id-1),qq_0(3,1:id-1))
    ALLOCATE(residual(1,3),resid(3))

    stageCoefficient(1)=0.1481_p2
    stageCoefficient(2)=0.4000_p2
    stageCoefficient(3)=1.0000_p2

    DO i = 1,id
        x_1D(i) = x(i,1)
    END DO

    DO i = 1,id-1
        DeltaX(i) = x_1D(i+1)-x_1D(i)
    END DO

    CALL Initialization_1D()
    
    IF (Riemann_Solver == 1 .OR. Riemann_Solver == 7 .OR. Riemann_Solver == 8) THEN
        IF (epsilon < 0.4) THEN
            factor = 1
        ELSE
            factor = 0
        END IF
    ELSE
        factor = 0
    END IF

    CALL Boundary()

    ncount=0
    DO WHILE (ncount <= maxCount)
        ncount = ncount+1

        CALL CalTimeStep()
        CALL Multistage_Runge_Kutta()

        resid=zero
        DO i=1,id-1
            DO n=1,3
                resid(n)=resid(n)+(dqq(n,i)/dt)**2
            END DO
        END DO

        IF (MOD(ncount-1,nst)==0)THEN
            WRITE(*,202) ncount-1
202         FORMAT(10("-"), " Step ", I7, 1X, " >>", 1X,10("-"))

            WRITE(*,203) LOG10(SQRT(SUM(resid(:))))
203         FORMAT(4X,"Residual=",F10.4)
        END IF

        IF(ISNAN(SUM(resid)))THEN
            PAUSE
            STOP
        END IF

        !IF(ncount .EQ. ((ncount-1)/nresidual)*nresidual+nst)then
        IF (MOD(ncount-1,nst)==0)THEN
            IF(ncount == 1)then
                OPEN(112,FILE="Results\Residual_1D.plt",FORM="FORMATTED")
            ELSE
                OPEN(112,FILE="Results\Residual_1D.plt",FORM="FORMATTED",POSITION="APPEND")
            ENDIF

            WRITE(112,"(1X,I8,4(E15.6))") ncount,log10(SQRT(resid(1)*resid(1)))
        END IF
        CLOSE(112)
    END DO

    CALL Output_1D()

    DEALLOCATE(x_1D,deltaT,DeltaX)
    DEALLOCATE(qq,dqq,rhs,qq_0)
    DEALLOCATE(residual,resid)
  
    WRITE(*,"(A80)")"=========================The 1D computation is finished========================="
    WRITE(*,*)

    END SUBROUTINE OneDShock