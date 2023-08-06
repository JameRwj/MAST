
SUBROUTINE IterationStep
	
	USE BasicData,ONLY:Riemann_Solver,zero
    USE BasicData_1D,ONLY:rhs
	
	IMPLICIT NONE
	
	rhs=zero
	
	SELECT CASE(Riemann_solver)
    CASE(1)
        CALL Roe_1D()
    CASE(2)
        CALL HLLC_1D()
    CASE(3)
        CALL HLL_1D()
    CASE(4)
        CALL van_Leer_1D()
    CASE(5)
        CALL AUSMplus_1D()
    CASE(6)
        CALL SLAU_1D()
    CASE(7)
        CALL HLLE_1D()
    CASE(8)
        CALL HLLEM_1D() 
    END SELECT
	
END SUBROUTINE IterationStep