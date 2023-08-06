    SUBROUTINE ReadSettings()

    USE BasicData,ONLY:Limiter,Riemann_Solver,epsilon,Mainf,Initialization_Method,TestCase,Reconstruction_Method
    USE BasicData_1D,ONLY:maxCount

    IMPLICIT NONE

    OPEN(103,FILE="Settings.dat")
    READ(103,*) Reconstruction_Method
    READ(103,*) Limiter
    READ(103,*) Riemann_Solver
    READ(103,*) TestCase
    READ(103,*) epsilon
    READ(103,*) Mainf
    READ(103,*) Initialization_Method
    READ(103,*) maxCount
    CLOSE(103)

    END SUBROUTINE ReadSettings