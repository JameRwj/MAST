    SUBROUTINE Information_Write()
    USE BasicData,ONLY : Limiter,Riemann_Solver,Mainf,epsilon,id,jd,Initialization_Method,TestCase,Reconstruction_Method

    IMPLICIT NONE

    CHARACTER*50 :: CharTemp1,CharTemp2
    CHARACTER*100:: CharGrid,CharLimiter,CharSolver,CharInitialization,CharTestCase,CharReconstruction
    INTEGER :: key,error

    WRITE(CharTemp1,*)id-1
    WRITE(CharTemp2,*)jd-1
    CharGrid = TRIM(ADJUSTL(TRIM(ADJUSTL(CharTemp1))//"*"//TRIM(ADJUSTL(CharTemp2))))

    error = 0
    
    SELECT CASE (Reconstruction_Method)   
    CASE(1)
        CharReconstruction = "MUSCL"
    CASE(2)
        CharReconstruction = "ROUND"
    CASE DEFAULT
        CharReconstruction = "Warning! The reconstruction method is wrong"
        error = 1
    END SELECT
        
    SELECT CASE (Limiter)
    CASE (0)
        CharLimiter = "No Limiter (First Order)"
    CASE (1)
        CharLimiter = "Superbee Limiter"
    CASE (2)
        CharLimiter = "van Leer Limiter"
    CASE (3)
        CharLimiter = "van Albada Limiter"
    CASE (4)
        CharLimiter = "minmod Limiter"
    CASE (5)
        CharLimiter = "Limiter of Deng Xi"
    CASE DEFAULT
        CharLimiter = "Warning! The limiter is wrong"
        error = 1
    END SELECT

    SELECT CASE (Riemann_Solver)
    CASE (1)
        CharSolver = "Roe Solver"
    CASE (2)
        CharSolver = "HLLC Solver"
    CASE (3)
        CharSolver = "HLL Solver"
    CASE (4)
        CharSolver = "van Leer Solver"
    CASE (5)
        CharSolver = "AUSMplus Solver"
    CASE (6)
        CharSolver = "SLAU Solver"
    CASE (7)
        CharSolver = "HLLE Solver"
    CASE (8)
        CharSolver = "HLLEM Solver"
    CASE DEFAULT
        CharSolver = "Warning! The solver is wrong"
        error = 1
    END SELECT

    SELECT CASE (TestCase)
    CASE(1)
        CharTestCase = "2D steady normal shock"
        SELECT CASE (Initialization_Method)
        CASE(1)
            CharInitialization = "the Rankine-Hugoniot conditions"
        CASE(2)
            CharInitialization = "the 1D computation"
        CASE DEFAULT
            CharInitialization = "Warning! The initialization method is wrong"
            error = 1
        END SELECT
    CASE(2)
        CharTestCase = "From the file in the InitialFlow floder"
    CASE DEFAULT
        CharTestCase = "Warning! The test case is wrong"
        error = 1
    END SELECT

    WRITE(*,*)
    WRITE(*,"(A80)")"===============================Computation Settings==============================="
    WRITE(*,"(A35,A60)")"The computational grid is:         " , CharGrid
    WRITE(*,*)
    WRITE(*,"(A35,A60)")"The reconstruction method is:      " , CharReconstruction
    WRITE(*,*)
    IF (Reconstruction_Method == 1) THEN
        WRITE(*,"(A35,A60)")"The used limiter is:               " , CharLimiter
        WRITE(*,*)
    END IF
    WRITE(*,"(A35,A60)")"The used Riemann solver is:        " , CharSolver
    WRITE(*,*)
    WRITE(*,"(A35,A60)")"The test case is:                  " , CharTestCase
    IF (TestCase == 1) THEN
        WRITE(*,*)
        WRITE(*,"(A35,F5.2)")"The Mach number is:                " , Mainf
        WRITE(*,*)
        WRITE(*,"(A35,F4.2)")"The numerical shock position is:   " , epsilon
        Write(*,*)
        WRITE(*,"(A35,A60)")"The flow field is initialized by:   " , CharInitialization
        WRITE(*,*)
    END IF
    
    IF (error==1) THEN
        WRITE(*,*)"There is at least one error in Setting.dat. Please check."
        PAUSE
        STOP
    END IF
    
    WRITE(*,"(A80)")"=====================Start the calculation?(Yes-1,No-else)======================"
    READ(*,*)key
    
    IF (key == 1) THEN
        WRITE(*,*)"The computation is ongoing, please wait......"
    ELSE 
        STOP
    END IF
        
    END SUBROUTINE Information_Write