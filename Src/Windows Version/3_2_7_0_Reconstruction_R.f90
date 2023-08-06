    SUBROUTINE Reconstruction_R(i,variableR,variable)

    USE BasicData ,ONLY: Reconstruction_Method,p2,ghostLayers,id

    IMPLICIT NONE

    INTEGER i

    REAL(p2) :: variableR,variable(1-ghostLayers:id+ghostLayers)

    SELECT CASE(Reconstruction_Method)
    CASE(1)
        CALL MUSCL_R_1D(i,variableR,variable)
    CASE(2)
        CALL ROUND_R_1D(i,variableR,variable)
    END SELECT

    END SUBROUTINE Reconstruction_R