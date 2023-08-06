    SUBROUTINE MUSCL_L_1D(i,variableL,variable)

    USE BasicData ,ONLY: p2,id,ghostLayers

    IMPLICIT NONE

    INTEGER i

    REAL(p2):: delta_plus,delta_minus,delta
    REAL(p2) :: variableL,fai
    REAL(p2) :: variable(1-ghostLayers:id+ghostLayers)

    delta_plus=variable(i)-variable(i-1)
    delta_minus=variable(i-1)-variable(i-2)
    CALL Limiter_function(delta_plus,delta_minus,delta,fai)
    variableL=variable(i-1)+delta
    
    END SUBROUTINE MUSCL_L_1D