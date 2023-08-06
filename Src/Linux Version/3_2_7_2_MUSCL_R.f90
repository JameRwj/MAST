    SUBROUTINE MUSCL_R_1D(i,variableR,variable)

    USE BasicData ,ONLY: p2,id,ghostLayers

    IMPLICIT NONE

    INTEGER i

    REAL(p2):: delta_plus,delta_minus,delta
    REAL(p2) :: variableR,fai
    REAL(p2) :: variable(1-ghostLayers:id+ghostLayers)

    delta_plus=variable(i+1)-variable(i)
    delta_minus=variable(i)-variable(i-1)
    CALL Limiter_function(delta_minus,delta_plus,delta,fai)
    variableR=variable(i)-delta

    END SUBROUTINE MUSCL_R_1D
