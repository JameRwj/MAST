    SUBROUTINE Reconstruct_R(i,j,U,Face,faiR,UR)
    USE BasicData,ONLY: p2,id,jd,ghostLayers,zero

    IMPLICIT NONE

    INTEGER i,j,Face

    REAL(p2) :: faiR,UR
    REAL(p2) :: delta_plus,delta_minus,delta
    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)

    SELECT CASE(Face)
    CASE(1)
        delta_plus = U(i+2,j)-U(i+1,j)
        delta_minus= U(i+1,j)-U(i,j)
    CASE(2)
        delta_plus = U(i,j+2)-U(i,j+1)
        delta_minus= U(i,j+1)-U(i,j)
    CASE(3)
        delta_plus = U(i+1,j)-U(i,j)
        delta_minus= U(i,j)  -U(i-1,j)
    CASE(4)
        delta_plus = U(i,j+1)-U(i,j)
        delta_minus= U(i,j)  -U(i,j-1)
    END SELECT

    CALL Limiter_function(delta_minus,delta_plus,delta,faiR)
    
   SELECT CASE(Face)
    CASE(1)
        UR = U(i+1,j)-delta
    CASE(2)
        UR = U(i,j+1)-delta
    CASE(3)
        UR = U(i,j)-delta
    CASE(4)
        UR = U(i,j)-delta
    END SELECT

    IF (delta_plus==zero)THEN 
        faiR = zero
    END IF
    
    END SUBROUTINE Reconstruct_R