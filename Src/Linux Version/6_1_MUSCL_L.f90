    SUBROUTINE MUSCL_L(i,j,U,Face,faiL,UL)
    USE BasicData,ONLY: p2,id,jd,ghostLayers,zero

    IMPLICIT NONE

    INTEGER i,j,Face

    REAL(p2) :: faiL,UL
    REAL(p2) :: delta_plus,delta_minus,delta
    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)

    SELECT CASE(Face)
    CASE(1)
        delta_plus = U(i+1,j)-U(i,j)
        delta_minus= U(i,j)  -U(i-1,j)
    CASE(2)
        delta_plus = U(i,j+1)-U(i,j)
        delta_minus= U(i,j)  -U(i,j-1)
    CASE(3)
        delta_plus = U(i,j)  -U(i-1,j)
        delta_minus= U(i-1,j)-U(i-2,j)
    CASE(4)
        delta_plus = U(i,j)  -U(i,j-1)
        delta_minus= U(i,j-1)-U(i,j-2)
    END SELECT

    CALL Limiter_function(delta_plus,delta_minus,delta,faiL)
    
    SELECT CASE(Face)
    CASE(1)
        UL = U(i,j)+delta
    CASE(2)
        UL = U(i,j)+delta
    CASE(3)
        UL = U(i-1,j)+delta
    CASE(4)
        UL = U(i,j-1)+delta
    END SELECT
    
    IF (delta_minus==zero)THEN 
        faiL = zero
    END IF
    
    END SUBROUTINE MUSCL_L