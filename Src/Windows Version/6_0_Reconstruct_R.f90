    SUBROUTINE Reconstruct_R(i,j,U,Face,faiR_1,faiR_2,faiR_3,UR)
    USE BasicData,ONLY: Reconstruction_Method,p2,id,jd,ghostLayers,zero
    
    IMPLICIT NONE
    
    INTEGER i,j,Face
    REAL(p2) :: faiR_1,faiR_2,faiR_3,UR
    REAL(p2) :: U(1-ghostLayers:id+ghostLayers-1,1-ghostLayers:jd+ghostLayers-1)
    
    faiR_1 = zero
    faiR_2 = zero
    faiR_3 = zero
    
    SELECT CASE(Reconstruction_Method)
    CASE(1)
        CALL MUSCL_R(i,j,U,Face,faiR_1,UR) 
    CASE(2)
        CALL ROUND_R(i,j,U,Face,faiR_1,faiR_2,faiR_3,UR)
    END SELECT
    
    END SUBROUTINE Reconstruct_R