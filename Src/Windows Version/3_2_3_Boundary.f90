    SUBROUTINE Boundary()

    USE BasicData,ONLY: id,gamma,mainf,p2,ghostlayers,one
    USE BasicData_1D,ONLY: rho_1D,u_1D,p_1D

    IMPLICIT NONE

    INTEGER :: i

    !    !left boundary
    DO i=0,-1,1-ghostLayers
        rho_1D(i) = one
        u_1D(i)   = one
        p_1D(i)   = one/gamma/mainf/mainf
    END DO        

    DO i=id,id+ghostlayers-1
        !right Boundary
        rho_1D(i) = rho_1D(id-1)
        u_1D(i)   = one/rho_1D(id)
        p_1D(i)   = p_1D(id-1)
    END DO

    END SUBROUTINE Boundary