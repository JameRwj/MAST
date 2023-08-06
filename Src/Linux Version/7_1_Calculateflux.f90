    SUBROUTINE Calculateflux(UL,UR,nx,ny,flux)
    USE BasicData,ONLY:p2,Riemann_solver

    IMPLICIT NONE

    REAL(p2),DIMENSION(4) :: UL,UR,flux
    REAL(p2) :: nx,ny

    SELECT CASE(Riemann_solver)
    CASE(1)
        CALL Roe(UL,UR,nx,ny,flux)
    CASE(2)
        CALL HLLC(UL,UR,nx,ny,flux)
    CASE(3)
        CALL HLL(UL,UR,nx,ny,flux)
    CASE(4)
        CALL van_Leer(UL,UR,nx,ny,flux)
    CASE(5)
        CALL AUSMplus(UL,UR,nx,ny,flux)
    CASE(6)
        CALL SLAU(UL,UR,nx,ny,flux)
    CASE(7)
        CALL HLLE(UL,UR,nx,ny,flux)
    CASE(8)
        CALL HLLEM(UL,UR,nx,ny,flux)
        
    END SELECT
    
    END SUBROUTINE Calculateflux