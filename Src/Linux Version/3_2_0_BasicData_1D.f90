    MODULE BasicData_1D

    INTEGER,PARAMETER :: p2 = SELECTED_REAL_KIND(P=15) !Double precision 15

    REAL(p2),PARAMETER:: cfl_1D=0.2_p2
    INTEGER,PARAMETER ::nst=1000
    
    INTEGER :: factor,ncount,maxCount

    REAL(p2) :: stageCoefficient(3),dt

    REAL(p2),ALLOCATABLE :: x_1D(:),rho_1D(:),u_1D(:),p_1D(:),DeltaX(:),deltaT(:)        
    REAL(p2),ALLOCATABLE :: qq(:,:),dqq(:,:),rhs(:,:),qq_0(:,:)
    REAL(p2),ALLOCATABLE :: residual(:,:),resid(:)

    END MODULE BasicData_1D