    MODULE BasicData

    INTEGER,PARAMETER :: p2 = SELECTED_REAL_KIND(P=15) !P=15 Double precision; p=6 is single precision

    CHARACTER*50 :: GridName
    
    INTEGER :: id,jd  
    INTEGER :: Limiter,Riemann_Solver,Initialization_Method,TestCase,Reconstruction_Method
    
    REAL(p2):: epsilon,Mainf
    
    INTEGER,PARAMETER :: ghostLayers = 2

    REAL(p2),PARAMETER ::  zero = 0.0_p2
    REAL(p2),PARAMETER ::   one = 1.0_p2
    REAL(p2),PARAMETER ::   two = 2.0_p2
    REAL(p2),PARAMETER ::  half = 0.5_p2
    REAL(p2),PARAMETER :: gamma = 1.4_p2
 
    REAL(p2),ALLOCATABLE :: x(:,:),y(:,:)
    REAL(p2),ALLOCATABLE :: sxx(:,:),sxy(:,:),syx(:,:),syy(:,:),vol(:,:)
		
    REAL(p2),ALLOCATABLE :: nxx(:,:),nxy(:,:),nyx(:,:),nyy(:,:)
    REAL(p2),ALLOCATABLE :: sxMod(:,:),syMod(:,:)

    REAL(p2),ALLOCATABLE :: rho(:,:),u(:,:),v(:,:),p(:,:)
    
    REAL(p2),DIMENSION(:,:),ALLOCATABLE :: VectorR,VectorL
    REAL(p2),DIMENSION(:)  ,ALLOCATABLE :: WR,WI

    END MODULE BasicData