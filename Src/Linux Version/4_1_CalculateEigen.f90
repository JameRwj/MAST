    SUBROUTINE CalculateEigen()
    USE BasicData
    !USE LAPACK95
    
    IMPLICIT NONE
    
    REAL(p2) :: len1,len2,len3,len4,s
    REAL(p2) :: nx1,ny1,nx2,ny2,nx3,ny3,nx4,ny4
    REAL(p2) :: rho_,u_,v_,q_
    
    INTEGER  :: col,row,n,i,j
    INTEGER  :: N_Matrix,INFO,LDA,LDVL,LDVR,LWORK

    REAL(p2),DIMENSION(4,4) :: faiL1_1,faiL2_1,faiL3_1,faiL4_1
    REAL(p2),DIMENSION(4,4) :: faiR1_1,faiR2_1,faiR3_1,faiR4_1
    REAL(p2),DIMENSION(4,4) :: faiL1_2,faiL2_2,faiL3_2,faiL4_2
    REAL(p2),DIMENSION(4,4) :: faiR1_2,faiR2_2,faiR3_2,faiR4_2
    REAL(p2),DIMENSION(4,4) :: faiL1_3,faiL2_3,faiL3_3,faiL4_3
    REAL(p2),DIMENSION(4,4) :: faiR1_3,faiR2_3,faiR3_3,faiR4_3
    REAL(p2),DIMENSION(4,4) :: A1,A2,A3,A4
    REAL(p2),DIMENSION(4,4) :: B1,B2,B3,B4
    REAL(p2),DIMENSION(4,4) :: C1,C2,C3,C4
    REAL(p2),DIMENSION(4,4) :: L1,L2,L3,L4
    REAL(p2),DIMENSION(4,4) :: R1,R2,R3,R4
    REAL(p2),DIMENSION(4,4) :: dW_dU,E

    
    REAL(p2),DIMENSION(:,:),ALLOCATABLE :: Matrix
    REAL(p2),DIMENSION(:)  ,ALLOCATABLE :: WORK
    
    row = jd-1
    col = id-1
    
    ALLOCATE(Matrix(4*row*col,4*col*row))
    
    Matrix = zero
    A1     = zero
    A2     = zero
    A3     = zero
    A4     = zero
    B1     = zero
    B2     = zero
    B3     = zero
    B4     = zero
    C1     = zero
    C2     = zero
    C3     = zero
    C4     = zero
    L1     = zero
    L2     = zero
    L3     = zero
    L4     = zero
    R1     = zero
    R2     = zero
    R3     = zero
    R4     = zero
    dW_dU  = zero
    
    DO i = 1,4
        E(i,i) = one
    END DO
        
    DO j = 1,row
        DO i = 1,col
            n = i+(j-1)*col
            
            nx1  = nxx(i+1,j)
            ny1  = nxy(i+1,j)
            nx2  = nyx(i,j+1)
            ny2  = nyy(i,j+1)
            nx3  = nxx(i,j)
            ny3  = nxy(i,j)
            nx4  = nyx(i,j)
            ny4  = nyy(i,j)
            len1 = SQRT(sxx(i+1,j+0)**2+sxy(i+1,j+0)**2)
            len2 = SQRT(syx(i+0,j+1)**2+syy(i+0,j+1)**2)
            len3 = SQRT(sxx(i+0,j+0)**2+sxy(i+0,j+0)**2)
            len4 = SQRT(syx(i+0,j+0)**2+syy(i+0,j+0)**2)
            s    = vol(i,j)
            
            CALL Gradientofflux_L(i,j,1,nx1,ny1,L1,faiL1_1,faiL1_2,faiL1_3)
            CALL Gradientofflux_L(i,j,2,nx2,ny2,L2,faiL2_1,faiL2_2,faiL2_3)
            CALL Gradientofflux_L(i,j,3,nx3,ny3,L3,faiL3_1,faiL3_2,faiL3_3)
            CALL Gradientofflux_L(i,j,4,nx4,ny4,L4,faiL4_1,faiL4_2,faiL4_3)
                                       
            CALL Gradientofflux_R(i,j,1,nx1,ny1,R1,faiR1_1,faiR1_2,faiR1_3)
            CALL Gradientofflux_R(i,j,2,nx2,ny2,R2,faiR2_1,faiR2_2,faiR2_3)
            CALL Gradientofflux_R(i,j,3,nx3,ny3,R3,faiR3_1,faiR3_2,faiR3_3)
            CALL Gradientofflux_R(i,j,4,nx4,ny4,R4,faiR4_1,faiR4_2,faiR4_3)
            
            rho_ = rho(i,j)
            u_   = u(i,j)
            v_   = v(i,j)
            q_   = u_**2+v_**2
            
            dW_dU = RESHAPE([one  , -u_/rho_ , -v_/rho_ , (gamma-one)*q_/two,&
                             zero , one/rho_ , zero     , -(gamma-one)*u_,&
                             zero , zero     , one/rho_ , -(gamma-one)*v_,&
                             zero , zero     , zero     , gamma-one] ,[4,4])

            SELECT CASE (Reconstruction_Method)
            CASE(1)
                A1 = MATMUL(L1,(E+faiL1_1)) * len1/s
                A2 = MATMUL(L2,(E+faiL2_1)) * len2/s
                A3 = MATMUL(R3,(E+faiR3_1)) * len3/s
                A4 = MATMUL(R4,(E+faiR4_1)) * len4/s
                
                B1 = (len1*MATMUL(R1,(E+faiR1_1)) - len3*MATMUL(R3,faiR3_1)) / s
                B2 = (len2*MATMUL(R2,(E+faiR2_1)) - len4*MATMUL(R4,faiR4_1)) / s
                B3 = (len3*MATMUL(L3,(E+faiL3_1)) - len1*MATMUL(L1,faiL1_1)) / s
                B4 = (len4*MATMUL(L4,(E+faiL4_1)) - len2*MATMUL(L2,faiL2_1)) / s
                
                C1 = MATMUL(R1,faiR1_1) * len1/s
                C2 = MATMUL(R2,faiR2_1) * len2/s
                C3 = MATMUL(L3,faiL3_1) * len3/s
                C4 = MATMUL(L4,faiL4_1) * len4/s
            CASE(2)
                A1 = (MATMUL(L1,faiL1_2) + (MATMUL(R1,faiR1_1))) * len1/s
                A2 = (MATMUL(L2,faiL2_2) + (MATMUL(R2,faiR2_1))) * len2/s
                A3 = (MATMUL(L3,faiL3_1) + (MATMUL(R3,faiR3_2))) * len3/s
                A4 = (MATMUL(L4,faiL4_1) + (MATMUL(R4,faiR4_2))) * len4/s
                
                B1 = (len1*(MATMUL(L1,(faiL1_1)) + MATMUL(R1,faiR1_2))&
                    + len3*MATMUL(R3,faiR3_3)) / s
                B2 = (len2*(MATMUL(L2,(faiL2_1)) + MATMUL(R2,faiR2_2))&
                    + len4*MATMUL(R4,faiR4_3)) / s
                B3 = (len3*(MATMUL(R3,(faiR3_1)) + MATMUL(L3,faiL3_2))&
                    + len1*MATMUL(L1,faiL1_3)) / s
                B4 = (len4*(MATMUL(R4,(faiR4_1)) + MATMUL(L4,faiL4_2))&
                    + len2*MATMUL(L2,faiL2_3)) / s
                
                C1 = -MATMUL(R1,faiR1_3) * len1/s
                C2 = -MATMUL(R2,faiR2_3) * len2/s
                C3 = -MATMUL(L3,faiL3_3) * len3/s
                C4 = -MATMUL(L4,faiL4_3) * len4/s
            END SELECT
            
            A1 = MATMUL(dW_DU,A1)
            A2 = MATMUL(dW_DU,A2)
            A3 = MATMUL(dW_DU,A3)
            A4 = MATMUL(dW_DU,A4)
            B1 = MATMUL(dW_DU,B1)
            B2 = MATMUL(dW_DU,B2)
            B3 = MATMUL(dW_DU,B3)
            B4 = MATMUL(dW_DU,B4)
            C1 = MATMUL(dW_DU,C1)
            C2 = MATMUL(dW_DU,C2)
            C3 = MATMUL(dW_DU,C3)
            C4 = MATMUL(dW_DU,C4)
            
            Matrix(4*n-3:4*n,4*n-3:4*n)= -(A1+A2+A3+A4)
            
            IF (i >= 2) THEN
                Matrix(4*n-3:4*n,4*n-7:4*n-4) = -B3
            END IF
            
            IF (i <= col-1) THEN
                Matrix(4*n-3:4*n,4*n+1:4*n+4) = -B1
            END IF
            
            IF (j >= 2) THEN
                Matrix(4*n-3:4*n,4*n-4*col-3:4*n-4*col) = -B4
            END IF
            
            IF (j <= row-1) THEN
                Matrix(4*n-3:4*n,4*n+4*col-3:4*n+4*col) = -B2
            END IF
            
            IF (i >= 3) THEN
                Matrix(4*n-3:4*n,4*n-11:4*n-8) = C3
            END IF
            
            IF (i <= col-2) THEN
                Matrix(4*n-3:4*n,4*n+5:4*n+8) = C1
            END IF
            
            IF (j >= 3) THEN
                Matrix(4*n-3:4*n,4*n-4*2*col-3:4*n-4*2*col) = C4
            END IF
            
            IF (j <= row-2) THEN
                Matrix(4*n-3:4*n,4*n+4*2*col-3:4*n+4*2*col) = C2
            END IF

        END DO
    END DO
    
    N_Matrix = 4*col*row
    
    LDA  = N_Matrix
    LDVL = N_Matrix
    LDVR = N_Matrix
    
    ALLOCATE(WR(N_Matrix),WI(N_Matrix),WORK(4*N_Matrix),&
        VectorL(LDVL,N_Matrix),VectorR(LDVR,N_Matrix))
    
    WR = 0.0
    WI = 0.0
    WORK = 0.0
    VectorL = 0.0
    VectorR = 0.0
    
    LWORK = 4*N_Matrix
    
    CALL DGEEV("Vectors","Vectors",N_Matrix,Matrix,LDA,WR,WI,VectorL,LDVL,VectorR,LDVR,WORK,LWORK,INFO)
    
    DEALLOCATE(WORK,Matrix)
    END SUBROUTINE CalculateEigen