    SUBROUTINE Output()

    USE BasicData,ONLY:p2,WR,WI,id,jd,VectorR,VectorL,zero

    IMPLICIT NONE

    INTEGER :: N_Matrix,i,j,k,n,m,ii,jj
    REAL(p2):: max_eigvalue

    CHARACTER*50 :: fileName
    CHARACTER*50 :: number

    REAL(p2),DIMENSION(:),ALLOCATABLE :: eigVector_rho,eigVector_u,eigVector_v,eigVector_p
    INTEGER ,DIMENSION(:),ALLOCATABLE :: position

    ALLOCATE(eigVector_rho((id-1)*(jd-1)),eigVector_u((id-1)*(jd-1)),eigVector_v((id-1)*(jd-1)),eigVector_p((id-1)*(jd-1)),position(4*(id-1)*(jd-1)))

    N_Matrix = SIZE(WR,1)

    OPEN(115,FILE="Results\Scatter.plt",FORM="FORMATTED")
    DO i = 1,N_Matrix
        WRITE(115,*)WR(i),WI(i)
    END DO
    CLOSE(115)

    Max_eigvalue = MAXVAL(WR)
    k = 0
    position = zero
    DO i = 1,N_Matrix
        IF (WR(i) == Max_eigvalue) THEN
            k = k + 1
            position(k) = i
        END IF
    END DO

    WRITE(*,"(A80)")"=====================================Result====================================="
    Write(*,"(A26)")"The maximal eigenvalue is:"
    DO i = 1,k
        WRITE(*,"(A26,F16.8,A4,F16.8,A2)")" ",WR(position(i)),"+",WI(position(i)),"i"
    END DO

    IF (Max_eigvalue>zero) THEN
        DO i = 1,k
            DO j = 1,N_Matrix,4
                n = 1 + FLOOR(j/4.0)
                eigVector_rho(n) = VectorR(j  ,position(i))
                eigVector_u(n)   = VectorR(j+1,position(i))
                eigVector_v(n)   = VectorR(j+2,position(i))
                eigVector_p(n)   = VectorR(j+3,position(i))

            END DO

            WRITE(number,*)i

            fileName = "Results\Eigenvector_rho_"//TRIM(ADJUSTL(number))//'.plt'
            OPEN(115,FILE=TRIM(ADJUSTL(fileName)),FORM="FORMATTED")
            WRITE(115,*)' variables="x","y","<greek>r</greek>" '
            WRITE(115,"(A10,A10,I2,A10,I2)")"zone","i=",id-1,"j=",jd-1

            fileName = "Results\Eigenvector_u_"//TRIM(ADJUSTL(number))//'.plt'
            OPEN(116,FILE=TRIM(ADJUSTL(fileName)),FORM="FORMATTED")
            WRITE(116,*)' variables="x","y","u" '
            WRITE(116,"(A10,A10,I2,A10,I2)")"zone","i=",id-1,"j=",jd-1

            fileName = "Results\Eigenvector_v_"//TRIM(ADJUSTL(number))//'.plt'
            OPEN(117,FILE=TRIM(ADJUSTL(fileName)),FORM="FORMATTED")
            WRITE(117,*)' variables="x","y","v" '
            WRITE(117,"(A10,A10,I2,A10,I2)")"zone","i=",id-1,"j=",jd-1

            fileName = "Results\Eigenvector_p_"//TRIM(ADJUSTL(number))//'.plt'
            OPEN(118,FILE=fileName,FORM="FORMATTED")
            WRITE(118,*)' variables="x","y","p" '
            WRITE(118,"(A10,A10,I2,A10,I2)")"zone","i=",id-1,"j=",jd-1

            m = 1
            DO ii = 1,jd-1
                DO jj = 1,id-1
                    WRITE(115,*)jj,ii,eigVector_rho(m)
                    WRITE(116,*)jj,ii,eigVector_u(m)
                    WRITE(117,*)jj,ii,eigVector_v(m)
                    WRITE(118,*)jj,ii,eigVector_p(m)
                    m = m+1
                END DO
            END DO

            CLOSE(115)
            CLOSE(116)
            CLOSE(117)
            CLOSE(118)
        END DO
    END IF

    DEALLOCATE(eigVector_rho,eigVector_u,eigVector_v,eigVector_p,position)
    END SUBROUTINE Output