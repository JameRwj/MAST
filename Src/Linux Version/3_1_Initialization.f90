    SUBROUTINE Initialization()

    USE BasicData,ONLY:x,y,rho,u,v,p,id,jd,epsilon,Mainf,zero,&
        one,two,id,jd,ghostLayers,p2,gamma,Initialization_Method,TestCase
    USE BasicData_1D,ONLY:rho_1D,u_1D,p_1D

    IMPLICIT NONE

    INTEGER :: i,j,ii

    REAL(p2) :: rhoL,uL,vL,pL
    REAL(p2) :: rhoR,uR,vR,pR
    REAL(p2) :: rhoM,uM,vM,pM

    REAL(p2) :: epsilonU,epsilonP

    SELECT CASE (TestCase)
    CASE(1)

        SELECT CASE (Initialization_Method)
        CASE(1)
            rhoL = one
            uL   = one
            vL   = zero
            pL   = one/gamma/mainf/mainf

            rhoR = one/(two/( (gamma+one)*mainf*mainf )+&
                (gamma-one)/(gamma+one))
            uR   = two/((gamma+one)*mainf*mainf)+(gamma-one)/(gamma+one)
            vR   = zero
            pR   = (two*gamma*mainf*mainf/(gamma+one)-(gamma-one)/&
                (gamma+one))/gamma/mainf/mainf

            epsilonU = one-(one-epsilon)/&
                sqrt(one+epsilon*(mainf*mainf-one)/(one+(gamma-one)*mainf*mainf/two))/&
                sqrt(one+epsilon*(mainf*mainf-one)/(one-two*gamma*mainf*mainf/(gamma-one)))
            epsilonP = epsilon/sqrt(one+(one-epsilon)*(gamma+one)/(gamma-one)*&
                (mainf*mainf-one)/mainf/mainf)

            rhoM = (one-epsilon)*rhoL+epsilon*rhoR
            uM   = (one-epsilonU)*uL+epsilonU*uR
            vM   = zero
            pM   = (one-epsilonP)*pL+epsilonP*pR

            ii = FLOOR(0.5*id)
            DO i = 1-ghostLayers,ii-1
                DO j = 1-ghostLayers,jd-1+ghostLayers
                    rho(i,j) = rhoL
                    u(i,j)   = uL
                    v(i,j)   = vL
                    p(i,j)   = pL
                END DO
            END DO

            DO j = 1-ghostLayers,jd-1+ghostLayers
                rho(ii,j) = rhoM
                u(ii,j)   = uM
                v(ii,j)   = vM
                p(ii,j)   = pM
            END DO
            
            DO i = ii+1,id-1+ghostLayers
                DO j = 1-ghostLayers,jd-1+ghostLayers
                    rho(i,j) = rhoR
                    u(i,j)   = uR
                    v(i,j)   = vR
                    p(i,j)   = pR
                END DO
            END DO

        CASE(2)

            CALL OneDShock()

            rhoL=one
            uL=one
            vL=zero
            pL=one/gamma/mainf/mainf

            rhoR=one/(two/( (gamma+one)*mainf*mainf )+(gamma-one)/(gamma+one))
            uR=two/((gamma+one)*mainf*mainf)+(gamma-one)/(gamma+one)
            vR=zero
            pR=(two*gamma*mainf*mainf/(gamma+one)-(gamma-one)/(gamma+one))/gamma/mainf/mainf
            
            DO i=1,id-1
                rho(i,1)= rho_1D(i)
                u(i,1)  = u_1D(i)
                p(i,1)  = p_1D(i)
            END DO

            ii=FLOOR(0.5*id)
            DO i=1,ii-1
                IF (ABS(rho(i,1)-rhoL)<1.0E-7) THEN
                    rho(i,1)=rhoL
                END IF
                IF (ABS(u(i,1)-uL)<1.0E-7) THEN
                    u(i,1)=uL
                END IF
                IF (ABS(p(i,1)-pL)<1.0E-7) THEN
                    p(i,1)=pL
                END IF
            END DO
            
            DO i=ii+1,id-1
                IF (ABS(rho(i,1)-rhoR)<1.0E-7) THEN
                    rho(i,1)=rhoR
                END IF
                IF (ABS(u(i,1)-uR)<1.0E-7) THEN
                    u(i,1)=uR
                END IF
                IF (ABS(p(i,1)-pR)<1.0E-7) THEN
                    p(i,1)=pR
                END IF
            END DO

            CLOSE(101)
            CLOSE(102)
            CLOSE(103)
            DO i=1-ghostLayers,0
                rho(i,1)=rho(1,1)
                u(i,1)=u(1,1)
                p(i,1)=p(1,1)
            END DO
            DO i=id,id-1+ghostLayers
                rho(i,1)=rho(id-1,1)
                u(i,1)=u(id-1,1)
                p(i,1)=p(id-1,1)
            END DO
            DO j=1-ghostLayers,jd-1+ghostLayers
                DO i=1-ghostLayers,id-1+ghostLayers
                    rho(i,j)=rho(i,1)
                    u(i,j)=  u(i,1)
                    p(i,j)=  p(i,1)
                END DO
            END DO

            DO j=1-ghostLayers,jd-1+ghostLayers
                DO i=1-ghostLayers,id-1+ghostLayers
                    rho(i,j) = rho(i,j)
                    u(i,j)   = u(i,j)
                    v(i,j)   = zero
                    p(i,j)   = p(i,j)
                END DO
            END DO
        END SELECT

        
    CASE(2)
        OPEN(113,FILE="InitialFlow_rho.dat")
        OPEN(114,FILE="InitialFlow_u.dat")
        OPEN(115,FILE="InitialFlow_v.dat")
        OPEN(116,FILE="InitialFlow_p.dat")
        DO i = 1,id-1
            DO j = 1,jd-1
                READ(113,*)rho(i,j)
                READ(114,*)u(i,j)
                READ(115,*)v(i,j)
                READ(116,*)p(i,j)
            END DO
        END DO

        DO i = 1,id-1
            DO j = 1-ghostLayers,0
                rho(i,j) = rho(i,1)
                u(i,j)   = u(i,1)
                v(i,j)   = v(i,1)
                p(i,j)   = p(i,1)
            END DO
        END DO
        DO i = 1,id-1
            DO j = jd,jd+ghostLayers-1
                rho(i,j) = rho(i,jd-1)
                u(i,j)   = u(i,jd-1)
                v(i,j)   = v(i,jd-1)
                p(i,j)   = p(i,jd-1)
            END DO
        END DO
        DO i = 1-ghostLayers,0
            DO j = 1,jd-1
                rho(i,j) = rho(1,j)
                u(i,j)   = u(1,j)
                v(i,j)   = v(1,j)
                p(i,j)   = p(1,j)
            END DO
        END DO
        DO i = id,id+ghostLayers-1
            DO j = 1,jd-1
                rho(i,j) = rho(id-1,j)
                u(i,j)   = u(id-1,j)
                v(i,j)   = v(id-1,j)
                p(i,j)   = p(id-1,j)
            END DO
        END DO

    END SELECT

    !OPEN(113,FILE="Results/FlowField.plt",FORM="FORMATTED")
    !WRITE(113,*) "TITLE             = result"
    !WRITE(113,*) "VARIABLES=   ", "x   ","y   ","rho   ","u   ","v   ","p   "
    !WRITE(113,*) "ZONE T=","ZONE1"
    !WRITE(113,*) "STRANDID=0, SOLUTIONTIME=0"
    !WRITE(113,*) "I= ",id," J= ",jd," K= ",1
    !WRITE(113,*) "DATAPACKING=BLOCK"
    !WRITE(113,*) "VARLOCATION=([3-6]=CELLCENTERED)"

    !DO j=1,jd
    !    DO i = 1,id
    !        WRITE(113,*)x(i,j)
    !    END DO
    !END DO
    !DO j=1,jd
    !    DO i = 1,id
    !        WRITE(113,*)y(i,j)
    !    END DO
    !END DO
    !DO j=1,jd-1
    !    DO i = 1,id-1
    !        WRITE(113,*)rho(i,j)
    !    END DO
    !END DO
    !DO j=1,jd-1
    !    DO i = 1,id-1
    !        WRITE(113,*)u(i,j)
    !    END DO
    !END DO
    !DO j=1,jd-1
    !    DO i = 1,id-1
    !        WRITE(113,*)v(i,j)
    !    END DO
    !END DO
    !DO j=1,jd-1
    !    DO i = 1,id-1
    !        WRITE(113,*)p(i,j)
    !    END DO
    !END DO
    !
    !CLOSE(113)

    
    END SUBROUTINE Initialization