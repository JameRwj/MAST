    SUBROUTINE ReadGrid()
    USE BasicData,ONLY:id,jd,x,y,sxx,sxy,sxMod,nxx,nxy,&
                       syx,syy,syMod,nyx,nyy,vol,half,GridName

    IMPLICIT NONE

    INTEGER :: ni,nj,i,j
    
    OPEN(104,FILE=TRIM(ADJUSTL(GridName)))
    READ(104,*) ni,nj
    DO i = 1,id
        DO j = 1,jd
            READ(104,*)x(i,j),y(i,j)
        END DO
    END DO
    CLOSE(104)
    
    DO i = 1,id
        DO j = 1,jd-1
            sxx(i,j) = y(i,j+1)-y(i,j)
            sxy(i,j) = x(i,j)-x(i,j+1)
            sxMod(i,j)=SQRT(sxx(i,j)**2+sxy(i,j)**2)
            nxx(i,j)=sxx(i,j)/sxMod(i,j)
            nxy(i,j)=sxy(i,j)/sxMod(i,j)
        END DO
    END DO

    DO j = 1,jd
        DO i = 1,id-1
            syx(i,j) = y(i,j)-y(i+1,j)
            syy(i,j) = x(i+1,j)-x(i,j)
            syMod(i,j)=SQRT(syx(i,j)**2+syy(i,j)**2)
            nyx(i,j)=syx(i,j)/syMod(i,j)
            nyy(i,j)=syy(i,j)/syMod(i,j)
        END DO
    END DO

    DO j=1,jd-1
        DO i=1,id-1
            vol(i,j)=half*((x(i+1,j+1)-x(i,j))*(y(i,j+1)-y(i+1,j))&
                          -(x(i,j+1)-x(i+1,j))*(y(i+1,j+1)-y(i,j)))
        END DO
    END DO

    END SUBROUTINE ReadGrid