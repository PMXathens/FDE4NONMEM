SUBROUTINE FDEGL(PARAM,Y0,TF,YF,N,Dim,AMT,Np)

! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                                                   ! PART ONE OF SUBROUTINE-NO MODIFICATIONS FROM USER
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

IMPLICIT NONE
! -----------------------------------------------------------------
                        !DIGITS
! -----------------------------------------------------------------

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
! -----------------------------------------------------------------
                        !SOLVER PARAMETERS
! -----------------------------------------------------------------

REAL (dp)            ::   TF,H
INTEGER            ::   i,ni,N,j,q,K,Dim,W,Np
! -----------------------------------------------------------------
                        !PROBLEM PARAMETERS
! -----------------------------------------------------------------

REAL (dp)           ::   ALPHA,k01
REAL (dp)           ::   OMEGA0 
REAL (dp)           ::   AMT
! -----------------------------------------------------------------
                         !SCALARS
! -----------------------------------------------------------------
REAL (dp)           ::   suma2,psi2,ffv2,ffvd2,y02,gf2,gdf2
! -----------------------------------------------------------------
                         !MATRICES
! -----------------------------------------------------------------
INTEGER             	:: unity(Dim,Dim)
REAL (dp)	:: y(N,Dim),y2(N),yt(Dim)
REAL (dp)           	:: PARAM(Np)
REAL (dp)           	:: suma(Dim)
REAL (dp)           	:: psi(Dim)
REAL (dp)     	:: tempt(Dim),temp2(5),temp(5,Dim)
REAL (dp)           	:: y0(Dim)
REAL (dp)           	:: t(N)
REAL (dp)           	:: omega(N)
REAL (dp)           	:: YF(Dim)
REAL (dp)	:: ffv(Dim),ffvd(Dim,Dim)
REAL (dp)           	:: tes(Dim)
REAL (dp)           	:: gf(Dim),gdf(dim,Dim),Inv_gdf(Dim,Dim)
! -----------------------------------------------------------------
                         !MATRICES MULTI-DOSE
! -----------------------------------------------------------------
INTEGER,save         :: counvar =1
REAL(dp),save     	:: doses(100)
REAL(dp),save         :: tdoses(100)

! -----------------------------------------------------------------
                         !MULTI-DOSE MATRICES
! -----------------------------------------------------------------
IF (tf==0) then

YF=Y0
ELSE

if (AMT/=0) THEN
counvar=counvar+1
doses(counvar)=AMT+doses(COUNVAR-1)
tdoses(counvar)=TF
endif
K=5
H=Tf/N



! -----------------------------------------------------------------
                         !PARAMETERS
! -----------------------------------------------------------------
OMEGA0=1
ALPHA=PARAM(1)

! -----------------------------------------------------------------
                         !UNITY MATRIX
! -----------------------------------------------------------------
unity = 0
do w = 1,Dim
unity(w,w)=1
end DO

! -----------------------------------------------------------------
                         !COEFFICIENTS
! -----------------------------------------------------------------
DO i=1, N
 IF (i==1) THEN
   OMEGA(i)=OMEGA0*(1-(1-ALPHA)/i)
 ELSE
   OMEGA(i)=OMEGA(i-1)*(1-(1-ALPHA)/i)
 ENDIF
END DO

! -----------------------------------------------------------------
                         !TIME MATRIX
! -----------------------------------------------------------------
DO i=1, N
 t(i)=H*i
END DO
! -----------------------------------------------------------------
                         !SOLVER
! -----------------------------------------------------------------

IF (DIM>1) THEN

  IF (counvar .ge. 2) THEN
       DO ni=1,N
	      suma=0.0
    	      IF (ni>1) THEN
            		DO q=1,  (counvar-1) 
                    		IF (tdoses(q)<=t(ni) .AND. t(ni)<tdoses(q+1)) THEN
                        		     temp(1,:)=y(ni-1,:)+doses(q)                    
                    		ELSE IF (tdoses(counvar)<=t(ni)) THEN
                        		     temp(1,:)=y(ni-1,:)+doses(counvar)
                    		ENDIF
            		END DO
            	          y(ni,:)= temp(1,:)                            
            		DO j=1,(ni-1)
                			ffv =  FFUN_ND(t(j),y(j,:),param)
                			suma(:)=suma(:)+omega(ni-j)*ffv(:)                   
            		END DO
    	     ELSE
            	         temp(1,:) = y0(:)
                           y(ni,:)=y0(:) 
    	     END IF
    	     DO q=1,(counvar-1) 
            		IF (tdoses(q)<=t(ni) .AND. t(ni)<tdoses(q+1)) THEN
                			psi(:)=y0(:)+h**alpha*suma(:)+doses(q)         
            		ELSE IF (tdoses(counvar)<=t(ni)) THEN
                			psi(:)=y0(:)+h**alpha*suma(:)+doses(counvar)    
            		ENDIF
     	     END DO

    	    DO i=2,k
        		DO q=1,(counvar-1) 
            			IF (tdoses(q)<=t(ni) .AND. t(ni)<tdoses(q+1)) THEN
                			      temp(i,:) = y(ni,:)+doses(q)                
            			ELSE IF (tdoses(counvar)<=t(ni)) THEN
                			      temp(i,:) = y(ni,:)+doses(counvar)    
            			END IF
        		END DO
        		temp(i,:) = y(ni,:);
        		ffv = FFUN_ND(t(ni),temp(i,:),param) 
        		gf(:) = temp(i,:)-h**(alpha)*ffv(:)-psi(:)  
        		ffvd = FFUNJ_ND(t(ni),temp(i,:),param)          
        		gdf = unity-h**(alpha)*ffvd
        		CALL inverse(gdf,Inv_gdf,Dim)
        		tes = MATMUL(Inv_gdf,gf)  
        		y(ni,:)=temp(i,:)-tes    
    	    END DO
       END DO
       YF=Y(N,:)
  ELSE
    
       DO ni=1,N
	suma=0.0;
        	IF (ni>1) THEN
            	       y(ni,:)= temp(1,:)                            
            		DO j=1,(ni-1)
            		yt=y(j,:)
            		ffv =  FFUN_ND(t(j),yt,param)
            		suma=suma+omega(ni-j)*ffv                 
            		END DO
        	ELSE
            	       temp(1,:) = y0(:)
                         y(ni,:)=y0(:)
        	END IF
    	psi=y0+h**alpha*suma
    
    	DO i=2,k
        
        		temp(i,:) = y(ni,:);
        		tempt=temp(i,:)
        		ffv =FFUN_ND(t(ni),tempt,param) 
        		gf = temp(i,:)-h**(alpha)*ffv-psi
        		ffvd = FFUNJ_ND(t(ni),tempt,param) 
        		gdf = unity-h**(alpha)*ffvd
       		CALL inverse(gdf,Inv_gdf,Dim)
        		tes = MATMUL(Inv_gdf,gf)  
        		y(ni,:)=temp(i,:)-tes
    	END DO
    
       END DO
       YF=Y(N,:)
  END IF
ELSE
       
  IF (counvar .ge. 2) THEN
       DO ni=1,N
                    suma2=0.0
    		IF (ni>1) THEN
            			DO q=1,  (counvar-1) 
                    			IF (tdoses(q)<=t(ni) .AND. t(ni)<tdoses(q+1)) THEN
                        				temp2(1)=y2(ni-1)+doses(q)                    
                    			ELSE IF (tdoses(counvar)<=t(ni)) THEN
                        				temp2(1)=y2(ni-1)+doses(counvar)
                    			ENDIF
            			END DO
            			y2(ni)= temp2(1)                            
            			DO j=1,(ni-1)
                				ffv2 =  FFUN_1D(t(j),y2(j),param)
                				suma2=suma2+omega(ni-j)*ffv2                
            			END DO
    		ELSE
            			temp2(1) = y0(1)
            			y2(ni)=y0(1) 
    		END IF
    		DO q=1,(counvar-1) 
            			IF (tdoses(q)<=t(ni) .AND. t(ni)<tdoses(q+1)) THEN
                				psi2=y0(1)+h**alpha*suma2+doses(q)         
            			ELSE IF (tdoses(counvar)<=t(ni)) THEN
                				psi2=y0(1)+h**alpha*suma2+doses(counvar)    
            			ENDIF
    		END DO

    		DO i=2,k
        			DO q=1,(counvar-1) 
            				IF (tdoses(q)<=t(ni) .AND. t(ni)<tdoses(q+1)) THEN
                					temp2(i) = y2(ni)+doses(q)                
           				ELSE IF (tdoses(counvar)<=t(ni)) THEN
                					temp2(i) = y2(ni)+doses(counvar)    
            				END IF
        			END DO
        			temp2(i) = y2(ni);
        			ffv2 =FFUN_1D(t(ni),temp2(i),param) 
        			gf2 = temp2(i)-h**(alpha)*ffv2-psi2  
    
        			ffvd2 =FFUNJ_1D(t(ni),temp2(i),param)          
        			gdf2 = 1-h**(alpha)*ffvd2      
        			y2(ni)=temp2(i)-gf2/gdf2  
    		END DO
       END DO
      YF=Y2(N)
ELSE

       DO ni=1,N
	suma2=0.0
    		IF (ni>1) THEN
    			y2(ni)= temp2(1)                            
    			DO j=1,(ni-1)
    			ffv2 =  FFUN_1D(t(j),y2(j),param)
    			suma2=suma2+omega(ni-j)*ffv2                  
    			END DO
    		ELSE
    			temp2(1) = y0(1)
    			y2(ni)=y0(1)
    		END IF
    		psi2=y0(1)+h**alpha*suma2
        		DO i=2,k
        
         			temp2(i) = y2(ni);
        			ffv2= FFUN_1D(t(ni),temp2(i),param) 
        			gf2= temp2(i)-h**(alpha)*ffv2-psi2 
        			ffvd2 = FFUNJ_1D(t(ni),temp2(i),param)          
        			gdf2 = 1-h**(alpha)*ffvd2
        			y2(ni)=temp2(i)-gf2/gdf2  
       		END DO

       END DO


      YF=Y2(N)

  END IF
ENDIF
ENDIF
! -----------------------------------------------------------------
! ----------------------------------------------------------------- 

contains

! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                                                   ! PART TWO OF SUBROUTINE-USER DECLARES FFUN AND FFUNJ
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



! ----------------------------------------------------------------- ----------------------------------------------------------------- 
  				!FUCNTIONS FOR SINGLE-COMPARTMENT!
! ----------------------------------------------------------------- ----------------------------------------------------------------- 

function FFUN_1D(t,y,param)
implicit none
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp)  :: ffun_1d
REAL (dp)   :: param(*) 
REAL (dp)   :: t
REAL (dp)   :: y
 
! -----------------------------------------------------------------
                              !local vars
! -----------------------------------------------------------------

REAL (dp)   :: k

! -----------------------------------------------------------------
! ----------------------------------------------------------------- 
k = param(2)
ffun_1d=-k*y

end function FFUN_1D

! -----------------------------------------------------------------
! -----------------------------------------------------------------

function FFUNJ_1D(t,y,param)
implicit none
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp)  :: ffunj_1d
REAL (dp)   :: t
REAL (dp)   :: y  
REAL (dp)   :: param(*)
! -----------------------------------------------------------------
                               !local vars
! -----------------------------------------------------------------

REAL (dp)   :: k

! -----------------------------------------------------------------
! -----------------------------------------------------------------
k=param(2)
ffunj_1d=-k

end function FFUNJ_1D


! ----------------------------------------------------------------- ----------------------------------------------------------------- 
  				!FUCNTIONS FOR MULTI-COMPARTMENT!
! ----------------------------------------------------------------- ----------------------------------------------------------------- 

function FFUN_ND(t,y,param)
implicit none
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


REAL (dp)   :: param(*)
REAL (dp)   :: t
REAL (dp)   :: y(*)
! -----------------------------------------------------------------
                              !local vars
! -----------------------------------------------------------------

REAL (dp)   :: k1,k2
REAL(dp)    :: ffun_nd(2)

! -----------------------------------------------------------------
! ----------------------------------------------------------------- 
k1 = param(2)
k2=param(3)

ffun_nd(1)=-k1*y(1)
ffun_nd(2)=k1*y(1)-k2*y(2)
end function FFUN_ND

! -----------------------------------------------------------------
! -----------------------------------------------------------------

function FFUNJ_ND (t,y,param)
implicit none
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp)   :: param(*)  
REAL (dp)   :: t
REAL (dp)   :: y(*)

! -----------------------------------------------------------------
                               !local vars
! -----------------------------------------------------------------

REAL (dp)   :: k1,k2
REAL (dp)  :: ffunj_nd(2,2)

! -----------------------------------------------------------------
! -----------------------------------------------------------------
k1=param(2)
k2=param(3)

ffunj_nd(1,1)=-k1
ffunj_nd(1,2)=0
ffunj_nd(2,1)=k1
ffunj_nd(2,2)=-k2

end function FFUNJ_ND
end subroutine FDEGL
! -----------------------------------------------------------------
                               ! SUBROUTINE FOR INVERSE MATRIX CALCULATION!
! -----------------------------------------------------------------
subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse