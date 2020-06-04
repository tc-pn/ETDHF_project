      subroutine FOURIER_TIME(Norb,TEMPSMAX,delta_t,loss_term) 

      IMPLICIT    NONE
      INTEGER     Norb,TEMPSMAX,Nstep
      PARAMETER   (Nstep = 91-1)
cc      PARAMETER   (Nstep = 45-1)
      real*8      delta_t,loss_term(0:TEMPSMAX,Norb)
      real*8      FCT(0:Nstep)
      INTEGER     I,L
      character*2 nb

 204  FORMAT(I2)

      DO L = 1,Norb
         DO I = 1,Nstep+1
            FCT(I-1) = loss_term(I-1,L)
         ENDDO
         call FOURIER(Nstep,FCT,FCT)
         write(nb,204) L
         OPEN(UNIT = 56,file='FOUR'//nb//'.dat',status='unknown')
         DO I = 1,((Nstep+1)/2)
            write(56,*) DFLOAT(I-1),FCT(I-1) 
         ENDDO
         close(56)
      ENDDO

      end

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine FOURIER(Nstep,FCT,X)     

      IMPLICIT NONE	
	
      INTEGER Nstep,I,IFAIL,Nmax
      PARAMETER (Nmax = 10000)
      real*8  FCT(0:Nstep),X(0:Nstep),Y(0:Nmax)

      EXTERNAL C06EAF

      DO I=0,Nstep
         Y(I) = FCT(I)         
      ENDDO

      call C06EAF(Y,Nstep+1,IFAIL)

      DO I = 0,(Nstep)/2
            IF( I .GT. 0) THEN
                X(I) = DSQRT(Y(I)**2.D0+Y(Nstep-I+1)**2.D0) 
            ELSE
                X(I) = DABS(Y(I))
            ENDIF


      ENDDO

      end


      subroutine diag_H_overlap(Iwave,pas,wave_temp
     &                         ,MF,E1,E2,OVERLAP)  
                           
      IMPLICIT NONE

      INTEGER     Noeud,Norbitalmax,Iwave

      PARAMETER   (Noeud= 2**6-1)
      PARAMETER   (Norbitalmax = 50)

      real*8      MF(0:Noeud),pas,HB,H1,E1,E2,OVERLAP
      real*8      Delta_E,NORM 
      COMPLEX*16  NORM1
      COMPLEX*16  wave_temp(0:Noeud,Norbitalmax) 
      COMPLEX*16  w_temp(0:Noeud,0:Noeud),PROD_SCAL
      INTEGER     N,K,I,J
cccccccccccccccccccccccc
c   NAG VARIABLES      c
cccccccccccccccccccccccc
      real*8 AR(Noeud+1,Noeud+1),AI(Noeud+1,Noeud+1)
      real*8 R(Noeud+1),VR(Noeud+1,Noeud+1),VI(Noeud+1,Noeud+1) 
      real*8 WK1(Noeud+1),WK3(Noeud+1),WK2(Noeud+1)
      
      INTEGER IFAIL,ICOUNT

      EXTERNAL F02AXF

      DATA   HB/20.7335D0/

      H1 = HB/pas**2

      Delta_E = 2*E2-2*E1

      IFAIL = 1

ctttttttttttttttttttteeeeeeeeeeeeeeeeessssssssssssssssssttttttttttttttt

      DO I = 1,Iwave
            NORM = 0.
            DO K = 0,Noeud
               NORM  = NORM + CDABS(WAVE_TEMP(K,I))**2.D0
            ENDDO

            IF((IDINT(1000*(NORM*pas)) .NE. 1000) .AND. 
     &         (IDINT(1000*(NORM*pas)) .NE. 999))  
     &         print *,'PB NORM HF',I,'=',NORM*pas
            IF(I .GT. 2) THEN
               NORM1 = DCMPLX(0.D0,0.D0)
               DO K = 0,Noeud
                 NORM1 = NORM1 + DCONJG(WAVE_TEMP(K,I))
     &                           *WAVE_TEMP(K,I-1)
               ENDDO
               IF(IDINT(1000*CDABS(NORM1)) .NE. 0)  
     &         print *,'PB ORTH HF',I,'=',CDABS(NORM1*pas)
           endif
      ENDDO

ctttttttttttttttttttteeeeeeeeeeeeeeeeessssssssssssssssssttttttttttttttt



      DO N = 1,Noeud+1
         DO K = 1,Noeud+1
            AR(K,N) = 0.D0
            AI(K,N) = 0.D0            
         ENDDO
      ENDDO
      DO K = 1,Noeud+1
         AR(K,K-1) = -H1
         AR(K,K+1) = -H1
         AR(K,K)   = 2*H1+MF(k-1)
      ENDDO

      call F02AXF(AR,Noeud+1,AI,Noeud+1,Noeud+1,R,VR,Noeud+1,VI
     &            ,Noeud+1,WK1,WK2,WK3,IFAIL)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   BE CAREFULL THE EIGENVALUES ARE IN ASCENDING ORDER
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
c  FILLING STATES+NORMALISATION
cccccccccccccccccccccccccccccccccccccc
      
      DO I = 0,Noeud
         NORM = 0.D0
         DO K = 0,Noeud
            W_temp(K,I) = DCMPLX(VR(K+1,I+1),VI(K+1,I+1))
            NORM = NORM + CDABS(W_temp(K,I))**2.D0
         ENDDO
         NORM = 1/DSQRT(NORM*pas)
         DO K = 0,Noeud
            W_temp(K,I) = W_temp(K,I)*NORM 
         ENDDO
      ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CALCULATION OF THE OVERLAP                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccc LOOPS ON THE WAVES OF THE BASE


      OVERLAP = 0.D0
      ICOUNT  = 0
      DO I = 0,Noeud
         IF(R(I+1) .LE. (E1+Delta_E)) THEN
            ICOUNT  = ICOUNT + 1
            DO J = 1,Iwave
               PROD_SCAL = DCMPLX(0.D0,0.D0)
               DO K = 0,Noeud
                  PROD_SCAL = PROD_SCAL + DCONJG(W_temp(K,I))*
     &                                    wave_temp(K,J) 
               ENDDO
               PROD_SCAL = pas*PROD_SCAL
               OVERLAP = OVERLAP + CDABS(PROD_SCAL)**2.D0
            ENDDO
         ENDIF
      ENDDO
      OVERLAP = DSQRT(OVERLAP/MIN0(ICOUNT,IWAVE))
      print *,'OVERLAP =',OVERLAP

      end


