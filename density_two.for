ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     program which gives many different two body
c     densities but with a fixed one-body density
c
c     "Exact Resolution"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine diagonalize2d(POT,pas,Nwave,FORCE)   

      IMPLICIT NONE

      INTEGER   Nmax,Noeud,N_ham,Nwave,Ntot
      parameter (Noeud=2**6-1,N_ham = (Noeud+1)**2)
      parameter (Ntot = N_ham*(N_ham+1)/2,Nmax=400)
      real*8    POT(0:Noeud,0:Noeud),pas,FORCE
      real*8    HC(Ntot),DENS_TEMP,Xmin
      real*8    LAMBDA,HB,MASS,H1,H2,BETA,P(N_ham),Ptot
      INTEGER   I,J,K,L,I_temp,J_temp,I_centre  
      character*2 nb     


      real*8    TOL,D(N_ham),E(N_ham),E2(N_ham)
      real*8    ALB,UB,EPS,EPS1,V(N_ham,Nmax),R(Nmax)     
      INTEGER   N_reel,ICOUNT(Nmax),IFAIL,M1,M2
      real*8    X(N_ham,7)
      LOGICAL   C(N_ham)
      real*8    X02AJF  

      EXTERNAL F01AYF
      EXTERNAL  F02BEF
      EXTERNAL F01AZF

      DATA HB/197.32705D0/,MASS/938.91897D0/

 110  FORMAT(I2)
      ALB = -200.
c      UB = (HB*3.1415D0/pas)**2.D0
c     &        /(2.D0*938.D0)/4.D0
      UB = +30.
      EPS = X02AJF()
      EPS1= 0.0E0
      IFAIL = 0
 
      H1 = 4*HB**2.D0/(2.D0*MASS)/pas**2.D0
      H2 = -H1/4
      I_centre = Noeud/2


ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONSRAINED FIELD - LAMBDA R^2 = - LAMBA(R1^2 + R2^2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Xmin=2.D0
      LAMBDA = FORCE*2.D0/pas**2.D0*Xmin**2.D0  
      LAMBDA=  LAMBDA*2.D0 + FORCE*4.D0/pas**2.D0


c      LAMBDA = 1.D0*pas
ccccccccccccccccccccccccccccccccccccccccc
c             TEMPERATURE
ccccccccccccccccccccccccccccccccccccccccc

      BETA = 1/5D0
c      BETA = 1.D0
c
c  TEST-ATTENTION
c
c       BETA = 1/0.1D0
c       LAMBDA = 0.
ccccccccccccccccccccccccccccccccccccccccc
c   FILLING THE CONSTRAINED HAMILTONIEN
ccccccccccccccccccccccccccccccccccccccccc

c      OPEN(UNIT=12,file='POTass.DAT',status='unknown')
c      DO I = 0,Noeud
c         IF(I .LE. I_centre) THEN
c            write(12,*) pas*DFLOAT(I),(POT(I,I)
c     &                  -lambda/2*(DFLOAT(I-I_centre))**2.D0)
c         ELSE
c            write(12,*) pas*DFLOAT(I),(POT(I,I)
c     &                  -lambda*(DFLOAT(I-I_centre))**2.D0)
c         ENDIF
cc      ENDDO
c      close(12)
c      PAUSE

      DO I = 0,Noeud
         DO J = 0,Noeud
            DO K = 0,Noeud
               DO L = 0,Noeud 
                  I_TEMP  = J+1+(Noeud+1)*I
                  J_TEMP  = L+1+(Noeud+1)*K
                  IF(J_TEMP .LE. I_TEMP) THEN
c                    HC(J+1+(Noeud+1)*I,L+1+(Noeud+1)*K) = 0.D0
                     HC(I_TEMP*(I_TEMP-1)/2+J_TEMP) = 0.D0
                  ENDIF
               ENDDO               
           ENDDO
        ENDDO
      ENDDO

      DO I = 0,Noeud
         DO J = 0,Noeud
            DO K = 0,Noeud
               DO L = 0,Noeud 
                  IF ((I .EQ. K) .AND. (J .EQ. L)) THEN
                      I_TEMP  = J+1+(Noeud+1)*I
                      J_TEMP  = L+1+(Noeud+1)*K
                      IF(J_TEMP .LE. I_TEMP) THEN
c     &                   HC(J+1+(Noeud+1)*I,L+1+(Noeud+1)*K) = 
                          HC(I_TEMP*(I_TEMP-1)/2+J_TEMP) = 
     &                       H1 + POT(I,J)-LAMBDA*
     &                       DFLOAT((I-I_centre)**2+(J-I_centre)**2)
cccccccccccccccccccccccccccccccccccccccccc
c    ASSYMETRIC CASE ONLY
cccccccccccccccccccccccccccccccccccccccccc
                          IF(I .LE. I_centre) THEN 
                             HC(I_TEMP*(I_TEMP-1)/2+J_TEMP)=
     &                          HC(I_TEMP*(I_TEMP-1)/2+J_TEMP)
     &                          +LAMBDA/2*
     &                          DFLOAT((I-I_centre)**2)
                          ENDIF
                          IF(J .LE. I_centre) THEN 
                             HC(I_TEMP*(I_TEMP-1)/2+J_TEMP)=
     &                          HC(I_TEMP*(I_TEMP-1)/2+J_TEMP)
     &                          +LAMBDA/2*
     &                          DFLOAT((J-I_centre)**2)
                          ENDIF
                      ENDIF
                  ENDIF

                  IF ((I .EQ. K) .AND. (IABS(J-L) .EQ. 1)) THEN
                      I_TEMP  = J+1+(Noeud+1)*I
                      J_TEMP  = L+1+(Noeud+1)*K
                      IF(J_TEMP .LE. I_TEMP) THEN  
                         HC(I_TEMP*(I_TEMP-1)/2+J_TEMP) = H2
c     &                  HC(J+1+(Noeud+1)*I,L+1+(Noeud+1)*K) = H2
                      ENDIF
                  ENDIF
                  IF ((J .EQ. L) .AND. (IABS(I-K) .EQ. 1)) THEN
                      I_TEMP  = J+1+(Noeud+1)*I
                      J_TEMP  = L+1+(Noeud+1)*K
                      IF(J_TEMP .LE. I_TEMP) THEN     
                         HC(I_TEMP*(I_TEMP-1)/2+J_TEMP) = H2 
c     &                  HC(J+1+(Noeud+1)*I,L+1+(Noeud+1)*K) = H2
                      ENDIF
                  ENDIF

               ENDDO
            ENDDO
         ENDDO
      ENDDO

cccccccccccccccccccccccccccccccccccccc
c     TRIANGULATION OF THE MATRICE
cccccccccccccccccccccccccccccccccccccc

      print *,'coucou0'

c      call F01AGF(N_ham,TOL,HC,N_ham,D,E,E2)
       call F01AYF(N_ham,TOL,HC,Ntot,D,E,E2)
 
      print *,'coucou1'


      call F02BEF(N_ham,D,ALB,UB,EPS,EPS1,E,E2,Nmax,N_reel,R,V
     &            ,N_ham,ICOUNT,X,C,IFAIL)

ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Calculation of the density probablities
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      print *,'Nreel =',N_reel

      Ptot = 0.D0
      DO I = 1,N_reel
          P(I) = DEXP(-BETA*R(I))
          Ptot = Ptot + P(I)
      ENDDO

      DO I = 1,N_reel
         PRINT*,'R(',I,')=',R(I)
      enddo

      DO I = 1,N_reel
         P(I) = P(I)/Ptot
c         IF (P(I) .GT. (0.00001D0)) Nwave = I
        Nwave = I  
      ENDDO

      PRINT*,'Nwave=',Nwave 

      M1 = 1
      M2 = Nwave 

      print *,'coucou2'

c      call F01AHF(N_ham,M1,M2,HC,N_ham,E,V,N_ham)
      call F01AZF(N_ham,M1,M2,HC,Ntot,V,N_ham)

      IF( Nwave .LE. Nmax ) THEN
         OPEN(UNIT=21,file='PIass.DAT',status='unknown')
          DO K = 1,Nwave
             write(21,*) P(k)
             write(nb,110) K
             OPEN(UNIT=22,file='PSIass'//nb//'.dat',status='unknown')
             DO I = 0,Noeud
                DO J = 0,Noeud
                   write(22,*) SNGL(V(J+1+(Noeud+1)*I,K)/pas)
                ENDDO
             ENDDO
             close(22)
          ENDDO
          close(21)
      ELSE
c         PRINT*,'PB-TO MUCH WAVE_PACKET!...'
      ENDIF
      print *,'coucou'
c
c  ONE-BODY DENSITY CALCULATION FOR TDHF
c
      OPEN(UNIT=22,file='DENSITY_ONEass.DAT',STATUS='UNKNOWN')
      OPEN(UNIT = 20,file='DENSITY_ONEbass.DAT',STATUS='UNKNOWN')
      DO I = 0,Noeud
         DO J = 0,Noeud
         DENS_TEMP = 0.D0
            DO K = 1,Nwave
               DO L = 0,Noeud
                  DENS_TEMP = DENS_TEMP + P(K)*V(L+1+(Noeud+1)*I,K)/pas
     &                                    *V(L+1+(Noeud+1)*J,K)/pas
               ENDDO
            ENDDO
            DENS_TEMP = DENS_TEMP*pas
            write(20,*) SNGL(DENS_TEMP)  
            IF(I .EQ. J) WRITE(22,*) SNGL(DFLOAT(J)*pas)
     &                  ,SNGL(DENS_TEMP)  
         ENDDO
      ENDDO
      CLOSE(20)
      CLOSE(22)
      end

