      subroutine INITIALISATION2(Norbitalmax,wave,pas,pot,K,t0,delta_t
     &                    ,I_TEMPSMAX,ni,Norbital)   
                                         
      IMPLICIT NONE

      INTEGER     Noeud,Nmax
      parameter   (Noeud = 2**6-1,Nmax = Noeud+1)

      integer     A,I_TEMPSMAX,Norbitalmax,Norbital   
      COMPLEX*16  WAVE(0:Noeud,Norbitalmax)     
      real*8      pas,K,t0,FORCE,HB,MASS,PI,ni(Norbitalmax)
      real*8      delta_t,Emax,delta_p
      real*8      POT(0:Noeud),Y,X
      real*8      DENS(Nmax,Nmax)
      real*8      FREQ,HB2,R02,FORCE1  
      INTEGER     I,J,I_centre

ccccccccccccccccccccccccccccccccccccccccccccccccc
c     NAG VARIABLES
ccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL*8   TOL,EPS
      REAL*8   D(Nmax),E(Nmax)
      REAL*8   V(Nmax,Nmax) 
      REAL*8   X02AJF
      INTEGER  IFAIL


      EXTERNAL X02AJF
      EXTERNAL F01AJF
      EXTERNAL F02AMF

      DATA HB/197.32705D0/,MASS/938.91897D0/
c
c   NB OF ITERATION WITH TIME EVOLUTION
c
c      I_TEMPSMAX = 1000
      I_TEMPSMAX = 100000

      EPS   =  X02AJF()
      IFAIL = 1
      TOL   = 0.0e0

      PI = DACOS(-1.D0)

      HB2  = HB**2

      A=40
      R02 =  0.
      FREQ = 41.D0*DFLOAT(A)**(-1.D0/3.D0)

      K = MASS*FREQ**2.D0/HB2

      I_centre = Noeud/2

      t0 = -20.0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   This gives a realistic NN cross section around 40 mb
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      A = 2


      
      FORCE  = pas**2.D0*MASS*FREQ**2.D0/2.D0/HB2
      FORCE  = FORCE/(((pas*DFLOAT(Noeud-I_centre))**2.D0)/2.D0-1)/2.D0

      FORCE = 4*FORCE

      FORCE1 = FORCE*2.D0

      FORCE  = FORCE*pas**2.D0



c
c   In MeV/C^2
c
      DO J = 0,Noeud
          POT(J) = FORCE*DFLOAT((J-I_centre)**4) 
     &             +FORCE1*DFLOAT((J-I_centre)**2)

      ENDDO


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DETERMINATION OF EIGENSTATES AND EIGENVALUES OF RHO^1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      OPEN(UNIT = 20,file='DENSITY_ONEbass.DAT',STATUS='UNKNOWN')
      Y = 0.D0
      DO I = 1,Nmax
         DO J = 1,Nmax
            read(20,*) X 
            IF(I .EQ. J) Y = Y + DBLE(X)  
            DENS(I,J) = DBLE(X)  
         ENDDO
      ENDDO
      DO I=1,Nmax
         DO J = 1,(I-1)
            DENS(I,J) = DENS(J,I)
         ENDDO
      ENDDO
      CLOSE(20)
c     print *,'Y=',Y*pas


       call F01AJF(Nmax,TOL,DENS,Nmax,D,E,V,Nmax)
       call F02AMF(Nmax,EPS,D,E,V,Nmax,IFAIL)

       Norbital = 0
       DO I=1,Nmax
          D(I) = D(I)*pas
          IF(D(I) .GT. 0.001) Norbital = Norbital+1
       ENDDO

       print *,'Norbital=',Norbital
       IF( Norbital. LE. Norbitalmax) THEN

ccccccccccccccccccccccccccccccccccccccccc
c     OCCUPATION NUMBERS
ccccccccccccccccccccccccccccccccccccccccc
          
          Y = 0.
          DO I = 1,Norbital
             NI(I) = D(Nmax-I+1)
             Y = Y +NI(I)
             print *,'NI(',I,')=',NI(I)
          ENDDO
c          print *,'Y=',Y
          DO I = 1,Norbital
             NI(I) = NI(I)/Y
          ENDDO

ccccccccccccccccccccccccccccccccccccccccc
c     SINGLE PARTICLE WAVE-PACKETS
ccccccccccccccccccccccccccccccccccccccccc


          DO I = 1,Norbital
             DO J = 1,Nmax
                WAVE(J-1,I) = DCMPLX(V(J,Nmax-I+1)/dsqrt(pas),0.D0)     
             ENDDO
          ENDDO
      ELSE
         PRINT *,'TO MUCH WAVES...'
      ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     NORMALISATION OF WAVE-PACKETS IS DONE IN NAG
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   INITIALISATION OF VARIABLES FOR EVOLUTION 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      delta_p = HB*2*PI/(pas*DFLOAT(Nmax))


      Emax = (HB*PI/pas)**2.D0
     &        /(2.D0*MASS)
 
      delta_t = HB*2.D0*PI/(Emax)


      print *,'delta_t=',delta_t

c      delta_t = 0.01

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   TEST RECALCULATION OF THE ONE-BODY DENSITY MATRICE 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      OPEN(UNIT= 23,file ='dens_rec.dat',status='unknown')
c      DO J = 0,Noeud
c         Y = 0.D0
c         DO I = 1,Norbital
c            Y = Y + NI(I)*CDABS(WAVE(J,I))**2.D0
c         ENDDO
c         WRITE(23,*) DFLOAT(J)*pas,Y
c      ENDDO
      end


