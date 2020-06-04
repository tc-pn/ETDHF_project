      subroutine INITIALISATION(noeud,wave,pas,pot,K,t0,delta_t
     &                          ,CIN_COS,CIN_SIN,I_TEMPSMAX
     &                          ,POT_COS,POT_SIN,FORCE)
                                
      IMPLICIT NONE

      integer     noeud,A,I_TEMPSMAX
      COMPLEX*16  WAVE(0:Noeud,0:Noeud)     
      REAL*8      POT(0:Noeud,0:Noeud)     
      real*8      pas,K,t0,FORCE,HB,MASS,PI
      real*8      delta_t,Emax,delta_p
      real*8      FREQ,ALPHA,P1,P0,NORM,HB2,R02  
      real*8      CIN_COS(0:Noeud,0:Noeud)
      real*8      CIN_SIN(0:Noeud,0:Noeud)
      real*8      POT_COS(0:Noeud,0:Noeud)
      real*8      POT_SIN(0:Noeud,0:Noeud)
      real*8      SIGMA,FORCE1

      INTEGER     Q0,Q1,I,J,I_centre

      DATA HB/197.32705D0/,MASS/938.91897D0/
c
c   NB OF ITERATION WITH TIME EVOLUTION
c
      I_TEMPSMAX = 9000
c      I_TEMPSMAX = 20

      PI = DACOS(-1.D0)
c
c    En MeV  FREQ = HB*OMEGA
c

      HB2  = HB**2

      A=40
c      R_02 = (1.2D0*DFLOAT(A)**(1.D0/3.D0))**2.D0
      R02 =  0.
      FREQ = 41.D0*DFLOAT(A)**(-1.D0/3.D0)
c      K=1.01*DFLOAT(A)**(1.D0/6.D0)

       K= MASS * FREQ**2.D0/HB2

      I_centre = Noeud/2

c
c  TWO BODY INTERACTION
c

      SIGMA = 2.

cccccccccccccccccccccc
c ATTENTION TEST     c
cc@$!)!^&$)&@%)ccccccc


         t0 = -20.D0

c        t0 = -300. 
c        t0 = 0. 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   This gives a realistic NN cross section around 40 mb
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c   TWO coherent states centered at P0,Q0,P1.Q1
c


      ALPHA = -pas**2.D0*MASS*FREQ/(2.D0*HB2)

      NORM  = 1.D0/DSQRT(2.D0)*(MASS*FREQ/(HB2*PI))**0.5D0

      FORCE  = pas**2.D0*MASS*FREQ**2.D0/2.D0/HB2
      FORCE  = FORCE/(((pas*DFLOAT(Noeud-I_centre))**2.D0)/2.D0-1)/2.D0

      FORCE = 4*FORCE

      FORCE1 = FORCE*2.D0

      FORCE  = FORCE*pas**2.D0

c
c  In unit of pas
c
      Q0 = +8.
      Q1 = +8.
c      Q0 = 0.
c      Q1 = 0.
c
c   In MeV/C^2
c
      P0 = 0.
      P1 = 0.

      DO J = 0,Noeud
         DO I = 0,Noeud

c                POT(I,J) = FORCE*DFLOAT((I-I_centre)**2
c     &                                  +(J-I_centre)**2)
                POT(I,J) = FORCE*DFLOAT((I-I_centre)**4
     &                                  +(J-I_centre)**4)
     &                     +FORCE1*DFLOAT((I-I_centre)**2
     &                                  +(J-I_centre)**2)
                           
               POT(I,J) = POT(I,J)+ t0/(DSQRT(2*PI)*SIGMA) 
     &         *DEXP(-pas**2.D0*DFLOAT(I-J)**2.D0/(2*sigma**2.D0))

c             POT(I,J) = POT(I,J)+ t0/2 *pas**2.D0*DFLOAT(I-J)**2.D0

c            IF( I.eq.J) THEN
c                POT(I,J) = POT(I,J) +t0/pas 
c            ENDIF
cccccccccccccccccccccccccccccccccccccccc
c    FIRST COHERENT STATE
cccccccccccccccccccccccccccccccccccccccc
            wave(I,J) = NORM*DCMPLX(DCOS(P0*pas*DFLOAT(I-I_centre)/HB
     &                              -DFLOAT(Q0)*pas*P0/(2.D0*HB))
     &                        ,DSIN(P0*pas*DFLOAT(I-I_centre)/HB
     &                              -DFLOAT(Q0)*pas*P0/(2.D0*HB))
     &                        )*DEXP(ALPHA*
     &                          DFLOAT(I-I_centre-Q0)**2.D0)
cccccccccccccccccccccccccccccccccccccccc
c    SECOND COHERENT STATE
cccccccccccccccccccccccccccccccccccccccc
     &                         *DCMPLX(DCOS(P1*pas*DFLOAT(J-I_centre)/HB
     &                                  -DFLOAT(Q1)*pas*P1/(2.D0*HB))
     &                         ,DSIN(P1*pas*DFLOAT(J-I_centre)/HB
     &                              -DFLOAT(Q1)*pas*P1/(2.D0*HB))
     &                        )*DEXP(ALPHA*
     &                         DFLOAT(J-I_centre-Q1)**2.D0)     
ccccccccccccccccccccccccccccccccccccccccc
c    SYMETRIZATION
ccccccccccccccccccccccccccccccccccccccccc 
c     &                        +NORM*DCMPLX(DCOS(P0*pas
c     &                              *DFLOAT(J-I_centre)/HB
c     &                              -pas*DFLOAT(Q0)*P0/(2.D0*HB))
c     &                        ,DSIN(P0*pas*DFLOAT(J-I_centre)/HB
c     &                              -DFLOAT(Q0)*pas*P0/(2.D0*HB))
c     &                        )*DEXP(ALPHA*
c     &                          DFLOAT(J-I_centre-Q0)**2.D0)
ccccccccccccccccccccccccccccccccccccccccc
cc    SECOND COHERENT STATE
ccccccccccccccccccccccccccccccccccccccccc
c     &                         *DCMPLX(DCOS(P1*DFLOAT(I-I_centre)/HB
c     &                                  -DFLOAT(Q1)*P1/(2.D0*HB))
cc    &                         ,DSIN(P1*DFLOAT(I-I_centre)/HB
c    &                              -DFLOAT(Q1)*P1/(2.D0*HB))
c    &                        )*DEXP(ALPHA*
c    &                         DFLOAT(I-I_centre-Q1)**2.D0)     
c                                     
         ENDDO
      ENDDO

cccccccccccccccccccccccccccccccccccccccccccccccccccc
c   INITIALISATION OF VARIABLES FOR EVOLUTION 
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      delta_p = HB*2*PI/(pas*DFLOAT(Noeud+1))


      Emax = (HB*PI/pas)**2.D0
     &        /(2.D0*MASS)
             
      delta_t = HB*2.D0*PI/(Emax)



      call COEF_CIN(Noeud,pas,delta_t
     &                ,CIN_COS,CIN_SIN)

      call COEF_POT(Noeud,delta_t
     &              ,POT,POT_COS,POT_SIN)   
  


      end

      subroutine normalise(Noeud,wave,pas)

      IMPLICIT NONE

      INTEGER Noeud
      COMPLEX*16 wave(0:Noeud,0:Noeud)
      real*8 pas,NORM
      INTEGER I,J

      NORM = 0.D0
      DO I = 0,Noeud
         DO J = 0,Noeud
            NORM = NORM + CDABS(WAVE(J,I))**2.D0
         ENDDO
      ENDDO
      NORM = 1/(sqrt(NORM)*pas)
      DO I = 0,Noeud
         DO J = 0,Noeud
           WAVE(J,I) = WAVE(J,I)*NORM
         ENDDO
      ENDDO

      end

      subroutine calcul_densite(pas,noeud,wave,densite)     

      IMPLICIT NONE

      integer noeud
      real*8  pas
      complex*16 WAVE(0:Noeud,0:Noeud)     
      real*8 densite(0:Noeud)
      integer I,J

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c BE CAREFULL THE WAVE-PACKET HAVE SPINS
c
c we have rho(up,up)   = rho(down,down)=densite
c         rho(up,down) = rho(down,up)  =-densite
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DO I = 0,noeud
         densite(I) = 0.D0 
         DO J = 0,noeud
            DENSITE(I) = DENSITE(I) + CDABS(wave(I,J))**2.D0
         ENDDO
         densite(I) = DENSITE(I)*pas	
      ENDDO
      end

      subroutine calcul_dens(pas,noeud,wave,dens)     

      IMPLICIT NONE

      integer noeud
      real*8  pas
      complex*16 WAVE(0:Noeud,0:Noeud)     
      complex*16 dens(0:Noeud,0:Noeud)
      integer I,J,K

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c BE CAREFULL THE WAVE-PACKET HAVE SPINS
c
c we have rho(up,up)   = rho(down,down)=densite
c         rho(up,down) = rho(down,up)  =-densite
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DO I = 0,noeud
         DO J = 0,noeud
            dens(J,I) = DCMPLX(0.D0,0.D0)
            DO K = 0,Noeud
               DENS(J,I) = DENS(J,I) + DCONJG(wave(J,K))*
     &                                 wave(I,K)
            ENDDO
            DENS(J,I) = DENS(J,I)*pas
         ENDDO
      ENDDO
      end

      subroutine COEF_CIN(network_size,delta_x,delta_t
     &                ,CIN_COS,CIN_SIN)    

      IMPLICIT NONE
      
      INTEGER network_size,I,J,I_TEMP,J_TEMP
      REAL*8  delta_x,delta_t 
      real*8 CIN_COS(0:network_size,0:network_size)
      real*8 CIN_SIN(0:network_size,0:network_size)

      REAL*8 PI,HB,MASS,NORM,VALUE

      DATA HB/197.32705D0/,MASS/938.91897D0/


      PI = DACOS(-1.D0)
       
      NORM = HB*PI**2*delta_t/((network_size+1)**2*delta_x**2*MASS)

      DO J = 0,network_size
         DO I = 0,network_size
c
c   essai
c
            IF(I .LE. ((network_size+1)/2)) then
               I_TEMP = I
            ELSE
               I_TEMP = (network_size+1)-I
            ENDIF
            IF(J .LE. ((network_size+1)/2)) then
               J_TEMP = J
            ELSE
               J_TEMP = (network_size+1)-J
            ENDIF
            VALUE  = NORM*(DFLOAT((I_TEMP)**2)
     &                     +DFLOAT((J_TEMP)**2))

               CIN_COS(I,J) = DCOS(VALUE)
               CIN_SIN(I,J) = DSIN(VALUE)
         ENDDO
      ENDDO

      end

      subroutine COEF_POT(network_size,delta_t
     &                ,POT,POT_COS,POT_SIN)     

      IMPLICIT NONE
      
      INTEGER network_size,I,J

      real*8 delta_t,HB,VALUE,NORM
      real*8 POT(0:network_size,0:network_size)
      real*8 POT_COS(0:network_size,0:network_size)
      real*8 POT_SIN(0:network_size,0:network_size)

      DATA HB/197.32705D0/

      NORM = delta_t/HB

      DO J = 0,network_size
         DO I = 0,network_size
            VALUE  =  NORM*POT(I,J)             
            POT_COS(I,J) = DCOS(VALUE)
            POT_SIN(I,J) = DSIN(VALUE)

         ENDDO
      ENDDO

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     EVOLUTION BY SPLITTING METHOD
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C
C    2-Dimensionnal SPLIT-OPERATOR Calculation: EVOLUTION
C    Use only real and imaginary densite of wave function
c

      subroutine SPLIT_OPERATOR(network_size,WF,WF1
     &                         ,REAL_PART,IMAG_PART
     &                         ,CIN_COS,CIN_SIN
     &                         ,POT_COS,POT_SIN)

      IMPLICIT NONE

      integer network_size

      real*8 REAL_PART(0:network_size,0:network_size)
      real*8 IMAG_PART(0:network_size,0:network_size)
      real*8 CIN_COS(0:network_size,0:network_size)
      real*8 CIN_SIN(0:network_size,0:network_size)
      real*8 POT_COS(0:network_size,0:network_size)
      real*8 POT_SIN(0:network_size,0:network_size)
      complex*16 WF(0:network_size,0:network_size)
      complex*16 WF1(0:network_size,0:network_size)

      INTEGER   N_SIZE,I,J

      N_SIZE = (network_size+1)**2

c      call COEF_POT(network_size,delta_t
c     &                ,POT,POT_COS,POT_SIN)   


       DO J = 0,network_size
            DO I = 0,network_size
               REAL_PART(I,J) = DREAL(WF(I,J))
               IMAG_PART(I,J) = DIMAG(WF(I,J))
            ENDDO
         ENDDO

         call FFT_CIN(network_size,REAL_PART,IMAG_PART 
     &                   ,N_SIZE,CIN_COS,CIN_SIN) 
        
         call POT_TERM(network_size,REAL_PART,IMAG_PART,
     &                 POT_COS,POT_SIN)  

         call FFT_CIN(network_size,REAL_PART,IMAG_PART 
     &                   ,N_SIZE,CIN_COS,CIN_SIN) 

         DO J = 0,network_size
            DO I = 0,network_size
               WF1(I,J) = DCMPLX(REAL_PART(I,J)
     &                           ,IMAG_PART(I,J))
            ENDDO
         ENDDO

      end

      subroutine FFT_CIN(network_size,REAL_PART,IMAG_PART
     &                   ,N_SIZE,CIN_COS,CIN_SIN) 
                 

      IMPLICIT NONE

      INTEGER network_size,N_SIZE
      REAL*8 REAL_PART(0:network_size,0:network_size)
      REAL*8 IMAG_PART(0:network_size,0:network_size)
      real*8 CIN_COS(0:network_size,0:network_size)
      real*8 CIN_SIN(0:network_size,0:network_size)

      INTEGER I,J

      INTEGER NDIM,DIM_MAX,LWORK
      PARAMETER (NDIM = 2,DIM_MAX = 260, LWORK=3*DIM_MAX)
      INTEGER ND(NDIM),IFAIL
      REAL*8  WORK(LWORK),X_TEMP

      EXTERNAL C06FJF,C06GCF

      N_size=1
      DO I = 1,NDIM
         ND(I) = network_size+1
         N_size=N_size*ND(I)
      ENDDO

c
c     FFT OF Wave Functions
c
      call C06FJF(NDIM,ND,N_SIZE,REAL_PART
     &            ,IMAG_PART,WORK,LWORK,IFAIL)

c
c    Multiplication by the kinetic factor
c

      DO J = 0,network_size
        DO I = 0,network_size

           X_TEMP = REAL_PART(I,J)  

           REAL_PART(I,J) = CIN_COS(I,J)*REAL_PART(I,J)  
     &                      +CIN_SIN(I,J)*IMAG_PART(I,J)
           IMAG_PART(I,J) = CIN_COS(I,J)*IMAG_PART(I,J)
     &                       -CIN_SIN(I,J)*X_TEMP

         ENDDO
       ENDDO

c
c  FFT-1 OF THE WAVE FUNCTION
c

      call C06GCF(IMAG_PART,N_SIZE,IFAIL)
      call C06FJF(NDIM,ND,N_SIZE,REAL_PART
     &            ,IMAG_PART,WORK,LWORK,IFAIL)
      call C06GCF(IMAG_PART,N_SIZE,IFAIL)

      end


      subroutine POT_TERM(network_size,REAL_PART,IMAG_PART,
     &                    POT_COS,POT_SIN)  

      IMPLICIT NONE

      INTEGER network_size,I,J
      real*8 REAL_PART(0:network_size,0:network_size)
      real*8 IMAG_PART(0:network_size,0:network_size)
      real*8 POT_COS(0:network_size,0:network_size)
      real*8 POT_SIN(0:network_size,0:network_size)
      real*8 X_TEMP

      DO J = 0,network_size
         DO I = 0,network_size

            X_TEMP = REAL_PART(I,J)  


            REAL_PART(I,J) = POT_COS(I,J)*REAL_PART(I,J)  
     &                       +POT_SIN(I,J)*IMAG_PART(I,J)  


            IMAG_PART(I,J) = POT_COS(I,J)*IMAG_PART(I,J)
     &                       -POT_SIN(I,J)*X_TEMP

         ENDDO
      ENDDO

      end


      subroutine POSITION(noeud,pas,densite,WAVE,
     &                   Xmoy,Pmoy,SIGMAX,SIGMAP,SIGMAA) 

      IMPLICIT NONE

      INTEGER Noeud
      COMPLEX*16  WAVE(0:Noeud,0:Noeud)     
      real*8 pas,Xmoy,Pmoy,SIGMAX,SIGMAP,SIGMAA     
      real*8 DENSITE(0:Noeud)
      INTEGER J

      Xmoy   = 0.D0
      Pmoy   = 0.D0
      SIGMAX = 0.D0
      SIGMAP = 0.D0
      SIGMAA = 0.D0
    
      DO J = 0,Noeud 
         Xmoy = Xmoy + DFLOAT(J)*DENSITE(J)
         SIGMAX = SIGMAX + (pas*DFLOAT(J))**2.D0*DENSITE(J)  
      ENDDO


      Xmoy = Xmoy*pas**2.D0
      SIGMAX = SIGMAX*pas
c
c  HERE SIGMAX IS X^2
c 

      end


      subroutine position_rel(pas,noeud,wave,POS)

      IMPLICIT NONE 

      INTEGER Noeud
      real*8 pas,POS
      complex*16 wave(0:Noeud,0:Noeud)
      INTEGER I,J

      POS = 0.D0
      DO I= 0,Noeud
         DO J= 0,Noeud
            IF( ((i+j) .LE. Noeud) .AND. (i-j) .GE. 0) THEN 
                POS = POS + CDABS(WAVE(I+J,I-J))**2.D0*J
            ENDIF
         ENDDO
      ENDDO

      POS = POS*pas**3.D0

      end


      subroutine ENERGY(noeud,Wave,Etot,pas,POT)   

      IMPLICIT NONE

      INTEGER    Noeud
      real*8     pas,Etot,POT(0:Noeud,0:Noeud)
      complex*16 wave(0:Noeud,0:Noeud)
      real*8     HB,HBd 
      INTEGER    I,J

      Etot = 0.D0
      HB  = -20.7335D0/pas**2.D0
      HBd = -HB*2.D0


      DO J = 0,Noeud
         Etot = Etot + HBd*CDABS(WAVE(0,J))**2.D0 
     &               + HB*DCONJG(Wave(0,J))*Wave(1,J)    
         DO I=1,Noeud-1
            Etot = Etot + HBd*CDABS(WAVE(I,J))**2.D0 +HB
     &                    *DCONJG(Wave(I,J))*(Wave(I-1,J)
     &                    +Wave(I+1,J))
         ENDDO
         Etot = Etot + HBd*CDABS(WAVE(Noeud,J))**2.D0 
     &                +HB*DCONJG(Wave(Noeud,J))*Wave(Noeud-1,J)
      ENDDO

      DO J = 0,Noeud
         Etot = Etot + HBd*CDABS(WAVE(J,0))**2.D0 
     &               + HB*DCONJG(Wave(J,0))*Wave(J,1)    
         DO I=1,Noeud-1
            Etot = Etot + HBd*CDABS(WAVE(J,I))**2.D0 +HB
     &                    *DCONJG(Wave(J,I))*(Wave(J,I-1)
     &                    +Wave(J,I+1))
         ENDDO
         Etot = Etot + HBd*CDABS(WAVE(J,Noeud))**2.D0 
     &                +HB*DCONJG(Wave(J,Noeud))*Wave(J,Noeud-1)
      ENDDO



cccccccccccccccccccccccccccccccccccccccccccccccccc
c   POTENTIAL PART
cccccccccccccccccccccccccccccccccccccccccccccccccc

      DO J = 0,Noeud
         DO I = 0,Noeud
            Etot = Etot + CDABS(WAVE(I,J))**2.D0*POT(I,J)
         ENDDO
      ENDDO

      Etot = Etot*pas**2.D0
c      print *,'E=',Etot
      end

c      subroutine DIAG_RO(pas,DENS,DD,Entropy)  
      subroutine DIAG_RO(pas,DENS,DD)

      IMPLICIT NONE

      INTEGER    noeud
      PARAMETER (Noeud=2**6-1)
      real*8     pas
      complex*16 DENS(0:Noeud,0:Noeud)
      complex*16 one_minus_DENS(0:Noeud,0:Noeud)
      complex*16 logDENS(0:Noeud,0:Noeud)
      complex*16 log_one_minus_DENS(0:Noeud,0:Noeud)
      complex*16 Entropy_intermediary(Noeud+1,Noeud+1)
      real*8     Entropy

c Variables for diag. the density

      real*8 AR(Noeud+1,Noeud+1),AI(Noeud+1,Noeud+1),DD(Noeud+1)
      real*8 WK1(Noeud+1),WK2(Noeud+1),WK3(Noeud+1)
      real*8 VR(Noeud+1,Noeud+1),VI(Noeud+1,Noeud+1)

c Variables for diag. 1. - density

      real*8 AR2(Noeud+1,Noeud+1),AI2(Noeud+1,Noeud+1),DD2(Noeud+1)
      real*8 WK12(Noeud+1),WK22(Noeud+1),WK32(Noeud+1)
      real*8 VR2(Noeud+1,Noeud+1),VI2(Noeud+1,Noeud+1)

      integer ifail,I,J,K      

      EXTERNAL F02AXF

c        open(unit=14,file='test.dat',status='unknown')
c        dO J=0,Noeud
c           write(14,*) pas*dfloat(J)
c     &                ,DREAL(DENS(J,J)) 
c        ENDDO
c        close(14)

c Calculation of 1 - density

c        do i = 0, Noeud
c            do j = 0, Noeud
c                one_minus_DENS(i,j) = - DENS(i,j)
c                if (i.eq.j) then
c                    one_minus_DENS(i,j) = 1. - DENS(i,j) 
c                endif
c            enddo
c        enddo

c Real and imaginary parts of the density

      do k = 0,Noeud
         DD(k+1) = 0.D0
         do j = 0,Noeud
            AI(j+1,k+1) = DIMAG(DENS(j,k))
            AR(j+1,k+1) = DREAL(DENS(j,k))
         enddo
         if(AI(k+1,k+1) .ne. 0) then
c            print*,'AI(',k+1,k+1,')=',AI(k+1,k+1)
            AI(k+1,k+1) = 0.0
         endif  
      enddo

c Real and imaginary parts of one minus the density

c      do k = 0,Noeud
c         DD2(k+1) = 0.D0
c         do j = 0,Noeud
c            AI2(j+1,k+1) = DIMAG(one_minus_DENS(j,k))
c            AR2(j+1,k+1) = DREAL(one_minus_DENS(j,k))
c         enddo
c         if(AI2(k+1,k+1) .ne. 0) then
c            print*,'AI(',k+1,k+1,')=',AI(k+1,k+1)
c            AI2(k+1,k+1) = 0.0
c         endif  
c      enddo

c Diagonalization of the density 

      call F02AXF(AR,Noeud+1,AI,Noeud+1,Noeud+1
     &            ,DD,VR,Noeud+1,
     &            VI,Noeud+1,WK1,WK2,WK3,IFAIL) 

      if(IFAIL.NE.0) then
         print *,'Error'       
      endif
      DO I=1,Noeud+1
         DD(I) = DD(I)*pas
      ENDDO

c Diagonalization of one minus the density 

c      call F02AXF(AR2,Noeud+1,AI2,Noeud+1,Noeud+1
c     &            ,DD2,VR2,Noeud+1,
c     &            VI2,Noeud+1,WK12,WK22,WK32,IFAIL) 

c      if(IFAIL.NE.0) then
c         print *,'Error'       
c      endif
c      DO I=1,Noeud+1
c         DD2(I) = DD2(I)*pas
c      ENDDO

c Calculation of the log of density and 

c      call logM(AR,AI,DD,VR,
c     &            VI,WK1,WK2,WK3,IFAIL,logDENS)

c      call logM(AR2,AI2,DD2,VR2,
c     &            VI2,WK12,WK22,WK32,IFAIL,log_one_minus_DENS)

c        Entropy_intermediary = - MATMUL(DENS, logDENS) - MATMUL(one_minus_DENS,log_one_minus_DENS) 

c        Entropy = 0.
c        do i = 1, Noeud
c            Entropy = Entropy + DBLE(Entropy_intermediary(i,i))
c        enddo

      end

c      subroutine logM(AR,AI,DD,VR,
c     &            VI,WK1,WK2,WK3,IFAIL,logDENS)

c        IMPLICIT NONE

c      INTEGER    noeud
c      parameter (noeud = 2**6 - 1)
c      real*8     pas
c      complex*16 DENS(0:Noeud,0:Noeud)
c      complex*16 logDENS(0:Noeud,0:Noeud)

c      real*8 AR(Noeud+1,Noeud+1),AI(Noeud+1,Noeud+1),DD(Noeud+1)
c      real*8 WK1(Noeud+1),WK2(Noeud+1),WK3(Noeud+1)
c      real*8 VR(Noeud+1,Noeud+1),VI(Noeud+1,Noeud+1)
c      integer ifail,I,J,K 

c      INTEGER n, LDA, LDVR, LDVL, LWORK, INFO
c      PARAMETER (n = noeud+1) 
c      INTEGER IPIV(n)
c      CHARACTER*1, JOBVL, JOBVR
c      COMPLEX*16 WORK(2 * n), invVR(n,n)
c      DOUBLE PRECISION RWORK(2 * n)

c      LDA = n
c      LDVL = n
c      LDVR = n
c      LWORK = 2 * n
c      INFO = 0 

c      JOBVL = 'N'
c      JOBVR = 'V'


c We already have the diagonalized matrix
c Calculation of the log of the density:
c        do i = 0, Noeud
c            do j = 0, Noeud
c            logDENS(i,j) = (0.,0.)

c            if (i .eq. j) then
c                logDENS(i,j) = log(DD(i))
c            endif
c            enddo
c        enddo

c Now we just have to compute the inverse of the matric de passage to get back into the original basis

c        invVR = VR

c      call zgetrf(n,n,invVR,LDA,IPIV,INFO)

c      call zgetri(n,invVR,LDA,IPIV,WORK,LWORK,INFO)

c      logDENS = MATMUL(MATMUL(VR,logDENS),invVR)

c      end
