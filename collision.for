ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      TO USE WITH TDHF-1D 
c     
c      INPUT: TDHF WAVE PACKETS, OCCUPATION NUMBERS
c
c      1) DETERMINE THE INSTANTANEOUS BASES (WITH PARTICLE STATES)
c
c      2) DETERMINE THE COLLISION TERM IN THIS BASE
c
c          -> TAKE THE INSTANTANEOUS BASES AND DO THE EVOLUTION WITH
c                  NEGATIVE TIME TO HAVE THE MEMORY EFFECT
c
c      3) REDIAGONALISE THE DENSITY
c
c      4) KEEP IMPORTANT OCCUPATION NUMBERS TO HAVE THE NEW TDHF BASE
c
c      OUTPUT: NEW TDHF STATES WITH ASSOCIATED OCCUPATIONS
c              NUMBERS
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    PROGRAM CALLED BY THE USUAL TDHF                               c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine COLLISION(Norbital,pas,delta_t
     &                     ,WAVE,NI,MF
     &                     ,CIN_COS,CIN_SIN
     &                     ,POT_COS,POT_SIN,GAU,T_MACRO,OVERLAP,it)
cccccccccccccccccccccccccccccccccccccccc
c  !!NOTE MF is the TOTAL FIELD!!
cccccccccccccccccccccccccccccccccccccccc                               
                          
      IMPLICIT NONE

      INTEGER     Norbitalmax,Norbital,Noeud,Ncoll_max
      INTEGER     Iwave,it
      PARAMETER   (Norbitalmax = 50)
      PARAMETER   (Noeud = 2**6-1)
ccccccccccccccccccccccccccccccccccccccccccccccccc
c   CONSTRAINT Ncoll_max*Norbital < Norbitalmax c
ccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccc
c   Ncoll_max DETERMINE THE NUMBER OF PARTICLE PER
c             HOLES n_p=Ncoll_max/N_hole-1 
ccccccccccccccccccccccccccccccccccccccccccccc

      real*8      pas,MF(0:Noeud),delta_t,T_MACRO
      real*8      CIN_COS(0:Noeud)
      real*8      CIN_SIN(0:Noeud)
      real*8      POT_COS(0:Noeud)
      real*8      POT_SIN(0:Noeud)
      real*8      ni(Norbitalmax) 
      real*8      OVERLAP
      real*8      GAU(0:Noeud,0:Noeud)
      real*8      C_EXP(0:Noeud),P_EXP(0:Noeud)
      COMPLEX*16  wave(0:Noeud,Norbitalmax)
      COMPLEX*16  wave_temp(0:Noeud,Norbitalmax) 

      real*8      ENIV(Norbitalmax)
      real*8      delta_tIMAG
      INTEGER     I,K

      real*8      NORM
      complex*16  NORM1

      character*2 nb

ccccccccccccccccccccc
c   INITIALISATION  c
ccccccccccccccccccccc

      delta_tIMAG   = delta_t*197./100


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   1)  INSTANTANEOUS BASE DETERMINATION
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 105  FORMAT(I2)
     
      call INSTANT(Norbitalmax,Ncoll_max 
     &            ,Norbital,delta_tIMAG
     &            ,pas,WAVE,wave_temp,MF,C_EXP,P_EXP
     &            ,Iwave,ENIV,ni,OVERLAP)
                        

      print *,'SORTIE DE INSTANT'

      DO I = Norbital+1,Iwave
         DO K = 0,Noeud
            WAVE(K,I) = wave_temp(K,I)
         ENDDO
      ENDDO

      print *,'Iwave=',Iwave

      call COLLISION_TERM(Norbital
     &                   ,Iwave,delta_t,wave_temp  
     &                   ,ni,CIN_COS,CIN_SIN,POT_COS
     &                   ,POT_SIN,GAU,T_MACRO,pas,wave,it,eniv)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   RETURN THE NEW NUMBER OF WAVES, TDHF WAVES 
C              AND OCCUPATIONS NUMBERS.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      print *,'SORTIE DE collision'

c      call DIAG_DENS(Norbital,Iwave,pas,wave,DENS,NI)

c      print *,'Norbital=',Norbital
             
c      DO I = 1,Iwave
c         write(nb,105) I
c         OPEN(UNIT=10,file='TEST'//nb//'.DAT',status='unknown')
c         DO K = 0,Noeud
c            write(10,*) pas*DFLOAT(K),Dreal(WAVE_temp(K,I))
c         ENDDO
c         CLOSE(10)
c      ENDDO

c      STOP


ctttttttttttttttttttteeeeeeeeeeeeeeeeessssssssssssssssssttttttttttttttt

c      DO I = 1,Iwave
c            NORM = 0.
c            DO K = 0,Noeud
c               NORM  = NORM + CDABS(WAVE_TEMP(K,I))**2.D0
c            ENDDO
c
c            IF((IDINT(1000*(NORM*pas)) .NE. 1000) .AND. 
c     &         (IDINT(1000*(NORM*pas)) .NE. 999))  
c     &         print *,'PB NORM HF',I,'=',NORM*pas
c            IF(I .GT. 2) THEN
c               NORM1 = DCMPLX(0.D0,0.D0)
c               DO K = 0,Noeud
c                 NORM1 = NORM1 + DCONJG(WAVE_TEMP(K,I))
c     &                           *WAVE_TEMP(K,I-1)
c               ENDDO
c               IF(IDINT(1000*CDABS(NORM1)) .NE. 0)  
c     &         print *,'PB ORTH HF',I,'=',CDABS(NORM1*pas)
c           endif
c      ENDDO
c
ctttttttttttttttttttteeeeeeeeeeeeeeeeessssssssssssssssssttttttttttttttt





      end

      subroutine INSTANT(Norbitalmax,Ncoll_max
     &                  ,Norbital,delta_tIMAG
     &                  ,pas,WAVE,wave_temp,MF,C_EXP,P_EXP
     &                  ,Iwave,ENIV,ni,OVERLAP)
                        
      IMPLICIT NONE

      INTEGER     Norbitalmax,Norbital,Noeud,Ncoll_max
      PARAMETER   (Noeud = 2**6-1)
      real*8      MF(0:Noeud),POT(0:Noeud)
      real*8      pas,NORM,delta_tIMAG
      real*8      ni(Norbitalmax)       
      COMPLEX*16  wave(0:Noeud,Norbitalmax)
      COMPLEX*16  wave_temp(0:Noeud,Norbitalmax) 
      INTEGER     Iwave,I,J,K
      real*8      C_EXP(0:Noeud),OVERLAP
      real*8      P_EXP(0:Noeud),NORM2
      complex*16  NORM1
      real*8      ENIV(Norbitalmax),X,Y,R2,FORCE
      real*8      E1,E2

      COMPLEX*16  WF(0:Noeud),WF1(0:Noeud)

      character*2 nb

c      real*8      G05CAF
c      EXTERNAL    G05CAF
c      EXTERNAL    G05CCF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     WE WILL PUT IN Iwave the counter of created states by
c     imaginary time method  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


 101  FORMAT(I2)

c      call G05CCF

      FORCE = 0.5D0*pas**2.D0*1.1D0*40**0.166667D0

      DO Iwave = 1,Norbital
         DO I = 0,Noeud
            wave_temp(I,Iwave) = wave(I,Iwave)
         ENDDO
      ENDDO

      Iwave = Norbital

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   USE OF SPLIT-OPERATOR METHOD IN IMAGINARY TIME        c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      call CIN_EXP(Noeud,pas,delta_timag,C_EXP) 

      DO I = 1,Norbital
         IF(NI(I) .GT. 0.01) THEN
            DO K = 0,Noeud
               WF(K) = Wave(K,I)
            ENDDO

cccccccccccccccccccccccccccccccccccc
c  Number of particles per holes
cccccccccccccccccccccccccccccccccccc
            Ncoll_max =Norbitalmax/Norbital 

            IF(Ncoll_max .GT. 2) Ncoll_max = 2
c            IF(Ncoll_max .GT. 3) Ncoll_max = 3

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cxxxx  ADDING STOCHASTICITY IN THE FIELD   xxxxxxxx
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

c           X = G05CAF(X)
c           Y = G05CAF(Y)
            X=rand()
            Y=rand()
           print *,'X=',X,'Y=',Y

           DO K = 0,Noeud
              R2 = FORCE*(DFLOAT(K))**2.D0
              POT(K) = Y*MF(K)+X*R2
           ENDDO
           call POT_EXP(Noeud,delta_timag,POT,P_EXP)  


            DO J = 1,(Ncoll_max-1)
               call SPLIT_IMAG(WF,WF1,C_EXP,P_EXP)   
               NORM = 0.
               DO K = 0,Noeud
                  NORM = NORM + CDABS(WF1(K))**2.D0  
               ENDDO

               NORM = 1/DSQRT(NORM*pas)

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c   FILLING OF THE ENTIRE BASE
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c
c     ORTHOGONALITY TEST
c
c     maybe NOT do IT 
ccccccccccccccccccccccccccccccccccccccccc
               NORM1 = DCMPLX(0.D0,0.D0)
               DO K = 0,Noeud
                  NORM1 = NORM1 + DCONJG(WF1(K))*WF(K)*NORM
               ENDDO
               NORM1 = NORM1*pas

               NORM2 = 0.
               DO K = 0,Noeud
                  WF1(K) = WF1(K)*NORM-NORM1*WF(K)
                  NORM2 = NORM2 + CDABS(WF1(K))**2.D0  
               ENDDO
               NORM2 = 1/DSQRT(NORM2*pas)

               Iwave = Iwave + 1

               DO K = 0,Noeud
                  WF(K) = WF1(K)*NORM2
                  WAVE_TEMP(K,Iwave) = WF(K) 
               ENDDO

            ENDDO
         ENDIF
      ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      SCHMIDT ORTHONORMALISATION OF THE BASE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c       print *,'coucou'

       call SCHMIDT(Norbitalmax,Iwave,pas,wave_temp)


c       STOP
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     MAYBE WE COULD TAKE THIS AS A BASE                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       call diag_ener(Norbital,Iwave,pas,wave_temp,MF,ENIV)
       
       call ENERGY_F(Norbitalmax,Iwave,pas,wave_temp,MF,ENIV)
    
       E1 = 10000.D0
       E2 = -10000.D0
ccccccccccccccccccccccccccccccccc
c  ROUGH SELECTION ON ENERGIES  c
ccccccccccccccccccccccccccccccccc

       DO I=1,Norbital
          IF(ENIV(I) .GE. E2) E2 = ENIV(I)
          IF(ENIV(I) .LE. E1) E1 = ENIV(I)
       ENDDO


cccccccccccccccccccccccccccccc
c CALCULATION OF THE OVERLAP c
cccccccccccccccccccccccccccccc

       call  diag_H_overlap(Iwave,pas,wave_temp
     &                      ,MF,E1,E2,OVERLAP)  

            

c       call CLAS_ENER(Noeud,Norbitalmax,
c     &                Norbital,Iwave,pas,wave_temp,ENIV,E1,E2)      

      end

      subroutine CLAS_ENER(Noeud,Norbitalmax,
     &                     Norbital,Iwave,pas,wave_temp,ENIV,E1,E2)      

      IMPLICIT    NONE

      INTEGER     Norbitalmax,Iwave,Noeud,Norbital
      REAL*8      pas,ENIV(Norbitalmax),E1,E2,Delta_E  
      COMPLEX*16  wave_temp(0:Noeud,Norbitalmax) 
      INTEGER     I,J,COUNT

      Delta_E = 2*E2-2*E1

      COUNT = 0
      DO I = Norbital+1,Iwave   
         IF(ENIV(I) .LE. (E1+Delta_E)) THEN
            COUNT = COUNT + 1 
            DO J = 0,Noeud
               wave_temp(J,Norbital+COUNT) = wave_temp(J,I)
            ENDDO
         ENDIF
      ENDDO
      Iwave = Norbital + COUNT
      end

      subroutine POT_EXP(Noeud,delta_timag,MF,P_EXP)  

      IMPLICIT NONE

      INTEGER   Noeud
      real*8    delta_timag
      real*8    MF(0:Noeud),P_EXP(0:Noeud)
      real*8    HB,NORM
      INTEGER   I

      DATA HB/197.32705D0/

      NORM = -delta_timag/HB

      DO I = 0,Noeud
            P_EXP(I) = DEXP(NORM*MF(I))
      ENDDO
      end

      subroutine CIN_EXP(Noeud,pas,delta_timag,C_EXP) 

      IMPLICIT NONE

      INTEGER   Noeud
      real*8    pas,delta_timag
      real*8    C_EXP(0:Noeud)
      
      REAL*8    PI,HB,MASS,NORM,VALUE
      INTEGER   I,I_TEMP

      DATA HB/197.32705D0/,MASS/938.91897D0/


      PI = DACOS(-1.D0)
       
      NORM = HB*PI**2*delta_timag/((Noeud+1)**2*pas**2*MASS)

      DO I = 0,Noeud
c
c   essai
c
         IF(I .LE. ((Noeud+1)/2)) then
            I_TEMP = I
         ELSE
            I_TEMP = (Noeud+1)-I
         ENDIF

         VALUE  = -NORM*DFLOAT((I_TEMP)**2)
         C_EXP(I) = DEXP(VALUE)
      ENDDO

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    IMAGINARY TIME SPLIT OPERATOR METHOD
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine SPLIT_IMAG(WF,WF1,C_EXP,P_EXP)                             

      IMPLICIT NONE

      INTEGER    Noeud
      parameter  (Noeud =  2**6-1)

      real*8     REAL_PART(0:Noeud)
      real*8     IMAG_PART(0:Noeud)
      real*8     C_EXP(0:Noeud)
      real*8     P_EXP(0:Noeud)
      complex*16 WF(0:Noeud)
      complex*16 WF1(0:Noeud)

      INTEGER   N_SIZE,I

      N_SIZE = (Noeud+1)
      DO I = 0,Noeud
         REAL_PART(I) = DREAL(WF(I))
         IMAG_PART(I) = DIMAG(WF(I))
      ENDDO

      call FFT_CIN_IMAG(Noeud,REAL_PART,IMAG_PART 
     &                  ,N_SIZE,C_EXP)  
        
      call POT_TERM_IMAG(Noeud,REAL_PART,IMAG_PART,
     &                   P_EXP)  

      call FFT_CIN_IMAG(Noeud,REAL_PART,IMAG_PART 
     &                  ,N_SIZE,C_EXP) 

      DO I = 0,Noeud
         WF1(I) = DCMPLX(REAL_PART(I)
     &                  ,IMAG_PART(I))
      ENDDO

      end

      subroutine FFT_CIN_IMAG(Noeud,REAL_PART,IMAG_PART 
     &                        ,N_SIZE,C_EXP)  

      IMPLICIT NONE

      INTEGER Noeud,N_SIZE
      REAL*8 REAL_PART(0:Noeud)
      REAL*8 IMAG_PART(0:Noeud)
      real*8 C_EXP(0:Noeud)

      INTEGER I

      INTEGER NDIM,DIM_MAX,LWORK
      PARAMETER (NDIM = 1,DIM_MAX = 520, LWORK=3*DIM_MAX)
      INTEGER ND(NDIM),IFAIL
      REAL*8  WORK(LWORK)


      EXTERNAL C06FJF,C06GCF

      N_size=1
      DO I = 1,NDIM
         ND(I) = Noeud+1
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

      DO I = 0,Noeud
         REAL_PART(I) = C_EXP(I)*REAL_PART(I)  
         IMAG_PART(I) = C_EXP(I)*IMAG_PART(I)
      ENDDO

c
c  FFT-1 OF THE WAVE FUNCTION
c

      call C06GCF(IMAG_PART,N_SIZE,IFAIL)
      call C06FJF(NDIM,ND,N_SIZE,REAL_PART
     &            ,IMAG_PART,WORK,LWORK,IFAIL)
      call C06GCF(IMAG_PART,N_SIZE,IFAIL)

      end


      subroutine POT_TERM_IMAG(Noeud,REAL_PART,IMAG_PART
     &              ,P_EXP) 

      IMPLICIT NONE

      INTEGER Noeud,I
      real*8 REAL_PART(0:Noeud)
      real*8 IMAG_PART(0:Noeud)
      real*8 P_EXP(0:Noeud)

      DO I = 0,Noeud
         REAL_PART(I) = P_EXP(I)*REAL_PART(I)  
         IMAG_PART(I) = P_EXP(I)*IMAG_PART(I)
      ENDDO

      end

cccccccccccccccccccccccccccccccccccccccccccccccccc
c   SCHMIDT FOR COMPLEX MATRICES
cccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine SCHMIDT(Norbitalmax,Iwave,pas,wave_temp)   

      IMPLICIT NONE 

      INTEGER     Noeud,Norbitalmax,Iwave,Nmax
      PARAMETER   (Noeud = 2**6-1)
      real*8      pas,NORM
      complex*16  NORM1
      COMPLEX*16  wave_temp(0:Noeud,Norbitalmax)
      COMPLEX*16  V(0:Noeud)
      INTEGER     I,J,K,ICOUNT

      ICOUNT = 0

      DO I = 1,Iwave
         IF (I .GT. 1) THEN
            NORM = 0.
            DO K = 0,Noeud
               NORM = NORM + CDABS(WAVE_TEMP(K,I))**2.D0      
            ENDDO
            NORM = 1/DSQRT(NORM*pas)
            DO K = 0,Noeud
              V(K) = WAVE_TEMP(K,I)*NORM
            ENDDO

cccccccccccccccccccccc
c   OVERLAP TEST     c
cccccccccccccccccccccc

            NORM1 = DCMPLX(0.D0,0.D0)
            DO J = 1,ICOUNT
               DO K = 0,Noeud
                  NORM1 = NORM1 + DCONJG(V(K))
     &                          *WAVE_TEMP(K,J)
               ENDDO
            ENDDO
            NORM = CDABS(NORM1 * pas)
c            print *,'OVERLAP=',I,NORM


            IF(NORM .LE. 0.95D0) THEN
c            IF(NORM .LE. 2.D0) THEN
ccccccccccccccccccccccc
c ORTHONORMALISATION  c
ccccccccccccccccccccccc
               DO J=1,ICOUNT
                  NORM1 = DCMPLX(0.D0,0.D0)
                  DO K = 0,Noeud
                     NORM1 = NORM1 + V(K)
     &                             *DCONJG(WAVE_TEMP(K,J))
                  ENDDO
                  NORM1 = NORM1*pas
                  DO K = 0,Noeud
                     V(K) = V(K) - NORM1 * WAVE_TEMP(K,J)
                  ENDDO
               ENDDO
               ICOUNT = ICOUNT + 1
               NORM = 0.
               DO K = 0,Noeud
                  WAVE_TEMP(K,ICOUNT) = V(K)
                  NORM  = NORM + CDABS(WAVE_TEMP(K,ICOUNT))**2.D0
               ENDDO
               NORM = 1/DSQRT(NORM*pas)
               DO K = 0,Noeud
                  WAVE_TEMP(K,ICOUNT) = wave_temp(K,ICOUNT)*NORM    
               ENDDO
ccccccccccccccccccccccccccccc
c END OF THE OVERLAP TEST   c
ccccccccccccccccccccccccccccc
            ENDIF
         ELSE
            ICOUNT = ICOUNT + 1
            NORM = 0.
            DO K = 0,Noeud
               NORM  = NORM + CDABS(WAVE_TEMP(K,I))**2.D0      
            ENDDO
            NORM = 1/DSQRT(NORM*pas)
            DO K = 0,Noeud
               WAVE_TEMP(K,I)   = WAVE_TEMP(K,I)*NORM  
            ENDDO
         ENDIF
      ENDDO 

      Iwave = ICOUNT

cccccccccccccccccccccccccccccccccccccc
c   TEST
cccccccccccccccccccccccccccccccccccccc
c      DO I = 1,Iwave
c            NORM = 0.
c            DO K = 0,Noeud
c               NORM  = NORM + CDABS(WAVE_TEMP(K,I))**2.D0
c            ENDDO
c            print *,'NORM',I,'=',NORM*pas
c            IF(I .GT. 2) THEN
c               NORM = 0.
c               DO K = 0,Noeud
c                 NORM = NORM + DCONJG(WAVE_TEMP(K,I))
c     &                         *WAVE_TEMP(K,I-1)
c               ENDDO
cc               print *,'ORTH',I,'=',NORM
c           endif
c      ENDDO
cc      OPEN(UNIT=10,file='TEST1.DAT',status='unknown')
c      DO K = 0,Noeud
c         write(10,*) pas*DFLOAT(K),CDABS(WAVE_TEMP(K,5))
c      ENDDO
c      CLOSE(10)
c      STOP
      end
ccccccccccccccccccccccccccccccccccccccccc
c CALCULATION OF THE ENERGY USING FFT   c
ccccccccccccccccccccccccccccccccccccccccc


      subroutine ENERGY_F(Norbitalmax,Iwave,pas,wave_temp,MF,ENIV)    
                 
      IMPLICIT NONE

      INTEGER     Noeud,Norbitalmax,Iwave

      PARAMETER   (Noeud= 2**6-1)

      real*8      MF(0:Noeud)
      real*8      ENIV(Norbitalmax),pas,EKIN,EPOT 
      COMPLEX*16  wave_temp(0:Noeud,Norbitalmax)
      real*8      REAL_PART(0:Noeud)
      real*8      IMAG_PART(0:Noeud)

      INTEGER     I_temp,I,N_SIZE,J
      REAL*8      PI,HB,MASS,NORM,VALUE(0:Noeud)

      DATA        HB/197.32705D0/,MASS/938.91897D0/

      PI = DACOS(-1.D0)
       
      NORM = 2*(HB*PI)**2/((Noeud+1)**2*pas**2*MASS)

cccccccccccccccccccccccccccccccccccccccc
c  CALCULATION OF P^2/2m               c
cccccccccccccccccccccccccccccccccccccccc
      DO I = 0,Noeud

         IF(I .LE. ((Noeud+1)/2)) then
            I_TEMP = I
         ELSE
            I_TEMP = (Noeud+1)-I
         ENDIF

         VALUE(I)  = NORM*DFLOAT((I_TEMP)**2)

      ENDDO

      N_SIZE = (Noeud+1)
      DO J=1,Iwave
         ENIV(J) =0
         DO I = 0,Noeud
            REAL_PART(I) = DREAL(WAVE_TEMP(I,J))
            IMAG_PART(I) = DIMAG(WAVE_TEMP(I,J))
         ENDDO
cccccccccccccc
c KINETIK PART
cccccccccccccc
         call EN_CIN(Noeud,REAL_PART,IMAG_PART,VALUE 
     &               ,N_SIZE,pas,EKIN)      
cccccccccccccccccc
c POTENTIAL PART c
cccccccccccccccccc
         EPOT = 0.D0
         DO I = 0,Noeud
            EPOT = EPOT + MF(I)*CDABS(WAVE_TEMP(I,J))**2.D0
         ENDDO

         EPOT  = EPOT*pas

         ENIV(J) = EKIN + EPOT
c         print *,'ENIV(',J,')=',ENIV(J)
      ENDDO 
      end

      subroutine EN_CIN(Noeud,REAL_PART,IMAG_PART,VALUE 
     &                 ,N_SIZE,pas,EKIN)  

      IMPLICIT  NONE

      INTEGER   Noeud,N_SIZE
      REAL*8    REAL_PART(0:Noeud)
      REAL*8    IMAG_PART(0:Noeud)
      REAL*8    VALUE(0:Noeud),EKIN,pas
      REAL*8    HB,PI,delta_p,Y

      INTEGER I

      INTEGER NDIM,DIM_MAX,LWORK
      PARAMETER (NDIM = 1,DIM_MAX = 520, LWORK=3*DIM_MAX)
      INTEGER ND(NDIM),IFAIL
      REAL*8  WORK(LWORK)


      EXTERNAL C06FJF,C06GCF

      DATA HB/197.32705D0/

      PI = DACOS(-1.D0)

      delta_p = HB*2*PI/(pas*DFLOAT(Noeud+1))


      N_size=1
      DO I = 1,NDIM
         ND(I) = Noeud+1
         N_size=N_size*ND(I)
      ENDDO
c
c     FFT OF Wave Functions
c
      call C06FJF(NDIM,ND,N_SIZE,REAL_PART
     &            ,IMAG_PART,WORK,LWORK,IFAIL)

c
c   Kinetic energy
c
      Y = 0.D0
      DO I = 0,Noeud
         Y = Y + REAL_PART(I)**2.D0
     &         + IMAG_PART(I)**2.D0   
      ENDDO
      Y = Y * delta_p

      EKIN = 0.D0

      DO I = 0,Noeud
         EKIN = EKIN + VALUE(I)*(REAL_PART(I)**2.D0
     &                           +IMAG_PART(I)**2.D0)   
      ENDDO

      EKIN = EKIN*delta_p/Y

      end

ccccccccccccccccccccccccccccccccccccccccc
c USUAL CALCULATION OF THE ENERGY       c
ccccccccccccccccccccccccccccccccccccccccc

      subroutine ENERGY(Norbitalmax,Iwave,pas,wave_temp,MF,ENIV)    

      IMPLICIT NONE

      INTEGER     Noeud,Norbitalmax,Iwave

      PARAMETER   (Noeud= 2**6-1)

      real*8      MF(0:Noeud),HB,H1
      real*8      ENIV(Norbitalmax),pas 
      COMPLEX*16  wave_temp(0:Noeud,Norbitalmax)
      complex*16  AX(0:noeud),BX(0:noeud),CX(0:noeud),DX(0:noeud)
      INTEGER     N,K

      DATA        HB/20.7335D0/

      H1 = HB/pas**2

      do k=0,noeud 
          AX(k)=DCMPLX(-H1,0.D0)
          BX(k)=DCMPLX(2*H1+MF(k),0.D0)
          CX(k)=DCMPLX(-H1,0.D0)
      enddo

      do n=1,Iwave
         ENIV(N) = 0.D0

         DX(0)=BX(0)*WAVE_TEMP(0,n)+CX(0)*WAVE_TEMP(1,n)
         do k=1,(noeud-1)
            DX(k)=AX(k)*WAVE_TEMP(k-1,n)+BX(k)*WAVE_TEMP(k,n)
     &      +CX(k)*WAVE_TEMP(k+1,n)
         enddo
         DX(noeud)=AX(noeud)*WAVE_TEMP(noeud-1,n)
     &            +BX(noeud)*WAVE_TEMP(noeud,n)

         do k=0,noeud
             ENIV(N) =ENIV(N)+DX(k)*DCONJG(WAVE_TEMP(k,n))
         enddo
         ENIV(N) = ENIV(N)*pas
c         print *,'ENIV(',N,')=',ENIV(N)
       enddo

      end 

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DIAGONALISATION OF THE ENERGY IN THE SUBSPACE
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine diag_ener(Norbital
     &                      ,Iwave,pas,wave_temp,MF,ENIV)
                           
      IMPLICIT NONE

      INTEGER     Noeud,Norbitalmax,Norbital,Iwave

      PARAMETER   (Noeud= 2**6-1)
      PARAMETER   (Norbitalmax = 50)

      real*8      MF(0:Noeud),pas,HB,H1
      real*8      ENIV(Norbitalmax),NORM 
      COMPLEX*16  wave_temp(0:Noeud,Norbitalmax)
      COMPLEX*16  w_temp(0:Noeud,Norbitalmax)
      complex*16  AX(0:noeud),BX(0:noeud),CX(0:noeud),DX(0:noeud)
      INTEGER     N,K,I,J
c      character*2 nb 
cccccccccccccccccccccccc
c   NAG VARIABLES      c
cccccccccccccccccccccccc
      real*8 AR(Noeud+1,Noeud+1),AI(Noeud+1,Noeud+1)
      real*8 R(Noeud+1),VR(Noeud+1,Noeud+1),VI(Noeud+1,Noeud+1) 
      real*8 WK1(Noeud+1),WK3(Noeud+1),WK2(Noeud+1)
      
      INTEGER IFAIL

      EXTERNAL F02AXF

      DATA   HB/20.7335D0/

c      print *,'Norbital =',Norbital
c      print *,'Iwave=',Iwave,'pas',pas

      H1 = HB/pas**2

      IFAIL = 1

      DO N = 1,Noeud+1
         DO K = 1,Noeud+1
            AR(K,N) = 0.D0
            AI(K,N) = 0.D0            
         ENDDO
      ENDDO
      do k=0,noeud 
          AX(k)=DCMPLX(-H1,0.D0)
          BX(k)=DCMPLX(2*H1+MF(k),0.D0)
          CX(k)=DCMPLX(-H1,0.D0)
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccc
c    DIAGONALISATION OF (1-P)H(1-P)           c
ccccccccccccccccccccccccccccccccccccccccccccccc

      DO I = Norbital+1,Iwave
         do n = Norbital+1,I
            DX(0)=BX(0)*WAVE_TEMP(0,n)+CX(0)*WAVE_TEMP(1,n)
            do k=1,(noeud-1)
               DX(k)=AX(k)*WAVE_TEMP(k-1,n)+BX(k)*WAVE_TEMP(k,n)
     &         +CX(k)*WAVE_TEMP(k+1,n)
            enddo
            DX(noeud)=AX(noeud)*WAVE_TEMP(noeud-1,n)
     &               +BX(noeud)*WAVE_TEMP(noeud,n)

            do k=0,noeud
               AR(I,N) = AR(I,N)+DREAL(DX(k)
     &                                 *DCONJG(WAVE_TEMP(k,I)))
               AI(I,N) = AI(I,N)+DIMAG(DX(k)
     &                                 *DCONJG(WAVE_TEMP(k,I)))
            enddo
            AR(I,N) = AR(I,N)*pas
            AI(I,N) = AI(I,N)*pas

            IF( I .EQ. N) AI(I,I) = 0.D0 

         enddo
      ENDDO


      call F02AXF(AR,Noeud+1,AI,Noeud+1,Iwave,R,VR,Noeud+1,VI
     &            ,Noeud+1,WK1,WK2,WK3,IFAIL)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   BE CAREFULL THE EIGENVALUES ARE IN ASCENDING ORDER
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      print *,'IFAIL=',IFAIL

      DO I = 1,Iwave
         DO K = 0,Noeud
            W_temp(K,I) = DCMPLX(0.D0,0.D0)
         ENDDO
      ENDDO


      DO I = 1,Iwave
         ENIV(I) = R(I)
c         print *,'R(',I,')=',R(I)
         IF ( R(I) .NE. 0) THEN
            DO K = 0,Noeud
               DO J = Norbital+1,Iwave
               W_temp(K,I) = W_temp(K,I)  
     &                       + DCMPLX(VR(J,I)*
     &                         DREAL(wave_temp(K,J))-
     &                         VI(J,I)*DIMAG(wave_temp(K,J)),
     &                         VR(J,I)*DIMAG(wave_temp(K,J))+
     &                         VI(J,I)*DREAL(wave_temp(K,J)))
               ENDDO
            ENDDO
         ENDIF
      ENDDO


      N = Norbital

      DO I = 1,Iwave
         IF ( R(I) .NE. 0) THEN
             N = N + 1
             NORM = 0. 
             DO K = 0,Noeud
                 wave_temp(K,N) = W_temp(K,I) 
                 NORM = NORM + CDABS(W_temp(K,I))**2.D0
             ENDDO
             NORM = 1/DSQRT(NORM*pas)
             DO K = 0,Noeud
                 wave_temp(K,N) = wave_temp(K,N)*NORM
             ENDDO
         ENDIF
      ENDDO

      end


                         
      subroutine COLLISION_TERM(Norbital
     &                         ,Iwave,delta_t,wave_temp  
     &                         ,ni,CIN_COS,CIN_SIN,POT_COS
     &                         ,POT_SIN,GAU,T_MACRO,pas,wave,it,eniv) 
                               
      IMPLICIT NONE

      INTEGER     Norbitalmax,Norbital,Iwave
      INTEGER     Noeud,ITIME                        

      PARAMETER   (Noeud        = 2**6-1)
      PARAMETER   (Norbitalmax  = 50)

c      PARAMETER   (ITIME = 600) 
c First test, tau = t_back
       PARAMETER   (ITIME = 240)

c Second test t_back = tau/2
c       PARAMETER   (ITIME = 120)

c third test t_back = tau/12
c        PARAMETER   (ITIME = 20)
     
c      PARAMETER   (ITIME = 300)      
c      PARAMETER   (ITIME = 91)      
c      PARAMETER   (ITIME = 45)      

ccccccccccccccccccccccccccccccccccccccccccccccc
c     TIME INTEGRATION PARAMETER ITIME        c
ccccccccccccccccccccccccccccccccccccccccccccccc

      real*8      ENIV(Norbitalmax)
      REAL*8      delta_t,T_MACRO,pas
      real*8      CIN_COS(0:Noeud)
      real*8      CIN_SIN(0:Noeud)
      real*8      POT_COS(0:Noeud)
      real*8      POT_SIN(0:Noeud)
      real*8      ni(Norbitalmax)
      real*8      GAU(0:Noeud,0:Noeud)
      COMPLEX*16  wave_temp(0:Noeud,Norbitalmax)
      COMPLEX*16  wave(0:Noeud,Norbitalmax)
      COMPLEX*16  DENS(Norbitalmax,Norbitalmax) 
      COMPLEX*16  F_coll(Norbitalmax,Norbitalmax)
      REAL*8      W_gain(Norbitalmax)
      REAL*8      W_loss(Norbitalmax)
      COMPLEX*16  V(Norbitalmax,Norbitalmax
     &             ,Norbitalmax,Norbitalmax)
      COMPLEX*16  V_0(Norbitalmax,Norbitalmax
     &             ,Norbitalmax,Norbitalmax) 
      COMPLEX*16  INTEGR(Norbitalmax,Norbitalmax
     &             ,Norbitalmax,Norbitalmax)

      INTEGER     TEMPSMAX,Norb   
      PARAMETER   (TEMPSMAX = 200, Norb = 3)   

       real*8      loss_term(0:TEMPSMAX,Norb)       
c      COMPLEX*16  gain_term(0:TEMPSMAX,Norb)
    


      INTEGER I,T0_OR_NOT,J,K,L,IT1,it
      INTEGER I_alpha,I_beta,I_delta
      real*8  CTE,CTE1
      real*8  HBAR,HB2

cccccccccccccccccccccccccccccccccccccc
c  AVERAGE RELAXATION TIME FROM TDHF
cccccccccccccccccccccccccccccccccccccc
      real*8  TAU

      DATA HBAR/197.625D0/

c      TAU  = 64.D0/4.D0
      TAU  = 64.D0
c      TAU  = 128.D0
c      TAU  = 64.D0/16.D0
c      TAU  = 64.D0/32.D0

      HB2  = HBAR**2.D0
      CTE  = -delta_t*T_MACRO*2/HB2 

      CTE1 = +delta_t*4/HB2
 
cccccccccccccccccccccccccccccccccccccc
c  BE CAREFULL TO THE FACTOR 4
cccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   OCCUPATION NUMBERS OF PARTICLE STATES              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO I = Norbital+1,Iwave
         NI(I) = 0.D0
      ENDDO
cccccccccccccccccccccccccccc
c    !! NEGATIVE TIME !!   c
cccccccccccccccccccccccccccc
      DO I = 0,Noeud
         CIN_SIN(I) = -CIN_SIN(I)
         POT_SIN(I) = -POT_SIN(I) 
      ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c   HIGH FREQUENCIES COMPONENT
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO I = 1,Norb
         DO J = 1,TEMPSMAX
            loss_term(J,I) = DCMPLX(0.D0,0.D0)
         ENDDO
      ENDDO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOOPS ON NEGATIVE TIME                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            if (it .eq. 10) then
            open(unit=160,file='ENIV_T.dat',status='unknown')
            endif
            write(160,*) it * 0.54, (eniv(i), i = 1,iwave)
            print*, it * 0.54, (eniv(i), i = 1,iwave)

     
      DO IT1 = 1,ITIME

c
c

c       print *,'I_neg=',IT1

         IF ( IT1 .EQ. 1) THEN
            T0_OR_NOT = 1
            call TWO_BODY_POT(Norbitalmax,Norbital,Iwave,Noeud
     &                       ,wave_temp,V_0,GAU,T0_OR_NOT,pas)
            DO L = 1,Iwave
               DO K = 1,Iwave
                  DO J =1,Iwave
                     DO I = 1,Iwave
                        IF((J .GT. Norbital .AND. K .GT. Norbital)
     &                                   .OR.
     &                  (J .GT. Norbital .AND. L .GT. Norbital)
     &                                   .OR.
     &                  (I .GT. Norbital .AND. K .GT. Norbital)
     &                                   .OR.
     &                  (I .GT. Norbital .AND. L .GT. Norbital))
     &                  GOTO 200
                        INTEGR(I,J,K,L) = V_0(I,J,K,L)/2.D0*DEXP(
     &                                   -DFLOAT(IT1)*delta_t/TAU)   
c                        IF( L .LE. Norb) THEN
c                            loss_term(IT1-1,L) = loss_term(IT1-1,L) +
c     &                                        DREAL(
c     &                                        V_0(I,J,L,K)/2.D0*CTE
c     &                                        *DEXP(-DFLOAT(IT1)
c     &                                              *delta_t/TAU)   
c     &                                        *V_0(L,K,I,J)  
c     &                                        *(NI(L)*NI(K)*(1-NI(I))
c     &                                        *(1-NI(J))
c     &                                        -(1-NI(L))*(1-NI(K))
c     &                                        *NI(I)*NI(J))
c     &                                             )
c
c                        ENDIF
  200                   CONTINUE
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO

         ELSE
            T0_OR_NOT = 0
            call TWO_BODY_POT(Norbitalmax,Norbital,Iwave,Noeud
     &                       ,wave_temp,V,GAU,T0_OR_NOT,pas)
            DO L = 1,Iwave
               DO K = 1,Iwave
                  DO J =1,Iwave
                     DO I = 1,Iwave
                        IF((J .GT. Norbital .AND. K .GT. Norbital)
     &                                   .OR.
     &                  (J .GT. Norbital .AND. L .GT. Norbital)
     &                                   .OR.
     &                  (I .GT. Norbital .AND. K .GT. Norbital)
     &                                   .OR.
     &                  (I .GT. Norbital .AND. L .GT. Norbital))
     &                  GOTO 240
                        INTEGR(I,J,K,L) = INTEGR(I,J,K,L)+
     &                                    V(I,J,K,L)
     &                                   *DEXP(
     &                                   -DFLOAT(IT1)*delta_t/TAU)
c                        IF( L .LE. Norb) THEN
c                            loss_term(IT1-1,L) = loss_term(IT1-1,L) +
c     &                                        DREAL(
c     &                                        V(I,J,L,K)*CTE
c     &                                        *DEXP(-DFLOAT(IT1)
c     &                                              *delta_t/TAU)   
c     &                                        *V_0(L,K,I,J)  
c     &                                        *(NI(L)*NI(K)*(1-NI(I))
c     &                                        *(1-NI(J))
c     &                                        -(1-NI(L))*(1-NI(K))
c     &                                        *NI(I)*NI(J))
c     &                                             )
c                        ENDIF
  240                   CONTINUE
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c  EVOLUTION IN NEGATIVE TIME x
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

         call EVOLUTION_2(Norbitalmax,WAVE_temp,wave_temp
     &                    ,CIN_COS,CIN_SIN
     &                    ,POT_COS,POT_SIN,Iwave)

cccccccccccccccccccccccccccccccccccc
c  END OF THE NEGATIVE TIME LOOPS  c
cccccccccccccccccccccccccccccccccccc

      ENDDO

cccccccccccccccccccccccccccccccccccccccccccccccccc
c   GAIN TERM IN THE DENSITY                                             
cccccccccccccccccccccccccccccccccccccccccccccccccc

      DO J = 1,Iwave
         W_gain(I) = 0.D0
         W_loss(I) = 0.D0
         DO I = 1,Iwave
            F_coll(I,J) = DCMPLX(0.D0,0.D0)
            DO I_alpha = 1,Iwave
               DO I_beta = 1,Iwave
                  DO I_delta = 1,Iwave
                     F_coll(I,J) = F_coll(I,J)
     &                   + V_0(I,I_delta,I_alpha,I_beta)
     &                   *INTEGR(I_alpha,I_beta,J,I_delta)  
     &                   *(NI(J)*NI(I_delta)*(1-NI(I_alpha))
     &                   *(1-NI(I_beta))
     &                   -(1-NI(J))*(1-NI(I_delta))
     &                   *NI(I_alpha)*NI(I_beta)
     &                   )
                     IF(I .EQ. J) THEN
                        W_gain(I) = W_gain(I)
     &                    + DREAL(V_0(I,I_delta,I_alpha,I_beta)
     &                    *INTEGR(I_alpha,I_beta,J,I_delta)  
     &                    *(1-NI(I_delta))
     &                    *NI(I_alpha)*NI(I_beta)
     &                             )
                        W_loss(I) = W_loss(I)
     &                    + DREAL(V_0(I,I_delta,I_alpha,I_beta)
     &                    *INTEGR(I_alpha,I_beta,J,I_delta)  
     &                    *NI(I_delta)*(1-NI(I_alpha))
     &                    *(1-NI(I_beta))
     &                    )
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   DENSITY MATRIX FILLING  dens(t+Dt)=dens(t)+ Dt... c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   DENS(I,J) = <i|delta \rho|j>
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DO I = 1,Iwave
         DO J = 1,Iwave
              DENS(J,I) = CTE*(F_coll(J,I)+DCONJG(F_coll(I,J)))
              IF(I .EQ. J) THEN
                W_loss(I) = CTE1*W_loss(I)
                W_gain(I) = CTE1*W_gain(I)
              ENDIF
         ENDDO
      ENDDO

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     FIRST ORDER PERTURBATION                        c
c     FOR OCCUPATIONS NUMBERS                         c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call PERT_one(Norbital,Iwave,DENS,W_gain
     &              ,W_loss,T_MACRO,pas,wave,NI)   
                   
cccccccccccccccccccccccccccccccccccccc
c    !! RETURN TO POSITIVE TIME !!   c
cccccccccccccccccccccccccccccccccccccc

      DO I = 0,Noeud
         CIN_SIN(I) = -CIN_SIN(I)
         POT_SIN(I) = -POT_SIN(I) 
      ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   TEST OF INTERVENTION OF HIGH FREQUENCIES         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      call FOURIER_TIME(Norb,TEMPSMAX,delta_t,loss_term)  
     

      end         

                   
      subroutine PERT_ONE(Norbital,Iwave,DENS,W_gain
     &                    ,W_loss,T_MACRO,pas,wave,NI)
                 
      IMPLICIT NONE

      INTEGER     Noeud,Norbitalmax,Norbital,Iwave

      PARAMETER   (Noeud= 2**6-1)
      PARAMETER   (Norbitalmax = 50)

      real*8      pas,NORM,Y,R(Norbitalmax)
      
      real*8      NI(Norbitalmax),R1(Norbitalmax) 
      COMPLEX*16  DENS(Norbitalmax,Norbitalmax) 
      COMPLEX*16  wave(0:Noeud,Norbitalmax)
      COMPLEX*16  w_temp(0:Noeud,Norbitalmax)
      REAL*8      W_gain(Norbitalmax),T_MACRO
      REAL*8      W_loss(Norbitalmax),TAU

      INTEGER     I,K,J,N

      character*2 nb

 106  FORMAT(I2)

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   NEW OCCUPATION NUMBERS/A KIND OF EXACT SOLUTION c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DO I=1,Iwave
         R1(I) = NI(I)+ DREAL(DENS(I,I))
c         print *,'R(',I,')=',R1(I)
      ENDDO

      DO I=1,Iwave
         TAU = W_gain(I)+W_loss(I)
c         print *,'TAU,I=',TAU,I
         IF(TAU .GT. 1.e-10) THEN 
            NORM = W_gain(I)/(W_gain(I)+W_loss(I))
c            print *,'EQU',I,'=',NORM
            NI(I) = (NI(I) - NORM)*DEXP(-T_MACRO*TAU) +NORM
         ELSE
            IF(NI(I) .LE. 0.D0 ) print *,'PROBLEM NI NEGATIVE'
         ENDIF
      ENDDO


      DO I = 1,Iwave
         DO K = 0,Noeud
            W_temp(K,I) = DCMPLX(0.D0,0.D0)
         ENDDO
      ENDDO

      Y =0.D0

      DO I = 1,Iwave
         IF(NI(I) .GT. 1.e-5) Y = Y +NI(I)  
      ENDDO
      DO I = 1,Iwave
         NI(I) = NI(I)/Y
      ENDDO

      N = 0

      Y = 0.D0
      DO I = 1,Iwave
         IF(NI(I) .GE. 0.01D0) THEN
            N = N + 1
            Y = Y +NI(I)  
            R(N) = NI(I)


            DO K = 0,Noeud
               DO J = 1,Iwave
                  IF( I .NE. J) THEN
                    IF(NI(I) .NE. NI(J)) THEN
                       W_temp(K,N) = W_temp(K,N) +DENS(J,I)
     &                             *wave(K,J)/(NI(I)-NI(J))
                    ELSE
                       print *,'DEGENERESCENCES'
                    ENDIF
                  ELSE
                    W_temp(K,N) = W_temp(K,N) + wave(K,I)
                  ENDIF
               ENDDO
            ENDDO

         ENDIF
      ENDDO


      Norbital = N
c      print *,'N=',N

cccccccccccccccccccccccccccccccccccc
c    NORMALISTATION OF STATES      c
cccccccccccccccccccccccccccccccccccc

      DO N = 1,Norbital
         NI(N) = R(N)/Y
         print *,'NI(',N,')=',NI(N)

         NORM = 0. 
         DO K = 0,Noeud
            NORM = NORM + CDABS(W_temp(K,N))**2.D0
         ENDDO
         NORM = 1/DSQRT(NORM*pas)
         DO K = 0,Noeud
            wave(K,N) = w_temp(K,N)*NORM
         ENDDO

      ENDDO
cccccccccccccccccccccccccccccccccccc
c    ORTHONORMALISATION OF STATES cc
cccccccccccccccccccccccccccccccccccc

      call SCHMIDT(Norbitalmax,Norbital,pas,wave)

ccccccccccccccccccccccccccccc
c    test
ccccccccccccccccccccccccccccc

c      DO I = 1,Norbital
c            NORM = 0.
c            DO J = 1,Norbital
c               NORM = 0.
c               DO K = 0,Noeud
c                  NORM = NORM + DCONJG(WAVE(K,I))
c     &                         *WAVE(K,J)
c               ENDDO
c               print *,'ORTH',I,J,'=',NORM*pas
c           ENDDO
c      ENDDO
c      DO I = 1,Norbital
c         write(nb,106) I
c         OPEN(UNIT=10,file='TEST'//nb//'.DAT',status='unknown')
c         DO K = 0,Noeud
c            write(10,*) pas*DFLOAT(K),Dreal(wave(K,I))
c         ENDDO
c         CLOSE(10)
c      ENDDO
c

c      STOP

      end
         

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculation of V(i,j,k,l)=<ij|V|kl>           x  
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c                                                   x
c IF T0_OR_NOT = 1 we calculate all but <ip|V|pj>   x
c                                       <ip|V|jp>   x
c                                                   x
c     where p = particle states                     x
c                                                   x
c IF T0_OR_NOT = 1 we calculate all but <ip|V|pl>   x
c                                       <ip|V|kp>   x
c                                       <pj|V|pl>   x
c                                       <pj|V|kp>   x
c                                                   x
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine TWO_BODY_POT(Norbitalmax,Norbital,Iwave,Noeud
     &                        ,wave_temp,V,GAU,T0_OR_NOT,pas)
                              
      IMPLICIT NONE

      INTEGER    Norbitalmax,Iwave,Noeud,T0_OR_NOT,Norbital
      real*8     GAU(0:Noeud,0:Noeud),pas
      COMPLEX*16  wave_temp(0:Noeud,Norbitalmax)
      COMPLEX*16 V(Norbitalmax,Norbitalmax
     &              ,Norbitalmax,Norbitalmax)
                    
      INTEGER I,J,K,L,I_1,I_2

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c  THIS COULD BE OPTIMIZED MAYBE MORE !
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      IF (T0_OR_NOT .EQ. 1) THEN
         DO L = 1,Iwave
c            print *,'L=',L
            DO K = 1,Iwave
               DO J =1,Iwave
                  DO I = 1,Iwave
                     V(I,J,K,L) = DCMPLX(0.D0,0.D0)
                     IF((J .GT. Norbital .AND. K .GT. Norbital)
     &                                   .OR.
     &                  (J .GT. Norbital .AND. L .GT. Norbital))
     &                  GOTO 210

                     DO I_1 = 0,Noeud
                        DO I_2 = 0,Noeud
                              V(I,J,K,L) = V(I,J,K,L) + GAU(I_1,I_2)
     &                                     *DCONJG(wave_temp(I_1,I))
     &                                     *DCONJG(wave_temp(I_2,J))
     &                                     *wave_temp(I_1,K)
     &                                     *wave_temp(I_2,L)

c**************************
c  ! GAU CONTAINS PAS**2  *
***************************
                        ENDDO
                     ENDDO 
 210                 CONTINUE
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
       ENDIF

      IF (T0_OR_NOT .EQ. 0) THEN
         DO L = 1,Iwave
c            print *,'L=',L
            DO K = 1,Iwave
               DO J =1,Iwave
                  DO I = 1,Iwave
                     V(I,J,K,L) = DCMPLX(0.D0,0.D0)
                     IF((J .GT. Norbital .AND. K .GT. Norbital)
     &                                   .OR.
     &                  (J .GT. Norbital .AND. L .GT. Norbital)
     &                                   .OR.
     &                  (I .GT. Norbital .AND. K .GT. Norbital)
     &                                   .OR.
     &                  (I .GT. Norbital .AND. L .GT. Norbital))
     &                  GOTO 220

                     DO I_1 = 0,Noeud
                        DO I_2 = 0,Noeud
                              V(I,J,K,L) = V(I,J,K,L) +GAU(I_1,I_2)
     &                                     *DCONJG(wave_temp(I_1,I))
     &                                     *DCONJG(wave_temp(I_2,J))
     &                                     *wave_temp(I_1,K)
     &                                     *wave_temp(I_2,L)
c**************************
c  ! GAU CONTAINS PAS**2  *
***************************
                        ENDDO
                     ENDDO
c                    PRINT*,'V(I,J,K,L',')=',V(I,J,K,L)
 220                 CONTINUE
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
       ENDIF
      end 
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TWO BODY COLLISION TERM IN R-REPRESENTATION
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine FILL_POT(Noeud,pas,t0,GAU)

      IMPLICIT NONE

      INTEGER    Noeud,I,J
      real*8     pas,t0
      real*8     GAU(0:Noeud,0:Noeud)
      real*8     SIGMA,PI,CTE1,CTE2

      SIGMA = 2.
      PI = DACOS(-1.D0)
      CTE1 = pas**2.D0*t0/(DSQRT(2*PI)*SIGMA) 
      CTE2 = -pas**2.D0/(2*sigma**2.D0)

      DO I = 0,Noeud
         DO J = 0,Noeud
            GAU(J,I) = CTE1*DEXP(CTE2*DFLOAT(J-I)**2.D0)  
         ENDDO
      ENDDO
      end

      subroutine diag_dens(Norbital,Iwave,pas,wave,DENS,NI)
                           
      IMPLICIT NONE

      INTEGER     Noeud,Norbitalmax,Norbital,Iwave

      PARAMETER   (Noeud= 2**6-1)
      PARAMETER   (Norbitalmax = 50)

      real*8      pas,Y
      
      real*8      NORM,NI(Norbitalmax) 
      COMPLEX*16  wave(0:Noeud,Norbitalmax)
      COMPLEX*16  DENS(Norbitalmax,Norbitalmax)
      COMPLEX*16  w_temp(0:Noeud,Norbitalmax)

      INTEGER     N,K,I,J

cccccccccccccccccccccccc
c   NAG VARIABLES      c
cccccccccccccccccccccccc

      real*8 AR(Norbitalmax,Norbitalmax)
      real*8 AI(Norbitalmax,Norbitalmax)
      real*8 R(Norbitalmax),VR(Norbitalmax,Norbitalmax)
      real*8 VI(Norbitalmax,Norbitalmax) 
      real*8 WK1(Norbitalmax),WK3(Norbitalmax),WK2(Norbitalmax)
      
      character*2 nb

      INTEGER IFAIL

      EXTERNAL F02AXF

c      print *,'Norbital =',Norbital
c      print *,'Iwave=',Iwave,'pas',pas

      IFAIL = 1
cccccccccccccccccccccccccccccccccccccc
c   TEST OF THE DIAGONALISATION      c
cccccccccccccccccccccccccccccccccccccc

      Y = 0.D0
      DO N = 1,Iwave
         Y = Y + DREAL(DENS(N,N))
         DO K = 1,Iwave
            AR(K,N) = DREAL(DENS(K,N))
            AI(K,N) = DIMAG(DENS(K,N))          
         ENDDO
      ENDDO

c      print *,'Y=',Y

      call F02AXF(AR,Norbitalmax,AI,Norbitalmax
     &            ,Iwave,R,VR,Norbitalmax,VI
     &            ,Norbitalmax,WK1,WK2,WK3,IFAIL)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   BE CAREFULL THE EIGENVALUES ARE IN ASCENDING ORDER
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO I = 1,Iwave
         DO K = 0,Noeud
            W_temp(K,I) = DCMPLX(0.D0,0.D0)
         ENDDO
      ENDDO


cccccccccccccccccccccccccccccccccc
c  !!   vp = vp*pas
c  !!   wave = wave/sqrt(pas)
cccccccccccccccccccccccccccccccccc

      DO I = 1,Iwave
         IF( R(I) .LT. 0) THEN
c            PRINT *,'ni values less than zero:'
            print *,'R(',I,')=',R(I)  
         ENDIF
c         print *,'R(',I,')=',R(I)  
      ENDDO
      N = 0
      Y = 0.
      DO I = 1,Iwave
         IF ( R(I) .GT. 0.005) THEN
            N = N + 1
            NI(N) = R(I)
            Y = Y + NI(N)
            DO K = 0,Noeud
               DO J = 1,Iwave
                  W_temp(K,I) = W_temp(K,I)  
     &                        + DCMPLX(VR(J,I)*
     &                          DREAL(wave(K,J))-
     &                          VI(J,I)*DIMAG(wave(K,J)),
     &                          VR(J,I)*DIMAG(wave(K,J))+
     &                          VI(J,I)*DREAL(wave(K,J)))
               ENDDO
            ENDDO
         ENDIF
      ENDDO

      Norbital = N
      N = 0

      DO I = 1,Iwave
         IF ( R(I) .GT. 0.005) THEN
             N = N + 1
             NI(N) = NI(N)/Y
             print *,'NI(',N,')=',NI(N)

             NORM = 0. 
             DO K = 0,Noeud
                 wave(K,N) = W_temp(K,I) 
                 NORM = NORM + CDABS(W_temp(K,I))**2.D0
             ENDDO
             NORM = 1/DSQRT(NORM*pas)
             DO K = 0,Noeud
                 wave(K,N) = wave(K,N)*NORM
             ENDDO
         ENDIF
      ENDDO
 102  FORMAT(I2)

      print *,'N=',Norbital

c      DO I = 1,Norbital
c         write(nb,102) I
c         OPEN(UNIT=10,file='TEST'//nb//'.DAT',status='unknown')
c         DO K = 0,Noeud
c            write(10,*) pas*DFLOAT(K),CDABS(WAVE(K,I))
c         ENDDO
c         CLOSE(10)
c      ENDDO
      end


