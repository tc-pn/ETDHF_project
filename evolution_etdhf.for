      subroutine evolution(wave,pot,I_TEMPSMAX,delta_t
     &                     ,t0,densite,ni,pas,Norbital)  
                             
      IMPLICIT NONE

      integer Noeud,Norbital,I_TEMPSMAX,INDEX,Norbitalmax
      integer N_coll
      parameter(Noeud =  2**6-1,Norbitalmax= 50)
      integer N_time
      parameter(N_time = 4,N_coll=10)
c      parameter(N_time = 4,N_coll=1)
      COMPLEX*16 wave(0:Noeud,Norbitalmax),FTEMP(0:Noeud,Norbitalmax)  
      real*8 POT(0:Noeud),densite(0:Noeud),t0,delta_t,pas
      real*8 HBAR,HB,DU,ni(Norbitalmax),XMU,X,U(0:Noeud)
      real*8 densite_real(0:Noeud)     
      real*8 GAU(0:Noeud,0:Noeud)
      real*8 H1,Y,TIME,XMOY,X2moy,Etot,Entropy
      real*8 CIN_COS(0:Noeud)
      real*8 CIN_SIN(0:Noeud)
      real*8 POT_COS(0:Noeud)
      real*8 POT_SIN(0:Noeud)
      real*8 MF(0:Noeud)
      real*8 T_MACRO,OVERLAP
      real*8 E_NIV(1+Noeud)

      character*2 nb

      INTEGER IT,n,k,j,i

      DATA HBAR/197.625D0/,HB/20.7335D0/

  100 FORMAT(I2)
      H1 = HB/pas**2
      IT=0
      DU=delta_t/HBAR

ccccccccccccccccccccccccccccccccccc
c  MACROSCOPIC COLLISION TERM     c
ccccccccccccccccccccccccccccccccccc

      T_MACRO = N_coll*delta_t

ccccccccccccccccccccccccccccccccccccccccccccccccccc
c    FILES SPECIFICATIONS
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  SPLIT OPERATOR METHOD EVOLUTION
c

      call FILL_POT(Noeud,pas,t0,GAU)

      call COEF_CIN(Noeud,pas,delta_t
     &                ,CIN_COS,CIN_SIN)



      open(unit=20,file='x_etdhfassnew.dat',status='unknown')
      open(unit=22,file='sx_etdhfassnew.dat',status='unknown')      
      open(unit=21,file='E_etdhfassnew.dat',status='unknown')
      open(unit=24,file='S_etdhfassnew.dat',status='unknown')
      DO I=1,10
         write(nb,100) I
         OPEN(UNIT=24+I,file='ENER'//nb//'assnew.dat',STATUS='UNKNOWN')
      ENDDO
      DO I=1,10
         write(nb,100) I
         OPEN(UNIT=35+I,file='ERHO'//nb//'assnew.dat',STATUS='UNKNOWN')   
         OPEN(UNIT=45+I,file='NI'//nb//'assnew.dat',STATUS='UNKNOWN')   
      ENDDO
      OPEN(UNIT=56,file='overlapassnew.dat',STATUS='UNKNOWN')   


      do n=1,Norbital
        do k=0,noeud
           FTEMP(k,n) = WAVE(k,n)
        enddo
      enddo



      INDEX=1
  101 INDEX=-INDEX
      XMU=float(1+INDEX)/4.



      do n=1,Norbital
        Y=0.
        do k=0,noeud
            Y=Y+CDABS(FTEMP(k,n))**2
        enddo      
        Y=1./sqrt(Y*pas)
        do k=0,noeud
          FTEMP(k,n)=FTEMP(k,n)*Y
        enddo
      enddo                                         

      IF (IT .GT. I_TEMPSMAX) GOTO 103
ccccccccccccccccccccccccccccccccccccccccccccccc
c    CALCULATION OF THE DENSITY
ccccccccccccccccccccccccccccccccccccccccccccccc
 

      call calcul_density(Noeud,Norbitalmax,FTEMP,ni,densite,pas,
     &                    densite_real,Norbital)   

cccccccccccccccccccccccccccccccccccccccccccccc
c    CALCULATION OF THE POTENTIAL
cccccccccccccccccccccccccccccccccccccccccccccc

      call MEAN_FIELD(Noeud,pas,densite,t0,MF)  


      do J = 0,Noeud
         X= MF(j) + POT(J)
         U(J)=(1.-XMU)*X+XMU*U(J)
      enddo

c      OPEN(unit=30,file='POT.DAT',status='unknown')
c      DO J = 0,Noeud
c         WRITE(30,*) PAS*DFLOAT(J),U(J)
c      ENDDO
c      close(30) 
c      PAUSE
      call COEF_POT(Noeud,delta_t
     &              ,U,POT_COS,POT_SIN)   


      IF(INDEX.eq.1) GO TO 112
         IF(IT .NE. 0) TIME = TIME + delta_t
         IT=IT+1

         IF(MOD(IT-1,20) .EQ.0) print *,'IT=',IT
cccccccccccccccccccccccccccccccccccccccccccccccccc
c   ENERGETIQUE CONSIDERATIONS
ccccccccccccccccccccccccccccccccccccccccccccccccccc

         IF( MOD(IT,N_coll) .EQ. 0) THEN
           call DIAG_EN(U,H1,E_NIV)  
           DO I=1,10
              WRITE(24+I,*) TIME,E_NIV(I)
           ENDDO
           call ENERGY_F(Noeud+1,Norbital,pas,wave,U,E_NIV)
           DO I=1,Norbital
              write(35+I,*) TIME,E_NIV(I)
              write(45+I,*) TIME,NI(I)
           ENDDO
         ENDIF


cccccccccccccccccccccccccccccccccccccccccccccccccccc
c    THE COLLISION TERM
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        TEST
c

         IF( MOD(IT,N_coll) .EQ. 0) THEN

            call COLLISION(Norbital,pas,delta_t
     &                     ,FTEMP,NI,U,CIN_COS,CIN_SIN
     &                     ,POT_COS,POT_SIN,GAU,T_MACRO
     &                     ,OVERLAP,it)
            write(56,*) TIME,OVERLAP
         ENDIF


         do i=1,Norbital
            do k=0,Noeud
              wave(k,i)=FTEMP(k,i)
            enddo
         enddo 
cccccccccccccccccccccccccccccccccccccccccccccc
c    calculation of the position
cccccccccccccccccccccccccccccccccccccccccccccc
        IF(MOD(IT,N_time) .EQ. 0) THEN
           Xmoy  = 0.D0
           X2moy = 0.D0
           DO J=0,Noeud
              Xmoy = Xmoy + DFLOAT(J)*densite_real(J)
              X2moy = X2moy + DFLOAT(J)**2.D0*densite_real(J)
           ENDDO
           Xmoy  = Xmoy*pas**2.D0
           X2moy = X2moy*pas**3.D0
           write(20,*) TIME,Xmoy
           write(22,*) TIME,DSQRT(X2moy - Xmoy**2.D0)  
c            write(22,*) TIME,DSQRT(X2moy)  
cccccccccccccc
c  ENTROPY 
cccccccccccccc
           CALL calcul_entropy(Norbital,Norbitalmax,NI,Entropy)
           write(24,*) TIME,entropy  

        ENDIF

        IF( (MOD(INT(IT-1),24*N_time) .EQ. 0)) THEN
            WRITE(nb,100) INT((IT-1)/(24*N_time))
            OPEN(UNIT=13,file='dens'//nb//'etdhfassnew.dat'
     &                  ,status='unknown')
            DO J=0,Noeud
               write(13,*) pas*DFLOAT(J),densite_real(J)
c               write(13,*) pas*DFLOAT(J),DREAL(WAve(J,1))
            ENDDO
            CLOSE(13)
        ENDIF



        call CALCUL_ENERGY(Norbitalmax,pas,POT,H1,WAVE,densite
     &                          ,NI,Etot,MF,Norbital)

c        IF(MOD(INT(TIME),100).EQ.0) print *,'Etot=',Etot

        write(21,*) TIME,Etot
  112 CONTINUE
   
      call EVOLUTION_2(Norbitalmax,WAVE,FTEMP
     &                 ,CIN_COS,CIN_SIN
     &                 ,POT_COS,POT_SIN,Norbital)

      GOTO 101

  103 CONTINUE
    
      close(20)
      close(21)
      close(22)
      close(24)
      DO I=1,10
         close(24+I)
         close(35+I)
         close(45+I)
      ENDDO

      end


      subroutine calcul_density(Noeud,Norbitalmax,FTEMP,ni,densite,pas,
     &                          densite_real,Norbital)     

      IMPLICIT NONE
      INTEGER Noeud,Norbitalmax ,Norbital

      COMPLEX*16 FTEMP(0:Noeud,Norbitalmax)
      real*8 NI(Norbitalmax),NORM,pas
      real*8 densite(0:Noeud),densite_real(0:Noeud)
      INTEGER I,J
      DO I = 0,Noeud
         densite(I) = 0.D0
         densite_real(I) = 0.D0
      ENDDO
      DO J=1,Norbital
         DO I = 0,Noeud
            densite(I) = densite(I) + CDABS(FTEMP(I,J))**2.D0*NI(J)
         ENDDO
      ENDDO

      NORM = 0.D0
      DO I = 0,Noeud
        densite(I) = densite(I)/2.D0
        NORM = NORM + densite(I)
      ENDDO
      NORM = NORM*pas 
      DO I = 0,Noeud
         densite_real(I) = densite(I)/NORM  
      ENDDO

      end

      
       subroutine CALCUL_ENERGY(Norbitalmax,pas,POT,H1,WAVE,densite
     &                          ,NI,Etot,MF,Norbital)   

       IMPLICIT NONE

       INTEGER   Noeud,Norbitalmax,Norbital
c       parameter(Noeud =  3*(2**8-1))
       parameter(Noeud =  2**6-1)
c       parameter(Noeud =  2**9-1)

       real*8    pas,POT(0:Noeud),H1,MF(0:Noeud)
       real*8    densite(0:Noeud),Etot,NI(Norbitalmax)
       complex*16 Wave(0:Noeud,Norbitalmax) 
       complex*16 AX(0:noeud),BX(0:noeud),CX(0:noeud),DX(0:noeud)
       INTEGER   N,K

       Etot = 0.D0
cccccccccccccccccccccccccccccccccccccccccccccccccc
c     KINETIK PART
cccccccccccccccccccccccccccccccccccccccccccccccccc


       do k=0,noeud 
          AX(k)=DCMPLX(-H1,0.D0)
          BX(k)=DCMPLX(2*H1+POT(k),0.D0)
          CX(k)=DCMPLX(-H1,0.D0)
       enddo

       do n=1,Norbital
          DX(0)=BX(0)*WAVE(0,n)+CX(0)*WAVE(1,n)
          do k=1,(noeud-1)
             DX(k)=AX(k)*WAVE(k-1,n)+BX(k)*WAVE(k,n)
     &       +CX(k)*WAVE(k+1,n)
          enddo
          DX(noeud)=AX(noeud)*WAVE(noeud-1,n)
     &             +BX(noeud)*WAVE(noeud,n)

          do k=0,noeud
             Etot=Etot+DX(k)*DCONJG(WAVE(k,n))*ni(n)
          enddo
       enddo

       Etot = Etot*2.D0
ccccccccccccccccccccccccccccccccccccccccccccccccc
c      POTENTIAL PART
ccccccccccccccccccccccccccccccccccccccccccccccccc

       DO k=0,Noeud
          Etot = Etot + 2.D0*MF(K)*DENSITE(K)
       ENDDO

       Etot = Etot*pas

c       print *,'Etot=',Etot

       end


      subroutine COEF_CIN(network_size,delta_x,delta_t
     &                ,CIN_COS,CIN_SIN)    

      IMPLICIT NONE
      
      INTEGER network_size,I,I_TEMP
      REAL*8  delta_x,delta_t 
      real*8 CIN_COS(0:network_size)
      real*8 CIN_SIN(0:network_size)

      REAL*8 PI,HB,MASS,NORM,VALUE

      DATA HB/197.32705D0/,MASS/938.91897D0/


      PI = DACOS(-1.D0)
       
      NORM = HB*PI**2*delta_t/((network_size+1)**2*delta_x**2*MASS)

      DO I = 0,network_size
c
c   essai
c
         IF(I .LE. ((network_size+1)/2)) then
            I_TEMP = I
         ELSE
            I_TEMP = (network_size+1)-I
         ENDIF

         VALUE  = NORM*DFLOAT((I_TEMP)**2)
         CIN_COS(I) = DCOS(VALUE)
         CIN_SIN(I) = DSIN(VALUE)
      ENDDO

      end

      subroutine COEF_POT(network_size,delta_t
     &                ,POT,POT_COS,POT_SIN)     

      IMPLICIT NONE
      
      INTEGER network_size,I

      real*8 delta_t,HB,VALUE,NORM
      real*8 POT(0:network_size)
      real*8 POT_COS(0:network_size)
      real*8 POT_SIN(0:network_size)

      DATA HB/197.32705D0/

      NORM = delta_t/HB

      DO I = 0,network_size
            VALUE  =  NORM*POT(I)             
            POT_COS(I) = DCOS(VALUE)
            POT_SIN(I) = DSIN(VALUE)
      ENDDO

      end
                       
C
C    2-Dimensionnal SPLIT-OPERATOR Calculation: EVOLUTION
C    Use only real and imaginary densite of wave function
c

      subroutine evolution_2(Norbitalmax,WF,WF1
     &                       ,CIN_COS,CIN_SIN
     &                       ,POT_COS,POT_SIN,Norbital)
                            
      IMPLICIT NONE

      INTEGER   network_size,Norbitalmax,Norbital
      parameter(network_size =  2**6-1)
c      parameter(network_size =  2**9-1)

c      parameter (network_size = 3*(2**8-1))


      real*8 REAL_PART(0:network_size)
      real*8 IMAG_PART(0:network_size)
      real*8 CIN_COS(0:network_size)
      real*8 CIN_SIN(0:network_size)
      real*8 POT_COS(0:network_size)
      real*8 POT_SIN(0:network_size)
      complex*16 WF(0:network_size,Norbitalmax)
      complex*16 WF1(0:network_size,Norbitalmax)

      INTEGER   N_SIZE,I,J

      N_SIZE = (network_size+1)
      DO J=1,Norbital
         DO I = 0,network_size
            REAL_PART(I) = DREAL(WF(I,J))
            IMAG_PART(I) = DIMAG(WF(I,J))
         ENDDO

         call FFT_CIN(network_size,REAL_PART,IMAG_PART 
     &               ,N_SIZE,CIN_COS,CIN_SIN) 
        
         call POT_TERM(network_size,REAL_PART,IMAG_PART,
     &                 POT_COS,POT_SIN)  

         call FFT_CIN(network_size,REAL_PART,IMAG_PART 
     &               ,N_SIZE,CIN_COS,CIN_SIN) 

         DO I = 0,network_size
            WF1(I,J) = DCMPLX(REAL_PART(I)
     &                       ,IMAG_PART(I))
         ENDDO
      ENDDO

      end

      subroutine FFT_CIN(network_size,REAL_PART,IMAG_PART
     &                   ,N_SIZE,CIN_COS,CIN_SIN) 
                 

      IMPLICIT NONE

      INTEGER network_size,N_SIZE
      REAL*8 REAL_PART(0:network_size)
      REAL*8 IMAG_PART(0:network_size)
      real*8 CIN_COS(0:network_size)
      real*8 CIN_SIN(0:network_size)

      INTEGER I

      INTEGER NDIM,DIM_MAX,LWORK
      PARAMETER (NDIM = 1,DIM_MAX = 520, LWORK=3*DIM_MAX)
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

      DO I = 0,network_size

         X_TEMP = REAL_PART(I)  

         REAL_PART(I) = CIN_COS(I)*REAL_PART(I)  
     &                 +CIN_SIN(I)*IMAG_PART(I)
         IMAG_PART(I) = CIN_COS(I)*IMAG_PART(I)
     &                 -CIN_SIN(I)*X_TEMP

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

      INTEGER network_size,I
      real*8 REAL_PART(0:network_size)
      real*8 IMAG_PART(0:network_size)
      real*8 POT_COS(0:network_size)
      real*8 POT_SIN(0:network_size)
      real*8 X_TEMP

      DO I = 0,network_size

         X_TEMP = REAL_PART(I)  

         REAL_PART(I) = POT_COS(I)*REAL_PART(I)  
     &                  +POT_SIN(I)*IMAG_PART(I)  


         IMAG_PART(I) = POT_COS(I)*IMAG_PART(I)
     &                    -POT_SIN(I)*X_TEMP

      ENDDO

      end

      subroutine MEAN_FIELD(Noeud,pas,densite,t0,MF)
  
      IMPLICIT NONE

      INTEGER   Noeud
      real*8    pas,densite(0:Noeud),t0,MF(0:Noeud)

      INTEGER I,J
      real*8 PI,SIGMA,CTE1,CTE2

      SIGMA = 2.
      PI = DACOS(-1.D0)
      CTE1 = 2.D0*pas*t0/(DSQRT(2*PI)*SIGMA) 
      CTE2 = -pas**2.D0/(2*sigma**2.D0)

      DO I = 0,Noeud
         MF(I) = 0.D0
      ENDDO

      DO I = 0,Noeud
         DO J = 0,Noeud
            MF(I) = MF(I) + CTE1*DEXP(CTE2*DFLOAT(I-J)**2.D0)
     &                      *densite(J)
         ENDDO
      ENDDO
      end


      subroutine calcul_entropy(Norbital,Norbitalmax,NI,Entropy)
                        
      IMPLICIT NONE

      integer  Norbitalmax,Norbital
      real*8   ni(Norbitalmax),Entropy
      integer  I

      Entropy = 0.
      DO I = 1,Norbital
         IF( NI(I) .NE. 1) THEN
            Entropy = Entropy 
     &             + ni(i)*DLOG(ni(i))
     &             +(1-ni(i))*DLOG(1-ni(i))  
         ENDIF
      ENDDO

      Entropy = -Entropy

      end


      subroutine DIAG_EN(MF,H1,E_NIV)     

      IMPLICIT NONE
      INTEGER noeud,J   
      PARAMETER (Noeud=2**6-1)
      real*8 pas,H1,X02AJF
      real*8 MF(0:Noeud)
      real*8 E_NIV(Noeud+1)
      real*8 EPS,D(Noeud+1),E(Noeud+1),IFAIL

      external F02AVF
      external X02AJF

      EPS    = X02AJF()

      IFAIL  = 1
    
      do j = 1,Noeud+1
         E(j) = -H1
         D(j) = 2*H1+MF(j-1)
      enddo
          
      call F02AVF(Noeud+1,EPS,D,E,IFAIL)

      do J = 1, Noeud+1
         E_NIV(J) = D(J)
      enddo
      end
