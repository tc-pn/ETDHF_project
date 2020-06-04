ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     program of Resolution of two coupled particles
c     1."Exact Resolution"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program principal

      IMPLICIT NONE

      integer     noeud,A,Nmax,Nwave
      real*8      pas,K,t0,delta_t,FORCE
      PARAMETER   (pas=0.3,noeud=2**6-1,Nmax=400)
      INTEGER    N_time,N_timemax,N_ref
      parameter  (N_time = 4,N_timemax=4000,N_ref=100)
ccccccccccccccccccccccccccccccccc
c  PUT Nmax = Nwave
ccccccccccccccccccccccccccccccccc
      COMPLEX*16  WAVE(0:Noeud,0:Noeud)
      COMPLEX*16  DENS(0:Noeud,0:Noeud)
      COMPLEX*16  DENS_T(0:Noeud,0:Noeud,N_ref) 
      real*4      X
      real*8      NI(Nmax),DD(Noeud+1),Y
      real*8      POT(0:Noeud,0:Noeud),densite(0:Noeud) 
      real*8      dens_moy(0:Noeud,N_ref) 
      real*8      CIN_COS(0:Noeud,0:Noeud)
      real*8      CIN_SIN(0:Noeud,0:Noeud)
      real*8      POT_COS(0:Noeud,0:Noeud)
      real*8      POT_SIN(0:Noeud,0:Noeud)
      real*8      REAL_PART(0:Noeud,0:Noeud)
      real*8      IMAG_PART(0:Noeud,0:Noeud)
      real*8      Xmoy,Pmoy,TIME,SIGMAX,SIGMAP,SIGMAA,POS,Etot
      integer     I,I_TEMPSMAX,J,Iwave,K1
      real*8      Entropy

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ANALYSING BEHAVIOUR
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
      real*8      X_average(N_timemax),SX_AVERAGE(N_timemax) 
                  
      character*2 nb,nb1


  100 FORMAT(I2)

      call INITIALISATION(noeud,wave,pas,pot,K,t0,delta_t
     &                    ,CIN_COS,CIN_SIN,I_TEMPSMAX
     &                    ,POT_COS,POT_SIN,FORCE)

c      GOTO 102
cccccccccccccccccccccccccccccccccccccc
c   IF WE WANT TO DIAGONALIZE D
cccccccccccccccccccccccccccccccccccccc

c      OPEN(UNIT=12,file='POT.DAT',status='unknown')
c      DO I = 0,Noeud
c         write(12,*) pas*DFLOAT(I),POT(I,I)
c      ENDDO
c      close(12)

      call diagonalize2d(POT,pas,Nwave,FORCE)   

c      STOP

ccccccccccccccccccccccccccccccccccccccccccccccccc
c  READ OCCUPATION PROBABILITIES
ccccccccccccccccccccccccccccccccccccccccccccccccc

 102  CONTINUE

      PRINT *,'ATTENTION-HERE PUT NWAVES'	

c      Nwave = 1
ccccccccccccccccccc
c  symetric
ccccccccccccccccccc
c      Nwave = 13
ccccccccccccccccccc
c  assymmetric
ccccccccccccccccccc

      Nwave = 12
 
      PRINT *,'Nwave=',Nwave
      DO J = 0,Noeud
         DO I = 0,Noeud
            DO K1 = 1,N_ref
               DENS_T(I,J,K1) = DCMPLX(0.D0,0.D0)
            ENDDO
         ENDDO
      ENDDO
      DO J = 1,N_ref
         DO I = 0,Noeud
            dens_moy(I,J) = 0.
         ENDDO
      ENDDO

      OPEN(UNIT=21,file='PIass.DAT',status='OLD')

      DO Iwave=1,Nwave
         read (21,*) NI(Iwave)
      ENDDO
      close(21)

      DO Iwave = 1,Nwave

         TIME = 0.

         write(nb,100) Iwave
         print *,'Wave',Iwave
ccccccccccccccccccccccccccccccccccccccccccccccccc
c  READ WAVE-PACKETS
ccccccccccccccccccccccccccccccccccccccccccccccccc
         OPEN(UNIT=10,file='PSIass'//nb//'.dat',status='unknown')
         DO I = 0,Noeud
            DO J = 0,Noeud           
               read(10,*) X
               WAVE(I,J) = DCMPLX(DBLE(X),0.D0)
            ENDDO
         ENDDO
         close(10)

c         OPEN(UNIT=12,file='x_twave'//nb//'ass.dat',status='unknown')
c         OPEN(UNIT=14,file='r_twave'//nb//'ass.dat',status='unknown')
c         OPEN(UNIT=15,file='e_twave'//nb//'ass.dat',status='unknown')

         DO I = 1,I_TEMPSMAX  
             IF(MOD(I-1,100) .EQ. 0)  print*,'I=',I

             call NORMALISE(noeud,wave,pas) 
   
             IF((MOD((I-1),N_time) .EQ. 0)) THEN 
                 call calcul_densite(pas,noeud,wave,densite)
                 call position_rel(pas,noeud,wave,POS)

                 IF (MOD(INT(I-1),(24*N_time)) .EQ. 0) THEN
c                 IF( (MOD(INT(I-1),(24*N_time)) .EQ. 0) .
c     &               AND. (TIME.LE.1001)) THEN

                     WRITE(nb1,100) INT((I-1)/(24*n_time)) 
c                     OPEN(UNIT=13,file='dens'//nb1//'wave'//nb//
c     &                                 'ass.dat',status='unknown')

                     call calcul_dens(pas,noeud,wave,dens)

                     DO J = 0,Noeud  
                        DO K1 = 0,Noeud
                           DENS_T(K1,J,(I-1)/(24*N_time)+1) = 
     &                        DENS_T(K1,J,(I-1)/(24*N_time)+1) 
     &                       + DENS(K1,J)*NI(Iwave)
                        ENDDO
                     ENDDO
  
                     DO J=0,Noeud  
c                        write(13,*) pas*DFLOAT(J),densite(J)
                        dens_moy(J,(I-1)/(24*N_time)+1) =
     &                    dens_moy(J,(I-1)/(24*n_time)+1) 
     &                    + densite(J)*NI(Iwave) 
                     ENDDO
c                     CLOSE(13)
                 ENDIF

                 call POSITION(noeud,pas,densite,WAVE,
     &                         Xmoy,Pmoy,SIGMAX,SIGMAP,SIGMAA)    

c                 write(12,*) TIME,Xmoy
c                 write(14,*) TIME,POS
                 X_average((I-1)/N_time+1) = 
     &                    X_average((I-1)/N_time+1) +NI(Iwave)*Xmoy
                 SX_AVERAGE((I-1)/N_time+1) = 
     &                    SX_AVERAGE((I-1)/N_time+1)+NI(Iwave)*SIGMAX 

                 call ENERGY(noeud,Wave,Etot,pas,POT)   
c                 write(15,*) TIME,Etot
             ENDIF

             call SPLIT_OPERATOR(noeud,Wave,Wave
     &                          ,REAL_PART,IMAG_PART
     &                         ,CIN_COS,CIN_SIN
     &                         ,POT_COS,POT_SIN)
             TIME = TIME + delta_t
         ENDDO
c         close(12)
c         close(14)
c         close(15)

      ENDDO
cccccccccccccccccccccccccccccccccccc
c   END OF THE LOOPS on WAVE_PACKET
cccccccccccccccccccccccccccccccccccc
      OPEN(UNIT = 16,file = 'Xexactass.dat',status='unknown')
      OPEN(UNIT = 17,file = 'SXexactass.dat',status='unknown')
      OPEN(UNIT = 18,file = 'NIexactass.dat',status='unknown')
      OPEN(UNIT = 200,file = 'entropyexactass.dat',status='unknown')
      TIME = 0.
      DO I = 1,I_TEMPSMAX  
         IF(MOD((I-1),N_time) .EQ. 0) THEN        
                                    print*, 'loul',TIME    
            write(16,*) TIME,X_average((I-1)/N_time+1) 
            write(17,*) TIME,DSQRT(SX_AVERAGE((I-1)/N_time+1)
     &                            -X_average((I-1)/N_time+1)**2.D0)
         ENDIF
         TIME = TIME + delta_t

c         IF (MOD(INT(I-1),(N_time)) .EQ. 0) THEN
         IF (MOD(INT(I-1),(24*N_time)) .EQ. 0) THEN
c         IF( (MOD(INT(I-1),(24*N_time)) .EQ. 0) .
c     &        AND. (TIME.LE.1001)) THEN

c              WRITE(18,*) '(TIME =',TIME

              write(nb1,100) INT((I-1)/(24*N_time)) 
              OPEN(unit=13,file='densm'//nb1//'ass.dat',status=
     &                                                 'unknown')
              DO J = 0,Noeud
                 write(13,*) pas*DFLOAT(j)
     &                       ,dens_moy(J,(I-1)/(24*n_time)+1)
              ENDDO
              CLOSE(13)
              DO K1 = 0,Noeud
                 DO J = 0,Noeud
                    DENS(J,K1) = DENS_T(J,K1,(I-1)/(24*N_time)+1) 
                 ENDDO
              ENDDO
c              call DIAG_RO(pas,DENS,DD,Entropy)
              call DIAG_RO(pas,DENS,DD)  
              Y = 0.
              DO J = 1,Noeud+1
                 IF( DD(Noeud-J+2) .GT. 0.005) THEN
                     Y = Y + DD(Noeud-J+2)
                 ENDIF
              ENDDO

              Entropy = 0.
              DO J = 1,Noeud+1
                 IF( DD(Noeud-J+2) .GT. 0.005) THEN
                    Entropy = Entropy - DD(Noeud-J+2)/Y * DLOG(DD(Noeud-J+2)/Y)
     &                      - (1. - DD(Noeud-J+2)/Y) * DLOG(1. - DD(Noeud-J+2)/Y) 
                    print*, 'lol', TIME,DD(Noeud-J+2)/Y , entropy
                    WRITE(18,*) TIME,DD(Noeud-J+2)/Y
                 ELSE
                    continue
                 ENDIF
              ENDDO
              write(200,*) TIME, Entropy
 123          CONTINUE
         ENDIF
         
      ENDDO

      CLOSE(16)
      CLOSE(17)
      CLOSE(18)
      CLOSE(200)

      end 

      include 'evolution_exact.for'

      include 'density_two.for'
