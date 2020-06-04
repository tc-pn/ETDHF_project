       subroutine DIAG_EN(Ncoupure,N0,pas,LL,NN,noeud
     &                        ,U,H1,A,NI,FONDE)

       IMPLICIT NONE
       INTEGER N0,noeud,il,J,K,IK,Ncoupure
       INTEGER LL(N0),NN(N0),L,IL_MAX,I,IVEC     
       real*8 pas,H1,X02AJF,A,mu,NI_H(N0),NI(N0)
       real*8 E_NIV(0:N0,Ncoupure),U(noeud,N0),POT(Ncoupure)
       real*8 EPS,D(Ncoupure),E(Ncoupure),IFAIL,E_ORD(N0)
       real*8 VECT_P(N0,Ncoupure),NORM
              
       real*8 ALB,UB,EPS1,E2(Ncoupure),R(N0),V(Ncoupure,N0) 
       real*8 X(Ncoupure,7)

       complex FONDE(noeud,N0)
       complex DENS(Ncoupure,Ncoupure)

       integer M,MM,II,JJ,KK

       INTEGER LL_ORD(N0),ICOUNT(N0)

       LOGICAL C(Ncoupure)

       external F02AVF
       external X02AJF
       external F02BEF

       EPS    = X02AJF()

       IFAIL  = 1
       il_MAX = 0
    
       do il=1,N0
          if(il_max .LE. LL(IL)) IL_MAX = LL(IL)
       enddo


       do IL = 1,IL_MAX+1
          L = IL-1
          K  = 1
          do while (L .NE. LL(K))
             K = K + 1
          enddo
          IK = K
          do j = 1,Ncoupure
             E(j) = -H1
             D(j) = 2*H1+U(j,IK)
          enddo
          
          call F02AVF(Ncoupure,EPS,D,E,IFAIL)
          do J = 1, Ncoupure
             E_NIV(L,J) = D(J)
          enddo      
       enddo
 
       call ORDONNE1(N0,Ncoupure,E_NIV,LL_ORD,E_ORD,IL_MAX) 

c       open(unit=56,file='el.dat',status='unknown')
c       do j=1,N0
c          write(56,*) E_ORD(J),LL_ORD(J)
c       enddo
c       close(56)



c
c   Calcul des vecteurs propres associes aux orbitales occupees
c

        ALB = E_ORD(1) - 1.
        UB  = E_ORD(N0) + 1.
        EPS = X02AJF()
        EPS1 = 0.0
        il_MAX = 0
                        
        do il=1,N0
             if(il_max .LE. LL_ORD(IL)) IL_MAX = LL_ORD(IL)
        enddo



        do IL = 1,IL_MAX+1
            L = IL-1

c
c     CALCULATION OF DENSITY ASSOCIATED TO THIS L 
c          
             do kk = 1,Ncoupure
                do jj = 1,Ncoupure
                  DENS(jj,kk) = 0.               
                enddo        
             enddo
             do kk = 1,Ncoupure  
                do jj = 1,Ncoupure
                   do ii = 1,N0
                      if(ll(ii) .eq. L) then
                        dens(jj,kk) = dens(jj,kk)+ ni(ii)
     &                     *FONDE(jj,ii)*CONJG(FONDE(kk,ii)) 
                      endif
                   enddo
                enddo
             enddo

            K  = 1
            do while (L .NE. LL(K))
                K = K + 1
            enddo
            IK = K

            do j = 1,Ncoupure
               E(j)  = -H1
               E2(j) = E(j)**2
               D(j) = 2*H1+U(j,IK)
            enddo

            E2(1) = 0.0
            E(1)  = 0.0
          
            call F02BEF(Ncoupure,D,ALB,UB,EPS,EPS1,E,E2,N0,
     &                  MM,R,V,Ncoupure,ICOUNT,X,C,IFAIL)

            do j = 1, N0
               do I = 1,MM
                  if(INT(100.*R(I)) .EQ. INT(100.*E_ORD(J))) THEN
                    NORM = 0.
                    do IVEC = 1, Ncoupure
                       VECT_P(J,IVEC) = V(IVEC,I)
c /sqrt(pas) 
                       NORM = NORM + VECT_P(J,IVEC)**2 
                    enddo
                    NORM = 1./sqrt(NORM*pas)
                    do IVEC = 1, Ncoupure
                       VECT_P(J,IVEC) = VECT_P(J,IVEC)*NORM 
                    enddo
                  endif
               enddo
            enddo
c
c   calcul des nombres d'occupation associes
c

           DO J=1,N0
              IF (LL_ORD(J) .EQ. L) THEN
                 NI_H(J) = 0.
                 DO II = 1,Ncoupure
                    DO JJ = 1,Ncoupure
                       NI_H(J) = NI_H(J) + VECT_P(J,II)   
     &                           *DENS(II,JJ)*VECT_P(J,JJ)   
                    ENDDO
                 ENDDO
                 NI_H(J) = NI_H(J)*pas**2
c                 PRINT *,'NI_H(',J,')=',NI_H(J),LL_ORD(J)
              ENDIF
           ENDDO
c
c   FIN de la boucle sur les L
c
       ENDDO
       OPEN (UNIt=53,file='EH_NIH.DAT',STATUS='UNKNOWN')
       DO J=1,N0
          WRITE(53,*) E_ORD(J),NI_H(J),LL_ORD(J)
       ENDDO
       CLOSE(53)
       end

       subroutine ORDONNE1(N0,Ncoupure,E_NIV,LL_ORD,E_ORD,IL_MAX)    


       IMPLICIT NONE

       INTEGER N0,Ncoupure,LL_ORD(N0),I,J,L,K_TEMP,IL_MAX
       INTEGER IL,K,IL_TEMP,JJ,IK
       real*8 E_NIV(0:N0,Ncoupure),E_ORD(N0),E

       J = 0

       do while (J.LT.N0)
          J = J + 1
          E = 1.0e12
          do JJ = 1, IL_MAX+1
             IL = JJ-1
             do IK = 1,Ncoupure
                K = Ncoupure+1-IK
                if((E .GT. E_NIV(IL,K)) .AND. 
     &             (INT(1000000.*E_NIV(IL,K)) .NE.0)) THEN
                 E = E_NIV(IL,K)
                 K_TEMP  = K
                 IL_TEMP = IL
               endif
             enddo
          enddo
          E_ORD(J)    = E
          LL_ORD(J)   = IL_TEMP
          E_NIV(IL_TEMP,K_TEMP) = 1.0e14
       enddo
       end

       subroutine DIAG_RO(Ncoupure,N0,pas,LL,NN,noeud
     &                        ,U,H1,NI,FONDE)

       IMPLICIT NONE

       INTEGER noeud,N0,Ncoupure,LLmax,LLL

       real*8 ni(N0),H1,K,XX
       real*8 pas,S,Entropy,U(noeud,N0),NORM      
       integer LL(N0),NN(N0)
       integer jj,kk,ii
       complex FONDE(noeud,N0)
       complex DENS(Ncoupure,Ncoupure)
       complex VECT_P(Ncoupure,N0)
c
c    parametre des routines nag
c            
       real*8 AR(Ncoupure,Ncoupure),AI(Ncoupure,Ncoupure),DD(Ncoupure)
       real*8 WK1(Ncoupure),WK2(Ncoupure),WK3(Ncoupure)
       real*8 VR(Ncoupure,Ncoupure),VI(Ncoupure,Ncoupure) 
       real*8 Hamilt(Ncoupure,Ncoupure)


       real*8 nni(Ncoupure),NI_RO(N0),E_RO(N0)
       integer ifail,COMPTEUR,LL_RO(N0)

      
       EXTERNAL F02AXF


       call LEVEL(N0,NN,LL)

       Entropy = 0
       llmax = 0

       do ii = 1,N0
          if(LL(ii) .ge. llmax) llmax=LL(ii)
       enddo             
                
       do kk = 1,Ncoupure
          nni(kk) = 0
       enddo

       COMPTEUR = 0


       do lll = 0,llmax
                      
          do kk = 1,Ncoupure
             do jj = 1,Ncoupure
                DENS(jj,kk) = 0.               
             enddo        
          enddo
          do kk = 1,Ncoupure  
            do jj = 1,Ncoupure
              do ii = 1,N0
                if(ll(ii).eq.lll) then

                  dens(jj,kk) = dens(jj,kk)+pas*ni(ii)
     &           *FONDE(jj,ii)*CONJG(FONDE(kk,ii)) 

                endif
              enddo
            enddo
          enddo

          do kk = 1,Ncoupure
            do jj = 1,Ncoupure
               AI(jj,kk) = AIMAG(DENS(jj,kk))
               AR(jj,kk) = REAL(DENS(jj,kk))
            enddo
               if(AI(kk,kk) .lt. 1.0e-10) then
                  AI(kk,kk) = 0.0
               endif  
          enddo


          call F02AXF(AR,Ncoupure,AI,Ncoupure,Ncoupure
     &                 ,DD,VR,Ncoupure,
     &                 VI,Ncoupure,WK1,WK2,WK3,IFAIL) 

           if(IFAIL.NE.0) then
               print *,'Error'       
           else
c
c     calcul du hamiltonien associe au L
c
               K = 1

               do while (LLL .NE. LL(K))
                 K = K + 1
               enddo
               do ii = 1,Ncoupure
                  do jj = 1,Ncoupure
                     Hamilt(ii,jj)=0
                  enddo
               enddo

               do ii = 1,Ncoupure
                  Hamilt(ii,ii) = 2*H1+U(II,K)
               enddo
               do ii = 1,Ncoupure-1
                  HAMILT(ii,ii+1) = -H1
               enddo
               do ii = 2,Ncoupure
                  HAMILT(ii,ii-1) = -H1
               enddo
  
               do ii=1,Ncoupure
                  nni(ii) = DD(Ncoupure-ii+1)
                  IF( nni(ii) .GT. 0.001D0) THEN
                      COMPTEUR = COMPTEUR + 1
                      NORM = 0.
                      DO JJ = 1,Ncoupure
                         VECT_P(JJ,COMPTEUR) = CMPLX(VR(JJ,
     &                                          Ncoupure-ii+1)
     &                                ,VI(JJ,Ncoupure-ii+1))     
                  
                         NORM = NORM + CABS(VECT_P(JJ,COMPTEUR))**2     
                         LL_RO(COMPTEUR)     = LLL
                         NI_RO(COMPTEUR)     = NNI(II)
                      ENDDO
c                         PRINT *,'ni,LL', NI_RO(COMPTEUR)   
c     &                                  ,LL_RO(COMPTEUR)      
                      NORM = 1./sqrt(NORM*pas)
                      DO JJ = 1,Ncoupure
                         XX = pas*float(JJ)
                         VECT_P(JJ,COMPTEUR) = 
     &                          VECT_P(JJ,COMPTEUR)*NORM
                      ENDDO
                      E_RO(COMPTEUR) = 0.
                      DO JJ = 1,Ncoupure
                         DO KK = 1,Ncoupure
                            E_RO(COMPTEUR) = E_RO(COMPTEUR) + 
     &                          CONJG(VECT_P(KK,COMPTEUR))
     &                          *Hamilt(KK,JJ)*VECT_P(JJ,COMPTEUR)
                         ENDDO
                      ENDDO
                      E_RO(COMPTEUR) = E_RO(COMPTEUR)*pas
c                      print *,'E_ro',E_RO(COMPTEUR)

                  ENDIF
               enddo
           ENDIF

c
c    Fin de la boucle sur les L
c
       ENDDO
       call ORDO_RO(COMPTEUR,N0,E_RO,NI_RO,LL_RO)    

       open(unit=55,file='NI_ERO.dat',status='unknown')

       do JJ=1,COMPTEUR
            write(55,*) E_RO(JJ),NI_RO(JJ),LL_RO(JJ)
       enddo
       CLOSE(55)

       END

       subroutine ORDO_RO(COMPTEUR,N0,E_RO,NI_RO,LL_RO) 
   
       IMPLICIT NONE

       INTEGER N0,COMPTEUR,llmax
       INTEGER LL_RO(N0),K,J,K_TEMP
       REAL*8  E_RO(N0),NI_RO(N0)
       REAL*8  E_ORD(COMPTEUR),NI_ORD(COMPTEUR),LL_ORD(COMPTEUR)
       REAL*8  E

       J = 0

       do while (J.LT.COMPTEUR)
          J = J + 1
          E = 1.0e12
          DO K = 1,COMPTEUR
             if(E .GT. E_RO(K)) THEN
                E = E_RO(K)
                K_TEMP = K
             ENDIF
          ENDDO
          E_ORD(J) = E_RO(K_TEMP)
          LL_ORD(J) = LL_RO(K_TEMP)
          NI_ORD(J) = NI_RO(K_TEMP)
          E_RO(K_TEMP) = 1.0e13
       ENDDO

       DO J = 1,COMPTEUR
          E_RO(J) = E_ORD(J)
          LL_RO(J) = LL_ORD(J)
          NI_RO(J) = NI_ORD(J)
          print *,E_ORD(J),J
       enddo

       end
