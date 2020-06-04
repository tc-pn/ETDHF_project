ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     program of Resolution of two coupled particles
c     1."Exact Resolution"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program principal

      IMPLICIT NONE

      integer     noeud,Norbitalmax,Norbital
      real*8      pas,K,t0,delta_t
c      PARAMETER   (pas=0.1,noeud=3*(2**8-1),Norbitalmax=1)
      PARAMETER   (pas=0.3,noeud=2**6-1,Norbitalmax=50)
c      PARAMETER   (pas=0.149706457,noeud=2**9-1,Norbitalmax=1)
      COMPLEX*16  WAVE(0:Noeud,Norbitalmax)
      real*8      POT(0:Noeud),densite(0:Noeud),NI(Norbitalmax) 
      integer     I_TEMPSMAX   
 
      call initialisation2(Norbitalmax,wave,pas,pot,K,t0,delta_t
     &                    ,I_TEMPSMAX,ni,Norbital)   
                               
      call evolution(wave,pot,I_TEMPSMAX,delta_t
     &               ,t0,densite,ni,pas,Norbital)  
                            

      end 

      include 'evolution_etdhf.for'

      include 'initial_c.for'

      include 'collision.for'

      include 'fourier.for'

