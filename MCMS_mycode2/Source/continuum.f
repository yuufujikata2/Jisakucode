      subroutine continuum(vmtz,imvz,kex,es,wnwj,ntas,refl,chi,lmoutx,
     &                     myid, numproces,allkex)

C
      INCLUDE 'cont.inc'
C
C
C
      COMMON/PARAM/VCON,XE,EV,E,IOUT,NOUT,NAT,
     1 NDAT,NSPINS,NAS,NT1,NTB,
     2 RS(AT_),XV(AT_),YV(AT_),ZV(AT_),Z(AT_),
     3 EXFACT(AT_),LMAXX(AT_),NZ(AT_),NSYMBL(AT_),
     4 NEQ(AT_),LCORE(AT_),KION,NAME(2),MLEQ(AT_)
      COMPLEX VCON,XE,EV

C
      COMMON /STATE/ CN(B_,N_),MN(B_,N_),
     1 IN(B_,N_),NATOM(B_,N_),LN(N_),
     2 NMS(N_),IMIN(N_,AT_),IMAX(N_,AT_),
     3 NLEQ(AT_),KTAU(AT_),NNS,ICORE,
     4 NUATOM,NDG,NLS(AT_),N0L(AT_),N0(AT_),
     5 NTERMS(AT_),LMAXN(AT_),NDIM,NDIMTR
      INTEGER*2 MN,IN,NATOM,IMIN,IMAX
C
C

c   peter's input variables    ( kex = #of energy-points <= NP_ )
      integer kex, lmoutx
      real vmtz, imvz, ES(NP_), wnwj(0:LMAX_,D_,NP_)
c
c   peter's output variables
      complex refl(NB_,NB_,NP_),chi(NB_,NB_,NP_)

c   MPI variables
      integer myid,numproces,allkex
      character chmyid
             
      write(chmyid,207) myid       
  207  format(I0)
       OPEN(4,FILE='sym'//chmyid//'.dat')
       OPEN(3,FILE='instr.dat')
       OPEN(14,FILE='fssym'//chmyid//'.out')
C
C***DO NOT CHANGE NNS=1. THIS VERSION USES SPIN INDEPENDENT POTENTIALS
      NNS=1
C
C***SET UP SYMMETRY INFORMATION FOR INITIAL STATE
C
      NAS=1
      CALL INPUT
C
      CALL SETUP

C
      CALL CONT(vmtz,imvz,kex,ES,wnwj,ntas,refl,chi,lmoutx,myid,
     &          numproces,allkex) 
C
      WRITE (6,160) 
  160 FORMAT(1X,'CONTINUUM MAIN EXIT ')
       CLOSE (4)
       CLOSE (3)
       CLOSE (14)

c     STOP
      END
C
C
C
      SUBROUTINE CONT(vmtz,imvz,kex,ES,wnwj,ntas,refl,chi,lmoutx,myid,
     &                numproces,allkex)
      INCLUDE 'cont.inc'
C
C
      COMMON/BESSEL/SBF(LTOT_),DSBF(LTOT_),SNF(LTOT_),DSNF(LTOT_)
      COMPLEX SBF,DSBF,SNF,DSNF
C
C
      COMMON/PARAM/VCON,XE,EV,E,IOUT,NOUT,NAT,
     1 NDAT,NSPINS,NAS,NT1,NTB,
     2 RS(AT_),XV(AT_),YV(AT_),ZV(AT_),Z(AT_),
     3 EXFACT(AT_),LMAXX(AT_),NZ(AT_),NSYMBL(AT_),
     4 NEQ(AT_),LCORE(AT_),KION,NAME(2),MLEQ(AT_)
      COMPLEX VCON,XE,EV
C
C
      COMMON /SECULR/ A(NTR_,NTR_),B(NB_,NB_),AE0(NTR_,NB_)
      COMPLEX*16 A,B,AE0
C
      COMMON /STATE/ CN(B_,N_),MN(B_,N_),
     1 IN(B_,N_),NATOM(B_,N_),LN(N_),
     2 NMS(N_),IMIN(N_,AT_),IMAX(N_,AT_),
     3 NLEQ(AT_),KTAU(AT_),NNS,ICORE,
     4 NUATOM,NDG,NLS(AT_),N0L(AT_),N0(AT_),
     5 NTERMS(AT_),LMAXN(AT_),NDIM,NDIMTR
      INTEGER*2 MN,IN,NATOM,IMIN,IMAX
C
      REAL*8 WKAREA
      COMPLEX CU
C
C.... DECLARATIONS FOR NAG SUBROUTINES
C
      CHARACTER UPLO
      INTEGER IPIV(NTR_)
      COMPLEX*16 WORK(LWORK)
      DATA UPLO/'L'/
C
c   peter's input variables
      integer kex, lmoutx
      real vmtz, imvz, ES(NP_), wnwj(0:LMAX_,D_,NP_)
c     
c   peter's output variables
      complex refl(NB_,NB_,NP_),chi(NB_,NB_,NP_)

c   MPI variables
      integer myid,numproces
      integer allkex
      integer numproces2
c
c   peter's local variables
      integer ndime,i,j,lm1,lm2
      complex dum
      complex*16 dum1,dum2,ucr(NB_,NB_),urc(NB_,NB_)
      complex*16 dumrefl(NB_,NB_), dumchi(NB_,NB_)
C....
      DATA THIRD/.3333333/, PAI/3.1415927/
C
C
cp    do iep=1,KEX  
cp       do iat=1,D_
cp          do il=0,LMAX_
cp             write(20,'("l=",I2," at=",I2," ep=",I3," wnwj=",F20.6)')
cp   &              il,iat,iep,wnwj(il,iat,iep)
cp          enddo
cp       enddo
cp    enddo
c
c  construct trafo matrices between real <-> cplx spherical harmonics bases
c
      call ucrurc(NB_,ucr,urc)
c
c

C
      CU = (0.0,1.0)
      PAI4=4.*PAI
      NTAS=NTERMS(NAS)
      NTB = NTAS
      NT1 = 0
C
C
Cp    DO 9 NE = 1, KEX, KP
Cp    ES(NE)=EMIN+(NE-1)*DE
C
      VCON = vmtz + (0.,1.) * imvz
      numproces2 = numproces - 1
c     for MPI    
      if ( myid .eq. numproces2 ) then         
      DO 9 NE = kex * myid + 1, allkex
C
      E=ES(NE)
C
      EV=E-VCON
C
C      IF(REAL(EV).LE.0.0) THEN
C      WRITE(6,137) E , VCON
C  137 FORMAT(1X, ' E-VCONR LES THAN ZERO FOR E =',F10.5,' AND VCON =',
C     1             2E15.7)
C        CALL EXIT 
C      ENDIF 
C
      XE=CSQRT(EV)
cp
cp      print*,XE, REAL(XE),AIMAG(XE)
      if ( AIMAG(XE).gt.0.) then
         XE = CMPLX( REAL(XE) , -AIMAG(XE) ) 
      endif
C
C CONSTRUCT SCATTERING MATRIX A(N,N) 
C
    5 CALL SMTX(wnwj(0,1,NE))
C


c  define parts of A=tau^-1 matrix: absorber block B,
c  and environment-absorber block AE0 
c  then shift back indices of A by
c  NTAS := NTERMS(NAS) = dim(absorber)
c
      do j=1,NTAS
       do i=1,NTAS
        B(i,j)=A(i,j)
       enddo
       do i=NTAS+1,NDIMTR
        AE0(i-NTAS,j)=A(i,j)
       enddo
      enddo
      do i=NTAS+1,NDIMTR
       do j=NTAS+1,NDIMTR
        A(i-NTAS,j-NTAS)=A(i,j)
       enddo
      enddo
      ndime=NDIMTR-NTAS

c      print*, "EV=", EV, " XE=", XE
c        print*
c        print*,"before inversion"
c        do i=1,5
c           write(*,*) ( A(i,j), j=1,5 )
c        enddo

c
c   calculate inverse of environment A-matrix, i.e. its tau matrix
c
      CALL ZSYTRF(UPLO,ndime,A,NTR_,IPIV,WORK,LWORK,INFO)
c     WRITE(6,60) WORK(1),INFO
  60  FORMAT(5X,' OPTIMAL LWORK = (',D14.7,',',D14.7,')','INFO=',I5)
      IF (INFO.NE.0) THEN
         PRINT *,'  THE FACTOR D IS SINGULAR'
         WRITE(6,61)
  61     FORMAT ('  THE FACTOR D IS SINGULAR')
         STOP
      ENDIF
C
      CALL ZSYTRI(UPLO,ndime,A,NTR_,IPIV,WORK,INFO)
c    
c   complete A^-1 matrix (zsytri only returns UPLO er part)
c
      do i=2,ndime
       do j=1,i-1
        if (UPLO.eq.'L'.or.UPLO.eq.'l') then
          A(j,i)=A(i,j)
         elseif (UPLO.eq.'U'.or.UPLO.eq.'u') then
          A(i,j)=A(j,i)
         else
          write(*,*)'UPLO not in L,l,U,u'
          stop
        endif
       enddo
      enddo
c        print*
c        print*,"after inversion"
c        do i=1,5
c           write(*,*) ( A(i,j), j=1,5 )
c        enddo
c   
c   now calculate reflectivity
c
c     print*,"A00=",B(1,1)   
c     print*,"AE0(i,1)=",(AE0(i,1), i = 1,5)
      do lm1=1,NTAS
       do lm2=1,NTAS
        dum1=(0.d0,0.d0)
        dum2=(0.d0,0.d0)
        do i=1,ndime
         do j=1,ndime
          dum1=dum1+AE0(i,lm1)*A(i,j)*AE0(j,lm2)
          dum2=dum2+DIMAG(AE0(i,lm1))*A(i,j)*AE0(j,lm2)
         enddo
        enddo
        refl(lm1,lm2,NE)=dum1
        chi(lm1,lm2,NE)=dum2
c       print out reflectivity of environment in REAL S.H. BASIS !
        IF (IOUT.EQ.5) THEN
          write(20,754)'lm1=',lm1,'lm2=',lm2,
     &      'refl=',refl(lm1,lm2,NE),'chi=',chi(lm1,lm2,NE)
        ENDIF
       enddo
      enddo
c 
c  now atomic t-matrix minus reflectivity
c
      do lm1=1,NTAS
       do lm2=1,NTAS
        B(lm1,lm2)=B(lm1,lm2)-refl(lm1,lm2,NE)
       enddo
      enddo
 754  format(1x,2(a4,i2,1x),2(2x,a5,2f10.7))
 755  format(3x,a19,2f10.7)
c
c  TRANSFORM REFLECTIVITY and chi FROM REAL TO COMPLEX representation 
c  using   O_cmplx(i,j) = \sum_{k,l} ucr(i,k) O_real(k,l) urc(l,j)  
c  if  (lmoutx==1), i.e. |_lm_> on _out_put uses comple_x_ spherical harmonics 
c
      if (lmoutx.eq.1) then
         do i=1,ntas
            do j=1,ntas
               dum1 = (0.,0.)
               dum2 = (0.,0.)
               do k=1,ntas
                  do l=1,ntas
                     dum1 = dum1 + ucr(i,k) * refl(k,l,NE) * urc(l,j)
                     dum2 = dum2 + ucr(i,k) *  chi(k,l,NE) * urc(l,j)
                  enddo
               enddo
               dumrefl(i,j) = dum1
               dumchi(i,j)  = dum2
            enddo
         enddo
         do i=1,ntas
            do j=1,ntas
               refl(i,j,NE) = dumrefl(i,j)
               chi(i,j,NE)  = dumchi(i,j)
            enddo
         enddo            
      endif 
C                
C
      CALL ZSYTRF(UPLO,NTAS,B,NB_,IPIV,WORK,LWORK,INFO) 
cp    WRITE(6,60) WORK(1),INFO
      IF (INFO.NE.0) THEN
         PRINT *,'  THE FACTOR D IS SINGULAR in B'
         STOP
      ENDIF
c
c
      CALL ZSYTRI(UPLO,NTAS,B,NB_,IPIV,WORK,INFO)
C
c   complete B^-1 matrix (zsytri only returns UPLO er part)
c
      do i=2,NTAS
       do j=1,i-1
        if (UPLO.eq.'L'.or.UPLO.eq.'l') then
          B(j,i)=B(i,j)
         elseif (UPLO.eq.'U'.or.UPLO.eq.'u') then
          B(i,j)=B(j,i)
         else
          write(*,*)'UPLO not in L,l,U,u'
          stop
        endif
       enddo
      enddo
c
c
      do lm1=1,NTAS
       do lm2=1,NTAS
        B(lm1,lm2)=B(lm1,lm2)*XE/PAI
       enddo
      enddo
c
C***WRITE OUT INVERSE OF MS MATRIX B(NTB,NTB) IF IOUT=5
C***AND RADIAL MATRIX ELEMENTS D(MG,J,II), J=NSTART,NLAST, II=1,NMI
C
C       IF(IOUT.EQ.5) THEN
C       DO 5739 NN=NSTART,NLAST
C       DO 5739 NM=1,NTB
C       WRITE (6,2821) NM,NN,B(NM,NN)
C 2821  FORMAT(' NM=',I3,' NN=',I3,' B(NM,NN):',2E14.6)
C 5739  CONTINUE
C       ENDIF
c
c  peter: other format to print out B-matrix
c     write(6,*)'Real(B(i,j))'
c     do i=1,NTAS
c       write(6,753)(REAL(B(i,j)),j=1,min(9,NTAS))
c     enddo
c     write(6,*)'Imag(B(i,j))'
c     do i=1,NTAS
c       write(6,753)(IMAG(B(i,j)),j=1,min(9,NTAS))
c     enddo
c753  format(f8.5,8(1x,f8.5))
C
    9 CONTINUE

      else 
      DO 11 NE = kex * myid + 1, kex * (myid + 1 )
C
      E=ES(NE)
C
      EV=E-VCON
C
C      IF(REAL(EV).LE.0.0) THEN
C      WRITE(6,137) E , VCON
C  137 FORMAT(1X, ' E-VCONR LES THAN ZERO FOR E =',F10.5,' AND VCON =',
C     1             2E15.7)
C        CALL EXIT 
C      ENDIF 
C
      XE=CSQRT(EV)
cp
cp      print*,XE, REAL(XE),AIMAG(XE)
      if ( AIMAG(XE).gt.0.) then
         XE = CMPLX( REAL(XE) , -AIMAG(XE) ) 
      endif
C
C CONSTRUCT SCATTERING MATRIX A(N,N) 
C
  10  CALL SMTX(wnwj(0,1,NE))
C


c  define parts of A=tau^-1 matrix: absorber block B,
c  and environment-absorber block AE0 
c  then shift back indices of A by
c  NTAS := NTERMS(NAS) = dim(absorber)
c
      do j=1,NTAS
       do i=1,NTAS
        B(i,j)=A(i,j)
       enddo
       do i=NTAS+1,NDIMTR
        AE0(i-NTAS,j)=A(i,j)
       enddo
      enddo
      do i=NTAS+1,NDIMTR
       do j=NTAS+1,NDIMTR
        A(i-NTAS,j-NTAS)=A(i,j)
       enddo
      enddo
      ndime=NDIMTR-NTAS

c      print*, "EV=", EV, " XE=", XE
c        print*
c        print*,"before inversion"
c        do i=1,5
c           write(*,*) ( A(i,j), j=1,5 )
c        enddo

c
c   calculate inverse of environment A-matrix, i.e. its tau matrix
c
      CALL ZSYTRF(UPLO,ndime,A,NTR_,IPIV,WORK,LWORK,INFO)
c     WRITE(6,60) WORK(1),INFO
      IF (INFO.NE.0) THEN
         PRINT *,'  THE FACTOR D IS SINGULAR'
         WRITE(6,61)
         STOP
      ENDIF
C
      CALL ZSYTRI(UPLO,ndime,A,NTR_,IPIV,WORK,INFO)
c    
c   complete A^-1 matrix (zsytri only returns UPLO er part)
c
      do i=2,ndime
       do j=1,i-1
        if (UPLO.eq.'L'.or.UPLO.eq.'l') then
          A(j,i)=A(i,j)
         elseif (UPLO.eq.'U'.or.UPLO.eq.'u') then
          A(i,j)=A(j,i)
         else
          write(*,*)'UPLO not in L,l,U,u'
          stop
        endif
       enddo
      enddo
c        print*
c        print*,"after inversion"
c        do i=1,5
c           write(*,*) ( A(i,j), j=1,5 )
c        enddo
c   
c   now calculate reflectivity
c
c     print*,"A00=",B(1,1)   
c     print*,"AE0(i,1)=",(AE0(i,1), i = 1,5)
      do lm1=1,NTAS
       do lm2=1,NTAS
        dum1=(0.d0,0.d0)
        dum2=(0.d0,0.d0)
        do i=1,ndime
         do j=1,ndime
          dum1=dum1+AE0(i,lm1)*A(i,j)*AE0(j,lm2)
          dum2=dum2+DIMAG(AE0(i,lm1))*A(i,j)*AE0(j,lm2)
         enddo
        enddo
        refl(lm1,lm2,NE)=dum1
        chi(lm1,lm2,NE)=dum2
c       print out reflectivity of environment in REAL S.H. BASIS !
        IF (IOUT.EQ.10) THEN
          write(20,754)'lm1=',lm1,'lm2=',lm2,
     &      'refl=',refl(lm1,lm2,NE),'chi=',chi(lm1,lm2,NE)
        ENDIF
       enddo
      enddo
c 
c  now atomic t-matrix minus reflectivity
c
      do lm1=1,NTAS
       do lm2=1,NTAS
        B(lm1,lm2)=B(lm1,lm2)-refl(lm1,lm2,NE)
       enddo
      enddo
c
c  TRANSFORM REFLECTIVITY and chi FROM REAL TO COMPLEX representation 
c  using   O_cmplx(i,j) = \sum_{k,l} ucr(i,k) O_real(k,l) urc(l,j)  
c  if  (lmoutx==1), i.e. |_lm_> on _out_put uses comple_x_ spherical harmonics 
c
      if (lmoutx.eq.1) then
         do i=1,ntas
            do j=1,ntas
               dum1 = (0.,0.)
               dum2 = (0.,0.)
               do k=1,ntas
                  do l=1,ntas
                     dum1 = dum1 + ucr(i,k) * refl(k,l,NE) * urc(l,j)
                     dum2 = dum2 + ucr(i,k) *  chi(k,l,NE) * urc(l,j)
                  enddo
               enddo
               dumrefl(i,j) = dum1
               dumchi(i,j)  = dum2
            enddo
         enddo
         do i=1,ntas
            do j=1,ntas
               refl(i,j,NE) = dumrefl(i,j)
               chi(i,j,NE)  = dumchi(i,j)
            enddo
         enddo            
      endif 
C                
C
      CALL ZSYTRF(UPLO,NTAS,B,NB_,IPIV,WORK,LWORK,INFO) 
cp    WRITE(6,60) WORK(1),INFO
      IF (INFO.NE.0) THEN
         PRINT *,'  THE FACTOR D IS SINGULAR in B'
         STOP
      ENDIF
c
c
      CALL ZSYTRI(UPLO,NTAS,B,NB_,IPIV,WORK,INFO)
C
c   complete B^-1 matrix (zsytri only returns UPLO er part)
c
      do i=2,NTAS
       do j=1,i-1
        if (UPLO.eq.'L'.or.UPLO.eq.'l') then
          B(j,i)=B(i,j)
         elseif (UPLO.eq.'U'.or.UPLO.eq.'u') then
          B(i,j)=B(j,i)
         else
          write(*,*)'UPLO not in L,l,U,u'
          stop
        endif
       enddo
      enddo
c
c
      do lm1=1,NTAS
       do lm2=1,NTAS
        B(lm1,lm2)=B(lm1,lm2)*XE/PAI
       enddo
      enddo
c
C***WRITE OUT INVERSE OF MS MATRIX B(NTB,NTB) IF IOUT=5
C***AND RADIAL MATRIX ELEMENTS D(MG,J,II), J=NSTART,NLAST, II=1,NMI
C
C       IF(IOUT.EQ.5) THEN
C       DO 5739 NN=NSTART,NLAST
C       DO 5739 NM=1,NTB
C       WRITE (6,2821) NM,NN,B(NM,NN)
C 2821  FORMAT(' NM=',I3,' NN=',I3,' B(NM,NN):',2E14.6)
C 5739  CONTINUE
C       ENDIF
c
c  peter: other format to print out B-matrix
c     write(6,*)'Real(B(i,j))'
c     do i=1,NTAS
c       write(6,753)(REAL(B(i,j)),j=1,min(9,NTAS))
c     enddo
c     write(6,*)'Imag(B(i,j))'
c     do i=1,NTAS
c       write(6,753)(IMAG(B(i,j)),j=1,min(9,NTAS))
c     enddo
c753  format(f8.5,8(1x,f8.5))
C
  11   CONTINUE

      end if
c
C  transpose reflectivity matrix in lm1 <-> lm2  in order to
c  annulate the effect of differing array orders in Fortran and C
c      

      if (myid.eq. numproces2 )then
      do iep = kex * myid +1,allkex
         do lm1 = 1, ntas
            do lm2 = 1, lm1-1
               dum = refl(lm1,lm2,iep) 
               refl(lm1,lm2,iep) = refl(lm2,lm1,iep) 
               refl(lm2,lm1,iep) = dum 
               dum = chi(lm1,lm2,iep) 
               chi(lm1,lm2,iep) = chi(lm2,lm1,iep) 
               chi(lm2,lm1,iep) = dum
            enddo
         enddo
      enddo
      write(6,*)'refl and chi: lm1 <-> lm2 transposition done for',
     &           allkex - myid * kex,' energy points'

      else
      do iep = kex * myid +1, kex * (myid + 1)
         do lm1 = 1, ntas
            do lm2 = 1, lm1-1
               dum = refl(lm1,lm2,iep) 
               refl(lm1,lm2,iep) = refl(lm2,lm1,iep) 
               refl(lm2,lm1,iep) = dum 
               dum = chi(lm1,lm2,iep) 
               chi(lm1,lm2,iep) = chi(lm2,lm1,iep) 
               chi(lm2,lm1,iep) = dum
            enddo
         enddo
      enddo
      write(6,*)'refl and chi: lm1 <-> lm2 transposition done for',
     &           kex,' energy points'

      end if


C  printout (for check) refl and chi at 1st energy point
c     write(20,'("1st energy point E =",F12.6)') ES(1)
c     do lm1 = 1, ntas
c        do lm2 = 1, ntas
c           write(20,335) lm1,lm2,refl(lm1,lm2,1),chi(lm1,lm2,1) 
c        enddo
c     enddo
 335  format('lm1=',I2,' lm2=',I2,' refl=',2(1x,F12.6),
     &       ' chi=',2(1x,F12.6))
C
      RETURN
C
      END
C 
      COMPLEX FUNCTION GMAT(L1,M1,L2,M2,YL,SBF,I)
      INCLUDE 'cont.inc'
C
C     G-MATRIX FOR POLYATOMIC MOLECULES USING REAL SPHERICAL HARMONICS
C
C
      COMMON/GAUNT/AI(WW_),INDEX(XX_)
C
      COMMON/PARAM/VCON,XE,EV,E
      COMPLEX XE,EV,VCON
C
      COMPLEX*16 GMATP
      COMPLEX*8 SBF
      DIMENSION SBF(LTOT_),YL(YL_)
C     TRUE DIMENSION YL(MYL),  MYL IS VARIABLE DEFINED IN SMTX
      LOGICAL MPHASE
      LOGICAL ENEG
      DATA SQR2 /1.414213562373/
      DATA PI4/12.56637061435916/
      DATA ZERO/0.0/
C
C
Corrected on 4nov07. Before:  
      ENEG=REAL(EV).LT.ZERO
C      ENEG=.FALSE.
      MM=IABS(M2-M1)
      LMIN=MAX0(IABS(L2-L1),MM)
      IF(MOD(LMIN+L2+L1,2).NE.0) LMIN=LMIN+1
      LMAX=L2+L1
      NP=MM+1+(LMIN*(LMIN+1))/2
      LD =2*LMIN+3
      IF(L2.GT.L1) GO TO 5
      MPHASE=.FALSE.
      LL=L1
      M=M1
      LP=L2
      MP=M2
      GO TO 6
    5 LL=L2
      M=M2
      LP=L1
      MP=M1
      MPHASE=.TRUE.
    6 IF(M.GE.0) GO TO 7
      M=-M
      MP=-MP
    7 ISUB=(LL*(LL+1)*(LL+2)*(3*LL+1))/24+((LL+1)*(LL+2)*M+LP*(LP+1))/2
     1  +IABS(MP)+1
      N=INDEX(ISUB)
      IF(MP.LT.0) N=N+MIN0(LP,(LL+LP-IABS(M+MP))/2)+1
      GMATP = (0.0D0,0.0D0) 
      NSGN=1
      LMIN1=LMIN+1
      LMAX1=LMAX+1
      DO 1 LP1=LMIN1,LMAX1,2
      L=LP1-1
      CLM=AI(N)
      N=N+1
      IF(ENEG) GO TO 2
      IF(NSGN.GT.0) GO TO 2
      CLM=-CLM
    2 GMATP = GMATP+SBF(L+1)*YL(NP)*CLM
      GMAT = GMATP
      NP=NP+LD
      LD=LD+4
    1 NSGN=-NSGN
      IF(MPHASE.AND.MOD(M+MP,2).NE.0) GMAT=-GMAT
      GMAT=GMAT*PI4
      IF(ENEG) GO TO 3
      IF(MOD(LMIN+L2-L1,4).NE.0) GMAT=-GMAT
      GO TO 4
    3 IF(MOD(L1,2).EQ.0) GMAT=-GMAT
4     IF (M1.EQ.0.OR.M2.EQ.0) GO TO 13
11    IF (M2.EQ.M1) GO TO 15
      GMAT=GMAT/SQR2
      GO TO 15
  13  GMAT=GMAT/(2.0,0.0)
   15 IF(M2.GE.M1) RETURN
      IF(MOD(MM,2).NE.0) GMAT=-GMAT
      IF(I.EQ.-1) GMAT=-GMAT
      RETURN
C                                                 
      END
C
      SUBROUTINE INPUT
C
      INCLUDE 'cont.inc'
C
C
      COMMON/PARAM/VCON,XE,EV,E,IOUT,NOUT,NAT,
     1 NDAT,NSPINS,NAS,NT1,NTB,
     2 RS(AT_),XV(AT_),YV(AT_),ZV(AT_),Z(AT_),
     3 EXFACT(AT_),LMAXX(AT_),NZ(AT_),NSYMBL(AT_),
     4 NEQ(AT_),LCORE(AT_),KION,NAME(2),MLEQ(AT_)
       COMPLEX VCON,XE,EV
C
C
      NOUT = 0
      NSPINS = 1
      INSTR = 3
      INSYM = 4
C
      READ(INSTR,*) NAT, NDAT, IOUT
C
      DO N=1,NAT
      READ (INSTR,*) NSYMBL(N),NEQ(N),NZ(N),XV(N),YV(N),ZV(N)
      Z(N)=NZ(N)
      ENDDO
c
cpeter  NEQ is lmax in instr.dat. In the end I set NEQ=0.
c
      NDIM = 0 
      DO N=1,NAT
        NDIM = NDIM + (NEQ(N)+1)*(NEQ(N)+1)
      ENDDO
      
      write(INSYM,'(2(I5))') NDIM, 2
      idim = 0
      do N = 1, NAT
        do l=0,NEQ(N)
          write(INSYM,'(3(I5))') l, 1, 0
          write(INSYM,'(3(I5),F15.10)') 0,1,N,1.
          idim = idim+1
          do m=1,l
            write(INSYM,'(3(I5))') l, 1, 0
            write(INSYM,'(3(I5),F15.10)') m,1,N,1.
            write(INSYM,'(3(I5))') l, 1, 0
            write(INSYM,'(3(I5),F15.10)') m,-1,N,1.
            idim = idim+2
          enddo
        enddo
      enddo
      if ( idim .ne. NDIM ) then
        print*,"PANIC: instr -> sym: idim=",idim,"!=",NDIM,"=NDIM"
        stop
      endif

      do N = 1, NAT
        NEQ(N) = 0
      enddo
      REWIND INSYM
C
      RETURN
      END
CC
C
      SUBROUTINE SETUP
      INCLUDE 'cont.inc'
C
      LOGICAL*4 PREV
C
      COMMON/GAUNT/AI(WW_),INDEX(XX_),YL(YL_),RAB(ZZ_)
C
      COMMON/PARAM/VCON,XE,EV,E,IOUT,NOUT,NAT,
     1 NDAT,NSPINS,NAS,NT1,NTB,
     2 RS(AT_),XV(AT_),YV(AT_),ZV(AT_),Z(AT_),
     3 EXFACT(AT_),LMAXX(AT_),NZ(AT_),NSYMBL(AT_),
     4 NEQ(AT_),LCORE(AT_),KION,NAME(2),MLEQ(AT_)
      COMPLEX VCON,XE,EV
C
      COMMON /STATE/ CN(B_,N_),MN(B_,N_),
     1 IN(B_,N_),NATOM(B_,N_),LN(N_),
     2 NMS(N_),IMIN(N_,AT_),IMAX(N_,AT_),
     3 NLEQ(AT_),KTAU(AT_),NNS,ICORE,
     4 NUATOM,NDG,NLS(AT_),N0L(AT_),N0(AT_),
     5 NTERMS(AT_),LMAXN(AT_),NDIM,NDIMTR
      INTEGER*2 MN,IN,NATOM,IMIN,IMAX
C
      LOGICAL*4 DOALL
      DATA SMALL,ZERO,ONE/1.E-5,0.0,1.0/
      DATA PI/3.14159265358979/,PI4/12.56637061435916/
C     INDEXX SHOULD BE THE SAME AS THE DIMENSION OF INDEX
      INTEGER*4 INDEXX/XX_/,JWW/WW_/,JYL/YL_/,NYL/0/
C
C
C
      PREV = .FALSE.
      INSYM = 4
      IOSYM = 14
C
      DOALL=.TRUE.
      GO TO 121
C
      ENTRY SYMM
      DOALL=.FALSE.
  121 IF(PREV ) GO TO 122
      DO 62 I=1,INDEXX
   62 INDEX(I)=0
      MAXSUB=0
      NCOEF=0
      MNYL = 0
      PREV=.TRUE.
  122 DO 59 N=1,NAT
      LMAXX(N)=0
      NLEQ(N)=0
      N0(N)=0
      N0L(N)=0
      LMAXN(N)=0
      NTERMS(N)=0
   59 NLS(N)=0
      NUATOM=0
      WRITE (6,327) INSYM,IOSYM
  327 FORMAT(' SYMMETRY INFORMATION READ IN FROM FILE',I3,/,
     X       ' SYMMETRY INFORMATION  WRITTEN  TO FILE',I3)
   60 CONTINUE
      READ (INSYM,1010) NDIM,NDG,NAME
1010  FORMAT (2I5,10X,10A4)
      WRITE (6,1011) NDG, NAME
1011  FORMAT (' DEGENERACY =',I3,'   REPRESENTATION =',10A4)
  101 FORMAT(16I5)
      IF (IOUT.GE.2) WRITE (IOSYM,103) NDIM,NDG,NAME
  103 FORMAT(' # BASIS FUNCTION INCLUDING O.S. =',I4,'   DEGENERACY=',
     1 I3,5X,2A4)
      WRITE(6,*) '******SYMMETRY WARNINGS********'
      DO 125 N=1,NDIM
      DO 131 NA=1,NAT
      IMIN(N,NA)=1
  131 IMAX(N,NA)=0
      READ (INSYM,101) LN(N),NMS(N),NONE
      IF (IOUT.GE.2) WRITE (IOSYM,104) N,LN(N),NMS(N),NONE
104   FORMAT ( 1X,'BASIS FUNCTION NO.',I3,'   L=',I3,'   NO. OF TERMS=',
     1         I3,10X,I3)
      NMN=NMS(N)
      IF (NONE.NE.0)
     $ READ (INSYM,1050) (MN(I,N),IN(I,N),NATOM(I,N),CN(I,N),I=1,NMN)
1050  FORMAT (4(2I2,I4,F12.9))
      DO 55 I=1,NMN
      IF (NONE.EQ.0)
     $ READ (INSYM,105) MN(I,N),IN(I,N),NATOM(I,N),CN(I,N)
  105 FORMAT(3I5,F15.10)
      IF (IOUT.GE.2) 
     * WRITE(IOSYM,106) MN(I,N),IN(I,N),NATOM(I,N),CN(I,N)
  106 FORMAT(30X,'M=',I3,'   I=',I3,'  ATOM NO.=',I3,'  CN=',F15.10)
      ININ=IN(I,N)
      IF (IABS(ININ).NE.1 .OR. MN(I,N).LT.0) WRITE(6,109)
      IF(MN(I,N).EQ.0.AND.IN(I,N).EQ.-1) WRITE(6,109)
  109 FORMAT('     WRONG VALUES FOR M AND I')
      NA=NATOM(I,N)
      IF(NLEQ(NA).EQ.0) NLEQ(NA)=NATOM(1,N)
      IF(NLEQ(NA).NE.NATOM(1,N)) WRITE(6,110)
  110 FORMAT(  '    INCONSISTENT NUMBERING OF ATOMS')
      IF(I.EQ.1) GO TO 58
      IF(NA.EQ.NATOM(I-1,N)) GO TO 56
      IF(NA.GT.NATOM(I-1,N)) GO TO 58
      WRITE(6,107)
  107 FORMAT('     SYMMETRY CARD OUT OF SEQUENCE')
      STOP
  58  IMIN(N,NA)=I
      I0=I
   56 IMAX(N,NA)=I
      IF(I0.EQ.I) GO TO 55
      I1=I-1
      DO 126 J=I0,I1
  126 IF(MN(I,N).EQ.MN(J,N).AND.IN(I,N).EQ.IN(J,N)) WRITE(6,111)
  111 FORMAT('    DUPLICATE TERMS IN BASIS FUNCTION   ')
   55 LMAXN(NA)=MAX0(LMAXN(NA),LN(N))
      NUATOM=MAX0(NUATOM,NLEQ(NA))
      NA=NATOM(1,N)
      NTERMS(NA)=NTERMS(NA)+1
      IF(N.EQ.1) GO TO 128
      IF(NA.LT.NATOM(1,N-1)) WRITE(6,107)
      IF(NA.NE.NATOM(1,N-1)) GO TO 128
      IF(LN(N).EQ.LN(N-1)) GO TO 125
      IF(LN(N).LT.LN(N-1)) WRITE(6,107)
  128 NLS(NA)=NLS(NA)+1
  125 CONTINUE
      NDIMTR=NDIM
      WRITE(6,999) NDIMTR
      IF (IOUT.GE.2) WRITE(IOSYM,999) NDIMTR
  999 FORMAT(' TRUE DIMENSION OF SECULAR MATRIX TO SOLVE =',I4)
      WRITE(6,112) NUATOM, NAME 
      WRITE(IOSYM,112) NUATOM, NAME 
  112 FORMAT(' NUMBER OF INEQUIVALENT ATOMS =',I4,
     *       ' FOR REPRESENTATION:',10A4) 
      N0(1)=1
      N0L(1)=1
      LMAXX(1)=MAX0(LMAXX(1),LMAXN(1))
      IF(NUATOM.EQ.1) GO TO 127
      DO 124 NA=2,NUATOM
      N0(NA)=N0(NA-1)+NTERMS(NA-1)
      N0L(NA)=N0L(NA-1)+NLS(NA-1)
  124 LMAXX(NA)=MAX0(LMAXN(NA),LMAXX(NA))
  127 DO 61 NN=1,NDIM
      NMN=NMS(NN)
      L=LN(NN)
      DO 63 I=1,NMN
      M=MN(I,NN)
      DO 64 NM=1,NDIM
      NMM=NMS(NM)
      LP=LN(NM)
      IF(L.LT.LP) GO TO 64
      LX=L+LP
      DO 65 J=1,NMM
      MP=MN(J,NM)
    7 ISUB=(L *(L +1)*(L +2)*(3*L +1))/24+((L +1)*(L +2)*M+LP*(LP+1))/2
     1  +MP+1
      IF (ISUB.LE.INDEXX) GO TO 113
      WRITE (6,102) ISUB
      WRITE (6,100) JWW,NCOEF,INDEXX,MAXSUB
      CALL MERR(151374)
  102 FORMAT('-JXX',I10,' CAN NOT CONTINUE'/'-INCOMPLETE')
  113 CONTINUE
      IF(INDEX(ISUB).NE.0) GO TO 65
   68 II2=1
      IF(MP.NE.0) II2=2
      INDEX(ISUB)=NCOEF+1
      MAXSUB=MAX0(MAXSUB,ISUB)
      DO 67 II=1,II2
      LMIN=MAX0(L-LP,IABS(M-MP))
      IF(MOD(LX-LMIN,2).NE.0) LMIN=LMIN+1
      IF(NCOEF+(LX-LMIN)/2 .LT. JWW) GO TO 3
      NCOEF = NCOEF + (LX-LMIN)/2+1
      GO TO 67
   3  DO 4 LL=LMIN,LX,2
      NCOEF=NCOEF+1
      ARG=(2*LL+1)*(2*L+1)/(PI4*(2*LP+1))
    4 AI(NCOEF)=CGC(LL,L,LP,0,0)*CGC(LL,L,LP,MP-M,M)*SQRT(ARG)
C     MODIFIED FOR DOUBLE PRECISION COMPILATION JUNE 27 1975
   67 MP=-MP
   65 CONTINUE
   64 CONTINUE
  63  CONTINUE
   61 CONTINUE
      WRITE (6,100) JWW,NCOEF,INDEXX,MAXSUB
      IF (IOUT.GE.2) WRITE (IOSYM,100) JWW,NCOEF,INDEXX,MAXSUB
  100 FORMAT(' DIMENSION JWW =',I6,' COULD BE',I6/
     1       ' DIMENSION JXX =',I6,' COULD BE',I6)
      IF (NCOEF.GT.JWW .OR. MAXSUB.GT.INDEXX) CALL MERR(151580)
      IF(.NOT.DOALL) RETURN
C
      ENTRY STRUCT
C      COMPUTE NUMBER OF ATOMS EQUIVALENT TO EACH DIFFERENT ATOM
      DO 130 NA=2,NUATOM
  130 KTAU(NA)=0
      KTAU(1)=1
      NYL = 1
      MLEQ(1) = NLEQ(1)
      DO 129 NA=2,NAT
      MLEQ(NA) = NLEQ(NA)
      IF(NLEQ(NA).EQ.0) GO TO 129
      KTAU(NLEQ(NA))=KTAU(NLEQ(NA))+1
      NA1=NA-1
      DO 415 NB=1,NA1
      IF(NLEQ(NB).EQ.0) GO TO 415
      NLAB=LMAXX(NLEQ(NA))+LMAXX(NLEQ(NB))
      NAB=((NA-1)*(NA-2))/2+NB
      MYL = (NLAB+1)*(NLAB+2)/2
      PHI=ZERO
      ZMU=ONE
      RAB(NAB)= (XV(NA)-XV(NB))**2+(YV(NA)-YV(NB))**2+(ZV(NA)-ZV(NB))**2
      IF(RAB(NAB).GT.SMALL) GO TO 42
      RAB(NAB)=ZERO
      GO TO 41
   42 RAB(NAB)=SQRT(RAB(NAB))
      ZMU=(ZV(NB)-ZV(NA))/RAB(NAB)
      RXY=      (XV(NA)-XV(NB))**2+(YV(NA)-YV(NB))**2
      IF(RXY.LT.SMALL) GO TO 41
      RXY=SQRT(RXY)
C before: PHI=ASIN((YV(NB)-YV(NA))/RXY)
      ARGXY=(YV(NB)-YV(NA))/RXY
      IF(ABS(ARGXY).GT.1.0) ARGXY=SIGN(1.0,ARGXY)
      PHI=ASIN(ARGXY)
      IF(XV(NB).GE.XV(NA)) GO TO 41
      PHI=PI-PHI
   41 IF((NYL+2*MYL-1).LE.JYL) CALL YLM1(NLAB,ZMU,PHI,YL(NYL),MYL)
      NYL = NYL+2*MYL
  415 CONTINUE
  129 CONTINUE
      NYL = NYL-1
      IF (MNYL.LT.NYL) MNYL=NYL
      WRITE (6,108) JYL,NYL
      IF (IOUT.GE.2) WRITE (IOSYM,108) JYL,NYL
  108 FORMAT(' DIMENSION JYL =',I8,' COULD BE',I8/1X)
      IF (NYL.GT.JYL) CALL MERR(152000)
      RETURN
      ENTRY INIT
      PREV=.FALSE.
      RETURN
C
      END
C
C
      SUBROUTINE SMTX(wnwj)
      INCLUDE 'cont.inc'
C
      COMMON/BESSEL/SBF(LTOT_),DSBF(LTOT_),SHF(LTOT_),DSHF(LTOT_)
      COMPLEX SBF,DSBF,SHF,DSHF
C
      COMMON/GAUNT/AI(WW_),INDEX(XX_),YL(YL_),RAB(ZZ_)
C
      COMMON/PARAM/VCON,XE,EV,E,IOUT,NOUT,NAT,
     1 NDAT,NSPINS,NAS,NT1,NTB,
     2 RS(AT_),XV(AT_),YV(AT_),ZV(AT_),Z(AT_),
     3 EXFACT(AT_),LMAXX(AT_),NZ(AT_),NSYMBL(AT_),
     4 NEQ(AT_),LCORE(AT_),KION,NAME(2),MLEQ(AT_)
      COMPLEX VCON,EV,XE
C
      COMMON /SECULR/ A(NTR_,NTR_),B(NB_,NB_),AE0(NTR_,NB_)
C
      COMPLEX*16 A,STMAT,B,AE0
C
      COMMON /STATE/ CN(B_,N_),MN(B_,N_),
     1 IN(B_,N_),NATOM(B_,N_),LN(N_),
     2 NMS(N_),IMIN(N_,AT_),IMAX(N_,AT_),
     3 NLEQ(AT_),KTAU(AT_),NNS,ICORE,
     4 NUATOM,NDG,NLS(AT_),N0L(AT_),N0(AT_),
     5 NTERMS(AT_),LMAXN(AT_),NDIM,NDIMTR
      INTEGER*2 MN,IN,NATOM,IMIN,IMAX
C
C  peter's variables
      real wnwj(0:LMAX_,D_)
c
      COMPLEX CSQRT,ARG,GM,GM1,GMAT,XEPI,XEPISR
C
      COMPLEX SHFP,DSHFP
C
      DATA ZERO,ONE,TWO/0.0,1.0,2.0/
      DATA PI/3.14159265358979/
C
C      XE= CSQRT(EV)
      XEPI=XE/PI
      XEPISR=CSQRT(XEPI)
      NL=0
      NS=(NNS-1)*NDAT
C
C INITIALIZE A TO ZERO AND B TO RHS UNIT VECTORS (NOUT=0) 
      do I=1,NDIMTR
       do J=1,NDIMTR
        A(I,J)=(0.0,0.0)
       enddo
      enddo
C
C           
cp    write(6,*)'NUATOM=',NUATOM
      write(6,*)'E = ', E
C
      DO 60 NA=1,NUATOM
      STMATP=0.E0
      NS=NS+1
      IF(NLEQ(NA).EQ.0) GO TO 60
      MOUT=1
   25 NT0A=N0(NA)            
      NTXA=NT0A+NTERMS(NA)-1
      IF(NEQ(NA).NE.0) GO TO 50
      L=-1
      NLP=-1
      ARG=XE*RS(NA)
      ML=LMAXN(NA)+1
      IF (ML.LT.3) ML = 3
      CALL CSBF(ARG,XE,ML,SBF,DSBF)
      CALL CSHF2(ARG,XE,ML,SHF,DSHF)
C
   43 DO 45 NN=NT0A,NTXA
      NNN=NN-NT1
      IF(LN(NN).EQ.L ) GO TO 35
      L=LN(NN)
      NL=NL+1
      NP=NL
Cp
Cp    write(20,'("from TMAT: STMAT=",2(E12.6))') STMAT 
      if ( REAL(EV) .GE. 0. ) then
        STMAT = wnwj(L,NA) + (0.,1.)      
      else
        STMAT = wnwj(L,NA) + (1.,0.)
      endif
              
Cp    write(20,'("from wnwj: STMAT=",2(E12.6)," L=",I2," NA=",I3)')
Cp   &     STMAT, L, NA 
  
  35  A(NNN,NNN)=STMAT
  45  CONTINUE
C
      GOTO  60
   50 NN0=N0(NEQ(NA))-NT1
      DO 55 NN=NT0A,NTXA
      NMN=NN-NT1
      A(NMN,NMN)=A(NN0,NN0)
   55 NN0=NN0+1
   60 CONTINUE
C
      IF(IOUT.EQ.5) THEN 
      DO 70 N=1,NDIMTR
      WRITE(6,1002) N, A(N,N)
 1002 FORMAT(' NN=',I3,' A(NN,NN) =',1P2E16.6)
   70 CONTINUE
      ENDIF
C

      IF (IOUT.EQ.5) WRITE (6,100) (MLEQ(I),I=1,NAT)
  100 FORMAT(1H0,'SMTX.MLEQ ',30I4)
      IF (IOUT.EQ.5) WRITE (6,101) (LMAXX(I),I=1,NAT)
  101 FORMAT (1X,'SMTX.LMAXX',30I4)
      IF (IOUT.EQ.5) WRITE (6,102) (NLEQ(I),I=1,NAT)
  102 FORMAT (1X,'SMTX.NLEQ ',30I4)
      IF (IOUT.EQ.5) WRITE (6,103) (LMAXN(I),I=1,NAT)
  103 FORMAT (1X,'SMTX.LMAXN',30I4)
      IF (IOUT.EQ.5) WRITE (6,104) (NTERMS(I),I=1,NAT)
  104 FORMAT (1X,'SMTXNTERMS',30I4)
      IF (IOUT.EQ.5) WRITE (6,105) (N0(I),I=1,NAT)
  105 FORMAT (1X,'SMTX.N0   ',30I4)
C
C 
C CALCULATE 'OFF DIAGONAL' BLOCKS OF G-MATS:
C         FOR ATOMIC SPHERES:
C                      GMAT:ARE CORRECTLY SYMMETRIZED THEN ADDED TO A
C                      B:STORED FOR LATTER CONSTRUCTION OF B (AX=B)
C
      NYLC=1
      DO 340 NA=2,NAT
      IF (MLEQ(NA).EQ.0) GO TO 340
      NT0A=N0(NLEQ(NA))
      NTXA=NT0A+NTERMS(NLEQ(NA))-1
      NA1=NA-1
      DO 330 NB=1,NA1
      IF (MLEQ(NB).EQ.0) GO TO 330
      NLAB = LMAXX(MLEQ(NA))+LMAXX(MLEQ(NB))
      MYL = (NLAB+1)*(NLAB+2)/2
      IF (NLEQ(NA)*NLEQ(NB).EQ.0) GO TO 320
      MOUT=1
      IF(NB.EQ.1.AND.NOUT.EQ.0.AND.MORLEX.NE.0) MOUT=2
      IF(NB.EQ.1.AND.NOUT.EQ.1) MOUT=2
      NAB=((NA-1)*(NA-2))/2+NB
      NT0B=N0(NLEQ(NB))
      NTXB=NT0B+NTERMS(NLEQ(NB))-1
      NYLS = NYLC+MYL
      ARG=XE*RAB(NAB)
      MLAB=LMAXN(NA)+LMAXN(NB)+1
      IF (MLAB.LT.3) MLAB = 3
      CALL CSHF2(ARG,(0.0,0.0),MLAB,SBF,DSBF)
C
      IF (IOUT.EQ.5) WRITE (6,180) NA,NB,( SBF(I),I=1,MLAB)
  180 FORMAT (1H0,'SMTX. SHF ',2I4,(1P8E12.4))
      IF (IOUT.EQ.5) WRITE (6,185) NA,NB,(DSBF(I),I=1,MLAB)
  185 FORMAT ( 1X,'SMTX.DSHF ',2I4,(1P8E12.4))
C
      DO 310  NN=NT0A,NTXA
      NNT1=NN-NT1
      LA=LN(NN)
      IMINA=IMIN(NN,NA)
      IMAXA=IMAX(NN,NA)
      IF(IMINA.GT.IMAXA) GO TO 310
      DO 300 NM=NT0B,NTXB
      NMT1=NM-NT1
      LB=LN(NM)
      IMINB=IMIN(NM,NB)
      IMAXB=IMAX(NM,NB)
      IF(IMINB.GT.IMAXB) GO TO 300
      DO 260 I=IMINA,IMAXA
      MA=MN(I,NN)
      IA=IN(I,NN)
      DO 260 J=IMINB,IMAXB
      IF(ABS(CN(I,NN)).LT.1.E-5 .OR. ABS(CN(J,NM)).LT.1.E-5) GOTO 260
      MB=MN(J,NM)
      IB=IN(J,NM)
      IF(IA.NE.IB) GOTO  200
      GM = GMAT(LA,MA,LB,MB,YL(NYLC),SBF,1)
      IF(MA.EQ.0.OR.MB.EQ.0) GO TO 220
      GM1 = GMAT(LA,MA,LB,-MB,YL(NYLC),SBF,1)
      IF(IA.EQ.-1) GM1=-GM1
      GOTO  215
  200 IF(MA.NE.MB) GOTO  205
      GM=(0.0,0.0)
      GOTO  210
  205 GM = GMAT(LA,MA,LB,MB,YL(NYLS),SBF,-1)
      IF(IA.EQ.-1) GM=-GM
      IF(MA.EQ.0.OR.MB.EQ.0) GO TO 220
  210 GM1 = -GMAT(LA,MA,LB,-MB,YL(NYLS),SBF,-1)
  215 IF(MOD(MB,2).NE.0) GM1=-GM1
      GOTO  225
  220 GM1=GM
  225 GM=GM+GM1
C      IF(MOUT.EQ.2) GOTO  235
      IF(NM.GT.NN) GOTO  230
      IF(NM.EQ.NN) GM=(2.0,0.0)*GM
      A(NNT1,NMT1)=A(NNT1,NMT1)+      GM*CN(I,NN)*CN(J,NM)
      GOTO  260
  230 A(NMT1,NNT1)=A(NMT1,NNT1)+      GM*CN(I,NN)*CN(J,NM)
      GOTO  260
C  235 CONTINUE
C     IF(NM.LT.NN) GOTO  255
C      WRITE(6,245) 
C  245 FORMAT(' SMTX EXIT: ERROR IN STORAGE FOR GMAT ELEMENTS')
C      CALL EXIT
C  255 BL(NNT1,NM) = BL(NNT1,NM)+    GM*CN(I,NN)*CN(J,NM)
  260 CONTINUE
  300 CONTINUE
  310 CONTINUE
  320 NYLC = NYLC+2*MYL
  330 CONTINUE
  340 CONTINUE
C
C
       IF(IOUT.EQ.5) THEN
       DO 5739 NM=2,NDIMTR
       NA1=NM-1
       NMP1=NM+NT1
       DO 5739 NN=1,NA1
       NNP1=NN+NT1
       WRITE (6,282) NMP1,NNP1,A(NM,NN)
 282   FORMAT(' NMP1=',I3,' NNP1=',I3,' A(NM,NN):',2E14.6)
 5739  CONTINUE
       ENDIF
C
C FILL IN REST OF SCATTERING MATRIX A
C
      DO 490  NN=2,NDIMTR
      NA1=NN-1
      DO 490  NM=1,NA1
  490 A(NM,NN)=A(NN,NM)
      RETURN
C
      END
C

      SUBROUTINE CSBF(X,Y,MAX,SBF,DSBF)
      REAL*8 XF1
      COMPLEX X,Y,SBF,DSBF,RAT,DSBF1,CSIN,Z,SBFJ,B,A,CCOS,CABS
      COMPLEX*16 SBFK,SBF1,SBF2,CDEXP
      INTEGER MAX,K,JMIN,KMAX
      DIMENSION SBF(MAX), DSBF(MAX)
C
C
C     GENERATES SPHERICAL BESSEL FUNCTIONS OF ORDER 0 - MAX-1 AND THEIR
C     FIRST DERIVATIVES WITH RESPECT TO R.  X=ARGUMENT= Y*R.
C     IF Y=0, NO DERIVATIVES ARE CALCULATED.  MAX MUST BE AT LEAST 3.
C     OSBF GENERATES ORDINARY SPHERICAL BESSEL FUNCTIONS.  MSBF - MODI-
C     FIED SPHERICAL BESSEL FUNCTIONS; OSNF - ORD. SPH. NEUMANN FCNS;
C     MSNF - MOD. SPH. NEUMANN FCNS; MSHF - MOD. SPH HANKEL FCNS
C
C
C
    1 IF (MAX.LT.1.OR.MAX.GT.2000) GO TO 99
      ABSX = ABS(X)
      IF(ABSX.LT.0.50 ) GO TO 18
C
C
C     BESSEL FUNCTIONS BY DOWNWARD RECURSION
      SBF2=(0.0,0.0)
      SBF1=1.0E-25*(0.5,0.5)
      IF(ABSX.LT.2.0) SBF1=1.0E-37*(0.5,0.5)
      JMIN=10+ABS(X)
      KMAX=MAX+JMIN-1
      K=MAX
      XF1=2*KMAX+1
      DO  10  J=1,KMAX
      SBFK=XF1*SBF1/X-SBF2
      SBF2=SBF1
      SBF1=SBFK
      XF1=XF1-2.0D0
      IF (J.LT.JMIN) GO TO 10
      SBF(K)=SBFK
      K=K-1
10    CONTINUE           
      RAT=CSIN(X)/(X*SBF(1))
   16 DO 17 K=1,MAX
   17 SBF(K)=RAT*SBF(K)
      DSBF1=-SBF(2)
      GO TO 26
C
C
C
C     SMALL ARGUMENTS
   18 Z=-(X*X*0.50)
      A=(1.0,0.0)           
      MMX=MAX
      IF (MAX.EQ.1.AND.Y.NE.(0.0,0.0)) MMX=2
      DO  30  J=1,MMX
      SBFJ=A
      B=A
      DO 31 I=1,20
      B=B*Z/(I*(2*(J+I)-1))
      SBFJ=SBFJ+B
      ABSB = ABS(B)
      ABSJ = ABS(SBFJ)
      IF (ABSB.LE.1.0E-07*ABSJ) GO TO 29
   31 CONTINUE
29    IF (J.EQ.2) DSBF1=-SBFJ
      IF (J.LE.MAX) SBF(J)=SBFJ
   30 A=A*X/ FLOAT(2*J+1)
      GO TO 26
C
C
C ENTRY TO CALCULATE SPHERICAL NEUMANN FUNCTIONS 
      ENTRY CSNF(X,Y,MAX,SBF,DSBF)
      SBF2=-CCOS(X)/X
      IF (MAX.EQ.1 .AND. Y.EQ.(0.0,0.0))  GO TO 2
      SBF1=(SBF2-CSIN(X))/X
      DSBF1=-SBF1
      GO TO 2
C
C
C ENTRY TO CALCULATE SPHERICAL HANKEL FUNCTIONS OF FIRST TYPE ('OUTGOING')
C************  NOTE :  RETURNS I [-(0.0,1.0)] TIMES HL1 ******************
      ENTRY CSHF1(X,Y,MAX,SBF,DSBF)
      SBF2=(0.0,1.0)*X
      SBF2=-(0.0,1.0)*CDEXP(SBF2)/SBF2
      IF (MAX.EQ.1 .AND. Y.EQ.(0.0,0.0))  GO TO 2
      SBF1=SBF2*((1.0/X)-(0.0,1.0))
      DSBF1=-SBF1
      GOTO 2
C
C
C ENTRY TO CALCULATE SPHERICAL HANKEL FUNCTIONS OF SECOND TYPE ('INGOING')
C************  NOTE :  RETURNS I [(0.0,1.0)] TIMES HL2 *******************
      ENTRY CSHF2(X,Y,MAX,SBF,DSBF)
      SBF2=-(0.0,1.0)*X
      SBF2=(0.0,1.0)*CDEXP(SBF2)/SBF2
      SBF1=SBF2*((1.0/X)+(0.0,1.0))
      DSBF1=-SBF1
2     SBF(1)=SBF2
      IF (MAX.LT.1.OR.MAX.GT.2000) GO TO 99
      IF (MAX.EQ.1) GO TO 26
      SBF(2)=SBF1
      IF (MAX.EQ.2) GO TO 26
      XF1=3.0D0
21    DO  22  I=3,MAX
      SBFK=XF1*SBF1/X-SBF2
      SBF(I)=SBFK
      SBF2=SBF1
      SBF1=SBFK
22    XF1=XF1+2.0D0
26    IF (Y.EQ.(0.0,0.0))  RETURN
      DSBF(1)=Y*DSBF1
      IF (MAX.EQ.1)  RETURN
      DO 9 I=2,MAX
    9 DSBF(I)=Y*(SBF(I-1)- FLOAT(I)*SBF(I)/X)
      RETURN
99    WRITE(6,100) MAX
100   FORMAT ('       SPHERICAL BESSEL FUNCTION ROUTINE - MAX=',I8)
 
      STOP   
      END


      SUBROUTINE YLM1(LMAX,Z,PHI,YL,MYL)
C     GENERATES REAL SPHERICAL HARMONICS, L = 0 TO LMAX; M =  0 TO L.
C     ARRANGED WITH M VARYING MOST RAPIDLY.  YL (I,1)=EVEN SPHERICAL
C     HARMONIC; YL(I,2)=ODD SPHERICAL HARMONIC.
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 PI4/12.56637061435917D0/,PI2/6.283185307179586D0/
      REAL*8 Y0/.282094791773878D0/,ZERO/0.0D0/,ONE/1.0D0/
      DIMENSION YL(MYL,2)
      REAL*4 Z,PHI,YL
      YL(1,2)=ZERO
      YL(1,1)=Y0
      IF(LMAX.EQ.0) RETURN
      X=DSQRT(ONE-Z*Z)
      SINPHI= SIN(PHI)
      COSPHI= COS(PHI)
      SINMP=ZERO
      COSMP=ONE
      PMM=ONE
      FAC=ONE
      ISUB=1
      MFAC2=0
      MFAC=1
C     P(M,M) BY RECURSION FORMULA
      LP1=LMAX+1
      DO 1 MP1=1,LP1
      M=MP1-1
      IF(M.EQ.0) GO TO 10
      FAC=FAC/(MFAC*MFAC2)
      PMM=-PMM *MFAC*X
      COSTP=COSMP*COSPHI-SINMP*SINPHI
      SINMP=SINMP*COSPHI+COSMP*SINPHI
      COSMP=COSTP
      FAC1=FAC/PI2
      MFAC=MFAC+2
      YL1=PMM*DSQRT(FAC1*MFAC)
      YL(ISUB,1)=YL1*COSMP
      YL(ISUB,2)=YL1*SINMP
      IF(M.EQ.LMAX) RETURN
      GO TO 11
   10 FAC1=FAC/PI4
   11 ISUB1=ISUB+M+1
      ISUB=ISUB1+1
      LFAC=MFAC
      PLM=PMM
      MFAC2=MFAC2+2
      M1=M+1
C     RECURSION FOR P(L,M), L = M+1 TO LMAX.
      DO 2 L=M1,LMAX
      LM=L-M
      LP=L+M
      PLM1=PLM*Z*LFAC
      LFAC=LFAC+2
      IF(L.EQ.M1) GO TO 20
      PLM1=PLM1-(LP-1)*PLM0
   20 PLM1=PLM1/LM
      IF(M.EQ.0) GO TO 21
      FAC1=LM*FAC1/LP
   21 YL1=DSQRT(LFAC*FAC1)*PLM1
      YL(ISUB1,1)=YL1*COSMP
      YL(ISUB1,2)=YL1*SINMP
      PLM0=PLM
      PLM=PLM1
    2 ISUB1=ISUB1+L+1
    1 CONTINUE
      RETURN
C                                                                    #
      END
C
      SUBROUTINE MERR(ISEQ)
      WRITE(6,1) ISEQ
    1 FORMAT(' STOP CAUSED BY MERR AT SEQUENCY NUMBER',I6)
      STOP
      END
C
      FUNCTION CGC(L1,L2,L3,M1,M2)
      IMPLICIT REAL*8(A-H,O-Z)
C CLEBSCH-GORDAN COEFFICIENT  EQ. 3.18, ROSE
      DIMENSION NUM(5),ND(5)
      REAL*4 CGC
      NFF=0
      M3=M1+M2
C     ARGUMENTS OF FACTORIALS
      NUM(1)=L3+L1-L2
      NUM(2)=L3-L1+L2
      NUM(3)=L1+L2-L3
      NUM(4)=L3+M3
      NUM(5)=L3-M3
      ND(1)=L1+L2+L3+1
      ND(2)=L1-M1
      ND(3)=L1+M1
      ND(4)=L2-M2
      ND(5)=L2+M2
C     CHECK TRIANGLE AND PROJECTION CONDITIONS
      DO 12 I=1,5
      IF(NUM(I)) 99,11,11
   11 IF(ND(I)) 99,12,12
   12 CONTINUE
      FF=1.
C     TWO SETS OF FACTORIAL PRODUCTS
      N=5
      DO 120 NFAC=1,2
      N1=N-1
C     ARRANGE ARGUMENTS IN DESCENDING ORDER
      DO 13 I=1,N1
      INUM=I
      ID=I
      I1=I+1
      DO 14 J=I1,N
      IF(NUM(J).LE.NUM(INUM)) GO TO 15
      INUM=J
   15 IF(ND(J).LE.ND(ID)) GO TO 14
      ID=J
   14 CONTINUE
      NTEMP=NUM(I)
      NUM(I)=NUM(INUM)
      NUM(INUM)=NTEMP
      NTEMP=ND(I)
      ND(I)=ND(ID)
   13 ND(ID)=NTEMP
C     COMPUTE FACTORIAL RATIOS
      DO 16 I=1,N
      IF(NUM(I)-ND(I)) 17,16,18
   17 JM=ND(I)
      IF(JM.EQ.1) GO TO 16
      J0=NUM(I)+1
      IF(NUM(I).EQ.0) J0=2
      DO 19 J=J0,JM
      IF(DABS(FF).GT.1.D-20) GO TO 19
      FF=FF*1.D20
      NFF=NFF-2
   19 FF=FF/DFLOAT(J)
      GO TO 16
   18 JM=NUM(I)
      IF(JM.EQ.1) GO TO 16
      J0=ND(I)+1
      IF(ND(I).EQ.0) J0=2
      DO 20 J=J0,JM
      IF(DABS(FF).LT.1.D 20) GO TO 20
      FF=FF/1.D20
      NFF=NFF+2
   20 FF=FF*DFLOAT(J)
   16 CONTINUE
      IF(NFAC.EQ.2) GO TO 21
      NFF=NFF/2
      FF=DSQRT((2*L3+1)*FF)
C     SECOND SET OF FACTORIAL ARGUMENTS
      NMIN=MAX0(0,L2+M3-L1)
      NUM(1)=L2+L3+M1-NMIN
      NUM(2)=L1-M1+NMIN
      NUM(3)=0
      ND(1)=NMIN
      IF(NMIN.EQ.0) ND(1)=L1-L2-M3
      ND(2)=L3-L1+L2-NMIN
      ND(3)=L3+M3-NMIN
  120 N=3
   21 IF(MOD(NMIN+L2+M2,2).EQ.0) GO TO 22
      FF=-FF
   22 FF=FF*1.D10**NFF
      CGCP = FF
      NMAX=MIN0(L3-L1+L2,L3+M3)
      CGC = CGCP
      IF(NMIN.GE.NMAX) RETURN
      NMIN=NMIN+1
      DO 23 NU=NMIN,NMAX
      FF= -(((L1-M1+NU)*(L3-L1+L2-NU+1)*(L3+M3-NU+1))/DFLOAT(NU*(NU+L1-L
     1  2-M3)*(L2+L3+M1-NU+1)))*FF
   23 CGCP = CGCP+FF
      CGC = CGCP
      RETURN
   99 CGC=0.0
      RETURN
      END
C
      subroutine ucrurc(NB_,ucr,urc)
c
c   construct L_real <-> L_complex transformation matrices 
c   ucr(i,j) = <L_c(i)|L_r(j))  urc(i,j) = (L_r(i)|L_c(j)> = urc(j,i)*
c   relies on ordering: l=0,lmax;  |m|=0,l; tau=1,2 (1/2 = +/- = even/odd) 
c   The non-trivial blocks of ucr are: 
c   [ < m|m+) < m|m-)  ]  = [   1       -i     ] = diag{ 1 , (-)^m ) * ucr22 
c   [ <-m|m+) <-m|m-)  ]    [ (-)^m    i*(-)^m ] 
c   Here,  m > 0 and  |>, |) denotes cplx, real state, respectively
c   and urc is the hermitian conjugate of ucr 
c
      implicit none
      integer NB_,lmax,ll,mm,taur,tauc,mfac,i,j
      complex*16 ucr22(2,2),urc22(2,2),ucr(NB_,NB_),urc(NB_,NB_)
      data ucr22/ (0.7071067811865475244,0.d0), 
     &            (0.7071067811865475244,0.d0),
     &            (0.d0,-0.7071067811865475244), 
     &            (0.d0, 0.7071067811865475244)/
      data urc22/ (0.7071067811865475244,0.d0),
     &            (0.d0,0.7071067811865475244),
     &            (0.7071067811865475244,0.d0), 
     &            (0.d0,-0.7071067811865475244)/
c
      lmax = sqrt( NB_ + 0.1 ) - 1
c
      do i = 1, NB_
        do j = 1, NB_
           ucr(i,j) = (0.d0, 0.d0) 
           urc(i,j) = (0.d0, 0.d0)
        enddo
      enddo
c
      i = 0
      do ll = 0, lmax
         i = i + 1
         ucr(i,i) = (1.d0, 0.d0)
         urc(i,i) = (1.d0, 0.d0)
         do mm = 1, ll
            do taur = 1, 2
               do tauc = 1, 2
                  mfac = 1
                  if ( tauc.eq.2.and.mod(mm,2).eq.1 ) mfac = -1
                  ucr( i + tauc, i + taur ) = mfac * ucr22(tauc,taur)
                  urc( i + taur, i + tauc ) = mfac * urc22(taur,tauc)
               enddo
            enddo
            i = i + 2
         enddo
      enddo
c
c     write(*,*)'ucr == < L_cplx | L_real > '
c     write(*,'(16(7x,I2)),/)') ( INT(sqrt(j-.5)), j=1,9 )
c     do i=1,9
c        write(*,'(16(1x,2(F4.1)),/)') ( ucr(i,j), j=1,9 )
c     enddo
c     print*
c     write(*,*)'urc == < L_real | L_cplx >'
c     write(*,'(16(7x,I2)),/)') ( INT(sqrt(j-.5)), j=1,9 )
c     do i=1,9
c        write(*,'(16(1x,2(F4.1)),/)') ( urc(i,j), j=1,9 )
c     enddo
c
      return
      END
c       (of ucrurc)
