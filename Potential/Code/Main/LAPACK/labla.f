*> \brief \b DAXPY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION DA
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION DX(*),DY(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    DAXPY constant times a vector plus a vector.
*>    uses unrolled loops for increments equal to one.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DY(I) + DA*DX(I)
            END DO
         END IF
         IF (N.LT.4) RETURN
         MP1 = M + 1
         DO I = MP1,N,4
            DY(I) = DY(I) + DA*DX(I)
            DY(I+1) = DY(I+1) + DA*DX(I+1)
            DY(I+2) = DY(I+2) + DA*DX(I+2)
            DY(I+3) = DY(I+3) + DA*DX(I+3)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
         END DO
      END IF
      RETURN
      END
*> \brief \b DCABS1
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DCABS1(Z)
* 
*       .. Scalar Arguments ..
*       COMPLEX*16 Z
*       ..
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DCABS1 computes |Re(.)| + |Im(.)| of a double complex number 
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup double_blas_level1
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION DCABS1(Z)
*
*  -- Reference BLAS level1 routine (version 3.6.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      COMPLEX*16 Z
*     ..
*     ..
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,DIMAG
*
      DCABS1 = ABS(DBLE(Z)) + ABS(DIMAG(Z))
      RETURN
      END
*> \brief \b DCOPY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
* 
*       .. Scalar Arguments ..
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION DX(*),DY(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    DCOPY copies a vector, x, to a vector, y.
*>    uses unrolled loops for increments equal to one.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,7)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DX(I)
            END DO
            IF (N.LT.7) RETURN
         END IF   
         MP1 = M + 1
         DO I = MP1,N,7
            DY(I) = DX(I)
            DY(I+1) = DX(I+1)
            DY(I+2) = DX(I+2)
            DY(I+3) = DX(I+3)
            DY(I+4) = DX(I+4)
            DY(I+5) = DX(I+5)
            DY(I+6) = DX(I+6)
         END DO
      ELSE      
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DY(IY) = DX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
*> \brief \b DDOT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
* 
*       .. Scalar Arguments ..
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION DX(*),DY(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    DDOT forms the dot product of two vectors.
*>    uses unrolled loops for increments equal to one.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      DDOT = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DTEMP + DX(I)*DY(I)
            END DO
            IF (N.LT.5) THEN
               DDOT=DTEMP
            RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     $            DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DTEMP + DX(IX)*DY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      DDOT = DTEMP
      RETURN
      END
*> \brief \b DGEMM
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER K,LDA,LDB,LDC,M,N
*       CHARACTER TRANSA,TRANSB
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGEMM  performs one of the matrix-matrix operations
*>
*>    C := alpha*op( A )*op( B ) + beta*C,
*>
*> where  op( X ) is one of
*>
*>    op( X ) = X   or   op( X ) = X**T,
*>
*> alpha and beta are scalars, and A, B and C are matrices, with op( A )
*> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n',  op( A ) = A.
*>
*>              TRANSA = 'T' or 't',  op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c',  op( A ) = A**T.
*> \endverbatim
*>
*> \param[in] TRANSB
*> \verbatim
*>          TRANSB is CHARACTER*1
*>           On entry, TRANSB specifies the form of op( B ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSB = 'N' or 'n',  op( B ) = B.
*>
*>              TRANSB = 'T' or 't',  op( B ) = B**T.
*>
*>              TRANSB = 'C' or 'c',  op( B ) = B**T.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry,  M  specifies  the number  of rows  of the  matrix
*>           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry,  N  specifies the number  of columns of the matrix
*>           op( B ) and the number of columns of the matrix C. N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>           On entry,  K  specifies  the number of columns of the matrix
*>           op( A ) and the number of rows of the matrix op( B ). K must
*>           be at least  zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*>           part of the array  A  must contain the matrix  A,  otherwise
*>           the leading  k by m  part of the array  A  must contain  the
*>           matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*>           least  max( 1, k ).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*>           part of the array  B  must contain the matrix  B,  otherwise
*>           the leading  n by k  part of the array  B  must contain  the
*>           matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*>           least  max( 1, n ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION.
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*>           supplied as zero then C need not be set on input.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*>           Before entry, the leading  m by n  part of the array  C must
*>           contain the matrix  C,  except when  beta  is zero, in which
*>           case C need not be set on entry.
*>           On exit, the array  C  is overwritten by the  m by n  matrix
*>           ( alpha*op( A )*op( B ) + beta*C ).
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>           On entry, LDC specifies the first dimension of C as declared
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup double_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine (version 3.6.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL NOTA,NOTB
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      IF (NOTA) THEN
          NROWA = M
          NCOLA = K
      ELSE
          NROWA = K
          NCOLA = M
      END IF
      IF (NOTB) THEN
          NROWB = K
      ELSE
          NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND.
     +    (.NOT.LSAME(TRANSA,'T'))) THEN
          INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND.
     +         (.NOT.LSAME(TRANSB,'T'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
          INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 13
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGEMM ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
*
*     And if  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (NOTB) THEN
          IF (NOTA) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
              DO 90 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  END IF
                  DO 80 L = 1,K
                      TEMP = ALPHA*B(L,J)
                      DO 70 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B + beta*C
*
              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  100                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  110             CONTINUE
  120         CONTINUE
          END IF
      ELSE
          IF (NOTA) THEN
*
*           Form  C := alpha*A*B**T + beta*C
*
              DO 170 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 130 I = 1,M
                          C(I,J) = ZERO
  130                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 140 I = 1,M
                          C(I,J) = BETA*C(I,J)
  140                 CONTINUE
                  END IF
                  DO 160 L = 1,K
                      TEMP = ALPHA*B(J,L)
                      DO 150 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  150                 CONTINUE
  160             CONTINUE
  170         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B**T + beta*C
*
              DO 200 J = 1,N
                  DO 190 I = 1,M
                      TEMP = ZERO
                      DO 180 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  180                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  190             CONTINUE
  200         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DGEMM .
*
      END
*> \brief \b DGEMV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER INCX,INCY,LDA,M,N
*       CHARACTER TRANS
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGEMV  performs one of the matrix-vector operations
*>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
*>
*> where alpha and beta are scalars, x and y are vectors and A is an
*> m by n matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry, TRANS specifies the operation to be performed as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*>
*>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
*>
*>              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of the matrix A.
*>           M must be at least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*>           Before entry, the leading m by n part of the array A must
*>           contain the matrix of coefficients.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, m ).
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is DOUBLE PRECISION array of DIMENSION at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*>           and at least
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*>           Before entry, the incremented array X must contain the
*>           vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION.
*>           On entry, BETA specifies the scalar beta. When BETA is
*>           supplied as zero then Y need not be set on input.
*> \endverbatim
*>
*> \param[in,out] Y
*> \verbatim
*>          Y is DOUBLE PRECISION array of DIMENSION at least
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*>           and at least
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*>           Before entry with BETA non-zero, the incremented array Y
*>           must contain the vector y. On exit, Y is overwritten by the
*>           updated vector y.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           Y. INCY must not be zero.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup double_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
*  -- Reference BLAS level2 routine (version 3.6.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +    .NOT.LSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGEMV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
          JX = KX
          IF (INCY.EQ.1) THEN
              DO 60 J = 1,N
                  TEMP = ALPHA*X(JX)
                  DO 50 I = 1,M
                      Y(I) = Y(I) + TEMP*A(I,J)
   50             CONTINUE
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  TEMP = ALPHA*X(JX)
                  IY = KY
                  DO 70 I = 1,M
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
   70             CONTINUE
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
*
*        Form  y := alpha*A**T*x + y.
*
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 100 J = 1,N
                  TEMP = ZERO
                  DO 90 I = 1,M
                      TEMP = TEMP + A(I,J)*X(I)
   90             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  100         CONTINUE
          ELSE
              DO 120 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  DO 110 I = 1,M
                      TEMP = TEMP + A(I,J)*X(IX)
                      IX = IX + INCX
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END
*> \brief \b DGER
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA
*       INTEGER INCX,INCY,LDA,M,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGER   performs the rank 1 operation
*>
*>    A := alpha*x*y**T + A,
*>
*> where alpha is a scalar, x is an m element vector, y is an n element
*> vector and A is an m by n matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of the matrix A.
*>           M must be at least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is DOUBLE PRECISION array of dimension at least
*>           ( 1 + ( m - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the m
*>           element vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is DOUBLE PRECISION array of dimension at least
*>           ( 1 + ( n - 1 )*abs( INCY ) ).
*>           Before entry, the incremented array Y must contain the n
*>           element vector y.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           Y. INCY must not be zero.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*>           Before entry, the leading m by n part of the array A must
*>           contain the matrix of coefficients. On exit, A is
*>           overwritten by the updated matrix.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*
*  -- Reference BLAS level2 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER INCX,INCY,LDA,M,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JY,KX
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGER  ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (INCY.GT.0) THEN
          JY = 1
      ELSE
          JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
   20     CONTINUE
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
          DO 40 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  IX = KX
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              END IF
              JY = JY + INCY
   40     CONTINUE
      END IF
*
      RETURN
*
*     End of DGER  .
*
      END
*> \brief \b DNRM2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
* 
*       .. Scalar Arguments ..
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION X(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DNRM2 returns the euclidean norm of a vector via the function
*> name, so that
*>
*>    DNRM2 := sqrt( x'*x )
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  -- This version written on 25-October-1982.
*>     Modified on 14-October-1993 to inline the call to DLASSQ.
*>     Sven Hammarling, Nag Ltd.
*> \endverbatim
*>
*  =====================================================================
      DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
*
*  -- Reference BLAS level1 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
      INTEGER IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
*     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE IF (N.EQ.1) THEN
          NORM = ABS(X(1))
      ELSE
          SCALE = ZERO
          SSQ = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
*
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (X(IX).NE.ZERO) THEN
                  ABSXI = ABS(X(IX))
                  IF (SCALE.LT.ABSXI) THEN
                      SSQ = ONE + SSQ* (SCALE/ABSXI)**2
                      SCALE = ABSXI
                  ELSE
                      SSQ = SSQ + (ABSXI/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
*
      DNRM2 = NORM
      RETURN
*
*     End of DNRM2.
*
      END
*> \brief \b DROT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION C,S
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION DX(*),DY(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    DROT applies a plane rotation.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
*
*  -- Reference BLAS level1 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION C,S
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*       code for both increments equal to 1
*
         DO I = 1,N
            DTEMP = C*DX(I) + S*DY(I)
            DY(I) = C*DY(I) - S*DX(I)
            DX(I) = DTEMP
         END DO
      ELSE
*
*       code for unequal increments or equal increments not equal
*         to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = C*DX(IX) + S*DY(IY)
            DY(IY) = C*DY(IY) - S*DX(IX)
            DX(IX) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
*> \brief \b DSCAL
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSCAL(N,DA,DX,INCX)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION DA
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION DX(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    DSCAL scales a vector by a constant.
*>    uses unrolled loops for increment equal to one.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 3/11/78.
*>     modified 3/93 to return if incx .le. 0.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DSCAL(N,DA,DX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DX(I) = DA*DX(I)
            END DO
            IF (N.LT.5) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            DX(I) = DA*DX(I)
            DX(I+1) = DA*DX(I+1)
            DX(I+2) = DA*DX(I+2)
            DX(I+3) = DA*DX(I+3)
            DX(I+4) = DA*DX(I+4)
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DX(I) = DA*DX(I)
         END DO
      END IF
      RETURN
      END
*> \brief \b DSWAP
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
* 
*       .. Scalar Arguments ..
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION DX(*),DY(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    interchanges two vectors.
*>    uses unrolled loops for increments equal one.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*       code for both increments equal to 1
*
*
*       clean-up loop
*
         M = MOD(N,3)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DX(I)
               DX(I) = DY(I)
               DY(I) = DTEMP
            END DO
            IF (N.LT.3) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,3
            DTEMP = DX(I)
            DX(I) = DY(I)
            DY(I) = DTEMP
            DTEMP = DX(I+1)
            DX(I+1) = DY(I+1)
            DY(I+1) = DTEMP
            DTEMP = DX(I+2)
            DX(I+2) = DY(I+2)
            DY(I+2) = DTEMP
         END DO
      ELSE
*
*       code for unequal increments or equal increments not equal
*         to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DX(IX)
            DX(IX) = DY(IY)
            DY(IY) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
*> \brief \b DSYMV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER INCX,INCY,LDA,N
*       CHARACTER UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYMV  performs the matrix-vector  operation
*>
*>    y := alpha*A*x + beta*y,
*>
*> where alpha and beta are scalars, x and y are n element vectors and
*> A is an n by n symmetric matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the upper or lower
*>           triangular part of the array A is to be referenced as
*>           follows:
*>
*>              UPLO = 'U' or 'u'   Only the upper triangular part of A
*>                                  is to be referenced.
*>
*>              UPLO = 'L' or 'l'   Only the lower triangular part of A
*>                                  is to be referenced.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the order of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n
*>           upper triangular part of the array A must contain the upper
*>           triangular part of the symmetric matrix and the strictly
*>           lower triangular part of A is not referenced.
*>           Before entry with UPLO = 'L' or 'l', the leading n by n
*>           lower triangular part of the array A must contain the lower
*>           triangular part of the symmetric matrix and the strictly
*>           upper triangular part of A is not referenced.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, n ).
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is DOUBLE PRECISION array of dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the n
*>           element vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION.
*>           On entry, BETA specifies the scalar beta. When BETA is
*>           supplied as zero then Y need not be set on input.
*> \endverbatim
*>
*> \param[in,out] Y
*> \verbatim
*>          Y is DOUBLE PRECISION array of dimension at least
*>           ( 1 + ( n - 1 )*abs( INCY ) ).
*>           Before entry, the incremented array Y must contain the n
*>           element vector y. On exit, Y is overwritten by the updated
*>           vector y.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           Y. INCY must not be zero.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
*  -- Reference BLAS level2 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 5
      ELSE IF (INCX.EQ.0) THEN
          INFO = 7
      ELSE IF (INCY.EQ.0) THEN
          INFO = 10
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYMV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
*
*     Set up the start points in  X  and  Y.
*
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (N-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (N-1)*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
*     First form  y := beta*y.
*
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,N
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,N
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,N
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,N
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  y  when A is stored in upper triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 60 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  DO 50 I = 1,J - 1
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(I)
   50             CONTINUE
                  Y(J) = Y(J) + TEMP1*A(J,J) + ALPHA*TEMP2
   60         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 80 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  IX = KX
                  IY = KY
                  DO 70 I = 1,J - 1
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(IX)
                      IX = IX + INCX
                      IY = IY + INCY
   70             CONTINUE
                  Y(JY) = Y(JY) + TEMP1*A(J,J) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
   80         CONTINUE
          END IF
      ELSE
*
*        Form  y  when A is stored in lower triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 100 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  Y(J) = Y(J) + TEMP1*A(J,J)
                  DO 90 I = J + 1,N
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(I)
   90             CONTINUE
                  Y(J) = Y(J) + ALPHA*TEMP2
  100         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 120 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  Y(JY) = Y(JY) + TEMP1*A(J,J)
                  IX = JX
                  IY = JY
                  DO 110 I = J + 1,N
                      IX = IX + INCX
                      IY = IY + INCY
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(IX)
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSYMV .
*
      END
*> \brief \b DSYR2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA
*       INTEGER INCX,INCY,LDA,N
*       CHARACTER UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYR2  performs the symmetric rank 2 operation
*>
*>    A := alpha*x*y**T + alpha*y*x**T + A,
*>
*> where alpha is a scalar, x and y are n element vectors and A is an n
*> by n symmetric matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the upper or lower
*>           triangular part of the array A is to be referenced as
*>           follows:
*>
*>              UPLO = 'U' or 'u'   Only the upper triangular part of A
*>                                  is to be referenced.
*>
*>              UPLO = 'L' or 'l'   Only the lower triangular part of A
*>                                  is to be referenced.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the order of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is DOUBLE PRECISION array of dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the n
*>           element vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is DOUBLE PRECISION array of dimension at least
*>           ( 1 + ( n - 1 )*abs( INCY ) ).
*>           Before entry, the incremented array Y must contain the n
*>           element vector y.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           Y. INCY must not be zero.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n
*>           upper triangular part of the array A must contain the upper
*>           triangular part of the symmetric matrix and the strictly
*>           lower triangular part of A is not referenced. On exit, the
*>           upper triangular part of the array A is overwritten by the
*>           upper triangular part of the updated matrix.
*>           Before entry with UPLO = 'L' or 'l', the leading n by n
*>           lower triangular part of the array A must contain the lower
*>           triangular part of the symmetric matrix and the strictly
*>           upper triangular part of A is not referenced. On exit, the
*>           lower triangular part of the array A is overwritten by the
*>           lower triangular part of the updated matrix.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, n ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*
*  -- Reference BLAS level2 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYR2 ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Set up the start points in X and Y if the increments are not both
*     unity.
*
      IF ((INCX.NE.1) .OR. (INCY.NE.1)) THEN
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
          IF (INCY.GT.0) THEN
              KY = 1
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
          JX = KX
          JY = KY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  A  when A is stored in the upper triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 20 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      DO 10 I = 1,J
                          A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
   10                 CONTINUE
                  END IF
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = KX
                      IY = KY
                      DO 30 I = 1,J
                          A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   30                 CONTINUE
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
   40         CONTINUE
          END IF
      ELSE
*
*        Form  A  when A is stored in the lower triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 60 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      DO 50 I = J,N
                          A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
   50                 CONTINUE
                  END IF
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = JX
                      IY = JY
                      DO 70 I = J,N
                          A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
   80         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSYR2 .
*
      END
*> \brief \b DSYR2K
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER K,LDA,LDB,LDC,N
*       CHARACTER TRANS,UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYR2K  performs one of the symmetric rank 2k operations
*>
*>    C := alpha*A*B**T + alpha*B*A**T + beta*C,
*>
*> or
*>
*>    C := alpha*A**T*B + alpha*B**T*A + beta*C,
*>
*> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
*> and  A and B  are  n by k  matrices  in the  first  case  and  k by n
*> matrices in the second case.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*>           triangular  part  of the  array  C  is to be  referenced  as
*>           follows:
*>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*>                                  is to be referenced.
*>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*>                                  is to be referenced.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry,  TRANS  specifies the operation to be performed as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   C := alpha*A*B**T + alpha*B*A**T +
*>                                        beta*C.
*>
*>              TRANS = 'T' or 't'   C := alpha*A**T*B + alpha*B**T*A +
*>                                        beta*C.
*>
*>              TRANS = 'C' or 'c'   C := alpha*A**T*B + alpha*B**T*A +
*>                                        beta*C.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry,  N specifies the order of the matrix C.  N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*>           of  columns  of the  matrices  A and B,  and on  entry  with
*>           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
*>           of rows of the matrices  A and B.  K must be at least  zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*>           part of the array  A  must contain the matrix  A,  otherwise
*>           the leading  k by n  part of the array  A  must contain  the
*>           matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*>           be at least  max( 1, k ).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*>           part of the array  B  must contain the matrix  B,  otherwise
*>           the leading  k by n  part of the array  B  must contain  the
*>           matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*>           then  LDB must be at least  max( 1, n ), otherwise  LDB must
*>           be at least  max( 1, k ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION.
*>           On entry, BETA specifies the scalar beta.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*>           upper triangular part of the array C must contain the upper
*>           triangular part  of the  symmetric matrix  and the strictly
*>           lower triangular part of C is not referenced.  On exit, the
*>           upper triangular part of the array  C is overwritten by the
*>           upper triangular part of the updated matrix.
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*>           lower triangular part of the array C must contain the lower
*>           triangular part  of the  symmetric matrix  and the strictly
*>           upper triangular part of C is not referenced.  On exit, the
*>           lower triangular part of the array  C is overwritten by the
*>           lower triangular part of the updated matrix.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>           On entry, LDC specifies the first dimension of C as declared
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least
*>           max( 1, n ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,N
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*
*     Test the input parameters.
*
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND.
     +         (.NOT.LSAME(TRANS,'T')) .AND.
     +         (.NOT.LSAME(TRANS,'C'))) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (K.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDB.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 12
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYR2K',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR.
     +    (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.ZERO) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 I = J,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  C := alpha*A*B**T + alpha*B*A**T + C.
*
          IF (UPPER) THEN
              DO 130 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 100 I = 1,J
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                  END IF
                  DO 120 L = 1,K
                      IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                          TEMP1 = ALPHA*B(J,L)
                          TEMP2 = ALPHA*A(J,L)
                          DO 110 I = 1,J
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 +
     +                                 B(I,L)*TEMP2
  110                     CONTINUE
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 150 I = J,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  END IF
                  DO 170 L = 1,K
                      IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                          TEMP1 = ALPHA*B(J,L)
                          TEMP2 = ALPHA*A(J,L)
                          DO 160 I = J,N
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 +
     +                                 B(I,L)*TEMP2
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
*
*        Form  C := alpha*A**T*B + alpha*B**T*A + C.
*
          IF (UPPER) THEN
              DO 210 J = 1,N
                  DO 200 I = 1,J
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 190 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  190                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     +                             ALPHA*TEMP2
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 J = 1,N
                  DO 230 I = J,N
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 220 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  220                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     +                             ALPHA*TEMP2
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSYR2K.
*
      END
*> \brief \b DSYRK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER K,LDA,LDC,N
*       CHARACTER TRANS,UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),C(LDC,*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYRK  performs one of the symmetric rank k operations
*>
*>    C := alpha*A*A**T + beta*C,
*>
*> or
*>
*>    C := alpha*A**T*A + beta*C,
*>
*> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
*> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
*> in the second case.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*>           triangular  part  of the  array  C  is to be  referenced  as
*>           follows:
*>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*>                                  is to be referenced.
*>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*>                                  is to be referenced.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry,  TRANS  specifies the operation to be performed as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.
*>
*>              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.
*>
*>              TRANS = 'C' or 'c'   C := alpha*A**T*A + beta*C.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry,  N specifies the order of the matrix C.  N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*>           of  columns   of  the   matrix   A,   and  on   entry   with
*>           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
*>           of rows of the matrix  A.  K must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*>           part of the array  A  must contain the matrix  A,  otherwise
*>           the leading  k by n  part of the array  A  must contain  the
*>           matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*>           be at least  max( 1, k ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION.
*>           On entry, BETA specifies the scalar beta.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*>           upper triangular part of the array C must contain the upper
*>           triangular part  of the  symmetric matrix  and the strictly
*>           lower triangular part of C is not referenced.  On exit, the
*>           upper triangular part of the array  C is overwritten by the
*>           upper triangular part of the updated matrix.
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*>           lower triangular part of the array C must contain the lower
*>           triangular part  of the  symmetric matrix  and the strictly
*>           upper triangular part of C is not referenced.  On exit, the
*>           lower triangular part of the array  C is overwritten by the
*>           lower triangular part of the updated matrix.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>           On entry, LDC specifies the first dimension of C as declared
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least
*>           max( 1, n ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDC,N
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*
*     Test the input parameters.
*
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND.
     +         (.NOT.LSAME(TRANS,'T')) .AND.
     +         (.NOT.LSAME(TRANS,'C'))) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (K.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 10
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYRK ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR.
     +    (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.ZERO) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 I = J,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  C := alpha*A*A**T + beta*C.
*
          IF (UPPER) THEN
              DO 130 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 100 I = 1,J
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                  END IF
                  DO 120 L = 1,K
                      IF (A(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*A(J,L)
                          DO 110 I = 1,J
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  110                     CONTINUE
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 150 I = J,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  END IF
                  DO 170 L = 1,K
                      IF (A(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*A(J,L)
                          DO 160 I = J,N
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
*
*        Form  C := alpha*A**T*A + beta*C.
*
          IF (UPPER) THEN
              DO 210 J = 1,N
                  DO 200 I = 1,J
                      TEMP = ZERO
                      DO 190 L = 1,K
                          TEMP = TEMP + A(L,I)*A(L,J)
  190                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 J = 1,N
                  DO 230 I = J,N
                      TEMP = ZERO
                      DO 220 L = 1,K
                          TEMP = TEMP + A(L,I)*A(L,J)
  220                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSYRK .
*
      END
*> \brief \b DTRMM
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA
*       INTEGER LDA,LDB,M,N
*       CHARACTER DIAG,SIDE,TRANSA,UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),B(LDB,*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTRMM  performs one of the matrix-matrix operations
*>
*>    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
*>
*> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*>
*>    op( A ) = A   or   op( A ) = A**T.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>           On entry,  SIDE specifies whether  op( A ) multiplies B from
*>           the left or right as follows:
*>
*>              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
*>
*>              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix A is an upper or
*>           lower triangular matrix as follows:
*>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*> \endverbatim
*>
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n'   op( A ) = A.
*>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c'   op( A ) = A**T.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of B. M must be at
*>           least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of B.  N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*>           zero then  A is not referenced and  B need not be set before
*>           entry.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>           A is DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*>           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*>           upper triangular part of the array  A must contain the upper
*>           triangular matrix  and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*>           lower triangular part of the array  A must contain the lower
*>           triangular matrix  and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*>           A  are not referenced either,  but are assumed to be  unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*>           then LDA must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*>           Before entry,  the leading  m by n part of the array  B must
*>           contain the matrix  B,  and  on exit  is overwritten  by the
*>           transformed matrix.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*  -- Reference BLAS level3 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*
*     Test the input parameters.
*
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRMM ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*A*B.
*
              IF (UPPER) THEN
                  DO 50 J = 1,N
                      DO 40 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              DO 30 I = 1,K - 1
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   30                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP*A(K,K)
                              B(K,J) = TEMP
                          END IF
   40                 CONTINUE
   50             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              B(K,J) = TEMP
                              IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                              DO 60 I = K + 1,M
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   60                         CONTINUE
                          END IF
   70                 CONTINUE
   80             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*A**T*B.
*
              IF (UPPER) THEN
                  DO 110 J = 1,N
                      DO 100 I = M,1,-1
                          TEMP = B(I,J)
                          IF (NOUNIT) TEMP = TEMP*A(I,I)
                          DO 90 K = 1,I - 1
                              TEMP = TEMP + A(K,I)*B(K,J)
   90                     CONTINUE
                          B(I,J) = ALPHA*TEMP
  100                 CONTINUE
  110             CONTINUE
              ELSE
                  DO 140 J = 1,N
                      DO 130 I = 1,M
                          TEMP = B(I,J)
                          IF (NOUNIT) TEMP = TEMP*A(I,I)
                          DO 120 K = I + 1,M
                              TEMP = TEMP + A(K,I)*B(K,J)
  120                     CONTINUE
                          B(I,J) = ALPHA*TEMP
  130                 CONTINUE
  140             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*B*A.
*
              IF (UPPER) THEN
                  DO 180 J = N,1,-1
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 150 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  150                 CONTINUE
                      DO 170 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 160 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  160                         CONTINUE
                          END IF
  170                 CONTINUE
  180             CONTINUE
              ELSE
                  DO 220 J = 1,N
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 190 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  190                 CONTINUE
                      DO 210 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 200 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  200                         CONTINUE
                          END IF
  210                 CONTINUE
  220             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*A**T.
*
              IF (UPPER) THEN
                  DO 260 K = 1,N
                      DO 240 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = ALPHA*A(J,K)
                              DO 230 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(K,K)
                      IF (TEMP.NE.ONE) THEN
                          DO 250 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              ELSE
                  DO 300 K = N,1,-1
                      DO 280 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = ALPHA*A(J,K)
                              DO 270 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  270                         CONTINUE
                          END IF
  280                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(K,K)
                      IF (TEMP.NE.ONE) THEN
                          DO 290 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  290                     CONTINUE
                      END IF
  300             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of DTRMM .
*
      END
*> \brief \b DTRMV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
* 
*       .. Scalar Arguments ..
*       INTEGER INCX,LDA,N
*       CHARACTER DIAG,TRANS,UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),X(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTRMV  performs one of the matrix-vector operations
*>
*>    x := A*x,   or   x := A**T*x,
*>
*> where x is an n element vector and  A is an n by n unit, or non-unit,
*> upper or lower triangular matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix is an upper or
*>           lower triangular matrix as follows:
*>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry, TRANS specifies the operation to be performed as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   x := A*x.
*>
*>              TRANS = 'T' or 't'   x := A**T*x.
*>
*>              TRANS = 'C' or 'c'   x := A**T*x.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit
*>           triangular as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the order of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n
*>           upper triangular part of the array A must contain the upper
*>           triangular matrix and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry with UPLO = 'L' or 'l', the leading n by n
*>           lower triangular part of the array A must contain the lower
*>           triangular matrix and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*>           A are not referenced either, but are assumed to be unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, n ).
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is DOUBLE PRECISION array of dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the n
*>           element vector x. On exit, X is overwritten with the
*>           tranformed vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
*
*  -- Reference BLAS level2 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOUNIT
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRMV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (N.EQ.0) RETURN
*
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  x := A*x.
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 10 I = 1,J - 1
                              X(I) = X(I) + TEMP*A(I,J)
   10                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   20             CONTINUE
              ELSE
                  JX = KX
                  DO 40 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 30 I = 1,J - 1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX + INCX
   30                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX + INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 50 I = N,J + 1,-1
                              X(I) = X(I) + TEMP*A(I,J)
   50                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   60             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 70 I = N,J + 1,-1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX - INCX
   70                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX - INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
*
*        Form  x := A**T*x.
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 100 J = N,1,-1
                      TEMP = X(J)
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 90 I = J - 1,1,-1
                          TEMP = TEMP + A(I,J)*X(I)
   90                 CONTINUE
                      X(J) = TEMP
  100             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 120 J = N,1,-1
                      TEMP = X(JX)
                      IX = JX
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 110 I = J - 1,1,-1
                          IX = IX - INCX
                          TEMP = TEMP + A(I,J)*X(IX)
  110                 CONTINUE
                      X(JX) = TEMP
                      JX = JX - INCX
  120             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140 J = 1,N
                      TEMP = X(J)
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 130 I = J + 1,N
                          TEMP = TEMP + A(I,J)*X(I)
  130                 CONTINUE
                      X(J) = TEMP
  140             CONTINUE
              ELSE
                  JX = KX
                  DO 160 J = 1,N
                      TEMP = X(JX)
                      IX = JX
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 150 I = J + 1,N
                          IX = IX + INCX
                          TEMP = TEMP + A(I,J)*X(IX)
  150                 CONTINUE
                      X(JX) = TEMP
                      JX = JX + INCX
  160             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of DTRMV .
*
      END
*> \brief \b DTRSM
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA
*       INTEGER LDA,LDB,M,N
*       CHARACTER DIAG,SIDE,TRANSA,UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),B(LDB,*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTRSM  solves one of the matrix equations
*>
*>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*>
*> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*>
*>    op( A ) = A   or   op( A ) = A**T.
*>
*> The matrix X is overwritten on B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>           On entry, SIDE specifies whether op( A ) appears on the left
*>           or right of X as follows:
*>
*>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*>
*>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix A is an upper or
*>           lower triangular matrix as follows:
*>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*> \endverbatim
*>
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n'   op( A ) = A.
*>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c'   op( A ) = A**T.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of B. M must be at
*>           least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of B.  N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*>           zero then  A is not referenced and  B need not be set before
*>           entry.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, k ),
*>           where k is m when SIDE = 'L' or 'l'  
*>             and k is n when SIDE = 'R' or 'r'.
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*>           upper triangular part of the array  A must contain the upper
*>           triangular matrix  and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*>           lower triangular part of the array  A must contain the lower
*>           triangular matrix  and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*>           A  are not referenced either,  but are assumed to be  unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*>           then LDA must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*>           Before entry,  the leading  m by n part of the array  B must
*>           contain  the  right-hand  side  matrix  B,  and  on exit  is
*>           overwritten by the solution matrix  X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*  -- Reference BLAS level3 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*
*     Test the input parameters.
*
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRSM ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*inv( A )*B.
*
              IF (UPPER) THEN
                  DO 60 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 30 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   30                     CONTINUE
                      END IF
                      DO 50 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 70 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   70                     CONTINUE
                      END IF
                      DO 90 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*inv( A**T )*B.
*
              IF (UPPER) THEN
                  DO 130 J = 1,N
                      DO 120 I = 1,M
                          TEMP = ALPHA*B(I,J)
                          DO 110 K = 1,I - 1
                              TEMP = TEMP - A(K,I)*B(K,J)
  110                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  120                 CONTINUE
  130             CONTINUE
              ELSE
                  DO 160 J = 1,N
                      DO 150 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          DO 140 K = I + 1,M
                              TEMP = TEMP - A(K,I)*B(K,J)
  140                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*B*inv( A ).
*
              IF (UPPER) THEN
                  DO 210 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 170 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  170                     CONTINUE
                      END IF
                      DO 190 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              DO 180 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 200 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  200                     CONTINUE
                      END IF
  210             CONTINUE
              ELSE
                  DO 260 J = N,1,-1
                      IF (ALPHA.NE.ONE) THEN
                          DO 220 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  220                     CONTINUE
                      END IF
                      DO 240 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              DO 230 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 250 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*inv( A**T ).
*
              IF (UPPER) THEN
                  DO 310 K = N,1,-1
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 270 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      END IF
                      DO 290 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 280 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  280                         CONTINUE
                          END IF
  290                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 300 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  300                     CONTINUE
                      END IF
  310             CONTINUE
              ELSE
                  DO 360 K = 1,N
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 320 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  320                     CONTINUE
                      END IF
                      DO 340 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 330 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  330                         CONTINUE
                          END IF
  340                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 350 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  350                     CONTINUE
                      END IF
  360             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of DTRSM .
*
      END
*> \brief \b IDAMAX
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION IDAMAX(N,DX,INCX)
* 
*       .. Scalar Arguments ..
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION DX(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    IDAMAX finds the index of the first element having maximum absolute value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup aux_blas
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 3/11/78.
*>     modified 3/93 to return if incx .le. 0.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.6.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DABS
*     ..
      IDAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
         DMAX = DABS(DX(1))
         DO I = 2,N
            IF (DABS(DX(I)).GT.DMAX) THEN
               IDAMAX = I
               DMAX = DABS(DX(I))
            END IF
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         IX = 1
         DMAX = DABS(DX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (DABS(DX(IX)).GT.DMAX) THEN
               IDAMAX = I
               DMAX = DABS(DX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
      END
*> \brief \b IZAMAX
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION IZAMAX(N,ZX,INCX)
* 
*       .. Scalar Arguments ..
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    IZAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)|
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup aux_blas
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, 1/15/85.
*>     modified 3/93 to return if incx .le. 0.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      INTEGER FUNCTION IZAMAX(N,ZX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.6.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
*     ..
*     .. External Functions ..
      DOUBLE PRECISION DCABS1
      EXTERNAL DCABS1
*     ..
      IZAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IZAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
         DMAX = DCABS1(ZX(1))
         DO I = 2,N
            IF (DCABS1(ZX(I)).GT.DMAX) THEN
               IZAMAX = I
               DMAX = DCABS1(ZX(I))
            END IF
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         IX = 1
         DMAX = DCABS1(ZX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (DCABS1(ZX(IX)).GT.DMAX) THEN
               IZAMAX = I
               DMAX = DCABS1(ZX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
      END
*> \brief \b ZCOPY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZCOPY(N,ZX,INCX,ZY,INCY)
* 
*       .. Scalar Arguments ..
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*),ZY(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZCOPY copies a vector, x, to a vector, y.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup complex16_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 4/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZCOPY(N,ZX,INCX,ZY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO I = 1,N
          ZY(I) = ZX(I)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            ZY(IY) = ZX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
*> \brief \b ZDOTU
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       COMPLEX*16 FUNCTION ZDOTU(N,ZX,INCX,ZY,INCY)
* 
*       .. Scalar Arguments ..
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*),ZY(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZDOTU forms the dot product of two complex vectors
*>      ZDOTU = X^T * Y
*>
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup complex16_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      COMPLEX*16 FUNCTION ZDOTU(N,ZX,INCX,ZY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.6.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      COMPLEX*16 ZTEMP
      INTEGER I,IX,IY
*     ..
      ZTEMP = (0.0d0,0.0d0)
      ZDOTU = (0.0d0,0.0d0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO I = 1,N
            ZTEMP = ZTEMP + ZX(I)*ZY(I)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            ZTEMP = ZTEMP + ZX(IX)*ZY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      ZDOTU = ZTEMP
      RETURN
      END
*> \brief \b ZGEMM
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
* 
*       .. Scalar Arguments ..
*       COMPLEX*16 ALPHA,BETA
*       INTEGER K,LDA,LDB,LDC,M,N
*       CHARACTER TRANSA,TRANSB
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGEMM  performs one of the matrix-matrix operations
*>
*>    C := alpha*op( A )*op( B ) + beta*C,
*>
*> where  op( X ) is one of
*>
*>    op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
*>
*> alpha and beta are scalars, and A, B and C are matrices, with op( A )
*> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n',  op( A ) = A.
*>
*>              TRANSA = 'T' or 't',  op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c',  op( A ) = A**H.
*> \endverbatim
*>
*> \param[in] TRANSB
*> \verbatim
*>          TRANSB is CHARACTER*1
*>           On entry, TRANSB specifies the form of op( B ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSB = 'N' or 'n',  op( B ) = B.
*>
*>              TRANSB = 'T' or 't',  op( B ) = B**T.
*>
*>              TRANSB = 'C' or 'c',  op( B ) = B**H.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry,  M  specifies  the number  of rows  of the  matrix
*>           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry,  N  specifies the number  of columns of the matrix
*>           op( B ) and the number of columns of the matrix C. N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>           On entry,  K  specifies  the number of columns of the matrix
*>           op( A ) and the number of rows of the matrix op( B ). K must
*>           be at least  zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is
*>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*>           part of the array  A  must contain the matrix  A,  otherwise
*>           the leading  k by m  part of the array  A  must contain  the
*>           matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*>           least  max( 1, k ).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is COMPLEX*16 array of DIMENSION ( LDB, kb ), where kb is
*>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*>           part of the array  B  must contain the matrix  B,  otherwise
*>           the leading  n by k  part of the array  B  must contain  the
*>           matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*>           least  max( 1, n ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is COMPLEX*16
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*>           supplied as zero then C need not be set on input.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array of DIMENSION ( LDC, n ).
*>           Before entry, the leading  m by n  part of the array  C must
*>           contain the matrix  C,  except when  beta  is zero, in which
*>           case C need not be set on entry.
*>           On exit, the array  C  is overwritten by the  m by n  matrix
*>           ( alpha*op( A )*op( B ) + beta*C ).
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>           On entry, LDC specifies the first dimension of C as declared
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup complex16_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine (version 3.6.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL CONJA,CONJB,NOTA,NOTB
*     ..
*     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
*     B  respectively are to be  transposed but  not conjugated  and set
*     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
*     and the number of rows of  B  respectively.
*
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      CONJA = LSAME(TRANSA,'C')
      CONJB = LSAME(TRANSB,'C')
      IF (NOTA) THEN
          NROWA = M
          NCOLA = K
      ELSE
          NROWA = K
          NCOLA = M
      END IF
      IF (NOTB) THEN
          NROWB = K
      ELSE
          NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.CONJA) .AND.
     +    (.NOT.LSAME(TRANSA,'T'))) THEN
          INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.CONJB) .AND.
     +         (.NOT.LSAME(TRANSB,'T'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
          INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 13
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGEMM ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (NOTB) THEN
          IF (NOTA) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
              DO 90 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  END IF
                  DO 80 L = 1,K
                      TEMP = ALPHA*B(L,J)
                      DO 70 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          ELSE IF (CONJA) THEN
*
*           Form  C := alpha*A**H*B + beta*C.
*
              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*B(L,J)
  100                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  110             CONTINUE
  120         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B + beta*C
*
              DO 150 J = 1,N
                  DO 140 I = 1,M
                      TEMP = ZERO
                      DO 130 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  130                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  140             CONTINUE
  150         CONTINUE
          END IF
      ELSE IF (NOTA) THEN
          IF (CONJB) THEN
*
*           Form  C := alpha*A*B**H + beta*C.
*
              DO 200 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 160 I = 1,M
                          C(I,J) = ZERO
  160                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 170 I = 1,M
                          C(I,J) = BETA*C(I,J)
  170                 CONTINUE
                  END IF
                  DO 190 L = 1,K
                      TEMP = ALPHA*DCONJG(B(J,L))
                      DO 180 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  180                 CONTINUE
  190             CONTINUE
  200         CONTINUE
          ELSE
*
*           Form  C := alpha*A*B**T + beta*C
*
              DO 250 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 210 I = 1,M
                          C(I,J) = ZERO
  210                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 220 I = 1,M
                          C(I,J) = BETA*C(I,J)
  220                 CONTINUE
                  END IF
                  DO 240 L = 1,K
                      TEMP = ALPHA*B(J,L)
                      DO 230 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  230                 CONTINUE
  240             CONTINUE
  250         CONTINUE
          END IF
      ELSE IF (CONJA) THEN
          IF (CONJB) THEN
*
*           Form  C := alpha*A**H*B**H + beta*C.
*
              DO 280 J = 1,N
                  DO 270 I = 1,M
                      TEMP = ZERO
                      DO 260 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*DCONJG(B(J,L))
  260                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  270             CONTINUE
  280         CONTINUE
          ELSE
*
*           Form  C := alpha*A**H*B**T + beta*C
*
              DO 310 J = 1,N
                  DO 300 I = 1,M
                      TEMP = ZERO
                      DO 290 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*B(J,L)
  290                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  300             CONTINUE
  310         CONTINUE
          END IF
      ELSE
          IF (CONJB) THEN
*
*           Form  C := alpha*A**T*B**H + beta*C
*
              DO 340 J = 1,N
                  DO 330 I = 1,M
                      TEMP = ZERO
                      DO 320 L = 1,K
                          TEMP = TEMP + A(L,I)*DCONJG(B(J,L))
  320                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  330             CONTINUE
  340         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B**T + beta*C
*
              DO 370 J = 1,N
                  DO 360 I = 1,M
                      TEMP = ZERO
                      DO 350 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  350                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  360             CONTINUE
  370         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of ZGEMM .
*
      END
*> \brief \b ZGEMV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
* 
*       .. Scalar Arguments ..
*       COMPLEX*16 ALPHA,BETA
*       INTEGER INCX,INCY,LDA,M,N
*       CHARACTER TRANS
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),X(*),Y(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGEMV  performs one of the matrix-vector operations
*>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
*>
*>    y := alpha*A**H*x + beta*y,
*>
*> where alpha and beta are scalars, x and y are vectors and A is an
*> m by n matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry, TRANS specifies the operation to be performed as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*>
*>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
*>
*>              TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of the matrix A.
*>           M must be at least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array of DIMENSION ( LDA, n ).
*>           Before entry, the leading m by n part of the array A must
*>           contain the matrix of coefficients.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, m ).
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX*16 array of DIMENSION at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*>           and at least
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*>           Before entry, the incremented array X must contain the
*>           vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is COMPLEX*16
*>           On entry, BETA specifies the scalar beta. When BETA is
*>           supplied as zero then Y need not be set on input.
*> \endverbatim
*>
*> \param[in,out] Y
*> \verbatim
*>          Y is COMPLEX*16 array of DIMENSION at least
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*>           and at least
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*>           Before entry with BETA non-zero, the incremented array Y
*>           must contain the vector y. On exit, Y is overwritten by the
*>           updated vector y.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           Y. INCY must not be zero.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup complex16_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
*  -- Reference BLAS level2 routine (version 3.6.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
      LOGICAL NOCONJ
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +    .NOT.LSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGEMV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
*
      NOCONJ = LSAME(TRANS,'T')
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
          JX = KX
          IF (INCY.EQ.1) THEN
              DO 60 J = 1,N
                  TEMP = ALPHA*X(JX)
                  DO 50 I = 1,M
                      Y(I) = Y(I) + TEMP*A(I,J)
   50             CONTINUE
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  TEMP = ALPHA*X(JX)
                  IY = KY
                  DO 70 I = 1,M
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
   70             CONTINUE
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
*
*        Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
*
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 110 J = 1,N
                  TEMP = ZERO
                  IF (NOCONJ) THEN
                      DO 90 I = 1,M
                          TEMP = TEMP + A(I,J)*X(I)
   90                 CONTINUE
                  ELSE
                      DO 100 I = 1,M
                          TEMP = TEMP + DCONJG(A(I,J))*X(I)
  100                 CONTINUE
                  END IF
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  110         CONTINUE
          ELSE
              DO 140 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  IF (NOCONJ) THEN
                      DO 120 I = 1,M
                          TEMP = TEMP + A(I,J)*X(IX)
                          IX = IX + INCX
  120                 CONTINUE
                  ELSE
                      DO 130 I = 1,M
                          TEMP = TEMP + DCONJG(A(I,J))*X(IX)
                          IX = IX + INCX
  130                 CONTINUE
                  END IF
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  140         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of ZGEMV .
*
      END
*> \brief \b ZGERU
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
* 
*       .. Scalar Arguments ..
*       COMPLEX*16 ALPHA
*       INTEGER INCX,INCY,LDA,M,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),X(*),Y(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGERU  performs the rank 1 operation
*>
*>    A := alpha*x*y**T + A,
*>
*> where alpha is a scalar, x is an m element vector, y is an n element
*> vector and A is an m by n matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of the matrix A.
*>           M must be at least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX*16 array of dimension at least
*>           ( 1 + ( m - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the m
*>           element vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is COMPLEX*16 array of dimension at least
*>           ( 1 + ( n - 1 )*abs( INCY ) ).
*>           Before entry, the incremented array Y must contain the n
*>           element vector y.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           Y. INCY must not be zero.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array of DIMENSION ( LDA, n ).
*>           Before entry, the leading m by n part of the array A must
*>           contain the matrix of coefficients. On exit, A is
*>           overwritten by the updated matrix.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup complex16_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*
*  -- Reference BLAS level2 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      INTEGER INCX,INCY,LDA,M,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JY,KX
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGERU ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (INCY.GT.0) THEN
          JY = 1
      ELSE
          JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
   20     CONTINUE
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
          DO 40 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  IX = KX
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              END IF
              JY = JY + INCY
   40     CONTINUE
      END IF
*
      RETURN
*
*     End of ZGERU .
*
      END
*> \brief \b ZSCAL
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZSCAL(N,ZA,ZX,INCX)
* 
*       .. Scalar Arguments ..
*       COMPLEX*16 ZA
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZSCAL scales a vector by a constant.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup complex16_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, 3/11/78.
*>     modified 3/93 to return if incx .le. 0.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZSCAL(N,ZA,ZX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      COMPLEX*16 ZA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,NINCX
*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
         DO I = 1,N
            ZX(I) = ZA*ZX(I)
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            ZX(I) = ZA*ZX(I)
         END DO
      END IF
      RETURN
      END
*> \brief \b ZSWAP
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZSWAP(N,ZX,INCX,ZY,INCY)
* 
*       .. Scalar Arguments ..
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*),ZY(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZSWAP interchanges two vectors.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup complex16_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZSWAP(N,ZX,INCX,ZY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      COMPLEX*16 ZTEMP
      INTEGER I,IX,IY
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*       code for both increments equal to 1
         DO I = 1,N
            ZTEMP = ZX(I)
            ZX(I) = ZY(I)
            ZY(I) = ZTEMP
         END DO
      ELSE
*
*       code for unequal increments or equal increments not equal
*         to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            ZTEMP = ZX(IX)
            ZX(IX) = ZY(IY)
            ZY(IY) = ZTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
*> \brief \b ZSYR2K
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
* 
*       .. Scalar Arguments ..
*       COMPLEX*16 ALPHA,BETA
*       INTEGER K,LDA,LDB,LDC,N
*       CHARACTER TRANS,UPLO
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZSYR2K  performs one of the symmetric rank 2k operations
*>
*>    C := alpha*A*B**T + alpha*B*A**T + beta*C,
*>
*> or
*>
*>    C := alpha*A**T*B + alpha*B**T*A + beta*C,
*>
*> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
*> and  A and B  are  n by k  matrices  in the  first  case  and  k by n
*> matrices in the second case.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*>           triangular  part  of the  array  C  is to be  referenced  as
*>           follows:
*>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*>                                  is to be referenced.
*>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*>                                  is to be referenced.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry,  TRANS  specifies the operation to be performed as
*>           follows:
*>
*>              TRANS = 'N' or 'n'    C := alpha*A*B**T + alpha*B*A**T +
*>                                         beta*C.
*>
*>              TRANS = 'T' or 't'    C := alpha*A**T*B + alpha*B**T*A +
*>                                         beta*C.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry,  N specifies the order of the matrix C.  N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*>           of  columns  of the  matrices  A and B,  and on  entry  with
*>           TRANS = 'T' or 't',  K  specifies  the number of rows of the
*>           matrices  A and B.  K must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*>           part of the array  A  must contain the matrix  A,  otherwise
*>           the leading  k by n  part of the array  A  must contain  the
*>           matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*>           be at least  max( 1, k ).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is COMPLEX*16 array of DIMENSION ( LDB, kb ), where kb is
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*>           part of the array  B  must contain the matrix  B,  otherwise
*>           the leading  k by n  part of the array  B  must contain  the
*>           matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*>           then  LDB must be at least  max( 1, n ), otherwise  LDB must
*>           be at least  max( 1, k ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is COMPLEX*16
*>           On entry, BETA specifies the scalar beta.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array of DIMENSION ( LDC, n ).
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*>           upper triangular part of the array C must contain the upper
*>           triangular part  of the  symmetric matrix  and the strictly
*>           lower triangular part of C is not referenced.  On exit, the
*>           upper triangular part of the array  C is overwritten by the
*>           upper triangular part of the updated matrix.
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*>           lower triangular part of the array C must contain the lower
*>           triangular part  of the  symmetric matrix  and the strictly
*>           upper triangular part of C is not referenced.  On exit, the
*>           lower triangular part of the array  C is overwritten by the
*>           lower triangular part of the updated matrix.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>           On entry, LDC specifies the first dimension of C as declared
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least
*>           max( 1, n ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup complex16_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,N
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP1,TEMP2
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
*     ..
*     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*
*     Test the input parameters.
*
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND.
     +         (.NOT.LSAME(TRANS,'T'))) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (K.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDB.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 12
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZSYR2K',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR.
     +    (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.ZERO) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 I = J,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  C := alpha*A*B**T + alpha*B*A**T + C.
*
          IF (UPPER) THEN
              DO 130 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 100 I = 1,J
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                  END IF
                  DO 120 L = 1,K
                      IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                          TEMP1 = ALPHA*B(J,L)
                          TEMP2 = ALPHA*A(J,L)
                          DO 110 I = 1,J
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 +
     +                                 B(I,L)*TEMP2
  110                     CONTINUE
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 150 I = J,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  END IF
                  DO 170 L = 1,K
                      IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                          TEMP1 = ALPHA*B(J,L)
                          TEMP2 = ALPHA*A(J,L)
                          DO 160 I = J,N
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 +
     +                                 B(I,L)*TEMP2
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
*
*        Form  C := alpha*A**T*B + alpha*B**T*A + C.
*
          IF (UPPER) THEN
              DO 210 J = 1,N
                  DO 200 I = 1,J
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 190 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  190                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     +                             ALPHA*TEMP2
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 J = 1,N
                  DO 230 I = J,N
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 220 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  220                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     +                             ALPHA*TEMP2
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of ZSYR2K.
*
      END
*> \brief \b ZTRMM
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
* 
*       .. Scalar Arguments ..
*       COMPLEX*16 ALPHA
*       INTEGER LDA,LDB,M,N
*       CHARACTER DIAG,SIDE,TRANSA,UPLO
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),B(LDB,*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTRMM  performs one of the matrix-matrix operations
*>
*>    B := alpha*op( A )*B,   or   B := alpha*B*op( A )
*>
*> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*>
*>    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>           On entry,  SIDE specifies whether  op( A ) multiplies B from
*>           the left or right as follows:
*>
*>              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
*>
*>              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix A is an upper or
*>           lower triangular matrix as follows:
*>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*> \endverbatim
*>
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n'   op( A ) = A.
*>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c'   op( A ) = A**H.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of B. M must be at
*>           least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of B.  N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*>           zero then  A is not referenced and  B need not be set before
*>           entry.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array of DIMENSION ( LDA, k ), where k is m
*>           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*>           upper triangular part of the array  A must contain the upper
*>           triangular matrix  and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*>           lower triangular part of the array  A must contain the lower
*>           triangular matrix  and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*>           A  are not referenced either,  but are assumed to be  unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*>           then LDA must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is (input/output) COMPLEX*16 array of DIMENSION ( LDB, n ).
*>           Before entry,  the leading  m by n part of the array  B must
*>           contain the matrix  B,  and  on exit  is overwritten  by the
*>           transformed matrix.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup complex16_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*  -- Reference BLAS level3 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*
*     Test the input parameters.
*
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOCONJ = LSAME(TRANSA,'T')
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRMM ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*A*B.
*
              IF (UPPER) THEN
                  DO 50 J = 1,N
                      DO 40 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              DO 30 I = 1,K - 1
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   30                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP*A(K,K)
                              B(K,J) = TEMP
                          END IF
   40                 CONTINUE
   50             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              B(K,J) = TEMP
                              IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                              DO 60 I = K + 1,M
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   60                         CONTINUE
                          END IF
   70                 CONTINUE
   80             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*A**T*B   or   B := alpha*A**H*B.
*
              IF (UPPER) THEN
                  DO 120 J = 1,N
                      DO 110 I = M,1,-1
                          TEMP = B(I,J)
                          IF (NOCONJ) THEN
                              IF (NOUNIT) TEMP = TEMP*A(I,I)
                              DO 90 K = 1,I - 1
                                  TEMP = TEMP + A(K,I)*B(K,J)
   90                         CONTINUE
                          ELSE
                              IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                              DO 100 K = 1,I - 1
                                  TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
  100                         CONTINUE
                          END IF
                          B(I,J) = ALPHA*TEMP
  110                 CONTINUE
  120             CONTINUE
              ELSE
                  DO 160 J = 1,N
                      DO 150 I = 1,M
                          TEMP = B(I,J)
                          IF (NOCONJ) THEN
                              IF (NOUNIT) TEMP = TEMP*A(I,I)
                              DO 130 K = I + 1,M
                                  TEMP = TEMP + A(K,I)*B(K,J)
  130                         CONTINUE
                          ELSE
                              IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                              DO 140 K = I + 1,M
                                  TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
  140                         CONTINUE
                          END IF
                          B(I,J) = ALPHA*TEMP
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*B*A.
*
              IF (UPPER) THEN
                  DO 200 J = N,1,-1
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 170 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  170                 CONTINUE
                      DO 190 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 180 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
  200             CONTINUE
              ELSE
                  DO 240 J = 1,N
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 210 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  210                 CONTINUE
                      DO 230 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 220 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  220                         CONTINUE
                          END IF
  230                 CONTINUE
  240             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*A**T   or   B := alpha*B*A**H.
*
              IF (UPPER) THEN
                  DO 280 K = 1,N
                      DO 260 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = ALPHA*A(J,K)
                              ELSE
                                  TEMP = ALPHA*DCONJG(A(J,K))
                              END IF
                              DO 250 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  250                         CONTINUE
                          END IF
  260                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = TEMP*A(K,K)
                          ELSE
                              TEMP = TEMP*DCONJG(A(K,K))
                          END IF
                      END IF
                      IF (TEMP.NE.ONE) THEN
                          DO 270 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      END IF
  280             CONTINUE
              ELSE
                  DO 320 K = N,1,-1
                      DO 300 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = ALPHA*A(J,K)
                              ELSE
                                  TEMP = ALPHA*DCONJG(A(J,K))
                              END IF
                              DO 290 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  290                         CONTINUE
                          END IF
  300                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = TEMP*A(K,K)
                          ELSE
                              TEMP = TEMP*DCONJG(A(K,K))
                          END IF
                      END IF
                      IF (TEMP.NE.ONE) THEN
                          DO 310 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  310                     CONTINUE
                      END IF
  320             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of ZTRMM .
*
      END
*> \brief \b ZTRMV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
* 
*       .. Scalar Arguments ..
*       INTEGER INCX,LDA,N
*       CHARACTER DIAG,TRANS,UPLO
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),X(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTRMV  performs one of the matrix-vector operations
*>
*>    x := A*x,   or   x := A**T*x,   or   x := A**H*x,
*>
*> where x is an n element vector and  A is an n by n unit, or non-unit,
*> upper or lower triangular matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix is an upper or
*>           lower triangular matrix as follows:
*>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry, TRANS specifies the operation to be performed as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   x := A*x.
*>
*>              TRANS = 'T' or 't'   x := A**T*x.
*>
*>              TRANS = 'C' or 'c'   x := A**H*x.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit
*>           triangular as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the order of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array of DIMENSION ( LDA, n ).
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n
*>           upper triangular part of the array A must contain the upper
*>           triangular matrix and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry with UPLO = 'L' or 'l', the leading n by n
*>           lower triangular part of the array A must contain the lower
*>           triangular matrix and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*>           A are not referenced either, but are assumed to be unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, n ).
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is (input/output) COMPLEX*16 array of dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the n
*>           element vector x. On exit, X is overwritten with the
*>           tranformed vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup complex16_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
*
*  -- Reference BLAS level2 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOCONJ,NOUNIT
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRMV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (N.EQ.0) RETURN
*
      NOCONJ = LSAME(TRANS,'T')
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  x := A*x.
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 10 I = 1,J - 1
                              X(I) = X(I) + TEMP*A(I,J)
   10                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   20             CONTINUE
              ELSE
                  JX = KX
                  DO 40 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 30 I = 1,J - 1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX + INCX
   30                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX + INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 50 I = N,J + 1,-1
                              X(I) = X(I) + TEMP*A(I,J)
   50                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   60             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 70 I = N,J + 1,-1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX - INCX
   70                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX - INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
*
*        Form  x := A**T*x  or  x := A**H*x.
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 110 J = N,1,-1
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 90 I = J - 1,1,-1
                              TEMP = TEMP + A(I,J)*X(I)
   90                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 100 I = J - 1,1,-1
                              TEMP = TEMP + DCONJG(A(I,J))*X(I)
  100                     CONTINUE
                      END IF
                      X(J) = TEMP
  110             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 140 J = N,1,-1
                      TEMP = X(JX)
                      IX = JX
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 120 I = J - 1,1,-1
                              IX = IX - INCX
                              TEMP = TEMP + A(I,J)*X(IX)
  120                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 130 I = J - 1,1,-1
                              IX = IX - INCX
                              TEMP = TEMP + DCONJG(A(I,J))*X(IX)
  130                     CONTINUE
                      END IF
                      X(JX) = TEMP
                      JX = JX - INCX
  140             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 170 J = 1,N
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 150 I = J + 1,N
                              TEMP = TEMP + A(I,J)*X(I)
  150                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 160 I = J + 1,N
                              TEMP = TEMP + DCONJG(A(I,J))*X(I)
  160                     CONTINUE
                      END IF
                      X(J) = TEMP
  170             CONTINUE
              ELSE
                  JX = KX
                  DO 200 J = 1,N
                      TEMP = X(JX)
                      IX = JX
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 180 I = J + 1,N
                              IX = IX + INCX
                              TEMP = TEMP + A(I,J)*X(IX)
  180                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 190 I = J + 1,N
                              IX = IX + INCX
                              TEMP = TEMP + DCONJG(A(I,J))*X(IX)
  190                     CONTINUE
                      END IF
                      X(JX) = TEMP
                      JX = JX + INCX
  200             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of ZTRMV .
*
      END
*> \brief \b ZTRSM
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
* 
*       .. Scalar Arguments ..
*       COMPLEX*16 ALPHA
*       INTEGER LDA,LDB,M,N
*       CHARACTER DIAG,SIDE,TRANSA,UPLO
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),B(LDB,*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTRSM  solves one of the matrix equations
*>
*>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*>
*> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*>
*>    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
*>
*> The matrix X is overwritten on B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>           On entry, SIDE specifies whether op( A ) appears on the left
*>           or right of X as follows:
*>
*>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*>
*>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix A is an upper or
*>           lower triangular matrix as follows:
*>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*> \endverbatim
*>
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n'   op( A ) = A.
*>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c'   op( A ) = A**H.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of B. M must be at
*>           least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of B.  N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*>           zero then  A is not referenced and  B need not be set before
*>           entry.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array of DIMENSION ( LDA, k ),
*>           where k is m when SIDE = 'L' or 'l'  
*>             and k is n when SIDE = 'R' or 'r'.
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*>           upper triangular part of the array  A must contain the upper
*>           triangular matrix  and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*>           lower triangular part of the array  A must contain the lower
*>           triangular matrix  and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*>           A  are not referenced either,  but are assumed to be  unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*>           then LDA must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX*16 array of DIMENSION ( LDB, n ).
*>           Before entry,  the leading  m by n part of the array  B must
*>           contain  the  right-hand  side  matrix  B,  and  on exit  is
*>           overwritten by the solution matrix  X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup complex16_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*  -- Reference BLAS level3 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(LDA,*),B(LDB,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*
*     Test the input parameters.
*
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOCONJ = LSAME(TRANSA,'T')
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRSM ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*inv( A )*B.
*
              IF (UPPER) THEN
                  DO 60 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 30 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   30                     CONTINUE
                      END IF
                      DO 50 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 70 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   70                     CONTINUE
                      END IF
                      DO 90 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*inv( A**T )*B
*           or    B := alpha*inv( A**H )*B.
*
              IF (UPPER) THEN
                  DO 140 J = 1,N
                      DO 130 I = 1,M
                          TEMP = ALPHA*B(I,J)
                          IF (NOCONJ) THEN
                              DO 110 K = 1,I - 1
                                  TEMP = TEMP - A(K,I)*B(K,J)
  110                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/A(I,I)
                          ELSE
                              DO 120 K = 1,I - 1
                                  TEMP = TEMP - DCONJG(A(K,I))*B(K,J)
  120                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
                          END IF
                          B(I,J) = TEMP
  130                 CONTINUE
  140             CONTINUE
              ELSE
                  DO 180 J = 1,N
                      DO 170 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          IF (NOCONJ) THEN
                              DO 150 K = I + 1,M
                                  TEMP = TEMP - A(K,I)*B(K,J)
  150                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/A(I,I)
                          ELSE
                              DO 160 K = I + 1,M
                                  TEMP = TEMP - DCONJG(A(K,I))*B(K,J)
  160                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
                          END IF
                          B(I,J) = TEMP
  170                 CONTINUE
  180             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*B*inv( A ).
*
              IF (UPPER) THEN
                  DO 230 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 190 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  190                     CONTINUE
                      END IF
                      DO 210 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              DO 200 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  200                         CONTINUE
                          END IF
  210                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 220 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  220                     CONTINUE
                      END IF
  230             CONTINUE
              ELSE
                  DO 280 J = N,1,-1
                      IF (ALPHA.NE.ONE) THEN
                          DO 240 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  240                     CONTINUE
                      END IF
                      DO 260 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              DO 250 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  250                         CONTINUE
                          END IF
  260                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 270 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  270                     CONTINUE
                      END IF
  280             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*inv( A**T )
*           or    B := alpha*B*inv( A**H ).
*
              IF (UPPER) THEN
                  DO 330 K = N,1,-1
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = ONE/A(K,K)
                          ELSE
                              TEMP = ONE/DCONJG(A(K,K))
                          END IF
                          DO 290 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  290                     CONTINUE
                      END IF
                      DO 310 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = A(J,K)
                              ELSE
                                  TEMP = DCONJG(A(J,K))
                              END IF
                              DO 300 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  300                         CONTINUE
                          END IF
  310                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 320 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  320                     CONTINUE
                      END IF
  330             CONTINUE
              ELSE
                  DO 380 K = 1,N
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = ONE/A(K,K)
                          ELSE
                              TEMP = ONE/DCONJG(A(K,K))
                          END IF
                          DO 340 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  340                     CONTINUE
                      END IF
                      DO 360 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = A(J,K)
                              ELSE
                                  TEMP = DCONJG(A(J,K))
                              END IF
                              DO 350 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  350                         CONTINUE
                          END IF
  360                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 370 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  370                     CONTINUE
                      END IF
  380             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of ZTRSM .
*
      END
*> \brief \b DGEQR2 computes the QR factorization of a general rectangular matrix using an unblocked algorithm.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DGEQR2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqr2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqr2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqr2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGEQR2 computes a QR factorization of a real m by n matrix A:
*> A = Q * R.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the m by n matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the min(m,n) by n upper trapezoidal matrix R (R is
*>          upper triangular if m >= n); the elements below the diagonal,
*>          with the array TAU, represent the orthogonal matrix Q as a
*>          product of elementary reflectors (see Further Details).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))
*>          The scalar factors of the elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleGEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix Q is represented as a product of elementary reflectors
*>
*>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**T
*>
*>  where tau is a real scalar, and v is a real vector with
*>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*>  and tau in TAU(i).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K
      DOUBLE PRECISION   AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQR2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
*
      DO 10 I = 1, K
*
*        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
*
         CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,
     $                TAU( I ) )
         IF( I.LT.N ) THEN
*
*           Apply H(i) to A(i:m,i+1:n) from the left
*
            AII = A( I, I )
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
            A( I, I ) = AII
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGEQR2
*
      END
*> \brief \b DGEQRF
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DGEQRF + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqrf.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqrf.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqrf.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGEQRF computes a QR factorization of a real M-by-N matrix A:
*> A = Q * R.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*>          upper triangular if m >= n); the elements below the diagonal,
*>          with the array TAU, represent the orthogonal matrix Q as a
*>          product of min(m,n) elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))
*>          The scalar factors of the elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,N).
*>          For optimum performance LWORK >= N*NB, where NB is
*>          the optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup doubleGEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix Q is represented as a product of elementary reflectors
*>
*>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**T
*>
*>  where tau is a real scalar, and v is a real vector with
*>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*>  and tau in TAU(i).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEQR2, DLARFB, DLARFT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DGEQRF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DGEQRF', ' ', M, N, -1,
     $                 -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code initially
*
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
*
*           Compute the QR factorization of the current block
*           A(i:m,i:i+ib-1)
*
            CALL DGEQR2( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H**T to A(i:m,i+ib:n) from the left
*
               CALL DLARFB( 'Left', 'Transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
*
*     Use unblocked code to factor the last or only block.
*
      IF( I.LE.K )
     $   CALL DGEQR2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,
     $                IINFO )
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DGEQRF
*
      END
*> \brief \b DGETF2 computes the LU factorization of a general m-by-n matrix using partial pivoting with row interchanges (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DGETF2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetf2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetf2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetf2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGETF2 computes an LU factorization of a general m-by-n matrix A
*> using partial pivoting with row interchanges.
*>
*> The factorization has the form
*>    A = P * L * U
*> where P is a permutation matrix, L is lower triangular with unit
*> diagonal elements (lower trapezoidal if m > n), and U is upper
*> triangular (upper trapezoidal if m < n).
*>
*> This is the right-looking Level 2 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the m by n matrix to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U; the unit diagonal elements of L are not stored.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (min(M,N))
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -k, the k-th argument had an illegal value
*>          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*>               has been completed, but the factor U is exactly
*>               singular, and division by zero will occur if it is used
*>               to solve a system of equations.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleGEcomputational
*
*  =====================================================================
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN 
      INTEGER            I, J, JP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH      
      INTEGER            IDAMAX
      EXTERNAL           DLAMCH, IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Compute machine safe minimum 
* 
      SFMIN = DLAMCH('S')  
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M ) THEN 
               IF( ABS(A( J, J )) .GE. SFMIN ) THEN 
                  CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 ) 
               ELSE 
                 DO 20 I = 1, M-J 
                    A( J+I, J ) = A( J+I, J ) / A( J, J ) 
   20            CONTINUE 
               END IF 
            END IF 
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGETF2
*
      END
*> \brief \b DGETRF2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       RECURSIVE SUBROUTINE DGETRF2( M, N, A, LDA, IPIV, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGETRF2 computes an LU factorization of a general M-by-N matrix A
*> using partial pivoting with row interchanges.
*>
*> The factorization has the form
*>    A = P * L * U
*> where P is a permutation matrix, L is lower triangular with unit
*> diagonal elements (lower trapezoidal if m > n), and U is upper
*> triangular (upper trapezoidal if m < n).
*>
*> This is the recursive version of the algorithm. It divides
*> the matrix into four submatrices:
*>            
*>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
*>    A = [ -----|----- ]  with n1 = min(m,n)
*>        [  A21 | A22  ]       n2 = n-n1
*>            
*>                                       [ A11 ]
*> The subroutine calls itself to factor [ --- ],
*>                                       [ A12 ]
*>                 [ A12 ]
*> do the swaps on [ --- ], solve A12, update A22,
*>                 [ A22 ]
*>
*> then calls itself to factor A22 and do the swaps on A21.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the M-by-N matrix to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U; the unit diagonal elements of L are not stored.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (min(M,N))
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*>                has been completed, but the factor U is exactly
*>                singular, and division by zero will occur if it is used
*>                to solve a system of equations.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup doubleGEcomputational
*
*  =====================================================================
      RECURSIVE SUBROUTINE dgetrf2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( lda, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN, TEMP
      INTEGER            I, IINFO, N1, N2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            IDAMAX
      EXTERNAL           dlamch, idamax
*     ..
*     .. External Subroutines ..
      EXTERNAL           dgemm, dscal, dlaswp, dtrsm, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      info = 0
      IF( m.LT.0 ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, m ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGETRF2', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 )
     $   RETURN

      IF ( m.EQ.1 ) THEN
*
*        Use unblocked code for one row case
*        Just need to handle IPIV and INFO
*
         ipiv( 1 ) = 1
         IF ( a(1,1).EQ.zero )
     $      info = 1
*
      ELSE IF( n.EQ.1 ) THEN
*
*        Use unblocked code for one column case
*
*
*        Compute machine safe minimum
*
         sfmin = dlamch('S')
*
*        Find pivot and test for singularity
*
         i = idamax( m, a( 1, 1 ), 1 )
         ipiv( 1 ) = i
         IF( a( i, 1 ).NE.zero ) THEN
*
*           Apply the interchange
*
            IF( i.NE.1 ) THEN
               temp = a( 1, 1 )
               a( 1, 1 ) = a( i, 1 )
               a( i, 1 ) = temp
            END IF
*
*           Compute elements 2:M of the column
*
            IF( abs(a( 1, 1 )) .GE. sfmin ) THEN
               CALL dscal( m-1, one / a( 1, 1 ), a( 2, 1 ), 1 )
            ELSE
               DO 10 i = 1, m-1
                  a( 1+i, 1 ) = a( 1+i, 1 ) / a( 1, 1 )
   10          CONTINUE
            END IF
*
         ELSE
            info = 1
         END IF
*
      ELSE
*
*        Use recursive code
*
         n1 = min( m, n ) / 2
         n2 = n-n1
*
*               [ A11 ]
*        Factor [ --- ]
*               [ A21 ]
*
         CALL dgetrf2( m, n1, a, lda, ipiv, iinfo )

         IF ( info.EQ.0 .AND. iinfo.GT.0 )
     $      info = iinfo
*
*                              [ A12 ]
*        Apply interchanges to [ --- ]
*                              [ A22 ]
*
         CALL dlaswp( n2, a( 1, n1+1 ), lda, 1, n1, ipiv, 1 )
*
*        Solve A12
*
         CALL dtrsm( 'L', 'L', 'N', 'U', n1, n2, one, a, lda, 
     $               a( 1, n1+1 ), lda )
*
*        Update A22
*
         CALL dgemm( 'N', 'N', m-n1, n2, n1, -one, a( n1+1, 1 ), lda, 
     $               a( 1, n1+1 ), lda, one, a( n1+1, n1+1 ), lda )
*
*        Factor A22
*
         CALL dgetrf2( m-n1, n2, a( n1+1, n1+1 ), lda, ipiv( n1+1 ),
     $                 iinfo )
*
*        Adjust INFO and the pivot indices
*
         IF ( info.EQ.0 .AND. iinfo.GT.0 )
     $      info = iinfo + n1
         DO 20 i = n1+1, min( m, n )
            ipiv( i ) = ipiv( i ) + n1
   20    CONTINUE
*
*        Apply interchanges to A21
*
         CALL dlaswp( n1, a( 1, 1 ), lda, n1+1, min( m, n), ipiv, 1 )
*
      END IF
      RETURN
*
*     End of DGETRF2
*
      END
*> \brief \b DGETRF
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DGETRF + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetrf.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetrf.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetrf.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGETRF computes an LU factorization of a general M-by-N matrix A
*> using partial pivoting with row interchanges.
*>
*> The factorization has the form
*>    A = P * L * U
*> where P is a permutation matrix, L is lower triangular with unit
*> diagonal elements (lower trapezoidal if m > n), and U is upper
*> triangular (upper trapezoidal if m < n).
*>
*> This is the right-looking Level 3 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the M-by-N matrix to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U; the unit diagonal elements of L are not stored.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (min(M,N))
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*>                has been completed, but the factor U is exactly
*>                singular, and division by zero will occur if it is used
*>                to solve a system of equations.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup doubleGEcomputational
*
*  =====================================================================
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETRF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL DGETRF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL DGETRF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DGETRF
*
      END
*> \brief \b DGETRI
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DGETRI + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetri.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetri.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetri.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LWORK, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGETRI computes the inverse of a matrix using the LU factorization
*> computed by DGETRF.
*>
*> This method inverts U and then computes inv(A) by solving the system
*> inv(A)*L = inv(U) for inv(A).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the factors L and U from the factorization
*>          A = P*L*U as computed by DGETRF.
*>          On exit, if INFO = 0, the inverse of the original matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,N).
*>          For optimal performance LWORK >= N*NB, where NB is
*>          the optimal blocksize returned by ILAENV.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
*>                singular and its inverse could not be computed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup doubleGEcomputational
*
*  =====================================================================
      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB,
     $                   NBMIN, NN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGEMV, DSWAP, DTRSM, DTRTRI, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NB = ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRI', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,
*     and the inverse is not computed.
*
      CALL DTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
      IF( INFO.GT.0 )
     $   RETURN
*
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = MAX( LDWORK*NB, 1 )
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DGETRI', ' ', N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = N
      END IF
*
*     Solve the equation inv(A)*L = inv(U) for inv(A).
*
      IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         DO 20 J = N, 1, -1
*
*           Copy current column of L to WORK and replace with zeros.
*
            DO 10 I = J + 1, N
               WORK( I ) = A( I, J )
               A( I, J ) = ZERO
   10       CONTINUE
*
*           Compute current column of inv(A).
*
            IF( J.LT.N )
     $         CALL DGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
   20    CONTINUE
      ELSE
*
*        Use blocked code.
*
         NN = ( ( N-1 ) / NB )*NB + 1
         DO 50 J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )
*
*           Copy current block column of L to WORK and replace with
*           zeros.
*
            DO 40 JJ = J, J + JB - 1
               DO 30 I = JJ + 1, N
                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
                  A( I, JJ ) = ZERO
   30          CONTINUE
   40       CONTINUE
*
*           Compute current block column of inv(A).
*
            IF( J+JB.LE.N )
     $         CALL DGEMM( 'No transpose', 'No transpose', N, JB,
     $                     N-J-JB+1, -ONE, A( 1, J+JB ), LDA,
     $                     WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
            CALL DTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
     $                  ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
   50    CONTINUE
      END IF
*
*     Apply column interchanges.
*
      DO 60 J = N - 1, 1, -1
         JP = IPIV( J )
         IF( JP.NE.J )
     $      CALL DSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
   60 CONTINUE
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DGETRI
*
      END
*> \brief \b DGGBAK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DGGBAK + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggbak.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggbak.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggbak.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V,
*                          LDV, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          JOB, SIDE
*       INTEGER            IHI, ILO, INFO, LDV, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   LSCALE( * ), RSCALE( * ), V( LDV, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGGBAK forms the right or left eigenvectors of a real generalized
*> eigenvalue problem A*x = lambda*B*x, by backward transformation on
*> the computed eigenvectors of the balanced pair of matrices output by
*> DGGBAL.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOB
*> \verbatim
*>          JOB is CHARACTER*1
*>          Specifies the type of backward transformation required:
*>          = 'N':  do nothing, return immediately;
*>          = 'P':  do backward transformation for permutation only;
*>          = 'S':  do backward transformation for scaling only;
*>          = 'B':  do backward transformations for both permutation and
*>                  scaling.
*>          JOB must be the same as the argument JOB supplied to DGGBAL.
*> \endverbatim
*>
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'R':  V contains right eigenvectors;
*>          = 'L':  V contains left eigenvectors.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of rows of the matrix V.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>          The integers ILO and IHI determined by DGGBAL.
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*> \endverbatim
*>
*> \param[in] LSCALE
*> \verbatim
*>          LSCALE is DOUBLE PRECISION array, dimension (N)
*>          Details of the permutations and/or scaling factors applied
*>          to the left side of A and B, as returned by DGGBAL.
*> \endverbatim
*>
*> \param[in] RSCALE
*> \verbatim
*>          RSCALE is DOUBLE PRECISION array, dimension (N)
*>          Details of the permutations and/or scaling factors applied
*>          to the right side of A and B, as returned by DGGBAL.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of columns of the matrix V.  M >= 0.
*> \endverbatim
*>
*> \param[in,out] V
*> \verbatim
*>          V is DOUBLE PRECISION array, dimension (LDV,M)
*>          On entry, the matrix of right or left eigenvectors to be
*>          transformed, as returned by DTGEVC.
*>          On exit, V is overwritten by the transformed eigenvectors.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the matrix V. LDV >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup doubleGBcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  See R.C. Ward, Balancing the generalized eigenvalue problem,
*>                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V,
     $                   LDV, INFO )
*
*  -- LAPACK computational routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   LSCALE( * ), RSCALE( * ), V( LDV, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, K
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, INT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      RIGHTV = LSAME( SIDE, 'R' )
      LEFTV = LSAME( SIDE, 'L' )
*
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND.
     $    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -4
      ELSE IF( N.EQ.0 .AND. IHI.EQ.0 .AND. ILO.NE.1 ) THEN
         INFO = -4
      ELSE IF( N.GT.0 .AND. ( IHI.LT.ILO .OR. IHI.GT.MAX( 1, N ) ) )
     $   THEN
         INFO = -5
      ELSE IF( N.EQ.0 .AND. ILO.EQ.1 .AND. IHI.NE.0 ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -8
      ELSE IF( LDV.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGBAK', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( M.EQ.0 )
     $   RETURN
      IF( LSAME( JOB, 'N' ) )
     $   RETURN
*
      IF( ILO.EQ.IHI )
     $   GO TO 30
*
*     Backward balance
*
      IF( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) THEN
*
*        Backward transformation on right eigenvectors
*
         IF( RIGHTV ) THEN
            DO 10 I = ILO, IHI
               CALL DSCAL( M, RSCALE( I ), V( I, 1 ), LDV )
   10       CONTINUE
         END IF
*
*        Backward transformation on left eigenvectors
*
         IF( LEFTV ) THEN
            DO 20 I = ILO, IHI
               CALL DSCAL( M, LSCALE( I ), V( I, 1 ), LDV )
   20       CONTINUE
         END IF
      END IF
*
*     Backward permutation
*
   30 CONTINUE
      IF( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) THEN
*
*        Backward permutation on right eigenvectors
*
         IF( RIGHTV ) THEN
            IF( ILO.EQ.1 )
     $         GO TO 50
*
            DO 40 I = ILO - 1, 1, -1
               K = INT(RSCALE( I ))
               IF( K.EQ.I )
     $            GO TO 40
               CALL DSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   40       CONTINUE
*
   50       CONTINUE
            IF( IHI.EQ.N )
     $         GO TO 70
            DO 60 I = IHI + 1, N
               K = INT(RSCALE( I ))
               IF( K.EQ.I )
     $            GO TO 60
               CALL DSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   60       CONTINUE
         END IF
*
*        Backward permutation on left eigenvectors
*
   70    CONTINUE
         IF( LEFTV ) THEN
            IF( ILO.EQ.1 )
     $         GO TO 90
            DO 80 I = ILO - 1, 1, -1
               K = INT(LSCALE( I ))
               IF( K.EQ.I )
     $            GO TO 80
               CALL DSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   80       CONTINUE
*
   90       CONTINUE
            IF( IHI.EQ.N )
     $         GO TO 110
            DO 100 I = IHI + 1, N
               K = INT(LSCALE( I ))
               IF( K.EQ.I )
     $            GO TO 100
               CALL DSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
  100       CONTINUE
         END IF
      END IF
*
  110 CONTINUE
*
      RETURN
*
*     End of DGGBAK
*
      END
*> \brief \b DGGBAL
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DGGBAL + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggbal.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggbal.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggbal.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE,
*                          RSCALE, WORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          JOB
*       INTEGER            IHI, ILO, INFO, LDA, LDB, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), LSCALE( * ),
*      $                   RSCALE( * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGGBAL balances a pair of general real matrices (A,B).  This
*> involves, first, permuting A and B by similarity transformations to
*> isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N
*> elements on the diagonal; and second, applying a diagonal similarity
*> transformation to rows and columns ILO to IHI to make the rows
*> and columns as close in norm as possible. Both steps are optional.
*>
*> Balancing may reduce the 1-norm of the matrices, and improve the
*> accuracy of the computed eigenvalues and/or eigenvectors in the
*> generalized eigenvalue problem A*x = lambda*B*x.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOB
*> \verbatim
*>          JOB is CHARACTER*1
*>          Specifies the operations to be performed on A and B:
*>          = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0
*>                  and RSCALE(I) = 1.0 for i = 1,...,N.
*>          = 'P':  permute only;
*>          = 'S':  scale only;
*>          = 'B':  both permute and scale.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the input matrix A.
*>          On exit,  A is overwritten by the balanced matrix.
*>          If JOB = 'N', A is not referenced.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,N)
*>          On entry, the input matrix B.
*>          On exit,  B is overwritten by the balanced matrix.
*>          If JOB = 'N', B is not referenced.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[out] IHI
*> \verbatim
*>          IHI is INTEGER
*>          ILO and IHI are set to integers such that on exit
*>          A(i,j) = 0 and B(i,j) = 0 if i > j and
*>          j = 1,...,ILO-1 or i = IHI+1,...,N.
*>          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
*> \endverbatim
*>
*> \param[out] LSCALE
*> \verbatim
*>          LSCALE is DOUBLE PRECISION array, dimension (N)
*>          Details of the permutations and scaling factors applied
*>          to the left side of A and B.  If P(j) is the index of the
*>          row interchanged with row j, and D(j)
*>          is the scaling factor applied to row j, then
*>            LSCALE(j) = P(j)    for J = 1,...,ILO-1
*>                      = D(j)    for J = ILO,...,IHI
*>                      = P(j)    for J = IHI+1,...,N.
*>          The order in which the interchanges are made is N to IHI+1,
*>          then 1 to ILO-1.
*> \endverbatim
*>
*> \param[out] RSCALE
*> \verbatim
*>          RSCALE is DOUBLE PRECISION array, dimension (N)
*>          Details of the permutations and scaling factors applied
*>          to the right side of A and B.  If P(j) is the index of the
*>          column interchanged with column j, and D(j)
*>          is the scaling factor applied to column j, then
*>            LSCALE(j) = P(j)    for J = 1,...,ILO-1
*>                      = D(j)    for J = ILO,...,IHI
*>                      = P(j)    for J = IHI+1,...,N.
*>          The order in which the interchanges are made is N to IHI+1,
*>          then 1 to ILO-1.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (lwork)
*>          lwork must be at least max(1,6*N) when JOB = 'S' or 'B', and
*>          at least 1 when JOB = 'N' or 'P'.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup doubleGBcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  See R.C. WARD, Balancing the generalized eigenvalue problem,
*>                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE,
     $                   RSCALE, WORK, INFO )
*
*  -- LAPACK computational routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, LDB, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), LSCALE( * ),
     $                   RSCALE( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   THREE, SCLFAC
      PARAMETER          ( THREE = 3.0D+0, SCLFAC = 1.0D+1 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICAB, IFLOW, IP1, IR, IRAB, IT, J, JC, JP1,
     $                   K, KOUNT, L, LCAB, LM1, LRAB, LSFMAX, LSFMIN,
     $                   M, NR, NRP2
      DOUBLE PRECISION   ALPHA, BASL, BETA, CAB, CMAX, COEF, COEF2,
     $                   COEF5, COR, EW, EWC, GAMMA, PGAMMA, RAB, SFMAX,
     $                   SFMIN, SUM, T, TA, TB, TC
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DDOT, DLAMCH
      EXTERNAL           LSAME, IDAMAX, DDOT, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG10, MAX, MIN, SIGN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND.
     $    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGBAL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         ILO = 1
         IHI = N
         RETURN
      END IF
*
      IF( N.EQ.1 ) THEN
         ILO = 1
         IHI = N
         LSCALE( 1 ) = ONE
         RSCALE( 1 ) = ONE
         RETURN
      END IF
*
      IF( LSAME( JOB, 'N' ) ) THEN
         ILO = 1
         IHI = N
         DO 10 I = 1, N
            LSCALE( I ) = ONE
            RSCALE( I ) = ONE
   10    CONTINUE
         RETURN
      END IF
*
      K = 1
      L = N
      IF( LSAME( JOB, 'S' ) )
     $   GO TO 190
*
      GO TO 30
*
*     Permute the matrices A and B to isolate the eigenvalues.
*
*     Find row with one nonzero in columns 1 through L
*
   20 CONTINUE
      L = LM1
      IF( L.NE.1 )
     $   GO TO 30
*
      RSCALE( 1 ) = ONE
      LSCALE( 1 ) = ONE
      GO TO 190
*
   30 CONTINUE
      LM1 = L - 1
      DO 80 I = L, 1, -1
         DO 40 J = 1, LM1
            JP1 = J + 1
            IF( A( I, J ).NE.ZERO .OR. B( I, J ).NE.ZERO )
     $         GO TO 50
   40    CONTINUE
         J = L
         GO TO 70
*
   50    CONTINUE
         DO 60 J = JP1, L
            IF( A( I, J ).NE.ZERO .OR. B( I, J ).NE.ZERO )
     $         GO TO 80
   60    CONTINUE
         J = JP1 - 1
*
   70    CONTINUE
         M = L
         IFLOW = 1
         GO TO 160
   80 CONTINUE
      GO TO 100
*
*     Find column with one nonzero in rows K through N
*
   90 CONTINUE
      K = K + 1
*
  100 CONTINUE
      DO 150 J = K, L
         DO 110 I = K, LM1
            IP1 = I + 1
            IF( A( I, J ).NE.ZERO .OR. B( I, J ).NE.ZERO )
     $         GO TO 120
  110    CONTINUE
         I = L
         GO TO 140
  120    CONTINUE
         DO 130 I = IP1, L
            IF( A( I, J ).NE.ZERO .OR. B( I, J ).NE.ZERO )
     $         GO TO 150
  130    CONTINUE
         I = IP1 - 1
  140    CONTINUE
         M = K
         IFLOW = 2
         GO TO 160
  150 CONTINUE
      GO TO 190
*
*     Permute rows M and I
*
  160 CONTINUE
      LSCALE( M ) = I
      IF( I.EQ.M )
     $   GO TO 170
      CALL DSWAP( N-K+1, A( I, K ), LDA, A( M, K ), LDA )
      CALL DSWAP( N-K+1, B( I, K ), LDB, B( M, K ), LDB )
*
*     Permute columns M and J
*
  170 CONTINUE
      RSCALE( M ) = J
      IF( J.EQ.M )
     $   GO TO 180
      CALL DSWAP( L, A( 1, J ), 1, A( 1, M ), 1 )
      CALL DSWAP( L, B( 1, J ), 1, B( 1, M ), 1 )
*
  180 CONTINUE
      GO TO ( 20, 90 )IFLOW
*
  190 CONTINUE
      ILO = K
      IHI = L
*
      IF( LSAME( JOB, 'P' ) ) THEN
         DO 195 I = ILO, IHI
            LSCALE( I ) = ONE
            RSCALE( I ) = ONE
  195    CONTINUE
         RETURN
      END IF
*
      IF( ILO.EQ.IHI )
     $   RETURN
*
*     Balance the submatrix in rows ILO to IHI.
*
      NR = IHI - ILO + 1
      DO 200 I = ILO, IHI
         RSCALE( I ) = ZERO
         LSCALE( I ) = ZERO
*
         WORK( I ) = ZERO
         WORK( I+N ) = ZERO
         WORK( I+2*N ) = ZERO
         WORK( I+3*N ) = ZERO
         WORK( I+4*N ) = ZERO
         WORK( I+5*N ) = ZERO
  200 CONTINUE
*
*     Compute right side vector in resulting linear equations
*
      BASL = LOG10( SCLFAC )
      DO 240 I = ILO, IHI
         DO 230 J = ILO, IHI
            TB = B( I, J )
            TA = A( I, J )
            IF( TA.EQ.ZERO )
     $         GO TO 210
            TA = LOG10( ABS( TA ) ) / BASL
  210       CONTINUE
            IF( TB.EQ.ZERO )
     $         GO TO 220
            TB = LOG10( ABS( TB ) ) / BASL
  220       CONTINUE
            WORK( I+4*N ) = WORK( I+4*N ) - TA - TB
            WORK( J+5*N ) = WORK( J+5*N ) - TA - TB
  230    CONTINUE
  240 CONTINUE
*
      COEF = ONE / DBLE( 2*NR )
      COEF2 = COEF*COEF
      COEF5 = HALF*COEF2
      NRP2 = NR + 2
      BETA = ZERO
      IT = 1
*
*     Start generalized conjugate gradient iteration
*
  250 CONTINUE
*
      GAMMA = DDOT( NR, WORK( ILO+4*N ), 1, WORK( ILO+4*N ), 1 ) +
     $        DDOT( NR, WORK( ILO+5*N ), 1, WORK( ILO+5*N ), 1 )
*
      EW = ZERO
      EWC = ZERO
      DO 260 I = ILO, IHI
         EW = EW + WORK( I+4*N )
         EWC = EWC + WORK( I+5*N )
  260 CONTINUE
*
      GAMMA = COEF*GAMMA - COEF2*( EW**2+EWC**2 ) - COEF5*( EW-EWC )**2
      IF( GAMMA.EQ.ZERO )
     $   GO TO 350
      IF( IT.NE.1 )
     $   BETA = GAMMA / PGAMMA
      T = COEF5*( EWC-THREE*EW )
      TC = COEF5*( EW-THREE*EWC )
*
      CALL DSCAL( NR, BETA, WORK( ILO ), 1 )
      CALL DSCAL( NR, BETA, WORK( ILO+N ), 1 )
*
      CALL DAXPY( NR, COEF, WORK( ILO+4*N ), 1, WORK( ILO+N ), 1 )
      CALL DAXPY( NR, COEF, WORK( ILO+5*N ), 1, WORK( ILO ), 1 )
*
      DO 270 I = ILO, IHI
         WORK( I ) = WORK( I ) + TC
         WORK( I+N ) = WORK( I+N ) + T
  270 CONTINUE
*
*     Apply matrix to vector
*
      DO 300 I = ILO, IHI
         KOUNT = 0
         SUM = ZERO
         DO 290 J = ILO, IHI
            IF( A( I, J ).EQ.ZERO )
     $         GO TO 280
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( J )
  280       CONTINUE
            IF( B( I, J ).EQ.ZERO )
     $         GO TO 290
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( J )
  290    CONTINUE
         WORK( I+2*N ) = DBLE( KOUNT )*WORK( I+N ) + SUM
  300 CONTINUE
*
      DO 330 J = ILO, IHI
         KOUNT = 0
         SUM = ZERO
         DO 320 I = ILO, IHI
            IF( A( I, J ).EQ.ZERO )
     $         GO TO 310
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( I+N )
  310       CONTINUE
            IF( B( I, J ).EQ.ZERO )
     $         GO TO 320
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( I+N )
  320    CONTINUE
         WORK( J+3*N ) = DBLE( KOUNT )*WORK( J ) + SUM
  330 CONTINUE
*
      SUM = DDOT( NR, WORK( ILO+N ), 1, WORK( ILO+2*N ), 1 ) +
     $      DDOT( NR, WORK( ILO ), 1, WORK( ILO+3*N ), 1 )
      ALPHA = GAMMA / SUM
*
*     Determine correction to current iteration
*
      CMAX = ZERO
      DO 340 I = ILO, IHI
         COR = ALPHA*WORK( I+N )
         IF( ABS( COR ).GT.CMAX )
     $      CMAX = ABS( COR )
         LSCALE( I ) = LSCALE( I ) + COR
         COR = ALPHA*WORK( I )
         IF( ABS( COR ).GT.CMAX )
     $      CMAX = ABS( COR )
         RSCALE( I ) = RSCALE( I ) + COR
  340 CONTINUE
      IF( CMAX.LT.HALF )
     $   GO TO 350
*
      CALL DAXPY( NR, -ALPHA, WORK( ILO+2*N ), 1, WORK( ILO+4*N ), 1 )
      CALL DAXPY( NR, -ALPHA, WORK( ILO+3*N ), 1, WORK( ILO+5*N ), 1 )
*
      PGAMMA = GAMMA
      IT = IT + 1
      IF( IT.LE.NRP2 )
     $   GO TO 250
*
*     End generalized conjugate gradient iteration
*
  350 CONTINUE
      SFMIN = DLAMCH( 'S' )
      SFMAX = ONE / SFMIN
      LSFMIN = INT( LOG10( SFMIN ) / BASL+ONE )
      LSFMAX = INT( LOG10( SFMAX ) / BASL )
      DO 360 I = ILO, IHI
         IRAB = IDAMAX( N-ILO+1, A( I, ILO ), LDA )
         RAB = ABS( A( I, IRAB+ILO-1 ) )
         IRAB = IDAMAX( N-ILO+1, B( I, ILO ), LDB )
         RAB = MAX( RAB, ABS( B( I, IRAB+ILO-1 ) ) )
         LRAB = INT( LOG10( RAB+SFMIN ) / BASL+ONE )
         IR = INT(LSCALE( I ) + SIGN( HALF, LSCALE( I ) ))
         IR = MIN( MAX( IR, LSFMIN ), LSFMAX, LSFMAX-LRAB )
         LSCALE( I ) = SCLFAC**IR
         ICAB = IDAMAX( IHI, A( 1, I ), 1 )
         CAB = ABS( A( ICAB, I ) )
         ICAB = IDAMAX( IHI, B( 1, I ), 1 )
         CAB = MAX( CAB, ABS( B( ICAB, I ) ) )
         LCAB = INT( LOG10( CAB+SFMIN ) / BASL+ONE )
         JC = INT(RSCALE( I ) + SIGN( HALF, RSCALE( I ) ))
         JC = MIN( MAX( JC, LSFMIN ), LSFMAX, LSFMAX-LCAB )
         RSCALE( I ) = SCLFAC**JC
  360 CONTINUE
*
*     Row scaling of matrices A and B
*
      DO 370 I = ILO, IHI
         CALL DSCAL( N-ILO+1, LSCALE( I ), A( I, ILO ), LDA )
         CALL DSCAL( N-ILO+1, LSCALE( I ), B( I, ILO ), LDB )
  370 CONTINUE
*
*     Column scaling of matrices A and B
*
      DO 380 J = ILO, IHI
         CALL DSCAL( IHI, RSCALE( J ), A( 1, J ), 1 )
         CALL DSCAL( IHI, RSCALE( J ), B( 1, J ), 1 )
  380 CONTINUE
*
      RETURN
*
*     End of DGGBAL
*
      END
*> \brief <b> DGGEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DGGEV + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggev.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggev.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggev.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,
*                         BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          JOBVL, JOBVR
*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
*      $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
*      $                   VR( LDVR, * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)
*> the generalized eigenvalues, and optionally, the left and/or right
*> generalized eigenvectors.
*>
*> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
*> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
*> singular. It is usually represented as the pair (alpha,beta), as
*> there is a reasonable interpretation for beta=0, and even for both
*> being zero.
*>
*> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
*> of (A,B) satisfies
*>
*>                  A * v(j) = lambda(j) * B * v(j).
*>
*> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
*> of (A,B) satisfies
*>
*>                  u(j)**H * A  = lambda(j) * u(j)**H * B .
*>
*> where u(j)**H is the conjugate-transpose of u(j).
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBVL
*> \verbatim
*>          JOBVL is CHARACTER*1
*>          = 'N':  do not compute the left generalized eigenvectors;
*>          = 'V':  compute the left generalized eigenvectors.
*> \endverbatim
*>
*> \param[in] JOBVR
*> \verbatim
*>          JOBVR is CHARACTER*1
*>          = 'N':  do not compute the right generalized eigenvectors;
*>          = 'V':  compute the right generalized eigenvectors.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A, B, VL, and VR.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA, N)
*>          On entry, the matrix A in the pair (A,B).
*>          On exit, A has been overwritten.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB, N)
*>          On entry, the matrix B in the pair (A,B).
*>          On exit, B has been overwritten.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] ALPHAR
*> \verbatim
*>          ALPHAR is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] ALPHAI
*> \verbatim
*>          ALPHAI is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION array, dimension (N)
*>          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
*>          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
*>          the j-th eigenvalue is real; if positive, then the j-th and
*>          (j+1)-st eigenvalues are a complex conjugate pair, with
*>          ALPHAI(j+1) negative.
*>
*>          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
*>          may easily over- or underflow, and BETA(j) may even be zero.
*>          Thus, the user should avoid naively computing the ratio
*>          alpha/beta.  However, ALPHAR and ALPHAI will be always less
*>          than and usually comparable with norm(A) in magnitude, and
*>          BETA always less than and usually comparable with norm(B).
*> \endverbatim
*>
*> \param[out] VL
*> \verbatim
*>          VL is DOUBLE PRECISION array, dimension (LDVL,N)
*>          If JOBVL = 'V', the left eigenvectors u(j) are stored one
*>          after another in the columns of VL, in the same order as
*>          their eigenvalues. If the j-th eigenvalue is real, then
*>          u(j) = VL(:,j), the j-th column of VL. If the j-th and
*>          (j+1)-th eigenvalues form a complex conjugate pair, then
*>          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
*>          Each eigenvector is scaled so the largest component has
*>          abs(real part)+abs(imag. part)=1.
*>          Not referenced if JOBVL = 'N'.
*> \endverbatim
*>
*> \param[in] LDVL
*> \verbatim
*>          LDVL is INTEGER
*>          The leading dimension of the matrix VL. LDVL >= 1, and
*>          if JOBVL = 'V', LDVL >= N.
*> \endverbatim
*>
*> \param[out] VR
*> \verbatim
*>          VR is DOUBLE PRECISION array, dimension (LDVR,N)
*>          If JOBVR = 'V', the right eigenvectors v(j) are stored one
*>          after another in the columns of VR, in the same order as
*>          their eigenvalues. If the j-th eigenvalue is real, then
*>          v(j) = VR(:,j), the j-th column of VR. If the j-th and
*>          (j+1)-th eigenvalues form a complex conjugate pair, then
*>          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
*>          Each eigenvector is scaled so the largest component has
*>          abs(real part)+abs(imag. part)=1.
*>          Not referenced if JOBVR = 'N'.
*> \endverbatim
*>
*> \param[in] LDVR
*> \verbatim
*>          LDVR is INTEGER
*>          The leading dimension of the matrix VR. LDVR >= 1, and
*>          if JOBVR = 'V', LDVR >= N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,8*N).
*>          For good performance, LWORK must generally be larger.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          = 1,...,N:
*>                The QZ iteration failed.  No eigenvectors have been
*>                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
*>                should be correct for j=INFO+1,...,N.
*>          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
*>                =N+2: error return from DTGEVC.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date April 2012
*
*> \ingroup doubleGEeigen
*
*  =====================================================================
      SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,
     $                  BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine (version 3.4.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     April 2012
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
     $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
     $                   VR( LDVR, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY
      CHARACTER          CHTEMP
      INTEGER            ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO,
     $                   IN, IRIGHT, IROWS, ITAU, IWRK, JC, JR, MAXWRK,
     $                   MINWRK
      DOUBLE PRECISION   ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS,
     $                   SMLNUM, TEMP
*     ..
*     .. Local Arrays ..
      LOGICAL            LDUMMA( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEQRF, DGGBAK, DGGBAL, DGGHRD, DHGEQZ, DLABAD,
     $                   DLACPY,DLASCL, DLASET, DORGQR, DORMQR, DTGEVC,
     $                   XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( LSAME( JOBVL, 'N' ) ) THEN
         IJOBVL = 1
         ILVL = .FALSE.
      ELSE IF( LSAME( JOBVL, 'V' ) ) THEN
         IJOBVL = 2
         ILVL = .TRUE.
      ELSE
         IJOBVL = -1
         ILVL = .FALSE.
      END IF
*
      IF( LSAME( JOBVR, 'N' ) ) THEN
         IJOBVR = 1
         ILVR = .FALSE.
      ELSE IF( LSAME( JOBVR, 'V' ) ) THEN
         IJOBVR = 2
         ILVR = .TRUE.
      ELSE
         IJOBVR = -1
         ILVR = .FALSE.
      END IF
      ILV = ILVL .OR. ILVR
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( IJOBVL.LE.0 ) THEN
         INFO = -1
      ELSE IF( IJOBVR.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVL.LT.1 .OR. ( ILVL .AND. LDVL.LT.N ) ) THEN
         INFO = -12
      ELSE IF( LDVR.LT.1 .OR. ( ILVR .AND. LDVR.LT.N ) ) THEN
         INFO = -14
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV. The workspace is
*       computed assuming ILO = 1 and IHI = N, the worst case.)
*
      IF( INFO.EQ.0 ) THEN
         MINWRK = MAX( 1, 8*N )
         MAXWRK = MAX( 1, N*( 7 +
     $                 ILAENV( 1, 'DGEQRF', ' ', N, 1, N, 0 ) ) )
         MAXWRK = MAX( MAXWRK, N*( 7 +
     $                 ILAENV( 1, 'DORMQR', ' ', N, 1, N, 0 ) ) )
         IF( ILVL ) THEN
            MAXWRK = MAX( MAXWRK, N*( 7 +
     $                 ILAENV( 1, 'DORGQR', ' ', N, 1, N, -1 ) ) )
         END IF
         WORK( 1 ) = MAXWRK
*
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY )
     $      INFO = -16
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = DLANGE( 'M', N, N, A, LDA, WORK )
      ILASCL = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      END IF
      IF( ILASCL )
     $   CALL DLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )
*
*     Scale B if max element outside range [SMLNUM,BIGNUM]
*
      BNRM = DLANGE( 'M', N, N, B, LDB, WORK )
      ILBSCL = .FALSE.
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      ELSE IF( BNRM.GT.BIGNUM ) THEN
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      END IF
      IF( ILBSCL )
     $   CALL DLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )
*
*     Permute the matrices A, B to isolate eigenvalues if possible
*     (Workspace: need 6*N)
*
      ILEFT = 1
      IRIGHT = N + 1
      IWRK = IRIGHT + N
      CALL DGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ),
     $             WORK( IRIGHT ), WORK( IWRK ), IERR )
*
*     Reduce B to triangular form (QR decomposition of B)
*     (Workspace: need N, prefer N*NB)
*
      IROWS = IHI + 1 - ILO
      IF( ILV ) THEN
         ICOLS = N + 1 - ILO
      ELSE
         ICOLS = IROWS
      END IF
      ITAU = IWRK
      IWRK = ITAU + IROWS
      CALL DGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ),
     $             WORK( IWRK ), LWORK+1-IWRK, IERR )
*
*     Apply the orthogonal transformation to matrix A
*     (Workspace: need N, prefer N*NB)
*
      CALL DORMQR( 'L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB,
     $             WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ),
     $             LWORK+1-IWRK, IERR )
*
*     Initialize VL
*     (Workspace: need N, prefer N*NB)
*
      IF( ILVL ) THEN
         CALL DLASET( 'Full', N, N, ZERO, ONE, VL, LDVL )
         IF( IROWS.GT.1 ) THEN
            CALL DLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB,
     $                   VL( ILO+1, ILO ), LDVL )
         END IF
         CALL DORGQR( IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL,
     $                WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      END IF
*
*     Initialize VR
*
      IF( ILVR )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, VR, LDVR )
*
*     Reduce to generalized Hessenberg form
*     (Workspace: none needed)
*
      IF( ILV ) THEN
*
*        Eigenvectors requested -- work on whole matrix.
*
         CALL DGGHRD( JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL,
     $                LDVL, VR, LDVR, IERR )
      ELSE
         CALL DGGHRD( 'N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA,
     $                B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IERR )
      END IF
*
*     Perform QZ algorithm (Compute eigenvalues, and optionally, the
*     Schur forms and Schur vectors)
*     (Workspace: need N)
*
      IWRK = ITAU
      IF( ILV ) THEN
         CHTEMP = 'S'
      ELSE
         CHTEMP = 'E'
      END IF
      CALL DHGEQZ( CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB,
     $             ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR,
     $             WORK( IWRK ), LWORK+1-IWRK, IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.GT.0 .AND. IERR.LE.N ) THEN
            INFO = IERR
         ELSE IF( IERR.GT.N .AND. IERR.LE.2*N ) THEN
            INFO = IERR - N
         ELSE
            INFO = N + 1
         END IF
         GO TO 110
      END IF
*
*     Compute Eigenvectors
*     (Workspace: need 6*N)
*
      IF( ILV ) THEN
         IF( ILVL ) THEN
            IF( ILVR ) THEN
               CHTEMP = 'B'
            ELSE
               CHTEMP = 'L'
            END IF
         ELSE
            CHTEMP = 'R'
         END IF
         CALL DTGEVC( CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL,
     $                VR, LDVR, N, IN, WORK( IWRK ), IERR )
         IF( IERR.NE.0 ) THEN
            INFO = N + 2
            GO TO 110
         END IF
*
*        Undo balancing on VL and VR and normalization
*        (Workspace: none needed)
*
         IF( ILVL ) THEN
            CALL DGGBAK( 'P', 'L', N, ILO, IHI, WORK( ILEFT ),
     $                   WORK( IRIGHT ), N, VL, LDVL, IERR )
            DO 50 JC = 1, N
               IF( ALPHAI( JC ).LT.ZERO )
     $            GO TO 50
               TEMP = ZERO
               IF( ALPHAI( JC ).EQ.ZERO ) THEN
                  DO 10 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) ) )
   10             CONTINUE
               ELSE
                  DO 20 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) )+
     $                      ABS( VL( JR, JC+1 ) ) )
   20             CONTINUE
               END IF
               IF( TEMP.LT.SMLNUM )
     $            GO TO 50
               TEMP = ONE / TEMP
               IF( ALPHAI( JC ).EQ.ZERO ) THEN
                  DO 30 JR = 1, N
                     VL( JR, JC ) = VL( JR, JC )*TEMP
   30             CONTINUE
               ELSE
                  DO 40 JR = 1, N
                     VL( JR, JC ) = VL( JR, JC )*TEMP
                     VL( JR, JC+1 ) = VL( JR, JC+1 )*TEMP
   40             CONTINUE
               END IF
   50       CONTINUE
         END IF
         IF( ILVR ) THEN
            CALL DGGBAK( 'P', 'R', N, ILO, IHI, WORK( ILEFT ),
     $                   WORK( IRIGHT ), N, VR, LDVR, IERR )
            DO 100 JC = 1, N
               IF( ALPHAI( JC ).LT.ZERO )
     $            GO TO 100
               TEMP = ZERO
               IF( ALPHAI( JC ).EQ.ZERO ) THEN
                  DO 60 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VR( JR, JC ) ) )
   60             CONTINUE
               ELSE
                  DO 70 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VR( JR, JC ) )+
     $                      ABS( VR( JR, JC+1 ) ) )
   70             CONTINUE
               END IF
               IF( TEMP.LT.SMLNUM )
     $            GO TO 100
               TEMP = ONE / TEMP
               IF( ALPHAI( JC ).EQ.ZERO ) THEN
                  DO 80 JR = 1, N
                     VR( JR, JC ) = VR( JR, JC )*TEMP
   80             CONTINUE
               ELSE
                  DO 90 JR = 1, N
                     VR( JR, JC ) = VR( JR, JC )*TEMP
                     VR( JR, JC+1 ) = VR( JR, JC+1 )*TEMP
   90             CONTINUE
               END IF
  100       CONTINUE
         END IF
*
*        End of eigenvector calculation
*
      END IF
*
*     Undo scaling if necessary
*
  110 CONTINUE
*
      IF( ILASCL ) THEN
         CALL DLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR )
         CALL DLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR )
      END IF
*
      IF( ILBSCL ) THEN
         CALL DLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
      END IF
*
      WORK( 1 ) = MAXWRK
      RETURN
*
*     End of DGGEV
*
      END
*> \brief \b DGGHRD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DGGHRD + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgghrd.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgghrd.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgghrd.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,
*                          LDQ, Z, LDZ, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          COMPQ, COMPZ
*       INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
*      $                   Z( LDZ, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGGHRD reduces a pair of real matrices (A,B) to generalized upper
*> Hessenberg form using orthogonal transformations, where A is a
*> general matrix and B is upper triangular.  The form of the
*> generalized eigenvalue problem is
*>    A*x = lambda*B*x,
*> and B is typically made upper triangular by computing its QR
*> factorization and moving the orthogonal matrix Q to the left side
*> of the equation.
*>
*> This subroutine simultaneously reduces A to a Hessenberg matrix H:
*>    Q**T*A*Z = H
*> and transforms B to another upper triangular matrix T:
*>    Q**T*B*Z = T
*> in order to reduce the problem to its standard form
*>    H*y = lambda*T*y
*> where y = Z**T*x.
*>
*> The orthogonal matrices Q and Z are determined as products of Givens
*> rotations.  They may either be formed explicitly, or they may be
*> postmultiplied into input matrices Q1 and Z1, so that
*>
*>      Q1 * A * Z1**T = (Q1*Q) * H * (Z1*Z)**T
*>
*>      Q1 * B * Z1**T = (Q1*Q) * T * (Z1*Z)**T
*>
*> If Q1 is the orthogonal matrix from the QR factorization of B in the
*> original equation A*x = lambda*B*x, then DGGHRD reduces the original
*> problem to generalized Hessenberg form.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] COMPQ
*> \verbatim
*>          COMPQ is CHARACTER*1
*>          = 'N': do not compute Q;
*>          = 'I': Q is initialized to the unit matrix, and the
*>                 orthogonal matrix Q is returned;
*>          = 'V': Q must contain an orthogonal matrix Q1 on entry,
*>                 and the product Q1*Q is returned.
*> \endverbatim
*>
*> \param[in] COMPZ
*> \verbatim
*>          COMPZ is CHARACTER*1
*>          = 'N': do not compute Z;
*>          = 'I': Z is initialized to the unit matrix, and the
*>                 orthogonal matrix Z is returned;
*>          = 'V': Z must contain an orthogonal matrix Z1 on entry,
*>                 and the product Z1*Z is returned.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>
*>          ILO and IHI mark the rows and columns of A which are to be
*>          reduced.  It is assumed that A is already upper triangular
*>          in rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI are
*>          normally set by a previous call to DGGBAL; otherwise they
*>          should be set to 1 and N respectively.
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA, N)
*>          On entry, the N-by-N general matrix to be reduced.
*>          On exit, the upper triangle and the first subdiagonal of A
*>          are overwritten with the upper Hessenberg matrix H, and the
*>          rest is set to zero.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB, N)
*>          On entry, the N-by-N upper triangular matrix B.
*>          On exit, the upper triangular matrix T = Q**T B Z.  The
*>          elements below the diagonal are set to zero.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
*>          On entry, if COMPQ = 'V', the orthogonal matrix Q1,
*>          typically from the QR factorization of B.
*>          On exit, if COMPQ='I', the orthogonal matrix Q, and if
*>          COMPQ = 'V', the product Q1*Q.
*>          Not referenced if COMPQ='N'.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q.
*>          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
*>          On entry, if COMPZ = 'V', the orthogonal matrix Z1.
*>          On exit, if COMPZ='I', the orthogonal matrix Z, and if
*>          COMPZ = 'V', the product Z1*Z.
*>          Not referenced if COMPZ='N'.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.
*>          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup doubleOTHERcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  This routine reduces A to Hessenberg and B to triangular form by
*>  an unblocked reduction, as described in _Matrix_Computations_,
*>  by Golub and Van Loan (Johns Hopkins Press.)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,
     $                   LDQ, Z, LDZ, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILQ, ILZ
      INTEGER            ICOMPQ, ICOMPZ, JCOL, JROW
      DOUBLE PRECISION   C, S, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARTG, DLASET, DROT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Decode COMPQ
*
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
*
*     Decode COMPZ
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF( ICOMPQ.LE.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPZ.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -4
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( ( ILQ .AND. LDQ.LT.N ) .OR. LDQ.LT.1 ) THEN
         INFO = -11
      ELSE IF( ( ILZ .AND. LDZ.LT.N ) .OR. LDZ.LT.1 ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGHRD', -INFO )
         RETURN
      END IF
*
*     Initialize Q and Z if desired.
*
      IF( ICOMPQ.EQ.3 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
      IF( ICOMPZ.EQ.3 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
*     Zero out lower triangle of B
*
      DO 20 JCOL = 1, N - 1
         DO 10 JROW = JCOL + 1, N
            B( JROW, JCOL ) = ZERO
   10    CONTINUE
   20 CONTINUE
*
*     Reduce A and B
*
      DO 40 JCOL = ILO, IHI - 2
*
         DO 30 JROW = IHI, JCOL + 2, -1
*
*           Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL)
*
            TEMP = A( JROW-1, JCOL )
            CALL DLARTG( TEMP, A( JROW, JCOL ), C, S,
     $                   A( JROW-1, JCOL ) )
            A( JROW, JCOL ) = ZERO
            CALL DROT( N-JCOL, A( JROW-1, JCOL+1 ), LDA,
     $                 A( JROW, JCOL+1 ), LDA, C, S )
            CALL DROT( N+2-JROW, B( JROW-1, JROW-1 ), LDB,
     $                 B( JROW, JROW-1 ), LDB, C, S )
            IF( ILQ )
     $         CALL DROT( N, Q( 1, JROW-1 ), 1, Q( 1, JROW ), 1, C, S )
*
*           Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1)
*
            TEMP = B( JROW, JROW )
            CALL DLARTG( TEMP, B( JROW, JROW-1 ), C, S,
     $                   B( JROW, JROW ) )
            B( JROW, JROW-1 ) = ZERO
            CALL DROT( IHI, A( 1, JROW ), 1, A( 1, JROW-1 ), 1, C, S )
            CALL DROT( JROW-1, B( 1, JROW ), 1, B( 1, JROW-1 ), 1, C,
     $                 S )
            IF( ILZ )
     $         CALL DROT( N, Z( 1, JROW ), 1, Z( 1, JROW-1 ), 1, C, S )
   30    CONTINUE
   40 CONTINUE
*
      RETURN
*
*     End of DGGHRD
*
      END
*> \brief \b DHGEQZ
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DHGEQZ + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dhgeqz.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dhgeqz.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dhgeqz.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT,
*                          ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK,
*                          LWORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          COMPQ, COMPZ, JOB
*       INTEGER            IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   ALPHAI( * ), ALPHAR( * ), BETA( * ),
*      $                   H( LDH, * ), Q( LDQ, * ), T( LDT, * ),
*      $                   WORK( * ), Z( LDZ, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DHGEQZ computes the eigenvalues of a real matrix pair (H,T),
*> where H is an upper Hessenberg matrix and T is upper triangular,
*> using the double-shift QZ method.
*> Matrix pairs of this type are produced by the reduction to
*> generalized upper Hessenberg form of a real matrix pair (A,B):
*>
*>    A = Q1*H*Z1**T,  B = Q1*T*Z1**T,
*>
*> as computed by DGGHRD.
*>
*> If JOB='S', then the Hessenberg-triangular pair (H,T) is
*> also reduced to generalized Schur form,
*> 
*>    H = Q*S*Z**T,  T = Q*P*Z**T,
*> 
*> where Q and Z are orthogonal matrices, P is an upper triangular
*> matrix, and S is a quasi-triangular matrix with 1-by-1 and 2-by-2
*> diagonal blocks.
*>
*> The 1-by-1 blocks correspond to real eigenvalues of the matrix pair
*> (H,T) and the 2-by-2 blocks correspond to complex conjugate pairs of
*> eigenvalues.
*>
*> Additionally, the 2-by-2 upper triangular diagonal blocks of P
*> corresponding to 2-by-2 blocks of S are reduced to positive diagonal
*> form, i.e., if S(j+1,j) is non-zero, then P(j+1,j) = P(j,j+1) = 0,
*> P(j,j) > 0, and P(j+1,j+1) > 0.
*>
*> Optionally, the orthogonal matrix Q from the generalized Schur
*> factorization may be postmultiplied into an input matrix Q1, and the
*> orthogonal matrix Z may be postmultiplied into an input matrix Z1.
*> If Q1 and Z1 are the orthogonal matrices from DGGHRD that reduced
*> the matrix pair (A,B) to generalized upper Hessenberg form, then the
*> output matrices Q1*Q and Z1*Z are the orthogonal factors from the
*> generalized Schur factorization of (A,B):
*>
*>    A = (Q1*Q)*S*(Z1*Z)**T,  B = (Q1*Q)*P*(Z1*Z)**T.
*> 
*> To avoid overflow, eigenvalues of the matrix pair (H,T) (equivalently,
*> of (A,B)) are computed as a pair of values (alpha,beta), where alpha is
*> complex and beta real.
*> If beta is nonzero, lambda = alpha / beta is an eigenvalue of the
*> generalized nonsymmetric eigenvalue problem (GNEP)
*>    A*x = lambda*B*x
*> and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the
*> alternate form of the GNEP
*>    mu*A*y = B*y.
*> Real eigenvalues can be read directly from the generalized Schur
*> form: 
*>   alpha = S(i,i), beta = P(i,i).
*>
*> Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
*>      Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
*>      pp. 241--256.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOB
*> \verbatim
*>          JOB is CHARACTER*1
*>          = 'E': Compute eigenvalues only;
*>          = 'S': Compute eigenvalues and the Schur form. 
*> \endverbatim
*>
*> \param[in] COMPQ
*> \verbatim
*>          COMPQ is CHARACTER*1
*>          = 'N': Left Schur vectors (Q) are not computed;
*>          = 'I': Q is initialized to the unit matrix and the matrix Q
*>                 of left Schur vectors of (H,T) is returned;
*>          = 'V': Q must contain an orthogonal matrix Q1 on entry and
*>                 the product Q1*Q is returned.
*> \endverbatim
*>
*> \param[in] COMPZ
*> \verbatim
*>          COMPZ is CHARACTER*1
*>          = 'N': Right Schur vectors (Z) are not computed;
*>          = 'I': Z is initialized to the unit matrix and the matrix Z
*>                 of right Schur vectors of (H,T) is returned;
*>          = 'V': Z must contain an orthogonal matrix Z1 on entry and
*>                 the product Z1*Z is returned.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices H, T, Q, and Z.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>          ILO and IHI mark the rows and columns of H which are in
*>          Hessenberg form.  It is assumed that A is already upper
*>          triangular in rows and columns 1:ILO-1 and IHI+1:N.
*>          If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0.
*> \endverbatim
*>
*> \param[in,out] H
*> \verbatim
*>          H is DOUBLE PRECISION array, dimension (LDH, N)
*>          On entry, the N-by-N upper Hessenberg matrix H.
*>          On exit, if JOB = 'S', H contains the upper quasi-triangular
*>          matrix S from the generalized Schur factorization.
*>          If JOB = 'E', the diagonal blocks of H match those of S, but
*>          the rest of H is unspecified.
*> \endverbatim
*>
*> \param[in] LDH
*> \verbatim
*>          LDH is INTEGER
*>          The leading dimension of the array H.  LDH >= max( 1, N ).
*> \endverbatim
*>
*> \param[in,out] T
*> \verbatim
*>          T is DOUBLE PRECISION array, dimension (LDT, N)
*>          On entry, the N-by-N upper triangular matrix T.
*>          On exit, if JOB = 'S', T contains the upper triangular
*>          matrix P from the generalized Schur factorization;
*>          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks of S
*>          are reduced to positive diagonal form, i.e., if H(j+1,j) is
*>          non-zero, then T(j+1,j) = T(j,j+1) = 0, T(j,j) > 0, and
*>          T(j+1,j+1) > 0.
*>          If JOB = 'E', the diagonal blocks of T match those of P, but
*>          the rest of T is unspecified.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= max( 1, N ).
*> \endverbatim
*>
*> \param[out] ALPHAR
*> \verbatim
*>          ALPHAR is DOUBLE PRECISION array, dimension (N)
*>          The real parts of each scalar alpha defining an eigenvalue
*>          of GNEP.
*> \endverbatim
*>
*> \param[out] ALPHAI
*> \verbatim
*>          ALPHAI is DOUBLE PRECISION array, dimension (N)
*>          The imaginary parts of each scalar alpha defining an
*>          eigenvalue of GNEP.
*>          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
*>          positive, then the j-th and (j+1)-st eigenvalues are a
*>          complex conjugate pair, with ALPHAI(j+1) = -ALPHAI(j).
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION array, dimension (N)
*>          The scalars beta that define the eigenvalues of GNEP.
*>          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
*>          beta = BETA(j) represent the j-th eigenvalue of the matrix
*>          pair (A,B), in one of the forms lambda = alpha/beta or
*>          mu = beta/alpha.  Since either lambda or mu may overflow,
*>          they should not, in general, be computed.
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
*>          On entry, if COMPZ = 'V', the orthogonal matrix Q1 used in
*>          the reduction of (A,B) to generalized Hessenberg form.
*>          On exit, if COMPZ = 'I', the orthogonal matrix of left Schur
*>          vectors of (H,T), and if COMPZ = 'V', the orthogonal matrix
*>          of left Schur vectors of (A,B).
*>          Not referenced if COMPZ = 'N'.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q.  LDQ >= 1.
*>          If COMPQ='V' or 'I', then LDQ >= N.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
*>          On entry, if COMPZ = 'V', the orthogonal matrix Z1 used in
*>          the reduction of (A,B) to generalized Hessenberg form.
*>          On exit, if COMPZ = 'I', the orthogonal matrix of
*>          right Schur vectors of (H,T), and if COMPZ = 'V', the
*>          orthogonal matrix of right Schur vectors of (A,B).
*>          Not referenced if COMPZ = 'N'.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.  LDZ >= 1.
*>          If COMPZ='V' or 'I', then LDZ >= N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,N).
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          = 1,...,N: the QZ iteration did not converge.  (H,T) is not
*>                     in Schur form, but ALPHAR(i), ALPHAI(i), and
*>                     BETA(i), i=INFO+1,...,N should be correct.
*>          = N+1,...,2*N: the shift calculation failed.  (H,T) is not
*>                     in Schur form, but ALPHAR(i), ALPHAI(i), and
*>                     BETA(i), i=INFO-N+1,...,N should be correct.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2013
*
*> \ingroup doubleGEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Iteration counters:
*>
*>  JITER  -- counts iterations.
*>  IITER  -- counts iterations run since ILAST was last
*>            changed.  This is therefore reset only when a 1-by-1 or
*>            2-by-2 block deflates off the bottom.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT,
     $                   ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK,
     $                   LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.5.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2013
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ, JOB
      INTEGER            IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   ALPHAI( * ), ALPHAR( * ), BETA( * ),
     $                   H( LDH, * ), Q( LDQ, * ), T( LDT, * ),
     $                   WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*    $                     SAFETY = 1.0E+0 )
      DOUBLE PRECISION   HALF, ZERO, ONE, SAFETY
      PARAMETER          ( HALF = 0.5D+0, ZERO = 0.0D+0, ONE = 1.0D+0,
     $                   SAFETY = 1.0D+2 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILAZR2, ILAZRO, ILPIVT, ILQ, ILSCHR, ILZ,
     $                   LQUERY
      INTEGER            ICOMPQ, ICOMPZ, IFIRST, IFRSTM, IITER, ILAST,
     $                   ILASTM, IN, ISCHUR, ISTART, J, JC, JCH, JITER,
     $                   JR, MAXIT
      DOUBLE PRECISION   A11, A12, A1I, A1R, A21, A22, A2I, A2R, AD11,
     $                   AD11L, AD12, AD12L, AD21, AD21L, AD22, AD22L,
     $                   AD32L, AN, ANORM, ASCALE, ATOL, B11, B1A, B1I,
     $                   B1R, B22, B2A, B2I, B2R, BN, BNORM, BSCALE,
     $                   BTOL, C, C11I, C11R, C12, C21, C22I, C22R, CL,
     $                   CQ, CR, CZ, ESHIFT, S, S1, S1INV, S2, SAFMAX,
     $                   SAFMIN, SCALE, SL, SQI, SQR, SR, SZI, SZR, T1,
     $                   TAU, TEMP, TEMP2, TEMPI, TEMPR, U1, U12, U12L,
     $                   U2, ULP, VS, W11, W12, W21, W22, WABS, WI, WR,
     $                   WR2
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   V( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANHS, DLAPY2, DLAPY3
      EXTERNAL           LSAME, DLAMCH, DLANHS, DLAPY2, DLAPY3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAG2, DLARFG, DLARTG, DLASET, DLASV2, DROT,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Decode JOB, COMPQ, COMPZ
*
      IF( LSAME( JOB, 'E' ) ) THEN
         ILSCHR = .FALSE.
         ISCHUR = 1
      ELSE IF( LSAME( JOB, 'S' ) ) THEN
         ILSCHR = .TRUE.
         ISCHUR = 2
      ELSE
         ISCHUR = 0
      END IF
*
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
*
*     Check Argument Values
*
      INFO = 0
      WORK( 1 ) = MAX( 1, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( ISCHUR.EQ.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPQ.EQ.0 ) THEN
         INFO = -2
      ELSE IF( ICOMPZ.EQ.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -5
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -6
      ELSE IF( LDH.LT.N ) THEN
         INFO = -8
      ELSE IF( LDT.LT.N ) THEN
         INFO = -10
      ELSE IF( LDQ.LT.1 .OR. ( ILQ .AND. LDQ.LT.N ) ) THEN
         INFO = -15
      ELSE IF( LDZ.LT.1 .OR. ( ILZ .AND. LDZ.LT.N ) ) THEN
         INFO = -17
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -19
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DHGEQZ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = DBLE( 1 )
         RETURN
      END IF
*
*     Initialize Q and Z
*
      IF( ICOMPQ.EQ.3 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
      IF( ICOMPZ.EQ.3 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
*
*     Machine Constants
*
      IN = IHI + 1 - ILO
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      ULP = DLAMCH( 'E' )*DLAMCH( 'B' )
      ANORM = DLANHS( 'F', IN, H( ILO, ILO ), LDH, WORK )
      BNORM = DLANHS( 'F', IN, T( ILO, ILO ), LDT, WORK )
      ATOL = MAX( SAFMIN, ULP*ANORM )
      BTOL = MAX( SAFMIN, ULP*BNORM )
      ASCALE = ONE / MAX( SAFMIN, ANORM )
      BSCALE = ONE / MAX( SAFMIN, BNORM )
*
*     Set Eigenvalues IHI+1:N
*
      DO 30 J = IHI + 1, N
         IF( T( J, J ).LT.ZERO ) THEN
            IF( ILSCHR ) THEN
               DO 10 JR = 1, J
                  H( JR, J ) = -H( JR, J )
                  T( JR, J ) = -T( JR, J )
   10          CONTINUE
            ELSE
               H( J, J ) = -H( J, J )
               T( J, J ) = -T( J, J )
            END IF
            IF( ILZ ) THEN
               DO 20 JR = 1, N
                  Z( JR, J ) = -Z( JR, J )
   20          CONTINUE
            END IF
         END IF
         ALPHAR( J ) = H( J, J )
         ALPHAI( J ) = ZERO
         BETA( J ) = T( J, J )
   30 CONTINUE
*
*     If IHI < ILO, skip QZ steps
*
      IF( IHI.LT.ILO )
     $   GO TO 380
*
*     MAIN QZ ITERATION LOOP
*
*     Initialize dynamic indices
*
*     Eigenvalues ILAST+1:N have been found.
*        Column operations modify rows IFRSTM:whatever.
*        Row operations modify columns whatever:ILASTM.
*
*     If only eigenvalues are being computed, then
*        IFRSTM is the row of the last splitting row above row ILAST;
*        this is always at least ILO.
*     IITER counts iterations since the last eigenvalue was found,
*        to tell when to use an extraordinary shift.
*     MAXIT is the maximum number of QZ sweeps allowed.
*
      ILAST = IHI
      IF( ILSCHR ) THEN
         IFRSTM = 1
         ILASTM = N
      ELSE
         IFRSTM = ILO
         ILASTM = IHI
      END IF
      IITER = 0
      ESHIFT = ZERO
      MAXIT = 30*( IHI-ILO+1 )
*
      DO 360 JITER = 1, MAXIT
*
*        Split the matrix if possible.
*
*        Two tests:
*           1: H(j,j-1)=0  or  j=ILO
*           2: T(j,j)=0
*
         IF( ILAST.EQ.ILO ) THEN
*
*           Special case: j=ILAST
*
            GO TO 80
         ELSE
            IF( ABS( H( ILAST, ILAST-1 ) ).LE.ATOL ) THEN
               H( ILAST, ILAST-1 ) = ZERO
               GO TO 80
            END IF
         END IF
*
         IF( ABS( T( ILAST, ILAST ) ).LE.BTOL ) THEN
            T( ILAST, ILAST ) = ZERO
            GO TO 70
         END IF
*
*        General case: j<ILAST
*
         DO 60 J = ILAST - 1, ILO, -1
*
*           Test 1: for H(j,j-1)=0 or j=ILO
*
            IF( J.EQ.ILO ) THEN
               ILAZRO = .TRUE.
            ELSE
               IF( ABS( H( J, J-1 ) ).LE.ATOL ) THEN
                  H( J, J-1 ) = ZERO
                  ILAZRO = .TRUE.
               ELSE
                  ILAZRO = .FALSE.
               END IF
            END IF
*
*           Test 2: for T(j,j)=0
*
            IF( ABS( T( J, J ) ).LT.BTOL ) THEN
               T( J, J ) = ZERO
*
*              Test 1a: Check for 2 consecutive small subdiagonals in A
*
               ILAZR2 = .FALSE.
               IF( .NOT.ILAZRO ) THEN
                  TEMP = ABS( H( J, J-1 ) )
                  TEMP2 = ABS( H( J, J ) )
                  TEMPR = MAX( TEMP, TEMP2 )
                  IF( TEMPR.LT.ONE .AND. TEMPR.NE.ZERO ) THEN
                     TEMP = TEMP / TEMPR
                     TEMP2 = TEMP2 / TEMPR
                  END IF
                  IF( TEMP*( ASCALE*ABS( H( J+1, J ) ) ).LE.TEMP2*
     $                ( ASCALE*ATOL ) )ILAZR2 = .TRUE.
               END IF
*
*              If both tests pass (1 & 2), i.e., the leading diagonal
*              element of B in the block is zero, split a 1x1 block off
*              at the top. (I.e., at the J-th row/column) The leading
*              diagonal element of the remainder can also be zero, so
*              this may have to be done repeatedly.
*
               IF( ILAZRO .OR. ILAZR2 ) THEN
                  DO 40 JCH = J, ILAST - 1
                     TEMP = H( JCH, JCH )
                     CALL DLARTG( TEMP, H( JCH+1, JCH ), C, S,
     $                            H( JCH, JCH ) )
                     H( JCH+1, JCH ) = ZERO
                     CALL DROT( ILASTM-JCH, H( JCH, JCH+1 ), LDH,
     $                          H( JCH+1, JCH+1 ), LDH, C, S )
                     CALL DROT( ILASTM-JCH, T( JCH, JCH+1 ), LDT,
     $                          T( JCH+1, JCH+1 ), LDT, C, S )
                     IF( ILQ )
     $                  CALL DROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1,
     $                             C, S )
                     IF( ILAZR2 )
     $                  H( JCH, JCH-1 ) = H( JCH, JCH-1 )*C
                     ILAZR2 = .FALSE.
                     IF( ABS( T( JCH+1, JCH+1 ) ).GE.BTOL ) THEN
                        IF( JCH+1.GE.ILAST ) THEN
                           GO TO 80
                        ELSE
                           IFIRST = JCH + 1
                           GO TO 110
                        END IF
                     END IF
                     T( JCH+1, JCH+1 ) = ZERO
   40             CONTINUE
                  GO TO 70
               ELSE
*
*                 Only test 2 passed -- chase the zero to T(ILAST,ILAST)
*                 Then process as in the case T(ILAST,ILAST)=0
*
                  DO 50 JCH = J, ILAST - 1
                     TEMP = T( JCH, JCH+1 )
                     CALL DLARTG( TEMP, T( JCH+1, JCH+1 ), C, S,
     $                            T( JCH, JCH+1 ) )
                     T( JCH+1, JCH+1 ) = ZERO
                     IF( JCH.LT.ILASTM-1 )
     $                  CALL DROT( ILASTM-JCH-1, T( JCH, JCH+2 ), LDT,
     $                             T( JCH+1, JCH+2 ), LDT, C, S )
                     CALL DROT( ILASTM-JCH+2, H( JCH, JCH-1 ), LDH,
     $                          H( JCH+1, JCH-1 ), LDH, C, S )
                     IF( ILQ )
     $                  CALL DROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1,
     $                             C, S )
                     TEMP = H( JCH+1, JCH )
                     CALL DLARTG( TEMP, H( JCH+1, JCH-1 ), C, S,
     $                            H( JCH+1, JCH ) )
                     H( JCH+1, JCH-1 ) = ZERO
                     CALL DROT( JCH+1-IFRSTM, H( IFRSTM, JCH ), 1,
     $                          H( IFRSTM, JCH-1 ), 1, C, S )
                     CALL DROT( JCH-IFRSTM, T( IFRSTM, JCH ), 1,
     $                          T( IFRSTM, JCH-1 ), 1, C, S )
                     IF( ILZ )
     $                  CALL DROT( N, Z( 1, JCH ), 1, Z( 1, JCH-1 ), 1,
     $                             C, S )
   50             CONTINUE
                  GO TO 70
               END IF
            ELSE IF( ILAZRO ) THEN
*
*              Only test 1 passed -- work on J:ILAST
*
               IFIRST = J
               GO TO 110
            END IF
*
*           Neither test passed -- try next J
*
   60    CONTINUE
*
*        (Drop-through is "impossible")
*
         INFO = N + 1
         GO TO 420
*
*        T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a
*        1x1 block.
*
   70    CONTINUE
         TEMP = H( ILAST, ILAST )
         CALL DLARTG( TEMP, H( ILAST, ILAST-1 ), C, S,
     $                H( ILAST, ILAST ) )
         H( ILAST, ILAST-1 ) = ZERO
         CALL DROT( ILAST-IFRSTM, H( IFRSTM, ILAST ), 1,
     $              H( IFRSTM, ILAST-1 ), 1, C, S )
         CALL DROT( ILAST-IFRSTM, T( IFRSTM, ILAST ), 1,
     $              T( IFRSTM, ILAST-1 ), 1, C, S )
         IF( ILZ )
     $      CALL DROT( N, Z( 1, ILAST ), 1, Z( 1, ILAST-1 ), 1, C, S )
*
*        H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHAR, ALPHAI,
*                              and BETA
*
   80    CONTINUE
         IF( T( ILAST, ILAST ).LT.ZERO ) THEN
            IF( ILSCHR ) THEN
               DO 90 J = IFRSTM, ILAST
                  H( J, ILAST ) = -H( J, ILAST )
                  T( J, ILAST ) = -T( J, ILAST )
   90          CONTINUE
            ELSE
               H( ILAST, ILAST ) = -H( ILAST, ILAST )
               T( ILAST, ILAST ) = -T( ILAST, ILAST )
            END IF
            IF( ILZ ) THEN
               DO 100 J = 1, N
                  Z( J, ILAST ) = -Z( J, ILAST )
  100          CONTINUE
            END IF
         END IF
         ALPHAR( ILAST ) = H( ILAST, ILAST )
         ALPHAI( ILAST ) = ZERO
         BETA( ILAST ) = T( ILAST, ILAST )
*
*        Go to next block -- exit if finished.
*
         ILAST = ILAST - 1
         IF( ILAST.LT.ILO )
     $      GO TO 380
*
*        Reset counters
*
         IITER = 0
         ESHIFT = ZERO
         IF( .NOT.ILSCHR ) THEN
            ILASTM = ILAST
            IF( IFRSTM.GT.ILAST )
     $         IFRSTM = ILO
         END IF
         GO TO 350
*
*        QZ step
*
*        This iteration only involves rows/columns IFIRST:ILAST. We
*        assume IFIRST < ILAST, and that the diagonal of B is non-zero.
*
  110    CONTINUE
         IITER = IITER + 1
         IF( .NOT.ILSCHR ) THEN
            IFRSTM = IFIRST
         END IF
*
*        Compute single shifts.
*
*        At this point, IFIRST < ILAST, and the diagonal elements of
*        T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
*        magnitude)
*
         IF( ( IITER / 10 )*10.EQ.IITER ) THEN
*
*           Exceptional shift.  Chosen for no particularly good reason.
*           (Single shift only.)
*
            IF( ( DBLE( MAXIT )*SAFMIN )*ABS( H( ILAST, ILAST-1 ) ).LT.
     $          ABS( T( ILAST-1, ILAST-1 ) ) ) THEN
               ESHIFT = H( ILAST, ILAST-1 ) /
     $                  T( ILAST-1, ILAST-1 )
            ELSE
               ESHIFT = ESHIFT + ONE / ( SAFMIN*DBLE( MAXIT ) )
            END IF
            S1 = ONE
            WR = ESHIFT
*
         ELSE
*
*           Shifts based on the generalized eigenvalues of the
*           bottom-right 2x2 block of A and B. The first eigenvalue
*           returned by DLAG2 is the Wilkinson shift (AEP p.512),
*
            CALL DLAG2( H( ILAST-1, ILAST-1 ), LDH,
     $                  T( ILAST-1, ILAST-1 ), LDT, SAFMIN*SAFETY, S1,
     $                  S2, WR, WR2, WI )
*
            IF ( ABS( (WR/S1)*T( ILAST, ILAST ) - H( ILAST, ILAST ) )
     $         .GT. ABS( (WR2/S2)*T( ILAST, ILAST ) 
     $         - H( ILAST, ILAST ) ) ) THEN
               TEMP = WR
               WR = WR2
               WR2 = TEMP
               TEMP = S1
               S1 = S2
               S2 = TEMP
            END IF
            TEMP = MAX( S1, SAFMIN*MAX( ONE, ABS( WR ), ABS( WI ) ) )
            IF( WI.NE.ZERO )
     $         GO TO 200
         END IF
*
*        Fiddle with shift to avoid overflow
*
         TEMP = MIN( ASCALE, ONE )*( HALF*SAFMAX )
         IF( S1.GT.TEMP ) THEN
            SCALE = TEMP / S1
         ELSE
            SCALE = ONE
         END IF
*
         TEMP = MIN( BSCALE, ONE )*( HALF*SAFMAX )
         IF( ABS( WR ).GT.TEMP )
     $      SCALE = MIN( SCALE, TEMP / ABS( WR ) )
         S1 = SCALE*S1
         WR = SCALE*WR
*
*        Now check for two consecutive small subdiagonals.
*
         DO 120 J = ILAST - 1, IFIRST + 1, -1
            ISTART = J
            TEMP = ABS( S1*H( J, J-1 ) )
            TEMP2 = ABS( S1*H( J, J )-WR*T( J, J ) )
            TEMPR = MAX( TEMP, TEMP2 )
            IF( TEMPR.LT.ONE .AND. TEMPR.NE.ZERO ) THEN
               TEMP = TEMP / TEMPR
               TEMP2 = TEMP2 / TEMPR
            END IF
            IF( ABS( ( ASCALE*H( J+1, J ) )*TEMP ).LE.( ASCALE*ATOL )*
     $          TEMP2 )GO TO 130
  120    CONTINUE
*
         ISTART = IFIRST
  130    CONTINUE
*
*        Do an implicit single-shift QZ sweep.
*
*        Initial Q
*
         TEMP = S1*H( ISTART, ISTART ) - WR*T( ISTART, ISTART )
         TEMP2 = S1*H( ISTART+1, ISTART )
         CALL DLARTG( TEMP, TEMP2, C, S, TEMPR )
*
*        Sweep
*
         DO 190 J = ISTART, ILAST - 1
            IF( J.GT.ISTART ) THEN
               TEMP = H( J, J-1 )
               CALL DLARTG( TEMP, H( J+1, J-1 ), C, S, H( J, J-1 ) )
               H( J+1, J-1 ) = ZERO
            END IF
*
            DO 140 JC = J, ILASTM
               TEMP = C*H( J, JC ) + S*H( J+1, JC )
               H( J+1, JC ) = -S*H( J, JC ) + C*H( J+1, JC )
               H( J, JC ) = TEMP
               TEMP2 = C*T( J, JC ) + S*T( J+1, JC )
               T( J+1, JC ) = -S*T( J, JC ) + C*T( J+1, JC )
               T( J, JC ) = TEMP2
  140       CONTINUE
            IF( ILQ ) THEN
               DO 150 JR = 1, N
                  TEMP = C*Q( JR, J ) + S*Q( JR, J+1 )
                  Q( JR, J+1 ) = -S*Q( JR, J ) + C*Q( JR, J+1 )
                  Q( JR, J ) = TEMP
  150          CONTINUE
            END IF
*
            TEMP = T( J+1, J+1 )
            CALL DLARTG( TEMP, T( J+1, J ), C, S, T( J+1, J+1 ) )
            T( J+1, J ) = ZERO
*
            DO 160 JR = IFRSTM, MIN( J+2, ILAST )
               TEMP = C*H( JR, J+1 ) + S*H( JR, J )
               H( JR, J ) = -S*H( JR, J+1 ) + C*H( JR, J )
               H( JR, J+1 ) = TEMP
  160       CONTINUE
            DO 170 JR = IFRSTM, J
               TEMP = C*T( JR, J+1 ) + S*T( JR, J )
               T( JR, J ) = -S*T( JR, J+1 ) + C*T( JR, J )
               T( JR, J+1 ) = TEMP
  170       CONTINUE
            IF( ILZ ) THEN
               DO 180 JR = 1, N
                  TEMP = C*Z( JR, J+1 ) + S*Z( JR, J )
                  Z( JR, J ) = -S*Z( JR, J+1 ) + C*Z( JR, J )
                  Z( JR, J+1 ) = TEMP
  180          CONTINUE
            END IF
  190    CONTINUE
*
         GO TO 350
*
*        Use Francis double-shift
*
*        Note: the Francis double-shift should work with real shifts,
*              but only if the block is at least 3x3.
*              This code may break if this point is reached with
*              a 2x2 block with real eigenvalues.
*
  200    CONTINUE
         IF( IFIRST+1.EQ.ILAST ) THEN
*
*           Special case -- 2x2 block with complex eigenvectors
*
*           Step 1: Standardize, that is, rotate so that
*
*                       ( B11  0  )
*                   B = (         )  with B11 non-negative.
*                       (  0  B22 )
*
            CALL DLASV2( T( ILAST-1, ILAST-1 ), T( ILAST-1, ILAST ),
     $                   T( ILAST, ILAST ), B22, B11, SR, CR, SL, CL )
*
            IF( B11.LT.ZERO ) THEN
               CR = -CR
               SR = -SR
               B11 = -B11
               B22 = -B22
            END IF
*
            CALL DROT( ILASTM+1-IFIRST, H( ILAST-1, ILAST-1 ), LDH,
     $                 H( ILAST, ILAST-1 ), LDH, CL, SL )
            CALL DROT( ILAST+1-IFRSTM, H( IFRSTM, ILAST-1 ), 1,
     $                 H( IFRSTM, ILAST ), 1, CR, SR )
*
            IF( ILAST.LT.ILASTM )
     $         CALL DROT( ILASTM-ILAST, T( ILAST-1, ILAST+1 ), LDT,
     $                    T( ILAST, ILAST+1 ), LDT, CL, SL )
            IF( IFRSTM.LT.ILAST-1 )
     $         CALL DROT( IFIRST-IFRSTM, T( IFRSTM, ILAST-1 ), 1,
     $                    T( IFRSTM, ILAST ), 1, CR, SR )
*
            IF( ILQ )
     $         CALL DROT( N, Q( 1, ILAST-1 ), 1, Q( 1, ILAST ), 1, CL,
     $                    SL )
            IF( ILZ )
     $         CALL DROT( N, Z( 1, ILAST-1 ), 1, Z( 1, ILAST ), 1, CR,
     $                    SR )
*
            T( ILAST-1, ILAST-1 ) = B11
            T( ILAST-1, ILAST ) = ZERO
            T( ILAST, ILAST-1 ) = ZERO
            T( ILAST, ILAST ) = B22
*
*           If B22 is negative, negate column ILAST
*
            IF( B22.LT.ZERO ) THEN
               DO 210 J = IFRSTM, ILAST
                  H( J, ILAST ) = -H( J, ILAST )
                  T( J, ILAST ) = -T( J, ILAST )
  210          CONTINUE
*
               IF( ILZ ) THEN
                  DO 220 J = 1, N
                     Z( J, ILAST ) = -Z( J, ILAST )
  220             CONTINUE
               END IF
               B22 = -B22
            END IF
*
*           Step 2: Compute ALPHAR, ALPHAI, and BETA (see refs.)
*
*           Recompute shift
*
            CALL DLAG2( H( ILAST-1, ILAST-1 ), LDH,
     $                  T( ILAST-1, ILAST-1 ), LDT, SAFMIN*SAFETY, S1,
     $                  TEMP, WR, TEMP2, WI )
*
*           If standardization has perturbed the shift onto real line,
*           do another (real single-shift) QR step.
*
            IF( WI.EQ.ZERO )
     $         GO TO 350
            S1INV = ONE / S1
*
*           Do EISPACK (QZVAL) computation of alpha and beta
*
            A11 = H( ILAST-1, ILAST-1 )
            A21 = H( ILAST, ILAST-1 )
            A12 = H( ILAST-1, ILAST )
            A22 = H( ILAST, ILAST )
*
*           Compute complex Givens rotation on right
*           (Assume some element of C = (sA - wB) > unfl )
*                            __
*           (sA - wB) ( CZ   -SZ )
*                     ( SZ    CZ )
*
            C11R = S1*A11 - WR*B11
            C11I = -WI*B11
            C12 = S1*A12
            C21 = S1*A21
            C22R = S1*A22 - WR*B22
            C22I = -WI*B22
*
            IF( ABS( C11R )+ABS( C11I )+ABS( C12 ).GT.ABS( C21 )+
     $          ABS( C22R )+ABS( C22I ) ) THEN
               T1 = DLAPY3( C12, C11R, C11I )
               CZ = C12 / T1
               SZR = -C11R / T1
               SZI = -C11I / T1
            ELSE
               CZ = DLAPY2( C22R, C22I )
               IF( CZ.LE.SAFMIN ) THEN
                  CZ = ZERO
                  SZR = ONE
                  SZI = ZERO
               ELSE
                  TEMPR = C22R / CZ
                  TEMPI = C22I / CZ
                  T1 = DLAPY2( CZ, C21 )
                  CZ = CZ / T1
                  SZR = -C21*TEMPR / T1
                  SZI = C21*TEMPI / T1
               END IF
            END IF
*
*           Compute Givens rotation on left
*
*           (  CQ   SQ )
*           (  __      )  A or B
*           ( -SQ   CQ )
*
            AN = ABS( A11 ) + ABS( A12 ) + ABS( A21 ) + ABS( A22 )
            BN = ABS( B11 ) + ABS( B22 )
            WABS = ABS( WR ) + ABS( WI )
            IF( S1*AN.GT.WABS*BN ) THEN
               CQ = CZ*B11
               SQR = SZR*B22
               SQI = -SZI*B22
            ELSE
               A1R = CZ*A11 + SZR*A12
               A1I = SZI*A12
               A2R = CZ*A21 + SZR*A22
               A2I = SZI*A22
               CQ = DLAPY2( A1R, A1I )
               IF( CQ.LE.SAFMIN ) THEN
                  CQ = ZERO
                  SQR = ONE
                  SQI = ZERO
               ELSE
                  TEMPR = A1R / CQ
                  TEMPI = A1I / CQ
                  SQR = TEMPR*A2R + TEMPI*A2I
                  SQI = TEMPI*A2R - TEMPR*A2I
               END IF
            END IF
            T1 = DLAPY3( CQ, SQR, SQI )
            CQ = CQ / T1
            SQR = SQR / T1
            SQI = SQI / T1
*
*           Compute diagonal elements of QBZ
*
            TEMPR = SQR*SZR - SQI*SZI
            TEMPI = SQR*SZI + SQI*SZR
            B1R = CQ*CZ*B11 + TEMPR*B22
            B1I = TEMPI*B22
            B1A = DLAPY2( B1R, B1I )
            B2R = CQ*CZ*B22 + TEMPR*B11
            B2I = -TEMPI*B11
            B2A = DLAPY2( B2R, B2I )
*
*           Normalize so beta > 0, and Im( alpha1 ) > 0
*
            BETA( ILAST-1 ) = B1A
            BETA( ILAST ) = B2A
            ALPHAR( ILAST-1 ) = ( WR*B1A )*S1INV
            ALPHAI( ILAST-1 ) = ( WI*B1A )*S1INV
            ALPHAR( ILAST ) = ( WR*B2A )*S1INV
            ALPHAI( ILAST ) = -( WI*B2A )*S1INV
*
*           Step 3: Go to next block -- exit if finished.
*
            ILAST = IFIRST - 1
            IF( ILAST.LT.ILO )
     $         GO TO 380
*
*           Reset counters
*
            IITER = 0
            ESHIFT = ZERO
            IF( .NOT.ILSCHR ) THEN
               ILASTM = ILAST
               IF( IFRSTM.GT.ILAST )
     $            IFRSTM = ILO
            END IF
            GO TO 350
         ELSE
*
*           Usual case: 3x3 or larger block, using Francis implicit
*                       double-shift
*
*                                    2
*           Eigenvalue equation is  w  - c w + d = 0,
*
*                                         -1 2        -1
*           so compute 1st column of  (A B  )  - c A B   + d
*           using the formula in QZIT (from EISPACK)
*
*           We assume that the block is at least 3x3
*
            AD11 = ( ASCALE*H( ILAST-1, ILAST-1 ) ) /
     $             ( BSCALE*T( ILAST-1, ILAST-1 ) )
            AD21 = ( ASCALE*H( ILAST, ILAST-1 ) ) /
     $             ( BSCALE*T( ILAST-1, ILAST-1 ) )
            AD12 = ( ASCALE*H( ILAST-1, ILAST ) ) /
     $             ( BSCALE*T( ILAST, ILAST ) )
            AD22 = ( ASCALE*H( ILAST, ILAST ) ) /
     $             ( BSCALE*T( ILAST, ILAST ) )
            U12 = T( ILAST-1, ILAST ) / T( ILAST, ILAST )
            AD11L = ( ASCALE*H( IFIRST, IFIRST ) ) /
     $              ( BSCALE*T( IFIRST, IFIRST ) )
            AD21L = ( ASCALE*H( IFIRST+1, IFIRST ) ) /
     $              ( BSCALE*T( IFIRST, IFIRST ) )
            AD12L = ( ASCALE*H( IFIRST, IFIRST+1 ) ) /
     $              ( BSCALE*T( IFIRST+1, IFIRST+1 ) )
            AD22L = ( ASCALE*H( IFIRST+1, IFIRST+1 ) ) /
     $              ( BSCALE*T( IFIRST+1, IFIRST+1 ) )
            AD32L = ( ASCALE*H( IFIRST+2, IFIRST+1 ) ) /
     $              ( BSCALE*T( IFIRST+1, IFIRST+1 ) )
            U12L = T( IFIRST, IFIRST+1 ) / T( IFIRST+1, IFIRST+1 )
*
            V( 1 ) = ( AD11-AD11L )*( AD22-AD11L ) - AD12*AD21 +
     $               AD21*U12*AD11L + ( AD12L-AD11L*U12L )*AD21L
            V( 2 ) = ( ( AD22L-AD11L )-AD21L*U12L-( AD11-AD11L )-
     $               ( AD22-AD11L )+AD21*U12 )*AD21L
            V( 3 ) = AD32L*AD21L
*
            ISTART = IFIRST
*
            CALL DLARFG( 3, V( 1 ), V( 2 ), 1, TAU )
            V( 1 ) = ONE
*
*           Sweep
*
            DO 290 J = ISTART, ILAST - 2
*
*              All but last elements: use 3x3 Householder transforms.
*
*              Zero (j-1)st column of A
*
               IF( J.GT.ISTART ) THEN
                  V( 1 ) = H( J, J-1 )
                  V( 2 ) = H( J+1, J-1 )
                  V( 3 ) = H( J+2, J-1 )
*
                  CALL DLARFG( 3, H( J, J-1 ), V( 2 ), 1, TAU )
                  V( 1 ) = ONE
                  H( J+1, J-1 ) = ZERO
                  H( J+2, J-1 ) = ZERO
               END IF
*
               DO 230 JC = J, ILASTM
                  TEMP = TAU*( H( J, JC )+V( 2 )*H( J+1, JC )+V( 3 )*
     $                   H( J+2, JC ) )
                  H( J, JC ) = H( J, JC ) - TEMP
                  H( J+1, JC ) = H( J+1, JC ) - TEMP*V( 2 )
                  H( J+2, JC ) = H( J+2, JC ) - TEMP*V( 3 )
                  TEMP2 = TAU*( T( J, JC )+V( 2 )*T( J+1, JC )+V( 3 )*
     $                    T( J+2, JC ) )
                  T( J, JC ) = T( J, JC ) - TEMP2
                  T( J+1, JC ) = T( J+1, JC ) - TEMP2*V( 2 )
                  T( J+2, JC ) = T( J+2, JC ) - TEMP2*V( 3 )
  230          CONTINUE
               IF( ILQ ) THEN
                  DO 240 JR = 1, N
                     TEMP = TAU*( Q( JR, J )+V( 2 )*Q( JR, J+1 )+V( 3 )*
     $                      Q( JR, J+2 ) )
                     Q( JR, J ) = Q( JR, J ) - TEMP
                     Q( JR, J+1 ) = Q( JR, J+1 ) - TEMP*V( 2 )
                     Q( JR, J+2 ) = Q( JR, J+2 ) - TEMP*V( 3 )
  240             CONTINUE
               END IF
*
*              Zero j-th column of B (see DLAGBC for details)
*
*              Swap rows to pivot
*
               ILPIVT = .FALSE.
               TEMP = MAX( ABS( T( J+1, J+1 ) ), ABS( T( J+1, J+2 ) ) )
               TEMP2 = MAX( ABS( T( J+2, J+1 ) ), ABS( T( J+2, J+2 ) ) )
               IF( MAX( TEMP, TEMP2 ).LT.SAFMIN ) THEN
                  SCALE = ZERO
                  U1 = ONE
                  U2 = ZERO
                  GO TO 250
               ELSE IF( TEMP.GE.TEMP2 ) THEN
                  W11 = T( J+1, J+1 )
                  W21 = T( J+2, J+1 )
                  W12 = T( J+1, J+2 )
                  W22 = T( J+2, J+2 )
                  U1 = T( J+1, J )
                  U2 = T( J+2, J )
               ELSE
                  W21 = T( J+1, J+1 )
                  W11 = T( J+2, J+1 )
                  W22 = T( J+1, J+2 )
                  W12 = T( J+2, J+2 )
                  U2 = T( J+1, J )
                  U1 = T( J+2, J )
               END IF
*
*              Swap columns if nec.
*
               IF( ABS( W12 ).GT.ABS( W11 ) ) THEN
                  ILPIVT = .TRUE.
                  TEMP = W12
                  TEMP2 = W22
                  W12 = W11
                  W22 = W21
                  W11 = TEMP
                  W21 = TEMP2
               END IF
*
*              LU-factor
*
               TEMP = W21 / W11
               U2 = U2 - TEMP*U1
               W22 = W22 - TEMP*W12
               W21 = ZERO
*
*              Compute SCALE
*
               SCALE = ONE
               IF( ABS( W22 ).LT.SAFMIN ) THEN
                  SCALE = ZERO
                  U2 = ONE
                  U1 = -W12 / W11
                  GO TO 250
               END IF
               IF( ABS( W22 ).LT.ABS( U2 ) )
     $            SCALE = ABS( W22 / U2 )
               IF( ABS( W11 ).LT.ABS( U1 ) )
     $            SCALE = MIN( SCALE, ABS( W11 / U1 ) )
*
*              Solve
*
               U2 = ( SCALE*U2 ) / W22
               U1 = ( SCALE*U1-W12*U2 ) / W11
*
  250          CONTINUE
               IF( ILPIVT ) THEN
                  TEMP = U2
                  U2 = U1
                  U1 = TEMP
               END IF
*
*              Compute Householder Vector
*
               T1 = SQRT( SCALE**2+U1**2+U2**2 )
               TAU = ONE + SCALE / T1
               VS = -ONE / ( SCALE+T1 )
               V( 1 ) = ONE
               V( 2 ) = VS*U1
               V( 3 ) = VS*U2
*
*              Apply transformations from the right.
*
               DO 260 JR = IFRSTM, MIN( J+3, ILAST )
                  TEMP = TAU*( H( JR, J )+V( 2 )*H( JR, J+1 )+V( 3 )*
     $                   H( JR, J+2 ) )
                  H( JR, J ) = H( JR, J ) - TEMP
                  H( JR, J+1 ) = H( JR, J+1 ) - TEMP*V( 2 )
                  H( JR, J+2 ) = H( JR, J+2 ) - TEMP*V( 3 )
  260          CONTINUE
               DO 270 JR = IFRSTM, J + 2
                  TEMP = TAU*( T( JR, J )+V( 2 )*T( JR, J+1 )+V( 3 )*
     $                   T( JR, J+2 ) )
                  T( JR, J ) = T( JR, J ) - TEMP
                  T( JR, J+1 ) = T( JR, J+1 ) - TEMP*V( 2 )
                  T( JR, J+2 ) = T( JR, J+2 ) - TEMP*V( 3 )
  270          CONTINUE
               IF( ILZ ) THEN
                  DO 280 JR = 1, N
                     TEMP = TAU*( Z( JR, J )+V( 2 )*Z( JR, J+1 )+V( 3 )*
     $                      Z( JR, J+2 ) )
                     Z( JR, J ) = Z( JR, J ) - TEMP
                     Z( JR, J+1 ) = Z( JR, J+1 ) - TEMP*V( 2 )
                     Z( JR, J+2 ) = Z( JR, J+2 ) - TEMP*V( 3 )
  280             CONTINUE
               END IF
               T( J+1, J ) = ZERO
               T( J+2, J ) = ZERO
  290       CONTINUE
*
*           Last elements: Use Givens rotations
*
*           Rotations from the left
*
            J = ILAST - 1
            TEMP = H( J, J-1 )
            CALL DLARTG( TEMP, H( J+1, J-1 ), C, S, H( J, J-1 ) )
            H( J+1, J-1 ) = ZERO
*
            DO 300 JC = J, ILASTM
               TEMP = C*H( J, JC ) + S*H( J+1, JC )
               H( J+1, JC ) = -S*H( J, JC ) + C*H( J+1, JC )
               H( J, JC ) = TEMP
               TEMP2 = C*T( J, JC ) + S*T( J+1, JC )
               T( J+1, JC ) = -S*T( J, JC ) + C*T( J+1, JC )
               T( J, JC ) = TEMP2
  300       CONTINUE
            IF( ILQ ) THEN
               DO 310 JR = 1, N
                  TEMP = C*Q( JR, J ) + S*Q( JR, J+1 )
                  Q( JR, J+1 ) = -S*Q( JR, J ) + C*Q( JR, J+1 )
                  Q( JR, J ) = TEMP
  310          CONTINUE
            END IF
*
*           Rotations from the right.
*
            TEMP = T( J+1, J+1 )
            CALL DLARTG( TEMP, T( J+1, J ), C, S, T( J+1, J+1 ) )
            T( J+1, J ) = ZERO
*
            DO 320 JR = IFRSTM, ILAST
               TEMP = C*H( JR, J+1 ) + S*H( JR, J )
               H( JR, J ) = -S*H( JR, J+1 ) + C*H( JR, J )
               H( JR, J+1 ) = TEMP
  320       CONTINUE
            DO 330 JR = IFRSTM, ILAST - 1
               TEMP = C*T( JR, J+1 ) + S*T( JR, J )
               T( JR, J ) = -S*T( JR, J+1 ) + C*T( JR, J )
               T( JR, J+1 ) = TEMP
  330       CONTINUE
            IF( ILZ ) THEN
               DO 340 JR = 1, N
                  TEMP = C*Z( JR, J+1 ) + S*Z( JR, J )
                  Z( JR, J ) = -S*Z( JR, J+1 ) + C*Z( JR, J )
                  Z( JR, J+1 ) = TEMP
  340          CONTINUE
            END IF
*
*           End of Double-Shift code
*
         END IF
*
         GO TO 350
*
*        End of iteration loop
*
  350    CONTINUE
  360 CONTINUE
*
*     Drop-through = non-convergence
*
      INFO = ILAST
      GO TO 420
*
*     Successful completion of all QZ steps
*
  380 CONTINUE
*
*     Set Eigenvalues 1:ILO-1
*
      DO 410 J = 1, ILO - 1
         IF( T( J, J ).LT.ZERO ) THEN
            IF( ILSCHR ) THEN
               DO 390 JR = 1, J
                  H( JR, J ) = -H( JR, J )
                  T( JR, J ) = -T( JR, J )
  390          CONTINUE
            ELSE
               H( J, J ) = -H( J, J )
               T( J, J ) = -T( J, J )
            END IF
            IF( ILZ ) THEN
               DO 400 JR = 1, N
                  Z( JR, J ) = -Z( JR, J )
  400          CONTINUE
            END IF
         END IF
         ALPHAR( J ) = H( J, J )
         ALPHAI( J ) = ZERO
         BETA( J ) = T( J, J )
  410 CONTINUE
*
*     Normal Termination
*
      INFO = 0
*
*     Exit (other than argument error) -- return optimal workspace size
*
  420 CONTINUE
      WORK( 1 ) = DBLE( N )
      RETURN
*
*     End of DHGEQZ
*
      END
*> \brief \b DISNAN tests input for NaN.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DISNAN + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/disnan.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/disnan.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/disnan.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       LOGICAL FUNCTION DISNAN( DIN )
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   DIN
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
*> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
*> future.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] DIN
*> \verbatim
*>          DIN is DOUBLE PRECISION
*>          Input to test for NaN.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      LOGICAL FUNCTION DISNAN( DIN )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   DIN
*     ..
*
*  =====================================================================
*
*  .. External Functions ..
      LOGICAL DLAISNAN
      EXTERNAL DLAISNAN
*  ..
*  .. Executable Statements ..
      DISNAN = DLAISNAN(DIN,DIN)
      RETURN
      END
*> \brief \b DLABAD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLABAD + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlabad.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlabad.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlabad.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLABAD( SMALL, LARGE )
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   LARGE, SMALL
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLABAD takes as input the values computed by DLAMCH for underflow and
*> overflow, and returns the square root of each of these values if the
*> log of LARGE is sufficiently large.  This subroutine is intended to
*> identify machines with a large exponent range, such as the Crays, and
*> redefine the underflow and overflow limits to be the square roots of
*> the values computed by DLAMCH.  This subroutine is needed because
*> DLAMCH does not compensate for poor arithmetic in the upper half of
*> the exponent range, as is found on a Cray.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in,out] SMALL
*> \verbatim
*>          SMALL is DOUBLE PRECISION
*>          On entry, the underflow threshold as computed by DLAMCH.
*>          On exit, if LOG10(LARGE) is sufficiently large, the square
*>          root of SMALL, otherwise unchanged.
*> \endverbatim
*>
*> \param[in,out] LARGE
*> \verbatim
*>          LARGE is DOUBLE PRECISION
*>          On entry, the overflow threshold as computed by DLAMCH.
*>          On exit, if LOG10(LARGE) is sufficiently large, the square
*>          root of LARGE, otherwise unchanged.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE DLABAD( SMALL, LARGE )
*
*  -- LAPACK auxiliary routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   LARGE, SMALL
*     ..
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          LOG10, SQRT
*     ..
*     .. Executable Statements ..
*
*     If it looks like we're on a Cray, take the square root of
*     SMALL and LARGE to avoid overflow and underflow problems.
*
      IF( LOG10( LARGE ).GT.2000.D0 ) THEN
         SMALL = SQRT( SMALL )
         LARGE = SQRT( LARGE )
      END IF
*
      RETURN
*
*     End of DLABAD
*
      END
*> \brief \b DLACPY copies all or part of one two-dimensional array to another.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLACPY + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlacpy.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlacpy.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlacpy.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
* 
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            LDA, LDB, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLACPY copies all or part of a two-dimensional matrix A to another
*> matrix B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies the part of the matrix A to be copied to B.
*>          = 'U':      Upper triangular part
*>          = 'L':      Lower triangular part
*>          Otherwise:  All of the matrix A
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          The m by n matrix A.  If UPLO = 'U', only the upper triangle
*>          or trapezoid is accessed; if UPLO = 'L', only the lower
*>          triangle or trapezoid is accessed.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,N)
*>          On exit, B = A in the locations specified by UPLO.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,M).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
*
*     End of DLACPY
*
      END
*> \brief \b DLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLADIV + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dladiv.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dladiv.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dladiv.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLADIV( A, B, C, D, P, Q )
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   A, B, C, D, P, Q
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLADIV performs complex division in  real arithmetic
*>
*>                       a + i*b
*>            p + i*q = ---------
*>                       c + i*d
*>
*> The algorithm is due to Michael Baudin and Robert L. Smith
*> and can be found in the paper
*> "A Robust Complex Division in Scilab"
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is DOUBLE PRECISION
*>          The scalars a, b, c, and d in the above expression.
*> \endverbatim
*>
*> \param[out] P
*> \verbatim
*>          P is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[out] Q
*> \verbatim
*>          Q is DOUBLE PRECISION
*>          The scalars p and q in the above expression.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date January 2013
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE DLADIV( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     January 2013
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   BS
      PARAMETER          ( BS = 2.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
*
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, DD, AB, CD, S, OV, UN, BE, EPS
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLADIV1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
      AA = A
      BB = B
      CC = C
      DD = D
      AB = MAX( ABS(A), ABS(B) )
      CD = MAX( ABS(C), ABS(D) )
      S = 1.0D0
      
      OV = DLAMCH( 'Overflow threshold' )
      UN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Epsilon' )
      BE = BS / (EPS*EPS)
      
      IF( AB >= HALF*OV ) THEN
         AA = HALF * AA
         BB = HALF * BB
         S  = TWO * S
      END IF
      IF( CD >= HALF*OV ) THEN
         CC = HALF * CC
         DD = HALF * DD
         S  = HALF * S
      END IF
      IF( AB <= UN*BS/EPS ) THEN
         AA = AA * BE
         BB = BB * BE
         S  = S / BE
      END IF
      IF( CD <= UN*BS/EPS ) THEN
         CC = CC * BE
         DD = DD * BE
         S  = S * BE
      END IF
      IF( ABS( D ).LE.ABS( C ) ) THEN
         CALL DLADIV1(AA, BB, CC, DD, P, Q)
      ELSE
         CALL DLADIV1(BB, AA, DD, CC, P, Q)
         Q = -Q
      END IF
      P = P * S
      Q = Q * S
*
      RETURN
*
*     End of DLADIV
*
      END

      

      SUBROUTINE DLADIV1( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     January 2013
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
*
*     .. Local Scalars ..
      DOUBLE PRECISION   R, T
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLADIV2
      EXTERNAL           DLADIV2
*     ..
*     .. Executable Statements ..
*
      R = D / C
      T = ONE / (C + D * R)
      P = DLADIV2(A, B, C, D, R, T)
      A = -A
      Q = DLADIV2(B, A, C, D, R, T)
*
      RETURN
*
*     End of DLADIV1
*
      END

      DOUBLE PRECISION FUNCTION DLADIV2( A, B, C, D, R, T )
*
*  -- LAPACK auxiliary routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     January 2013
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, R, T
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*
*     .. Local Scalars ..
      DOUBLE PRECISION   BR
*     ..
*     .. Executable Statements ..
*
      IF( R.NE.ZERO ) THEN
         BR = B * R
         IF( BR.NE.ZERO ) THEN
            DLADIV2 = (A + BR) * T
         ELSE
            DLADIV2 = A * T + (B * T) * R
         END IF
      ELSE
         DLADIV2 = (A + D * (B / C)) * T
      END IF
*
      RETURN
*
*     End of DLADIV12
*
      END
*> \brief \b DLAE2 computes the eigenvalues of a 2-by-2 symmetric matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLAE2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlae2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlae2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlae2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   A, B, C, RT1, RT2
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
*>    [  A   B  ]
*>    [  B   C  ].
*> On return, RT1 is the eigenvalue of larger absolute value, and RT2
*> is the eigenvalue of smaller absolute value.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION
*>          The (1,1) element of the 2-by-2 matrix.
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION
*>          The (1,2) and (2,1) elements of the 2-by-2 matrix.
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is DOUBLE PRECISION
*>          The (2,2) element of the 2-by-2 matrix.
*> \endverbatim
*>
*> \param[out] RT1
*> \verbatim
*>          RT1 is DOUBLE PRECISION
*>          The eigenvalue of larger absolute value.
*> \endverbatim
*>
*> \param[out] RT2
*> \verbatim
*>          RT2 is DOUBLE PRECISION
*>          The eigenvalue of smaller absolute value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  RT1 is accurate to a few ulps barring over/underflow.
*>
*>  RT2 may be inaccurate if there is massive cancellation in the
*>  determinant A*C-B*B; higher precision or correctly rounded or
*>  correctly truncated arithmetic would be needed to compute RT2
*>  accurately in all cases.
*>
*>  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*>  Underflow is harmless if the input data is 0 or exceeds
*>     underflow_threshold / macheps.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, RT1, RT2
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AB, ACMN, ACMX, ADF, DF, RT, SM, TB
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Compute the eigenvalues
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        Includes case AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        Includes case RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN
*
*     End of DLAE2
*
      END
*> \brief \b DLAEV2 computes the eigenvalues and eigenvectors of a 2-by-2 symmetric/Hermitian matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLAEV2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaev2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaev2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaev2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
*>    [  A   B  ]
*>    [  B   C  ].
*> On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
*> eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
*> eigenvector for RT1, giving the decomposition
*>
*>    [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
*>    [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION
*>          The (1,1) element of the 2-by-2 matrix.
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION
*>          The (1,2) element and the conjugate of the (2,1) element of
*>          the 2-by-2 matrix.
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is DOUBLE PRECISION
*>          The (2,2) element of the 2-by-2 matrix.
*> \endverbatim
*>
*> \param[out] RT1
*> \verbatim
*>          RT1 is DOUBLE PRECISION
*>          The eigenvalue of larger absolute value.
*> \endverbatim
*>
*> \param[out] RT2
*> \verbatim
*>          RT2 is DOUBLE PRECISION
*>          The eigenvalue of smaller absolute value.
*> \endverbatim
*>
*> \param[out] CS1
*> \verbatim
*>          CS1 is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[out] SN1
*> \verbatim
*>          SN1 is DOUBLE PRECISION
*>          The vector (CS1, SN1) is a unit right eigenvector for RT1.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  RT1 is accurate to a few ulps barring over/underflow.
*>
*>  RT2 may be inaccurate if there is massive cancellation in the
*>  determinant A*C-B*B; higher precision or correctly rounded or
*>  correctly truncated arithmetic would be needed to compute RT2
*>  accurately in all cases.
*>
*>  CS1 and SN1 are accurate to a few ulps barring over/underflow.
*>
*>  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*>  Underflow is harmless if the input data is 0 or exceeds
*>     underflow_threshold / macheps.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            SGN1, SGN2
      DOUBLE PRECISION   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     $                   TB, TN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Compute the eigenvalues
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        Includes case AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        Includes case RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
*
*     Compute the eigenvector
*
      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
*
*     End of DLAEV2
*
      END
*> \brief \b DLAG2 computes the eigenvalues of a 2-by-2 generalized eigenvalue problem, with scaling as necessary to avoid over-/underflow.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLAG2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlag2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlag2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlag2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAG2( A, LDA, B, LDB, SAFMIN, SCALE1, SCALE2, WR1,
*                         WR2, WI )
* 
*       .. Scalar Arguments ..
*       INTEGER            LDA, LDB
*       DOUBLE PRECISION   SAFMIN, SCALE1, SCALE2, WI, WR1, WR2
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAG2 computes the eigenvalues of a 2 x 2 generalized eigenvalue
*> problem  A - w B, with scaling as necessary to avoid over-/underflow.
*>
*> The scaling factor "s" results in a modified eigenvalue equation
*>
*>     s A - w B
*>
*> where  s  is a non-negative scaling factor chosen so that  w,  w B,
*> and  s A  do not overflow and, if possible, do not underflow, either.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA, 2)
*>          On entry, the 2 x 2 matrix A.  It is assumed that its 1-norm
*>          is less than 1/SAFMIN.  Entries less than
*>          sqrt(SAFMIN)*norm(A) are subject to being treated as zero.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= 2.
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB, 2)
*>          On entry, the 2 x 2 upper triangular matrix B.  It is
*>          assumed that the one-norm of B is less than 1/SAFMIN.  The
*>          diagonals should be at least sqrt(SAFMIN) times the largest
*>          element of B (in absolute value); if a diagonal is smaller
*>          than that, then  +/- sqrt(SAFMIN) will be used instead of
*>          that diagonal.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= 2.
*> \endverbatim
*>
*> \param[in] SAFMIN
*> \verbatim
*>          SAFMIN is DOUBLE PRECISION
*>          The smallest positive number s.t. 1/SAFMIN does not
*>          overflow.  (This should always be DLAMCH('S') -- it is an
*>          argument in order to avoid having to call DLAMCH frequently.)
*> \endverbatim
*>
*> \param[out] SCALE1
*> \verbatim
*>          SCALE1 is DOUBLE PRECISION
*>          A scaling factor used to avoid over-/underflow in the
*>          eigenvalue equation which defines the first eigenvalue.  If
*>          the eigenvalues are complex, then the eigenvalues are
*>          ( WR1  +/-  WI i ) / SCALE1  (which may lie outside the
*>          exponent range of the machine), SCALE1=SCALE2, and SCALE1
*>          will always be positive.  If the eigenvalues are real, then
*>          the first (real) eigenvalue is  WR1 / SCALE1 , but this may
*>          overflow or underflow, and in fact, SCALE1 may be zero or
*>          less than the underflow threshhold if the exact eigenvalue
*>          is sufficiently large.
*> \endverbatim
*>
*> \param[out] SCALE2
*> \verbatim
*>          SCALE2 is DOUBLE PRECISION
*>          A scaling factor used to avoid over-/underflow in the
*>          eigenvalue equation which defines the second eigenvalue.  If
*>          the eigenvalues are complex, then SCALE2=SCALE1.  If the
*>          eigenvalues are real, then the second (real) eigenvalue is
*>          WR2 / SCALE2 , but this may overflow or underflow, and in
*>          fact, SCALE2 may be zero or less than the underflow
*>          threshhold if the exact eigenvalue is sufficiently large.
*> \endverbatim
*>
*> \param[out] WR1
*> \verbatim
*>          WR1 is DOUBLE PRECISION
*>          If the eigenvalue is real, then WR1 is SCALE1 times the
*>          eigenvalue closest to the (2,2) element of A B**(-1).  If the
*>          eigenvalue is complex, then WR1=WR2 is SCALE1 times the real
*>          part of the eigenvalues.
*> \endverbatim
*>
*> \param[out] WR2
*> \verbatim
*>          WR2 is DOUBLE PRECISION
*>          If the eigenvalue is real, then WR2 is SCALE2 times the
*>          other eigenvalue.  If the eigenvalue is complex, then
*>          WR1=WR2 is SCALE1 times the real part of the eigenvalues.
*> \endverbatim
*>
*> \param[out] WI
*> \verbatim
*>          WI is DOUBLE PRECISION
*>          If the eigenvalue is real, then WI is zero.  If the
*>          eigenvalue is complex, then WI is SCALE1 times the imaginary
*>          part of the eigenvalues.  WI will always be non-negative.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE DLAG2( A, LDA, B, LDB, SAFMIN, SCALE1, SCALE2, WR1,
     $                  WR2, WI )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB
      DOUBLE PRECISION   SAFMIN, SCALE1, SCALE2, WI, WR1, WR2
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = ONE / TWO )
      DOUBLE PRECISION   FUZZY1
      PARAMETER          ( FUZZY1 = ONE+1.0D-5 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   A11, A12, A21, A22, ABI22, ANORM, AS11, AS12,
     $                   AS22, ASCALE, B11, B12, B22, BINV11, BINV22,
     $                   BMIN, BNORM, BSCALE, BSIZE, C1, C2, C3, C4, C5,
     $                   DIFF, DISCR, PP, QQ, R, RTMAX, RTMIN, S1, S2,
     $                   SAFMAX, SHIFT, SS, SUM, WABS, WBIG, WDET,
     $                   WSCALE, WSIZE, WSMALL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
      RTMIN = SQRT( SAFMIN )
      RTMAX = ONE / RTMIN
      SAFMAX = ONE / SAFMIN
*
*     Scale A
*
      ANORM = MAX( ABS( A( 1, 1 ) )+ABS( A( 2, 1 ) ),
     $        ABS( A( 1, 2 ) )+ABS( A( 2, 2 ) ), SAFMIN )
      ASCALE = ONE / ANORM
      A11 = ASCALE*A( 1, 1 )
      A21 = ASCALE*A( 2, 1 )
      A12 = ASCALE*A( 1, 2 )
      A22 = ASCALE*A( 2, 2 )
*
*     Perturb B if necessary to insure non-singularity
*
      B11 = B( 1, 1 )
      B12 = B( 1, 2 )
      B22 = B( 2, 2 )
      BMIN = RTMIN*MAX( ABS( B11 ), ABS( B12 ), ABS( B22 ), RTMIN )
      IF( ABS( B11 ).LT.BMIN )
     $   B11 = SIGN( BMIN, B11 )
      IF( ABS( B22 ).LT.BMIN )
     $   B22 = SIGN( BMIN, B22 )
*
*     Scale B
*
      BNORM = MAX( ABS( B11 ), ABS( B12 )+ABS( B22 ), SAFMIN )
      BSIZE = MAX( ABS( B11 ), ABS( B22 ) )
      BSCALE = ONE / BSIZE
      B11 = B11*BSCALE
      B12 = B12*BSCALE
      B22 = B22*BSCALE
*
*     Compute larger eigenvalue by method described by C. van Loan
*
*     ( AS is A shifted by -SHIFT*B )
*
      BINV11 = ONE / B11
      BINV22 = ONE / B22
      S1 = A11*BINV11
      S2 = A22*BINV22
      IF( ABS( S1 ).LE.ABS( S2 ) ) THEN
         AS12 = A12 - S1*B12
         AS22 = A22 - S1*B22
         SS = A21*( BINV11*BINV22 )
         ABI22 = AS22*BINV22 - SS*B12
         PP = HALF*ABI22
         SHIFT = S1
      ELSE
         AS12 = A12 - S2*B12
         AS11 = A11 - S2*B11
         SS = A21*( BINV11*BINV22 )
         ABI22 = -SS*B12
         PP = HALF*( AS11*BINV11+ABI22 )
         SHIFT = S2
      END IF
      QQ = SS*AS12
      IF( ABS( PP*RTMIN ).GE.ONE ) THEN
         DISCR = ( RTMIN*PP )**2 + QQ*SAFMIN
         R = SQRT( ABS( DISCR ) )*RTMAX
      ELSE
         IF( PP**2+ABS( QQ ).LE.SAFMIN ) THEN
            DISCR = ( RTMAX*PP )**2 + QQ*SAFMAX
            R = SQRT( ABS( DISCR ) )*RTMIN
         ELSE
            DISCR = PP**2 + QQ
            R = SQRT( ABS( DISCR ) )
         END IF
      END IF
*
*     Note: the test of R in the following IF is to cover the case when
*           DISCR is small and negative and is flushed to zero during
*           the calculation of R.  On machines which have a consistent
*           flush-to-zero threshhold and handle numbers above that
*           threshhold correctly, it would not be necessary.
*
      IF( DISCR.GE.ZERO .OR. R.EQ.ZERO ) THEN
         SUM = PP + SIGN( R, PP )
         DIFF = PP - SIGN( R, PP )
         WBIG = SHIFT + SUM
*
*        Compute smaller eigenvalue
*
         WSMALL = SHIFT + DIFF
         IF( HALF*ABS( WBIG ).GT.MAX( ABS( WSMALL ), SAFMIN ) ) THEN
            WDET = ( A11*A22-A12*A21 )*( BINV11*BINV22 )
            WSMALL = WDET / WBIG
         END IF
*
*        Choose (real) eigenvalue closest to 2,2 element of A*B**(-1)
*        for WR1.
*
         IF( PP.GT.ABI22 ) THEN
            WR1 = MIN( WBIG, WSMALL )
            WR2 = MAX( WBIG, WSMALL )
         ELSE
            WR1 = MAX( WBIG, WSMALL )
            WR2 = MIN( WBIG, WSMALL )
         END IF
         WI = ZERO
      ELSE
*
*        Complex eigenvalues
*
         WR1 = SHIFT + PP
         WR2 = WR1
         WI = R
      END IF
*
*     Further scaling to avoid underflow and overflow in computing
*     SCALE1 and overflow in computing w*B.
*
*     This scale factor (WSCALE) is bounded from above using C1 and C2,
*     and from below using C3 and C4.
*        C1 implements the condition  s A  must never overflow.
*        C2 implements the condition  w B  must never overflow.
*        C3, with C2,
*           implement the condition that s A - w B must never overflow.
*        C4 implements the condition  s    should not underflow.
*        C5 implements the condition  max(s,|w|) should be at least 2.
*
      C1 = BSIZE*( SAFMIN*MAX( ONE, ASCALE ) )
      C2 = SAFMIN*MAX( ONE, BNORM )
      C3 = BSIZE*SAFMIN
      IF( ASCALE.LE.ONE .AND. BSIZE.LE.ONE ) THEN
         C4 = MIN( ONE, ( ASCALE / SAFMIN )*BSIZE )
      ELSE
         C4 = ONE
      END IF
      IF( ASCALE.LE.ONE .OR. BSIZE.LE.ONE ) THEN
         C5 = MIN( ONE, ASCALE*BSIZE )
      ELSE
         C5 = ONE
      END IF
*
*     Scale first eigenvalue
*
      WABS = ABS( WR1 ) + ABS( WI )
      WSIZE = MAX( SAFMIN, C1, FUZZY1*( WABS*C2+C3 ),
     $        MIN( C4, HALF*MAX( WABS, C5 ) ) )
      IF( WSIZE.NE.ONE ) THEN
         WSCALE = ONE / WSIZE
         IF( WSIZE.GT.ONE ) THEN
            SCALE1 = ( MAX( ASCALE, BSIZE )*WSCALE )*
     $               MIN( ASCALE, BSIZE )
         ELSE
            SCALE1 = ( MIN( ASCALE, BSIZE )*WSCALE )*
     $               MAX( ASCALE, BSIZE )
         END IF
         WR1 = WR1*WSCALE
         IF( WI.NE.ZERO ) THEN
            WI = WI*WSCALE
            WR2 = WR1
            SCALE2 = SCALE1
         END IF
      ELSE
         SCALE1 = ASCALE*BSIZE
         SCALE2 = SCALE1
      END IF
*
*     Scale second eigenvalue (if real)
*
      IF( WI.EQ.ZERO ) THEN
         WSIZE = MAX( SAFMIN, C1, FUZZY1*( ABS( WR2 )*C2+C3 ),
     $           MIN( C4, HALF*MAX( ABS( WR2 ), C5 ) ) )
         IF( WSIZE.NE.ONE ) THEN
            WSCALE = ONE / WSIZE
            IF( WSIZE.GT.ONE ) THEN
               SCALE2 = ( MAX( ASCALE, BSIZE )*WSCALE )*
     $                  MIN( ASCALE, BSIZE )
            ELSE
               SCALE2 = ( MIN( ASCALE, BSIZE )*WSCALE )*
     $                  MAX( ASCALE, BSIZE )
            END IF
            WR2 = WR2*WSCALE
         ELSE
            SCALE2 = ASCALE*BSIZE
         END IF
      END IF
*
*     End of DLAG2
*
      RETURN
      END
*> \brief \b DLAISNAN tests input for NaN by comparing two arguments for inequality.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLAISNAN + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaisnan.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaisnan.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaisnan.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   DIN1, DIN2
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> This routine is not for general use.  It exists solely to avoid
*> over-optimization in DISNAN.
*>
*> DLAISNAN checks for NaNs by comparing its two arguments for
*> inequality.  NaN is the only floating-point value where NaN != NaN
*> returns .TRUE.  To check for NaNs, pass the same variable as both
*> arguments.
*>
*> A compiler must assume that the two arguments are
*> not the same variable, and the test will not be optimized away.
*> Interprocedural or whole-program optimization may delete this
*> test.  The ISNAN functions will be replaced by the correct
*> Fortran 03 intrinsic once the intrinsic is widely available.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] DIN1
*> \verbatim
*>          DIN1 is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] DIN2
*> \verbatim
*>          DIN2 is DOUBLE PRECISION
*>          Two numbers to compare for inequality.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   DIN1, DIN2
*     ..
*
*  =====================================================================
*
*  .. Executable Statements ..
      DLAISNAN = (DIN1.NE.DIN2)
      RETURN
      END
*> \brief \b DLALN2 solves a 1-by-1 or 2-by-2 linear system of equations of the specified form.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLALN2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaln2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaln2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaln2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B,
*                          LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )
* 
*       .. Scalar Arguments ..
*       LOGICAL            LTRANS
*       INTEGER            INFO, LDA, LDB, LDX, NA, NW
*       DOUBLE PRECISION   CA, D1, D2, SCALE, SMIN, WI, WR, XNORM
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), X( LDX, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLALN2 solves a system of the form  (ca A - w D ) X = s B
*> or (ca A**T - w D) X = s B   with possible scaling ("s") and
*> perturbation of A.  (A**T means A-transpose.)
*>
*> A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
*> real diagonal matrix, w is a real or complex value, and X and B are
*> NA x 1 matrices -- real if w is real, complex if w is complex.  NA
*> may be 1 or 2.
*>
*> If w is complex, X and B are represented as NA x 2 matrices,
*> the first column of each being the real part and the second
*> being the imaginary part.
*>
*> "s" is a scaling factor (.LE. 1), computed by DLALN2, which is
*> so chosen that X can be computed without overflow.  X is further
*> scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
*> than overflow.
*>
*> If both singular values of (ca A - w D) are less than SMIN,
*> SMIN*identity will be used instead of (ca A - w D).  If only one
*> singular value is less than SMIN, one element of (ca A - w D) will be
*> perturbed enough to make the smallest singular value roughly SMIN.
*> If both singular values are at least SMIN, (ca A - w D) will not be
*> perturbed.  In any case, the perturbation will be at most some small
*> multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
*> are computed by infinity-norm approximations, and thus will only be
*> correct to a factor of 2 or so.
*>
*> Note: all input quantities are assumed to be smaller than overflow
*> by a reasonable factor.  (See BIGNUM.)
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] LTRANS
*> \verbatim
*>          LTRANS is LOGICAL
*>          =.TRUE.:  A-transpose will be used.
*>          =.FALSE.: A will be used (not transposed.)
*> \endverbatim
*>
*> \param[in] NA
*> \verbatim
*>          NA is INTEGER
*>          The size of the matrix A.  It may (only) be 1 or 2.
*> \endverbatim
*>
*> \param[in] NW
*> \verbatim
*>          NW is INTEGER
*>          1 if "w" is real, 2 if "w" is complex.  It may only be 1
*>          or 2.
*> \endverbatim
*>
*> \param[in] SMIN
*> \verbatim
*>          SMIN is DOUBLE PRECISION
*>          The desired lower bound on the singular values of A.  This
*>          should be a safe distance away from underflow or overflow,
*>          say, between (underflow/machine precision) and  (machine
*>          precision * overflow ).  (See BIGNUM and ULP.)
*> \endverbatim
*>
*> \param[in] CA
*> \verbatim
*>          CA is DOUBLE PRECISION
*>          The coefficient c, which A is multiplied by.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,NA)
*>          The NA x NA matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of A.  It must be at least NA.
*> \endverbatim
*>
*> \param[in] D1
*> \verbatim
*>          D1 is DOUBLE PRECISION
*>          The 1,1 element in the diagonal matrix D.
*> \endverbatim
*>
*> \param[in] D2
*> \verbatim
*>          D2 is DOUBLE PRECISION
*>          The 2,2 element in the diagonal matrix D.  Not used if NW=1.
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,NW)
*>          The NA x NW matrix B (right-hand side).  If NW=2 ("w" is
*>          complex), column 1 contains the real part of B and column 2
*>          contains the imaginary part.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of B.  It must be at least NA.
*> \endverbatim
*>
*> \param[in] WR
*> \verbatim
*>          WR is DOUBLE PRECISION
*>          The real part of the scalar "w".
*> \endverbatim
*>
*> \param[in] WI
*> \verbatim
*>          WI is DOUBLE PRECISION
*>          The imaginary part of the scalar "w".  Not used if NW=1.
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is DOUBLE PRECISION array, dimension (LDX,NW)
*>          The NA x NW matrix X (unknowns), as computed by DLALN2.
*>          If NW=2 ("w" is complex), on exit, column 1 will contain
*>          the real part of X and column 2 will contain the imaginary
*>          part.
*> \endverbatim
*>
*> \param[in] LDX
*> \verbatim
*>          LDX is INTEGER
*>          The leading dimension of X.  It must be at least NA.
*> \endverbatim
*>
*> \param[out] SCALE
*> \verbatim
*>          SCALE is DOUBLE PRECISION
*>          The scale factor that B must be multiplied by to insure
*>          that overflow does not occur when computing X.  Thus,
*>          (ca A - w D) X  will be SCALE*B, not B (ignoring
*>          perturbations of A.)  It will be at most 1.
*> \endverbatim
*>
*> \param[out] XNORM
*> \verbatim
*>          XNORM is DOUBLE PRECISION
*>          The infinity-norm of X, when X is regarded as an NA x NW
*>          real matrix.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          An error flag.  It will be set to zero if no error occurs,
*>          a negative number if an argument is in error, or a positive
*>          number if  ca A - w D  had to be perturbed.
*>          The possible values are:
*>          = 0: No error occurred, and (ca A - w D) did not have to be
*>                 perturbed.
*>          = 1: (ca A - w D) had to be perturbed to make its smallest
*>               (or only) singular value greater than SMIN.
*>          NOTE: In the interests of speed, this routine does not
*>                check the inputs for errors.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE DLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B,
     $                   LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      LOGICAL            LTRANS
      INTEGER            INFO, LDA, LDB, LDX, NA, NW
      DOUBLE PRECISION   CA, D1, D2, SCALE, SMIN, WI, WR, XNORM
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), X( LDX, * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            ICMAX, J
      DOUBLE PRECISION   BBND, BI1, BI2, BIGNUM, BNORM, BR1, BR2, CI21,
     $                   CI22, CMAX, CNORM, CR21, CR22, CSI, CSR, LI21,
     $                   LR21, SMINI, SMLNUM, TEMP, U22ABS, UI11, UI11R,
     $                   UI12, UI12S, UI22, UR11, UR11R, UR12, UR12S,
     $                   UR22, XI1, XI2, XR1, XR2
*     ..
*     .. Local Arrays ..
      LOGICAL            RSWAP( 4 ), ZSWAP( 4 )
      INTEGER            IPIVOT( 4, 4 )
      DOUBLE PRECISION   CI( 2, 2 ), CIV( 4 ), CR( 2, 2 ), CRV( 4 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLADIV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Equivalences ..
      EQUIVALENCE        ( CI( 1, 1 ), CIV( 1 ) ),
     $                   ( CR( 1, 1 ), CRV( 1 ) )
*     ..
*     .. Data statements ..
      DATA               ZSWAP / .FALSE., .FALSE., .TRUE., .TRUE. /
      DATA               RSWAP / .FALSE., .TRUE., .FALSE., .TRUE. /
      DATA               IPIVOT / 1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4,
     $                   3, 2, 1 /
*     ..
*     .. Executable Statements ..
*
*     Compute BIGNUM
*
      SMLNUM = TWO*DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      SMINI = MAX( SMIN, SMLNUM )
*
*     Don't check for input errors
*
      INFO = 0
*
*     Standard Initializations
*
      SCALE = ONE
*
      IF( NA.EQ.1 ) THEN
*
*        1 x 1  (i.e., scalar) system   C X = B
*
         IF( NW.EQ.1 ) THEN
*
*           Real 1x1 system.
*
*           C = ca A - w D
*
            CSR = CA*A( 1, 1 ) - WR*D1
            CNORM = ABS( CSR )
*
*           If | C | < SMINI, use C = SMINI
*
            IF( CNORM.LT.SMINI ) THEN
               CSR = SMINI
               CNORM = SMINI
               INFO = 1
            END IF
*
*           Check scaling for  X = B / C
*
            BNORM = ABS( B( 1, 1 ) )
            IF( CNORM.LT.ONE .AND. BNORM.GT.ONE ) THEN
               IF( BNORM.GT.BIGNUM*CNORM )
     $            SCALE = ONE / BNORM
            END IF
*
*           Compute X
*
            X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / CSR
            XNORM = ABS( X( 1, 1 ) )
         ELSE
*
*           Complex 1x1 system (w is complex)
*
*           C = ca A - w D
*
            CSR = CA*A( 1, 1 ) - WR*D1
            CSI = -WI*D1
            CNORM = ABS( CSR ) + ABS( CSI )
*
*           If | C | < SMINI, use C = SMINI
*
            IF( CNORM.LT.SMINI ) THEN
               CSR = SMINI
               CSI = ZERO
               CNORM = SMINI
               INFO = 1
            END IF
*
*           Check scaling for  X = B / C
*
            BNORM = ABS( B( 1, 1 ) ) + ABS( B( 1, 2 ) )
            IF( CNORM.LT.ONE .AND. BNORM.GT.ONE ) THEN
               IF( BNORM.GT.BIGNUM*CNORM )
     $            SCALE = ONE / BNORM
            END IF
*
*           Compute X
*
            CALL DLADIV( SCALE*B( 1, 1 ), SCALE*B( 1, 2 ), CSR, CSI,
     $                   X( 1, 1 ), X( 1, 2 ) )
            XNORM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) )
         END IF
*
      ELSE
*
*        2x2 System
*
*        Compute the real part of  C = ca A - w D  (or  ca A**T - w D )
*
         CR( 1, 1 ) = CA*A( 1, 1 ) - WR*D1
         CR( 2, 2 ) = CA*A( 2, 2 ) - WR*D2
         IF( LTRANS ) THEN
            CR( 1, 2 ) = CA*A( 2, 1 )
            CR( 2, 1 ) = CA*A( 1, 2 )
         ELSE
            CR( 2, 1 ) = CA*A( 2, 1 )
            CR( 1, 2 ) = CA*A( 1, 2 )
         END IF
*
         IF( NW.EQ.1 ) THEN
*
*           Real 2x2 system  (w is real)
*
*           Find the largest element in C
*
            CMAX = ZERO
            ICMAX = 0
*
            DO 10 J = 1, 4
               IF( ABS( CRV( J ) ).GT.CMAX ) THEN
                  CMAX = ABS( CRV( J ) )
                  ICMAX = J
               END IF
   10       CONTINUE
*
*           If norm(C) < SMINI, use SMINI*identity.
*
            IF( CMAX.LT.SMINI ) THEN
               BNORM = MAX( ABS( B( 1, 1 ) ), ABS( B( 2, 1 ) ) )
               IF( SMINI.LT.ONE .AND. BNORM.GT.ONE ) THEN
                  IF( BNORM.GT.BIGNUM*SMINI )
     $               SCALE = ONE / BNORM
               END IF
               TEMP = SCALE / SMINI
               X( 1, 1 ) = TEMP*B( 1, 1 )
               X( 2, 1 ) = TEMP*B( 2, 1 )
               XNORM = TEMP*BNORM
               INFO = 1
               RETURN
            END IF
*
*           Gaussian elimination with complete pivoting.
*
            UR11 = CRV( ICMAX )
            CR21 = CRV( IPIVOT( 2, ICMAX ) )
            UR12 = CRV( IPIVOT( 3, ICMAX ) )
            CR22 = CRV( IPIVOT( 4, ICMAX ) )
            UR11R = ONE / UR11
            LR21 = UR11R*CR21
            UR22 = CR22 - UR12*LR21
*
*           If smaller pivot < SMINI, use SMINI
*
            IF( ABS( UR22 ).LT.SMINI ) THEN
               UR22 = SMINI
               INFO = 1
            END IF
            IF( RSWAP( ICMAX ) ) THEN
               BR1 = B( 2, 1 )
               BR2 = B( 1, 1 )
            ELSE
               BR1 = B( 1, 1 )
               BR2 = B( 2, 1 )
            END IF
            BR2 = BR2 - LR21*BR1
            BBND = MAX( ABS( BR1*( UR22*UR11R ) ), ABS( BR2 ) )
            IF( BBND.GT.ONE .AND. ABS( UR22 ).LT.ONE ) THEN
               IF( BBND.GE.BIGNUM*ABS( UR22 ) )
     $            SCALE = ONE / BBND
            END IF
*
            XR2 = ( BR2*SCALE ) / UR22
            XR1 = ( SCALE*BR1 )*UR11R - XR2*( UR11R*UR12 )
            IF( ZSWAP( ICMAX ) ) THEN
               X( 1, 1 ) = XR2
               X( 2, 1 ) = XR1
            ELSE
               X( 1, 1 ) = XR1
               X( 2, 1 ) = XR2
            END IF
            XNORM = MAX( ABS( XR1 ), ABS( XR2 ) )
*
*           Further scaling if  norm(A) norm(X) > overflow
*
            IF( XNORM.GT.ONE .AND. CMAX.GT.ONE ) THEN
               IF( XNORM.GT.BIGNUM / CMAX ) THEN
                  TEMP = CMAX / BIGNUM
                  X( 1, 1 ) = TEMP*X( 1, 1 )
                  X( 2, 1 ) = TEMP*X( 2, 1 )
                  XNORM = TEMP*XNORM
                  SCALE = TEMP*SCALE
               END IF
            END IF
         ELSE
*
*           Complex 2x2 system  (w is complex)
*
*           Find the largest element in C
*
            CI( 1, 1 ) = -WI*D1
            CI( 2, 1 ) = ZERO
            CI( 1, 2 ) = ZERO
            CI( 2, 2 ) = -WI*D2
            CMAX = ZERO
            ICMAX = 0
*
            DO 20 J = 1, 4
               IF( ABS( CRV( J ) )+ABS( CIV( J ) ).GT.CMAX ) THEN
                  CMAX = ABS( CRV( J ) ) + ABS( CIV( J ) )
                  ICMAX = J
               END IF
   20       CONTINUE
*
*           If norm(C) < SMINI, use SMINI*identity.
*
            IF( CMAX.LT.SMINI ) THEN
               BNORM = MAX( ABS( B( 1, 1 ) )+ABS( B( 1, 2 ) ),
     $                 ABS( B( 2, 1 ) )+ABS( B( 2, 2 ) ) )
               IF( SMINI.LT.ONE .AND. BNORM.GT.ONE ) THEN
                  IF( BNORM.GT.BIGNUM*SMINI )
     $               SCALE = ONE / BNORM
               END IF
               TEMP = SCALE / SMINI
               X( 1, 1 ) = TEMP*B( 1, 1 )
               X( 2, 1 ) = TEMP*B( 2, 1 )
               X( 1, 2 ) = TEMP*B( 1, 2 )
               X( 2, 2 ) = TEMP*B( 2, 2 )
               XNORM = TEMP*BNORM
               INFO = 1
               RETURN
            END IF
*
*           Gaussian elimination with complete pivoting.
*
            UR11 = CRV( ICMAX )
            UI11 = CIV( ICMAX )
            CR21 = CRV( IPIVOT( 2, ICMAX ) )
            CI21 = CIV( IPIVOT( 2, ICMAX ) )
            UR12 = CRV( IPIVOT( 3, ICMAX ) )
            UI12 = CIV( IPIVOT( 3, ICMAX ) )
            CR22 = CRV( IPIVOT( 4, ICMAX ) )
            CI22 = CIV( IPIVOT( 4, ICMAX ) )
            IF( ICMAX.EQ.1 .OR. ICMAX.EQ.4 ) THEN
*
*              Code when off-diagonals of pivoted C are real
*
               IF( ABS( UR11 ).GT.ABS( UI11 ) ) THEN
                  TEMP = UI11 / UR11
                  UR11R = ONE / ( UR11*( ONE+TEMP**2 ) )
                  UI11R = -TEMP*UR11R
               ELSE
                  TEMP = UR11 / UI11
                  UI11R = -ONE / ( UI11*( ONE+TEMP**2 ) )
                  UR11R = -TEMP*UI11R
               END IF
               LR21 = CR21*UR11R
               LI21 = CR21*UI11R
               UR12S = UR12*UR11R
               UI12S = UR12*UI11R
               UR22 = CR22 - UR12*LR21
               UI22 = CI22 - UR12*LI21
            ELSE
*
*              Code when diagonals of pivoted C are real
*
               UR11R = ONE / UR11
               UI11R = ZERO
               LR21 = CR21*UR11R
               LI21 = CI21*UR11R
               UR12S = UR12*UR11R
               UI12S = UI12*UR11R
               UR22 = CR22 - UR12*LR21 + UI12*LI21
               UI22 = -UR12*LI21 - UI12*LR21
            END IF
            U22ABS = ABS( UR22 ) + ABS( UI22 )
*
*           If smaller pivot < SMINI, use SMINI
*
            IF( U22ABS.LT.SMINI ) THEN
               UR22 = SMINI
               UI22 = ZERO
               INFO = 1
            END IF
            IF( RSWAP( ICMAX ) ) THEN
               BR2 = B( 1, 1 )
               BR1 = B( 2, 1 )
               BI2 = B( 1, 2 )
               BI1 = B( 2, 2 )
            ELSE
               BR1 = B( 1, 1 )
               BR2 = B( 2, 1 )
               BI1 = B( 1, 2 )
               BI2 = B( 2, 2 )
            END IF
            BR2 = BR2 - LR21*BR1 + LI21*BI1
            BI2 = BI2 - LI21*BR1 - LR21*BI1
            BBND = MAX( ( ABS( BR1 )+ABS( BI1 ) )*
     $             ( U22ABS*( ABS( UR11R )+ABS( UI11R ) ) ),
     $             ABS( BR2 )+ABS( BI2 ) )
            IF( BBND.GT.ONE .AND. U22ABS.LT.ONE ) THEN
               IF( BBND.GE.BIGNUM*U22ABS ) THEN
                  SCALE = ONE / BBND
                  BR1 = SCALE*BR1
                  BI1 = SCALE*BI1
                  BR2 = SCALE*BR2
                  BI2 = SCALE*BI2
               END IF
            END IF
*
            CALL DLADIV( BR2, BI2, UR22, UI22, XR2, XI2 )
            XR1 = UR11R*BR1 - UI11R*BI1 - UR12S*XR2 + UI12S*XI2
            XI1 = UI11R*BR1 + UR11R*BI1 - UI12S*XR2 - UR12S*XI2
            IF( ZSWAP( ICMAX ) ) THEN
               X( 1, 1 ) = XR2
               X( 2, 1 ) = XR1
               X( 1, 2 ) = XI2
               X( 2, 2 ) = XI1
            ELSE
               X( 1, 1 ) = XR1
               X( 2, 1 ) = XR2
               X( 1, 2 ) = XI1
               X( 2, 2 ) = XI2
            END IF
            XNORM = MAX( ABS( XR1 )+ABS( XI1 ), ABS( XR2 )+ABS( XI2 ) )
*
*           Further scaling if  norm(A) norm(X) > overflow
*
            IF( XNORM.GT.ONE .AND. CMAX.GT.ONE ) THEN
               IF( XNORM.GT.BIGNUM / CMAX ) THEN
                  TEMP = CMAX / BIGNUM
                  X( 1, 1 ) = TEMP*X( 1, 1 )
                  X( 2, 1 ) = TEMP*X( 2, 1 )
                  X( 1, 2 ) = TEMP*X( 1, 2 )
                  X( 2, 2 ) = TEMP*X( 2, 2 )
                  XNORM = TEMP*XNORM
                  SCALE = TEMP*SCALE
               END IF
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DLALN2
*
      END
*> \brief \b DLAMCH
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAMCH determines double precision machine parameters.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] CMACH
*> \verbatim
*>          Specifies the value to be returned by DLAMCH:
*>          = 'E' or 'e',   DLAMCH := eps
*>          = 'S' or 's ,   DLAMCH := sfmin
*>          = 'B' or 'b',   DLAMCH := base
*>          = 'P' or 'p',   DLAMCH := eps*base
*>          = 'N' or 'n',   DLAMCH := t
*>          = 'R' or 'r',   DLAMCH := rnd
*>          = 'M' or 'm',   DLAMCH := emin
*>          = 'U' or 'u',   DLAMCH := rmin
*>          = 'L' or 'l',   DLAMCH := emax
*>          = 'O' or 'o',   DLAMCH := rmax
*>          where
*>          eps   = relative machine precision
*>          sfmin = safe minimum, such that 1/sfmin does not overflow
*>          base  = base of the machine
*>          prec  = eps*base
*>          t     = number of (base) digits in the mantissa
*>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*>          emin  = minimum exponent before (gradual) underflow
*>          rmin  = underflow threshold - base**(emin-1)
*>          emax  = largest exponent before overflow
*>          rmax  = overflow threshold  - (base**emax)*(1-eps)
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   RND, EPS, SFMIN, SMALL, RMACH
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DIGITS, EPSILON, HUGE, MAXEXPONENT,
     $                   MINEXPONENT, RADIX, TINY
*     ..
*     .. Executable Statements ..
*
*
*     Assume rounding, not chopping. Always.
*
      RND = ONE
*
      IF( ONE.EQ.RND ) THEN
         EPS = EPSILON(ZERO) * 0.5
      ELSE
         EPS = EPSILON(ZERO)
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         SFMIN = TINY(ZERO)
         SMALL = ONE / HUGE(ZERO)
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = RADIX(ZERO)
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = EPS * RADIX(ZERO)
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = DIGITS(ZERO)
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = MINEXPONENT(ZERO)
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = tiny(zero)
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = MAXEXPONENT(ZERO)
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = HUGE(ZERO)
      ELSE
         RMACH = ZERO
      END IF
*
      DLAMCH = RMACH
      RETURN
*
*     End of DLAMCH
*
      END
************************************************************************
*> \brief \b DLAMC3
*> \details
*> \b Purpose:
*> \verbatim
*> DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*> the addition of  A  and  B ,  for use in situations where optimizers
*> might hold one of these in a register.
*> \endverbatim
*> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*> \date November 2015
*> \ingroup auxOTHERauxiliary
*>
*> \param[in] A
*> \verbatim
*>          A is a DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is a DOUBLE PRECISION
*>          The values A and B.
*> \endverbatim
*>
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.6.0) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2010
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
*> \brief \b DLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general rectangular matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLANGE + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlange.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlange.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlange.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
* 
*       .. Scalar Arguments ..
*       CHARACTER          NORM
*       INTEGER            LDA, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLANGE  returns the value of the one norm,  or the Frobenius norm, or
*> the  infinity norm,  or the  element of  largest absolute value  of a
*> real matrix A.
*> \endverbatim
*>
*> \return DLANGE
*> \verbatim
*>
*>    DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*>             (
*>             ( norm1(A),         NORM = '1', 'O' or 'o'
*>             (
*>             ( normI(A),         NORM = 'I' or 'i'
*>             (
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] NORM
*> \verbatim
*>          NORM is CHARACTER*1
*>          Specifies the value to be returned in DLANGE as described
*>          above.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.  When M = 0,
*>          DLANGE is set to zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.  When N = 0,
*>          DLANGE is set to zero.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          The m by n matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(M,1).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
*>          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
*>          referenced.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleGEauxiliary
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE, TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASSQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( MIN( M, N ).EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, M
               TEMP = ABS( A( I, J ) )
               IF( VALUE.LT.TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
*
*        Find norm1(A).
*
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, M
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            IF( VALUE.LT.SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
*
*        Find normI(A).
*
         DO 50 I = 1, M
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, M
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, M
            TEMP = WORK( I )
            IF( VALUE.LT.TEMP .OR. DISNAN( TEMP ) ) VALUE = TEMP
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL DLASSQ( M, A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      DLANGE = VALUE
      RETURN
*
*     End of DLANGE
*
      END
*> \brief \b DLANHS returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of an upper Hessenberg matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLANHS + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanhs.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanhs.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanhs.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DLANHS( NORM, N, A, LDA, WORK )
* 
*       .. Scalar Arguments ..
*       CHARACTER          NORM
*       INTEGER            LDA, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLANHS  returns the value of the one norm,  or the Frobenius norm, or
*> the  infinity norm,  or the  element of  largest absolute value  of a
*> Hessenberg matrix A.
*> \endverbatim
*>
*> \return DLANHS
*> \verbatim
*>
*>    DLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*>             (
*>             ( norm1(A),         NORM = '1', 'O' or 'o'
*>             (
*>             ( normI(A),         NORM = 'I' or 'i'
*>             (
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] NORM
*> \verbatim
*>          NORM is CHARACTER*1
*>          Specifies the value to be returned in DLANHS as described
*>          above.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.  When N = 0, DLANHS is
*>          set to zero.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          The n by n upper Hessenberg matrix A; the part of A below the
*>          first sub-diagonal is not referenced.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(N,1).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
*>          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
*>          referenced.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleOTHERauxiliary
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION DLANHS( NORM, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASSQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, MIN( N, J+1 )
               SUM = ABS( A( I, J ) )
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
*
*        Find norm1(A).
*
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, MIN( N, J+1 )
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
*
*        Find normI(A).
*
         DO 50 I = 1, N
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, MIN( N, J+1 )
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, N
            SUM = WORK( I )
            IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL DLASSQ( MIN( N, J+1 ), A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      DLANHS = VALUE
      RETURN
*
*     End of DLANHS
*
      END
*> \brief \b DLANST returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a real symmetric tridiagonal matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLANST + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanst.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanst.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanst.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )
* 
*       .. Scalar Arguments ..
*       CHARACTER          NORM
*       INTEGER            N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   D( * ), E( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLANST  returns the value of the one norm,  or the Frobenius norm, or
*> the  infinity norm,  or the  element of  largest absolute value  of a
*> real symmetric tridiagonal matrix A.
*> \endverbatim
*>
*> \return DLANST
*> \verbatim
*>
*>    DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*>             (
*>             ( norm1(A),         NORM = '1', 'O' or 'o'
*>             (
*>             ( normI(A),         NORM = 'I' or 'i'
*>             (
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] NORM
*> \verbatim
*>          NORM is CHARACTER*1
*>          Specifies the value to be returned in DLANST as described
*>          above.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.  When N = 0, DLANST is
*>          set to zero.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (N)
*>          The diagonal elements of A.
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (N-1)
*>          The (n-1) sub-diagonal or super-diagonal elements of A.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   ANORM, SCALE, SUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASSQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         ANORM = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         ANORM = ABS( D( N ) )
         DO 10 I = 1, N - 1
            SUM = ABS( D( I ) )
            IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
            SUM = ABS( E( I ) )
            IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
   10    CONTINUE
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' .OR.
     $         LSAME( NORM, 'I' ) ) THEN
*
*        Find norm1(A).
*
         IF( N.EQ.1 ) THEN
            ANORM = ABS( D( 1 ) )
         ELSE
            ANORM = ABS( D( 1 ) )+ABS( E( 1 ) )
            SUM = ABS( E( N-1 ) )+ABS( D( N ) )
            IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
            DO 20 I = 2, N - 1
               SUM = ABS( D( I ) )+ABS( E( I ) )+ABS( E( I-1 ) )
               IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
   20       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         IF( N.GT.1 ) THEN
            CALL DLASSQ( N-1, E, 1, SCALE, SUM )
            SUM = 2*SUM
         END IF
         CALL DLASSQ( N, D, 1, SCALE, SUM )
         ANORM = SCALE*SQRT( SUM )
      END IF
*
      DLANST = ANORM
      RETURN
*
*     End of DLANST
*
      END
*> \brief \b DLANSY returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a real symmetric matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLANSY + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlansy.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlansy.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlansy.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK )
* 
*       .. Scalar Arguments ..
*       CHARACTER          NORM, UPLO
*       INTEGER            LDA, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLANSY  returns the value of the one norm,  or the Frobenius norm, or
*> the  infinity norm,  or the  element of  largest absolute value  of a
*> real symmetric matrix A.
*> \endverbatim
*>
*> \return DLANSY
*> \verbatim
*>
*>    DLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*>             (
*>             ( norm1(A),         NORM = '1', 'O' or 'o'
*>             (
*>             ( normI(A),         NORM = 'I' or 'i'
*>             (
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] NORM
*> \verbatim
*>          NORM is CHARACTER*1
*>          Specifies the value to be returned in DLANSY as described
*>          above.
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          symmetric matrix A is to be referenced.
*>          = 'U':  Upper triangular part of A is referenced
*>          = 'L':  Lower triangular part of A is referenced
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.  When N = 0, DLANSY is
*>          set to zero.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          The symmetric matrix A.  If UPLO = 'U', the leading n by n
*>          upper triangular part of A contains the upper triangular part
*>          of the matrix A, and the strictly lower triangular part of A
*>          is not referenced.  If UPLO = 'L', the leading n by n lower
*>          triangular part of A contains the lower triangular part of
*>          the matrix A, and the strictly upper triangular part of A is
*>          not referenced.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(N,1).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
*>          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
*>          WORK is not referenced.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup doubleSYauxiliary
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASSQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           LSAME, DISNAN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, N
               DO 10 I = 1, J
                  SUM = ABS( A( I, J ) )
                  IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               DO 30 I = J, N
                  SUM = ABS( A( I, J ) )
                  IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR.
     $         ( NORM.EQ.'1' ) ) THEN
*
*        Find normI(A) ( = norm1(A), since A is symmetric).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( A( J, J ) )
   60       CONTINUE
            DO 70 I = 1, N
               SUM = WORK( I )
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( A( J, J ) )
               DO 90 I = J + 1, N
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL DLASSQ( J-1, A( 1, J ), 1, SCALE, SUM )
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL DLASSQ( N-J, A( J+1, J ), 1, SCALE, SUM )
  120       CONTINUE
         END IF
         SUM = 2*SUM
         CALL DLASSQ( N, A, LDA+1, SCALE, SUM )
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      DLANSY = VALUE
      RETURN
*
*     End of DLANSY
*
      END
*> \brief \b DLAPY2 returns sqrt(x2+y2).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLAPY2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   X, Y
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*> overflow.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] X
*> \verbatim
*>          X is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is DOUBLE PRECISION
*>          X and Y specify the values x and y.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         DLAPY2 = W
      ELSE
         DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
*
*     End of DLAPY2
*
      END
*> \brief \b DLAPY3 returns sqrt(x2+y2+z2).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLAPY3 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy3.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy3.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy3.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z )
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   X, Y, Z
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
*> unnecessary overflow.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] X
*> \verbatim
*>          X is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] Z
*> \verbatim
*>          Z is DOUBLE PRECISION
*>          X, Y and Z specify the values x, y and z.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y, Z
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, ZABS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO ) THEN
*     W can be zero for max(0,nan,0)
*     adding all three entries together will make sure
*     NaN will not disappear.
         DLAPY3 =  XABS + YABS + ZABS
      ELSE
         DLAPY3 = W*SQRT( ( XABS / W )**2+( YABS / W )**2+
     $            ( ZABS / W )**2 )
      END IF
      RETURN
*
*     End of DLAPY3
*
      END
*> \brief \b DLARFB applies a block reflector or its transpose to a general rectangular matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLARFB + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfb.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfb.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfb.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
*                          T, LDT, C, LDC, WORK, LDWORK )
* 
*       .. Scalar Arguments ..
*       CHARACTER          DIRECT, SIDE, STOREV, TRANS
*       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),
*      $                   WORK( LDWORK, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLARFB applies a real block reflector H or its transpose H**T to a
*> real m by n matrix C, from either the left or the right.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply H or H**T from the Left
*>          = 'R': apply H or H**T from the Right
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N': apply H (No transpose)
*>          = 'T': apply H**T (Transpose)
*> \endverbatim
*>
*> \param[in] DIRECT
*> \verbatim
*>          DIRECT is CHARACTER*1
*>          Indicates how H is formed from a product of elementary
*>          reflectors
*>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*> \endverbatim
*>
*> \param[in] STOREV
*> \verbatim
*>          STOREV is CHARACTER*1
*>          Indicates how the vectors which define the elementary
*>          reflectors are stored:
*>          = 'C': Columnwise
*>          = 'R': Rowwise
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The order of the matrix T (= the number of elementary
*>          reflectors whose product defines the block reflector).
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is DOUBLE PRECISION array, dimension
*>                                (LDV,K) if STOREV = 'C'
*>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
*>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
*>          The matrix V. See Further Details.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
*>          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
*>          if STOREV = 'R', LDV >= K.
*> \endverbatim
*>
*> \param[in] T
*> \verbatim
*>          T is DOUBLE PRECISION array, dimension (LDT,K)
*>          The triangular k by k matrix T in the representation of the
*>          block reflector.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T. LDT >= K.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension (LDC,N)
*>          On entry, the m by n matrix C.
*>          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (LDWORK,K)
*> \endverbatim
*>
*> \param[in] LDWORK
*> \verbatim
*>          LDWORK is INTEGER
*>          The leading dimension of the array WORK.
*>          If SIDE = 'L', LDWORK >= max(1,N);
*>          if SIDE = 'R', LDWORK >= max(1,M).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date June 2013
*
*> \ingroup doubleOTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The shape of the matrix V and the storage of the vectors which define
*>  the H(i) is best illustrated by the following example with n = 5 and
*>  k = 3. The elements equal to 1 are not stored; the corresponding
*>  array elements are modified but restored on exit. The rest of the
*>  array is not used.
*>
*>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*>
*>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*>                   ( v1  1    )                     (     1 v2 v2 v2 )
*>                   ( v1 v2  1 )                     (        1 v3 v3 )
*>                   ( v1 v2 v3 )
*>                   ( v1 v2 v3 )
*>
*>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*>
*>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
*>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
*>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*>                   (     1 v3 )
*>                   (        1 )
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
*
*  -- LAPACK auxiliary routine (version 3.5.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2013
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DTRMM
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
*
      IF( LSAME( STOREV, 'C' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1 )    (first K rows)
*                     ( V2 )
*           where  V1  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H**T * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
*
*              W := C1**T
*
               DO 10 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10          CONTINUE
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2**T * V2
*
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,
     $                        ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T**T  or  W * T
*
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W**T
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2 * W**T
*
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,
     $                        -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1**T
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W**T
*
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
   30          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C1
*
               DO 40 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T**T
*
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V**T
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2**T
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1**T
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
*
         ELSE
*
*           Let  V =  ( V1 )
*                     ( V2 )    (last K rows)
*           where  V2  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H**T * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
*
*              W := C2**T
*
               DO 70 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
   70          CONTINUE
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1**T * V1
*
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T**T  or  W * T
*
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W**T
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1 * W**T
*
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,
     $                        -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
*
*              W := W * V2**T
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $                     ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W**T
*
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
   80             CONTINUE
   90          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C2
*
               DO 100 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T**T
*
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V**T
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1**T
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
*
*              W := W * V2**T
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
     $                     ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W
*
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
*
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1  V2 )    (V1: first K columns)
*           where  V1  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H**T * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
*
*              W := C1**T
*
               DO 130 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
*
*              W := W * V1**T
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2**T * V2**T
*
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                        C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE,
     $                        WORK, LDWORK )
               END IF
*
*              W := W * T**T  or  W * T
*
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V**T * W**T
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2**T * W**T
*
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
     $                        V( 1, K+1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W**T
*
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
*
*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
*
*              W := C1
*
               DO 160 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
*
*              W := W * V1**T
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2**T
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T**T
*
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
*
            END IF
*
         ELSE
*
*           Let  V =  ( V1  V2 )    (V2: last K columns)
*           where  V2  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H**T * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
*
*              W := C2**T
*
               DO 190 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
  190          CONTINUE
*
*              W := W * V2**T
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
     $                     ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1**T * V1**T
*
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                        C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T**T  or  W * T
*
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V**T * W**T
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1**T * W**T
*
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
     $                        V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W**T
*
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
  200             CONTINUE
  210          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
*
*              W := C2
*
               DO 220 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
*
*              W := W * V2**T
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
     $                     ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1**T
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T**T
*
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
*
            END IF
*
         END IF
      END IF
*
      RETURN
*
*     End of DLARFB
*
      END
*> \brief \b DLARF applies an elementary reflector to a general rectangular matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLARF + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
* 
*       .. Scalar Arguments ..
*       CHARACTER          SIDE
*       INTEGER            INCV, LDC, M, N
*       DOUBLE PRECISION   TAU
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLARF applies a real elementary reflector H to a real m by n matrix
*> C, from either the left or the right. H is represented in the form
*>
*>       H = I - tau * v * v**T
*>
*> where tau is a real scalar and v is a real vector.
*>
*> If tau = 0, then H is taken to be the unit matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': form  H * C
*>          = 'R': form  C * H
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is DOUBLE PRECISION array, dimension
*>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*>          The vector v in the representation of H. V is not used if
*>          TAU = 0.
*> \endverbatim
*>
*> \param[in] INCV
*> \verbatim
*>          INCV is INTEGER
*>          The increment between elements of v. INCV <> 0.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION
*>          The value tau in the representation of H.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension (LDC,N)
*>          On entry, the m by n matrix C.
*>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*>          or C * H if SIDE = 'R'.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension
*>                         (N) if SIDE = 'L'
*>                      or (M) if SIDE = 'R'
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILADLR, ILADLC
      EXTERNAL           LSAME, ILADLR, ILADLC
*     ..
*     .. Executable Statements ..
*
      APPLYLEFT = LSAME( SIDE, 'L' )
      LASTV = 0
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
         IF( APPLYLEFT ) THEN
            LASTV = M
         ELSE
            LASTV = N
         END IF
         IF( INCV.GT.0 ) THEN
            I = 1 + (LASTV-1) * INCV
         ELSE
            I = 1
         END IF
!     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
            LASTV = LASTV - 1
            I = I - INCV
         END DO
         IF( APPLYLEFT ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILADLC(LASTV, N, C, LDC)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILADLR(M, LASTV, C, LDC)
         END IF
      END IF
!     Note that lastc.eq.0 renders the BLAS operations null; no special
!     case is needed at this level.
      IF( APPLYLEFT ) THEN
*
*        Form  H * C
*
         IF( LASTV.GT.0 ) THEN
*
*           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
*
            CALL DGEMV( 'Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV,
     $           ZERO, WORK, 1 )
*
*           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
*
            CALL DGER( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
*
*        Form  C * H
*
         IF( LASTV.GT.0 ) THEN
*
*           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
*
            CALL DGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC,
     $           V, INCV, ZERO, WORK, 1 )
*
*           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
*
            CALL DGER( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
*
*     End of DLARF
*
      END
*> \brief \b DLARFG generates an elementary reflector (Householder matrix).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLARFG + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfg.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfg.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfg.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
* 
*       .. Scalar Arguments ..
*       INTEGER            INCX, N
*       DOUBLE PRECISION   ALPHA, TAU
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   X( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLARFG generates a real elementary reflector H of order n, such
*> that
*>
*>       H * ( alpha ) = ( beta ),   H**T * H = I.
*>           (   x   )   (   0  )
*>
*> where alpha and beta are scalars, and x is an (n-1)-element real
*> vector. H is represented in the form
*>
*>       H = I - tau * ( 1 ) * ( 1 v**T ) ,
*>                     ( v )
*>
*> where tau is a real scalar and v is a real (n-1)-element
*> vector.
*>
*> If the elements of x are all zero, then tau = 0 and H is taken to be
*> the unit matrix.
*>
*> Otherwise  1 <= tau <= 2.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the elementary reflector.
*> \endverbatim
*>
*> \param[in,out] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION
*>          On entry, the value alpha.
*>          On exit, it is overwritten with the value beta.
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is DOUBLE PRECISION array, dimension
*>                         (1+(N-2)*abs(INCX))
*>          On entry, the vector x.
*>          On exit, it is overwritten with the vector v.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between elements of X. INCX > 0.
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION
*>          The value tau.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
      EXTERNAL           DLAMCH, DLAPY2, DNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = DNRM2( N-1, X, INCX )
*
      IF( XNORM.EQ.ZERO ) THEN
*
*        H  =  I
*
         TAU = ZERO
      ELSE
*
*        general case
*
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
            RSAFMN = ONE / SAFMIN
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            XNORM = DNRM2( N-1, X, INCX )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         END IF
         TAU = ( BETA-ALPHA ) / BETA
         CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
*
*        If ALPHA is subnormal, it may lose relative accuracy
*
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
*
      RETURN
*
*     End of DLARFG
*
      END
*> \brief \b DLARFT forms the triangular factor T of a block reflector H = I - vtvH
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLARFT + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarft.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarft.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarft.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
* 
*       .. Scalar Arguments ..
*       CHARACTER          DIRECT, STOREV
*       INTEGER            K, LDT, LDV, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLARFT forms the triangular factor T of a real block reflector H
*> of order n, which is defined as a product of k elementary reflectors.
*>
*> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
*>
*> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
*>
*> If STOREV = 'C', the vector which defines the elementary reflector
*> H(i) is stored in the i-th column of the array V, and
*>
*>    H  =  I - V * T * V**T
*>
*> If STOREV = 'R', the vector which defines the elementary reflector
*> H(i) is stored in the i-th row of the array V, and
*>
*>    H  =  I - V**T * T * V
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] DIRECT
*> \verbatim
*>          DIRECT is CHARACTER*1
*>          Specifies the order in which the elementary reflectors are
*>          multiplied to form the block reflector:
*>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*> \endverbatim
*>
*> \param[in] STOREV
*> \verbatim
*>          STOREV is CHARACTER*1
*>          Specifies how the vectors which define the elementary
*>          reflectors are stored (see also Further Details):
*>          = 'C': columnwise
*>          = 'R': rowwise
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the block reflector H. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The order of the triangular factor T (= the number of
*>          elementary reflectors). K >= 1.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is DOUBLE PRECISION array, dimension
*>                               (LDV,K) if STOREV = 'C'
*>                               (LDV,N) if STOREV = 'R'
*>          The matrix V. See further details.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i).
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is DOUBLE PRECISION array, dimension (LDT,K)
*>          The k by k triangular factor T of the block reflector.
*>          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
*>          lower triangular. The rest of the array is not used.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T. LDT >= K.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleOTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The shape of the matrix V and the storage of the vectors which define
*>  the H(i) is best illustrated by the following example with n = 5 and
*>  k = 3. The elements equal to 1 are not stored.
*>
*>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*>
*>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*>                   ( v1  1    )                     (     1 v2 v2 v2 )
*>                   ( v1 v2  1 )                     (        1 v3 v3 )
*>                   ( v1 v2 v3 )
*>                   ( v1 v2 v3 )
*>
*>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*>
*>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
*>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
*>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*>                   (     1 v3 )
*>                   (        1 )
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, PREVLASTV, LASTV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DTRMV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( LSAME( DIRECT, 'F' ) ) THEN
         PREVLASTV = N
         DO I = 1, K
            PREVLASTV = MAX( I, PREVLASTV )
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO J = 1, I
                  T( J, I ) = ZERO
               END DO
            ELSE
*
*              general case
*
               IF( LSAME( STOREV, 'C' ) ) THEN
*                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( LASTV, I ).NE.ZERO ) EXIT
                  END DO
                  DO J = 1, I-1
                     T( J, I ) = -TAU( I ) * V( I , J )
                  END DO   
                  J = MIN( LASTV, PREVLASTV )
*
*                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**T * V(i:j,i)
*
                  CALL DGEMV( 'Transpose', J-I, I-1, -TAU( I ), 
     $                        V( I+1, 1 ), LDV, V( I+1, I ), 1, ONE, 
     $                        T( 1, I ), 1 )
               ELSE
*                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( I, LASTV ).NE.ZERO ) EXIT
                  END DO
                  DO J = 1, I-1
                     T( J, I ) = -TAU( I ) * V( J , I )
                  END DO   
                  J = MIN( LASTV, PREVLASTV )
*
*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**T
*
                  CALL DGEMV( 'No transpose', I-1, J-I, -TAU( I ),
     $                        V( 1, I+1 ), LDV, V( I, I+1 ), LDV, ONE,
     $                        T( 1, I ), 1 )
               END IF
*
*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
*
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T,
     $                     LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
               IF( I.GT.1 ) THEN
                  PREVLASTV = MAX( PREVLASTV, LASTV )
               ELSE
                  PREVLASTV = LASTV
               END IF
            END IF
         END DO
      ELSE
         PREVLASTV = 1
         DO I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO J = I, K
                  T( J, I ) = ZERO
               END DO
            ELSE
*
*              general case
*
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
*                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( LASTV, I ).NE.ZERO ) EXIT
                     END DO
                     DO J = I+1, K
                        T( J, I ) = -TAU( I ) * V( N-K+I , J )
                     END DO   
                     J = MAX( LASTV, PREVLASTV )
*
*                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**T * V(j:n-k+i,i)
*
                     CALL DGEMV( 'Transpose', N-K+I-J, K-I, -TAU( I ),
     $                           V( J, I+1 ), LDV, V( J, I ), 1, ONE,
     $                           T( I+1, I ), 1 )
                  ELSE
*                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( I, LASTV ).NE.ZERO ) EXIT
                     END DO
                     DO J = I+1, K
                        T( J, I ) = -TAU( I ) * V( J, N-K+I )
                     END DO   
                     J = MAX( LASTV, PREVLASTV )
*
*                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**T
*
                     CALL DGEMV( 'No transpose', K-I, N-K+I-J,
     $                    -TAU( I ), V( I+1, J ), LDV, V( I, J ), LDV,
     $                    ONE, T( I+1, I ), 1 )
                  END IF
*
*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
*
                  CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I,
     $                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
                  IF( I.GT.1 ) THEN
                     PREVLASTV = MIN( PREVLASTV, LASTV )
                  ELSE
                     PREVLASTV = LASTV
                  END IF
               END IF
               T( I, I ) = TAU( I )
            END IF
         END DO
      END IF
      RETURN
*
*     End of DLARFT
*
      END
*> \brief \b DLARTG generates a plane rotation with real cosine and real sine.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLARTG + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlartg.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlartg.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlartg.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLARTG( F, G, CS, SN, R )
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   CS, F, G, R, SN
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLARTG generate a plane rotation so that
*>
*>    [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
*>    [ -SN  CS  ]     [ G ]     [ 0 ]
*>
*> This is a slower, more accurate version of the BLAS1 routine DROTG,
*> with the following other differences:
*>    F and G are unchanged on return.
*>    If G=0, then CS=1 and SN=0.
*>    If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
*>       floating point operations (saves work in DBDSQR when
*>       there are zeros on the diagonal).
*>
*> If F exceeds G in magnitude, CS will be positive.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] F
*> \verbatim
*>          F is DOUBLE PRECISION
*>          The first component of vector to be rotated.
*> \endverbatim
*>
*> \param[in] G
*> \verbatim
*>          G is DOUBLE PRECISION
*>          The second component of vector to be rotated.
*> \endverbatim
*>
*> \param[out] CS
*> \verbatim
*>          CS is DOUBLE PRECISION
*>          The cosine of the rotation.
*> \endverbatim
*>
*> \param[out] SN
*> \verbatim
*>          SN is DOUBLE PRECISION
*>          The sine of the rotation.
*> \endverbatim
*>
*> \param[out] R
*> \verbatim
*>          R is DOUBLE PRECISION
*>          The nonzero component of the rotated vector.
*>
*>  This version has a few statements commented out for thread safety
*>  (machine parameters are computed on each entry). 10 feb 03, SJH.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE DLARTG( F, G, CS, SN, R )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CS, F, G, R, SN
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
*     LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
*     ..
*     .. Save statement ..
*     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
*     ..
*     .. Data statements ..
*     DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
*     IF( FIRST ) THEN
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $            LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
*        FIRST = .FALSE.
*     END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 )
     $         GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 )
     $         GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
            CS = -CS
            SN = -SN
            R = -R
         END IF
      END IF
      RETURN
*
*     End of DLARTG
*
      END
*> \brief \b DLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLASCL + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlascl.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlascl.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlascl.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          TYPE
*       INTEGER            INFO, KL, KU, LDA, M, N
*       DOUBLE PRECISION   CFROM, CTO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLASCL multiplies the M by N real matrix A by the real scalar
*> CTO/CFROM.  This is done without over/underflow as long as the final
*> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
*> A may be full, upper triangular, lower triangular, upper Hessenberg,
*> or banded.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TYPE
*> \verbatim
*>          TYPE is CHARACTER*1
*>          TYPE indices the storage type of the input matrix.
*>          = 'G':  A is a full matrix.
*>          = 'L':  A is a lower triangular matrix.
*>          = 'U':  A is an upper triangular matrix.
*>          = 'H':  A is an upper Hessenberg matrix.
*>          = 'B':  A is a symmetric band matrix with lower bandwidth KL
*>                  and upper bandwidth KU and with the only the lower
*>                  half stored.
*>          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
*>                  and upper bandwidth KU and with the only the upper
*>                  half stored.
*>          = 'Z':  A is a band matrix with lower bandwidth KL and upper
*>                  bandwidth KU. See DGBTRF for storage details.
*> \endverbatim
*>
*> \param[in] KL
*> \verbatim
*>          KL is INTEGER
*>          The lower bandwidth of A.  Referenced only if TYPE = 'B',
*>          'Q' or 'Z'.
*> \endverbatim
*>
*> \param[in] KU
*> \verbatim
*>          KU is INTEGER
*>          The upper bandwidth of A.  Referenced only if TYPE = 'B',
*>          'Q' or 'Z'.
*> \endverbatim
*>
*> \param[in] CFROM
*> \verbatim
*>          CFROM is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] CTO
*> \verbatim
*>          CTO is DOUBLE PRECISION
*>
*>          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
*>          without over/underflow if the final result CTO*A(I,J)/CFROM
*>          can be represented without over/underflow.  CFROM must be
*>          nonzero.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
*>          storage type.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          0  - successful exit
*>          <0 - if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH, DISNAN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
*
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
*
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO .OR. DISNAN(CFROM) ) THEN
         INFO = -4
      ELSE IF( DISNAN(CTO) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.
     $         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.
     $            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )
     $             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.
     $            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.
     $            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASCL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
*
      CFROMC = CFROM
      CTOC = CTO
*
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      IF( CFROM1.EQ.CFROMC ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
         MUL = CTOC / CFROMC
         DONE = .TRUE.
         CTO1 = CTOC
      ELSE
         CTO1 = CTOC / BIGNUM
         IF( CTO1.EQ.CTOC ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
            MUL = CTOC
            DONE = .TRUE.
            CFROMC = ONE
         ELSE IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
            MUL = SMLNUM
            DONE = .FALSE.
            CFROMC = CFROM1
         ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
            MUL = BIGNUM
            DONE = .FALSE.
            CTOC = CTO1
         ELSE
            MUL = CTOC / CFROMC
            DONE = .TRUE.
         END IF
      END IF
*
      IF( ITYPE.EQ.0 ) THEN
*
*        Full matrix
*
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
*
      ELSE IF( ITYPE.EQ.1 ) THEN
*
*        Lower triangular matrix
*
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        Upper triangular matrix
*
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        Upper Hessenberg matrix
*
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
*
      ELSE IF( ITYPE.EQ.4 ) THEN
*
*        Lower half of a symmetric band matrix
*
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
*
      ELSE IF( ITYPE.EQ.5 ) THEN
*
*        Upper half of a symmetric band matrix
*
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
*
      ELSE IF( ITYPE.EQ.6 ) THEN
*
*        Band matrix
*
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
*
      END IF
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of DLASCL
*
      END
*> \brief \b DLASET initializes the off-diagonal elements and the diagonal elements of a matrix to given values.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLASET + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaset.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaset.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaset.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
* 
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            LDA, M, N
*       DOUBLE PRECISION   ALPHA, BETA
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLASET initializes an m-by-n matrix A to BETA on the diagonal and
*> ALPHA on the offdiagonals.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies the part of the matrix A to be set.
*>          = 'U':      Upper triangular part is set; the strictly lower
*>                      triangular part of A is not changed.
*>          = 'L':      Lower triangular part is set; the strictly upper
*>                      triangular part of A is not changed.
*>          Otherwise:  All of the matrix A is set.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION
*>          The constant to which the offdiagonal elements are to be set.
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION
*>          The constant to which the diagonal elements are to be set.
*> \endverbatim
*>
*> \param[out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On exit, the leading m-by-n submatrix of A is set as follows:
*>
*>          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
*>          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
*>          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
*>
*>          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Set the strictly upper triangular or trapezoidal part of the
*        array to ALPHA.
*
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
*
*        Set the strictly lower triangular or trapezoidal part of the
*        array to ALPHA.
*
         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
*
      ELSE
*
*        Set the leading m-by-n submatrix to ALPHA.
*
         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      END IF
*
*     Set the first min(M,N) diagonal elements to BETA.
*
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE
*
      RETURN
*
*     End of DLASET
*
      END
*> \brief \b DLASR applies a sequence of plane rotations to a general rectangular matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLASR + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasr.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasr.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasr.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
* 
*       .. Scalar Arguments ..
*       CHARACTER          DIRECT, PIVOT, SIDE
*       INTEGER            LDA, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLASR applies a sequence of plane rotations to a real matrix A,
*> from either the left or the right.
*> 
*> When SIDE = 'L', the transformation takes the form
*> 
*>    A := P*A
*> 
*> and when SIDE = 'R', the transformation takes the form
*> 
*>    A := A*P**T
*> 
*> where P is an orthogonal matrix consisting of a sequence of z plane
*> rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
*> and P**T is the transpose of P.
*> 
*> When DIRECT = 'F' (Forward sequence), then
*> 
*>    P = P(z-1) * ... * P(2) * P(1)
*> 
*> and when DIRECT = 'B' (Backward sequence), then
*> 
*>    P = P(1) * P(2) * ... * P(z-1)
*> 
*> where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
*> 
*>    R(k) = (  c(k)  s(k) )
*>         = ( -s(k)  c(k) ).
*> 
*> When PIVOT = 'V' (Variable pivot), the rotation is performed
*> for the plane (k,k+1), i.e., P(k) has the form
*> 
*>    P(k) = (  1                                            )
*>           (       ...                                     )
*>           (              1                                )
*>           (                   c(k)  s(k)                  )
*>           (                  -s(k)  c(k)                  )
*>           (                                1              )
*>           (                                     ...       )
*>           (                                            1  )
*> 
*> where R(k) appears as a rank-2 modification to the identity matrix in
*> rows and columns k and k+1.
*> 
*> When PIVOT = 'T' (Top pivot), the rotation is performed for the
*> plane (1,k+1), so P(k) has the form
*> 
*>    P(k) = (  c(k)                    s(k)                 )
*>           (         1                                     )
*>           (              ...                              )
*>           (                     1                         )
*>           ( -s(k)                    c(k)                 )
*>           (                                 1             )
*>           (                                      ...      )
*>           (                                             1 )
*> 
*> where R(k) appears in rows and columns 1 and k+1.
*> 
*> Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
*> performed for the plane (k,z), giving P(k) the form
*> 
*>    P(k) = ( 1                                             )
*>           (      ...                                      )
*>           (             1                                 )
*>           (                  c(k)                    s(k) )
*>           (                         1                     )
*>           (                              ...              )
*>           (                                     1         )
*>           (                 -s(k)                    c(k) )
*> 
*> where R(k) appears in rows and columns k and z.  The rotations are
*> performed without ever forming P(k) explicitly.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          Specifies whether the plane rotation matrix P is applied to
*>          A on the left or the right.
*>          = 'L':  Left, compute A := P*A
*>          = 'R':  Right, compute A:= A*P**T
*> \endverbatim
*>
*> \param[in] PIVOT
*> \verbatim
*>          PIVOT is CHARACTER*1
*>          Specifies the plane for which P(k) is a plane rotation
*>          matrix.
*>          = 'V':  Variable pivot, the plane (k,k+1)
*>          = 'T':  Top pivot, the plane (1,k+1)
*>          = 'B':  Bottom pivot, the plane (k,z)
*> \endverbatim
*>
*> \param[in] DIRECT
*> \verbatim
*>          DIRECT is CHARACTER*1
*>          Specifies whether P is a forward or backward sequence of
*>          plane rotations.
*>          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1)
*>          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1)
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  If m <= 1, an immediate
*>          return is effected.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  If n <= 1, an
*>          immediate return is effected.
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension
*>                  (M-1) if SIDE = 'L'
*>                  (N-1) if SIDE = 'R'
*>          The cosines c(k) of the plane rotations.
*> \endverbatim
*>
*> \param[in] S
*> \verbatim
*>          S is DOUBLE PRECISION array, dimension
*>                  (M-1) if SIDE = 'L'
*>                  (N-1) if SIDE = 'R'
*>          The sines s(k) of the plane rotations.  The 2-by-2 plane
*>          rotation part of the matrix P(k), R(k), has the form
*>          R(k) = (  c(k)  s(k) )
*>                 ( -s(k)  c(k) ).
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          The M-by-N matrix A.  On exit, A is overwritten by P*A if
*>          SIDE = 'R' or by A*P**T if SIDE = 'L'.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASR ', INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  P * A
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*        Form A * P**T
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DLASR
*
      END
*> \brief \b DLASRT sorts numbers in increasing or decreasing order.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLASRT + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasrt.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasrt.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasrt.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLASRT( ID, N, D, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          ID
*       INTEGER            INFO, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   D( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> Sort the numbers in D in increasing order (if ID = 'I') or
*> in decreasing order (if ID = 'D' ).
*>
*> Use Quick Sort, reverting to Insertion sort on arrays of
*> size <= 20. Dimension of STACK limits N to about 2**32.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ID
*> \verbatim
*>          ID is CHARACTER*1
*>          = 'I': sort D in increasing order;
*>          = 'D': sort D in decreasing order.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The length of the array D.
*> \endverbatim
*>
*> \param[in,out] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (N)
*>          On entry, the array to be sorted.
*>          On exit, D has been sorted into increasing order
*>          (D(1) <= ... <= D(N) ) or into decreasing order
*>          (D(1) >= ... >= D(N) ), depending on ID.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERcomputational
*
*  =====================================================================
      SUBROUTINE DLASRT( ID, N, D, INFO )
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
*     ..
*     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT
      DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
*     ..
*     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input paramters.
*
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
*
*        Do Insertion sort on D( START:ENDD )
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( D( J ).GT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE
*
         ELSE
*
*           Sort into increasing order
*
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( D( J ).LT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE
*
         END IF
*
      ELSE IF( ENDD-START.GT.SELECT ) THEN
*
*        Partition D( START:ENDD ) and stack parts, largest one first
*
*        Choose partition entry as median of 3
*
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF( D1.LT.D2 ) THEN
            IF( D3.LT.D1 ) THEN
               DMNMX = D1
            ELSE IF( D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3.LT.D2 ) THEN
               DMNMX = D2
            ELSE IF( D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( D( J ).LT.DMNMX )
     $         GO TO 70
   80       CONTINUE
            I = I + 1
            IF( D( I ).GT.DMNMX )
     $         GO TO 80
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 60
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
*
*           Sort into increasing order
*
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( D( J ).GT.DMNMX )
     $         GO TO 100
  110       CONTINUE
            I = I + 1
            IF( D( I ).LT.DMNMX )
     $         GO TO 110
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 90
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 )
     $   GO TO 10
      RETURN
*
*     End of DLASRT
*
      END
*> \brief \b DLASSQ updates a sum of squares represented in scaled form.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLASSQ + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlassq.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlassq.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlassq.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
* 
*       .. Scalar Arguments ..
*       INTEGER            INCX, N
*       DOUBLE PRECISION   SCALE, SUMSQ
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   X( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLASSQ  returns the values  scl  and  smsq  such that
*>
*>    ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*>
*> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
*> assumed to be non-negative and  scl  returns the value
*>
*>    scl = max( scale, abs( x( i ) ) ).
*>
*> scale and sumsq must be supplied in SCALE and SUMSQ and
*> scl and smsq are overwritten on SCALE and SUMSQ respectively.
*>
*> The routine makes only one pass through the vector x.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of elements to be used from the vector X.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is DOUBLE PRECISION array, dimension (N)
*>          The vector for which a scaled sum of squares is computed.
*>             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between successive values of the vector X.
*>          INCX > 0.
*> \endverbatim
*>
*> \param[in,out] SCALE
*> \verbatim
*>          SCALE is DOUBLE PRECISION
*>          On entry, the value  scale  in the equation above.
*>          On exit, SCALE is overwritten with  scl , the scaling factor
*>          for the sum of squares.
*> \endverbatim
*>
*> \param[in,out] SUMSQ
*> \verbatim
*>          SUMSQ is DOUBLE PRECISION
*>          On entry, the value  sumsq  in the equation above.
*>          On exit, SUMSQ is overwritten with  smsq , the basic sum of
*>          squares from which  scl  has been factored out.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
*     ..
*     .. External Functions ..
      LOGICAL            DISNAN
      EXTERNAL           DISNAN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            ABSXI = ABS( X( IX ) )
            IF( ABSXI.GT.ZERO.OR.DISNAN( ABSXI ) ) THEN
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     End of DLASSQ
*
      END
*> \brief \b DLASV2 computes the singular value decomposition of a 2-by-2 triangular matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLASV2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasv2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasv2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasv2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL )
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLASV2 computes the singular value decomposition of a 2-by-2
*> triangular matrix
*>    [  F   G  ]
*>    [  0   H  ].
*> On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the
*> smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and
*> right singular vectors for abs(SSMAX), giving the decomposition
*>
*>    [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]
*>    [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] F
*> \verbatim
*>          F is DOUBLE PRECISION
*>          The (1,1) element of the 2-by-2 matrix.
*> \endverbatim
*>
*> \param[in] G
*> \verbatim
*>          G is DOUBLE PRECISION
*>          The (1,2) element of the 2-by-2 matrix.
*> \endverbatim
*>
*> \param[in] H
*> \verbatim
*>          H is DOUBLE PRECISION
*>          The (2,2) element of the 2-by-2 matrix.
*> \endverbatim
*>
*> \param[out] SSMIN
*> \verbatim
*>          SSMIN is DOUBLE PRECISION
*>          abs(SSMIN) is the smaller singular value.
*> \endverbatim
*>
*> \param[out] SSMAX
*> \verbatim
*>          SSMAX is DOUBLE PRECISION
*>          abs(SSMAX) is the larger singular value.
*> \endverbatim
*>
*> \param[out] SNL
*> \verbatim
*>          SNL is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[out] CSL
*> \verbatim
*>          CSL is DOUBLE PRECISION
*>          The vector (CSL, SNL) is a unit left singular vector for the
*>          singular value abs(SSMAX).
*> \endverbatim
*>
*> \param[out] SNR
*> \verbatim
*>          SNR is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[out] CSR
*> \verbatim
*>          CSR is DOUBLE PRECISION
*>          The vector (CSR, SNR) is a unit right singular vector for the
*>          singular value abs(SSMAX).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Any input parameter may be aliased with any output parameter.
*>
*>  Barring over/underflow and assuming a guard digit in subtraction, all
*>  output quantities are correct to within a few units in the last
*>  place (ulps).
*>
*>  In IEEE arithmetic, the code works correctly if one matrix element is
*>  infinite.
*>
*>  Overflow will not occur unless the largest singular value itself
*>  overflows or is within a few ulps of overflow. (On machines with
*>  partial overflow, like the Cray, overflow may occur if the largest
*>  singular value is within a factor of 2 of overflow.)
*>
*>  Underflow is harmless if underflow is gradual. Otherwise, results
*>  may correspond to a matrix modified by perturbations of size near
*>  the underflow threshold.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   FOUR
      PARAMETER          ( FOUR = 4.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            GASMAL, SWAP
      INTEGER            PMAX
      DOUBLE PRECISION   A, CLT, CRT, D, FA, FT, GA, GT, HA, HT, L, M,
     $                   MM, R, S, SLT, SRT, T, TEMP, TSIGN, TT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Executable Statements ..
*
      FT = F
      FA = ABS( FT )
      HT = H
      HA = ABS( H )
*
*     PMAX points to the maximum absolute element of matrix
*       PMAX = 1 if F largest in absolute values
*       PMAX = 2 if G largest in absolute values
*       PMAX = 3 if H largest in absolute values
*
      PMAX = 1
      SWAP = ( HA.GT.FA )
      IF( SWAP ) THEN
         PMAX = 3
         TEMP = FT
         FT = HT
         HT = TEMP
         TEMP = FA
         FA = HA
         HA = TEMP
*
*        Now FA .ge. HA
*
      END IF
      GT = G
      GA = ABS( GT )
      IF( GA.EQ.ZERO ) THEN
*
*        Diagonal matrix
*
         SSMIN = HA
         SSMAX = FA
         CLT = ONE
         CRT = ONE
         SLT = ZERO
         SRT = ZERO
      ELSE
         GASMAL = .TRUE.
         IF( GA.GT.FA ) THEN
            PMAX = 2
            IF( ( FA / GA ).LT.DLAMCH( 'EPS' ) ) THEN
*
*              Case of very large GA
*
               GASMAL = .FALSE.
               SSMAX = GA
               IF( HA.GT.ONE ) THEN
                  SSMIN = FA / ( GA / HA )
               ELSE
                  SSMIN = ( FA / GA )*HA
               END IF
               CLT = ONE
               SLT = HT / GT
               SRT = ONE
               CRT = FT / GT
            END IF
         END IF
         IF( GASMAL ) THEN
*
*           Normal case
*
            D = FA - HA
            IF( D.EQ.FA ) THEN
*
*              Copes with infinite F or H
*
               L = ONE
            ELSE
               L = D / FA
            END IF
*
*           Note that 0 .le. L .le. 1
*
            M = GT / FT
*
*           Note that abs(M) .le. 1/macheps
*
            T = TWO - L
*
*           Note that T .ge. 1
*
            MM = M*M
            TT = T*T
            S = SQRT( TT+MM )
*
*           Note that 1 .le. S .le. 1 + 1/macheps
*
            IF( L.EQ.ZERO ) THEN
               R = ABS( M )
            ELSE
               R = SQRT( L*L+MM )
            END IF
*
*           Note that 0 .le. R .le. 1 + 1/macheps
*
            A = HALF*( S+R )
*
*           Note that 1 .le. A .le. 1 + abs(M)
*
            SSMIN = HA / A
            SSMAX = FA*A
            IF( MM.EQ.ZERO ) THEN
*
*              Note that M is very tiny
*
               IF( L.EQ.ZERO ) THEN
                  T = SIGN( TWO, FT )*SIGN( ONE, GT )
               ELSE
                  T = GT / SIGN( D, FT ) + M / T
               END IF
            ELSE
               T = ( M / ( S+T )+M / ( R+L ) )*( ONE+A )
            END IF
            L = SQRT( T*T+FOUR )
            CRT = TWO / L
            SRT = T / L
            CLT = ( CRT+SRT*M ) / A
            SLT = ( HT / FT )*SRT / A
         END IF
      END IF
      IF( SWAP ) THEN
         CSL = SRT
         SNL = CRT
         CSR = SLT
         SNR = CLT
      ELSE
         CSL = CLT
         SNL = SLT
         CSR = CRT
         SNR = SRT
      END IF
*
*     Correct signs of SSMAX and SSMIN
*
      IF( PMAX.EQ.1 )
     $   TSIGN = SIGN( ONE, CSR )*SIGN( ONE, CSL )*SIGN( ONE, F )
      IF( PMAX.EQ.2 )
     $   TSIGN = SIGN( ONE, SNR )*SIGN( ONE, CSL )*SIGN( ONE, G )
      IF( PMAX.EQ.3 )
     $   TSIGN = SIGN( ONE, SNR )*SIGN( ONE, SNL )*SIGN( ONE, H )
      SSMAX = SIGN( SSMAX, TSIGN )
      SSMIN = SIGN( SSMIN, TSIGN*SIGN( ONE, F )*SIGN( ONE, H ) )
      RETURN
*
*     End of DLASV2
*
      END
*> \brief \b DLASWP performs a series of row interchanges on a general rectangular matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLASWP + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaswp.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaswp.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaswp.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
* 
*       .. Scalar Arguments ..
*       INTEGER            INCX, K1, K2, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLASWP performs a series of row interchanges on the matrix A.
*> One row interchange is initiated for each of rows K1 through K2 of A.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the matrix of column dimension N to which the row
*>          interchanges will be applied.
*>          On exit, the permuted matrix.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*> \endverbatim
*>
*> \param[in] K1
*> \verbatim
*>          K1 is INTEGER
*>          The first element of IPIV for which a row interchange will
*>          be done.
*> \endverbatim
*>
*> \param[in] K2
*> \verbatim
*>          K2 is INTEGER
*>          The last element of IPIV for which a row interchange will
*>          be done.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (K2*abs(INCX))
*>          The vector of pivot indices.  Only the elements in positions
*>          K1 through K2 of IPIV are accessed.
*>          IPIV(K) = L implies rows K and L are to be interchanged.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between successive values of IPIV.  If IPIV
*>          is negative, the pivots are applied in reverse order.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleOTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Modified by
*>   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
*
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASWP
*
      END
*> \brief \b DLATRD reduces the first nb rows and columns of a symmetric/Hermitian matrix A to real tridiagonal form by an orthogonal similarity transformation.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DLATRD + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlatrd.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlatrd.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlatrd.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
* 
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            LDA, LDW, N, NB
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), E( * ), TAU( * ), W( LDW, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLATRD reduces NB rows and columns of a real symmetric matrix A to
*> symmetric tridiagonal form by an orthogonal similarity
*> transformation Q**T * A * Q, and returns the matrices V and W which are
*> needed to apply the transformation to the unreduced part of A.
*>
*> If UPLO = 'U', DLATRD reduces the last NB rows and columns of a
*> matrix, of which the upper triangle is supplied;
*> if UPLO = 'L', DLATRD reduces the first NB rows and columns of a
*> matrix, of which the lower triangle is supplied.
*>
*> This is an auxiliary routine called by DSYTRD.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          symmetric matrix A is stored:
*>          = 'U': Upper triangular
*>          = 'L': Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The number of rows and columns to be reduced.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*>          n-by-n upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading n-by-n lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the strictly upper
*>          triangular part of A is not referenced.
*>          On exit:
*>          if UPLO = 'U', the last NB columns have been reduced to
*>            tridiagonal form, with the diagonal elements overwriting
*>            the diagonal elements of A; the elements above the diagonal
*>            with the array TAU, represent the orthogonal matrix Q as a
*>            product of elementary reflectors;
*>          if UPLO = 'L', the first NB columns have been reduced to
*>            tridiagonal form, with the diagonal elements overwriting
*>            the diagonal elements of A; the elements below the diagonal
*>            with the array TAU, represent the  orthogonal matrix Q as a
*>            product of elementary reflectors.
*>          See Further Details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= (1,N).
*> \endverbatim
*>
*> \param[out] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (N-1)
*>          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
*>          elements of the last NB columns of the reduced matrix;
*>          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
*>          the first NB columns of the reduced matrix.
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (N-1)
*>          The scalar factors of the elementary reflectors, stored in
*>          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
*>          See Further Details.
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is DOUBLE PRECISION array, dimension (LDW,NB)
*>          The n-by-nb matrix W required to update the unreduced part
*>          of A.
*> \endverbatim
*>
*> \param[in] LDW
*> \verbatim
*>          LDW is INTEGER
*>          The leading dimension of the array W. LDW >= max(1,N).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleOTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  If UPLO = 'U', the matrix Q is represented as a product of elementary
*>  reflectors
*>
*>     Q = H(n) H(n-1) . . . H(n-nb+1).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**T
*>
*>  where tau is a real scalar, and v is a real vector with
*>  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
*>  and tau in TAU(i-1).
*>
*>  If UPLO = 'L', the matrix Q is represented as a product of elementary
*>  reflectors
*>
*>     Q = H(1) H(2) . . . H(nb).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**T
*>
*>  where tau is a real scalar, and v is a real vector with
*>  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
*>  and tau in TAU(i).
*>
*>  The elements of the vectors v together form the n-by-nb matrix V
*>  which is needed, with W, to apply the transformation to the unreduced
*>  part of the matrix, using a symmetric rank-2k update of the form:
*>  A := A - V*W**T - W*V**T.
*>
*>  The contents of A on exit are illustrated by the following examples
*>  with n = 5 and nb = 2:
*>
*>  if UPLO = 'U':                       if UPLO = 'L':
*>
*>    (  a   a   a   v4  v5 )              (  d                  )
*>    (      a   a   v4  v5 )              (  1   d              )
*>    (          a   1   v5 )              (  v1  1   a          )
*>    (              d   1  )              (  v1  v2  a   a      )
*>    (                  d  )              (  v1  v2  a   a   a  )
*>
*>  where d denotes a diagonal element of the reduced matrix, a denotes
*>  an element of the original matrix that is unchanged, and vi denotes
*>  an element of the vector defining H(i).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), E( * ), TAU( * ), W( LDW, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IW
      DOUBLE PRECISION   ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DGEMV, DLARFG, DSCAL, DSYMV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Reduce last NB columns of upper triangle
*
         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF( I.LT.N ) THEN
*
*              Update A(1:i,i)
*
               CALL DGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ),
     $                     LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL DGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ),
     $                     LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
            END IF
            IF( I.GT.1 ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(1:i-2,i)
*
               CALL DLARFG( I-1, A( I-1, I ), A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = A( I-1, I )
               A( I-1, I ) = ONE
*
*              Compute W(1:i-1,i)
*
               CALL DSYMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1,
     $                     ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
                  CALL DGEMV( 'Transpose', I-1, N-I, ONE, W( 1, IW+1 ),
     $                        LDW, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
                  CALL DGEMV( 'Transpose', I-1, N-I, ONE, A( 1, I+1 ),
     $                        LDA, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
               END IF
               CALL DSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
               ALPHA = -HALF*TAU( I-1 )*DDOT( I-1, W( 1, IW ), 1,
     $                 A( 1, I ), 1 )
               CALL DAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF
*
   10    CONTINUE
      ELSE
*
*        Reduce first NB columns of lower triangle
*
         DO 20 I = 1, NB
*
*           Update A(i:n,i)
*
            CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ),
     $                  LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ),
     $                  LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            IF( I.LT.N ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(i+2:n,i)
*
               CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                      TAU( I ) )
               E( I ) = A( I+1, I )
               A( I+1, I ) = ONE
*
*              Compute W(i+1:n,i)
*
               CALL DSYMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', N-I, I-1, ONE, W( I+1, 1 ), LDW,
     $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ),
     $                     LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ),
     $                     LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
               ALPHA = -HALF*TAU( I )*DDOT( N-I, W( I+1, I ), 1,
     $                 A( I+1, I ), 1 )
               CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF
*
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of DLATRD
*
      END
*> \brief \b DORG2L generates all or part of the orthogonal matrix Q from a QL factorization determined by sgeqlf (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DORG2L + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorg2l.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorg2l.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorg2l.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, K, LDA, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DORG2L generates an m by n real matrix Q with orthonormal columns,
*> which is defined as the last n columns of a product of k elementary
*> reflectors of order m
*>
*>       Q  =  H(k) . . . H(2) H(1)
*>
*> as returned by DGEQLF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix Q. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix Q. M >= N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines the
*>          matrix Q. N >= K >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the (n-k+i)-th column must contain the vector which
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as
*>          returned by DGEQLF in the last k columns of its array
*>          argument A.
*>          On exit, the m by n matrix Q.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The first dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by DGEQLF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument has an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleOTHERcomputational
*
*  =====================================================================
      SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, II, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2L', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns 1:n-k to columns of the unit matrix
*
      DO 20 J = 1, N - K
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( M-N+J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = 1, K
         II = N - K + I
*
*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
*
         A( M-N+II, II ) = ONE
         CALL DLARF( 'Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A,
     $               LDA, WORK )
         CALL DSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
         A( M-N+II, II ) = ONE - TAU( I )
*
*        Set A(m-k+i+1:m,n-k+i) to zero
*
         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of DORG2L
*
      END
*> \brief \b DORG2R generates all or part of the orthogonal matrix Q from a QR factorization determined by sgeqrf (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DORG2R + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorg2r.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorg2r.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorg2r.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, K, LDA, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DORG2R generates an m by n real matrix Q with orthonormal columns,
*> which is defined as the first n columns of a product of k elementary
*> reflectors of order m
*>
*>       Q  =  H(1) H(2) . . . H(k)
*>
*> as returned by DGEQRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix Q. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix Q. M >= N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines the
*>          matrix Q. N >= K >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the i-th column must contain the vector which
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as
*>          returned by DGEQRF in the first k columns of its array
*>          argument A.
*>          On exit, the m-by-n matrix Q.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The first dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by DGEQRF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument has an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleOTHERcomputational
*
*  =====================================================================
      SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns k+1:n to columns of the unit matrix
*
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = K, 1, -1
*
*        Apply H(i) to A(i:m,i:n) from the left
*
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M )
     $      CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
*
*        Set A(1:i-1,i) to zero
*
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of DORG2R
*
      END
*> \brief \b DORGQL
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DORGQL + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgql.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgql.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgql.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, K, LDA, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DORGQL generates an M-by-N real matrix Q with orthonormal columns,
*> which is defined as the last N columns of a product of K elementary
*> reflectors of order M
*>
*>       Q  =  H(k) . . . H(2) H(1)
*>
*> as returned by DGEQLF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix Q. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix Q. M >= N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines the
*>          matrix Q. N >= K >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the (n-k+i)-th column must contain the vector which
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as
*>          returned by DGEQLF in the last k columns of its array
*>          argument A.
*>          On exit, the M-by-N matrix Q.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The first dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by DGEQLF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= max(1,N).
*>          For optimum performance LWORK >= N*NB, where NB is the
*>          optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument has an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup doubleOTHERcomputational
*
*  =====================================================================
      SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KK, L, LDWORK, LWKOPT,
     $                   NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORG2L, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
            NB = ILAENV( 1, 'DORGQL', ' ', M, N, K, -1 )
            LWKOPT = N*NB
         END IF
         WORK( 1 ) = LWKOPT
*
         IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DORGQL', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQL', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the first block.
*        The last kk columns are handled by the block method.
*
         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
*
*        Set A(m-kk+1:m,1:n-kk) to zero.
*
         DO 20 J = 1, N - KK
            DO 10 I = M - KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the first or only block.
*
      CALL DORG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            IF( N-K+I.GT.1 ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i+ib-1) . . . H(i+1) H(i)
*
               CALL DLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB,
     $                      A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
*
               CALL DLARFB( 'Left', 'No transpose', 'Backward',
     $                      'Columnwise', M-K+I+IB-1, N-K+I-1, IB,
     $                      A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA,
     $                      WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H to rows 1:m-k+i+ib-1 of current block
*
            CALL DORG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA,
     $                   TAU( I ), WORK, IINFO )
*
*           Set rows m-k+i+ib:m of current block to zero
*
            DO 40 J = N - K + I, N - K + I + IB - 1
               DO 30 L = M - K + I + IB, M
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGQL
*
      END
*> \brief \b DORGQR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DORGQR + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgqr.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgqr.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgqr.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, K, LDA, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DORGQR generates an M-by-N real matrix Q with orthonormal columns,
*> which is defined as the first N columns of a product of K elementary
*> reflectors of order M
*>
*>       Q  =  H(1) H(2) . . . H(k)
*>
*> as returned by DGEQRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix Q. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix Q. M >= N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines the
*>          matrix Q. N >= K >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the i-th column must contain the vector which
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as
*>          returned by DGEQRF in the first k columns of its array
*>          argument A.
*>          On exit, the M-by-N matrix Q.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The first dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by DGEQRF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= max(1,N).
*>          For optimum performance LWORK >= N*NB, where NB is the
*>          optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument has an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup doubleOTHERcomputational
*
*  =====================================================================
      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the last block.
*        The first kk columns are handled by the block method.
*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
*
*        Set A(1:kk,kk+1:n) to zero.
*
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the last or only block.
*
      IF( KK.LT.N )
     $   CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
     $                TAU( KK+1 ), WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(i:m,i+ib:n) from the left
*
               CALL DLARFB( 'Left', 'No transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H to rows i:m of current block
*
            CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
*
*           Set rows 1:i-1 of current block to zero
*
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGQR
*
      END
*> \brief \b DORGTR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DORGTR + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgtr.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgtr.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgtr.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DORGTR generates a real orthogonal matrix Q which is defined as the
*> product of n-1 elementary reflectors of order N, as returned by
*> DSYTRD:
*>
*> if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
*>
*> if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U': Upper triangle of A contains elementary reflectors
*>                 from DSYTRD;
*>          = 'L': Lower triangle of A contains elementary reflectors
*>                 from DSYTRD.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix Q. N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the vectors which define the elementary reflectors,
*>          as returned by DSYTRD.
*>          On exit, the N-by-N orthogonal matrix Q.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (N-1)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by DSYTRD.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= max(1,N-1).
*>          For optimum performance LWORK >= (N-1)*NB, where NB is
*>          the optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup doubleOTHERcomputational
*
*  =====================================================================
      SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, J, LWKOPT, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DORGQL, DORGQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N-1 ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( UPPER ) THEN
            NB = ILAENV( 1, 'DORGQL', ' ', N-1, N-1, N-1, -1 )
         ELSE
            NB = ILAENV( 1, 'DORGQR', ' ', N-1, N-1, N-1, -1 )
         END IF
         LWKOPT = MAX( 1, N-1 )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGTR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Q was determined by a call to DSYTRD with UPLO = 'U'
*
*        Shift the vectors which define the elementary reflectors one
*        column to the left, and set the last row and column of Q to
*        those of the unit matrix
*
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               A( I, J ) = A( I, J+1 )
   10       CONTINUE
            A( N, J ) = ZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            A( I, N ) = ZERO
   30    CONTINUE
         A( N, N ) = ONE
*
*        Generate Q(1:n-1,1:n-1)
*
         CALL DORGQL( N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO )
*
      ELSE
*
*        Q was determined by a call to DSYTRD with UPLO = 'L'.
*
*        Shift the vectors which define the elementary reflectors one
*        column to the right, and set the first row and column of Q to
*        those of the unit matrix
*
         DO 50 J = N, 2, -1
            A( 1, J ) = ZERO
            DO 40 I = J + 1, N
               A( I, J ) = A( I, J-1 )
   40       CONTINUE
   50    CONTINUE
         A( 1, 1 ) = ONE
         DO 60 I = 2, N
            A( I, 1 ) = ZERO
   60    CONTINUE
         IF( N.GT.1 ) THEN
*
*           Generate Q(2:n,2:n)
*
            CALL DORGQR( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK,
     $                   LWORK, IINFO )
         END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DORGTR
*
      END
*> \brief \b DORM2R multiplies a general matrix by the orthogonal matrix from a QR factorization determined by sgeqrf (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DORM2R + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorm2r.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorm2r.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorm2r.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
*                          WORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          SIDE, TRANS
*       INTEGER            INFO, K, LDA, LDC, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DORM2R overwrites the general real m by n matrix C with
*>
*>       Q * C  if SIDE = 'L' and TRANS = 'N', or
*>
*>       Q**T* C  if SIDE = 'L' and TRANS = 'T', or
*>
*>       C * Q  if SIDE = 'R' and TRANS = 'N', or
*>
*>       C * Q**T if SIDE = 'R' and TRANS = 'T',
*>
*> where Q is a real orthogonal matrix defined as the product of k
*> elementary reflectors
*>
*>       Q = H(1) H(2) . . . H(k)
*>
*> as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
*> if SIDE = 'R'.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**T from the Left
*>          = 'R': apply Q or Q**T from the Right
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N': apply Q  (No transpose)
*>          = 'T': apply Q**T (Transpose)
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines
*>          the matrix Q.
*>          If SIDE = 'L', M >= K >= 0;
*>          if SIDE = 'R', N >= K >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,K)
*>          The i-th column must contain the vector which defines the
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by
*>          DGEQRF in the first k columns of its array argument A.
*>          A is modified by the routine but restored on exit.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*>          If SIDE = 'L', LDA >= max(1,M);
*>          if SIDE = 'R', LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by DGEQRF.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension (LDC,N)
*>          On entry, the m by n matrix C.
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension
*>                                   (N) if SIDE = 'L',
*>                                   (M) if SIDE = 'R'
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleOTHERcomputational
*
*  =====================================================================
      SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      DOUBLE PRECISION   AII
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ is the order of Q
*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORM2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) )
     $     THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
*
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
*
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
*
*           H(i) is applied to C(i:m,1:n)
*
            MI = M - I + 1
            IC = I
         ELSE
*
*           H(i) is applied to C(1:m,i:n)
*
            NI = N - I + 1
            JC = I
         END IF
*
*        Apply H(i)
*
         AII = A( I, I )
         A( I, I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ),
     $               LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
*
*     End of DORM2R
*
      END
*> \brief \b DORMQR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DORMQR + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormqr.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormqr.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormqr.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
*                          WORK, LWORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          SIDE, TRANS
*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DORMQR overwrites the general real M-by-N matrix C with
*>
*>                 SIDE = 'L'     SIDE = 'R'
*> TRANS = 'N':      Q * C          C * Q
*> TRANS = 'T':      Q**T * C       C * Q**T
*>
*> where Q is a real orthogonal matrix defined as the product of k
*> elementary reflectors
*>
*>       Q = H(1) H(2) . . . H(k)
*>
*> as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
*> if SIDE = 'R'.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**T from the Left;
*>          = 'R': apply Q or Q**T from the Right.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N':  No transpose, apply Q;
*>          = 'T':  Transpose, apply Q**T.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines
*>          the matrix Q.
*>          If SIDE = 'L', M >= K >= 0;
*>          if SIDE = 'R', N >= K >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,K)
*>          The i-th column must contain the vector which defines the
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by
*>          DGEQRF in the first k columns of its array argument A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*>          If SIDE = 'L', LDA >= max(1,M);
*>          if SIDE = 'R', LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by DGEQRF.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension (LDC,N)
*>          On entry, the M-by-N matrix C.
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>          If SIDE = 'L', LWORK >= max(1,N);
*>          if SIDE = 'R', LWORK >= max(1,M).
*>          For good performance, LWORK should generally be larger.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup doubleOTHERcomputational
*
*  =====================================================================
      SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1,
     $                     TSIZE = LDT*NBMAX )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWT, JC, LDWORK,
     $                   LWKOPT, MI, NB, NBMIN, NI, NQ, NW
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORM2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Compute the workspace requirements
*
         NB = MIN( NBMAX, ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N, K,
     $        -1 ) )
         LWKOPT = MAX( 1, NW )*NB + TSIZE
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IF( LWORK.LT.NW*NB+TSIZE ) THEN
            NB = (LWORK-TSIZE) / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMQR', SIDE // TRANS, M, N, K,
     $              -1 ) )
         END IF
      END IF
*
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
*
*        Use unblocked code
*
         CALL DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,
     $                IINFO )
      ELSE
*
*        Use blocked code
*
         IWT = 1 + NW*NB
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR.
     $       ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
*
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
*
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL DLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ),
     $                   LDA, TAU( I ), WORK( IWT ), LDT )
            IF( LEFT ) THEN
*
*              H or H**T is applied to C(i:m,1:n)
*
               MI = M - I + 1
               IC = I
            ELSE
*
*              H or H**T is applied to C(1:m,i:n)
*
               NI = N - I + 1
               JC = I
            END IF
*
*           Apply H or H**T
*
            CALL DLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI,
     $                   IB, A( I, I ), LDA, WORK( IWT ), LDT,
     $                   C( IC, JC ), LDC, WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DORMQR
*
      END
*> \brief \b DSTEQR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DSTEQR + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsteqr.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsteqr.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsteqr.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          COMPZ
*       INTEGER            INFO, LDZ, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
*> symmetric tridiagonal matrix using the implicit QL or QR method.
*> The eigenvectors of a full or band symmetric matrix can also be found
*> if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
*> tridiagonal form.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] COMPZ
*> \verbatim
*>          COMPZ is CHARACTER*1
*>          = 'N':  Compute eigenvalues only.
*>          = 'V':  Compute eigenvalues and eigenvectors of the original
*>                  symmetric matrix.  On entry, Z must contain the
*>                  orthogonal matrix used to reduce the original matrix
*>                  to tridiagonal form.
*>          = 'I':  Compute eigenvalues and eigenvectors of the
*>                  tridiagonal matrix.  Z is initialized to the identity
*>                  matrix.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (N)
*>          On entry, the diagonal elements of the tridiagonal matrix.
*>          On exit, if INFO = 0, the eigenvalues in ascending order.
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (N-1)
*>          On entry, the (n-1) subdiagonal elements of the tridiagonal
*>          matrix.
*>          On exit, E has been destroyed.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
*>          On entry, if  COMPZ = 'V', then Z contains the orthogonal
*>          matrix used in the reduction to tridiagonal form.
*>          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
*>          orthonormal eigenvectors of the original symmetric matrix,
*>          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
*>          of the symmetric tridiagonal matrix.
*>          If COMPZ = 'N', then Z is not referenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.  LDZ >= 1, and if
*>          eigenvectors are desired, then  LDZ >= max(1,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (max(1,2*N-2))
*>          If COMPZ = 'N', then WORK is not referenced.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  the algorithm has failed to find all the eigenvalues in
*>                a total of 30*N iterations; if INFO = i, then i
*>                elements of E have not converged to zero; on exit, D
*>                and E contain the elements of a symmetric tridiagonal
*>                matrix which is orthogonally similar to the original
*>                matrix.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup auxOTHERcomputational
*
*  =====================================================================
      SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND,
     $                   LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1,
     $                   NM1, NMAXIT
      DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2,
     $                   S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASET, DLASR,
     $                   DLASRT, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEQR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Determine the unit roundoff and over/underflow thresholds.
*
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
*
*     Compute the eigenvalues and eigenvectors of the tridiagonal
*     matrix.
*
      IF( ICOMPZ.EQ.2 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
*
      NMAXIT = N*MAXIT
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
      NM1 = N - 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 160
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $          1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = DLANST( 'M', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
*
*     Choose between QL and QR iteration
*
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
*
      IF( LEND.GT.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+
     $             SAFMIN )GO TO 60
   50       CONTINUE
         END IF
*
         M = LEND
*
   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 80
*
*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL DLASR( 'R', 'V', 'B', N, 2, WORK( L ),
     $                     WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 40
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )
     $         E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
*
   70    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL DLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ),
     $                  Z( 1, L ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
*
*        Eigenvalue found.
*
   80    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 40
         GO TO 140
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+
     $             SAFMIN )GO TO 110
  100       CONTINUE
         END IF
*
         M = LEND
*
  110    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 130
*
*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL DLASR( 'R', 'V', 'F', N, 2, WORK( M ),
     $                     WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 90
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M )
     $         E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
*
  120    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL DLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ),
     $                  Z( 1, M ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
*
*        Eigenvalue found.
*
  130    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 90
         GO TO 140
*
      END IF
*
*     Undo scaling if necessary
*
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      END IF
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.LT.NMAXIT )
     $   GO TO 10
      DO 150 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  150 CONTINUE
      GO TO 190
*
*     Order eigenvalues and eigenvectors.
*
  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN
*
*        Use Quick Sort
*
         CALL DLASRT( 'I', N, D, INFO )
*
      ELSE
*
*        Use Selection Sort to minimize swaps of eigenvectors
*
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
*
  190 CONTINUE
      RETURN
*
*     End of DSTEQR
*
      END
*> \brief \b DSTERF
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DSTERF + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsterf.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsterf.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsterf.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSTERF( N, D, E, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   D( * ), E( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
*> using the Pal-Walker-Kahan variant of the QL or QR algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (N)
*>          On entry, the n diagonal elements of the tridiagonal matrix.
*>          On exit, if INFO = 0, the eigenvalues in ascending order.
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (N-1)
*>          On entry, the (n-1) subdiagonal elements of the tridiagonal
*>          matrix.
*>          On exit, E has been destroyed.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  the algorithm failed to find all of the eigenvalues in
*>                a total of 30*N iterations; if INFO = i, then i
*>                elements of E have not converged to zero.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup auxOTHERcomputational
*
*  =====================================================================
      SUBROUTINE DSTERF( N, D, E, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M,
     $                   NMAXIT
      DOUBLE PRECISION   ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC,
     $                   OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN,
     $                   SIGMA, SSFMAX, SSFMIN, RMAX
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           DLAMCH, DLANST, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAE2, DLASCL, DLASRT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DSTERF', -INFO )
         RETURN
      END IF
      IF( N.LE.1 )
     $   RETURN
*
*     Determine the unit roundoff for this environment.
*
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
      RMAX = DLAMCH( 'O' )
*
*     Compute the eigenvalues of the tridiagonal matrix.
*
      NMAXIT = N*MAXIT
      SIGMA = ZERO
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 170
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      DO 20 M = L1, N - 1
         IF( ABS( E( M ) ).LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $       1 ) ) ) )*EPS ) THEN
            E( M ) = ZERO
            GO TO 30
         END IF
   20 CONTINUE
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = DLANST( 'M', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10      
      IF( (ANORM.GT.SSFMAX) ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
*
      DO 40 I = L, LEND - 1
         E( I ) = E( I )**2
   40 CONTINUE
*
*     Choose between QL and QR iteration
*
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
*
      IF( LEND.GE.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   50    CONTINUE
         IF( L.NE.LEND ) THEN
            DO 60 M = L, LEND - 1
               IF( ABS( E( M ) ).LE.EPS2*ABS( D( M )*D( M+1 ) ) )
     $            GO TO 70
   60       CONTINUE
         END IF
         M = LEND
*
   70    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 90
*
*        If remaining matrix is 2 by 2, use DLAE2 to compute its
*        eigenvalues.
*
         IF( M.EQ.L+1 ) THEN
            RTE = SQRT( E( L ) )
            CALL DLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 50
            GO TO 150
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        Form shift.
*
         RTE = SQRT( E( L ) )
         SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        Inner loop
*
         DO 80 I = M - 1, L, -1
            BB = E( I )
            R = P + BB
            IF( I.NE.M-1 )
     $         E( I+1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
   80    CONTINUE
*
         E( L ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 50
*
*        Eigenvalue found.
*
   90    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 50
         GO TO 150
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
  100    CONTINUE
         DO 110 M = L, LEND + 1, -1
            IF( ABS( E( M-1 ) ).LE.EPS2*ABS( D( M )*D( M-1 ) ) )
     $         GO TO 120
  110    CONTINUE
         M = LEND
*
  120    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 140
*
*        If remaining matrix is 2 by 2, use DLAE2 to compute its
*        eigenvalues.
*
         IF( M.EQ.L-1 ) THEN
            RTE = SQRT( E( L-1 ) )
            CALL DLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
            D( L ) = RT1
            D( L-1 ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 100
            GO TO 150
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        Form shift.
*
         RTE = SQRT( E( L-1 ) )
         SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        Inner loop
*
         DO 130 I = M, L - 1
            BB = E( I )
            R = P + BB
            IF( I.NE.M )
     $         E( I-1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I+1 )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
  130    CONTINUE
*
         E( L-1 ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 100
*
*        Eigenvalue found.
*
  140    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 100
         GO TO 150
*
      END IF
*
*     Undo scaling if necessary
*
  150 CONTINUE
      IF( ISCALE.EQ.1 )
     $   CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
      IF( ISCALE.EQ.2 )
     $   CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.LT.NMAXIT )
     $   GO TO 10
      DO 160 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  160 CONTINUE
      GO TO 180
*
*     Sort eigenvalues in increasing order.
*
  170 CONTINUE
      CALL DLASRT( 'I', N, D, INFO )
*
  180 CONTINUE
      RETURN
*
*     End of DSTERF
*
      END
*> \brief <b> DSYEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DSYEV + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyev.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyev.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyev.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, UPLO
*       INTEGER            INFO, LDA, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYEV computes all eigenvalues and, optionally, eigenvectors of a
*> real symmetric matrix A.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBZ
*> \verbatim
*>          JOBZ is CHARACTER*1
*>          = 'N':  Compute eigenvalues only;
*>          = 'V':  Compute eigenvalues and eigenvectors.
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangle of A is stored;
*>          = 'L':  Lower triangle of A is stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA, N)
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the
*>          leading N-by-N upper triangular part of A contains the
*>          upper triangular part of the matrix A.  If UPLO = 'L',
*>          the leading N-by-N lower triangular part of A contains
*>          the lower triangular part of the matrix A.
*>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*>          orthonormal eigenvectors of the matrix A.
*>          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*>          or the upper triangle (if UPLO='U') of A, including the
*>          diagonal, is destroyed.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is DOUBLE PRECISION array, dimension (N)
*>          If INFO = 0, the eigenvalues in ascending order.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of the array WORK.  LWORK >= max(1,3*N-1).
*>          For optimal efficiency, LWORK >= (NB+2)*N,
*>          where NB is the blocksize for DSYTRD returned by ILAENV.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, the algorithm failed to converge; i
*>                off-diagonal elements of an intermediate tridiagonal
*>                form did not converge to zero.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup doubleSYeigen
*
*  =====================================================================
      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE,
     $                   LLWORK, LWKOPT, NB
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANSY
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASCL, DORGTR, DSCAL, DSTEQR, DSTERF, DSYTRD,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
*
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( 1, ( NB+2 )*N )
         WORK( 1 ) = LWKOPT
*
         IF( LWORK.LT.MAX( 1, 3*N-1 ) .AND. .NOT.LQUERY )
     $      INFO = -8
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         RETURN
      END IF
*
      IF( N.EQ.1 ) THEN
         W( 1 ) = A( 1, 1 )
         WORK( 1 ) = 2
         IF( WANTZ )
     $      A( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 )
     $   CALL DLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
*
*     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
*
      INDE = 1
      INDTAU = INDE + N
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      CALL DSYTRD( UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ),
     $             WORK( INDWRK ), LLWORK, IINFO )
*
*     For eigenvalues only, call DSTERF.  For eigenvectors, first call
*     DORGTR to generate the orthogonal matrix, then call DSTEQR.
*
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, WORK( INDE ), INFO )
      ELSE
         CALL DORGTR( UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ),
     $                LLWORK, IINFO )
         CALL DSTEQR( JOBZ, N, W, WORK( INDE ), A, LDA, WORK( INDTAU ),
     $                INFO )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
*     Set WORK(1) to optimal workspace size.
*
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of DSYEV
*
      END
*> \brief \b DSYTD2 reduces a symmetric matrix to real symmetric tridiagonal form by an orthogonal similarity transformation (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DSYTD2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytd2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytd2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytd2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal
*> form T by an orthogonal similarity transformation: Q**T * A * Q = T.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          symmetric matrix A is stored:
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*>          n-by-n upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading n-by-n lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the strictly upper
*>          triangular part of A is not referenced.
*>          On exit, if UPLO = 'U', the diagonal and first superdiagonal
*>          of A are overwritten by the corresponding elements of the
*>          tridiagonal matrix T, and the elements above the first
*>          superdiagonal, with the array TAU, represent the orthogonal
*>          matrix Q as a product of elementary reflectors; if UPLO
*>          = 'L', the diagonal and first subdiagonal of A are over-
*>          written by the corresponding elements of the tridiagonal
*>          matrix T, and the elements below the first subdiagonal, with
*>          the array TAU, represent the orthogonal matrix Q as a product
*>          of elementary reflectors. See Further Details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (N)
*>          The diagonal elements of the tridiagonal matrix T:
*>          D(i) = A(i,i).
*> \endverbatim
*>
*> \param[out] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (N-1)
*>          The off-diagonal elements of the tridiagonal matrix T:
*>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (N-1)
*>          The scalar factors of the elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleSYcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  If UPLO = 'U', the matrix Q is represented as a product of elementary
*>  reflectors
*>
*>     Q = H(n-1) . . . H(2) H(1).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**T
*>
*>  where tau is a real scalar, and v is a real vector with
*>  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*>  A(1:i-1,i+1), and tau in TAU(i).
*>
*>  If UPLO = 'L', the matrix Q is represented as a product of elementary
*>  reflectors
*>
*>     Q = H(1) H(2) . . . H(n-1).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**T
*>
*>  where tau is a real scalar, and v is a real vector with
*>  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*>  and tau in TAU(i).
*>
*>  The contents of A on exit are illustrated by the following examples
*>  with n = 5:
*>
*>  if UPLO = 'U':                       if UPLO = 'L':
*>
*>    (  d   e   v2  v3  v4 )              (  d                  )
*>    (      d   e   v3  v4 )              (  e   d              )
*>    (          d   e   v4 )              (  v1  e   d          )
*>    (              d   e  )              (  v1  v2  e   d      )
*>    (                  d  )              (  v1  v2  v3  e   d  )
*>
*>  where d and e denote diagonal and off-diagonal elements of T, and vi
*>  denotes an element of the vector defining H(i).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO, HALF
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0,
     $                   HALF = 1.0D0 / 2.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      DOUBLE PRECISION   ALPHA, TAUI
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DLARFG, DSYMV, DSYR2, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTD2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A
*
         DO 10 I = N - 1, 1, -1
*
*           Generate elementary reflector H(i) = I - tau * v * v**T
*           to annihilate A(1:i-1,i+1)
*
            CALL DLARFG( I, A( I, I+1 ), A( 1, I+1 ), 1, TAUI )
            E( I ) = A( I, I+1 )
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(1:i,1:i)
*
               A( I, I+1 ) = ONE
*
*              Compute  x := tau * A * v  storing x in TAU(1:i)
*
               CALL DSYMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO,
     $                     TAU, 1 )
*
*              Compute  w := x - 1/2 * tau * (x**T * v) * v
*
               ALPHA = -HALF*TAUI*DDOT( I, TAU, 1, A( 1, I+1 ), 1 )
               CALL DAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w**T - w * v**T
*
               CALL DSYR2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A,
     $                     LDA )
*
               A( I, I+1 ) = E( I )
            END IF
            D( I+1 ) = A( I+1, I+1 )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = A( 1, 1 )
      ELSE
*
*        Reduce the lower triangle of A
*
         DO 20 I = 1, N - 1
*
*           Generate elementary reflector H(i) = I - tau * v * v**T
*           to annihilate A(i+2:n,i)
*
            CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                   TAUI )
            E( I ) = A( I+1, I )
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(i+1:n,i+1:n)
*
               A( I+1, I ) = ONE
*
*              Compute  x := tau * A * v  storing y in TAU(i:n-1)
*
               CALL DSYMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, TAU( I ), 1 )
*
*              Compute  w := x - 1/2 * tau * (x**T * v) * v
*
               ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), 1, A( I+1, I ),
     $                 1 )
               CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w**T - w * v**T
*
               CALL DSYR2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1,
     $                     A( I+1, I+1 ), LDA )
*
               A( I+1, I ) = E( I )
            END IF
            D( I ) = A( I, I )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = A( N, N )
      END IF
*
      RETURN
*
*     End of DSYTD2
*
      END
*> \brief \b DSYTRD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DSYTRD + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrd.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrd.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrd.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ),
*      $                   WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYTRD reduces a real symmetric matrix A to real symmetric
*> tridiagonal form T by an orthogonal similarity transformation:
*> Q**T * A * Q = T.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangle of A is stored;
*>          = 'L':  Lower triangle of A is stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*>          N-by-N upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading N-by-N lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the strictly upper
*>          triangular part of A is not referenced.
*>          On exit, if UPLO = 'U', the diagonal and first superdiagonal
*>          of A are overwritten by the corresponding elements of the
*>          tridiagonal matrix T, and the elements above the first
*>          superdiagonal, with the array TAU, represent the orthogonal
*>          matrix Q as a product of elementary reflectors; if UPLO
*>          = 'L', the diagonal and first subdiagonal of A are over-
*>          written by the corresponding elements of the tridiagonal
*>          matrix T, and the elements below the first subdiagonal, with
*>          the array TAU, represent the orthogonal matrix Q as a product
*>          of elementary reflectors. See Further Details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (N)
*>          The diagonal elements of the tridiagonal matrix T:
*>          D(i) = A(i,i).
*> \endverbatim
*>
*> \param[out] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (N-1)
*>          The off-diagonal elements of the tridiagonal matrix T:
*>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (N-1)
*>          The scalar factors of the elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= 1.
*>          For optimum performance LWORK >= N*NB, where NB is the
*>          optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup doubleSYcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  If UPLO = 'U', the matrix Q is represented as a product of elementary
*>  reflectors
*>
*>     Q = H(n-1) . . . H(2) H(1).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**T
*>
*>  where tau is a real scalar, and v is a real vector with
*>  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*>  A(1:i-1,i+1), and tau in TAU(i).
*>
*>  If UPLO = 'L', the matrix Q is represented as a product of elementary
*>  reflectors
*>
*>     Q = H(1) H(2) . . . H(n-1).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**T
*>
*>  where tau is a real scalar, and v is a real vector with
*>  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*>  and tau in TAU(i).
*>
*>  The contents of A on exit are illustrated by the following examples
*>  with n = 5:
*>
*>  if UPLO = 'U':                       if UPLO = 'L':
*>
*>    (  d   e   v2  v3  v4 )              (  d                  )
*>    (      d   e   v3  v4 )              (  e   d              )
*>    (          d   e   v4 )              (  v1  e   d          )
*>    (              d   e  )              (  v1  v2  e   d      )
*>    (                  d  )              (  v1  v2  v3  e   d  )
*>
*>  where d and e denote diagonal and off-diagonal elements of T, and vi
*>  denotes an element of the vector defining H(i).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ),
     $                   WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLATRD, DSYR2K, DSYTD2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size.
*
         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NX = N
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
*
*        Determine when to cross over from blocked to unblocked code
*        (last block is always handled by unblocked code).
*
         NX = MAX( NB, ILAENV( 3, 'DSYTRD', UPLO, N, -1, -1, -1 ) )
         IF( NX.LT.N ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  determine the
*              minimum value of NB, and reduce NB or force use of
*              unblocked code by setting NX = N.
*
               NB = MAX( LWORK / LDWORK, 1 )
               NBMIN = ILAENV( 2, 'DSYTRD', UPLO, N, -1, -1, -1 )
               IF( NB.LT.NBMIN )
     $            NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
*
      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A.
*        Columns 1:kk are handled by the unblocked method.
*
         KK = N - ( ( N-NX+NB-1 ) / NB )*NB
         DO 20 I = N - NB + 1, KK + 1, -NB
*
*           Reduce columns i:i+nb-1 to tridiagonal form and form the
*           matrix W which is needed to update the unreduced part of
*           the matrix
*
            CALL DLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK,
     $                   LDWORK )
*
*           Update the unreduced submatrix A(1:i-1,1:i-1), using an
*           update of the form:  A := A - V*W**T - W*V**T
*
            CALL DSYR2K( UPLO, 'No transpose', I-1, NB, -ONE, A( 1, I ),
     $                   LDA, WORK, LDWORK, ONE, A, LDA )
*
*           Copy superdiagonal elements back into A, and diagonal
*           elements into D
*
            DO 10 J = I, I + NB - 1
               A( J-1, J ) = E( J-1 )
               D( J ) = A( J, J )
   10       CONTINUE
   20    CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
         CALL DSYTD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
      ELSE
*
*        Reduce the lower triangle of A
*
         DO 40 I = 1, N - NX, NB
*
*           Reduce columns i:i+nb-1 to tridiagonal form and form the
*           matrix W which is needed to update the unreduced part of
*           the matrix
*
            CALL DLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ),
     $                   TAU( I ), WORK, LDWORK )
*
*           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
*           an update of the form:  A := A - V*W**T - W*V**T
*
            CALL DSYR2K( UPLO, 'No transpose', N-I-NB+1, NB, -ONE,
     $                   A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE,
     $                   A( I+NB, I+NB ), LDA )
*
*           Copy subdiagonal elements back into A, and diagonal
*           elements into D
*
            DO 30 J = I, I + NB - 1
               A( J+1, J ) = E( J )
               D( J ) = A( J, J )
   30       CONTINUE
   40    CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
         CALL DSYTD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ),
     $                TAU( I ), IINFO )
      END IF
*
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DSYTRD
*
      END
*> \brief \b DTGEVC
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DTGEVC + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtgevc.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtgevc.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtgevc.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL,
*                          LDVL, VR, LDVR, MM, M, WORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          HOWMNY, SIDE
*       INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N
*       ..
*       .. Array Arguments ..
*       LOGICAL            SELECT( * )
*       DOUBLE PRECISION   P( LDP, * ), S( LDS, * ), VL( LDVL, * ),
*      $                   VR( LDVR, * ), WORK( * )
*       ..
*  
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTGEVC computes some or all of the right and/or left eigenvectors of
*> a pair of real matrices (S,P), where S is a quasi-triangular matrix
*> and P is upper triangular.  Matrix pairs of this type are produced by
*> the generalized Schur factorization of a matrix pair (A,B):
*>
*>    A = Q*S*Z**T,  B = Q*P*Z**T
*>
*> as computed by DGGHRD + DHGEQZ.
*>
*> The right eigenvector x and the left eigenvector y of (S,P)
*> corresponding to an eigenvalue w are defined by:
*> 
*>    S*x = w*P*x,  (y**H)*S = w*(y**H)*P,
*> 
*> where y**H denotes the conjugate tranpose of y.
*> The eigenvalues are not input to this routine, but are computed
*> directly from the diagonal blocks of S and P.
*> 
*> This routine returns the matrices X and/or Y of right and left
*> eigenvectors of (S,P), or the products Z*X and/or Q*Y,
*> where Z and Q are input matrices.
*> If Q and Z are the orthogonal factors from the generalized Schur
*> factorization of a matrix pair (A,B), then Z*X and Q*Y
*> are the matrices of right and left eigenvectors of (A,B).
*> 
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'R': compute right eigenvectors only;
*>          = 'L': compute left eigenvectors only;
*>          = 'B': compute both right and left eigenvectors.
*> \endverbatim
*>
*> \param[in] HOWMNY
*> \verbatim
*>          HOWMNY is CHARACTER*1
*>          = 'A': compute all right and/or left eigenvectors;
*>          = 'B': compute all right and/or left eigenvectors,
*>                 backtransformed by the matrices in VR and/or VL;
*>          = 'S': compute selected right and/or left eigenvectors,
*>                 specified by the logical array SELECT.
*> \endverbatim
*>
*> \param[in] SELECT
*> \verbatim
*>          SELECT is LOGICAL array, dimension (N)
*>          If HOWMNY='S', SELECT specifies the eigenvectors to be
*>          computed.  If w(j) is a real eigenvalue, the corresponding
*>          real eigenvector is computed if SELECT(j) is .TRUE..
*>          If w(j) and w(j+1) are the real and imaginary parts of a
*>          complex eigenvalue, the corresponding complex eigenvector
*>          is computed if either SELECT(j) or SELECT(j+1) is .TRUE.,
*>          and on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is
*>          set to .FALSE..
*>          Not referenced if HOWMNY = 'A' or 'B'.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices S and P.  N >= 0.
*> \endverbatim
*>
*> \param[in] S
*> \verbatim
*>          S is DOUBLE PRECISION array, dimension (LDS,N)
*>          The upper quasi-triangular matrix S from a generalized Schur
*>          factorization, as computed by DHGEQZ.
*> \endverbatim
*>
*> \param[in] LDS
*> \verbatim
*>          LDS is INTEGER
*>          The leading dimension of array S.  LDS >= max(1,N).
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is DOUBLE PRECISION array, dimension (LDP,N)
*>          The upper triangular matrix P from a generalized Schur
*>          factorization, as computed by DHGEQZ.
*>          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks
*>          of S must be in positive diagonal form.
*> \endverbatim
*>
*> \param[in] LDP
*> \verbatim
*>          LDP is INTEGER
*>          The leading dimension of array P.  LDP >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] VL
*> \verbatim
*>          VL is DOUBLE PRECISION array, dimension (LDVL,MM)
*>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
*>          contain an N-by-N matrix Q (usually the orthogonal matrix Q
*>          of left Schur vectors returned by DHGEQZ).
*>          On exit, if SIDE = 'L' or 'B', VL contains:
*>          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P);
*>          if HOWMNY = 'B', the matrix Q*Y;
*>          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by
*>                      SELECT, stored consecutively in the columns of
*>                      VL, in the same order as their eigenvalues.
*>
*>          A complex eigenvector corresponding to a complex eigenvalue
*>          is stored in two consecutive columns, the first holding the
*>          real part, and the second the imaginary part.
*>
*>          Not referenced if SIDE = 'R'.
*> \endverbatim
*>
*> \param[in] LDVL
*> \verbatim
*>          LDVL is INTEGER
*>          The leading dimension of array VL.  LDVL >= 1, and if
*>          SIDE = 'L' or 'B', LDVL >= N.
*> \endverbatim
*>
*> \param[in,out] VR
*> \verbatim
*>          VR is DOUBLE PRECISION array, dimension (LDVR,MM)
*>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
*>          contain an N-by-N matrix Z (usually the orthogonal matrix Z
*>          of right Schur vectors returned by DHGEQZ).
*>
*>          On exit, if SIDE = 'R' or 'B', VR contains:
*>          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P);
*>          if HOWMNY = 'B' or 'b', the matrix Z*X;
*>          if HOWMNY = 'S' or 's', the right eigenvectors of (S,P)
*>                      specified by SELECT, stored consecutively in the
*>                      columns of VR, in the same order as their
*>                      eigenvalues.
*>
*>          A complex eigenvector corresponding to a complex eigenvalue
*>          is stored in two consecutive columns, the first holding the
*>          real part and the second the imaginary part.
*>          
*>          Not referenced if SIDE = 'L'.
*> \endverbatim
*>
*> \param[in] LDVR
*> \verbatim
*>          LDVR is INTEGER
*>          The leading dimension of the array VR.  LDVR >= 1, and if
*>          SIDE = 'R' or 'B', LDVR >= N.
*> \endverbatim
*>
*> \param[in] MM
*> \verbatim
*>          MM is INTEGER
*>          The number of columns in the arrays VL and/or VR. MM >= M.
*> \endverbatim
*>
*> \param[out] M
*> \verbatim
*>          M is INTEGER
*>          The number of columns in the arrays VL and/or VR actually
*>          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M
*>          is set to N.  Each selected real eigenvector occupies one
*>          column and each selected complex eigenvector occupies two
*>          columns.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (6*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  the 2-by-2 block (INFO:INFO+1) does not have a complex
*>                eigenvalue.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup doubleGEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Allocation of workspace:
*>  ---------- -- ---------
*>
*>     WORK( j ) = 1-norm of j-th column of A, above the diagonal
*>     WORK( N+j ) = 1-norm of j-th column of B, above the diagonal
*>     WORK( 2*N+1:3*N ) = real part of eigenvector
*>     WORK( 3*N+1:4*N ) = imaginary part of eigenvector
*>     WORK( 4*N+1:5*N ) = real part of back-transformed eigenvector
*>     WORK( 5*N+1:6*N ) = imaginary part of back-transformed eigenvector
*>
*>  Rowwise vs. columnwise solution methods:
*>  ------- --  ---------- -------- -------
*>
*>  Finding a generalized eigenvector consists basically of solving the
*>  singular triangular system
*>
*>   (A - w B) x = 0     (for right) or:   (A - w B)**H y = 0  (for left)
*>
*>  Consider finding the i-th right eigenvector (assume all eigenvalues
*>  are real). The equation to be solved is:
*>       n                   i
*>  0 = sum  C(j,k) v(k)  = sum  C(j,k) v(k)     for j = i,. . .,1
*>      k=j                 k=j
*>
*>  where  C = (A - w B)  (The components v(i+1:n) are 0.)
*>
*>  The "rowwise" method is:
*>
*>  (1)  v(i) := 1
*>  for j = i-1,. . .,1:
*>                          i
*>      (2) compute  s = - sum C(j,k) v(k)   and
*>                        k=j+1
*>
*>      (3) v(j) := s / C(j,j)
*>
*>  Step 2 is sometimes called the "dot product" step, since it is an
*>  inner product between the j-th row and the portion of the eigenvector
*>  that has been computed so far.
*>
*>  The "columnwise" method consists basically in doing the sums
*>  for all the rows in parallel.  As each v(j) is computed, the
*>  contribution of v(j) times the j-th column of C is added to the
*>  partial sums.  Since FORTRAN arrays are stored columnwise, this has
*>  the advantage that at each step, the elements of C that are accessed
*>  are adjacent to one another, whereas with the rowwise method, the
*>  elements accessed at a step are spaced LDS (and LDP) words apart.
*>
*>  When finding left eigenvectors, the matrix in question is the
*>  transpose of the one in storage, so the rowwise method then
*>  actually accesses columns of A and B at each step, and so is the
*>  preferred method.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL,
     $                   LDVL, VR, LDVR, MM, M, WORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   P( LDP, * ), S( LDS, * ), VL( LDVL, * ),
     $                   VR( LDVR, * ), WORK( * )
*     ..
*
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, SAFETY
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0,
     $                   SAFETY = 1.0D+2 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COMPL, COMPR, IL2BY2, ILABAD, ILALL, ILBACK,
     $                   ILBBAD, ILCOMP, ILCPLX, LSA, LSB
      INTEGER            I, IBEG, IEIG, IEND, IHWMNY, IINFO, IM, ISIDE,
     $                   J, JA, JC, JE, JR, JW, NA, NW
      DOUBLE PRECISION   ACOEF, ACOEFA, ANORM, ASCALE, BCOEFA, BCOEFI,
     $                   BCOEFR, BIG, BIGNUM, BNORM, BSCALE, CIM2A,
     $                   CIM2B, CIMAGA, CIMAGB, CRE2A, CRE2B, CREALA,
     $                   CREALB, DMIN, SAFMIN, SALFAR, SBETA, SCALE,
     $                   SMALL, TEMP, TEMP2, TEMP2I, TEMP2R, ULP, XMAX,
     $                   XSCALE
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   BDIAG( 2 ), SUM( 2, 2 ), SUMS( 2, 2 ),
     $                   SUMP( 2, 2 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DLABAD, DLACPY, DLAG2, DLALN2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Decode and Test the input parameters
*
      IF( LSAME( HOWMNY, 'A' ) ) THEN
         IHWMNY = 1
         ILALL = .TRUE.
         ILBACK = .FALSE.
      ELSE IF( LSAME( HOWMNY, 'S' ) ) THEN
         IHWMNY = 2
         ILALL = .FALSE.
         ILBACK = .FALSE.
      ELSE IF( LSAME( HOWMNY, 'B' ) ) THEN
         IHWMNY = 3
         ILALL = .TRUE.
         ILBACK = .TRUE.
      ELSE
         IHWMNY = -1
         ILALL = .TRUE.
      END IF
*
      IF( LSAME( SIDE, 'R' ) ) THEN
         ISIDE = 1
         COMPL = .FALSE.
         COMPR = .TRUE.
      ELSE IF( LSAME( SIDE, 'L' ) ) THEN
         ISIDE = 2
         COMPL = .TRUE.
         COMPR = .FALSE.
      ELSE IF( LSAME( SIDE, 'B' ) ) THEN
         ISIDE = 3
         COMPL = .TRUE.
         COMPR = .TRUE.
      ELSE
         ISIDE = -1
      END IF
*
      INFO = 0
      IF( ISIDE.LT.0 ) THEN
         INFO = -1
      ELSE IF( IHWMNY.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDS.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDP.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTGEVC', -INFO )
         RETURN
      END IF
*
*     Count the number of eigenvectors to be computed
*
      IF( .NOT.ILALL ) THEN
         IM = 0
         ILCPLX = .FALSE.
         DO 10 J = 1, N
            IF( ILCPLX ) THEN
               ILCPLX = .FALSE.
               GO TO 10
            END IF
            IF( J.LT.N ) THEN
               IF( S( J+1, J ).NE.ZERO )
     $            ILCPLX = .TRUE.
            END IF
            IF( ILCPLX ) THEN
               IF( SELECT( J ) .OR. SELECT( J+1 ) )
     $            IM = IM + 2
            ELSE
               IF( SELECT( J ) )
     $            IM = IM + 1
            END IF
   10    CONTINUE
      ELSE
         IM = N
      END IF
*
*     Check 2-by-2 diagonal blocks of A, B
*
      ILABAD = .FALSE.
      ILBBAD = .FALSE.
      DO 20 J = 1, N - 1
         IF( S( J+1, J ).NE.ZERO ) THEN
            IF( P( J, J ).EQ.ZERO .OR. P( J+1, J+1 ).EQ.ZERO .OR.
     $          P( J, J+1 ).NE.ZERO )ILBBAD = .TRUE.
            IF( J.LT.N-1 ) THEN
               IF( S( J+2, J+1 ).NE.ZERO )
     $            ILABAD = .TRUE.
            END IF
         END IF
   20 CONTINUE
*
      IF( ILABAD ) THEN
         INFO = -5
      ELSE IF( ILBBAD ) THEN
         INFO = -7
      ELSE IF( COMPL .AND. LDVL.LT.N .OR. LDVL.LT.1 ) THEN
         INFO = -10
      ELSE IF( COMPR .AND. LDVR.LT.N .OR. LDVR.LT.1 ) THEN
         INFO = -12
      ELSE IF( MM.LT.IM ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTGEVC', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      M = IM
      IF( N.EQ.0 )
     $   RETURN
*
*     Machine Constants
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      BIG = ONE / SAFMIN
      CALL DLABAD( SAFMIN, BIG )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      SMALL = SAFMIN*N / ULP
      BIG = ONE / SMALL
      BIGNUM = ONE / ( SAFMIN*N )
*
*     Compute the 1-norm of each column of the strictly upper triangular
*     part (i.e., excluding all elements belonging to the diagonal
*     blocks) of A and B to check for possible overflow in the
*     triangular solver.
*
      ANORM = ABS( S( 1, 1 ) )
      IF( N.GT.1 )
     $   ANORM = ANORM + ABS( S( 2, 1 ) )
      BNORM = ABS( P( 1, 1 ) )
      WORK( 1 ) = ZERO
      WORK( N+1 ) = ZERO
*
      DO 50 J = 2, N
         TEMP = ZERO
         TEMP2 = ZERO
         IF( S( J, J-1 ).EQ.ZERO ) THEN
            IEND = J - 1
         ELSE
            IEND = J - 2
         END IF
         DO 30 I = 1, IEND
            TEMP = TEMP + ABS( S( I, J ) )
            TEMP2 = TEMP2 + ABS( P( I, J ) )
   30    CONTINUE
         WORK( J ) = TEMP
         WORK( N+J ) = TEMP2
         DO 40 I = IEND + 1, MIN( J+1, N )
            TEMP = TEMP + ABS( S( I, J ) )
            TEMP2 = TEMP2 + ABS( P( I, J ) )
   40    CONTINUE
         ANORM = MAX( ANORM, TEMP )
         BNORM = MAX( BNORM, TEMP2 )
   50 CONTINUE
*
      ASCALE = ONE / MAX( ANORM, SAFMIN )
      BSCALE = ONE / MAX( BNORM, SAFMIN )
*
*     Left eigenvectors
*
      IF( COMPL ) THEN
         IEIG = 0
*
*        Main loop over eigenvalues
*
         ILCPLX = .FALSE.
         DO 220 JE = 1, N
*
*           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or
*           (b) this would be the second of a complex pair.
*           Check for complex eigenvalue, so as to be sure of which
*           entry(-ies) of SELECT to look at.
*
            IF( ILCPLX ) THEN
               ILCPLX = .FALSE.
               GO TO 220
            END IF
            NW = 1
            IF( JE.LT.N ) THEN
               IF( S( JE+1, JE ).NE.ZERO ) THEN
                  ILCPLX = .TRUE.
                  NW = 2
               END IF
            END IF
            IF( ILALL ) THEN
               ILCOMP = .TRUE.
            ELSE IF( ILCPLX ) THEN
               ILCOMP = SELECT( JE ) .OR. SELECT( JE+1 )
            ELSE
               ILCOMP = SELECT( JE )
            END IF
            IF( .NOT.ILCOMP )
     $         GO TO 220
*
*           Decide if (a) singular pencil, (b) real eigenvalue, or
*           (c) complex eigenvalue.
*
            IF( .NOT.ILCPLX ) THEN
               IF( ABS( S( JE, JE ) ).LE.SAFMIN .AND.
     $             ABS( P( JE, JE ) ).LE.SAFMIN ) THEN
*
*                 Singular matrix pencil -- return unit eigenvector
*
                  IEIG = IEIG + 1
                  DO 60 JR = 1, N
                     VL( JR, IEIG ) = ZERO
   60             CONTINUE
                  VL( IEIG, IEIG ) = ONE
                  GO TO 220
               END IF
            END IF
*
*           Clear vector
*
            DO 70 JR = 1, NW*N
               WORK( 2*N+JR ) = ZERO
   70       CONTINUE
*                                                 T
*           Compute coefficients in  ( a A - b B )  y = 0
*              a  is  ACOEF
*              b  is  BCOEFR + i*BCOEFI
*
            IF( .NOT.ILCPLX ) THEN
*
*              Real eigenvalue
*
               TEMP = ONE / MAX( ABS( S( JE, JE ) )*ASCALE,
     $                ABS( P( JE, JE ) )*BSCALE, SAFMIN )
               SALFAR = ( TEMP*S( JE, JE ) )*ASCALE
               SBETA = ( TEMP*P( JE, JE ) )*BSCALE
               ACOEF = SBETA*ASCALE
               BCOEFR = SALFAR*BSCALE
               BCOEFI = ZERO
*
*              Scale to avoid underflow
*
               SCALE = ONE
               LSA = ABS( SBETA ).GE.SAFMIN .AND. ABS( ACOEF ).LT.SMALL
               LSB = ABS( SALFAR ).GE.SAFMIN .AND. ABS( BCOEFR ).LT.
     $               SMALL
               IF( LSA )
     $            SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )
               IF( LSB )
     $            SCALE = MAX( SCALE, ( SMALL / ABS( SALFAR ) )*
     $                    MIN( BNORM, BIG ) )
               IF( LSA .OR. LSB ) THEN
                  SCALE = MIN( SCALE, ONE /
     $                    ( SAFMIN*MAX( ONE, ABS( ACOEF ),
     $                    ABS( BCOEFR ) ) ) )
                  IF( LSA ) THEN
                     ACOEF = ASCALE*( SCALE*SBETA )
                  ELSE
                     ACOEF = SCALE*ACOEF
                  END IF
                  IF( LSB ) THEN
                     BCOEFR = BSCALE*( SCALE*SALFAR )
                  ELSE
                     BCOEFR = SCALE*BCOEFR
                  END IF
               END IF
               ACOEFA = ABS( ACOEF )
               BCOEFA = ABS( BCOEFR )
*
*              First component is 1
*
               WORK( 2*N+JE ) = ONE
               XMAX = ONE
            ELSE
*
*              Complex eigenvalue
*
               CALL DLAG2( S( JE, JE ), LDS, P( JE, JE ), LDP,
     $                     SAFMIN*SAFETY, ACOEF, TEMP, BCOEFR, TEMP2,
     $                     BCOEFI )
               BCOEFI = -BCOEFI
               IF( BCOEFI.EQ.ZERO ) THEN
                  INFO = JE
                  RETURN
               END IF
*
*              Scale to avoid over/underflow
*
               ACOEFA = ABS( ACOEF )
               BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
               SCALE = ONE
               IF( ACOEFA*ULP.LT.SAFMIN .AND. ACOEFA.GE.SAFMIN )
     $            SCALE = ( SAFMIN / ULP ) / ACOEFA
               IF( BCOEFA*ULP.LT.SAFMIN .AND. BCOEFA.GE.SAFMIN )
     $            SCALE = MAX( SCALE, ( SAFMIN / ULP ) / BCOEFA )
               IF( SAFMIN*ACOEFA.GT.ASCALE )
     $            SCALE = ASCALE / ( SAFMIN*ACOEFA )
               IF( SAFMIN*BCOEFA.GT.BSCALE )
     $            SCALE = MIN( SCALE, BSCALE / ( SAFMIN*BCOEFA ) )
               IF( SCALE.NE.ONE ) THEN
                  ACOEF = SCALE*ACOEF
                  ACOEFA = ABS( ACOEF )
                  BCOEFR = SCALE*BCOEFR
                  BCOEFI = SCALE*BCOEFI
                  BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
               END IF
*
*              Compute first two components of eigenvector
*
               TEMP = ACOEF*S( JE+1, JE )
               TEMP2R = ACOEF*S( JE, JE ) - BCOEFR*P( JE, JE )
               TEMP2I = -BCOEFI*P( JE, JE )
               IF( ABS( TEMP ).GT.ABS( TEMP2R )+ABS( TEMP2I ) ) THEN
                  WORK( 2*N+JE ) = ONE
                  WORK( 3*N+JE ) = ZERO
                  WORK( 2*N+JE+1 ) = -TEMP2R / TEMP
                  WORK( 3*N+JE+1 ) = -TEMP2I / TEMP
               ELSE
                  WORK( 2*N+JE+1 ) = ONE
                  WORK( 3*N+JE+1 ) = ZERO
                  TEMP = ACOEF*S( JE, JE+1 )
                  WORK( 2*N+JE ) = ( BCOEFR*P( JE+1, JE+1 )-ACOEF*
     $                             S( JE+1, JE+1 ) ) / TEMP
                  WORK( 3*N+JE ) = BCOEFI*P( JE+1, JE+1 ) / TEMP
               END IF
               XMAX = MAX( ABS( WORK( 2*N+JE ) )+ABS( WORK( 3*N+JE ) ),
     $                ABS( WORK( 2*N+JE+1 ) )+ABS( WORK( 3*N+JE+1 ) ) )
            END IF
*
            DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )
*
*                                           T
*           Triangular solve of  (a A - b B)  y = 0
*
*                                   T
*           (rowwise in  (a A - b B) , or columnwise in (a A - b B) )
*
            IL2BY2 = .FALSE.
*
            DO 160 J = JE + NW, N
               IF( IL2BY2 ) THEN
                  IL2BY2 = .FALSE.
                  GO TO 160
               END IF
*
               NA = 1
               BDIAG( 1 ) = P( J, J )
               IF( J.LT.N ) THEN
                  IF( S( J+1, J ).NE.ZERO ) THEN
                     IL2BY2 = .TRUE.
                     BDIAG( 2 ) = P( J+1, J+1 )
                     NA = 2
                  END IF
               END IF
*
*              Check whether scaling is necessary for dot products
*
               XSCALE = ONE / MAX( ONE, XMAX )
               TEMP = MAX( WORK( J ), WORK( N+J ),
     $                ACOEFA*WORK( J )+BCOEFA*WORK( N+J ) )
               IF( IL2BY2 )
     $            TEMP = MAX( TEMP, WORK( J+1 ), WORK( N+J+1 ),
     $                   ACOEFA*WORK( J+1 )+BCOEFA*WORK( N+J+1 ) )
               IF( TEMP.GT.BIGNUM*XSCALE ) THEN
                  DO 90 JW = 0, NW - 1
                     DO 80 JR = JE, J - 1
                        WORK( ( JW+2 )*N+JR ) = XSCALE*
     $                     WORK( ( JW+2 )*N+JR )
   80                CONTINUE
   90             CONTINUE
                  XMAX = XMAX*XSCALE
               END IF
*
*              Compute dot products
*
*                    j-1
*              SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
*                    k=je
*
*              To reduce the op count, this is done as
*
*              _        j-1                  _        j-1
*              a*conjg( sum  S(k,j)*x(k) ) - b*conjg( sum  P(k,j)*x(k) )
*                       k=je                          k=je
*
*              which may cause underflow problems if A or B are close
*              to underflow.  (E.g., less than SMALL.)
*
*
               DO 120 JW = 1, NW
                  DO 110 JA = 1, NA
                     SUMS( JA, JW ) = ZERO
                     SUMP( JA, JW ) = ZERO
*
                     DO 100 JR = JE, J - 1
                        SUMS( JA, JW ) = SUMS( JA, JW ) +
     $                                   S( JR, J+JA-1 )*
     $                                   WORK( ( JW+1 )*N+JR )
                        SUMP( JA, JW ) = SUMP( JA, JW ) +
     $                                   P( JR, J+JA-1 )*
     $                                   WORK( ( JW+1 )*N+JR )
  100                CONTINUE
  110             CONTINUE
  120          CONTINUE
*
               DO 130 JA = 1, NA
                  IF( ILCPLX ) THEN
                     SUM( JA, 1 ) = -ACOEF*SUMS( JA, 1 ) +
     $                              BCOEFR*SUMP( JA, 1 ) -
     $                              BCOEFI*SUMP( JA, 2 )
                     SUM( JA, 2 ) = -ACOEF*SUMS( JA, 2 ) +
     $                              BCOEFR*SUMP( JA, 2 ) +
     $                              BCOEFI*SUMP( JA, 1 )
                  ELSE
                     SUM( JA, 1 ) = -ACOEF*SUMS( JA, 1 ) +
     $                              BCOEFR*SUMP( JA, 1 )
                  END IF
  130          CONTINUE
*
*                                  T
*              Solve  ( a A - b B )  y = SUM(,)
*              with scaling and perturbation of the denominator
*
               CALL DLALN2( .TRUE., NA, NW, DMIN, ACOEF, S( J, J ), LDS,
     $                      BDIAG( 1 ), BDIAG( 2 ), SUM, 2, BCOEFR,
     $                      BCOEFI, WORK( 2*N+J ), N, SCALE, TEMP,
     $                      IINFO )
               IF( SCALE.LT.ONE ) THEN
                  DO 150 JW = 0, NW - 1
                     DO 140 JR = JE, J - 1
                        WORK( ( JW+2 )*N+JR ) = SCALE*
     $                     WORK( ( JW+2 )*N+JR )
  140                CONTINUE
  150             CONTINUE
                  XMAX = SCALE*XMAX
               END IF
               XMAX = MAX( XMAX, TEMP )
  160       CONTINUE
*
*           Copy eigenvector to VL, back transforming if
*           HOWMNY='B'.
*
            IEIG = IEIG + 1
            IF( ILBACK ) THEN
               DO 170 JW = 0, NW - 1
                  CALL DGEMV( 'N', N, N+1-JE, ONE, VL( 1, JE ), LDVL,
     $                        WORK( ( JW+2 )*N+JE ), 1, ZERO,
     $                        WORK( ( JW+4 )*N+1 ), 1 )
  170          CONTINUE
               CALL DLACPY( ' ', N, NW, WORK( 4*N+1 ), N, VL( 1, JE ),
     $                      LDVL )
               IBEG = 1
            ELSE
               CALL DLACPY( ' ', N, NW, WORK( 2*N+1 ), N, VL( 1, IEIG ),
     $                      LDVL )
               IBEG = JE
            END IF
*
*           Scale eigenvector
*
            XMAX = ZERO
            IF( ILCPLX ) THEN
               DO 180 J = IBEG, N
                  XMAX = MAX( XMAX, ABS( VL( J, IEIG ) )+
     $                   ABS( VL( J, IEIG+1 ) ) )
  180          CONTINUE
            ELSE
               DO 190 J = IBEG, N
                  XMAX = MAX( XMAX, ABS( VL( J, IEIG ) ) )
  190          CONTINUE
            END IF
*
            IF( XMAX.GT.SAFMIN ) THEN
               XSCALE = ONE / XMAX
*
               DO 210 JW = 0, NW - 1
                  DO 200 JR = IBEG, N
                     VL( JR, IEIG+JW ) = XSCALE*VL( JR, IEIG+JW )
  200             CONTINUE
  210          CONTINUE
            END IF
            IEIG = IEIG + NW - 1
*
  220    CONTINUE
      END IF
*
*     Right eigenvectors
*
      IF( COMPR ) THEN
         IEIG = IM + 1
*
*        Main loop over eigenvalues
*
         ILCPLX = .FALSE.
         DO 500 JE = N, 1, -1
*
*           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or
*           (b) this would be the second of a complex pair.
*           Check for complex eigenvalue, so as to be sure of which
*           entry(-ies) of SELECT to look at -- if complex, SELECT(JE)
*           or SELECT(JE-1).
*           If this is a complex pair, the 2-by-2 diagonal block
*           corresponding to the eigenvalue is in rows/columns JE-1:JE
*
            IF( ILCPLX ) THEN
               ILCPLX = .FALSE.
               GO TO 500
            END IF
            NW = 1
            IF( JE.GT.1 ) THEN
               IF( S( JE, JE-1 ).NE.ZERO ) THEN
                  ILCPLX = .TRUE.
                  NW = 2
               END IF
            END IF
            IF( ILALL ) THEN
               ILCOMP = .TRUE.
            ELSE IF( ILCPLX ) THEN
               ILCOMP = SELECT( JE ) .OR. SELECT( JE-1 )
            ELSE
               ILCOMP = SELECT( JE )
            END IF
            IF( .NOT.ILCOMP )
     $         GO TO 500
*
*           Decide if (a) singular pencil, (b) real eigenvalue, or
*           (c) complex eigenvalue.
*
            IF( .NOT.ILCPLX ) THEN
               IF( ABS( S( JE, JE ) ).LE.SAFMIN .AND.
     $             ABS( P( JE, JE ) ).LE.SAFMIN ) THEN
*
*                 Singular matrix pencil -- unit eigenvector
*
                  IEIG = IEIG - 1
                  DO 230 JR = 1, N
                     VR( JR, IEIG ) = ZERO
  230             CONTINUE
                  VR( IEIG, IEIG ) = ONE
                  GO TO 500
               END IF
            END IF
*
*           Clear vector
*
            DO 250 JW = 0, NW - 1
               DO 240 JR = 1, N
                  WORK( ( JW+2 )*N+JR ) = ZERO
  240          CONTINUE
  250       CONTINUE
*
*           Compute coefficients in  ( a A - b B ) x = 0
*              a  is  ACOEF
*              b  is  BCOEFR + i*BCOEFI
*
            IF( .NOT.ILCPLX ) THEN
*
*              Real eigenvalue
*
               TEMP = ONE / MAX( ABS( S( JE, JE ) )*ASCALE,
     $                ABS( P( JE, JE ) )*BSCALE, SAFMIN )
               SALFAR = ( TEMP*S( JE, JE ) )*ASCALE
               SBETA = ( TEMP*P( JE, JE ) )*BSCALE
               ACOEF = SBETA*ASCALE
               BCOEFR = SALFAR*BSCALE
               BCOEFI = ZERO
*
*              Scale to avoid underflow
*
               SCALE = ONE
               LSA = ABS( SBETA ).GE.SAFMIN .AND. ABS( ACOEF ).LT.SMALL
               LSB = ABS( SALFAR ).GE.SAFMIN .AND. ABS( BCOEFR ).LT.
     $               SMALL
               IF( LSA )
     $            SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )
               IF( LSB )
     $            SCALE = MAX( SCALE, ( SMALL / ABS( SALFAR ) )*
     $                    MIN( BNORM, BIG ) )
               IF( LSA .OR. LSB ) THEN
                  SCALE = MIN( SCALE, ONE /
     $                    ( SAFMIN*MAX( ONE, ABS( ACOEF ),
     $                    ABS( BCOEFR ) ) ) )
                  IF( LSA ) THEN
                     ACOEF = ASCALE*( SCALE*SBETA )
                  ELSE
                     ACOEF = SCALE*ACOEF
                  END IF
                  IF( LSB ) THEN
                     BCOEFR = BSCALE*( SCALE*SALFAR )
                  ELSE
                     BCOEFR = SCALE*BCOEFR
                  END IF
               END IF
               ACOEFA = ABS( ACOEF )
               BCOEFA = ABS( BCOEFR )
*
*              First component is 1
*
               WORK( 2*N+JE ) = ONE
               XMAX = ONE
*
*              Compute contribution from column JE of A and B to sum
*              (See "Further Details", above.)
*
               DO 260 JR = 1, JE - 1
                  WORK( 2*N+JR ) = BCOEFR*P( JR, JE ) -
     $                             ACOEF*S( JR, JE )
  260          CONTINUE
            ELSE
*
*              Complex eigenvalue
*
               CALL DLAG2( S( JE-1, JE-1 ), LDS, P( JE-1, JE-1 ), LDP,
     $                     SAFMIN*SAFETY, ACOEF, TEMP, BCOEFR, TEMP2,
     $                     BCOEFI )
               IF( BCOEFI.EQ.ZERO ) THEN
                  INFO = JE - 1
                  RETURN
               END IF
*
*              Scale to avoid over/underflow
*
               ACOEFA = ABS( ACOEF )
               BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
               SCALE = ONE
               IF( ACOEFA*ULP.LT.SAFMIN .AND. ACOEFA.GE.SAFMIN )
     $            SCALE = ( SAFMIN / ULP ) / ACOEFA
               IF( BCOEFA*ULP.LT.SAFMIN .AND. BCOEFA.GE.SAFMIN )
     $            SCALE = MAX( SCALE, ( SAFMIN / ULP ) / BCOEFA )
               IF( SAFMIN*ACOEFA.GT.ASCALE )
     $            SCALE = ASCALE / ( SAFMIN*ACOEFA )
               IF( SAFMIN*BCOEFA.GT.BSCALE )
     $            SCALE = MIN( SCALE, BSCALE / ( SAFMIN*BCOEFA ) )
               IF( SCALE.NE.ONE ) THEN
                  ACOEF = SCALE*ACOEF
                  ACOEFA = ABS( ACOEF )
                  BCOEFR = SCALE*BCOEFR
                  BCOEFI = SCALE*BCOEFI
                  BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
               END IF
*
*              Compute first two components of eigenvector
*              and contribution to sums
*
               TEMP = ACOEF*S( JE, JE-1 )
               TEMP2R = ACOEF*S( JE, JE ) - BCOEFR*P( JE, JE )
               TEMP2I = -BCOEFI*P( JE, JE )
               IF( ABS( TEMP ).GE.ABS( TEMP2R )+ABS( TEMP2I ) ) THEN
                  WORK( 2*N+JE ) = ONE
                  WORK( 3*N+JE ) = ZERO
                  WORK( 2*N+JE-1 ) = -TEMP2R / TEMP
                  WORK( 3*N+JE-1 ) = -TEMP2I / TEMP
               ELSE
                  WORK( 2*N+JE-1 ) = ONE
                  WORK( 3*N+JE-1 ) = ZERO
                  TEMP = ACOEF*S( JE-1, JE )
                  WORK( 2*N+JE ) = ( BCOEFR*P( JE-1, JE-1 )-ACOEF*
     $                             S( JE-1, JE-1 ) ) / TEMP
                  WORK( 3*N+JE ) = BCOEFI*P( JE-1, JE-1 ) / TEMP
               END IF
*
               XMAX = MAX( ABS( WORK( 2*N+JE ) )+ABS( WORK( 3*N+JE ) ),
     $                ABS( WORK( 2*N+JE-1 ) )+ABS( WORK( 3*N+JE-1 ) ) )
*
*              Compute contribution from columns JE and JE-1
*              of A and B to the sums.
*
               CREALA = ACOEF*WORK( 2*N+JE-1 )
               CIMAGA = ACOEF*WORK( 3*N+JE-1 )
               CREALB = BCOEFR*WORK( 2*N+JE-1 ) -
     $                  BCOEFI*WORK( 3*N+JE-1 )
               CIMAGB = BCOEFI*WORK( 2*N+JE-1 ) +
     $                  BCOEFR*WORK( 3*N+JE-1 )
               CRE2A = ACOEF*WORK( 2*N+JE )
               CIM2A = ACOEF*WORK( 3*N+JE )
               CRE2B = BCOEFR*WORK( 2*N+JE ) - BCOEFI*WORK( 3*N+JE )
               CIM2B = BCOEFI*WORK( 2*N+JE ) + BCOEFR*WORK( 3*N+JE )
               DO 270 JR = 1, JE - 2
                  WORK( 2*N+JR ) = -CREALA*S( JR, JE-1 ) +
     $                             CREALB*P( JR, JE-1 ) -
     $                             CRE2A*S( JR, JE ) + CRE2B*P( JR, JE )
                  WORK( 3*N+JR ) = -CIMAGA*S( JR, JE-1 ) +
     $                             CIMAGB*P( JR, JE-1 ) -
     $                             CIM2A*S( JR, JE ) + CIM2B*P( JR, JE )
  270          CONTINUE
            END IF
*
            DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )
*
*           Columnwise triangular solve of  (a A - b B)  x = 0
*
            IL2BY2 = .FALSE.
            DO 370 J = JE - NW, 1, -1
*
*              If a 2-by-2 block, is in position j-1:j, wait until
*              next iteration to process it (when it will be j:j+1)
*
               IF( .NOT.IL2BY2 .AND. J.GT.1 ) THEN
                  IF( S( J, J-1 ).NE.ZERO ) THEN
                     IL2BY2 = .TRUE.
                     GO TO 370
                  END IF
               END IF
               BDIAG( 1 ) = P( J, J )
               IF( IL2BY2 ) THEN
                  NA = 2
                  BDIAG( 2 ) = P( J+1, J+1 )
               ELSE
                  NA = 1
               END IF
*
*              Compute x(j) (and x(j+1), if 2-by-2 block)
*
               CALL DLALN2( .FALSE., NA, NW, DMIN, ACOEF, S( J, J ),
     $                      LDS, BDIAG( 1 ), BDIAG( 2 ), WORK( 2*N+J ),
     $                      N, BCOEFR, BCOEFI, SUM, 2, SCALE, TEMP,
     $                      IINFO )
               IF( SCALE.LT.ONE ) THEN
*
                  DO 290 JW = 0, NW - 1
                     DO 280 JR = 1, JE
                        WORK( ( JW+2 )*N+JR ) = SCALE*
     $                     WORK( ( JW+2 )*N+JR )
  280                CONTINUE
  290             CONTINUE
               END IF
               XMAX = MAX( SCALE*XMAX, TEMP )
*
               DO 310 JW = 1, NW
                  DO 300 JA = 1, NA
                     WORK( ( JW+1 )*N+J+JA-1 ) = SUM( JA, JW )
  300             CONTINUE
  310          CONTINUE
*
*              w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling
*
               IF( J.GT.1 ) THEN
*
*                 Check whether scaling is necessary for sum.
*
                  XSCALE = ONE / MAX( ONE, XMAX )
                  TEMP = ACOEFA*WORK( J ) + BCOEFA*WORK( N+J )
                  IF( IL2BY2 )
     $               TEMP = MAX( TEMP, ACOEFA*WORK( J+1 )+BCOEFA*
     $                      WORK( N+J+1 ) )
                  TEMP = MAX( TEMP, ACOEFA, BCOEFA )
                  IF( TEMP.GT.BIGNUM*XSCALE ) THEN
*
                     DO 330 JW = 0, NW - 1
                        DO 320 JR = 1, JE
                           WORK( ( JW+2 )*N+JR ) = XSCALE*
     $                        WORK( ( JW+2 )*N+JR )
  320                   CONTINUE
  330                CONTINUE
                     XMAX = XMAX*XSCALE
                  END IF
*
*                 Compute the contributions of the off-diagonals of
*                 column j (and j+1, if 2-by-2 block) of A and B to the
*                 sums.
*
*
                  DO 360 JA = 1, NA
                     IF( ILCPLX ) THEN
                        CREALA = ACOEF*WORK( 2*N+J+JA-1 )
                        CIMAGA = ACOEF*WORK( 3*N+J+JA-1 )
                        CREALB = BCOEFR*WORK( 2*N+J+JA-1 ) -
     $                           BCOEFI*WORK( 3*N+J+JA-1 )
                        CIMAGB = BCOEFI*WORK( 2*N+J+JA-1 ) +
     $                           BCOEFR*WORK( 3*N+J+JA-1 )
                        DO 340 JR = 1, J - 1
                           WORK( 2*N+JR ) = WORK( 2*N+JR ) -
     $                                      CREALA*S( JR, J+JA-1 ) +
     $                                      CREALB*P( JR, J+JA-1 )
                           WORK( 3*N+JR ) = WORK( 3*N+JR ) -
     $                                      CIMAGA*S( JR, J+JA-1 ) +
     $                                      CIMAGB*P( JR, J+JA-1 )
  340                   CONTINUE
                     ELSE
                        CREALA = ACOEF*WORK( 2*N+J+JA-1 )
                        CREALB = BCOEFR*WORK( 2*N+J+JA-1 )
                        DO 350 JR = 1, J - 1
                           WORK( 2*N+JR ) = WORK( 2*N+JR ) -
     $                                      CREALA*S( JR, J+JA-1 ) +
     $                                      CREALB*P( JR, J+JA-1 )
  350                   CONTINUE
                     END IF
  360             CONTINUE
               END IF
*
               IL2BY2 = .FALSE.
  370       CONTINUE
*
*           Copy eigenvector to VR, back transforming if
*           HOWMNY='B'.
*
            IEIG = IEIG - NW
            IF( ILBACK ) THEN
*
               DO 410 JW = 0, NW - 1
                  DO 380 JR = 1, N
                     WORK( ( JW+4 )*N+JR ) = WORK( ( JW+2 )*N+1 )*
     $                                       VR( JR, 1 )
  380             CONTINUE
*
*                 A series of compiler directives to defeat
*                 vectorization for the next loop
*
*
                  DO 400 JC = 2, JE
                     DO 390 JR = 1, N
                        WORK( ( JW+4 )*N+JR ) = WORK( ( JW+4 )*N+JR ) +
     $                     WORK( ( JW+2 )*N+JC )*VR( JR, JC )
  390                CONTINUE
  400             CONTINUE
  410          CONTINUE
*
               DO 430 JW = 0, NW - 1
                  DO 420 JR = 1, N
                     VR( JR, IEIG+JW ) = WORK( ( JW+4 )*N+JR )
  420             CONTINUE
  430          CONTINUE
*
               IEND = N
            ELSE
               DO 450 JW = 0, NW - 1
                  DO 440 JR = 1, N
                     VR( JR, IEIG+JW ) = WORK( ( JW+2 )*N+JR )
  440             CONTINUE
  450          CONTINUE
*
               IEND = JE
            END IF
*
*           Scale eigenvector
*
            XMAX = ZERO
            IF( ILCPLX ) THEN
               DO 460 J = 1, IEND
                  XMAX = MAX( XMAX, ABS( VR( J, IEIG ) )+
     $                   ABS( VR( J, IEIG+1 ) ) )
  460          CONTINUE
            ELSE
               DO 470 J = 1, IEND
                  XMAX = MAX( XMAX, ABS( VR( J, IEIG ) ) )
  470          CONTINUE
            END IF
*
            IF( XMAX.GT.SAFMIN ) THEN
               XSCALE = ONE / XMAX
               DO 490 JW = 0, NW - 1
                  DO 480 JR = 1, IEND
                     VR( JR, IEIG+JW ) = XSCALE*VR( JR, IEIG+JW )
  480             CONTINUE
  490          CONTINUE
            END IF
  500    CONTINUE
      END IF
*
      RETURN
*
*     End of DTGEVC
*
      END
*> \brief \b DTRTI2 computes the inverse of a triangular matrix (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DTRTI2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrti2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrti2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrti2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          DIAG, UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTRTI2 computes the inverse of a real upper or lower triangular
*> matrix.
*>
*> This is the Level 2 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the matrix A is upper or lower triangular.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>          Specifies whether or not the matrix A is unit triangular.
*>          = 'N':  Non-unit triangular
*>          = 'U':  Unit triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the triangular matrix A.  If UPLO = 'U', the
*>          leading n by n upper triangular part of the array A contains
*>          the upper triangular matrix, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading n by n lower triangular part of the array A contains
*>          the lower triangular matrix, and the strictly upper
*>          triangular part of A is not referenced.  If DIAG = 'U', the
*>          diagonal elements of A are also not referenced and are
*>          assumed to be 1.
*>
*>          On exit, the (triangular) inverse of the original matrix, in
*>          the same storage format.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -k, the k-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup doubleOTHERcomputational
*
*  =====================================================================
      SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DTRMV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTI2', -INFO )
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Compute inverse of upper triangular matrix.
*
         DO 10 J = 1, N
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
*
*           Compute elements 1:j-1 of j-th column.
*
            CALL DTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA,
     $                  A( 1, J ), 1 )
            CALL DSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE
      ELSE
*
*        Compute inverse of lower triangular matrix.
*
         DO 20 J = N, 1, -1
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
            IF( J.LT.N ) THEN
*
*              Compute elements j+1:n of j-th column.
*
               CALL DTRMV( 'Lower', 'No transpose', DIAG, N-J,
     $                     A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL DSCAL( N-J, AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of DTRTI2
*
      END
*> \brief \b DTRTRI
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download DTRTRI + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrtri.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrtri.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrtri.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          DIAG, UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTRTRI computes the inverse of a real upper or lower triangular
*> matrix A.
*>
*> This is the Level 3 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  A is upper triangular;
*>          = 'L':  A is lower triangular.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>          = 'N':  A is non-unit triangular;
*>          = 'U':  A is unit triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the triangular matrix A.  If UPLO = 'U', the
*>          leading N-by-N upper triangular part of the array A contains
*>          the upper triangular matrix, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading N-by-N lower triangular part of the array A contains
*>          the lower triangular matrix, and the strictly upper
*>          triangular part of A is not referenced.  If DIAG = 'U', the
*>          diagonal elements of A are also not referenced and are
*>          assumed to be 1.
*>          On exit, the (triangular) inverse of the original matrix, in
*>          the same storage format.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
*>               matrix is singular and its inverse can not be computed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup doubleOTHERcomputational
*
*  =====================================================================
      SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J, JB, NB, NN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTRMM, DTRSM, DTRTI2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check for singularity if non-unit.
*
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
         INFO = 0
      END IF
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DTRTRI', UPLO // DIAG, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code
*
         IF( UPPER ) THEN
*
*           Compute inverse of upper triangular matrix
*
            DO 20 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
*
*              Compute rows 1:j-1 of current block column
*
               CALL DTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, ONE, A, LDA, A( 1, J ), LDA )
               CALL DTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
*
*              Compute inverse of current diagonal block
*
               CALL DTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
   20       CONTINUE
         ELSE
*
*           Compute inverse of lower triangular matrix
*
            NN = ( ( N-1 ) / NB )*NB + 1
            DO 30 J = NN, 1, -NB
               JB = MIN( NB, N-J+1 )
               IF( J+JB.LE.N ) THEN
*
*                 Compute rows j+jb:n of current block column
*
                  CALL DTRMM( 'Left', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,
     $                        A( J+JB, J ), LDA )
                  CALL DTRSM( 'Right', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, -ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
*
*              Compute inverse of current diagonal block
*
               CALL DTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
   30       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DTRTRI
*
      END
*> \brief \b IEEECK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download IEEECK + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ieeeck.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ieeeck.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ieeeck.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
* 
*       .. Scalar Arguments ..
*       INTEGER            ISPEC
*       REAL               ONE, ZERO
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> IEEECK is called from the ILAENV to verify that Infinity and
*> possibly NaN arithmetic is safe (i.e. will not trap).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ISPEC
*> \verbatim
*>          ISPEC is INTEGER
*>          Specifies whether to test just for inifinity arithmetic
*>          or whether to test for infinity and NaN arithmetic.
*>          = 0: Verify infinity arithmetic only.
*>          = 1: Verify infinity and NaN arithmetic.
*> \endverbatim
*>
*> \param[in] ZERO
*> \verbatim
*>          ZERO is REAL
*>          Must contain the value 0.0
*>          This is passed to prevent the compiler from optimizing
*>          away this code.
*> \endverbatim
*>
*> \param[in] ONE
*> \verbatim
*>          ONE is REAL
*>          Must contain the value 1.0
*>          This is passed to prevent the compiler from optimizing
*>          away this code.
*>
*>  RETURN VALUE:  INTEGER
*>          = 0:  Arithmetic failed to produce the correct answers
*>          = 1:  Arithmetic produced the correct answers
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
*
*  -- LAPACK auxiliary routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
*     ..
*     .. Executable Statements ..
      IEEECK = 1
*
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
*
*
*
*     Return if we were only asked to check infinity arithmetic
*
      IF( ISPEC.EQ.0 )
     $   RETURN
*
      NAN1 = POSINF + NEGINF
*
      NAN2 = POSINF / NEGINF
*
      NAN3 = POSINF / POSINF
*
      NAN4 = POSINF*ZERO
*
      NAN5 = NEGINF*NEGZRO
*
      NAN6 = NAN5*ZERO
*
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      RETURN
      END
*> \brief \b ILADLC scans a matrix for its last non-zero column.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ILADLC + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladlc.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iladlc.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iladlc.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION ILADLC( M, N, A, LDA )
* 
*       .. Scalar Arguments ..
*       INTEGER            M, N, LDA
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ILADLC scans A for its last non-zero column.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          The m by n matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      INTEGER FUNCTION ILADLC( M, N, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER I
*     ..
*     .. Executable Statements ..
*
*     Quick test for the common case where one corner is non-zero.
      IF( N.EQ.0 ) THEN
         ILADLC = N
      ELSE IF( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILADLC = N
      ELSE
*     Now scan each column from the end, returning with the first non-zero.
         DO ILADLC = N, 1, -1
            DO I = 1, M
               IF( A(I, ILADLC).NE.ZERO ) RETURN
            END DO
         END DO
      END IF
      RETURN
      END
*> \brief \b ILADLR scans a matrix for its last non-zero row.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ILADLR + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladlr.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iladlr.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iladlr.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION ILADLR( M, N, A, LDA )
* 
*       .. Scalar Arguments ..
*       INTEGER            M, N, LDA
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ILADLR scans A for its last non-zero row.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          The m by n matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      INTEGER FUNCTION ILADLR( M, N, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER I, J
*     ..
*     .. Executable Statements ..
*
*     Quick test for the common case where one corner is non-zero.
      IF( M.EQ.0 ) THEN
         ILADLR = M
      ELSE IF( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILADLR = M
      ELSE
*     Scan up each column tracking the last zero row seen.
         ILADLR = 0
         DO J = 1, N
            I=M
            DO WHILE((A(MAX(I,1),J).EQ.ZERO).AND.(I.GE.1))
               I=I-1
            ENDDO
            ILADLR = MAX( ILADLR, I )
         END DO
      END IF
      RETURN
      END
*> \brief \b ILAENV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ILAENV + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaenv.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaenv.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaenv.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
* 
*       .. Scalar Arguments ..
*       CHARACTER*( * )    NAME, OPTS
*       INTEGER            ISPEC, N1, N2, N3, N4
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ILAENV is called from the LAPACK routines to choose problem-dependent
*> parameters for the local environment.  See ISPEC for a description of
*> the parameters.
*>
*> ILAENV returns an INTEGER
*> if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
*> if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*>
*> This version provides a set of parameters which should give good,
*> but not optimal, performance on many of the currently available
*> computers.  Users are encouraged to modify this subroutine to set
*> the tuning parameters for their particular machine using the option
*> and problem size information in the arguments.
*>
*> This routine will not function correctly if it is converted to all
*> lower case.  Converting it to all upper case is allowed.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ISPEC
*> \verbatim
*>          ISPEC is INTEGER
*>          Specifies the parameter to be returned as the value of
*>          ILAENV.
*>          = 1: the optimal blocksize; if this value is 1, an unblocked
*>               algorithm will give the best performance.
*>          = 2: the minimum block size for which the block routine
*>               should be used; if the usable block size is less than
*>               this value, an unblocked routine should be used.
*>          = 3: the crossover point (in a block routine, for N less
*>               than this value, an unblocked routine should be used)
*>          = 4: the number of shifts, used in the nonsymmetric
*>               eigenvalue routines (DEPRECATED)
*>          = 5: the minimum column dimension for blocking to be used;
*>               rectangular blocks must have dimension at least k by m,
*>               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*>          = 6: the crossover point for the SVD (when reducing an m by n
*>               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*>               this value, a QR factorization is used first to reduce
*>               the matrix to a triangular form.)
*>          = 7: the number of processors
*>          = 8: the crossover point for the multishift QR method
*>               for nonsymmetric eigenvalue problems (DEPRECATED)
*>          = 9: maximum size of the subproblems at the bottom of the
*>               computation tree in the divide-and-conquer algorithm
*>               (used by xGELSD and xGESDD)
*>          =10: ieee NaN arithmetic can be trusted not to trap
*>          =11: infinity arithmetic can be trusted not to trap
*>          12 <= ISPEC <= 16:
*>               xHSEQR or related subroutines,
*>               see IPARMQ for detailed explanation
*> \endverbatim
*>
*> \param[in] NAME
*> \verbatim
*>          NAME is CHARACTER*(*)
*>          The name of the calling subroutine, in either upper case or
*>          lower case.
*> \endverbatim
*>
*> \param[in] OPTS
*> \verbatim
*>          OPTS is CHARACTER*(*)
*>          The character options to the subroutine NAME, concatenated
*>          into a single character string.  For example, UPLO = 'U',
*>          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*>          be specified as OPTS = 'UTN'.
*> \endverbatim
*>
*> \param[in] N1
*> \verbatim
*>          N1 is INTEGER
*> \endverbatim
*>
*> \param[in] N2
*> \verbatim
*>          N2 is INTEGER
*> \endverbatim
*>
*> \param[in] N3
*> \verbatim
*>          N3 is INTEGER
*> \endverbatim
*>
*> \param[in] N4
*> \verbatim
*>          N4 is INTEGER
*>          Problem dimensions for the subroutine NAME; these may not all
*>          be required.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup auxOTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The following conventions have been used when calling ILAENV from the
*>  LAPACK routines:
*>  1)  OPTS is a concatenation of all of the character options to
*>      subroutine NAME, in the same order that they appear in the
*>      argument list for NAME, even if they are not used in determining
*>      the value of the parameter specified by ISPEC.
*>  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*>      that they appear in the argument list for NAME.  N1 is used
*>      first, N2 second, and so on, and unused problem dimensions are
*>      passed a value of -1.
*>  3)  The parameter value returned by ILAENV is checked for validity in
*>      the calling subroutine.  For example, ILAENV is used to retrieve
*>      the optimal blocksize for STRTRI as follows:
*>
*>      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*>      IF( NB.LE.1 ) NB = MAX( 1, N )
*> \endverbatim
*>
*  =====================================================================
      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*
*  -- LAPACK auxiliary routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. External Functions ..
      INTEGER            IEEECK, IPARMQ
      EXTERNAL           IEEECK, IPARMQ
*     ..
*     .. Executable Statements ..
*
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,
     $        130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
   10 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $             I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
*
      GO TO ( 50, 60, 70 )ISPEC
*
   50 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      ELSE IF( C2.EQ.'GG' ) THEN
         NB = 32
         IF( C3.EQ.'HD3' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         END IF
      END IF
      ILAENV = NB
      RETURN
*
   60 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'GG' ) THEN
         NBMIN = 2
         IF( C3.EQ.'HD3' ) THEN
            NBMIN = 2
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
   70 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'GG' ) THEN
         NX = 128
         IF( C3.EQ.'HD3' ) THEN
            NX = 128
         END IF
      END IF
      ILAENV = NX
      RETURN
*
   80 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
   90 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  100 CONTINUE
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  110 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  120 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
  130 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ILAENV = 25
      RETURN
*
  140 CONTINUE
*
*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
*
  150 CONTINUE
*
*     ISPEC = 11: infinity arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
*
  160 CONTINUE
*
*     12 <= ISPEC <= 16: xHSEQR or related subroutines.
*
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
*
*     End of ILAENV
*
      END
*> \brief \b IPARMQ
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download IPARMQ + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iparmq.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iparmq.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iparmq.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
* 
*       .. Scalar Arguments ..
*       INTEGER            IHI, ILO, ISPEC, LWORK, N
*       CHARACTER          NAME*( * ), OPTS*( * )
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>      This program sets problem and machine dependent parameters
*>      useful for xHSEQR and related subroutines for eigenvalue
*>      problems. It is called whenever
*>      IPARMQ is called with 12 <= ISPEC <= 16
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ISPEC
*> \verbatim
*>          ISPEC is integer scalar
*>              ISPEC specifies which tunable parameter IPARMQ should
*>              return.
*>
*>              ISPEC=12: (INMIN)  Matrices of order nmin or less
*>                        are sent directly to xLAHQR, the implicit
*>                        double shift QR algorithm.  NMIN must be
*>                        at least 11.
*>
*>              ISPEC=13: (INWIN)  Size of the deflation window.
*>                        This is best set greater than or equal to
*>                        the number of simultaneous shifts NS.
*>                        Larger matrices benefit from larger deflation
*>                        windows.
*>
*>              ISPEC=14: (INIBL) Determines when to stop nibbling and
*>                        invest in an (expensive) multi-shift QR sweep.
*>                        If the aggressive early deflation subroutine
*>                        finds LD converged eigenvalues from an order
*>                        NW deflation window and LD.GT.(NW*NIBBLE)/100,
*>                        then the next QR sweep is skipped and early
*>                        deflation is applied immediately to the
*>                        remaining active diagonal block.  Setting
*>                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
*>                        multi-shift QR sweep whenever early deflation
*>                        finds a converged eigenvalue.  Setting
*>                        IPARMQ(ISPEC=14) greater than or equal to 100
*>                        prevents TTQRE from skipping a multi-shift
*>                        QR sweep.
*>
*>              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
*>                        a multi-shift QR iteration.
*>
*>              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
*>                        following meanings.
*>                        0:  During the multi-shift QR/QZ sweep,
*>                            blocked eigenvalue reordering, blocked
*>                            Hessenberg-triangular reduction,
*>                            reflections and/or rotations are not
*>                            accumulated when updating the
*>                            far-from-diagonal matrix entries.
*>                        1:  During the multi-shift QR/QZ sweep,
*>                            blocked eigenvalue reordering, blocked
*>                            Hessenberg-triangular reduction,
*>                            reflections and/or rotations are
*>                            accumulated, and matrix-matrix
*>                            multiplication is used to update the
*>                            far-from-diagonal matrix entries.
*>                        2:  During the multi-shift QR/QZ sweep,
*>                            blocked eigenvalue reordering, blocked
*>                            Hessenberg-triangular reduction,
*>                            reflections and/or rotations are
*>                            accumulated, and 2-by-2 block structure
*>                            is exploited during matrix-matrix
*>                            multiplies.
*>                        (If xTRMM is slower than xGEMM, then
*>                        IPARMQ(ISPEC=16)=1 may be more efficient than
*>                        IPARMQ(ISPEC=16)=2 despite the greater level of
*>                        arithmetic work implied by the latter choice.)
*> \endverbatim
*>
*> \param[in] NAME
*> \verbatim
*>          NAME is character string
*>               Name of the calling subroutine
*> \endverbatim
*>
*> \param[in] OPTS
*> \verbatim
*>          OPTS is character string
*>               This is a concatenation of the string arguments to
*>               TTQRE.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is integer scalar
*>               N is the order of the Hessenberg matrix H.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>               It is assumed that H is already upper triangular
*>               in rows and columns 1:ILO-1 and IHI+1:N.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is integer scalar
*>               The amount of workspace available.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup auxOTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>       Little is known about how best to choose these parameters.
*>       It is possible to use different values of the parameters
*>       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
*>
*>       It is probably best to choose different parameters for
*>       different matrices and different parameters at different
*>       times during the iteration, but this has not been
*>       implemented --- yet.
*>
*>
*>       The best choices of most of the parameters depend
*>       in an ill-understood way on the relative execution
*>       rate of xLAQR3 and xLAQR5 and on the nature of each
*>       particular eigenvalue problem.  Experiment may be the
*>       only practical way to determine which choices are most
*>       effective.
*>
*>       Following is a list of default values supplied by IPARMQ.
*>       These defaults may be adjusted in order to attain better
*>       performance in any particular computational environment.
*>
*>       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
*>                        Default: 75. (Must be at least 11.)
*>
*>       IPARMQ(ISPEC=13) Recommended deflation window size.
*>                        This depends on ILO, IHI and NS, the
*>                        number of simultaneous shifts returned
*>                        by IPARMQ(ISPEC=15).  The default for
*>                        (IHI-ILO+1).LE.500 is NS.  The default
*>                        for (IHI-ILO+1).GT.500 is 3*NS/2.
*>
*>       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
*>
*>       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
*>                        a multi-shift QR iteration.
*>
*>                        If IHI-ILO+1 is ...
*>
*>                        greater than      ...but less    ... the
*>                        or equal to ...      than        default is
*>
*>                                0               30       NS =   2+
*>                               30               60       NS =   4+
*>                               60              150       NS =  10
*>                              150              590       NS =  **
*>                              590             3000       NS =  64
*>                             3000             6000       NS = 128
*>                             6000             infinity   NS = 256
*>
*>                    (+)  By default matrices of this order are
*>                         passed to the implicit double shift routine
*>                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
*>                         values of NS are used only in case of a rare
*>                         xLAHQR failure.
*>
*>                    (**) The asterisks (**) indicate an ad-hoc
*>                         function increasing from 10 to 64.
*>
*>       IPARMQ(ISPEC=16) Select structured matrix multiply.
*>                        (See ISPEC=16 above for details.)
*>                        Default: 3.
*> \endverbatim
*>
*  =====================================================================
      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
*
*  -- LAPACK auxiliary routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
*
*  ================================================================
*     .. Parameters ..
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14,
     $                   ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14,
     $                   NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
*     ..
*     .. Local Scalars ..
      INTEGER            NH, NS
      INTEGER            I, IC, IZ
      CHARACTER          SUBNAM*6
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
*     ..
*     .. Executable Statements ..
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR.
     $    ( ISPEC.EQ.IACC22 ) ) THEN
*
*        ==== Set the number simultaneous shifts ====
*
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 )
     $      NS = 4
         IF( NH.GE.60 )
     $      NS = 10
         IF( NH.GE.150 )
     $      NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
         IF( NH.GE.590 )
     $      NS = 64
         IF( NH.GE.3000 )
     $      NS = 128
         IF( NH.GE.6000 )
     $      NS = 256
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
*
      IF( ISPEC.EQ.INMIN ) THEN
*
*
*        ===== Matrices of order smaller than NMIN get sent
*        .     to xLAHQR, the classic double shift algorithm.
*        .     This must be at least 11. ====
*
         IPARMQ = NMIN
*
      ELSE IF( ISPEC.EQ.INIBL ) THEN
*
*        ==== INIBL: skip a multi-shift qr iteration and
*        .    whenever aggressive early deflation finds
*        .    at least (NIBBLE*(window size)/100) deflations. ====
*
         IPARMQ = NIBBLE
*
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
*
*        ==== NSHFTS: The number of simultaneous shifts =====
*
         IPARMQ = NS
*
      ELSE IF( ISPEC.EQ.INWIN ) THEN
*
*        ==== NW: deflation window size.  ====
*
         IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
         ELSE
            IPARMQ = 3*NS / 2
         END IF
*
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
*
*        ==== IACC22: Whether to accumulate reflections
*        .     before updating the far-from-diagonal elements
*        .     and whether to use 2-by-2 block structure while
*        .     doing it.  A small amount of work could be saved
*        .     by making this choice dependent also upon the
*        .     NH=IHI-ILO+1.
*
*
*        Convert NAME to upper case if the first character is lower case.
*
         IPARMQ = 0
         SUBNAM = NAME
         IC = ICHAR( SUBNAM( 1: 1 ) )
         IZ = ICHAR( 'Z' )
         IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*           ASCII character set
*
            IF( IC.GE.97 .AND. IC.LE.122 ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( IC.GE.97 .AND. IC.LE.122 )
     $               SUBNAM( I: I ) = CHAR( IC-32 )
               END DO
            END IF
*
         ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*           EBCDIC character set
*
            IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $          ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $          ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC+64 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $                ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $                ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $                I ) = CHAR( IC+64 )
               END DO
            END IF
*
         ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*           Prime machines:  ASCII+128
*
            IF( IC.GE.225 .AND. IC.LE.250 ) THEN
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( IC.GE.225 .AND. IC.LE.250 )
     $               SUBNAM( I: I ) = CHAR( IC-32 )
               END DO
            END IF
         END IF
*
         IF( SUBNAM( 2:6 ).EQ.'GGHRD' .OR.
     $       SUBNAM( 2:6 ).EQ.'GGHD3' ) THEN
            IPARMQ = 1
            IF( NH.GE.K22MIN )
     $         IPARMQ = 2
         ELSE IF ( SUBNAM( 4:6 ).EQ.'EXC' ) THEN
            IF( NH.GE.KACMIN )
     $         IPARMQ = 1
            IF( NH.GE.K22MIN )
     $         IPARMQ = 2
         ELSE IF ( SUBNAM( 2:6 ).EQ.'HSEQR' .OR.
     $             SUBNAM( 2:5 ).EQ.'LAQR' ) THEN
            IF( NS.GE.KACMIN )
     $         IPARMQ = 1
            IF( NS.GE.K22MIN )
     $         IPARMQ = 2
         END IF
*
      ELSE
*        ===== invalid value of ispec =====
         IPARMQ = -1
*
      END IF
*
*     ==== End of IPARMQ ====
*
      END
*> \brief \b LSAME
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*      LOGICAL FUNCTION LSAME( CA, CB )
*
*     .. Scalar Arguments ..
*      CHARACTER          CA, CB
*     ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> LSAME returns .TRUE. if CA is the same letter as CB regardless of
*> case.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] CA
*> \verbatim
*> \endverbatim
*>
*> \param[in] CB
*> \verbatim
*>          CA and CB specify the single characters to be compared.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      LOGICAL FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
*> \brief \b XERBLA
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download XERBLA + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/xerbla.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/xerbla.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/xerbla.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE XERBLA( SRNAME, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER*(*)      SRNAME
*       INTEGER            INFO
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> XERBLA  is an error handler for the LAPACK routines.
*> It is called by an LAPACK routine if an input parameter has an
*> invalid value.  A message is printed and execution stops.
*>
*> Installers may consider modifying the STOP statement in order to
*> call system-specific exception-handling facilities.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SRNAME
*> \verbatim
*>          SRNAME is CHARACTER*(*)
*>          The name of the routine which called XERBLA.
*> \endverbatim
*>
*> \param[in] INFO
*> \verbatim
*>          INFO is INTEGER
*>          The position of the invalid parameter in the parameter list
*>          of the calling routine.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
*     ..
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          LEN_TRIM
*     ..
*     .. Executable Statements ..
*
      WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END
*> \brief \b ZGETF2 computes the LU factorization of a general m-by-n matrix using partial pivoting with row interchanges (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZGETF2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgetf2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgetf2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgetf2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGETF2( M, N, A, LDA, IPIV, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGETF2 computes an LU factorization of a general m-by-n matrix A
*> using partial pivoting with row interchanges.
*>
*> The factorization has the form
*>    A = P * L * U
*> where P is a permutation matrix, L is lower triangular with unit
*> diagonal elements (lower trapezoidal if m > n), and U is upper
*> triangular (upper trapezoidal if m < n).
*>
*> This is the right-looking Level 2 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the m by n matrix to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U; the unit diagonal elements of L are not stored.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (min(M,N))
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -k, the k-th argument had an illegal value
*>          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*>               has been completed, but the factor U is exactly
*>               singular, and division by zero will occur if it is used
*>               to solve a system of equations.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup complex16GEcomputational
*
*  =====================================================================
      SUBROUTINE ZGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN
      INTEGER            I, J, JP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            IZAMAX
      EXTERNAL           DLAMCH, IZAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGERU, ZSCAL, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Compute machine safe minimum
*
      SFMIN = DLAMCH('S') 
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IZAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL ZSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M ) THEN
               IF( ABS(A( J, J )) .GE. SFMIN ) THEN
                  CALL ZSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
               ELSE
                  DO 20 I = 1, M-J
                     A( J+I, J ) = A( J+I, J ) / A( J, J )
   20             CONTINUE
               END IF
            END IF
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL ZGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ),
     $                  LDA, A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of ZGETF2
*
      END
*> \brief \b ZGETRF2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       RECURSIVE SUBROUTINE ZGETRF2( M, N, A, LDA, IPIV, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGETRF2 computes an LU factorization of a general M-by-N matrix A
*> using partial pivoting with row interchanges.
*>
*> The factorization has the form
*>    A = P * L * U
*> where P is a permutation matrix, L is lower triangular with unit
*> diagonal elements (lower trapezoidal if m > n), and U is upper
*> triangular (upper trapezoidal if m < n).
*>
*> This is the recursive version of the algorithm. It divides
*> the matrix into four submatrices:
*>            
*>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
*>    A = [ -----|----- ]  with n1 = min(m,n)
*>        [  A21 | A22  ]       n2 = n-n1
*>            
*>                                       [ A11 ]
*> The subroutine calls itself to factor [ --- ],
*>                                       [ A12 ]
*>                 [ A12 ]
*> do the swaps on [ --- ], solve A12, update A22,
*>                 [ A22 ]
*>
*> then calls itself to factor A22 and do the swaps on A21.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the M-by-N matrix to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U; the unit diagonal elements of L are not stored.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (min(M,N))
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*>                has been completed, but the factor U is exactly
*>                singular, and division by zero will occur if it is used
*>                to solve a system of equations.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup complex16GEcomputational
*
*  =====================================================================
      RECURSIVE SUBROUTINE zgetrf2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( lda, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      parameter( one = ( 1.0d+0, 0.0d+0 ),
     $                     zero = ( 0.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN
      COMPLEX*16         TEMP
      INTEGER            I, IINFO, N1, N2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            IZAMAX
      EXTERNAL           dlamch, izamax
*     ..
*     .. External Subroutines ..
      EXTERNAL           zgemm, zscal, zlaswp, ztrsm, zerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      info = 0
      IF( m.LT.0 ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, m ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGETRF2', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 )
     $   RETURN

      IF ( m.EQ.1 ) THEN
*
*        Use unblocked code for one row case
*        Just need to handle IPIV and INFO
*
         ipiv( 1 ) = 1
         IF ( a(1,1).EQ.zero )
     $      info = 1
*
      ELSE IF( n.EQ.1 ) THEN
*
*        Use unblocked code for one column case
*
*
*        Compute machine safe minimum
*
         sfmin = dlamch('S')
*
*        Find pivot and test for singularity
*
         i = izamax( m, a( 1, 1 ), 1 )
         ipiv( 1 ) = i
         IF( a( i, 1 ).NE.zero ) THEN
*
*           Apply the interchange
*
            IF( i.NE.1 ) THEN
               temp = a( 1, 1 )
               a( 1, 1 ) = a( i, 1 )
               a( i, 1 ) = temp
            END IF
*
*           Compute elements 2:M of the column
*
            IF( abs(a( 1, 1 )) .GE. sfmin ) THEN
               CALL zscal( m-1, one / a( 1, 1 ), a( 2, 1 ), 1 )
            ELSE
               DO 10 i = 1, m-1
                  a( 1+i, 1 ) = a( 1+i, 1 ) / a( 1, 1 )
   10          CONTINUE
            END IF
*
         ELSE
            info = 1
         END IF

      ELSE
*
*        Use recursive code
*
         n1 = min( m, n ) / 2
         n2 = n-n1
*
*               [ A11 ]
*        Factor [ --- ]
*               [ A21 ]
*
         CALL zgetrf2( m, n1, a, lda, ipiv, iinfo )

         IF ( info.EQ.0 .AND. iinfo.GT.0 )
     $      info = iinfo
*
*                              [ A12 ]
*        Apply interchanges to [ --- ]
*                              [ A22 ]
*
         CALL zlaswp( n2, a( 1, n1+1 ), lda, 1, n1, ipiv, 1 )
*
*        Solve A12
*
         CALL ztrsm( 'L', 'L', 'N', 'U', n1, n2, one, a, lda, 
     $               a( 1, n1+1 ), lda )
*
*        Update A22
*
         CALL zgemm( 'N', 'N', m-n1, n2, n1, -one, a( n1+1, 1 ), lda, 
     $               a( 1, n1+1 ), lda, one, a( n1+1, n1+1 ), lda )
*
*        Factor A22
*
         CALL zgetrf2( m-n1, n2, a( n1+1, n1+1 ), lda, ipiv( n1+1 ),
     $                 iinfo )
*
*        Adjust INFO and the pivot indices
*
         IF ( info.EQ.0 .AND. iinfo.GT.0 )
     $      info = iinfo + n1
         DO 20 i = n1+1, min( m, n )
            ipiv( i ) = ipiv( i ) + n1
   20    CONTINUE
*
*        Apply interchanges to A21
*
         CALL zlaswp( n1, a( 1, 1 ), lda, n1+1, min( m, n), ipiv, 1 )
*
      END IF
      RETURN
*
*     End of ZGETRF2
*
      END
*> \brief \b ZGETRF
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZGETRF + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgetrf.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgetrf.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgetrf.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGETRF computes an LU factorization of a general M-by-N matrix A
*> using partial pivoting with row interchanges.
*>
*> The factorization has the form
*>    A = P * L * U
*> where P is a permutation matrix, L is lower triangular with unit
*> diagonal elements (lower trapezoidal if m > n), and U is upper
*> triangular (upper trapezoidal if m < n).
*>
*> This is the right-looking Level 3 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the M-by-N matrix to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U; the unit diagonal elements of L are not stored.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (min(M,N))
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*>                has been completed, but the factor U is exactly
*>                singular, and division by zero will occur if it is used
*>                to solve a system of equations.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup complex16GEcomputational
*
*  =====================================================================
      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZGETRF2, ZLASWP, ZTRSM
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'ZGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL ZGETRF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL ZGETRF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL ZLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL ZLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL ZGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of ZGETRF
*
      END
*> \brief \b ZGETRI
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZGETRI + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgetri.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgetri.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgetri.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LWORK, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         A( LDA, * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGETRI computes the inverse of a matrix using the LU factorization
*> computed by ZGETRF.
*>
*> This method inverts U and then computes inv(A) by solving the system
*> inv(A)*L = inv(U) for inv(A).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the factors L and U from the factorization
*>          A = P*L*U as computed by ZGETRF.
*>          On exit, if INFO = 0, the inverse of the original matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          The pivot indices from ZGETRF; for 1<=i<=N, row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,N).
*>          For optimal performance LWORK >= N*NB, where NB is
*>          the optimal blocksize returned by ILAENV.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
*>                singular and its inverse could not be computed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup complex16GEcomputational
*
*  =====================================================================
      SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB,
     $                   NBMIN, NN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZGEMV, ZSWAP, ZTRSM, ZTRTRI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NB = ILAENV( 1, 'ZGETRI', ' ', N, -1, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRI', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Form inv(U).  If INFO > 0 from ZTRTRI, then U is singular,
*     and the inverse is not computed.
*
      CALL ZTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
      IF( INFO.GT.0 )
     $   RETURN
*
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = MAX( LDWORK*NB, 1 )
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'ZGETRI', ' ', N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = N
      END IF
*
*     Solve the equation inv(A)*L = inv(U) for inv(A).
*
      IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         DO 20 J = N, 1, -1
*
*           Copy current column of L to WORK and replace with zeros.
*
            DO 10 I = J + 1, N
               WORK( I ) = A( I, J )
               A( I, J ) = ZERO
   10       CONTINUE
*
*           Compute current column of inv(A).
*
            IF( J.LT.N )
     $         CALL ZGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
   20    CONTINUE
      ELSE
*
*        Use blocked code.
*
         NN = ( ( N-1 ) / NB )*NB + 1
         DO 50 J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )
*
*           Copy current block column of L to WORK and replace with
*           zeros.
*
            DO 40 JJ = J, J + JB - 1
               DO 30 I = JJ + 1, N
                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
                  A( I, JJ ) = ZERO
   30          CONTINUE
   40       CONTINUE
*
*           Compute current block column of inv(A).
*
            IF( J+JB.LE.N )
     $         CALL ZGEMM( 'No transpose', 'No transpose', N, JB,
     $                     N-J-JB+1, -ONE, A( 1, J+JB ), LDA,
     $                     WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
            CALL ZTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
     $                  ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
   50    CONTINUE
      END IF
*
*     Apply column interchanges.
*
      DO 60 J = N - 1, 1, -1
         JP = IPIV( J )
         IF( JP.NE.J )
     $      CALL ZSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
   60 CONTINUE
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of ZGETRI
*
      END
*> \brief \b ZLASWP performs a series of row interchanges on a general rectangular matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZLASWP + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaswp.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaswp.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaswp.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
* 
*       .. Scalar Arguments ..
*       INTEGER            INCX, K1, K2, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLASWP performs a series of row interchanges on the matrix A.
*> One row interchange is initiated for each of rows K1 through K2 of A.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the matrix of column dimension N to which the row
*>          interchanges will be applied.
*>          On exit, the permuted matrix.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*> \endverbatim
*>
*> \param[in] K1
*> \verbatim
*>          K1 is INTEGER
*>          The first element of IPIV for which a row interchange will
*>          be done.
*> \endverbatim
*>
*> \param[in] K2
*> \verbatim
*>          K2 is INTEGER
*>          The last element of IPIV for which a row interchange will
*>          be done.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (K2*abs(INCX))
*>          The vector of pivot indices.  Only the elements in positions
*>          K1 through K2 of IPIV are accessed.
*>          IPIV(K) = L implies rows K and L are to be interchanged.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between successive values of IPIV.  If IPIV
*>          is negative, the pivots are applied in reverse order.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup complex16OTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Modified by
*>   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      COMPLEX*16         TEMP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
*
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLASWP
*
      END
*> \brief \b ZLASYF computes a partial factorization of a complex symmetric matrix using the Bunch-Kaufman diagonal pivoting method.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLASYF + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlasyf.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlasyf.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlasyf.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, KB, LDA, LDW, N, NB
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         A( LDA, * ), W( LDW, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLASYF computes a partial factorization of a complex symmetric matrix
*> A using the Bunch-Kaufman diagonal pivoting method. The partial
*> factorization has the form:
*>
*> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
*>       ( 0  U22 ) (  0   D  ) ( U12**T U22**T )
*>
*> A  =  ( L11  0 ) ( D    0  ) ( L11**T L21**T )  if UPLO = 'L'
*>       ( L21  I ) ( 0   A22 ) (  0       I    )
*>
*> where the order of D is at most NB. The actual order is returned in
*> the argument KB, and is either NB or NB-1, or N if N <= NB.
*> Note that U**T denotes the transpose of U.
*>
*> ZLASYF is an auxiliary routine called by ZSYTRF. It uses blocked code
*> (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
*> A22 (if UPLO = 'L').
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          symmetric matrix A is stored:
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The maximum number of columns of the matrix A that should be
*>          factored.  NB should be at least 2 to allow for 2-by-2 pivot
*>          blocks.
*> \endverbatim
*>
*> \param[out] KB
*> \verbatim
*>          KB is INTEGER
*>          The number of columns of A that were actually factored.
*>          KB is either NB-1 or NB, or N if N <= NB.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*>          n-by-n upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading n-by-n lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the strictly upper
*>          triangular part of A is not referenced.
*>          On exit, A contains details of the partial factorization.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D.
*>
*>          If UPLO = 'U':
*>             Only the last KB elements of IPIV are set.
*>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
*>             interchanged and D(k,k) is a 1-by-1 diagonal block.
*>
*>             If IPIV(k) = IPIV(k-1) < 0, then rows and columns
*>             k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
*>             is a 2-by-2 diagonal block.
*>
*>          If UPLO = 'L':
*>             Only the first KB elements of IPIV are set.
*>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
*>             interchanged and D(k,k) is a 1-by-1 diagonal block.
*>
*>             If IPIV(k) = IPIV(k+1) < 0, then rows and columns
*>             k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1)
*>             is a 2-by-2 diagonal block.
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is COMPLEX*16 array, dimension (LDW,NB)
*> \endverbatim
*>
*> \param[in] LDW
*> \verbatim
*>          LDW is INTEGER
*>          The leading dimension of the array W.  LDW >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
*>               has been completed, but the block diagonal matrix D is
*>               exactly singular.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2013
*
*> \ingroup complex16SYcomputational
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>  November 2013,  Igor Kozachenko,
*>                  Computer Science Division,
*>                  University of California, Berkeley
*> \endverbatim
*
*  =====================================================================
      SUBROUTINE ZLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
*
*  -- LAPACK computational routine (version 3.5.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2013
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KB, LDA, LDW, N, NB
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), W( LDW, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            IMAX, J, JB, JJ, JMAX, JP, K, KK, KKW, KP,
     $                   KSTEP, KW
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, ROWMAX
      COMPLEX*16         D11, D21, D22, R1, T, Z
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IZAMAX
      EXTERNAL           LSAME, IZAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZCOPY, ZGEMM, ZGEMV, ZSCAL, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX, MIN, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( Z ) = ABS( DBLE( Z ) ) + ABS( DIMAG( Z ) )
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Initialize ALPHA for use in choosing pivot block size.
*
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Factorize the trailing columns of A using the upper triangle
*        of A and working backwards, and compute the matrix W = U12*D
*        for use in updating A11
*
*        K is the main loop index, decreasing from N in steps of 1 or 2
*
*        KW is the column of W which corresponds to column K of A
*
         K = N
   10    CONTINUE
         KW = NB + K - N
*
*        Exit from loop
*
         IF( ( K.LE.N-NB+1 .AND. NB.LT.N ) .OR. K.LT.1 )
     $      GO TO 30
*
*        Copy column K of A to column KW of W and update it
*
         CALL ZCOPY( K, A( 1, K ), 1, W( 1, KW ), 1 )
         IF( K.LT.N )
     $      CALL ZGEMV( 'No transpose', K, N-K, -CONE, A( 1, K+1 ), LDA,
     $                  W( K, KW+1 ), LDW, CONE, W( 1, KW ), 1 )
*
         KSTEP = 1
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = CABS1( W( K, KW ) )
*
*        IMAX is the row-index of the largest off-diagonal element in

*
         IF( K.GT.1 ) THEN
            IMAX = IZAMAX( K-1, W( 1, KW ), 1 )
            COLMAX = CABS1( W( IMAX, KW ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*
*           Column K is zero or underflow: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              no interchange, use 1-by-1 pivot block
*
               KP = K
            ELSE
*
*              Copy column IMAX to column KW-1 of W and update it
*
               CALL ZCOPY( IMAX, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 )
               CALL ZCOPY( K-IMAX, A( IMAX, IMAX+1 ), LDA,
     $                     W( IMAX+1, KW-1 ), 1 )
               IF( K.LT.N )
     $            CALL ZGEMV( 'No transpose', K, N-K, -CONE,
     $                        A( 1, K+1 ), LDA, W( IMAX, KW+1 ), LDW,
     $                        CONE, W( 1, KW-1 ), 1 )
*
*              JMAX is the column-index of the largest off-diagonal
*              element in row IMAX, and ROWMAX is its absolute value
*
               JMAX = IMAX + IZAMAX( K-IMAX, W( IMAX+1, KW-1 ), 1 )
               ROWMAX = CABS1( W( JMAX, KW-1 ) )
               IF( IMAX.GT.1 ) THEN
                  JMAX = IZAMAX( IMAX-1, W( 1, KW-1 ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( W( JMAX, KW-1 ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 no interchange, use 1-by-1 pivot block
*
                  KP = K
               ELSE IF( CABS1( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 interchange rows and columns K and IMAX, use 1-by-1
*                 pivot block
*
                  KP = IMAX
*
*                 copy column KW-1 of W to column KW of W
*
                  CALL ZCOPY( K, W( 1, KW-1 ), 1, W( 1, KW ), 1 )
               ELSE
*
*                 interchange rows and columns K-1 and IMAX, use 2-by-2
*                 pivot block
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
*           ============================================================
*
*           KK is the column of A where pivoting step stopped
*
            KK = K - KSTEP + 1
*
*           KKW is the column of W which corresponds to column KK of A
*
            KKW = NB + KK - N
*
*           Interchange rows and columns KP and KK.
*           Updated column KP is already stored in column KKW of W.
*
            IF( KP.NE.KK ) THEN
*
*              Copy non-updated column KK to column KP of submatrix A
*              at step K. No need to copy element into column K
*              (or K and K-1 for 2-by-2 pivot) of A, since these columns
*              will be later overwritten.
*
               A( KP, KP ) = A( KK, KK )
               CALL ZCOPY( KK-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ),
     $                     LDA )
               IF( KP.GT.1 )
     $            CALL ZCOPY( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )
*
*              Interchange rows KK and KP in last K+1 to N columns of A
*              (columns K (or K and K-1 for 2-by-2 pivot) of A will be
*              later overwritten). Interchange rows KK and KP
*              in last KKW to NB columns of W.
*
               IF( K.LT.N )
     $            CALL ZSWAP( N-K, A( KK, K+1 ), LDA, A( KP, K+1 ),
     $                        LDA )
               CALL ZSWAP( N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ),
     $                     LDW )
            END IF
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column kw of W now holds
*
*              W(kw) = U(k)*D(k),
*
*              where U(k) is the k-th column of U
*
*              Store subdiag. elements of column U(k)
*              and 1-by-1 block D(k) in column k of A.
*              NOTE: Diagonal element U(k,k) is a UNIT element
*              and not stored.
*                 A(k,k) := D(k,k) = W(k,kw)
*                 A(1:k-1,k) := U(1:k-1,k) = W(1:k-1,kw)/D(k,k)
*
               CALL ZCOPY( K, W( 1, KW ), 1, A( 1, K ), 1 )
               R1 = CONE / A( K, K )
               CALL ZSCAL( K-1, R1, A( 1, K ), 1 )
*
            ELSE
*
*              2-by-2 pivot block D(k): columns kw and kw-1 of W now hold
*
*              ( W(kw-1) W(kw) ) = ( U(k-1) U(k) )*D(k)
*
*              where U(k) and U(k-1) are the k-th and (k-1)-th columns
*              of U
*
*              Store U(1:k-2,k-1) and U(1:k-2,k) and 2-by-2
*              block D(k-1:k,k-1:k) in columns k-1 and k of A.
*              NOTE: 2-by-2 diagonal block U(k-1:k,k-1:k) is a UNIT
*              block and not stored.
*                 A(k-1:k,k-1:k) := D(k-1:k,k-1:k) = W(k-1:k,kw-1:kw)
*                 A(1:k-2,k-1:k) := U(1:k-2,k:k-1:k) =
*                 = W(1:k-2,kw-1:kw) * ( D(k-1:k,k-1:k)**(-1) )
*
               IF( K.GT.2 ) THEN
*
*                 Compose the columns of the inverse of 2-by-2 pivot
*                 block D in the following way to reduce the number
*                 of FLOPS when we myltiply panel ( W(kw-1) W(kw) ) by
*                 this inverse
*
*                 D**(-1) = ( d11 d21 )**(-1) =
*                           ( d21 d22 )
*
*                 = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
*                                        ( (-d21 ) ( d11 ) )
*
*                 = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
*
*                   * ( ( d22/d21 ) (      -1 ) ) =
*                     ( (      -1 ) ( d11/d21 ) )
*
*                 = 1/d21 * 1/(D22*D11-1) * ( ( D11 ) (  -1 ) ) =
*                                           ( ( -1  ) ( D22 ) )
*
*                 = 1/d21 * T * ( ( D11 ) (  -1 ) )
*                               ( (  -1 ) ( D22 ) )
*
*                 = D21 * ( ( D11 ) (  -1 ) )
*                         ( (  -1 ) ( D22 ) )
*
                  D21 = W( K-1, KW )
                  D11 = W( K, KW ) / D21
                  D22 = W( K-1, KW-1 ) / D21
                  T = CONE / ( D11*D22-CONE )
                  D21 = T / D21
*
*                 Update elements in columns A(k-1) and A(k) as
*                 dot products of rows of ( W(kw-1) W(kw) ) and columns
*                 of D**(-1)
*
                  DO 20 J = 1, K - 2
                     A( J, K-1 ) = D21*( D11*W( J, KW-1 )-W( J, KW ) )
                     A( J, K ) = D21*( D22*W( J, KW )-W( J, KW-1 ) )
   20             CONTINUE
               END IF
*
*              Copy D(k) to A
*
               A( K-1, K-1 ) = W( K-1, KW-1 )
               A( K-1, K ) = W( K-1, KW )
               A( K, K ) = W( K, KW )
*
            END IF
*
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         END IF
*
*        Decrease K and return to the start of the main loop
*
         K = K - KSTEP
         GO TO 10
*
   30    CONTINUE
*
*        Update the upper triangle of A11 (= A(1:k,1:k)) as
*
*        A11 := A11 - U12*D*U12**T = A11 - U12*W**T
*
*        computing blocks of NB columns at a time
*
         DO 50 J = ( ( K-1 ) / NB )*NB + 1, 1, -NB
            JB = MIN( NB, K-J+1 )
*
*           Update the upper triangle of the diagonal block
*
            DO 40 JJ = J, J + JB - 1
               CALL ZGEMV( 'No transpose', JJ-J+1, N-K, -CONE,
     $                     A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, CONE,
     $                     A( J, JJ ), 1 )
   40       CONTINUE
*
*           Update the rectangular superdiagonal block
*
            CALL ZGEMM( 'No transpose', 'Transpose', J-1, JB, N-K,
     $                  -CONE, A( 1, K+1 ), LDA, W( J, KW+1 ), LDW,
     $                  CONE, A( 1, J ), LDA )
   50    CONTINUE
*
*        Put U12 in standard form by partially undoing the interchanges
*        in columns k+1:n looping backwards from k+1 to n
*
         J = K + 1
   60    CONTINUE
*
*           Undo the interchanges (if any) of rows JJ and JP at each
*           step J
*
*           (Here, J is a diagonal index)
            JJ = J
            JP = IPIV( J )
            IF( JP.LT.0 ) THEN
               JP = -JP
*              (Here, J is a diagonal index)
               J = J + 1
            END IF
*           (NOTE: Here, J is used to determine row length. Length N-J+1
*           of the rows to swap back doesn't include diagonal element)
            J = J + 1
            IF( JP.NE.JJ .AND. J.LE.N )
     $         CALL ZSWAP( N-J+1, A( JP, J ), LDA, A( JJ, J ), LDA )
         IF( J.LT.N )
     $      GO TO 60
*
*        Set KB to the number of columns factorized
*
         KB = N - K
*
      ELSE
*
*        Factorize the leading columns of A using the lower triangle
*        of A and working forwards, and compute the matrix W = L21*D
*        for use in updating A22
*
*        K is the main loop index, increasing from 1 in steps of 1 or 2
*
         K = 1
   70    CONTINUE
*
*        Exit from loop
*
         IF( ( K.GE.NB .AND. NB.LT.N ) .OR. K.GT.N )
     $      GO TO 90
*
*        Copy column K of A to column K of W and update it
*
         CALL ZCOPY( N-K+1, A( K, K ), 1, W( K, K ), 1 )
         CALL ZGEMV( 'No transpose', N-K+1, K-1, -CONE, A( K, 1 ), LDA,
     $               W( K, 1 ), LDW, CONE, W( K, K ), 1 )
*
         KSTEP = 1
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = CABS1( W( K, K ) )
*
*        IMAX is the row-index of the largest off-diagonal element in

*
         IF( K.LT.N ) THEN
            IMAX = K + IZAMAX( N-K, W( K+1, K ), 1 )
            COLMAX = CABS1( W( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*
*           Column K is zero or underflow: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              no interchange, use 1-by-1 pivot block
*
               KP = K
            ELSE
*
*              Copy column IMAX to column K+1 of W and update it
*
               CALL ZCOPY( IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1 )
               CALL ZCOPY( N-IMAX+1, A( IMAX, IMAX ), 1, W( IMAX, K+1 ),
     $                     1 )
               CALL ZGEMV( 'No transpose', N-K+1, K-1, -CONE, A( K, 1 ),
     $                     LDA, W( IMAX, 1 ), LDW, CONE, W( K, K+1 ),
     $                     1 )
*
*              JMAX is the column-index of the largest off-diagonal
*              element in row IMAX, and ROWMAX is its absolute value
*
               JMAX = K - 1 + IZAMAX( IMAX-K, W( K, K+1 ), 1 )
               ROWMAX = CABS1( W( JMAX, K+1 ) )
               IF( IMAX.LT.N ) THEN
                  JMAX = IMAX + IZAMAX( N-IMAX, W( IMAX+1, K+1 ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( W( JMAX, K+1 ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 no interchange, use 1-by-1 pivot block
*
                  KP = K
               ELSE IF( CABS1( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 interchange rows and columns K and IMAX, use 1-by-1
*                 pivot block
*
                  KP = IMAX
*
*                 copy column K+1 of W to column K of W
*
                  CALL ZCOPY( N-K+1, W( K, K+1 ), 1, W( K, K ), 1 )
               ELSE
*
*                 interchange rows and columns K+1 and IMAX, use 2-by-2
*                 pivot block
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
*           ============================================================
*
*           KK is the column of A where pivoting step stopped
*
            KK = K + KSTEP - 1
*
*           Interchange rows and columns KP and KK.
*           Updated column KP is already stored in column KK of W.
*
            IF( KP.NE.KK ) THEN
*
*              Copy non-updated column KK to column KP of submatrix A
*              at step K. No need to copy element into column K
*              (or K and K+1 for 2-by-2 pivot) of A, since these columns
*              will be later overwritten.
*
               A( KP, KP ) = A( KK, KK )
               CALL ZCOPY( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ),
     $                     LDA )
               IF( KP.LT.N )
     $            CALL ZCOPY( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
*
*              Interchange rows KK and KP in first K-1 columns of A
*              (columns K (or K and K+1 for 2-by-2 pivot) of A will be
*              later overwritten). Interchange rows KK and KP
*              in first KK columns of W.
*
               IF( K.GT.1 )
     $            CALL ZSWAP( K-1, A( KK, 1 ), LDA, A( KP, 1 ), LDA )
               CALL ZSWAP( KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW )
            END IF
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column k of W now holds
*
*              W(k) = L(k)*D(k),
*
*              where L(k) is the k-th column of L
*
*              Store subdiag. elements of column L(k)
*              and 1-by-1 block D(k) in column k of A.
*              (NOTE: Diagonal element L(k,k) is a UNIT element
*              and not stored)
*                 A(k,k) := D(k,k) = W(k,k)
*                 A(k+1:N,k) := L(k+1:N,k) = W(k+1:N,k)/D(k,k)
*
               CALL ZCOPY( N-K+1, W( K, K ), 1, A( K, K ), 1 )
               IF( K.LT.N ) THEN
                  R1 = CONE / A( K, K )
                  CALL ZSCAL( N-K, R1, A( K+1, K ), 1 )
               END IF
*
            ELSE
*
*              2-by-2 pivot block D(k): columns k and k+1 of W now hold
*
*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
*
*              where L(k) and L(k+1) are the k-th and (k+1)-th columns
*              of L
*
*              Store L(k+2:N,k) and L(k+2:N,k+1) and 2-by-2
*              block D(k:k+1,k:k+1) in columns k and k+1 of A.
*              (NOTE: 2-by-2 diagonal block L(k:k+1,k:k+1) is a UNIT
*              block and not stored)
*                 A(k:k+1,k:k+1) := D(k:k+1,k:k+1) = W(k:k+1,k:k+1)
*                 A(k+2:N,k:k+1) := L(k+2:N,k:k+1) =
*                 = W(k+2:N,k:k+1) * ( D(k:k+1,k:k+1)**(-1) )
*
               IF( K.LT.N-1 ) THEN
*
*                 Compose the columns of the inverse of 2-by-2 pivot
*                 block D in the following way to reduce the number
*                 of FLOPS when we myltiply panel ( W(k) W(k+1) ) by
*                 this inverse
*
*                 D**(-1) = ( d11 d21 )**(-1) =
*                           ( d21 d22 )
*
*                 = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
*                                        ( (-d21 ) ( d11 ) )
*
*                 = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
*
*                   * ( ( d22/d21 ) (      -1 ) ) =
*                     ( (      -1 ) ( d11/d21 ) )
*
*                 = 1/d21 * 1/(D22*D11-1) * ( ( D11 ) (  -1 ) ) =
*                                           ( ( -1  ) ( D22 ) )
*
*                 = 1/d21 * T * ( ( D11 ) (  -1 ) )
*                               ( (  -1 ) ( D22 ) )
*
*                 = D21 * ( ( D11 ) (  -1 ) )
*                         ( (  -1 ) ( D22 ) )
*
                  D21 = W( K+1, K )
                  D11 = W( K+1, K+1 ) / D21
                  D22 = W( K, K ) / D21
                  T = CONE / ( D11*D22-CONE )
                  D21 = T / D21
*
*                 Update elements in columns A(k) and A(k+1) as
*                 dot products of rows of ( W(k) W(k+1) ) and columns
*                 of D**(-1)
*
                  DO 80 J = K + 2, N
                     A( J, K ) = D21*( D11*W( J, K )-W( J, K+1 ) )
                     A( J, K+1 ) = D21*( D22*W( J, K+1 )-W( J, K ) )
   80             CONTINUE
               END IF
*
*              Copy D(k) to A
*
               A( K, K ) = W( K, K )
               A( K+1, K ) = W( K+1, K )
               A( K+1, K+1 ) = W( K+1, K+1 )
*
            END IF
*
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         END IF
*
*        Increase K and return to the start of the main loop
*
         K = K + KSTEP
         GO TO 70
*
   90    CONTINUE
*
*        Update the lower triangle of A22 (= A(k:n,k:n)) as
*
*        A22 := A22 - L21*D*L21**T = A22 - L21*W**T
*
*        computing blocks of NB columns at a time
*
         DO 110 J = K, N, NB
            JB = MIN( NB, N-J+1 )
*
*           Update the lower triangle of the diagonal block
*
            DO 100 JJ = J, J + JB - 1
               CALL ZGEMV( 'No transpose', J+JB-JJ, K-1, -CONE,
     $                     A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, CONE,
     $                     A( JJ, JJ ), 1 )
  100       CONTINUE
*
*           Update the rectangular subdiagonal block
*
            IF( J+JB.LE.N )
     $         CALL ZGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB,
     $                     K-1, -CONE, A( J+JB, 1 ), LDA, W( J, 1 ),
     $                     LDW, CONE, A( J+JB, J ), LDA )
  110    CONTINUE
*
*        Put L21 in standard form by partially undoing the interchanges
*        of rows in columns 1:k-1 looping backwards from k-1 to 1
*
         J = K - 1
  120    CONTINUE
*
*           Undo the interchanges (if any) of rows JJ and JP at each
*           step J
*
*           (Here, J is a diagonal index)
            JJ = J
            JP = IPIV( J )
            IF( JP.LT.0 ) THEN
               JP = -JP
*              (Here, J is a diagonal index)
               J = J - 1
            END IF
*           (NOTE: Here, J is used to determine row length. Length J
*           of the rows to swap back doesn't include diagonal element)
            J = J - 1
            IF( JP.NE.JJ .AND. J.GE.1 )
     $         CALL ZSWAP( J, A( JP, 1 ), LDA, A( JJ, 1 ), LDA )
         IF( J.GT.1 )
     $      GO TO 120
*
*        Set KB to the number of columns factorized
*
         KB = K - 1
*
      END IF
      RETURN
*
*     End of ZLASYF
*
      END
*> \brief \b ZSYMV computes a matrix-vector product for a complex symmetric matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZSYMV + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsymv.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsymv.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsymv.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZSYMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
* 
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INCX, INCY, LDA, N
*       COMPLEX*16         ALPHA, BETA
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), X( * ), Y( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZSYMV  performs the matrix-vector  operation
*>
*>    y := alpha*A*x + beta*y,
*>
*> where alpha and beta are scalars, x and y are n element vectors and
*> A is an n by n symmetric matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the upper or lower
*>           triangular part of the array A is to be referenced as
*>           follows:
*>
*>              UPLO = 'U' or 'u'   Only the upper triangular part of A
*>                                  is to be referenced.
*>
*>              UPLO = 'L' or 'l'   Only the lower triangular part of A
*>                                  is to be referenced.
*>
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the order of the matrix A.
*>           N must be at least zero.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>           On entry, ALPHA specifies the scalar alpha.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension ( LDA, N )
*>           Before entry, with  UPLO = 'U' or 'u', the leading n by n
*>           upper triangular part of the array A must contain the upper
*>           triangular part of the symmetric matrix and the strictly
*>           lower triangular part of A is not referenced.
*>           Before entry, with UPLO = 'L' or 'l', the leading n by n
*>           lower triangular part of the array A must contain the lower
*>           triangular part of the symmetric matrix and the strictly
*>           upper triangular part of A is not referenced.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, N ).
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension at least
*>           ( 1 + ( N - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the N-
*>           element vector x.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is COMPLEX*16
*>           On entry, BETA specifies the scalar beta. When BETA is
*>           supplied as zero then Y need not be set on input.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in,out] Y
*> \verbatim
*>          Y is COMPLEX*16 array, dimension at least
*>           ( 1 + ( N - 1 )*abs( INCY ) ).
*>           Before entry, the incremented array Y must contain the n
*>           element vector y. On exit, Y is overwritten by the updated
*>           vector y.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           Y. INCY must not be zero.
*>           Unchanged on exit.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup complex16SYauxiliary
*
*  =====================================================================
      SUBROUTINE ZSYMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INCX, INCY, LDA, N
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
      COMPLEX*16         TEMP1, TEMP2
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = 1
      ELSE IF( N.LT.0 ) THEN
         INFO = 2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = 5
      ELSE IF( INCX.EQ.0 ) THEN
         INFO = 7
      ELSE IF( INCY.EQ.0 ) THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSYMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ) .OR. ( ( ALPHA.EQ.ZERO ) .AND. ( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set up the start points in  X  and  Y.
*
      IF( INCX.GT.0 ) THEN
         KX = 1
      ELSE
         KX = 1 - ( N-1 )*INCX
      END IF
      IF( INCY.GT.0 ) THEN
         KY = 1
      ELSE
         KY = 1 - ( N-1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE ) THEN
         IF( INCY.EQ.1 ) THEN
            IF( BETA.EQ.ZERO ) THEN
               DO 10 I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20 I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO ) THEN
               DO 30 I = 1, N
                  Y( IY ) = ZERO
                  IY = IY + INCY
   30          CONTINUE
            ELSE
               DO 40 I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY = IY + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Form  y  when A is stored in upper triangle.
*
         IF( ( INCX.EQ.1 ) .AND. ( INCY.EQ.1 ) ) THEN
            DO 60 J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               DO 50 I = 1, J - 1
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2 = TEMP2 + A( I, J )*X( I )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*A( J, J ) + ALPHA*TEMP2
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80 J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX = KX
               IY = KY
               DO 70 I = 1, J - 1
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2 = TEMP2 + A( I, J )*X( IX )
                  IX = IX + INCX
                  IY = IY + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*A( J, J ) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y  when A is stored in lower triangle.
*
         IF( ( INCX.EQ.1 ) .AND. ( INCY.EQ.1 ) ) THEN
            DO 100 J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               Y( J ) = Y( J ) + TEMP1*A( J, J )
               DO 90 I = J + 1, N
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2 = TEMP2 + A( I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120 J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               Y( JY ) = Y( JY ) + TEMP1*A( J, J )
               IX = JX
               IY = JY
               DO 110 I = J + 1, N
                  IX = IX + INCX
                  IY = IY + INCY
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2 = TEMP2 + A( I, J )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZSYMV
*
      END
*> \brief \b ZSYR performs the symmetric rank-1 update of a complex symmetric matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZSYR + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyr.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyr.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyr.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZSYR( UPLO, N, ALPHA, X, INCX, A, LDA )
* 
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INCX, LDA, N
*       COMPLEX*16         ALPHA
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), X( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZSYR   performs the symmetric rank 1 operation
*>
*>    A := alpha*x*x**H + A,
*>
*> where alpha is a complex scalar, x is an n element vector and A is an
*> n by n symmetric matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the upper or lower
*>           triangular part of the array A is to be referenced as
*>           follows:
*>
*>              UPLO = 'U' or 'u'   Only the upper triangular part of A
*>                                  is to be referenced.
*>
*>              UPLO = 'L' or 'l'   Only the lower triangular part of A
*>                                  is to be referenced.
*>
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the order of the matrix A.
*>           N must be at least zero.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>           On entry, ALPHA specifies the scalar alpha.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension at least
*>           ( 1 + ( N - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the N-
*>           element vector x.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension ( LDA, N )
*>           Before entry, with  UPLO = 'U' or 'u', the leading n by n
*>           upper triangular part of the array A must contain the upper
*>           triangular part of the symmetric matrix and the strictly
*>           lower triangular part of A is not referenced. On exit, the
*>           upper triangular part of the array A is overwritten by the
*>           upper triangular part of the updated matrix.
*>           Before entry, with UPLO = 'L' or 'l', the leading n by n
*>           lower triangular part of the array A must contain the lower
*>           triangular part of the symmetric matrix and the strictly
*>           upper triangular part of A is not referenced. On exit, the
*>           lower triangular part of the array A is overwritten by the
*>           lower triangular part of the updated matrix.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, N ).
*>           Unchanged on exit.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup complex16SYauxiliary
*
*  =====================================================================
      SUBROUTINE ZSYR( UPLO, N, ALPHA, X, INCX, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INCX, LDA, N
      COMPLEX*16         ALPHA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, IX, J, JX, KX
      COMPLEX*16         TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = 1
      ELSE IF( N.LT.0 ) THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSYR  ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ) .OR. ( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Set the start point in X if the increment is not unity.
*
      IF( INCX.LE.0 ) THEN
         KX = 1 - ( N-1 )*INCX
      ELSE IF( INCX.NE.1 ) THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Form  A  when A is stored in upper triangle.
*
         IF( INCX.EQ.1 ) THEN
            DO 20 J = 1, N
               IF( X( J ).NE.ZERO ) THEN
                  TEMP = ALPHA*X( J )
                  DO 10 I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE
            JX = KX
            DO 40 J = 1, N
               IF( X( JX ).NE.ZERO ) THEN
                  TEMP = ALPHA*X( JX )
                  IX = KX
                  DO 30 I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX = IX + INCX
   30             CONTINUE
               END IF
               JX = JX + INCX
   40       CONTINUE
         END IF
      ELSE
*
*        Form  A  when A is stored in lower triangle.
*
         IF( INCX.EQ.1 ) THEN
            DO 60 J = 1, N
               IF( X( J ).NE.ZERO ) THEN
                  TEMP = ALPHA*X( J )
                  DO 50 I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            DO 80 J = 1, N
               IF( X( JX ).NE.ZERO ) THEN
                  TEMP = ALPHA*X( JX )
                  IX = JX
                  DO 70 I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX = IX + INCX
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZSYR
*
      END
*> \brief \b ZSYTF2 computes the factorization of a real symmetric indefinite matrix, using the diagonal pivoting method (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZSYTF2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytf2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytf2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytf2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZSYTF2( UPLO, N, A, LDA, IPIV, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZSYTF2 computes the factorization of a complex symmetric matrix A
*> using the Bunch-Kaufman diagonal pivoting method:
*>
*>    A = U*D*U**T  or  A = L*D*L**T
*>
*> where U (or L) is a product of permutation and unit upper (lower)
*> triangular matrices, U**T is the transpose of U, and D is symmetric and
*> block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
*>
*> This is the unblocked version of the algorithm, calling Level 2 BLAS.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          symmetric matrix A is stored:
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*>          n-by-n upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading n-by-n lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the strictly upper
*>          triangular part of A is not referenced.
*>
*>          On exit, the block diagonal matrix D and the multipliers used
*>          to obtain the factor U or L (see below for further details).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D.
*>
*>          If UPLO = 'U':
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
*>             interchanged and D(k,k) is a 1-by-1 diagonal block.
*>
*>             If IPIV(k) = IPIV(k-1) < 0, then rows and columns
*>             k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
*>             is a 2-by-2 diagonal block.
*>
*>          If UPLO = 'L':
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
*>             interchanged and D(k,k) is a 1-by-1 diagonal block.
*>
*>             If IPIV(k) = IPIV(k+1) < 0, then rows and columns
*>             k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1)
*>             is a 2-by-2 diagonal block.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -k, the k-th argument had an illegal value
*>          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
*>               has been completed, but the block diagonal matrix D is
*>               exactly singular, and division by zero will occur if it
*>               is used to solve a system of equations.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2013
*
*> \ingroup complex16SYcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  If UPLO = 'U', then A = U*D*U**T, where
*>     U = P(n)*U(n)* ... *P(k)U(k)* ...,
*>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
*>  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
*>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
*>  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
*>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
*>
*>             (   I    v    0   )   k-s
*>     U(k) =  (   0    I    0   )   s
*>             (   0    0    I   )   n-k
*>                k-s   s   n-k
*>
*>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
*>  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
*>  and A(k,k), and v overwrites A(1:k-2,k-1:k).
*>
*>  If UPLO = 'L', then A = L*D*L**T, where
*>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
*>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
*>  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
*>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
*>  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
*>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
*>
*>             (   I    0     0   )  k-1
*>     L(k) =  (   0    I     0   )  s
*>             (   0    v     I   )  n-k-s+1
*>                k-1   s  n-k-s+1
*>
*>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
*>  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
*>  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
*> \endverbatim
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>  09-29-06 - patch from
*>    Bobby Cheng, MathWorks
*>
*>    Replace l.209 and l.377
*>         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*>    by
*>         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN
*>
*>  1-96 - Based on modifications by J. Lewis, Boeing Computer Services
*>         Company
*> \endverbatim
*
*  =====================================================================
      SUBROUTINE ZSYTF2( UPLO, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine (version 3.5.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2013
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IMAX, J, JMAX, K, KK, KP, KSTEP
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, ROWMAX
      COMPLEX*16         D11, D12, D21, D22, R1, T, WK, WKM1, WKP1, Z
*     ..
*     .. External Functions ..
      LOGICAL            DISNAN, LSAME
      INTEGER            IZAMAX
      EXTERNAL           DISNAN, LSAME, IZAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZSCAL, ZSWAP, ZSYR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( Z ) = ABS( DBLE( Z ) ) + ABS( DIMAG( Z ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSYTF2', -INFO )
         RETURN
      END IF
*
*     Initialize ALPHA for use in choosing pivot block size.
*
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
*
      IF( UPPER ) THEN
*
*        Factorize A as U*D*U**T using the upper triangle of A
*
*        K is the main loop index, decreasing from N to 1 in steps of
*        1 or 2
*
         K = N
   10    CONTINUE
*
*        If K < 1, exit from loop
*
         IF( K.LT.1 )
     $      GO TO 70
         KSTEP = 1
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = CABS1( A( K, K ) )
*
*        IMAX is the row-index of the largest off-diagonal element in
*        column K, and COLMAX is its absolute value.
*        Determine both COLMAX and IMAX.
*
         IF( K.GT.1 ) THEN
            IMAX = IZAMAX( K-1, A( 1, K ), 1 )
            COLMAX = CABS1( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO .OR. DISNAN(ABSAKK) ) THEN
*
*           Column K is zero or underflow, or contains a NaN:
*           set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              no interchange, use 1-by-1 pivot block
*
               KP = K
            ELSE
*
*              JMAX is the column-index of the largest off-diagonal
*              element in row IMAX, and ROWMAX is its absolute value
*
               JMAX = IMAX + IZAMAX( K-IMAX, A( IMAX, IMAX+1 ), LDA )
               ROWMAX = CABS1( A( IMAX, JMAX ) )
               IF( IMAX.GT.1 ) THEN
                  JMAX = IZAMAX( IMAX-1, A( 1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( A( JMAX, IMAX ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 no interchange, use 1-by-1 pivot block
*
                  KP = K
               ELSE IF( CABS1( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 interchange rows and columns K and IMAX, use 1-by-1
*                 pivot block
*
                  KP = IMAX
               ELSE
*
*                 interchange rows and columns K-1 and IMAX, use 2-by-2
*                 pivot block
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
            KK = K - KSTEP + 1
            IF( KP.NE.KK ) THEN
*
*              Interchange rows and columns KK and KP in the leading
*              submatrix A(1:k,1:k)
*
               CALL ZSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )
               CALL ZSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ),
     $                     LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = A( K-1, K )
                  A( K-1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF
            END IF
*
*           Update the leading submatrix
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column k now holds
*
*              W(k) = U(k)*D(k)
*
*              where U(k) is the k-th column of U
*
*              Perform a rank-1 update of A(1:k-1,1:k-1) as
*
*              A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T
*
               R1 = CONE / A( K, K )
               CALL ZSYR( UPLO, K-1, -R1, A( 1, K ), 1, A, LDA )
*
*              Store U(k) in column k
*
               CALL ZSCAL( K-1, R1, A( 1, K ), 1 )
            ELSE
*
*              2-by-2 pivot block D(k): columns k and k-1 now hold
*
*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
*
*              where U(k) and U(k-1) are the k-th and (k-1)-th columns
*              of U
*
*              Perform a rank-2 update of A(1:k-2,1:k-2) as
*
*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T
*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T
*
               IF( K.GT.2 ) THEN
*
                  D12 = A( K-1, K )
                  D22 = A( K-1, K-1 ) / D12
                  D11 = A( K, K ) / D12
                  T = CONE / ( D11*D22-CONE )
                  D12 = T / D12
*
                  DO 30 J = K - 2, 1, -1
                     WKM1 = D12*( D11*A( J, K-1 )-A( J, K ) )
                     WK = D12*( D22*A( J, K )-A( J, K-1 ) )
                     DO 20 I = J, 1, -1
                        A( I, J ) = A( I, J ) - A( I, K )*WK -
     $                              A( I, K-1 )*WKM1
   20                CONTINUE
                     A( J, K ) = WK
                     A( J, K-1 ) = WKM1
   30             CONTINUE
*
               END IF
*
            END IF
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         END IF
*
*        Decrease K and return to the start of the main loop
*
         K = K - KSTEP
         GO TO 10
*
      ELSE
*
*        Factorize A as L*D*L**T using the lower triangle of A
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2
*
         K = 1
   40    CONTINUE
*
*        If K > N, exit from loop
*
         IF( K.GT.N )
     $      GO TO 70
         KSTEP = 1
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = CABS1( A( K, K ) )
*
*        IMAX is the row-index of the largest off-diagonal element in
*        column K, and COLMAX is its absolute value.
*        Determine both COLMAX and IMAX.
*
         IF( K.LT.N ) THEN
            IMAX = K + IZAMAX( N-K, A( K+1, K ), 1 )
            COLMAX = CABS1( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO .OR. DISNAN(ABSAKK) ) THEN
*
*           Column K is zero or underflow, or contains a NaN:
*           set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              no interchange, use 1-by-1 pivot block
*
               KP = K
            ELSE
*
*              JMAX is the column-index of the largest off-diagonal
*              element in row IMAX, and ROWMAX is its absolute value
*
               JMAX = K - 1 + IZAMAX( IMAX-K, A( IMAX, K ), LDA )
               ROWMAX = CABS1( A( IMAX, JMAX ) )
               IF( IMAX.LT.N ) THEN
                  JMAX = IMAX + IZAMAX( N-IMAX, A( IMAX+1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( A( JMAX, IMAX ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 no interchange, use 1-by-1 pivot block
*
                  KP = K
               ELSE IF( CABS1( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 interchange rows and columns K and IMAX, use 1-by-1
*                 pivot block
*
                  KP = IMAX
               ELSE
*
*                 interchange rows and columns K+1 and IMAX, use 2-by-2
*                 pivot block
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
            KK = K + KSTEP - 1
            IF( KP.NE.KK ) THEN
*
*              Interchange rows and columns KK and KP in the trailing
*              submatrix A(k:n,k:n)
*
               IF( KP.LT.N )
     $            CALL ZSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
               CALL ZSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ),
     $                     LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = A( K+1, K )
                  A( K+1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF
            END IF
*
*           Update the trailing submatrix
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column k now holds
*
*              W(k) = L(k)*D(k)
*
*              where L(k) is the k-th column of L
*
               IF( K.LT.N ) THEN
*
*                 Perform a rank-1 update of A(k+1:n,k+1:n) as
*
*                 A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T
*
                  R1 = CONE / A( K, K )
                  CALL ZSYR( UPLO, N-K, -R1, A( K+1, K ), 1,
     $                       A( K+1, K+1 ), LDA )
*
*                 Store L(k) in column K
*
                  CALL ZSCAL( N-K, R1, A( K+1, K ), 1 )
               END IF
            ELSE
*
*              2-by-2 pivot block D(k)
*
               IF( K.LT.N-1 ) THEN
*
*                 Perform a rank-2 update of A(k+2:n,k+2:n) as
*
*                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**T
*                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**T
*
*                 where L(k) and L(k+1) are the k-th and (k+1)-th
*                 columns of L
*
                  D21 = A( K+1, K )
                  D11 = A( K+1, K+1 ) / D21
                  D22 = A( K, K ) / D21
                  T = CONE / ( D11*D22-CONE )
                  D21 = T / D21
*
                  DO 60 J = K + 2, N
                     WK = D21*( D11*A( J, K )-A( J, K+1 ) )
                     WKP1 = D21*( D22*A( J, K+1 )-A( J, K ) )
                     DO 50 I = J, N
                        A( I, J ) = A( I, J ) - A( I, K )*WK -
     $                              A( I, K+1 )*WKP1
   50                CONTINUE
                     A( J, K ) = WK
                     A( J, K+1 ) = WKP1
   60             CONTINUE
               END IF
            END IF
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         END IF
*
*        Increase K and return to the start of the main loop
*
         K = K + KSTEP
         GO TO 40
*
      END IF
*
   70 CONTINUE
      RETURN
*
*     End of ZSYTF2
*
      END
*> \brief \b ZSYTRF
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZSYTRF + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytrf.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytrf.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytrf.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, LWORK, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         A( LDA, * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZSYTRF computes the factorization of a complex symmetric matrix A
*> using the Bunch-Kaufman diagonal pivoting method.  The form of the
*> factorization is
*>
*>    A = U*D*U**T  or  A = L*D*L**T
*>
*> where U (or L) is a product of permutation and unit upper (lower)
*> triangular matrices, and D is symmetric and block diagonal with
*> with 1-by-1 and 2-by-2 diagonal blocks.
*>
*> This is the blocked version of the algorithm, calling Level 3 BLAS.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangle of A is stored;
*>          = 'L':  Lower triangle of A is stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*>          N-by-N upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading N-by-N lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the strictly upper
*>          triangular part of A is not referenced.
*>
*>          On exit, the block diagonal matrix D and the multipliers used
*>          to obtain the factor U or L (see below for further details).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D.
*>          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
*>          interchanged and D(k,k) is a 1-by-1 diagonal block.
*>          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
*>          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
*>          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
*>          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
*>          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of WORK.  LWORK >=1.  For best performance
*>          LWORK >= N*NB, where NB is the block size returned by ILAENV.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
*>                has been completed, but the block diagonal matrix D is
*>                exactly singular, and division by zero will occur if it
*>                is used to solve a system of equations.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup complex16SYcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  If UPLO = 'U', then A = U*D*U**T, where
*>     U = P(n)*U(n)* ... *P(k)U(k)* ...,
*>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
*>  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
*>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
*>  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
*>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
*>
*>             (   I    v    0   )   k-s
*>     U(k) =  (   0    I    0   )   s
*>             (   0    0    I   )   n-k
*>                k-s   s   n-k
*>
*>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
*>  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
*>  and A(k,k), and v overwrites A(1:k-2,k-1:k).
*>
*>  If UPLO = 'L', then A = L*D*L**T, where
*>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
*>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
*>  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
*>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
*>  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
*>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
*>
*>             (   I    0     0   )  k-1
*>     L(k) =  (   0    I     0   )  s
*>             (   0    v     I   )  n-k-s+1
*>                k-1   s  n-k-s+1
*>
*>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
*>  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
*>  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            IINFO, IWS, J, K, KB, LDWORK, LWKOPT, NB, NBMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLASYF, ZSYTF2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size
*
         NB = ILAENV( 1, 'ZSYTRF', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSYTRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = LDWORK*NB
         IF( LWORK.LT.IWS ) THEN
            NB = MAX( LWORK / LDWORK, 1 )
            NBMIN = MAX( 2, ILAENV( 2, 'ZSYTRF', UPLO, N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = 1
      END IF
      IF( NB.LT.NBMIN )
     $   NB = N
*
      IF( UPPER ) THEN
*
*        Factorize A as U*D*U**T using the upper triangle of A
*
*        K is the main loop index, decreasing from N to 1 in steps of
*        KB, where KB is the number of columns factorized by ZLASYF;
*        KB is either NB or NB-1, or K for the last block
*
         K = N
   10    CONTINUE
*
*        If K < 1, exit from loop
*
         IF( K.LT.1 )
     $      GO TO 40
*
         IF( K.GT.NB ) THEN
*
*           Factorize columns k-kb+1:k of A and use blocked code to
*           update columns 1:k-kb
*
            CALL ZLASYF( UPLO, K, NB, KB, A, LDA, IPIV, WORK, N, IINFO )
         ELSE
*
*           Use unblocked code to factorize columns 1:k of A
*
            CALL ZSYTF2( UPLO, K, A, LDA, IPIV, IINFO )
            KB = K
         END IF
*
*        Set INFO on the first occurrence of a zero pivot
*
         IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO
*
*        Decrease K and return to the start of the main loop
*
         K = K - KB
         GO TO 10
*
      ELSE
*
*        Factorize A as L*D*L**T using the lower triangle of A
*
*        K is the main loop index, increasing from 1 to N in steps of
*        KB, where KB is the number of columns factorized by ZLASYF;
*        KB is either NB or NB-1, or N-K+1 for the last block
*
         K = 1
   20    CONTINUE
*
*        If K > N, exit from loop
*
         IF( K.GT.N )
     $      GO TO 40
*
         IF( K.LE.N-NB ) THEN
*
*           Factorize columns k:k+kb-1 of A and use blocked code to
*           update columns k+kb:n
*
            CALL ZLASYF( UPLO, N-K+1, NB, KB, A( K, K ), LDA, IPIV( K ),
     $                   WORK, N, IINFO )
         ELSE
*
*           Use unblocked code to factorize columns k:n of A
*
            CALL ZSYTF2( UPLO, N-K+1, A( K, K ), LDA, IPIV( K ), IINFO )
            KB = N - K + 1
         END IF
*
*        Set INFO on the first occurrence of a zero pivot
*
         IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO + K - 1
*
*        Adjust IPIV
*
         DO 30 J = K, K + KB - 1
            IF( IPIV( J ).GT.0 ) THEN
               IPIV( J ) = IPIV( J ) + K - 1
            ELSE
               IPIV( J ) = IPIV( J ) - K + 1
            END IF
   30    CONTINUE
*
*        Increase K and return to the start of the main loop
*
         K = K + KB
         GO TO 20
*
      END IF
*
   40 CONTINUE
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of ZSYTRF
*
      END
*> \brief \b ZSYTRI
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZSYTRI + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytri.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytri.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytri.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         A( LDA, * ), WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZSYTRI computes the inverse of a complex symmetric indefinite matrix
*> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
*> ZSYTRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are stored
*>          as an upper or lower triangular matrix.
*>          = 'U':  Upper triangular, form is A = U*D*U**T;
*>          = 'L':  Lower triangular, form is A = L*D*L**T.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the block diagonal matrix D and the multipliers
*>          used to obtain the factor U or L as computed by ZSYTRF.
*>
*>          On exit, if INFO = 0, the (symmetric) inverse of the original
*>          matrix.  If UPLO = 'U', the upper triangular part of the
*>          inverse is formed and the part of A below the diagonal is not
*>          referenced; if UPLO = 'L' the lower triangular part of the
*>          inverse is formed and the part of A above the diagonal is
*>          not referenced.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D
*>          as determined by ZSYTRF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (2*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
*>               inverse could not be computed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup complex16SYcomputational
*
*  =====================================================================
      SUBROUTINE ZSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K, KP, KSTEP
      COMPLEX*16         AK, AKKP1, AKP1, D, T, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTU
      EXTERNAL           LSAME, ZDOTU
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZCOPY, ZSWAP, ZSYMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSYTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check that the diagonal matrix D is nonsingular.
*
      IF( UPPER ) THEN
*
*        Upper triangular storage: examine D from bottom to top
*
         DO 10 INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
      ELSE
*
*        Lower triangular storage: examine D from top to bottom.
*
         DO 20 INFO = 1, N
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   20    CONTINUE
      END IF
      INFO = 0
*
      IF( UPPER ) THEN
*
*        Compute inv(A) from the factorization A = U*D*U**T.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = 1
   30    CONTINUE
*
*        If K > N, exit from loop.
*
         IF( K.GT.N )
     $      GO TO 40
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Invert the diagonal block.
*
            A( K, K ) = ONE / A( K, K )
*
*           Compute column K of the inverse.
*
            IF( K.GT.1 ) THEN
               CALL ZCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL ZSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO,
     $                     A( 1, K ), 1 )
               A( K, K ) = A( K, K ) - ZDOTU( K-1, WORK, 1, A( 1, K ),
     $                     1 )
            END IF
            KSTEP = 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Invert the diagonal block.
*
            T = A( K, K+1 )
            AK = A( K, K ) / T
            AKP1 = A( K+1, K+1 ) / T
            AKKP1 = A( K, K+1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K, K ) = AKP1 / D
            A( K+1, K+1 ) = AK / D
            A( K, K+1 ) = -AKKP1 / D
*
*           Compute columns K and K+1 of the inverse.
*
            IF( K.GT.1 ) THEN
               CALL ZCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL ZSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO,
     $                     A( 1, K ), 1 )
               A( K, K ) = A( K, K ) - ZDOTU( K-1, WORK, 1, A( 1, K ),
     $                     1 )
               A( K, K+1 ) = A( K, K+1 ) -
     $                       ZDOTU( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 )
               CALL ZCOPY( K-1, A( 1, K+1 ), 1, WORK, 1 )
               CALL ZSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO,
     $                     A( 1, K+1 ), 1 )
               A( K+1, K+1 ) = A( K+1, K+1 ) -
     $                         ZDOTU( K-1, WORK, 1, A( 1, K+1 ), 1 )
            END IF
            KSTEP = 2
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           Interchange rows and columns K and KP in the leading
*           submatrix A(1:k+1,1:k+1)
*
            CALL ZSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 )
            CALL ZSWAP( K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = A( K, K+1 )
               A( K, K+1 ) = A( KP, K+1 )
               A( KP, K+1 ) = TEMP
            END IF
         END IF
*
         K = K + KSTEP
         GO TO 30
   40    CONTINUE
*
      ELSE
*
*        Compute inv(A) from the factorization A = L*D*L**T.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = N
   50    CONTINUE
*
*        If K < 1, exit from loop.
*
         IF( K.LT.1 )
     $      GO TO 60
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Invert the diagonal block.
*
            A( K, K ) = ONE / A( K, K )
*
*           Compute column K of the inverse.
*
            IF( K.LT.N ) THEN
               CALL ZCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL ZSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1,
     $                     ZERO, A( K+1, K ), 1 )
               A( K, K ) = A( K, K ) - ZDOTU( N-K, WORK, 1, A( K+1, K ),
     $                     1 )
            END IF
            KSTEP = 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Invert the diagonal block.
*
            T = A( K, K-1 )
            AK = A( K-1, K-1 ) / T
            AKP1 = A( K, K ) / T
            AKKP1 = A( K, K-1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K-1, K-1 ) = AKP1 / D
            A( K, K ) = AK / D
            A( K, K-1 ) = -AKKP1 / D
*
*           Compute columns K-1 and K of the inverse.
*
            IF( K.LT.N ) THEN
               CALL ZCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL ZSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1,
     $                     ZERO, A( K+1, K ), 1 )
               A( K, K ) = A( K, K ) - ZDOTU( N-K, WORK, 1, A( K+1, K ),
     $                     1 )
               A( K, K-1 ) = A( K, K-1 ) -
     $                       ZDOTU( N-K, A( K+1, K ), 1, A( K+1, K-1 ),
     $                       1 )
               CALL ZCOPY( N-K, A( K+1, K-1 ), 1, WORK, 1 )
               CALL ZSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1,
     $                     ZERO, A( K+1, K-1 ), 1 )
               A( K-1, K-1 ) = A( K-1, K-1 ) -
     $                         ZDOTU( N-K, WORK, 1, A( K+1, K-1 ), 1 )
            END IF
            KSTEP = 2
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           Interchange rows and columns K and KP in the trailing
*           submatrix A(k-1:n,k-1:n)
*
            IF( KP.LT.N )
     $         CALL ZSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 )
            CALL ZSWAP( KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = A( K, K-1 )
               A( K, K-1 ) = A( KP, K-1 )
               A( KP, K-1 ) = TEMP
            END IF
         END IF
*
         K = K - KSTEP
         GO TO 50
   60    CONTINUE
      END IF
*
      RETURN
*
*     End of ZSYTRI
*
      END
*> \brief \b ZTRTI2 computes the inverse of a triangular matrix (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZTRTI2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrti2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrti2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrti2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTRTI2( UPLO, DIAG, N, A, LDA, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          DIAG, UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTRTI2 computes the inverse of a complex upper or lower triangular
*> matrix.
*>
*> This is the Level 2 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the matrix A is upper or lower triangular.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>          Specifies whether or not the matrix A is unit triangular.
*>          = 'N':  Non-unit triangular
*>          = 'U':  Unit triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the triangular matrix A.  If UPLO = 'U', the
*>          leading n by n upper triangular part of the array A contains
*>          the upper triangular matrix, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading n by n lower triangular part of the array A contains
*>          the lower triangular matrix, and the strictly upper
*>          triangular part of A is not referenced.  If DIAG = 'U', the
*>          diagonal elements of A are also not referenced and are
*>          assumed to be 1.
*>
*>          On exit, the (triangular) inverse of the original matrix, in
*>          the same storage format.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -k, the k-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup complex16OTHERcomputational
*
*  =====================================================================
      SUBROUTINE ZTRTI2( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J
      COMPLEX*16         AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZSCAL, ZTRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTRTI2', -INFO )
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Compute inverse of upper triangular matrix.
*
         DO 10 J = 1, N
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
*
*           Compute elements 1:j-1 of j-th column.
*
            CALL ZTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA,
     $                  A( 1, J ), 1 )
            CALL ZSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE
      ELSE
*
*        Compute inverse of lower triangular matrix.
*
         DO 20 J = N, 1, -1
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
            IF( J.LT.N ) THEN
*
*              Compute elements j+1:n of j-th column.
*
               CALL ZTRMV( 'Lower', 'No transpose', DIAG, N-J,
     $                     A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL ZSCAL( N-J, AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of ZTRTI2
*
      END
*> \brief \b ZTRTRI
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZTRTRI + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrtri.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrtri.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrtri.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTRTRI( UPLO, DIAG, N, A, LDA, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          DIAG, UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTRTRI computes the inverse of a complex upper or lower triangular
*> matrix A.
*>
*> This is the Level 3 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  A is upper triangular;
*>          = 'L':  A is lower triangular.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>          = 'N':  A is non-unit triangular;
*>          = 'U':  A is unit triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the triangular matrix A.  If UPLO = 'U', the
*>          leading N-by-N upper triangular part of the array A contains
*>          the upper triangular matrix, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading N-by-N lower triangular part of the array A contains
*>          the lower triangular matrix, and the strictly upper
*>          triangular part of A is not referenced.  If DIAG = 'U', the
*>          diagonal elements of A are also not referenced and are
*>          assumed to be 1.
*>          On exit, the (triangular) inverse of the original matrix, in
*>          the same storage format.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
*>               matrix is singular and its inverse can not be computed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup complex16OTHERcomputational
*
*  =====================================================================
      SUBROUTINE ZTRTRI( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J, JB, NB, NN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZTRMM, ZTRSM, ZTRTI2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTRTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check for singularity if non-unit.
*
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
         INFO = 0
      END IF
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'ZTRTRI', UPLO // DIAG, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL ZTRTI2( UPLO, DIAG, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code
*
         IF( UPPER ) THEN
*
*           Compute inverse of upper triangular matrix
*
            DO 20 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
*
*              Compute rows 1:j-1 of current block column
*
               CALL ZTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, ONE, A, LDA, A( 1, J ), LDA )
               CALL ZTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
*
*              Compute inverse of current diagonal block
*
               CALL ZTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
   20       CONTINUE
         ELSE
*
*           Compute inverse of lower triangular matrix
*
            NN = ( ( N-1 ) / NB )*NB + 1
            DO 30 J = NN, 1, -NB
               JB = MIN( NB, N-J+1 )
               IF( J+JB.LE.N ) THEN
*
*                 Compute rows j+jb:n of current block column
*
                  CALL ZTRMM( 'Left', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,
     $                        A( J+JB, J ), LDA )
                  CALL ZTRSM( 'Right', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, -ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
*
*              Compute inverse of current diagonal block
*
               CALL ZTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
   30       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZTRTRI
*
      END
