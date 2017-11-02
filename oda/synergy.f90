! ===================================================================================
PROGRAM Synergy

! ===================================================================================
! Program to investigate the potential synergistic interaction of observations
! Ross Bannister
! National Centre for Earth Observation
! January 2011
! FORTRAN 90
!
! On Unix
! ------------------
! setup studio12
! To compile on Unix
! f90 synergy.f90 -o synergy.out
! To run on Unix
! synergy.out
! ------------------
! Input files: ObsNetwork, DiagsPos
! Output file: B_Pa_Errs.dat


! ===================================================================================

  ! Declare variables

  IMPLICIT NONE

  INTEGER,          PARAMETER :: MaxNoOD = 10         ! Max No. of obs/diags
  INTEGER                     :: p                    ! No. of observations
  INTEGER,          PARAMETER :: Nx = 20              ! No. of x-points
  INTEGER,          PARAMETER :: Ny = 20              ! No. of y-points
  INTEGER,          PARAMETER :: n = 3*Nx*Ny          ! Size of model space
  DOUBLE PRECISION, PARAMETER :: d = 100000.0         ! Grid length
  DOUBLE PRECISION, PARAMETER :: Density = 0.5        ! Density of air
  DOUBLE PRECISION, PARAMETER :: f = 0.0001           ! Coriolis parameter
  DOUBLE PRECISION, PARAMETER :: L = 500000.0         ! p-p correlation lengthscale
  DOUBLE PRECISION, PARAMETER :: sigma_p = 1.0        ! b/g pressure errors
  INTEGER                     :: ObsType(1:MaxNoOD)   ! Obs type 1=p, 2=u, 3=v
  INTEGER                     :: Obs_xpos(1:MaxNoOD)  ! x positions of obs
  INTEGER                     :: Obs_ypos(1:MaxNoOD)  ! y positions of obs
  DOUBLE PRECISION            :: ObsErrs(1:MaxNoOD)   ! Observation std err
  DOUBLE PRECISION            :: B(1:n, 1:n)          ! B-matrix
  DOUBLE PRECISION            :: H(1:MaxNoOD,1:n)     ! Obs operator
  DOUBLE PRECISION            :: Pa(1:n, 1:n)         ! Anal error covariance matrix
  INTEGER                     :: NDiagPos             ! No. of diagnostics positions
  LOGICAL,          PARAMETER :: Multivariate =.true. ! True for multivariate B
  LOGICAL,          PARAMETER :: PaDiagOnly = .true.  ! True to calc diag of Pa only



  ! 1. Read-in data about the observation network
  ! ---------------------------------------------
  PRINT *, 'Reading-in the observations'
  CALL ReadObsInfo ( MaxNoOD, p, Nx, Ny, ObsType(1:MaxNoOD), Obs_xpos(1:MaxNoOD), &
                     Obs_ypos(1:MaxNoOD), ObsErrs(1:MaxNoOD), 'ObsNetwork' )

  ! 2. Compute the B-matrix
  ! -----------------------
  PRINT *, 'Computing the B-matrix'
  CALL CalcB ( Nx, Ny, n, B(1:n,1:n), d, Density, f, L, sigma_p, Multivariate )

  ! 3. Calculate the H-matrix
  ! -------------------------
  PRINT *, 'Computing the H-matrix'
  CALL CalcH (p, Nx, Ny, n, H(1:p,1:n), ObsType(1:p), Obs_xpos(1:p), Obs_ypos(1:p) )

  ! 4. Calculate the analysis error covariance matrix
  ! -------------------------------------------------
  PRINT *, 'Calculating the analysis error covariance matrix'
  CALL CalcAnalError ( p, Nx, Ny, n, H(1:p,1:n), B(1:n,1:n), ObsErrs(1:p), &
                       Pa(1:n,1:n), PaDiagOnly )

  ! 5. Read-in diagnostics positions and output diagnostics at observation locations
  ! --------------------------------------------------------------------------------
  PRINT *, 'Outputting B and Pa diagnostics at observation locations'
  CALL B_Pa_Diags ( MaxNoOD, Nx, Ny, n, B(1:n,1:n), Pa(1:n,1:n), &
                    'DiagsPos','B_Pa_Errs.dat' )
 
  PRINT *, 'Program completed'

CONTAINS

! ===================================================================================
SUBROUTINE CalcB ( Nx,      &      ! No. of x-points
                   Ny,      &      ! No. of y-points
                   n,       &      ! Size of model space
                   B,       &      ! B-matrix
                   d,       &      ! Grid length
                   Density, &      ! Density
                   f,       &      ! Coriolis
                   L,       &      ! Lengthscale (grid points)
                   sigma_p, &      ! Background error of pressure
                   Multivariate )  ! .true. or .false.

! Subroutine to determine the B-matrix according to geostrophy in either a uni- or
! multi-variate fashion
! Ross Bannister
! National Centre for Earth Observation
! January 2011
! FORTRAN 90

  IMPLICIT NONE

  ! Subroutine arguments
  INTEGER,          INTENT(IN)    :: Nx
  INTEGER,          INTENT(IN)    :: Ny
  INTEGER,          INTENT(IN)    :: n
  DOUBLE PRECISION, INTENT(INOUT) :: B(1:n, 1:n)
  DOUBLE PRECISION, INTENT(IN)    :: d
  DOUBLE PRECISION, INTENT(IN)    :: Density
  DOUBLE PRECISION, INTENT(IN)    :: f
  DOUBLE PRECISION, INTENT(IN)    :: L
  DOUBLE PRECISION, INTENT(IN)    :: sigma_p
  LOGICAL,          INTENT(IN)    :: Multivariate

  ! Local variables
  INTEGER                         :: i, j, ip, jp, npoints
  INTEGER                         :: element_p, element_u, element_v
  INTEGER                         :: element_pp, element_up, element_vp
  DOUBLE PRECISION                :: dis2, mu, dx, dy, dx2, dy2, minushalfinvL2
  DOUBLE PRECISION                :: sigma_p2, Density2, f2, invL2

  npoints          = Nx*Ny
  invL2            = 1.0 / (L*L)
  minushalfinvL2   = -0.5 * invL2
  sigma_p2         = sigma_p * sigma_p
  Density2         = Density * Density
  f2               = f * f

  DO i = 1, Nx
    DO j = 1, Ny
      element_p = (i-1)*Ny + j
      element_u = element_p + npoints
      element_v = element_u + npoints
      DO ip = 1, Nx
        DO jp = 1, Ny
          element_pp = (ip-1)*Ny + jp
          element_up = element_pp + npoints
          element_vp = element_up + npoints

          ! Distance squared between (i,j) and (ip,jp)
          dx   = DBLE(i-ip) * d
          dx2  = dx*dx
          dy   = DBLE(j-jp) * d
          dy2  = dy*dy
          dis2 = dx2 + dy2

          ! Pressure correlation value
          mu = EXP(dis2*minushalfinvL2)

          ! p-p covariances
          B(element_p,element_pp) = sigma_p2 * mu

          ! p-u covariances
          IF (Multivariate) THEN
            B(element_p,element_up) = -1.0 * &
                                      invL2 * sigma_p2 * mu * dy / (Density * f)
          ELSE
            B(element_p,element_up) = 0.0
          END IF

          ! p-v covariances
          IF (Multivariate) THEN
            B(element_p,element_vp) = invL2 * sigma_p2 * mu * dx / (Density * f)
          ELSE
            B(element_p,element_vp) = 0.0
          END IF

          ! u-p covariances
          IF (Multivariate) THEN
            B(element_u,element_pp) = invL2 * sigma_p2 * mu * dy / (Density * f)
          ELSE
            B(element_u,element_pp) = 0.0
          END IF

          ! u-u covariances
          B(element_u,element_up) = invL2 * sigma_p2 * mu * (1.0 - invL2 * dy2) / &
                                    (Density2 * f2)

          ! u-v covariances
          IF (Multivariate) THEN
            B(element_u,element_vp) = invL2 * invL2 * sigma_p2 * mu * dx * dy / &
                                      (Density2 * f2)
          ELSE
            B(element_u,element_vp) = 0.0
          END IF

          ! v-p covariances
          IF (Multivariate) THEN
            B(element_v,element_pp) = -1.0 * &
                                      invL2 * sigma_p2 * mu * dx / (Density * f)
          ELSE
            B(element_v,element_pp) = 0.0
          END IF

          ! v-u covariances
          IF (Multivariate) THEN
            B(element_v,element_up) = invL2 * invL2 * sigma_p2 * mu * dx * dy / &
                                      (Density2 * f2)
          ELSE
            B(element_v,element_up) = 0.0
          END IF

          ! v-v covariances
          B(element_v,element_vp) = invL2 * sigma_p2 * mu * (1.0 - invL2 * dx2) / &
                                    (Density2 * f2)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE CalcB


! ===================================================================================
SUBROUTINE CalcH ( p,        &      ! No. of observations
                   Nx,       &      ! No. of x-points
                   Ny,       &      ! No. of y-points
                   n,        &      ! Size of model space
                   H,        &      ! H-matrix
                   ObsType,  &      ! Obs type 1=p, 2=u, 3=v
                   Obs_xpos, &      ! x positions of obs
                   Obs_ypos )       ! y positions of obs

! Subroutine to determine the H-matrix according to the observations specified
! Ross Bannister
! National Centre for Earth Observation
! January 2011
! FORTRAN 90

  IMPLICIT NONE

  ! Subroutine arguments
  INTEGER,          INTENT(IN)    :: p
  INTEGER,          INTENT(IN)    :: Nx
  INTEGER,          INTENT(IN)    :: Ny
  INTEGER,          INTENT(IN)    :: n
  DOUBLE PRECISION, INTENT(INOUT) :: H(1:p,1:n)
  INTEGER,          INTENT(IN)    :: ObsType(1:p)
  INTEGER,          INTENT(IN)    :: Obs_xpos(1:p)
  INTEGER,          INTENT(IN)    :: Obs_ypos(1:p)

  ! Local variables
  INTEGER                :: k, element, npoints

  npoints    = Nx*Ny
  H(1:p,1:n) = 0.0

  DO k = 1, p
    ! Find the element in the model space of this observation
    element      = (Obs_xpos(k)-1)*Ny + Obs_ypos(k) + npoints * (ObsType(k)-1)
    H(k,element) = 1.0
  END DO
END SUBROUTINE CalcH


! ===================================================================================
SUBROUTINE CalcAnalError ( p,        &      ! No. of observations
                           Nx,       &      ! No. of x-points
                           Ny,       &      ! No. of y-points
                           n,        &      ! Size of model space
                           H,        &      ! H-matrix
                           B,        &      ! B-matrix
                           ObsErrs,  &      ! Observation std err
                           Pa,       &      ! Pa-matrix
                           DiagOnly )       ! true to calc diagonal elements only

! Subroutine to determine the analysis error covariance matrix for a given B, H, R
! Ross Bannister
! National Centre for Earth Observation
! January 2011
! FORTRAN 90

  IMPLICIT NONE

  ! Subroutine arguments
  INTEGER,          INTENT(IN)    :: p
  INTEGER,          INTENT(IN)    :: Nx
  INTEGER,          INTENT(IN)    :: Ny
  INTEGER,          INTENT(IN)    :: n
  DOUBLE PRECISION, INTENT(IN)    :: H(1:p,1:n)
  DOUBLE PRECISION, INTENT(IN)    :: B(1:n,1:n)
  DOUBLE PRECISION, INTENT(IN)    :: ObsErrs(1:p)
  DOUBLE PRECISION, INTENT(INOUT) :: Pa(1:n,1:n)
  LOGICAL,          INTENT(IN)    :: DiagOnly

  ! Local variables
  INTEGER                         :: i, j, u, v, from, to, npoints
  DOUBLE PRECISION                :: HBHT_R(1:p,1:p)
  DOUBLE PRECISION                :: HBHT_Rinv(1:p,1:p)
  DOUBLE PRECISION                :: InvCheck(1:p,1:p)
  DOUBLE PRECISION                :: K(1:n,1:p)
  DOUBLE PRECISION                :: I_KH(1:n,1:n)
  DOUBLE PRECISION                :: total

  npoints    = Nx*Ny

  ! Use formulae
  ! Pa = (I - KH) B
  ! K = BH^T (R + HBH^T)^-1

  ! 1. Calculate HBH^T + R
  ! ----------------------
  PRINT *, '  Calculating HBH^T + R'
  DO i = 1, p
    DO j = 1, p
      ! Calculate HBH^T
      total= 0.0
      DO u = 1, n
        DO v = 1, n
          total = total + H(i,u) * H(j,v) * B(u,v)
        END DO
      END DO
      HBHT_R(i,j) = total
    END DO
    ! Add on obs errors (diagonal elements only)
    HBHT_R(i,i) = HBHT_R(i,i) + ObsErrs(i) * ObsErrs(i)
  END DO

  ! 2. Calculate (HBH^T + R)^-1
  ! ---------------------------
  PRINT *, '  Inverting this'
  CALL MatrixInverse (p, HBHT_R(1:p,1:p), HBHT_Rinv(1:p,1:p))

  ! 2.a - Check that this is the inverse
  ! ------------------------------------
  DO i = 1, p
    DO j = 1, p
      total= 0.0
      DO u = 1, p
        total = total + HBHT_Rinv(i,u) * HBHT_R(u,j)
      END DO
      InvCheck(i,j) = total
    END DO
  END DO
  PRINT *,'Report from subroutine CalcAnalError - inverse check - this should be I'
  DO i = 1, p
    PRINT '(20F8.4)', (InvCheck(i,j), j=1,p)
  END DO

  ! 3. Calculate K
  ! --------------
  PRINT *, '  Calculating K'
  DO i = 1, n
    DO j = 1, p
      total= 0.0
      DO u = 1, n
        DO v = 1, p
          total = total + B(i,u) * H(v,u) * HBHT_Rinv(v,j)
        END DO
      END DO
      K(i,j) = total
    END DO
  END DO

  ! 4. Calculate I-KH
  ! -----------------
  PRINT *, '  Calculating I-KH'
  DO i = 1, n
    DO j = 1, n
      ! Calculate -KH
      total= 0.0
      DO u = 1, p
        total = total + K(i,u) * H(u,j)
      END DO
      I_KH(i,j) = -1.0 * total
    END DO
    ! Calculate I_KH
    I_KH(i,i) = I_KH(i,i) + 1.0
  END DO

  ! 5. Calculate Pa
  ! ---------------
  PRINT *, '  Calculating Pa'
  DO i = 1, n
    IF (DiagOnly) THEN
      from = i
      to   = i
    ELSE
      from = 1
      to   = n
    END IF
    IF (.NOT.DiagOnly) PRINT *, '  Dealing with row ', i, 'of ', n
    DO j = from, to
      total= 0.0
      DO u = 1, n
        total = total + I_KH(i,u) * B(u,j)
      END DO
      Pa(i,j) = total
    END DO
  END DO

END SUBROUTINE CalcAnalError



! ===================================================================================
SUBROUTINE ReadObsInfo ( Max,      &      ! Maximum No. of obs
                         p,        &      ! No. of observations
                         Nx,       &      ! No. of x-points
                         Ny,       &      ! No. of y-points
                         ObsType,  &      ! Obs type 1=p, 2=u, 3=v
                         Obs_xpos, &      ! x positions of obs
                         Obs_ypos, &      ! y positions of obs
                         ObsErrs,  &      ! Observation std err
                         Filename )       ! Input filename

! Subroutine to read-in from file observation information
! Ross Bannister
! National Centre for Earth Observation
! January 2011
! FORTRAN 90

  IMPLICIT NONE

  ! Subroutine arguments
  INTEGER,           INTENT(IN)    :: Max
  INTEGER,           INTENT(INOUT) :: p
  INTEGER,           INTENT(IN)    :: Nx
  INTEGER,           INTENT(IN)    :: Ny
  INTEGER,           INTENT(INOUT) :: ObsType(1:Max)
  INTEGER,           INTENT(INOUT) :: Obs_xpos(1:Max)
  INTEGER,           INTENT(INOUT) :: Obs_ypos(1:Max)
  DOUBLE PRECISION,  INTENT(INOUT) :: ObsErrs(1:Max)
  CHARACTER (LEN=*), INTENT(IN)    :: Filename

  ! Local variables
  INTEGER                          :: k

  ! Read-in the data
  OPEN (12, file=Filename)
  READ (12,*) p
  IF (p > Max) THEN
    PRINT *,'Too many observations'
    STOP
  END IF
  PRINT *, 'Reading ', p, ' observations'
  PRINT *, '--------------------------------------------------------------'
  PRINT *, 'Observations'
  DO k = 1, p
    READ (12,*) ObsType(k), Obs_xpos(k), Obs_ypos(k), ObsErrs(k)
    PRINT '(4I4,F8.4)', k, ObsType(k), Obs_xpos(k), Obs_ypos(k), ObsErrs(k)
  END DO
  PRINT *, '--------------------------------------------------------------'
  CLOSE(12)

  ! Do some checks
  DO k = 1, p
    IF ((ObsType(k) < 1) .OR. (ObsType(k) > 3)) THEN
      PRINT *, 'Problem with obs ', k, '(type out of range)'
      STOP
    END IF
    IF ((Obs_xpos(k) < 1) .OR. (Obs_xpos(k) > Nx)) THEN
      PRINT *, 'Problem with obs ', k, '(x position out of range)'
      STOP
    END IF
    IF ((Obs_ypos(k) < 1) .OR. (Obs_ypos(k) > Ny)) THEN
      PRINT *, 'Problem with obs ', k, '(y position out of range)'
      STOP
    END IF
  END DO

END SUBROUTINE ReadObsInfo



! ===================================================================================
SUBROUTINE MatrixInverse ( n,        &      ! Size of matrix
                           A,        &      ! Matrix to be inverted
                           Ainv )           ! Inverse

! Subroutine to calculate the inverse of a matrix
! (very inefficient method, so suitable only for very small n)
! Ross Bannister
! National Centre for Earth Observation
! January 2011
! FORTRAN 90

  IMPLICIT NONE

  ! Subroutine arguments
  INTEGER,          INTENT(IN)    :: n
  DOUBLE PRECISION, INTENT(IN)    :: A(1:n,1:n)
  DOUBLE PRECISION, INTENT(INOUT) :: Ainv(1:n,1:n)

  ! Local variables
  INTEGER                         :: i, j, u
  DOUBLE PRECISION                :: total
  DOUBLE PRECISION                :: Evals(1:n,1:n)
  DOUBLE PRECISION                :: Evecs(1:n,1:n)  !Note - rows not columns
  DOUBLE PRECISION                :: lim

  IF (n == 1) THEN
    Ainv(1,1) = 1.0 / A(1,1)

  ELSE

    lim = 0.0001

    ! Diagonalize the matrix A
    CALL Jacobi (A(1:n,1:n), Evals(1:n,1:n), Evecs(1:n,1:n), n, lim)

    ! Construct the inverse
    DO i = 1, n
      DO j = 1, n
        total = 0.0
        DO u = 1, n
          total = total + Evecs(u,i) * Evecs(u,j) / Evals(u,u)
        END DO
        Ainv(i,j) = total
      END DO
    END DO
  END IF
END SUBROUTINE MatrixInverse



! ===================================================================================
  SUBROUTINE Jacobi (A,B,T,N,lim)
!     To diagonalize a square symmetric matrix by a sequence of Jacobi rotations
!     Ross Bannister, January 2001
!     A (input) matrix in original representation
!     B (output) matrix in diagonal representation
!     T (output) matrix of eigenvectors (rows)
!     N (input) the order of the matrix
!     lim (input) convergence criteria, e.g. lim=0.01

      implicit none

      DOUBLE PRECISION, INTENT(IN)    :: A(1:N,1:N)
      DOUBLE PRECISION, INTENT(INOUT) :: B(1:N,1:N)
      DOUBLE PRECISION, INTENT(INOUT) :: T(1:N,1:N)
      INTEGER,          INTENT(IN)    :: N
      DOUBLE PRECISION, INTENT(IN)    :: lim


      INTEGER                         :: i,j,p,q,k,lp1,lp2,sweeps
      DOUBLE PRECISION                :: off,diag,theta,t1,t2,tt
      DOUBLE PRECISION                :: cc,ss,sc,Bip,Bqj,Tpj,Bpp,Bqq,dis,s,c,norm


!     The matrix B is initially equal to A
      do j=1,N
        do i=1,N
          B(i,j)=A(i,j)
        enddo
      enddo

!     The matrix T is initially the identity matrix
      do j=1,N
        do i=1,N
          if (i.eq.j) then
            T(i,j)=1.0
            else
            T(i,j)=0.0
          endif
        enddo
      enddo

!      print 101,'------------------------------------------------------'
!      print *,'The original matrix and transformation N=',N
!      do lp1=1,N
!        print 100,(B(lp1,lp2),lp2=1,N),(T(lp1,lp2),lp2=1,N)
!      enddo

!     Calculate the sum of the square of off-diagonal elements of B
      off=0.0
      do j=2,N
        do i=1,j-1
          off=off+B(i,j)*B(i,j)
        enddo
      enddo
      off=2.0*off

!     Calculate the sum of the square of diagonal elements of B
      diag=0.0
      do i=1,N
        diag=diag+B(i,i)*B(i,i)
      enddo
!      print 102,'off/diag =',off/diag

!     The number of sweeps completed
      sweeps=0
10    continue

!     Loop round each upper-right element (a 'sweep')
      do p=2,N
        do q=1,p-1
!          print*
!          print*,'Working with position ',q,p
          if (abs(B(q,p)).gt.(0.0001)) then
!           In order to find the new representation of the matrix,
!           work out some details regarding the rotation, U to eliminate
!           element q,p of B (p>q).
            theta=(B(p,p)-B(q,q))/(2.0*B(q,p))
            dis=sqrt(theta*theta+1)
!           Two possible values of t
            t1=-theta+dis
            t2=-theta-dis
!           Choose the smaller value
            if (abs(t1).lt.abs(t2)) then
              tt=t1
              else
              tt=t2
            endif
!           Find the values of s and c
            c=1/sqrt(tt*tt+1)
            s=tt*c
!            print*,theta,tt,c,s
!           Modify the operator UBU(trans) (piecewise)
!           Update only upper-right part of matrix
            cc=c*c
            ss=s*s
            sc=s*c
!           Columns p and q:
            do i=1,p-1
              if (i.ne.q) then
                Bip=B(i,p)
                B(i,p)=c*Bip+s*B(i,q)
                if (i.le.q) B(i,q)=c*B(i,q)-s*Bip
              endif
            enddo
!           Rows p and q:
            do j=q+1,N
              if (j.ne.p) then
                Bqj=B(q,j)
                B(q,j)=c*Bqj-s*B(p,j)
                if (j.gt.p) B(p,j)=c*B(p,j)+s*Bqj
              endif
            enddo
!           Mixed elements
!            B(q,p)=(cc-ss)*B(q,p)+sc*(B(q,q)-B(p,p))
!            print*,'This should be zero ',B(q,p)
            B(q,p)=0.0

!           New diagonal elements
            Bpp=B(p,p)
            Bqq=B(q,q)
            B(p,p)=cc*Bpp+ss*B(q,q)+2.0*sc*B(p,q)
            B(q,q)=cc*B(q,q)+ss*Bpp-2.0*sc*B(p,q)

!           Modified value of sum of square of off-diagonal elements...
            off=off-2.0*B(p,q)*B(p,q)
!           ...and diagonal elements
            diag=diag-Bpp*Bpp-Bqq*Bqq+B(p,p)*B(p,p)+B(q,q)*B(q,q)

!           Update lower-left of matrix
            do k=1,N
              if (k.lt.p) then
                B(p,k)=B(k,p)
                if (k.lt.q) B(q,k)=B(k,q)
              endif
              if (k.gt.q) then
                B(k,q)=B(q,k)
                if (k.gt.p) B(k,p)=B(p,k)
              endif
            enddo

!           Modify the transformation T(new)=U T(old)
!           Only rows p and q are affected
            do j=1,N
              Tpj=T(p,j)
              T(p,j)=c*Tpj+s*T(q,j)
              T(q,j)=c*T(q,j)-s*Tpj
            enddo
!            print 101,'------------------------------------------------'
!            print 101,'The new matrix and new transformation'
!            do lp1=1,N
!              print 100,(B(lp1,lp2),lp2=1,N),(T(lp1,lp2),lp2=1,N)
!            enddo
!            print 102,'off/diag =',off/diag
            else
!            print*,'This element is small anyway!'
          endif
        enddo
      enddo

      sweeps=sweeps+1

      if (((off/diag).gt.lim).or.(sweeps.lt.3)) goto 10

!     Normalize the eigenvectors
      do i=1,N
        norm=0.0
        do j=1,N
          norm=norm+T(i,j)*T(i,j)
        enddo
        norm=sqrt(norm)
!        print*,'This should be 1 ',norm
        do j=1,N
          T(i,j)=T(i,j)/norm
        enddo
      enddo


100   format (7F8.3,6F7.1)
101   format ((A))
102   format ((A),F12.6)

      return
  END SUBROUTINE Jacobi



! ===================================================================================
SUBROUTINE B_Pa_Diags ( Max,         &   ! Max No. of diagnostic positions
                        Nx,          &   ! No. of x-points
                        Ny,          &   ! No. of y-points
                        n,           &   ! Size of model space
                        B,           &   ! B-matrix
                        Pa,          &   ! Pa-matrix
                        FilenamePos, &   ! Input filename
                        FilenameDiags )  ! Output filename

! Subroutine to show b/g and analysis error standard deviations at obs locations
! Ross Bannister
! National Centre for Earth Observation
! January 2011
! FORTRAN 90

  IMPLICIT NONE

  ! Subroutine arguments
  INTEGER,           INTENT(IN)    :: Max
  INTEGER,           INTENT(IN)    :: Nx
  INTEGER,           INTENT(IN)    :: Ny
  INTEGER,           INTENT(IN)    :: n
  DOUBLE PRECISION,  INTENT(IN)    :: B(1:n,1:n)
  DOUBLE PRECISION,  INTENT(IN)    :: Pa(1:n,1:n)
  CHARACTER (LEN=*), INTENT(IN)    :: FilenamePos
  CHARACTER (LEN=*), INTENT(IN)    :: FilenameDiags

  ! Local variables
  INTEGER                          :: k, element_p, element_u, element_v, npoints
  INTEGER                          :: Diag_x(1:Max), Diag_y(1:Max)
  INTEGER                          :: NDiagPos

  npoints    = Nx*Ny

  ! Read-in the diagnostic positions
  PRINT *, '  Reading-in the diagnostic positions'
  OPEN (12, file=FilenamePos)
  READ (12,*) NDiagPos
  IF (NDiagPos > Max) THEN
    PRINT *, 'Too many diagnostic positions'
    STOP
  END IF
  PRINT *, 'There are ', NDiagPos, ' diagnostic positions'
  DO k = 1, NDiagPos
    READ (12,*) Diag_x(k), Diag_y(k)
  END DO
  CLOSE(12)

  ! Do some checks
  DO k = 1, NDiagPos
    IF ((Diag_x(k) < 1) .OR. (Diag_x(k) > Nx)) THEN
      PRINT *, 'Problem with diagnostic ', k, '(x position out of range)'
      STOP
    END IF
    IF ((Diag_y(k) < 1) .OR. (Diag_y(k) > Ny)) THEN
      PRINT *, 'Problem with diagnostic ', k, '(y position out of range)'
      STOP
    END IF
  END DO

  ! Calculate and output the diagnostics
  PRINT *, '  Calculating and outputting the diagnostics'
  OPEN (12, file=FilenameDiags)
  WRITE (12,*) '=========================================================='
  DO k = 1, NDiagPos
    WRITE (12,'(A,I4,A,2I4)') 'Location of diagnostic ', k, ' = ', Diag_x(k), Diag_y(k)
    element_p = (Diag_x(k)-1)*Ny + Diag_y(k)
    element_u = element_p + npoints
    element_v = element_u + npoints
    WRITE (12,*) 'Background error variances for p, u, v'
    WRITE (12,'(4F16.8)') B(element_p,element_p), &
                          B(element_u,element_u), &
                          B(element_v,element_v)
    WRITE (12,*) 'Analysis error variances for p, u, v'
    WRITE (12,'(4F16.8)') Pa(element_p,element_p), &
                          Pa(element_u,element_u), &
                          Pa(element_v,element_v)
 END DO
  WRITE (12,*) '=========================================================='
  CLOSE (12)

END SUBROUTINE B_Pa_Diags


END PROGRAM Synergy
! ===================================================================================
