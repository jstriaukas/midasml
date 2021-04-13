! ------------------------------------------------------------------------
! panelsglfitF.f90: block coordinate descent for sg-LASSO panel regression.
! ------------------------------------------------------------------------
SUBROUTINE panelsglfitF(nf, T, gamma, ngroups, gindex, nobs, nvars, x, y, pf, dfmax, pmax, &
& nlam, flmin, ulam, eps, peps, isd, intr, maxit, nalam, b0, a0, beta, ibeta, nbeta, &
& alam, npass, jerr)
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: nobs, nvars, dfmax, pmax, nlam, nalam, isd, intr, ngroups, nf, T
    INTEGER :: npass, jerr, maxit, gindex(ngroups)
    INTEGER :: ibeta(pmax), nbeta(nlam)
    DOUBLE PRECISION :: flmin, eps, peps, gamma, maxlam
    DOUBLE PRECISION :: x(nobs, nvars), y(nobs)
    DOUBLE PRECISION :: pf(nvars)
    DOUBLE PRECISION :: beta(pmax, nlam), b0(nlam), a0(nf, nlam)
    DOUBLE PRECISION :: ulam(nlam), alam(nlam), tmp
    DOUBLE PRECISION :: ymb(nf), xmb(nf, nvars), yn(nobs), xn(nobs,nvars)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: j, l, nk, ierr, k, festart, feend
    INTEGER, DIMENSION(:), ALLOCATABLE :: ju
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xmean, xnorm, maj
    ! -------- ALLOCATE VARIABLES -------- !
    ALLOCATE(ju(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE(xmean(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE(xnorm(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE(maj(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    IF (jerr /= 0) RETURN
    CALL chkvars(nobs, nvars, x, ju)
    IF (MAXVAL(pf) <= 0.0D0) THEN
    jerr = 10000
    RETURN
    END IF
    pf = MAX(0.0D0, pf)

    ! -------------------- FIXED EFFECTS --------------------- !
    DO k = 1, nf
        festart = (k-1)*T+1
        feend = k*T
        ymb(k) = SUM(y(festart:feend))/T
        yn(festart:feend) = y(festart:feend) - ymb(k)
        DO j = 1, nvars
            xmb(k,j) = SUM(x(festart:feend,j))/T
            xn(festart:feend,j) = x(festart:feend,j) - xmb(k,j)
        END DO

    END DO
    ! -------------------- COMPUTE MAJ --------------------- !
    isd = 0
    intr = 0
    CALL standard(nobs, nvars, xn, ju, isd, intr, xmean, xnorm, maj)
    ! -------------------- COMPUTE LAMBDA --------------------- !
    IF (ulam(1) .EQ. -1.0D0) THEN
        CALL maxlambda(nvars, nobs, xn, yn, gamma, gindex, ngroups, pf, maxlam)
        ulam(1) = maxlam
        DO j = 2, nlam
            tmp = LOG(maxlam) + (LOG(maxlam*flmin) - LOG(maxlam)) * (j - 1) / (nlam - 1)
            ulam(j) = EXP(tmp)
        END DO
    END IF

    IF (gamma == 1.0D0) THEN
        ! -------------------- CALL lassofitpathF --------------------- !
        CALL lassofitpathF(maj, nobs, nvars, xn, yn, ju, pf, dfmax, &
                           & pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, ibeta, nbeta, alam, &
                           & npass, jerr, intr)
    ELSE
        ! -------------------- CALL sglfitpathF --------------------- !
        CALL sglfitpathF(maj, gamma, ngroups, gindex, nobs, nvars, xn, yn, ju, pf, dfmax, &
                         & pmax, nlam, flmin, ulam, eps, peps, maxit, nalam, b0, beta, ibeta, &
                         & nbeta, alam, npass, jerr, intr)
    END IF
    IF (jerr > 0) RETURN ! CHECK ERROR AFTER CALLING FUNCTION
    ! ----------- TRANSFORM BETA BACK TO THE ORIGINAL SCALE ----------- !
    DO l = 1, nalam
        nk = nbeta(l)
        IF (isd == 1) THEN
            DO j = 1, nk
                beta(j,l) = beta(j,l)/xnorm(ibeta(j))
            END DO
        END IF
    END DO
    ! ----------- COMPUTE FIXED EFFECTS ----------- !
    DO l = 1, nlam
        nk = nbeta(l)
        a0(:,l) = ymb -  MATMUL(xmb(:,ibeta(1:nk)),beta(1:nk,l))
    END DO
    RETURN
END SUBROUTINE panelsglfitF
