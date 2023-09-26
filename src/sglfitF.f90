! ------------------------------------------------------------------
! sglfitF.f90: block coordinate descent for sg-LASSO MSE regression.
! ------------------------------------------------------------------
SUBROUTINE sglfitF(gamma, ngroups, gindex, nobs, nvars, x, y, pf, dfmax, pmax, &
& nlam, flmin, ulam, eps, peps, isd, intr, maxit, nalam, b0, beta, ibeta, nbeta, &
& alam, npass, jerr)

  IMPLICIT NONE
  ! -------- INPUT VARIABLES -------- !
  INTEGER :: nobs, nvars, dfmax, pmax, nlam, nalam, isd, intr, ngroups
  INTEGER :: npass, jerr, maxit, gindex(ngroups)
  INTEGER :: ibeta(pmax), nbeta(nlam)
  DOUBLE PRECISION :: flmin, eps, peps, gamma, maxlam
  DOUBLE PRECISION :: x(nobs, nvars), y(nobs)
  DOUBLE PRECISION :: pf(nvars)
  DOUBLE PRECISION :: beta(pmax, nlam), b0(nlam)
  DOUBLE PRECISION :: ulam(nlam), alam(nlam), tmp
  ! -------- LOCAL DECLARATIONS -------- !
  INTEGER :: j, l, nk, ierr
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
  ! -------------------- STANDARDIZE & COMPUTE MAJ --------------------- !
  CALL standard(nobs, nvars, x, ju, isd, intr, xmean, xnorm, maj)
  ! -------------------- COMPUTE LAMBDA --------------------- !
  IF (ulam(1) .EQ. -1.0D0) THEN
    CALL maxlambda(nvars, nobs, x, y, gamma, gindex, ngroups, pf, maxlam)
    ulam(1) = maxlam
    DO j = 2, nlam
        tmp = LOG(maxlam) + (LOG(maxlam*flmin) - LOG(maxlam)) * (j - 1) / (nlam - 1)
        ulam(j) = EXP(tmp)
    END DO
  END IF

  IF (gamma == 1.0D0) THEN
    ! -------------------- CALL lassofitpathF --------------------- !
    CALL lassofitpathF(maj, nobs, nvars, x, y, ju, pf, dfmax, &
      & pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, ibeta, nbeta, alam, &
      & npass, jerr, intr)
  ELSE
    ! -------------------- CALL sglfitpathF --------------------- !
    CALL sglfitpathF(maj, gamma, ngroups, gindex, nobs, nvars, x, y, ju, pf, dfmax, &
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
    IF (intr == 1) THEN
        b0(l) = b0(l) - DOT_PRODUCT(beta(1:nk,l),xmean(ibeta(1:nk)))
    END IF
  END DO
  DEALLOCATE(ju,xmean,xnorm,maj)
  RETURN
END SUBROUTINE sglfitF


! --------------------------------- sglfitpathF --------------------------------- !
SUBROUTINE sglfitpathF(maj, gamma, ngroups, gindex, nobs, nvars, x, y, ju, pf, dfmax, &
& pmax, nlam, flmin, ulam, eps, peps, maxit, nalam, b0, beta, m, nbeta, alam, &
& npass, jerr, intr)

  IMPLICIT NONE
  ! -------- INPUT VARIABLES -------- !
  INTEGER :: mnl, nobs, nvars, dfmax, pmax, nlam, maxit, nalam, npass, jerr, intr, ngroups
  INTEGER :: ju(nvars), m(pmax), nbeta(nlam), gindex(ngroups)
  DOUBLE PRECISION :: eps, gamma, peps, steps(ngroups)
  DOUBLE PRECISION :: x(nobs, nvars), y(nobs), maj(nvars)
  DOUBLE PRECISION :: pf(nvars)
  DOUBLE PRECISION :: beta(pmax, nlam), b0(nlam)
  DOUBLE PRECISION :: ulam(nlam), alam(nlam)
  ! -------- LOCAL DECLARATIONS -------- !
  INTEGER,  PARAMETER :: mnlam = 6
  DOUBLE PRECISION :: d, dif, oldb, u, al, flmin
  DOUBLE PRECISION,  DIMENSION(:), ALLOCATABLE :: b, oldbeta, r
  DOUBLE PRECISION :: gw, tmp
  INTEGER :: gstart, gend
  INTEGER :: k, j, l, g, vrg, ctr, ierr, ni, me, pln
  INTEGER :: gs, gj, skip
  INTEGER,  DIMENSION(:),  ALLOCATABLE :: mm
  ! -------- ALLOCATE VARIABLES -------- !
  ALLOCATE(b(0:nvars), STAT=jerr)
  ALLOCATE(oldbeta(0:nvars), STAT=ierr)
  jerr = jerr + ierr
  ALLOCATE(mm(1:nvars), STAT=ierr)
  jerr = jerr + ierr
  ALLOCATE(r(1:nobs), STAT=ierr)
  jerr = jerr + ierr
  IF (jerr /= 0) RETURN
  ! ---------------- INITIALIZATION ---------------- !
  b = 0.0D0
  oldbeta = 0.0D0
  m = 0
  mm = 0
  npass = 0
  ni = npass
  mnl = MIN(mnlam,nlam)
  ju = 0
  r = y
  DO k = 1, ngroups
    gend = gindex(k)
    IF (k == 1) THEN
        gstart = 1
    ELSE
        gstart = gindex(k-1) + 1
    END IF
    tmp = SUM(MATMUL(TRANSPOSE(x(:,gstart:gend)),x(:,gstart:gend))/nobs * MATMUL(TRANSPOSE(x(:,gstart:gend)),x(:,gstart:gend))/nobs)
    steps(k) = 1/SQRT(tmp)
  END DO

  ! ----------------- LAMBDA LOOP (OUTMOST LOOP) ------------------- !
  DO l = 1, nlam
        al = ulam(l)
        ctr = 0
        pln = 0
        ! ------------------ OUTER LOOP -------------------- !
        DO
            IF (intr == 1) oldbeta(0) = b(0)
            IF (ni > 0) oldbeta(m(1:ni)) = b(m(1:ni))
            ! ----------------- MIDDLE LOOP -------------------- !
            DO
                npass = npass + 1
                dif = 0.0D0
                IF (intr == 1) oldbeta(0) = b(0)
                IF (ni > 0) oldbeta(m(1:ni)) = b(m(1:ni))
                pln = pln + 1
                ! ----------------- GROUP LOOP -------------------- !
                DO k = 1, ngroups
                    gend = gindex(k)
                    IF (k == 1) THEN
                        gstart = 1
                    ELSE
                        gstart = gindex(k-1) + 1
                    END IF
                    gs = gend - gstart + 1
                    gw = 0.0D0
                    DO gj = gstart, gend
                        gw = gw + pf(gj)
                    END DO
                    gw = SQRT(gw)
                    skip = 1
                    IF (pln == 1) THEN
                        skip = 0
                    END IF
                    IF (ju(k) == 1) THEN
                        skip = 0
                    END IF
                    IF (skip == 0) THEN
                        ! --- sg-LASSO PROXIMAL MAP --- !
                        CALL prox_sgl(gstart, gend, nvars, nobs, x, r, b(1:nvars), al, gamma, pf, peps, gw, steps(k))
                        ! UPDATE REMAINING VARIABLES
                        DO g = gstart, gend
                            !IF (ABS(b(g))<eps) b(g) = 0.0D0
                            IF (ABS(b(g))>0.0D0) THEN
                                IF (pln == 1) THEN
                                    ju(k) = 1
                                END IF
                                d = oldbeta(g) - b(g)
                                dif = MAX(dif, maj(g) * d**2)
                                IF (mm(g) == 0) THEN
                                    ni = ni + 1
                                    IF (ni > pmax) EXIT
                                    mm(g) = ni
                                    m(ni) = g ! RECORD ACTIVE VARIABLES
                                END IF
                            END IF
                        END DO
                    END IF ! ----------> END GROUP UPDATES
                END DO ! ----------> END GROUP LOOP
                IF (ni > pmax) EXIT
                IF (intr == 1) THEN
                    oldb = oldbeta(0)
                    u = 0.0D0
                    DO ! BEGIN GRADIENT DESCENT
                      d = SUM(r)/nobs
                      IF (d**2 < eps) EXIT
                      b(0) = b(0) + d
                      r = r - d
                    END DO ! END GRADIENT DESCENT
                    d = b(0) - oldb
                    IF (ABS(d) > 0.0D0) dif = MAX(dif, d**2)
                END IF
                IF (dif < eps) EXIT
            END DO ! ----------> END MIDDLE LOOP
            IF (ni > pmax) EXIT
            ! -------------- FINAL CHECK ---------------- !
            vrg = 1
            IF (intr == 1) THEN
            IF ((b(0) - oldbeta(0))**2 >= eps) vrg = 0
            END IF
            DO j = 1, ni
                IF ((b(m(j)) - oldbeta(m(j)))**2 >= eps) THEN
                    vrg = 0
                    EXIT
                END IF
            END DO
            IF (vrg == 1) EXIT
            ctr = ctr + 1
            IF (ctr > maxit) THEN
             jerr = - l
             RETURN
            END IF
        END DO ! -------> END OUTER LOOP
        ! ----------- FINAL UPDATE & SAVE RESULTS ------------ !
        IF (ni > pmax) THEN
         jerr = - 10000 - l
         EXIT
        END IF
        IF (ni > 0) beta(1:ni,l) = b(m(1:ni))
        nbeta(l) = ni
        b0(l) = b(0)
        alam(l) = al
        nalam = l
        IF (l < mnl) CYCLE
        IF (flmin >= 1.0D0) CYCLE
        me = COUNT(ABS(beta(1:ni,l)) > 0.0D0)
        IF (me > dfmax) EXIT
  END DO ! -------> END LAMBDA LOOP
  DEALLOCATE(b,oldbeta,mm,r)
  RETURN
END SUBROUTINE sglfitpathF


! --------------------------------- lassofitpathF --------------------------------- !
SUBROUTINE lassofitpathF(maj, nobs, nvars, x, y, ju, pf, dfmax, &
& pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, m, nbeta, alam, &
& npass, jerr, intr)

    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: mnl, nobs, nvars, dfmax, pmax, nlam, maxit, nalam, npass, jerr, intr
    INTEGER :: ju(nvars), m(pmax), nbeta(nlam)
    DOUBLE PRECISION :: eps
    DOUBLE PRECISION :: x(nobs, nvars), y(nobs), maj(nvars)
    DOUBLE PRECISION :: pf(nvars)
    DOUBLE PRECISION :: beta(pmax, nlam), b0(nlam)
    DOUBLE PRECISION :: ulam(nlam), alam(nlam)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER,  PARAMETER :: mnlam = 6
    DOUBLE PRECISION :: tmp, d, dif, oldb, u, v, al, flmin
    DOUBLE PRECISION,  DIMENSION(:),  ALLOCATABLE :: b, oldbeta, r

    INTEGER :: k, j, l, vrg, ctr, ierr, ni, me, pln, iter
    INTEGER,  DIMENSION(:),  ALLOCATABLE :: mm
    ! -------- ALLOCATE VARIABLES -------- !
    ALLOCATE(b(0:nvars), STAT=jerr)
    ALLOCATE(oldbeta(0:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE(mm(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE(r(1:nobs), STAT=ierr)
    jerr = jerr + ierr
    IF (jerr /= 0) RETURN
    ! ---------------- INITIALIZATION ---------------- !
    r = y
    b = 0.0D0
    oldbeta = 0.0D0
    m = 0
    mm = 0
    npass = 0
    ni = npass
    mnl = MIN(mnlam,nlam)
    ju = 0
    ! ----------------- LAMBDA LOOP (OUTMOST LOOP) ------------------- !
    DO l = 1, nlam
        al = ulam(l)
        ctr = 0
        ! ------------------ OUTER LOOP -------------------- !
        DO
            IF (intr == 1) oldbeta(0) = b(0)
            IF (ni > 0) oldbeta(m(1:ni)) = b(m(1:ni))
            pln = 0
            ! ----------------- MIDDLE LOOP -------------------- !
            DO
                npass = npass + 1
                dif = 0.0D0
                pln = pln + 1
                DO k = 1, nvars
                    IF (pln == 1) THEN
                        oldb = b(k)
                        iter = 1
                        DO ! BEGIN PROXIMAL COORDINATE DESCENT
                            !u = 0.0D0
                            u = maj(k) * b(k) + DOT_PRODUCT(r,x(:,k))/nobs
                            v = ABS(u) - al * pf(k)
                            IF (v > 0.0D0) THEN
                                tmp = SIGN(v,u)/maj(k)
                            ELSE
                                tmp = 0.0D0
                            END IF
                            d = tmp - b(k)
                            IF (d**2 < eps) EXIT
                            IF (iter > maxit) EXIT
                            b(k) = tmp
                            r = r - x(:,k) * d
                            iter = iter + 1
                        END DO ! END PROXIMAL GRADIENT DESCENT
                        d = b(k) - oldb
                        IF (ABS(d) > 0.0D0) THEN
                            dif = dif + maj(k) * d**2
                            IF (mm(k) == 0) THEN
                                ni = ni + 1
                                IF (ni > pmax) EXIT
                                mm(k) = ni
                                m(ni) = k ! RECORD ACTIVE VARIABLES
                            END IF
                        END IF
                        IF (ABS(b(k))>0.0D0) ju(k) = 1
                    ELSE
                            IF (ju(k) == 1) THEN
                                oldb = b(k)
                                DO ! BEGIN PROXIMAL GRADIENT DESCENT
                                    !u = 0.0D0
                                    u = DOT_PRODUCT(r,x(:,k))
                                    u = maj(k) * b(k) + u/nobs
                                    v = ABS(u) - al * pf(k)
                                    IF (v > 0.0D0) THEN
                                        tmp = SIGN(v,u)/maj(k)
                                    ELSE
                                        tmp = 0.0D0
                                    END IF
                                    d = tmp - b(k)
                                    IF (d**2 < eps) EXIT
                                    b(k) = tmp
                                    r = r - x(:,k) * d
                                END DO ! END PROXIMAL GRADIENT DESCENT
                                d = b(k) - oldb
                                IF (ABS(d) > 0.0D0) THEN
                                        dif = MAX(dif, maj(k) * d**2)
                                END IF
                            END IF
                    END IF
                END DO
                IF (ni > pmax) EXIT
                IF (intr == 1) THEN
                    oldb = b(0)
                    iter = 1
                    DO ! BEGIN GRADIENT DESCENT
                        d = SUM(r)/nobs
                        IF (d**2 < eps) EXIT
                        IF (iter > maxit) EXIT
                        b(0) = b(0) + d
                        r = r - d
                        iter = iter + 1
                    END DO ! END GRADIENT DESCENT
                    d = b(0) - oldb
                    IF (ABS(d) > 0.0D0) dif = MAX(dif, d**2)
                END IF
                IF (dif < eps) EXIT
            END DO ! ----------> END MIDDLE LOOP
            IF (ni > pmax) EXIT
            ! -------------- FINAL CHECK ---------------- !
            vrg = 1
            IF (intr == 1) THEN
                IF ((b(0) - oldbeta(0))**2 >= eps) vrg = 0
            END IF
            DO j = 1, ni
                IF ((b(m(j)) - oldbeta(m(j)))**2 >= eps) THEN
                    vrg = 0
                    EXIT
                END IF
            END DO
            IF (vrg == 1) EXIT
            ctr = ctr + 1
            IF (ctr > maxit) THEN
                jerr = - l
                RETURN
            END IF
        END DO ! -------> END OUTER LOOP
        ! ----------- FINAL UPDATE & SAVE RESULTS ------------ !
        IF (ni > pmax) THEN
            jerr = - 10000 - l
            EXIT
        END IF
        IF (ni > 0) beta(1:ni,l) = b(m(1:ni))
        nbeta(l) = ni
        b0(l) = b(0)
        alam(l) = al
        nalam = l
        IF (l < mnl) CYCLE
        IF (flmin >= 1.0D0) CYCLE
        me = COUNT(ABS(beta(1:ni,l)) > 0.0D0)
        IF (me > dfmax) EXIT
    END DO ! -------> END LAMBDA LOOP
    DEALLOCATE(b,oldbeta,r,mm)
    RETURN
END SUBROUTINE lassofitpathF



SUBROUTINE prox_sgl(gstart, gend, nvars, nobs, x, r, b, al, gamma, pf, peps, gw, step)

    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: gstart, gend, nvars, nobs
    DOUBLE PRECISION :: x(nobs,nvars), r(nobs), b(nvars), al, gamma, pf(nvars), peps, gw, step
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: g
    DOUBLE PRECISION :: u, v, scl, tmp, maxg, normg, d, bold(nvars), vg, big = 9.9D30, s
    s = step
    ! -------- BEGIN PROGRAM -------- !
    DO
        bold(gstart:gend) = b(gstart:gend)
        !--------- LASSO PART ----------!
        DO g = gstart, gend
            u = b(g) + s * DOT_PRODUCT(x(:,g),r)/nobs
            !S(.) map
            !tmp = max(u - al * gamma * pf(g), 0.0D0) - max(-u - al * gamma * pf(g), 0.0D0)
            v = ABS(u) - s * al * gamma * pf(g)
            IF (v > 0.0D0) THEN
                tmp = SIGN(v,u)
            ELSE
                tmp = 0.0D0
            END IF
            b(g) = tmp
        END DO
        !--------- g-LASSO PART ----------!
        ! L2 norm of b_g
        !normg = NORM2(b(gstart:gend))
        normg = SQRT(DOT_PRODUCT(b(gstart:gend), b(gstart:gend)))
        ! Initialize storage vars
        maxg = 0.0D0
        vg = s * gw * al * (1.0D0-gamma)/normg
        IF (normg .EQ. 0.0D0) THEN
            vg = big
        END IF

        DO g = gstart, gend
            scl = 1.0D0 - pf(g) * vg
            scl = MAX(scl, 0.0D0)
            !l_2,1 norm map
            tmp = scl*b(g)
            d = tmp - bold(g)
            r = r - x(:,g)*d
            maxg = MAX(maxg, ABS(d))
            b(g) = tmp
        END DO
        !--------- CHECK CONVERGENCE ----------!
        IF (maxg < peps) EXIT
    END DO
END SUBROUTINE prox_sgl
