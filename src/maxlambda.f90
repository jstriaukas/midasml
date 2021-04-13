! ------------------------------------------------------------------
! maxlambda.f90: compute maxlam such that for a given gamma
!
!              min_b |y-Xb|_T^2 + maxlam Omega_gamma(b)
!
!             gives b = 0
!
! ------------------------------------------------------------------
SUBROUTINE maxlambda(nvars, nobs, x, y, gamma, gindex, ngroups, pf, maxlam)

    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: nvars, nobs, ngroups, gindex(ngroups)
    DOUBLE PRECISION :: x(nobs,nvars), y(nobs), pf(nvars), gamma, maxlam
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: k, c, nzvars
    INTEGER :: gstart, gend, gs, gj
    DOUBLE PRECISION :: gw, xy(nvars), r(nobs)
    DOUBLE PRECISION :: wmaxg(ngroups), lb, rb

    ! -------- BEGIN PROGRAM -------- !
    c = 0
    r = y
    nzvars = 0
    DO k = 1, nvars
        IF (pf(k) .EQ. 0.0D0) THEN
        nzvars = nzvars + 1
        END IF
    END DO
    IF (nzvars .NE. 0) THEN
        !CALL rnz(nvars, nobs, nzvars, y, x, r, pf)
    END IF
    xy = MATMUL(TRANSPOSE(x),r)/nobs


    IF (gamma .EQ. 1.0D0) THEN
        maxlam = MAXVAL(ABS(xy))
    ELSE
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
            IF (gw == 0.0D0) THEN
                wmaxg(k) = 0.0D0
            ELSE
                IF (gamma .EQ. 0.0D0) THEN
                    rb = NORM2(xy(gstart:gend))
                    wmaxg(k) = rb/gw
                ELSE
                    lb = 0.0D0
                    rb = MAXVAL(ABS(xy(gstart:gend)))/gamma
                    CALL solvewmaxg(gstart, gend, gamma, lb, rb, gw, pf, xy, nvars)
                    wmaxg(k) = rb
                END IF
            END IF
        END DO
        maxlam = MAXVAL(wmaxg)
    END IF
    !--- ADD SMALL NUMBER TO ENSURE b = 0 @ maxlam (DUE TO NUMERICAL IMPRESSION)
    maxlam = maxlam + 1E-5

END SUBROUTINE maxlambda


SUBROUTINE  solvewmaxg(gstart, gend, gamma, lb, rb, gw, pf, xy, nvars)
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: gstart, gend, nvars
    DOUBLE PRECISION :: gamma, lb, rb, gw, pf(nvars), xy(nvars)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: stopflag, indexi
    DOUBLE PRECISION ::  tol = 1E-13, mp, fl, fm, fr, tmpl, tmpm, tmpr

    stopflag = 0
    DO
        mp = 0.5 * (lb + rb)
        fl = 0.0D0
        fm = 0.0D0
        fr = 0.0D0
        tmpl = 0.0D0
        tmpm = 0.0D0
        tmpr = 0.0D0
        DO indexi =  gstart, gend
            tmpl = ABS(xy(indexi)) - gamma * lb * pf(indexi)
            tmpm = ABS(xy(indexi)) - gamma * mp * pf(indexi)
            tmpr = ABS(xy(indexi)) - gamma * rb * pf(indexi)
            IF (tmpl > 0.0D0) THEN
                fl = fl + tmpl * tmpl
            END IF
            IF (tmpm > 0.0D0) THEN
                fm = fm + tmpm * tmpm
            END IF
            IF (tmpr > 0.0D0) THEN
                fr = fr + tmpr * tmpr
            END IF
        END DO
        fl = fl - (1.0D0 - gamma) * (1.0D0 - gamma) * lb * lb * gw * gw
        fm = fm - (1.0D0 - gamma) * (1.0D0 - gamma) * mp * mp * gw * gw
        fr = fr - (1.0D0 - gamma) * (1.0D0 - gamma) * rb * rb * gw * gw
        IF (fl * fm < 0.0D0) THEN
            IF (ABS(lb - mp) > tol) THEN
                rb = mp
            ELSE
                stopflag = 1
            END IF
        ELSE
            IF (fm * fr < 0.0D0) THEN
                IF (ABS(mp - rb) > tol) THEN
                    lb = mp
                ELSE
                    stopflag = 1
                END IF
            ELSE
                stopflag = 1
            END IF
        END IF
       IF (stopflag .EQ. 1) EXIT
    END DO
    rb = mp

END SUBROUTINE solvewmaxg


SUBROUTINE rnz(nvars, nobs, nzvars, y, x, r, pf)

    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: nvars, nobs, nzvars
    DOUBLE PRECISION :: x(nobs,nvars), y(nobs), pf(nvars), r(nobs), xact(nobs,nzvars), yn(nobs), xactn(nobs,nzvars)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: i, j, c, lwork
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: work
    ALLOCATE(work(100 * nobs * nzvars))
    ! -------- BEGIN PROGRAM -------- !
    c = 0
    r = y
    yn = y
    lwork = 100 * nobs * nzvars
    DO i = 1, nvars
        IF (pf(i) .EQ. 0.0D0) THEN
            c = c + 1
            xact(:,c) = x(:,i)
        END IF
    END DO
    xactn = xact
    !CALL DGELS('No transpose', nobs, nzvars, 1, xactn, nobs, yn, nobs, work, lwork, info)
    DO j = 1, nzvars
        r = r - xact(:,j)*yn(j)
    END DO

END SUBROUTINE rnz
