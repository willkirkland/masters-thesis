!> @file generalSolver.f90

!> Contains the interface and calls to LAPACK routines for general solutions of
!! linear equations that are not handled by the direct solvers.
!! @author{Will Kirkland; wkirklan@vols.utk.edu}

module generalSolver_mod

    use kind_parameters, only : ikind, rSgle, rDble
    
    implicit none

    !> General interface for solving linear systems using LAPACK routines.
    interface generalSolver
        module procedure lapack_single, lapack_double, lapack_extended, &
                         lapack_quad
    end interface

    integer, parameter :: rExt = 10  ! 10 byte extended precision kind
    integer, parameter :: rQuad = 16 ! 16 byte quad precision kind

contains

    !> Solves linear system AX = B using LAPACK single-precision routines.
    !! @param A Matrix of linear coefficients
    !! @param B Known vector
    !! @param X Unknown vector
    !! @param ndim Number of dimensions of system
    subroutine lapack_single(A, B, X, ndim)

        ! Arguments
        real(kind=rSgle), dimension(ndim,ndim), intent(in) :: A
        real(kind=rSgle), dimension(ndim),intent(in) :: B
        real(kind=rSgle), dimension(ndim),intent(out) :: X
        integer(ikind) :: ndim

        ! Locals
        integer(ikind),dimension(ndim) :: ipiv ! Pivot indices
        integer(ikind) :: info ! LAPACK status

        call sgetrf(ndim, ndim, A, ndim, ipiv, info)
        if (info .ne. 0) then
            if (info .lt. 0) then
                print *,"LAPACK SGETRF error: ",-info,"th argument has illegal value."
            else
                print *,"LAPACK SGETRF error: A(",info,",",info,") has zero value." 
            endif
            STOP
        endif

        X = B  ! Function will overwrite X with solution vector
        call sgetrs('N', ndim, 1, A, ndim, ipiv, X, ndim, info)
        if (info .ne. 0) then
            print *,"LAPACK SGETRS error: ",-info,"the argument has illegal value."
            STOP
        endif

    end subroutine lapack_single

    !> Solves linear system AX = B using LAPACK double-precision routines.
    !! @param A Matrix of linear coefficients
    !! @param B Known vector
    !! @param X Unknown vector
    !! @param ndim Number of dimensions of system
    subroutine lapack_double(A, B, X, ndim)

        ! Arguments
        real(kind=rDble), dimension(ndim,ndim), intent(in) :: A
        real(kind=rDble), dimension(ndim),intent(in) :: B
        real(kind=rDble), dimension(ndim),intent(out) :: X
        integer(ikind) :: ndim

        ! Locals
        real(kind=rDble), dimension(ndim,ndim) :: A_orig
        integer(ikind),dimension(ndim) :: ipiv ! Pivot indices
        integer(ikind) :: info ! LAPACK status
        real(kind=rDble), dimension(1) :: ferr, berr
        real(kind=rDble), dimension(3*ndim) :: work
        integer(ikind), dimension(ndim) :: iwork

        A_orig = A
        call dgetrf(ndim, ndim, A, ndim, ipiv, info)
        if (info .ne. 0) then
            if (info .lt. 0) then
                print *,"LAPACK DGETRF error: ",-info,"th argument has illegal value."
            else
                print *,"LAPACK DGETRF error: A(",info,",",info,") has zero value." 
            endif
            STOP
        endif

        X = B  ! Function will overwrite X with solution vector
        call dgetrs('N', ndim, 1, A, ndim, ipiv, X, ndim, info)
        if (info .ne. 0) then
            print *,"LAPACK DGETRS error: ",-info,"the argument has illegal value."
            STOP
        endif

        !call dgerfs('N', ndim, 1, A_orig, ndim, A, ndim, ipiv, B, ndim, X, &
        !            ndim, ferr, berr, work, iwork, info)
        !if (info .ne. 0) then
        !    print *,"LAPACK DGERFS error: ",-info,"the argument has illegal value."
        !    STOP
        !endif

    end subroutine lapack_double

    !> Solves linear system AX = B using LAPACK double-precision routines when
    !! arguments are extended precision (only accurate to double precision)
    !! @param A Matrix of linear coefficients
    !! @param B Known vector
    !! @param X Unknown vector
    !! @param ndim Number of dimensions of system
    subroutine lapack_extended(A, B, X, ndim)

        ! Arguments
        real(kind=rExt), dimension(ndim,ndim), intent(in) :: A
        real(kind=rExt), dimension(ndim),intent(in) :: B
        real(kind=rExt), dimension(ndim),intent(out) :: X
        integer(ikind) :: ndim

        ! Locals
        integer(ikind),dimension(ndim) :: ipiv ! Pivot indices
        integer(ikind) :: info ! LAPACK status

        call dgetrf(ndim, ndim, real(A,rDble), ndim, ipiv, info)
        if (info .ne. 0) then
            if (info .lt. 0) then
                print *,"LAPACK DGETRF error: ",-info,"th argument has illegal value."
            else
                print *,"LAPACK DGETRF error: A(",info,",",info,") has zero value." 
            endif
            STOP
        endif

        X = B  ! Function will overwrite X with solution vector
        call dgetrs('N', ndim, 1, real(A,rDble), ndim, ipiv, real(X,rDble), ndim, info)
        if (info .ne. 0) then
            print *,"LAPACK DGETRS error: ",-info,"the argument has illegal value."
            STOP
        endif

    end subroutine lapack_extended

    !> Solves linear system AX = B using LAPACK double-precision routines when
    !! arguments are quad precision (only accurate to double precision)
    !! @param A Matrix of linear coefficients
    !! @param B Known vector
    !! @param X Unknown vector
    !! @param ndim Number of dimensions of system
    subroutine lapack_quad(A, B, X, ndim)

        ! Arguments
        real(kind=rQuad), dimension(ndim,ndim), intent(in) :: A
        real(kind=rQuad), dimension(ndim),intent(in) :: B
        real(kind=rQuad), dimension(ndim),intent(out) :: X
        integer(ikind) :: ndim

        ! Locals
        integer(ikind),dimension(ndim) :: ipiv ! Pivot indices
        integer(ikind) :: info ! LAPACK status

        call dgetrf(ndim, ndim, real(A,rDble), ndim, ipiv, info)
        if (info .ne. 0) then
            if (info .lt. 0) then
                print *,"LAPACK DGETRF error: ",-info,"th argument has illegal value."
            else
                print *,"LAPACK DGETRF error: A(",info,",",info,") has zero value." 
            endif
            STOP
        endif

        X = B  ! Function will overwrite X with solution vector
        call dgetrs('N', ndim, 1, real(A,rDble), ndim, ipiv, real(X,rDble), ndim, info)
        if (info .ne. 0) then
            print *,"LAPACK DGETRS error: ",-info,"the argument has illegal value."
            STOP
        endif

    end subroutine lapack_quad

end module generalSolver_mod

