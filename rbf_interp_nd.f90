module rbf_interp_nd_mod
    use kind_parameters, only: ikind, rkind
  implicit none

contains

subroutine rbf_interp_nd ( m, nd, xd, r0, phi, w, xi, fi )

!*****************************************************************************80
!
!! RBF_INTERP_ND evaluates a radial basis function interpolant.
!
!  The original version of RBF_INTERP_ND by John Burkardt is available at 
!  http://people.sc.fsu.edu/~jburkardt/f_src/rbf_interp_nd/rbf_interp_nd.f90
!
!  This version of RBF_INTERP_ND has been modified to include constant and 
!  linear terms in the radial basis function fitting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    5 Aug 2015
!
!  Author:
!
!    John Burkardt, modified by Will Kirkland, Univ. of Tennessee
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer ( kind = ikind ) M, the spatial dimension.
!
!    Input, integer ( kind = ikind ) ND, the number of data points.
!
!    Input, real ( kind = rkind ) XD(M,ND), the data points.
!
!    Input, real ( kind = rkind ) R0, a scale factor.  R0 should be larger than 
!    the typical separation between points, but smaller than the maximum 
!    separation.  The value of R0 has a significant effect on the resulting 
!    interpolant.
!
!    Input, subroutine PHI ( N, R, R0, V ), a subroutine to evaluate the radial
!    basis functions.
!
!    Input, real ( kind = rkind ) W(ND), the weights, as computed by RBF_WEIGHTS.
!
!    Input, real ( kind = rkind ) XI(M), the interpolation point.
!
!    Output, real ( kind = rkind ) FI, the interpolated value.
!
  implicit none

  integer ( kind = ikind ) m
  integer ( kind = ikind ) nd

  real ( kind = rkind ) fi
  integer ( kind = ikind ) j
  external phi
  real ( kind = rkind ) r(nd)
  real ( kind = rkind ) r0
  real ( kind = rkind ) v(nd)
  real ( kind = rkind ) w(nd+m+1)
  real ( kind = rkind ) xd(m,nd)
  real ( kind = rkind ) xi(m)
  real ( kind = rkind ) xi_aug(m+1)

    do j = 1, nd
    r(j) = sqrt ( sum ( ( xi(1:m) - xd(1:m,j) )**2 ) )
    end do

    call phi ( nd, r, r0, v )

  fi = dot_product ( v, w(1:nd) )

    xi_aug(1) = 1
    do j = 1,m
    xi_aug(j+1) = xi(j)
    enddo
  fi = fi + dot_product(w(nd+1:nd+m+1),xi_aug)

  return
end subroutine rbf_interp_nd

subroutine rbf_weight ( m, nd, nrhs, xd, r0, phi, fd, w)

!*****************************************************************************80
!
!! RBF_WEIGHT computes weights for radial basis function interpolation.
!
!  The original version of RBF_INTERP_ND by John Burkardt is available at 
!  http://people.sc.fsu.edu/~jburkardt/f_src/rbf_interp_nd/rbf_interp_nd.f90
!
!  This version of RBF_INTERP_ND has been modified to include constant and 
!  linear terms in the radial basis function fitting and to replace LINPACK 
!  routines with calls to LAPACK (assumed to be included as library at compile
!  time).
!
!  Discussion:
!
!    We assume that there are N (nonsingular) equations in N unknowns.
!
!    However, it should be clear that, if we are willing to do some kind
!    of least squares calculation, we could allow for singularity,
!    inconsistency, or underdetermine systems.  This could be associated
!    with data points that are very close or repeated, a smaller number
!    of data points than function values, or some other ill-conditioning
!    of the system arising from a peculiarity in the point spacing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    5 Aug 2015
!
!  Author:
!
!    John Burkardt, modified by Will Kirkland, Univ. of Tennessee
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer ( kind = ikind ) M, the spatial dimension.
!
!    Input, integer ( kind = ikind ) ND, the number of data points.
!
!    Input, integer ( kind = ikind ) NRHS, the number of right-hand sides 
!
!    Input, real ( kind = rkind ) XD(M,ND), the data points.
!
!    Input, real ( kind = rkind ) R0, a scale factor.  R0 should be larger than 
!    the typical separation between points, but smaller than the maximum 
!    separation.  The value of R0 has a significant effect on the resulting 
!    interpolant.
!
!    Input, subroutine PHI ( N, R, R0, V ), a subroutine to evaluate the radial
!    basis functions.
!
!    Input, real ( kind = rkind ) FD(ND), the function values at the data points.
!
!    Output, real ( kind = rkind ) W(ND), the weights.
!
  implicit none

  integer ( kind = ikind ) m
  integer ( kind = ikind ) nd
  integer ( kind = ikind ) nrhs

  real ( kind = rkind ) a(nd+m+1, nd+m+1) 
  real ( kind = rkind ) fd(nd, nrhs)
  integer ( kind = ikind ) i
  integer ( kind = ikind ) j
  external phi
  real ( kind = rkind ) r(nd)
  real ( kind = rkind ) r0
  real ( kind = rkind ) v(nd)
  real ( kind = rkind ) w(nd+m+1, nrhs)
  real ( kind = rkind ) xd(m,nd)
  real ( kind = rkind ) fd_aug(nd+m+1, nrhs)
  integer ( kind = ikind ) ipiv(nd+m+1)
  integer ( kind = ikind ) info

  a = 0. 
  fd_aug = 0.
  ipiv = 0
  info = 0

  do i = 1, nd

    do j = 1, nd
      r(j) = sqrt ( sum ( ( xd(1:m,i) - xd(1:m,j) )**2 ) )
    end do

    call phi ( nd, r, r0, v )

    a(i,1:nd) = v(1:nd)

    a(i,nd+1) = 1.
    do j = 1, m
      a(i,nd+1+j) = xd(j,i)
    enddo
    a(nd+1,i) = 1.
   
  end do
  
  do i=1,m
    do j=1,nd
      a(nd+1+i,j) = xd(i,j)
    enddo
  enddo
  fd_aug(1:nd, :) = fd(1:nd, :)

!
!  Solve for the weights.
!
  w = fd_aug
  call dgetrf ( nd+m+1, nd+m+1, a, nd+m+1, ipiv, info )
  call dgetrs ( 'N', nd+m+1, nrhs, a, nd+m+1, ipiv, w, nd+m+1, info)

  return
end subroutine rbf_weight

subroutine phi1 ( n, r, r0, v )

!*****************************************************************************80
!
!! PHI1 evaluates the multiquadric radial basis function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer ( kind = ikind ) N, the number of points.
!
!    Input, real ( kind = rkind ) R(N), the radial separation.
!    0 < R.
!
!    Input, real ( kind = rkind ) R0, a scale factor.
!
!    Output, real ( kind = rkind ) V(N), the value of the radial basis function.
!
  implicit none

  integer ( kind = ikind ) n

  real ( kind = rkind ) r(n)
  real ( kind = rkind ) r0
  real ( kind = rkind ) v(n)

  v(1:n) = sqrt ( r(1:n)**2 + r0**2 )

  return
end subroutine phi1

subroutine phi2 ( n, r, r0, v )

!*****************************************************************************80
!
!! PHI2 evaluates the inverse multiquadric radial basis function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer ( kind = ikind ) N, the number of points.
!
!    Input, real ( kind = rkind ) R(N), the radial separation.
!    0 < R.
!
!    Input, real ( kind = rkind ) R0, a scale factor.
!
!    Output, real ( kind = rkind ) V(N), the value of the radial basis function.
!
  implicit none

  integer ( kind = ikind ) n

  real ( kind = rkind ) r(n)
  real ( kind = rkind ) r0
  real ( kind = rkind ) v(n)

  v = 1.0D+00 / sqrt ( r**2 + r0**2 )

  return
end subroutine phi2

subroutine phi3 ( n, r, r0, v )

!*****************************************************************************80
!
!! PHI3 evaluates the thin-plate spline radial basis function.
!
!  Discussion:
!
!    Note that PHI3(R,R0) is negative if R < R0.  Thus, for this basis function,
!    it may be desirable to choose a value of R0 smaller than any possible R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer ( kind = ikind ) N, the number of points.
!
!    Input, real ( kind = rkind ) R(N), the radial separation.
!    0 < R.
!
!    Input, real ( kind = rkind ) R0, a scale factor.
!
!    Output, real ( kind = rkind ) V(N), the value of the radial basis function.
!
  implicit none

  integer ( kind = ikind ) n

  integer ( kind = ikind ) i
  real ( kind = rkind ) r(n)
  real ( kind = rkind ) r0
  real ( kind = rkind ) v(n)

  do i = 1, n
    if ( r(i) .le. 0.0D+00 ) then
      v(i) = 0.0D+00
    else
      v(i) = r(i)**2 * log ( r(i) / r0 )
    end if
  end do

  return
end subroutine phi3

subroutine phi4 ( n, r, r0, v )

!*****************************************************************************80
!
!! PHI4 evaluates the gaussian radial basis function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer ( kind = ikind ) N, the number of points.
!
!    Input, real ( kind = rkind ) R(N), the radial separation.
!    0 < R.
!
!    Input, real ( kind = rkind ) R0, a scale factor.
!
!    Output, real ( kind = rkind ) V(N), the value of the radial basis function.
!
  implicit none

  integer ( kind = ikind ) n

  real ( kind = rkind ) r(n)
  real ( kind = rkind ) r0
  real ( kind = rkind ) v(n)

  v(1:n) = exp ( - 0.5D+00 * r(1:n)**2 / r0**2 )

  return
end subroutine phi4

end module rbf_interp_nd_mod
