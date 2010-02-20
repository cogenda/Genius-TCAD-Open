subroutine getnp2 ( px, py, x, y, nr, lcell, lnext, xmin, ymin, &
  dx, dy, np, dsq )

!*****************************************************************************80
!
!! GETNP2 seeks the closest unmarked node to a point.
!
!  Discussion:
!
!    GETNP2 uses the cell method to find the closest unmarked node NP
!    to a specified point P, given a set of N nodes and the data structure 
!    defined by STORE2.
!
!    NP is then marked by negating LNEXT(NP).  Thus, the closest M nodes to
!    P may be determined by a sequence of M calls to this routine.  
!
!    If the point P is itself actually a node K, and you want to find the
!    nearest point to P that is not node K, then you must be sure to mark
!    node K before calling.
!
!    The search is begun in the cell containing or closest to P and proceeds 
!    outward in rectangular layers until all cells which contain points 
!    within distance R of P have been searched.  R is the distance from P to 
!    the first unmarked node encountered, or infinite if no unmarked nodes
!    are present.
!
!    Input parameters other than LNEXT are not altered by this routine.  
!    With the exception of ( PX, PY ) and the signs of LNEXT elements, 
!    these parameters should be unaltered from their values on output 
!    from subroutine STORE2.
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    Robert Renka,
!    University of North Texas
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 660: QSHEP2D, Quadratic Shepard method for bivariate
!    interpolation of scattered data,
!    ACM Transactions on Mathematical Software,
!    Volume 14, 1988, pages 149-150.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) PX, PY, the (X,Y) coordinates of the point P
!    whose nearest unmarked neighbor is to be found.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes at which
!    data has been supplied.
!
!    Input, integer NR, the number of rows and columns in the cell grid.
!    NR must be at least 1.
!
!    Input, integer LCELL(NR,NR), array of nodal indices associated
!    with cells.
!
!    Input/output, integer LNEXT(N), contains next-node indices ( or their 
!    negatives ).  On return, if the output value of NP is nonzero, then
!    LNEXT(NP) will be negative.
!
!    Input, real ( kind = 8 ) XMIN, YMIN, DX, DY, the minimum nodal X, Y
!    coordinates, and the X, Y dimensions of a cell.  DX and DY must be
!    positive.
!
!    Output, integer NP, the index into the vectors X and Y of the nearest
!    unmarked node to the point P.  NP will be 0 if all nodes are marked 
!    or if the values of NR, DX, DY are illegal.  LNEXT(NP) will be less
!    than 0 if NP is nonzero (this marks node NP as being used now).
!
!    Output, real ( kind = 8 ) DSQ, if NP is nonzero, then DSQ is the 
!    squared distance between P and node NP.
!
!  Local Parameters:
!
!    first = true iff the first unmarked node has yet to be encountered,
!
!    imin, imax, jmin, jmax = cell indices defining the range of the search,
!
!    delx, dely = px-xmin and py-ymin,
!
!    i0, j0 = cell containing or closest to P,
!
!    i1, i2, j1, j2 = cell indices of the layer whose intersection with 
!    the range defined by imin,...,jmax is currently being searched.
!
  implicit none

  integer nr

  real ( kind = 8 ) delx
  real ( kind = 8 ) dely
  real ( kind = 8 ) dsq
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  logical first
  integer i
  integer i0
  integer i1
  integer i2
  integer imax
  integer imin
  integer j
  integer j0
  integer j1
  integer j2
  integer jmax
  integer jmin
  integer l
  integer lcell(nr,nr)
  integer lmin
  integer ln
  integer lnext(*)
  integer np
  real ( kind = 8 ) px
  real ( kind = 8 ) py
  real ( kind = 8 ) r
  real ( kind = 8 ) rsmin
  real ( kind = 8 ) rsq
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xp
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yp

  xp = px
  yp = py
!
!  Test for invalid input parameters.
!
  if ( nr < 1 .or. dx <= 0.0D+00 .or. dy <= 0.0D+00 ) then
    np = 0
    dsq = 0.0D+00
  end if
!
!  Initialize parameters:
!
  first = .true.
  imin = 1
  imax = nr
  jmin = 1
  jmax = nr
  delx = xp - xmin
  dely = yp - ymin

  i0 = int ( delx / dx ) + 1
  i0 = max ( i0, 1 )
  i0 = min ( i0, nr )

  j0 = int ( dely / dy ) + 1
  j0 = max ( j0, 1 )
  j0 = min ( j0, nr )

  i1 = i0
  i2 = i0
  j1 = j0
  j2 = j0
!
!  Outer loop on layers, inner loop on layer cells, excluding
!  those outside the range (imin,imax) x (jmin,jmax).
!
1 continue

  do j = j1, j2

    if ( jmax < j ) then
      exit
    end if

    if ( j < jmin ) then
      cycle
    end if

    do i = i1, i2

      if ( imax < i ) then
        exit
      end if

      if ( i < imin ) then
        cycle
      end if

      if ( j /= j1 .and. j /= j2 .and. i /= i1 .and. i /= i2 ) then
        cycle
      end if
!
!  Search cell (i,j) for unmarked nodes l.
!
      l = lcell(i,j)

      if ( 0 < l ) then
!
!  Loop on nodes in cell (i,j).
!
2       continue

        ln = lnext(l)
!
!  Node L is the first unmarked neighbor of P encountered.
!
!  Initialize LMIN to the current candidate for np, and
!  rsmin to the squared distance from p to lmin.  imin,
!  imax, jmin, and jmax are updated to define the smal-
!  lest rectangle containing a circle of radius r =
!  sqrt(rsmin) centered at p, and contained in (1,nr) x
!  (1,nr) (except that, if p is outside the rectangle
!  defined by the nodes, it is possible that imin .gt.
!  nr, imax < 1, nr < jmin, or jmax < 1).
!
        if ( 0 <= ln ) then

          rsq = ( x(l) - xp )**2 + ( y(l) - yp )**2

          if ( first ) then

            lmin = l
            rsmin = rsq
            r = sqrt ( rsmin )

            imin = int ( ( delx - r ) / dx ) + 1
            imin = max ( imin, 1 )

            imax = int ( ( delx + r ) / dx ) + 1
            imax = min ( imax, nr )

            jmin = int ( ( dely - r ) / dy ) + 1
            jmin = max ( jmin, 1 )

            jmax = int ( ( dely + r ) / dy ) + 1
            jmax = min ( jmax, nr )

            first = .false.

          else

            if ( rsq < rsmin ) then
              lmin = l
              rsmin = rsq
            end if
 
          end if

        end if

        if ( abs ( ln ) /= l ) then
          l = abs ( ln )
          go to 2
        end if

      end if

    end do

  end do
!
!  Test for termination of loop on cell layers.
!
  if ( imin < i1 .or. i2 < imax .or. jmin < j1 .or. j2 < jmax ) then

    i1 = i1 - 1
    i2 = i2 + 1
    j1 = j1 - 1
    j2 = j2 + 1
    go to 1

  end if

  if ( first ) then
    np = 0
    dsq = 0.0D+00
  else
    np = lmin
    dsq = rsmin
    lnext(lmin) = -lnext(lmin)
  end if

  return
end
subroutine givens ( a, b, c, s )

!*****************************************************************************80
!
!! GIVENS constructs a Givens plane rotation.
!
!  Discussion:
!
!    The transformation has the form of a 2 by 2 matrix G(C,S):
!
!      (   C  S )
!      ( - S  C )
!
!    where C*C + S*S = 1, which zeroes the second entry of the
!    the column vector ( A, B ) when C and S are properly chosen.
!    A call to GIVENS is normally followed by a call to ROTATE
!    which computes the product of G(C,S) with a 2 by N matrix.
!
!  Modified:
!
!    10 July 1999
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A, B.
!
!    On input, A and B define the 2-vector whose second entry (B) is
!    to be annihilated by a Givens rotation.
!
!    On output, A has been overwritten by a value
!      R = +/- SQRT ( A*A + B*B )
!    and B has been overwritten by a value Z which allows C
!    and S to be recovered as:
!
!      if | Z | <= 1, then
!        C = SQRT ( 1 - Z*Z ), 
!        S = Z
!      else if 1 < | Z | then
!        C = 1 / Z, 
!        S = SQRT ( 1 - C*C ).
!
!    Output, real ( kind = 8 ) C, S, the components of the Givens
!    transformation, which may be computed by:
!
!      C = +/- A / SQRT ( A*A + B*B )
!      S = +/- B / SQRT ( A*A + B*B )
!
!  Local parameters:
!
!  r =        c*a + s*b = +/-sqrt(a*a+b*b)
!  u,v =   variables used to scale a and b for computing r
!
!  abs(b) < abs(a)
!
!  Note that r has the sign of a, 0 < c, and s has
!  sign(a)*sign(b).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) u
  real ( kind = 8 ) v

  if ( abs ( b ) < abs ( a ) ) then

    u = 2.0D+00 * a
    v = b / u
    r = sqrt ( 0.25D+00 + v * v ) * u
    c = a / r
    s = 2.0D+00 * v * c
    b = s
    a = r
!
!  abs(a) <= abs(b)
!
!  Store r in a.
!  Note that r has the sign of b, 0 < s, and c has sign(a)*sign(b).
!
  else if ( b /= 0.0D+00 ) then

    u = 2.0D+00 * b
    v = a / u
    a = sqrt ( 0.25D+00 + v * v ) * u
    s = b / a
    c = 2.0D+00 * v * s

    if ( c /= 0.0D+00 ) then
      b = 1.0D+00 / c
    else
      b = 1.0D+00
    end if
!
!  a = b = 0.
!
  else

    c = 1.0D+00
    s = 0.0D+00

  end if

  return
end
subroutine qs2grd ( px, py, n, x, y, f, nr, lcell, lnext, xmin, &
  ymin, dx, dy, rmax, rsq, a, q, qx, qy, ier )

!*****************************************************************************80
!
!! QS2GRD evaluates the interpolant and its first spatial derivatives.
!
!  Discussion:
!
!    QS2GRD computes the value and the gradient at the point (PX,PY) 
!    of the interpolatory function Q, defined by QSHEP2 for a given set
!    of scattered data.  Q(X,Y) is a weighted sum of quadratic
!    nodal functions.
!
!    Input parameters are not altered by this subroutine.  The parameters 
!    other than PX and PY should be input unaltered from their values 
!    on output from QSHEP2.  This subroutine should not be called if a 
!    nonzero error flag was returned by QSHEP2.
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    Robert Renka,
!    University of North Texas
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 660: QSHEP2D, Quadratic Shepard method for bivariate
!    interpolation of scattered data,
!    ACM Transactions on Mathematical Software,
!    Volume 14, 1988, pages 149-150.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) PX, PY, the coordinates of the point at which the
!    interpolant and its derivatives are to be evaluated.
!
!    Input, integer N, the number of nodes and data values which
!    are to be interpolated.  N must be at least 6. 
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes at which
!    data has been supplied.
!
!    Input, real ( kind = 8 ) F(N), the data values at the nodes.
!
!    Input, integer NR, the number of rows and columns in the cell 
!    grid.  Refer to subroutine STORE2 for details.  NR must be at least 1.
!
!    Input, integer LCELL(NR,NR), array of nodal indices associated
!    with cells.
!
!    Input, integer LNEXT(N), contains next-node indices.
!
!    Input, real ( kind = 8 ) XMIN, YMIN, DX, DY, the minimum nodal X, Y
!    coordinates, and the X, Y dimensions of a cell.  Computed by QSHEP2.
!
!    Input, real ( kind = 8 ) RMAX, the square root of the largest element
!    in RSQ, the maximum radius of influence.  Computed by QSHEP2.
!
!    Input, real ( kind = 8 ) RSQ(N), the squared radii which enter into the
!    weights defining the interpolant Q.  Computed by QSHEP2.
!
!    Input, real ( kind = 8 ) A(5,N), the coefficients for the nodal functions 
!    defining the interpolant Q.  Computed by QSHEP2.
!
!    Output, real ( kind = 8 ) Q, QX, QY, the value of the interpolant, and
!    its derivatives with respect to X and Y, at (PX,PY).
!
!    Output, integer IER, error indicator.
!    0, if no errors were encountered.
!    1, if N, NR, DX, DY or RMAX is invalid.
!    2, if no errors were encountered but (PX,PY) is not within the 
!       radius R(K) for any node K and thus Q = QX = QY = 0.
!
  implicit none

  integer n
  integer nr

  real ( kind = 8 ) a(5,n)
  real ( kind = 8 ) delx
  real ( kind = 8 ) dely
  real ( kind = 8 ) ds
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ) f(n)
  integer i
  integer ier
  integer imax
  integer imin
  integer j
  integer jmax
  integer jmin
  integer k
  integer kp
  integer lcell(nr,nr)
  integer lnext(n)
  real ( kind = 8 ) px
  real ( kind = 8 ) py
  real ( kind = 8 ) q
  real ( kind = 8 ) qk
  real ( kind = 8 ) qkx
  real ( kind = 8 ) qky
  real ( kind = 8 ) qx
  real ( kind = 8 ) qy
  real ( kind = 8 ) rd
  real ( kind = 8 ) rds
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rs
  real ( kind = 8 ) rsq(n)
  real ( kind = 8 ) sw
  real ( kind = 8 ) swq
  real ( kind = 8 ) swqx
  real ( kind = 8 ) swqy
  real ( kind = 8 ) sws
  real ( kind = 8 ) swx
  real ( kind = 8 ) swy
  real ( kind = 8 ) t
  real ( kind = 8 ) w
  real ( kind = 8 ) wx
  real ( kind = 8 ) wy
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xp
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yp

  xp = px
  yp = py

  if ( n < 6 ) then
    ier = 1
    return
  else if ( nr < 1 ) then
    ier = 1
    return
  else if ( dx <= 0.0D+00 ) then
    ier = 1
    return
  else if ( dy <= 0.0D+00 ) then
    ier = 1
    return
  else if ( rmax < 0.0D+00 ) then
    ier = 1
    return
  end if
!
!  Set imin, imax, jmin, and jmax to cell indices defining
!  the range of the search for nodes whose radii include P.
!  The cells which must be searched are those inter-
!  sected by (or contained in) a circle of radius rmax
!  centered at p.
!
  imin = int ( ( xp - xmin - rmax ) / dx ) + 1
  imin = max ( imin, 1 )

  imax = int ( ( xp - xmin + rmax ) / dx ) + 1
  imax = min ( imax, nr )

  jmin = int ( ( yp - ymin - rmax ) / dy ) + 1
  jmin = max ( jmin, 1 )

  jmax = int ( ( yp - ymin + rmax ) / dy ) + 1
  jmax = min ( jmax, nr )
!
!  Test for no cells within the circle of radius RMAX.
!
  if ( imax < imin .or. jmax < jmin ) then
    q = 0.0D+00
    qx = 0.0D+00
    qy = 0.0D+00
    ier = 2
    return
  end if
!
!  Q = swq/sw = sum(w(k)*q(k))/sum(w(k)) where the sum is
!  from k = 1 to n, q(k) is the quadratic nodal function,
!  and w(k) = ((r-d)+/(r*d))**2 for radius r(k) and distance d(k).  Thus
!
!    qx = (swqx*sw - swq*swx)/sw**2  and
!    qy = (swqy*sw - swq*swy)/sw**2
!
!  where swqx and swx are partial derivatives with respect
!  to x of swq and sw, respectively.  swqy and swy are 
!  defined similarly.
!
  sw = 0.0D+00
  swx = 0.0D+00
  swy = 0.0D+00
  swq = 0.0D+00
  swqx = 0.0D+00
  swqy = 0.0D+00
!
!  Outer loop on cells (I,J).
!
  do j = jmin, jmax

    do i = imin, imax

      k = lcell(i,j)
!
!  Inner loop on nodes K.
!
      if ( k /= 0 ) then

        do

          delx = xp - x(k)
          dely = yp - y(k)
          ds = delx * delx + dely * dely
          rs = rsq(k)

          if ( ds == 0.0D+00 ) then
            q = f(k)
            qx = a(4,k)
            qy = a(5,k)
            ier = 0
            return
          end if

          if ( ds < rs ) then

            rds = rs * ds
            rd = sqrt ( rds )
            w = ( rs + ds - rd - rd ) / rds
            t = 2.0D+00 * ( rd - rs ) / ( ds * rds )
            wx = delx * t
            wy = dely * t
            qkx = 2.0D+00 * a(1,k) * delx + a(2,k) * dely
            qky = a(2,k) * delx + 2.0D+00 * a(3,k) * dely
            qk = ( qkx * delx + qky * dely ) / 2.0D+00
            qkx = qkx + a(4,k)
            qky = qky + a(5,k)
            qk = qk + a(4,k) * delx + a(5,k) * dely + f(k)
            sw = sw + w
            swx = swx + wx
            swy = swy + wy
            swq = swq + w * qk
            swqx = swqx + wx * qk + w * qkx
            swqy = swqy + wy * qk + w * qky

          end if

          kp = k
          k = lnext(kp)

          if ( k == kp ) then
            exit
          end if

        end do

      end if

    end do

  end do
!
!  SW = 0 if and only if P is not within the radius R(K) for any node K.
!
  if ( sw /= 0.0D+00 ) then

    q = swq / sw
    sws = sw * sw
    qx = ( swqx * sw - swq * swx ) / sws
    qy = ( swqy * sw - swq * swy ) / sws
    ier = 0

  else

    q = 0.0D+00
    qx = 0.0D+00
    qy = 0.0D+00
    ier = 2

  end if

  return
end
subroutine qshep2 ( n, x, y, f, nq, nw, nr, lcell, lnext, xmin, &
  ymin, dx, dy, rmax, rsq, a, ier )

!*****************************************************************************80
!
!! QSHEP2 computes an interpolant to scattered data in the plane.
!
!  Discussion:
!
!    QSHEP2 computes a set of parameters A and RSQ defining a smooth, 
!    once continuously differentiable, bi-variate function Q(X,Y) which 
!    interpolates given data values F at scattered nodes (X,Y).  
!
!    The interpolant function Q(X,Y) may be evaluated at an arbitrary point 
!    by passing the parameters A and RSQ to the function QS2VAL.  The
!    first derivatives dQdX(X,Y) and dQdY(X,Y) may be evaluated by 
!    subroutine QS2GRD.
!
!    The interpolation scheme is a modified quadratic Shepard method:
!
!      Q = ( W(1) * Q(1) + W(2) * Q(2) + .. + W(N) * Q(N) ) 
!        / ( W(1)        + W(2)        + .. + W(N) )
!
!    for bivariate functions W(K) and Q(K).  The nodal functions are given by
!
!      Q(K)(X,Y) = 
!          F(K)
!        + A(1,K) * ( X - X(K) )**2 
!        + A(2,K) * ( X - X(K) ) * ( Y - Y(K) )
!        + A(3,K) * ( Y - Y(K) )**2 
!        + A(4,K) * ( X - X(K) )
!        + A(5,K) * ( Y - Y(K) ).
!
!    Thus, Q(K) is a quadratic function which interpolates the
!    data value at node K.  Its coefficients A(*,K) are obtained
!    by a weighted least squares fit to the closest NQ data
!    points with weights similar to W(K).  Note that the radius
!    of influence for the least squares fit is fixed for each
!    K, but varies with K.
!
!    The weights are taken to be
!
!      W(K)(X,Y) = ( (R(K)-D(K))+ / R(K) * D(K) )**2
!
!    where (R(K)-D(K))+ = 0 if R(K) <= D(K) and D(K)(X,Y) is
!    the euclidean distance between (X,Y) and (X(K),Y(K)).  The
!    radius of influence R(K) varies with K and is chosen so
!    that NW nodes are within the radius.  Note that W(K) is
!    not defined at node (X(K),Y(K)), but Q(X,Y) has limit F(K)
!    as (X,Y) approaches (X(K),Y(K)).
!
!  Author:
!
!    Robert Renka,
!    University of North Texas
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 660: QSHEP2D, Quadratic Shepard method for bivariate
!    interpolation of scattered data,
!    ACM Transactions on Mathematical Software,
!    Volume 14, 1988, pages 149-150.
!
!  Parameters:
!
!    Input, integer N, the number of nodes (X,Y) at which data values
!    are given.  N must be at least 6.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes at which
!    data has been supplied.
!
!    Input, real ( kind = 8 ) F(N), the data values.
!
!    Input, integer NQ, the number of data points to be used in the least
!    squares fit for coefficients defining the nodal functions Q(K).  
!    A highly recommended value is NQ = 13.  
!    NQ must be at least 5, and no greater than the minimum of 40 and N-1.
!
!    Input, integer NW, the number of nodes within (and defining) the radii
!    of influence R(K) which enter into the weights W(K).  For N 
!    sufficiently large, a recommended value is NW = 19.   NW must be
!    at least 1, and no greater than the minimum of 40 and N-1.
!
!    Input, integer NR, the number of rows and columns in the cell grid 
!    defined in subroutine STORE2.  A rectangle containing the nodes 
!    is partitioned into cells in order to increase search efficiency.  
!    NR = SQRT(N/3) is recommended.  NR must be at least 1.
!
!    Output, integer LCELL(NR,NR), array of nodal indices associated
!    with cells.
!
!    Output, integer LNEXT(N), contains next-node indices ( or their 
!    negatives ).
!
!    Output, real ( kind = 8 ) XMIN, YMIN, DX, DY, the minimum nodal X, Y
!    coordinates, and the X, Y dimensions of a cell.
!
!    Output, real ( kind = 8 ) RMAX, the square root of the largest element
!    in RSQ, the maximum radius of influence.
!
!    Output, real ( kind = 8 ) RSQ(N), the squared radii which enter into 
!    the weights defining the interpolant Q.
!
!    Output, real ( kind = 8 ) A(5,N), the coefficients for the nodal functions 
!    defining the interpolant Q.
!
!    Output, integer IER, error indicator.
!    0, if no errors were encountered.
!    1, if N, NQ, NW, or NR is out of range.
!    2, if duplicate nodes were encountered.
!    3, if all nodes are collinear.
!
!  Local parameters:
!
! av =        root-mean-square distance between k and the
!             nodes in the least squares system (unless
!             additional nodes are introduced for stabil-
!             ity).      the first 3 columns of the matrix
!             are scaled by 1/avsq, the last 2 by 1/av
! avsq =      av*av
! b =         transpose of the augmented regression matrix
! c =         first component of the plane rotation used to
!             zero the lower triangle of b**t -- computed
!             by subroutine givens
! ddx,ddy =   local variables for dx and dy
! dmin =      minimum of the magnitudes of the diagonal
!             elements of the regression matrix after
!             zeros are introduced below the diagonal
! DTOL =      tolerance for detecting an ill-conditioned
!             system.  the system is accepted when DTOL <= DMIN.
! fk =        data value at node k -- f(k)
! i =         index for a, b, and npts
! ib =        do-loop index for back solve
! ierr =      error flag for the call to store2
! irow =      row index for b
! j =         index for a and b
! jp1 =       j+1
! k =         nodal function index and column index for a
! lmax =      maximum number of npts elements (must be con-
!             sistent with the dimension statement above)
! lnp =       current length of npts
! neq =       number of equations in the least squares fit
! nn,nnr =    local copies of n and nr
! np =        npts element
! npts =      array containing the indices of a sequence of
!             nodes to be used in the least squares fit
!             or to compute rsq.  the nodes are ordered
!             by distance from k and the last element
!             (usually indexed by lnp) is used only to
!             determine rq, or rsq(k) if NQ < NW.
! nqwmax =    max(nq,nw)
! rq =        radius of influence which enters into the
!             weights for q(k) (see subroutine setup2)
! rs =        squared distance between k and npts(lnp) --
!             used to compute rq and rsq(k)
! rsmx =      maximum rsq element encountered
! rsold =     squared distance between k and npts(lnp-1) --
!             used to compute a relative change in rs
!             between succeeding npts elements
! RTOL =      tolerance for detecting a sufficiently large
!             relative change in rs.  if the change is
!             not greater than RTOL, the nodes are
!             treated as being the same distance from k
! rws =       current value of rsq(k)
! s =         second component of the plane givens rotation
! SF =        marquardt stabilization factor used to damp
!             out the first 3 solution components (second
!             partials of the quadratic) when the system
!             is ill-conditioned.  as SF increases, the
!             fitting function approaches a linear
! sum2 =      sum of squared euclidean distances between
!             node k and the nodes used in the least
!             squares fit (unless additional nodes are
!             added for stability)
! t =         temporary variable for accumulating a scalar
!             product in the back solve
! xk,yk =     coordinates of node k -- x(k), y(k)
! xmn,ymn =   local variables for xmin and ymin
!
  implicit none

  integer n
  integer nr

  real ( kind = 8 ) a(5,n)
  real ( kind = 8 ) av
  real ( kind = 8 ) avsq
  real ( kind = 8 ) b(6,6)
  real ( kind = 8 ) c
  real ( kind = 8 ) ddx
  real ( kind = 8 ) ddy
  real ( kind = 8 ) dmin
  real ( kind = 8 ), parameter :: dtol = 0.01D+00
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) fk
  integer i
  integer ier
  integer ierr
  integer irow
  integer j
  integer jp1
  integer k
  integer lcell(nr,nr)
  integer lmax
  integer lnext(n)
  integer lnp
  integer neq
  integer nn
  integer nnr
  integer np
  integer npts(40)
  integer nq
  integer nqwmax
  integer nw
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rq
  real ( kind = 8 ) rs
  real ( kind = 8 ) rsmx
  real ( kind = 8 ) rsold
  real ( kind = 8 ) rsq(n)
  real ( kind = 8 ), parameter :: rtol = 1.0D-05
  real ( kind = 8 ) rws
  real ( kind = 8 ) s
  real ( kind = 8 ), parameter :: SF = 1.0D+00
  real ( kind = 8 ) sum2
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xk
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xmn
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yk
  real ( kind = 8 ) ymin
  real ( kind = 8 ) ymn

  nn = n
  nnr = nr
  nqwmax = max ( nq, nw )
  lmax = min ( 40, n-1 )

  if ( nq < 5 ) then
    ier = 1
    return
  else if ( nw < 1 ) then
    ier = 1
    return
  else if ( lmax < nqwmax ) then
    ier = 1
    return
  else if ( nr < 1 ) then
    ier = 1
    return
  end if
!
!  Create the cell data structure, and initialize RSMX.
!
  call store2 ( nn, x, y, nnr, lcell, lnext, xmn, ymn, ddx, ddy, ierr )

  if ( ierr /= 0 ) then
    xmin = xmn
    ymin = ymn
    dx = ddx
    dy = ddy
    ier = 3
    return
  end if

  rsmx = 0.0D+00
!
!  Outer loop on node K.
!
  do k = 1, nn

    xk = x(k)
    yk = y(k)
    fk = f(k)
!
!  Mark node K to exclude it from the search for nearest neighbors.
!
    lnext(k) = -lnext(k)
!
!  Initialize for loop on NPTS.
!
    rs = 0.0D+00
    sum2 = 0.0D+00
    rws = 0.0D+00
    rq = 0.0D+00
    lnp = 0
!
!  Compute NPTS, LNP, RWS, NEQ, RQ, and AVSQ.
!
1   continue

    sum2 = sum2 + rs

    if ( lnp == lmax ) then
      go to 3
    end if

    lnp = lnp + 1
    rsold = rs

    call getnp2 ( xk, yk, x, y, nnr, lcell, lnext, xmn, ymn, ddx, ddy, np, rs )

    if ( rs == 0.0D+00 ) then
      ier = 2
      return
    end if

    npts(lnp) = np

    if ( ( rs - rsold ) / rs < RTOL ) then
      go to 1
    end if

    if ( rws == 0.0D+00 .and. nw < lnp ) then
      rws = rs
    end if
!
!  RQ = 0 (not yet computed) and NQ < lnp.     
!
!  RQ = sqrt(RS) is sufficiently large to (strictly) include NQ nodes.  
!
!  The least squares fit will include NEQ = LNP - 1 equations for 
!  5 <= NQ <= NEQ < LMAX <= N-1.
!
    if ( rq == 0.0D+00 .and. nq < lnp ) then
      neq = lnp - 1
      rq = sqrt ( rs )
      avsq = sum2 / real ( neq, kind = 8 )
    end if

    if ( nqwmax < lnp ) then
      go to 4
    else
      go to 1
    end if
!
!  All LMAX nodes are included in NPTS.   RWS and/or RQ**2 is
!  (arbitrarily) taken to be 10 percent larger than the
!  distance RS to the last node included.
!
3   continue

    if ( rws == 0.0D+00 ) then
      rws = 1.1D+00 * rs
    end if

    if ( rq == 0.0D+00 ) then
      neq = lmax
      rq = sqrt ( 1.1D+00 * rs )
      avsq = sum2 / real ( neq, kind = 8 )
    end if

4   continue
!
!  Store RSQ(K), update RSMX if necessary, and compute AV.
!
    rsq(k) = rws
    rsmx = max ( rsmx, rws )
    av = sqrt ( avsq )
!
!  Set up the augmented regression matrix (transposed) as the
!  columns of B, and zero out the lower triangle (upper
!  triangle of B) with Givens rotations -- QR decomposition
!  with orthogonal matrix Q not stored.
!
    i = 0

5   continue

    i = i + 1
    np = npts(i)
    irow = min ( i, 6 )

    call setup2 ( xk, yk, fk, x(np), y(np), f(np), av, avsq, rq, b(1,irow) )

    if ( i == 1 ) then
      go to 5
    end if

    do j = 1, irow-1
      jp1 = j + 1
      call givens ( b(j,j), b(j,irow), c, s )
      call rotate ( 6-j, c, s, b(jp1,j), b(jp1,irow) )
    end do

    if ( i < neq ) then
      go to 5
    end if
!
!  Test the system for ill-conditioning.
!
    dmin =  min ( abs ( b(1,1) ), abs ( b(2,2) ), abs ( b(3,3) ), &
      abs ( b(4,4) ), abs ( b(5,5) ) )

    if ( DTOL <= dmin * rq ) then
      go to 13
    end if

    if ( neq == lmax ) then
      go to 10
    end if
!
!  Increase RQ and add another equation to the system to improve conditioning.  
!  The number of NPTS elements is also increased if necessary.
!
7   continue

    rsold = rs
    neq = neq + 1

    if ( neq == lmax ) then
      go to 9
    end if
!
!   NEQ < LNP.
!
    if ( neq /= lnp ) then
      np = npts(neq+1)
      rs = ( x(np) - xk )**2 + ( y(np) - yk )**2
      if ( ( rs - rsold ) / rs < rtol ) then
        go to 7
      end if
      rq = sqrt(rs)
      go to 5
    end if
!
!  Add an element to NPTS.
!
    lnp = lnp + 1
    call getnp2 ( xk, yk, x, y, nnr, lcell, lnext, xmn, ymn, ddx, ddy, np, rs )

    if ( np == 0 ) then
      ier = 2
      return
    end if

    npts(lnp) = np

    if ( ( rs - rsold ) / rs < rtol ) then
      go to 7
    end if

    rq = sqrt ( rs )
    go to 5

9   continue

    rq = sqrt ( 1.1D+00 * rs )
    go to 5
!
!  Stabilize the system by damping second partials.  Add multiples of the 
!  first three unit vectors to the first three equations.
!
10  continue

    do i = 1, 3

      b(i,6) = sf
      b(i+1:6,6) = 0.0D+00

      do j = i, 5
        jp1 = j + 1
        call givens ( b(j,j), b(j,6), c, s )
        call rotate ( 6-j, c, s, b(jp1,j), b(jp1,6) )
      end do

    end do
!
!  Test the stabilized system for ill-conditioning.
!
    dmin = min ( abs ( b(1,1) ), abs ( b(2,2) ), abs ( b(3,3) ), &
      abs ( b(4,4) ), abs ( b(5,5) ) )

    if ( dmin * rq < dtol ) then
      xmin = xmn
      ymin = ymn
      dx = ddx
      dy = ddy
      ier = 3
      return
    end if
!
!  Solve the 5 by 5 triangular system for the coefficients.
!
13  continue

    do i = 5, 1, -1

      t = 0.0D+00

      do j = i+1, 5
        t = t + b(j,i) * a(j,k)
      end do

      a(i,k) = ( b(6,i) - t ) / b(i,i)

    end do
!
!  Scale the coefficients to adjust for the column scaling.
!
    a(1:3,k) = a(1:3,k) / avsq
    a(4,k) = a(4,k) / av
    a(5,k) = a(5,k) / av
!
!  Unmark K and the elements of NPTS.
!
    lnext(k) = - lnext(k)

    do i = 1, lnp
      np = npts(i)
      lnext(np) = - lnext(np)
    end do

  end do
!
!  No errors encountered.
!
  xmin = xmn
  ymin = ymn
  dx = ddx
  dy = ddy
  rmax = sqrt ( rsmx )
  ier = 0

  return
end
function qs2val ( px, py, n, x, y, f, nr, lcell, lnext, xmin, &
  ymin, dx, dy, rmax, rsq, a )

!*****************************************************************************80
!
!! QS2VAL evaluates the interpolant function at a point.
!
!  Discussion:
!
!    QS2VAL returns the value Q(PX,PY) where Q is the weighted sum of 
!    quadratic nodal functions defined by QSHEP2.  If the spatial 
!    derivatives of Q are also desired, call QS2GRD instead.
!
!    Input parameters are not altered by this function.  The
!    parameters other than PX and PY should be input unaltered
!    from their values on output from QSHEP2.  This function
!    should not be called if a nonzero error flag was returned
!    by QSHEP2.
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    Robert Renka,
!    University of North Texas
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 660: QSHEP2D, Quadratic Shepard method for bivariate
!    interpolation of scattered data,
!    ACM Transactions on Mathematical Software,
!    Volume 14, 1988, pages 149-150.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) PX, PY, the (X,Y) coordinates of the point P at
!    which Q is to be evaluated.
!
!    Input, integer N, the number of nodes and data values to be 
!    interpolated.  N must be at least 6.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes at which
!    data has been supplied.
!
!    Input, real ( kind = 8 ) F(N), the data values at the nodes.
!
!    Input, integer NR, the number of rows and columns in the cell grid.
!    Refer to subroutine STORE2.  NR must be at least 1.
!
!    Input, integer LCELL(NR,NR), the array of nodal indices associated
!    with cells.  Refer to STORE2.
!
!    Input, integer LNEXT(N), the next-node indices.  Refer to STORE2.
!
!    Input, real ( kind = 8 ) XMIN, YMIN, DX, DY, the minimum nodal X, Y
!    coordinates, and the X, Y dimensions of a cell.  Computed by QSHEP2.
!
!    Input, real ( kind = 8 ) RMAX, the square root of the largest element
!    in RSQ, the maximum radius of influence.  Computed by QSHEP2.
!
!    Input, real ( kind = 8 ) RSQ(N), the squared radii which enter into the
!    weights defining the interpolant Q.  Computed by QSHEP2.
!
!    Input, real ( kind = 8 ) A(5,N), the coefficients for the nodal functions 
!    defining the interpolant Q.  Computed by QSHEP2.
!
!    Output, real ( kind = 8 ) QS2VAL, the interpolated function value
!    at (PX,PY).
!
  implicit none

  integer n
  integer nr

  real ( kind = 8 ) a(5,n)
  real ( kind = 8 ) delx
  real ( kind = 8 ) dely
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ) f(n)
  integer i
  integer imax
  integer imin
  integer j
  integer jmax
  integer jmin
  real ( kind = 8 ) ds
  integer k
  integer kp
  integer lcell(nr,nr)
  integer lnext(n)
  real ( kind = 8 ) px
  real ( kind = 8 ) py
  real ( kind = 8 ) qs2val
  real ( kind = 8 ) rd
  real ( kind = 8 ) rds
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rs
  real ( kind = 8 ) rsq(n)
  real ( kind = 8 ) sw
  real ( kind = 8 ) swq
  real ( kind = 8 ) w
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xp
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yp

  xp = px
  yp = py
  qs2val = 0.0D+00

  if ( n < 6  ) then
    return  
  else if ( nr < 1  ) then
    return
  else if ( dx <= 0.0D+00 ) then
    return
  else if ( dy <= 0.0D+00 ) then
    return
  else if ( rmax < 0.0D+00 ) then
    return
  end if
!
!  Set imin, imax, jmin, and jmax to cell indices defining
!  the range of the search for nodes whose radii include
!  p.  The cells which must be searched are those intersected
!  by (or contained in) a circle of radius rmax
!  centered at P.
!
  imin = int ( ( xp - xmin - rmax ) / dx ) + 1
  imin = max ( imin, 1 )

  imax = int ( ( xp - xmin + rmax ) / dx ) + 1
  imax = min ( imax, nr )

  jmin = int ( ( yp - ymin - rmax ) / dy ) + 1
  jmin = max ( jmin, 1 )

  jmax = int ( ( yp - ymin + rmax ) / dy ) + 1
  jmax = min ( jmax, nr )
!
!  Test for no cells within the circle of radius RMAX.
!
  if ( imax < imin .or. jmax < jmin ) then
    qs2val = 0.0D+00
    return
  end if
!
!  Accumulate weight values in SW and weighted nodal function
!  values in swq.  the weights are w(k) = ((r-d)+/(r*d))**2
!  for r**2 = rsq(k) and d = distance between p and node K.
!
  sw = 0.0D+00
  swq = 0.0D+00

  do j = jmin, jmax

    do i = imin, imax

      k = lcell(i,j)

      if ( k /= 0 ) then

        do

          delx = xp - x(k)
          dely = yp - y(k)
          ds = delx * delx + dely * dely
          rs = rsq(k)

          if ( ds < rs ) then

            if ( ds == 0.0D+00 ) then
              qs2val = f(k)
              return
            end if

            rds = rs * ds
            rd = sqrt ( rds )
            w = ( rs + ds - rd - rd ) / rds
            sw = sw + w

            swq = swq + w * ( f(k) + a(1,k) * delx * delx &
              + a(2,k) * delx * dely + a(3,k) * dely * dely &
              + a(4,k) * delx + a(5,k) * dely )

          end if

          kp = k
          k = lnext(kp)

          if ( k == kp ) then
            exit
          end if

        end do

      end if

    end do

  end do
!
!  SW = 0 if and only if P is not within the radius R(K) for any node K.
!
  if ( sw == 0.0D+00 ) then
    qs2val = 0.0D+00
  else
    qs2val = swq / sw
  end if

  return
end
subroutine rotate ( n, c, s, x, y )

!*****************************************************************************80
!
!! ROTATE applies a Givens rotation.
!
!  Discussion:
!
!    The rotation has the form:
!
!      (   C  S )
!      ( - S  C )
!
!    and is essentially applied to a 2 by N matrix:
!
!      ( X(1) X(2) ... X(N) )
!      ( Y(1) Y(2) ... Y(N) )
!
!  Modified:
!
!    28 June 1999
!
!  Parameters:
!
!    Input, integer N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) C, S, the cosine and sine entries of the Givens
!    rotation matrix.  These may be determined by subroutine GIVENS.
!
!    Input/output, real ( kind = 8 ) X(N), Y(N), the rotated vectors. 
!
  implicit none

  integer n

  real ( kind = 8 ) c
  integer i
  real ( kind = 8 ) s
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xi
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yi

  if ( n <= 0 ) then
    return
  else if ( c == 1.0D+00 .and. s == 0.0D+00 ) then
    return
  end if

  do i = 1, n
    xi = x(i)
    yi = y(i)
    x(i) =   c * xi + s * yi
    y(i) = - s * xi + c * yi
  end do

  return
end
subroutine setup2 ( xk, yk, fk, xi, yi, fi, s1, s2, r, row )

!*****************************************************************************80
!
!! SETUP2 sets up a row of the least squares regression matrix.
!
!  Discussion:
!
!    SETUP2 sets up the I-th row of an augmented regression matrix for 
!    a weighted least-squares fit of a quadratic function Q(X,Y) to a set 
!    of data values F, where Q(XK,YK) = FK.  
!
!    The first 3 columns are quadratic terms, and are scaled by 1/S2.
!    The fourth and fifth columns represent linear terms, and are scaled 
!    by 1/S1.  
!
!    If D = 0, or R <= D, the weight is
!      0,
!    else if D < R, the weight is 
!      (R-D)/(R*D), 
!    where D is the distance between nodes I and K, and R is a maximum
!    influence distance.
!
!  Modified:
!
!    05 July 1999
!
!  Author:
!
!    Robert Renka,
!    University of North Texas
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 660: QSHEP2D, Quadratic Shepard method for bivariate
!    interpolation of scattered data,
!    ACM Transactions on Mathematical Software,
!    Volume 14, 1988, pages 149-150.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XK, YK, FK, the coordinates and value of the data
!    at data node K.
!
!    Input, real ( kind = 8 ) XI, YI, FI, the coorindates and value of the data
!    at data node I.
!
!    Input, real ( kind = 8 ) S1, S2, reciprocals of the scale factors.
!
!    Input, real ( kind = 8 ) R, the maximum radius of influence about node K.
!
!    Output, real ( kind = 8 ) ROW(6), a row of the augmented regression matrix.
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ) fi
  real ( kind = 8 ) fk
  integer i
  real ( kind = 8 ) r
  real ( kind = 8 ) row(6)
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) w
  real ( kind = 8 ) xi
  real ( kind = 8 ) xk
  real ( kind = 8 ) yi
  real ( kind = 8 ) yk

  dx = xi - xk
  dy = yi - yk

  d = sqrt ( dx * dx + dy * dy )

  if ( d <= 0.0D+00 .or. r <= d ) then

    row(1:6) = 0.0D+00

  else

    w = ( r - d ) / r / d

    row(1) = dx * dx * w / s2
    row(2) = dx * dy * w / s2
    row(3) = dy * dy * w / s2
    row(4) = dx * w / s1
    row(5) = dy * w / s1
    row(6) = ( fi - fk ) * w

  end if

  return
end
subroutine store2 ( n, x, y, nr, lcell, lnext, xmin, ymin, dx, dy, ier )

!*****************************************************************************80
!
!! STORE2 creates a cell data structure for the scattered data.
!
!  Discussion:
!
!    STORE2 is given a set of N arbitrarily distributed nodes in the 
!    plane and creates a data structure for a cell-based method of 
!    solving closest-point problems.  The smallest rectangle containing 
!    all the nodes is partitioned into an NR by NR uniform grid of cells, 
!    and nodes are associated with cells.      
!
!    In particular, the data structure stores the indices of the nodes 
!    contained in each cell.  For a uniform random distribution of nodes, 
!    the nearest node to an arbitrary point can be determined in constant
!    expected time.
!
!  Modified:
!
!    05 July 1999
!
!  Author:
!
!    Robert Renka
!    University of North Texas
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 660: QSHEP2D, Quadratic Shepard method for bivariate
!    interpolation of scattered data,
!    ACM Transactions on Mathematical Software,
!    Volume 14, 1988, pages 149-150.
!
!  Parameters:
!
!    Input, integer N, the number of data nodes.  N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the data nodes.
!
!    Input, integer NR, the number of rows and columns in the grid.  The
!    cell density, or average number of data nodes per cell, is
!      D = N / ( NR * NR ).
!    A recommended value, based on empirical evidence, is 
!      D = 3. 
!    Hence, the corresponding value of NR is recommended to be about
!      NR = SQRT ( N / 3 ).  
!    NR must be at least 1.
!
!    Output, integer LCELL(NR,NR), an array set up so that LCELL(I,J)
!    contains the index (for X and Y) of the first data node (that is, the
!    data node with smallest index) in the (I,J) cell.  LCELL(I,J) will be 0 if 
!    no data nodes are contained in the (I,J) cell.  The upper right corner of 
!    the (I,J) cell has coordinates 
!      ( XMIN + I * DX, YMIN + J * DY ).
!
!    Output, integer LNEXT(N), an array of next-node indices.  LNEXT(K)
!    contains the index of the next node in the cell which contains node K, 
!    or LNEXT(K) = K if K is the last node in the cell.
!    The data nodes contained in a cell are ordered by their indices.
!    If, for example, cell (I,J) contains nodes 2, 3, and 5 and no others, 
!    then:
!
!      LCELL(I,J) = 2, (index of the first data node)
!
!      LNEXT(2) = 3, 
!      LNEXT(3) = 5,
!      LNEXT(5) = 5.
!
!    Output, real ( kind = 8 ) XMIN, YMIN, the X, Y coordinates of the lower
!    left corner of the rectangle defined by the data nodes.  The upper right 
!    corner is ( XMAX, YMAX ), where
!      XMAX = XMIN + NR * DX,
!      YMAX = YMIN + NR * DY.
!
!    Output, real ( kind = 8 ) DX, DY, the X and Y dimensions of the 
!    individual cells.
!      DX = ( XMAX - XMIN ) / NR
!      DY = ( YMAX - YMIN ) / NR,
!    where XMIN, XMAX, YMIN and YMAX are the extrema of X and Y.
!
!    Output, integer IER, an error indicator.
!    0, if no errors were encountered.
!    1, if N < 2 or NR < 1.
!    2, if DX = 0 or DY = 0.
!
  implicit none

  integer n
  integer nr

  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer i
  integer ier
  integer j
  integer k
  integer l
  integer lcell(nr,nr)
  integer lnext(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  ier = 0

  if ( n < 2 ) then
    ier = 1
    return
  end if

  if ( nr < 1 ) then
    ier = 1
    return
  end if
!
!  Compute the dimensions of the (X,Y) rectangle containing all the data nodes.
!
  xmin = minval ( x(1:n) )
  xmax = maxval ( x(1:n) )
  ymin = minval ( y(1:n) )
  ymax = maxval ( y(1:n) )
!
!  Compute the dimensions of a single cell.
!
  dx = ( xmax - xmin ) / real ( nr, kind = 8 )
  dy = ( ymax - ymin ) / real ( nr, kind = 8 )
!
!  Test for zero area.
!
  if ( dx == 0.0D+00 .or. dy == 0.0D+00 ) then
    ier = 2
    return
  end if
!
!  Initialize LCELL.
!
  lcell(1:nr,1:nr) = 0
!
!  Loop on nodes, storing indices in LCELL and LNEXT.
!
  do k = n, 1, -1

    i = int ( ( x(k) - xmin ) / dx ) + 1
    i = min ( i, nr )

    j = int ( ( y(k) - ymin ) / dy ) + 1
    j = min ( j, nr )

    l = lcell(i,j)

    if ( l /= 0 ) then
      lnext(k) = l
    else
      lnext(k) = k
    end if

    lcell(i,j) = k

  end do

  return
end
