module legendre_gauss_lobatto
   use, intrinsic :: iso_fortran_env, only : wp => real64
   real(wp), parameter :: pi = acos(-1._wp)
   contains
   !! elemental functions can take instead of a floating point : x, a vector x(:) of a desired size
   !! without explicitly invoking the size of the vector x

   !! Jacobi calculates the jacobi hypergeometric orthogonal function
   real(wp) pure elemental function Jacobi(N, alpha, beta, x)
      implicit none
      integer,  intent(in) :: N
      real(wp), intent(in) :: alpha, beta, x
      real(wp) :: Pk, Pk_1, Pk_2, reg
      integer  :: k

      Pk_2 = 1.0_wp
      if (N.eq.0) then
         Jacobi = Pk_2
         return
      endif
      Pk_1 = ((alpha - beta) + (alpha + beta + 2.0_wp) * x) / 2.0_wp
      if (N.eq.1) then
         Jacobi = Pk_1
         return
      endif

      do k = 2, N
         reg = 2._wp*real(k,wp) + alpha + beta
         Pk = (reg - 1._wp) * ( reg*(reg-2._wp)*x + alpha**2 - beta**2 )*Pk_1 - &
              2._wp*(real(k,wp)+alpha-1._wp)*(real(k,wp)+beta-1._wp)*reg*Pk_2
         Pk = Pk / ( 2._wp*k*(reg-k)*(reg-2._wp) )
         Pk_2 = Pk_1
         Pk_1 = Pk
      end do
      Jacobi = Pk
   end function Jacobi

   !! jacobi_normalized calculates the orthonormal polynomial basis
   pure elemental function jacobi_normalized(N, alpha, beta, x)
      implicit none
      integer,  intent(in) :: N
      real(wp), intent(in) :: alpha, beta, x
      real(wp) :: gamma_n, jacobi_normalized
      integer  :: i
      gamma_n = 2**(alpha+beta+1._wp) / (2*N+alpha+beta+1._wp)
      gamma_n = gamma_n * gamma(N+alpha+1._wp) * gamma(N+beta+1._wp)
      gamma_n = gamma_n / gamma(N+alpha+beta+1._wp)
      do i=N,1,-1
         gamma_n = gamma_n / real(i,wp)
      enddo
      jacobi_normalized = jacobi(N, alpha, beta, x) / sqrt(gamma_n)
   end function jacobi_normalized

   real(wp) pure elemental function grad_jacobi(N, alpha, beta, x)
      implicit none
      integer , intent(in) :: N
      real(wp), intent(in) :: alpha, beta, x
      if (N.eq.0) then
         grad_jacobi = 0._wp
      else
         grad_jacobi = gamma(alpha+beta+N+2) / (2*gamma(alpha+beta+N+1)) * &
                       Jacobi(N-1,alpha+1._wp, beta+1._wp, x)
      endif
   end function grad_jacobi

   real(wp) pure elemental function grad_jacobi_normalized(N, alpha, beta, x)
      implicit none
      integer , intent(in) :: N
      real(wp), intent(in) :: alpha, beta, x
      if (N.eq.0) then
         grad_jacobi_normalized = 0._wp
      else
         grad_jacobi_normalized = sqrt(N*(N+alpha+beta+1))* &
                                  jacobi_normalized(N-1, alpha+1._wp,beta+1._wp,x)
      endif
   end function grad_jacobi_normalized

   !! calculates the Legendre-Gauss-Lobatto cites, used in the quadrature
   !! input is the order of polynomial, the output is of size = order + 1 = N + 1
   pure function lgl_nodes(N)
      real(wp), parameter :: eps = 1.e-15_wp
      integer, intent(in) :: N
      real(wp) :: lgl_nodes(N+1)
      real(wp) :: P(N+1,N+1), x_old(N+1), x_new(N+1)
      integer  :: i
      x_new = cos(pi*linspace(N+1, 0._wp, real(N,wp)) / real(N,wp))
      P = 0._wp; x_old = 2._wp
      do while (maxval(abs(x_new-x_old)).gt.eps)
         x_old = x_new
         P(:,1) = 1._wp
         P(:,2) = x_new
         do i=3,N+1
            P(:,i) = ( (2*i-3)*x_new*P(:,i-1) - (i-2)*P(:,i-2) ) / (i-1)
         enddo
         x_new = x_old - ( x_new*P(:,N+1)-P(:,N) ) / ( (N+1)*P(:,N+1) );
      enddo
      lgl_nodes = x_new(N+1:1:-1)
   end function lgl_nodes

   !function gauss_legendre_nodes(N)
   !   real(wp), parameter :: eps = 1.e-15_wp
   !   integer, intent(in) :: N
   !   real(wp) :: gauss_legendre_nodes(N+1)
   !   real(wp) :: P, x_old, x_new(N+2), x, dP
   !   integer :: i
   !   do i=1,N+2
   !      x_new(i) = cos((4.d0*i-1.d0)*pi/(4.d0*(n+2)+2.d0))
   !      print *, x_new(i)
   !   enddo
   !   do i=1,N+1
   !      x = 0.5_wp*(x_new(i)+x_new(i+1))
   !      x_old = x_new(i)
   !      do while(abs(x-x_old).gt.eps)
   !         P  = jacobi(N,0._wp, 0._wp, x)
   !         dP = grad_jacobi(N,0._wp, 0._wp, x)
   !         x_old = x
   !         x = x - P / dP
   !      enddo
   !      gauss_legendre_nodes(i) = x
   !   enddo
   !   print *, ''
   !   do i=1,N+1
   !      print *, gauss_legendre_nodes(i)
   !   enddo
   !end function gauss_legendre_nodes

   !! calculates the Legendre-Gauss-Lobatto weights, used in the quadrature
   !! N :: is the degree of the corresponding polynomial
   !! For example : 
   !! N = 3 are the interpolation points for a polynomial of order = 2, needs N = 3 coefficients
   pure function lgl_weights(N)
      implicit none
      integer, intent(in) :: N
      integer  :: j
      real(wp) :: lgl_weights(N), x_lgl(N)
      x_lgl = lgl_nodes(N-1)
      lgl_weights(1) = -1._wp
      lgl_weights(N) = 1._wp
      do j=2,N-1
         lgl_weights(j) = 2._wp / (real(n*(n-1),wp) * Jacobi(N-1,0._wp, 0._wp, x_lgl(j))**2)
      enddo
   end function lgl_weights

   !! calculates numerical integration from N number of samples of f(x)
   !! integrates exactly a polynomial of order 2*N-3
   pure function lgl_quadrature(N, f, a, b)
      integer, intent(in) :: N
      real(wp), intent(in) :: f(N), a, b
      real(wp) :: lgl_quadrature, weights(N)
      weights = lgl_weights(N)
      lgl_quadrature = dot_product(weights(2:N-1), f(2:N-1)) + 2._wp/real(n,wp)/real(n-1,wp)*(f(1)+f(N))
      lgl_quadrature = 0.5_wp*(b-a)*lgl_quadrature
   end function lgl_quadrature

   !! auxiliary routine, needed for the initial guess of lgl nodes
   pure function linspace(N, a, b)
      implicit none
      integer , intent(in)  :: N
      real(wp), intent(in)  :: a, b
      real(wp) :: linspace(N), step, junk
      integer  :: i
      step = (b-a)/real((N-1),wp)
      junk = a-step
      do i=1,N
          junk = junk + step
          linspace(i) = junk
      enddo
   end function linspace

   function pade_reconstruction(N, u)
      implicit none
      integer, intent(in) :: N !! N >= 8 must
      real(wp), intent(in):: u(N+1)
      real(wp) :: f(N+1), pade_reconstruction(N+1), cites(N+1), gauss_nodes(N+1)
      real(wp), allocatable :: A(:,:), rhs(:), lhs(:,:), q_tilde(:), Q(:), p_tilde(:)
      integer :: M, L, i, j, cut
      pade_reconstruction = 0._wp
      cites = lgl_nodes(N)
      M = N - 5
      L = N - 6
      print *, N
      gauss_nodes = lgl_nodes(N)
      print *, gauss_nodes
      stop
      print *, M, L, N
      allocate(A(L,L+1))
      do i=M+1,M+L
         do j=0,L
            f = u * cites**j * jacobi_normalized(i, 0._wp, 0._wp, cites)
            A(i-M,j+1) = lgl_quadrature(N, f, -1._wp, 1._wp)
         enddo
      enddo
      allocate(rhs(L))
      lhs = -1e+2-1e+23
      cut = 3
      rhs = -A(:,cut)
      allocate(lhs(L,L))
      lhs(:,:cut-1) = A(:,:cut-1)
      lhs(:,cut:) = A(:,cut+1:)
      deallocate(A)
      call my_dgesv(L, lhs, rhs)
      allocate(q_tilde(L+1))
      q_tilde(1) = 1._wp
      q_tilde(2:L+1) = rhs
      deallocate(lhs, rhs)
      !print *, q_tilde
      allocate(Q(N+1))
      ! calculate Q
      do i=1,N ! for every calculation point
         Q(i) = 0._wp
         do j=0,L ! add every monomial of the expansion
            Q(i) = Q(i) + q_tilde(j+1) * cites(i)**j
         enddo
      enddo
      print *, Q
      allocate(p_tilde(M+1))
      print *, M, L, N
      error stop 'here'
      !do i=0,M
      !   f = Q * u * jacobi_normalized(i, 0._wp, 0._wp, cites)
      !   p_dilde(i+1) = lgl_quadrature(
      !enddo





   end function pade_reconstruction

   subroutine my_dgesv(N,A, B)
      implicit none
      external dgesv
      integer, intent(in)             :: N
      double precision, intent(in)    :: A(N,N)
      double precision, intent(inout) :: B(N)
      integer :: ipiv(n), info
      call DGESV( N, 1, A, N, ipiv, B, N, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The diagonal element of the triangular factor of A,'
         WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
         WRITE(*,*)'A is singular; the solution could not be computed.'
         STOP
      END IF
   endsubroutine my_dgesv
end module legendre_gauss_lobatto


