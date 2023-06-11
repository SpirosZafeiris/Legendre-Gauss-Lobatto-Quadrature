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

   !! calculates numerical integration from N number of samples of f(x)
   !! integrates exactly a polynomial of order 2*N-3
   pure function LGL_quadrature(N, f, a, b)
      integer, intent(in) :: N
      real(wp), intent(in) :: f(N), a, b
      real(wp) :: lgl_quadrature, weights(N)
      real(wp) :: buf(N,2)
      buf = LGL_nodes_weights(N)
      weights = buf(:,2)
      lgl_quadrature = dot_product(weights(2:N-1), f(2:N-1)) + 2._wp/real(n,wp)/real(n-1,wp)*(f(1)+f(N))
      lgl_quadrature = 0.5_wp*(b-a)*lgl_quadrature
   end function LGL_quadrature

   pure function GQ_nodes_weights(N)
      real(wp), parameter :: eps = 1.e-15_wp
      integer, intent(in) :: N
      real(wp) :: GQ_nodes_weights(N,2)
      integer  :: size_k, counter, kk
      integer, allocatable :: k(:)
      real(wp), allocatable :: theta(:), x(:), Pm2(:), Pm1(:), P(:), PP(:), PPm1(:), PPm2(:), dx(:)
      size_k = ceiling(real(N,wp)/real(2,wp))
      allocate(k(size_k), theta(size_k))
      k = [(n+mod(N,2))/2:1:-1]
      theta = pi*real(4*k-1,wp)/real(4*n+2,wp)
      allocate(x(size_k), Pm2(size_k), Pm1(size_k), P(size_k), PP(size_k), PPm1(size_k), PPm2(size_k), dx(size_k))
      x = (1._wp - real(N-1,wp)/real(8*N**3,wp) - 1._wp/real(384*N**4,wp)*(39._wp-28._wp/sin(theta)**2) ) * cos(theta)
      ! Initialise
      Pm2 = 1._wp
      Pm1 = x
      PPm2 = 0._wp
      PPm1 = 1._wp
      dx = 1.e+14_wp
      counter = 0

      ! Loop until convergence
      do while ( maxval(abs(dx)) > eps .and. counter < 10 )
         counter = counter + 1
         do kk = 1,N-1
            P = ((2*kk+1)*Pm1*x-kk*Pm2)/(kk+1)
            Pm2 = Pm1
            Pm1 = P
            PP = ((2*kk+1)*(Pm2+x*PPm1)-kk*PPm2)/(kk+1)
            PPm2 = PPm1
            PPm1 = PP
         enddo
         ! Newton step
         dx = -P/PP;
         ! Newton update
         x = x + dx;
         ! Reinitialise
         Pm2 = 1._wp
         Pm1 = x
         PPm2 = 0._wp
         PPm1 = 1._wp
      end do
   !
   ! Once more for derivatives
   do kk = 1,N-1
      P = ( (2*kk+1)*Pm1*x - kk*Pm2 ) / real(kk+1,wp)
      Pm2 = Pm1
      Pm1 = P
      PP = ( (2*kk+1)*(Pm2+x*PPm1) - kk*PPm2 ) / real(kk+1,wp)
      PPm2 = PPm1
      PPm1 = PP
   enddo

   GQ_nodes_weights(:,1) = [-x(size_k:1+mod(N,2):-1) , x] !! quadrature nodes
   GQ_nodes_weights(:,2) = [PP(size(PP):1+mod(N,2):-1) , PP]
   GQ_nodes_weights(:,2) = 2._wp/((1._wp-GQ_nodes_weights(:,1)**2)*GQ_nodes_weights(:,2)**2) !! quadrature weights
   end function GQ_nodes_weights

   pure function jacpts(N, alpha, beta)
      implicit none
      integer, intent(in) :: N
      real(wp), intent(in) :: alpha, beta
      real(wp) :: xw1(N,2), xw2(N,2), jacpts(N,2)
      integer :: idx(N)
      xw1 = jacpts_main(N, alpha, beta, .true.)
      xw2 = jacpts_main(N, beta, alpha, .false.)
      xw1(:,1) = xw1(N:1:-1,1)
      xw1(:,2) = xw1(N:1:-1,2)
      jacpts = 0._wp
      jacpts(:,1) = xw1(:,1) - xw2(:,1)
      jacpts(:,2) = xw1(:,2) + xw2(:,2)
      idx = argsort(N, jacpts(1:N,1))
      jacpts(:,1) = jacpts(idx, 1)
      jacpts(:,2) = jacpts(idx, 2)
      jacpts(:,2) = 1._wp/((1._wp-jacpts(:,1)**2)*jacpts(:,2)**2) !Quadrature weights
   end function jacpts

   pure function jacpts_main(N, a, b, flag)
      implicit none
      real(wp), parameter   :: eps = 1.e-15_wp
      integer, intent(in)   :: N
      real(wp), intent(in)  :: a, b
      logical, intent(in)   :: flag
      integer, allocatable  :: r(:)
      real(wp), allocatable :: C(:), T(:), x(:), dx(:), P(:), PP(:)
      real(wp) :: jacpts_main(N,2)
      integer :: size_r, l
      jacpts_main = 0._wp

      ! Asymptotic formula (WKB) - only positive x.
      if (flag) then
         size_r = ceiling(real(N,wp)/real(2,wp))
      else
         size_r = floor(real(n,wp)/real(2,wp))
      endif
      allocate(x(size_r), r(size_r), c(size_r), T(size_r), dx(size_r), P(size_r), PP(size_r))
      r = [size_r:1:-1]
      C = (real(2*r,wp)+a-0.5_wp)*pi/(real(2*n,wp)+a+b+1._wp);
      T = C + 1._wp/(real(2*n,wp)+a+b+1._wp)**2 * ((0.25_wp-a**2)*cotan(0.5_wp*C) - (0.25_wp-b**2)*tan(0.5_wp*C))
      x = cos(T)
   
      ! Initialise
      dx = 1.e+12_wp; 
      l = 0;
      ! Loop until convergence
      do while ( (maxval(abs(dx)) > eps) .and. (l < 10) )
         l = l + 1;
         P  = jacobi(N, a, b, x)
         PP = grad_jacobi(N, a, b, x)
         dx = -P/PP
         x = x + dx
      enddo
      jacpts_main(1:size_r,1) = x
      jacpts_main(1:size_r,2) = grad_jacobi(N, a, b, x)
   end function jacpts_main

   pure function LGL_nodes_weights(N)
      integer, intent(in) :: N
      real(wp) :: LGL_nodes_weights(N,2), buf(N-2,2), tmp
      buf = jacpts(N-2, 1._wp, 1._wp)
      LGL_nodes_weights(1,1) = -1._wp
      LGL_nodes_weights(2:N-1,1) = buf(:,1)
      LGL_nodes_weights(N,1) = 1._wp
      LGL_nodes_weights(:,2) = buf(:,2)
      tmp = 2._wp / real(n*(n-1),wp)
      LGL_nodes_weights([1,N],2) = tmp
      do j=2,N-1
         LGL_nodes_weights(j,2) = tmp / (Jacobi(N-1,0._wp, 0._wp, LGL_nodes_weights(j,1))**2)
      enddo
   end function LGL_nodes_weights

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

   pure function argsort(N, x, limit)
      implicit none
         !! input size of array, number of elements in the array to sort thEnd <= N+1
      real(wp), parameter :: tol = 1e-15
      integer , intent(in):: N
      integer , intent(in), optional :: limit
      real(wp), intent(in):: x(N) !! input array
      integer  :: argsort(N) !! output array of sorted argsortices
      real(wp) :: tmp, y(N)
      integer  :: i,j,tmpi, thEnd
      if (.not.present(limit)) then
         thEnd = N+1
      else
         thEnd = limit
      endif
      Y = X
      do i=1,N
         argsort(i) = i
      enddo
      do i=1,N-1
         do j=i+1,N
            if ((Y(j)-Y(i)).lt.(-tol)) then
               tmpi = argsort(j)  !
               argsort(j) = argsort(i)! swap indices
               argsort(i) = tmpi  !
               tmp  = Y(j)    !
               Y(j) = Y(i)    ! swap values
               Y(i) = tmp     !
            endif
         enddo
         if (i.eq.thEnd) exit
      enddo
   end function argsort
end module legendre_gauss_lobatto
