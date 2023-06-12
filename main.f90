program main
   use legendre_gauss_lobatto
   implicit none
   integer, parameter :: N = 101
   integer, parameter :: degree = 52
   double precision :: x(N), y(N), cites(degree),weights(degree),fa(degree), num_integral,integral
   double precision :: cites_physical(degree), x_a, x_b, pade(degree+1), buffer(degree,2)
   integer :: i
   !! test functions
   x = linspace(N=N, a=0.d0, b=1.d0)
   y = jacobi(N=degree, alpha=1.d0, beta=1.d0, x=x)
   y = jacobi_normalized(N=degree, alpha=1.d0, beta=1.d0, x=x)
   y = grad_jacobi(N=degree, alpha=1.d0, beta=1.d0, x=x)
   y = grad_jacobi_normalized(N=degree, alpha=1.d0, beta=1.d0, x=x)

   buffer = GQ_nodes_weights(degree)
   cites = buffer(:,1)
   weights = buffer(:,2)
   
   fa = -1._wp
   do i=1,degree
      if (cites(i) > 0._wp) fa(i) = 1._wp
   enddo
   do i=1,degree
      write(22,*) cites(i), fa(i)
   enddo
   pade = pade_reconstruction(degree, fa)
   do i=1,degree
      write(33,*) cites(i), pade(i)
   enddo
   !! test numerical integration 
   x_a = -1.d0; x_b = 2.d0
   cites_physical = x_a + 0.5d0*(x_b-x_a) * (cites+1.d0)
   num_integral = lgl_quadrature(N=degree, f=fa, a=-1.d0, b=2.d0)
   integral = exp(2.d0)-exp(-1.d0)
   !write(*,'(a,E15.5)') 'Numerical Integration error : ', abs(num_integral - integral)
end program main
