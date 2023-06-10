program main
   use legendre_gauss_lobatto
   implicit none
   integer, parameter :: N = 101
   integer, parameter :: degree = 6
   double precision :: x(N), y(101), cites(degree),weights(degree),fa(degree), num_integral,integral
   double precision :: cites_physical(degree), x_a, x_b, buffer(degree,2)

   !! test functions
   x = linspace(N=N, a=0.d0, b=1.d0)
   y = jacobi(N=degree, alpha=1.d0, beta=1.d0, x=x)
   y = jacobi_normalized(N=degree, alpha=1.d0, beta=1.d0, x=x)
   y = grad_jacobi(N=degree, alpha=1.d0, beta=1.d0, x=x)
   y = grad_jacobi_normalized(N=degree, alpha=1.d0, beta=1.d0, x=x)

   buffer = LGL_nodes_weights(degree)
   cites = buffer(:,1)
   weights = buffer(:,2)
   
   !! test numerical integration 
   x_a = -1.d0; x_b = 2.d0
   cites_physical = x_a + 0.5d0*(x_b-x_a) * (cites+1.d0)
   fa = exp(cites_physical)
   num_integral = lgl_quadrature(N=degree, f=fa, a=-1.d0, b=2.d0)
   integral = exp(2.d0)-exp(-1.d0)
   write(*,'(a,E15.5)') 'Numerical Integration error : ', abs(num_integral - integral)
end program main
