program main
  
  use types, only: wp
  use chill, only: reading, init_value, solver, writing
  
  implicit none
  
  integer :: i, n
  real(wp) :: H, time, nu, A, u_0, delta_y
  real(wp), dimension(:), allocatable :: u, y
  
  call reading(n, H, time, nu, A, u_0)
  
  delta_y = H / (n - 1)
  
  allocate(u(n))
  
  call init_value(u, u_0)
  
  call solver(u, delta_y, time, nu, A)
  
  allocate(y(n))
  
  y = [(delta_y * (i - 1), i = 1, n)]
  
  call writing(y, u)
  
  deallocate(u, y)
  
  call system("gnuplot -p plot.plt")
  
end program main
