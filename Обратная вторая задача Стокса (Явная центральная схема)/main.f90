program main
  
  use types, only: wp
  use chill, only: reading, solver, writing
  
  implicit none
  
  integer :: i, n
  real(wp) :: time, nu, A, omega, delta, delta_y, nperiods
  real(wp), dimension(:), allocatable :: u, y
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  
  call reading(n, nu, A, omega)
  
  delta = sqrt(2.0_wp * nu / omega)
  
  delta_y = delta / (n - 1)
  
  read*, nperiods
  time = nperiods * 2.0_wp * pi / omega
  
  allocate(u(n), source=0.0_wp)
  
  call solver(u, time, delta_y, nu, A, omega)
  
  allocate(y(n))
  
  y = [(delta_y * (i - 1), i = 1, n)]
  
  call writing(y, u)
  
  deallocate(u, y)
  
  call system("gnuplot -p plot.plt")
  
end program main
