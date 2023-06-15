module chill
  
  use types, only: wp
  
  implicit none
  
contains
  
  subroutine solver(u, time, delta_y, nu, A, omega)
    real(wp), dimension(:), intent(inout) :: u
    real(wp), intent(in) :: time, delta_y, nu, A, omega
    real(wp), dimension(:), allocatable :: u_next, aa, bb, dd
    real(wp) :: t, delta_t, Ca, Cb, C1, C2, C3, C4
    integer :: i, n
    
    n = size(u, dim=1)
    
    delta_t = delta_y ** 2.0_wp / (6.0_wp * nu)
    
    allocate(aa(n-2), bb(n-2), dd(n-2))
    
    Ca = nu / delta_y ** 2.0_wp
    Cb = - 2.0_wp / delta_t - 2.0_wp * nu / delta_y ** 2.0_wp
    
    aa = Ca
    bb = Cb
    
    allocate(u_next(n))
    
    C1 = - 2.0_wp / delta_t
    C3 = nu * delta_t / (2.0_wp * delta_y ** 2.0_wp)
    t = delta_t
    do while (t <= time)
      C2 = A * cos(omega * t)
      C4 = -C2 * delta_t / 2.0_wp
      dd(1) = C1 * u(2) - aa(1) * u(1) + C2
      do i = 2, n - 3
        dd(i) = C1 * u(i+1) + C2
      end do
      dd(n-2) = C1 * u(n-1) - aa(n-2) * u(n) + C2
      call sweep_method(aa, bb, aa, dd, u(2:n-1))
      do i = 2, n - 1
        u_next(i) = u(i) + C3 * (u(i-1) - 2.0_wp * u(i) + u(i+1)) + C4
      end do
      u(2:n-1) = u_next(2:n-1)
      u(n) = u(n-1)
      t = t + delta_t
    end do
    
    deallocate(aa, bb, dd, u_next)
    
  end subroutine solver
  
  subroutine sweep_method(a, b, c, d, x)
    real(wp), dimension(:), intent(in) :: a, b, c, d
    real(wp), dimension(:), intent(out) :: x
    real(wp), dimension(:), allocatable :: alpha, beta, gama
    integer :: i, n
    
    n = size(x, dim=1)
    
    allocate(alpha(n-1), beta(n), gama(n))
    
    gama(1) = b(1)
    alpha(1) = - c(1) / gama(1)
    beta(1) = d(1) / gama(1)
    
    do i = 2, n - 1
      gama(i) = b(i) + a(i) * alpha(i-1)
      alpha(i) = - c(i) / gama(i)
      beta(i) = (d(i) - a(i) * beta(i-1)) / gama(i)
    end do
    
    gama(n) = b(n) + a(n) * alpha(n-1)
    beta(n) = (d(n) - a(n) * beta(n-1)) / gama(n)
    
    x(n) = beta(n)
    
    do i = n - 1, 1, -1
      x(i) = beta(i) + alpha(i) * x(i+1)
    end do
    
    deallocate(alpha, beta, gama)
    
  end subroutine sweep_method
  
  subroutine reading(n_y, nu, A, omega)
    integer, intent(out) :: n_y
    real(wp), intent(out) :: nu, A, omega
    integer, parameter :: io = 101
    integer :: ios
    character(len=256) :: str
    
    open(unit=io, file='input.dat', iostat=ios, iomsg=str, status="old", action="read")
    if (ios /= 0) stop trim(str)
    read(io, *) n_y
    read(io, *) nu
    read(io, *) A
    read(io, *) omega
    close(io)
    
  end subroutine reading
  
  subroutine writing(y, u)
    real(wp), dimension(:), intent(in) :: y, u
    integer, parameter :: io = 102
    integer :: ios, i, n
    character(len=256) :: str
    
    n = size(y, dim=1)
    open(unit=io, file='result.dat', iostat=ios, iomsg=str, status="replace", action="write")
    if (ios /= 0) stop trim(str)
    do i = 1, n
      write(io, *) u(i), y(i)
    end do
    close(io)
    
  end subroutine writing
  
end module chill
