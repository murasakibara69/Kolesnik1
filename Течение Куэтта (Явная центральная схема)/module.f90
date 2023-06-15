module chill
  
  use types, only: wp
  
  implicit none
  
contains
  
  subroutine solver(u, delta_y, time, nu, A)
    real(wp), dimension(:), intent(inout) :: u
    real(wp), intent(in) :: delta_y, time, nu, A
    real(wp), dimension(:), allocatable :: u_next
    real(wp) :: t, delta_t, C1, C2
    integer :: i, n
    
    n = size(u, dim=1)
    
    delta_t = delta_y ** 2.0_wp / (6.0_wp * nu)
    
    allocate(u_next(n))
    
    t = delta_t
    C1 = nu * delta_t / delta_y ** 2.0_wp 
    C2 = A * delta_t 
    do while (t <= time)
      do i = 2, n - 1
        u_next(i) = u(i) + C1 * (u(i-1) - 2.0_wp * u(i) + u(i+1)) + C2
      end do
      u(2:n-1) = u_next(2:n-1)
      t = t + delta_t
    end do
    
    deallocate(u_next)
    
  end subroutine solver
  
  subroutine init_value(u, u_0)
    real(wp), dimension(:), intent(out) ::u
    real(wp), intent(in) :: u_0
    integer :: n
    
    n = size(u, dim=1)
    
    u(1:n-1) = 0.0_wp
    u(n) = u_0
    
  end subroutine init_value
  
  subroutine reading(n_y, H, time, nu, A, u_0)
    integer, intent(out) :: n_y
    real(wp), intent(out) :: H, time, nu, A, u_0
    integer, parameter :: io = 101
    integer :: ios
    character(len=256) :: str
    
    open(unit=io, file='input.dat', iostat=ios, iomsg=str, status="old", action="read")
    if (ios /= 0) stop trim(str)
    read(io, *) n_y
    read(io, *) h
    read(io, *) time
    read(io, *) nu
    read(io, *) A
    read(io, *) u_0
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
