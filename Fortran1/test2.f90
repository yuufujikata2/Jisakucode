program main
  implicit none
  real,dimension(3,5) :: a=1
  integer i,j
  real f
  a(2:3,1:2)=5
  !do i=1,3
  !  do j=1,5
  !    print*,a(i,j)
  !  end do
  !end do
 ! print*,a
   f=myfunc(a)
contains
  real function myfunc(x)
    real,intent(in) :: x(:,:)
    integer i
    myfunc=0
    do i=1,ubound(x,1)
      do j=1,ubound(x,2)
        print*,x(i,j)
      end do
    end do
  end function myfunc
end program main
