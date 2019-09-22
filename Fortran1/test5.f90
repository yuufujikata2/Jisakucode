program main
  use my_modu
  implicit none
  real S
  real a,b

  a=5 
  b=1
  S=myinteg(myfunc,a,b)
  print*,S
contains 
  real function myfunc(x)
    real,intent(in) :: x
    myfunc=sqrt(x)
  end function myfunc
end program main
