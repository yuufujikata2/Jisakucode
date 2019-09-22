program main
  integer i,j,k
  integer x
  integer a,b
  i=1
  j=2
  k=3
  a=1
  b=100
  x=func(i,j,k)
  print*,real(1/100)
  print*,x,i,j,k
contains  
  integer function func(a,b,c)
    integer,intent(in) :: a
    integer,intent(inout) :: b,c
    b=a*10
    c=b+10
    func=b+c
  end function 

end program main
