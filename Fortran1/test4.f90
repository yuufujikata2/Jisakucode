program main
  use my_modu
  implicit none
  real,dimension(3,4):: a
  real,dimension(4,2):: b
  real,allocatable,dimension(:,:)::x
  a=reshape((/1,2,3,4,5,6,7,8,9,10,11,12/),(/3,4/))
  b=reshape((/1,2,3,4,5,6,7,8/),(/4,2/))
  x=myproduct(a,b)
  print*,x

end program main
