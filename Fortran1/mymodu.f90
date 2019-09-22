module my_modu
  implicit none
  !real,dimension(3,4):: a
  !real,dimension(4,2):: b
  !real,allocatable,dimension(:,:)::x
  !a=reshape((/1,2,3,4,5,6,7,8,9,10,11,12/),(/3,4/))
  !b=reshape((/1,2,3,4,5,6,7,8/),(/4,2/))
  !x=myproduct(a,b)
  !print*,x
contains
  function myproduct(a,b)
    implicit none
    real,intent(in):: a(:,:),b(:,:)
    real,dimension(ubound(a,1),ubound(b,2)) :: myproduct
    integer i,j,k
    if(ubound(a,2)/=ubound(b,1))then
      print*,'error'
      stop
    end if
    do i=1,ubound(a,1)
      do j=1,ubound(b,2)
        myproduct(i,j)=0
        do k=1,ubound(a,2)
          myproduct(i,j)=myproduct(i,j)+a(i,k)*b(k,j) 
        end do
      end do
    end do
  end function myproduct
  
  function myinteg(f,a,b)
    implicit none
    interface
      real function f(x)
        real,intent(in) :: x
      end function f
    end interface
    real,intent(in)::a,b
    real myinteg
    integer i, n
    myinteg=0
    n=10000
    do i=0,n
      myinteg=myinteg+f(a+(b-a)*real(i)/real(n))*real((b-a)/real(n))
    end do
  end function myinteg
end module my_modu
