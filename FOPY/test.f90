subroutine test(n,A)
  implicit none
  integer,intent(in) :: n
  integer,intent(in) :: A(15,8)

  print*,A
end subroutine test
subroutine hyouji(haichi)
  implicit none
  integer,intent(in) :: haichi(15,8)
  integer i,j
  do i=15,1,-1
    do j=1,8
      if (j<=7) then
        write(*,fmt='(I3)',advance='no') haichi(i,j)
      else
        write(*,fmt='(I3)') haichi(i,j)
      end if 
    end do
  end do
end subroutine hyouji

