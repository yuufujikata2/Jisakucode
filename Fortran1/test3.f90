program main
  implicit none
  integer,dimension(15,8) :: haichi=0
  integer i,j,k
  double precision t1,t2
  integer m
  integer :: c1=0,c2=0,c3=0,c4=0

  haichi(1,:)=20
  haichi(:,1)=20
  haichi(:,8)=20
  haichi(14,2:7)=(/0,0,0,0,0,0/)
  haichi(13,2:7)=(/0,0,0,0,0,0/)
  haichi(12,2:7)=(/0,0,0,0,0,0/)
  haichi(11,2:7)=(/0,0,0,0,0,0/)
  haichi(10,2:7)=(/0,0,0,0,0,0/)
  haichi(9,2:7)= (/0,0,0,0,0,0/)
  haichi(8,2:7)= (/0,0,0,0,0,0/)
  haichi(7,2:7)= (/0,0,0,0,0,0/)
  haichi(6,2:7)= (/0,0,0,0,0,0/)
  haichi(5,2:7)= (/2,0,0,0,0,4/)
  haichi(4,2:7)= (/2,1,3,2,4,4/)
  haichi(3,2:7)= (/2,2,1,3,2,2/)
  haichi(2,2:7)= (/1,1,3,3,2,4/)
  call hyouji(haichi)
  m=rensa(haichi)
  call hyouji(haichi)
contains
  integer function rakka(haichi)
    implicit none
    integer,intent(inout) :: haichi(:,:)
    integer i,j
    integer :: rakka_c
    rakka_c=1
    rakka=0
    do while(rakka_c/=0)
      rakka_c=0
      do i=2,7 !ubound(haichi,2)
        do j=2,14 !ubound(haichi,1)
          if (haichi(j,i)/=0 .and. haichi(j-1,i)==0) then
            haichi(j-1,i)=haichi(j,i)
            haichi(j,i)=0
            rakka_c=rakka_c+1
          end if
        end do
      end do
      rakka=rakka+1
    end do
  end function rakka

  subroutine hyouji(haichi)
    implicit none
    integer,intent(in) :: haichi(:,:)
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
  
  recursive function renketu(haichi,i,j,count1,count2,count3,count4) result(m)
    implicit none
    integer,intent(in) :: i,j
    integer,intent(inout)  :: haichi(:,:)
    integer a
    integer :: m
    integer,intent(inout) :: count1,count2,count3,count4

    a=haichi(j,i)
    if (a==0 .or. a==5 .or. a==11 .or. a==12 .or. a==13 .or. a==14 .or.a==15 .or. a==20) then 
      m=1
    else
      if (haichi(j,i)==1) then
        count1=count1+1
        haichi(j,i)=11
      else if (haichi(j,i)==2) then 
        count2=count2+1
        haichi(j,i)=12
      else if (haichi(j,i)==3) then
        count3=count3+1
        haichi(j,i)=13
      else if (haichi(j,i)==4) then
        count4=count4+1
        haichi(j,i)=14
      end if
      if (a==haichi(j,i+1)) then
        m=renketu(haichi,i+1,j,count1,count2,count3,count4)
      else if (haichi(j,i+1)==5) then
        haichi(j,i+1)=15
      end if
      if (a==haichi(j,i-1)) then
        m=renketu(haichi,i-1,j,count1,count2,count3,count4)
      else if (haichi(j,i-1)==5) then
        haichi(j,i-1)=15
      end if
      if (j/=12 .and. j/=13 .and. j/=14)then
        if (a==haichi(j+1,i)) then
          m=renketu(haichi,i,j+1,count1,count2,count3,count4)
        else if (haichi(j,i+1)==5) then
          haichi(j,i+1)=15
        end if
      end if
      if (a==haichi(j-1,i)) then
        m=renketu(haichi,i,j-1,count1,count2,count3,count4)
      else if (haichi(j-1,i)==5) then
        haichi(j-1,i)=15
      end if
    end if
    m=1

  end function renketu

  recursive function renketu2(haichi,i,j) result(m)
    implicit none
    integer,intent(in) :: i,j
    integer,intent(inout)  :: haichi(:,:)
    integer a
    integer :: m

    a=haichi(j,i)
    if (a==0 .or. a==5 .or. a==1 .or. a==2 .or. a==3 .or. a==4 .or.a==15 .or. a==20) then 
      m=1
    else
      if (haichi(j,i)==11) then
        haichi(j,i)=0
      else if (haichi(j,i)==12) then 
        haichi(j,i)=0
      else if (haichi(j,i)==13) then
        haichi(j,i)=0
      else if (haichi(j,i)==14) then
        haichi(j,i)=0
      end if
      if (a==haichi(j,i+1)) then
        m=renketu2(haichi,i+1,j)
      else if (haichi(j,i+1)==15) then
        haichi(j,i+1)=0
      end if
      if (a==haichi(j,i-1)) then
        m=renketu2(haichi,i-1,j)
      else if (haichi(j,i-1)==15) then
        haichi(j,i-1)=0
      end if
      if (j/=12 .and. j/=13 .and. j/=14)then
        if (a==haichi(j+1,i)) then
          m=renketu2(haichi,i,j+1)
        else if (haichi(j,i+1)==15) then
          haichi(j,i+1)=0
        end if
      end if
      if (a==haichi(j-1,i)) then
        m=renketu2(haichi,i,j-1)
      else if (haichi(j-1,i)==15) then
        haichi(j-1,i)=0
      end if
    end if
    m=1

  end function renketu2


  integer function renketukeshi(haichi)
    implicit none
    integer,intent(inout) :: haichi(:,:)
    integer count1,count2,count3,count4
    integer i,j
    integer m

    renketukeshi=0
    do i=2,7
      do j=2,13
        if (haichi(j,i)==0 .or. haichi(j,i)==11 .or. haichi(j,i)==12 .or. haichi(j,i)==13 .or. haichi(j,i)==14)cycle
        count1=0
        count2=0
        count3=0
        count4=0
        m=renketu(haichi,i,j,count1,count2,count3,count4)
        if (count1>=4 .or. count2>=4 .or. count3>=4 .or. count4>=4)then
          renketukeshi=renketukeshi+count1+count2+count3+count4
          m=renketu2(haichi,i,j)
        end if
      end do
    end do
    do i=2,7
      do j=2,13
        if (haichi(j,i)==0)cycle
        if (haichi(j,i)==11)then
          haichi(j,i)=1
        else if (haichi(j,i)==12)then
          haichi(j,i)=2
        else if (haichi(j,i)==13)then
          haichi(j,i)=3
        else if (haichi(j,i)==14)then
          haichi(j,i)=4
        else if (haichi(j,i)==15)then
          haichi(j,i)=5
        end if
      end do
    enddo

  end function renketukeshi  

  recursive function rensa(haichi) result(m)
    implicit none
    integer,intent(inout) :: haichi(:,:)
    integer m,j,k
    
    j=renketukeshi(haichi) 
    k=rakka(haichi)
    if (k==1)then
      m=1
    else
      call hyouji(haichi)
      m=rensa(haichi)
    end if
  end function rensa
end program main
