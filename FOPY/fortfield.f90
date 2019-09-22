module fortfield
  implicit none
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

  recursive function rensa(haichi,rensa_c) result(m)
    implicit none
    integer,intent(inout) :: haichi(:,:)
    integer m,j,k
    integer,intent(inout) ::  rensa_c
    j=renketukeshi(haichi) 
    if(j/=0) rensa_c=rensa_c+1
    k=rakka(haichi)
    if (k==1)then
      m=1
    else
      call hyouji(haichi)
      m=rensa(haichi,rensa_c)
    end if
  end function rensa

  integer function renketushirabe(haichi,i,j)
    implicit none
    integer,intent(inout) :: haichi(:,:)
    integer,intent(in) :: i,j
    integer :: count1=0,count2=0,count3=0,count4=0
    integer,dimension(15,8) :: karihaichi
    integer m
    integer k,l

    karihaichi=haichi
    m=renketu(haichi,i,j,count1,count2,count3,count4)
    if (count1>=4 .or. count2>=4 .or. count3>=4 .or. count4>=4) then
      m=renketu2(haichi,i,j)
      renketushirabe=1
    else
      do k=2,7
        do l=2,13
          if (haichi(l,k)==0) exit
          if (haichi(l,k)==11) then
            haichi(l,k)=1
          else if (haichi(l,k)==12)then
            haichi(l,k)=2
          else if (haichi(l,k)==13)then
            haichi(l,k)=3
          else if (haichi(l,k)==14)then
            haichi(l,k)=4
          end if
        end do
      end do 
      renketushirabe=0
    end if
  end function renketushirabe

end module fortfield


integer function  rensashirabe(haichi)
  use fortfield
  implicit none
  integer,dimension(6,12) :: kekka_rensasuu=0
  integer :: rensa_c=0
  integer,intent(inout) :: haichi(15,8)
  integer  karihaichi(15,8)
  integer i,j,k
  integer m
  integer :: rensamax=0

  karihaichi=haichi
  do i=2,7
    do j=13,2,-1
      if (haichi(j,i)/=0) exit
      if (haichi(j-2,i)==0 .or. haichi(j-2,i)==20) cycle
      if (haichi(j-1,i)==0 .and. (haichi(j,i-1)==0 .or. haichi(j,i-1)==20) .and. (haichi(j,i+1)==0 .or. haichi(j,i+1)==20)) cycle
      rensa_c=0
      do k=1,4
        haichi(j,i)=k
        if (renketushirabe(haichi,i,j)==1)then
          m=rakka(haichi)
          m=rensa(haichi,rensa_c)
          kekka_rensasuu(i-1,j-1)=rensa_c
          haichi=karihaichi
        else
          haichi(j,i)=0
        end if
      end do
    end do
  end do
  do j=1,12
    do i=1,6
      if (kekka_rensasuu(i,j)>rensamax) rensamax=kekka_rensasuu(i,j)
    end do
  end do
  rensashirabe=rensamax

end function rensashirabe







