!-------------------------
!
! Calcurate cluster when call from python script 
!
!-------------------------
module fortcluster
  use ISO_C_binding
  implicit none
  type atom
    character(5) :: atomname
    integer iz
    double precision x,y,z
  end type atom
contains

  double precision function veclength(vec)
    double precision,intent(in) :: vec(3)
    integer j
    veclength=0.
    do j=1,3
      veclength = veclength + vec(j)**2
    end do
    veclength=sqrt(veclength)
  end function veclength
  
  subroutine baburusort(allatom,natoms)
    implicit none
    type(atom),intent(inout) :: allatom(:)
    integer,intent(in) :: natoms
    logical :: baburu=.true.
    integer :: baburu_c=0
    integer i
    double precision karix ,kariy ,kariz
    integer kariiz
    character(5) ::  kariname

    do while(baburu)
      baburu_c=baburu_c+1
      baburu=.false.
      do i=1,natoms-baburu_c
        if (allatom(i)%x**2+allatom(i)%y**2+allatom(i)%z**2 > &
           & allatom(i+1)%x**2+allatom(i+1)%y**2+allatom(i+1)%z**2) then
           kariname=allatom(i+1)%atomname
           kariiz=allatom(i+1)%iz
           karix=allatom(i+1)%x
           kariy=allatom(i+1)%y
           kariz=allatom(i+1)%z
           allatom(i+1)%atomname=allatom(i)%atomname
           allatom(i+1)%iz=allatom(i)%iz
           allatom(i+1)%x=allatom(i)%x
           allatom(i+1)%y=allatom(i)%y
           allatom(i+1)%z=allatom(i)%z
           allatom(i)%atomname=kariname
           allatom(i)%iz=kariiz
           allatom(i)%x=karix
           allatom(i)%y=kariy
           allatom(i)%z=kariz
           baburu=.true.
        end if
      end do
    end do
  end subroutine baburusort         
end module fortcluster

subroutine cluster(Plat,nbas,cart,alat,alliz,sa,sb,sc)
  use fortcluster
  implicit none
  double precision,parameter :: FACTOR=5
  double precision,parameter :: BOHR=0.529177
  double precision,intent(in) :: Plat(3,3)
  integer,intent(in) :: nbas
  double precision,intent(in) :: cart(nbas,3)
  double precision,intent(in) :: alat
  integer,intent(in) :: alliz(nbas)
  integer,intent(in) :: sa
  integer,intent(in) :: sb
  integer,intent(in) :: sc
  integer i0 , i1 , i2 , k ,j
  integer :: natoms=0
  integer :: n=1
  double precision y_xyz(3)
  double precision y_ctrl(3)
  integer nbas2
  integer lmax
  type(atom),allocatable :: allclass(:)
  type(atom),allocatable :: allatom_xyz(:)
  type(atom),allocatable :: allatom_ctrl(:)
  integer :: count24 = 0, count8 = 0, count0 = 0
  
  
  open(17,file="cluster.in",status="old")
  do j=1,5
   read(17,*)
  end do

  read(17,*) nbas2

  if (nbas/=nbas2)then
    print *,"error1"
    stop 
  end if

 allocate(allclass(nbas))
 
  do j=1,nbas
    read (17,*) allclass(j)%atomname ,&
               &  allclass(j)%x , allclass(j)%y , allclass(j)%z
  end do
  
  close(17)

  print*,""
  print*,"allsite"
  do j=1,nbas
    write (*,40) allclass(j)%atomname,allclass(j)%x,allclass(j)%y,allclass(j)%z
  end do
  40 format(A4," ",F10.6," ",F10.6," "F10.6) 
  
  natoms=sa * sb * sc * nbas

  write (*,51) "natoms=",natoms 
  51 format(A7,I4)

  allocate(allatom_ctrl(natoms))
  allocate(allatom_xyz(natoms))

  do i0=0,sa-1
    do i1=0,sb-1
      do i2=0,sc-1
        do k=1,nbas
          do j=1,3
            y_ctrl(j)=i0*Plat(1,j)+i1*Plat(2,j)+i2*plat(3,j) &
                  & +cart(k,j) 
            y_xyz(j) = y_ctrl(j)* alat * BOHR 
          end do
          allatom_ctrl(n)%atomname=allclass(k)%atomname
          allatom_ctrl(n)%iz=alliz(k)
          allatom_ctrl(n)%x=y_ctrl(1)
          allatom_ctrl(n)%y=y_ctrl(2)
          allatom_ctrl(n)%z=y_ctrl(3)
          allatom_xyz(n)%atomname=allclass(k)%atomname
          allatom_xyz(n)%iz=alliz(k)
          allatom_xyz(n)%x=y_xyz(1)
          allatom_xyz(n)%y=y_xyz(2)
          allatom_xyz(n)%z=y_xyz(3)
          n=n+1
        end do
      end do
    end do
  end do
  
  
  call baburusort(allatom_xyz,natoms)
  call baburusort(allatom_ctrl,natoms)
  
         
  open(18,file="susel.xyz",status="replace")
  
  write (18,100) natoms
  write (18,*)
  do j=1,natoms
    write (18,101) allatom_xyz(j)%atomname,allatom_xyz(j)%x,allatom_xyz(j)%y,allatom_xyz(j)%z 
  end do
  100 format(I4)
  101 format(A4," ",F10.6," ",F10.6," "F10.6) 
  close(18)

  open(19,file="susel.ctrl",status="replace")
  write(19,102) natoms 
  do j=1,natoms
    if (allatom_ctrl(j)%iz == 24) then
      write (19,103) allatom_ctrl(j)%atomname,"Cr",count24,allatom_ctrl(j)%x,allatom_ctrl(j)%y,allatom_ctrl(j)%z 
      count24 = count24 + 1
    else if (allatom_ctrl(j)%iz == 8) then
      write (19,103) allatom_ctrl(j)%atomname,"O",count8,allatom_ctrl(j)%x,allatom_ctrl(j)%y,allatom_ctrl(j)%z 
      count8 = count8 + 1
    else if (allatom_ctrl(j)%iz == 0) then
      write (19,103) allatom_ctrl(j)%atomname,"E",count0,allatom_ctrl(j)%x,allatom_ctrl(j)%y,allatom_ctrl(j)%z 
      count0 = count0 + 1
    endif
  end do  
  102 format(I5,I5,I5)
  103 format(A4,"a",A4,I0,"  ",F10.6," ",F10.6," ",F10.6)

  deallocate(allclass)
  deallocate(allatom_xyz)
  deallocate(allatom_ctrl)
end subroutine cluster

