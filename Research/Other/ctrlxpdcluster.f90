!-------------------------
!
! Calcurate cluster(surface)  when call from python script 
!
!-------------------------
module fortcluster
  use ISO_C_binding
  implicit none
  type atom
    character(5) :: atomname
    integer iz
    real x,y,z
  end type atom
contains

  real function veclength(vec)
    real,intent(in) :: vec(3)
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
    real karix ,kariy ,kariz
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

  subroutine rotate(cart,Plat,nbas,miller,newcart,newPlat,neworigincart)

    implicit none
    real,parameter :: PI=3.14159

    real,intent(inout) :: cart(nbas,3)
    real,intent(inout) :: Plat(3,3)
    integer,intent(in) :: nbas
    real,intent(in) :: miller(3)
    real,intent(out) :: newcart(nbas,3)
    real,intent(out) :: newplat(3,3)
    real,intent(out) :: neworigincart(3)
    
    real cart2(nbas,3)
    real Plat2(3,3),tPlat(3,3)
    real origincart(3),origincart2(3)

    real millerlength
    real millervec(3)
    real alpha,sita
    real Ra(3,3),Rs(3,3)

    integer i,j,k

 
    do i=1,3
      millervec(i) = (miller(1) * Plat(1,i) + miller(2) * Plat(2,i) + miller(3) * Plat(3,i))
    end do

    millerlength = veclength(millervec)


    if (millerlength ==0.) then
      newPlat = Plat
      newcart = cart
      neworigincart = 0.
      return
    end if 

    millervec = millervec / millerlength

    neworigincart = (/0.,0.,millerlength/)

    if (millervec(1)== 0.) then
      alpha = PI / 2
    else
      alpha = atan(millervec(2) / millervec(1))
    end if

    if (millervec(3) == 0.) then
      sita = PI / 2
    else
      sita = atan(sqrt(millervec(1)**2 + millervec(2)**2) / millervec(3))
    end if

    Ra(1,:) = (/cos(alpha),-sin(alpha),0./)
    Ra(2,:) = (/sin(alpha),cos(alpha),0./)
    Ra(3,:) = (/0., 0., 1./)

    Rs(1,:) = (/cos(-sita), 0., sin(-sita)/)
    Rs(2,:) = (/0., 1., 0./)
    Rs(3,:) = (/-sin(-sita), 0., cos(-sita)/)


    do i =1,3
      do j =1,3
        tPlat(i,j)= Plat(i,j)
      end do
    end do

    do i = 1,3
      do j = 1,3
        Plat(j,i) = tPlat(i,j)
      end do
    end do

    do i = 1,nbas
      do j =1,3
        cart2(i,j) = Ra(j,1) * cart(i,1) + Ra(j,2) * cart(i,2) + Ra(j,3) * cart(i,3)
      end do
    end do

    do i = 1,nbas
      do j =1,3
        newcart(i,j) = Rs(j,1) * cart2(i,1) + Rs(j,2) * cart2(i,2) + Rs(j,3) * cart2(i,3)
      end do
    end do

    do i = 1,3
      do j =1,3
        Plat2(i,j) = Ra(i,1) * Plat(1,j) + Ra(i,2) * Plat(2,j) + Ra(i,3) * Plat(3,j)
      end do
    end do

    do i = 1,3
      do j =1,3
        newPlat(i,j) = Rs(i,1) * Plat2(1,j) + Rs(i,2) * Plat2(2,j) + Rs(i,3) * Plat2(3,j)
      end do
    end do

    do i = 1,3
      do j = 1,3
        tPlat(i,j) = newPlat(i,j)
      end do
    end do

    do i = 1,3
      do j = 1,3
        newPlat(i,j) = tPlat(j,i)
      end do
    end do
    
  end subroutine rotate
end module fortcluster

subroutine cluster(clradius,Plat,nbas,cart,alat,alliz,miller)
  use fortcluster
  implicit none
  real,parameter :: FACTOR=5
  real,parameter :: BOHR=0.529177
  real,intent(in) :: clradius
  real,intent(inout) :: Plat(3,3)
  integer,intent(in) :: nbas
  real,intent(inout) :: cart(nbas,3)
  real,intent(in) :: alat
  integer,intent(in) :: alliz(nbas)
  real,intent(inout) :: miller(3)
  real cutoff
  integer i0max , i1max , i2max
  integer i0 , i1 , i2 , k ,j
  integer :: natoms=0
  integer :: n=1
  real y(3)
  integer nbas2
  integer lmax
  type(atom),allocatable :: allclass(:)
  type(atom),allocatable :: allatom(:)
  real newcart(nbas,3)
  real neworigincart(3)
  real newPlat(3,3)
  
  
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

  call rotate(cart,Plat,nbas,miller,newcart,newPlat,neworigincart)


  cutoff = clradius
  i0max = cutoff * FACTOR / veclength(Plat(1,:))
  i1max = cutoff * FACTOR / veclength(Plat(2,:))
  i2max = cutoff * FACTOR / veclength(Plat(3,:))
  
  write (*,50) "i0max=",i0max,"i1max=",i1max,"i2max=",i2max
  50 format(A6,I3," ",A6,I3," ",A6,I3)
  do i0=-i0max,i0max
    do i1=-i1max,i1max
      do i2=-i2max,i2max
        do k=1,nbas
          do j=1,3
            y(j)=i0*newPlat(1,j)+i1*newPlat(2,j)+i2*newplat(3,j) &
                 &  +newcart(k,j) - neworigincart(j)
            y(j) = y(j)* alat * BOHR 
          end do
          if (veclength(y)< clradius .and. y(3) < 0.001) then
            natoms=natoms+1
          end if
        end do
      end do
    end do
  end do

  write (*,51) "natoms=",natoms 
  51 format(A7,I4)

  allocate(allatom(natoms))

  do i0=-i0max,i0max
    do i1=-i1max,i1max
      do i2=-i2max,i2max
        do k=1,nbas
          do j=1,3
            y(j)=i0*newPlat(1,j)+i1*newPlat(2,j)+i2*newplat(3,j) &
                  & +newcart(k,j) - neworigincart(j)
            y(j) = y(j)* alat * BOHR 
          end do
          if (veclength(y)< clradius .and. y(3) < 0.001) then
            allatom(n)%atomname=allclass(k)%atomname
            allatom(n)%iz=alliz(k)
            allatom(n)%x=y(1)
            allatom(n)%y=y(2)
            allatom(n)%z=y(3)
            n=n+1
          end if
        end do
      end do
    end do
  end do
  
  print* ,"coordinate" 
  do j=1,natoms
    write (*,101) allatom(j)%atomname,allatom(j)%x,allatom(j)%y,allatom(j)%z 
  end do
  
  call baburusort(allatom,natoms)

  print*,"sort coordinate"
  do j=1,natoms
    write (*,101) allatom(j)%atomname,allatom(j)%x,allatom(j)%y,allatom(j)%z 
  end do
         
  open(18,file="cluster.xyz",status="replace")
  
  write (18,100) natoms
  write (18,*)
  do j=1,natoms
    write (18,101) allatom(j)%atomname,allatom(j)%x,allatom(j)%y,allatom(j)%z 
  end do
  100 format(I4)
  101 format(A4," ",F10.6," ",F10.6," "F10.6) 
  close(18)

  open(19,file="instr",status="replace")
  write(19,102) natoms , natoms , 2
  do j=1,natoms
    if (allatom(j)%iz==0) then
      lmax=0 
    else if (allatom(j)%iz<=18) then 
      lmax=1
    else 
      lmax=2
    end if 
    write(19,103) allatom(j)%atomname,lmax,allatom(j)%iz, &
                  & allatom(j)%x/BOHR,allatom(j)%y/BOHR,allatom(j)%z/BOHR
  
  end do  
  102 format(I5,I5,I5)
  103 format(A4," ",I1," ",I2,"  ",F10.6," ",F10.6," ",F10.6)

  deallocate(allclass)
  deallocate(allatom)
end subroutine cluster

