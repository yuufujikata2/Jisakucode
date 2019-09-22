      program instr2sym
      implicit none
      integer nat, idum, iz, i, l, m, lmax, nstates 
      double precision x, y, z 
      character*10 type
      open(4, FILE ='sym.dat', status='unknown')
      open(3, FILE ='instr.dat') 
      read(3,*) nat, idum, idum
      nstates = 0
      do i = 1, nat 
        read(3,*) type, lmax, iz, x,y,z 
        nstates = nstates + (lmax+1)*(lmax+1) 
      enddo
      write(4,'(2(I5))') nstates, 2
      close(3)
      open(3, FILE ='instr.dat')
      read(3,*) nat, idum, idum
      idum = 0 
      do i = 1, nat
        read(3,*) type, lmax, iz, x,y,z
        do l=0,lmax 
          write(4,'(3(I5))') l, 1, 0 
          write(4,'(3(I5),F15.10)') 0,1,i,1. 
          idum = idum+1
          do m=1,l
            write(4,'(3(I5))') l, 1, 0
            write(4,'(3(I5),F15.10)') m,1,i,1.
            write(4,'(3(I5))') l, 1, 0
            write(4,'(3(I5),F15.10)') m,-1,i,1.
            idum = idum+2
          enddo
        enddo
      enddo
      close(3) 
      close(4) 
      if ( idum .ne. nstates ) then
        print*,"PANIC: instr -> sym: idum=",idum,"!=",nstates,"=nstates"
        stop
      endif

      END
