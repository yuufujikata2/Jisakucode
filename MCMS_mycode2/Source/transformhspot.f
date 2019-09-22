      implicit none      
      character*5      nsymbl
      character*1      line(75)
      integer          neq, nz, lcore, kmax, kplace, nc, ichg(10)
      integer          i, k, icard, nineqat, nat
      double precision x, y, z, rs, exfact, r(5), v(5), vcon
c
      nineqat = 0
      nat = 0
      open(10,file='.scratch',form='formatted',status='new')
 1    read(5,'(a5)') nsymbl
      backspace(5)
      if ( nsymbl(4:4).ne.'.') then
         nat = nat + 1
         read(5,'(a5,3i2,2i4,5f11.6,t76,i5)') 
     +        nsymbl, neq, nz, lcore, kmax, kplace, 
     +        x, y, z, rs, exfact, nc
         write(10,'(a5,2(1x,i2),2(1x,i4),4(1x,f11.6))')
     +        nsymbl, neq, nz, kmax, kplace, x, y, z, rs
c           
         if (neq .eq. 0) then
            nineqat = nineqat + 1
            read(5,'(10i5,t76,i5)') (ichg(i),i=1,10),nc        
            write(10,'(10(1x,i5))')  (ichg(i),i=1,10)
            read(5,'(t76,i5,t2,1p5e14.7)') nc, (r(i),i=1,5)
            write(10,'(5(1x,1pe14.7))') (r(i),i=1,5)
            do k = 1, kmax, 5
               icard = min0(kmax,k+4)-k+1
               read(5,'(t76,i5,t2,1p5e14.7)') nc, (v(i),i=1,icard)
               write(10,'(5(1x,1pe14.7))') (v(i),i=1,icard)
            enddo
         endif
         go to 1
      else
         read(5,'(t76,i5,t2,1p5e14.7)') nc, vcon
      endif
      close(10)
c
      open(10,file='.scratch',status='old')
      write(6,'(2(1x,i5),1x,1pe14.7)') nineqat, nat, vcon     
 2    read(10,'(80a1)',end=99) line 
      write(6,'(80a1)') line
      go to 2
 99   close(10, status='delete')
      write(6,'(a5)') 'XXXXX'
c
      End
