! Erweitertes Programm zur Berechnung der Primzahlen

      program prime3
      
      implicit none
      integer i,j,k,l
      integer xx
      parameter(xx=100000)
      integer r(xx)

      do i=1,xx
         r(i)=i
      enddo

      do j=2,sqrt(real(xx))
         do k=j,xx,j
            if (k.gt.j) r(k)=0
         enddo
      enddo

      do l=1,xx
         if (r(l).ne.0) write(*,*) r(l)
      enddo
      
      end
