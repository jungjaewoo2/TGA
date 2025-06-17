c     ============================================
      subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
c     ============================================
c
c     # set auxiliary arrays 
c     # dummy routine when no auxiliary arrays
c
c     
      implicit double precision (a-h,o-z)
c      dimension aux(1-mbc:maxmx+mbc, *)
      dimension aux(1-mbc:maxmx+mbc, maux)
c
       return
       end
