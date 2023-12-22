       subroutine bh_out(xs,vs,bodys)
*
*      Additional output of bh neighbour particles 
*
       include "common6.h"
       integer nb, cnt 
     &, ind_sort(nmax),var_sort(nmax)
*
*    use single precision
      real*4 tmp_r(nmax),tmp_v,tmp_a,tmp_adot,fs(3,nmax),
     &       fdots(3,nmax)
      REAL*4  XS(3,NMAX),VS(3,NMAX),BODYS(NMAX),RHOS(NMAX),AS(20)
      CHARACTER*30 OUTFILE
      CHARACTER*20 TCHAR
*
      do i=ifirst,ntot 
*    create array of distances:
           tmp_r(i-ifirst+1) = sqrt(xs(1,i)**2+xs(2,i)**2+xs(3,i)**2)
*    Convert force & derivative to single precision:
        do j=1,3
           fs(j,i) = real( f(j,i) )
           fdots(j,i) = real( fdot(j,i) )        
        enddo
*
      enddo
*
*    Make sure time is updated
      ttot = time+toff
*    Number of BH neighbour particles
      nb = 100
*
*    Sort distances and create array of indexes
*    using indexx routine from num receipts library
      call indexx(ntot-ifirst+1, tmp_r, ind_sort)
*
*     Split the bh_dat.31 files by time
         call string_left(TCHAR,TTOT,DELTAT)
         write(OUTFILE,118) TCHAR
 118     format('bh_dat.31_',A20)
*
       OPEN (UNIT=31,STATUS='UNKNOWN',FORM='FORMATTED',FILE=OUTFILE)
*
      do i=1,nb
         ii = ind_sort(i)
         tmp_v = sqrt(vs(1,ii)**2+vs(2,ii)**2+vs(3,ii)**2)
         tmp_a = sqrt(fs(1,ii)**2+fs(2,ii)**2+fs(3,ii)**2)
         tmp_adot = sqrt(fdots(1,ii)**2+fdots(2,ii)**2+fdots(3,ii)**2)  
*          
       write (31, 98) 
     &    ttot,cmbh,i,name(ii),kstar(ii),bodys(ii),(xs(j,ii),j=1,3),
     &    tmp_r(ii),(vs(j,ii), j=1,3), tmp_v, (fs(j,ii), j=1,3), tmp_a,
     &    (fdots(j,ii), j=1,3), tmp_adot, step(ii), stepr(ii)
   98  format (1x,1p,2e13.5,2i10,i4,19(e15.6))   
      enddo 
*    Put newline character after time interval
          write (31,99) 
   99     format ()     
*
      call flush(31)
*
      return
*
      end
