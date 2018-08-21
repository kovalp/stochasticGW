!                     
!                     
!    The routine(s) in this file are a part of the  
!                     StochasticGW                  
!    suite, developed 2012-2018, and copyrighted    
!    to the authors of StochasticGW , see:          
!                                                   
!                 www.stochasticgw.com              
!                                                   
!   If you use or modify any part of this routine   
!   the header should be kept, unmodified.          
!                                                   
!                                                   
!                                                   
subroutine exchange_stoch
  use gwm, only : ne, ns, ge, n, ns_blk, pt, dv, normpt, exce
  implicit none
  integer ie, is, st
  complex*16  c1
  complex*16, allocatable :: crho(:), cpot(:)
  allocate(crho(n),cpot(n),stat=st); if(st/=0) stop' crho problems'

  exce = 0d0
  ieloop : do ie=1,ne
     isloop: do is=1,ns
        
        crho(:) = 1/sqrt(dv)*ge(:,ie)*pt(:,is)  ! 1/sqrt(dv) since somewhere in gw_core we 
                                                !   defined  ge = sqrt(dv)*orbital
        call cvh_sub( crho, cpot, n)
        
        c1 = dv*sum(conjg(crho)*cpot)*normpt(is)**2
        call check_real_1(c1,' exce ')
        
        exce(ie) = exce(ie) - dble(c1)/ dble(ns_blk)
     enddo isloop
  enddo ieloop
  deallocate(crho,cpot)
end subroutine exchange_stoch

