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
subroutine exch_sum
  use simple_mpi, only : allsum_r8
  use gwm
  implicit none  
  if(trace.and..not.rdexch) stop ' trace requries rdexch '
  if(rdexch) then
     call check(ne,1,' ne 1 ' )
     exce(:) = exchange_value
  else if(det_dft) then
     call allsum_r8(exce, size(exce))
  else
     exce = 0d0
  end if
end subroutine exch_sum
