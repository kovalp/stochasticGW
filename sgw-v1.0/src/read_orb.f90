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
subroutine read_orb
  use gwm
  implicit none
  if(orb_indx==-1.or.(.not.det_dft)) then
     call read_orb_1
  else
     call read_orb_fromwf
  endif
end subroutine read_orb



