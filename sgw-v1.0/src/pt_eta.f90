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
subroutine pt_eta
  use gwm
  implicit none

  if(det_tddft) then
     call  pt_det_tddft
     call eta_det_tddft
     return
  end if

  call pt_stoch  ! prepare stochastic pt, to be projected.

  if(det_dft) then
     call pt_eta_effic
  else
     stop ' det_dft problem '
  endif

  
end subroutine pt_eta

