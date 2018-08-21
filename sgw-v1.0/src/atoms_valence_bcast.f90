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
subroutine bcast_valence
  use atoms, only : valch_a, ch_a, valence_list
  use gwm,   only : na
  use simple_mpi, only : bcast_i
  implicit none
  integer i
  call bcast_i(valence_list,size(valence_list),0)
  do i=1,size(valch_a)
     valch_a(i) = valence_list(ch_a(i))
  enddo
end subroutine bcast_valence
