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
subroutine get_vxc_spn(dens, n, nsp, vxc)  ! include core_correction
  use gwm, only : xc_type, vxc0, rho_core
  implicit none
  integer n,nsp,sp, st
  real*8  dens(n, nsp)
  real*8   vxc(n, nsp)
  real*8, allocatable :: vx(:,:), vc(:,:), dens_tot(:,:)
  if(xc_type==1) then
     vxc = vxc0
  else
     allocate(vx(n, nsp), vc(n, nsp), dens_tot(n, nsp),stat=st); if(st/=0) stop ' vx vc  '
     vx = 0d0
     vc = 0d0
     do sp=1,nsp
        dens_tot(:,sp) = dens(:,sp) + rho_core(:)/dble(nsp)
     enddo
     call get_vx_spin(    dens_tot, n, nsp, vx);
     call get_vcn_spin(   dens_tot, n, nsp, vc);
     vxc = vx + vc
     deallocate(vx,vc,dens_tot)
  end if
end subroutine get_vxc_spn
