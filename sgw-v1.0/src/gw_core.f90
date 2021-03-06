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
subroutine gw_core
  use simple_mpi, only : color_reduce_sum_c16, rank
  use gwm
  implicit none
  integer it, st
  integer,     save   :: n2=2
  integer             :: map_sp_etaxi(2)
  integer, external   :: is_it_power2
  real*8              :: wg, r, nrm_tmp
  real*8, external    :: ran_ps

  call alloc_ge_gge
  if(eta_flg) call set_seed_zeta()
  call set_ear
  call set_ge
  allocate(g(n),stat=st); call check0(st,' g allocation problem ')
  call reset_g
  call calc_gge
  call write_ge
  call allocate_calc_del_zeta  ! del, zeta: calc. in eta_rank, distributed.
  
  !
  ! all cores in a given color will project for psi representing W
  ! one core per color for eta(and xi) representing G(t<0).
  !  
  allocate(eta(n),  stat=st); call check0(st,' eta ')
  allocate(pt(n,ns),stat=st); call check0(st,' pt  ')
  call pt_eta
  call write_pt

  !
  ! finished the projection stage
  !
  !
  ! exchange from pt 
  !
  call xc_expect
  if(gamflg==1.or.gamflg==2) call gam_prep
  !
  ! now propagation
  !

  call vo_make

  if(allocated(vh0))   deallocate(vh0)
  if(allocated(pt))    deallocate(pt)
  if(allocated(del))   deallocate(del)

  !
  ! now eta_xi
  !
  if(eta_flg) then
     allocate(xi(n), stat=st); call check0(st,' xi ')
     xi=zeta-eta
  end if
  if(allocated(zeta)) deallocate(zeta)

  allocate(etaxi(n, n2), stat=st); call check0(st,' etaxi ')
  call prep_norm_etaxi
  map_sp_etaxi(:)=sp0
  if(allocated(eta)) deallocate(eta)
  if(allocated(xi))  deallocate(xi)

  !
  !
  ! now ct=<phi eta_or_xi(t) | W(t) | phi zeta >
  ! Need to sum within color
  !

  call prnt(1,'DEBUG: pre G0 propagation')

  ct=0d0
  allocate(vo(n,nspv),stat=st); call check0(st,' vo ')
  do it=0,nt
     call vo1(it)    
     wg = wgf(it); if(it==0)wg=wg/2d0
     call prepcv
     ct( it,:) = ct( it,:) + cvxi*normxi*wg
     ct(-it,:) = ct(-it,:) - cveta*normeta*wg 
     call propdt(etaxi,n,n2,map_sp_etaxi,0)

     if(mod(it,100)==0) then
       call prnt(1,'DEBUG: G0 prop after : 0,100,200,... steps') 
     end if
  enddo
  deallocate(vo)
  if(gamflg.eq.1.or.gamflg.eq.2) call color_reduce_sum_c16(ct, size(ct), ptheader_rank) ! no need in direct

  call gwplt

  if(allocated(gam))    deallocate(gam)
  if(allocated(etaxi))  deallocate(etaxi)
  deallocate(gge, ge, g, stat=st); call check0(st,' gge dealloc. ')
contains
  subroutine alloc_ge_gge
    implicit none
    if(.not.allocated(ge)) then
       allocate(ge(n,ne), stat=st); if(st/=0) stop ' ge  allocation problem '
    endif
    if(.not.allocated(gge)) then
       allocate(gge( ne), stat=st); if(st/=0) stop ' gge allocation problem '
    endif
  end subroutine alloc_ge_gge

  subroutine prepcv
    use gwm
    implicit none
    integer ie
    do ie=1,ne
       cveta(ie) = sum(conjg(etaxi(:,1))*vo(:,sp0)*ge(:,ie))*dv
       cvxi( ie) = sum(      etaxi(:,2) *vo(:,sp0)*ge(:,ie))*dv
    enddo
  end subroutine prepcv
  
  subroutine reset_g
    implicit none
    integer ie
    g = 0d0
    do ie=1,ne
       g = g + ge(:,ie)
    enddo
  end subroutine reset_g

  subroutine calc_gge
    implicit none
    integer ie
    do ie=1,ne
       gge(ie) = sum(ge(:,ie)**2.d0)*dv
    enddo
  end subroutine calc_gge

  subroutine set_ear
    implicit none
    if(ne/=1) then
       write(6,*)' problem, not set for   ne ',ne,' stopping '
       stop
    end if
    ear(:)= eorb
  end subroutine set_ear

end subroutine gw_core

