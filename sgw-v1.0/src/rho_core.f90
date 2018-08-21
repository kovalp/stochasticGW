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
subroutine make_rho_core
  use kb_mod
  use simple_mpi, only : bcast_r8
  use gwm, only : rho_core
  use ppm, only : core_correction, rho_core_a, nrppmx
  implicit none
  integer, save :: i1=1
  integer       :: ia, ma,  mpd, ig, i
  integer       :: ix0,iy0,iz0
  integer       :: ix, iy, iz
  integer       :: ipx, ipy, ipz
  integer       :: mmx, mmy, mmz
  real*8        :: d1,d2,r0(3),rg(3),rn(3),rtop 
  real*8, parameter :: toll_rho = 1d-12
  
  if(i1/=1) return
  if(i1==1) i1=-1
  r0 = (/ -dble(nx)/2d0*dx, -dble(ny)/2d0*dy, -dble(nz)/2d0*dz /)
  mpd=0
  if(periodic)mpd=2  ! see mpd below to unerstand.
  call check_cell(mpd)
  
  if(.not.allocated(rho_core)) then
     allocate(rho_core(ng), stat=st)
     call check0(st,' rho_core ')
  else
     call check(size(rho_core),ng,' rho_core ')
  end if
  rho_core = 0d0
  
  ialoop : do ia=1,na
     ma = mapai(ia)
     rtop = -99d0
     core_check : if(core_correction(ma)) then
        
        do i=1,nrppmx
            if(rho_core_a(i,ma)>toll_rho) then
              rtop = max(rtop,rrpp(i,ma))
           end if
        enddo
        call check_r_le(0d0,rtop,' 0-rtop ')
       
        mmx = rtop/dx + 1;    
        mmy = rtop/dy + 1;    
        mmz = rtop/dz + 1;    
        if(rnk==0)write(17,*)' ia,ch,rtop                        ',ia,ch_a(ia),rtop
        if(rnk==0)write(17,*)' ia,ch,mmxyz incusuion for rho_cut ',ia,ch_a(ia),mmx,mmy,mmz
        if(periodic.and. 2d0*rtop>min(nx*dx, ny*dy, nz*dz)) then
           stop ' ERROR: stopping. periodic, & cell too small, images of nuclei in NL overlap '
        endif
        
        pzloop : do ipz=-mpd,mpd
        pyloop : do ipy=-mpd,mpd
        pxloop : do ipx=-mpd,mpd
           
           rn = cnt(:,ia) + (/ ipx*nx*dx, ipy*ny*dy, ipz*nz*dz /)
           
           if(nx>1) then; ix0 = (rn(1) - r0(1))/dx; else ; ix0=1; endif
           if(ny>1) then; iy0 = (rn(2) - r0(2))/dy; else ; iy0=1; endif
           if(nz>1) then; iz0 = (rn(3) - r0(3))/dz; else ; iz0=1; endif
              
           izloop : do iz=max(1,iz0-mmz),min(nz,iz0+mmz)
           iyloop : do iy=max(1,iy0-mmy),min(ny,iy0+mmy)
           ixloop : do ix=max(1,ix0-mmx),min(nx,ix0+mmx)
              ig = ix+(iy-1)*nx+(iz-1)*nx*ny
              rg = r0 + (/ (ix-1)*dx, (iy-1)*dy, (iz-1)*dz /) 
              d2 = sum((rg-rn)**2)
              nearby: if(d2.le.rtop**2) then
                 d1 = sqrt(d2)
                 call interp_rho_add(d1, rho_core(ig), rrpp(1,ma), rho_core_a(1,ma), nrpp(ma))
              end if nearby
           end do ixloop
           end do iyloop
           end do izloop

        end do pxloop
        end do pyloop
        end do pzloop
     endif core_check
  enddo ialoop

  if(rnk==0) then
     write(17,*)' max rhocore, integral rhocore ',maxval(rho_core), sum(rho_core)*dv
     call flush(17)
  endif
  call bcast_r8(rho_core, size(rho_core),0)
contains
  subroutine check_cell(mpd)
    use simple_mpi, only : sync_mpi, rank
    implicit none
    integer mpd, i
    real*8  rL(3),amn(3),amx(3),fctr
    rL = (/ nx*dx, ny*dy, nz*dz /)
    amn = (/ minval(cnt(1,:)),minval(cnt(2,:)),minval(cnt(3,:)) /)
    amx = (/ maxval(cnt(1,:)),maxval(cnt(2,:)),maxval(cnt(3,:)) /)
    
    rnk0: if(rank==0) then
       write(17,*)' Checking atoms vs. cell '
       write(17,*)' atoms_minmax, cell length '
       do i=1,3
          write(17,*)' direction: ',i,' atoms_minmax ',real(amn(i)),real(amx(i)),&
               ' cell-length ',real(rl(i))
       enddo
       prdc: if(periodic) then
          call check(mpd,2,' mpd, 2 ')
          write(17,*)' periodic,     atoms cell should be within -fctr*cell,fctr*cell, fctr=1.5 ' 
          fctr = 1.5d0
       else
          call check(mpd,0,' mpd, 0 ')
          write(17,*)' non-periodic, atoms cell should be within -fctr*cell,fctr*cell, fctr=0.5 ' 
          fctr = 0.5d0
       end if prdc
       do i=1,3
          call check_r_le(-fctr*rL(i),amn(i), ' -fctr*rL, minval_cnt ')
          call check_r_le(amx(i), fctr*rL(i), '  maxval_cnt, fctr*rL ')
       enddo
    end if rnk0
    call sync_mpi
  end subroutine check_cell
end subroutine 
