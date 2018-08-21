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
subroutine calc_cvkb(nxb,nyb,nzb,nrp,rp,vr,cvkb)
  use kb_mod, only : dx,dy,dz,dv,periodic
  implicit none
  integer nxb, nyb,nzb,nrp,ix,iy,iz,ir,st
  real*8  rp(nrp),vr(nrp), rv(3),vrlg,rscalar,rlg
  real*8  aa_here,charge_here
  real*8, allocatable, dimension(:) :: vr2, rplog
  complex*16 cvkb( nxb, nyb, nzb )
  
  if(periodic) then
     aa_here = max(7d0/min(nxb*dx,nyb*dy,nzb*dz), 5d0/rp(nrp))
     write(17,*)' aa_here ',aa_here
  end if

  !
  ! OK, now we have vr on a long grid.  Now we need to put in on a big 3-d grid.
  !  First, fit to spline

  allocate(vr2(nrp), rplog(nrp), stat=st); call check0(st,' vr2 ')

  do ir=1,nrp
     rplog(ir) = log(rp(ir))
  enddo

  call SPLINE_dn(rplog,vr,nrp,vr2)

  charge_here = vr(nrp)*rp(nrp)

  do iz=1,nzb
     do iy=1,nyb
        do ix=1,nxb
           rv = (/ (ix-1-nxb/2)*dx, (iy-1-nyb/2)*dy, (iz-1-nzb/2)*dz /)
           rscalar = max(1d-10,sqrt(sum(rv**2)))
           rlg = log(rscalar)

           if(rscalar<rp(1)) then
              vrlg = vr(1)
           elseif(rscalar>rp(nrp)) then
               vrlg = charge_here/rscalar
            else
               call SPLINT_dn(rplog,vr,vr2,nrp,rlg,vrlg) 
            end if
            
            if(periodic) vrlg = vrlg - charge_here * erf(aa_here*rscalar)/rscalar 
            cvkb(ix,iy,iz) = vrlg
        end do
     end do
  end do
  ! now we have the 3d potential in r space. convert to k:
  call fft3d_forward_many( nxb,nyb, nzb,1,cvkb)
  cvkb = cvkb*dv
  call shift_cvkb
  if(periodic) call add_long_range_periodic
  call check_real(cvkb, size(cvkb))
  deallocate(vr2,rplog)
contains
  subroutine interp_lin(rplog,vr,nrp,rlg,vrlg)
    use kb_mod, only : dv
    implicit none
    integer nrp, i, j,k
    real*8    vr(nrp)
    real*8 rplog(nrp)
    real*8 rlg, vrlg
    i=1
    j=nrp 
    do while(j>i+1)
       k=(i+j)/2
       if(rlg<rplog(k)) then
          j=k
       else
          i=k
       endif
    end do
    vrlg = (vr(i)*(rplog(j)-rlg)+vr(j)*(rlg-rplog(i)))/(rplog(j)-rplog(i)) 
  end subroutine interp_lin

  subroutine shift_cvkb
    use kb_mod, only : dx, dy, dz
    implicit none
    integer ikx,iky,ikz
    real*8  dkx,dky,dkz,kx,ky,kz,pi,r0(3)
    complex*16, parameter :: ci = (0d0,1d0)
    pi=dacos(-1d0)
    r0 = (/ -nxb/2d0*dx, -nyb/2d0*dy, -nzb/2d0*dz /)
    do ikz=1,nzb
       do iky=1,nyb
          do ikx=1,nxb
             dkx = 0d0 ; if(dx>1e-8) dkx = 2.d0*pi/(dble(Nxb)*dx)
             dky = 0d0 ; if(dy>1e-8) dky = 2.d0*pi/(dble(Nyb)*dy)
             dkz = 0d0 ; if(dz>1e-8) dkz = 2.d0*pi/(dble(Nzb)*dz)
             
             kx = dble(ikx-1)* dkx
             ky = dble(iky-1)* dky
             kz = dble(ikz-1)* dkz
             
             if(kx>pi/dx) kx = kx-2d0*pi/dx  
             if(ky>pi/dy) ky = ky-2d0*pi/dy
             if(kz>pi/dz) kz = kz-2d0*pi/dz
             
             cvkb(ikx,iky,ikz) = &
             cvkb(ikx,iky,ikz) * exp(-ci*( r0(1)*kx+r0(2)*ky+r0(3)*kz))
          enddo
       enddo
    end do
  end subroutine shift_cvkb

  subroutine add_long_range_periodic
    use kb_mod, only : dx, dy, dz
    implicit none
    integer ikx,iky,ikz
    real*8  dkx,dky,dkz,kx,ky,kz,pi,r0(3),k2
    complex*16, parameter :: ci = (0d0,1d0)
    pi=dacos(-1d0)
    r0 = (/ -nxb/2d0*dx, -nyb/2d0*dy, -nzb/2d0*dz /)
    do ikz=1,nzb
       do iky=1,nyb
          do ikx=1,nxb
             dkx = 0d0 ; if(dx>1e-8) dkx = 2.d0*pi/(dble(nxb)*dx)
             dky = 0d0 ; if(dy>1e-8) dky = 2.d0*pi/(dble(nyb)*dy)
             dkz = 0d0 ; if(dz>1e-8) dkz = 2.d0*pi/(dble(nzb)*dz)
             
             kx = dble(ikx-1)* dkx
             ky = dble(iky-1)* dky
             kz = dble(ikz-1)* dkz
             
             if(kx>pi/dx) kx = kx-2d0*pi/dx  
             if(ky>pi/dy) ky = ky-2d0*pi/dy
             if(kz>pi/dz) kz = kz-2d0*pi/dz
             
             k2 = kx**2+ky**2+kz**2
             if(k2>1d-4) then
                cvkb(ikx,iky,ikz) = &
                cvkb(ikx,iky,ikz) + &
                                  charge_here * 4d0*pi/k2*exp(-k2/4d0/aa_here**2)
             else
                cvkb(ikx,iky,ikz) = 0d0
             endif
          enddo
       enddo
    enddo


  end subroutine add_long_range_periodic
end subroutine calc_cvkb

