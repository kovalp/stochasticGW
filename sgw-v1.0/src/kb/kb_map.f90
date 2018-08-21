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
subroutine kb_map
  use kb_mod
  implicit none
  integer, save :: i1=1
  integer       :: ia, ma, nga, mpd, ig, rep
  integer       :: ix0,iy0,iz0
  integer       :: ix, iy, iz
  integer       :: ipx, ipy, ipz
  integer       :: mmx, mmy, mmz
  real*8        :: d2,r0(3),rg(3),rn(3)
  logical       :: found
  
  if(i1/=1) return
  if(i1==1) i1=-1

  r0 = (/ -dble(nx)/2d0*dx, -dble(ny)/2d0*dy, -dble(nz)/2d0*dz /)

  mpd=0
  if(periodic)mpd=2  ! see mpd below to unerstand.
  call check_cell(mpd)
  
  ngsml = 0
  reploop : do rep = 1,2
     ialoop : do ia=1,na
        ma = mapai(ia)
        mmx = racut_a(ma)/dx + 1;    
        mmy = racut_a(ma)/dy + 1;    
        mmz = racut_a(ma)/dz + 1;    
        if(rnk==0)write(17,*)' ia,ch,mmxyz--integer radius of checking ',ia,ch_a(ia),mmx,mmy,mmz
        if(periodic.and. 2d0*racut_a(ma)>min(nx*dx, ny*dy, nz*dz)) then
           stop ' ERROR: stopping. periodic, & cell too small, images of nuclei in NL overlap '
        endif
        nga = 0
        pzloop : do ipz=-mpd,mpd
           pyloop : do ipy=-mpd,mpd
              pxloop : do ipx=-mpd,mpd
                 
                 rn = cnt(:,ia) + (/ ipx*nx*dx, ipy*ny*dy, ipz*nz*dz /)
                 
                 if(nx>1) then; ix0 = (rn(1) - r0(1))/dx; else ; ix0=1; endif
                 if(ny>1) then; iy0 = (rn(2) - r0(2))/dy; else ; iy0=1; endif
                 if(nz>1) then; iz0 = (rn(3) - r0(3))/dz; else ; iz0=1; endif
                          
                 izloop : do iz=max(1,iz0-mmz),min(nz,iz0+mmz)
                    iyloop : do iy=max(1,iy0-mmy),min(ny,iy0+mmy)
                       ixloop: do ix=max(1,ix0-mmx),min(nx,ix0+mmx)
                          found = .false.
                          ig = ix+(iy-1)*nx+(iz-1)*nx*ny
                          rg = r0        + (/ (ix-1)*dx, (iy-1)*dy, (iz-1)*dz /) 
                          d2 = sum((rg-rn)**2)
                          nearby: if(d2<racut_a(ma)**2) then
                             if(found)stop ' ERROR: 2 atomic images overlapping in NL '
                             found=.true.
                             nga = nga+1
                             if(rep==2) then
                                mapkbg(  nga,ia) = ig
                                rgn(:,nga,ia) = rg(:)-rn(:) 
                             end if
                          end if nearby
                       end do ixloop
                    end do iyloop
                 end do izloop
              end do pxloop
           end do pyloop
        end do pzloop
        ngs(ia)=nga
     end do ialoop
     
     rep1 : if(rep==1) then
        ngsml = maxval(ngs(1:na))         
        if(rnk==0) write(17,*)' minval ngs, maxval(=ngsml) ',minval(ngs(1:na)),ngsml; 
        allocate(mapkbg(    ngsml, na), stat=st); call check(st,0,' map_kb_g ' )
        allocate(rgn(     3,ngsml, na), stat=st); call check(st,0,' rgn      ' )
        mapkbg = -999999
     end if rep1

  end do reploop

  call flush(17)
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
end subroutine kb_map
