!
! Program gazdag.f90
! Fortran MEX fucntion callable by gazdagmig.m to perform 
! image construction for Gazdag phase-shift migration.
!
!
! Copyright (C) 2005,  Andreas Tzanis
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.!
!
!
      subroutine mexfunction(nlhs, plhs, nrhs, prhs)
!
!     Gateway routine for subroutine gazdag.f90
!     Author   : Andreas Tzanis, 
!                Department of Geophysics, University of Athens 
!                atzanis@geol.uoa.gr
!                November 2005.
!
!
 integer*4 plhs(*), prhs(*)
 integer*4 nlhs, nrhs
 integer*4 mxCreateDoubleMatrix, mxGetPr, mxGetPi
 integer*4 mxGetM, mxGetN, mxGetString
 integer*4 dt_pr, dx_pr, vmig_pr, fkr_pr, fki_pr
 integer*4 imgr_pr, imgi_pr, ns_pr, ntr_pr
 integer*4 M, N, strstatus
 real*8    imgr,imgi,fkr,fki,dt,dx,vmig,dns,dntr
 character*8 mode
!
! CHECK FOR PROPER NUMBER OF ARGUMENTS
 if (nrhs .ne. 8) then
     call mexErrMsgTxt('GAZDAG > 8 input arguments required')
 elseif (nlhs .gt. 2) then
     call mexErrMsgTxt('GAZDAG > 2 output arguments required')
 endif
! Assign pointers to the input parameters 
 dt_pr     = mxGetPr(prhs(1))
 dx_pr     = mxGetPr(prhs(2))
 vmig_pr   = mxGetPr(prhs(3))
 fkr_pr    = mxGetPr(prhs(4))
 fki_pr    = mxGetPr(prhs(5))
 ns_pr     = mxGetPr(prhs(6))
 ntr_pr    = mxGetPr(prhs(7))
 strstatus = mxGetString(prhs(8), mode, 8)
! Take the dimensions of the input data data f-k transform  
 M = mxGetM(prhs(4))
 N = mxGetN(prhs(4))
! create matrices for the return arguments
 plhs(1) = mxCreateDoubleMatrix(M,N,0)
 plhs(2) = mxCreateDoubleMatrix(M,N,0)
 imgr_pr = mxGetPr(plhs(1))
 imgi_pr = mxGetPr(plhs(2))
! Compute the image in SUBROUTINE GAZDAG
 call gazdag(%val(imgr_pr),%val(imgi_pr),%val(dt_pr),&
      &%val(dx_pr),%val(vmig_pr),%val(fkr_pr),%val(fki_pr),&
      &%val(ns_pr),%val(ntr_pr), mode) 
 return
 end
!
 subroutine gazdag(imgr,imgi,dt,dx,vmig,fkr,fki,dns,dntr,mode) 
!
!  Gazdag Phase-shifting migration of zero-offset GPR data in 
!  homogeneous halfspaces (migration velocity vmig is constant),
!  layered halfspaces (migration velocity vmig is a function of time
!  and length(vmig) == ns
! 
!
!  NOTE : The operator "cc" is conjugated, (opposite to what books 
!         usually say), to account for the engineering definition 
!         of the Fourier kernel in MATLAB 
! Author: Andreas Tzanis, 
!         Department of Geophysics, 
!         University of Athens 
!         atzanis@geol.uoa.gr
!         November 2005
!
!
 character mode*8
 character(64) command
 integer*4  mexEvalString
 integer    ns, ntr, nx, nz, iw, ikx, ikz, nw, ntau, err_alloc
 real*8     dt, dx, dz, w0, kz0, dw, dkx, dkz, pi, vkx2, w2
 real*8      w, kx, kz, kx0, phase, dtau, dns, dntr
 real*8     ft, ftau, tmax, coss, vmigc
 real*8     fkr(*), fki(*), imgr(*), imgi(*), vmig(*)
 complex*16 cc
 complex*16, allocatable :: fk(:,:), img(:,:)
!
 ns  = idint(dns)
 ntr = idint(dntr)
! Set up physical constants 
 pi   = 3.14159268d0
 nw   = ns
 w0   = -pi/dt
 dw   = 2*pi/(dfloat(ns)*dt)
 nx   = ntr
 kx0  = -pi/dx
 dkx  = 2*pi/(dfloat(nx)*dx)        
 ntau = ns   
 dtau = dt
 ft   = 0
 ftau = ft
 tmax  = ft  + (dfloat(ntau)-1.0d0)*dtau
! allocate memory to the work variables
 allocate(fk(ns,ntr),img(ns,ntr),STAT=err_alloc)
 if (err_alloc .ne. 0) then 
     call mexErrMsgTxt('GAZDAG > Memory allocation error - Aborting')
     return
 endif

! Load the data onto the work arrays
 do i=1,ntr
    ik = (i-1)*nw +1 
    do j = 1,nw
       jk = ik + j -1
       fk(j,i)  = dcmplx(fkr(jk), fki(jk))
       img(j,i) = dcmplx(0.0D0,0.0D0)
    enddo
 enddo
! Preare the progress bar
 istat =  mexEvalString('h=waitbar(0,''Imaging FK_x -> TK_x'');')

if (mode .eq. 'uniformv' .or. mode .eq. 'UNIFORMV')  then

! Compute image - constant velocity structures
 vmigc = vmig(1)
 do iw = 1,nw			
    w  = w0 + (dfloat(iw)-1.0d0)*dw
    if (w.eq.0.0) then 
        w = 1.0d-10/dt
    endif
! Report progress  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    progress = float(iw)/float(nw)
    write(command,'(a8, f6.4, a3)') 'waitbar(', progress, ',h);'
    istat =  mexEvalString(command(1:len_trim(command)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ikx = 1,nx  
        kx  = kx0 + (dfloat(ikx)-1.0d0)*dkx;
        w2   = w*w
        vkx2 = (vmigc*vmigc * kx*kx)/4.0d0
        if (w2 .gt. vkx2) then
           phase = real(-w * dtau * dsqrt(1.0d0 - vkx2/w2))
           cc    = dconjg(dcmplx(dcos(phase),dsin(phase)))
! Accumulate image summed over all frequencies
           do itau = 1,ntau
               fk(iw,ikx)   = fk(iw,ikx) * cc
               img(itau,ikx)= img(itau,ikx) + fk(iw,ikx)
           enddo  
        else
              fk(iw,ikx) = dcmplx(0.0,0.0)
	    endif  
    enddo                                            ! ikx loop
 enddo                                               ! iw loop

elseif (mode .eq. 'layeredv' .or. mode  .eq. 'LAYEREDV') then

! Compute image - layered velocity structures
 do ikx = 1,nx                                     ! Loop over wavenumbers
     kx  = kx0 + (dfloat(ikx)-1.0d0)*dkx
! Report progress  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     prgrs = float(ikx)/float(nx)
     write(command,'(a8, f6.4, a3)') 'waitbar(', prgrs, ',h);'
     istat =  mexEvalString(command(1:len_trim(command)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do itau = 1, ntau                             ! Loop over time
         tau = ftau + (dfloat(itau)-1.0d0)*dtau
         do iw = 1,nw                              ! Loop over frequencies
             w  = w0 + (dfloat(iw)-1.0d0)*dw
             if (w .eq. 0.0d0) then 
                 w = 1d-10/dt
	         endif
             coss = 1.0 - (0.5 * vmig(itau)*kx/w)**2
             if (coss .ge. (tau/tmax)**2) then 
                phase = real(-w*dt*dsqrt(coss))
                cc    = dconjg(dcmplx(dcos(phase),dsin(phase)))
                fk(iw,ikx)   = fk(iw,ikx) * cc
             else
                fk(iw,ikx) = dcmplx(0.0d0,0.0d0)
             endif                                
             img(itau,ikx)= img(itau,ikx) + fk(iw,ikx)
         enddo                                     ! iw loop
     enddo                                         ! itau loop
 enddo                                    

endif                                              ! if uniform or layered

! Load the image onto the output array
 do i=1,ntr
  	ik = (i-1)*nw +1 
    do j = 1,nw
       jk = ik + j -1 
       imgr(jk) = dreal(img(j,i))/dfloat(nw)
       imgi(jk) = dimag(img(j,i))/dfloat(nw)
    enddo
 enddo
! Clear the work arrays
 deallocate(fk, img)
! Close the progress bar
 istat =  mexEvalString('close(h); clear h')
 istat =  mexEvalString('disp(''GAZDAG > Returning with success'');')
 return
 end
