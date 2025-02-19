!
!     Migration of zero-offset GPR data 
!
! *** The program assumes that data has already been preprocessed and 
!     transformed to the F-K domain by the calling program (MATLAB or OCTAVE).
!     MIGRATE only performs the integrations, contructs the image and exports
!     the image for post-processing by MATLAB or OCTAVE.
! *** Migration method/procedure is determined by the (imported) control 
!     keyword MOPTION. Available options are:
!
! 1.  GAZDAG PHASE-SHIFTING MIGRATION: 
! ==> For Constant Velocity Structures MOPTION = 'CVMIGG' and the migration 
!     velocity "vmigc" is read as a scalar. 
! ==> For Layered Velocity Structures  MOPTION = 'LVMIGG' and migration 
!     velocities are imported as a vector "vmigv" of velocity vs time.
!     => The code (nesting of loops) is different between the two options and
!        CASE CVMIGG is considerably faster than CASE LVMIGG
!     ** The upward continuation operator "cc" is conjugated (opposite to what 
!        books usually say), to account for (engineering) definition of the 
!        Fourier kernel in MATLAB / OCTAVE
!
! 2.  STOLT F-K MIGRATION: 
! ==> For Constant Velocity Structures MOPTION = 'CSTOLT' and the migration 
!     velocity "vmigc" is read as a scalar. 
!     ** The "moved out" frequency w is taken with sign opposite to that 
!        reported in most books, to account for the conjugate definition of 
!        the Fourier kernel in MATLAB 
! ==> Acknowledgement: Stolt code was adapted from John Claerbout's RATFOR 
!     routine found in "Imaging the Earth's Interior", Chapter 3
!
!     Author   : Andreas Tzanis, 
!                Department of Geophysics, University of Athens 
!     Created  : 9 October 2003
!

$debug
character moption*6
character buffer*255, infile*255, outfile*255
integer*2 status
integer   ns, ntr, nx, nz, iw, ikx, ikz, nw, err_alloc
real      dt, dx, dz, w0, kz0, dw, dkx, dkz, pi, vmigc, w, kx, kz, signum
complex   cc
! Uses dynamically allocated arrays : Feature will not work for FORTRAN 77 
real, allocatable :: vmigv(:), rr(:,:), ii(:,:)
complex, allocatable :: up(:,:), img(:,:)
pi   = 3.14159268

! Begin main
! input and output file names are passed as arguments
call getarg(int2(1),buffer,status)
infile = adjustl(buffer(1:status))
call getarg(int2(2),buffer,status)
outfile = adjustl(buffer(1:status))
! Input file has been prepared by calling program (MATLAB or OCTAVE) 
open(1,file=infile,form='binary')
read(1) moption
! ns   : No of samples in the time / frequency axis
! ntr  : No of traces in profile
if (moption.eq.'cvmigg'.or.moption.eq.'CVMIGG') then
   read(1) ns, ntr, dt, dx, vmigc
   allocate(rr(ns,ntr),ii(ns,ntr),up(ns,ntr),img(ns,ntr),STAT=err_alloc)
   if (err_alloc .ne. 0) then 
        print *, " Memory allocation error"
	    stop
   endif

else if (moption.eq.'lvmigg'.or.moption.eq.'LVMIGG') then
   read(1) ns, ntr, dt, dx
   allocate(rr(ns,ntr),ii(ns,ntr),vmigv(ns),up(ns,ntr),img(ns,ntr),STAT=err_alloc)
   if (err_alloc .ne. 0) then 
        print *, " Memory allocation error"
	    stop
   endif
! get Velocity vs time for layered structures
   do i=1,ns
       read(1) vmigv(i)
   enddo

else if (moption.eq.'cstolt'.or.moption.eq.'CSTOLT') then
   read(1) ns, ntr, dt, dx, dz, vmigc
   allocate(rr(ns,ntr),ii(ns,ntr),up(ns,ntr),img(ns,ntr),STAT=err_alloc)
   if (err_alloc .ne. 0) then 
        print *, " Memory allocation error"
	    stop
   endif

endif
! Import data Get FFT'ed to the F-K domain by MATLAB or OCTAVE 
do i=1,ns
read(1) (rr(i,j), j=1,ntr)
enddo
do i=1,ns
read(1) (ii(i,j), j=1,ntr)
enddo
close(1,STATUS='delete')
do i=1,ns
do j=1,ntr
up(i,j) = cmplx(rr(i,j),ii(i,j))
enddo
enddo
! Free some memory
deallocate( rr, ii )

! Set up physical constants 
nw    = ns;     w0  = -pi/dt;      dw  = 2*pi/(ns*dt);
nx    = ntr;   kx0  = -pi/dx;     dkx  = 2*pi/(nx*dx);        

percent = 0.
ddec    = 0.
!!!!!!!!!!  PHASE-SHIFTING METHOD - constant velocity structures  !!!!!!!!!!
if (moption.eq.'cvmigg'.or.moption.eq.'CVMIGG') then

! Set upadditional physical constants
ntau  = ns   
 dtau  = dt;   

 do iw = 1,nw										  ! Loop over frequencies
     w  = w0 + (iw-1)*dw
     if (w==0.0) w = 1e-10/dt
!!!!! Report progress
      percent= 100.0*(nw-(nw-iw+1))/nw
      if (percent > ddec) then
          write(*,'(''  Done '',i3,'' % '')')	int(percent)
          ddec = ddec + 10.0
      endif
     do ikx = 1,nx                                    ! Loop over wavenumbers
         kx  = kx0 + (ikx-1)*dkx;
         w2   = w*w
         vkx2 = (vmigc*vmigc * kx*kx)/4.
         if (w2 > vkx2) then
             phase = real(-w * dtau * sqrt(1.- vkx2/w2))
             cc    = conjg(cmplx(cos(phase),sin(phase)))
! Accumulate image summed over all frequencies
             do itau = 1,ntau
                 up(iw,ikx)   = up(iw,ikx) * cc
                 img(itau,ikx)= img(itau,ikx) + up(iw,ikx)
             enddo                                    ! itau loop
         else
                 up(iw,ikx) = cmplx(0.0,0.0)
         endif                                        ! if loop
     enddo                                            ! ikx loop
 enddo                                                ! iw loop

!!!!!!!!!!  PHASE-SHIFTING METHOD - layered  velocity structures  !!!!!!!!!!
else if (moption.eq.'lvmigg'.or.moption.eq.'LVMIGG') then

! Set up additional physical constants 
ntau  = ns;   
 dtau  = dt;   
  ft    = 0;     ftau = ft;        
   tmax  = ft  + (ntau-1)*dtau;

do ikx = 1,nx                                         ! Loop over wavenumbers
    kx  = kx0 + (ikx-1)*dkx;
!!!!! Report progress
    percent= 100.0*(nx-(nx-ikx+1))/nx
    if (percent > ddec) then
        write(*,'(''  Done '',i3,'' % '')')	int(percent)
        ddec = ddec + 10.0
    endif
! Loop over migrated times  - accumulate image by summing over all frequencies
    do itau = 1, ntau
        tau = ftau + (itau-1)*dtau;
       do iw = 1,nw                                  ! Loop over frequencies
           w  = w0 + (iw-1)*dw
		   if (w==0.0) w = 1e-10/dt
           coss = 1.0 - (0.5 * vmigv(itau) * kx / w )**2
           if (coss >= (tau/tmax)**2) then 
               phase = real(-w*dt*sqrt(coss))
               cc    = conjg(cmplx(cos(phase),sin(phase)))
               up(iw,ikx)   = up(iw,ikx) * cc
           else
               up(iw,ikx) = cmplx(0.0,0.0)
           endif                                     ! if loop
           img(itau,ikx)= img(itau,ikx) + up(iw,ikx)
       enddo                                         ! iw loop
    enddo                                            ! itau loop
enddo                                                ! ikx loop

!!!!!!!!!!  STOLT F-K METHOD - constant velocity     !!!!!!!!!!!!!!
else if (moption.eq.'cstolt'.or.moption.eq.'CSTOLT') then

! Set up additional physical constants 
nz    = ns    
 kz0  = -pi/dz     
  dkz  = 2*pi/(nz*dz)         

do ikz = 1, nz                   ! Zero the first wavenumbers
    img(ikz,1) = 0.
enddo
do ikx = 1, nx
    img(1,ikx) = 0.
enddo
do ikx = 2, nx                    ! Loop over horizontal wavenumbers 
    kx = kx0 + dkx*(ikx-1)
!!!!! Report progress
      percent= 100.0*(nx-(nx-ikx+1))/nx
      if (percent > ddec) then
          write(*,'(''  Done '',i3,'' % '')')	int(percent)
          ddec = ddec + 10.0
      endif
    do ikz = 2, nz;               ! Loop over vertical (depth) wavenumbers
        kz = kz0 + dkz*(ikz-1)
!        w = -signum(kz) * sqrt( kx*kx + kz*kz)*vmigc /2
        w = signum(kz) * sqrt( kx*kx + kz*kz)*vmigc /2 
! Change variable and interpolate F --> Kz
        call cinterp1(w, nw, w0, dw, up(1,ikx), img(ikz,ikx) )
		img(ikz,ikx) = img(ikz,ikx)*(vmigc*abs(kz)/sqrt(kx*kx + kz*kz))
        
    enddo                         ! ikz loop
enddo                             ! ikx loop

endif                                                ! MOPTION loop

! Done - export the image for post-processing 
! In case of the phase-shifting method scale by nw as if it were a proper FT
! In case of the Stolt method do not scale - will be done in the calling program 
if (moption.eq.'cstolt'.or.moption.eq.'CSTOLT')    nw = 1
open(1,file=outfile,form='binary')
do i=1,ntr
write(1) (real(img(j,i)/float(nw)), j=1,ns)
enddo
do i=1,ntr
write(1) (aimag(img(j,i)/float(nw)), j=1,ns)
enddo
close(1)

end

subroutine cinterp1(x , nx, x0, dx, cb, cbx)
! Interpolates from 'dip moved out' frequency to Kz 
! *** Code adapted from John Claerbout's RATFOR routine 
!     found in "Imaging the Earth's Interior", Chapter 3
integer   ix, ixc, nx
real      x, xc, x0, dx, fraction
complex   cb(nx), cbx
xc        = (x-x0) / dx
ixc       = xc 
fraction  = xc - ixc
ix        = 1 + ixc
if (ix.lt.1) then
    cbx = cb(1)
elseif (ix+1.gt.nx) then
    cbx = cb(nx)
else
    cbx = (1. - fraction) * cb(ix) + fraction * cb(ix+1);
endif
return
end

real function signum(x)
! *** Code adapted from John Claerbout's RATFOR routine 
!     found in "Imaging the Earth's Interior", Chapter 3
real x
if       (x > 0) then 
    signum =  1.
elseif   (x < 0) then 
    signum = -1.
else
    signum = 0.
endif
return
end
    
