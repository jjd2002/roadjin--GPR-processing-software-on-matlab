!
! Program radar2d.f90
! Fortran MEX fucntion callable by splitstep2dmodel.m, to 
! do synthetic GPR modelling with the split-step method.
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
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
!
! The gateway routine for radar2d.f90
!
! Author : Andreas Tzanis, 
!          Department of Geophysics, 
!          University of Athens 
!          atzanis@geol.uoa.gr
!          November 2005
!
!
integer*4 plhs(*), prhs(*)
integer*4 nlhs, nrhs
integer*4 mxCreateDoubleMatrix, mxGetPr
integer*4 mxGetM, mxGetN
integer*4 dt_pr, dx_pr, K_pr, M_pr, q_pr, ntr_pr
integer*4 nt_pr, nz_pr, dz_pr, fc_pr 
integer*4 band_pr, dm_pr, dum_pr
real*8    dntr,dnt,dx,dnz,dz,dt,dfc,band 

! CHECK FOR PROPER NUMBER OF ARGUMENTS
if (nrhs .ne. 12) then
   call mexErrMsgTxt('RADAR2D requires 12 input arguments')
elseif (nlhs .gt. 1 .or. nlhs .eq. 0) then
   call mexErrMsgTxt('RADAR2D requires 1 output argument')
endif
! Assign pointers to the input parameters 
K_pr    = mxGetPr(prhs(1))
M_pr    = mxGetPr(prhs(2))
Q_pr    = mxGetPr(prhs(3))
ntr_pr  = mxGetPr(prhs(4))
nt_pr   = mxGetPr(prhs(5))
dx_pr   = mxGetPr(prhs(6))
nz_pr   = mxGetPr(prhs(7))
dz_pr   = mxGetPr(prhs(8))
dt_pr   = mxGetPr(prhs(9))
fc_pr   = mxGetPr(prhs(10))
band_pr = mxGetPr(prhs(11))
dum_pr  = mxGetPr(prhs(12))
! create matrices for the return arguments
M = mxgetm(prhs(12))
N = mxgetn(prhs(12))
plhs(1) = mxCreateDoubleMatrix(M,N,0)
dm_pr= mxGetPr(plhs(1))     ! the model radargram
!
!Call radar2d
!Note that the call to radar2d.f90 contains a dummy variable that goes 
!unused. This is a dirty and memory consuming, but working method by which  
!to pass dynamically the dimensions of the model radargram from the caller 
!routine. For hitherto unspecified reasons, more elegant approaches 
!involving "mxgetscalar" and "mxcopyptrtointeger4" among others, could
!not extract the value of ntr (# traces) and nt (# time samples)
!from their respective pointers ...
call radar2d(%val(dm_pr), %val(K_pr), %val(M_pr), %val(Q_pr),&
   &%val(ntr_pr), %val(nt_pr), %val(dx_pr), %val(nz_pr), %val(dz_pr),&
   &%val(dt_pr), %val(fc_pr), %val(band_pr), %val(dum_pr))
!call mxdestroyarray(K_pr)
!call mxdestroyarray(M_pr)
!call mxdestroyarray(Q_pr)
!call mxdestroyarray(dum_pr)
return
end
!
subroutine radar2d(model,per,mu,q,dntr,dnt,dx,dnz,dz,dt,dfc,band,dummy)
!
!********************************************************************
!*                SPLIT-STEP 2D GPR MODELLING                       *	
!*         After Grandjean, G.and Durand, H., 1999,                 *
!*         Computers and Geosciences25 141-149.                     *
!*         Based on program radar2d4.c by G. Grandjean (1997)       *
!*         Coded by : Andreas Tzanis, Department of Geophysics,     *
!*                    University of Athens, atzanis@geol.uoa.gr     *
!********************************************************************
!
implicit none
real(8),parameter    :: pi  = 3.141592653589793d0
real(8),parameter    :: mu0 =4d-7*pi
real(8),parameter    :: eps0=8.8592d-12
integer(4),parameter :: mhz = 1000000
! - - - local declarations - - -
integer(2) :: status
integer(4) :: mexEvalString
integer(4) :: npfft, alloc_err, m2, istat
integer(4) :: n,ikx,iw,iw2,iwc,ix,it,iz,iw1, i, j, k
integer(4) :: f1,f2,nx2,nt2,nt,nz,ntr,fc,fc2 
integer(4) :: tet1,tet2
real*4     :: progress
real*8     :: kz,ctr,radpattern,difc,tet,nq,wr,expq,vomegaxz,bomegaxz
real*8     :: aomegaxz,cref,sigma,kx,den,rshift,den1,den2,source  
real*8     :: dx, dz, dt, dntr, dnt, dnz, dfc, band, bomega, aomega 
real*8     :: per(*), mu(*), q(*), model(*), dummy(*)
complex*16 :: shift,phase,denkc,denc,denkz,zdif,komegaxz
complex*16 :: shift1,shift2,cref1,cref2
real*8,allocatable     :: v(:,:),qtmp(:,:),vmean(:),qm(:),omega(:)
real*8,allocatable     :: g(:),waven(:),rwork(:),iwork(:)
complex*16,allocatable :: cp(:,:),rc(:,:),komega(:)
character(128)         :: matcommand
!
ntr  = idint(dntr)                   ! # traces
nt   = idint(dnt)                    ! # samples (time)
nz   = idint(dnz)                    ! # depth (depth)
fc   = idint(dfc)                 ! Antenna central frequency
nx2  = npfft(int(ntr*1.1),2*ntr) 
nt2  = npfft(int(nt*1.1),2*nt)  
!**** Determine frequency range  
tet1 = 60
tet2 = 80
wr   = 2.0*pi*fc*1.0e6
dt   = dt/1.0e9
f1   = int(fc*mhz*0.25)
f2   = int(fc*mhz*2)
fc2  = int(fc*mhz)
iw1  = int(nt2*dt*f1)
iw2  = int(nt2*dt*f2)
iwc  = int(nt2*dt*fc2)
if (band < 0.0 .or. band > 0.0 ) then
   iw1 = 1
   iw2 = nt2/2+1   
end if
!report status
matcommand = "disp('> Computed parameters, OK.')"
istat =  mexEvalString(matcommand(1:len_trim(matcommand)))

!**** ALLOCATE WORK ARRAYS ****
allocate(cp(nt2,nx2),rc(nz,nx2),komega(nt2),stat=alloc_err)
allocate(v(nz,ntr),qtmp(nz,ntr),stat=alloc_err)
allocate(vmean(nz),qm(nz),g(nt2),waven(nx2),omega(nt2),stat=alloc_err)
m2 = max(nx2,nt2)
allocate(rwork(m2),iwork(m2),stat=alloc_err) 
!**** initialize complex arrays
do ix = 1,nx2
!   rwork(ix)   = 0.0
!   iwork(ix)   = 0.0
   do iz = 1,nz
      rc(iz,ix) = cmplx(0.0,0.0)
   end do
   do iz = 1,nt2
     cp(iz,ix) = cmplx(0.0,0.0)
   end do
end do
!report status
matcommand = "disp('> Memory allocation OK.')"
istat =  mexEvalString(matcommand(1:len_trim(matcommand)))

!**** COMPUTE VELOCITY STRUCTURE AND REFLECTION COEFFICIENTS ****
do ix = 1,ntr
   i = (ix-1)*nz +1 
   do iz = 1,nz
      j = i + iz -1
! velocity structure
      v(iz,ix) = 1.0/(dsqrt(mu0*mu(j)*eps0*per(j))*  &
                   dcos((pi/4.0d0)*(1.0d0-((2.0d0/pi)*datan(q(j))))))
	  qtmp(iz,ix) = sngl(q(j))
! reflection coefficients
      if (iz > 1) then 
         cref1 = cdsqrt( dcmplx(mu(j), 0.0) / &
                 dcmplx( per(j)*dsin(datan(q(j))), &
	             dcos(per(j)*datan(q(j))) ))
         cref2 = cdsqrt( dcmplx(mu(j-1),0.0) / &
                 dcmplx( per(j-1)*dsin(datan(q(j-1))), &
                 dcos(per(j-1)*datan(q(j-1))) ))
         rc(iz,ix) = (cref1 - cref2)/(cref1 + cref2)
      endif
   enddo
enddo
!**** compute mean velocity and quality factor
call mean_velocity(vmean, v, ntr, nz)
call mean_q(qm, qtmp, ntr, nz)
!report status
matcommand = "disp('> Computed velocity structure, OK.')"
istat =  mexEvalString(matcommand(1:len_trim(matcommand)))

!**** Compute angular frequencies and wavenumbers
do ikx = 1,nx2
   kx = dfloat(ikx)*2.0d0*pi/dfloat(nx2)
   if (kx > pi) kx = 2.0d0*pi-kx
   waven(ikx) = kx/dx
end do
do iw = iw1,iw2
   omega(iw) = (dfloat(iw)*2.0d0*pi/dfloat(nt2)) / dt
end do
if (band > 0.0) then
   sigma = dfloat(iw2-iw1)*band
   do iw = iw1,iw2
       g(iw) = dexp(-1.0*(dfloat(iw-iwc)/sigma)*(dfloat(iw-iwc)/sigma))
   end do
end if
!report status
matcommand = "disp('> Computed frequencies and wavenumbers, OK.')"
istat =  mexEvalString(matcommand(1:len_trim(matcommand)))

!**** BEGIN FORWARD COMPUTTIONS
!**** Transform X -> kx
do iz = 1,nz
   do ix = 1,ntr
	  rwork(ix) = dreal(rc(iz,ix)) 
	  iwork(ix) = dimag(rc(iz,ix)) 
   end do
   do ix = ntr+1,nx2
	  rwork(ix) = 0.0d0 
	  iwork(ix) = 0.0d0 
   end do
   call fft(rwork,iwork,nx2,nx2,nx2,-1)
   do ix = 1,nx2
      rc(iz,ix) = dcmplx(rwork(ix), iwork(ix)) 
   end do
end do
matcommand = "disp('> Initial FT X -> Kx, OK.')"
istat =  mexEvalString(matcommand(1:len_trim(matcommand)))
matcommand = "disp('> Beginning forward computations.')"
istat =  mexEvalString(matcommand(1:len_trim(matcommand)))
!*******************************
!*     LOOP OVER DEPTHS
!*******************************
istat =  mexEvalString('h=waitbar(0,''Computing... please wait'');')
do iz = nz,1,-1
!**** Report progress  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	progress = float(nz-iz+1)/float(nz)
    write(matcommand,'(a8, f6.4, a3)') 'waitbar(', progress, ',h);'
    istat =  mexEvalString(matcommand(1:len_trim(matcommand)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   nq = (2.0d0/pi)*datan(qm(iz))
   expq = (1.0d0 - nq)/2.0d0
!*******************************
!*    LOOP OVER WAVENUMBERS
!*******************************
   do ikx = 1,nx2
      den = (waven(ikx)*waven(ikx))
      denkc = dcmplx(den,0.0)
!*******************************
!*    LOOP OVER FREQUENCIES
!*******************************
      do iw = iw1,iw2
         ctr = (omega(iw)/wr)**expq
         bomega = omega(iw)/( vmean(iz) * ctr )
         aomega = bomega*dtan((1.0d0 - nq)*(pi/4.0d0))
         komega(iw) = dcmplx(bomega,aomega)
!*******************************
!*  INCLUDE RADIATION PATTERN
!*******************************
         difc = (omega(iw)/vmean(1))**2 - den
         if (difc > 0.0) then
            kz = dsqrt(difc)
            tet = datan(waven(ikx)/kz)
            radpattern = source(tet1,tet2,tet)
          else
            radpattern = 1.0d0
          end if
          denc = komega(iw)*komega(iw)
          zdif = denc - denkc
          denkz = cdsqrt(zdif)
          den1 = dreal(denkz)
          den2 = dimag(denkz)
          if ( den2 > 0.0) then
             rshift = dexp(-1.0d0*dz*den2)
          else
             rshift = dexp(1.0d0*dz*den2)
          end if
          shift1 = cdexp(dcmplx(0.0, dz*den1))
          shift2 = shift1*rshift
          cp(iw,ikx) = cp(iw,ikx)*shift2
          if (band > 0.0) cp(iw,ikx) = cp(iw,ikx)*g(iw)
          cp(iw,ikx) = cp(iw,ikx) + rc(iz,ikx)
          cp(iw,ikx) = cp(iw,ikx)*radpattern
      end do
   end do
!**** Inverse transform Kx -> X
   do iw = iw1,iw2
      do ikx = 1,nx2
  	     rwork(ikx) = dreal(cp(iw,ikx)*(1.0d0/dfloat(nx2)))
	     iwork(ikx) = dimag(cp(iw,ikx)*(1.0d0/dfloat(nx2)))
      end do
      call fft(rwork,iwork,nx2,nx2,nx2,1)
      do ikx = 1,nx2
         cp(iw,ikx) = dcmplx(rwork(ikx), iwork(ikx))
     end do
   end do
!**************************************************
!*    CORRECTION FOR LATERAL VELOVITY VARIATIONS
!**************************************************
   do ix = 1,ntr
      nq = (2.0d0/pi)*datan(qtmp(iz,ix))
      expq = (1.0d0 - nq)/2.0d0
      do iw = iw1,iw2
        ctr = (omega(iw)/wr)**expq
        bomegaxz = omega(iw)/( (v(iz,ix) * ctr)/2.0d0 )
        aomegaxz = bomegaxz*dtan((1.0d0-nq)*(pi/4.0d0))
        komegaxz = dcmplx(bomegaxz,aomegaxz)
        zdif = komega(iw) - komegaxz
        phase = dcmplx(0.0,-1.0d0)*zdif
        shift = cdexp(phase*dz)
        cp(iw,ix) = cp(iw,ix)*shift
      end do
   end do
!**** Transform X -> Kx
   do iw = iw1,iw2
      do ikx = 1,nx2
  	     rwork(ikx) = dreal(cp(iw,ikx))
	     iwork(ikx) = dimag(cp(iw,ikx))
      end do
      call fft(rwork,iwork,nx2,nx2,nx2,-1)
      do ikx = 1,nx2
         cp(iw,ikx) = dcmplx(rwork(ikx), iwork(ikx))
      end do
   end do
end do
!*************************
!* END LOOP OVER DEPTHS
!*************************
matcommand = "disp('> Finished forward computations, OK.')"
istat =  mexEvalString(matcommand(1:len_trim(matcommand)))
!**** inverse transform kx -> x
do iw = iw1,iw2
   do ikx = 1,nx2
      rwork(ikx) = dreal(cp(iw,ikx)*(1.0d0/dfloat(nx2)))
      iwork(ikx) = dimag(cp(iw,ikx)*(1.0d0/dfloat(nx2)))
   end do
   call fft(rwork,iwork,nx2,nx2,nx2,1)
   do ikx = 1,nx2
      cp(iw,ikx) = dcmplx(rwork(ikx), iwork(ikx))
   end do
end do
istat =  mexEvalString('close(h);')
matcommand = "disp('> Modeling: Inverse FT Kx -> X, OK.')"
istat =  mexEvalString(matcommand(1:len_trim(matcommand)))
!**** DONE ...
do ix = 1,ntr
   do iw = 1,nt2
  	  rwork(iw) = 0.0d0
	  iwork(iw) = 0.0d0
   end do
   do iw = iw1,iw2
     rwork(iw) = dreal(cp(iw,ix))
     rwork(nt2-iw+1) = dreal(cp(iw,ix))
     iwork(iw) = dimag(cp(iw,ix))
	 iwork(nt2-iw+1) = -dimag(cp(iw,ix))
   end do
   call fft(rwork,iwork,nt2,nt2,nt2,-1)
   i = (ix-1)*nt +1 
   do it = 1,nt
      j = i + it -1
      model(j) = rwork(it)
   enddo
end do
matcommand = "disp('> Done inverse FT Omega -> T OK : OVER and OUT.')"
istat =  mexEvalString(matcommand(1:len_trim(matcommand)))

!**** Free memory ****
deallocate(cp,rc,komega)
deallocate(qtmp,v,vmean,waven,omega,qm,g,rwork,iwork)
end
!
! --------------------------------------------------
subroutine mean_velocity(vmean, v, ntr, nz)
!**** computes mean velocity
implicit none
integer :: nz,ntr,i,j
real*8  :: vmean(nz), v(nz,ntr)
! - - - begin - - -
  do j = 1,nz
    vmean(j) = 0.0
    do i = 1,ntr
      vmean(j) = vmean(j)+ v(j,i)/2.
    end do
    vmean(j) = vmean(j)/float(ntr)
  end do
end subroutine
! --------------------------------------------------
subroutine mean_q(qm, q, ntr, nz)
!*****  computes mean quality factor
implicit none
integer :: nz,ntr,i,j
real*8  :: qm(nz), q(nz, ntr)
! - - - begin - - -
do j = 1,nz
  qm(j) = 0.0
  do i = 1,ntr
    qm(j) = qm(j)+ q(j,i)
  end do
  qm(j) = qm(j)/float(ntr)
end do
end subroutine
! --------------------------------------------------
function source(tet1,tet2,tet)  result (output_8)
!*****  Computes Antenna radiation pattern
implicit none
real(8) :: output_8                                       
real(8),parameter :: pi = 3.141592653589793d0
integer :: tet1,tet2
real(8) :: tet,src,alf1,alf2
! - - - begin - - -
alf1 = pi*tet1/180.0d0
alf2 = pi*tet2/180.0d0
if (tet <= alf1) then
   output_8 = (1.0)
   return
else
   if (tet <= alf2) then
      src = (0.42d0-0.5d0*dcos((1.0d0+(tet-alf1)/(alf2-alf1))*pi)+  &
             0.08d0*dcos((1.0d0+(tet-alf1)/(alf2-alf1))*2.0d0*pi))
      output_8 = src
      return
   else
      output_8 = 0.0d0
      return
   end if
end if
return
end function
! --------------------------------------------------
function npfft(nmin,nmax) result (output_4)
!*************************************************************************
!Return optimal n between nmin and nmax for prime factor fft.
!*************************************************************************
!Input:
!nmin		lower bound on returned value (see notes below)
!nmax		desired (but not guaranteed) upper bound on returned value
!!Returned:	valid n for prime factor fft
!The returned n will be composed of mutually prime factors from
!the set {2,3,4,5,7,8,9,11,13,16}.  Because n cannot exceed
!720720 = 5*7*9*11*13*16, 720720 is returned if nmin exceeds 720720.
!If nmin does not exceed 720720, then the returned n will not be 
!less than nmin.  The optimal n is chosen to minimize the estimated
!cost of performing the fft, while satisfying the constraint, if
!possible, that n not exceed nmax.
!*************************************************************************
!uthor:   Dave Hale, Colorado School of Mines
!F90 code: Andreas Tzanis, University of Athens
!*************************************************************************
implicit none
! - - - arg types - - -
integer :: output_4                                       
integer :: nmin,nmax
! - - - local declarations - - -
integer :: i,j
integer*4,parameter :: ntab=240
integer*4 :: nn(ntab)
real*4    :: c(ntab)
data nn/  &
     1,     2,     3,     4,     5,     6,     7,     8, &
     9,    10,    11,    12,    13,    14,    15,    16, &
    18,    20,    21,    22,    24,    26,    28,    30, &
    33,    35,    36,    39,    40,    42,    44,    45, &
    48,    52,    55,    56,    60,    63,    65,    66, &
    70,    72,    77,    78,    80,    84,    88,    90, &
    91,    99,   104,   105,   110,   112,   117,   120, &
   126,   130,   132,   140,   143,   144,   154,   156, &
   165,   168,   176,   180,   182,   195,   198,   208, &
   210,   220,   231,   234,   240,   252,   260,   264, &
   273,   280,   286,   308,   312,   315,   330,   336, &
   360,   364,   385,   390,   396,   420,   429,   440, &
   455,   462,   468,   495,   504,   520,   528,   546, &
   560,   572,   585,   616,   624,   630,   660,   693, &
   715,   720,   728,   770,   780,   792,   819,   840, &
   858,   880,   910,   924,   936,   990,  1001,  1008, &
  1040,  1092,  1144,  1155,  1170,  1232,  1260,  1287, &
  1320,  1365,  1386,  1430,  1456,  1540,  1560,  1584, &
  1638,  1680,  1716,  1820,  1848,  1872,  1980,  2002, &
  2145,  2184,  2288,  2310,  2340,  2520,  2574,  2640, &
  2730,  2772,  2860,  3003,  3080,  3120,  3276,  3432, &
  3465,  3640,  3696,  3960,  4004,  4095,  4290,  4368, &
  4620,  4680,  5005,  5040,  5148,  5460,  5544,  5720, &
  6006,  6160,  6435,  6552,  6864,  6930,  7280,  7920, &
  8008,  8190,  8580,  9009,  9240,  9360, 10010, 10296, &
 10920, 11088, 11440, 12012, 12870, 13104, 13860, 15015, &
 16016, 16380, 17160, 18018, 18480, 20020, 20592, 21840, &
 24024, 25740, 27720, 30030, 32760, 34320, 36036, 40040, &
 45045, 48048, 51480, 55440, 60060, 65520, 72072, 80080, &
 90090,102960,120120,144144,180180,240240,360360,720720/
data c/  &
 0.000052,0.000061,0.000030,0.000053,0.000066,0.000067,0.000071,0.000062, &
 0.000079,0.000080,0.000052,0.000069,0.000103,0.000123,0.000050,0.000086, &
 0.000108,0.000101,0.000098,0.000135,0.000090,0.000165,0.000084,0.000132, &
 0.000158,0.000138,0.000147,0.000207,0.000156,0.000158,0.000176,0.000171, &
 0.000185,0.000227,0.000242,0.000194,0.000215,0.000233,0.000288,0.000271, &
 0.000248,0.000247,0.000285,0.000395,0.000285,0.000209,0.000332,0.000321, &
 0.000372,0.000400,0.000391,0.000358,0.000440,0.000367,0.000494,0.000413, &
 0.000424,0.000549,0.000480,0.000450,0.000637,0.000497,0.000590,0.000626, &
 0.000654,0.000536,0.000656,0.000611,0.000730,0.000839,0.000786,0.000835, &
 0.000751,0.000826,0.000926,0.000991,0.000852,0.000820,0.001053,0.000987, &
 0.001152,0.000952,0.001299,0.001155,0.001270,0.001156,0.001397,0.001173, &
 0.001259,0.001471,0.001569,0.001767,0.001552,0.001516,0.002015,0.001748, &
 0.001988,0.001921,0.001956,0.002106,0.001769,0.002196,0.002127,0.002454, &
 0.002099,0.002632,0.002665,0.002397,0.002711,0.002496,0.002812,0.002949, &
 0.003571,0.002783,0.003060,0.003392,0.003553,0.003198,0.003726,0.003234, &
 0.004354,0.003800,0.004304,0.003975,0.004123,0.004517,0.005066,0.003902, &
 0.004785,0.005017,0.005599,0.005380,0.005730,0.005323,0.005112,0.006658, &
 0.005974,0.006781,0.006413,0.007622,0.006679,0.007032,0.007538,0.007126, &
 0.007979,0.007225,0.008961,0.008818,0.008427,0.009004,0.009398,0.010830, &
 0.012010,0.010586,0.012058,0.011673,0.011700,0.011062,0.014313,0.013021, &
 0.014606,0.013216,0.015789,0.016988,0.014911,0.016393,0.016741,0.018821, &
 0.018138,0.018892,0.018634,0.020216,0.022455,0.022523,0.026087,0.023474, &
 0.024590,0.025641,0.030303,0.025253,0.030364,0.031250,0.029412,0.034404, &
 0.037500,0.034091,0.040214,0.037221,0.042735,0.040214,0.042980,0.045872, &
 0.049505,0.049834,0.055762,0.057034,0.054945,0.056818,0.066667,0.065502, &
 0.068182,0.065217,0.075000,0.078534,0.087719,0.081081,0.084270,0.102740, &
 0.106383,0.105634,0.119048,0.123967,0.119048,0.137615,0.140187,0.154639, &
 0.168539,0.180723,0.180723,0.220588,0.241935,0.254237,0.254237,0.288462, &
 0.357143,0.357143,0.384615,0.384615,0.454545,0.517241,0.576923,0.625000, &
 0.833333,0.789474,1.153846,1.153846,1.875000,2.500000,3.750000,7.5000000/
! - - - begin - - -
i = 1
do while (i <= ntab .AND. nn(i) <= nmin)
   i = i+1
   j = i
   do while (j<=ntab .AND. nn(j) <= nmax)
      j = j+1
	  IF (c(j) < c(i) ) i = j
   enddo
enddo
output_4 = nn(i)
return
end function
! --------------------------------------------------
      subroutine fft(a,b,ntot,n,nspan,isn)
!  multivariate complex fourier transform, computed in place
!    using mixed-radix fast fourier transform algorithm.
!  by r. c. singleton, stanford research institute, sept. 1968
!  arrays a and b originally hold the real and imaginary
!    components of the data, and return the real and
!    imaginary components of the resulting fourier coefficients.
!  multivariate data is indexed according to the fortran
!    array element successor function, without limit
!    on the number of implied multiple subscripts.
!    the subroutine is called once for each variate.
!    the calls for a multivariate transform may be in any order.
!  ntot is the total number of complex data values.
!  n is the dimension of the current variable.
!  nspan/n is the spacing of consecutive data values
!    while indexing the current variable.
!  the sign of isn determines the sign of the complex
!    exponential, and the magnitude of isn is normally one.
!  a tri-variate transform with a(n1,n2,n3), b(n1,n2,n3)
!    is computed by
!      call fft(a,b,n1*n2*n3,n1,n1,1)
!      call fft(a,b,n1*n2*n3,n2,n1*n2,1)
!      call fft(a,b,n1*n2*n3,n3,n1*n2*n3,1)
!  for a single-variate transform,
!    ntot = n = nspan = (number of complex data values), e.g.
!      call fft(a,b,n,n,n,1)
!  the data can alternatively be stored in a single complex array c
!    in standard fortran fashion, i.e. alternating real and imaginary
!    parts. then with most fortran compilers, the complex array c can
!    be equivalenced to a real array a, the magnitude of isn changed
!    to two to give correct indexing increment, and a(1) and a(2) used
!    to pass the initial addresses for the sequences of real and
!    imaginary values, e.g.
!       complex c(ntot)
!       real    a(2*ntot)
!       equivalence (c(1),a(1))
!       call fft(a(1),a(2),ntot,n,nspan,2)
!  arrays at(maxf), ck(maxf), bt(maxf), sk(maxf), and np(maxp)
!    are used for temporary storage.  if the available storage
!    is insufficient, the program is terminated by a stop.
!    maxf must be .ge. the maximum prime factor of n.
!    maxp must be .gt. the number of prime factors of n.
!    in addition, if the square-free portion k of n has two or
!    more prime factors, then maxp must be .ge. k-1.
      real*8 a(1),b(1)
!  array storage in nfac for a maximum of 15 prime factors of n.
!  if n has more than one square-free factor, the product of the
!    square-free factors must be .le. 210
      dimension nfac(11),np(209)
!  array storage for maximum prime factor of 23
      dimension at(23),ck(23),bt(23),sk(23)
      equivalence (i,ii)
!  the following two constants should agree with the array dimensions.
      maxp=209
!
! Date: Wed, 9 Aug 1995 09:38:49 -0400
! From: ldm@apollo.numis.nwu.edu
      maxf=23
!
      if(n .lt. 2) return
      inc=isn
      c72=0.30901699437494742
      s72=0.95105651629515357
      s120=0.86602540378443865
      rad=6.2831853071796
      if(isn .ge. 0) go to 10
      s72=-s72
      s120=-s120
      rad=-rad
      inc=-inc
   10 nt=inc*ntot
      ks=inc*nspan
      kspan=ks
      nn=nt-inc
      jc=ks/n
      radf=rad*float(jc)*0.5
      i=0
      jf=0
!  determine the factors of n
      m=0
      k=n
      go to 20
   15 m=m+1
      nfac(m)=4
      k=k/16
   20 if(k-(k/16)*16 .eq. 0) go to 15
      j=3
      jj=9
      go to 30
   25 m=m+1
      nfac(m)=j
      k=k/jj
   30 if(mod(k,jj) .eq. 0) go to 25
      j=j+2
      jj=j**2
      if(jj .le. k) go to 30
      if(k .gt. 4) go to 40
      kt=m
      nfac(m+1)=k
      if(k .ne. 1) m=m+1
      go to 80
   40 if(k-(k/4)*4 .ne. 0) go to 50
      m=m+1
      nfac(m)=2
      k=k/4
   50 kt=m
      j=2
   60 if(mod(k,j) .ne. 0) go to 70
      m=m+1
      nfac(m)=j
      k=k/j
   70 j=((j+1)/2)*2+1
      if(j .le. k) go to 60
   80 if(kt .eq. 0) go to 100
      j=kt
   90 m=m+1
      nfac(m)=nfac(j)
      j=j-1
      if(j .ne. 0) go to 90
!  compute fourier transform
  100 sd=radf/float(kspan)
      cd=2.0*sin(sd)**2
      sd=sin(sd+sd)
      kk=1
      i=i+1
      if(nfac(i) .ne. 2) go to 400
!  transform for factor of 2 (including rotation factor)
      kspan=kspan/2
      k1=kspan+2
  210 k2=kk+kspan
      ak=a(k2)
      bk=b(k2)
      a(k2)=a(kk)-ak
      b(k2)=b(kk)-bk
      a(kk)=a(kk)+ak
      b(kk)=b(kk)+bk
      kk=k2+kspan
      if(kk .le. nn) go to 210
      kk=kk-nn
      if(kk .le. jc) go to 210
      if(kk .gt. kspan) go to 800
  220 c1=1.0-cd
      s1=sd
  230 k2=kk+kspan
      ak=a(kk)-a(k2)
      bk=b(kk)-b(k2)
      a(kk)=a(kk)+a(k2)
      b(kk)=b(kk)+b(k2)
      a(k2)=c1*ak-s1*bk
      b(k2)=s1*ak+c1*bk
      kk=k2+kspan
      if(kk .lt. nt) go to 230
      k2=kk-nt
      c1=-c1
      kk=k1-k2
      if(kk .gt. k2) go to 230
      ak=c1-(cd*c1+sd*s1)
      s1=(sd*c1-cd*s1)+s1
      c1=2.0-(ak**2+s1**2)
      s1=c1*s1
      c1=c1*ak
      kk=kk+jc
      if(kk .lt. k2) go to 230
      k1=k1+inc+inc
      kk=(k1-kspan)/2+jc
      if(kk .le. jc+jc) go to 220
      go to 100
!  transform for factor of 3 (optional code)
  320 k1=kk+kspan
      k2=k1+kspan
      ak=a(kk)
      bk=b(kk)
      aj=a(k1)+a(k2)
      bj=b(k1)+b(k2)
      a(kk)=ak+aj
      b(kk)=bk+bj
      ak=-0.5*aj+ak
      bk=-0.5*bj+bk
      aj=(a(k1)-a(k2))*s120
      bj=(b(k1)-b(k2))*s120
      a(k1)=ak-bj
      b(k1)=bk+aj
      a(k2)=ak+bj
      b(k2)=bk-aj
      kk=k2+kspan
      if(kk .lt. nn) go to 320
      kk=kk-nn
      if(kk .le. kspan) go to 320
      go to 700
!  transform for factor of 4
  400 if(nfac(i) .ne. 4) go to 600
      kspnn=kspan
      kspan=kspan/4
  410 c1=1.0
      s1=0
  420 k1=kk+kspan
      k2=k1+kspan
      k3=k2+kspan
      akp=a(kk)+a(k2)
      akm=a(kk)-a(k2)
      ajp=a(k1)+a(k3)
      ajm=a(k1)-a(k3)
      a(kk)=akp+ajp
      ajp=akp-ajp
      bkp=b(kk)+b(k2)
      bkm=b(kk)-b(k2)
      bjp=b(k1)+b(k3)
      bjm=b(k1)-b(k3)
      b(kk)=bkp+bjp
      bjp=bkp-bjp
      if(isn .lt. 0) go to 450
      akp=akm-bjm
      akm=akm+bjm
      bkp=bkm+ajm
      bkm=bkm-ajm
      if(s1 .eq. 0) go to 460
  430 a(k1)=akp*c1-bkp*s1
      b(k1)=akp*s1+bkp*c1
      a(k2)=ajp*c2-bjp*s2
      b(k2)=ajp*s2+bjp*c2
      a(k3)=akm*c3-bkm*s3
      b(k3)=akm*s3+bkm*c3
      kk=k3+kspan
      if(kk .le. nt) go to 420
  440 c2=c1-(cd*c1+sd*s1)
      s1=(sd*c1-cd*s1)+s1
      c1=2.0-(c2**2+s1**2)
      s1=c1*s1
      c1=c1*c2
      c2=c1**2-s1**2
      s2=2.0*c1*s1
      c3=c2*c1-s2*s1
      s3=c2*s1+s2*c1
      kk=kk-nt+jc
      if(kk .le. kspan) go to 420
      kk=kk-kspan+inc
      if(kk .le. jc) go to 410
      if(kspan .eq. jc) go to 800
      go to 100
  450 akp=akm+bjm
      akm=akm-bjm
      bkp=bkm-ajm
      bkm=bkm+ajm
      if(s1 .ne. 0) go to 430
  460 a(k1)=akp
      b(k1)=bkp
      a(k2)=ajp
      b(k2)=bjp
      a(k3)=akm
      b(k3)=bkm
      kk=k3+kspan
      if(kk .le. nt) go to 420
      go to 440
!  transform for factor of 5 (optional code)
  510 c2=c72**2-s72**2
      s2=2.0*c72*s72
  520 k1=kk+kspan
      k2=k1+kspan
      k3=k2+kspan
      k4=k3+kspan
      akp=a(k1)+a(k4)
      akm=a(k1)-a(k4)
      bkp=b(k1)+b(k4)
      bkm=b(k1)-b(k4)
      ajp=a(k2)+a(k3)
      ajm=a(k2)-a(k3)
      bjp=b(k2)+b(k3)
      bjm=b(k2)-b(k3)
      aa=a(kk)
      bb=b(kk)
      a(kk)=aa+akp+ajp
      b(kk)=bb+bkp+bjp
      ak=akp*c72+ajp*c2+aa
      bk=bkp*c72+bjp*c2+bb
      aj=akm*s72+ajm*s2
      bj=bkm*s72+bjm*s2
      a(k1)=ak-bj
      a(k4)=ak+bj
      b(k1)=bk+aj
      b(k4)=bk-aj
      ak=akp*c2+ajp*c72+aa
      bk=bkp*c2+bjp*c72+bb
      aj=akm*s2-ajm*s72
      bj=bkm*s2-bjm*s72
      a(k2)=ak-bj
      a(k3)=ak+bj
      b(k2)=bk+aj
      b(k3)=bk-aj
      kk=k4+kspan
      if(kk .lt. nn) go to 520
      kk=kk-nn
      if(kk .le. kspan) go to 520
      go to 700
!  transform for odd factors
  600 k=nfac(i)
      kspnn=kspan
      kspan=kspan/k
      if(k .eq. 3) go to 320
      if(k .eq. 5) go to 510
      if(k .eq. jf) go to 640
      jf=k
      s1=rad/float(k)
      c1=cos(s1)
      s1=sin(s1)
      if(jf .gt. maxf) go to 998
      ck(jf)=1.0
      sk(jf)=0.0
      j=1
  630 ck(j)=ck(k)*c1+sk(k)*s1
      sk(j)=ck(k)*s1-sk(k)*c1
      k=k-1
      ck(k)=ck(j)
      sk(k)=-sk(j)
      j=j+1
      if(j .lt. k) go to 630
  640 k1=kk
      k2=kk+kspnn
      aa=a(kk)
      bb=b(kk)
      ak=aa
      bk=bb
      j=1
      k1=k1+kspan
  650 k2=k2-kspan
      j=j+1
      at(j)=a(k1)+a(k2)
      ak=at(j)+ak
      bt(j)=b(k1)+b(k2)
      bk=bt(j)+bk
      j=j+1
      at(j)=a(k1)-a(k2)
      bt(j)=b(k1)-b(k2)
      k1=k1+kspan
      if(k1 .lt. k2) go to 650
      a(kk)=ak
      b(kk)=bk
      k1=kk
      k2=kk+kspnn
      j=1
  660 k1=k1+kspan
      k2=k2-kspan
      jj=j
      ak=aa
      bk=bb
      aj=0.0
      bj=0.0
      k=1
  670 k=k+1
      ak=at(k)*ck(jj)+ak
      bk=bt(k)*ck(jj)+bk
      k=k+1
      aj=at(k)*sk(jj)+aj
      bj=bt(k)*sk(jj)+bj
      jj=jj+j
      if(jj .gt. jf) jj=jj-jf
      if(k .lt. jf) go to 670
      k=jf-j
      a(k1)=ak-bj
      b(k1)=bk+aj
      a(k2)=ak+bj
      b(k2)=bk-aj
      j=j+1
      if(j .lt. k) go to 660
      kk=kk+kspnn
      if(kk .le. nn) go to 640
      kk=kk-nn
      if(kk .le. kspan) go to 640
!  multiply by rotation factor (except for factors of 2 and 4)
  700 if(i .eq. m) go to 800
      kk=jc+1
  710 c2=1.0-cd
      s1=sd
  720 c1=c2
      s2=s1
      kk=kk+kspan
  730 ak=a(kk)
      a(kk)=c2*ak-s2*b(kk)
      b(kk)=s2*ak+c2*b(kk)
      kk=kk+kspnn
      if(kk .le. nt) go to 730
      ak=s1*s2
      s2=s1*c2+c1*s2
      c2=c1*c2-ak
      kk=kk-nt+kspan
      if(kk .le. kspnn) go to 730
      c2=c1-(cd*c1+sd*s1)
      s1=s1+(sd*c1-cd*s1)
      c1=2.0-(c2**2+s1**2)
      s1=c1*s1
      c2=c1*c2
      kk=kk-kspnn+jc
      if(kk .le. kspan) go to 720
      kk=kk-kspan+jc+inc
      if(kk .le. jc+jc) go to 710
      go to 100
!  permute the results to normal order---done in two stages
!  permutation for square factors of n
  800 np(1)=ks
      if(kt .eq. 0) go to 890
      k=kt+kt+1
      if(m .lt. k) k=k-1
      j=1
      np(k+1)=jc
  810 np(j+1)=np(j)/nfac(j)
      np(k)=np(k+1)*nfac(j)
      j=j+1
      k=k-1
      if(j .lt. k) go to 810
      k3=np(k+1)
      kspan=np(2)
      kk=jc+1
      k2=kspan+1
      j=1
      if(n .ne. ntot) go to 850
!  permutation for single-variate transform (optional code)
  820 ak=a(kk)
      a(kk)=a(k2)
      a(k2)=ak
      bk=b(kk)
      b(kk)=b(k2)
      b(k2)=bk
      kk=kk+inc
      k2=kspan+k2
      if(k2 .lt. ks) go to 820
  830 k2=k2-np(j)
      j=j+1
      k2=np(j+1)+k2
      if(k2 .gt. np(j)) go to 830
      j=1
  840 if(kk .lt. k2) go to 820
      kk=kk+inc
      k2=kspan+k2
      if(k2 .lt. ks) go to 840
      if(kk .lt. ks) go to 830
      jc=k3
      go to 890
!  permutation for multivariate transform
  850 k=kk+jc
  860 ak=a(kk)
      a(kk)=a(k2)
      a(k2)=ak
      bk=b(kk)
      b(kk)=b(k2)
      b(k2)=bk
      kk=kk+inc
      k2=k2+inc
      if(kk .lt. k) go to 860
      kk=kk+ks-jc
      k2=k2+ks-jc
      if(kk .lt. nt) go to 850
      k2=k2-nt+kspan
      kk=kk-nt+jc
      if(k2 .lt. ks) go to 850
  870 k2=k2-np(j)
      j=j+1
      k2=np(j+1)+k2
      if(k2 .gt. np(j)) go to 870
      j=1
  880 if(kk .lt. k2) go to 850
      kk=kk+jc
      k2=kspan+k2
      if(k2 .lt. ks) go to 880
      if(kk .lt. ks) go to 870
      jc=k3
  890 if(2*kt+1 .ge. m) return
      kspnn=np(kt+1)
!  permutation for square-free factors of n
      j=m-kt
      nfac(j+1)=1
  900 nfac(j)=nfac(j)*nfac(j+1)
      j=j-1
      if(j .ne. kt) go to 900
      kt=kt+1
      nn=nfac(kt)-1
      if(nn .gt. maxp) go to 998
      jj=0
      j=0
      go to 906
  902 jj=jj-k2
      k2=kk
      k=k+1
      kk=nfac(k)
  904 jj=kk+jj
      if(jj .ge. k2) go to 902
      np(j)=jj
  906 k2=nfac(kt)
      k=kt+1
      kk=nfac(k)
      j=j+1
      if(j .le. nn) go to 904
!  determine the permutation cycles of length greater than 1
      j=0
      go to 914
  910 k=kk
      kk=np(k)
      np(k)=-kk
      if(kk .ne. j) go to 910
      k3=kk
  914 j=j+1
      kk=np(j)
      if(kk .lt. 0) go to 914
      if(kk .ne. j) go to 910
      np(j)=-j
      if(j .ne. nn) go to 914
      maxf=inc*maxf
!  reorder a and b, following the permutation cycles
      go to 950
  924 j=j-1
      if(np(j) .lt. 0) go to 924
      jj=jc
  926 kspan=jj
      if(jj .gt. maxf) kspan=maxf
      jj=jj-kspan
      k=np(j)
      kk=jc*k+ii+jj
      k1=kk+kspan
      k2=0
  928 k2=k2+1
      at(k2)=a(k1)
      bt(k2)=b(k1)
      k1=k1-inc
      if(k1 .ne. kk) go to 928
  932 k1=kk+kspan
      k2=k1-jc*(k+np(k))
      k=-np(k)
  936 a(k1)=a(k2)
      b(k1)=b(k2)
      k1=k1-inc
      k2=k2-inc
      if(k1 .ne. kk) go to 936
      kk=k2
      if(k .ne. j) go to 932
      k1=kk+kspan
      k2=0
  940 k2=k2+1
      a(k1)=at(k2)
      b(k1)=bt(k2)
      k1=k1-inc
      if(k1 .ne. kk) go to 940
      if(jj .ne. 0) go to 926
      if(j .ne. 1) go to 924
  950 j=k3+1
      nt=nt-kspnn
      ii=nt-inc+1
      if(nt .ge. 0) go to 924
      return
!  error finish, insufficient array storage
  998 isn=0
      istat=mexEvalString('ERROR : array bounds exceeded within subroutine fft')
      stop
      end
