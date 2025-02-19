function [dmig,zmig,nz] = splitstepmig(d, dt, dx, v2d, zv, dz)
%
% SPLITSTEPMIG - Split-step Fourier depth migration.
%              * Depending on the user's input, this routine may, account
%                for the frequency dependence of the phase velocity
%                (dispersion). This is done with the approach of Bano
%                (1996, Geophys. J. Int., 279 - 288), whereby the phase
%                velocity depends on a non-dispersive term Vo and an
%                exponential correction as:
%                V(w) = Vo*(w / wc)^p,  p = (1 - n)/2,  n = (2/pi)*atan(Q),
%                where w is the ang. frequency, wc is the antenna central
%                ang. frequency, and Q is the quality factor. 
%              * Both the layer reference velocity and layer perturbation
%                velocities are reduced for frequency dependence in the
%                F-domain.
%             => To account for dispersion effects, introduce the velocity
%                structure "v2d" as a {1 x 3} cell array, with sub-array
%                elements {Vo Q antenna_freqquency} in this order. 
%             => If the velocity structure "v2d" is introduced as a cell
%                array with fewer than 3 elements, or as a single matrix,
%                the routine will use only the non-dispersive term Vo. 
%
%   Usage : [dmig, zmig, nz] = splitstepmig(d, dt, dx, v2d, zv, dz)
%
%  Inputs :
%       d : The [ns x ntr] array of the equally-spaced GPR section
%      dt : The sampling rate in the ns dimension (in nanosecs)
%      dx : Trace spacing in the ntr dimension (in m).
%     v2d : The velocity structure.
%           ** If it is a {1 x 3} cell array it must be structured as: 
%              1. v2d{1) : [ns x ntr] array, the non-dispersive term of the
%                 velocity model.  
%              2. v2d{2} : [ns x ntr] array, the quality factor of the 2-D
%                 model. 
%              3. v2d{3) : Scalar, float, antenna central frequency.
%           ** If v2d is a {1 x 1} or a {1 x 2} cell array, then v2d{1)
%              must be the [ns x ntr] non-dispersive term of the velocity
%              model. 
%           ** If v2d is a single matrix, then it must be the [ns x ntr]
%              non-dispersive term of the velocity model. 
%      zv : [ns x 1] vector, the z-coordinates of the velocity model.
%      dz : Depth spacing used for integration (stepping size).
%
% Outputs : 
%    dmig : [nz x ntr] array, the depth migrated GPR section.
%    zmig : [nz x 1] vector, the depth coordinates of the migrated section.
%      nz : The vertical dimension of dmig == dimension of zmig
%
%    Uses : Attached subfunction "tofk.m"   : Performs F-K transform
%           Attached subfunction "cosbel.m" : Cosine-tapered boxcar window      
%
%Requires : checkcomma.m
%
% Credits : Some programming tips were taken from the split-step migration
%           routine by G.F. Margrave, CREWES Project, U. of Calgary, 1996. 
% 
%  Author : Andreas Tzanis, 
%           Department of Geophysics, 
%           University of Athens, 
%           atzanis@geol.uoa.gr
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%

%%% Empty data array
if isempty(d),
    info = {'The input data array is empty! Aborting!'}; 
    erh = errordlg( info,'SPLITSTEPMIG: ERROR');
    uiwait(erh)
    dmig = [];        zmig = [];    nz=[];
    return
end
%%% Velocity model and its depth coordinates not supplied
if nargin < 5 || nargin < 4,
    info = {'Velocity information not supplied or incomplete.' ...
          'Please generate a 2D velocity model using GET2DVELOCITY'}; 
    erh = errordlg( info,'SPLITSTEPMIG: ERROR');
    uiwait(erh)
    dmig = [];        zmig = [];    nz=[];
    return
end
if nargin < 6,           % Depth stepping size may not have been furnished
    dz = 0.02; 
    cldef=cellstr(num2str(dz));
    answer = inputdlg('Give depth spacing in m', ...
        'SPLITSTEPMIG : Requesting dz',1 , cldef); 
    if isempty(answer)
        dmig = [];   zmig = [];    nz = [];
        return      
    end
    lb = checkcomma(answer);     
    dz = str2num(lb(1,:));
end

% Extract the velocity model and quality factor model from input cell array.
Dispersion = 1;
if ~iscell(v2d),
    Dispersion = 0;
    velmod = v2d;
elseif size(v2d,2) < 3,
    Dispersion = 0;
    velmod = v2d{1};
else
    velmod = v2d{1};
    Q      = v2d{2};
    antenna_freq = v2d{3};
    wc = 2*pi*antenna_freq*1e6;     % antenna central angular frequency
end
clear v2d;

[ns,ntr] = size(d);
% Trap an obvious error ...
ncv = size(velmod,2);
if ntr ~= ncv, 
    erh = errordlg('Mismatch between the sizes of data and velocity model',...
        'SPLITSTEPMIG: ERROR');
    uiwait(erh)
    dmig = [];        zmig = [];    nz = [];
    return
end

% Get zmax and fmax interactively
zmax = zv(ns);                            
fmax = 0.6*(1/(2*dt));  
fnyquist = 1/(2*dt);
clb = cell(2,1);                                        
clb(1)   = cellstr('Give maximum depth to step to in m');
clb(2)   = cellstr(['Give maximum frequency to start from. '...
                    'Default=0.6*Nyquist,  ', ...
                    ['Nyquist=' num2str(fnyquist,4) 'GHz']]);
cldef(1) = cellstr(num2str(zmax));
cldef(2) = cellstr(num2str(fmax));
answer = inputdlg(clb,'SPLITSTEPMIG : Give ZMAX and FMAX',1 , cldef); 
if isempty(answer)
    dmig = [];    zmig = [];    nz=[];
    return      
end
lb = checkcomma(answer);    % ensure that decimals will be correct even if 
zmax = str2num(lb(1,:));    % ',' typed instead of '.'
fmax = str2num(lb(2,:)); 

% Begin
tic                                % time the process
z         = 0:dz:zmax;
nz        = length(z);
%forward F-K tranform
hwbar      = waitbar(0,'Transforming TX - FK_x ...');
[fkdata,f,kx] = tofk(d,dt,dx,10);
if isempty( fkdata ),
    dmig = [];
    zmig = [];
    return
end
kx2     = kx.^2;              % kx spectrum will be wrapped
df      = f(2)-f(1); 
nfmax   = round(fmax/df)+1;
f1      = f(2:nfmax);
f2      = f1.^2;
fkdata  = fkdata(2:nfmax,:);  % don't bother with dc or f > fmax
nfmax   = nfmax-1; 

dmig    = zeros(nz,ntr);
hwbar   = waitbar(0,hwbar, ['Imaging FK_x -> ZX in ' int2str(nz) ...
        ' depth steps']);
for iz = 1:nz
% get velocity at iz'th depth
    izv = find(zv <= z(iz));
    izv = izv(length(izv));
    if Dispersion,
% frequency dependence of phase velocity after Bano, 1996.
        mq = mean(Q(izv,:));     % mean quality factor at iz'th depth
        expq = 0.5*(1 - (2/pi)*atan(mq));  
        disperse = (2*pi*f1*1e9/wc).^expq; 
    else
        disperse = ones(size(f1));
    end
	v   = 0.5*velmod(izv,:);     % mean (reference) velocity at iz'th depth 
    vm  = mean(v);
    vmd = vm * disperse;         % frequency dependence of ref. velocity
	for jk = 1:length(kx)
		fev  = kx(jk)*vm;        % first non-evanescent frequency (approx.)
		nfev = max([round(fev/df),1]);
%        nfev = 1;
		if ( nfev <= nfmax )
            % compute phase shift
            vmd2 = vmd(nfev:nfmax).*vmd(nfev:nfmax);
			phase = (2*pi*f1(nfev:nfmax)./vmd(nfev:nfmax)).*(sqrt(1 - ...
					vmd2.*kx2(jk)./f2(nfev:nfmax) )-1);
            %apply phase shift
			fkdata((nfev:nfmax),jk)=fkdata((nfev:nfmax),jk).*exp(i*phase*dz);
			if ( nfev > 1 )
				fkdata(1:nfev-1,jk) = 0;
			end
		else
			fkdata(:,jk) = 0;
		end
	end
% inverse fft over Kx and apply static phase shift
% also re-zero the zero pad, also image
    if Dispersion,
        nq = (2/pi)*atan(Q(izv,:));
        expq = 0.5*(1 - nq);
    end
	for jf = 1:nfmax
		tmp        = fft(fkdata(jf,1:ntr));
% frequency dependence of phase velocity along scan axis        
        if Dispersion,
            disperse = (2*pi*f1(jf)*1e9/wc).^expq;
        else
            disperse = ones(1,ntr);
        end
        vd  = disperse.*v;
   		tmp      = tmp.*exp(i*2*pi*f1(jf)*dz./vd);
		dmig(iz,:) = dmig(iz,:)+2*real(tmp);
		fkdata(jf,:)  = [ifft(tmp) zeros(1,length(kx)- ntr)];
	end
    waitbar(iz/nz, hwbar);
end
dmig  = dmig/(2*(nfmax+1));
zmig  = z;
time  = toc;
disp(['SPLITSTEPMIG finished in ' num2str(time) ' seconds'])
close(hwbar)
return
%
function [fks, f, kx] = tofk(d, dt, dx, percent)
%
%   TOFK :   Performs a 2-D F-K transform on a real valued data matrix.
%            Only the positive frequencies are kept, while all the
%            wavenumbers are returned.  
%
%  USAGE : [fks, f, kx] = tofk(d, dt, dx, ntpad, nxpad, percent)
%
% Inputs :  
%      d : The [ns x ntr] array of equally-spaced GPR section
%     dt : The sampling rate in the ns dimension (in nanosecs).
%     dx : Trace spacing in the ntr dimension (in m).
% percent: Apply a cosine taper to both t and x dimensions. Length of taper
%          is this percentage of the length of x and t axes. Default = 0. 
%
%Outputs : 
%    fks : The complex valued F-K transform of d.
%      f : The vector of frequency coordinates for the rows of fks
%     kx : The vector of wavenumber coordinates for the columns of fks
% 
% Author:  Andreas Tzanis, 
%          Department of Geophysics,
%          University of Athens, 
%          atzanis@geol.uoa.gr
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%
[ns,ntr] = size(d);
if ( nargin < 4 ) 
    percent = 0.; 
end
%nxpad = 2^nextpow2(ntr);         % pad rows with zeros until this size
nxpad = ntr;
%ntpad = 2^nextpow2(ns);          % pad columns with zeros until this size
ntpad = 2*ns;
if ( percent > 0 )
	 cw = cosbell(ns, percent);
	 cw = cw * ones(1,ntr);
	 d  = d.*cw;
     clear cw
end
% pad time axis if needed
if ( ns < ntpad ),
	d   = [d;  zeros(ntpad-ns, ntr)];
 	ns  = ntpad; 
end
% Transform to frequency domain and keep only the positive frequencies
fxs  = fft(d);
fxs  = fxs(1:round(ns/2+1), :);
fnyq = 1./( 2 * dt );
nf   = size(fxs,1);
f    = linspace(0.,fnyq,nf)';
% ok taper and pad in x
if ( percent > 0 )
    cw  = cosbell(ntr,percent)';
	cw  = ones(nf,1) * cw;
	fxs = fxs.*cw;
	clear cw;
end
if( ntr < nxpad )
	fxs = [fxs zeros(nf, nxpad-ntr)];
	ntr =  nxpad;
end
%fft on rows
fks   = ifft(fxs.',ntr).';
% compute kx to be used in calling function
kxnyq = 1./(2. * dx );
dkx   = 2.*kxnyq/ntr;
if rem(ntr,2),
    kx    = [0:dkx:kxnyq-dkx -kxnyq:dkx:0]; 
else
    kx    = [0:dkx:kxnyq-dkx -kxnyq:dkx:-dkx];
end
return
%
function w = cosbell(n,percent)
%
% COSBEL : Returns an N-point boxcar window tapered at both ends with half
%          cosine bells. The flat region extends over the
%          (100-2*percent)*n/100 central samples, with each ascending and
%          descending half cosine bell tapering the first and last percent
%          points.  
%
%  Usage : w = cosbell(n,percent)
%          w = cosbell(n)
% 
% Inputs  :
%      n  : The length of the window. 
% percent : Percent taper on the ends of the window.
%
% Output  : w -> The length N cosine tapered window
%
%  Author   : Andreas Tzanis,
%             Dept. of Geophysics, 
%             University of Athens
%             atzanis@geol.uoa.gr
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

if ( percent > 50 ) || ( percent < 0 )
   error(' invalid percent for COSBELL')
end
% set defaults
if nargin < 2
    percent = 10;
end
if length(n) > 1
    n = length(n);
end
m  = 2.*percent*n/100.;
m  = 2*floor(m/2);
m2 = m/2;
del= 3.1415927/(m/2);
w  = [ 0.5*(1.0+cos((m2-[1:1:m2])*del))'; ...
            ones(n-m,1); ...
            0.5*(1.0+cos(([(n-m2+1):1:n]-(n-m2+1))*del))'];
return
