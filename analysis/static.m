function dstcorr = static(d,dt,xyz_Tx,xyz_Rx,sdel,wv,swv,direction)   
%
% STATIC  : Apply topographic corrections to single or zero-offset GPR data
%
% Usage   :  dstcorr = static(d,dt,xyz_Tx,xyz_Rx,sdel,wv,swv,direction)
%
%  Inputs :
%       d : The input [ns x ntr] GPR section        
%      dt : Sampling interval
%  xyz_Tx : [ntr x 3] array with the X, Y and Z coordinates of the source
%           antenna in a local frame of reference. For monostatic GPR
%           systems (zero-offset data) xyz_Tx = xyz_Rx
%  xyz_Rx : [ntr x 3] array with the X, Y and Z coordinates of the receiver
%           antenna in a local frame of reference For monostatic GPR 
%           systems (zero-offset data) xyz_Rx = xyz_Tx
%      wv : Weathering layer velocity 
%     swv : Subweathering layer velocity 
%   sdel  : Datum elevation for the source (Tx) antenna. 
%direction: The sense of correction. = -1 shift down in time (default)
%                                    =  1 shift up in time
%
%  Output : 
% dstcorr : [ns, ntr] GPR section reduced for topography.
%
% Author  :  Andreas Tzanis,
%            Department of Geophysics, 
%            University of Athens,
%            atzanis@geol.uoa.gr
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%

% Trap a very-easy-to-do error
if wv <= 0 || wv > 0.2998 || swv <= 0 || swv > 0.2998,
    erh = errordlg('Impossible velocity value! Please try again!',...
            'STATIC: ERROR');
    uiwait(erh);
    return
end
sdepth = 0;             % source depth - assume that it is on surface 
tmin   = 0;             % start correction from time zero
[ns , ntr] = size(d);
dstcorr = zeros(ns,ntr);
traveltime = [tmin:dt:(ns-1)*dt]';
hw = waitbar(0,'Applying Correction ...');
for ix = 1:ntr
    relev = xyz_Rx(ix,3);  
    selev = xyz_Tx(ix,3);  
    tsd = (-selev + sdel + sdepth)/swv;
    tstat = tsd + (selev - relev)/wv;
% Compute output times 
    ts = zeros(ns,1);
    for itime=1:ns
        ts(itime) = tmin + (itime-1)*dt + direction*tstat;
    end
% Interpolate new data / set out-of-range values equal to zero
    tr1 = interp1(traveltime,d(:,ix),ts,'cubic',0);
% Load corrected data array
    dstcorr(:,ix) = tr1;
    waitbar(ix/ntr);
end
close(hw)
return