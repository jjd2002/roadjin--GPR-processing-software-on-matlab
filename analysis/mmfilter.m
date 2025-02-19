function  dmf = mmfilter( d, nv, nx, fmode )
%
%MMFILTER : Applies a 1-D or 2-D spatial smoothing filter by sliding a user
%           defined 2-D window over a 2-D data wall. The data corresponding
%           to the central element of the window is substituted by the mean
%           or median of the window data. 
%       ==> Zero padding equal to 1/2 the window length is applied,
%           therefore the edges of the output data are expected to be
%           distorted. 
%       ==> For meaningful results, the maximum size of the smoothing
%           window should be kept reasonably small, say smaller than 20%
%           the size of the data matrix.
%
%   Usage : dmf  = mmfilter( d, nv, nx, fmode )
%
%  Inputs : 
%       d : The 2-D data matrix (GPR section)
%      nv : The vertical (column) dimension of the filter window 
%      nx : The horizontal (row) dimension of the filter window
%   fmode : "medifilt" : Keyword invoking the median filter mode. 
%           "meanfilt" : Keyword invoking the mean filter mode
%
% Outputs : 
%     dmf : The smoothed (filtered) data matrix
%
%   Author: Andreas Tzanis,
%           Department of Geophysics,
%           University of Athens
%           atzanis@geol.uoa.gr
%
% Copyright (C) 2008, Andreas Tzanis. All rights reserved.
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

%%%%%   Get data size
[ns,ntr]=size(d);
     
%%%%%   Check filter parameters
if nv < 1 || nx < 1,
    erh = errordlg(['Window size cannot be less than 1 element long! ' ...
            'Please try again!'],'MFILTER : ERROR');
    uiwait(erh);
    dmf = []; 
    return
end
if nv > ns/4 || nx > ntr/4,
    erh = errordlg(['Window size cannot be more than 25% of the data!' ...
            ' Please try again!'],'MFILTER : ERROR');
    uiwait(erh);
    dmf = []; 
    return
end

% Now proceed
switch fmode

    %%% Apply MEAN FILTER
    case 'meanfilt'

        h = ones(nv,nx);
        dmf = filter2(h,d,'same');
        dmf = dmf/(nv*nx);

    %%% Apply MEDIAN FILTER
    case 'medifilt'

        %%%%%    Define window size around the center
        if rem(nv,2) ~=0,
            nslow  = ceil(nv/2);
            nshigh = floor(nv/2);
        else
            nslow  = nv(1)/2;
            nshigh = nv/2;
        end
        if rem(nx,2) ~=0,
            ntrlow  = ceil(nx/2);
            ntrhigh = floor(nx/2);
        else
            ntrlow  = nx/2;
            ntrhigh = nx/2;
        end
        %%%%%   Pad data with zeros
        dpad = [ zeros(ns,ntrlow)    d    zeros(ns,ntrhigh)];
        dpad = [ zeros(nslow,ntrlow+ntr+ntrhigh); dpad ;  ...
            zeros(nshigh,ntrlow+ntr+ntrhigh)];

        dmf = dpad;
        tic
        nn = nv*nx;
        ii = 1:1:nn;
        h=waitbar(0,'Median Filtering in progress ... Please wait!');
        for ix = ntrlow:ntr+ntrlow
            for iy = nslow:ns+nslow
                s = dpad(iy-nslow+1:iy+nshigh,ix-ntrlow+1:ix+ntrhigh);
                dmf(iy,ix) = median(s(ii));
            end
            waitbar(ix*iy/(ns*ntr),h)
        end
        dmf = dmf(nslow+1:ns+nslow,ntrlow+1:ntr+ntrlow);   
        close(h)
        time = toc;
        disp(['MEDIAN FILTERING finished in ' num2str(time) ' seconds'])

    %%% Wrong filter mode given
    otherwise

        erh = errordlg('Wrong filter mode - Aborting!','MFILTER : ERROR');
        uiwait(erh);
        dmf = [];
        return

end
return

