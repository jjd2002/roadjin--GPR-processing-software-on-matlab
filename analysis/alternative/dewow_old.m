function dwow = dewow( d );
%
%  DEWOW : Dewows each trace (column) of the 2-D GPR data matrix d by
%          removing the mean and approximately 2 cycles of very low
%          frequency content.  
%
%   Usage : dwow = dewow( d );
%
%   Input :    d,  the 2-D GPR section
%
%  Output : dwow,  the reduced GPR section
%
%Requires : baseline.m operated in the high-pass filter mode.
%
%  Author : Andreas Tzanis,
%           Department of Geophysics, 
%           University of Athens
%           atzanis@geol.uoa.gr
%          (C) 2005, Andreas Tzanis, all rights reserved
%

[ns,ntr] = size(d);
dwow     = zeros(ns,ntr); 
ncycles  = floor(ns/3);
h = waitbar(0,'Processing');
for itr=1:ntr;
    dwow(:,itr) = baseline(ncycles,d(:,itr),'high','cubic');
    dwow(:,itr) = dwow(:,itr) - mean(dwow(:,itr));
%    dwow(:,itr) = d(:,itr) - mean(d(:,itr));
%    dwow(:,itr) = baseline(ncycles,dwow(:,itr),'high','cubic');
    waitbar(itr/ntr,h);
end;
close(h);
return