function [dsfg, ww] = rmswforegr( d )
%
% RMSWFOREGR : Removes a sliding-window foreground trace from the 2-D
%              GPR section passed in d, to reduce dipping features.
%
%      Usage : [dsfg, ww] = rmswforegr( d )
%
%      Input :  d is the 2-D GPR data matrix 
%
%    Outputs :
%       dsfg : The filtered data matrix
%         ww : The width of the sliding window
%
%   Author   :  Andreas Tzanis,
%               Dept. of Geophysics, 
%               University of Athens
%               atzanis@geol.uoa.gr
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
%% 

[Samples_per_scan, No_traces] = size(d);
%%%%%  Inquire window size
answer = inputdlg('以道数形式给出时窗宽度',' ',1);
if isempty(answer), 
    dsfg = []; ww=[]
    return; 
end;
lb = char(answer{1});
comma = findstr(lb,',');
if ~isempty(comma),  
    lb(comma) = '.'; 
end
ww = floor(str2num(lb));    ww=ww(1);
dsfg        = d;
hw = waitbar(0,'Filtering in Progress');
 %%% Do the first "ww/2" traces
ii           = 1 : floor(ww/2);
bg           = mean(d(:,ii)')';
dsfg(:,ii)  = bg*ones(1,length(ii));
%%% Do the middle (ww/2 - (No_traces - ww/2)) traces
for i = floor(ww/2)+1 : No_traces - floor(ww/2)
    ii          = i - floor(ww/2) : i + floor(ww/2);
    bg          = mean(d(:,ii)')'; 
    dsfg(:,i)  = bg;
    waitbar(i/No_traces);
end
 %%% Do the last "ww/2" traces
ii          = No_traces - floor(ww/2) + 1 : No_traces;
bg          = mean(d(:,ii)')';
fg          =  - dsfg(:,ii);
dsfg(:,ii) = bg*ones(1,length(ii));
close(hw);
return           