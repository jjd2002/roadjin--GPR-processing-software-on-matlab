function [dsbg, ww] = rmswbackgr( d )
%
% RMSWBACKGR :  Removes a sliding-window background trace to reduce
%               subhorizontal features
%
%      Usage :  [dsbg, ww] = rmswbackgr( d )
%
%      Input :  d is the input GPR data matrix 
%
%    Outputs :
%       dsbg : The filtered data matrix
%         ww : The width of the sliding window
%
%   Author   : Andreas Tzanis,
%              Dept. of Geophysics, 
%              University of Athens
%              atzanis@geol.uoa.gr
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

[Samples_per_scan, No_traces] = size(d);
%%%%%  Inquire window size
%disp([ 'No of traces in data is ' num2str(No_traces) '.']);
%ww = input('Give window width in # traces : ');
answer = inputdlg('以道数形式给出时窗宽度',' ',1);
if isempty(answer), 
    dsbg = []; ww=[];
    return; 
end;
lb = char(answer{1});
    comma = findstr(lb,',');
    if ~isempty(comma),  
        lb(comma) = '.'; 
    end
    ww = floor(str2num(lb));    ww=ww(1);

% Set up surrogate data
dsbg = d;

hw = waitbar(0,'Filtering in Progress');
 %%% Do the first "ww/2" traces
    ii          = 1 : floor(ww/2);
    bg          = mean(d(:,ii)')';
for i = 1 : floor(ww/2)
    dsbg(:,i)  = dsbg(:,i) - bg;
    waitbar(i/No_traces);
end
%%% Do the middle (full) range (ww/2 - (No_traces - ww/2)) 
for i = floor(ww/2)+1 : No_traces - floor(ww/2)
    ii         = i - floor(ww/2) : i + floor(ww/2);
    bg         = mean(d(:,ii)')'; 
    dsbg(:,i) = dsbg(:,i) - bg;
    waitbar(i/No_traces);
end
 %%% Do the last "ww/2" traces
    ii          = No_traces - floor(ww/2) + 1 : No_traces;
    bg          = mean(d(:,ii)')';
for i = No_traces - floor(ww/2) + 1 : No_traces
    dsbg(:,i)  = dsbg(:,i) - bg;
    waitbar(i/No_traces);
end
close(hw);
           