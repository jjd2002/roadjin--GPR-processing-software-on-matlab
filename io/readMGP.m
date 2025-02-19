function DATA = readMGP(MGPfile)
%
% READMGP : Import the IPD (current input data) structure from a file
%           written with MATGPR's MGP binary format. For details about the
%           MGP format refer to function "save2MGP.m" 
%
% Usage   : DATA = readMGP(MGPfile)
%
% Input   : MPGfile, the complete path and file name of the input MGP file 
%
% Output  : DATA, the MATGPR data structure
%
% Requires: initdatastr.m
%
% Author : Andreas Tzanis
%          Department of Geophysics, 
%          University of Athens
%          atzanis@geol.uoa.gr
%
% Copyright (C) 2008 Andreas Tzanis. All rights reserved.
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

global ENVAR

DATA = initdatastr;
fid=fopen(MGPfile,'r',ENVAR.endian);

% Data origin
ls = fread(fid,1,'int16');
s  = fread(fid,ls,'uchar');
DATA.origin = char(s');
% Original data path 
ls = fread(fid,1,'int16');
s  = fread(fid,ls,'uchar');
DATA.pname = char(s');
% Original data file name
ls = fread(fid,1,'int16');
s  = fread(fid,ls,'uchar');
DATA.fname = char(s');
% Number of samples 
DATA.ns = fread(fid,1,'int16');
% Number of traces
DATA.ntr = fread(fid,1,'int16');
% Data matrix
DATA.d = fread(fid,[DATA.ns DATA.ntr],'float');
% Sampling rate
s = fread(fid,1,'float');
if s ~=0,
    DATA.dt = s;
end
% Make traveltime vector
if ~isempty(DATA.dt),
    DATA.tt2w = 0 : DATA.dt : (DATA.ns-1)*DATA.dt;
end
% Signal position
DATA.sigpos = fread(fid,1,'float');
% First depth estimate (if converted)
d1 = fread(fid,1,'float');
% Depth spacing (if converted)
s = fread(fid,1,'float');
if s ~=0,
    DATA.dz = s;
end
% Make vector of depths 
if ~isempty(DATA.dz),
    DATA.z = d1 : DATA.dz : d1 + (DATA.ns-1)*DATA.dz;
end
% Vertical axis label
ls = fread(fid,1,'int16');
s = fread(fid,ls,'uchar');
DATA.zlab = char(s');
% Fist trace location
t1 = fread(fid,1,'float');
% Trace spacing
s = fread(fid,1,'float');
if s ~= 0;
    DATA.dx = s;
    % Make scan axis vector
    DATA.x = t1 : DATA.dx : t1 + (DATA.ntr-1)*DATA.dx;
else
    DATA.x = 1 : 1 : DATA.ntr;     % data collected at equal time intervals
end
% Scan axis label
ls = fread(fid,1,'int16');
s = fread(fid,ls,'uchar');
DATA.xlab = char(s');
% Marker traces, if any
nm = fread(fid,1,'int16');
mm = fread(fid,1,'int16');
if nm~=0 && mm~=0,
    DATA.markertr = fread(fid,[nm mm],'float');
end
% Transmitter coordinates in global (survey) reference frame
nm = fread(fid,1,'int16');
mm = fread(fid,1,'int16');
if nm~=0 && mm~=0,
    DATA.xyz.Tx = fread(fid,[nm mm],'float');
end
% Receiver coordinates in global (survey) reference frame
nm = fread(fid,1,'int16');
mm = fread(fid,1,'int16');
if nm~=0 && mm~=0,
    DATA.xyz.Rx = fread(fid,[nm mm],'float');
end
% Transmitter - Receiver separation
DATA.TxRx = fread(fid,1,'float');
% Antenna name if written in header files
ls = fread(fid,1,'int16');
s = fread(fid,ls,'uchar');
DATA.Antenna = char(s');
% Header gain for DZT data
nm = fread(fid,1,'int16');
mm = fread(fid,1,'int16');
if nm ~=0 && mm ~=0,
    DATA.DZThdgain = fread(fid,[nm mm],'float');
end
% Number of time saved 
DATA.TimesSaved = fread(fid,1,'int16');
% Header comments
nm = fread(fid,1,'int16');
mm = fread(fid,1,'int16');
if nm~=0 && mm~=0,
    s = fread(fid,[nm mm],'uchar');
end
DATA.comments = char(s); 
% Processing history
nm = fread(fid,1,'int16'); 
DATA.history = cell(nm,1);
for i=1:nm,
    ls = fread(fid,1,'int16');
    s  = fread(fid,ls,'uchar');
    DATA.history(i) = cellstr(char(s'));
end
% Done
fclose(fid);
return