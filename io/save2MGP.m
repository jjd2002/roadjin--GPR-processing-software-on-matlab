function save2MGP(IPD,MGPfile)
%
% SAVE2MGP : Save the IPD (current input data) structure using MATGPR's
%            MGP binary format. The MGP format is as follows:
%--------------------------------------------------------------------------
%     Field             Type       Description 
%  or Variable
%--------------------------------------------------------------------------
%     lo                int16      Length of IPD.origin
%     IPD.origin        uchar      Data origin 
%     lp                int16      Length of IPD.pname
%     IPD.pname         uchar      Data path 
%     lf                int16      Length of IPD.fname
%     IPD.fname         uchar      Data file name   
%     IPD.ns            int16      No of Samples per trace
%     IPD.ntr           int16      Number of traces
%     IPD.d             float      [IPD.ns x IPD.ntr} data matrix
%     IPD.dt            float      Sampling rate
%     IPD.sigpos        float      Signal position
%     d1                float      First depth estimate (if converted)
%     IPD.dz            float      Depth spacing (if converted)
%     lz                int16      Length of IPD.zlab
%     IPD.zlab          uchar      Vertical axis label
%     t1                float      First trace location
%     IPD.dx            float      Trace spacing
%     lx                int16      Length of IPD.xlab
%     IPD.xlab          uchar      Scan axis label
%     mr1               int16      No of rows in IPD.markertr
%     mc1               int16      No of columns in IPD.markertr
%     IPD.markertr      float      Marker traces and coordinates
%     mr2               int16      No of rows in IPD.xyz.Tx
%     mc2               int16      No of columns in IPD.xyz.Tx
%     IPD.xyz.Tx        float      Tx coordinates in survey reference frame
%     mr3               int16      No of rows in IPD.xyz.Rx
%     mc3               int16      No of columns in IPD.xyz.Rx
%     IPD.xyz.Rx        float      Rx coordinates in survey reference frame
%     IPD.TxRx          float      Transmitter - Receiver separation
%     la                int16      Length of IPD.Antenna
%     IPD.Antenna       uchar      Antenna name if in header files
%     mr4               int16      No of rows in IPD.DZThdgain
%     mc4               int16      No of columns in IPD.DZThdgain
%     IPD.DZThdgain     float      Header gain for DZT data
%     IPD.TimesSaved    int16      Number of time saved
%     mr5               int16      No of rows in IPD.comments
%     mc5               int16      No of columns in IPD.comments
%     IPD.comments      uchar      Header comments
%     mr6               int16      No of rows in IPD.history
%  -> for i = 1:mr6  - - - - - - - - - - - - - - - - - - - - - - - - -
% |   lh                int16      Length of i'th row in IPD.history  |
% |   IPD.history{i}    uchar      i'th item in IPD.history           |
% |                                                                   |
% |_ _ _ _ _ _ _ _ _ Loop over mr6 entries _ _ _ _ _ _ _ _ _ _ _ _ _ _|
%
% Usage: save2MGP(IPD, MGPfile)
%
% Input: 
%        IPD : the MATGPR current input data structure
%    MPGfile : The complete path and file name of the output MGP file 
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
fid=fopen(MGPfile,'w',ENVAR.endian);

% Data origin
lo = length(IPD.origin);
fwrite(fid,lo,'int16');
fwrite(fid,IPD.origin,'uchar');
% Original data path 
lp = length(IPD.pname);
fwrite(fid,lp,'int16');
fwrite(fid,IPD.pname,'uchar');
% Original data file name
lf = length(IPD.fname);
fwrite(fid,lf,'int16');
fwrite(fid,IPD.fname,'uchar');
% Number of samples 
fwrite(fid,IPD.ns,'int16');
% Number of traces
fwrite(fid,IPD.ntr,'int16');
% Data matrix
fwrite(fid,IPD.d,'float');
% Sampling rate
if ~isempty(IPD.dt);
    fwrite(fid,IPD.dt,'float');
else
    s = 0;
    fwrite(fid,s,'float');
end
% Signal position
fwrite(fid,IPD.sigpos,'float');
% First depth estimate
if ~isempty(IPD.z);
    fwrite(fid,IPD.z(1),'float');
else
    s = 0;
    fwrite(fid,s,'float');
end
% Depth spacing 
if ~isempty(IPD.dz);
    fwrite(fid,IPD.dz,'float');
else
    s = 0;
    fwrite(fid,s,'float');
end
% Vertical axis label
lz = length(IPD.zlab);
fwrite(fid,lz,'int16');
fwrite(fid,IPD.zlab,'uchar');
% First trace location
fwrite(fid,IPD.x(1),'float');
% Trace spacing
if ~isempty(IPD.dx);
    fwrite(fid,IPD.dx,'float');
else
    s = 0;
    fwrite(fid,s,'float');
end
% Scan axis label
lx = length(IPD.xlab);
fwrite(fid,lx,'int16');
fwrite(fid,IPD.xlab,'uchar');
% Marker traces, if any
[mr, mc] = size(IPD.markertr);
fwrite(fid,mr,'int16');
fwrite(fid,mc,'int16');
if mr~=0 && mc~=0,
    fwrite(fid,IPD.markertr,'float');
end
% Transmitter coordinates in global (survey) reference frame
[mr, mc] = size(IPD.xyz.Tx);
fwrite(fid,mr,'int16');
fwrite(fid,mc,'int16');
if mr~=0 && mc~=0,
    fwrite(fid,IPD.xyz.Tx,'float');
end
% Receiver coordinates in global (survey) reference frame
[mr, mc] = size(IPD.xyz.Rx);
fwrite(fid,mr,'int16');
fwrite(fid,mc,'int16');
if mr~=0 && mc~=0,
    fwrite(fid,IPD.xyz.Rx,'float');
end
% Transmitter - Receiver separation
fwrite(fid,IPD.TxRx,'float');
% Antenna name if written in header files
la = length(IPD.Antenna);
fwrite(fid,la,'int16');
fwrite(fid,IPD.Antenna,'uchar');
% Header gain for DZT data
[mr, mc] = size(IPD.DZThdgain);
fwrite(fid,mr,'int16');
fwrite(fid,mc,'int16');
if mr ~=0 && mc ~=0,
    fwrite(fid,IPD.DZThdgain,'float');
end
% Number of time saved
fwrite(fid,IPD.TimesSaved,'int16');
% Header comments
[mr, mc] = size(IPD.comments);
fwrite(fid,mr,'int16');
fwrite(fid,mc,'int16');
if mr~=0 && mc~=0,
    fwrite(fid,IPD.comments,'uchar');
end
% Processing history
mr = size(IPD.history,1);
fwrite(fid,mr,'int16');
for i=1:mr
    s  = IPD.history{i};
    lh = length(s);
    fwrite(fid,lh,'int16');
    fwrite(fid,s,'uchar');
end

% Done
fclose(fid);
return

