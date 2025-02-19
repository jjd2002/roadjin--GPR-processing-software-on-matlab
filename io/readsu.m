function DATA = readsu()
%
%  READSU : Imports GPR data and related information from file stored in
%           the Seismix Unix format (.SU)
%
%   Usage : DATA = readsu;
%
% RETURNS : A MATGPR data structure DATA containg the GPR data and other
%           necessary parameters. 
%
%REQUIRES : initdatastr.m, initSUheader.m, readsuheader.m
%
%Author   : Andreas Tzanis
%            Department of Geophysics, 
%            University of Athens
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

% Get global environmental variables
global ENVAR
% If the function is running alone, ENVAR was not initialized. Must define
% the endian!  
if ~isstruct(ENVAR) || ~isfield(ENVAR,'endian'),
    [computer_system, maxsize, endian] = computer;
    endian = lower(endian);
    currentworkdir = [pwd '/'];
    clear computer_system maxsize;
else
    endian = ENVAR.endian;
    currentworkdir = ENVAR.currentworkdir;
end

% Initialize the input data structure and header.
DATA = initdatastr;
TRACEHDR = initSUheader(DATA);

%%% Get data file name
[DATA.fname, DATA.pname]= uigetfile('*.su; *.SU',...
    'Give SU data file name', currentworkdir);
if DATA.fname == 0, 
    erh = errordlg('Error in file specification!','READSU: ERROR');
    uiwait(erh)
    return, 
end;
%%% Data origin
DATA.origin = 'SU File';
%%% Open file and import data
fid=fopen([DATA.pname DATA.fname],'r',endian);
hw = waitbar(0,'Importing from SU format');
% Read the first header to get basic information
TRACEHDR = readsuheader(fid);
DATA.sigpos = TRACEHDR.delrt/1000;  % signal position (but be careful!)
DATA.dt     = TRACEHDR.d1;
DATA.ns     = TRACEHDR.ns;
DATA.ntr    = TRACEHDR.ntr;
if DATA.ntr == 0,                     % This is possible for some data sets
    fseek(fid,0,'eof');
    lastbyte = ftell(fid);
    bytespertrace = 240 + DATA.ns*4;
    DATA.ntr = lastbyte/bytespertrace;
end
if TRACEHDR.d2 == 0,
    DATA.dx = [];
else
    DATA.dx     = TRACEHDR.d2;
end
% Test scaling factor - if not set, assign the default value of -3
if TRACEHDR.scalco == 0,
    TRACEHDR.scalco = -3;
end
DATA.TxRx   = TRACEHDR.offset * 10^TRACEHDR.scalco; 
First_trace_location =  TRACEHDR.f2;

% Initialize work arrays
DATA.d = zeros(DATA.ns,DATA.ntr);
XT = zeros(DATA.ntr,1);  YT = XT;   ZT = XT;   XR = XT;   YR = XT;   ZR = ZT;   
fseek(fid,0,'bof');
for i = 1:DATA.ntr,
    TRACEHDR = readsuheader(fid);
    % Test scaling factors - if not set, assign the default value = -3
    if TRACEHDR.scalel == 0,
        TRACEHDR.scalel = -3;
    end
    if TRACEHDR.scalco == 0,
        TRACEHDR.scalco = -3;
    end
    ZT(i) = TRACEHDR.selev * 10^TRACEHDR.scalel;
    XT(i) = TRACEHDR.sx * 10^TRACEHDR.scalco;
    YT(i) = TRACEHDR.sy * 10^TRACEHDR.scalco;
    ZR(i) = TRACEHDR.gelev * 10^TRACEHDR.scalel;
    XR(i) = TRACEHDR.gx * 10^TRACEHDR.scalco;
    YR(i) = TRACEHDR.gy * 10^TRACEHDR.scalco;
    if ZT(i) ~= 0 && ZR == 0,
        ZR(i) = ZT(i);
    end
    if ZT(i) == 0 && ZR(i) ~= 0,
        ZT(i) = ZR(i);
    end
    if (XT(i)+YT(i)) ~= 0 && (XR(i)+YR(i)) == 0,
        XR(i) = XT(i);
        YR(i) = YT(i);
    end
    if (XT(i)+YT(i)) == 0 && (XR(i)+YR(i)) ~= 0,
        XT(i) = XR(i);
        YT(i) = YR(i);
    end
    if TRACEHDR.mark ~= 0,
        DATA.markertr = [DATA.markertr; [i (XT(i)+XR(i))/2 ...
            (YT(i)+YR(i))/2  (ZT(i)+ZR(i))/2] ];
    end
    tr   = fread(fid,DATA.ns,'float');
    DATA.d(:,i) = tr;
    waitbar(i/DATA.ntr,hw);
end
fclose(fid);
close(hw);

% if No of samples is odd, reomove last to make it even - easier for FFT's,
% filtering and such ...
if mod(DATA.ns,2) ~=0,
    DATA.d  = DATA.d(1:DATA.ns - 1,:);
    DATA.ns = DATA.ns -1;
end

%%%% Two-way traveltime
DATA.tt2w   = 0 : DATA.dt : DATA.dt*(DATA.ns - 1);
DATA.zlab   = 'Traveltime (ns)';
%%%% Scan Line axis
if isempty(DATA.dx),
    DATA.x    = 1:1:DATA.ntr;
    DATA.xlab = 'Scan Axis (# Traces)';
else 
    DATA.x    = First_trace_location : DATA.dx : First_trace_location + ...
        ((DATA.ntr-1)*DATA.dx);
    DATA.xlab = 'Scan Axis (meters)';
end

%%% Check for marked traces and report
if isempty(DATA.markertr),      
    disp('READSU > THERE ARE NO MARKER TRACES');
end
if ~isempty(DATA.markertr),
    if isempty(find(DATA.markertr(:,2) + DATA.markertr(:,3))),
        disp('READSU > THERE ARE MARKER TRACES WITHOUT TOPOGRAPHIC INFO')
        DATA.markertr = Data.markertr(:,1);
    else
        disp('READSU > THERE ARE MARKER TRACES WITH TOPOGRAPHIC INFO ')
    end
end
if isempty(find(XT + YT)),
    disp('READSU > XYZ TRACE COORDINATES NOT INCLUDED IN THIS SATA SET');
end
if ~isempty(find(XT + YT )),
    disp('READSU > XYZ TRACE COORDINATES ARE INCLUDED IN THIS DATA SET');
    DATA.xyz.Tx = [XT YT ZT];
    DATA.xyz.Rx = [XR YR ZR];
end

% Antenna name will generally not be stored in SU format
DATA.Antenna = 'NA (SU format)';        
%%% And that's all folks ...
return
