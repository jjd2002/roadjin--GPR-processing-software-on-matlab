function [DATA, HEADER] = readrd3()
%
% READRD3 : Imports SINGLE channel RAMAC (.RD3) georadar data.
%
%   Usage : DATA = readrd3;     [DATA, HEADER] = readrd3;
%
% RETURNS : A structure DATA containg the GPR data and all other necessary 
%           parameters and, optionally, a structure HEADER containg the 
%           header of the imported RD3 data file taken from the associated
%           file *.RAD
%
% REQUIRES: initdatastr.m
%
%  Author : Andreas Tzanis, 
%           Department of Geophysics, 
%           University of Athens
%           atzanis@geol.uoa.gr
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

global ENVAR
% If the function is running alone, ENVAR was not initialized. Must define
% the endian! Use the system's default
if isempty(ENVAR) || ~isfield(ENVAR,'endian'),
    [computer_system, maxsize, endian] = computer;
    endian = lower(endian);
    clear computer_system maxsize;
else
    endian = ENVAR.endian;
end

% Initialize the input data structure.
DATA = initdatastr;

% First tell the origin and type of the original data
DATA.origin = 'RAMAC Monostatic';
%%% Get data file name 
[DATA.fname, DATA.pname]= uigetfile('*.rd3; *.RD3',...
    'Give RAMAC data file name', ENVAR.currentworkdir);
if DATA.fname == 0, 
    return, 
end;
%%% Make header and marker file names 
Inhdfname = [DATA.fname(1:findstr(DATA.fname,'.')) 'rad'];
Inmkfname = [DATA.fname(1:findstr(DATA.fname,'.')) 'mkn'];

%%% Now open and read header file
disp('READRD3 > Reading and decoding header information')
fid  = fopen([DATA.pname Inhdfname],'r');
HEADER  = readRD3header(fid);
%%% Assign values to the DATA structure
DATA.ns      = HEADER.Samples_per_scan;
DATA.dt      = 1000/HEADER.Frequency;
DATA.sigpos  = HEADER.Signal_position;
DATA.dx      = HEADER.Distance_interval;
if DATA.dx == 0,
    DATA.dx = [];
end
DATA.Antenna = HEADER.Antennas;
DATA.TxRx    = HEADER.Antenna_separation;
DATA.ntr     = HEADER.Last_trace;
%%%% Assign header comments (but first wrap to 80 column text)
hc=[]; 
for i=1:ceil(length(HEADER.Comment)/80), 
    if i*80 >length(HEADER.Comment), 
        hc=strvcat(hc,HEADER.Comment((i-1)*80+1:length(HEADER.Comment))); 
    else 
        hc=strvcat(hc,HEADER.Comment((i-1)*80+1:i*80)); 
    end; 
end;
DATA.comments = str2mat('ORIGINAL HEADER FILE COMMENTS:',hc);
fclose(fid);                           
disp('        ==> OK')

%%% Now open and read data file
disp('READRD3 > Reading data')
fid    = fopen([DATA.pname DATA.fname],'r',endian);
DATA.d = fread(fid,[DATA.ns,DATA.ntr],'short');
fclose(fid);
disp('        ==> OK')

%%%% Two-way traveltime
DATA.tt2w       = 0 : DATA.dt : DATA.dt*(DATA.ns - 1);
DATA.zlab       = 'Traveltime (ns)';
%%%% Scan Line axis
if isempty(DATA.dx),
    DATA.x      = 1 : 1 : DATA.ntr;
    DATA.xlab   = 'Scan Axis (# Traces)';
else 
    DATA.x      = 0 : DATA.dx : (DATA.ntr - 1)*DATA.dx;
    DATA.xlab   = 'Scan Axis (meters)';
end

%%%% Check and import marker information 
if ~exist([DATA.pname Inmkfname],'file'), 
    return
end
disp('READRD3 > Checking for MARKER TRACES')
fid = fopen([DATA.pname Inmkfname],'r');
for i=1:11,
    hdline=fgetl(fid);
end
i=0;
while ~feof(fid),
    i = i + 1;
    hdline=fgetl(fid);
    dummy = str2num(hdline);
    DATA.markertr(i) = dummy(1);
end
fclose(fid);
if isempty(DATA.markertr),
    disp('        ==> THERE ARE NO MARKED TRACES')
else
    DATA.markertr = DATA.markertr(:);
    disp(['        ==> FOUND ' num2str(length(DATA.markertr)) ...
        ' MARKER TRACES'])
end;    
%%% And that's all folks ...
return

function HDR = readRD3header(fid)
%%%% Import the RAMAC header structure
hdline         = fgetl(fid);
HDR.Samples_per_scan   = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline         = fgetl(fid);
HDR.Frequency          = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline         = fgetl(fid);
HDR.Frequency_steps    = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline         = fgetl(fid);
HDR.Signal_position    = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
HDR.Raw_Signal_position= str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
HDR.Distance_flag      = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
HDR.Time_flag          = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
HDR.Program_flag       = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
HDR.External_flag      = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
HDR.Time_interval      = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
HDR.Distance_interval  = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
HDR.Operator           = fgetl(fid);
HDR.Customer           = fgetl(fid);
HDR.Site               = fgetl(fid);
hdline = fgetl(fid);
HDR.Antennas           = hdline(findstr(hdline,':')+1:length(hdline));
hdline                 = fgetl(fid); % Antenna orientation is not yet valid field
%HDR.Antenna_orientation= str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline                 = fgetl(fid);
HDR.Antenna_separation = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
HDR.Comment            = fgetl(fid);
hdline = fgetl(fid);
HDR.Timewindow         = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
HDR.Stacks             = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
HDR.Stack_exponent     = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
HDR.Stacking_time      = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
HDR.Last_trace         = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
HDR.Stop_position      = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
%%% Merge all text into a single char array.
HDR.Comment = cat(2,HDR.Operator,'  ',HDR.Customer,'  ',HDR.Site,'  ',...
    HDR.Comment);
while ~feof(fid),
    hdline = fgetl(fid);
    HDR.System_calibration = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    HDR.Start_position     = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    HDR.Short_flag         = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    HDR.Intermediate_flag  = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    HDR.Long_flag          = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    HDR.Preprocessing      = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    HDR.High               = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    HDR.Low                = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    HDR.Fixed_increment    = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    HDR.Fixed_moves_up     = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    HDR.Fixed_moves_down   = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    HDR.Fixed_position     = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    if feof(fid), 
        HDR.Wheel_calibration  = [];
        HDR.Positive_direction = [];
        break, 
    end;
    HDR.Wheel_calibration  = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    HDR.Positive_direction = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
end                 % end of file loop
return
