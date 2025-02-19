function [DATA, HEADER] = readdt1()
%
% READDT1 : Imports single channel PULSE EKKO (.DT1) GPR data.
%
%   Usage : DATA = readdt1;     [DATA, HEADER] = readdt1;
%
% RETURNS : A structure DATA containing the GPR data and other necessary 
%           parameters and, optionally, a structure HEADER containing the 
%           information about the imported DT1 data file, as obtained from
%           the associated header file *.HD 
%
% Requires: initdatastr.m, dt1read.m 
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
% Initialize the input data structure.
DATA = initdatastr;

% Assign data origin
DATA.origin = 'PULSE EKKO Bistatic';

%%% Get data file name 
[DATA.fname, DATA.pname]= uigetfile('*.dt1; *.DT1',...
    'Give PULSE EKKO data file name', ENVAR.currentworkdir);
if DATA.fname == 0, 
    return, 
end;

%%% Make header file name 
Inhdfname = [DATA.fname(1:findstr(DATA.fname,'.')) 'hd'];

%%% Now open and read header file
disp('READDT1 : Reading and decoding header information')
fid  = fopen([DATA.pname Inhdfname],'r');
% Read header file
HEADER = readDT1hdfile(fid);
DATA.ntr    = HEADER.Num_traces;
DATA.ns     = HEADER.Samples_per_scan;
DATA.dt     = HEADER.Timewindow/HEADER.Samples_per_scan;
DATA.sigpos = HEADER.Signal_position * DATA.dt;
DATA.dx     = HEADER.Step_size_used;
DATA.Antenna= HEADER.Nominal_Frequency;
DATA.TxRx   = HEADER.Antenna_separation;
fclose(fid);                           
%%%% Assign header comments (but first wrap to 80 column text)
hc=[]; 
for i=1:ceil(length(HEADER.comment)/80), 
    if i*80 >length(HEADER.comment), 
        hc=strvcat(hc,HEADER.comment((i-1)*80+1:length(HEADER.comment))); 
    else 
        hc=strvcat(hc,HEADER.comment((i-1)*80+1:i*80)); 
    end; 
end;
DATA.comments = str2mat('ORIGINAL HEADER FILE COMMENTS:',hc);
disp('        ==> OK')

%%% Now open and read data file
disp('READDT1 : Reading data')
[DATA.d, trheaders] = dt1read([DATA.pname DATA.fname]);
disp('        ==> OK')

%%%% Two-way traveltime
DATA.tt2w    = 0 : DATA.dt : DATA.dt*(DATA.ns - 1);
DATA.zlab    = 'Traveltime (ns)';
%%%% Scan Line axis
if isempty(DATA.dx),
    DATA.x   = 1 : 1 : DATA.ntr;
    DATA.xlab= 'Scan Axis (# Traces)';
else 
    DATA.x   = HEADER.Starting_position : DATA.dx : (DATA.ntr - 1)*DATA.dx;
    DATA.xlab= ['Scan Axis (' HEADER.position_units ')'];
end

%%%% Source locations
DATA.xyz.Tx = zeros(DATA.ntr,3);
DATA.xyz.Rx = zeros(DATA.ntr,3);
for i=1:DATA.ntr
    DATA.xyz.Tx(i,:)=[trheaders(i).x_tra trheaders(i).y_tra ...
        trheaders(i).z_tra];
    DATA.xyz.Rx(i,:)=[trheaders(i).x_rec trheaders(i).y_rec ...
        trheaders(i).z_rec];
end

disp('        ==> THERE ARE NO MARKER TRACES IN PULSE EKKO DATA')
%%% And that's all folks ...
return

function HDR = readDT1hdfile(fid)
%%%% Import the Pulse Ekko header file in structure HDR
tmpcell = cell(1);
j=0;
while ~feof(fid),
    hdline = fgetl(fid);
    if ~isempty(hdline),
        j = j+1;
        tmpcell(j) = cellstr(hdline);
    end
end
HDR.comment = ' ';
for j=1:5,
    hdline = [tmpcell{j} '          '];
    if ~(strcmp(hdline(1:6),'NUMBER')), 
        HDR.comment = cat(2,HDR.comment, ' ', hdline);
    else
        break
    end
end;
HDR.Num_traces       =  str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
i = j+1;    hdline   = tmpcell{i};
HDR.Samples_per_scan = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
i = i+1;    hdline   = tmpcell{i};
HDR.Signal_position  = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
i = i+1;    hdline   = tmpcell{i};
HDR.Timewindow       = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
i = i+1;    hdline   = tmpcell{i};
HDR.Starting_position= str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
i = i+1;    hdline   = tmpcell{i};
HDR.Final_position   = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
i = i+1;    hdline   = tmpcell{i};
HDR.Step_size_used   = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
i = i+1;    hdline   = tmpcell{i};
HDR.position_units   = hdline(findstr(hdline,'=')+1:length(hdline));
i = i+1;    hdline   = tmpcell{i};
HDR.Nominal_Frequency= hdline(findstr(hdline,'=')+1:length(hdline));
i = i+1;    hdline   = tmpcell{i};
HDR.Antenna_separation= str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
i = i+1;    hdline   = tmpcell{i};
HDR.Pulser_voltage   = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
i = i+1;    hdline   = tmpcell{i};
HDR.Number_of_stacks = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
i = i+1;    hdline   = tmpcell{i};
HDR.Survey_mode      = hdline(findstr(hdline,'=')+1:length(hdline));
for j = i+1:length(tmpcell),
    HDR.comment = cat(2, HDR.comment, tmpcell{j});
end
return

%function HDR = readDT1hdfile(fid)
%%%% Import the Pulse Ekko header file in structure HDR
%HDR.comment = ' ';
%for i=1:5,
%    hdline = fgetl(fid);
%    hdline = [hdline '          '];
%    if ~(strcmp(hdline(1:6),'NUMBER')), 
%        HDR.comment = cat(2,HDR.comment, ' ', hdline);
%    else
%        break
%    end
%end;
%HDR.Num_traces       =  str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
%hdline      = fgetl(fid);
%HDR.Samples_per_scan = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
%hdline      = fgetl(fid);
%HDR.Signal_position  = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
%hdline         = fgetl(fid);
%HDR.Timewindow       = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
%hdline      = fgetl(fid);
%HDR.Starting_position= str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
%hdline = fgetl(fid);
%HDR.Final_position   = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
%hdline      = fgetl(fid);
%HDR.Step_size_used   = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
%hdline      = fgetl(fid);
%HDR.position_units   = hdline(findstr(hdline,'=')+1:length(hdline));
%hdline      = fgetl(fid);
%HDR.Nominal_Frequency= hdline(findstr(hdline,'=')+1:length(hdline));
%hdline      = fgetl(fid);
%HDR.Antenna_separation= str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
%hdline      = fgetl(fid);
%HDR.Pulser_voltage   = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
%hdline      = fgetl(fid);
%HDR.Number_of_stacks = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
%hdline      = fgetl(fid);
%HDR.Survey_mode      = hdline(findstr(hdline,'=')+1:length(hdline));
%while ~feof(fid),
%    hdline = fgetl(fid);
%    HDR.comment = cat(2, HDR.comment, hdline);
%end
%return
