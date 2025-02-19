function [DATA, SEGYreelhdr] = readsegy()  
%
% READSEGY : Import GPR data from file stored in the SEG-Y format, 
%            Revision 0 or Revision 1.
%
%    Usage : DATA = readsegy;     [DATA, SEGYreelhdr] = readsegy;
%
%  RETURNS : A MATGPR data structure DATA containing the GPR data and the
%            other necessary parameters and, optionally, a structure
%            SEGYreelhdr with the reel header of the SEGY file
%
% REQUIRES : initdatastr.m, readsuheader.m, ibm2num.m
%
%   CAVEAT : It is not recommended to use this program as is, in order to
%            read seismic data. For seismic data, either modify this
%            program or use something else. 
%
% Author : Andreas Tzanis
%          Department of Geophysics, 
%          University of Athens
%          atzanis@geol.uoa.gr
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
% the endian! Use the system's default
if isempty(ENVAR) || ~isfield(ENVAR,'endian'),
    [computer_system, maxsize, endian] = computer;
    endian = lower(endian);
    clear computer_system maxsize;
else
    endian = ENVAR.endian;
end

% Initialize the input data structure and header.
DATA          = initdatastr;
SEGYtracehdr  = initSEGYtraceheader(DATA);
SEGYdataformat= initSEGYdataformats;

%%% Get data file name 
[DATA.fname, DATA.pname]= uigetfile('*.segy; *.sgy; *.SEGY; *.SGY',...
    'Give SEG-Y data file name',ENVAR.currentworkdir);
if DATA.fname == 0, 
    erh = errordlg('Error in file specification!','READSEGY : ERROR');
    uiwait(erh)
    return, 
end;
%%% Data origin
DATA.origin = 'SEG-Y file';
%%% Open file and import data
fid=fopen([DATA.pname DATA.fname],'r',endian);
%%% ====== Read SEG-Y reel header =========================================
% Read 3200-byte textual header block (ASCII coded) 
SEGYreelhdr.headertext= fread(fid,3200,'uchar'); 
SEGYreelhdr.headertext= char(reshape(SEGYreelhdr.headertext, 40, 80));
% Read binary coded block  
SEGYreelhdr.jobid = fread(fid,1,'int32');  % Job identification number  
SEGYreelhdr.lino  = fread(fid,1,'int32');  % Line number (one line per reel)  
SEGYreelhdr.reno  = fread(fid,1,'int32');  % Reel number  
SEGYreelhdr.ntrpr = fread(fid,1,'int16');  % Number of data traces / record  
SEGYreelhdr.nart  = fread(fid,1,'int16');  % Number of auxiliary traces per 
                                           % record   
SEGYreelhdr.hdt   = fread(fid,1,'uint16'); % Sample interval in microseconds
                                           % for this reel  
SEGYreelhdr.dto   = fread(fid,1,'uint16'); % Same for original field record  
SEGYreelhdr.hns   = fread(fid,1,'uint16'); % Number of samples per trace 
                                           % for this reel  
SEGYreelhdr.nso   = fread(fid,1,'uint16'); % Number of samples per trace 
                                           % for original field recording  
SEGYreelhdr.format= fread(fid,1,'int16');  % Data sample format code: 
                                           % 1, IBM floating point (4 bytes) 
                                           % 2, Fixed-Point (4 bytes) 
                                           % 3, Fixed-Point (2 bytes) 
                                           % 5, IEEE floating point (4 bytes)
                                           % 8, Fixed-point two's complement
                                           %    (1 byte)
%%% If data format old or not supported, return
if SEGYreelhdr.format==4 || SEGYreelhdr.format==6 || SEGYreelhdr.format==7,
    erh = errordlg('Data format OBSOLETE or NOT SUPPORTED in MATGPR', ...
        'READ SEG-Y: Error');
    uiwait(erh);
    DATA = initdatastr;
    SEGYreelhdr = initSEGYreelheader(DATA);
    return
end
SEGYreelhdr.fold  = fread(fid,1,'int16');  % CDP fold expected per CDP 
                                           % ensemble 
SEGYreelhdr.tshort= fread(fid,1,'int16');  % Trace sorting code: 
                                           % 1 = as recorded (no sorting) 
                                           % 2 = CDP ensemble 
                                           % 3 = single fold continuous profile 
                                           % 4 = horizontally stacked 
SEGYreelhdr.vscode= fread(fid,1,'int16');  % Vertical sum code: 
                                           % 1 = no sum 
                                           % 2 = two sum ... 
                                           % N = N sum (N = 32,767) 
SEGYreelhdr.hsfs  = fread(fid,1,'int16');  % Sweep frequency at start  
SEGYreelhdr.hsfe  = fread(fid,1,'int16');  % Sweep frequency at end  
SEGYreelhdr.hslen = fread(fid,1,'int16');  % Sweep length (ms)  
SEGYreelhdr.hstyp = fread(fid,1,'int16');  % Sweep type code: 
                                           % 1 = linear,     2 = parabolic, 
                                           % 3 = exponential,    4 = other 
SEGYreelhdr.schn  = fread(fid,1,'int16');  % Trace number of sweep channel 
SEGYreelhdr.hstas = fread(fid,1,'int16');  % Sweep trace taper length (msec) 
                                           % at start if tapered (the taper
                                           % starts at zero time and is
                                           % effective for this length) 
SEGYreelhdr.hstae = fread(fid,1,'int16');  % Sweep trace taper length (msec) 
                                           % at end (the ending taper
                                           % starts at sweep length minus
                                           % the taper length at end)   
SEGYreelhdr.htatyp= fread(fid,1,'int16');  % Sweep trace taper type code: 
                                           % 1 = linear,   2 = cos-squared, 
                                           % 3 = other 
SEGYreelhdr.hcorr = fread(fid,1,'int16');  % Correlated data traces code: 
                                           % 1 2 3 4 5 6 7 8 95,
                                           % 1 = no 2 = yes 
SEGYreelhdr.bgrcv = fread(fid,1,'int16');  % Binary gain recovered code: 
                                           % 1 = yes, 2 = no 
SEGYreelhdr.rcvm  = fread(fid,1,'int16');  % Amplitude recovery method code: 
                                           % 1 = none, 
                                           % 2 = spherical divergence, 
                                           % 3 = AGC,         % 4 = other 
SEGYreelhdr.mfeet = fread(fid,1,'int16');  % Measurement system code: 
                                           % 1 = meters;    2 = feet 
SEGYreelhdr.polyt = fread(fid,1,'int16');  % Impulse signal polarity code: 
                                           % 1 = increase in pressure or
                                           %     upward geophone case
                                           %     movement gives negative
                                           %     number on tape  
                                           % 2 = increase in pressure or
                                           %     upward geophone case
                                           %     movement gives positive
                                           %     number on tape  
SEGYreelhdr.vpol  = fread(fid,1,'int16');  % Vibratory polarity code 
SEGYreelhdr.dt    = fread(fid,1,'float32');% Sampling rate (MATGPR specific 
                                           % field)
SEGYreelhdr.t1    = fread(fid,1,'float32');% First sample location (MATGPR 
                                           % specific field)
SEGYreelhdr.dx    = fread(fid,1,'float32');% Trace spacing (MATGPR specific 
                                           % field)
SEGYreelhdr.x1    = fread(fid,1,'float32');% First trace location (MATGPR 
                                           % specific field)
SEGYreelhdr.unass1= fread(fid,112,'int16');% Unassigned 224 bytes
SEGYreelhdr.revision=fread(fid,1,'uint16');% SEG-Y revision number
SEGYreelhdr.fixedlen=fread(fid,1,'int16'); % Fixed length trace flag
SEGYreelhdr.ntexthdrs=fread(fid,1,'int16');% Number of extended textual headers
SEGYreelhdr.unass2 = fread(fid,47,'int16');% Uassigned 94 bytes.
%%% ====== Done reading the reel header - Proceed with trace data =========

%%% Check fixed length flag - variable length taces not supported
if SEGYreelhdr.revision == 1 && SEGYreelhdr.fixedlen ~= 1,
    erh = errordlg('Variable trace length NOT SUPPORTED in MATGPR!',...
        'READSEGY : ERROR');
    uiwait(erh);
    DATA = initdatastr;
    SEGYreelhdr = initSEGYreelheader(DATA);
    return
end
%%% Check Revision Number against Data Sample Format. 
%%% Fixed point formats of revision 0 are not supported
if SEGYreelhdr.revision == 0 && SEGYreelhdr.format > 1,
    erh = errordlg('Revision = 0 and Format > 1 NOT SUPPORTED in MATGPR!',...
        'READSEGY : ERROR');
    uiwait(erh);
    DATA = initdatastr;
    SEGYreelhdr = initSEGYreelheader(DATA);
    return
end

%%% Begin taking data
hw = waitbar(0,['Importing from SEG-Y Rev.' ...
    num2str(SEGYreelhdr.revision) ', ' ...
    SEGYdataformat(SEGYreelhdr.format).name]);

%%% ========= REVISION 0 ==================================================
if SEGYreelhdr.revision == 0,
    %%% Read the first trace header to get basic information
    SEGYtracehdr  = readsuheader(fid);
    %%% Find if file is SU compliant by trying to verify that trace header
    %%% bytes 181-184 and 189-193 contain dt and dx respectively
    if SEGYtracehdr.d2 ~= 0 && SEGYtracehdr.d1 ~=0 ...
            && (round(SEGYtracehdr.d1 * 1000) == SEGYreelhdr.hdt), 
        disp('READSEGY > Importing from SU compliant Rev.0 data file')
        DATA.dt  =  SEGYtracehdr.d1;         % File is SU compliant
        DATA.ns  =  SEGYtracehdr.ns;
        DATA.ntr =  SEGYtracehdr.ntr;
        DATA.dx  =  SEGYtracehdr.d2;
        First_trace_location =  SEGYtracehdr.f2;
        DATA.sigpos  = SEGYtracehdr.delrt/1000;
    else                                     
        disp('READSEGY > Importing from unknown Rev.0 data file')
        DATA.dt = SEGYreelhdr.hdt/1000;      % Other SEG-Y revision 0 file  
        DATA.dx   = 0;                       % will resolve this later
        First_trace_location = 0;                
        fseek(fid,0,'eof');                  % find number of traces  
        lastbyte  = ftell(fid);                  
        databytes = lastbyte - 3600;
        tracebytes= 240 + SEGYdataformat(SEGYreelhdr.format).bytes;
        DATA.ntr  = databytes/tracebytes;
        DATA.sigpos  = 0; 
    end
    % Get transmitter - receiver offset
    DATA.TxRx =  SEGYtracehdr.offset/1000; 
    % Initialize work arrays
    DATA.d = zeros(DATA.ns,DATA.ntr);
    XT = [];   YT = [];   ZT = [];   XR = [];   YR = [];   ZR = [];   
    fseek(fid,3600,'bof');
    for i=1:DATA.ntr,
        SEGYtracehdr = readsuheader(fid);
        % Test scaling factors - if not set assign the default value = -3
        if SEGYtracehdr.scalel == 0,
            SEGYtracehdr.scalel = -3;
        end
        if SEGYtracehdr.scalco == 0,
            SEGYtracehdr.scalco = -3;
        end
        ZT(i) = SEGYtracehdr.selev * 10^SEGYtracehdr.scalel;
        XT(i) = SEGYtracehdr.sx * 10^SEGYtracehdr.scalco;
        YT(i) = SEGYtracehdr.sy * 10^SEGYtracehdr.scalco;
        ZR(i) = SEGYtracehdr.gelev * 10^SEGYtracehdr.scalel;
        XR(i) = SEGYtracehdr.gx * 10^SEGYtracehdr.scalco;
        YR(i) = SEGYtracehdr.gy * 10^SEGYtracehdr.scalco;
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
        if SEGYtracehdr.mark ~= 0,
        DATA.markertr = [DATA.markertr; i (XT(i)+XR(i))/2 ...
            (YT(i)+YR(i))/2  (ZT(i)+ZR(i))/2 ];
        end
        % Read trace data in IBM floating point format
        tr  = fread(fid,DATA.ns,'uint32');
        tr1 = ibm2num(uint32(tr));
        DATA.d(:,i) = tr1;
        waitbar(i/DATA.ntr,hw);
    end
end

%%% ========= REVISION 1 ==================================================
if SEGYreelhdr.revision == 1, 
    %%% Take into account the possibility of extended textual headers,
    %%% although not supported in current version of MATGPR
    DataBlock = 3600 + 3200 * SEGYreelhdr.ntexthdrs;
    %%% Read the first trace header to get basic information
    fseek(fid, DataBlock, 'bof');
    SEGYtracehdr = readsegytraceheader(fid);
    % Get number of samples per trace
    DATA.ns   = SEGYreelhdr.hns;
    %%% Find if file is MATGPR compliant - use several tests
    if (SEGYreelhdr.dx ~=0 && SEGYtracehdr.dx~=0) ...     
            && SEGYreelhdr.dx == SEGYtracehdr.dx ...   
               && round(SEGYreelhdr.dt * 1000) == SEGYtracehdr.dt,  
        disp('READSEGY > Importing from MATGPR compliant Rev.1 data file')
        DATA.dt = SEGYreelhdr.dt;                % File is MATGPR compliant
        DATA.dx   = SEGYreelhdr.dx;              
        First_trace_location = SEGYreelhdr.x1;
        DATA.ntr  = SEGYtracehdr.ntr; 
        DATA.sigpos  = SEGYtracehdr.delrt/1000;
    else                                         % File is other SEG-Y 
        disp('READSEGY > Importing from unknown Rev.1 data file')
        DATA.dt = SEGYreelhdr.hdt/1000;          
        DATA.dx   = 0;                           % will resolve this later
        First_trace_location = 0;                
        fseek(fid,0,'eof');                      % find number of traces  
        lastbyte  = ftell(fid);                  
        databytes = lastbyte - DataBlock;
        tracebytes= 240 + SEGYdataformat(SEGYreelhdr.format).bytes;
        DATA.ntr  = databytes/tracebytes;
        DATA.sigpos  = 0; 
    end
    % Get transmitter - receiver offset
    DATA.TxRx = SEGYtracehdr.offset/1000;
    % Initialize work arrays
    DATA.d = zeros(DATA.ns,DATA.ntr);
    XT = [];   YT = [];   ZT = [];   XR = [];   YR = [];   ZR = [];   
    % Rewind to end of reel header block and start reading traces
    fseek(fid, DataBlock, 'bof');
    for i=1:DATA.ntr
        % Read trace header
        SEGYtracehdr = readsegytraceheader(fid); 
        % Test scaling factors - if not set, assign the default value = -3
        if SEGYtracehdr.scalel == 0,
            SEGYtracehdr.scalel = -3;
        end
        if SEGYtracehdr.scalco == 0,
            SEGYtracehdr.scalco = -3;
        end
        ZT(i) = SEGYtracehdr.selev * 10^SEGYtracehdr.scalel;
        XT(i) = SEGYtracehdr.sx * 10^SEGYtracehdr.scalco;
        YT(i) = SEGYtracehdr.sy * 10^SEGYtracehdr.scalco;
        ZR(i) = SEGYtracehdr.gelev * 10^SEGYtracehdr.scalel;
        XR(i) = SEGYtracehdr.gx * 10^SEGYtracehdr.scalco;
        YR(i) = SEGYtracehdr.gy * 10^SEGYtracehdr.scalco;
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
        if SEGYtracehdr.mark ~= 0,
        DATA.markertr = [DATA.markertr; i (XT(i)+XR(i))/2 ...
            (YT(i)+YR(i))/2  (ZT(i)+ZR(i))/2 ];
        end
        % Read trace data 
        if SEGYreelhdr.format == 1,         % IBM floating point format
            tr  = fread(fid,DATA.ns,'uint32');
            tr1 = ibm2num(uint32(tr));
        elseif SEGYreelhdr.format > 1,      % All other formats
            tr1 = fread(fid,DATA.ns,SEGYdataformat(SEGYreelhdr.format).format);
        end
        DATA.d(:,i)  =  tr1;
        waitbar(i/DATA.ntr,hw);
    end
end
close(hw)
fclose(fid);
% if No of samples is odd, remove last to make it even - easier for FFT's,
% filtering and such ...
if mod(DATA.ns,2) ~=0,
    DATA.d  = DATA.d(1:DATA.ns - 1,:);
    DATA.ns = DATA.ns -1;
end

%%% Two-way traveltime
DATA.tt2w   = 0 : DATA.dt : DATA.dt*(DATA.ns - 1);
DATA.zlab   = 'Traveltime (ns)';

%%% If trace spacing dx = 0 (read error or non-MATGPR compliant SEG-Y file)
%%% try to work out a solution using the trace coordinates. If this doesn't
%%% work either, dx will remain naught 
if DATA.dx == 0 && ~isempty(find(XT + YT)),
    ds = sqrt((XT-XT(1)).^2 + (YT-YT(1)).^2);
    DATA.dx = mean(diff(ds));
end
%%% Create the scan axis vector
if DATA.dx == 0,
    DATA.x    = 1:1:DATA.ntr;
    DATA.xlab = 'Scan Axis (# Traces)';
else 
    DATA.x    = First_trace_location : DATA.dx : First_trace_location + ...
        ((DATA.ntr-1)*DATA.dx);
    DATA.xlab = 'Scan Axis (meters)';
end

%%% Check for marker traces and report
if isempty(DATA.markertr),      
    disp('READSEGY > THERE ARE NO MARKER TRACES');
end
if ~isempty(DATA.markertr),
    if isempty(find(DATA.markertr(:,2) + DATA.markertr(:,3))),
        disp('READSEGY > THERE ARE MARKER TRACES WITHOUT TOPOGRAPHIC INFO')
        DATA.markertr = DATA.markertr(:,1);
    else
        disp('READSEGY > THERE ARE MARKER TRACES WITH TOPOGRAPHIC INFO ')
    end
end
if isempty(find(XT + YT)),
    disp('READSEGY > XYZ TRACE COORDINATES NOT INCLUDED IN THIS SATA SET');
end
if ~isempty(find(XT + YT )),
    disp('READSEGY > XYZ TRACE COORDINATES ARE INCLUDED IN THIS DATA SET');
    DATA.xyz.Tx = [XT' YT' ZT'];
    DATA.xyz.Rx = [XR' YR' ZR'];
end

% Antenna name will generally not be stored in SU format
DATA.Antenna = 'NA (imported from SEG-Y format)';        
DATA.comments= SEGYreelhdr.headertext;
%%% And that's all folks ...
return

function SEGYdataformat = initSEGYdataformats()
%%% Initialize the data sample format structure 
SEGYdataformat(1).format='uint32'; 
SEGYdataformat(2).format='int32';  
SEGYdataformat(3).format='int16'; 
SEGYdataformat(4).format='';  
SEGYdataformat(5).format='float32'; 
SEGYdataformat(6).format='';  
SEGYdataformat(7).format='';  
SEGYdataformat(8).format='int8';  

SEGYdataformat(1).bytes = 4;
SEGYdataformat(2).bytes = 4;  
SEGYdataformat(3).bytes = 2; 
SEGYdataformat(4).bytes = 0;  
SEGYdataformat(5).bytes = 4; 
SEGYdataformat(6).bytes = 0;  
SEGYdataformat(7).bytes = 0;  
SEGYdataformat(8).bytes = 1;  

SEGYdataformat(1).name='4-byte IBM Floating Point';
SEGYdataformat(2).name='4-byte two''s complement';
SEGYdataformat(3).name='2-byte two''s complement';
SEGYdataformat(4).name='4-byte fixed-point with gain - NOT SUPPORTED';
SEGYdataformat(5).name='4-byte IEEE floating-point';
SEGYdataformat(6).name='DSF Not currently used';
SEGYdataformat(7).name='DSF Not currently used';
SEGYdataformat(8).name='1-byte, two''s complement integer';
return
