function [DATA, SEGYreelhdr] = readZONDsegydata()  
%
% READZONDSEGYDATA: 
%            Import ZOND GPR data from the particular SEG-Y format created, 
%            by ZOND's own Prism v2.5 analysis pakage.
%            SEG-Y Revision 0
%
%    Usage : DATA = readZONDsegydata;
%            [DATA, SEGYreelhdr] = readZONDsegydata;
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
DATA.origin = 'ZOND SEG-Y file';
%%% Open file and import data
fid=fopen([DATA.pname DATA.fname],'r',endian);
%%% ====== Read SEG-Y reel header =========================================
% Read 3200-byte textual header block (ASCII coded) 
SEGYreelhdr.headertext= fread(fid,3200,'uchar');
SEGYreelhdr.headertext= char(SEGYreelhdr.headertext');
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
    erh = errordlg('Data format OBSOLETE and NOT SUPPORTED in MATGPR', ...
        'READ ZOND SEG-Y: Error');
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
SEGYreelhdr.unass1= fread(fid,120,'int16');% Unassigned 224 bytes
SEGYreelhdr.revision=fread(fid,1,'uint16');% SEG-Y revision number
SEGYreelhdr.fixedlen=fread(fid,1,'int16'); % Fixed length trace flag
SEGYreelhdr.ntexthdrs=fread(fid,1,'int16');% Number of extended textual headers
SEGYreelhdr.unass2 = fread(fid,47,'int16');% Uassigned 94 bytes.
%%% ====== Done reading the reel header - Proceed with trace data =========

%%% Check fixed length flag - variable length taces not supported
if SEGYreelhdr.revision == 1 && SEGYreelhdr.fixedlen ~= 1,
    erh = errordlg('Variable trace length NOT SUPPORTED in MATGPR!',...
        'READ ZOND SEG-Y : ERROR');
    uiwait(erh);
    DATA = initdatastr;
    SEGYreelhdr = initSEGYreelheader(DATA);
    return
end
%%% Check Revision Number against Data Sample Format. 
%%% Fixed point formats of revision 0 are not supported
%if SEGYreelhdr.revision == 0 && SEGYreelhdr.format > 1,
%    erh = errordlg('Revision = 0 and Format > 1 NOT SUPPORTED in MATGPR!',...
%        'READSEGY : ERROR');
%    uiwait(erh);
%    DATA = initdatastr;
%    SEGYreelhdr = initSEGYreelheader(DATA);
%    return
%end

%%% Begin taking data
hw = waitbar(0,['Importing from SEG-Y Rev.' ...
    num2str(SEGYreelhdr.revision) ', ' ...
    SEGYdataformat(SEGYreelhdr.format).name]);
%%% Read the first trace header to get basic information
SEGYtracehdr  = readZONDsuheader(fid);
DATA.dt     =  SEGYtracehdr.dtsu/1000;              % File is SU compliant
DATA.ns     =  SEGYtracehdr.ns;
DATA.sigpos = SEGYtracehdr.delrt/1000;
fseek(fid,0,'eof');                  % find number of traces  
lastbyte  = ftell(fid);                  
databytes = lastbyte - 3600;
tracebytes= 240 + DATA.ns*SEGYdataformat(SEGYreelhdr.format).bytes;
DATA.ntr  = databytes/tracebytes;
% Get transmitter - receiver offset
DATA.TxRx =  SEGYtracehdr.offset/1000; 
% Initialize work arrays
DATA.d = zeros(DATA.ns,DATA.ntr);
DATA.x = zeros(1,DATA.ntr);
fseek(fid,3600,'bof');
for i=1:DATA.ntr,
    SEGYtracehdr = readZONDsuheader(fid);
% %      disp([ SEGYtracehdr.sx  SEGYtracehdr.sy  SEGYtracehdr.scalco])
% %      disp([ SEGYtracehdr.gx  SEGYtracehdr.gy  SEGYtracehdr.scalco])
% disp([ SEGYtracehdr.gelev  SEGYtracehdr.scalco])
    % Test scaling factors - if not set assign the default value = -3
    if abs(SEGYtracehdr.scalco) > 9,
        SEGYtracehdr.scalco = ...
            sign(SEGYtracehdr.scalco)*log10(abs(SEGYtracehdr.scalco));
    end
    if SEGYtracehdr.scalco == 0,
        SEGYtracehdr.scalco = -3;
    end
    % ZOND stores trace coordinates in the SEG-Y field reserved for
    % receiver group x-coordinates. All the other fields pertaining to
    % source and receiver coordinates are naught
    DATA.x(i) = SEGYtracehdr.gx * 10^SEGYtracehdr.scalco;
    %
    if i > 1 && SEGYtracehdr.mark ~= 0,
        DATA.markertr = [DATA.markertr; i DATA.x(i) 0 0];
    end
    % Read trace data 
    if SEGYreelhdr.format == 1,         % IBM floating point format
        tr  = fread(fid,DATA.ns,'uint32');
        tr1 = ibm2num(uint32(tr));
    elseif SEGYreelhdr.format > 1,      % All other formats
        tr1 = fread(fid,DATA.ns,SEGYdataformat(SEGYreelhdr.format).format);
    end
    DATA.d(:,i) = tr1;
    waitbar(i/DATA.ntr,hw);
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
%%% Trace spacing
DATA.dx = mean(diff(DATA.x));
DATA.xlab = 'Scan Axis (meters)';

%%% Check for marker traces and report
if isempty(DATA.markertr),      
    disp('READ ZOND SEG-Y > THERE ARE NO MARKER TRACES');
end
if ~isempty(DATA.markertr),
    disp('READ ZOND SEG-Y > THERE ARE MARKER TRACES ')
end

% Antenna info is stored in the Reel Header
str2 = findstr(SEGYreelhdr.headertext,'MHz');
if isempty(str2)
    str2 = findstr(SEGYreelhdr.headertext,'Mhz');
end
ws = find(isspace(SEGYreelhdr.headertext(str2-6:str2)));
str1 = str2 - 6 + ws;
DATA.Antenna = SEGYreelhdr.headertext(str1:str2+2);
% Dump the rest of the reel header to the "comments filed of IPD
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

function SUHDR = readZONDsuheader(fid)
%
% READSUHEADER : Read SU trace headers and returns their content in the
%                structure SUHDR 
%
%       Author : Andreas Tzanis, 
%                Department of Geophysics, 
%                University of Athens
%                (C) 2005, Andreas Tzanis, allrights reserved
%
%
SUHDR.tracl   = fread(fid,1,'int32');   % Trace sequence number within line
SUHDR.tracr   = fread(fid,1,'int32');   % Trace sequence number within reel
SUHDR.fldr    = fread(fid,1,'int32');   % Field record number
SUHDR.tracf   = fread(fid,1,'int32');   % Trace number within field record
SUHDR.ep      = fread(fid,1,'int32');   % Energy source point number
SUHDR.cdp     = fread(fid,1,'int32');   % CDP ensemble number 
SUHDR.cdpt    = fread(fid,1,'int32');   % Trace number within CDP ensemble 
SUHDR.trid    = fread(fid,1,'int16');   % Trace identification code:
                                        % MATGPR ID code = 7182 (ASCII codes
                                        % [71 82] == 'GR' for GeoRadar) 
SUHDR.nvs     = fread(fid,1,'int16');   % Number of vertically summed traces 
                                        % (see vscode in reel header structure)
SUHDR.nhs     = fread(fid,1,'int16');   % Number of horizontally summed traces 
                                        % (see vscode in reel header structure) 
SUHDR.duse    = fread(fid,1,'int16');   % Data use: 1  = production 2  = test 

% UNSUSED BY ZOND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUHDR.offset  = fread(fid,1,'int32');   % Distance from TX to RX group
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUHDR.gelev   = fread(fid,1,'int32');   % Receiver group elevation 

% UNSUSED BY ZOND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUHDR.selev   = fread(fid,1,'int32');   % Source elevation from sea level 
SUHDR.sdepth  = fread(fid,1,'int32');   % Source depth below surface 
SUHDR.gdel    = fread(fid,1,'int32');   % Datum elevation at receiver group
SUHDR.sdel    = fread(fid,1,'int32');   % Datum elevation at source 
SUHDR.swdep   = fread(fid,1,'int32');   % Water depth at source 
SUHDR.gwdep   = fread(fid,1,'int32');   % Water depth at receiver group
SUHDR.scalel  = fread(fid,1,'int16');   % Scale factor for previous 7 entries 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        
SUHDR.scalco  = fread(fid,1,'int16');   % Scale factor for next 4 entries 
                                        % with value +/- 10 to the power
                                        % 0, 1, 2, 3, or 4 (if positive, 
                                        % multiply, if negative divide) 
SUHDR.sx      = fread(fid,1,'int32');   % X source coordinate
SUHDR.sy      = fread(fid,1,'int32');   % Y source coordinate 
SUHDR.gx      = fread(fid,1,'int32');   % X group coordinate
SUHDR.gy      = fread(fid,1,'int32');   % Y group coordinate
SUHDR.counit  = fread(fid,1,'int16');   % Coordinate units code: for previous four entries 
                                        % 1  = length (meters or feet) 
                                        % 2  = seconds of arc (the X
                                        % -longitude Y -latitude, positive to 
                                        % the east of Greenwich or north of
                                        % the equator
                                        
% UNSUSED BY ZOND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUHDR.wevel   = fread(fid,1,'int16');   % Weathering velocity
SUHDR.swevel  = fread(fid,1,'int16');   % Subweathering velocity
SUHDR.sut     = fread(fid,1,'int16');   % Uphole time at source
SUHDR.gut     = fread(fid,1,'int16');   % Uphole time at receiver group
SUHDR.sstat   = fread(fid,1,'int16');   % Source static correction
SUHDR.gstat   = fread(fid,1,'int16');   % Group static correction
SUHDR.tstat   = fread(fid,1,'int16');   % Total static applied
SUHDR.laga    = fread(fid,1,'int16');   % Lag time A
SUHDR.lagb    = fread(fid,1,'int16');   % Lag time B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUHDR.delrt   = fread(fid,1,'int16');   % Delay recording time, time in ms between initiation time of energy source 
                                        % and time when recording of data samples begins (for deep water work if 
                                        % recording does not start at zero time)
                                       
% UNSUSED BY ZOND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUHDR.muts    = fread(fid,1,'int16');   % Mute time--start
SUHDR.mute    = fread(fid,1,'int16');   % Mute time--end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUHDR.ns      = fread(fid,1,'uint16');  % Number of samples in this trace
SUHDR.dtsu    = fread(fid,1,'uint16');  % sample interval; in micro-seconds

% UNSUSED BY ZOND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUHDR.gain    = fread(fid,1,'int16');   % Gain type of field instruments code:
SUHDR.igc     = fread(fid,1,'int16');   % Instrument gain constant
SUHDR.igi     = fread(fid,1,'int16');   % Instrument early or initial gain
SUHDR.corrsu  = fread(fid,1,'int16');   % Correlated:  1  = no;  2  = yes 
SUHDR.sfs     = fread(fid,1,'int16');   % Sweep frequency at start
SUHDR.sfe     = fread(fid,1,'int16');   % Sweep frequency at end
SUHDR.slen    = fread(fid,1,'int16');   % Sweep length in ms 
SUHDR.styp    = fread(fid,1,'int16');   % Sweep type code
SUHDR.stas    = fread(fid,1,'int16');   % Sweep trace taper length at start 
SUHDR.stae    = fread(fid,1,'int16');   % Sweep trace taper length at end 
SUHDR.tatyp   = fread(fid,1,'int16');   % Taper type
SUHDR.afilf   = fread(fid,1,'int16');   % Alias filter frequency if used
SUHDR.afils   = fread(fid,1,'int16');   % Alias filter slope 
SUHDR.nofilf  = fread(fid,1,'int16');   % Notch filter frequency if used 
SUHDR.nofils  = fread(fid,1,'int16');   % Notch filter slope 
SUHDR.lcf     = fread(fid,1,'int16');   % Low cut frequency if used
SUHDR.hcf     = fread(fid,1,'int16');   % High cut frequncy if used
SUHDR.lcs     = fread(fid,1,'int16');   % Low cut slope
SUHDR.hcs     = fread(fid,1,'int16');   % High cut slope
SUHDR.year    = fread(fid,1,'int16');   % Year data recorded 
SUHDR.day     = fread(fid,1,'int16');   % Day of year 
SUHDR.hour    = fread(fid,1,'int16');   % Hour of day (24 hour clock)
SUHDR.minute  = fread(fid,1,'int16');   % Minute of hour 
SUHDR.sec     = fread(fid,1,'int16');   % Second of minute 
SUHDR.timbas  = fread(fid,1,'int16');   % Time basis code: 1  = local; 2  = GMT; 3  = other 
SUHDR.trwf    = fread(fid,1,'int16');   % Trace weighting factor
SUHDR.grnors  = fread(fid,1,'int16');   % 
SUHDR.grnofr  = fread(fid,1,'int16');   % 
SUHDR.grnlof  = fread(fid,1,'int16');   % 
SUHDR.gaps    = fread(fid,1,'int16');   % 
SUHDR.otrav   = fread(fid,1,'int16');   % Overtravel taper code
SUHDR.d1      = fread(fid,1,'float32'); % Sample spacing for non-seismic data 
SUHDR.f1      = fread(fid,1,'float32'); % First sample location for non-seismic data 
SUHDR.d2      = fread(fid,1,'float32'); % Sample spacing between traces 
SUHDR.f2      = fread(fid,1,'float32'); % First trace location 
SUHDR.ungpow  = fread(fid,1,'float32'); % Negative of power used for dynamic range compression 
SUHDR.unscale = fread(fid,1,'float32'); % Reciprocal of scaling factor to normalize range
SUHDR.ntr     = fread(fid,1,'int32');   % Number of traces 
SUHDR.mark    = fread(fid,1,'int16');   % Mark selected traces 
SUHDR.shortpad  = fread(fid,1,'int16'); % Alignment padding 
SUHDR.unass   = fread(fid,13,'int16');  % Unassigned
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUHDR.mark    = fread(fid,1,'uint16');   % Marker #

% END FUNCTION READSUHEADER