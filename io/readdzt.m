function [DATA, HDR ] = readdzt()
%
% READDZT : Imports SINGLE channel GSSI (.DZT) georadar data.
%
%   Usage : DATA = readdzt;    [DATA, HEADER] = readdzt;  
%
% RETURNS : A structure DATA containg the GPR data and all other necessary 
%           parameters and, optionally, a structure HEADER containg the 
%           header of the imported DZT data file
%
% REQUIRES: initdatastr.m
%
% Caveat  : At present accepts only NEW style 1024 byte headers. 
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
% Initialize the data ahe header structures.
DATA = initdatastr;

% Get data file name and open
[DATA.fname, DATA.pname]= uigetfile('*.dzt; *.DZT', ...
    'Give DZT data file name',ENVAR.currentworkdir);
if DATA.fname == 0,  
    return
end;

% Open file and proceed
fid=fopen([DATA.pname DATA.fname],'r');

% First tell which kind of data the header refers to
DATA.origin = 'GSSI Mosostatic';
%%% Determine No of headers and abort if > 1
Number_of_headers = 0;
tag = fread(fid,1,'ushort'); 
if tag==hex2dec('00ff'), Number_of_headers = 1; end;  
if tag==hex2dec('01ff'), Number_of_headers = 2; end;  
if tag==hex2dec('02ff'), Number_of_headers = 3; end;  
if tag==hex2dec('03ff'), Number_of_headers = 4; end;  
if Number_of_headers > 1, 
    errordlg(['Number of headers is ' num2str(Number_of_headers) ... 
        '. The program does not handle more than one header'],'READDZT: ERROR')
    fclose(fid); 
    return
end 
%%% Determine Header size and abort if not 1024 bytes
Header_size = fread(fid,1,'ushort');
if Header_size == 512, 
    errordlg(['Header Size is ' num2str(Header_size) ... 
        '. Only NEW STYLE data with 1024 byte headers allowed'],...
        'READDZT : ERROR')
    fclose(fid); 
    return
end

%%% Read the header 
disp('READDZT > Reading header inforrmation')
HDR = readDZTheader(fid);
%%% Assign the DATA structure
DATA.ns = HDR.Samples_per_scan;
%%%% If Scans_per_meter == 0, then data were taken without survey wheel, 
%%%% at equal time intervals - must be transformed to equal spacing using
%%%% marker information. If ~=0, the data traces are equally spaced 
if HDR.Scans_per_meter ~= 0.0,
    DATA.dx = 1/HDR.Scans_per_meter;
end
DATA.sigpos  = HDR.Position;
DATA.Antenna = char(HDR.Antenna_name');
%%%% For GSSI monostatic antennae, separation TxRx =is naught 
DATA.TxRx = 0.0;
%%%% Determine range-gain function
if HDR.Size_of_range_gain > 0,
    drgbreak  = floor(HDR.Samples_per_scan/ (HDR.No_rg_breaks - 1));
    rgbreaks  = [0:1:(HDR.No_rg_breaks-1)]*drgbreak;     
    RangeGain = HDR.Range_gain;
%%%% interpolate from range gain to gain per sample (decibels)
    dummy     = 1 : 1 : HDR.Samples_per_scan; 
    DATA.DZThdgain  = [];
    DATA.DZThdgain  = interp1(rgbreaks(1:2),RangeGain(1:2),...
        dummy(1:rgbreaks(2)),'linear',0);
    for i=2:HDR.No_rg_breaks-1,
        DATA.DZThdgain = [DATA.DZThdgain ...
            interp1(rgbreaks(i:i+1),RangeGain(i:i+1),...
            dummy(rgbreaks(i)+1:rgbreaks(i+1)),'linear',0)];
    end
    if rgbreaks(HDR.No_rg_breaks) < HDR.Samples_per_scan,
        DATA.DZThdgain = [DATA.DZThdgain ...
            zeros(1,HDR.Samples_per_scan-rgbreaks(HDR.No_rg_breaks))+...
            RangeGain(HDR.No_rg_breaks)];
    end
else
    HDR.No_rg_breaks = 0;
    HDR.Range_gain   = [];
    DATA.DZThdgain   = [];
end
%%%% Assign header comments (but first wrap comments to 80 column text)
hc=[]; 
for i=1:ceil(length(HDR.Comment)/80), 
    if i*80 >length(HDR.Comment), 
        hc=strvcat(hc,HDR.Comment((i-1)*80+1:length(HDR.Comment))); 
    else 
        hc=strvcat(hc,HDR.Comment((i-1)*80+1:i*80)); 
    end; 
end;
DATA.comments = str2mat('ORIGINAL FILE HEADER COMMENTS:',hc);

disp('        ==> OK')

%%%% Determine Number of Traces 
disp('READDZT : Reading data')
fseek(fid, 0, 'eof');   
last_byte  = ftell(fid);  
databytes  = last_byte - HDR.Header_size;  
if HDR.Bits_per_word == 8, 
    DATA.ntr  = databytes/HDR.Samples_per_scan;
elseif HDR.Bits_per_word == 16, 
    DATA.ntr  = (databytes/2)/HDR.Samples_per_scan;
elseif Bits_perword == 32, 
    DATA.ntr  = (databytes/4)/HDR.Samples_per_scan;
elseif HDR.Bits_per_word == 64, 
    DATA.ntr  = (databytes/8)/HDR.Samples_per_scan;
end
%%%% Check the type of the data  
    if HDR.Bits_per_word == 8,   datatype = 'uint8';            % 'uchar'
elseif HDR.Bits_per_word ==16,   datatype = 'uint16';           % 'ushort'
elseif HDR.Bits_per_word ==32,   datatype = 'uint32';           % 'uint'
   end
%%%% Import the data  
fseek(fid, HDR.Header_size, 'bof');       % move to beginning of the file 
datasize=[HDR.Samples_per_scan,inf];      % and skip header                           
DATA.d=fread(fid,datasize,datatype);
disp('        ==> OK')
fclose(fid);                              % close the data file

%%%% Samplig interval and two-way traveltime
DATA.dt     = HDR.Range / (HDR.Samples_per_scan - 1);
DATA.tt2w   = 0 : DATA.dt : DATA.dt*(HDR.Samples_per_scan - 1);
DATA.zlab   = 'Traveltime (ns)';

%%%% Scan Line axis
if HDR.Scans_per_meter == 0,
    DATA.x    = 1 : 1 : DATA.ntr;
    DATA.xlab = 'Scan Axis (# Traces)';
else 
    DATA.x    = 0 : DATA.dx : (DATA.ntr - 1)*DATA.dx;
    DATA.xlab = 'Scan Axis (meters)';
end

%%%% The first two samples of each trace contain info on the health status
%%%% of traces and marker locations 
first2samples = DATA.d(1:2,:);                                  
%%%% Check the health status of traces - first sample should have all 
%%%% bits set                                                           
disp(['READDZT > Checking traces of file ' DATA.fname ' for CONSISTENCY'])
if HDR.Bits_per_word == 8, 
    temp = find(first2samples(1,:)==hex2dec('FF') | ...
        first2samples(1,:)==hex2dec('80') | ...
        first2samples(1,:)==hex2dec('7F'));    
elseif HDR.Bits_per_word == 16, 
    temp = find(first2samples(1,:)==hex2dec('FFFF') | ...
        first2samples(1,:)==hex2dec('8000') | ...
        first2samples(1,:)==hex2dec('7FFF'));
end
if length(first2samples(1,:)) ~= length(temp), 
    disp(['        ==> Possible inconsistency : ' ...
     'All bits of the first scan sample have not been set in all traces']), 
else
    disp('        ==> OK')
end

%%%% Check for MARKER traces in a dumb brute force manner !!!
disp('READDZT > Checking for MARKER TRACES')
if isempty(DATA.markertr),
    smarker = hex2dec('64');    
    DATA.markertr = find(first2samples(2,:)==smarker); 
end;
if isempty(DATA.markertr),
    smarker = hex2dec('E8');    
    DATA.markertr = find(first2samples(2,:)==smarker); 
end;
if isempty(DATA.markertr),
    smarker = hex2dec('E800');  
    DATA.markertr = find(first2samples(2,:)==smarker); 
end;
if isempty(DATA.markertr),
    smarker = hex2dec('E1');    
    DATA.markertr = find(first2samples(2,:)==smarker); 
end;
if isempty(DATA.markertr),
    smarker = hex2dec('E100');  
    DATA.markertr = find(first2samples(2,:)==smarker); 
end;
if isempty(DATA.markertr),
    smarker = hex2dec('EC');    
    DATA.markertr = find(first2samples(2,:)==smarker); 
end;
if isempty(DATA.markertr),
    smarker = hex2dec('EC00');  
    DATA.markertr = find(first2samples(2,:)==smarker); 
end;
if isempty(DATA.markertr),
    smarker = hex2dec('F1');    
    DATA.markertr = find(first2samples(2,:)==smarker); 
end;
if isempty(DATA.markertr),
    smarker = hex2dec('F100');  
    DATA.markertr = find(first2samples(2,:)==smarker); 
end;
if isempty(DATA.markertr),
    disp('        ==> THERE ARE NO MARKER TRACES'),
else
    disp(['        ==> FOUND ' num2str(length(DATA.markertr)) ...
        ' MARKER TRACES'])
end;
%%% Make sure this is a column vector
DATA.markertr = DATA.markertr(:);

%%% The first two samples of each trace are used for consistency testing
%%% and marker information - they can be discarded ...
%%% but only if the signal position has not been adjusted before ...
if DATA.sigpos ~= 0,
    DATA.d         = DATA.d(3:DATA.ns,:);
    DATA.tt2w      = 0 : DATA.dt : DATA.dt*(DATA.ns - 3);
    if ~isempty(DATA.DZThdgain),
        DATA.DZThdgain = DATA.DZThdgain(3:DATA.ns);
    end
    DATA.ns        = size(DATA.d,1);
    DATA.sigpos    = DATA.sigpos - 2*DATA.dt;
end

%%%% In a final step, convert unsigned integers to signed 
for i=1:DATA.ntr
    if HDR.Bits_per_word == 8,   
        DATA.d(:,i) = DATA.d(:,i) - 128; 
    elseif HDR.Bits_per_word ==16,  
        DATA.d(:,i) = DATA.d(:,i) - 32768; 
    elseif HDR.Bits_per_word ==32,   
        DATA.d(:,i) = DATA.d(:,i) - 2147483647; 
    end
end
%%% And that's all folks %%% And that's all folks ...
return

function HDR = readDZTheader(fid)
%%%% Returns a DZT data file header in structure HDR
fseek(fid,0,'bof');
HDR.tag                  = fread(fid,1,'ushort');                
HDR.Header_size          = fread(fid,1,'ushort');
HDR.Samples_per_scan     = fread(fid,1,'ushort');      
HDR.Bits_per_word        = fread(fid,1,'ushort');
HDR.Binary_offset        = fread(fid,1,'short');
HDR.Scans_per_second     = fread(fid,1,'float');
HDR.Scans_per_meter      = fread(fid,1,'float');
HDR.Meters_per_mark      = fread(fid,1,'float');
HDR.Position             = fread(fid,1,'float');
HDR.Range                = fread(fid,1,'float');
HDR.Scans_per_pass       = fread(fid,1,'ushort');
HDR.CreateDate.sec       = fread(fid,1,'ubit5')*2;     
HDR.CreateDate.min       = fread(fid,1,'ubit6');
HDR.CreateDate.hour      = fread(fid,1,'ubit5');
HDR.CreateDate.day       = fread(fid,1,'ubit5');
HDR.CreateDate.month     = fread(fid,1,'ubit4');
HDR.CreateDate.year      = fread(fid,1,'ubit7')+1980;
HDR.ModifyDate.sec       = fread(fid,1,'ubit5');
HDR.ModifyDate.min       = fread(fid,1,'ubit6');
HDR.ModifyDate.hour      = fread(fid,1,'ubit5');
HDR.ModifyDate.day       = fread(fid,1,'ubit5');
HDR.ModifyDate.month     = fread(fid,1,'ubit4');
HDR.ModifyDate.year      = fread(fid,1,'ubit7');
if HDR.ModifyDate.year ~= 0, 
    HDR.ModifyDate.year  = HDR.ModifyDate.year + 1980; 
end
HDR.Offset_to_range_gain = fread(fid,1,'ushort');
HDR.Size_of_range_gain   = fread(fid,1,'ushort');
HDR.Offset_to_text       = fread(fid,1,'ushort');
HDR.Size_of_text         = fread(fid,1,'ushort');
HDR.Offset_to_proc_hist  = fread(fid,1,'ushort');
HDR.Size_of_proc_hist    = fread(fid,1,'ushort');
HDR.Number_of_channels   = fread(fid,1,'ushort');
HDR.Dielectric_constant  = fread(fid,1,'float');
HDR.Top_pos_in_m         = fread(fid,1,'float');
HDR.Range_in_m           = fread(fid,1,'float');
HDR.reserved             = fread(fid,31,'uchar');
HDR.Data_type            = fread(fid,1,'uchar');
HDR.Antenna_name         = fread(fid,14,'uchar'); 
HDR.Channel_mask         = fread(fid,1,'ushort');
HDR.This_file_name       = fread(fid,12,'uchar'); 
HDR.Checksum             = fread(fid,1,'ushort');
%%%% Read header gain data 
if HDR.Size_of_range_gain > 0,
    HDR.No_rg_breaks     = fread(fid,1,'ushort'); % uint16
    HDR.Range_gain       = fread(fid,(HDR.Size_of_range_gain-2)/4,'float')'; 
end
%%%% Read comments 
fseek(fid,HDR.Offset_to_text,'bof');      
if HDR.Size_of_text > 0, 
    HDR.Comment = fread(fid,HDR.Size_of_text,'char'); 
    if size(HDR.Comment == 1)
        HDR.Comment = char(transpose(HDR.Comment));
        ii = find(HDR.Comment < 32);
        HDR.Comment(ii) = char(32);
    end
else
    HDR.Comment = ' ';
end
%%%% Read Processing history 
if HDR.Size_of_proc_hist == 0, 
    HDR.Proc_hist = [];
    return
end
fseek(fid,HDR.Offset_to_proc_hist,'bof');
current_byte = ftell(fid) - HDR.Offset_to_proc_hist;
while HDR.Offset_to_proc_hist + current_byte < HDR.Offset_to_proc_hist + ...
        HDR.Size_of_proc_hist,
    ph = fread(fid,2,'char');
    switch ph(1)
        case 1                                               % PR_VIIRL
            HDR.Proc_hist.viirl  = fread(fid,1,'float');
        case 2                                               % PR_VIIRH
            HDR.Proc_hist.viirh  = fread(fid,1,'float');
        case 3                                                % PR_VTCL
            HDR.Proc_hist.vtcl   = fread(fid,1,'float');
        case 4                                               % PR_VTCH
            HDR.Proc_hist.vtch   = fread(fid,1,'float');
        case 5                                               % PR_VBOXL
            HDR.Proc_hist.vboxl  = fread(fid,1,'float');                   
        case 6                                               % PR_VBOXH    
            HDR.Proc_hist.vboxh  = fread(fid,1,'float');
        case 7                                               % PR_VTRIL
            HDR.Proc_hist.vtril  = fread(fid,1,'float');
        case 8                                               % PR_VTRIH    
            HDR.Proc_hist.vtrih  = fread(fid,1,'float');
        case 9                                               % PR_VFIRL 
            HDR.Proc_hist.vfirl  = fread(fid,1,'float');
        case 10                                             % PR_VFIRH 
            HDR.Proc_hist.vfirh  = fread(fid,1,'float');
        case 11                                              % PR_HIIRL 
            HDR.Proc_hist.hiirl  = fread(fid,1,'float');
        case 12                                              % PR_HIIRH
            HDR.Proc_hist.hiirh  = fread(fid,1,'float');
        case 13                                               % PR_HTCL
            HDR.Proc_hist.htcl   = fread(fid,1,'float');
        case 14                                               % PR_HTCH 
            HDR.Proc_hist.htch   = fread(fid,1,'float');
        case 15                                              % PR_HBOXL 
            HDR.Proc_hist.hboxl  = fread(fid,1,'float');
        case 16                                              % PR_HBOXH
            HDR.Proc_hist.hboxh  = fread(fid,1,'float');
        case 17                                              % PR_HTRIL
            HDR.Proc_hist.htril  = fread(fid,1,'float');
        case 18                                              % PR_HTRIH
            HDR.Proc_hist.htrih  = fread(fid,1,'float');
        case 19                                              % PR_HFIRL
            HDR.Proc_hist.hfirl  = fread(fid,1,'float');
        case 20                                              % PR_HFIRH
            HDR.Proc_hist.hfirh  = fread(fid,1,'float');
        case 21                                                % PR_LPF
            HDR.Proc_hist.lp     = fread(fid,1,'float');
        case 22                                                % PR_HPF
            HDR.Proc_hist.hp     = fread(fid,1,'float');
        case 23                                                % PR_STS
            HDR.Proc_hist.sts    = fread(fid,1,'float');
    end        % switch loop
    current_byte = ftell(fid) - HDR.Offset_to_proc_hist;
end        % while loop
return
