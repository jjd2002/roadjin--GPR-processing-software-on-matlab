function HDR = reviewdztheader()
%
% REVIEWDZTHEADER : Import single channel GSSI (.DZT) georadar data and
% return a structure HDR containing the file header  
%
% Caveat  : At present accepts only NEW style 1024 byte headers. 
%
%  Author : Andreas Tzanis
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
%%% Initialize header structure
HDR = initheaderstr;
%%%% Get DZT data file name
[Infname, Inpname]= uigetfile('*.dzt','Give DZT data file name', ...
    ENVAR.currentworkdir);
if isequal(Infname,0) | isequal(Inpname,0), 
    return
end;
%%% Open data file and check for compliance with matGPR
fid=fopen([Inpname Infname],'r');
HDR.tag = fread(fid,1,'ushort');             % Determine No of headers
if HDR.tag==hex2dec('00ff'), Number_of_headers = 1; end;  
if HDR.tag==hex2dec('01ff'), Number_of_headers = 2; end;  
if HDR.tag==hex2dec('02ff'), Number_of_headers = 3; end;  
if HDR.tag==hex2dec('03ff'), Number_of_headers = 4; end;  
if Number_of_headers > 1, 
    errordlg(['Number of headers is ' num2str(Number_of_headers) ... 
        '. Cannot handle more than one header'],'ERROR')
    fclose(fid); 
    HDR = initheaderstr;
    return
end 
HDR.Header_size = fread(fid,1,'ushort');
if HDR.Header_size == 512, 
    errordlg(['Header Size is ' num2str(HDR.Header_size) ... 
        '. Only NEW STYLE data headers allowed'],...
        'ERROR')
    fclose(fid); 
    HDR = initheaderstr;
    return
end
%%%%% Since all is OK, pprepare to display the header %%%%%%%%%%%%%%%%%%%%
HDRfig = findobj('tag','headerinfofigure');
if ~isempty(HDRfig),
    figure(HDRfig);
    clf;
    matgprwindows('setup'); 
    matgprwindows('updateadd');
else
    HDRfig = figure('name','Header Information',...
        'tag','headerinfofigure', ...
        'numbertitle','off', ...
        'position',[500   200   500   500], ...
        'menubar', 'none', ...
        'CreateFcn',['matgprwindows(''setup''); ' ...
                     'matgprwindows(''updateadd'');'],...
        'DeleteFcn',['matgprwindows(''updateremove'', ' ...
                     'findobj(''tag'',''headerinfofigure''));']);
end
datainfotext = 'DZT Header to be displayed here';
datatitlestring =[Inpname Infname ];
rand('seed',sum(100*clock));
%whitebg(HDRfig,[rand rand rand]); 
c1 = 0.5;    c2 = 0.5;    c3 = 0.5; 
whitebg(HDRfig,[c1 c2 c3]); 
set(gca,'visible','off','drawmode','fast');
text1 = text(0.45,1.05,'DZT HEADER VIEWER') ;
set(text1,'FontSize',13,'Color','k','FontWeight','bold','horizontalal','center')
% prepare to display header info on a GPR data file
top=0.9;
left=0.05;
right=0.95;
bottom=0.05;
labelheight=0.07;
spacing=0.02;
% Draw the text window frame
frameBorder=0.02;
framePosition=[left-frameBorder bottom-frameBorder ...
        (right-left)+2*frameBorder (top-bottom)+2*frameBorder];
uicontrol( 'Style','frame', ...
    'Units','normalized', 'Position',framePosition, ...
    'BackgroundColor',[0.0 0.5 0.5]);
% Draw the text label
labelPosition=[left top-labelheight (right-left) labelheight];
uicontrol( 'Style','text', ...
    'Units','normalized', 'Position',labelPosition, ...
    'BackgroundColor',[0.0 0.5 0.5], 'ForegroundColor',[1 1 1], ...
    'String',datatitlestring, 'fontsize', 10, 'fontweight', 'demi');
% Display the info box and display the info text
textPosition=[left bottom (right-left) top-bottom-labelheight-2*spacing];
uicontrol( 'Style','edit', 'tag', 'headerinfobox', ...
    'Units','normalized', 'Max',18, 'String',datainfotext, ...
    'BackgroundColor',[1 1 1], ...
    'Position',textPosition);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Now start decoding the header 
HDR.Samples_per_scan     = fread(fid,1,'ushort');  % Samples per trace
HDR.Bits_per_word        = fread(fid,1,'ushort');
HDR.Binary_offset        = fread(fid,1,'short');
HDR.Scans_per_second     = fread(fid,1,'float');
HDR.Scans_per_meter      = fread(fid,1,'float');
%%%% If Scans_per_meter == 0, then data were taken without survey wheel, 
%%%% at equal time intervals - irregular trace spacing. If ~=0, the data 
%%%% traces are equally spaced 
if HDR.Scans_per_meter ~= 0.0,
    HDR.dx = 1/HDR.Scans_per_meter;
end
HDR.Meters_per_mark      = fread(fid,1,'float');
HDR.Position             = fread(fid,1,'float');
HDR.Range                = fread(fid,1,'float');
HDR.Scans_per_pass       = fread(fid,1,'ushort');
HDR.CreateDate.sec       = fread(fid,1,'ubit5')*2;            % Date info 
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
    HDR.ModifyDate.year = HDR.ModifyDate.year + 1980; 
end
Offset_to_range_gain     = fread(fid,1,'ushort');            % Header info 
Size_of_range_gain       = fread(fid,1,'ushort');
Offset_to_text           = fread(fid,1,'ushort');
Size_of_text             = fread(fid,1,'ushort');
Offset_to_proc_hist      = fread(fid,1,'ushort');
Size_of_proc_hist        = fread(fid,1,'ushort');
HDR.Number_of_channels   = fread(fid,1,'ushort');
HDR.Dielectric_constant  = fread(fid,1,'float');
HDR.Top_pos_in_m         = fread(fid,1,'float');
HDR.Range_in_m           = fread(fid,1,'float');
Reserved                 = fread(fid,31,'uchar');
HDR.Data_type            = fread(fid,1,'uchar');
HDR.Antenna_name         = fread(fid,14,'uchar'); 
%%%% For GSSI monostatic instruments antenna separation (TxRx) is naught 
HDR.Channel_mask         = fread(fid,1,'ushort');
HDR.This_file_name       = fread(fid,12,'uchar'); 
HDR.Checksum             = fread(fid,1,'ushort');
%%%% Read gain data and determine range-gain function
if Size_of_range_gain > 0,
    HDR.No_rg_breaks     = fread(fid,1,'ushort');                 % uint16
    Rg_break_delta       = floor(HDR.Samples_per_scan/...
                           (HDR.No_rg_breaks - 1));
    HDR.Rg_breaks        = [0:1:(HDR.No_rg_breaks-1)]*Rg_break_delta;     
    HDR.Range_gain       = fread(fid,(Size_of_range_gain-2)/4,'float')'; 
end
%%%% Samplig interval 
HDR.dt     = HDR.Range / (HDR.Samples_per_scan - 1);
%%%% Determine Number of Traces 
fseek(fid, 0, 'eof');   
last_byte  = ftell(fid);  
databytes  = last_byte - HDR.Header_size;  
if HDR.Bits_per_word == 8, 
    HDR.ntr  = databytes/HDR.Samples_per_scan;
elseif HDR.Bits_per_word == 16, 
    HDR.ntr  = (databytes/2)/HDR.Samples_per_scan;
elseif Bits_perword == 32, 
    HDR.ntr  = (databytes/4)/HDR.Samples_per_scan;
elseif HDR.Bits_per_word == 64, 
    HTR.ntr  = (databytes/8)/HDR.Samples_per_scan;
end
%%%%%% Generate basic report %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clb=cell(1);
txt = ['Bits per word          : ' num2str(HDR.Bits_per_word)];
clb(1) = cellstr(txt);
txt = ['Scans per second       : ' num2str(floor(HDR.Scans_per_second))];
clb(length(clb)+1) = cellstr(txt);
txt = ['Scans per pass         : ' num2str(HDR.Scans_per_pass)];
clb(length(clb)+1) = cellstr(txt);
txt = ['Samples per trace      : ' num2str(HDR.Samples_per_scan)];
clb(length(clb)+1) = cellstr(txt);
txt = ['Sampling interval      : ' num2str(HDR.dt) ' ns'];
clb(length(clb)+1) = cellstr(txt);
txt = ['Time window (range)    : ' num2str(HDR.Range) ' ns'];
clb(length(clb)+1) = cellstr(txt);
txt = ['Signal Position        : ' num2str(HDR.Position) ' ns'];
clb(length(clb)+1) = cellstr(txt);
txt = ['Number of traces       : ' num2str(HDR.ntr)];
clb(length(clb)+1) = cellstr(txt);
if HDR.Scans_per_meter ~=0, 
    txt = ['Traces per meter       : ' num2str(HDR.Scans_per_meter) ' m^-^1'];
    clb(length(clb)+1) = cellstr(txt);
    txt = ['Trace spacing          : ' num2str(HDR.dx) ' m'];
    clb(length(clb)+1) = cellstr(txt);
else
    txt = 'Traces per meter       : Irregular';
    clb(length(clb)+1) = cellstr(txt);
    txt = 'Trace spacing          : Unequally spaced';
    clb(length(clb)+1) = cellstr(txt);
end
txt = ['Meters per mark        : ' num2str(HDR.Meters_per_mark)];
clb(length(clb)+1) = cellstr(txt);
txt = ['Dielectric Constant    : ' num2str(HDR.Dielectric_constant)];
clb(length(clb)+1) = cellstr(txt);
txt = ['Range in meters        : ' num2str(HDR.Range_in_m) ' m'];
clb(length(clb)+1) = cellstr(txt);
txt = ['Creation Date          : '  num2str(HDR.CreateDate.day) ...
   '-' num2str(HDR.CreateDate.month) '-' num2str(HDR.CreateDate.year)...
   ', ' num2str(HDR.CreateDate.hour) ':' num2str(HDR.CreateDate.min) ':'...
   num2str(HDR.CreateDate.sec)];
clb(length(clb)+1) = cellstr(txt);
if HDR.ModifyDate.year~=0,
    txt = ['Modification Date      : '  num2str(HDR.ModifyDate.day) ...
             '-' num2str(HDR.ModifyDate.month) '-' num2str(HDR.ModifyDate.year)];
    clb(length(clb)+1) = cellstr(txt);
else
    txt = 'Modification Date       : NOT MODIFIED';
    clb(length(clb)+1) = cellstr(txt);
end
txt = ['Antenna                : ' char(HDR.Antenna_name')];
clb(length(clb)+1) = cellstr(txt);
txt = ['File name              : ' char(HDR.This_file_name')];
clb(length(clb)+1) = cellstr(txt);
txt = ['Checksum               : ' num2str(HDR.Checksum)];
clb(length(clb)+1) = cellstr(txt);
if Size_of_range_gain > 0,
    txt = 'Range Gain             : Sample            db';
    clb(length(clb)+1) = cellstr(txt);
    for i=1:HDR.No_rg_breaks,
        txt = sprintf('                           %-3d           %-7.3f\n',...
            [HDR.Rg_breaks(i) HDR.Range_gain(i)]);
        clb(length(clb)+1) = cellstr(txt);
    end
else
    txt = 'Range Gain             : NOT ASSIGNED';
    clb(length(clb)+1) = cellstr(txt);
end
%%%% OK, now import comments and processing history 
fseek(fid,Offset_to_text,'bof');      
if Size_of_text > 0, 
    txt = 'Comments               : ';
    clb(length(clb)+1) = cellstr(txt);
    HDR.Comment = fread(fid,Size_of_text,'char'); 
    txt = char(HDR.Comment');
    clb(length(clb)+1) = cellstr(txt);
else
    txt = 'Comments               : NO COMMENTS';   
    clb(length(clb)+1) = cellstr(txt);
end                
%%%% Read Processing history 
txt = 'Processing history     : ';
    clb(length(clb)+1) = cellstr(txt);
if Size_of_proc_hist ==0, 
    txt = sprintf('\t\t\t\t\t\t  NO PROCESSING HISTORY');
    clb(length(clb)+1) = cellstr(txt);
else
    fseek(fid,Offset_to_proc_hist,'bof');
    current_byte = ftell(fid) - Offset_to_proc_hist;
    while Offset_to_proc_hist + current_byte < Offset_to_proc_hist + Size_of_proc_hist,
        ph = fread(fid,2,'char');
        switch ph(1)
            case 1                                               % PR_VIIRL
                HDR.Proc_hist.viirl  = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tVert. IIR low pass  N=%d F=%g', ...
                    ph(2),HDR.Proc_hist.viirl));
            case 2                                               % PR_VIIRH
                HDR.Proc_hist.viirh  = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tVert. IIR high pass N=%d F=%g',...
                    ph(2),HDR.Proc_hist.viirh));
            case 3                                                % PR_VTCL
                HDR.Proc_hist.vtcl   = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tVert. low pass  TC=%g',...
                    HDR.Proc_hist.vtcl));
            case 4                                               % PR_VTCH
                HDR.Proc_hist.vtch   = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tVert. high pass TC=%g',...
                    HDR.Proc_hist.vtch));
            case 5                                               % PR_VBOXL
                HDR.Proc_hist.vboxl  = fread(fid,1,'float');                   
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tVert. boxcar low pass    N=%d',...
                    HDR.Proc_hist.vboxl));
            case 6                                               % PR_VBOXH    
                HDR.Proc_hist.vboxh  = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tVert. boxcar high pass   N=%d',...
                    HDR.Proc_hist.vboxh));
            case 7                                               % PR_VTRIL
                HDR.Proc_hist.vtril  = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tVert. triangle low pass  N=%d',...
                    HDR.Proc_hist.vtril));
            case 8                                               % PR_VTRIH    
                HDR.Proc_hist.vtrih  = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tVert. triangle high pass N=%d',...
                    HDR.Proc_hist.vtrih));
            case 9                                               % PR_VFIRL 
                HDR.Proc_hist.vfirl  = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tVert. FIR low pass  N=%d F=%g',...
                    ph(2),HDR.Proc_hist.vfirl));
            case 10                                             % PR_VFIRH 
                HDR.Proc_hist.vfirh  = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tVert. FIR high pass N=%d F=%g',...
                    ph(2),HDR.Proc_hist.vfirh));
            case 11                                              % PR_HIIRL 
                HDR.Proc_hist.hiirl  = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tHoriz. IIR low pass  N=%d F=%g',...
                    ph(2),HDR.Proc_hist.hiirl));
            case 12                                              % PR_HIIRH
                HDR.Proc_hist.hiirh  = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tHoriz. IIR high pass N=%d F=%g',...
                    ph(2),HDR.Proc_hist.hiirh));
            case 13                                               % PR_HTCL
                HDR.Proc_hist.htcl   = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tHoriz. IIR stacking TC=%g',...
                    HDR.Proc_hist.htcl));
            case 14                                               % PR_HTCH 
                HDR.Proc_hist.htch   = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tHoriz. IIR background removal TC=%g',...
                    HDR.Proc_hist.htch));
            case 15                                              % PR_HBOXL 
                HDR.Proc_hist.hboxl  = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tHoriz. boxcar stacking   N=%d',...
                    HDR.Proc_hist.hboxl));
            case 16                                              % PR_HBOXH
                HDR.Proc_hist.hboxh  = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tHoriz. boxcar background removal N=%d',...
                    HDR.Proc_hist.hboxh));
            case 17                                              % PR_HTRIL
                HDR.Proc_hist.htril  = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tHoriz. triangle stacking N=%d',...
                    HDR.Proc_hist.htril));
            case 18                                              % PR_HTRIH
                HDR.Proc_hist.htrih  = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tHoriz. triangle background removal N=%d',...
                    HDR.Proc_hist.htrih));
            case 19                                              % PR_HFIRL
                HDR.Proc_hist.hfirl  = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tHoriz. FIR low pass  N=%d F=%g',...
                    ph(2),HDR.Proc_hist.hfirl));
            case 20                                              % PR_HFIRH
                HDR.Proc_hist.hfirh  = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tHoriz. FIR high pass N=%d F=%g',...
                ph(2),HDR.Proc_hist.hfirh));
            case 21                                                % PR_LPF
                HDR.Proc_hist.lp     = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tOld style low pass filter N=%d F=%g',...
                    ph(2),HDR.Proc_hist.lp));
            case 22                                                % PR_HPF
                HDR.Proc_hist.hp     = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tOld style high pass filter N=%d F=%g',...
                    ph(2),HDR.Proc_hist.hp));
            case 23                                                % PR_STS
                HDR.Proc_hist.sts    = fread(fid,1,'float');
                clb(length(clb)+1) = ...
                    cellstr(sprintf('\tStatic stacking N=%d',...
                    HDR.Proc_hist.sts));
        end        % switch loop
        current_byte = ftell(fid) - Offset_to_proc_hist;
    end        % while loop
end        % if loop
fclose(fid);                       % Nothing more to do - close data file

%%% Print report
fontn = get(0,'FixedWidthFontName');
do = ['set(findobj(''tag'',''headerinfobox''),''String'',clb,' ...
      '''horizontalal'','...
      '''left'',''fontname'',fontn);' ];
er = 'reviewDZTheader'; 
eval(do,er);   
return

function HDR = initheaderstr()
%%%% Initializes the header structure
HDR.datasource           = [];
HDR.tag                  = [];           % No of headers in the data file
HDR.Header_size          = []; 
HDR.Samples_per_scan     = [];           % Data sampling rate
HDR.Bits_per_word        = [];
HDR.Binary_offset        = [];
HDR.Scans_per_second     = [];
HDR.Scans_per_meter      = [];
HDR.dx                   = 0;
HDR.ntr                  = [];
HDR.Meters_per_mark      = [];
HDR.Position             = [];           % signal position
HDR.Range                = [];
HDR.Scans_per_pass       = [];
HDR.CreateDate.sec       = [];            % Date info 
HDR.CreateDate.min       = [];
HDR.CreateDate.hour      = [];
HDR.CreateDate.day       = [];
HDR.CreateDate.month     = [];
HDR.CreateDate.year      = [];
HDR.ModifyDate.sec       = [];
HDR.ModifyDate.min       = [];
HDR.ModifyDate.hour      = [];
HDR.ModifyDate.day       = [];
HDR.ModifyDate.month     = [];
HDR.ModifyDate.year      = [];
HDR.Number_of_channels   = [];
HDR.Dielectric_constant  = [];
HDR.Top_pos_in_m         = [];
HDR.Range_in_m           = [];
HDR.Antenna_name         = [];
HDR.Channel_mask         = [];
HDR.This_file_name       = [];
HDR.Checksum             = [];
HDR.No_rg_breaks         = [];
HDR.Rg_breaks            = [];
HDR.Range_gain           = [];
HDR.Comment              = [];
HDR.Proc_hist            = [];
return
