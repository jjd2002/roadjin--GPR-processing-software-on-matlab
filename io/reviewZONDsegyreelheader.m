function HDR = reviewZONDsegyreelheader()
%
% REVIEWSEGYREELHEADER: Import and display the reel header of a SEG-Y 
% data file. Optionally return the header in a data structure HDR. 
%
% Author : Andreas Tzanis, 
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
SEGYdataformat= initSEGYdataformats;

%%% Get data file name 
[Infname, Inpname]= uigetfile('*.segy; *.sgy; *.SEGY; *.SGY',...
    'Give SEG-Y data file name', ENVAR.currentworkdir);
if Infname == 0, 
    erh = errordlg('Error in file specification!',...
        'REVIEWSEGYREELHEADER : ERROR');
    uiwait(erh)
    return, 
end;

%%%%% Prepare to display the header %%%%%%%%%%%%%%%%%%%%
HDRfig = findobj('tag','headerinfofigure');
if ~isempty(HDRfig),
    figure(HDRfig);
    clf;
    matgprwindows('setup');
    matgprwindows('updateadd');
else
    HDRfig = figure('name','Header Information',...
        'tag','headerinfofigure','numbertitle','off', ...
        'position',[500   300   500   600], ...
        'menubar', 'none', ...
        'CreateFcn',['matgprwindows(''setup''); ' ...
                     'matgprwindows(''updateadd'');'],...
        'DeleteFcn',['matgprwindows(''updateremove'', ' ...
                     'findobj(''tag'',''headerinfofigure''));']);
end
datainfotext = 'SEG-Y Header File to be displayed here';
datatitlestring =[Inpname Infname ];
rand('seed',sum(100*clock));
%whitebg(HDRfig,[rand rand rand]); 
c1 = 0.5;    c2 = 0.5;    c3 = 0.5; 
whitebg(HDRfig,[c1 c2 c3]); 
set(gca,'visible','off','drawmode','fast');
text1 = text(0.45,1.05,'SEG-Y HEADER FILE VIEWER') ;
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
%%% Read SEG-Y reel header 

fid=fopen([Inpname Infname],'r',endian);
% Read 3200-byte textual header block (ASCII coded) 
HDR.headertext= fread(fid,3200,'uchar'); 
% Read binary coded block  
HDR.jobid = fread(fid,1,'int32');  % Job identification number  
HDR.lino  = fread(fid,1,'int32');  % Line number (one line per reel)  
HDR.reno  = fread(fid,1,'int32');  % Reel number  
HDR.ntrpr = fread(fid,1,'int16');  % Number of data traces / record  
HDR.nart  = fread(fid,1,'int16');  % Auxiliary traces per record   
HDR.hdt   = fread(fid,1,'uint16'); % Sample interval in microseconds
HDR.dto   = fread(fid,1,'uint16'); % Same for original field record  
HDR.hns   = fread(fid,1,'uint16'); % Number of samples per trace 
HDR.nso   = fread(fid,1,'uint16'); % Same for original field record  
HDR.format= fread(fid,1,'int16');  % Data sample format code 
HDR.fold  = fread(fid,1,'int16');  % CDP fold 
HDR.tshort= fread(fid,1,'int16');  % Trace sorting code: 
HDR.vscode= fread(fid,1,'int16');  % Vertical sum code: 
HDR.hsfs  = fread(fid,1,'int16');  % Sweep frequency at start  
HDR.hsfe  = fread(fid,1,'int16');  % Sweep frequency at end  
HDR.hslen = fread(fid,1,'int16');  % Sweep length (ms)  
HDR.hstyp = fread(fid,1,'int16');  % Sweep type code 
HDR.schn  = fread(fid,1,'int16');  % Trace number of sweep channel 
HDR.hstas = fread(fid,1,'int16');  % Sweep trace taper length at start 
HDR.hstae = fread(fid,1,'int16');  % Sweep trace taper length at end 
HDR.htatyp= fread(fid,1,'int16');  % Sweep trace taper type code 
HDR.hcorr = fread(fid,1,'int16');  % Correlated data traces code 
HDR.bgrcv = fread(fid,1,'int16');  % Binary gain recovered code 
HDR.rcvm  = fread(fid,1,'int16');  % Amplitude recovery method code 
HDR.mfeet = fread(fid,1,'int16');  % Measurement system code: 
HDR.polyt = fread(fid,1,'int16');  % Impulse signal polarity code: 
HDR.vpol  = fread(fid,1,'int16');  % Vibratory polarity code 
HDR.dt    = fread(fid,1,'float32');% Sampling rate (MATGPR specific) 
HDR.t1    = fread(fid,1,'float32');% First sample location (MATGPR specific)
HDR.dx    = fread(fid,1,'float32');% Trace spacing (MATGPR specific)
HDR.x1    = fread(fid,1,'float32');% First trace location (MATGPR specific)
HDR.unass1= fread(fid,112,'int16');% Unassigned 224 bytes
HDR.revision=fread(fid,1,'uint16');% SEG-Y revision number
HDR.fixedlen=fread(fid,1,'int16'); % Fixed length trace flag
HDR.ntexthdrs=fread(fid,1,'int16');% Number of extended textual headers
HDR.unass2 = fread(fid,47,'int16');% Uassigned 94 bytes.
fclose(fid);

%%%% Generate report
clb = cell(1);
txt = ['SEG-Y Revision            : ' num2str(HDR.revision)];
clb(1) = cellstr(txt);
txt = ['Job ID number             : ' num2str(HDR.jobid)];
clb(length(clb)+1) = cellstr(txt);
txt = ['Line number               : ' num2str(HDR.lino)];
clb(length(clb)+1) = cellstr(txt);
txt = ['Reel number               : ' num2str(HDR.reno)];
clb(length(clb)+1) = cellstr(txt);
txt = ['Sampling interval (SEGY)  : ' num2str(HDR.hdt/1000)];
clb(length(clb)+1) = cellstr(txt);
txt = ['Samples per trace         : ' num2str(HDR.hns)];
clb(length(clb)+1) = cellstr(txt);
txt = ['Number of traces/record   : ' num2str(HDR.ntrpr)];
clb(length(clb)+1) = cellstr(txt);
txt = ['Fixed trace length flag   : ' num2str(HDR.fixedlen)];
clb(length(clb)+1) = cellstr(txt);
txt = ['Data Sample Format        : ' SEGYdataformat(HDR.format).name];
clb(length(clb)+1) = cellstr(txt);
txt = ' ';        
clb(length(clb)+1) = cellstr(txt);
txt = 'SEG-Y TEXTUAL HEADER BLOCK';
clb(length(clb)+1) = cellstr(txt);
txt = ' ';        
clb(length(clb)+1) = cellstr(txt);
clb(length(clb)+1)=cellstr( char(HDR.headertext')); 

%%% Print report
fontn = get(0,'FixedWidthFontName');
do = ['set(findobj(''tag'',''headerinfobox''),''String'',clb, ' ...
      '''horizontalal'', ''left'',''fontname'',fontn);' ];
eval(do);   
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
