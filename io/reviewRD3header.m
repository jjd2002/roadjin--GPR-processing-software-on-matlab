function HDR = reviewdrd3header()
%
% REVIEWRD3HEADER : Import and display the .RAD header file associated
% with a .RD3 file of RAMAC georadar (Mala Geophysics). Optionally return
% the header in a data structure HDR. 
%
%Author : Andreas Tzanis, 
%         Department of Geophysics, 
%         University of Athens
%         atzanis@geol.uoa.gr
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
%%% Get data file name 
[Infname, Inpname]= uigetfile('*.rd3; *.RD3',...
    'Give RAMAC data file name', ENVAR.currentworkdir);
if Infname == 0, 
    return, 
end;
%%% Make header and marker file names 
Inhdfname = [Infname(1:findstr(Infname,'.')) 'rad'];
Inmkfname = [Infname(1:findstr(Infname,'.')) 'mkn'];

%%% Now open and read header file
fid  = fopen([Inpname Inhdfname],'r');
% First tell which kind of data the header refers to
HDR.datasource = 'RAMAC';
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
datainfotext = 'RAMAC Header to be displayed here';
datatitlestring =[Inpname Infname ];
rand('seed',sum(100*clock));
%whitebg(HDRfig,[rand rand rand]); 
c1 = 0.5;    c2 = 0.5;    c3 = 0.5; 
whitebg(HDRfig,[c1 c2 c3]); 
set(gca,'visible','off','drawmode','fast');
text1 = text(0.45,1.05,'RAMAC HEADER VIEWER') ;
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
%%%%  Now start decoding the header and generating report
clb=cell(1);
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Samples_per_scan = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Frequency = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Frequency_steps = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Signal_position = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Raw_Signal_position= str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Distance_flag = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Time_flag = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Program_flag = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.External_flag = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Time_interval = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Distance_interval = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Operator = hdline(findstr(hdline,':')+1:length(hdline));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Customer = hdline(findstr(hdline,':')+1:length(hdline));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Site = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Antennas = hdline(findstr(hdline,':')+1:length(hdline));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
%HDR.Antenna_orientation=str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Antenna_separation = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Comment = hdline(findstr(hdline,':')+1:length(hdline));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Timewindow = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Stacks = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Stack_exponent = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Stacking_time = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Last_trace = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
hdline = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Stop_position = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
if ~feof(fid),
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    HDR.System_calibration = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    HDR.Start_position = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    HDR.Short_flag = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    HDR.Intermediate_flag  = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    HDR.Long_flag = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    HDR.Preprocessing = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    HDR.High = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    HDR.Low = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    HDR.Fixed_increment = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    HDR.Fixed_moves_up = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    HDR.Fixed_moves_down = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    HDR.Fixed_position = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    HDR.Wheel_calibration  = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    HDR.Positive_direction = str2num(hdline(findstr(hdline,':')+1:length(hdline)));
end                 % end of file loop
fclose(fid);                           

%%% Print report
fontn = get(0,'FixedWidthFontName');
do = ['set(findobj(''tag'',''headerinfobox''),''String'',clb,''horizontalal'','...
        '''left'',''fontname'',fontn);' ];
er = 'viewrd3header'; 
eval(do,er);   
return

function HDR = initheaderstr()
%%%% Initializes the header structure
HDR.datasource           = [];
HDR.Samples_per_scan     = [];
HDR.Frequency            = [];
HDR.Frequency_steps      = [];
HDR.Signal_position      = [];
HDR.Raw_Signal_position  = [];
HDR.Distance_flag        = [];
HDR.Time_flag            = []; 
HDR.Program_flag         = []; 
HDR.External_flag        = [];
HDR.Time_interval        = [];
HDR.Distance_interval    = [];
HDR.Operator             = [];
HDR.Customer             = []; 
HDR.Site                 = [];
HDR.Antennas             = [];
HDR.Antenna_orientation  = [];
HDR.Antenna_separation   = [];
HDR.Comment              = [];
HDR.Timewindow           = [];
HDR.Stacks               = [];
HDR.Stack_exponent       = [];
HDR.Stacking_time        = [];
HDR.Last_trace           = [];
HDR.Stop_position        = [];
HDR.System_calibration   = [];
HDR.Start_position       = [];
HDR.Short_flag           = [];
HDR.Intermediate_flag    = [];
HDR.Long_flag            = [];
HDR.Preprocessing        = [];
HDR.High                 = [];
HDR.Low                  = [];
HDR.Fixed_increment      = [];
HDR.Fixed_moves_up       = [];
HDR.Fixed_moves_down     = [];
HDR.Fixed_position       = [];
HDR.Wheel_calibration    = [];
HDR.Positive_direction   = [];
return