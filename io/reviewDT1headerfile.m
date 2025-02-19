function HDR = reviewDT1headerfile()
%
% REVIEWDT1HEADERFILE: Import and display the .HD header file associated
% with a .DT1 file of PULSE EKKO georadar data. Optionally return the
% header in a data structure HDR. 
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

global ENVAR
%%% Get data file name 
[Infname, Inpname]= uigetfile('*.dt1; *.DT1', ...
    'Give PULSE EKKO data file name', ENVAR.currentworkdir);
if Infname == 0, 
    return, 
end;
%%% Make header file name 
Inhdfname = [Infname(1:findstr(Infname,'.')) 'hd'];

%%% Now open and read header file
fid  = fopen([Inpname Inhdfname],'r');
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
datainfotext = 'PULSE EKKO Header File to be displayed here';
datatitlestring =[Inpname Infname ];
rand('seed',sum(100*clock));
%whitebg(HDRfig,[rand rand rand]); 
c1 = 0.5;    c2 = 0.5;    c3 = 0.5; 
whitebg(HDRfig,[c1 c2 c3]); 
set(gca,'visible','off','drawmode','fast');
text1 = text(0.45,1.05,'PULSE EKKO HEADER FILE VIEWER') ;
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
HDR.comment = ' ';
for i=1:5,
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    hdline = [hdline '          '];
    if ~(strcmp(hdline(1:6),'NUMBER')), 
        HDR.comment = cat(2,HDR.comment, ' ', hdline);
    else
        break
    end
end;
HDR.Num_traces = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
hdline  = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Samples_per_scan = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
hdline  = fgetl(fid);
clb(length(clb)+1)= cellstr(hdline);
HDR.Signal_position = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
hdline  = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Timewindow = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
hdline  = fgetl(fid);
clb(length(clb)+1)= cellstr(hdline);
HDR.Starting_position= str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
hdline  = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Final_position = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
hdline  = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Step_size_used = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
hdline  = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.position_units = hdline(findstr(hdline,'=')+1:length(hdline));
hdline  = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Nominal_Frequency = hdline(findstr(hdline,'=')+1:length(hdline));
hdline  = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Antenna_separation= str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
hdline  = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Pulser_voltage = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
hdline  = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Number_of_stacks = str2num(hdline(findstr(hdline,'=')+1:length(hdline)));
hdline  = fgetl(fid);
clb(length(clb)+1) = cellstr(hdline);
HDR.Survey_mode = hdline(findstr(hdline,'=')+1:length(hdline));
while ~feof(fid),
    hdline = fgetl(fid);
    clb(length(clb)+1) = cellstr(hdline);
    HDR.comment = cat(2, HDR.comment, hdline);
end
fclose(fid);

%%% Print report
fontn = get(0,'FixedWidthFontName');
do = ['set(findobj(''tag'',''headerinfobox''),''String'',clb,''horizontalal'','...
        '''left'',''fontname'',fontn);' ];
eval(do);   
return
