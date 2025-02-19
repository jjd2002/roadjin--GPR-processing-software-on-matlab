function [d_out,x_out,ntr_out,markertr_out,xyz_out] = ...
    editscanaxis(d,x,dx,markertr,xyz)
%
% EDITSCANAXIS : Reduces the size of the GPR data array by discarding
% groups of traces at the beginning or at the end of the section抯 scan
% axis (scanline). It is also used to extract groups of traces or even
% large parts of the section to a separate data set. The routine
% initializes a GUI, in which the user specifies the desired operation
% (e.g. trim, or extract, input mode). Then, it will process the data via 
% the function handle callback "edt". The reduced data set returns from
% "edt" by means of interim global variables that are kept strictly within
% the scope of "editscanaxis". The interim variables are copied onto the
% output variables just before "editscanaxis" returns.  
%
% USAGE : [d_out,x_out,ntr_out,markertr_out,xyz_out] = ...
%                        editscanaxis(d, x, dx, markertr,xyz)
%
%   Inputs :
% d        : The [ns x ntr] GPR data array
% x        : The [1 x ntr] scan axis coordinate vector
% dx       : The trace spacing in meters
% markertr : Array holding ID numbers and coordinates of the control
%           (marker) traces.
% xyz      : Structure with fields:
%            xyz.Tx : [ntr x 3] array with the X, Y and Z coordinates of
%                     the source antenna in a local frame of reference. For
%                     monostatic GPR systems (zero-offset data) xyz.Tx =
%                     xyz.Rx  
%            xyz.Rx : [ntr x 3] array with the X, Y and Z coordinates of
%                     the receiver antenna in a local frame of reference
%                     For monostatic GPR systems (zero-offset data)xyz.Rx =
%                     xyz.Tx  
%
%    Outputs : 
%      d_out : The trimmed GPR data array
%      x_out : The trimmed scan axis coordinate vector
%    ntr_out : The new dimension of the scan axis and horizontal
%             dimension of "d_out".
%markertr_out: Array with updated data on marker traces 
%    xyz_out : The structure or trimed source and receiver coordinates 
%
% Requires : checkcomma.m
%
%   Author   : Andreas Tzanis,
%              Dept. of Geophysics,   
%              University of Athens
%              atzanis@geol.uoa.gr
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

% First Make sure we'll not be shown the door because of a leftover gui
edscfig = findobj('Tag','edscanaxgui');
if ~isempty(edscfig),
    delete(edscfig);
end

global dummyd dummyx dummymark dummyxyz
% The (locally) global variables dummyd will return the edited section
% from the function handle callback "edt"
dummyd      = []; 
dummyx      = []; 
dummymark   = []; 
dummyxyz.Tx = [];
dummyxyz.Rx = [];

ntr = size(d,2);
% Default outputs are : NO ACTION
d_out   = d;
x_out   = x;
ntr_out = ntr;
markertr_out = markertr;
xyz_out = xyz;

% Create GUI interface
figure('Name','编辑扫描轴', 'Tag', 'edscanaxgui', ...
    'NumberTitle','off', 'Position',[300 200 350 350 ], ...
    'MenuBar','none', ...
    'Color',[0.76 0.86 0.69]);         % 'Color',[0 0.5 0.5]);
%    'CreateFcn',['matgprwindows(''setup''); ' ...
%                 'matgprwindows(''updateadd'');'],...
%    'DeleteFcn',['matgprwindows(''updateremove'', ' ...
%                 'findobj(''tag'',''edscanaxgui''));']);
axis off
fcolor = [0.76 0.87 0.78];             % fcolor = get(gcf,'color');

 % Input mode
uicontrol('Style','Radiobutton','Tag','usecursor', 'Units','normalized',...
    'Position',[.05 .8 .35 .05], 'String','采用光标','Value',1,...
    'BackgroundColor', fcolor,...
    'CallBack',@rd1)
uicontrol('Style','Radiobutton','Tag','usefingers','Units','normalized',...
    'Position',[.05 .7 .35 .05], 'String','手动输入','Value',0,...
    'BackgroundColor', fcolor,...
    'CallBack',@rd2)

% Traces groups
if isempty(dx),
    text('Color',[0 0 0 ],  'Units','normalized', ...
        'Position', [-0.125 0.63 0 ],...
        'FontSize',9 , 'FontWeight','bold' , ...
        'String','Trim from START to Trace # ');
else
    text('Color',[0 0 0 ],  'Units','normalized', ...
        'Position', [-0.125 0.63 0 ],...
        'FontSize',9 , 'FontWeight','bold' , ...
        'String','裁掉起始点到此点的内容 (m)');
end
uicontrol('Style','edit', 'Tag', 'startposition', 'enable', 'off', ... 
     'Units','normalized', 'Position',[.05 .54 .35 .06], ...
     'BackgroundColor', fcolor,...
     'String', '');

if isempty(dx),
    text('Color',[0 0 0 ],  'Units','normalized', ...
        'Position', [-0.125 0.5 0 ],...
        'FontSize',9 , 'FontWeight','bold' , ...
        'String','Trim from Trace # to END');
else
    text('Color',[0 0 0 ],  'Units','normalized', ...
        'Position', [-0.125 0.5 0 ],...
        'FontSize',9 , 'FontWeight','bold' , ...
        'String','裁掉此点到终点的内容 (m)');
end  
uicontrol('Style','edit', 'Tag', 'endposition', 'enable', 'off', ...
     'Units', 'normalized', 'Position',[.05 .43 .35 .06], ...
     'BackgroundColor', fcolor,...
     'String', '');
 
% Choose action 
text('Color',[0 0 0 ],  'Units','normalized',  'Position',[0.61 0.93 0 ],...
     'FontSize',9 , 'FontWeight','bold' , 'String','功能');
uicontrol('Style','Checkbox', 'Tag','trimstart', 'Units','normalized',...
    'Position',[0.6 0.80 0.34 0.05],  ...
    'BackgroundColor', fcolor,...
    'String','裁剪前侧', 'Value',0,'CallBack',@cb1);
uicontrol('Style','Checkbox', 'Tag','trimend', 'Units','normalized',...
    'Position',[0.6 0.70 0.34 0.05],  ...
    'BackgroundColor', fcolor,...
    'String','裁剪后侧','Value',0,'CallBack',@cb2)
uicontrol('Style','Checkbox', 'Tag','extractpart', 'Units','normalized',...
    'Position',[0.6 0.60 0.34 0.05],  ...
    'BackgroundColor', fcolor,...
    'String','两侧剪裁','Value',0,'CallBack',@cb3)

%  Proceed ...
hgo = uicontrol('Style','Pushbutton', 'Units','normalized', ...
    'Position',[.1 .15 .25 .1 ], 'String','开始', ...
    'backgroundcolor',  [0.6 1 0], ...
    'Fontweight', 'bold', 'Fontsize',12, ...
    'Callback', {@edt,d,x,dx,markertr,xyz} );
%  Changed my mind ...
uicontrol('Style','Pushbutton','Units','normalized', ...
    'Position',[.65 .15 .25 .1 ],  'String','取消', ...
    'backgroundcolor', [1 0.2 0], ...
    'Fontweight', 'bold', 'Fontsize',12, ...
    'Callback', 'delete(findobj(''tag'',''edscanaxgui'')); ');
% Wait until object hgo is deleted, otherwise EDITSCANAXIS will return
% and  when hgo's callback is executed there will be no result passed to
% the caller (matGPR) 
waitfor(hgo)
% Assign the output variables 
d_out        = dummyd;
x_out        = dummyx;
ntr_out      = size(d_out,2);
markertr_out = dummymark;
xyz_out      = dummyxyz; 
clear global dummyd dummyx dummymark dummyxyz
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edt(hgo, eventdata, d, x, dx, markertr, xyz)
% Function callback of 'GO' pushbutton 
global dummyd dummyx dummymark dummyxyz
handles = guihandles(gcf);
ntr = size(d,2);
nfrom = 1;
nto   = ntr;
%%%%%%% Fingertip input mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (get(handles.usefingers,'Value') == get(handles.usefingers,'Max')),

    if (get(handles.trimstart,'Value') == get(handles.trimstart,'Max')) || ...
            (get(handles.extractpart,'Value') == get(handles.extractpart,'Max')),
        str = get(findobj('tag','startposition'),'string');
        str = checkcomma(str);
        if isempty(num2str(str)),
            dummyd = [];
            delete(findobj('tag','edscanaxgui'));  
            return
        end
        nfrom = str2num(str);
        if ~isempty(dx),
            nfrom = floor(nfrom/dx);
        end
    end
        
    if (get(handles.trimend,'Value') == get(handles.trimend,'Max')) || ...
            (get(handles.extractpart,'Value') == get(handles.extractpart,'Max')),
        str = get(findobj('tag','endposition'),'string');
        str = checkcomma(str);
        if isempty(num2str(str)),
            dummyd = [];
            delete(findobj('tag','edscanaxgui'));  
            return
        end
        nto  = str2num(str);
        if ~isempty(dx),
            nto = floor(nto/dx);
        end
    end
end

%%%%%%% Cursor input mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (get(handles.usecursor,'Value') == get(handles.usecursor,'Max')),

    if (get(handles.trimstart,'Value') == get(handles.trimstart,'Max')) || ...
            (get(handles.extractpart,'Value') == get(handles.extractpart,'Max')),
%%%%%   Get the to trim from
        datafig = findobj('tag','datafigure');
        if ~ishandle(datafig),
            erh = errordlg('No data figure for the filter design interface',...
                'EDIT SCAN AXIS : ERROR');
            uiwait(erh);
            dummyd = [];
            delete(findobj('tag','edscanaxgui'));  
            return
         end
         message     = msgbox('Please define trace to trim from ',...
             'EDIT SCAN AXIS : HELP');
         messagepos  = get(message,'position');
         set(message,'pos',[20 60 150 messagepos(4)]);
         pause(0.5)
         figure(datafig)
         %pause(0.5)
         %[xtr,ytr]=ginput(1);
         if isempty(dx),
             [xtr,ytr]=getpoint('', ' (ns)'); 
             nfrom = floor(xtr);
         else
             [xtr,ytr]=getpoint(' (m)', ' (ns)'); 
             nfrom = floor((xtr-x(1))/dx);
         end
         close(message);
    end
    
    if (get(handles.trimend,'Value') == get(handles.trimend,'Max')) || ...
            (get(handles.extractpart,'Value') == get(handles.extractpart,'Max')),
        
        datafig = findobj('tag','datafigure');
        if ~ishandle(datafig),
            erh = errordlg('No data figure for the filter design interface',...
                'EDIT SCAN AXIS : ERROR');
            uiwait(erh);
            dummyd = [];
            delete(findobj('tag','edscanaxgui'));  
            return
         end
         message     = msgbox('Please define trace to trim to ',...
             'EDIT SCAN AXIS : HELP');
         messagepos  = get(message,'position');
         set(message,'pos',[20 60 150 messagepos(4)]);
         pause(0.5)
         figure(datafig)
         %pause(0.5)
         %[xtr,ytr]=ginput(1);
         if isempty(dx),
             [xtr,ytr]=getpoint('', ' (ns)');   
             nto = floor(xtr);
         else
             [xtr,ytr]=getpoint(' (m)', ' (ns)');  
             nto = floor((xtr-x(1))/dx);
         end
         close(message);
    end
    
end
% Trap a quite usual error
if nfrom >= nto,
    disp('EDIT SCAN AXIS > Last trace >= First trace! Sorting ...');
    tmp   = nto;
    nto   = nfrom;
    nfrom = tmp;
end    
dummyd = d(:,nfrom:nto);
ntrnew = size(dummyd,2);
if isempty(dx),
    dummyx = 1 : 1: ntrnew;
else
    dummyx = x(nfrom) + [0 : dx : (ntrnew -1)*dx];
end
% Update marker trace information
if ~isempty(markertr),
    marks     = find(markertr(:,1) >= nfrom & markertr(:,1) <= nto);
    dummymark = markertr(marks,:);
% Only the marker IDs modified - the coordinates should not change, as they
% correspond to the reference frame
    dummymark(:,1) = dummymark(:,1)- nfrom +1;
end
if ~isempty(xyz.Tx),
    dummyxyz.Tx = xyz.Tx(nfrom : nto,:);
    dummyxyz.Rx = xyz.Rx(nfrom : nto,:);
end
delete(findobj('tag','edscanaxgui'));  
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rd1(usecursor,eventdata)
handles = guihandles(gcf);
if (get(usecursor,'Value') == get(usecursor,'Max')),
    set(handles.usefingers,'value',0)
    set(handles.startposition,'String','',...
        'enable','off')
    set(handles.endposition,'String','',...
        'enable','off')
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rd2(usefingers,eventdata)
handles = guihandles(gcf);
if (get(usefingers,'Value') == get(usefingers,'Max')),
    set(handles.usecursor,'value',0)
end
if (get(handles.trimstart,'Value') == get(handles.trimstart,'Max')),
    set(handles.startposition,'string', '', 'enable','on')
    set(handles.endposition,'String','',...
        'enable','off')
end
if (get(handles.trimend,'Value') == get(handles.trimend,'Max')),
    set(handles.startposition,'String','',...
        'enable','off')
    set(handles.endposition,'string', '', 'enable','on')
end
if (get(handles.extractpart,'Value') == get(handles.extractpart,'Max')),
    set(handles.startposition,'string', '', 'enable','on')
    set(handles.endposition,'string', '', 'enable','on')
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb1(trimstart, eventdata)
% Function callback of checkbox 'trimstart'
handles = guihandles(gcf);
if (get(trimstart,'Value') == get(trimstart,'Max')),
    set(handles.trimend,'value',0)
    set(handles.extractpart,'value',0)
end
if (get(handles.usecursor,'Value') == get(handles.usecursor,'Max')),
    set(handles.startposition,'String','',...
        'enable','off')
    set(handles.endposition,'String','',...
        'enable','off')
end
if (get(handles.usefingers,'Value') == get(handles.usefingers,'Max')),
    set(handles.startposition,'string', '', 'enable','on')
    set(handles.endposition,'String','',...
        'enable','off')
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb2(trimend, eventdata)
% Function callback of checkbox 'trimend'
handles = guihandles(gcf);
if (get(trimend,'Value') == get(trimend,'Max')),
    set(handles.trimstart,'value',0)
    set(handles.extractpart,'value',0)
end
if (get(handles.usecursor,'Value') == get(handles.usecursor,'Max')),
    set(handles.startposition,'String','',...
        'enable','off')
    set(handles.endposition,'String','',...
        'enable','off')
end
if (get(handles.usefingers,'Value') == get(handles.usefingers,'Max')),
    set(handles.startposition,'String','',...
        'enable','off')
    set(handles.endposition,'string', '', 'enable','on')
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb3(extractpart, eventdata)
% Function callback of checkbox 'extractpart'
handles = guihandles(gcf);
if (get(extractpart,'Value') == get(extractpart,'Max')),
    set(handles.trimstart,'value',0)
    set(handles.trimend,'value',0)
end
if (get(handles.usecursor,'Value') == get(handles.usecursor,'Max')),
    set(handles.startposition,'String','',...
        'enable','off')
    set(handles.endposition,'String','',...
        'enable','off')
end
if (get(handles.usefingers,'Value') == get(handles.usefingers,'Max')),
    set(handles.startposition,'string', '', 'enable','on')
    set(handles.endposition,'string', '', 'enable','on')
end
return
