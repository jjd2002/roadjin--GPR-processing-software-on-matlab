function dstc =staticsgui(d, dt, xyz_Tx, xyz_Rx)
%
% STATICSGUI: Interface and driver for program "static.m", to do
% topographic corrections. The routine initializes a GUI, in which the user
% will specify the velocity of the weathering and subweathering layers to
% be used in the topographic correction, the sense of the correction (up or
% down), and the elevation datum in which to refer the reduced GPR section.
% Then, it will drive the routine "static.m" via the function handle
% callback "stc". The statics-corrected section returns from "stc" by means
% of an interim global variable that is kept strictly within the scope of
% "staticsgui". The interim variable is copied onto the output variable
% just before "staticsgui" returns.  
%
%  Inputs :
%      d  : (m x n) matrix of the observed data (GPR section);
%     dt  : Sampling rate (in ns)
%  xyz_Tx : [n x 3] array with the X, Y and Z coordinates of the source
%           antenna in a local frame of reference. For monostatic GPR
%           systems (zero-offset data) xyz_Tx = xyz_Rx
%  xyz_Rx : [n x 3] array with the X, Y and Z coordinates of the receiver
%           antenna in a local frame of reference For monostatic GPR 
%           systems (zero-offset data) xyz_Rx = xyz_Tx
%
%  Output : 
%    dstc : (m x n) matrix of the statics-corrected data
%
%Requires : static.m
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
%

global dummystc
% The (locally) global variable dummystc will return the corrected section
% from the function handle callback "stc"
dummystc = []; 
dstc = [];

% Trap errors
No_traces = size(d,2);
if isempty(xyz_Tx),
    erh = errordlg('Elevation data does not exist! Aborting!',...
        'STATICSGUI : ERROR');
    uiwait(erh);  
    return
end
if size(xyz_Tx,1) ~= No_traces,
    erh = errordlg(['Elevation data and radargram are INCONSISTENT '...
            '(different number of traces)!'], 'STATICS : ERROR');
    uiwait(erh);  
    return
end
if size(xyz_Tx,2) ~= 3,
    erh = errordlg(['The topographic data matrix has non-standard structure!'...
            'Do not know how to handle it!'], 'STATICS : ERROR');
    uiwait(erh);  
    return
end

% Default sub-weathering layer velocity
swv = 0.1;
% The weathering layer velocity is useful only for bistatic data.
% For monostatic data it doesn't matter what it will be!
wv = 0.1;

% Create GUI interface
figure('Name','Static correction Parametres', 'Tag', 'staticgui', ...
    'NumberTitle','off', 'Position',[300 200 350 350 ], ...
    'MenuBar','none', ...
    'Color',[0.76 0.86 0.69]);         % 'Color',[0 0.5 0.5]);
%    'CreateFcn',['matgprwindows(''setup''); ' ...
%                 'matgprwindows(''updateadd'');'],...
%    'DeleteFcn',['matgprwindows(''updateremove'', ' ...
%                 'findobj(''tag'',''staticgui''));']);
axis off
fcolor = [0.76 0.87 0.78];             % fcolor = get(gcf,'color');

% Migration Velocities for static corrections
text('Color',[0 0 0 ],  'Units','normalized', ...
    'Position', [-0.125 0.93 0 ], 'FontSize',9 , 'FontWeight','bold' , ...
    'String','Weathering Layer Velocity');
uicontrol('Style','edit', 'Tag', 'wvgui', 'enable', 'off', ... 
    'Units','normalized', 'Position',[.05 .80 .15 .05], ...
    'BackgroundColor', fcolor,...
    'String',num2str(wv));
     
text('Color',[0 0 0 ],  'Units','normalized', ...
    'Position',[-0.125 0.80 0 ],'FontSize',9 , 'FontWeight','bold', ...
    'String','Subweathering Velocity');
uicontrol('Style','edit', 'Tag', 'swvgui', ...
    'Units', 'normalized', 'Position',[.05 .69 .15 .05], ...
    'BackgroundColor', fcolor,...
    'String',num2str(swv));

 % Direction of static corrections    
uicontrol('Style','Radiobutton','Tag','shiftdn', 'Units','normalized',...
    'Position',[.05 .60 .35 .05], ...
    'String','Shift down (default)','Value',1,...
    'BackgroundColor', fcolor,...
    'CallBack',@rd1)
uicontrol('Style','Radiobutton','Tag','shiftup','Units','normalized',...
    'Position',[.05 .50 .25 .05], 'String','Shift up','Value',0, ...
    'BackgroundColor', fcolor,...
    'CallBack',@rd2)

% Choose elevation datum
text('Color',[0 0 0 ],  'Units','normalized', ...
    'Position',[0.61 0.93 0 ], 'FontSize',9 , 'FontWeight','bold' , ...
    'String','Elevation Datum');
uicontrol('Style','Checkbox', 'Tag','datummax', 'Units','normalized',...
    'Position',[0.6 0.80 0.34 0.05],  ...
    'BackgroundColor', fcolor,...
    'String','Maximum elevation', 'Value',1,'CallBack',@cb1);
uicontrol('Style','Checkbox', 'Tag','datummin', 'Units','normalized',...
    'Position',[0.6 0.70 0.34 0.05],  ...
    'BackgroundColor', fcolor,...
    'String','Minimum elevation','Value',0,'CallBack',@cb2)
uicontrol('Style','Checkbox', 'Tag','datumave', 'Units','normalized',...
    'Position',[0.6 0.60 0.34 0.05],  ...
    'BackgroundColor', fcolor,...
    'String','Mean elevation','Value',0,'CallBack',@cb3)
uicontrol('Style','Checkbox', 'Tag', 'datum1tr', 'Units','normalized',...
    'Position',[0.6 0.50 0.34 0.05],  ...
    'BackgroundColor', fcolor,...
    'String','First Trace','Value',0,'CallBack',@cb4)
text('Color',[0 0 0 ], 'Position',[0.925 0.385 0 ], 'FontSize',10 , ...
    'FontWeight','bold' , 'String','User');
uicontrol('Style','edit', 'Tag', 'datumusr', 'Units','normalized',...
    'BackgroundColor', fcolor,...
    'Position',[.6 .40 .225 .05],  'String','', 'CallBack', @edtusr)

%  Proceed ...
hgo = uicontrol('Style','Pushbutton', 'Units','normalized', ...
    'Position',[.1 .15 .25 .1 ], 'String','Go', ...
    'backgroundcolor', [0.6 1 0], ...     
    'Fontweight', 'bold', 'Fontsize',12, ...
    'Callback', {@stc,d,dt,xyz_Tx,xyz_Rx} );
%  Changed my mind ...
uicontrol('Style','Pushbutton','Units','normalized', ...
    'Position',[.65 .15 .25 .1 ],  'String','Cancel', ...
    'backgroundcolor', [1 0.2 0], ...
    'Fontweight', 'bold', 'Fontsize',12, ...
    'Callback', 'delete(findobj(''tag'',''staticgui'')); ');
% Wait until object hgo is deleted, otherwise STATICSGUI will return and 
% when hgo's callback is executed there will be no result to pass to the 
% caller (do_statics) 
waitfor(hgo)
% Assign the output variable dstc (statically corrected GPR section)
dstc = dummystc;
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stc(hgo, eventdata, d, dt, xyz_Tx, xyz_Rx)
% Function callback of 'GO' pushbutton to invoke function 'static' and do
% static computations
global dummystc
% The weathering layer velocity is useful only for bistatic data.
% For monostatic data it doesn't matter what it will be!
strwv = get(findobj('tag','wvgui'),'String');  
comma = findstr(strwv,',');
if ~isempty(comma), 
    strwv(comma) = '.'; 
end
wv = str2num(strwv);
strswv= get(findobj('tag','swvgui'),'String');
comma = findstr(strswv,',');
if ~isempty(comma), 
    strswv(comma) = '.'; 
end
swv = str2num(strswv);
% Trap a very easy to do error
if wv <= 0 || wv > 0.2998 || swv <= 0 || swv > 0.2998,
    erh = errordlg('Impossible velocity value! Please try again!',...
            'STATICSGUI: ERROR');
    uiwait(erh);
    return
end
if get(findobj('tag','shiftdn'),'value')== 1;
    direction = -1;
end; 
if get(findobj('tag','shiftup'),'value') == 1;
    direction =  1;
end; 
if get(findobj('tag','datummax'),'value') == 1;
    datum= max(xyz_Tx(:,3));
end; 
if get(findobj('tag','datummin'),'value') == 1; 
    datum= min(xyz_Tx(:,3)); 
end; 
if get(findobj('tag','datumave'),'value') == 1; 
    datum= mean(xyz_Tx(:,3));
end;
if get(findobj('tag','datum1tr'),'value') == 1; 
    datum= xyz_Tx(1,3); 
end;
if ~isempty(str2num(get(findobj('tag','datumusr'),'string'))); 
    datum = str2num(get(findobj('tag','datumusr'),'String')); 
end; 
dummystc = static(d,dt,xyz_Tx,xyz_Rx,datum,wv,swv,direction); 
delete(findobj('tag','staticgui'));  
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rd1(shiftdn,eventdata)
handles = guihandles(gcf);
if (get(shiftdn,'Value') == get(shiftdn,'Max')),
    set(handles.shiftup,'value',0)
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rd2(shiftup,eventdata)
handles = guihandles(gcf);
if (get(shiftup,'Value') == get(shiftup,'Max')),
    set(handles.shiftdn,'value',0)
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb1(datummax, eventdata)
% Function callback of checkbox 'datummax'
handles = guihandles(gcf);
if (get(datummax,'Value') == get(datummax,'Max')),
    set(handles.datummin,'value',0)
    set(handles.datumave,'value',0)
    set(handles.datum1tr,'value',0)
    set(handles.datumusr,'string','')
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb2(datummin, eventdata)
% Function callback of checkbox 'datummin'
handles = guihandles(gcf);
if (get(datummin,'Value') == get(datummin,'Max')),
    set(handles.datummax,'value',0)
    set(handles.datumave,'value',0)
    set(handles.datum1tr,'value',0)
    set(handles.datumusr,'string','')
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb3(datumave, eventdata)
% Function callback of checkbox 'datumave'
handles = guihandles(gcf);
if (get(datumave,'Value') == get(datumave,'Max')),
    set(handles.datummax,'value',0)
    set(handles.datummin,'value',0)
    set(handles.datum1tr,'value',0)
    set(handles.datumusr,'string','')
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb4(datum1tr, eventdata)
% Function callback of checkbox 'datum1tr'
handles = guihandles(gcf);
if (get(datum1tr,'Value') == get(datum1tr,'Max')),
    set(handles.datummax,'value',0)
    set(handles.datumave,'value',0)
    set(handles.datummin,'value',0)
    set(handles.datumusr,'string','')
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edtusr(datumusr, eventdata)
% Function callback of checkbox 'datumusr'
handles = guihandles(gcf);
if ~isempty(get(datumusr,'string')),
    set(handles.datummin,'value',0)
    set(handles.datummax,'value',0)
    set(handles.datumave,'value',0)
    set(handles.datum1tr,'value',0)
end
return
