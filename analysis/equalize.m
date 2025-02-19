function deq = equalize( d )
%
% EQUALIZE : Equalizes traces by making the sum of the absolute values of
% all samples in a trace the same for all traces. For instance, if the
% trace that is selected to be used for equalization has a sum of all
% absolute values of 1000 (base = 1000), then the sum of absolute values of
% every trace in the data set will be set to 1000 using a multiplying
% factor. The routine initializes a GUI, in which the user specifies the
% trace which will be used as the BASE. This can be:
% a) the first trace
% b) the mean (average) trace
% c) the median trace
% d) the n'th trace (i.e. any trace) specified by the user
% e) Any positive number (base value) provided by the user.
% Then, it will process the data via the function handle callback "eqlz".
% The reduced data set returns from "eqlz" by means of an interim global
% variable that are kept strictly within the scope of "equalize". The
% interim variable are copied onto the output variables just before
% "equalize" returns.   
%
% USAGE : deq = equalize( d )
%
%  Inputs :
%       d : The [ns x ntr] GPR data array
%
% Outputs : 
%     deq : The equalized-trace GPR data array
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

% First make sure we'll not be shown the door out because of a leftover gui
eqfig = findobj('Tag','equalizegui');
if ~isempty(eqfig),
    delete(eqfig);
end

global dummyd 
% The (locally) global variables dummyd will return the edited section
% from the function handle callback "edt"
dummyd  = []; 

% Default outputs are : VOID
deq     = [];

% Create GUI interface
figure('Name','Equalize Traces', 'Tag', 'equalizegui', ...
    'NumberTitle','off', 'Position',[300 200 350 350 ], ...
    'MenuBar','none', 'Color',[0.76 0.86 0.69]);
%    'CreateFcn',['matgprwindows(''setup''); ' ...
%                 'matgprwindows(''updateadd'');'],...
%    'DeleteFcn',['matgprwindows(''updateremove'', ' ...
%                 'findobj(''tag'',''edscanaxgui''));']);
axis off
fcolor = [0.76 0.87 0.78];             % fcolor = get(gcf,'color');

% Equalize to: 
text('Color',[0 0 0 ],  'Units','normalized',  'Position',[0.11 0.93 0 ],...
     'FontSize',9 , 'FontWeight','bold' , 'String','Equalize to:');
uicontrol('Style','Checkbox', 'Tag','eqfirst', 'Units','normalized',...
    'Position',[0.2 0.80 0.34 0.05],  ...
    'backgroundcolor',fcolor, ...
    'String','First Trace', 'Value',0,'CallBack',@cb1);
uicontrol('Style','Checkbox', 'Tag','eqmean', 'Units','normalized',...
    'backgroundcolor',fcolor, ...
    'Position',[0.2 0.70 0.34 0.05],  ...
    'String','Mean Trace','Value',0,'CallBack',@cb2)
uicontrol('Style','Checkbox', 'Tag','eqmed', 'Units','normalized',...
    'Position',[0.2 0.60 0.34 0.05],  ...
    'backgroundcolor',fcolor, ...
    'String','Median Trace','Value',0,'CallBack',@cb3)
text('Color',[0 0 0 ], 'Position',[0.425 0.505 0 ], 'FontSize',10 , ...
    'FontWeight','bold' , 'String','Give Trace #');
uicontrol('Style','edit', 'Tag', 'equser2', 'Units','normalized',...
    'Position',[.2 .48 .225 .075], ...
    'backgroundcolor',fcolor, ...
    'String','', 'CallBack', @equsr2)
text('Color',[0 0 0 ], 'Position',[0.425 0.385 0 ], 'FontSize',10 , ...
    'FontWeight','bold' , 'String','Give Base');
uicontrol('Style','edit', 'Tag', 'equser1', 'Units','normalized',...
    'Position',[.2 .38 .225 .075],  ...
    'backgroundcolor',fcolor, ...
    'String','', 'CallBack', @equsr1)

%  Proceed ...
hgo = uicontrol('Style','Pushbutton', 'Units','normalized', ...
    'Position',[.1 .15 .2 .1 ], 'String','Go', ...
    'backgroundcolor', [0.6 1 0], ...
    'Fontweight', 'bold', 'Fontsize',12, ...
    'Callback', {@eqlz,d} );
%  Changed my mind ...
uicontrol('Style','Pushbutton','Units','normalized', ...
    'Position',[.4 .15 .2 .1 ],  'String','Cancel', ...
    'backgroundcolor', [1 0.2 0], ...
    'Fontweight', 'bold', 'Fontsize',12, ...
    'Callback', 'delete(findobj(''tag'',''equalizegui'')); ');
% Wait until object hgo is deleted, otherwise EQUALIZE will return
% and  when hgo's callback is executed there will be no result passed to
% the caller (matGPR) 
waitfor(hgo)
% Assign the output variables 
deq  = dummyd;
clear global dummyd 
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eqlz(hgo, eventdata, d)
% Function callback of 'GO' pushbutton 
global dummyd 
handles   = guihandles(gcf);
[ns, ntr] = size(d);
dummyd    = zeros(ns,ntr);
base      = [];

% Get base value for equalization 
if (get(handles.eqfirst,'Value') == get(handles.eqfirst,'Max')),
    base = sum(abs(d(:,1)));
end
if (get(handles.eqmean,'Value') == get(handles.eqmean,'Max')),
    base = sum(abs(mean(d')));
end
if (get(handles.eqmed,'Value') == get(handles.eqmed,'Max')) ,
    base = sum(abs(median(d')));
end
if ~isempty(get(findobj('tag','equser2'),'string')),
    numtrace = floor(str2num(get(findobj('tag','equser2'),'string')));
    base = sum(abs(d(:,numtrace)));
end
if ~isempty(get(findobj('tag','equser1'),'string')),
    base = str2num(get(findobj('tag','equser1'),'string'));
end
if isempty(base),
    erh = errordlg('Equalization base NOT specified! Aborting!',...
            'EQUALIZEGUI: ERROR');
    uiwait(erh);
    dummyd = [];
    delete(findobj('tag','equalizegui'));
    return
end

% Equalize to base value
for i=1:ntr; 
    ss = sum(abs(d(:,i))); 
    factor = base/ss; 
    dummyd(:,i) = d(:,i)*factor; 
end;

% Exit
delete(findobj('tag','equalizegui'));
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb1(eqfirst, eventdata)
% Function callback of checkbox 'eqfirst'
handles = guihandles(gcf);
if (get(eqfirst,'Value') == get(eqfirst,'Max')),
    set(handles.eqmean,'value',0)
    set(handles.eqmed,'value',0)
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb2(eqmean, eventdata)
% Function callback of checkbox 'eqmean'
handles = guihandles(gcf);
if (get(eqmean,'Value') == get(eqmean,'Max')),
    set(handles.eqfirst,'value',0)
    set(handles.eqmed,'value',0)
    set(handles.equser1,'string','')
    set(handles.equser2,'string','')
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb3(eqmed, eventdata)
% Function callback of checkbox 'eqmed'
handles = guihandles(gcf);
if (get(eqmed,'Value') == get(eqmed,'Max')),
    set(handles.eqfirst,'value',0)
    set(handles.eqmean,'value',0)
    set(handles.equser1,'string','')
    set(handles.equser2,'string','')
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function equsr2(equser2, eventdata)
% Function callback of checkbox 'equsr2'
handles = guihandles(gcf);
if ~isempty(get(equser2,'string')),
    set(handles.eqfirst,'value',0)
    set(handles.eqmean,'value',0)
    set(handles.eqmed,'value',0)
    set(handles.equser1,'string','')
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function equsr1(equser1, eventdata)
% Function callback of checkbox 'equsr1'
handles = guihandles(gcf);
if ~isempty(get(equser1,'string')),
    set(handles.eqfirst,'value',0)
    set(handles.eqmean,'value',0)
    set(handles.eqmed,'value',0)
    set(handles.equser2,'string','')
end
return
