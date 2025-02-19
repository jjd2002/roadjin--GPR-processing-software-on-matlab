function [dsp, ns_new, tt_new] = sigposition(action, dt, dx, d);
%
% SIGPOSITION : Determine and adjust the true time zero of the GPR data
%               passed with the array "d".
%
%     Usage   : [dsp, ns_new, tt_new ] = sigposition('begin',dt,dx,d);
%               This is a multiple-entry program, recursively calling upon
%               itself. Always launch the function with keyword action =
%               'begin'. 
%
%     Inputs  :
%       d     : The (ns x n) GPR section 
%       dt    : The sampling rate 
%       dx    : The trace spacing in m (empty if traces are unequally
%               spaced). 
%      action : is a keyword. This is a multiple-entry recursive program.
%               Each entry is activated by a keyword. The keywords are
%               pre-programmed. For manual operation, always launch with
%               action = 'begin'. 
%
%     Outputs :
%       dsp   : (ns_new x n) matrix of the adjusted radargram
%      ns_new : is the new value of ns (samples per scan)
%      tt_new : The updated vector of time coordinates (2-way traveltime)
%
% REQUIRES: viewtraces.m
%
% Author  : Andreas Tzanis
%           Department of Geophysics, 
%           University of Athens
%           atzanis@geol.uoa.gr
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%

global dsp
% Note: The result is passed from the callback subfunction GETSP as a
% global variable "dsp", strictly kept within the scope of SIGPOSITION

if ~strcmp(action,'begin')   
    trace_fig   = findobj('Tag', 'viewindatatraces');
    trace_axes  = findobj('Tag','viewintraceaxes');
    xrange = get(trace_axes,'XLim');
    yrange = get(trace_axes,'YLim');
end

if isequal(action,'begin'),
    viewtraces(dt,dx,d,'indata');
    trace_fig   = findobj('Tag', 'viewindatatraces');
    trace_axes  = findobj('Tag','viewintraceaxes');
    set(trace_fig,'Backingstore','off','WindowButtonMotionFcn',...
        'sigposition(''motion'');');
    set(trace_fig,'WindowButtonDownFcn',{@getsp,dt,dx,d});
%%% Wait until the object "viewindatatraces" is deleted, otherwise 
%%% SIGPOTITION will return and when the WindowButtonDown function callback
%%% is executed, there will be no result passed to the caller (matGPR) 
    waitfor(findobj('tag','viewindatatraces'));
    ns_new = size(dsp,1);
    tt_new = 0 : dt : (ns_new - 1)*dt;
end

if strcmp(action,'motion'),
    cursorstate = get(trace_fig,'Pointer');
    cp = get(trace_axes,'CurrentPoint');
    cx = cp(1,1);
    cy = cp(1,2);
    online = cy > yrange(1) & cy < yrange(2) & ...
        cx > xrange(1) & cx < xrange(2) ;
    if online && strcmp(cursorstate,'arrow'),
        set(trace_fig,'Pointer','crosshair');
    elseif ~online && strcmp(cursorstate,'crosshair'),
        set(trace_fig,'Pointer','arrow');
    end
    if online, 
        cpos = sprintf('t = %7.3f ns',cy);
        set(get(trace_axes,'title'),'string',cpos);
    else
        set(get(trace_axes,'title'),'string','');
    end
end

if strcmp(action,'down'),
    set(trace_fig,'Pointer','crosshair');
    cp = get(trace_axes,'CurrentPoint');
    newx=cp(1,1);
    if newx > xrange(2)
       newx = xrange(2);
    end
    if newx < xrange(1)
       newx = xrange(1);
    end
    newt=cp(1,2);
    if newt > yrange(2)
       newt = yrange(2);
    end
    if newt < yrange(1)
       newt = yrange(1);
    end

    new_zero = floor(newt/dt);
    picked = ['现在选择的位置为 t = ' num2str(dt*new_zero,4) ' ns'];
    set(get(trace_axes,'title'),'string',picked);
    picked = [picked ', Sample = ' num2str(new_zero)];
    clb = { picked '是否重新选择 ？' };
    qbut = questdlg(clb,'SIGPOSITION: Request','是','否','取消','是');
    switch qbut
    case '否'
        answer = inputdlg('Give Signal (Sample #):','Please confirm',...
            1,cellstr(num2str(new_zero)));
        if isempty(answer),
            dsp = [];
            delete(trace_fig);
            return
        end
        new_zero = str2num(answer{1});
        [ns, ntr] = size(d);
        newns = ns - new_zero;
        if mod(newns,2) ~=0,
            new_zero = new_zero - 1;
        end
        dsp    = d(new_zero+1:ns,:);
        delete(trace_fig);
        return
    case '是'
        dsp = [];
    case '取消'
        dsp = [];
        delete(trace_fig);
        return
    end
    set(trace_fig,'WindowButtonMotionFcn','sigposition(''motion'');');
    set(trace_fig,'WindowButtonUpFcn','sigposition(''up'');');
end

if strcmp(action,'up'),
    set(trace_fig,'WindowButtonMotionFcn','sigposition(''motion'');');
    set(trace_fig,'WindowButtonUpFcn','');
end

function getsp(trace_fig,eventdata,dt,dx,d)
% Function handle callback of the 'WindowButtonDown function - now the
% true time zero will be defined and the data matrix will be adjusted
global dsp
dsp = sigposition('down',dt,dx,d);
return