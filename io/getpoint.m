function [x, y, button] = getpoint(s1,s2)
%
%GETPOINT : Interface for function "getxy.m". Grab from current axes, a
%           pair of points [ x y ] by pointing and clicking. Accurate
%           pointing is facilitated by a streaming display at the bottom of
%           the figure, of the current pointer position (coordinates) in
%           data units. The display is set up by "getxy" and is erased as
%           soon as the program returns. 
%       ==> Returns the number of mouse button pressed. 
%       ==> In effect, this is a surrogate "ginput" routine, without
%           keyboard support, but with the added flexibility of letting you 
%           know the exact location of the pointer.
%
%   Usage : [x, y, button] = getpoint(s1,s2);
%           [x, y ] = getpoint(s1,s2);
%           [x, y ] = getpoint;
%
%  Inputs : 
%  s1, s2 : Respectively the physical units of the data in the X- and Y-
%           axes, e.g. 'meters', 'nanoseconds', etc.
%
% Outputs : 
%    x, y : The coordinates of the clicked point in data units.
%  button : The mouse button pressed: 1 = LEFT   button
%                                     2 = MIDDLE button
%                                     3 = RIGHT  button
%
%    Note : The routine communicates with "getxy" through global variables.
%           These are kept strictly within the scope of "getpoint". For
%           details see the self-documentation of  "getpoint" below.
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
global xypoint xstr ystr buttonpressed
x    = [];    y    = [];   button=[]; 
xstr = '';    ystr = '';  
if nargin >= 1,
    xstr = s1;
end
if nargin == 2,
    ystr = s2;
end
getxy('begin');
x      = [x; xypoint(1)];
y      = [y; xypoint(2)];
button = buttonpressed;
clear global xypoint xstr ystr buttonpressed
pause(0.25)
return
%%% End function GETPOINT

function getxy(action)
%
%   GETXY : Get from the current axes, a pair of points [ x, y ], by
%           pointing and clicking. Accurate pointing is facilitated by a 
%           streaming display at the bottom of the figure, of the current
%           pointer position (coordinates) in data units. The display is
%           set up by "getxy" and is erased as  soon as it exits returns.  
%       ==> Also returns the number of mouse button pressed. 
%
%   Usage : getxy('begin');
%
%  Inputs : 
%  action : Keyword determining the runtime behaviour of "getxy". 
%           ==> This is a multiple-entry recursive program. Each entry is
%               activated by a keyword. The keywords are pre-programmed.
%               The only option available (and permisible) to the a user is
%               action = 'begin', which initializes the program, setting up
%               its interfaces on the current axes.
%
% Outputs : None. The routine communicates its runtime parameters and the
%           results to the caller (driver) function "getpoint.m" (or any
%           other caller function thereof), via global variables. These
%           variables are:
%           xypoint : [1 x 2] vector with the x, y coordinates of the
%                     clicked point in data units.
%           button  : The mouse button pressed: 1 = LEFT   button
%                                               2 = MIDDLE button
%                                               3 = RIGHT  button
%        xstr, ystr : Respectively strings indicating the physical units of
%                     the X- and Y- coordinates of the data plotted in the
%                     current axes, e.g. 'meters', 'nanoseconds', etc.
%
%  Author : Andreas Tzanis
%           Department of Geophysics, 
%           University of Athens
%           atzanis@geol.uoa.gr
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%
global xypoint xstr ystr buttonpressed
%%% Setup the stage for grabbing a point from the axes
if strcmp(action,'begin'),
    uicontrol('style','edit','tag','newx','parent',gcf,...
        'units','normalized','string','X-COORDINATES' ,...
        'position',[0.05 0.01 0.3 0.05],'horizontalal','center',...
        'fontsize',9);
    uicontrol('style','edit','tag','newy','parent',gcf,...
        'units','normalized','string','Y-COORDINATES' ,...
        'position',[0.65 0.01 0.3 0.05],'horizontalal','center',...
        'fontsize',9);
    set(gcf,'Backingstore','off',...
        'WindowButtonMotionFcn','getxy(''motion'');',...
        'WindowButtonDownFcn', ...
        'getxy(''down'');');
% Interrupt execution of the caller function until getxy exits
% gracefully
    waitfor(findobj('tag','newx'))
end

if strcmp(action,'motion'),
    xl = get(gca,'xlim');               % x-axis limits
    yl = get(gca,'ylim');               % y-axis limits
    cp = get(gca,'currentpoint');
    cx = cp(1,1);
    cy = cp(1,2);
    online = cy > yl(1) && cy < yl(2) && cx > xl(1) && cx < xl(2) ;
    if online,
        set(gcf,'Pointer','crosshair');
        str = ['X = ' num2str(cx) xstr];
        set(findobj('tag','newx'),'string',str)
        str = ['Y = ' num2str(cy) ystr];
        set(findobj('tag','newy'),'string',str)
    elseif ~online,
        set(gcf,'Pointer','arrow');
        set(findobj('tag','newx'),'string','Out of range' )
        set(findobj('tag','newy'),'string','Out of range' )
    end
end

%%% Set cursor shape and display cursor position 
if strcmp(action,'down'),
%%% Left mouse button clicked: Set a point 
        cp = get(gca ,'currentpoint');       % Get coordinates
        xypoint = [ cp(1,1) cp(1,2) ];
        if strcmp(get(gcf,'selectiontype'),'normal'),
            buttonpressed = 1;
        elseif strcmp(get(gcf,'selectiontype'),'alt'),
            buttonpressed = 3;
        elseif strcmp(get(gcf,'selectiontype'),'extend'),
            buttonpressed = 2;
        end
        set(gcf,'Backingstore','off','WindowButtonMotionFcn',''); 
        set(gcf,'WindowButtonDownFcn', '');
        set(gcf,'pointer','arrow')
        delete(findobj('tag','newx'))
        delete(findobj('tag','newy'))
end
return
%%% End function GETXY    
