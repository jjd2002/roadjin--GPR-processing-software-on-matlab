function [x, y ] = setpolygon(str1, str2)
%
% SETPOLYGON : Interface for function "setpoly". Define, on the current
%              axes, a polygonal line comprising a set of points [ x(n)
%              y(n) ] by pointing and clicking. 
%          ==> The polygonal line is clipped with respect to the axes
%              limits.  
%          ==> Accurate pointing is facilitated by a streaming display at
%              the bottom of the figure, of the current pointer position
%              (coordinates) in data units. The display is set up by
%              "setpoly" and is erased as soon as the program returns. 
%          ==> Work flow is controlled exclusively with the mouse: 
%              ** To set an x,y pair of coordinates click the LEFT mouse
%                 button.  
%              ** To correct mistakes click the MIDDLE button. Backstepping
%                 is possible until the polyline runs out of points. 
%              ** To finish, click the RIGHT mouse button.
%
%      Usage : [x, y] = setpolygon;
%
%     Inputs : Optional!
% str1, str2 : Respectively strings indicating the physical units of
%                     the X- and Y- coordinates of the data plotted in the
%                     current axes, e.g. 'meters', 'nanoseconds', etc.
%
%    Outputs : 
%       x, y : The coordinates of the polygonal line in data units.
%
%       Note : The routine communicates with "setpoly" through global
%              variables. These are kept strictly within the scope of
%              "setpolygon.m". For details see the self-documentation of 
%              "setpoly" below.
%
%     Author : Andreas Tzanis
%              Department of Geophysics, 
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
global poly xo yo marks p3 xstr ystr
xstr = '';    ystr = '';     
% Check input arguments
if nargin >= 1,
    xstr = str1;
end
if nargin == 2,
    ystr = str2;
end
% Initialize output
x    = []; 
y    = []; 
poly = [];
% Get polygon
setpoly('begin');
% Assign output
if ~isempty(poly),
    x    = poly(:,1);
    y    = poly(:,2);
end
clear global xo yo marks p3 xstr ystr
return
%%% End function SETPOLYGON    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setpoly(action);
%
% SETPOLY : Define, on the current axes, a set of points [ x(n) y(n) ]
%           by pointing and clicking. 
%           ==> The polygonal line is clipped with respect to the axes
%               limits. 
%           ==> Accurate pointing is facilitated by a streaming display at
%               the bottom of the figure, of the current pointer position  
%               (coordinates) in data units. The display is erased as soon
%               as the program returns. 
%           ==> Work flow is controlled exclusively with the mouse: 
%               ** To set an x,y pair of coordinates click the LEFT mouse
%                  button.  
%               ** To reshape the polyline or correct mistakes backstep by
%                  clicking the MIDDLE button. Backstepping is possible
%                  until the polyline runs out of points. 
%               ** To finish, click the RIGHT mouse button.
%
%   Usage : setpoly('begin');
%
%  Inputs : 
%  action : Keyword determining the runtime behaviour of "setpoly". 
%           ==> This is a multiple-entry recursive program. Each entry is
%               activated by a keyword. The keywords are pre-programmed.
%               The only option available (and permisible) to the a user is
%               action = 'begin', which initializes the program, setting up
%               its interfaces on the current axes. 
%
% Outputs : None. The routine communicates its runtime parameters and the
%           results to the caller (driver) function "setpolygon.m" (or any
%           other caller function thereof), via global variables. These
%           variables are:
%           poly    : The polygonal line
%           xo, yo  : auxiliary variables
%           marks   : array of handles to the marker objects demarkating
%                     the vertices of the polygon (stars)
%              p3   : handle to the semi-transparent patch object outlining
%                     the area enclosed by the polygon.
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

global poly xo yo marks lin2 p3 xstr ystr

%%% -----------------------------------------------------------------------
%%% ENTRY "BEGIN": Setup the runtime environment 
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
        'WindowButtonMotionFcn','setpoly(''motion'');',...
        'WindowButtonDownFcn', ...
        'setpoly(''down''); ');
%
        set(gca,'nextplot','add')
% Interrupt execution of the caller function until setpoly exits
% gracefully
    waitfor(findobj('tag','newx'))
end

%%% -----------------------------------------------------------------------
%%% ENTRY "MOTION": The WindowbuttonMotionFcn - Set cursor shape and
%%% streaming display of cursor position 
if strcmp(action,'motion'),
    nmark = size(poly,1);
    if exist('p3') & ishandle(p3),
        set(p3,'visible','off')
    end
    xl = get(gca,'xlim');               % x-axis limits
    yl = get(gca,'ylim');               % y-axis limits
    cp = get(gca,'currentpoint');
    cx = cp(1,1);
    cy = cp(1,2);
%    online = cy > yl(1) && cy < yl(2) && cx > xl(1) && cx < xl(2) ;
%    if online,
        set(gcf,'Pointer','crosshair');
        str = ['X = ' num2str(cx) ' ' xstr ];
        set(findobj('tag','newx'),'string',str)
        str = ['Y = ' num2str(cy) ' ' ystr];
        set(findobj('tag','newy'),'string',str)
%    elseif ~online,
%        set(gcf,'Pointer','arrow');
%        set(findobj('tag','newx'),'string','Out of range' )
%        set(findobj('tag','newy'),'string','Out of range' )
%    end
    if nmark ==1,
        p3 = line([poly(:,1); cp(1,1)],[poly(:,2); cp(1,2)],'color','k');
    elseif nmark == 2,
        p3 = patch([poly(:,1); cp(1,1)],[poly(:,2); cp(1,2)],'k',...
            'edgealpha',0.5, 'facealpha',0.1);
    elseif nmark > 2,
        set(p3,'xdata',[poly(:,1); cp(1,1)],'ydata',[poly(:,2); cp(1,2)],...
            'visible','on');
    end
end

%%% -----------------------------------------------------------------------
%%% ENTRY "DOWN": The WindowbuttonDownFcn - Define the polygonal line 
if strcmp(action,'down'),
%    if exist('p3') & ishandle(p3),
%        set(p3,'visible','off')
%    end
    nmark = size(poly,1);
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    dx = diff(xl);
    dy = diff(yl);
% Set threshold for deciding proximity to axes limits
    dtest = 0.01;
    
%%% Left mouse button clicked: Set a point 
    if strcmp(get(gcf,'selectiontype'),'normal'),
        
        cp = get(gca ,'currentpoint');       % Get coordinates
        xp = cp(1,1);
        yp = cp(1,2);
%%% If xp or yp are very close to, or beyond axes limits, snap them to 
%%% axes limits by linear interpolation.
%        slope = 0;
        if nmark > 0,
% check left        
        if xp < xl(1) && xo > xl(1) ,
            slope = (yo - yp)/(xo - xp);
            yp = yp + slope*(xl(1)-xp);
            xp = xl(1);
        end
        if xp > xl(1) && xo < xl(1) ,
            slope = (yp - yo)/(xp - xo);
            yp = yp + slope*(xl(1)-xp);
            xp = xl(1);
        end
        d = abs((xl(1)-xp)/dx);
        if d < dtest,
            xp = xl(1);
        end
% check right        
        if xp > xl(2) && xo < xl(2),
            slope = (yo - yp)/(xo - xp);
            yp = yp + slope*(xl(2)-xp);
            xp = xl(2);
        end
        if xp < xl(2) && xo > xl(2),
            slope = (yp - yo)/(xp - xo);
            yp = yp + slope*(xl(2)-xp);
            xp = xl(2);
        end
        d = abs((xl(2)-xp)/dx);
        if d < dtest,
            xp = xl(2);
        end
% check top
        if yp < yl(1) && yo > yl(1),
            den = (xo - xp);
            if den~=0,
                slope = (yo - yp)/den;
                y0 = yp + slope*(xl(1)-xp);
                xp = (yl(1) - y0)/slope;
                yp = yl(1);
            else
                yp = yl(1);
            end
        end
        if yp > yl(1) && yo < yl(1),
            den = (xp - xo);
            if den~=0,
                slope = (yp - yo)/den;
                y0 = yp + slope*(xl(1)-xp);
                xp = (yl(1) - y0)/slope;
                yp = yl(1);
            else
                yp = yl(1);
            end
        end
        d = abs((yl(1)-yp)/dy);
        if d < dtest,
            yp = yl(1);
        end
% check bottom
        if yp > yl(2) && yo < yl(2),
            den = (xo - xp);
            if den~=0,
                slope = (yo - yp)/den;
                y0 = yp + slope*(xl(1)-xp);
                xp = (yl(2) - y0)/slope;
                yp = yl(2);
            else
                yp = yl(2);
            end
        end
        if yp < yl(2) && yo > yl(2),
            den = (xp - xo);
            if den ~= 0,
                slope = (yp - yo)/den;
                y0 = yp + slope*(xl(1)-xp);
                xp = (yl(2) - y0)/slope;
                yp = yl(2);
            else
                yp = yl(2);
            end
        end
        d = abs((yl(2)-yp)/dy);
        if d < dtest,
            yp = yl(2);
        end
  end
%%% Set and draw the new vertex
        xo = cp(1,1);
        yo = cp(1,2);
        poly = [poly;  ...
                xp yp];
        nmark = nmark + 1;
        marks(nmark) = plot(poly(nmark,1),poly(nmark,2),'hk',...
            'markerfacecolor','w', 'MarkerSize',10,'LineWidth',1.0);
    end

%%% Right mouse button clicked: Clean up and Return to Caller 
    if strcmp(get(gcf,'selectiontype'),'alt'),
%%% Clip the plygon before returning
         if size(poly,1) >= 1,
             ii    = find(poly(:,1) < xl(1));
             poly(ii,1) = xl(1);
             ii    = find(poly(:,1) > xl(2));
             poly(ii,1) = xl(2);
             ii    = find(poly(:,2) < yl(1));
             poly(ii,2) = yl(1);
             ii    = find(poly(:,2) > yl(2));
             poly(ii,2) = yl(2);
         end
         set(gcf,'Backingstore','off','WindowButtonMotionFcn',''); 
         set(gcf,'WindowButtonDownFcn', '');
         set(gcf,'pointer','arrow')
         delete(findobj('tag','newx'))
         delete(findobj('tag','newy'))
         delete(marks(1:nmark))
         delete(p3)
    end
    
%%% Middle button of Left & Right clicked: Mistake; Clear last vertex
    if strcmp(get(gcf,'selectiontype'),'extend'),
         if nmark <= 0,
             return
         elseif nmark == 1,
             poly = [];
             delete(marks(1))
         else
             poly = poly(1:nmark-1,:);
             xo = poly(size(poly,1),1);
             yo = poly(size(poly,1),2);
             delete(marks(nmark))
             set(p3,'xdata',poly(:,1),'ydata',poly(:,2),'visible','on');
             drawnow
         end
    end
end
return
%%% End function SETPOLY    