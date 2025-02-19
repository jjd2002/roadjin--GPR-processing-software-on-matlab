function matgprwindows(action,object)
%
% MATGPRWINDOWS : The MATGPR fast window switching service. 
%                 This service generates a custom "Windows" menu, which
%                 provides for fast swithcing between MATGPR windows. The
%                 "matgprwindows" function will create and display windows
%                 EXCLUSIVE to MATGPR and, moreover, only those registered
%                 with the service (see the 'setup' block below). These
%                 would normally be windows with information deserving to
%                 be left on the screen for future reference (and consume
%                 memory), and not trivial stuff, as for instance dialog
%                 boxes. 
%             ==> This program is part of the MATGPR suite.
%
%         Usage : The function will normally not be used through the
%                 command line. Although it is not neccessary, it should
%                 rather be used in the CreateFcn and DeleteFcn callbacks
%                 (on creation or destruction of a window). For example:
%                 ....
%                 figure('name','My Figure','Tag','MyNewFigure', ...
%                        'menubar','none', ... , etc. , ...
%                        'CreateFcn',['figuretools; imagecolors; ' ... 
%                                     'matgprwindows(''setup''); ' ...
%                                     'matgprwindows(''updateadd'');'],...
%                        'DeleteFcn',['matgprwindows(''updateremove'', '...  
%                                     'findobj('Tag','MyNewFigure'));']);
%                  ....
%                  Keywords : 
%                     setup : Initialize the service in the current figure
%                 updateadd : After creation, add a figure's name to the
%                             "Windows" menus of all open figures (windows) 
%              updateremove : Upon deletion, remove a figure's name from
%                             the "Windows" menus of the remaining open
%                             figures   
%
%        Author : Andreas Tzanis
%                 Department of Geophysics, 
%                 University of Athens
%                 atzanis@geol.uoa.gr
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

global winarray
%%% Setup the Windows menu at the Current Figure (usually upon creation of
%%% the figure).
if strcmp(action,'setup'),
    if nargin == 1,
        object = gcf;
    end
    winarray = [];
    if isempty(findobj(object,'Tag','winmenu')),
        wmenu = uimenu(object,'Label', '窗口跳转', 'Tag', 'winmenu');
    else
        wmenu = findobj(object,'Tag','winmenu');
    end
%%% Search for open matGPR windows and add them to the list    
    if ~isempty(findobj('Tag','fi0')),
        winarray = [winarray; findobj('Tag','fi0')];
        uimenu(wmenu,'Label','软件主页面',...
            'Callback','figure(findobj(''Tag'',''fi0''));');
    end
    if ~isempty(findobj('tag','headerinfofigure')),
        winarray = [winarray; findobj('tag','headerinfofigure')];
        uimenu(wmenu,'Label','Header Info',...
            'Callback','figure(findobj(''tag'',''headerinfofigure''));');
    end
    if ~isempty(findobj('tag','datafigure')),
        winarray = [winarray; findobj('tag','datafigure')];
        uimenu(wmenu,'Label','GPR 数据',...
            'Callback','figure(findobj(''tag'',''datafigure''));');
    end
    if ~isempty(findobj('tag','procdatafigure')),
        winarray = [winarray; findobj('tag','procdatafigure')];
        uimenu(wmenu,'Label','处理过的 GPR 数据',...
            'Callback','figure(findobj(''tag'',''procdatafigure''));');
    end
    if ~isempty(findobj('tag','viewindatatraces')),
        winarray = [winarray; findobj('tag','viewindatatraces')];
        uimenu(wmenu,'Label','单道查看',...
            'Callback','figure(findobj(''tag'',''viewindatatraces''));');
    end
    if ~isempty(findobj('tag','viewoutdatatraces')),
        winarray = [winarray; findobj('tag','viewoutdatatraces')];
        uimenu(wmenu,'Label','处理过的单道查看',...
            'Callback','figure(findobj(''tag'',''viewoutdatatraces''));');
    end
    % if ~isempty(findobj('tag','viewindataspectra')),
    %     winarray = [winarray; findobj('Tag','viewindataspectra')];
    %     uimenu(wmenu,'Label','Spectra Viewer',...
    %         'Callback','figure(findobj(''tag'',''viewindataspectra''));');
    % end
    % if ~isempty(findobj('tag','viewoutdataspectra')),
    %     winarray = [winarray; findobj('Tag','viewoutdataspectra')];
    %     uimenu(wmenu,'Label','Processed Spectra Viewer',...
    %         'Callback','figure(findobj(''tag'',''viewoutdataspectra''));');
    % end
    % if ~isempty(findobj('tag','attenuationfigure')),
    %     winarray = [winarray; findobj('Tag','attenuationfigure')];
    %     uimenu(wmenu,'Label','Attenuation Characteristics',...
    %         'Callback','figure(findobj(''tag'',''attenuationfigure''));');
    % end
% This is a dialog box that could be forgotten open on occasion
%    if ~isempty(findobj('tag','edscanaxgui')),
%        winarray = [winarray; findobj('Tag','edscanaxgui')];
%        uimenu(wmenu,'Label','Edit Scan Axis GUI',...
%            'Callback','figure(findobj(''tag'',''edscanaxgui''));');
%    end
    % if ~isempty(findobj('tag','markercheckfigure')),
    %     winarray = [winarray; findobj('Tag','markercheckfigure')];
    %     uimenu(wmenu,'Label','Marker Trace Editor',...
    %         'Callback','figure(findobj(''tag'',''markercheckfigure''));');
    % end
    if ~isempty(findobj('tag','velanfigure')),
        winarray = [winarray; findobj('Tag','velanfigure')];
        uimenu(wmenu,'Label','速度拟合',...
            'Callback','figure(findobj(''tag'',''velanfigure''));');
    end
% This is a dialog box that could be forgotten open on occasion
%    if ~isempty(findobj('tag','staticgui')),
%        winarray = [winarray; findobj('Tag','staticgui')];
%        uimenu(wmenu,'Label','Statics GUI',...
%            'Callback','figure(findobj(''tag'',''staticgui''));');
%    end
    % if ~isempty(findobj('tag','modelfigure')),
    %     winarray = [winarray; findobj('Tag','modelfigure')];
    %     uimenu(wmenu,'Label','Model Builder',...
    %         'Callback','figure(findobj(''tag'',''modelfigure''));');
    % end
    % if ~isempty(findobj('tag','badtrfigure')),
    %     winarray = [winarray; findobj('Tag','badtrfigure')];
    %     uimenu(wmenu,'Label','去除坏道',...
    %         'Callback','figure(findobj(''tag'',''badtrfigure''));');
    % end
    % if ~isempty(findobj('tag','isofigure')),
    %     winarray = [winarray; findobj('Tag','isofigure')];
    %     uimenu(wmenu,'Label','Isosurface Display',...
    %         'Callback','figure(findobj(''tag'',''isofigure''));');
    % end
    % if ~isempty(findobj('tag','iso3dcontrols')),
    %     winarray = [winarray; findobj('Tag','iso3dcontrols')];
    %     uimenu(wmenu,'Label','Iso-Surface Display Controls',...
    %         'Callback','figure(findobj(''tag'',''iso3dcontrols''));');
    % end
    % if ~isempty(findobj('tag','slice3dfigure')),
    %     winarray = [winarray; findobj('Tag','slice3dfigure')];
    %     uimenu(wmenu,'Label','3D Slice Display',...
    %         'Callback','figure(findobj(''tag'',''slice3dfigure''));');
    % end
    % if ~isempty(findobj('tag','slice3dcontrols')),
    %     winarray = [winarray; findobj('Tag','slice3dcontrols')];
    %     uimenu(wmenu,'3-D Slice Display Controls',...
    %         'Callback','figure(findobj(''tag'',''slice3dcontrols''));');
    % end
%%%    
%%% IF YOU WANT TO REGISTER A NEW WINDOW WITH THE SERVICE, 
%%% INSERT AFTER THIS POINT AND BEFORE RETURN!     
%%%    
    return
end

%%% After creation, add a figure's name to the Windows menu of all open
%%% figures  
if strcmp(action,'updateadd'),
    winarray2 = winarray;
    for i=1:length(winarray2),
        if ishandle(winarray2(i)),
            cc = get(winarray2(i),'children');
            wm = findobj(cc,'tag','winmenu');
%            delete(wm);
            cc2 = get(wm, 'children');
            delete(cc2)
            matgprwindows('setup',winarray2(i))
        end
    end
end

%%% Upon deletion, remove a figure's name from the Windows menu of all open
%%% figures  
if strcmp(action,'updateremove'),
    ii = find(winarray == object);
    winarray2 = winarray(  setxor([1:1:length(winarray)], ii  ) );
    winarray1 = flipud(winarray);
    ii = find(winarray1 == object);
    for i=1:length(winarray1),
        if ishandle(winarray1(i)),
            cc = get(winarray1(i),'children');
            wm = findobj(cc,'tag','winmenu');
            cc2 = get(wm, 'children');
            delete(cc2(ii))
        end
    end
    winarray = winarray2;
end
return