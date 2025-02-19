function [IPD, OPD] = undoxx(action)
%
%UNDOXX : UNDO previous processing steps / Restore data to previous states
%
% Usage : undo('initialize') to set up the GUI to display the Processing
%         History and allow sellection of an undo point. Afterwards, the
%         routine recurses on itself ([IPD, OPD] = undo('process') 
%         to perform the operation. 
%
%Author : Andreas Tzanis,
%         Department of Geophysics, 
%         University of Athens
%         atzanis@geol.uoa.gr
%
% Copyright (C) 2008, Andreas Tzanis. All rights reserved.
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

% Startup the GUI to assemble profiles for concatenation
if strcmp(action,'initialize'),
    
    % Get Processing History
    UNDO = get(findobj('tag','fi0'),'userdata');
    if isempty(UNDO),
        disp('UNDO Warning > The Undo buffer is empty! No action taken!') 
        delete(findobj('tag','undo_gui'));
        return
    end
    % Create GUI interface
    figure('Name','Undo', 'Tag', 'undo_gui', ...
        'NumberTitle','off', 'Position',[200 450 500 300 ], ...
        'MenuBar','none', ...
        'Color',[0.76 0.86 0.69]);         % 'Color',[0 0.5 0.5]);
     axis off
     fcolor = [0.76 0.87 0.78];            % fcolor = get(gcf,'color');
    % Display Processing History
    text('Color',[0 0 0 ],  'Units','normalized', ...
        'Position', [0.0 0.95 0 ], 'FontSize',10 , ...
        'FontWeight','bold', ...
        'String','Choose an Undo Point in the Processing History:');
    uicontrol('Style', 'listbox','Tag','undobox', ...
        'Units','normalized', ...
        'Position',[0.05 0.1 0.75 0.73], ...
        'backgroundcolor', fcolor, ...
        'string', strvcat(UNDO{:,2}), ...
        'fontsize', 10); %#ok<VCAT>
    % Proceed ...
    uicontrol('Style','Pushbutton','tag','gobutn', ...
        'Units','normalized', 'Position',[.82 .475 .15 .1 ], ...
        'BackGroundColor', [0.6 1 0], ...
        'String','GO', ...
        'Callback', '[IPD, OPD] = undoxx(''process''); ');
    % Cancel ...
    uicontrol('Style','Pushbutton','Units','normalized', ... 
        'Position',[.82 .2 .15 .1 ], ...
        'BackGroundColor',[1 0.2 0], ...
        'String','Cancel', 'Callback', ...
        'delete(findobj(''tag'',''undo_gui'')); return; ');
    return
end

% File assembly collected - now concatenate
if strcmp(action,'process'),

    UNDO  = get(findobj('tag','fi0'),'userdata');
    nundo = get(findobj('tag','undobox'),'value');
    IPD = UNDO{nundo,1};
    if nundo > 1,
        UNDO = UNDO(1:nundo-1,:);
    elseif nundo == 1,
        UNDO = [];
    end
    set(findobj('tag','fi0'),'userdata',UNDO);

    % Output data no longer needed
    OPD = discardprocdata;
    %%% Clear data related figures
    if ishandle(findobj('tag','procdatafigure')),         % out data figure
        delete(findobj('tag','procdatafigure'))
    end
    if ishandle(findobj('tag','viewoutdatatraces')),     % out trace viewer
        delete(findobj('tag','viewoutdatatraces'))
    end
    if ishandle(findobj('tag','viewoutdataspectra')),  % out spectra viewer
        delete(findobj('tag','viewoutdataspectra'))
    end
    if ishandle(findobj('tag','viewindatatraces')),       % in trace viewer
        delete(findobj('tag','viewindatatraces'))
    end
    if ishandle(findobj('tag','viewindataspectra')),    % in spectra viewer
        delete(findobj('tag','viewindataspectra'))
    end
    %%% Display information of updated data
    showinfo(IPD)
    %%% Update the "Current Input Data" figure
    viewdata(IPD.x,IPD.tt2w,IPD.d,'indata',IPD.xlab,IPD.zlab);

% Wrap up and exit
    delete(findobj('tag','undo_gui')); 
    return
end


