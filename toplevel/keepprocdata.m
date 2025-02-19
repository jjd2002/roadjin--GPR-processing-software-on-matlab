function [IPD, OPD] = keepprocdata(IPD, OPD)
%
% KEEPPROCDATA : Utility to replace the current input data (before a
% processing step)with the current output data (after the processing step).  
%
% Usage    : [IPD, OPD] = keepprocdata(IPD, OPD)
%
% Requires : showinfo.m
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
if isempty(OPD.d),        % Check for processed data
    erh = errordlg('There is NO OUTPUT DATA - Nothing to do!'...
        ,'KEEPPROCDATA: ERROR');
    uiwait(erh); 
    return
end
%%% Check if data has been depth migrated - they cannot be held, only saved
if ~isempty(findstr(OPD.zlab,'Depth')),
    txt = cell(1);
    txt(1) = cellstr('This is the end of a line! Depth migrated data cannot be held!');
    txt(2) = cellstr('Depth migrated data can only be Saved, Printed or Discarded!'  );
    msg = msgbox(txt,'KEEPPROCDATA: NOTIFICATION','help');
    uiwait(msg); 
    return
end

%%% Load IPD onto the UNDO buffer and store
UNDO   = get(findobj('tag','fi0'),'userdata');
nundos = size(UNDO,1);
%%% To conserve memory, allow only 4 levels of undo/restore actions
if nundos >= ENVAR.undolevels,
    UNDO = [UNDO(2:nundos,:); {IPD, IPD.history{size(IPD.history,1)}}];
else
    UNDO = [UNDO; {IPD, IPD.history{size(IPD.history,1)}}];
end
set(findobj('tag','fi0'),'userdata',UNDO);
set(findobj('tag','undoxx'),'Enable','on');
set(findobj('tag','undoundo'),'Enable','on');

%%% Replace Current Input data
IPD  =  OPD;

%%% Re-initialize the Output data structure
OPD = initdatastr;

%%% Clear processed data figures
if ishandle(findobj('tag','procdatafigure')),    % Clear the out data figure
    delete(findobj('tag','procdatafigure'))
end
if ishandle(findobj('tag','viewoutdatatraces')), % Clear out trace viewer 
    delete(findobj('tag','viewoutdatatraces'))
end
if ishandle(findobj('tag','viewoutdataspectra')),% Clear out spectra viewer
    delete(findobj('tag','viewoutdataspectra'))
end
if ishandle(findobj('tag','viewindatatraces')),  % Clear in trace viewer 
    delete(findobj('tag','viewindatatraces'))
end
if ishandle(findobj('tag','viewindataspectra')), % Clear in spectra viewer 
    delete(findobj('tag','viewindataspectra'))
end

%%% Display information of updated data 
showinfo(IPD)

%%% Update the "Current Input Data" figure    
viewdata(IPD.x,IPD.tt2w,IPD.d,'indata',IPD.xlab,IPD.zlab); 
return
