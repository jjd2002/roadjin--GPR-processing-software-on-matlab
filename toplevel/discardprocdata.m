function OUTDATA = discardprocdata
%
% DISCARDPROCDATA: Utility to erase the current output data from memory. 
%
% Usage    : OUTDATA = discardprocdata
%
% Requires : initdatastr.m
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
%  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA
%

%%% Re-initialize the Output Data structure
OUTDATA = initdatastr;
%%% Clear figures
if ishandle(findobj('tag','procdatafigure')),    % clear out data figure
    delete(findobj('tag','procdatafigure'))
end
if ishandle(findobj('tag','viewoutdatatraces')), % clear out trace viewer 
    delete(findobj('tag','viewoutdatatraces'))
end
if ishandle(findobj('tag','viewoutdataspectra')),% clear out spectra viewer
    delete(findobj('tag','viewoutdataspectra'))
end
return