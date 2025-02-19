function showinfo(DATA)
%
% SHOWINFO : 1. Assembles a column cell array of strings with information
%               about the input data structure and its processing history.
%            2. Assembles information about the velocity models used for
%               interpretation, or derived from data analysis
%            3. Displays this information on the MATGPR information window
%               using the "welcome.m" utility.   
%
%    Usage : showinfo(DATA)
%
%   Inputs : DATA - the IPD or OPD data structures of MATGPR
%
% Requires : welcome.m
%
%   Author : Andreas Tzanis,
%            Department of Geophysics, 
%            University of Athens
%            atzanis@geol.uoa.gr
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

% Initialize global variables
global VS

clb=cell(1);
txt = ['Data Origin        = ' DATA.origin];
clb(1) = cellstr(txt);
%%% Data file path and name
txt = ['Input file name    = ' [DATA.pname DATA.fname]];
clb(length(clb)+1) = cellstr(txt);
%%% DZT header gain
if ~isempty(DATA.DZThdgain),
    txt = ['DZT header gain    = Loaded'];
else
    txt = ['DZT header gain    = NA'];
end
clb(length(clb)+1) = cellstr(txt);
%%% # traces
txt = ['Number of traces   = ' num2str(DATA.ntr)];
clb(length(clb)+1) = cellstr(txt);
%%% trace spacing
if ~isempty(DATA.dx),
    txt = ['Trace spacing      = ' num2str(DATA.dx,4) ... 
        'm;  (Section length = ' num2str(abs((DATA.ntr-1)*DATA.dx),4) 'm)'];
else
    txt = ['Trace spacing      = Unequally spaced'];
end
clb(length(clb)+1) = cellstr(txt);
%%% trace coordinates
if ~isempty(DATA.xyz.Tx) && ...
        isempty(find(size(DATA.xyz.Tx)~=size(DATA.xyz.Rx))),
    txt = ['x,y,z data         = Tx:[' num2str(size(DATA.xyz.Tx,1)) ...
          ' x ' num2str(size(DATA.xyz.Tx,2)) '], Rx:[' ...
          num2str(size(DATA.xyz.Rx,1)) ' x ' ...
          num2str(size(DATA.xyz.Rx,2)) '], OK and loaded'];
else
    txt = ['x,y,z data         = NA'];
end
clb(length(clb)+1) = cellstr(txt);
%%% Marker information
if ~isempty(DATA.markertr),
    [nmrk, nxyz] = size(DATA.markertr);
    if nxyz == 1,
        txt = ['Marker Information = ' num2str(nmrk) ' x ' ...
            num2str(nxyz) ', Marker ID numbers only'];
    elseif nxyz == 4,
        txt = ['Marker Information = ' num2str(nmrk) ' x ' ...
            num2str(nxyz) ', Marker ID + Coordinates'];
    else
        txt = ['Marker Information = Unknown data structure, please check'];
    end
else
    txt = ['Marker Information = NA'];
end
clb(length(clb)+1) = cellstr(txt);
%%% # time samples 
txt = ['Samples per trace  = ' num2str(DATA.ns)];
clb(length(clb)+1) = cellstr(txt);
%%% Sampling rate
txt = ['Sampling interval  = ' num2str(DATA.dt,4) 'ns;  (Time window = ' ... 
        num2str((DATA.ns-1)*DATA.dt,5) 'ns)'];
clb(length(clb)+1) = cellstr(txt);
%%% Signal position
if ~isempty(DATA.sigpos),
    txt = ['Signal position    = ' num2str(DATA.sigpos,5) ' ns' ];
else
    txt = ['Signal position    = NA']; 
end;
clb(length(clb)+1) = cellstr(txt);
%%% Antenna used
txt = ['Antenna            = ' DATA.Antenna];
clb(length(clb)+1) = cellstr(txt);
%%% Tx-Rx separation
txt = ['Antenna offset     = ' num2str(DATA.TxRx) 'm'];
clb(length(clb)+1) = cellstr(txt);
%%% Halfspace velocity if estimated
if isfield(VS,'v1d') && ~isempty(VS.v1d) && size(VS.v1d,1)==1,
    txt = ['Halfspace velocity = ' num2str(VS.v1d(1),4) 'm/ns'];
    clb(length(clb)+1) = cellstr(txt);
elseif  isfield(VS,'v1d') && ~isempty(VS.v1d) && size(VS.v1d,1) > 1,
%%% Current 1-D velocity model, if existing
    txt = ['1-D Velocity model = Loaded (' num2str(size(VS.v1d,1)-1) ... 
            ' layers + basal halfspace)'];
    clb(length(clb)+1) = cellstr(txt);
end
%%% Current 2-D velocity model, if existing
if isfield(VS,'v2d') && ~isempty(VS.v2d),
    if iscell(VS.v2d),
            txt = ['2-D Velocity model = {' num2str(size(VS.v2d,1)) ' x ' ... 
                 num2str(size(VS.v2d,2)) '} Cell Array with contents:'];
            clb(length(clb)+1) = cellstr(txt);
        if size(VS.v2d,2) >= 1,
            txt = ['                     V0 = [' num2str(size(VS.v2d{1},1)) ...
                ' x '  num2str(size(VS.v2d{1},2)) '] Non-Dispersive term']; 
            clb(length(clb)+1) = cellstr(txt);
        end
        if size(VS.v2d,2) >= 2,
            txt = ['                      Q = [' num2str(size(VS.v2d{2},1)) ...
                ' x '  num2str(size(VS.v2d{2},2)) '] Quality factor']; 
            clb(length(clb)+1) = cellstr(txt);
        end
        if size(VS.v2d,2) == 3,
            txt = ['                     fc = ' num2str(VS.v2d{3}) ...
                ' MHz Antenna Frequency']; 
            clb(length(clb)+1) = cellstr(txt);
        end
    else
        txt = ['2-D Velocity model = Loaded [' num2str(size(VS.v2d,1)) ... 
            ' x ' num2str(size(VS.v2d,2)) '] Non-Dispersive term'];
       clb(length(clb)+1) = cellstr(txt);
    end
end
%%% Processing history of current data set
if size(DATA.history,1) > 1,
    txt = ['Processing History = ' char(DATA.history{2})];
    clb(length(clb)+1) = cellstr(txt);
    for j=3:size(DATA.history,1),
        txt = ['                     ' char(DATA.history{j})];
        clb(length(clb)+1)  = cellstr(txt);
    end
else
    txt = ['Processing History = ' char(DATA.history) ];
    clb(length(clb)+1)  =   cellstr(txt);
end

%%%%%   Display the information
welcome(clb);    

return
