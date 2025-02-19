function DATA = initdatastr()

% INITDATASTR: Initializes the Input and Output matGPR data structures
%
% Author : Andreas Tzanis
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
DATA.origin     = [];         % String, specifies the type of system that 
                              % collected the original data set
                              % (manufacturer and monostatic of bistatic).  
DATA.pname      = [];         % String, path of the data file
DATA.fname      = [];         % String, name of the data file
DATA.d          = [];         % 2-D array, matrix of input data 
DATA.ns         = [];         % Scalar, No of samples per column 
                              % (trace / scan) in DATA.d 
DATA.dt         = [];         % Scalar, temporal sampling rate of the 
                              % columns (traces) of DATA.d
DATA.tt2w       = [];         % Vector, 2-way traveltime
DATA.sigpos     = [];         % Scalar, signal position (time-zero 
                              % determined during data acquisition)
DATA.dz         = [];         % Scalar, spatial (depth) sampling rate of 
                              % the columns (traces) of DATA.d - applies to
                              % depth migrated data
DATA.z          = [];         % Vector, depth 
DATA.zlab       = [];         % String, label of the vertical axis, used  
                              % for display 
DATA.ntr        = [];         % Scalar, No of traces (samples per row) 
                              % in DATA.d
DATA.dx         = [];         % Scalar, trace spacing
DATA.x          = [];         % Vector, horizontal coordinates of the rows  
                              % of DATA.d (scan line). Can be either number 
                              % of traces for unequally spaced data, or 
                              % distances from the start of the scanline
DATA.xlab       = [];         % String, label of the horizontal axis, used
                              % for display
DATA.markertr   = [];         % Vector, ID numbers of marker traces 
DATA.xyz.Tx     = [];         % Vector, coordinates of the transmitter 
                              % location in a local frame of reference 
DATA.xyz.Rx     = [];         % Vector, coordinates of the receiver 
                              % location in a local frame of reference 
DATA.TxRx       = [];         % Scalar, antenna offset (Tx - Rx distance)
DATA.Antenna    = [];         % String, antenna name / designation 
DATA.DZThdgain  = [];         % Vector, variable gain settings used for
                              % data acquisition with GSSI systems (DZT)
DATA.TimesSaved = 0;          % Scalar, the number of times that a "save data"
                              % operation has been done during the current
                              % session. 
DATA.comments   = [];         % Header comments of original data set
DATA.history = cellstr('Raw Data');     % Cell array of strings to hold the processing history
                                        % of the current data set
if ~isstruct(ENVAR) || ~isfield(ENVAR,'CURRENT_VERSION'),
    DATA.matgprversion = 'Release 2, 2008';
else
    DATA.matgprversion  = ENVAR.CURRENT_VERSION;  % MATGPR version
end
return
