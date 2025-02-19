function SUHDR = initSUheader(DATA)
%
% INITSUHEADER : Initializes the SU trace header structure SUHDR and 
% assigns some standard parameters passed with the MATGPR data structure
% DATA (either IPD or OPD).  
%
% Usage  : initSUheader(DATA)
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

% Date data are being created or exported
epoch     = datevec(now);
dayofyear = datenum(epoch(1:3)) - datenum([epoch(1),1,1]);
%Alternatively use decyear.m
%dayofyear = decyear(epoch(1:3)) - decyear([epoch(1),1,1]);
% Initialize header
SUHDR.tracl   = 0;                      % Trace sequence number within line
SUHDR.tracr   = 0;                      % Trace sequence number within reel
SUHDR.fldr    = 1;                      % Field record number
SUHDR.tracf   = 0;                      % Trace number within field record
SUHDR.ep      = 0;                      % Energy source point number
SUHDR.cdp     = 0;                      % CDP ensemble number 
SUHDR.cdpt    = 0;                      % Trace number within CDP ensemble 
SUHDR.trid    = 7182;                   % Trace identification code:
                                        % MATGPR ID code = 7182 (ASCII codes
                                        % [71 82] == 'GR' for GeoRadar) 
SUHDR.nvs     = 0;                      % Number of vertically summed traces 
                                        % (see vscode in reel header 
                                        % Structure - PULSE EKKO standard)
SUHDR.nhs     = 0;                      % number of horizontally summed traces 
                                        % (see vscode in reel header 
                                        % structure - PULSE EKKO standard) 
SUHDR.duse    = 1;                      % Data use: 1 = production 2 = test 
if ~isempty(DATA.TxRx),
    SUHDR.offset = DATA.TxRx*1000;      % Distance from source point to 
else                                    % receiver group. 
    SUHDR.offset = 0;                   % (negative if opposite to direction 
end                                     % in which the line was shot)
SUHDR.gelev   = 0;                      % Receiver group elevation from sea 
                                        % level (above sea level is positive)
SUHDR.selev   = 0;                      % Source elevation from sea level 
                                        % (above sea level is positive) 
SUHDR.sdepth  = 0;                      % Source depth below surface (positive) 
SUHDR.gdel    = 0;                      % Datum elevation at receiver group
SUHDR.sdel    = 0;                      % Datum elevation at source 
SUHDR.swdep   = 0;                      % Water depth at source 
SUHDR.gwdep   = 0;                      % Water depth at receiver group
SUHDR.scalel  = -3;                     % Scale factor for previous 7 
                                        % entries with value +/- 10 to the power 
                                        % power 0, 1, 2, 3, or 4 (if  
                                        % posituve multiply, if negative divide) 
SUHDR.scalco  = -3;                     % Scale factor for next 4 entries 
                                        % with value +/- 10 to the power 
                                        % 0, 1, 2, 3, or 4 (if positive,
                                        % multiply, if negative divide)  
SUHDR.sx      = 0;                      % X source coordinate
SUHDR.sy      = 0;                      % Y source coordinate 
SUHDR.gx      = 0;                      % X group coordinate
SUHDR.gy      = 0;                      % Y group coordinate
SUHDR.counit  = 1;                      % Coordinate units code: for 
                                        % previous four entries 
                                        % 1  = length (meters or feet) 
                                        % 2  = seconds of arc (the 
                                        % X -longitude Y -latitude, 
                                        % positive to the east of Greenwich
                                        % or north of the equator.  
SUHDR.wevel   = 0;                      % Weathering velocity
SUHDR.swevel  = 0;                      % Subweathering velocity
SUHDR.sut     = 0;                      % Uphole time at source
SUHDR.gut     = 0;                      % Uphole time at receiver group
SUHDR.sstat   = 0;                      % Source static correction
SUHDR.gstat   = 0;                      % Group static correction
SUHDR.tstat   = 0;                      % Total static applied
SUHDR.laga    = 0;                      % Lag time A, time in ms between 
                                        % end of 240-byte trace
                                        % identification header and time
                                        % break, positive if time break
                                        % occurs after end of header,  
                                        % time break is defined as the
                                        % initiation pulse which maybe
                                        % recorded  on an auxiliary trace
                                        % or as otherwise specified by the
                                        % recording system.
SUHDR.lagb = 0;                         % Lag time B, time in ms between 
                                        % the time break and the initiation
                                        % time of the energy source, may be
                                        % positive or negative. 
SUHDR.delrt   = DATA.sigpos*1000;       % Delay recording time, time in ms 
                                        % between initiation time of energy
                                        % source and time when recording of
                                        % data samples begins (for deep
                                        % water work if recording does not
                                        % start at zero time).  
SUHDR.muts    = 0;                      % Mute time--start
SUHDR.mute    = 0;                      % Mute time--end
if ~isempty(DATA.ns),
    SUHDR.ns  = DATA.ns;                % Number of samples in this trace
else
    SUHDR.ns  = 0;
end
if ~isempty(DATA.dt),
    SUHDR.dt  = DATA.dt*1000;           % Sampling interval; in micro-secs
else
    SUHDR.dt  = 0;
end
SUHDR.gain    = 0;                      % Gain type of field instruments 
                                        % code: 1 = fixed, 2 = binary, 
                                        % 3 = floating point, 4 - N  =
                                        % optional use.   
SUHDR.igc     = 0;                      % Instrument gain constant
SUHDR.igi     = 0;                      % Instrument early or initial gain
SUHDR.corrsu  = 0;                      % Correlated:  1 = no;  2 = yes 
SUHDR.sfs     = 0;                      % Sweep frequency at start
SUHDR.sfe     = 0;                      % Sweep frequency at end
SUHDR.slen    = 0;                      % Sweep length in ms 
SUHDR.styp    = 0;                      % Sweep type code: 1 = linear; 
                                        % 2 = parabolic; 3 = exponential; 
                                        % 4 = other 
SUHDR.stas    = 0;                      % Sweep trace taper length at start in ms
SUHDR.stae    = 0;                      % Sweep trace taper length at end in ms
SUHDR.tatyp   = 0;                      % Taper type: 1 =linear, 2 =cos^2, 3 =other
SUHDR.afilf   = 0;                      % Alias filter frequency if used
SUHDR.afils   = 0;                      % Alias filter slope 
SUHDR.nofilf  = 0;                      % Notch filter frequency if used 
SUHDR.nofils  = 0;                      % Notch filter slope 
SUHDR.lcf     = 0;                      % Low cut frequency if used
SUHDR.hcf     = 0;                      % High cut frequncy if used
SUHDR.lcs     = 0;                      % Low cut slope
SUHDR.hcs     = 0;                      % High cut slope
SUHDR.year    = epoch(1);               % Year data recorded 
SUHDR.day     = dayofyear;              % Day of year 
SUHDR.hour    = epoch(4);               % Hour of day (24 hour clock)
SUHDR.minute  = epoch(5);               % Minute of hour 
SUHDR.sec     = epoch(6);               % Second of minute 
SUHDR.timbas  = 1;                      % Time basis code: 1 = local; 2 = GMT; 
                                        % 3 = other. 
SUHDR.trwf    = 0;                      % Trace weighting factor, defined 
                                        % as 1/2^N volts for the least
                                        % sigificant bit.
SUHDR.grnors  = 0;                      % Geophone group number of roll 
                                        % switch position one
SUHDR.grnofr  = 0;                      % Geophone group number of trace 
                                        % one within original field record
SUHDR.grnlof  = 0;                      % Geophone group number of last 
                                        % trace within original field record 
SUHDR.gaps    = 0;                      % Size (total number of groups dropped)
SUHDR.otrav   = 0;                      % Overtravel taper code: 1 = down 
                                        % (or behind) 2 = up (or ahead) 
% Local assignments, SU version of SEGY headers 
if ~isempty(DATA.dt),
    SUHDR.d1  = DATA.dt;                % Sample spacing for non-seismic data 
else
    SUHDR.d1  = 0;
end
if ~isempty(DATA.tt2w),
    SUHDR.f1      = DATA.tt2w(1);       % First sample location for 
else                                    % non-seismic data 
    SUHDR.f1  = 0;
end
if ~isempty(DATA.dx),
    SUHDR.d2  = DATA.dx;                % Sample spacing between traces 
else
    SUHDR.d2  = 0;
end
if ~isempty(DATA.x),
    SUHDR.f2  = DATA.x(1);              % First trace location 
else
    SUHDR.f2  = 0;
end
SUHDR.ungpow  = 0;                      % Negative of power used for 
                                        % dynamic range compression 
SUHDR.unscale = 0;                      % Reciprocal of scaling factor 
                                        % to normalize range
if ~isempty(DATA.ntr),
    SUHDR.ntr = DATA.ntr;               % Number of traces 
else
    SUHDR.ntr = 0;
end
SUHDR.mark    = 0;                      % Mark selected traces 
SUHDR.shortpad= 0;                      % Alignment padding 
SUHDR.unass   = zeros(1,14);            % Unassigned

return
% END FUNCTION INITSUHEADER