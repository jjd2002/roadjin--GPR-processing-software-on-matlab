function SEGYtracehdr = initSEGYtraceheader(DATA)
%
% INITSEGYTRACEHEADER : Initializes the SEG-Y U trace header structure
% "SEGYtracehdr" and assigns some standard parameters. These are passed by
% the input structure DATA (which can be either IPD or OPD).  
%
% Usage  : initSEGYtraceheader(DATA)
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
SEGYtracehdr.tracl  = 0;            % Trace sequence number within line
SEGYtracehdr.tracr  = 0;            % Trace sequence number within reel
SEGYtracehdr.fldr   = 1;            % Field record number
SEGYtracehdr.tracf  = 0;            % Trace number within field record
SEGYtracehdr.ep     = 0;            % Energy source point number
SEGYtracehdr.cdp    = 0;            % CDP ensemble number 
SEGYtracehdr.cdpt   = 0;            % Trace number within CDP ensemble 
SEGYtracehdr.trid   = 7182;         % Trace identification code:
                                    % MATGPR ID code = 7182 (ASCII codes
                                    % [71 82] == 'GR' for GeoRadar) 
SEGYtracehdr.nvs    = 0;            % Number of vertically summed traces 
                                    % (see vscode in reel header structure)
SEGYtracehdr.nhs    = 0;            % number of horizontally summed traces 
                                    % (see vscode in reel header structure) 
SEGYtracehdr.duse   = 1;            % Data use: 1 = production 2 = test 
if ~isempty(DATA.TxRx),
    SEGYtracehdr.offset = DATA.TxRx*1000;  % Distance from source point to 
else                                       % receiver group (negative if 
    SEGYtracehdr.offset = 0;               % opposite to directionin which 
end                                        % the line was shot)
SEGYtracehdr.gelev  = 0;            % Receiver group elevation from sea 
                                    % level (above sea level is positive)
SEGYtracehdr.selev  = 0;            % Source elevation from sea level 
                                    % (above sea level is positive) 
SEGYtracehdr.sdepth = 0;            % Source depth below surface (positive) 
SEGYtracehdr.gdel   = 0;            % Datum elevation at receiver group
SEGYtracehdr.sdel   = 0;            % Datum elevation at source 
SEGYtracehdr.swdep  = 0;            % Water depth at source 
SEGYtracehdr.gwdep  = 0;            % Water depth at receiver group
SEGYtracehdr.scalel = -3;           % Scale factor for previous 7 
                                    % entries with value +/- 10 to the 
                                    % power 0, 1, 2, 3, or 4 (if positive  
                                    % multiply, if negative divide) 
SEGYtracehdr.scalco = -3;           % Scale factor for next 4 entries 
                                    % with value +/- 10 to the power 
                                    % 0, 1, 2, 3, or 4 (if positive,
                                    % multiply, if negative divide)  
SEGYtracehdr.sx     = 0;            % Source X-coordinate
SEGYtracehdr.sy     = 0;            % Source Y-coordinate 
SEGYtracehdr.gx     = 0;            % Group X-coordinate
SEGYtracehdr.gy     = 0;            % Group Y-coordinate
SEGYtracehdr.counit = 1;            % Coordinate units code for previous
                                    % four entries 
                                    % 1  = length (meters or feet) 
                                    % 2  = seconds of arc (the X-longitude
                                    % Y -latitude, positive to the east
                                    % of Greenwich)                                   % or north of the equator.  
SEGYtracehdr.wevel  = 0;            % Weathering velocity
SEGYtracehdr.swevel = 0;            % Subweathering velocity
SEGYtracehdr.sut    = 0;            % Uphole time at source
SEGYtracehdr.gut    = 0;            % Uphole time at receiver group
SEGYtracehdr.sstat  = 0;            % Source static correction
SEGYtracehdr.gstat  = 0;            % Group static correction
SEGYtracehdr.tstat  = 0;            % Total static applied
SEGYtracehdr.laga   = 0;            % Lag time A, time in ms between 
                                    % end of 240-byte trace identification
                                    % header and time break, positive if 
                                    % time break occurs after end of
                                    % header. Time break is defined as the
                                    % initiation pulse which maybe recorded
                                    % on an auxiliary trace or as otherwise 
                                    % specified by the recording system.
SEGYtracehdr.lagb   = 0;            % Lag time B, time in ms between 
                                    % the time break and the initiation
                                    % time of the energy source, may be
                                    % positive or negative. 
SEGYtracehdr.delrt  = DATA.sigpos*1000; % Delay recording time, time in ms 
                                        % between initiation time of energy
                                        % source and time when recording of
                                        % data samples begins (for deep
                                        % water work if recording does not
                                        % start at zero time).  
SEGYtracehdr.muts   = 0;            % Mute time--start
SEGYtracehdr.mute   = 0;            % Mute time--end
if ~isempty(DATA.ns),
    SEGYtracehdr.ns = DATA.ns;      % Number of samples in this trace
else
    SEGYtracehdr.ns = 0;
end
if ~isempty(DATA.dt),
    SEGYtracehdr.dt = DATA.dt*1000; % Sampling interval; in micro-secs
else
    SEGYtracehdr.dt = 0;
end
SEGYtracehdr.gain   = 0;            % Gain type of field instruments 
                                    % code: 1 = fixed, 2 = binary, 
                                    % 3 = floating point, 4 - N  =
                                    % optional use.   
SEGYtracehdr.igc    = 0;            % Instrument gain constant
SEGYtracehdr.igi    = 0;            % Instrument early or initial gain
SEGYtracehdr.corrsu = 0;            % Correlated:  1 = no;  2 = yes 
SEGYtracehdr.sfs    = 0;            % Sweep frequency at start
SEGYtracehdr.sfe    = 0;            % Sweep frequency at end
SEGYtracehdr.slen   = 0;            % Sweep length in ms 
SEGYtracehdr.styp   = 0;            % Sweep type code: 1 = linear; 
                                    % 2 = parabolic; 3 = exponential; 
                                    % 4 = other 
SEGYtracehdr.stas   = 0;            % Sweep trace taper length at start in ms
SEGYtracehdr.stae   = 0;            % Sweep trace taper length at end in ms
SEGYtracehdr.tatyp  = 0;            % Taper type: 1 =linear, 2 =cos^2, 3 =other
SEGYtracehdr.afilf  = 0;            % Alias filter frequency if used
SEGYtracehdr.afils  = 0;            % Alias filter slope 
SEGYtracehdr.nofilf = 0;            % Notch filter frequency if used 
SEGYtracehdr.nofils = 0;            % Notch filter slope 
SEGYtracehdr.lcf    = 0;            % Low cut frequency if used
SEGYtracehdr.hcf    = 0;            % High cut frequncy if used
SEGYtracehdr.lcs    = 0;            % Low cut slope
SEGYtracehdr.hcs    = 0;            % High cut slope
SEGYtracehdr.year   = epoch(1);     % Year data recorded 
SEGYtracehdr.day    = dayofyear;    % Day of year 
SEGYtracehdr.hour   = epoch(4);     % Hour of day (24 hour clock)
SEGYtracehdr.minute = epoch(5);     % Minute of hour 
SEGYtracehdr.sec    = epoch(6);     % Second of minute 
SEGYtracehdr.timbas = 1;            % Time basis code: 1 = local; 
                                    % 2 = GMT; 3 = other. 
SEGYtracehdr.trwf   = 0;            % Trace weighting factor, defined 
                                    % as 1/2^N volts for the least
                                    % sigificant bit.
SEGYtracehdr.grnors = 0;            % Geophone group number of roll 
                                    % switch position one
SEGYtracehdr.grnofr = 0;            % Geophone group number of trace 
                                    % one within original field record
SEGYtracehdr.grnlof = 0;            % Geophone group number of last 
                                    % trace within original field record 
SEGYtracehdr.gaps   = 0;            % Size (total number of groups dropped)
SEGYtracehdr.otrav  = 0;            % Overtravel taper code: 1 = down 
                                    % (or behind) 2 = up (or ahead) 
SEGYtracehdr.cdpx   = 0;            % X coordinate of CDP position of this trace
SEGYtracehdr.cdpy   = 0;            % Y coordinate of CDP position of this trace
SEGYtracehdr.inline3D= 0;           % In-line number for 3-D poststack data
SEGYtracehdr.crossline3D = 0;       % Cross-line number for 3-D poststack data
SEGYtracehdr.shotpt = 0;            % Shot Point
SEGYtracehdr.shotptscalar =0;       % Shot Point Scalar
SEGYtracehdr.tvmu   = -1;           % Trace Value Measurement Unit
SEGYtracehdr.tcm    = 0;            % Transduction Constant Mantissa
SEGYtracehdr.tcp    = 0;            % Transduction Constant Power
SEGYtracehdr.tdu    = 0;            % Transduction Units
SEGYtracehdr.traceid= 0;            % Trace Identifier
SEGYtracehdr.sthdr  = 0;            % Scalar Trace Header
SEGYtracehdr.srctype= 0;            % Source Type
SEGYtracehdr.sedm   = 0;            % Source Energy Direction Mantissa
SEGYtracehdr.sede   = 0;            % Source Energy Direction Exponent
SEGYtracehdr.smm    = 0;            % Source Measurement Mantissa
SEGYtracehdr.sme    = 0;            % Source Measurement Exponent
SEGYtracehdr.smu    = 0;            % Source Measurement Unit     
if ~isempty(DATA.dx),
    SEGYtracehdr.dx = DATA.dx;      % Sample spacing between traces 
else
    SEGYtracehdr.dx = 0;
end
if ~isempty(DATA.ntr),
    SEGYtracehdr.ntr= DATA.ntr;     % Number of traces 
else
    SEGYtracehdr.ntr= 0;
end
SEGYtracehdr.mark   = 0;            % Mark selected traces 

return
% END FUNCTION INITSEGYTRACEHEADER