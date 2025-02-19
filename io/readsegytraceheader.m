function SEGYtracehdr = readsegytraceheader(fid)
%
% READSEGYTRACEHDR: Read SEG-Y trace headers and return their content in
% the structure "SEGYtracehdr"
%
% Usage : SEGYtracehdr = readSEGYtracehdr(fid)
%
%Author : Andreas Tzanis, 
%         Department of Geophysics, 
%         University of Athens
%         atzanis@geol.uoa.gr
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

SEGYtracehdr.tracl = fread(fid,1,'int32');   % Trace sequence number within line
SEGYtracehdr.tracr = fread(fid,1,'int32');   % Trace sequence number within reel
SEGYtracehdr.fldr  = fread(fid,1,'int32');   % Field record number
SEGYtracehdr.tracf = fread(fid,1,'int32');   % Trace number within field record
SEGYtracehdr.ep    = fread(fid,1,'int32');   % Energy source point number
SEGYtracehdr.cdp   = fread(fid,1,'int32');   % CDP ensemble number 
SEGYtracehdr.cdpt  = fread(fid,1,'int32');   % Trace number within CDP ensemble 
SEGYtracehdr.trid  = fread(fid,1,'int16');   % Trace identification code:
                                             % MATGPR ID code = 7182 (ASCII
                                             % codes [71 82] == 'GR' for
                                             % GeoRadar) 
SEGYtracehdr.nvs   = fread(fid,1,'int16');   % Number of vertically summed  
                                             % traces 
SEGYtracehdr.nhs   = fread(fid,1,'int16');   % Number of horizontally summed 
                                             % traces 
SEGYtracehdr.duse  = fread(fid,1,'int16');   % Data use: 1 = production 2 = test 
SEGYtracehdr.offset= fread(fid,1,'int32');   % Distance from source point 
                                             % to receiver group (negative
                                             % if opposite to direction in
                                             % which the line was shot)   
SEGYtracehdr.gelev = fread(fid,1,'int32');   % Receiver group elevation 
                                             % from sea level (above sea
                                             % level is positive) 
SEGYtracehdr.selev = fread(fid,1,'int32');   % Source elevation from sea 
                                             % level (above sea level is
                                             % positive)  
SEGYtracehdr.sdepth= fread(fid,1,'int32');   % Source depth below surface 
                                             % (positive) 
SEGYtracehdr.gdel  = fread(fid,1,'int32');   % Datum elevation at receiver group
SEGYtracehdr.sdel  = fread(fid,1,'int32');   % Datum elevation at source 
SEGYtracehdr.swdep = fread(fid,1,'int32');   % Water depth at source 
SEGYtracehdr.gwdep = fread(fid,1,'int32');   % Water depth at receiver group
SEGYtracehdr.scalel= fread(fid,1,'int16');   % Scale factor for previous 7 
                                             % entries with value +/- 10 to
                                             % the power 0, 1, 2, 3, or 4
                                             % (if positive, multiply, if
                                             % negative divide)   
SEGYtracehdr.scalco= fread(fid,1,'int16');   % scale factor for next 4 
                                             % entries with value +/- 10 to
                                             % the power   0, 1, 2, 3, or 4
                                             % (if positive, multiply, if
                                             % negative divide)   
SEGYtracehdr.sx    = fread(fid,1,'int32');   % Source X-coordinate
SEGYtracehdr.sy    = fread(fid,1,'int32');   % Source Y-coordinate 
SEGYtracehdr.gx    = fread(fid,1,'int32');   % Group X-coordinate
SEGYtracehdr.gy    = fread(fid,1,'int32');   % Group Y-coordinate
SEGYtracehdr.counit= fread(fid,1,'int16');   % Coordinate units code for 
                                             % previous four entries 
                                             % 1  = length (meters or feet) 
                                             % 2  = seconds of arc (the X
                                             % -longitude Y -latitude,
                                             % positive to the east of
                                             % Greenwich or north of the
                                             % equator   
SEGYtracehdr.wevel = fread(fid,1,'int16');   % Weathering velocity
SEGYtracehdr.swevel= fread(fid,1,'int16');   % Subweathering velocity
SEGYtracehdr.sut   = fread(fid,1,'int16');   % Uphole time at source
SEGYtracehdr.gut   = fread(fid,1,'int16');   % Uphole time at receiver group
SEGYtracehdr.sstat = fread(fid,1,'int16');   % Source static correction
SEGYtracehdr.gstat = fread(fid,1,'int16');   % Group static correction
SEGYtracehdr.tstat = fread(fid,1,'int16');   % Total static applied
SEGYtracehdr.laga  = fread(fid,1,'int16');   % Lag time A, time in ms 
                                             % between end of 240-byte
                                             % trace identification header
                                             % and time break, positive if
                                             % time break occurs after end
                                             % of header.
SEGYtracehdr.lagb  = fread(fid,1,'int16');   % Lag time B, time in ms 
                                             % between the time break and
                                             % the initiation time  
                                             % of the energy source, may be
                                             % positive or negative  
SEGYtracehdr.delrt = fread(fid,1,'int16');   % Delay recording time. Time 
                                             % in ms between initiation
                                             % time of energy source and
                                             % time when recording of data
                                             % samples begins 
SEGYtracehdr.muts  = fread(fid,1,'int16');   % Mute time--start
SEGYtracehdr.mute  = fread(fid,1,'int16');   % Mute time--end
SEGYtracehdr.ns    = fread(fid,1,'uint16');  % Number of samples in this trace
SEGYtracehdr.dt    = fread(fid,1,'uint16');  % Sample interval; in micro-seconds
SEGYtracehdr.gain  = fread(fid,1,'int16');   % Gain type of field instruments code: 
                                             % 1 = fixed;   2 = binary; 
                                             % 3 = floating point; 
                                             % 4 - N  = optional use 
SEGYtracehdr.igc   = fread(fid,1,'int16');   % Instrument gain constant
SEGYtracehdr.igi   = fread(fid,1,'int16');   % Instrument early or initial gain
SEGYtracehdr.corrsu= fread(fid,1,'int16');   % Correlated: 1 = no;  2 = yes 
SEGYtracehdr.sfs   = fread(fid,1,'int16');   % Sweep frequency at start
SEGYtracehdr.sfe   = fread(fid,1,'int16');   % Sweep frequency at end
SEGYtracehdr.slen  = fread(fid,1,'int16');   % Sweep length in ms 
SEGYtracehdr.styp  = fread(fid,1,'int16');   % Sweep type code: 1  = linear; 
                                             % 2 = parabolic; 3 = exponential;
                                             % 4 = other 
SEGYtracehdr.stas  = fread(fid,1,'int16');   % Sweep trace taper length at 
                                             % start in ms
SEGYtracehdr.stae  = fread(fid,1,'int16');   % Sweep trace taper length at 
                                             % end in ms
SEGYtracehdr.tatyp = fread(fid,1,'int16');   % Taper type: 1 =linear, 
                                             % 2 =cos^2, 3 =other
SEGYtracehdr.afilf = fread(fid,1,'int16');   % Alias filter frequency if used
SEGYtracehdr.afils = fread(fid,1,'int16');   % Alias filter slope 
SEGYtracehdr.nofilf= fread(fid,1,'int16');   % Notch filter frequency if used 
SEGYtracehdr.nofils= fread(fid,1,'int16');   % Notch filter slope 
SEGYtracehdr.lcf   = fread(fid,1,'int16');   % Low cut frequency if used
SEGYtracehdr.hcf   = fread(fid,1,'int16');   % High cut frequncy if used
SEGYtracehdr.lcs   = fread(fid,1,'int16');   % Low cut slope
SEGYtracehdr.hcs   = fread(fid,1,'int16');   % High cut slope
SEGYtracehdr.year  = fread(fid,1,'int16');   % Year data recorded 
SEGYtracehdr.day   = fread(fid,1,'int16');   % Day of year 
SEGYtracehdr.hour  = fread(fid,1,'int16');   % Hour of day (24 hour clock)
SEGYtracehdr.minute= fread(fid,1,'int16');   % Minute of hour 
SEGYtracehdr.sec   = fread(fid,1,'int16');   % Second of minute 
SEGYtracehdr.timbas= fread(fid,1,'int16');   % Time basis code: 
                                             % 1 = local; 2 = GMT; 3 = other 
SEGYtracehdr.trwf  = fread(fid,1,'int16');   % Trace weighting factor
SEGYtracehdr.grnors= fread(fid,1,'int16');   % Geophone group number of 
                                             % roll switch position one
SEGYtracehdr.grnofr= fread(fid,1,'int16');   % Geophone group number of 
                                             % trace one within original
                                             % field record 
SEGYtracehdr.grnlof= fread(fid,1,'int16');   % Geophone group number of 
                                             % last trace within original
                                             % field record  
SEGYtracehdr.gaps  = fread(fid,1,'int16');   % Gap Size (total number of 
                                             % groups dropped)
SEGYtracehdr.otrav = fread(fid,1,'int16');   % Overtravel taper code: 
                                             % 1 = down (or behind); 2 = up
                                             % (or ahead)  
SEGYtracehdr.cdpx  = fread(fid,1,'int32');   % X coordinate of CDP position 
SEGYtracehdr.cdpy  = fread(fid,1,'int32');   % Y coordinate of CDP position 
SEGYtracehdr.inline3D =fread(fid,1,'int32'); % In-line number for 3-D 
                                             % poststack data
SEGYtracehdr.crossline3D =fread(fid,1,'int32'); % Cross-line number for 3-D 
                                                % poststack data
SEGYtracehdr.shotpt= fread(fid,1,'int32');      % Shot Point
SEGYtracehdr.shotptscalar=fread(fid,1,'int16'); % Shot Point Scalar
SEGYtracehdr.tvmu  = fread(fid,1,'int16');   % Trace Value Measurement Unit
SEGYtracehdr.tcm   = fread(fid,1,'int32');   % Transduction Constant Mantissa
SEGYtracehdr.tcp   = fread(fid,1,'int16');   % Transduction Constant Power
SEGYtracehdr.tdu   = fread(fid,1,'int16');   % Transduction Units
SEGYtracehdr.traceid  = fread(fid,1,'int16'); % Trace Identifier
SEGYtracehdr.sthdr = fread(fid,1,'int16');    % Scalar Trace Header
SEGYtracehdr.srctype  = fread(fid,1,'int16'); % Source Type
SEGYtracehdr.sedm  = fread(fid,1,'int32');   % Source Energy Direction Mantissa
SEGYtracehdr.sede  = fread(fid,1,'int16');   % Source Energy Direction Exponent
SEGYtracehdr.smm   = fread(fid,1,'int32');   % Source Measurement Mantissa
SEGYtracehdr.sme   = fread(fid,1,'int16');   % Source Measurement Exponent
SEGYtracehdr.smu   = fread(fid,1,'int16');   % Source Measurement Unit 
SEGYtracehdr.dx    = fread(fid,1,'float32'); % Sample spacing between traces 
SEGYtracehdr.ntr   = fread(fid,1,'int16');   % Number of traces 
SEGYtracehdr.mark  = fread(fid,1,'int16');   % Mark selected traces 

% END FUNCTION READSEGYTRACEHEADER