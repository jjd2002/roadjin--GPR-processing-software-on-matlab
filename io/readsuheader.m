function SUHDR = readsuheader(fid)
%
% READSUHEADER : Read SU trace headers and returns their content in the
%                structure SUHDR 
%
%       Author : Andreas Tzanis, 
%                Department of Geophysics, 
%                University of Athens
%                (C) 2005, Andreas Tzanis, allrights reserved
%
%

SUHDR.tracl   = fread(fid,1,'int32');   % Trace sequence number within line
SUHDR.tracr   = fread(fid,1,'int32');   % Trace sequence number within reel
SUHDR.fldr    = fread(fid,1,'int32');   % Field record number
SUHDR.tracf   = fread(fid,1,'int32');   % Trace number within field record
SUHDR.ep      = fread(fid,1,'int32');   % Energy source point number
SUHDR.cdp     = fread(fid,1,'int32');   % CDP ensemble number 
SUHDR.cdpt    = fread(fid,1,'int32');   % Trace number within CDP ensemble 
SUHDR.trid    = fread(fid,1,'int16');   % Trace identification code:
                                        % MATGPR ID code = 7182 (ASCII codes
                                        % [71 82] == 'GR' for GeoRadar) 
SUHDR.nvs     = fread(fid,1,'int16');   % Number of vertically summed traces 
                                        % (see vscode in reel header structure)
SUHDR.nhs     = fread(fid,1,'int16');   % Number of horizontally summed traces 
                                        % (see vscode in reel header structure) 
SUHDR.duse    = fread(fid,1,'int16');   % Data use: 1  = production 2  = test 
SUHDR.offset  = fread(fid,1,'int32');   % Distance from source point to receiver group 
                                        % (negative if opposite to direction in which the line was shot) 
SUHDR.gelev   = fread(fid,1,'int32');   % Receiver group elevation from sea level (above sea level is positive)
SUHDR.selev   = fread(fid,1,'int32');   % Source elevation from sea level (above sea level is positive) 
SUHDR.sdepth  = fread(fid,1,'int32');   % Source depth below surface (positive) 
SUHDR.gdel    = fread(fid,1,'int32');   % Datum elevation at receiver group
SUHDR.sdel    = fread(fid,1,'int32');   % Datum elevation at source 
SUHDR.swdep   = fread(fid,1,'int32');   % Water depth at source 
SUHDR.gwdep   = fread(fid,1,'int32');   % Water depth at receiver group
SUHDR.scalel  = fread(fid,1,'int16');   % Scale factor for previous 7 entries with value +/- 10 to the power 
                                        % 0, 1, 2, 3, or 4 (if positive, multiply, if negative divide) 
SUHDR.scalco  = fread(fid,1,'int16');   % Scale factor for next 4 entries with value +/- 10 to the power 
                                        % 0, 1, 2, 3, or 4 (if positive, multiply, if negative divide) 
SUHDR.sx      = fread(fid,1,'int32');   % X source coordinate
SUHDR.sy      = fread(fid,1,'int32');   % Y source coordinate 
SUHDR.gx      = fread(fid,1,'int32');   % X group coordinate
SUHDR.gy      = fread(fid,1,'int32');   % Y group coordinate
SUHDR.counit  = fread(fid,1,'int16');   % Coordinate units code: for previous four entries 
                                        % 1  = length (meters or feet) 
                                        % 2  = seconds of arc (the X -longitude Y -latitude, positive to
                                        % the east of Greenwich or north of the equator 
SUHDR.wevel   = fread(fid,1,'int16');   % Weathering velocity
SUHDR.swevel  = fread(fid,1,'int16');   % Subweathering velocity
SUHDR.sut     = fread(fid,1,'int16');   % Uphole time at source
SUHDR.gut     = fread(fid,1,'int16');   % Uphole time at receiver group
SUHDR.sstat   = fread(fid,1,'int16');   % Source static correction
SUHDR.gstat   = fread(fid,1,'int16');   % Group static correction
SUHDR.tstat   = fread(fid,1,'int16');   % Total static applied
SUHDR.laga    = fread(fid,1,'int16');   % Lag time A, time in ms between end of 240-byte trace identification 
                                        % header and time break, positive if time break occurs after end of header, 
                                        % time break is defined as the initiation pulse which maybe recorded 
                                        % on an auxiliary trace or as otherwise specified by the recording system 
SUHDR.lagb    = fread(fid,1,'int16');   % Lag time B, time in ms between the time break and the initiation time 
                                        % of the energy source, may be positive or negative 
SUHDR.delrt   = fread(fid,1,'int16');   % Delay recording time, time in ms between initiation time of energy source 
                                        % and time when recording of data samples begins (for deep water work if 
                                        % recording does not start at zero time) 
SUHDR.muts    = fread(fid,1,'int16');   % Mute time--start
SUHDR.mute    = fread(fid,1,'int16');   % Mute time--end
SUHDR.ns      = fread(fid,1,'uint16');  % Number of samples in this trace
SUHDR.dtsu    = fread(fid,1,'uint16');  % sample interval; in micro-seconds
SUHDR.gain    = fread(fid,1,'int16');   % Gain type of field instruments code:
                                        % 1  = fixed 2  = binary 3  = floating point 4 - N  = optional use 
SUHDR.igc     = fread(fid,1,'int16');   % Instrument gain constant
SUHDR.igi     = fread(fid,1,'int16');   % Instrument early or initial gain
SUHDR.corrsu  = fread(fid,1,'int16');   % Correlated:  1  = no;  2  = yes 
SUHDR.sfs     = fread(fid,1,'int16');   % Sweep frequency at start
SUHDR.sfe     = fread(fid,1,'int16');   % Sweep frequency at end
SUHDR.slen    = fread(fid,1,'int16');   % Sweep length in ms 
SUHDR.styp    = fread(fid,1,'int16');   % Sweep type code: 1  = linear; 2  = parabolic; 3  = exponential; 4  = other 
SUHDR.stas    = fread(fid,1,'int16');   % Sweep trace taper length at start in ms
SUHDR.stae    = fread(fid,1,'int16');   % Sweep trace taper length at end in ms
SUHDR.tatyp   = fread(fid,1,'int16');   % Taper type: 1 =linear, 2 =cos^2, 3 =other
SUHDR.afilf   = fread(fid,1,'int16');   % Alias filter frequency if used
SUHDR.afils   = fread(fid,1,'int16');   % Alias filter slope 
SUHDR.nofilf  = fread(fid,1,'int16');   % Notch filter frequency if used 
SUHDR.nofils  = fread(fid,1,'int16');   % Notch filter slope 
SUHDR.lcf     = fread(fid,1,'int16');   % Low cut frequency if used
SUHDR.hcf     = fread(fid,1,'int16');   % High cut frequncy if used
SUHDR.lcs     = fread(fid,1,'int16');   % Low cut slope
SUHDR.hcs     = fread(fid,1,'int16');   % High cut slope
SUHDR.year    = fread(fid,1,'int16');   % Year data recorded 
SUHDR.day     = fread(fid,1,'int16');   % Day of year 
SUHDR.hour    = fread(fid,1,'int16');   % Hour of day (24 hour clock)
SUHDR.minute  = fread(fid,1,'int16');   % Minute of hour 
SUHDR.sec     = fread(fid,1,'int16');   % Second of minute 
SUHDR.timbas  = fread(fid,1,'int16');   % Time basis code: 1  = local; 2  = GMT; 3  = other 
SUHDR.trwf    = fread(fid,1,'int16');   % Trace weighting factor
SUHDR.grnors  = fread(fid,1,'int16');   % Geophone group number of roll switch position one
SUHDR.grnofr  = fread(fid,1,'int16');   % Geophone group number of trace one within original field record
SUHDR.grnlof  = fread(fid,1,'int16');   % Geophone group number of last trace within original field record 
SUHDR.gaps    = fread(fid,1,'int16');   % Size (total number of groups dropped)
SUHDR.otrav   = fread(fid,1,'int16');   % Overtravel taper code: 1  = down (or behind) 2  = up (or ahead) 
% Local assignments, SU version of SEG-Y headers 
SUHDR.d1      = fread(fid,1,'float32'); % Sample spacing for non-seismic data 
SUHDR.f1      = fread(fid,1,'float32'); % First sample location for non-seismic data 
SUHDR.d2      = fread(fid,1,'float32'); % Sample spacing between traces 
SUHDR.f2      = fread(fid,1,'float32'); % First trace location 
SUHDR.ungpow  = fread(fid,1,'float32'); % Negative of power used for dynamic range compression 
SUHDR.unscale = fread(fid,1,'float32'); % Reciprocal of scaling factor to normalize range
SUHDR.ntr     = fread(fid,1,'int32');   % Number of traces 
SUHDR.mark    = fread(fid,1,'int16');   % Mark selected traces 
SUHDR.shortpad  = fread(fid,1,'int16'); % Alignment padding 
SUHDR.unass   = fread(fid,14,'int16');  % Unassigned
%SUHDR.unass   = fread(fid,28,'char');  % Unassigned

% END FUNCTION READSUHEADER