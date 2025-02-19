function [SEGYreelhdr, SEGYdataformat] = initSEGYreelheader(DATA)
%
% INITSEGYREELHEADER : Initialize parameters needed for exporting data in
% the SEG-Y format. These include: 
% 1. The SEG-Y reel header structure. Field values are assigned with
%    information passed through the input data structure DATA (== IPD), and 
%    through the ENVAR collection of global runtime variables
% 2. The "SEGYdataformat" structure of available data sample formats.
%
% Usage: initSEGYreelheader(DATA)
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

% Get the MATGPR runtime variables suite
global ENVAR
DataFormat = ENVAR.SEGY_DataSampleFormat;
segyrev    = 1;

%%% ====== PART 1: Initialize the Reel header structure ===================
%%% Step 1: Prepare dialog box to get job ID and textual header
clb = cell(4,1);                            % prompts
clb(1,1) = cellstr('Enter Job ID Number');
clb(2,1) = cellstr('Enter Line Number');
clb(3,1) = cellstr('Enter Reel Number');
clb(4,1) = cellstr('Enter Comments');
cldef = cell(4,1);                          % default answers
cldef(1,1) = cellstr('1');
cldef(2,1) = cellstr('1');
cldef(3,1) = cellstr('1');
cldef(4,1) = cellstr(' ');
num_lines = [1 80; 1 80; 1 80; 10 80];      % number of lines in inpdlg
answer = inputdlg(clb,'SEG-Y WRITE: Reel ID Parameters', num_lines, ...
    cldef,'on');
if isempty(answer),
    SEGYreelhdr = [];
    SEGYdataformat = [];
    return
end

%%% Step 2: Initialize the 3200 byte text block
blankline = char(ones(1,80)*32);
SEGYreelhdr.headertext = blankline;
for i=1:39,
    SEGYreelhdr.headertext =[SEGYreelhdr.headertext; blankline]; 
end
%%% Lead line
rightnow = ['Created  with MATGPR on ' datestr(now) ...
            '. Format: SEG-Y Revision 1'];
%%% User Comments
usertext = answer{4,1};
%%% Create character array with lead line, user comments, header comments
%%% and processing history and insert to reel header.
htext = strvcat(rightnow, ...
                 usertext, ...
                  DATA.comments, ...
                   'PROCESSING HISTORY FOLLOWS:', ...
                     str2mat(DATA.history));
ht = size(htext);
if ht(1) > 40,
    ht(1) = 40;
end
if ht(2) > 80,
    ht(1) = 80;
end
SEGYreelhdr.headertext(1:ht(1),1:ht(2)) = htext(1:ht(1),1:ht(2));

%%% Step 3: Initialize the binary block
jobid = str2num(answer{1});
SEGYreelhdr.jobid = jobid;         % Job identification number
linenum = str2num(answer{2});
SEGYreelhdr.lino  = linenum;       % Line number (only one line per reel)  
reelnum = str2num(answer{3});
SEGYreelhdr.reno  = reelnum;       % Reel number  
SEGYreelhdr.ntrpr = DATA.ntr;      % Number of traces per record  
SEGYreelhdr.nart  = 0;             % Number of auxiliary traces per record  
SEGYreelhdr.hdt   = DATA.dt*1000;  % Sample interval in micro secs for 
                                   % this reel 
SEGYreelhdr.dto   = 0;             % Same for original field recording  
SEGYreelhdr.hns   = DATA.ns;       % Number of samples per trace for this
                                   % reel 
SEGYreelhdr.nso   = 0;             % Number of samples per trace for 
                                   % original field recording      
SEGYreelhdr.format= DataFormat;    % Data sample format code: 
                                   % 1 = IBM floating point (4 byte) 
                                   % 2 = Fixed point, two's complement (4 byte) 
                                   % 3 = Fixed point, two's complement (2 byte) 
                                   % 4 = Fixed point w/gain code (UNSUPPORTED)
                                   % 5 = IEEE floating point (4 byte)
                                   % 6 = Not currently used
                                   % 7 = Not currently used
                                   % 8 = Fixed-point two's complement (1 byte)
SEGYreelhdr.fold  = 0;             % CDP fold expected per CDP ensemble
SEGYreelhdr.tsort = 1;             % Trace sorting code: 
                                   % 1 = as recorded (no sorting) 
                                   % 2 = CDP ensemble 
                                   % 3 = single fold continuous profile 
                                   % 4 = horizontally stacked 
SEGYreelhdr.vscode= 1;             % Vertical sum code: 
                                   % 1 = no sum 
                                   % 2 = two sum ... N = N sum (N =
                                   % 32,767)
SEGYreelhdr.hsfs   = 0;            % Sweep frequency at start  
SEGYreelhdr.hsfe   = 0;            % Sweep frequency at end 
SEGYreelhdr.hslen  = 0;            % Sweep length (ms)  
SEGYreelhdr.hstyp  = 4;            % Sweep type code: 
                                   % 1 = linear, 
                                   % 2 = parabolic, 
                                   % 3 = exponential, 
                                   % 4 = other    
SEGYreelhdr.schn   = 0;            % Trace number of sweep channel 
SEGYreelhdr.hstas  = 0;            % Trace number of sweep channel 
SEGYreelhdr.hstae  = 0;            % Sweep trace taper length (msec) at
                                   % end (the ending taper starts at sweep
                                   % Length minus the taper length at end)  
SEGYreelhdr.htatyp = 3;            % Sweep trace taper type code: 
                                   % 1 = linear, 
                                   % 2 = cos-squared, 
                                   % 3 = other  
SEGYreelhdr.hcorr  = 1;            % Correlated data traces code: 
                                   % 1 2 3 4 5 6 7 8 95;  1 = no 2 = yes 
SEGYreelhdr.bgrcv  = 0;            % Binary gain recovered code: 
                                   % 1 = yes, 2 = no 
SEGYreelhdr.rcvm   = 4;            % Amplitude recovery method code: 
                                   % 1 = none, 
                                   % 2 = spherical divergence, 
                                   % 3 = AGC, 
                                   % 4 = other 
SEGYreelhdr.mfeet  = 1;            % Measurement system code: 
                                   % 1 = meters 
                                   % 2 = feet 
SEGYreelhdr.polyt  = 0;            % Impulse signal polarity code: 
                                   % 1 = increase in pressure or upward
                                   % geophone case movement gives negative
                                   % number on tape.  
                                   % 2 = increase in pressure or upward
                                   % geophone case movement gives positive
                                   % number on tape  
SEGYreelhdr.vpol   = 0;            % Vibratory polarity code:
SEGYreelhdr.dt     = DATA.dt;      % Trace spacing (MATGPR specific 
                                   % field)
SEGYreelhdr.t1     = DATA.tt2w(1); % First trace location (MATGPR 
                                   % specific field)
SEGYreelhdr.dx     = DATA.dx;      % Trace spacing (MATGPR specific 
                                   % field)
SEGYreelhdr.x1     = DATA.x(1);    % First trace location (MATGPR 
                                   % specific field)
SEGYreelhdr.unass1= zeros(112,1);  % Unassigned, 224 bytes to export as int16
SEGYreelhdr.revision = segyrev;    % SEG-Y revision number
SEGYreelhdr.fixedlen = 1;          % Fixed length trace flag
SEGYreelhdr.ntexthdrs= 0;          % Number of extended textual headers
SEGYreelhdr.unass2= zeros(47,1);   % Unassigned, 94 bytes to export as int16

%%% ====== PART 2: Initialize the data sample format structure ============
SEGYdataformat(1).format='uint32'; 
SEGYdataformat(2).format='int32';  
SEGYdataformat(3).format='int16'; 
SEGYdataformat(4).format='';  
SEGYdataformat(5).format='float32'; 
SEGYdataformat(6).format='';  
SEGYdataformat(7).format='';  
SEGYdataformat(8).format='int8';  

SEGYdataformat(1).name='4-byte IBM Floating Point';
SEGYdataformat(2).name='4-byte two''s complement';
SEGYdataformat(3).name='2-byte two''s complement';
SEGYdataformat(4).name='4-byte fixed-point with gain - NOT SUPPORTED';
SEGYdataformat(5).name='4-byte IEEE floating-point';
SEGYdataformat(6).name='DSF Not currently used';
SEGYdataformat(7).name='DSF Not currently used';
SEGYdataformat(8).name='1-byte, two''s complement integer';
