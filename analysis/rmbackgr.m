function dbg = rmbackgr( d )
%
%   RMBACKGR : Remove a global background trace from the GPR section
%              passed in "d". The background trace is the average trace
%              determined by adding all traces together and dividing by the
%              number of traces. 
%
%      Usage : dbg = rmbackgr( d )
%  
%      Input :   d, the input GPR section
%    
%     Output : dbg, the reduced GPR section
%
%     Author : Andreas Tzanis,
%              Department of Geophysics, 
%              University of Athens
%              atzanis@geol.uoa.gr
%             (C) 2005, Andreas Tzanis, all rights reserved

[m,n] = size( d );
backgr = mean( d' )';
dbg = d - backgr * ones( 1, n );
return
