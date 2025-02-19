function OPD = do_SOMETHING(IPD, other_argument)
%
% Callback function to drive the routine "SOMETHING.m"
% TEMPLATE FOR THE FAST CREATION OF MATGPR CALLBACK FUNCTIONS 
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata;                           % discard the current OPD
% Trap common errors ...
% A most common mistake 
if isempty(IPD.d),                               
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end
% Another common mistake - delete if not applicable
if isempty(IPD.dx),                              
    erh = errordlg('Operation possible only for equally-spaced data!', ...
        'MATGPR : ERROR');
    uiwait(erh); 
    return;
end;
% Test CONDITIONS for other_argument - delete if not applicable
if  other_argument,                     
    erh = warndlg('Insert your ERROR or WARNING message!', ...
        'MATGPR : ERROR!'); 
    uiwait(erh); 
    return; 
end; 

% Now Proceed ...
OPD = IPD;                                       % Copy IPD to OPD

% Run SOMETHING.m
OPD = SOMETHING(IPD, other_argument);

% Check if Operation has been canceled or aborted and issue a message
if isempty(OPD.d),
    disp('SOMETHING > Operation aborted - No O/P data returned!');
    return
end

% Display data
viewdata(OPD.x,OPD.z,OPD.d,'outdata',OPD.xlab,OPD.zlab);
title('SOMETHING has been done')

% Update processing history
iss = size(OPD.history,1); 
text = 'SOMETHING has indeed been done';
OPD.history(iss+1,1) = cellstr(text); 

return
