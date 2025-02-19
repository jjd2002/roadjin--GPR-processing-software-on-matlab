function OPD = do_mmfilter(IPD,filtertype)
%
% Callback function to drive the routine "mmfilter.m"
% Apply a 1-D or 2-D, Mean or Median spatial filter
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%

OPD = discardprocdata; 
% Trap common errors
if isempty(IPD.d), 
   erh = errordlg('No data to process!', 'MATGPR: ERROR'); 
   uiwait(erh);  
   return
end

%%%%%  Get size of filter window 
clb    = cell(2,1);                                        
clb(1) = cellstr(['沿时间轴的样本数, >=1']);
clb(2) = cellstr(['沿坐标轴的样本数,  >=1']);
answer = inputdlg(clb,'Give Filter Dimensions',1 );          
if isempty(answer),                                   % operation canceled
    return
end
fdim = str2num(char(answer));                          

% Proceed ...
OPD   = IPD;
OPD.d = mmfilter(IPD.d, fdim(1), fdim(2), filtertype);   
if isempty(OPD.d),
    disp('MMFILTER > Operation aborted - No data returned!');
    return
end
viewdata(OPD.x,OPD.tt2w,OPD.d,'outdata',OPD.xlab,OPD.zlab); 
% Update history
iss = size(OPD.history,1); 
fsize = ['[ ' num2str(fdim(1)) ' x ' num2str(fdim(2)) ' ]'];
if strcmp(lower(filtertype),'medifilt'), %#ok<STCI>
    title([fsize ' Median filtered data']);
    text = ['Applied ' fsize ' Median spatial filter'];
elseif strcmp(lower(filtertype),'meanfilt'), %#ok<STCI>
    title([fsize ' Mean filtered data']);
    text = ['Applied ' fsize ' Mean spatial filter'];
end
OPD.history(iss+1,1) = cellstr(text); 
return
