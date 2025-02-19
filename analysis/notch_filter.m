function dnf = notch_filter(scanline, d, dt)
%
% NOTCH_FILTER : Design and apply a notch frequency filter 
%
%        Usage : dnf = notch_filter(scanline, d, dt) 
%   
%          dnf : output, 2-D array of filtered data
%    scanline  : input, vector of trace distances along the scanline (MUST BE EQUALLY SPACED)
%            d : input, 2-D data array 
%           dt : input, Sampling interval 
%
%       Author : Andreas Tzanis,
%                Department of Geophysics, 
%                University of Athens
%    Finalized : 17 November 2003
%
%

[ns, ntr] = size(d);
fn=1/(2*dt*0.001);                          % Nyquist multiplied by 1000 to give frequency axis in fake MHz
df=fn/((ns)/2);
f=df:df:fn;
%%%%%%%%   Get a test trace to design filter 
datafig = findobj('tag','datafigure');
if ~ishandle(datafig),
    erh = errordlg('No data figure for the filter design interface','NOTCH_FILTER : ERROR');
    uiwait(erh);
    return
end
%%%% display the help message
message     = msgbox('Please click on a test trace','NOTCH_FILTER : HELP');
messagepos  = get(message,'position');
set(message,'pos',[20 60 150 messagepos(4)]);
pause(0.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(datafig)
pause(0.5)
[xtr,ytr]=ginput(1);
trnum = round((xtr-scanline(1))/(scanline(2)-scanline(1)));
tr=d(:,trnum);
tr=tr-mean(tr);
fttr=fft(tr); 
%%%%%%%%  Design the filter
reply = 'No';
while strcmp(reply,'Yes')~=1,
  if exist('trspfig') & ishandle(trspfig),
      figure(trspfig)
  else
     trspfig = figure('name','Select Frequecy to Remove');
  end
  plot(f,abs(fttr(1:ns/2))); set(gca,'nextplot','add');
  xlabel('Frequency (MHz)');   ylabel(['Amplitude spectrum of trace ' num2str(trnum) ]); 
  axis tight
  [fcutoff,dummy]=ginput(1);
% This design algorithm assumes that Nyquist is 0.5, not 1 - fool it accordingly!
  wo = 2*pi*fcutoff/(2*fn);
% Set a zero on the unit circle at the angular frequency wo. 
  rez = cos(wo);          imz = sin(wo);
% Add a pole close to the zero, but inside the unit circle.  
  rez1 = .99*cos(wo);     imz1 = .99*sin(wo);
% Define numerator and denominator polynomials of the b and a
  b = conv([1 -rez-sqrt(-1)*imz],[1 -rez+sqrt(-1)*imz]);
  a = conv([1 -rez1-sqrt(-1)*imz1],[1 -rez1+sqrt(-1)*imz1]);
% Test the filter
  b = [b'; zeros(ns-length(b),1)];
  B = fft(b);
  a = [a'; zeros(ns-length(a),1)];
  A = fft(a);
  H = B./A; 
  H = H/max(abs(H));                 % normalize the filter
  fttrf  = fttr.*H;
  fttrf  = fttrf.*conj(H);
  figure(trspfig); 
  cla
  plot(f,abs(fttr(1:ns/2))/max(abs(fttr)),f,abs(fttrf(1:ns/2))/max(abs(fttr)),f,abs(H(1:ns/2)),'r');
  xlabel('Frequency (MHz)'); ylabel(['Amplitude spectrum of trace ' num2str(trnum) ]); 
  title(['Notch frequency = ' num2str(fcutoff,5)]);
  axis tight
  drawnow
  reply = questdlg('Is this Satisfactory ? ');
  if strcmp(reply,'Cancel')==1,
      delete(trspfig); clear trspfig
      dnf=[];    
      return
  end
end
delete(trspfig); 
if exist('message') & ishandle(message),
    delete(message)
end

%%%% Filter in the frequency domain twice to preserve phase
hw = helpdlg('Processing... Please wait ...');
drawnow

H     = H * ones(1, ntr);         % filter already in frequency domain - create ns x ntr matrix operator
fttr  = fft(d);                   % transform data to frequecy domain
fttr  = fttr.*H;                  % filter forwards
fttr  = fttr.*conj(H);            % filter backwards
dnf   = real(ifft(fttr));         % invert to time domain

close(hw);
return

