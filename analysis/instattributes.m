function dinst = instattributes(din, dt, attribute)
%
% INSTATTRIBUTES : Computes instantaneous attributes of GPR sections. These
%                  are the instantaneous amplitude, phase in [-90, 90],
%                  phase in [-180, 180], unwrapped (continuous) phase and
%                  instantaneous frequency. 
%
%         Usage  : dinst = instattributes(din, dt, attribute)
%
%         INPUTS :
%            din : the 2-D GPR section 
%             dt : the sampling rate 
%      attribute : Keyword. Sellects attribute with one of: 
%                 'amplitude' for instantaneous amplitude (envelope)
%                 'atan'  for instantaneous phase in [-90, 90]
%                 'atan2' for instantaneous phase in [-180, 180]
%                 'unwraped' for continuous instantaneous phase 
%                 'ifreq' for instantaneous frequency
%
%        OUTPUTS : 
%          dinst : the resulting instantaneous attribute
%
%         Author : Andreas Tzanis,
%                  Dept. of Geophysics,
%                  University of Athens
%                  atzanis@geol.uoa.gr
%                  (C) 2005, Andreas Tzanis, all rights reserved.
%
%

[ns, ntr] = size(din);
dinst     = zeros(ns,ntr);
ns3       = 3*ns; 
hh = waitbar(0,'Working ...');
for i=1:ntr,
    tr = din(:,i);
% Hilbert transform %%%%%%%%%%%%%%%%%%%%%%%
    ftr = fft(tr,ns3);
    h  = zeros(ns3,1);
    if ns3 >0 && 2*fix(ns3/2)==ns3   % even
        h([1 ns3/2+1]) = 1;
        h(2:ns3/2) = 2;
    elseif ns3 > 0                 % odd 
        h(1) = 1;
        h(2:(ns3+1)/2) = 2;
    end
    ht = ifft(ftr.*h);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    ht = ht(1:ns,1);
% instantaneous amplitude (envelope)
    if strcmp(lower(attribute),'amplitude'),
        dinst(:,i) = abs(ht);
    end
% instantaneous phase in [-90 90]
    if strcmp(lower(attribute),'atan'),
        jj = find(real(ht));
        dinst(jj,i) = atan(imag(ht(jj))./real(ht(jj)))*180/pi;
        kk = setxor(jj,[1:1:ns]');
        dinst(kk,i) = 90;
    end
% instantaneous phase in [-180, 180]
    if strcmp(lower(attribute),'atan2'),
        dinst(:,i) = angle(ht)*180/pi;
    end
% unwraped instantaneous phase
    if strcmp(lower(attribute),'unwraped'),
        dinst(:,i) = unwrap(angle(ht))*180/pi;
    end
% instantaneous frequency
    if strcmp(lower(attribute),'ifreq'),
        dinst(:,i) = [diff(unwrap(angle(ht))/dt); 0];
    end
%    
    waitbar(i/ntr,hh);
end
close(hh); 
return