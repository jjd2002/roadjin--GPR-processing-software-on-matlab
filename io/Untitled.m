
% tsample=1:200000;
% dt=10^-9;
% t=tsample*dt;
% fa=100*10^6
% fb=300*10^6
% fc=500*10^6
% data=sin(fa*2*pi*t)+4*cos(fb*2*pi*t)+2*cos(fc*2*pi*t);
% data=data';
traceNO=3
data=w9m_500e6;
[NS,NT]=size(data);
d=data(:,traceNO);
d = [d; zeros(size(d))];    % Double length of traces by padding with zeros 
[Samples_per_scan, No_traces] = size(d);
nf = Samples_per_scan/2;    % Resulting number of frequencies 
trvalue = fix( No_traces /2);
tfvalue = nf;
fn=1/(2*dt)      % Nyquist*1000 gives frequency axis in fake MHz  
df=fn/nf;
f=df : df : nf*df;
for traceNO=1:NT

fttr(:,traceNO) = abs(fft(d));         % Get the spectra of all traces
end 

setappdata(gcf,'fttr',fttr)
figure(1)
fttr = fttr(1:nf,:);
spaxsc = 4;
subplot(2,1,1)
pls = plot(f*10^-6,(fttr(1:nf)));
set(gca, 'xlim', [f(1)*10^-6,f(nf)*10^-6]);
y0 = [min(min(fttr)) max(max(fttr))];
set(gca, 'ylim', y0);
xlabel('Frequency (MHz)');
ylb = ylabel('Amplitude');
subplot(2,1,2)
plot(t,data(:,traceNO))