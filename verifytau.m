ampp = 0.5;
N=512;
taut = 30;
z_0 = (N/2 - N/3)*2*pi/N;
freqss = 2:0.5:10;
freqs = 2*pi/1000.0*freqss;
Amp = 100*sqrt(2*a)*pi;
pp = zeros(length(taut),1);
theta = 0.02;
a = 0.5;
k = 20.0;
kc = (Amp*(8*sqrt(2*pi)*a));
U_0 = (1+ sqrt(1 - k/kc))*Amp/(4*sqrt(pi)*a*k);
%r = w00(20)*(8*sqrt(2*pi)*a)
for i = 1:length(freqs)
    freq = freqs(i)
tau = taut
amp = ampp
%tcn = @t (z_0 - theta)*amp/tau - exp(t) - sqrt()/tau^2*sin(freq*t)/omega - 1/tau^2*sin(freq*t)/omega
T = tau/amp*log(z_0/theta) ;                          
lamb = 1-(k/kc);
fcn = @(phi) sin(freq*T + phi)/freq^2 + cos(freq*T + phi)/freq  - sin(phi)/freq^2-cos(phi)/freq + amp*tau/(1-lamb) *(exp(-(1-lamb)*T/tau)-1)*(cos(phi)/2*freq - sin(phi)/2);
%fcn = @(phi) cos(freq*T + phi) - cos(phi)
pp(i) = fzero(fcn, 1);
pp(i) = mod(pp(i),2*pi)
end

plot(freqss,pp)
hold on;
xlabel('coupling strength')
ylabel('Optimal phase difference')
set(gca,'fontsize',14)