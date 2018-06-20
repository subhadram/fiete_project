amp = 50;
N=512;
tau = 30;
z_0 = 150*2*pi/N
freqss = 1:1:100;
freqs = 2*pi/1000.0*freqss;
w00 = 12.5;
pp = zeros(length(tau),1);
theta = 0.2
a = pi;
k = 200.0;
w = a
%r = w00(20)*(8*sqrt(2*pi)*a)
for i = 1:length(freqs)
    freq = freqs(i)
T = tau/amp*log(z_0/theta)
lamb = 1-(0.1/w00*(8*sqrt(2*pi)*a));
fcn = @(phi) sin(freq*T + phi)/freq^2 + cos(freq*T + phi)/freq  - sin(phi)/freq^2-cos(phi)/freq + amp*tau/(1-lamb) *(exp(-(1-lamb)*T/tau)-1)*(cos(phi)/2*freq - sin(phi)/2)
pp(i) = fzero(fcn, 1)
pp(i) = mod(pp(i),2*pi)
end

plot(freqss,pp)
xlabel('Frequency (Hz)')
ylabel('Optimal phase difference')
set(gca,'fontsize',14)