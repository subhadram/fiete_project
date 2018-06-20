amp = 2.0
N=512;
tau = 30;
z_0 = 150*2*pi/N
freqs = 10;
freq = 2*pi/1000.0*freqs;
w00 = 12.5;
pp = zeros(length(tau),1);
theta = 0.1
a = pi/2;
kk = 1:1:125
w = a
%r = w00(20)*(8*sqrt(2*pi)*a)
for i = 1:length(kk)
    k = kk(i)
T = tau/amp*log(z_0/theta)
lamb = 1-(k/(w00*(8*sqrt(2*pi)*a)));
fcn = @(phi) sin(freq*T + phi)/freq^2 + cos(freq*T + phi)/freq  - sin(phi)/freq^2-cos(phi)/freq + amp*tau/(1-lamb) *(exp(-(1-lamb)*T/tau)-1)*(cos(phi)/2*freq - sin(phi)/2)
pp(i) = fzero(fcn, 1);
ss(i) = T;
end

plot(kk,pp)
hold on;
%plot(amp,ss)
xlabel('Coupling strength')
ylabel('Optimal phase difference')
set(gca,'fontsize',14)