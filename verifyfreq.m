amp = 0.2
N=512;
tau = 30;
z_0 = 150*2*pi/N
freqs = 1:1:100;
freq = 2*pi/1000.0*freqs;
w00 = 12.5;
pp = zeros(length(tau),1);
theta = 0.1
a = pi/2;
k = 0.1:0.1:10
w = a
%r = w00(20)*(8*sqrt(2*pi)*a)
for i = 1:length(k)
T = tau/amp(i)*log(z_0/theta)
lamb = 1-(0.1/(w00*(8*sqrt(2*pi)*a)));
fcn = @(phi) sin(freq*T + phi)/freq^2 + cos(freq*T + phi)/freq  - sin(phi)/freq^2-cos(phi)/freq + amp(i)*tau/(1-lamb) *(exp(-(1-lamb)*T/tau)-1)*(cos(phi)/2*freq - sin(phi)/2)
pp(i) = fzero(fcn, 1);
ss(i) = T;
end

plot(amp,pp)
hold on;
%plot(amp,ss)
xlabel('Coupling strength')
ylabel('Optimal phase difference')
set(gca,'fontsize',14)