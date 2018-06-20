ampp = 0.1:0.1:1;
ttt = zeros(length(ampp),1);
N=512;
taut = 30;
z_0 = (N/2 - N/3)*2*pi/N;
freqss = 10;
freqs = 2*pi/1000.0*freqss;
Amp = 100*sqrt(2*a)*pi;
pp = zeros(length(taut),1);
theta = 0.02;
a = 0.5;
k = 20.0;
kc = (Amp*(8*sqrt(2*pi)*a));
U_0 = (1+ sqrt(1 - k/kc))*Amp/(4*sqrt(pi)*a*k);
lamb = 1-(k/(kc*(8*sqrt(2*pi)*a)));
%r = w00(20)*(8*sqrt(2*pi)*a)
for i = 1:length(ampp)
    freq = freqs;
tau = taut;
amp = ampp(i);
%tcn = @t (z_0 - theta)*amp/tau - exp(t) - sqrt()/tau^2*sin(freq*t)/omega - 1/tau^2*sin(freq*t)/omega
%T = tau/amp*erf(z_0)*exp()()
fcn = @(phi) (-T + (amp*log(z_0/theta) -amp*tau/((1-lamb)^2)*(exp(-(1-lamb)*phi(2)/tau) -1) + amp/tau/(1-lamb + tau)*(exp(-theta^2/8*a^2) - exp(-z_0^2/8*a^2) + 2/U_0*(1+amp)*phi(2)/(amp*(1-lamb)) - 2/(amp*U_0*(1-lamb)^2)*tau*(1+amp)*(exp(-phi(2)*(1-lamb)/tau) -1) + 2/U_0*(1+amp)/(amp*(1-lamb))*(-cos(omega*phi(2))/(2*freq^2) + sin(omega*phi(2))/(2*freq) - cos(omega*phi(2) + phi(1))/(2*freq^2) + sin(omega*phi(2) + phi(1))/(2*freq) + 1/(2*freq^2) - amp*cos(phi(1))/(2*freq^2) - 2/(amp*U_0*(1-lamb)^2)*tau*(exp(-phi(2)*(1-lamb)/tau)) -1)*(1/2 + amp*sin(phi(1))/(2*freq) + amp*cos(phi(1))/2)))/tau)
%sin(freq*phi(2) + phi(1))/freq^2 - cos(freq*phi(2) + phi(1))/freq  + sin(phi(1))/freq^2-cos(phi(1))/freq + amp*tau/(1-lamb) *(exp(-(1-lamb)*(phi(2))/tau)-1)*(cos(phi(1))/2*freq - sin(phi(1))/2))/tau;
%fcn = @(phi) cos(freq*T + phi) - cos(phi)
x0 = [1,100];
xxx = fminsearch(fcn,x0);
pp(i) = xxx(1);
pp(i) = mod(pp(i),2*pi)
ttt(i) = xxx(2)
end

plot(ampp,pp)

hold on;
xlabel('coupling strength')
ylabel('Optimal phase difference')
set(gca,'fontsize',14)