tic
ampp = 0.1:0.1:1;
psi = 0:0.05:2.0 ;
psi = pi*psi;
trial = 10;
taut = 30;
freqs = 10;
freq = 2*pi*freqs/1000;
tpsi = zeros(length(psi),1)
minangamp = zeros(length(ampp),1) 
for tr=1:length(ampp)

for g = 1:length(psi)

a = 0.5;
amp = ampp(tr);
k = 20;
N = 512;
Ampk = 100*sqrt(2*a)*pi;
p = N/(2*pi);
x = 1:N;
%z = (-N/2:1:N/2-1);
xpos = N/2;
xpos2 = N/3;
z_0 = (xpos - xpos2)*2*pi/N
pssi = psi(g)
theta = 0.02;
a = 0.5;
k = 20.0;
kc = (Ampk*(8*sqrt(2*pi)*a));
U_0 = (1+ sqrt(1 - k/kc))*Amp/(4*sqrt(pi)*a*k);
lamb =1 - sqrt(1-k/kc)
dt = 1;
T = 1500;

tau = 30;
timet = 0:dt:T;

omega = freq;

z = zeros(length(timet),1);

check = 0;
for t = 2:length(timet)
    Rt = 1  + amp/(1-lamb)*exp(-(1-lamb)*t/tau) + 1/(2-lamb)*(exp(-(z(t-1)-z_0)^2/(8*a^2))  - exp(-z_0^2/(8*a^2))) + sqrt(2)/(U_0*(1-lamb))*(1+amp)*(1 - exp(-(1-lamb)*t/tau))+ 1/(sqrt(2)*U_0*tau)*(sin(omega*t)/(omega) + cos(omega*t) + amp*sin(omega*t + pssi)/(omega)  + amp*cos(omega*t+pssi) - (exp(-(1-lamb)*t/tau)*(1 + sin(pssi)/(omega) + cos(pssi)/2))); 
    z(t) = z(t-1) + ((amp/tau)*(z_0-z(t-1))*exp(-(z(t-1) - z_0)^2/(8*a^2))/Rt)*dt;
    
   
if check == 0        
sdd = z_0 - theta;
if round(z(t),4) == round(sdd,4);
    tpsi(g) =  t;
check =1;
end
end

end


figure(1)
plot(timet,z)
hold on;
end
[loc,mn] = min(tpsi);
minangamp(tr) = psi(mn) ;

end

