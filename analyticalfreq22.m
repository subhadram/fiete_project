tic
ampp = 0.01:0.01:0.10;
psi = 0:0.05:2.0 ;
psi = pi*psi;
trial = 10:10:100;
taut = 30;
freqs = 10;
freq = 2*pi*freqs/1000;
tpsi = zeros(length(psi),1)
amptim = zeros(length(ampp),1) 
minangamp = zeros(length(ampp),1) 
for tr=1:length(ampp)

for g = 1:length(psi)

a = 0.5;
amp = ampp(tr);
k = 0.20;
N = 512;
Ampk = sqrt(2*pi*a);
p = N/(2*pi);
x = 1:N;
%z = (-N/2:1:N/2-1);
xpos = 4*N/20;
xpos2 = 3*N/20;
z_0 = (xpos - xpos2)*2*pi/N
pssi = psi(g)
theta = 0.002;
a = 0.5;
k = 0.20;
kc = (Ampk*(8*sqrt(2*pi)*a));
U_0 = (1+ sqrt(1 - k/kc))*Ampk/(4*sqrt(pi)*a*k);
lamb =1 - sqrt(1-k/kc);
dt = 1;
T = 10000;

tau = taut;
timet = 0:dt:T;

omega = freq;

z = zeros(length(timet),1);
a0 = zeros(length(timet),1);
a1 = zeros(length(timet),1);
a2 = zeros(length(timet),1);
check = 0;
for t = 3:length(timet)
    %dz = z(t-1) - z(t-2);
    I2 = amp*U_0*exp(-(z_0 - z(t-1))^2/(8*a^2))*(z_0 - z(t-1))^2/(2*a)*(sqrt(sqrt(2*pi)*a))+ 1/(sqrt(sqrt(2*pi)*a*8)*a^2)*(1+amp + amp*cos(omega*t + pssi) + cos(omega*t))*(4*a^2 - sqrt(pi)*2*a);
    I1 = amp*U_0*exp(-(z_0 - z(t-1))^2/(8*a^2))*(z_0 - z(t-1))/(2*a)*(sqrt(sqrt(2*pi)*a));
   dz = ((2*a/tau)*(I1 + a1(t-1))/(sqrt(sqrt(2*pi)*a)*U_0 + a0(t-1) + sqrt(0.5)*a2(t-1)));
    I0 = amp*U_0*exp(-(z_0 - z(t-1))^2/(8*a^2))*(sqrt(sqrt(2*pi)*a)) + (sqrt(sqrt(8*pi)*a))*(1+amp + amp*cos(omega*t + pssi) + cos(omega*t));
    a0(t) = a0(t-1) + (I0/tau + a1(t-1)*dz- sqrt(2)/(4*tau)*a2(t-1)- (1-lamb)/tau*a0(t-1))*dt;
    a1(t) = a1(t-1) + (I1/tau - (sqrt(sqrt(2*pi)*a)*U_0 + a0(t-1) - sqrt(2)*a2(t-1))/(2*a)*dz)*dt;
    a2(t) = a2(t-1) + (I2/tau - sqrt(2)*a1(t-1)/(2*a)*dz - 1/(2*tau)*a2(t-1))*dt;
    %Rt = 1  + amp/(1-lamb)*exp(-(1-lamb)*t/tau) + amp/(2-lamb)*(exp(-(z(t-1)-z_0)^2/(8*a^2))  - exp(-z_0^2/(8*a^2))) + sqrt(2)/(U_0*(1-lamb))*(1+amp)*(1 - exp(-(1-lamb)*t/tau))+ 1/(sqrt(2)*U_0*tau)*(sin(omega*t)/(omega) + cos(omega*t) + amp*sin(omega*t + pssi)/(omega)  + amp*cos(omega*t+pssi) - (exp(-(1-lamb)*t/tau)*(1 + sin(pssi)/(omega) + cos(pssi)/2))); 
    %z(t) = z(t-1) + ((amp/tau)*(z_0-z(t-1))*exp(-(z(t-1) - z_0)^2/(8*a^2))+ a1(t-1)/())*dt;
    z(t) = z(t-1) + ((2*a/tau)*(I1 + a1(t-1))/(sqrt(sqrt(2*pi)*a)*U_0 + a0(t-1) + sqrt(0.5)*a2(t-1)))*dt;
    
   
if check == 0        
sdd = z_0 - theta;
if abs(z(t)) >= abs(z_0-theta)
    tpsi(g) =  t;
check =1;
end
end

end


figure(1)
plot(timet,z)
hold on;
end
figure(tr)
plot(psi,tpsi)
hold on;
[loc,mn] = min(tpsi);
minangamp(tr) = psi(mn) ;
amptim(tr) = tpsi(1)
end

