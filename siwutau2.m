tic
ampp = 0.1;
psi = 0;%:0.1:2.0 ;
psi = pi*psi;
trial = 10;
tdd = zeros(length(psi),1);
taut = 10:10:100;
minang  = zeros(length(ampp),1);
mintim  = zeros(length(ampp),1);
freqs = 10;
freq = 2*pi*freqs/1000;
for tr=1:length(taut)

for g = 1:length(psi)
check = 0;
a = 0.5;
am = ampp;
k = 20;
N = 512;
Amp = 100*sqrt(2*a)*pi;
p = N/(2*pi);
x = 1:N;
z = (-N/2:1:N/2-1);
xpos = N/2;
xpos2 = N/3;
zz = round(xpos2);
cs_i =Amp/(sqrt(2*pi)*a)*exp(-(z*2*pi/N).^2/(2*a^2));
cs_i = circshift(cs_i, [0 N/2 - 1]);
J = zeros(N,N);
%psi = 0:0.1:7 ;
%amp = ampp(d);
%lat = zeros(length(psi),trials);
%nospike = zeros(length(psi),trials);


for i = 1:N
    
    J(i,:) = circshift(cs_i,[0 i - 0]);
    
end

figure(1)
surf(J','LineStyle','none'),view([0,90])
set(gca,'fontsize',14), colormap(jet), colorbar
dt = 1;
T = 2000;

tau = taut(tr);
timet = 0:dt:T;

omega = freq;
sw = zeros(length(timet),1);
sw(500:T) = 1;
%r_inti = (0.010*cos(2*pi/N * (x-xpos)));
r_inti = 0.1*exp(-((x-xpos)*2*pi/N).^2/(2*a^2));
r_inti2 = 0.1*exp(-((x-xpos2)*2*pi/N).^2/(0.2*a^2));
r1 = zeros(N,length(timet));
r2 = zeros(N,length(timet));
u1 = zeros(N,length(timet));
u2 = zeros(N,length(timet));
cc = zeros(length(timet),1);
u1(:,1) = r_inti;
u2(:,1) = r_inti2;
%r1(N/2,1) = 0.01;

for t = 2:length(timet)
    tim = timet(t);
    
    u1(:,t) = u1(:,t-1) + (1+cos(omega*t) + sw(t)*am*u2(:,t-1)+ p*(J(:,:)*r1(:,t-1)) - u1(:,t-1))*dt/tau;
    r1(:,t) = (u1(:,t).^2)/(1 + k*p*(sum(u1(:,t).^2)));
    u2(:,t) = u2(:,t-1) + ( 1+cos(omega*t + psi(g))+ p*(J(:,:)*r2(:,t-1)) - u2(:,t-1))*dt/tau;
    r2(:,t) = (u2(:,t).^2)/(1 + k*p*(sum(u2(:,t).^2)));
    
   
  if check == 0     
  [pks,locs] = findpeaks(u1(:,t));
  if locs ==zz
tdd(g) = t;
check =1;
  end
  end
end


surf(u1','LineStyle','none'),view([0,90])
set(gca,'fontsize',14), colormap(jet), colorbar

end

figure(7)
plot(psi,tdd,'o')
hold on;
[mn,ind] = min(tdd);
minang(tr) = psi(ind);
mintim(tr) = mn;
end
toc

figure(8)
plot(ampp,minang)
