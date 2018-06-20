
clear all;
tic
%ampp = [0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25];
%ampp = 0.1:0.1:10.0;
ampp = [0.5];
psi = 0:0.1:2.0 ;
psi = pi*psi;
%freqs = 1:1:100;
freqs = [10];
freq = 2*pi/1000.0*freqs;
%w00 = [0.2,0.4,0.8,1.6,2.5];
w00 = [2.5];
taut = [30];
%taut = 1:1:100;
trials = 10;

spikedata = zeros(length(freq),length(taut),length(w00),length(ampp),length(psi));
spikesem = zeros(length(freq),length(taut),length(w00),length(ampp),length(psi));
angless = zeros(length(freq),length(taut),length(w00),length(ampp),trials);
for fr =  1:length(freq)
for tt = 1:length(taut)
for ff = 1:length(w00)
for d = 1:length(ampp)
%network parameters
N = 512;
x = 0:511;
b = 0.5; 
w0 = w00(ff);
xpos = N/2;
xpos2 = N;
dt      = 1;                   
tau     = taut(tt);
tau2 = 30;
T       = 1000; 
taus = 20;
a = zeros(T,1);
a(501:T) = 1;

sss = sqrt(trials);


%recurrent connection
cs_i = w0*(-5 + 5*cos(2*pi/N * x));
cs_j = 2.5*(-5 + 5*cos(2*pi/N * x));
cs_i = circshift(cs_i, [0, -1]);
cs_j = circshift(cs_j, [0, -1]);
J = zeros(N,N);
J2 = zeros(N,N);
%psi = 0:0.1:7 ;
amp = ampp(d);
lat = zeros(length(psi),trials);
nospike = zeros(length(psi),trials);


for i = 1:N
    J(i,:) = circshift(cs_i,[0 i + 0]);
    J2(i,:) = circshift(cs_j,[0 i + 0]);
end


for q=1:trials
for j=1:length(psi)
    check = 0;
m_inti = 0*ones(N,1) + (1e-1*cos(2*pi/N * (x-xpos)))';
m = zeros(N,T/dt);
s1 = zeros(N,T/dt);
s2 = zeros(N,T/dt);
sp1 = zeros(N,T/dt);
sp2 = zeros(N,T/dt);
m(:,1) = m_inti;
n_inti = 0*ones(N,1) + (1e-1*cos(2*pi/N * (x-xpos2)))';
n = zeros(N,T/dt);
n(:,1) = n_inti;
%omega = 0.0628;
omega = freq(fr);
b0 = 1.0;
for t = 2 : T/dt
    %I = 1/N * J2 * m(:,t-1) + b + (cos((omega*t)));
    I = 1/N * J2 * m(:,t-1) + b + (cos((omega*t)+psi(j)));
    I2 = 1/N * J * n(:,t-1) + b + (cos(omega*t))+ a(t)*amp*s1(:,t-1);
    Phi = I.*(I>=0);
    Phi2 = I2.*(I2>=0);
   
    %update population activity
    m(:,t) = m(:,t-1) + (-m(:,t-1)+ Phi )*dt/tau2;
    n(:,t) = n(:,t-1) + (-n(:,t-1) +  Phi2 )*dt/tau;
    prob1 = rand(N,1);
    prob2 = rand(N,1);
    index1 = (prob1 <=m(:,t));
    sp1(index1,t) = 1; %spikes
    index2 = (prob2 <=n(:,t));
    sp2(index2,t) = 1; %spikes 
    s1(:,t) = s1(:,t-1) +(- s1(:,t-1) + sp1(:,t-1))*dt/taus;
    s2(:,t) = s2(:,t-1) +(- s2(:,t-1) + sp2(:,t-1))*dt/taus; %not really necessary
    if check == 0
        if t > 500
        if sp2(xpos,t) == 1
            lat(j,q) = t-500;
            nospike(j,q) = sum(sp1(xpos,500:t));
            check = 1;
        end
    end
    end
    
    
        
end
end
[spikemin,index] = min(nospike(:,q));
angless(fr,tt,ff,d,q) = psi(index);
end
spikedata(fr,tt,ff,d,:) = mean(nospike(:,:),2);
spikesem(fr,tt,ff,d,:) = std(nospike(:,:),0,2)/sss;
spikeanglem(fr,tt,ff,d)= mean(angless(fr,tt,ff,d,:));
spikeanglesem(fr,tt,ff,d) = std(angless(fr,tt,ff,d,:))/sss;


end

end
tt
end

end
%only the last 


ang = squeeze(spikeanglem);
angsem = squeeze(spikeanglesem);

figure(1)
errorbar(taut,ang,angsem,'-o')
xlabel('tau (msec)')
ylabel('Optimal phase difference')
set(gca,'fontsize',14)

% for t=1:length(taut)
%     figure(t)
%     title(sprintf('tau = %d msec 10 Hz',taut(t)))
%     for w = 1:length(w00)
%         errorbar(ampp,ang(t,w,:),angsem(t,w,:),'-o')
%         hold off;
%         legend('recurrent weight = 0.2','recurrent weight = 0.4','recurrent weight = 0.8','recurrent weight = 1.6','recurrent weight = 2.5')
%     end
% end
    
% figure('Position',[100,200,1000,400]);
% %suptitle('For coupling strength = 0.1')
% subplot(2,2,1)
% surf(s1','LineStyle','none'),view([0,90]),ylim([1,T/dt])
% xlim([1,512])%,xticks([1,16,32,48,64])
% xlabel('neuron'), ylabel('time (msec)')
% set(gca,'fontsize',14), colormap(jet), colorbar, caxis([0,1])
% 
% subplot(2,2,2)
% surf(s2','LineStyle','none'),view([0,90]),ylim([1,T/dt])
% xlim([1,512])%,xticks([1,16,32,48,64])
% xlabel('neuron'), ylabel('time (msec)')
% set(gca,'fontsize',14), colormap(jet), colorbar, caxis([0,1])
% 
% subplot(2,2,3)
% plot(s2(:,100),'LineWidth',2);
% hold on;
% plot(s2(:,800),'LineWidth',2);
% legend('At time 100 msec','At time 800 msec')
% xlabel('neuron'), ylabel('s')
% set(gca,'fontsize',14)
% 
% subplot(2,2,4)
% plot(psi,mean(lat(:,:),2),'LineWidth',2)
% xlabel('Phase difference')
% ylabel('latency in msec')
% set(gca,'fontsize',14)
% 
% figure()
% %suptitle('For coupling strength = 0.1')
% plot(psi,mean(nospike(:,:),2),'LineWidth',2)
% xlabel('Phase difference')
% ylabel('number of pre-spike for the post-spike')
% set(gca,'fontsize',14)
 toc

