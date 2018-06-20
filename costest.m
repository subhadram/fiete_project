dt = 0.1;
T = 1000;
time = 0:dt:T
x = zeros(length(time),1);
x(1) = 1.0;
tau = 30.0;
freqs = 10;
freq = 2*pi*freqs/1000;
omega = freq*dt
for t=2:length(time)
    x(t) = x(t-1) + (cos(omega*t))*dt/tau;
    
end
plot(time,x)
hold on;