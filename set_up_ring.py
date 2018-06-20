from numpy import *
from pylab import *
import random
import matplotlib.image as mpimg

N= 360


X = arange(0,360,1)
Y = X*pi/360.0
Z = cos(Y)

#plot(Y,Z)
#show()



W = zeros((N,N))
#print W
w0 = 22.5
sig =1.50
for i in range(N):
    for j in range(N):
        W[i,j] = -w0 + (w0+10)*exp(-((Y[i] - Y[j])**2.0) / 2.0 * sig*sig)


plot(Y,W[:,0])
show()
w1 = zeros((N,N))
w2 = zeros((N,N))
for i in range(N):
    w1[:,i] = roll(W[:,180],i+180)
    w2[:,i] = roll(W[:,180],i+180)
imshow(w1,aspect='auto',origin='lower', cmap = cm.rainbow)
#imgplot.set_cmap('nipy_spectral')
xlabel("Neuron Index")
ylabel("Neuron Index")
#imgplot.ylim(0,360)
colorbar
show()

t = 100.0
dt = 1.00
tau = 60.0
time = arange(0,t,dt)

m_inti = (100*exp(cos(2*pi/N * (X-0))))
n_inti = (100*exp(cos(2*pi/N * (X-232))))
s1 = zeros((len(time),N))
s2 = zeros((len(time),N))
s2[0,:] = m_inti
s1[0,:] = n_inti
a1 = 0.0
#s1[0,180] = 0.1
omega = 30.0
plot(m_inti)
plot(m_inti)
show()

plot(s1[0,:])
plot(s2[0,:])
show()

phi = 10
b = ones(N)
#b = 1000.0*b
for k in range(0,len(time)-1):
    tt = k*0.1
    att = 0
    for i in range(0,N):
        xx = 0.0
        xx1 = 0.0
        for j in range(0,N):
            xx1 =  + (w1[i,j]*s1[k,j])
            xx2 =  + (w2[i,j]*s2[k,j])
        yy1 = xx1 + b[i] 
        #yy = xx + b[i] cos(omega*tt +phi)
        yy2 = xx2 + b[i] 
        if (yy1< 0.0):
            yy1 = 0.0
        if (yy2 < 0.0):
            yy2 = 0.0
	if (k >50):
            att = 1    
        s1[k+1,i] = s1[k,i]+(-s1[k,i] + yy1 + att*a1*s2[k,i] + cos(omega*tt) )*dt/tau
        s2[k+1,i] = s2[k,i]+(-s2[k,i] + yy2 + cos(omega*tt +phi))*dt/tau

imshow(s1,aspect='auto',origin='lower', cmap = cm.rainbow,norm=matplotlib.colors.LogNorm())

xlabel("Neuron Index")
ylabel("time msec")

colorbar()
show()

imshow(s2,aspect='auto',origin='lower', cmap = cm.rainbow,norm=matplotlib.colors.LogNorm())

xlabel("Neuron Index")
ylabel("time msec")

colorbar()
show()
