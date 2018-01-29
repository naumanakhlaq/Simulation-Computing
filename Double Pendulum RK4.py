import scipy
import numpy as np
import pylab as plt

h = 0.1  # step size
D = 0.3     # Damping Co-Efficient
m = 0.01    #mass of bob 1
M = 1.0    #mass of bob 2
R = .1#M/m    #ratio of masses
g = 9.81   #abs value for acceleration due to gravity
l = 1.5    #length
G = 0.0#D/(m* np.sqrt(g*l))   #natural unit
  
final_time = 100 #final time  

I = np.matrix([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1, 0],
               [0, 0, 0, 1]])

L = np.matrix([[0, 0, 1, 0],
               [0, 0, 0, 1],
               [-(R + 1), R, -G, 0],
               [(R + 1), -(R + 1), G*(1 - R**-1), -G/R]])

y0 = np.matrix([[0.01], 
               [0], [0], [0]])

pltyth = []
pltp = []
pltyph = []
E = []

t = 0

while t < final_time:
    try :
        k1 = L*h*y0
        k2 = L*h*(y0 + 0.5*k1)
        k3 = L*h*(y0 + 0.5*k2)
        k4 = L*h*(y0 + k3)
        
        y = y0 + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        y0 = y
        E.append(2*y[0,0]**2 + 2*y[2,0]**2 + y[3, 0]**2 + y[1,0]**2)   #k=1
        pltp.append(t)
        pltyth.append(y[0, 0])
        pltyph.append(y[1, 0])
        t += h
        #print t
    except: print 'no point'
    
   
    #print pltp
    #print plty


plt.subplots_adjust(hspace=.35)

plt.figure(1)

plt.suptitle('RK4 Method - Double Pendulum', fontsize=18, weight='bold')
plt.autoscale(True)

plt.subplot(211)
plt.title('Angular displacement')
plt.plot(pltp, pltyth, label=r'$\theta$')
plt.plot(pltp, pltyph, label=r'$\phi$')
plt.xlabel('Time')
plt.ylabel(r'$\theta$')
plt.grid(True)
plt.legend()

plt.subplot(212)
plt.title('Total energy')
plt.plot(pltp, E, 'r-')
plt.xlabel('Time')
plt.ylabel('Energy')
plt.grid(True)

plt.show()