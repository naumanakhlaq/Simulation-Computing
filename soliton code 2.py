import numpy as np              #RELEVANT MODULES
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib import animation   #to animate the plot

diff = False     #set true for diffusion fde (Task 3.3)
frames = 2000     #number of frames to animate

#DEAFULT: ONE SOLITON PROPOGATING
#SEE LINE 114-115 FOR FUNCTION DECLARATION.
#FOR WAVE BREAKING, CHANGE U0 AND MAKE SURE ONLY ONE SOLITON IN def soliton(a, x):

def soliton(a, x): #initilizes the soliton
    result = []
    
    for i in x:
        pt = 12 * (a ** 2) * (1 / (np.cosh(a * i) ** 2)) #+ (12 * ((a + da) ** 2) / (np.cosh((a + da) * (i + 5.0)) ** 2))
        result.append(pt)                                #^^ uncomment for two solitons ^^
	
    result = np.array(result)
    return result
    
    
def f2(u, h): # derivative for space descritisation
    
    if diff:                     #forward function calling
        return f2diff(u, h)     #due to the boolean declared above, different fde used
    else:                       #saves on repeating code for another rk4 method
        
        i = 0         #STANDARD METHOD FOR ALL OTHER TASKS (EX. 3.3)
        ut = []
    
        try:     #for error debugging purposes mainly
            while i < len(u): #sets up the loop
                
                if i < (len(u) - 2):       # boundary conditions
                    
                    uj = ((-(0.25 / h)) * (((u[i+1]) ** 2) - (u[i-1]) ** 2) - 
                        (0.5 / (h ** 3)) * (u[i+2] - (2 * u[i+1]) + (2 * u[i-1]) - u[i-2]))
                    
                    ut.append(uj)
                
                elif i == (len(u) - 2):   #boundary conditions       
                    
                    uj = ((-(0.25 / h)) * (((u[i+1]) ** 2) - (u[i-1]) ** 2) - 
                        (0.5 / (h ** 3)) * (u[0] - (2 * u[i+1]) + (2 * u[i-1]) - u[i-2]))
                    
                    ut.append(uj)
                
                else:
                    
                    uj = ((-(0.25 / h)) * (((u[0]) ** 2) - (u[i-1]) ** 2) - 
                        (0.5 / (h ** 3)) * (u[1] - (2 * u[0]) + (2 * u[i-1]) - u[i-2]))  
                        
                    ut.append(uj)
                
                i += 1
            
            ut = np.array(ut)
            return ut
            
        except: 
            print 'array error, investigate'
        
def f2diff(u, h):  #different Finite difference method for the diffusion equation
    i = 0
    ut = []
    
    while i < len(u):
        
        if i <= (len(u) - 2):  #boundary conditions
        
            uj = ((-(0.25 / h)) * (((u[i+1]) ** 2) - (u[i-1]) ** 2) 
                    + (1 / (h ** 2)) * (u[i-1] - (2 * u[i]) + u[i+1]))
            ut.append(uj)
                
        else:
            
            uj = ((-(0.25 / h)) * (((u[0]) ** 2) - (u[i-1]) ** 2) 
                    + (1 / (h ** 2)) * (u[i-1] - (2 * u[i]) + u[0]))
            ut.append(uj)
        
        i += 1
    
    ut = np.array(ut)
    return ut    

def rk4(u, h, dt): # For the time discretisation of the KdV equation
    #RK4 method
    k1 = dt * f2(u, h)
    k2 = dt * f2((u + (0.5 * k1)), h)
    k3 = dt * f2((u + (0.5 * k2)), h)
    k4 = dt * f2((u + k3), h)
    
    uj1 = u + ((1.0 / 6.0) * (k1 + (2 * k2) + (2 * k3) + k4))  #propogates to the next point
    return uj1

    
fig = plt.figure()                               #setting up the figure
ax = plt.axes(xlim=(-15, 15), ylim=(-4, 50))     # setting the scale
plt.xlabel('x')
plt.ylabel('u')
line, = ax.plot([], [], lw=2)                    #setting up for animation

#constants
a = 1.
da = 0.3

No = 200.
xl = 15
x = np.linspace(-xl, xl, No)   #outputs an array of x values


#u0 = soliton(a,x)                 #INITIAL FUNCTION
u0 = 25*np.exp(-(0.3*x)**2)# for wavebreaking

h = 2*xl/No         #picked to ensure stability

def prop(h):    #propogation function used to store arrays which then pass to animate
    result = [u0]
    for i in range(frames):
        result.append(rk4(result[i], h, 0.001))    #0.001 is the dt used
    
    return result
        
print(x)

uj = prop(h)    #array which contains the array of points to animate
massarr = []
timearr = []

def init2():
    line.set_data([], [])    #clear the line
    return line,
    
def animate(i):
    line.set_data(x, uj[i])  #replot the line
    ax.set_title('Time of Soliton Propogation: '+str(0.001*i)+ 's')
    mass = np.trapz(uj[i], dx=h) + uj[i][0]*h + uj[i][-1]*h          #to account for full area
    
    timearr.append(0.001*i)  
    massarr.append(mass)
    
    mom = np.trapz(uj[i]**2, dx=h)
    print mass, mom     #debugging the mass to check for stability
    return line,
    
anim = animation.FuncAnimation(fig, animate, init_func=init2,
                               frames=frames, interval=1, blit = True)  #blit True for smoother and faster animation, title will not be dynamic though

    
plt.grid()
plt.show()  #show the plot