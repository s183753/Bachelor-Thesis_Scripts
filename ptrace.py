#---------------------------------------------------------------------------#
# Victor Lorentzen, s183753                                                 #
# BSc. Physics and Nanotechnology                                           #
# Bachelor's Thesis: Ray Tracing of O-, X- and EBW-waves in Cold Plasma     #
#---------------------------------------------------------------------------#

from subprocess import Popen, PIPE
import numpy as np
import math 
import matplotlib.pyplot as plt

# Open C++ file to use for calculations using file named "trace"
p = Popen(['./trace'], shell=True, stdout=PIPE, stdin=PIPE)

# Physical constants/variables
omega = 2*math.pi*140*10**(9)           # Input frequency
m_e = 9.1093837015e-31                  # Mass of electrons
q = 1.602176634e-19                     # Charge of electrons
epsilon_0 = 8.8541878128e-12            # Vacuum permitivity
c = 299792458;                          # Speed of light in vacuum

# Define variables specific to system

# FOR EBW
N = 400000                               # Number of Iterations
delta_t = 0.00001                         # Time step size

#N = 30000
#delta_t = 0.0001

# FOR X-mode
#N = 300
#delta_t = 0.01

#N = 10000
#delta_t = 0.0002
# O-mode
N = 1000
delta_t = 0.002

#alphas = np.arange(math.pi, 3/2*math.pi,step = 1/32*math.pi)

alphas = np.linspace(3/4*math.pi,5/4*math.pi,2)
#alphas = alphas + math.pi/2
#alphas = alphas + math.pi
#alphas = np.append(alphas,alphas)
alphas = alphas.tolist()

#alphas = [math.pi,math.pi + 1/4*math.pi, math.pi + 1/8*math.pi]      # First angle
#alphas = [math.pi + 1/4*math.pi,math.pi + 1/3*math.pi, math.pi]  
#alphas = [math.pi,math.pi]
#alphas = [math.pi/2]
alphas = [math.pi]
#alphas = [0,math.pi/4,math.pi/2,math.pi*3/4,math.pi,math.pi*5/4,math.pi*6/4,math.pi*7/4]
#alphas = [0]    
#alphas = [math.pi/2,math.pi*5/4]

runs = len(alphas)                                # Number of runs
betas = [(math.pi)/2]*runs                # Second Angle

L_x_min = 2                             # Lower Limit of radius
L_x = 5                                 # Upper Limit of radius
L_y = 1.5                               # Upper Limit of Z axis
L_y_min = -L_y                          # Lower Limit of Z axis

#startCoord = 2.01                      # Start from left
#startCoord = 4.4                        # Start from center
startCoord = 4.999                     # Start from right

# MODES
# 1 = O-mode, -1 = X-mode, 0 = EBW
mode = [-1]*runs
mode[0] = 1 
#onos = [1]*int(runs/2)
#mode[0:1] = onos



r_init = []
for i in range(runs):
    r_init.append([startCoord - 0.2*i,0.1*i,0])    # Initial conditions for r-vector
dim = len(r_init[0])                    # Dimensions. dim-Dimensional space


# notes index of end of ray
OBIndex = np.full(runs,N-1)

parameters = [N, delta_t, L_x_min, L_x, L_y_min, L_y, runs, dim, betas[0]]

# Input parameters to C++ program
for i in range(len(parameters)):
    current = bytes(str(parameters[i]) + "\n", 'UTF-8')     # Remember + "\n"
    p.stdin.write(current)                                  # Write to std::cin
    p.stdin.flush()                                     

for i in range(runs):
    current = bytes(str(mode[i]) + "\n", 'UTF-8')
    print(mode[i])
    p.stdin.write(current)
    p.stdin.flush()                                     

nGrid = np.zeros([100,100])
BGrid = np.zeros([100,100])
TGrid = np.zeros([100,100])
for i in range(100):
    for j in range(100):
        nGrid[i][j] = float(p.stdout.readline().strip())
        BGrid[i][j] = float(p.stdout.readline().strip())
        TGrid[i][j] = float(p.stdout.readline().strip())
        
for i in range(runs):
    for j in range(dim):
        current = bytes(str(r_init[i][j]) + "\n", 'UTF-8')
        p.stdin.write(current)
        p.stdin.flush()
   
k0 = [0]*runs
for i in range(runs):
    k0[i] = float(p.stdout.readline().strip())
    print(k0[i])

    
#k0 = math.sqrt((omega**2 - (2.43e19*q**2)/(m_e*epsilon_0))/c**2)      # k0
k_init = [0]*runs
for i in range(runs):
    if dim == 3:
        #k_init[i] = [k0[i]*math.cos(alphas[i]),k0[i]*math.sin(alphas[i]),0]
        k_init[i] = [k0[i]*np.cos(alphas[i])*np.sin(betas[i]),k0[i]*np.sin(alphas[i])*np.sin(betas[i]),np.cos(betas[i])]
    
    else:
        k_init[i] = [k0[i]*math.cos(alphas[i]),k0[i]*math.sin(alphas[i])]       # Initial conditions for k-vector. Maybe change to input angle
    
for i in range(runs):
    for j in range(dim):
        current = bytes(str(k_init[i][j]) + "\n", 'UTF-8')
        p.stdin.write(current)
        p.stdin.flush()


n_e = np.zeros((runs,N))
theta = np.zeros((runs,N))
modes = np.zeros((runs,N))
ms = np.zeros((runs,N))

for i in range(N):
    for j in range(runs):
        #print("Got here: " + str(i + 1))
        ms[j][i] = float(p.stdout.readline().strip())
        modes[j][i] = float(p.stdout.readline().strip())
        theta[j][i] = float(p.stdout.readline().strip())
        if (i > 10000 - 1 and i % 10000 == 0):
            print("Reached N = " + str(i))
        #n_e [j][i] = float(p.stdout.readline().strip())


ECRsteps = 1000
RECR = np.zeros(ECRsteps)
RUH = np.zeros(ECRsteps)
RLH = np.zeros(ECRsteps)
for i in range(ECRsteps):
    RECR[i] = float(p.stdout.readline().strip())
    RUH[i] = float(p.stdout.readline().strip())
    RLH[i] = float(p.stdout.readline().strip())
    #n_e [j][i] = float(p.stdout.readline().strip())



# Cyclotron frequency and plasma frequency output from C++
omega_c = np.zeros((runs,N))
omega_p_2 = np.zeros((runs,N))
for i in range(N):
    for j in range(runs):
        omega_c[j][i] = float(p.stdout.readline().strip())
        omega_p_2[j][i] = float(p.stdout.readline().strip())

# Wave vector: k
k1 = np.zeros((runs,N))
k2 = np.zeros((runs,N))
k3 = np.zeros((runs,N))
for i in range(N):
    for j in range(runs):
        k1[j][i] = float(p.stdout.readline().strip())
        k2[j][i] = float(p.stdout.readline().strip())
        if dim == 3:
            k3[j][i] = float(p.stdout.readline().strip())
        #print("k-vector: (" + str(k1[i]) + ", " + str(k2[i]) + ")")

# Position vector: r
x = np.zeros((runs,N))
y = np.zeros((runs,N))
z = np.zeros((runs,N))
outOfBounds = [0]*runs
for i in range(N):
    for j in range(runs):
        x[j][i] = float(p.stdout.readline().strip())
        y[j][i] = float(p.stdout.readline().strip())
        z[j][i] = float(p.stdout.readline().strip())
        #if (y[j][i] >= L_y or y[j][i] <= L_y_min):  # Out of Bounds Z-axis
        #    outOfBounds[j] = 1
        #if (math.sqrt(x[j][i]**2 + z[j][i]**2) >= L_x or math.sqrt(x[j][i]**2 + z[j][i]**2) <= L_x_min):
        #    outOfBounds[j] = 1              # Out of Bounds radially
        if (math.sqrt( (math.sqrt(x[j][i]**2 + z[j][i]**2)- (L_x+L_x_min)/2)**2 + (y[j][i] - (L_y + L_y_min)/2)**2)) >= (L_x - L_x_min)/2: 
            outOfBounds[j] = 1
        if outOfBounds[j]:
            
            if OBIndex[j] == N - 1:
                print(math.sqrt( (math.sqrt(x[j][i]**2 + z[j][i]**2)- (L_x-L_x_min)/2)**2 + (y[j][i] - (L_y - L_y_min)/2)**2))
                OBIndex[j] = i
            x[j][i] = float("NaN")
            y[j][i] = float("NaN")
            z[j][i] = float("NaN")
            k1[j][i] = float("NaN")
            k2[j][i] = float("NaN")
            omega_c[j][i] = float("NaN")
            omega_p_2[j][i] = float("NaN")
            #n_e[j][i] = float("NaN")


k_theoX = np.zeros((runs,N))
k_theoB = np.zeros((runs,N))
k_theoO = np.zeros((runs,N))
for i in range(N):
    for j in range(runs):
        k_theoX[j][i] = float(p.stdout.readline().strip())
        k_theoB[j][i] = float(p.stdout.readline().strip())
        k_theoO[j][i] = float(p.stdout.readline().strip())


#--------------------------------------------------------------------------------#
# Plotting

# From grid coordinates to real space. Used in Contour
N_nGrid = 100
slopex = (N_nGrid)/(L_x - L_x_min);
interceptx = (N_nGrid) - slopex*(L_x + 0.01*L_x);
slopey = (N_nGrid)/(L_y - L_y_min);
intercepty = (N_nGrid) - slopey*(L_y + 0.01*L_x);
xs = [0]*N_nGrid
ys = [0]*N_nGrid
for i in range(N_nGrid):
    xs[i] = (i - interceptx)/slopex
    ys[i] = (i - intercepty)/slopey



Rsteps1 = np.linspace(L_x_min,L_x,ECRsteps)
Rsteps = np.ones(N)*3.5
Zsteps = np.linspace(L_y_min,L_y,ECRsteps)
R = np.sqrt(x**2 + z**2)
for i in range(runs):
    if (mode[i] == 1):
        plt.plot(R[i],y[i],'--',color = 'tab:red')       # PLOTTING RAY FOR EACH RUN
    if (mode[i] == -1):
        plt.plot(R[i],y[i],'--',color = 'tab:orange')       # PLOTTING RAY FOR EACH RUN
        #plt.plot(R[i],y[i],'--')
    if (mode[i] == 0):
        plt.plot(R[i],y[i],'--',color = 'tab:orange')
        

    #plt.plot(R[i],y[i], label = "r" + str(i + 1))       # PLOTTING RAY FOR EACH RUN
#plt.plot(Rsteps,Zsteps, '--', color = "c",label = '$\omega_{ce} = \omega$')                                # PLOTTING ECR
#plt.plot(RECR,Zsteps,'--', color = 'c', label = '$\omega_{ce} = \omega$')
#plt.plot(RUH,Zsteps,'--', color = 'tab:blue', label = '$\omega_{pe}^2 + \omega_{ce}^2 = \omega^2$')
#plt.plot(RLH,Zsteps,'--', color = 'm', label = '$\omega_{pe}^2 = \omega(\omega - \omega_{ce}$)')

plt.plot(RECR,Zsteps,'--', color = 'c', label = 'ECR')
plt.plot(RUH,Zsteps,'--', color = 'tab:blue', label = 'UHR')
plt.plot(RLH,Zsteps,'--', color = 'm', label = '$R_{Stix} = 0$')

BC = plt.Circle(((L_x+L_x_min)/2,(L_y+L_y_min)/2),(L_x-L_x_min)/2,color = 'k', fill = False)
plt.gca().add_patch(BC)

plt.axis([L_x_min - 0.1, L_x + 0.1,L_y_min - 0.1, L_y + 0.1])
plt.legend(bbox_to_anchor = (0.0, 1), loc = 'upper left')
plt.contour(xs, ys, nGrid.T)

plt.xlabel("$R [m]$")
plt.ylabel("$Z [m]$")
plt.title("Cross Section Trajectory in Overdense Plasma")
plt.show()

# PLOTTING RAY FOR EACH RUN
for i in range(runs):
    if (mode[i] == 1):
        plt.plot(R[i],y[i],'--',color = 'tab:red')
    if (mode[i] == -1):
        plt.plot(R[i],y[i],'--',color = 'tab:orange')
        #plt.plot(R[i],y[i],'--')
    if (mode[i] == 0):
        plt.plot(R[i],y[i],'-',color = 'tab:orange')
        

plt.plot(RECR,Zsteps,'--', color = 'c', label = 'ECR')
plt.plot(RUH,Zsteps,'--', color = 'tab:blue', label = 'UHR')
plt.plot(RLH,Zsteps,'--', color = 'm', label = '$R_{Stix} = 0$')

BC = plt.Circle(((L_x+L_x_min)/2,(L_y+L_y_min)/2),(L_x-L_x_min)/2,color = 'k', fill = False)
plt.gca().add_patch(BC)

plt.axis([3.9,4.7,-0.5,0.5])
plt.legend(bbox_to_anchor = (0.0, 1), loc = 'upper left')
plt.contour(xs, ys, nGrid.T)

plt.xlabel("$R [m]$")
plt.ylabel("$Z [m]$")
plt.title("Cross Section Trajectory in Overdense Plasma")
plt.show()


steps = np.linspace(0,1,100)
steps2 = np.linspace(0,1.5,100)
# R = 0
K1 = 1 - steps

# L = 0?
K2 = np.ones(len(steps))
for i in range(len(steps)):
    K2[i] = math.sqrt(1 - steps[i])
    
# L = 0
Kx = (steps**2 - 1)/(steps-1.000001)
ones = np.ones(len(steps))


for i in range(runs):
    if (mode[i] == 1):
        plt.plot(omega_p_2[i]/omega**2, omega_c[i]/omega,'--',color = 'tab:red')
    if (mode[i] == 0):
        plt.plot(omega_p_2[i]/omega**2, omega_c[i]/omega,'--',color = 'tab:orange')
    if (mode[i] == -1):
        plt.plot(omega_p_2[i]/omega**2, omega_c[i]/omega,'--',color = 'tab:orange')
    #ax.plot(omega_p_2[i][:OBIndex[i]]/omega**2, omega_c[i][:OBIndex[i]]/omega, label = "Ray angle: " + str(round(alphas[i],2)))
plt.plot(steps,K1, label = 'R = 0')
plt.plot(steps,K2, label = 'S = 0')
plt.plot(Kx,5*steps, label = 'L = 0')
plt.plot(ones,5*steps)
plt.plot(steps2,ones)
#ax.plot(steps,5*ones, label = r'L $\rightarrow \infty$')

plt.title(r'CMA Diagram')
plt.xlabel('X')
plt.ylabel('Y',rotation = 0)
plt.axis([0,1.2,0,2])
#legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')

plt.show()

'''
plt.plot(R.T,omega_p_2.T*m_e*epsilon_0/q**2,color = "tab:blue")
plt.title("Electron Density vs Radius")
plt.xlabel("$R\, [m]$")
plt.ylabel("$n_e\, [m^{-3}]$")
#plt.plot(x.T,omega_p_2.T/omega**2)
#plt.title("X vs. $r_x [m]$")
plt.show()


plt.plot(R.T,omega_c.T*m_e/q, color = 'g')
plt.title("Magnetic Field Strength vs Radius")
plt.xlabel("$R\, [m]$")
plt.ylabel("$B\, [T]$")
#plt.plot(x.T,omega_c.T/omega)
#plt.title("Y vs. $r_x [m]$")
plt.show()
'''
plt.plot(xs,nGrid)
plt.title("Electron Density vs Radius")
plt.xlabel("$R\, [m]$")
plt.ylabel("$n_e\, [m^{-3}]$")
plt.show()

plt.plot(xs,TGrid)
plt.title("Temperature vs Radius")
plt.xlabel("$R\, [m]$")
plt.ylabel("$T\, [eV]$")
plt.show()
'''
plt.plot(xs,nGrid/1e20)
plt.title("Electron Density and Temperature with Island")
plt.xlabel("$R\, [m]$")
plt.ylabel("H")
plt.legend(["$n_e = H\,10^{20}\, [m^{-3}]$","$T = 1.602H\,10^{-19}\, [eV]$"])
plt.show()
'''

k = np.zeros((runs,N))
#ktheo = np.zeros((runs,N))
for i in range(runs):
    k[i] = np.sqrt(k1[i]**2 + k2[i]**2 + k3[i]**2)
#    ktheo[i] = np.sqrt((omega**2 - omega_p_2[i])/c**2)

# K THEORETICAL USED FROM 
for i in range(runs):
    #plt.plot(R[i],k_theoX[i],lw = 2,label = "Theoretical X: " + str(i + 1))
    #plt.plot(R[i],k_theoB[i],lw = 2,label = "Theoretical EBW: " + str(i + 1))
    
    if (mode[i] == 1):
        plt.plot(R[i],k[i], color = 'tab:red')
        plt.plot(R[i],k_theoO[i], '--', color = 'k')
    if (mode[i] == 0):
        plt.plot(R[i],k[i], color = 'y')
        plt.plot(R[i],k_theoB[i], '--', color = 'tab:blue')
        plt.plot(R[i],k_theoX[i], '--', color = 'm')
    if (mode[i] == -1):
        plt.plot(R[i],k[i], color = 'tab:orange')
        #plt.plot(R[i],k_theoB[i], '--', color = 'tab:blue')
        plt.plot(R[i],k_theoX[i], '--', color = 'k')
        


plt.title("Wavenumber vs Radius in Overdense Plasma")
plt.legend(['Simulation','Theoretical EBW','Theoretical X'], loc = 'upper right')
plt.legend(['Simulation O','Theoretical O', 'Simulation X','Theoretical X'], loc = 'upper right')
#plt.legend(['Simulation','Theoretical X'], loc = 'upper right')
plt.xlabel("$R\, [m]$")
plt.ylabel("$k \, [m^{-1}]$")
plt.show()

#for i in range(runs):
#    plt.plot(x[i],n_e[i],label = "$n_e$ " + str(i + 1 ))

#plt.title("$n_e$ vs $R$")
#plt.legend(bbox_to_anchor = (0.0, 1), loc = 'upper left')
#plt.xlabel("$x\, [m]$")
#plt.ylabel("$n_e\, [m^{-3}]$")
#plt.show()

'''
plt.contour(xs, ys, nGrid)
plt.title("Electron Density")
plt.xlabel("$R\, [m]$")
plt.ylabel("$Z\, [m]$")
plt.show()

plt.contour(xs, ys, TGrid)
plt.title("Temperature Distribution")
plt.xlabel("$R\, [m]$")
plt.ylabel("$Z\, [m]$")
plt.show()
'''

'''
plt.contour(xs, ys, BGrid)
plt.title("Contour plot $B$ Grid")
plt.xlabel("$R\, [m]$")
plt.ylabel("$Z\, [m]$")
plt.show()
'''


plt.plot(x.T,z.T,color = 'tab:red')
BC = plt.Circle((0,0),L_x_min,color = 'k', fill = False)
plt.gca().add_patch(BC)
BC2 = plt.Circle((0,0),L_x,color = 'k', fill = False)
plt.gca().add_patch(BC2)
#BC3 = plt.Circle(((L_x+L_x_min)/2,(L_y+L_y_min)/2),(L_x-L_x_min)/2,color = 'k', fill = False)
#plt.gca().add_patch(BC3)
plt.axis([-L_x,L_x,-L_x,L_x])
plt.title("Toroidal Position")
plt.xlabel("$x\, [m]$")
plt.ylabel("$z\, [m]$")
plt.show()

plt.plot(R.T)
plt.title("Radius R vs Iterations N")
plt.xlabel("$N$")
plt.ylabel("$R\, [m]$")
plt.show()



plt.plot(R.T,k_theoB.T/1e3,color = 'b', label = 'EBW')
plt.plot(R.T,k_theoX.T/1e3,color = 'r',label = "Slow X-mode")
plt.title("$k\, [mm^{-1}]$ vs $R\, [m]$")
plt.xlabel("$R\, [m]$")
plt.ylabel("$k\, [mm^{-1}]$")
plt.legend(['EBW','Slow X-mode'], loc = 'upper right')
plt.show()
