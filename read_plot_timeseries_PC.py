import pencil_old as pc
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#import pylab as plt
from scipy.signal import savgol_filter
matplotlib.use('TkAgg')

param=pc.read_param()

ts=pc.read_ts() #pencil_old command
#df = pd.read_csv("latest_run-1446487.csv")

#import constants for orbiters

#gamma=param.gamma                 #adiabatic index
gamma=1                            #isothermal index
alpha=param.density_power_law       #surface density index
beta=param.temperature_power_law    #temp index 
b=param.r_smooth[1]                 #Softening parameter - distance btwn two objects and role of gravity there
cs0=param.cs0                       #initial sound speed. cs = sound speed! 
rho0=param.rho0                     #initial density
gravity=param.g0                    #Gravity constant - Gsmbhmass_live=param_live.pmass[0]          #SMBH Mass
smbhmass=param.pmass[0]             #SMBH Mass
xi=beta-((gamma-1)*alpha)           #initial entropy index - chaos!! 
orbiter1mass=param.pmass[1]         #mass ratio of the first orbiter divided by smbh
cp=param.cp                         #specific heat at constant pressure/thermodynamic constant w/constant pressure
Sigma=rho0                          #surface density near an orbiter is equal to initial density bc its a 2D model
r=ts.xq2       
r2=ts.xq3                           #radial position of second orbiter
Omega=ts.vyq2
vrad = ts.vxq2
v2=Omega**2+vrad**2
Omega3=ts.vyq3
vrad3=ts.vxq3
v2_outerBH=Omega3**2+vrad3**3

cs=cs0                              #sound speeed
vk=np.sqrt(gravity*smbhmass/r)      #Keplerian orbit of orbiter 1
h=cs/vk                             #aspect ratio - how thick the disk is
b_h=b/h                             #smoothing param divided by aspect ratio

Gamma0=np.power(orbiter1mass/h,2)*Sigma*np.power(r,4)*np.power(Omega,2)    #torque normalization constant

innerBH_torque=ts.torqint_2+ts.torqext_2
extBH_torque=ts.torqint_3+ts.torqext_3

#df = pd.DataFrame(innerBH_torque)
#tau_filtered = df[["innerBH_torque"]].apply(savgol_filter,  window_length=len(innerBH_torque), polyorder=2)

#total_innerBH_torque=gamma*tau_filtered/Gamma0
total_innerBH_torque=gamma*innerBH_torque/Gamma0
total_extBH_torque=gamma*extBH_torque/Gamma0

#semi-major axis
a=1./(2./r-v2)
omdot = Omega/r

a2=1./(2./r2-v2_outerBH)
omdot3=Omega3/r2

#eccentricity
h2=r**2*omdot
xtmp=h2**2/a
e2=np.sqrt(1 - xtmp)

h3=r2**2*omdot3
xtmp3=h3**2/a2
e3=np.sqrt(1 - xtmp3)

time=np.array(ts.t)/(2*np.pi)

#fig = plt.figure(figsize=(20,10))
#ax1 = fig.add_subplot(1,2,1)

ax1=plt.subplot(311)
plt.plot(time,a, label='semimajor axis: innerBH')
plt.plot(time,a2, label='semimajor axis: outerBH')
plt.legend(loc=0)
plt.ylabel('a', fontsize=14)
plt.title('noSGfree 2 body run')
plt.ylim(0,3)
plt.tick_params('x', labelbottom=False)

ax2=plt.subplot(312, sharex=ax1)
plt.plot(time,e2,label='inner_BH(free), eccentricity')
plt.plot(time,e3,label='outer_BH(free), eccentricity')
plt.legend()
plt.ylabel('e', fontsize=14)
plt.tick_params('x', labelbottom=False)

ax3=plt.subplot(313, sharex=ax1)
plt.plot(time,total_innerBH_torque,label='inner_BH(LIVE),q=1e-5')  
plt.plot(time,total_extBH_torque,label='inner_BH(LIVE),q=1e-5')  
plt.ylabel('$\gamma\Gamma/\Gamma_0$', fontsize=14)
plt.xlabel ('Time (Orbits)')
plt.legend()

plt.show()
