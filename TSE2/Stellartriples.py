from numpy.linalg import norm
from numpy import sqrt
from scipy import integrate,stats
from tqdm import tqdm
import numpy as np
from numba import jit
import matplotlib.pyplot as plt
import uuid
from isolatedstar import IsolatedCO, IsolatedStar
from binarystar import BinaryStar
from events import *
from units import *
from evolution import evolve
import my_lib
from common import power_law
import sys

# Set to False on the cluster
show_pbar = True

### Random output file name ###
file = str(uuid.uuid4())

### Metallicity ###
z = 0.0002

### Final integration time ###
t_bound = 13e9*yr

### Draw inner binary from Sana's distribution ###
m1 = power_law(-2.3,22.,100.)*Msol
m2 = power_law(-.1,.1,1.)*m1
m12 = m1+m2
ein = power_law(-.45,0.01,0.99)
Pin = power_law(-.55,1e0,1e4)*day
ain = (np.sqrt(Pin/2./np.pi)*m12*G)**(1/3)

### Read tertiary parameters from command line ###
#m3,eout,aout = float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3])
m3 = np.random.uniform(0,100)
eout = power_law(1.,0.01,0.99)
aout = 10.**np.random.uniform(1,3)*ain
m3 *= Msol
#aout *= AU
print(m3/Msol,eout,aout/AU)

f,Om,om,i = np.random.uniform(0,2*np.pi,size=4)
i = np.arccos(np.random.uniform(-1,1))
#i = 91./360.*2*np.pi

initial = np.array([z,m1/Msol,m2/Msol,m3/Msol,ain/AU,aout/AU,ein,eout,f,Om,om,i],dtype='f')

print("###  Initial parameters  ###")
print("%17s" % "Metallicity:","%8.3f" % z)
print("%17s" % "Masses:","%8.3f" % float(m1/Msol),"%8.3f" % float(m2/Msol),"%8.3f" % float(m3/Msol),"Solar masses")
print("%17s" % "SMAs:","%8.3f" % float(ain/AU),"%8.3f" % float(aout/AU),"AU")
print("%17s" % "Eccentricities:","%8.3f" % ein,"%8.3f" % eout)
print("%17s" % "Rel. inclination:","%8.3f" % i)
print("%17s" % "Max. integration:","%8.3f" % float(t_bound/Myr),"Myr")
print("")

globals.star1 = IsolatedStar(z,14,m1/Msol,m1/Msol)
globals.star2 = IsolatedStar(z,14,m2/Msol,m2/Msol)
globals.star3 = IsolatedStar(z,14,m3/Msol,m3/Msol)

y = np.zeros(25)
y[0] = ein
y[5] = sqrt(1.-ein**2)
y[6] = ain
#y[15] = 1.*m1**2*G/c
#y[18] = 1.*m2**2*G/c

m12 = m1+m2
mu = m1*m2/m12
m123 = m12+m3
muout = m12*m3/m123
p = (1.-eout**2)*aout
r = p/(1.+eout*np.cos(f))
u1 = np.array([np.cos(Om),np.sin(Om),0.])
u2 = np.array([-np.cos(i)*np.sin(Om),np.cos(i)*np.cos(Om),np.sin(i)])
y[7:10] = r*(u1*np.cos(f+om)+u2*np.sin(f+om))
y[10:13] = sqrt(G*m123/p)*(-u1*(eout*np.sin(om)+np.sin(f+om))+u2*(eout*np.cos(om)+np.cos(f+om)))
Pout = 2.*np.pi*sqrt(aout**3/G/m123)

y[21] = globals.star1.bcm[0,12+14]*1./yr
y[24] = globals.star2.bcm[0,12+14]*1./yr


pbar = tqdm(total=t_bound)

def evolve_pbar(t,y):
    globals.m1 = globals.star1.get_m(t/Myr)*Msol
    globals.m2 = globals.star2.get_m(t/Myr)*Msol
    globals.m3 = globals.star3.get_m(t/Myr)*Msol
    globals.k1 = globals.star1.get_k(t/Myr)
    globals.k2 = globals.star2.get_k(t/Myr)
    globals.k3 = globals.star3.get_k(t/Myr)
    globals.dm1 = globals.star1.get_dm(t/Myr)*Msol/Myr
    globals.dm2 = globals.star2.get_dm(t/Myr)*Msol/Myr
    globals.dm3 = globals.star3.get_dm(t/Myr)*Msol/Myr
    globals.r1 = globals.star1.get_r(t/Myr)*Rsol
    globals.r2 = globals.star2.get_r(t/Myr)*Rsol
    globals.r3 = globals.star3.get_r(t/Myr)*Rsol
    globals.tco1 = globals.star1.tco*Myr
    globals.tco2 = globals.star2.tco*Myr
    globals.tco3 = globals.star3.tco*Myr
    
    if(show_pbar & (t>pbar.n)):
        pbar.n = int(t)
        pbar.refresh()
    return evolve(t,y)

def evolve_Kepler(t,y,inner_binary):
    m1 = inner_binary.get_m1(t/Myr)*Msol
    m2 = inner_binary.get_m2(t/Myr)*Msol
    m3 = globals.star3.get_m(t/Myr)*Msol
    m123 = m1+m2+m3
    routv = y[0:3]
    rout = norm(routv)
    voutv = y[3:6]
    yp = np.zeros_like(y)
    yp[0:3] += voutv
    yp[3:6] += -G*m123*routv/rout**3 # Keplerian
    return yp

event_list = [merger,stability,PRLO,SRLO,TRLO,PSN,SSN,TSN,TUnbound] # Events to track
event_names = ["DCO merger","Dynamical instability","Primary RLO","Secondary RLO","Tertiary RLO","Primary SN","Secondary SN","Tertiary SN","Unbound outer orbit"]

t_init = 0.
t_eval = np.arange(t_init,t_bound,step_t)
sol = integrate.solve_ivp(evolve_pbar, [t_init,t_bound], y, method="Radau", events=event_list, t_eval=t_eval, rtol=rtol, atol=atol)#, max_step=Pout/100.)
t_final = sol.t
y_final = sol.y

print(sol.message)

# Check if and which event terminated integration
while(sol.status==1):
    maxs = [a.size if a.size==0 else np.max(a) for a in sol.t_events]
    event_ind = np.argmax(maxs)
    event_t = np.max(maxs)
    print("Terminating event:",event_names[event_ind])
    print("Terminating time (Myr):",event_t/Myr)

    if((event_ind==0) | (event_ind==1) | (event_ind==4) | (event_ind==8)):
        print("Simulation finished.")
        break
    ### Temporary to check the donor type ###
    #elif(event_ind==2):
    #    ev = (sol.y).T[-1][0:3]
    #    e = norm(ev)
    #    kw1,m1,rad1,_ = globals.star1.get_basics(event_t/Myr)
    #    kw2,m2,rad2,_ = globals.star2.get_basics(event_t/Myr)
    #    #print(event_t/Myr,kw1,m1,rad1*Rsol/AU,kw2,m2,rad2*Rsol/AU,(sol.y).T[-1][6]/AU,e,initial)
    #    output = np.array([event_t/Myr,kw1,m1,rad1*Rsol/AU,kw2,m2,rad2*Rsol/AU,(sol.y).T[-1][6]/AU,e],dtype='f')
    #    output = np.append(output,initial)
    #    print(output)
    #    with open('RLOs.txt','a') as f: np.savetxt(f, [output])   
    #    print("RLO. Simulation finished.")
    #    break

    elif((event_ind==2) | (event_ind==3)):
        print(sol.y_events[event_ind][-1][6]/Rsol)
        print("SMA (Rsol):",(sol.y).T[-1][6]/Rsol)
        print("Eccentricity:",norm((sol.y).T[-1][0:3]))
        print("LogR:",np.log10(globals.r1/Rsol),np.log10(globals.r2/Rsol))
        
        print("Initiate Roche lobe overflow.")
        
        ev = (sol.y).T[-1][0:3]
        e = norm(ev)
        a = (sol.y).T[-1][6]

        print("LogR:",(globals.r1)/(a*(1.-e)*Roche(globals.m1,globals.m2)),np.log10(globals.r2/Rsol))

        kw1,m1,rad1,dm1 = globals.star1.get_basics(event_t/Myr) # Leave BSE units
        kw2,m2,rad2,dm2 = globals.star2.get_basics(event_t/Myr)
        print(kw1,kw2)
        mass01,lum1,massc1,radc1,menv1,renv1,ospin1,epoch1 = globals.star1.get_full(event_t/Myr) # Leave BSE units
        mass02,lum2,massc2,radc2,menv2,renv2,ospin2,epoch2 = globals.star2.get_full(event_t/Myr)

        kstar = (kw1,kw2)
        mass = (m1,m2)
        rad = (rad1,rad2)
        mass0 = (mass01,mass02)
        lum = (lum1,lum2)
        massc = (massc1,massc2)
        radc = (radc1,radc2)
        menv = (menv1,menv2)
        renv = (renv1,renv2)
        ospin = (ospin1,ospin2)
        epoch = (epoch1,epoch2)
        tms = (0.,0.)

        tphys = event_t/Myr
        dtp = 0.
        tb = 2.*np.pi*np.sqrt(a**3/G/m12)/day

        inner_binary = BinaryStar(z,kstar,mass0,mass,rad,lum,massc,radc,menv,renv,ospin,epoch,tms,tphys,100.,dtp,tb,e)
        tphys,outcome,i  = inner_binary.get_outcome(interaction=True)
        if(outcome!="END_RCHE"): break

        kw1,m1,r1,mass01,lum1,massc1,radc1,menv1,renv1,ospin1,epoch1 = inner_binary.get_star1(tphys,bpp=True,i=i)
        kw2,m2,r2,mass02,lum2,massc2,radc2,menv2,renv2,ospin2,epoch2 = inner_binary.get_star2(tphys,bpp=True,i=i)

        # Update orbital elements
        (sol.y).T[-1][0:3] *= inner_binary.bcm[i,31]/norm((sol.y).T[-1][0:3])
        (sol.y).T[-1][3:6] *= np.sqrt(1.-inner_binary.bcm[i,31]**2)/norm((sol.y).T[-1][3:6])
        (sol.y).T[-1][19:22] = (sol.y).T[-1][3:6]*ospin1/yr/norm((sol.y).T[-1][3:6])
        (sol.y).T[-1][22:25] = (sol.y).T[-1][3:6]*ospin2/yr/norm((sol.y).T[-1][3:6])
        #sol_Kepler = integrate.solve_ivp(evolve_Kepler, [event_t,tphys*Myr], (sol.y).T[-1][7:13], args=(inner_binary), method="Radau", rtol=rtol, atol=atol)#, max_step=Pout/100.)
        #(sol.y).T[-1][7:10] = (sol_Kepler.y).T[-1][0:3]
        #(sol.y).T[-1][10:13] = (sol_Kepler.y).T[-1][3:6]
        if(kw1<10):
            globals.star1 = IsolatedStar(z,kw1,mass01,m1,r1,lum1,massc1,radc1,menv1,renv1,ospin1,epoch1,tphys)
        else:
            globals.star1 = IsolatedCO(kw1,m1)
            # No longer check for Primary SN & RLO
            PSN.terminal = False
            PRLO.terminal = False
        if(kw2<10):
            globals.star2 = IsolatedStar(z,kw2,mass02,m2,r2,lum2,massc2,radc2,menv2,renv2,ospin2,epoch2,tphys)
        else:
            globals.star2 = IsolatedCO(kw2,m2)
            # No longer check for Primary SN & RLO
            SSN.terminal = False
            SRLO.terminal = False

        t_init = tphys*Myr
    elif((event_ind==5) | (event_ind==6)):
        if(event_ind==5):
            subject_star = globals.star1
            companion_mass = m2
            print("Apply kick to the primary.")
        if(event_ind==6):
            subject_star = globals.star2
            companion_mass = m1
            print("Apply kick to the secondary.")

        print("Fallback fraction:",subject_star.ffb) 
        print("CO mass (Msol):",subject_star.mpost)
        print("CO type:",subject_star.kco)
        e = np.array(norm((sol.y).T[-1][0:3]),dtype='f8')
        a = np.array(y[6]/Rsol,dtype='f8')
        print("Pre-kick SMA (AU):",y[6]/AU)
        print("Pre-kick eccentricity:",e)
        kw = np.array(subject_star.kco,dtype='i8')
        mpre = np.array(subject_star.mpre/Msol,dtype='f8')
        mpost = np.array(subject_star.mpost/Msol,dtype='f8')
        m2 = np.array(companion_mass/Msol,dtype='f8')
        my_lib.kicksn.ffb,my_lib.kicksn.directcollapse,my_lib.kicksn.ecs,my_lib.kicksn.mfin = subject_star.ffb,subject_star.directcollapse,subject_star.ecs,subject_star.mfin
        jorb = np.zeros(1,dtype='f8')
        vs = np.zeros(3,dtype='f8')
        my_lib.kicksn.ev = (sol.y).T[-1][0:3]
        my_lib.kicksn.jv = (sol.y).T[-1][3:6]
        my_lib.kick(kw,mpre,mpost,m2,e,a,jorb,vs)
        vs *= kms
        a *= Rsol
        print("Post-kick SMA (AU):",a/AU)
        print("Post-kick eccentricity:",e)
        print("COM kick (km/s):",norm(vs))
        
        (sol.y).T[-1][0:3] = my_lib.kicksn.ev
        (sol.y).T[-1][3:6] = my_lib.kicksn.jv
        (sol.y).T[-1][6] = a
        print("vs",vs)
        print((sol.y).T[-1][10:13])
        (sol.y).T[-1][10:13] += vs
        t_init = event_t

        if((a<0.) | (e<0.) | (e>1.)):
            print("Triple disrupted.")
            break

        if(event_ind==5): globals.star1 = IsolatedCO(globals.star1.kco,globals.star1.mpost)
        if(event_ind==6): globals.star2 = IsolatedCO(globals.star2.kco,globals.star2.mpost)

    elif(event_ind==7):
        print("Apply kick to the tertiary.")
        print("Fallback fraction:",globals.star3.ffb) 
        print("CO mass (Msol):",globals.star3.mpost)
        print("CO type:",globals.star3.kco)
        vk = np.random.uniform(size=3) # random direction
        vk /= norm(vk) # make unit vector
        vk *= stats.maxwell.rvs(scale=vk_sigma) # Maxwellian kick veloctiy
        vk *= (1.-globals.star3.ffb) # apply fallback
        print("Kick velocity (km/s):",norm(vk))
        (sol.y).T[-1][10:13] += vk # apply kick
        globals.star3 = IsolatedCO(globals.star3.kco,globals.star3.mpost)
        t_init = event_t
    else:
        print("Error. Which case are you looking for?",event_ind)
        break

    t_eval = np.arange(t_init,t_bound,step_t)
    sol = integrate.solve_ivp(evolve_pbar, [t_init,t_bound], (sol.y).T[-1], method="Radau", events=event_list, t_eval=t_eval, rtol=rtol, atol=atol)#, max_step=Pout/100.)
    print(sol.message)
    t_final = np.append(t_final,sol.t)
    y_final = np.append(y_final,sol.y,axis=1)

#np.save("output/"+file+".npy",np.concatenate((y_final,[t_final]),axis=0))
np.savetxt("output/"+file+".in",np.append(initial,sol.status))

fig, axs = plt.subplots(1,3,figsize=(12,4))

axs[0].plot(t_final,1.-np.sqrt(y_final[0]**2+y_final[1]**2+y_final[2]**2))
axs[0].set_yscale("log")
axs[0].set_ylabel(r"$1-e$")
axs[0].set_xlabel(r"$t$ [yr]")
axs[0].set_xlim(0.,None)

axs[1].plot(t_final,y_final[6]/AU)
axs[1].set_yscale("log")
axs[1].set_ylabel(r"$a$ [AU]")
axs[1].set_xlabel(r"$t$ [yr]")
axs[1].set_xlim(0.,None)

plt.tight_layout()
plt.savefig("output/"+file+".png")
plt.close()

sys.exit()

axs[0,1].plot(t_final,y_final[6]/AU,"ro")
axs[0,1].set_yscale("log")
axs[0,1].set_ylabel(r"$a$ [AU]")
axs[0,1].set_xlabel(r"$t$ [yr]")
axs[0,1].set_xlim(0.,None)

E = (y_final[10]**2+y_final[11]**2+y_final[12]**2)/2-G*m123/np.sqrt(y_final[7]**2+y_final[8]**2+y_final[9]**2)

axs[0,2].plot(t_final,E/E[0])
axs[0,2].set_yscale("log")
axs[0,2].set_ylabel(r"$E_{out}/E_{out}(0)$")
axs[0,2].set_xlabel(r"$t$ [yr]")
axs[0,2].set_xlim(0.,None)

axs[1,0].plot(y_final[7]/AU,y_final[8]/AU,"r.")
axs[1,1].plot(y_final[7]/AU,y_final[9]/AU,"r.")
axs[1,2].plot(y_final[8]/AU,y_final[9]/AU,"r.")
axs[1,0].set_xlim(-10*aout/AU,10*aout/AU)
axs[1,0].set_ylim(-10*aout/AU,10*aout/AU)
axs[1,1].set_xlim(-10*aout/AU,10*aout/AU)
axs[1,1].set_ylim(-10*aout/AU,10*aout/AU)
axs[1,2].set_xlim(-10*aout/AU,10*aout/AU)
axs[1,2].set_ylim(-10*aout/AU,10*aout/AU)
axs[1,0].set_xlabel(r"$x_{\rm out}$ [AU]")
axs[1,0].set_ylabel(r"$y_{\rm out}$ [AU]")
axs[1,1].set_xlabel(r"$x_{\rm out}$ [AU]")
axs[1,1].set_ylabel(r"$z_{\rm out}$ [AU]")
axs[1,2].set_xlabel(r"$y_{\rm out}$ [AU]")
axs[1,2].set_ylabel(r"$z_{\rm out}$ [AU]")

plt.tight_layout()
plt.savefig("output/"+file+".png")
plt.close()

