import numpy as np
from scipy.integrate import solve_ivp
import sys
import astropy.units as u
from astropy.constants import G,c


# Import termination events from python script located in a different directory
sys.path.append('../src/')
from main import plot, store, evolve_remaining_inner_binary, primary_RL, secondary_RL, tertiary_RL, unstable, primary_SN, secondary_SN, tertiary_SN, DCO_merger, Unphysical, effectively_isolated_DCO, CustomEvent, model_RLO, apply_inner_SN, apply_outer_SN, evolve
import orbittools as ot
from stars import SingleStar, InteractingBinaryStar
from initialconditions import ic

####################################################################################################
# Main routine
####################################################################################################
if __name__ == '__main__':

    # Initial conditions
    m1 = ic.m1
    m2 = ic.m2
    m3 = ic.m3
    a_in = ic.a1
    a_out = ic.a2
    e_in = ic.e1
    e_out = ic.e2
    cos_i_in = ic.cosi1
    omega_in = ic.omega1
    Omega_in = ic.Omega1
    cos_i_out = ic.cosi2
    omega_out = ic.omega2
    Omega_out = ic.Omega2

    lvec_in,evec_in,nvec_in = ot.orbital_elements_to_vectorial_elements(cos_i_in, Omega_in, omega_in)
    lvec_out,evec_out,nvec_out = ot.orbital_elements_to_vectorial_elements(cos_i_out, Omega_out, omega_out)

    print(ic)

    # Print how the inner binary would have evolve if it were isolated
    print('Evolve inner binary as if it was isolated',end='\n\n')
    P_in = ot.orbital_period(ic.a1,m=ic.m1+ic.m2,units=(u.Rsun,u.day,u.Msun))
    _ = InteractingBinaryStar(mass_1=ic.m1,mass_2=ic.m2,mass0_1=ic.m1,mass0_2=ic.m2,
                                period=P_in,ecc=ic.e1,type1=1,type2=1,just_print=True)

    # Initial state vector
    y0 = np.zeros(26)

    y0[0:3] = evec_in*e_in
    y0[3:6] = lvec_in*np.sqrt(1-e_in**2)
    y0[6] = a_in

    y0[7:10] = evec_out*e_out
    y0[10:13] = lvec_out*np.sqrt(1-e_out**2)
    y0[13] = a_out

    # Star objects
    print(' ')
    print('Evolve triple system',end='\n\n')
    print('Evolve single stars (ignore entries for dummy secondary)',end='\n\n')
    star1 = SingleStar(mass0_1=m1)
    star2 = SingleStar(mass0_1=m2)
    star3 = SingleStar(mass0_1=m3)

    # Stellar spins
    y0[14:17] = lvec_in*star1.interpolators['ospin'](0)*1e6
    y0[17:20] = lvec_in*star2.interpolators['ospin'](0)*1e6

    # Store initial conditions
    store(t=0,y=y0,star1=star1,star2=star2,star3=star3,status='Initial')

    # Terminating events
    CustomEvent.terminal = True
    primary_RL.terminal = True
    secondary_RL.terminal = True
    tertiary_RL.terminal = True
    unstable.terminal = True
    primary_SN.terminal = True
    secondary_SN.terminal = True
    tertiary_SN.terminal = True
    DCO_merger.terminal = True
    Unphysical.terminal = True
    effectively_isolated_DCO.terminal = True
    events = [CustomEvent,primary_RL,secondary_RL,tertiary_RL,unstable,primary_SN,secondary_SN,tertiary_SN,DCO_merger,Unphysical,effectively_isolated_DCO]
    event_label = ['Custom event','Primary Roche lobe overflow','Secondary Roche lobe overflow','Tertiary Roche lobe overflow','Unstable','Primary supernova','Secondary supernova','Tertiary supernova','DCO merger','Unphysical','Effectively isolated DCO']

    # Do a short burn-in integration
    sol = solve_ivp(evolve, [0,1e-3], y0, args=(star1,star2,star3), method=ic.method, events=events, atol=ic.atol, rtol=ic.rtol, max_step=ic.max_step)

    if sol.status == -1:
        print('Integration failed at burn-in. Likely unrealistic initial conditions')
        store(t=0,y=y0,star1=star1,star2=star2,star3=star3,status='Failed at burn-in')
        sys.exit()
    elif sol.status == 1:
        print('Termination event at burn-in. Likely unrealistic initial conditions')
        print(sol.t_events)
        store(t=0,y=y0,star1=star1,star2=star2,star3=star3,status='Failed at burn-in')
        sys.exit()
    else:
        sol.status = 2 # Set dummy status to initiate the while loop below

    t_sol = sol.t
    y_sol = sol.y

    m1_sol = star1.interpolators['m'](t_sol)
    m2_sol = star2.interpolators['m'](t_sol)
    m3_sol = star3.interpolators['m'](t_sol)
    logr1_sol = star1.interpolators['logr'](t_sol)
    logr2_sol = star2.interpolators['logr'](t_sol)
    logr3_sol = star3.interpolators['logr'](t_sol)

    while sol.status != 0:
        # Get the last solution
        t0 = t_sol[-1]
        y0 = y_sol[:,-1]

        # Do a short integration
        sol = solve_ivp(evolve, [t0,ic.max_time], y0, args=(star1,star2,star3), method=ic.method, events=events, atol=ic.atol, rtol=ic.rtol, max_step=ic.max_step)

        # Append the solution
        t_sol = np.concatenate((t_sol,sol.t))
        y_sol = np.concatenate((y_sol,sol.y),axis=1)
        m1_sol = np.append(m1_sol,star1.interpolators['m'](sol.t))
        m2_sol = np.append(m2_sol,star2.interpolators['m'](sol.t))
        m3_sol = np.append(m3_sol,star3.interpolators['m'](sol.t))
        logr1_sol = np.append(logr1_sol,star1.interpolators['logr'](sol.t))
        logr2_sol = np.append(logr2_sol,star2.interpolators['logr'](sol.t))
        logr3_sol = np.append(logr3_sol,star3.interpolators['logr'](sol.t))

        # If integration terminated, check which termination event
        if sol.status == 1:
            # Check which entry of sol.t_events is not empty
            for i in range(len(sol.t_events)):
                if sol.t_events[i].size != 0:
                    break

            print('A termination event occured.',end='\n\n')

            if i == 0:
                print('Custom event at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Custom event')
                store(t=sol.t_events[i][0],y=sol.y_events[i][-1],star1=star1,star2=star2,star3=star3,status='Custom event')
                sys.exit()
            elif i == 1:
                print('Primary Roche lobe overflow at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Primary Roche lobe overflow')
                store(t=sol.t_events[i][0],y=sol.y_events[i][-1],star1=star1,star2=star2,star3=star3,status='Begin primary RLO')
                t_new,y_new,event_status,star1,star2 = model_RLO(sol.t_events[i][0],sol.y[:,-1],star1,star2)
                if event_status == -1:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End primary RLO (no bpp)')
                    sys.exit()
                elif event_status == -2:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End primary RLO (merger)')
                    sys.exit()
                elif event_status == -3:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End primary RLO (not intact)')
                    sys.exit()
                else:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End primary RLO')
            elif i == 2:
                print('Secondary Roche lobe overflow at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Secondary Roche lobe overflow')
                t_new,y_new,event_status,star1,star2 = model_RLO(sol.t_events[i][0],sol.y[:,-1],star1,star2)
                if event_status == -1:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End secondary RLO (no bpp)')
                    sys.exit()
                elif event_status == -2:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End secondary RLO (merger)')
                    sys.exit()
                elif event_status == -3:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End secondary RLO (not intact)')
                    sys.exit()
                else:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End secondary RLO')
            elif i == 3:
                print('Tertiary Roche lobe overflow at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Tertiary Roche lobe overflow')
                store(t=sol.t_events[i][0],y=sol.y_events[i][-1],star1=star1,star2=star2,star3=star3,status='Tertiary RLO')
                sys.exit()
            elif i == 4:
                print('Unstable at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Unstable')
                store(t=sol.t_events[i][0],y=sol.y_events[i][-1],star1=star1,star2=star2,star3=star3,status='Unstable')
                sys.exit()
            elif i == 5:
                ic.primary_SN = True
                print('Primary supernova at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Primary supernova')
                t_new,y_new,event_status = apply_inner_SN(sol.t_events[i][0],sol.y[:,-1],star1,star2,star3)
                if star1.interpolators['k'](t_new) == 14:
                    y_new[20:23] = y_new[14:17]/np.linalg.norm(y_new[14:17])*G*star1.interpolators['m'](t_new)**2/c # BH spins
                if event_status == -1:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Primary SN (inner binary unbound)')
                    sys.exit()
                elif event_status == -2:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Primary SN (outer binary unbound)')
                    evolve_remaining_inner_binary(t_new,y_new,star1,star2,star3)
                    sys.exit()
                else:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Primary SN')
            elif i == 6:
                ic.secondary_SN = True
                print('Secondary supernova at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Secondary supernova')
                t_new,y_new,event_status = apply_inner_SN(sol.t_events[i][0],sol.y[:,-1],star2,star1,star3)
                if star2.interpolators['k'](t_new) == 14:
                    y_new[23:26] = y_new[17:20]/np.linalg.norm(y_new[17:20])*G*star2.interpolators['m'](t_new)**2/c # BH spins
                if event_status == -1:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Secondary SN (inner binary unbound)')
                    sys.exit()
                elif event_status == -2:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Secondary SN (outer binary unbound)')
                    evolve_remaining_inner_binary(t_new,y_new,star1,star2,star3)
                    sys.exit()
                else:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Secondary SN')
            elif i == 7:
                ic.tertiary_SN = True
                print('Tertiary supernova at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Tertiary supernova')
                t_new,y_new,event_status = apply_outer_SN(sol.t_events[i][0],sol.y[:,-1],star1,star2,star3)
                if event_status == -1:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Tertiary SN (outer binary unbound)')
                    evolve_remaining_inner_binary(t_new,y_new,star1,star2,star3)
                    sys.exit()
                else:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Tertiary SN')
            elif i == 8:
                print('DCO merger at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='DCO merger')
                store(t=sol.t_events[i][0],y=sol.y_events[i][-1],star1=star1,star2=star2,star3=star3,status='DCO merger')
                sys.exit()
            elif i == 9:
                print('Unphysical solution at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Unphysical solution')
                store(t=sol.t_events[i][0],y=sol.y_events[i][-1],star1=star1,star2=star2,star3=star3,status='Unphysical solution')
                sys.exit()
            elif i == 10:
                print('Effectively isolated DCO at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Effectively isolated DCO')
                store(t=sol.t_events[i][0],y=sol.y_events[i][-1],star1=star1,star2=star2,star3=star3,status='Effectively isolated DCO')
                sys.exit()

            t_sol = np.append(t_sol,t_new)
            y_sol = np.append(y_sol,y_new.reshape(-1,1),axis=1)
            m1_sol = np.append(m1_sol,star1.interpolators['m'](t_new))
            m2_sol = np.append(m2_sol,star2.interpolators['m'](t_new))
            m3_sol = np.append(m3_sol,star3.interpolators['m'](t_new))
            logr1_sol = np.append(logr1_sol,star1.interpolators['logr'](t_new))
            logr2_sol = np.append(logr2_sol,star2.interpolators['logr'](t_new))
            logr3_sol = np.append(logr3_sol,star3.interpolators['logr'](t_new))

        elif sol.status == -1:
            print('Integration step failed.')
            plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Integration step failed.')
            store(t=sol.t[-1],y=sol.y[:,-1],star1=star1,star2=star2,star3=star3,status='Integration step failed')
            sys.exit()

    print(sol.message)
    plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Finished integration')
    store(t=t_sol[-1],y=y_sol[:,-1],star1=star1,star2=star2,star3=star3,status='Finished integration')