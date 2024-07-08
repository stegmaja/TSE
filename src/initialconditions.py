import argparse
import numpy as np
import orbittools as ot
import astropy.units as u
import os

class InitialConditions:
    def __init__(self):
        # Using argparse to parse the command line arguments
        parser = argparse.ArgumentParser(description='Initial conditions for the simulation')

        # Masses
        parser.add_argument('--m1', type=float, default=1., help='Initial mass (Msun) of the primary star (default: 1 Msun)')
        parser.add_argument('--m2', type=float, default=1., help='Initial mass (Msun) of the secondary star (default: 1 Msun)')
        parser.add_argument('--m3', type=float, default=1., help='Initial mass (Msun) of the tertiary star (default: 1 Msun)')

        # Semi-major axes
        parser.add_argument('--a1', type=float, default=2e3, help='Initial semi-major axis (Rsun) of the inner binary (default: 2e3 Rsun)')
        parser.add_argument('--a2', type=float, default=2e5, help='Initial semi-major axis (Rsun) of the outer binary (default: 2e5 Rsun)')

        # Eccentricities
        parser.add_argument('--e1', type=float, default=0.1, help='Initial eccentricity of the inner binary (default: 0.1)')
        parser.add_argument('--e2', type=float, default=0.5, help='Initial eccentricity of the outer binary (default: 0.5)')

        # Inclinations
        parser.add_argument('--cosi1', type=float, default=1., help='Initial cosine of the inclination of the inner binary (default: 1)')
        parser.add_argument('--cosi2', type=float, default=1., help='Initial cosine of the inclination of the outer binary (default: 1)')

        # Arguments of pericenter
        parser.add_argument('--omega1', type=float, default=0., help='Initial argument of pericenter of the inner binary (default: 0)')
        parser.add_argument('--omega2', type=float, default=0., help='Initial argument of pericenter of the outer binary (default: 0)')

        # Longitudes of the ascending nodes
        parser.add_argument('--Omega1', type=float, default=np.pi, help='Initial longitude of the ascending node of the inner binary (default: pi)')
        parser.add_argument('--Omega2', type=float, default=0., help='Initial longitude of the ascending node of the outer binary (default: 0)')

        # Integration method
        parser.add_argument('--method', type=str, default='RK45', help='Integration method (see scipy solve_ivp) (default: RK45)')

        # Integration tolerances
        parser.add_argument('--rtol', type=float, default=1e-8, help='Relative tolerance for the integration (default: 1e-8)')
        parser.add_argument('--atol', type=float, default=1e-8, help='Absolute tolerance for the integration (default: 1e-8)')

        # Maximum integration time
        parser.add_argument('--max_time', type=float, default=1e2, help='Maximum integration time (Myr) (default: 1e2 Myr)')

        # Maximum step size
        parser.add_argument('--max_step', type=float, default=1e-2, help='Maximum step size (Myr) (default: 1e-2 Myr)')

        # Metallicity
        parser.add_argument('--Z', type=float, default=0.02, help='Metallicity of the system (default: 0.02)')

        # neta
        parser.add_argument('--neta', type=float, default=0.5, help='Rejuvenation factor (default: 0.5)')

        # bwind
        parser.add_argument('--bwind', type=float, default=0.0, help='Wind mass loss parameter (default: 0.0)')

        # hewind
        parser.add_argument('--hewind', type=float, default=1.0, help='Helium star wind mass loss parameter (default: 1.0)')

        # alpha1
        parser.add_argument('--alpha1', type=float, default=3.0, help='Common envelope efficiency parameter (default: 3.0)')

        # lambda
        parser.add_argument('--lamb', type=float, default=0.1, help='Common envelope parameter (default: 0.1)')

        # ceflag
        parser.add_argument('--ceflag', type=int, default=0, help='Common envelope flag (default: 0)')

        # tflag
        parser.add_argument('--tflag', type=int, default=1, help='Thermally unstable flag (default: 1)')

        # ifflag
        parser.add_argument('--ifflag', type=int, default=0, help='Magnetic braking flag (default: 0)')

        # wdflag
        parser.add_argument('--wdflag', type=int, default=1, help='White dwarf flag (default: 1)')

        # bhflag
        parser.add_argument('--bhflag', type=int, default=2, help='Black hole flag (default: 2)')

        # nsflag
        parser.add_argument('--nsflag', type=int, default=2, help='Neutron star flag (default: 2)')

        # piflag
        parser.add_argument('--piflag', type=int, default=1, help='Pulsar flag (default: 1)')

        # mxns
        parser.add_argument('--mxns', type=int, default=3, help='Maximum number of neutron stars (default: 3)')

        # idum
        parser.add_argument('--idum', type=int, default=678555, help='Random MOBSE seed (default: 678555)')

        # pts1
        parser.add_argument('--pts1', type=float, default=0.05, help='Pulsar term 1 (default: 0.05)')

        # pts2
        parser.add_argument('--pts2', type=float, default=0.01, help='Pulsar term 2 (default: 0.01)')

        # pts3
        parser.add_argument('--pts3', type=float, default=0.02, help='Pulsar term 3 (default: 0.02)')

        # sigma
        parser.add_argument('--sigma', type=float, default=265.0, help='Sigma (default: 265.0)')

        # beta
        parser.add_argument('--beta', type=float, default=0.125, help='Beta (default: 0.125)')

        # xi
        parser.add_argument('--xi', type=float, default=1.0, help='Xi (default: 1.0)')

        # acc2
        parser.add_argument('--acc2', type=float, default=1.5, help='Accretion 2 (default: 1.5)')

        # epsnov
        parser.add_argument('--epsnov', type=float, default=0.001, help='Epsilon novae (default: 0.001)')

        # eddfac
        parser.add_argument('--eddfac', type=float, default=1.0, help='Eddington factor (default: 1.0)')

        # gamma
        parser.add_argument('--gamma', type=float, default=-1.0, help='Gamma (default: -1.0)')

        # seed 
        parser.add_argument('--seed', type=int, default=42, help='Random seed (default: 42)')

        # Create random initial conditions
        parser.add_argument('--random', type=bool, default=False, help='Create random initial conditions from Sana+ 2012? (default: False)')  

        # Should Galactic tides be included?
        parser.add_argument('--galactic_tides', type=bool, default=False, help='Include Galactic tides? (default: False)')

        # Initial x coordinate in the galactic frame
        parser.add_argument('--x_MW', type=float, default=8., help='Initial x coordinate in the galactic frame (default: 10 kpc)')     

        # Initial y coordinate in the galactic frame
        parser.add_argument('--y_MW', type=float, default=0., help='Initial y coordinate in the galactic frame (default: 0 kpc)')

        # Initial z coordinate in the galactic frame
        parser.add_argument('--z_MW', type=float, default=0., help='Initial z coordinate in the galactic frame (default: 0 kpc)')

        # Initial x velocity in the galactic frame
        parser.add_argument('--vx_MW', type=float, default=0., help='Initial x velocity in the galactic frame (default: 0 km/s)')

        # Initial y velocity in the galactic frame
        parser.add_argument('--vy_MW', type=float, default=175., help='Initial y velocity in the galactic frame (default: 175 km/s)')

        # Initial z velocity in the galactic frame
        parser.add_argument('--vz_MW', type=float, default=0., help='Initial z velocity in the galactic frame (default: 0 km/s)') 

        # Should stellar tides be included?
        parser.add_argument('--stellar_tides', type=bool, default=False, help='Include stellar tides? (default: False)')

        args = parser.parse_args()

        np.random.seed(args.seed)

        if args.random:

            args.m1 = ot.power_law_sample(-2.3, 8, 20)
            args.m2 = args.m1*ot.power_law_sample(-0.1, 0.1, 1.0)
            args.m3 = ot.power_law_sample(-2., 0.1, 20)

            P_in = 10**ot.power_law_sample(-0.55, 0.15, 5.5)
            args.a1 = ot.semi_major_axis(P_in,m=args.m1+args.m2,units=(u.Rsun,u.day,u.Msun))
            args.e1 = ot.power_law_sample(-0.42, 0.0, 0.9)

            args.cosi1 = np.random.uniform(-1,1)
            args.omega1 = np.random.uniform(0,2*np.pi)
            args.Omega1 = np.pi

            initial_stability = False

            while not initial_stability:
                args.a2 = 10**np.random.uniform(np.log10(args.a1),
                                            np.log10((1e4*u.AU).to(u.Rsun).value))
                args.e2 = np.sqrt(np.random.uniform(0,1))

                initial_stability = ot.check_triple_stability(args.a1,args.a2,args.e2,args.m1+args.m2,args.m3)

            args.cosi2 = np.random.uniform(-1,1)
            args.omega2 = np.random.uniform(0,2*np.pi)
            args.Omega2 = 0

        self.m1 = args.m1
        self.m2 = args.m2
        self.m3 = args.m3
        self.a1 = args.a1
        self.a2 = args.a2
        self.e1 = max(args.e1,1e-3)
        self.e2 = max(args.e2,1e-3)
        self.cosi1 = args.cosi1
        self.cosi2 = args.cosi2
        self.omega1 = args.omega1
        self.omega2 = args.omega2
        self.Omega1 = args.Omega1
        self.Omega2 = args.Omega2
        self.method = args.method
        self.rtol = args.rtol
        self.atol = args.atol
        self.max_time = args.max_time
        self.max_step = args.max_step
        self.Z = args.Z
        self.neta = args.neta
        self.bwind = args.bwind
        self.hewind = args.hewind
        self.alpha1 = args.alpha1
        self.lamb = args.lamb
        self.ceflag = args.ceflag
        self.tflag = args.tflag
        self.ifflag = args.ifflag
        self.wdflag = args.wdflag
        self.bhflag = args.bhflag
        self.nsflag = args.nsflag
        self.piflag = args.piflag
        self.mxns = args.mxns
        self.idum = args.idum
        self.pts1 = args.pts1
        self.pts2 = args.pts2
        self.pts3 = args.pts3
        self.sigma1 = args.sigma
        self.sigma2 = args.sigma
        self.beta = args.beta
        self.xi = args.xi
        self.acc2 = args.acc2
        self.epsnov = args.epsnov
        self.eddfac = args.eddfac
        self.gamma = args.gamma
        self.seed = args.seed
        self.galactic_tides = args.galactic_tides
        self.x_MW = args.x_MW
        self.y_MW = args.y_MW
        self.z_MW = args.z_MW
        self.vx_MW = args.vx_MW
        self.vy_MW = args.vy_MW
        self.vz_MW = args.vz_MW
        self.stellar_tides = args.stellar_tides

        self.SRC_DIR = os.getcwd()
        self.MOBSE_DIR = self.SRC_DIR+'/../mobse/src'

        self.mobse_input = '../input/' + str(self.seed)+'_'+str(self.Z).replace('.','') + '.input'
        self.mobse_output = '../' + str(self.seed)+'_'+str(self.Z).replace('.','') + '.out'
        self.mobse_log = '../' + str(self.seed)+'_'+str(self.Z).replace('.','') + '_bpp.out'
        
    def __str__(self):
        return f'Initial conditions:\n Seed={self.seed},\n Z={self.Z:.6f},\n m1={self.m1:.6f},\n m2={self.m2:.6f},\n m3={self.m3:.6f},\n a_in={self.a1:.6e},\n e_in={self.e1:.6f},\n a_out={self.a2:.6e},\n e_out={self.e2:.6f},\n cos_i_in={self.cosi1:.6f},\n omega_in={self.omega1:.6f},\n Omega_in={self.Omega1:.6f},\n cos_i_out={self.cosi2:.6f},\n omega_out={self.omega2:.6f},\n Omega_out={self.Omega2:.6f}\n'
    
ic = InitialConditions()