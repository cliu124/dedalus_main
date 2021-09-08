import numpy as np
from dedalus import public as de
import time


class flag(object):
    
    def __init__(self):
        self.Lx=np.pi
        self.Lz=np.pi
        self.Ly=np.pi
        self.Nx=16
        self.Ny=16
        self.Nz=16
        self.spectral_z='Fourier'
        self.name='test'
        self.flow='not_defined'
        
        self.Ra_ratio=1# the parameter for IFSC
        
        self.post_store_dt=1
        self.stop_sim_time=1
        self.ks=1# parameter for the large scale shear in IFSC with shear
        self.F_sin=1# amplitude for the large scale shear in IFSC with shear
        self.F_sin_2ks=0
        self.F_sin_3ks=0
        self.F_sin_4ks=0
        
        self.phase_2ks=0
        self.phase_3ks=0
        self.phase_4ks=0
        
        self.current_path='./'#This is the current folder path that might need to be specified is run on cluster
    
    def print_screen(self,logger):
        flag_attrs=vars(self)
        #print(', '.join("%s: %s, \n" % item for item in flag_attrs.items()))
        logger.info(', Attributes: Value,\n,')
        logger.info(', '.join("%s: %s, \n" % item for item in flag_attrs.items()))

    def print_file(self):
        flag_text=open(self.current_path+'/analysis'+'/flag.txt','w+')
        flag_attrs=vars(self)
        print(', Attributes: 123,\n ------\n-------\n------',file=flag_text)
        print(', test: 123,',file=flag_text)
        print(', '+', '.join("%s: %s, \n" % item for item in flag_attrs.items()),file=flag_text)
        flag_text.close()
        
    def build_domain(self):
        if self.flow in ['IFSC_2D_without_shear', 'IFSC_2D_with_shear']:
            x_basis = de.Fourier('x', self.Nx, interval=(0,self.Lx), dealias=3/2)
            z_basis = de.Fourier('z', self.Nz, interval=(0,self.Lz), dealias=3/2)
            domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)
        return domain

    def governing_equation(self,domain):
        if self.flow in ['IFSC_2D_without_shear','IFSC_2D_with_shear']:
            problem = de.IVP(domain,variables=['p','u','w','S','T'])
            problem.parameters['Ra_ratio']=self.Ra_ratio
            if self.flow == 'IFSC_2D_without_shear':
                problem.add_equation("- (dx(dx(u))+dz(dz(u)) ) +dx(p) = 0", condition="(nx!=0) or (nz!=0)")
            elif self.flow == 'IFSC_2D_with_shear':
                ##specify the background shear... this is kolmogorov type shear... 
                ##This is the amplitude and wavenumber of the fundamental frequency forcing
                problem.parameters['ks']=self.ks
                problem.parameters['F_sin']=self.F_sin
                
                ##Amplitude of the other frequency
                problem.parameters['F_sin_2ks']=self.F_sin_2ks
                problem.parameters['F_sin_3ks']=self.F_sin_3ks
                problem.parameters['F_sin_4ks']=self.F_sin_4ks
                
                ##phase of other frequency
                problem.parameters['phase_2ks']=self.phase_2ks
                problem.parameters['phase_3ks']=self.phase_3ks
                problem.parameters['phase_4ks']=self.phase_4ks
                problem.add_equation("- (dx(dx(u))+dz(dz(u)) ) +dx(p) = F_sin*sin(ks*z)+F_sin_2ks*sin(2*ks*z+phase_2ks)+F_sin_3ks*sin(3*ks*z+phase_3ks)+F_sin_4ks*sin(4*ks*z+phase_4ks)", condition="(nx!=0) or (nz!=0)")

            problem.add_equation("p=0",condition="(nx==0) and (nz==0)")
            problem.add_equation("u=0",condition="(nx==0) and (nz==0)")
            problem.add_equation("dx(u)+dz(w)=0",condition="(nx!=0) or (nz!=0)")
            problem.add_equation("dt(S) - (dx(dx(S)) + dz(dz(S))) +w =-u*dx(S)-w*dz(S) ")
            problem.add_equation(" - ( dx(dx(w)) + dz(dz(w)) ) + dz(p) -(T-S*Ra_ratio)  =0")
            problem.add_equation(" - ( dx(dx(T)) + dz(dz(T)) ) + w =0")
        # This is assumed to be doubly periodic, no boundary conditions
        #elif self.flow == "IFSC_2D_with_shear":
        #    problem = de.IVP(domain,variables=['p','u','w','S','T'])
            
        #    problem.add_equation("p=0",condition="(nx==0) and (nz==0)")
        #    problem.add_equation("u=0",condition="(nx==0) and (nz==0)")
        #    problem.add_equation("dx(u)+dz(w)=0",condition="(nx!=0) or (nz!=0)")
        #    problem.add_equation("dt(S) - (dx(dx(S)) + dz(dz(S))) +w =-u*dx(S)-w*dz(S) ")
        #    problem.add_equation(" - ( dx(dx(w)) + dz(dz(w)) ) + dz(p) -(T-S*Ra_ratio)  =0")
        #    problem.add_equation(" - ( dx(dx(T)) + dz(dz(T)) ) + w =0")
            #this is triple periodic, no boundary conditions.
        elif self.flow == "channel":
            problem.parameters['Re']=self.Re
            problem.add_equation("dx(u) + dy(v)+wz=0")
            problem.add_equation("dt(u) - 1/Re*(dx(dx(u)) + dy(dy(u)) + dz(uz) ) + dx(p) = -u*dx(u) - v*dy(u) - w*uz + 2/Re")
            problem.add_equation("dt(v) - 1/Re*(dx(dx(v)) + dy(dy(v)) + dz(vz) ) + dy(p) = -u*dx(v) - v*dy(v) - w*vz")
            problem.add_equation("dt(w) - 1/Re*(dx(dx(w)) + dy(dy(w)) + dz(wz) ) + dz(p) = -u*dx(w) v*dy(w) - w*wz")
            problem.add_equation("uz - dz(u) =0")
            problem.add_equation("vz - dz(v) =0")
            problem.add_equation("wz - dz(w) =0")
            problem.add_bc("left(u) = 0")
            problem.add_bc("left(v) = 0")
            problem.add_bc("left(w) = 0")
            problem.add_bc("right(u) = 0")
            problem.add_bc("right(v) = 0")
            problem.add_bc("right(w) = 0", condition="(nx !=0) or (ny !=0)")
            problem.add_bc("right(p) = 0", condition="(nx == 0) and (ny == 0)")
        elif self.flow=='not_defined':
            raise TypeError('flag.flow is not defined yet') 
        return problem

    def initial_condition(self,domain,solver):
        #This initial condition also need to be modified
        if self.flow in ['IFSC_2D_without_shear', 'IFSC_2D_with_shear']:
    
            x = domain.grid(0)
            z = domain.grid(1)
            u = solver.state['u']
            w = solver.state['w']
            S = solver.state['S']
            p = solver.state['p']
            T = solver.state['T']
            
            gshape = domain.dist.grid_layout.global_shape(scales=1)
            slices = domain.dist.grid_layout.slices(scales=1)
            rand = np.random.RandomState(seed=23)
            noise = rand.standard_normal(gshape)[slices]
            
            k_opt=(1/2*(-2-self.Ra_ratio+np.sqrt(self.Ra_ratio**2+8*self.Ra_ratio)))**(1/4)
            w['g'] = self.A_elevator*np.sin(k_opt*x) + self.A_noise*noise
            u['g'] = self.A_noise*noise \
                + self.A_shear*self.F_sin/self.ks**2*np.sin(self.ks*z)\
                + self.A_shear*self.F_sin_2ks/(2*self.ks)**2*np.sin(2*self.ks*z+self.phase_2ks) \
                + self.A_shear*self.F_sin_3ks/(3*self.ks)**2*np.sin(3*self.ks*z+self.phase_3ks) \
                + self.A_shear*self.F_sin_4ks/(4*self.ks)**2*np.sin(4*self.ks*z+self.phase_4ks)
                
            S['g'] = -1/self.Ra_ratio*(k_opt**2+1/k_opt**2)*self.A_elevator*np.sin(k_opt*x) + self.A_noise*noise
            T['g'] = -1/(k_opt**2)*self.A_elevator*np.sin(k_opt*x) + self.A_noise*noise
            p['g'] = self.A_noise*noise

    def run(self,solver,cfl,domain,logger):
        ##This CFL condition need to be modified for different simulation configuration.
        if self.flow in ['IFSC_2D_without_shear', 'IFSC_2D_with_shear']:
            cfl.add_velocities(('u','w'))

        logger.info('Starting loop')
        solver.stop_wall_time = np.inf
        solver.stop_iteration = np.inf
        start_time=time.time()
        while solver.ok:
            dt = cfl.compute_dt()    
            solver.step(dt)
            if solver.iteration % 10 == 0:
                logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
        
        end_time = time.time()
        # Print statistics
        logger.info('Run time: %f' %(end_time-start_time))
        logger.info('Iterations: %i' %solver.iteration)
        logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))

        
    def post_store(self,solver):
        #This post-processing variable need to be modified for different flow configuration
        if self.flow in ['IFSC_2D_without_shear','IFSC_2D_with_shear']:
            analysis = solver.evaluator.add_file_handler('analysis',sim_dt=self.post_store_dt)
            analysis.add_task('S',layout='g',name='S')
            analysis.add_task('T',layout='g',name='T')
            analysis.add_task('u',layout='g',name='u')
            analysis.add_task('w',layout='g',name='w')
            
            analysis.add_task("S",layout='c',name='S_coeff')
            analysis.add_task("T",layout='c',name='T_coeff')
            analysis.add_task("u",layout='c',name='u_coeff')
            analysis.add_task("w",layout='c',name='w_coeff')

            #analysis.add_system(solver.state,layout = 'c')

