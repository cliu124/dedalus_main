import numpy as np
from dedalus import public as de
import time
import pathlib
from scipy import linalg


class flag(object):
    
    def __init__(self):
        self.Lx=np.pi
        self.Lz=np.pi
        self.Ly=np.pi
        self.Nx=0 ##these default values are zero.. this can avoid some error if I forget to setup the grid point numbers
        self.Ny=0
        self.Nz=0
        
        ##option for self.flow that has been developed:
        #delete the old IFSC_2D_with_shear and IFSC_2D_without_shear option.. these two has been unified together..
        #Update 2021/09/13: also add the option of double_diffusive_2D
        #
        
        #IFSC_2D
        #double_diffusive_2D
        self.flow='not_defined'
        
        self.Ra_ratio=1.1# the parameter for IFSC... fix this default value as 1.1... although not used for double diffusive.. Ra_ratio=1 will cause division by zero error.
        
        self.post_store_dt=1
        self.stop_sim_time=1
        self.ks=1# parameter for the large scale shear in IFSC with shear
        self.F_sin=0# amplitude for the large scale shear in IFSC with shear
        self.F_sin_2ks=0
        self.F_sin_3ks=0
        self.F_sin_4ks=0
        
        self.phase_2ks=0
        self.phase_3ks=0
        self.phase_4ks=0
        
        self.current_path='./'#This is the current folder path that might need to be specified is run on cluster
    
        #These two values as one corresponds to the salf finger
        #If both of them are set as -1, then we have the diffusive regime
        ##IFSC should also apply in that regime...
        self.dy_T_mean=1
        self.dy_S_mean=1
    
        ##These three parameters are used for the double diffusive convection in primitive variable...
        self.R_rho_T2S=1
        self.tau=1
        self.Pr=1
        
        ##These six parameters are used for the double_diffusive_shear_2D. This is the most general formulation that has six parameters
        self.Re=1 #The Reynolds number appearing in front of the inertial term in momentum
        self.Pe_T=1 #The Peclet number appearing in front of the inertial term in temperature
        self.Pe_S=1 #The Peclet number appearing in front of the inertial term in salinity
        self.Ra_T=1 #The Rayleigh number appearing in front of the temperature term, defined as Ra_T=g\alpha T_z L^4/\nu \kappa_T
        self.Ra_S2T=1 #The Rayleigh number appearing in front of the salinity term, this is defined based salintiy over temperature, thus Ra_T=g\beta S_z L^4/\nu \kappa_T
        #self.tau=1 #This tau is not necessary as it has been defined... #This is the diffusivity ratio, \kappa_S/\kappa_T 
        
        self.A_elevator=0
        self.k_elevator=0.5
        self.A_noise=0
        self.A_shear=0        
        
        self.A_secondary_T=0
        self.A_secondary_S=0
        self.k_secondary=0
        
        self.flow_sub_double_diffusive_shear_2D='double_diffusive_2D'
        self.shear_Radko2016_reduced='primitive'
        
        self.lambda_elevator=0
                
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
        if self.flow in ['IFSC_2D','double_diffusive_2D','double_diffusive_shear_2D','porous_media_2D']:
            x_basis = de.Fourier('x', self.Nx, interval=(0,self.Lx), dealias=3/2)
            z_basis = de.Fourier('z', self.Nz, interval=(0,self.Lz), dealias=3/2)
            domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)
        return domain

    def governing_equation(self,domain):
        if self.flow in ['IFSC_2D']:
            problem = de.IVP(domain,variables=['p','u','w','S','T'])
            problem.parameters['Ra_ratio']=self.Ra_ratio
            problem.parameters['dy_T_mean']=self.dy_T_mean
            problem.parameters['dy_S_mean']=self.dy_S_mean
            
            #Update 2021/09/12, change the language to specify the background shear
            #test whether these amplitude of shear is zero....
            if self.F_sin == 0:
                print('without shear')
                problem.add_equation("- (dx(dx(u))+dz(dz(u)) ) +dx(p) = 0", condition="(nx!=0) or (nz!=0)")
            else:
                print('with shear')
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

            problem.add_equation(" - ( dx(dx(w)) + dz(dz(w)) ) + dz(p) -(T-S*Ra_ratio)  =0")            
            problem.add_equation("p=0",condition="(nx==0) and (nz==0)")
            problem.add_equation("u=0",condition="(nx==0) and (nz==0)")
            problem.add_equation("dx(u)+dz(w)=0",condition="(nx!=0) or (nz!=0)")
            problem.add_equation(" - ( dx(dx(T)) + dz(dz(T)) ) + dy_T_mean*w =0")
            problem.add_equation("dt(S) - (dx(dx(S)) + dz(dz(S))) + dy_S_mean*w =-u*dx(S)-w*dz(S) ")

            
        elif self.flow in ['double_diffusive_2D']:
            problem = de.IVP(domain,variables=['p','u','w','S','T'])
            problem.parameters['R_rho_T2S']=self.R_rho_T2S
            problem.parameters['tau']=self.tau
            problem.parameters['Pr']=self.Pr
            problem.parameters['dy_T_mean']=self.dy_T_mean
            problem.parameters['dy_S_mean']=self.dy_S_mean
            
            #Update 2021/09/12, change the language to specify the background shear
            #test whether these amplitude of shear is zero....
            
            ##Note that this is different from the IFSC,,, here I do not need to constraint that (nx!=0) or (nz!=0) because at nx=nz=0, it is just dt(u)=0, a valid equation.. 
            if self.F_sin == 0:
                print('without shear')
                problem.add_equation("dt(u)- Pr*(dx(dx(u))+dz(dz(u)) ) + Pr*dx(p) = -u*dx(u)-w*dz(u)")
            else:
                print('with shear')
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
                problem.add_equation("dt(u) - Pr*(dx(dx(u))+dz(dz(u)) ) +Pr*dx(p) = -u*dx(u)-w*dz(u)+ Pr*(F_sin*sin(ks*z)+F_sin_2ks*sin(2*ks*z+phase_2ks)+F_sin_3ks*sin(3*ks*z+phase_3ks)+F_sin_4ks*sin(4*ks*z+phase_4ks))")

            #problem.add_equation("u=0",condition="(nx==0) and (nz==0)") #Note that for the primitive equation,,, this singularity for u momentum is not there...
            problem.add_equation("dx(u)+dz(w)=0",condition="(nx!=0) or (nz!=0)")
            problem.add_equation("p=0",condition="(nx==0) and (nz==0)")
            problem.add_equation(" dt(w) - Pr*( dx(dx(w)) + dz(dz(w)) ) + Pr*dz(p) -Pr*(T-S/R_rho_T2S)  =-u*dx(w)-w*dz(w)")
            problem.add_equation("dt(S) - tau*(dx(dx(S)) + dz(dz(S))) + dy_S_mean*w =-u*dx(S)-w*dz(S) ")
            problem.add_equation(" dt(T) - ( dx(dx(T)) + dz(dz(T)) ) + dy_T_mean*w =-u*dx(T)-w*dz(T)")
       
        elif self.flow in ['double_diffusive_shear_2D']:
            problem = de.IVP(domain,variables=['p','u','w','S','T'])
            problem.parameters['Re']=self.Re
            problem.parameters['Pe_T']=self.Pe_T
            problem.parameters['Pe_S']=self.Pe_S
            problem.parameters['Ra_T']=self.Ra_T
            problem.parameters['Ra_S2T']=self.Ra_S2T
            problem.parameters['tau']=self.tau
            
            problem.parameters['dy_T_mean']=self.dy_T_mean
            problem.parameters['dy_S_mean']=self.dy_S_mean
            
            #Update 2021/09/12, change the language to specify the background shear
            #test whether these amplitude of shear is zero....
            
            ##Note that this is different from the IFSC,,, here I do not need to constraint that (nx!=0) or (nz!=0) because at nx=nz=0, it is just dt(u)=0, a valid equation.. 
            if self.F_sin == 0:
                print('without shear')
                if self.Re == 0:
                    problem.add_equation("- (dx(dx(u))+dz(dz(u)) ) + dx(p) = 0",condition="(nx!=0) or (nz!=0)")
                    problem.add_equation("u=0",condition="(nx==0) and (nz==0)")
                else:
                    problem.add_equation("Re*dt(u)- (dx(dx(u))+dz(dz(u)) ) + dx(p) = Re*( -u*dx(u)-w*dz(u))")
            else:
                print('with shear')
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
                
                if self.Re == 0:
                    problem.add_equation("- (dx(dx(u))+dz(dz(u)) ) +dx(p) =  (F_sin*sin(ks*z)+F_sin_2ks*sin(2*ks*z+phase_2ks)+F_sin_3ks*sin(3*ks*z+phase_3ks)+F_sin_4ks*sin(4*ks*z+phase_4ks))",condition="(nx!=0) or (nz!=0)")
                    problem.add_equation("u=0",condition="(nx==0) and (nz==0)")
                else:
                    problem.add_equation("Re*dt(u) - (dx(dx(u))+dz(dz(u)) ) +dx(p) = Re*( -u*dx(u)-w*dz(u) )+ (F_sin*sin(ks*z)+F_sin_2ks*sin(2*ks*z+phase_2ks)+F_sin_3ks*sin(3*ks*z+phase_3ks)+F_sin_4ks*sin(4*ks*z+phase_4ks))")

            if self.Re ==0:
                problem.add_equation("- ( dx(dx(w)) + dz(dz(w)) ) + dz(p) -(Ra_T*T-Ra_S2T*S)  =0")
            else:
                problem.add_equation("Re*dt(w)- ( dx(dx(w)) + dz(dz(w)) ) + dz(p) -(Ra_T*T-Ra_S2T*S)  = Re*( -u*dx(w)-w*dz(w) )")

            problem.add_equation("dx(u)+dz(w)=0",condition="(nx!=0) or (nz!=0)")
            problem.add_equation("p=0",condition="(nx==0) and (nz==0)")
            if self.Pe_T == 0:
                problem.add_equation(" - ( dx(dx(T)) + dz(dz(T)) ) + dy_T_mean*w =0")
            else:
                problem.add_equation(" Pe_T*dt(T) - ( dx(dx(T)) + dz(dz(T)) ) + dy_T_mean*w =Pe_T*( -u*dx(T)-w*dz(T) )")

            problem.add_equation("Pe_S*dt(S) - tau*(dx(dx(S)) + dz(dz(S))) + dy_S_mean*w =Pe_S*( -u*dx(S)-w*dz(S) ) ")

        elif self.flow =='porous_media_2D':
            problem = de.IVP(domain,variables=['p','u','w','T'])
            problem.parameters['Ra_T']=self.Ra_T
            problem.add_equation(" u + dx(p) = 0", condition="(nx!=0) or (nz!=0)")
            problem.add_equation(" w + dz(p) - Ra_T*T  =0")            
            problem.add_equation("p=0",condition="(nx==0) and (nz==0)")
            problem.add_equation("u=0",condition="(nx==0) and (nz==0)")
            problem.add_equation("dx(u)+dz(w)=0",condition="(nx!=0) or (nz!=0)")
            problem.add_equation("dt(T) - (dx(dx(T)) + dz(dz(T)))  - w =-u*dx(T)-w*dz(T) ")

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
        else:
            raise TypeError('flag.flow is not defined yet') 
        return problem

    def initial_condition(self,domain,solver):
        if not pathlib.Path('restart.h5').exists():
            print('setup initial condition')
            #This initial condition also need to be modified
            if self.flow in ['IFSC_2D','double_diffusive_2D','double_diffusive_shear_2D']:
        
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
                
                ##Add the random noise
                u0=self.A_noise*noise
                w0=self.A_noise*noise
                S0=self.A_noise*noise
                T0=self.A_noise*noise
                p0=self.A_noise*noise
                #Add the background shear
                u0 = u0+ self.A_shear*self.F_sin/self.ks**2*np.sin(self.ks*z)\
                    + self.A_shear*self.F_sin_2ks/(2*self.ks)**2*np.sin(2*self.ks*z+self.phase_2ks) \
                    + self.A_shear*self.F_sin_3ks/(3*self.ks)**2*np.sin(3*self.ks*z+self.phase_3ks) \
                    + self.A_shear*self.F_sin_4ks/(4*self.ks)**2*np.sin(4*self.ks*z+self.phase_4ks)
                
                #Compute the eigenvalue problem of elevator mode and add the elevator mode into initial condition
                if self.flow=='IFSC_2D':
                    k_opt=(1/2*(-2-self.Ra_ratio+np.sqrt(self.Ra_ratio**2+8*self.Ra_ratio)))**(1/4)
        
                    w0 =w0 +self.A_elevator*np.sin(k_opt*x)
          
                    S0 =S0 -1/self.Ra_ratio*(k_opt**2+1/k_opt**2)*self.A_elevator*np.sin(k_opt*x)
                    T0 =T0 -1/(k_opt**2)*self.A_elevator*np.sin(k_opt*x)
                elif self.flow in ['double_diffusive_2D']:
                    k2=self.k_elevator**2
                    
                    #This A is the linear system for the eigenvalue problem... We already take the elevator mode, so that k_y=0 (vertical wavenumber)
                    #This should work for both finger regime and diffusive regime..
                    A=[[-k2*self.Pr, self.Pr, -self.Pr/self.R_rho_T2S],
                        [-self.dy_T_mean, -k2, 0],
                        [-self.dy_S_mean, 0, -self.tau*k2]];
                    #Compute eigenvalue and eigenvector of 
                    
                    eig_val,eig_vec=linalg.eig(A) #use linear algebra package to compute eigenvalue
                    eig_val_max_ind=np.argmax(np.real(eig_val)) #compute the index of the eigenvalue
                    eig_vec_max=eig_vec[:,eig_val_max_ind] #get the corresponding eigen vector
                    
                    self.lambda_elevator=eig_val[eig_val_max_ind]
                    w0 =w0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[0]) #set the results weighted by the corresponding eigenvector 
                    T0 =T0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[1])
                    S0 =S0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[2])
                    #print(w0)
                    ##This is sample code to setup the initial condition as whatever I want...
                    #S0=np.array([1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16])
                    #S['g']=S0[slices]
                elif self.flow in ['double_diffusive_shear_2D']:
                    k2=self.k_elevator**2
                    
                    #This A is the linear system for the eigenvalue problem... We already take the elevator mode, so that k_y=0 (vertical wavenumber)
                    #This should work for both finger regime and diffusive regime..
                    #This is the eigenvalue problem only for no shear case....
                    A=[[-k2, self.Ra_T, -self.Ra_S2T],
                        [-self.dy_T_mean, -k2, 0],
                        [-self.dy_S_mean, 0, -self.tau*k2]]
                    B=[[self.Re, 0,0],
                       [0,self.Pe_T, 0],
                       [0,0,self.Pe_S]];
                    #Compute eigenvalue and eigenvector of 
                    
                    #Update 2021/10/05, use the scipy package of the linalg. This can also solve the generalized eigenvalue problem
                    eig_val,eig_vec=linalg.eig(A,B) #use linear algebra package to compute eigenvalue
                    eig_val[np.isinf(eig_val)]=-np.inf
                    eig_val_max_ind=np.argmax(np.real(eig_val)) #compute the index of the eigenvalue
                    eig_vec_max=eig_vec[:,eig_val_max_ind] #get the corresponding eigen vector
                    
                    self.lambda_elevator=eig_val[eig_val_max_ind]
                    w0 =w0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[0]) #set the results weighted by the corresponding eigenvector 
                    T0 =T0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[1])
                    S0 =S0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[2])
                    #print(w0)
                    ##This is sample code to setup the initial condition as whatever I want...
                    #S0=np.array([1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16])
                    #S['g']=S0[slices]
                
                #superpose the field of secondary mode that has variation in vertical direction
                #Right now just superpose a simple one...
                
                #u0 =u0 + self.A_secondary*np.real(np.exp(1j*self.k_secondary*z))
                #w0 =w0 + self.A_secondary*np.real(np.exp(1j*self.k_secondary*z)) #set the results weighted by the corresponding eigenvector 
                T0 =T0 + self.A_secondary_T*np.real(np.exp(1j*self.k_secondary*z))
                S0 =S0 + self.A_secondary_S*np.real(np.exp(1j*self.k_secondary*z))
                    
                
                u['g']=u0
                w['g']=w0
                T['g']=T0
                S['g']=S0
                p['g']=p0
            
            elif self.flow in ['porous_media_2D']:
                x = domain.grid(0)
                z = domain.grid(1)
                u = solver.state['u']
                w = solver.state['w']
                p = solver.state['p']
                T = solver.state['T']
                
                gshape = domain.dist.grid_layout.global_shape(scales=1)
                slices = domain.dist.grid_layout.slices(scales=1)
                rand = np.random.RandomState(seed=23)
                noise = rand.standard_normal(gshape)[slices]
                
                ##Add the random noise
                u0=self.A_noise*noise
                w0=self.A_noise*noise
                T0=self.A_noise*noise
                p0=self.A_noise*noise
                
                k2=self.k_elevator**2
                    
                #The matrix to compute the eigenvalue of the porous media convection...
                A=[[-1, self.Ra_T],
                    [1, -k2] ]
                B=[[0, 0],
                   [0, 1]];
                
                eig_val,eig_vec=linalg.eig(A,B) #use linear algebra package to compute eigenvalue
                eig_val[np.isinf(eig_val)]=-np.inf
                eig_val_max_ind=np.argmax(np.real(eig_val)) #compute the index of the eigenvalue
                eig_vec_max=eig_vec[:,eig_val_max_ind] #get the corresponding eigen vector
                
                self.lambda_elevator=eig_val[eig_val_max_ind]
                w0 =w0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[0]) #set the results weighted by the corresponding eigenvector 
                T0 =T0 + self.A_elevator*np.real(np.exp(1j*self.k_elevator*x))*np.real(eig_vec_max[1])
                
                T0 =T0 + self.A_secondary_T*np.real(np.exp(1j*self.k_secondary*z))

                u['g']=u0
                w['g']=w0
                T['g']=T0
                p['g']=p0
        else:
            #Restart
            print('restart')
            write, last_dt = solver.load_state('restart.h5', -1)

        
        
    def run(self,solver,cfl,domain,logger):
        ##This CFL condition need to be modified for different simulation configuration.
        ##Note this line is very important and it needs to be added!!!!!
        if self.flow in ['IFSC_2D','double_diffusive_2D','double_diffusive_shear_2D','porous_media_2D']:
            cfl.add_velocities(('u','w'))

        logger.info('Starting loop')
        solver.stop_wall_time = np.inf
        solver.stop_iteration = np.inf
        start_time=time.time()
        while solver.ok:
            dt = cfl.compute_dt()    
            solver.step(dt)
            if solver.iteration % 100 == 0:
                logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
        
        end_time = time.time()
        # Print statistics
        logger.info('Run time: %f' %(end_time-start_time))
        logger.info('Iterations: %i' %solver.iteration)
        logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))

        
    def post_store(self,solver):
        #This post-processing variable need to be modified for different flow configuration
        if self.flow in ['IFSC_2D','double_diffusive_2D','double_diffusive_shear_2D']:
            analysis = solver.evaluator.add_file_handler('analysis',sim_dt=self.post_store_dt)
            analysis.add_task('S',layout='g',name='S')
            analysis.add_task('T',layout='g',name='T')
            analysis.add_task('u',layout='g',name='u')
            analysis.add_task('w',layout='g',name='w')
            analysis.add_task('p',layout='g',name='p')
            
            analysis.add_task("S",layout='c',name='S_coeff')
            analysis.add_task("T",layout='c',name='T_coeff')
            analysis.add_task("u",layout='c',name='u_coeff')
            analysis.add_task("w",layout='c',name='w_coeff')
            analysis.add_task("p",layout='c',name='p_coeff')
            
        elif self.flow in ['porous_media_2D']:
            analysis = solver.evaluator.add_file_handler('analysis',sim_dt=self.post_store_dt)
            analysis.add_task('T',layout='g',name='T')
            analysis.add_task('u',layout='g',name='u')
            analysis.add_task('w',layout='g',name='w')
            analysis.add_task('p',layout='g',name='p')
            
            analysis.add_task("T",layout='c',name='T_coeff')
            analysis.add_task("u",layout='c',name='u_coeff')
            analysis.add_task("w",layout='c',name='w_coeff')
            analysis.add_task("p",layout='c',name='p_coeff')
