import numpy as np
from mpi4py import MPI
import pathlib
#from IPython import display
from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.tools import post
import logging

##These lib are used for post-processing
import h5py
import pathlib
import subprocess
import shutil
import dedalus_setup
shutil.rmtree('analysis',ignore_errors=True)

logger = logging.getLogger(__name__)
flag=dedalus_setup.flag()


flag.flow='porous_media_2D'

#-----------------setup storing for post-processing
flag.post_store_dt=0.5;
flag.stop_sim_time=200;

#------------ print these parameters in the screen
flag.print_screen(logger)

#---------main loop to run the dedalus 
domain=flag.build_domain()
problem=flag.governing_equation(domain)
ts = de.timesteppers.RK443
solver =  problem.build_solver(ts)
flag.initial_condition(domain,solver)
solver.stop_sim_time = flag.stop_sim_time
flag.post_store(solver)
flag.print_file() #move print file to here.
flag.run(solver,cfl,domain,logger)
