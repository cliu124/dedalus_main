#
# Real Ginzburg-Landau Equation on a Time-Dependent Domain using Dedalus
#

import os
import argparse

import numpy as np
from scipy import optimize
from scipy import integrate
from scipy import linalg
import matplotlib.pyplot as plt
from matplotlib import ticker

import dedalus.public as d3
import h5py

# Extra packages not in default dedalus3 install: moviepy, pandas (dependencies: xlrd, openpyxl)
import moviepy.editor as mpy
from moviepy.video.io.bindings import mplfig_to_npimage
import pandas as pd



### CONFIG ###

# File paths
regimes_sheet_name = "RGLE_parameter_regimes.xlsx"
sims_folder = "sims"
analysis_folder = "analysis"
if not os.path.isdir(sims_folder):
    os.mkdir(sims_folder)
if not os.path.isdir(analysis_folder):
    os.mkdir(analysis_folder)

# Set basic matplotlib settings
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = '12'



### PARSE COMMAND LINE ARGUMENTS ###

parser = argparse.ArgumentParser()
parser.add_argument('id', type=int, help='id of parameter regime to run')
parser.add_argument('-s', '--simulate', default=False, action='store_true', help="run Dedalus simulation and initial meshes")
parser.add_argument('-a', '--analysis', nargs='*', type=int, default=[], help="Create basic analysis plots. 0 for interactive. Input numbers for specific plots.")
parser.add_argument('-m', '--movieplayback', type=np.float64, default=0, help="playback speed of movies, set to 0 if no movie")

args = parser.parse_args()

id = args.id
simulate = args.simulate
analysis = args.analysis
movieplayback = args.movieplayback



### LOAD PARAMETER REGIME DATAFRAME FROM EXCEL FILE ###

regimes_df = pd.read_excel(regimes_sheet_name)

n = int(regimes_df['n'][id])
t_total = regimes_df['t_total'][id]
t_step = regimes_df['t_step'][id]
L_type = regimes_df['L'][id]
a = regimes_df['a'][id]
b = regimes_df['b'][id]
if isinstance(b, str) and 'pi' in b:
    # of the form c*pi
    b = float(b[:-3]) * np.pi
Q = int(regimes_df['Q'][id])
mu = regimes_df['mu'][id]
R0 = regimes_df['R0'][id]
k = int(regimes_df['k'][id])
ak_plus0 = regimes_df['akpR'][id] + 1j*regimes_df['akpI'][id]
ak_minus0 = regimes_df['akmR'][id] + 1j*regimes_df['akmI'][id]
dilution = bool(regimes_df['dilution'][id])
try:
    times = [np.float64(s) for s in regimes_df['times'][id].split(",")]
except:
    times = []




### PARAMETER SETUP ###


def get_L_for_type(L_type):
    if L_type == "exp":
        L = lambda t: np.exp(a*t)
        L_dot = lambda t: a*np.exp(a*t)
        L_str = r'e^{' + rf'{a}' r't}'
    elif L_type == "sin":
        L = lambda t: 1 + a*np.sin(b*t)
        L_dot = lambda t: a*b*np.cos(b*t)
        L_str = rf'1 + {a}\sin({b}t)'
        L_str = r'1 + \frac{1}{2}\sin(\frac{\pi}{10}t)'
    elif L_type == "expdec":
        L = lambda t: np.exp(-a*t)
        L_dot = lambda t: -a*np.exp(-a*t)
        L_str = r'e^{-' + rf'{a}' r't}'
    elif L_type == "rational":
        L = lambda t: a - (t + 1)**b
        L_dot = lambda t: -b*(t + 1)**(b-1)
        L_str = rf'{a}' + r'- (t+1)^{' + rf'{b}' + r'}'
    elif L_type == "sqrtsin":
        L = lambda t: np.sqrt(1+(1/25)*np.sin(2*np.pi*t))
        L_dot = lambda t: (np.pi*np.cos(2*np.pi*t)) / (25*np.sqrt((1+(1/25)*np.sin(2*np.pi*t))))
        L_str = r'\sqrt{1+\frac{1}{25}\sin(2\pi t)}'
    elif L_type == "const": 
        L = lambda t: 1
        L_dot = lambda t: 0
        L_str = r'1'
    else:
        raise ValueError(f'{L_type} not a valid L(t)')
    return L, L_dot, L_str


L, L_dot, L_str = get_L_for_type(L_type)
t_output = (1 / t_step) * 1e-2 # how many time steps before outputting a data point: allow for 1e-2 time coarsening

initial_condition = lambda xi: R0 * np.exp(1j*Q*xi) + ak_plus0 * np.exp(1j*(Q+k)*xi) + np.conjugate(ak_minus0) * np.exp(1j*(Q-k)*xi)

a_desc = r"a_{k+} = " + f"{ak_plus0}" + r", a_{k-} = " + f"{ak_minus0}"
description = fr'$Q={Q}$, $\mu={mu}$, $k = {k}$, ${a_desc}$, $L(t) = {L_str}$'

name = f"RGLE_{str(id).zfill(3)}"
sim_base_name = f'sim_{name}'
analysis_path = analysis_folder + '/' + name






### SIMULATION ###

def run_simulation():
    if os.path.isdir(sims_folder + '/' + sim_base_name):
        yn = input(f"Data for parameter regime {id} already exists. Do you wish to proceed (y/n)? ")
        if yn != "y":
            print("Simulation aborted.")
            return

    # Basis
    xicoord = d3.Coordinate('xi')
    dist = d3.Distributor(xicoord, dtype=np.complex128)
    xibasis = d3.ComplexFourier(xicoord, n, bounds=(0, 2*np.pi), dealias=2)

    # Fields
    A = dist.Field(name='A', bases=xibasis)

    # Substitutions
    dxi = lambda xi: d3.Differentiate(xi, xicoord)
    magsq_A = A * np.conj(A)

    # Problem
    problem = d3.IVP([A], namespace = globals() | locals())

    # fix Dedalus 3 bug: see namespace_additions and namespace in problems.py
    additions = {}
    additions[problem.time] = problem.sim_time_field
    problem.namespace.update(additions)

    if L_type == "const":
        problem.add_equation("dt(A) - mu*A - dxi(dxi(A)) = -magsq_A * A")
    elif dilution:
        problem.add_equation("dt(A) - mu*A = (L_dot(t)/L(t))*dxi(A) + (1/(L(t)**2))*dxi(dxi(A)) - magsq_A * A")
    else:
        problem.add_equation("dt(A) - mu*A = (1/(L(t)**2))*dxi(dxi(A)) - magsq_A * A")

    # Solver
    solver = problem.build_solver(d3.RK222)
    solver.stop_sim_time = t_total

    # Initial conditions
    xi = dist.local_grid(xibasis)
    A['g'] = initial_condition(xi)

    # Setup analysis
    A.change_scales(1)
    analysis = solver.evaluator.add_file_handler(sims_folder + '/' + sim_base_name, iter=t_output)
    analysis.add_tasks(solver.state, layout='g')
    analysis.add_task(A, layout='c', name='A_coeff')

    # Main loop
    print(f'Total iterations: {int(round(t_total / t_step))}')
    while solver.proceed:
        solver.step(t_step)
        if solver.iteration % 10000 == 0:
            print('Completed iteration {}'.format(solver.iteration))




### ANALYSIS ###


# Convenience method to save plot to file
def save_fig(suf):
    if not os.path.isdir(analysis_path):
        os.mkdir(analysis_path)
        print(f"Created folder: {analysis_path}")
    file_name = f"{analysis_path}/{name}_{suf}.png"
    plt.savefig(fname=file_name)
    print(f"Saved to file: {file_name}")

# Convenience method to convert real time to index
def time_to_ind(time):
    return int(time / (t_step * t_output))



# Plot A over t vs. x
def plot_meshes(xi, t, A, A_coeff):
    plt.figure(figsize=(14, 11), dpi=150)

    plt.subplot(2,2,1)

    plt.pcolormesh(xi, t, np.real(A), shading='nearest')
    plt.colorbar()
    plt.xlabel(r'$\xi$')
    plt.ylabel('t')
    plt.title(description + fr': Re $A$')
    plt.tight_layout()

    plt.subplot(2,2,2)

    plt.pcolormesh(xi, t, np.imag(A), shading='nearest')
    plt.colorbar()
    plt.xlabel(r'$\xi$')
    plt.ylabel('t')
    plt.title(description + fr': Im $A$')
    plt.tight_layout()

    R = np.abs(A)

    plt.subplot(2,2,3)

    plt.pcolormesh(xi, t, R, shading='nearest')
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('t')
    plt.title(description + fr': R')
    plt.tight_layout()

    plt.subplot(2,2,4)

    freq_bound = 4*Q
    freq = np.arange(0, freq_bound, 1)

    plt.pcolormesh(freq, t, np.abs(A_coeff[:,0:freq_bound]), shading='nearest')
    plt.colorbar()
    plt.xlabel('wavenumber')
    plt.ylabel('t')
    plt.title(description + fr': Fourier')
    plt.tight_layout()

    save_fig("mesh")


def make_snapshots(xi, t, A):
    """Makes snapshots for each time in times"""
    for time in times:
        ind = time_to_ind(time)
        if ind >= t.size:
            ind = -1 + t.size
        plt.figure(figsize=(16, 12), dpi=150)

        plt.subplot(2,2,1)

        plt.plot(xi, np.real(A[ind,:]))
        plt.xlabel(r'$\xi$')
        plt.ylabel(r'Re $A$')

        plt.subplot(2,2,2)

        plt.plot(xi, np.imag(A[ind,:]))
        plt.xlabel(r'$\xi$')
        plt.ylabel(r'Im $A$')
        
        R = np.abs(A)

        plt.subplot(2,2,3)

        plt.plot(xi, R[ind,:])
        plt.xlabel(r'$\xi$')
        plt.ylabel(r'$|A|$')

        plt.subplot(2,2,4)

        freq_bound = 4*Q
        freq = np.arange(0, freq_bound, 1)

        plt.stem(freq, np.abs(A_coeff[ind,0:freq_bound]), linefmt='C0', markerfmt=" ", basefmt="C0-")
        plt.xlabel(r'Wavenumber')
        plt.ylabel(r'Fourier Coefficient, Complex Norm')

        plt.suptitle(description + fr": $t = {time}$")

        save_fig(f"{str(ind).zfill(5)}_snapshot") 
    
    print("Snapshots created.")


def plot_dilution(xi, t, A):
    t = np.array(t)
    A = np.array(A)

    ind = 0 # np.random.randint(0, n)
    A_theory = np.sqrt(mu - (Q**2)/(L(t)**2))*np.exp(1j*Q*(xi[ind] + np.log(L(t))))

    plt.figure(figsize=(8,8), dpi=150)

    plt.plot(t, np.real(A)[:,ind], label=r"Re $A(x=0)$, numerical", color="royalblue")
    plt.plot(t, np.real(A_theory), label=r"Re $A(x=0)$, theory", color="lightsteelblue", linestyle="dashed")
    plt.plot(t, np.imag(A)[:,ind], label=r"Im $A(x=0)$, numerical", color="orange")
    plt.plot(t, np.imag(A_theory), label=r"Im $A(x=0)$, theory", color="navajowhite", linestyle="dashed")
    plt.xlabel(r"$t$")
    plt.title(description + fr': Dilution')
    plt.legend()

    save_fig("dilution")


def plot_ptb(xi, t, A_coeff):
    # ONLY WORKS FOR TIME-INDEPENDENT!!
    amp2 = mu - Q**2
    disc = np.sqrt((2*Q*k)**2 + amp2**2)
    lambda_plus = -amp2 - k**2 + disc
    lambda_minus = -amp2 - k**2 - disc

    v11 = (2*Q*k - disc) / amp2 # vk+
    v12 = 1
    v21 = (2*Q*k + disc) / amp2 # vk-
    v22 = 1

    c1 = (ak_plus0 - (v21/v22)*ak_minus0) / (v11 - (v12*v21/v22))
    c2 = (ak_minus0 - c1*v12) / v22

    ak_plus = c1*v11*np.exp(lambda_plus*t) + c2*v21*np.exp(lambda_minus*t)
    ak_minus = c1*v12*np.exp(lambda_plus*t) + c2*v22*np.exp(lambda_minus*t)

    fig, ax = plt.subplots(1, 1, figsize=(10,8), dpi=150)

    ax.plot(t, np.abs(A_coeff[:,Q+k]), color="orange", label=r"numerical $|a_{k+}|$")
    ax.plot(t, np.abs(A_coeff[:,Q-k]), color="royalblue", label=r"numerical $|a_{k-}|$")
    ax.plot(t, np.abs(ak_plus), linestyle="dashed", color="navajowhite", label=r"asymptotic $|a_{k+}|$")
    ax.plot(t, np.abs(ak_minus), linestyle="dashed", color="lightsteelblue", label=r"asymptotic $|a_{k-}|$")
    ax.hlines(y=R0, xmin=0, xmax=np.max(t), linestyle="dashed", color="gray", label="base state")
    ax.set_xlabel(r'$t$')
    ax.set_title(description + fr': ptb')
    ax.set_ylim(0, 1.1*R0)
    ax.legend(loc="lower right")

    axin = ax.inset_axes([0.05, 0.55, 0.4, 0.4])
    axin.plot(t, np.abs(A_coeff[:,Q+k]), color="orange")
    axin.plot(t, np.abs(A_coeff[:,Q-k]), color="royalblue")
    axin.plot(t, np.abs(ak_plus), linestyle="dashed", color="navajowhite")
    axin.plot(t, np.abs(ak_minus), linestyle="dashed", color="lightsteelblue")
    if len(times) != 0:
        t_max = times[0] * 2 / 3
    else:
        t_max = t_total / 4
    axin.set_xlim(0, t_max)
    t_max_ind = time_to_ind(t_max)
    axin.set_ylim(0, 1.1*np.max([np.abs(ak_plus)[:t_max_ind], np.abs(ak_minus)[:t_max_ind]]))
    axin.tick_params(axis='x', labelsize=5)
    axin.tick_params(axis='y', labelsize=5)

    ax.indicate_inset_zoom(axin, edgecolor="black")
    
    save_fig("ptb")


def plot_relaxation(xi, t, A, A_coeff):
    amp_diff = np.max(np.abs(A), axis=1) - np.min(np.abs(A), axis=1)

    plt.figure(figsize=(12, 8), dpi=150)
    plt.suptitle(description + fr': relaxation')

    plt.subplot(1,2,1)

    plt.plot(t, amp_diff)
    plt.xlabel(r'$t$')
    plt.ylabel(r'max - min amplitude')
    plt.xlim(left=times[0])

    plt.subplot(1,2,2)

    def model_lin(x, m, c):
        return m*x + c

    log_Q_coeff = np.log(np.abs(A_coeff[:,Q]))
    
    plt.plot(t, log_Q_coeff, label="numerical") # only works for simple phase slips
    plt.xlabel(r'$t$')
    plt.ylabel(r"numerical $\ln |a_Q|$")
    plt.xlim(left=times[0])

    fit_start_time = times[-2] if len(times) >= 3 else times[-1]
    fit_start_ind = time_to_ind(fit_start_time)
    t_trimmed = t[fit_start_ind:]

    lin_opt, lin_cov = optimize.curve_fit(model_lin, t_trimmed, log_Q_coeff[fit_start_ind:])
    m, c = lin_opt[0], lin_opt[1]
    predicted = model_lin(t_trimmed, m, c)

    plt.plot(t_trimmed, predicted, linestyle="dashed", label=f"linear fit, slope: {m:.3f}")
    plt.legend()
    
    save_fig("relax")

def plot_delay(xi, t, A_coeff):
    # Analogous to plot_ptb but for time-dependent case

    _, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,8), dpi=150)

    amp2 = mu - (Q**2 / L(t)**2)
    A = np.array([
        [-amp2 - ((2*Q*k + k**2) / L(t)**2), -amp2],
        [-amp2, -amp2 - ((-2*Q*k + k**2) / L(t)**2)]
    ]).transpose(2,0,1)
    energy_mat = A + A.conj().transpose(0,2,1)
    max_eig = []
    for i in range(0, t.size):
        max_eig.append(np.max(np.real(linalg.eigvals(energy_mat[i,:,:]))))

    energy_initial = np.abs(ak_plus0)**2 + np.abs(ak_minus0)**2
    energy_integrated = integrate.cumulative_trapezoid(np.array(max_eig), x=t, initial=0)

    ak_plus_abs = np.abs(A_coeff[:,Q+k])
    ak_minus_abs = np.abs(A_coeff[:,Q-k])
    energy_actual = ak_plus_abs**2 + ak_minus_abs**2

    energy_log = np.log(energy_actual / energy_initial)
    x_max = times[2] # can change depending on what I'm looking for
    x_max_ind = time_to_ind(x_max)
    y_min = -1.2*np.abs(np.min(energy_log))
    y_max = 1.2*np.abs(np.max(energy_log[:x_max_ind]))

    ax1.plot(t, energy_log, label=r"energy actual")
    ax1.plot(t, energy_integrated, linestyle="dashed", label=r"energy integrated")
    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"energy, log scale")
    ax1.legend()
    ax1.set_xlim(0, x_max)
    ax1.set_ylim(y_min, y_max)

    ax2.plot(t, energy_actual, label="energy actual")
    ax2.plot(t, np.exp(energy_integrated), linestyle="dashed", label="energy integrated")
    ax2.set_xlabel(r"$t$")
    ax2.set_ylabel(r"energy")
    ax2.legend()
    ax2.set_xlim(0,x_max)
    ax2.set_ylim(0,np.exp(y_max))

    plt.suptitle(description + fr": bifurcation delay")

    save_fig("delay")




### MOVIES ###

def create_movie(xi, t, A, speed=1, lagrangian=False, gif=False):
    if len(times) < 2:
        print("Must input at least two `times` to create movie.")
        return
    
    start_time = times[0]
    start_ind = time_to_ind(start_time)
    end_time = times[-1]

    x_max = 2*np.pi*np.max(L(t)) if not lagrangian else 2*np.pi
    y_max = np.max(np.abs(A)) + 2

    fig, (ax1, axempty, ax2) = plt.subplots(3, 1, figsize=(8,8), facecolor="white", gridspec_kw={'height_ratios': [5, 0.1, 1]})

    axempty.set_visible(False)

    if lagrangian:
        x = xi
    else:
        x = L(start_time)*xi
    A_cut = A[start_ind,:]
    ax1_real, = ax1.plot(x, np.real(A_cut), linewidth=3, color="red", solid_capstyle='round', label=fr"Re $A$")
    ax1_imag, = ax1.plot(x, np.imag(A_cut), linewidth=3, color="orange", solid_capstyle='round', label=fr"Im $A$")
    ax1_abs, = ax1.plot(x, np.abs(A_cut), linewidth=1, color="green", solid_capstyle='round', label=fr"$|A|$")
    
    if not lagrangian:
        vert_line_x = np.array([2*np.pi, 2*np.pi]) * L(start_time)
        vert_line_y = np.array([-y_max, y_max])
        ax1_vert, = ax1.plot(vert_line_x, vert_line_y, color="gray", linestyle="dashed")

    ax1_time = ax1.text(0.85*x_max, -0.9*y_max, rf"$t={(start_time):.2f}$", color="blue")

    ax1.legend(loc='upper right')

    ax1.set_xlim(0, x_max)
    ax1.set_ylim(-y_max, y_max)
    
    title = description + fr': $A$'
    if lagrangian:
        ax1.set_xlabel(r"$\xi$")
        title += " Lagrangian"
    else:
        ax1.set_xlabel(r"$x$")
        title += " Eulerian"
    ax1.set_title(title)

    ax2.yaxis.set_major_locator(ticker.NullLocator())
    ax2.spines.right.set_color('none')
    ax2.spines.left.set_color('none')
    ax2.xaxis.set_ticks_position('bottom')
    ax2.set_xlim(0.8*Q**2, 3.2 * Q**2)
    ax2.set_ylim(-1, 1)

    k = np.arange(1, 2*Q + 1, 1)
    eckhaus_pts = 3*Q**2 - 0.5*k**2
    ax2.vlines(eckhaus_pts, -1, 1, linestyles='dashed', colors='gray')

    stringify = np.vectorize(lambda x, k: f"{x} ({k})")
    labels = stringify(eckhaus_pts, k)

    ax2.tick_params(which='major', labelsize=9)
    ax2.xaxis.set_major_locator(ticker.FixedLocator(eckhaus_pts))
    ax2.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
    ax2.set_xlabel(r"Dashed lines: $3Q^2 - \frac{1}{2}k^2$ ($k$). Dot: $L^2(t)\mu$.")
    ax2.set_title("Eckhaus instability")

    ax2_point, = ax2.plot(L(start_time)**2 * mu, 0, 'o')

    def make_frame(t):
        real_t = start_time + t*speed
        if lagrangian:
            x = xi
        else:
            x = L(real_t)*xi
        A_cut = A[time_to_ind(real_t),:]

        ax1_real.set_data(x, np.real(A_cut))
        ax1_imag.set_data(x, np.imag(A_cut))
        ax1_abs.set_data(x, np.abs(A_cut))
        if not lagrangian:
            ax1_vert.set_xdata(np.array([2*np.pi, 2*np.pi]) * L(real_t))
        ax1_time.set_text(rf"$t={(real_t):.2f}$")

        ax2_point.set_xdata(L(real_t)**2 * mu)
        
        return mplfig_to_npimage(fig)

    file_name = f"{analysis_path}/{name}_movie.mp4" if not lagrangian else f"{analysis_path}/{name}_movie_lagrangian.mp4"
    animation = mpy.VideoClip(make_frame, duration=int((end_time - start_time)/speed))
    animation.write_videofile(filename=file_name, fps=24)



### RUN FUNCTIONS BASED ON SELECTED COMMAND LINE ARGUMENTS ###

if simulate:
    run_simulation()

with h5py.File(f'{sims_folder}/{sim_base_name}/{sim_base_name}_s1/{sim_base_name}_s1_p0.h5', mode='r') as file:
    A = file['tasks']['A']
    t = A.dims[0]['sim_time']
    xi = A.dims[1][0]

    xi, t, A = np.array(xi), np.array(t), np.array(A)

    A_coeff = np.array(file['tasks']['A_coeff'])
    if simulate:
        plot_meshes(xi, t, A, A_coeff)
    if len(analysis) > 0:
        if 0 in analysis:
            # Interactive mode
            yn = input("Create snapshots? ")
            if yn == "y":
                make_snapshots(xi, t, A)
            yn = input("Create ptb/delay plot? ")
            if yn == "y":
                if L_type == "const":
                    plot_ptb(xi, t, A_coeff)
                else:
                    plot_delay(xi, t, A_coeff)
            yn = input("Create relaxation plot? ")
            if yn == "y":
                plot_relaxation(xi, t, A, A_coeff)
            if L_type != "const":
                yn = input("Create dilution plot? ")
                if yn == "y":
                    plot_dilution(xi, t, A)
        else:
            if 1 in analysis:
                make_snapshots(xi, t, A)
            if 2 in analysis:
                if L_type == "const":
                    plot_ptb(xi, t, A_coeff)
                else:
                    plot_delay(xi, t, A_coeff)
            if 3 in analysis:
                plot_relaxation(xi, t, A, A_coeff)
            if 4 in analysis:
                plot_dilution(xi, t, A)
    if movieplayback != 0:
        create_movie(xi, t, A, speed=movieplayback, lagrangian=False)
        if L_type != "const":
            create_movie(xi, t, A, speed=movieplayback, lagrangian=True)
