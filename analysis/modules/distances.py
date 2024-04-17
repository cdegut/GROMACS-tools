import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from .plot import plt_smooth, ss_subplot, nicey_ax
import math
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from typing import Tuple

#MARK: RG 
#################

def Rg(atomistic_system, solvent_value = True):
    fig, ax = plt.subplots()
    fig.set(figwidth=12)
    Rgyr = []
    protein = atomistic_system.select_atoms("protein")
    for ts in atomistic_system.trajectory:
        Rgyr.append((atomistic_system.trajectory.time, protein.radius_of_gyration()))

    radius = np.array(Rgyr)
    radius = radius.T

    ax.plot(radius[0]/1000, radius[1], label='Rg', alpha = 0.4)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel(r'Rg  ($\AA$)')


    # plt median and + 1 stdev
    if len(radius[0]) > 1000:
        factor = 1000
    else:
        factor = 10

    plt_smooth(ax,radius[1], radius[0],factor)

    if solvent_value:
        theta_Rg, expanded_Rg, collapsed_Rg = random_walk_Rgs(len(atomistic_system.residues.resids))
        ax.axhline(y = expanded_Rg, color = 'g', linestyle = '-', alpha = 0.2, label=f"good solvent {expanded_Rg:.3}$\\AA$")
        ax.axhline(y = theta_Rg, color = 'g', linestyle = '-', alpha = 0.5, label=f"theta solvent {theta_Rg:.3}$\\AA$")
        ax.axhline(y = collapsed_Rg, color = 'g', linestyle = '-', alpha = 0.2, label=f"bad solvent {collapsed_Rg:.3}$\\AA$")

    rg_average = np.mean(radius[1,-1000:])
    ax.text(x=1.03, y=0.55, ha="left", va="top", s=f"Rg last ns: \n{rg_average:.3}$\\AA$", transform = ax.transAxes) 
    #plt_median(ax, radius[1], label=True)

    ax.set_xlim(0, len(radius.T)/100 )
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, 10))

    start, end = ax.get_ylim()
    yticks = np.arange(int(start), int(end), 5)
    ax.yaxis.set_ticks(yticks)

    ax.grid(visible = True, linestyle = '--', alpha=0.4)

    plt.legend(loc=(1.01, 0.6))
    plt.show()

    return Rgyr

def random_walk_Rgs(n_resids, Khun_lenght=False):
    
    typical_amino_acid = 3.6 # Angstrom
    #Khun_lenght = 8.8 # Angstrom

    theta_Rg =  typical_amino_acid * math.sqrt(n_resids/6)
    expanded_Rg = (typical_amino_acid * n_resids**(3/5) )/ math.sqrt(6)
    collapsed_Rg = (typical_amino_acid * n_resids**(1/3) )/ math.sqrt(6)

    if Khun_lenght:
        theta_Rg = math.sqrt((Khun_lenght / typical_amino_acid)) * theta_Rg 
        expanded_Rg = ((Khun_lenght / typical_amino_acid)**(1-3/5)) * expanded_Rg
        collapsed_Rg = ((Khun_lenght / typical_amino_acid)**(1-1/3)) * collapsed_Rg

    return theta_Rg, expanded_Rg, collapsed_Rg

def get_rolling_Rgs(atomistic_system, start, stop, residues_window):
    rolling_Rgs = np.zeros((stop - start, len(atomistic_system.residues.resids) - residues_window - 1))

    for k, ts in enumerate(range(start, stop)):
        atomistic_system.trajectory[ts]

        for i in range(1, len(atomistic_system.residues.resids) - residues_window):
            start = i
            finish = i + residues_window
            sele =  atomistic_system.select_atoms(f"resid {start} to {finish}")
            rolling_Rgs[k][i-1]= sele.radius_of_gyration()
        

    center_resids = np.arange(residues_window/2, len(atomistic_system.residues.resids) - residues_window/2 -1, 1)
    rolling_Rgs_avg = np.mean(rolling_Rgs, axis = 0)

    return np.column_stack((center_resids, rolling_Rgs_avg)).T

###############################
### 3D Distances plots #######
###############################

def calculate_3D_distance(atomistic_system):
    from MDAnalysis.analysis import distances

    n_res = len(atomistic_system.residues.resids)
    frames = len(atomistic_system.trajectory)
    distances_3Darray = np.zeros((frames,n_res,n_res))

    for ts in atomistic_system.trajectory:
        res_com = atomistic_system.atoms.center_of_mass(compound='residues')
        self_distances = distances.self_distance_array(res_com)
        sq_dist_arr = np.zeros((n_res, n_res))
        triu = np.triu_indices_from(sq_dist_arr, k=1)
        sq_dist_arr[triu] = self_distances
        sq_dist_arr.T[triu] = self_distances
        distances_3Darray[ts.frame, :, :] = sq_dist_arr
    
    return distances_3Darray

def plot_distances_slices(distances_3Darray, contact_start, contact_finish, dist_max):
    import matplotlib as mpl

    fig, axes = plt.subplots(1,3)
    fig.set(figwidth=15)
    fig.suptitle('Distance between centers-of-mass of residues', fontsize=18)

    initial_10frames = np.mean(distances_3Darray[contact_start *100 : contact_start *100 + 10], axis=0)
    final_10frames = np.mean(distances_3Darray[contact_finish *100 -10 :contact_finish *100], axis=0)

    im0 = axes[0].pcolormesh(initial_10frames,vmin=1, vmax=dist_max, cmap = plt.cm.inferno_r)
    im1 = axes[1].pcolormesh(final_10frames,vmin=1, vmax=dist_max, cmap = plt.cm.inferno_r)

    distances_average = np.mean(distances_3Darray[contact_start *100 : contact_finish *100], axis=0)                                                 
    im2 = axes[2].pcolormesh(distances_average,vmin=1, vmax=dist_max, cmap = plt.cm.inferno_r)
    # plt.pcolor gives a rectangular grid by default
    # so we need to make our heatmap square
    for ax in axes:
        ax.set_aspect('equal')
        ax.set_ylabel('Residue IDs')
        ax.set_xlabel('Residue IDs')
        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(start, end, 10))
        start, end = ax.get_ylim()
        ax.yaxis.set_ticks(np.arange(start, end, 10))
        ax.grid(visible = True, linestyle = '--', alpha=0.4)
        
    axes[0].set_title('Initial 10 frames (0.1ns)')
    axes[1].set_title('Final 10  frames (0.1ns)')
    axes[2].set_title(f'Average over {contact_finish - contact_start}ns')

    # colorbar
    cax,kw = mpl.colorbar.make_axes([ax for ax in axes.flat])
    plt.colorbar(im2, cax=cax, **kw)


#MARK: plot_distances


def multi_plot_distances(Rgyr,  distances_3Darray, contact_start, contact_finish, dist_max, dsspline_start, Rg_start = 0,
                      rolling_Rgs=False, Khun_lengh=False, residues_window=10,
                      atomistic_system=None) -> tuple[plt.Figure, tuple[plt.Axes, plt.Axes, plt.Axes | None, plt.Axes]]:

    n_resids =  len(atomistic_system.residues.resids)
    #Replot the rgyr with AOI
    radius = np.array(Rgyr)
    if Rg_start != 0:
        radius = radius[Rg_start*100:]
    radius = radius.T
    theta_Rg, expanded_Rg, collapsed_Rg = random_walk_Rgs(n_resids)
    
    fig= plt.figure(constrained_layout=True)
    if rolling_Rgs.any():
        spec = gridspec.GridSpec(4,1,height_ratios=[1,4,0.8,0.5] ,figure=fig)
    else:
        spec = gridspec.GridSpec(3,1,height_ratios=[1,4,0.5], figure=fig)

#### plot 1 Rg
    fig.set(figwidth=10, figheight=12)
    ax0 = fig.add_subplot(spec[0,0])

    ax0.set_xlabel('Time (ns)')
    ax0.set_ylabel(r'Rg  ($\AA$)')


    if len(radius[0]) > 1000:
        factor = 1000
    else:
        factor = 10
    plt_smooth(ax0,radius[1], radius[0],factor)

    ax0.axhline(y = expanded_Rg, color = 'g', linestyle = '-', alpha = 0.2, label=f"good solvent {expanded_Rg:.3}$\\AA$")
    ax0.axhline(y = theta_Rg, color = 'g', linestyle = '-', alpha = 0.4, label=f"theta solvent {theta_Rg:.3}$\\AA$")
    ax0.axhline(y = collapsed_Rg, color = 'g', linestyle = '-', alpha = 0.2, label=f"bad solvent {collapsed_Rg:.3}$\\AA$")
    ax0.axvspan(contact_start, contact_finish, color='red', alpha=0.3, label="Selected timeframe")

    ax0.set_xlim(Rg_start, len(radius.T)/100 +  Rg_start)
    start, end = ax0.get_xlim()
    ax0.xaxis.set_ticks(np.arange(start, end, 10))

    start, end = ax0.get_ylim()
    yticks = np.arange(int(start), int(end), 10)
    ax0.yaxis.set_ticks(yticks)
    ax0.grid(visible = True, linestyle = '--', alpha=0.4)

    rg_average = np.mean(radius[1,contact_start*100:contact_finish*100])
    ax0.text(x=1.03, y=0.1, ha="left", va="top", s=f"Rg (red): {rg_average:.3}($\\AA$)", transform = ax0.transAxes) 
    ax0.legend(loc=(1.01, 0.15))

####### plot colormesh

    ax1 = fig.add_subplot(spec[1,0])
    distances_average = np.mean(distances_3Darray[contact_start *100 : contact_finish *100], axis=0)  
    im2 = ax1.pcolormesh(distances_average,vmin=1, vmax=dist_max, cmap = plt.cm.inferno_r) 
    start, end = ax1.get_xlim()
    ax1.xaxis.set_ticks(np.arange(start, end, 5))
    start, end = ax1.get_ylim()
    ax1.yaxis.set_ticks(np.arange(start, end, 5))
    ax1.grid(visible = True, linestyle = '--', alpha=0.4)

    axins = inset_axes(ax1,width="5%", height="50%",  loc='center left', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax1.transAxes, borderpad=0) 
    fig.colorbar(im2, cax = axins, ax=ax1, label="CoM to CoM distance  ($\\AA$)") 


    ######## plotrolling Rg
    if rolling_Rgs.any():
        ax2 = fig.add_subplot(spec[2,0], sharex=ax1)
        ax2.set_ylabel(r'Rg  ($\AA$)')
        ax2.set_xlabel(f"Center of {residues_window} residues window")

        ax2.xaxis.set_ticks(np.arange(0,n_resids , 5))
        ax2.grid(visible = True, linestyle = '--', alpha=0.4)
        ax2.plot(rolling_Rgs[0], rolling_Rgs[1])

        theta_Rg, expanded_Rg, collapsed_Rg = random_walk_Rgs(residues_window, Khun_lengh)
        #ax2.scatter([],[], alpha=0,  label=f"With lK ={Khun_lengh}$\\AA$:") ## fake scatter to add a legend handle
        ax2.axhline(y = expanded_Rg, color = 'g', linestyle = '-', alpha = 0.2, label=f"good solvent {expanded_Rg:.3}$\\AA$")
        ax2.axhline(y = theta_Rg, color = 'g', linestyle = '-', alpha = 0.5, label=f"theta solvent {theta_Rg:.3}$\\AA$")
        #ax2.axhline(y = collapsed_Rg, color = 'g', linestyle = '-', alpha = 0.2, label="bad solvent") # not applicable 
        ax2.legend(loc=(1.01, 0.1), title= f"With lK ={Khun_lengh}$\\AA$:")
    
    ss = fig.add_subplot(spec[-1,0], sharex=ax1)
    ss_subplot(fig, ss, dsspline_start)
    ss.legend(loc=(1.01, 0.2))

    
    ax0.set_title("Protein Rg during simulation", loc='left')
    ax1.set_title("Residues to residues distance matrix", loc='left')
    ax2.set_title("Rolling Rgs", loc='left')
    ss.set_title("Secondary structure", loc='left')

    return fig , (ax0, ax1, ax2, ss)

def matrix_zoom(ax: plt.Axes, roi: Tuple[Tuple], dist_max: int, fig, distances_3Darray, contact_start, contact_finish) -> plt.Axes:
    
    distances_average = np.mean(distances_3Darray[contact_start *100 : contact_finish *100], axis=0)  
    im = ax.pcolormesh(distances_average,vmin=1, vmax=dist_max, cmap = plt.cm.inferno_r) 

    ax.set_xlim(roi[0])
    ax.set_ylim(roi[1])
    ax.xaxis.set_ticks(np.arange(roi[0][0], roi[0][1], 5))
    ax.yaxis.set_ticks(np.arange(roi[1][0], roi[1][1], 5))
    ax.grid(visible = True, linestyle = '--', alpha=0.4)
    #axins = inset_axes(ax,width="5%", height="50%",  loc='center left', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax.transAxes, borderpad=0) 
    fig.colorbar(im, ax=ax, label="CoM to CoM distance  ($\\AA$)")

    return ax

def plot_every_diagonal(matrix, contact_start,contact_finish, cutoff):
    distances_average = np.mean(matrix[contact_start *100 : contact_finish *100], axis=0)
    diagonals = []

    for i in range(0,len(distances_average[0])):
        diagonals.append([distances_average[i,i]])
        for j in range(1,i+1):
            k = i-j
            l = i+j
            if l <= len(distances_average[0]) -1:
                diagonals[i].append(distances_average[k,l])

    conatac_source = []
    for diag in diagonals:
        sum = 0
        for element in diag:
            if element < cutoff/4:
                sum += 4
            elif element < cutoff/2:
                sum += 2
            elif element < cutoff:
                sum += 1
        conatac_source.append(sum)

    fig = plt.figure()
    ax = fig.add_subplot()
    fig.set(figwidth=12)
    ax.plot( np.arange(1, len(conatac_source)+1 ), conatac_source, )

    ax.xaxis.set_ticks(np.arange(0, len(conatac_source), 5))
    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(np.arange(start, end, 2))
    ax.grid(visible = True, linestyle = '--', alpha=0.4)

def sele_distance(atomistic_system, seleA, seleB, random_walk_step, Khun_lengh):

    distances = []
    for ts in atomistic_system.trajectory:
        CoM_A = atomistic_system.select_atoms(seleA).center_of_mass()
        CoM_B = atomistic_system.select_atoms(seleB).center_of_mass()
        dist = math.sqrt((CoM_B[0] - CoM_A[0])**2 + (CoM_B[1] - CoM_A[1])**2 + (CoM_B[2] - CoM_A[2])**2)
        distances.append([dist, atomistic_system.trajectory.time])

    distances = np.array(distances).T
    fig = plt.figure()
    ax = fig.add_subplot()
    fig.set(figwidth=10)
    ax.plot(distances[1]/1000, distances[0], alpha=0.4)
    nicey_ax(ax)
    plt_smooth(ax,distances[0],distances[1],1000)

    theta_Rg, expanded_Rg, collapsed_Rg = random_walk_Rgs(random_walk_step, Khun_lengh)
    ax.axhline(y = expanded_Rg, color = 'g', linestyle = '-', alpha = 0.2, label="good Solvent")
    ax.axhline(y = theta_Rg, color = 'g', linestyle = '-', alpha = 0.5, label="theta Solvent")
    #ax.axhline(y = collapsed_Rg, color = 'g', linestyle = '-', alpha = 0.2, label="bad Solvent")
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel(r'Rg  ($\AA$)')

    #plt_median(ax,distances[0], label=True)