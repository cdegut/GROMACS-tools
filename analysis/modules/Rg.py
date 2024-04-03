import numpy as np
import matplotlib.pyplot as plt
from .plot import plt_smooth, ss_subplot
import math

### RG ##########
#################

def Rg(atomistic_system):
    fig, ax = plt.subplots()
    fig.set(figwidth=12)
    Rgyr = []
    protein = atomistic_system.select_atoms("protein")
    for ts in atomistic_system.trajectory:
        Rgyr.append((atomistic_system.trajectory.time, protein.radius_of_gyration()))

    radius = np.array(Rgyr)
    radius = radius.T
    # Plot the data - RMSF against residue index

    ax.plot(radius[0]/1000, radius[1], label='Rg', alpha = 0.4)
    # Add axis labels
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Rg  ($\AA$)')

    # Set axis limits
    ax.set_xlim(0, len(radius.T)/100 )
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, 10))
    ax.grid(visible = True, linestyle = '--', alpha=0.4)

    # plt median and + 1 stdev
    if len(radius[0]) > 1000:
        factor = 1000
    else:
        factor = 10

    plt_smooth(ax,radius[1], radius[0],factor)

    theta_Rg, expanded_Rg, collapsed_Rg = random_walk_Rgs(len(atomistic_system.residues.resids))
    ax.axhline(y = expanded_Rg, color = 'g', linestyle = '-', alpha = 0.2, label="good Solvent")
    ax.axhline(y = theta_Rg, color = 'g', linestyle = '-', alpha = 0.5, label="theta Solvent")
    ax.axhline(y = collapsed_Rg, color = 'g', linestyle = '-', alpha = 0.2, label="bad Solvent")
    #plt_median(ax, radius[1], label=True)
    # Show legend
    plt.legend()
    # Show plot
    plt.show()

    return Rgyr

def random_walk_Rgs(n_resids):
    
    typical_step_lengh = 3.6 # Angstrom

    theta_Rg =  typical_step_lengh * math.sqrt(n_resids/6)
    expanded_Rg = (typical_step_lengh * n_resids**(3/5) )/ math.sqrt(6)
    collapsed_Rg = (typical_step_lengh * n_resids**(1/3) )/ math.sqrt(6)
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

def plot_distances_HD(Rgyr, distances_3Darray, contact_start, contact_finish, dist_max, dsspline_start, do_rolling_Rgs=False, atomistic_system=None,  residues_window=10):
    import matplotlib.gridspec as gridspec

    #Replot the rgyr with AOI
    radius = np.array(Rgyr)
    radius = radius.T
    # Plot the data - RMSF against residue index
    fig= plt.figure(constrained_layout=True)
    if do_rolling_Rgs:
        spec = gridspec.GridSpec(4,1,height_ratios=[1,4,1,0.5], figure=fig)
    else:
        spec = gridspec.GridSpec(3,1,height_ratios=[1,4,0.5], figure=fig)

    fig.set(figwidth=10, figheight=10)
    ax0 = fig.add_subplot(spec[0,0])
    ax1 = fig.add_subplot(spec[1,0])


    ax0.set_xlabel('Time (ns)')
    ax0.set_ylabel('Rg  ($\AA$)')

    # Set axis limits
    ax0.set_xlim(0, len(radius.T)/100 )
    start, end = ax0.get_xlim()
    ax0.xaxis.set_ticks(np.arange(start, end, 10))
    ax0.grid(visible = True, linestyle = '--', alpha=0.4)


    #ax.vlines(10, start, end, colors='k', linestyles='dotted')

    # plt median and + 1 stdev
    if len(radius[0]) > 1000:
        factor = 1000
    else:
        factor = 10

    plt_smooth(ax0,radius[1], radius[0],factor)
    ax0.axvspan(contact_start, contact_finish, color='red', alpha=0.5)

    distances_average = np.mean(distances_3Darray[contact_start *100 : contact_finish *100], axis=0)  
    #distances_average = np.rot90(distances_average,3) 
    im2 = ax1.pcolormesh(distances_average,vmin=1, vmax=dist_max, cmap = plt.cm.inferno_r) 

    start, end = ax1.get_xlim()
    ax1.xaxis.set_ticks(np.arange(start, end, 5))
    start, end = ax1.get_ylim()
    ax1.yaxis.set_ticks(np.arange(start, end, 5))
    ax1.grid(visible = True, linestyle = '--', alpha=0.4)

    plt.colorbar(im2)

    if do_rolling_Rgs:
        ax2 = fig.add_subplot(spec[2,0], sharex=ax1)
        rolling_Rgs = get_rolling_Rgs(atomistic_system, contact_start*100, contact_finish*100, residues_window)

        n_resids =  len(atomistic_system.residues.resids)
        ax2.xaxis.set_ticks(np.arange(0,n_resids , 5))
        ax2.grid(visible = True, linestyle = '--', alpha=0.4)
        ax2.plot(rolling_Rgs[0], rolling_Rgs[1])
    
    ss = fig.add_subplot(spec[-1,0], sharex=ax1)
    ss_subplot(fig, ss, dsspline_start)

    plt.legend()


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