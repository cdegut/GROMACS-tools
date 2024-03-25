
import numpy as np
# Import matplotlib - for plotting data
import matplotlib.pyplot as plt
# Import statistics Library
import statistics
import MDAnalysis as mda
import pandas as pd 
import MDAnalysis.analysis.align as align
import math



def do_trajectory_CAalignement(atomistic_system, sim_path, trajectory_file_name):
    
    average = align.AverageStructure(atomistic_system, atomistic_system, select='protein and name CA', ref_frame=0)
    print("Calculating averaged structure:")
    average.run(verbose=True)
    averaged_ref = average.results.universe
    # Align all structure on the averaged one
    aligner = align.AlignTraj(atomistic_system, 
                                  averaged_ref, 
                                  select='protein and name CA',
                                  in_memory=False,
                                  filename=sim_path + trajectory_file_name.split(".")[0] + "_aligned.xtc")
    
    print("align trajectory on the averaged one, save as " + sim_path + trajectory_file_name.split(".")[0] + "_aligned.xtc")
    aligner.run(verbose=True)



## Def some plotting functions
##################################
def nicey_ax(ax):
    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(np.arange(0, end, 10))
    ax.grid(visible = True, linestyle = '--', alpha=0.4)

def plt_median(ax,array, n=1, positive=True, negative=True, label=False, inside=False):
    # plt median and +- stdev
    median = statistics.median(array)
    stdev = statistics.stdev(array[2:]) #skip 2 1st frame
    
    if label:
        if inside:
            ax.text(x=0.05, y=0.15, ha="left", va="top", s=f"Median = {median:.2f}\nStdev = {stdev:.2f}", transform = ax.transAxes)        
        else:
            ax.text(x=0, y=-0.07, ha="left", va="top", s=f"Median = {median:.2f}\nStdev = {stdev:.2f}", transform = ax.transAxes)


    ax.axhline(y = median, color = 'r', linestyle = '-', alpha = 0.5)
    if positive:
        ax.axhline(y = median + n*stdev, color = 'r', linestyle = '-', alpha=0.2)
    if negative:
        ax.axhline(y = median - n*stdev, color = 'r', linestyle = '-', alpha =0.2)
    
    return median, stdev


def plt_smooth(ax,data_array,time_array,window):
    ## Use convolution to smooth an array on a sliding window anď plot it on the ax plot
    avg = np.convolve(data_array, np.ones(window)/window, mode='valid')
    slice_start = int(window/2-1)
    slice_finish = int(window/2)
    ax.plot(time_array[slice_start:-slice_finish]/1000, avg, label= str(window) + "ps avg")


def edr_plot(energy_like_terms, ax1, data_label, unit=False, divide=False, edr_only = False):
    #get data
    edr_data = energy_like_terms.get_data(data_label)
    
    if divide:
        result = []
        for val in edr_data[data_label]:
            result.append(val/divide)
        edr_data[data_label] = result
        
    #plot data
    ax1.plot(edr_data['Time']/1000, edr_data[data_label], label = "Instant")

    # Plt axis
    ax1.set_xlabel('Time (ns)')
    if unit:
        ylabel = f"{data_label} ({unit})"
    else:
        ylabel = data_label
    ax1.set_ylabel(ylabel)

    ax1.set_title(data_label)

    plt_smooth(ax1,edr_data[data_label],edr_data['Time'],20)
    median, stdev = plt_median(ax1,edr_data[data_label], label = True, inside = True)

    # Set axis limits
    ax1.set_ylim(median - 5 * stdev, median + 5 * stdev)
    if not edr_only:
        ax1.set_xlim(0, len(edr_data['Time'])/100)  


def box_plot(energy_like_terms):
    fig, axs = plt.subplots(3, sharex=True)
    box_X_data = energy_like_terms.get_data('Box-X')
    box_Y_data = energy_like_terms.get_data('Box-Y')
    box_Z_data = energy_like_terms.get_data('Box-Z')
    axs[0].plot(box_X_data['Time']/1000, box_X_data['Box-X'])
    plt_smooth(axs[0],box_X_data['Box-X'],box_X_data['Time'],20)
    plt_median(axs[0],box_X_data['Box-X'], positive=False, negative=False)
    axs[0].set_title("X dimension")

    axs[1].plot(box_Y_data['Time']/1000, box_Y_data['Box-Y'])
    plt_smooth(axs[1],box_Y_data['Box-Y'],box_Y_data['Time'],20)
    plt_median(axs[1],box_Y_data['Box-Y'], positive=False, negative=False)
    axs[1].set_title("Y dimension")

    axs[2].plot(box_Z_data['Time']/1000, box_Z_data['Box-Z'])
    plt_smooth(axs[2],box_Z_data['Box-Z'],box_X_data['Time'],20)
    plt_median(axs[2],box_Z_data['Box-Z'], positive=False, negative=False)
    axs[2].set_title("Z dimension")
    axs[2].set(xlabel='Time (ps)')

    for ax in axs:
        ax.set(ylabel='Size ($\AA$)', xlim=(0, len(box_Z_data['Time'])/100))
    
    #return 
    xlim = len(box_Z_data['Time'])/100
    fig.set(figheight=9)
    # Set axis limits
    plt.xlim(0, xlim)
    plt.show()


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
    ax.axhline(y = expanded_Rg, color = 'g', linestyle = '-', alpha = 0.2, label="good solvant")
    ax.axhline(y = theta_Rg, color = 'g', linestyle = '-', alpha = 0.5, label="theta solvant")
    ax.axhline(y = collapsed_Rg, color = 'g', linestyle = '-', alpha = 0.2, label="bad solvant")
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
########## 
# Simple distance between selecion center of mass

def sele_distance(atomistic_system, seleA, seleB):

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
    plt_median(ax,distances[0], label=True)

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

def plot_distances_HD(Rgyr, distances_3Darray, contact_start, contact_finish, dist_max):
    import matplotlib.gridspec as gridspec

    #Replot the rgyr with AOI
    radius = np.array(Rgyr)
    radius = radius.T
    # Plot the data - RMSF against residue index
    fig= plt.figure(constrained_layout=True)
    spec = gridspec.GridSpec(2,1,height_ratios=[1,4], figure=fig)
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

################ Salt Bridges ###############
#############################################
def salt_bridges_detection(atomistic_system: mda.Universe, acidic, basic, salt_bridge_cutoff=4.5):
    
    from MDAnalysis.analysis import contacts
    import pandas as pd

    frames = len(atomistic_system.trajectory)
    distances_3Darray = np.zeros((frames,len(acidic), len(basic)))
    byresidues_charges_distances_3Darray = np.zeros( (frames, len(np.unique(acidic.resids)) , len(np.unique(basic.resids))) )

    total_contacts = []

    for ts in atomistic_system.trajectory:
        distances = contacts.distance_array(acidic.positions, basic.positions)    
        #dist_arr = np.zeros((len(acidic), len(basic)))
        df = pd.DataFrame(distances, index = acidic.resids, columns=basic.resids )
        reduced_df = df.groupby(level=0).min().T.groupby(level=0).min().T
        distances_3Darray[ts.frame, :, :] = distances
        byresidues_charges_distances_3Darray[ts.frame, :, :] = reduced_df
        
        current_contacts = reduced_df <= salt_bridge_cutoff
        total_contacts.append([current_contacts.sum().sum(), atomistic_system.trajectory.time])
        
    total_contacts = np.array(total_contacts).T

    return distances_3Darray, byresidues_charges_distances_3Darray, total_contacts

def salt_bridges_plot(total_contacts, salt_bridge_cutoff):
    fig, ax = plt.subplots()
    fig.set(figwidth=12)
    ax.plot(total_contacts[1]/1000, total_contacts[0], label='# Salt bridges', alpha = 0.4)
    # Add axis labels
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel(f'# Charge within {salt_bridge_cutoff} ($\AA$)')

    # Set axis limits and grid
    ax.set_xlim(0, len(total_contacts.T)/100 )
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, 10))
    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(np.arange(int(start), int(end)))
    ax.grid(visible = True, linestyle = '--', alpha=0.4)


    plt_smooth(ax,total_contacts[0], total_contacts[1], 1000)
    #plt_median(ax, total_contacts, label=True)
    # Show legend
    plt.legend()
    # Show plot
    plt.show()

def get_AB_resnmid(basic: mda.AtomGroup, acidic: mda.AtomGroup):
    basic_resname = []
    acidic_resname = []
    for i, resids in enumerate(basic.resids):
        if basic.resnames[i] == "ARG":
            resn = "R"
        else:
            resn = "K"

        name = f"{resn}{basic.resids[i]}"
        if name not in basic_resname:
            basic_resname.append(name)

    for i, resids in enumerate(acidic.resids):
        if acidic.resnames[i] == "GLU":
            resn = "E"
        else:
            resn = "D"
        
        name = f"{resn}{acidic.resids[i]}"
        if name not in acidic_resname:
            acidic_resname.append(name)
    return acidic_resname, basic_resname

def salt_brdige_overtime(acidic: mda.AtomGroup, basic: mda.AtomGroup, charges_distances_3Darray, atomistic_system,contact_start, contact_finish, dist_max):
    import pandas as pd
    unique_acidic = np.unique(acidic.resids)
    unique_basic = np.unique(basic.resids)
    acidic_resname, basic_resname = get_AB_resnmid(basic, acidic)

    fig, axes = plt.subplots(3,1)
    fig.set(figwidth=12, figheight=30)
    fig.suptitle('Distance between charged atoms', fontsize=18)

    final_10frames = np.mean(charges_distances_3Darray[contact_finish *100 -10 :contact_finish *100], axis=0)
    initial_10frames = np.mean(charges_distances_3Darray[contact_start *100 : contact_start *100 + 10], axis=0)
    distances_average = np.mean(charges_distances_3Darray[contact_start *100 : contact_finish *100], axis=0)  


    def process_frames(frames):
        # generate full size array
        df = pd.DataFrame(frames, index = acidic.resids, columns=basic.resids ).groupby(level=0).min().T.groupby(level=0).min().T # regroup atoms pair
        full_df = pd.DataFrame(np.zeros((len(atomistic_system.residues.resids),len(atomistic_system.residues.resids))))
        merged = pd.concat([df, full_df]).groupby(level=0).max().T.groupby(level=0).max().T

        # restore diagonal symetry
        triu = np.triu(merged)
        tril = np.tril(merged)
        merged = triu.T + tril.T +merged
        merged[merged == 0] = 100 
        return merged


    im0 = axes[0].pcolormesh(process_frames(initial_10frames),vmin=3, vmax=dist_max, cmap = plt.cm.inferno_r)
    im1 = axes[1].pcolormesh(process_frames(final_10frames),vmin=3, vmax=dist_max, cmap = plt.cm.inferno_r)                                              
    im2 = axes[2].pcolormesh(process_frames(distances_average),vmin=3, vmax=dist_max, cmap = plt.cm.inferno_r)

    #ticks
    start, end = axes[0].get_xlim()
    tl =np.arange(start, end, 10)

    for ax in axes:
        ax.plot([0, 1], [0, 1], color="grey", transform=ax.transAxes)
        ax.plot(unique_acidic, unique_acidic, 'bo',ms = 4, alpha = 0.5)
        ax.plot(unique_basic, unique_basic, 'ro',ms = 4, alpha = 0.5)

        ax.xaxis.set_ticks(tl, labels= [""] * len(tl))
        ax.set_xticks(unique_acidic,acidic_resname,  rotation=45, minor=True)
        #ax.xticks(rotation=45, ha="right", minor=True)

        ax.yaxis.set_ticks(tl, labels= [""] * len(tl))
        ax.set_yticks(unique_basic, basic_resname, rotation=45, minor=True)

        ax.grid(visible = True, linestyle = '--', alpha=0.4)

    axes[0].set_title('Initial 10 frames (0.1ns)')
    axes[1].set_title('Final 10  frames (0.1ns)')
    axes[2].set_title(f'Average over {contact_finish - contact_start}ns')

    #cax,kw = mpl.colorbar.make_axes([ax for ax in axes.flat])
    #plt.colorbar(im2, cax=cax, **kw)
    plt.show()

def salt_brdige_singlepoint(acidic: mda.AtomGroup, basic: mda.AtomGroup, charges_distances_3Darray, atomistic_system, contact_start, dist_max):  
    
    def process_as_dataframes(frames):
        # generate full size array
        df = pd.DataFrame(frames, index = acidic.resids, columns=basic.resids ).groupby(level=0).min().T.groupby(level=0).min().T # regroup atoms pair
        full_df = pd.DataFrame(np.zeros((len(atomistic_system.residues.resids),len(atomistic_system.residues.resids))))
        merged = pd.concat([df, full_df]).groupby(level=0).max().T.groupby(level=0).max().T

        # restore diagonal symetry
        merged =merged + merged.T
        merged[merged == 0] = 100 
        return merged
    
    unique_acidic = np.unique(acidic.resids)
    unique_basic = np.unique(basic.resids)
    acidic_resname, basic_resname = get_AB_resnmid(basic, acidic)

    fig, ax = plt.subplots()
    fig.set(figwidth=12, figheight=10)
    fig.suptitle('Distance between charged atoms', fontsize=18)

    initial_10frames = np.mean(charges_distances_3Darray[contact_start *100 : contact_start *100 + 10], axis=0)

    im0 = ax.pcolormesh(process_as_dataframes(initial_10frames),vmin=3, vmax=dist_max, cmap = plt.cm.inferno_r)

    #ticks
    start, end = ax.get_xlim()
    tl =np.arange(start, end, 10)

    ax.plot([0, 1], [0, 1], color="grey", transform=ax.transAxes)
    ax.plot(unique_acidic, unique_acidic, 'bo',ms = 4, alpha = 0.5)
    ax.plot(unique_basic, unique_basic, 'ro',ms = 4, alpha = 0.5)

    ax.xaxis.set_ticks(tl, labels= [""] * len(tl), minor=True)
    ax.set_xticks(unique_acidic,acidic_resname,  rotation=45)

    ax.yaxis.set_ticks(tl, labels= [""] * len(tl), minor=True)
    ax.set_yticks(unique_basic, basic_resname, rotation=45)

    ax.grid(visible = True, linestyle = '--', alpha=0.4, which='minor')

    ax.set_title('Initial 10 frames (0.1ns)')

    plt.colorbar(im0)
    plt.show()

## RMSD ############
####################

def plot_RMSD(rmsd, time, RMSD_groups, RMSD_groups_name, plot_0):
    fig, axs = plt.subplots(len(RMSD_groups)+1)
    n = 2   
    for ax in axs:
        try:
            ax.plot(time, rmsd[n],  label="RMSD backbone", alpha=0.4)

            if len(rmsd[1]) > 10000:
                plt_smooth(ax,rmsd[n],rmsd[1],1000)
            else:
                plt_smooth(ax,rmsd[n],rmsd[1],20)

            #smoothed and median + stddev view
            plt_median(ax,rmsd[n],label=True)

            # Add axis labels
            ax.set_ylabel("RMSD  ($\AA$)")

            ax.set_title(RMSD_groups_name[n-2])

            # Set axis limits
            ax.set_xlim(plot_0,len(rmsd[n])/100)
            ax.set_ylim(min(rmsd[n][plot_0+2:]),round(max(rmsd[n]),1)+0.1)
            n = n + 1
        except:
            print("Group definition or name probably wrong")
            break

    axs[-1].set_xlabel("Time (ns)")
    # Add legend
    plt.legend()
    fig.set(figheight=14)
    # Show the plot
    plt.show()