
import numpy as np
# Import matplotlib - for plotting data
import matplotlib.pyplot as plt
# Import statistics Library
import statistics



def do_trajectory_CAalignement(atomistic_system, sim_path, trajectory_file_name):
    import MDAnalysis.analysis.align as align
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
    plt_median(ax, radius[1], label=True)
    # Show legend
    plt.legend()
    # Show plot
    plt.show()

    return Rgyr

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

def plot_distances(Rgyr, distances_3Darray, contact_start, contact_finish, dist_max):
    import matplotlib as mpl

    #Replot the rgyr with AOI
    radius = np.array(Rgyr)
    radius = radius.T
    # Plot the data - RMSF against residue index
    fig, ax = plt.subplots()
    fig.set(figwidth=12)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Rg  ($\AA$)')

    # Set axis limits
    ax.set_xlim(0, len(radius.T)/100 )
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, 10))
    ax.grid(visible = True, linestyle = '--', alpha=0.4)


    #ax.vlines(10, start, end, colors='k', linestyles='dotted')

    # plt median and + 1 stdev
    if len(radius[0]) > 1000:
        factor = 1000
    else:
        factor = 10

    plt_smooth(ax,radius[1], radius[0],factor)
    plt.axvspan(contact_start, contact_finish, color='red', alpha=0.5)
    # Show legend
    plt.legend()
    # Show plot

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
    axes[0].set_title('Initial 10 frames (0.1ns)')
    axes[1].set_title('Final 10  frames (0.1ns)')
    axes[2].set_title(f'Average over {contact_finish - contact_start}ns')

    # colorbar
    cax,kw = mpl.colorbar.make_axes([ax for ax in axes.flat])
    plt.colorbar(im2, cax=cax, **kw)

def plot_distances_HD(distances_3Darray, contact_start, contact_finish, dist_max):
    import matplotlib as mpl
    fig, ax = plt.subplots()
    fig.set(figwidth=12, figheight=10)

    distances_average = np.mean(distances_3Darray[contact_start *100 : contact_finish *100], axis=0)   
    im2 = ax.pcolormesh(distances_average,vmin=1, vmax=dist_max, cmap = plt.cm.inferno_r) 

    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, 5))
    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(np.arange(start, end, 5))
    ax.grid(visible = True, linestyle = '--', alpha=0.4)
    plt.colorbar(im2)


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