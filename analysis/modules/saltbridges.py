import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from .plot import plt_smooth

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