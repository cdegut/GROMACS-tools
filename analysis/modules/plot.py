import matplotlib.pyplot as plt
import statistics
import numpy as np

###############################
#####
def ss_subplot(fig, ss: plt.Axes, dsspline_str) -> plt.Axes:
    ss.axhline(y = 4, color = 'y', linestyle = '-', label="Random Coil")
    ss.axis("off")
    ss.grid(visible = True, linestyle = '--', alpha=0.4)

    alpha_markers = []
    beta_markers = []

    for i,  res in  enumerate(dsspline_str):
        if res == "H":
            alpha_markers.append(i)
        if res == "E":
            beta_markers.append(i)

    alpha_markers = find_consecutive_numbers(alpha_markers, 3)
    beta_markers = find_consecutive_numbers(beta_markers, 4)
    ss.scatter(alpha_markers,[4 for x in alpha_markers], 
               zorder=10, color="r", marker="8", linewidths=8, label="Alpha helix (>3)")
    ss.scatter(beta_markers,[4 for x in beta_markers], 
               zorder=9, color="limegreen", marker=">", linewidths=5, label="Beta Strand (>4)")
    
    return ss

def find_consecutive_numbers(numbers, lenght):
    consecutive_numbers = []
    current_series = []

    for num in numbers:
        if not current_series or num == current_series[-1] + 1:
            current_series.append(num)
        else:
            if len(current_series) >= 4:
                consecutive_numbers.extend(current_series)
            current_series = [num]

    # Check if the last series is also consecutive
    if len(current_series) >= lenght:
        consecutive_numbers.extend(current_series)

    return consecutive_numbers


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


def edr_plot(energy_like_terms, ax1, data_label, unit=False, divide=False):
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
