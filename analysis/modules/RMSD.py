
import matplotlib.pyplot as plt
from .plot import plt_smooth, plt_median


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