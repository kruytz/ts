import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import time

pd.read_csv("valid_serves.csv").to_parquet("valid_serves.parquet")
from scipy.stats import gaussian_kde

try:
    from numba import njit, prange
    numba_available = True
except ImportError:
    numba_available = False
    print("Info: Numba not installed. The fastest CPU option will be disabled. Run 'pip install numba'.")


@njit(parallel=True)
def numba_kde(data_points, grid_points, bandwidth):
    n_data, d = data_points.shape; n_grid, _ = grid_points.shape
    density = np.zeros(n_grid)
    for i in prange(n_grid):
        diff = (grid_points[i] - data_points) / bandwidth
        energy = np.sum(diff**2, axis=1) / 2.0
        density[i] = np.sum(np.exp(-energy))
    norm_factor = n_data * np.prod(bandwidth) * (2 * np.pi)**(d / 2.0)
    density /= norm_factor
    return density


def draw_court_matplotlib(ax):
    net_position = 9.0; court_length = 18.0; court_width = 9.0; attack_line_dist = 3.0
    opponent_court = patches.Rectangle((net_position, -court_width / 2), net_position, court_width, linewidth=2, edgecolor='black', facecolor='none', zorder=2)
    ax.add_patch(opponent_court)
    server_court = patches.Rectangle((0, -court_width / 2), net_position, court_width, linewidth=2, edgecolor='black', facecolor='whitesmoke', zorder=0)
    ax.add_patch(server_court)
    ax.plot([net_position, net_position], [-court_width / 2, court_width / 2], 'k-', linewidth=3, zorder=3, label='Net')
    ax.plot([net_position - attack_line_dist, net_position - attack_line_dist], [-court_width/2, court_width/2], 'k--', linewidth=1, zorder=2)
    ax.plot([net_position + attack_line_dist, net_position + attack_line_dist], [-court_width/2, court_width/2], 'k--', linewidth=1, zorder=2)
    ax.text(13.5, 0, "Opponent's Side", ha='center', va='center', fontsize=14, alpha=0.5, zorder=1)
    ax.set_xticks(np.arange(0, 19, 1))
    ax.set_yticks(np.arange(-5, 6, 1))


def plot_heatmap(data, method='hist'):
    """Plots a heatmap using various CPU-based methods."""
    final_x = data['final_x']; final_z = data['final_z']
    fig, ax = plt.subplots(figsize=(13, 8))
    draw_court_matplotlib(ax)
    start_time = time.time()
    
    cmap_choice = 'inferno'

    if method == 'hist_highres':
        print("Generating High-Resolution Smoothed Histogram...")
        bins_x = 72; bins_z = 36
        counts, xedges, yedges = np.histogram2d(final_x, final_z, bins=(bins_x, bins_z), range=[[9, 18], [-4.5, 4.5]])
        
        im = ax.imshow(counts.T, 
                       interpolation='gaussian',
                       origin='lower',
                       cmap=cmap_choice,
                       extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                       aspect='auto')
        
        cbar = fig.colorbar(im, ax=ax); cbar.set_label('Density of Serve Landings (Smoothed)')
        ax.set_title('Heatmap of Valid Serve Landings (HD Smoothed Histogram)', fontsize=16, fontweight='bold')

    elif method == 'kde_numba':
        if not numba_available:
            print("Numba is not available. Cannot run this option."); return
        print("Generating Kernel Density Estimate (Numba JIT on all CPU cores)...")
        data_points = data[['final_x', 'final_z']].values
        n_data, d = data_points.shape

        if n_data < 2:
            print(f"Error: Cannot perform KDE with {n_data} data points. Need at least 2.")
            return
        
        exponent = -1. / (d + 4)
        silverman_factor = np.exp(exponent * np.log(n_data))
        
        std_devs = data_points.std(axis=0)

        base_bandwidth = silverman_factor * std_devs
        
        if n_data >= 100000:
            adaptive_scaling_factor =1.9 - np.log10(n_data)/1.6 + np.log10(max(1, n_data / 1000)) #OG: 1
        elif n_data >= 10000 and n_data < 100000:
            adaptive_scaling_factor =1.5 - np.log10(n_data)/1.6 + np.log10(max(1, n_data / 1000)) #OG: 1
        else:
            adaptive_scaling_factor =0.9 + np.log10(max(1, n_data / 1000)) #OG: 1
        bandwidth = base_bandwidth * adaptive_scaling_factor
        
        print(f"Using all {n_data} data points.")
        print(f"Numerically stable base bandwidth: ({base_bandwidth[0]:.4f}, {base_bandwidth[1]:.4f})")
        print(f"Adaptive scaling factor: {adaptive_scaling_factor:.4f}")
        print(f"Final applied bandwidth: ({bandwidth[0]:.4f}, {bandwidth[1]:.4f})")

        #xx, zz = np.mgrid[9:18:200j, -4.5:4.5:100j]

        if n_data > 90000:
            xx, zz = np.mgrid[9:18:200j, -4.5:4.5:100j]
        else:
            xx, zz = np.mgrid[9:18:400j, -4.5:4.5:200j]
            
        grid_points = np.vstack([xx.ravel(), zz.ravel()]).T
        f = numba_kde(data_points, grid_points, bandwidth)
        f_reshaped = f.reshape(xx.shape)
        
        contour = ax.contourf(xx, zz, f_reshaped, levels=100, cmap=cmap_choice, alpha=0.9)
        cbar = fig.colorbar(contour, ax=ax); cbar.set_label('Kernel Density Estimate (Numba)')
        ax.set_title(f'Heatmap (N-KDE on {len(data)} serves)', fontsize=16, fontweight='bold')

    print(f"Plot generation took: {time.time() - start_time:.4f} seconds")
    ax.set_xlabel('Horizontal Distance from Serve Line (meters)', fontsize=12)
    ax.set_ylabel('Sideways Distance from Center (meters)', fontsize=12)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(-1, 19)
    ax.set_ylim(-5.5, 5.5)
    ax.grid(True, linestyle=':', alpha=0.3)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # filename = "valid_serves.csv"  // OLD
    filename = "valid_serves.parquet" 

    try:
        # data = pd.read_csv(filename, engine='pyarrow') // OLD
        data = pd.read_parquet(filename)
        print(f"Successfully loaded '{filename}'. Found {len(data)} valid serves.\n")
    except Exception as e:
        print(f"Could not process file: {e}")
        exit()

    data.dropna(subset=['final_x', 'final_z'], inplace=True)
    if len(data) == 0:
        print("No valid data points found in the file to plot.")
        exit()

    print("--- Choose a visualization method ---")
    print("1. High-Res Smoothed Histogram (Fast, Recommended)")
    if numba_available:
        print("2. KDE with Numba (N-KDE)")
    
    choice_options = "1"
    if numba_available:
        choice_options += ", 2"
    
    choice = input(f"Enter your choice ({choice_options.replace(' ', '')}): ")
    
    if choice == '1':
        plot_heatmap(data, method='hist_highres')
    elif choice == '2' and numba_available:
        plot_heatmap(data, method='kde_numba')
    else:
        print("Invalid choice.")