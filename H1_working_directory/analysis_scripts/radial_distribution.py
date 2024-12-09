import numpy as np

def calculate_radial_distribution(positions, box_length, num_bins, particle_density):
    # To calculate radial distribution function
    num_particles = len(positions)
    max_distance = box_length / 2  # Maximum distance due to periodic boundary conditions
    bin_edges = np.linspace(0, max_distance, num_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    g_r = np.zeros(num_bins)

    # Compute pair distances
    for i in range(num_particles):
        for j in range(i + 1, num_particles):
            # Minimum image convention
            delta = np.abs(positions[i] - positions[j])
            delta = np.where(delta > box_length / 2, box_length - delta, delta)
            distance = np.linalg.norm(delta)

            if distance < max_distance:
                bin_index = int(distance / max_distance * num_bins)
                g_r[bin_index] += 2  # Count both (i, j) and (j, i)

    # Normalize g(r)
    bin_volumes = 4 / 3 * np.pi * (bin_edges[1:] ** 3 - bin_edges[:-1] ** 3)
    normalization = num_particles * particle_density * bin_volumes
    g_r /= normalization

    return bin_centers, g_r


# Inputs
if __name__ == "__main__":
    # Load particle positions (Nx3 array)
    positions = np.loadtxt("anal_files/positions_liquid.txt")  # Use sim_data
    lp = 4.236854996336203
    box_length = 4 * lp  # Simulation box length (in units of distance)
    num_bins = 100  # Number of bins for the histogram
    num_particles = 256  # Number of particles
    volume = box_length ** 3  # Volume of the simulation box
    particle_density = num_particles / volume

    # Calculate g(r)
    radii, g_r = calculate_radial_distribution(positions, box_length, num_bins, particle_density)

    # Save and plot the results
    np.savetxt("radial_distribution.txt", np.column_stack((radii, g_r)))

    import matplotlib.pyplot as plt
    plt.plot(radii, g_r, label="Radial Distribution Function")
    plt.xlabel("r (distance)")
    plt.ylabel("g(r)")
    plt.title("Radial Distribution Function for Liquid State")
    plt.legend()
    plt.savefig("radial_dist.png")
    plt.show()
