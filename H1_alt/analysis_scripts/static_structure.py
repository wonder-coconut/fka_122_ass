import numpy as np
import matplotlib.pyplot as plt

def compute_static_structure_factor(positions, box_length, num_particles, q_vectors):
# To compute static structure factor
    num_q_vectors = len(q_vectors)
    S_q = np.zeros(num_q_vectors)

    for idx, q in enumerate(q_vectors):
        # Compute exp(i * q · r) for all particles
        phase_factors = np.exp(1j * np.dot(positions, q))
        S_q[idx] = np.abs(np.sum(phase_factors)) ** 2

    # Normalize S(q) by the number of particles
    S_q /= num_particles

    # Calculate magnitudes of q-vectors
    q_magnitudes = np.linalg.norm(q_vectors, axis=1)
    return q_magnitudes, S_q


def generate_q_vectors(box_length, max_q, grid_size):
# To Generate a 3D grid of q-vectors consistent with periodic boundary conditions.
    dq = 2 * np.pi / box_length  # Smallest spacing in q-space
    q_max_idx = int(max_q / dq)
    indices = np.arange(-q_max_idx, q_max_idx + 1)
    q_vectors = np.array([[dq * nx, dq * ny, dq * nz]
                          for nx in indices for ny in indices for nz in indices])
    q_vectors = q_vectors[np.linalg.norm(q_vectors, axis=1) <= max_q]
    return q_vectors


# Inputs
if __name__ == "__main__":
    # Load particle positions
    positions = np.loadtxt("positions_liquid.txt")  # Use sim_data
    box_length = 10.0  # Simulation box length (in units of distance)
    num_particles = len(positions)  # Number of particles
    max_q = 2 * np.pi / (box_length / 2)  # Maximum q value
    grid_size = 10  # Adjust grid size for resolution in q-space

    # Generate q-vectors
    q_vectors = generate_q_vectors(box_length, max_q, grid_size)

    # Compute static structure factor
    q_magnitudes, S_q = compute_static_structure_factor(positions, box_length, num_particles, q_vectors)

    # Save and plot results
    np.savetxt("static_structure_factor.txt", np.column_stack((q_magnitudes, S_q)))

    plt.plot(q_magnitudes, S_q, 'o', label="Static Structure Factor")
    plt.xlabel("|q| (wavevector magnitude)")
    plt.ylabel("S(q)")
    plt.title("Static Structure Factor for Liquid State")
    plt.legend()
    plt.show()
