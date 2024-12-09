import numpy as np

def calculate_heat_capacity(energy_data, temperature, num_particles, boltzmann_constant=8.617341273302701E-05):
    
    # Calculate the mean energy and energy fluctuation
    energy_mean = np.mean(energy_data)
    energy_squared_mean = np.mean(energy_data ** 2)
    
    # Variance of energy (fluctuation)
    energy_variance = energy_squared_mean - energy_mean ** 2
    
    # Heat capacity
    c_v = (energy_variance / (boltzmann_constant * temperature ** 2)) * num_particles
    return c_v

# Inputs
energy_solid = np.loadtxt("anal_files/energy_solid.txt")  # Sim_data for the solid state
energy_liquid = np.loadtxt("anal_files/energy_liquid.txt")  # Sim_data for the liquid state

#Unit conversion if required otherwise comment it 
temperature_solid = 500 + 273.15  # Convert from Celsius to Kelvin
temperature_liquid = 700 + 273.15  # Convert from Celsius to Kelvin

num_particles = 256  # Number of atoms in the system

# Calculate C_V for solid and liquid states
c_v_solid = calculate_heat_capacity(energy_solid, temperature_solid, num_particles)
c_v_liquid = calculate_heat_capacity(energy_liquid, temperature_liquid, num_particles)

print(f"Heat Capacity at Constant Volume (Solid): {c_v_solid:.3e} eV/K")
print(f"Heat Capacity at Constant Volume (Liquid): {c_v_liquid:.3e} E/K")
