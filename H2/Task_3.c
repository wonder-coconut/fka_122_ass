void task3(int argc, char *argv[], gsl_rng *r) 
{
    FILE *energy_alpha_t3 = fopen("op_text/energy_evol_alpha_t3.txt", "w");

    int N = atoi(argv[2]);    // Number of Monte Carlo samples
    double d = atof(argv[3]); // Step size for MCMC moves

    // Energy and statistical inefficiency parameters
    double r1[] = {0, 0, 0};
    double r2[] = {0, 0, 0};
    double alpha, energy_sum, mean_energy, variance, sd, e_err;
    double alpha_ll = 0.05, alpha_ul = 0.25, alpha_step = 0.0025;

    double *energies = (double *)malloc(N * sizeof(double));
    double n_s = 19; // Hardcoded statistical inefficiency from Problem 2

    // Loop over \alpha values
    for (alpha = alpha_ll; alpha <= alpha_ul; alpha += alpha_step) {
        energy_sum = 0.0;

        // Sampling phase
        for (int i = 0; i < N; i++) {
            mcmc_step iteration = mcmc_displace_all(r1, r2, d, alpha, r);
            energies[i] = iteration.energy;
            energy_sum += iteration.energy;
        }

        // Compute mean energy
        mean_energy = energy_sum / N;

        // Compute variance
        variance = 0.0;
        for (int i = 0; i < N; i++) {
            variance += pow(energies[i] - mean_energy, 2);
        }
        variance /= (N - 1); // Sample variance

        // Standard deviation and error calculation
        sd = sqrt(variance);
        e_err = sd / sqrt(N / n_s); // Propagate error using statistical inefficiency

        // Write results to file and console
        fprintf(energy_alpha_t3, "%f\t%f\t%.10f\n", alpha, mean_energy, e_err);
        printf("%f\t%f\t%.10f\n", alpha, mean_energy, e_err);
    }

    free(energies);
    fclose(energy_alpha_t3);
}
