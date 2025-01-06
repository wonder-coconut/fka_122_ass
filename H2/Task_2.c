void task2_sampling(int argc, char *argv[], gsl_rng *r) {
    // File outputs for energy and sampled configurations
    FILE *op_energy = fopen("op_text/mcmc_energy_t2.txt", "w");
    FILE *op_samples = fopen("op_text/samples.txt", "w");

    // Input parameters
    int N = atoi(argv[2]);    // Total number of samples
    double d = atof(argv[3]); // Step size for MCMC moves

    // Initial positions and trial wavefunction parameter
    double r1[] = {-100.0, 100.0, -100.0}; 
    double r2[] = {100.0, -100.0, 100.0};
    double alpha = 0.1;

    mcmc_step iteration;

    // Variables to track energy and equilibration
    double energy = 0.0, previous_energy = 0.0;
    int eq_steps = 0;
    double tolerance = 1e-3; // Tolerance for equilibration
    bool equilibrated = false;

    // Equilibration phase
    for (int i = 0; i < 1E5; i++) {
        iteration = mcmc_displace_all(r1, r2, d, alpha, r);
        if (i > 0 && fabs(iteration.energy - previous_energy) < tolerance) {
            eq_steps = i;
            equilibrated = true;
            break;
        }
        previous_energy = iteration.energy;
    }

    if (!equilibrated) {
        printf("Equilibration not reached within %d steps. Using max equilibration.\n", 1E5);
        eq_steps = 1E5; // Default equilibration steps
    }
    printf("Equilibration steps: %d\n", eq_steps);

    // Sampling phase
    for (int i = 0; i < N; i++) {
        iteration = mcmc_displace_all(r1, r2, d, alpha, r);
        fprintf(op_energy, "%f\n", iteration.energy);
        fprintf(op_samples, "%f\t%f\t%f\t%f\t%f\t%f\n", r1[0], r1[1], r1[2], r2[0], r2[1], r2[2]);
        energy += iteration.energy;
    }

    printf("Acceptance ratio: %f\n", iteration.accepted * 1.0 / N);
    printf("Average energy: %f\n", energy / N);

    fclose(op_energy);
    fclose(op_samples);
}

double autocorrelation(double *data, int data_len, int lag) {
    double mean = average(data, data_len);
    double numerator = 0.0, denominator = 0.0;

    for (int i = 0; i < data_len - lag; i++) {
        numerator += (data[i] - mean) * (data[i + lag] - mean);
    }

    for (int i = 0; i < data_len; i++) {
        denominator += (data[i] - mean) * (data[i] - mean);
    }

    return numerator / denominator;
}

double calculate_statistical_inefficiency_autocorr(double *data, int data_len) {
    double inefficiency = 1.0;
    for (int lag = 1; lag < data_len; lag++) {
        double autocorr = autocorrelation(data, data_len, lag);
        if (autocorr < 0.01) break; // Threshold for insignificant correlation
        inefficiency += 2.0 * autocorr;
    }
    return inefficiency;
}

double block_average(double *data, int data_len, int block_size) {
    int m_b = data_len / block_size; // Number of blocks
    double *block_data = (double *)malloc(m_b * sizeof(double));

    for (int j = 0; j < m_b; j++) {
        block_data[j] = 0.0;
        for (int i = 0; i < block_size; i++) {
            block_data[j] += data[i + j * block_size];
        }
        block_data[j] /= block_size;
    }

    double inefficiency = block_size * variance(block_data, m_b) / variance(data, data_len);
    free(block_data);
    return inefficiency;
}

void task2_analysis() {
    FILE *ip_energy = fopen("op_text/mcmc_energy_t2.txt", "r");
    int N;
    fscanf(ip_energy, "%d\n", &N);

    double *energy_data = (double *)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        fscanf(ip_energy, "%lf\n", &energy_data[i]);
    }

    // Compute statistical inefficiency using autocorrelation
    double ineff_autocorr = calculate_statistical_inefficiency_autocorr(energy_data, N);
    printf("Autocorrelation-based inefficiency: %f\n", ineff_autocorr);

    // Compute statistical inefficiency using block averaging
    int block_size = sqrt(N);
    double ineff_block = block_average(energy_data, N, block_size);
    printf("Block averaging-based inefficiency: %f\n", ineff_block);

    free(energy_data);
    fclose(ip_energy);
}
