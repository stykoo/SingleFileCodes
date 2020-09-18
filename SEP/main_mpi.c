#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include "simul.h"

int parseArgs(int argc, char **argv, long *n_sites, long *n_parts,
		      double *prob, long *n_simuls, char *fname,
			  int *determ, int *n_times, double **tfins) {
	if (argc < 8) {
		fprintf(stderr,
				"Usage: %s n_sites n_parts prob n_simuls fname determ tfin1"
			    " [tfin2 ...]\n",
		        argv[0]);
		return 1;
	}
	*n_sites = atol(argv[1]);
	*n_parts = atol(argv[2]);
	*prob = atof(argv[3]);
	*n_simuls = atol(argv[4]);
	strcpy(fname, argv[5]);
	// Deterministic initial conditions or not
	*determ = atoi(argv[6]);

	*n_times = argc - 7;
	*tfins = malloc(*n_times * sizeof(double));
	for (long k = 0 ; k < *n_times ; ++k) {
		(*tfins)[k] = atof(argv[7 + k]);
	}
	return 0;
}

void export(long * observables[N_OBS], char *fname, long n_sites,
		    long n_simul_tot, char *header) {
	FILE *file = NULL;
	file = fopen(fname, "w+");
	//assert(file);
	fprintf(file, "%s", header);
	fprintf(file, "# r er e+*er e-*er");
	fprintf(file, " X*er X*e+*er X*e-*er X");
#ifdef LARGE_NOBS
	fprintf(file, " X^2*er X^2*e+*er X^2*e-*er X^2");
	fprintf(file, " X^3*er X^3*e+*er X^3*e-*er X^3");
#endif
	fprintf(file, "\n");

	double fac = 1.0 / n_simul_tot;

	for (long i = 0 ; i < n_sites ; ++i) {
		fprintf(file, "%ld ", i);
		for (int k = 0 ; k < N_OBS ; ++k) {
			fprintf(file, " %.10e", fac * observables[k][i]);
		}
		fprintf(file, "\n");
	}
	fclose(file);
}

void exportCums(long * obs_sum[N_OBS], char *fname, long n_simul_tot,
                long n_times, double *tfins, char *header) {
	FILE *file = NULL;
	file = fopen(fname, "w+");

	fprintf(file, "%s\n", header);
	fprintf(file, "# t X X^2 X^3 X^4 X^5 X^6\n");
	double fac = 1.0 / n_simul_tot;

	for (long i = 0 ; i < n_times ; ++i) {
		fprintf(file, "%lf ", tfins[i]);
		for (int k = 0 ; k < N_OBS ; ++k) {
			fprintf(file, " %.10e", fac * obs_sum[i * N_OBS + k][0]);
		}
		fprintf(file, "\n");
	}

	fclose(file);
}

// From https://arxiv.org/abs/1005.4117
/*unsigned long seedgen1(void)  {
    unsigned long s, seed, pid;
    pid = getpid();
    s = time(NULL);
    seed = abs(((s*181)*((pid-83)*359))%104729);
    return seed;
}*/

// F*** it's hard to find reliable advice on how to seed parallel generators
// Let's read from /dev/urandom...
unsigned long seedgen2(void)  {
	int fd = open("/dev/urandom", O_RDONLY);
    unsigned long seed;
	read(fd, &seed, sizeof(unsigned long));
	close(fd);
    return seed;
}

int main(int argc, char **argv) {
	long n_sites, n_parts, n_simuls;
	int n_times;
	double prob;
	double *tfins;
	char fname[200];
	int determ, status;

	// Initializations
	status = parseArgs(argc, argv, &n_sites, &n_parts, &prob, &n_simuls,
			           fname, &determ, &n_times, &tfins);
	if (status) {
		return status;
	}

	int world_rank, world_size;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	char header[500];
	if (world_rank == 0) {
		sprintf(header,
				"# n_sites=%ld, n_parts=%ld, prob=%lf, "
				"n_simuls=%dx%ld, fname=%s, determ=%d",
			    n_sites, n_parts, prob, world_size, n_simuls,
				fname, determ);
		printf("%s\nTimes: ", header);
		for (int k = 0 ; k < n_times ; ++k) {
			printf("%g ", tfins[k]);
		}
		printf("\nNumber of processors: %d\n", world_size);
	}

	const int n_obs_tot = N_OBS * n_times;
	long **observables = (long **) malloc(n_obs_tot * sizeof(long *));
	// It is crucial to allocate the array of pointers on all processes
	long **obs_sum = (long **) malloc(n_obs_tot * sizeof(long *));
	assert(observables);
	assert(obs_sum);
#ifdef ONLY_CUMS
	long n_sites_eff = 1;
#else
	long n_sites_eff = n_sites;
#endif
	for (int k = 0 ; k < n_obs_tot ; ++k) {
		observables[k] = (long *) calloc(n_sites_eff, sizeof(long));
		assert(observables[k]);
		if (world_rank == 0) {
			obs_sum[k] = (long *) calloc(n_sites_eff, sizeof(long));
			assert(obs_sum[k]);
		}
	}
	// End of initializations

	// Run the simulation
	// Generate the seed
	unsigned long local_seed = seedgen2();

	run(n_sites, n_parts, prob, n_simuls, n_times, tfins, observables,
		local_seed, determ);

	// We sum all the observables from all the processes
	for (int k = 0 ; k < n_obs_tot ; ++k) {
		MPI_Reduce(observables[k], obs_sum[k], n_sites_eff, MPI_LONG, MPI_SUM,
		           0, MPI_COMM_WORLD);
	}

	// Export
	if (world_rank == 0) {
#ifdef ONLY_CUMS
		exportCums(obs_sum, fname, n_simuls * world_size, n_times,
		           tfins, header);
#else
		for (int k = 0 ; k < n_times ; ++k) {
			char hd[500], ff[100];
			sprintf(hd, "%s, tfin=%g\n", header, tfins[k]); // Header
			sprintf(ff, "%s_%g.dat", fname, tfins[k]); // Filename
			export(&obs_sum[k * N_OBS], ff, n_sites, n_simuls * world_size, hd);
		}
#endif
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	// Free memory
	free(tfins);
	for (int k = 0 ; k < n_obs_tot ; ++k) {
		free(observables[k]);
		if (world_rank == 0)
			free(obs_sum[k]);
	}
	free(observables);
	free(obs_sum);

	return 0;
}
