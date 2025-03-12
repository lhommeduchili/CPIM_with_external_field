/* C code to expand the Contact Process (CP) code to incorporate 
   a basic Ising Model (IM) into the occupancy states of the CP */

#include <stdlib.h>
#include "mt64.h"    /* Pseudo-random number generation MT library (64 bit) */
#include <math.h>    /* Math to transform random numbers from continuous to discrete */
#include <time.h>    /* Used to seed pseudo-random number generator */
#include <stdio.h>

/* Lattice Size */
#define X_SIZE 256
#define Y_SIZE 256

/* Defaults */
/* Number of replicas for each parameter set  */
#define REPLICAS 1
/* Default birth/colonization rate/probability and scale ranges */
#define BETA  1.0

/* Default mortality/extinction rate/probability and scale ranges */
#define DELTA  0.0

/* Default differentiation rate/probability and scale ranges */
#define ALPHA  1.0

/* Strength of the coupling in positive terms (J = -1*COUPLING kBT units)
   we should have 1/2 if we do not want to double count pairs,
   therefore use a positive number! */
#define COUPLING (1)

/* Default temperature and scale ranges */
#define TEMPERATURE 2.27

/* Default interaction radius */
#define RADIUS 1

/* Default initial condition chosen */
#define INIT 6

/* Default external field (B) and scale ranges */
#define EXTERNAL_FIELD 0.0


/* Structure with the simulation data */
struct simulation {
    int lattice_configuration[X_SIZE][Y_SIZE]; /* Store lattice configuration */
    int init_option;            /* Choice of initial condition */
    int generation_time;        /* Generations simulated */
    int Ising_neighborhood;     /* Ising Neighborhood: r=1 (NN), r=2 (NNN), etc */
    int num_neighbors;          /* Effective number of neighbors a la Von Neumann 
                                   in the Ising-like interactions */
    int occupancy;              /* Lattice occupancy */
    int vacancy;                /* Lattice vacancy */
    int up;                     /* Number of spins in the up (+1) state */
    int down;                   /* Number of spins in the down (-1) state */
    double birth_rate;          /* Contact Process' birth rate */
    double death_rate;          /* Contact Process' death rate */
    double differentiation_rate;/* Differentiation into spin state */
    double T;                   /* Ising's temperature */
    double J;                   /* Ising's coupling: ferro (-kB) or 
                                   anti-ferro (+kB) */
    double B;                   /* Ising's external field (B) */
    double lambda_rate;         /* Contact-Ising Monte Carlo bias */
    double *transition_probs;  /* Pointer for transition probabilities' array */
} s; /* Instance s of the structure to hold the simulation */


/* Other functions */

/* Function to set the number of neighbors */
void set_num_neighbors(void)
{
    s.num_neighbors = pow(s.Ising_neighborhood, 2) + pow(s.Ising_neighborhood + 1, 2) - 1;
}


void init_transition_probs(void)
{
    int number_of_configs = 2 * s.num_neighbors + 1;
    s.transition_probs = (double *) malloc(number_of_configs * sizeof(double));
}


void set_transition_probs(void)
{   
    int number_of_configs = 2 * s.num_neighbors + 1;
    double spin_energy, spin_energy_diff;
    int i = 0;
    double *temp_trans_probs = realloc(s.transition_probs, number_of_configs * sizeof(double));
    
    s.transition_probs = temp_trans_probs;
    // g_print("Transition probabilities changed\n");
    for (int n = -s.num_neighbors; n <= s.num_neighbors; n++)
    {
        spin_energy = (n * s.J) - s.B;
        spin_energy_diff = -(2) * spin_energy;
        s.transition_probs[i] = exp (-spin_energy_diff/s.T);
        // g_print("Transition probability for neighborhood configuration %d = %f \n", n, s.transition_probs[i]);
        i++;
    }
}


double get_transition_prob(int spin_energy_difference)
{
    int idx =  (int) (spin_energy_difference/2) + s.num_neighbors;
    // g_print("Transition probability retrieved for spin E diff of %d = %f \n", spin_energy_difference, s.transition_probs[idx]);
    return s.transition_probs[idx];
}


/* Function to get the closest neighbors (separated by r sites) of a given site 
   in the lattice */
void get_neighborhood(int x, int y, int r, int* neighbors)
{
    int n = 0;
    /* Iterate over a square of r^2 sites centered at (i, j) = (0,0) */
    for (int i = -r; i <= r; i++) 
    {
        for (int j = -r; j <= r; j++) 
        {
            /* Skip the focal site at the origin */
            if (i == 0 && j == 0) 
            {
                continue;
            }
            /* Skip sites that are beyond (|i|+|j|) = r to match 
               von Neumann's neighborhood connectivity */
            else if (abs(i) + abs(j) > r) 
            {
                continue;
            }
            neighbors[n] = s.lattice_configuration[(X_SIZE + x + i) % X_SIZE]
                                                    [(Y_SIZE + y + j) % Y_SIZE];
            n++;
        }
    }
}


/* Function to compute the energy value of the site located at (x, y) */
double compute_energy(int x, int y)
{
    double energy;
    int spin = 0;
    double neighborhood_configuration = 0;
    int neighborhood[s.num_neighbors];

    if (s.lattice_configuration[x][y] == 1) 
    {
        spin = 1;
    }
    else if (s.lattice_configuration[x][y] == -1) 
    {
        spin = -1;
    }

    get_neighborhood(x, y, s.Ising_neighborhood, neighborhood);
  
    for (int n = 0; n < s.num_neighbors; n++) 
    {
        /* If neighbor has no spin, go to the next one */
        if (neighborhood[n] == 0 || neighborhood[n] == 2) 
        {
            continue;
        }
        else if (neighborhood[n] == 1) 
        {
            neighborhood_configuration += 1;
        }
        else if (neighborhood[n] == -1) 
        {
            neighborhood_configuration -= 1;
        }
    }
    energy = (double) spin * (s.J * neighborhood_configuration - s.B);
    
    return energy;
}


/* Update function */
void update_lattice(void)
{
    int random_neighbor_state, random_neighbor;
    double random_spin;
    // Energies
    double spin_energy, spin_energy_diff;
    // Probability of reactions
    double transition_probability;
    int random_x_coor, random_y_coor;

    /* For the Contact Process, we always consider NN interactions */
    for (int site = 0; site < (int)(Y_SIZE * X_SIZE); site++)
    {
        /* Pick a random focal site */
        random_x_coor = (int)floor(genrand64_real1() * X_SIZE);
        random_y_coor = (int)floor(genrand64_real1() * Y_SIZE);

        switch (s.lattice_configuration[random_x_coor][random_y_coor])
        {
            case 0: /* Site is empty */
                /* Choose a random neighbor from the num_neighbors possible ones */
                random_neighbor = (int)floor(genrand64_real3() * 4);
                switch (random_neighbor)
                {
                    case 0: /* South */
                        random_neighbor_state = s.lattice_configuration[random_x_coor][(int)((Y_SIZE + random_y_coor - 1) % Y_SIZE)];
                        break;
                    case 1: /* North */
                        random_neighbor_state = s.lattice_configuration[random_x_coor][(int)((Y_SIZE + random_y_coor + 1) % Y_SIZE)];
                        break;
                    case 2: /* East */
                        random_neighbor_state = s.lattice_configuration[(int)((X_SIZE + random_x_coor - 1) % X_SIZE)][random_y_coor];
                        break;
                    case 3: /* West */
                        random_neighbor_state = s.lattice_configuration[(int)((X_SIZE + random_x_coor + 1) % X_SIZE)][random_y_coor];
                        break;
                }
                /* If its random neighbor is occupied: put a copy at the focal site
                   with probability birth_rate * dt */
                if (genrand64_real2() < s.birth_rate)
                {
                    switch (random_neighbor_state)
                    {
                        case 2:
                            s.lattice_configuration[random_x_coor][random_y_coor] = 2;
                            s.occupancy++;
                            s.vacancy--;
                            break;
                        case 1:
                            s.lattice_configuration[random_x_coor][random_y_coor] = 1;
                            s.occupancy++;
                            s.vacancy--;
                            s.up++;
                            break;
                        case -1:
                            s.lattice_configuration[random_x_coor][random_y_coor] = -1;
                            s.occupancy++;
                            s.vacancy--;
                            s.down++;
                            break;
                        case 0:
                            s.lattice_configuration[random_x_coor][random_y_coor] = 0;
                            break;
                    }
                }
                break; /* break case 0 */

            case 2: /* Focal point is in the occupied, undifferentiated state */
                /* First, we check if the site survives */
                if (genrand64_real2() < s.death_rate)
                {
                    s.lattice_configuration[random_x_coor][random_y_coor] = 0;
                    s.occupancy--;
                    s.vacancy++;
                }
                else if (genrand64_real2() < s.differentiation_rate)
                {
                    /* Set an occupied site in the middle of the lattice */
                    random_spin = (int)((genrand64_int64() % 2) * 2) - 1;
                    if (random_spin == 1)
                    {
                        s.lattice_configuration[random_x_coor][random_y_coor] = random_spin;
                        s.up++;
                    }
                    else if (random_spin == -1)
                    {
                        s.lattice_configuration[random_x_coor][random_y_coor] = random_spin;
                        s.down++;
                    }
                }
                break;

            case 1: /* Focal point is in the up (+1) state */
                spin_energy = compute_energy(random_x_coor, random_y_coor);
                spin_energy_diff = -(2) * spin_energy;
                transition_probability = get_transition_prob(spin_energy_diff);
                if (genrand64_real2() < s.death_rate)
                {
                    s.lattice_configuration[random_x_coor][random_y_coor] = 0;
                    s.occupancy--;
                    s.vacancy++;
                    s.up--;
                }
                else if (spin_energy_diff < 0 || genrand64_real2() < transition_probability)
                {
                    s.lattice_configuration[random_x_coor][random_y_coor] = -1;
                    s.up--;
                    s.down++;
                }
                break;

            case -1: /* Focal point is in the down (-1) state */
                spin_energy = compute_energy(random_x_coor, random_y_coor);
                spin_energy_diff = -(2) * spin_energy;
                transition_probability = get_transition_prob(spin_energy_diff);
                if (genrand64_real2() < s.death_rate)
                {
                    s.lattice_configuration[random_x_coor][random_y_coor] = 0;
                    s.occupancy--;
                    s.vacancy++;
                    s.down--;
                }
                else if (spin_energy_diff < 0 || genrand64_real2() < transition_probability)
                {
                    s.lattice_configuration[random_x_coor][random_y_coor] = 1;
                    s.up++;
                    s.down--;
                }
                break;
        }
    }

    s.generation_time++;
}


/* Callback to initialize the lattice */
static void init_lattice(void)
{
    int random_spin;
    int x, y;
    /* Fill the lattice with 0s (unoccupied state) */
    for (x = 0; x < X_SIZE; x++)
    {
        for (y = 0; y < Y_SIZE; y++)
        {
            s.lattice_configuration[x][y] = 0;
        }
    }
    s.occupancy = 0;
    s.up = 0;
    s.down = 0;
    s.vacancy = X_SIZE * Y_SIZE;

    switch (s.init_option)
    {
        case 1:
            /* Set an occupied site in the middle of the lattice */
            random_spin = ((genrand64_int64() % 2) * 2) - 1;
            if (random_spin == 1)
            {
                s.lattice_configuration[X_SIZE / 2][Y_SIZE / 2] = random_spin;
                s.up++;
                s.vacancy--;
                s.occupancy++;
            }
            else if (random_spin == -1)
            {
                s.lattice_configuration[X_SIZE / 2][Y_SIZE / 2] = random_spin;
                s.down++;
                s.vacancy--;
                s.occupancy++;
            }
            break;

        case 2:
            /* Set an undifferentiated site in the middle of the lattice */
            s.lattice_configuration[X_SIZE / 2][Y_SIZE / 2] = 2;
            s.vacancy--;
            s.occupancy++;
            break;

        case 3:
            /* Set a small (r=2) cluster with undifferentiated sites in the middle of the lattice */
            for (x = X_SIZE / 2 - 2; x < X_SIZE / 2 + 2; x++)
            {
                for (y = X_SIZE / 2 - 2; y < X_SIZE / 2 + 2; y++)
                {
                    s.lattice_configuration[x][y] = 2;
                    s.occupancy++;
                    s.vacancy--;
                }
            }
            break;

        case 4:
            /* Set a small (r=2) cluster with random spins in the middle of the lattice */
            for (x = X_SIZE / 2 - 2; x < X_SIZE / 2 + 2; x++)
            {
                for (y = X_SIZE / 2 - 2; y < X_SIZE / 2 + 2; y++)
                {
                    random_spin = ((genrand64_int64() % 2) * 2) - 1;
                    if (random_spin == 1)
                    {
                        s.lattice_configuration[x][y] = random_spin;
                        s.up++;
                        s.vacancy--;
                        s.occupancy++;
                    }
                    else if (random_spin == -1)
                    {
                        s.lattice_configuration[x][y] = random_spin;
                        s.down++;
                        s.vacancy--;
                        s.occupancy++;
                    }
                }
            }
            break;

        case 5:
            /* Set up a lattice fully occupied with undifferentiated particles */
            for (x = 0; x < X_SIZE; x++)
            {
                for (y = 0; y < Y_SIZE; y++)
                {
                    s.lattice_configuration[x][y] = 2;
                    s.occupancy++;
                    s.vacancy--;
                }
            }
            break;
        
        case 6:
            /* Set up a lattice fully occupied with differentiated particles */
            random_spin = ((genrand64_int64() % 2) * 2) - 1;
            for (x = 0; x < X_SIZE; x++)
            {
                for (y = 0; y < Y_SIZE; y++)
                {
                    s.lattice_configuration[x][y] = random_spin;
                    s.occupancy++;
                    s.vacancy--;
                    if (random_spin == 1)
                    {
                        s.up++;
                    }
                    else if (random_spin == -1)
                    {
                        s.down++;
                    }
                }
            }
            break;
    }
    s.generation_time = 0;
}


static void initialize_simulation(void)
{
    /* Initialize Mersenne Twister algorithm for random number generation */
    unsigned int seed = (unsigned int)time(NULL);
    init_genrand64(seed);

    /* Set default parameters of the simulation */
    /* Initial condition option */
    s.init_option = (int)INIT;

    /* Contact Process */
    s.birth_rate = (double)BETA;
    // s.death_rate = (double)DELTA;

    /* Cell differentiation */
    s.differentiation_rate = (double)ALPHA;

    /* Ising Model */
    /* Interaction radius */
    s.Ising_neighborhood = (int)RADIUS;

    set_num_neighbors();

    /* Temperature */
    // s.T = (double)TEMPERATURE;

    /* Spin coupling */
    s.J = -1 * (double)COUPLING;

    /* External field */
    s.B = (double)EXTERNAL_FIELD;
    
    /* Initialize double array that holds transition probability values */
    init_transition_probs();
    /* Compute transition probabilities and store them in the array */
    set_transition_probs();
}


/* Main function spanning a Gtk Application object */
int main (int argc, char **argv)
    {
    /* Create a new text file and populate with the data headers */
    FILE *datafile;
    datafile = fopen(argv[1], "w");
    fprintf(datafile, "replica,lattice_size,birth_rate,death_rate,temperature,external_field,generation_time,occupancy,up,down\n");
    fclose(datafile);

    /* The following variables are used for calculating the variance of 
    num_gens_averaged (e.g. 10) consecutive occupancy values */
    int num_gens_hotstart = 100000;
    int num_gens_averaged = 1000000;

    /* The general strategy is to perform REPLICAS times the whole swap of the 
    the parameters*/
	for (int N = 0; N < REPLICAS; N++)
	    {
		for (double temp = 3.0; temp >= 0.1; temp -= 0.1)
            {
            s.T = temp;
            for (double death_rate = 0.0; death_rate <= 0.50; death_rate += 0.1)
                {
                s.death_rate = death_rate;
                initialize_simulation ();
                init_lattice();
                /* We run the simulation for some generations... */
                while (s.generation_time < num_gens_hotstart)
                    {
                    update_lattice ();
                    }
                /* We now save the next num_gens_averaged generation times' data. */
                for (int t = 0; t < num_gens_averaged + 1; t++)
                    {
                    if (t % 1000 == 0)
                        {
                        datafile = fopen(argv[1], "a+");
                        fprintf(datafile, "%d,%d,%f,%f,%f,%f,%d,%f,%d,%d\n", 
                                N, 
                                X_SIZE,
                                s.birth_rate,
                                s.death_rate,
                                s.T,
                                s.B,
                                s.generation_time, 
                                (double) s.occupancy / (double) (X_SIZE * Y_SIZE),
                                s.up,
                                s.down);
                        fclose(datafile);
                        }
                    update_lattice();
                    }
                }
            }
	    }
    return 0;
    }