/* C code to expand the Contact Process (CP) code to incorporate 
   a basic Ising Model (IM) into the occupancy states of the CP */

#include <stdlib.h>
#include <gtk/gtk.h> /* GUI, Gtk library */
#include "mt64.h"    /* Pseudo-random number generation MT library (64 bit) */
#include <math.h>    /* Math to transform random numbers from continuous to discrete */
#include <time.h>    /* Used to seed pseudo-random number generator */
#include <stdio.h>

/* Lattice Size */
#define X_SIZE 250
#define Y_SIZE 250

/* Defaults */
#define SAMPLE_RATE 1
/* Default birth/colonization rate/probability and scale ranges */
#define BETA  1.0
#define BETA_STEP 0.00001
#define BETA_MIN 0.000000
#define BETA_MAX 1.0
/* Default mortality/extinction rate/probability and scale ranges */
#define DELTA  0.00001
#define DELTA_STEP  0.00001
#define DELTA_MIN 0.000000
#define DELTA_MAX 1.0
/* Default differentiation rate/probability and scale ranges */
#define ALPHA  1.0
#define ALPHA_STEP 0.01
#define ALPHA_MIN  0.00
#define ALPHA_MAX  1.0
/* Strength of the coupling in positive terms (J = -1*COUPLING kBT units)
   we should have 1/2 if we do not want to double count pairs,
   therefore use a positive number! */
#define COUPLING (1)
/* Default temperature and scale ranges */
#define TEMPERATURE 2.27
#define TEMPERATURE_STEP 0.0000001
#define TEMPERATURE_MIN  0.0000001
#define TEMPERATURE_MAX  15 
/* Default interaction radius */
#define RADIUS 1
/* Default initial condition chosen */
#define INIT 1
/* Default external field (B) and scale ranges */
#define EXTERNAL_FIELD 0.0
#define EXTERNAL_FIELD_STEP 0.0000001
#define EXTERNAL_FIELD_MIN -1.0
#define EXTERNAL_FIELD_MAX 1.0

/* Structure with the simulation data */
struct simulation {
    int lattice_configuration[X_SIZE][Y_SIZE]; /* Store lattice configuration */
    gint run;                   /* Time handler tag */
    gboolean running;           /* Are we running? */
    int init_option;            /* Choice of initial condition */
    gboolean initialized;       /* Have we been initialized? */
    int generation_time;        /* Generations simulated */
    int Ising_neighborhood;     /* Ising Neighborhood: r=1 (NN), r=2 (NNN), etc */
    int num_neighbors;          /* Effective number of neighbors a la Von Neumann 
                                   in the Ising-like interactions */
    int occupancy;              /* Lattice occupancy */
    int vacancy;                /* Lattice vacancy */
    int up;                     /* Number of spins in the up (+1) state */
    int down;                   /* Number of spins in the down (-1) state */
    int display_rate;           /* Display rate: to paint the lattice */
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


/* GDK's pixel buffer functions implemented at the end of document */
void put_pixel(GdkPixbuf *pixbuf, int x, int y, 
               guchar red, guchar green, guchar blue, guchar alpha);
static void paint_a_background (gpointer data);
static void paint_lattice (gpointer data);


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
    g_print("Transition probabilities changed\n");
    for (int n = -s.num_neighbors; n <= s.num_neighbors; n++)
    {
        spin_energy = (n * s.J) - s.B;
        spin_energy_diff = -(2) * spin_energy;
        s.transition_probs[i] = exp (-spin_energy_diff/s.T);
        g_print("Transition probability for neighborhood configuration %d = %f \n", n, s.transition_probs[i]);
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
int update_lattice(gpointer data)
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
    if (s.generation_time % s.display_rate == 0)
    {
        paint_lattice(data);
        g_print("Gen: %d \t Vacancy: %f \t Occupancy: %f \t Up: %f \t Down: %f\n",
                s.generation_time, (double)s.vacancy / (double)(Y_SIZE * X_SIZE), (double)s.occupancy / (double)(Y_SIZE * X_SIZE), (double)s.up / (double)(s.occupancy), (double)s.down / (double)s.occupancy);
    }

    return 0;
}


/* Time handler to connect update function to the gtk loop */
gboolean time_handler (gpointer data)
   { 
    update_lattice (data);
    return TRUE;
    }


/* Callback to initialize the lattice */
static void init_lattice(GtkWidget *widget, gpointer data)
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

    s.initialized = TRUE;
    s.generation_time = 0;
    paint_lattice(data);
    g_print("Lattice initialized\n");
}


/* Stop simulation control */
static void stop_simulation(gpointer data)
{
    if (s.running)
    {
        g_source_remove(s.run);
        s.running = FALSE;
        g_print("Simulation stopped\n");
    }
}


/* Callback to launch dialog with info (GitHub's wiki) */
static void on_button_show_about(GtkWidget *widget, gpointer data)
{
    GdkPixbuf *pixbuf = gdk_pixbuf_new_from_file("X-Institute_logo_small.tif", NULL);
    GtkWidget *dialog = gtk_about_dialog_new();
    gtk_about_dialog_set_program_name(GTK_ABOUT_DIALOG(dialog), "Contact Process Ising Model App");
    gtk_about_dialog_set_version(GTK_ABOUT_DIALOG(dialog), "version 3.0, 2024");
    gtk_about_dialog_set_copyright(GTK_ABOUT_DIALOG(dialog), "Open-Source Software");
    gtk_about_dialog_set_comments(GTK_ABOUT_DIALOG(dialog), "The Contact Process Ising Model (CPIM).");
    gtk_about_dialog_set_website(GTK_ABOUT_DIALOG(dialog), "https://github.com/jekeymer/Contact-Process-Ising-Model/wiki");
    gtk_about_dialog_set_logo(GTK_ABOUT_DIALOG(dialog), pixbuf);
    g_object_unref(pixbuf);
    gtk_dialog_run(GTK_DIALOG(dialog));
    gtk_widget_destroy(dialog);
}


/* Callback to start simulation */
static void on_button_start_simulation(GtkWidget *button, gpointer data)
{
    if (!s.running && s.initialized)
    {
        s.run = g_idle_add((GSourceFunc)time_handler, GTK_IMAGE(data));
        s.running = TRUE;
        g_print("Simulation started\n");
    }
}


/* Callback to stop simulation */
static void on_button_stop_simulation(GtkWidget *button, gpointer data)
{
    stop_simulation(data);
}


/* Callback to change initial conditions -- dirty
   get_active() method is cleaner as I could use only one handler */
/* Init 1 */
static void on_radio_initial_condition_1(GtkWidget *button, gpointer data)
{
    char *id_radio = (char *)data;
    g_print("%s\n", id_radio);
    s.init_option = 1;
}


/* Init 2 */
static void on_radio_initial_condition_2(GtkWidget *button, gpointer data)
{
    char *id_radio = (char *)data;
    g_print("%s\n", id_radio);
    s.init_option = 2;
}


/* Init 3 */
static void on_radio_initial_condition_3(GtkWidget *button, gpointer data)
{
    char *id_radio = (char *)data;
    g_print("%s\n", id_radio);
    s.init_option = 3;
}


/* Init 4 */
static void on_radio_initial_condition_4(GtkWidget *button, gpointer data)
{
    char *id_radio = (char *)data;
    g_print("%s\n", id_radio);
    s.init_option = 4;
}


/* Init 5 */
static void on_radio_initial_condition_5(GtkWidget *button, gpointer data)
{
    char *id_radio = (char *)data;
    g_print("%s\n", id_radio);
    s.init_option = 5;
}


/* Init 6 */
static void on_radio_initial_condition_6(GtkWidget *button, gpointer data)
{
    char *id_radio = (char *)data;
    g_print("%s\n", id_radio);
    s.init_option = 6;
}



/* Callback to change Ising NN (r=1) vs NNN (r=2) conditions -- dirty
   get_active() method is cleaner as we could use only one handler */
/* NN; r = 1 */
static void on_radio_NN(GtkWidget *button, gpointer data)
{
    char *id_radio = (char *)data;
    g_print("%s\n", id_radio);
    s.Ising_neighborhood = 1;
    set_num_neighbors();
    set_transition_probs();
}


/* NNN; r = 2 */
static void on_radio_NNN(GtkWidget *button, gpointer data)
{
    char *id_radio = (char *)data;
    g_print("%s\n", id_radio);
    s.Ising_neighborhood = 2;
    set_num_neighbors();
    set_transition_probs();
}


/* NNNN; r = 3 */
static void on_radio_NNNN(GtkWidget *button, gpointer data)
{
    char *id_radio = (char *)data;
    g_print("%s\n", id_radio);
    s.Ising_neighborhood = 3;
    set_num_neighbors();
    set_transition_probs();
}


/* Callback to change Ising J = -kB (ferro) vs J = +kB (anti-ferro) -- dirty
   get_active() method is cleaner as I could use only one handler */
/* Ferro: J = -kB */
static void on_radio_ferro(GtkWidget *button, gpointer data)
{
    char *id_radio = (char *)data;
    g_print("%s\n", id_radio);
    s.J = -1 * (float)COUPLING;
    set_transition_probs();
}


/* Anti-Ferro: J = +kB */
static void on_radio_anti_ferro(GtkWidget *button, gpointer data)
{
    char *id_radio = (char *)data;
    g_print("%s\n", id_radio);
    s.J = 1 * (float)COUPLING;
    set_transition_probs();
}


/* Callback to respond to Gtk scale slide move event */
static void birth_rate_scale_moved(GtkRange *range, gpointer user_data)
{
    gdouble pos = gtk_range_get_value(range);
    s.birth_rate = (float)pos;
}


/* Callback to respond to Gtk scale slide move event */
static void differentiation_rate_scale_moved(GtkRange *range, gpointer user_data)
{
    gdouble pos = gtk_range_get_value(range);
    s.differentiation_rate = (float)pos;
}


/* Callback to respond to Gtk scale slide move event */
static void death_rate_scale_moved(GtkRange *range, gpointer user_data)
{
    gdouble pos = gtk_range_get_value(range);
    s.death_rate = (float)pos;
}


/* Callback to respond to Gtk scale slide move event */
static void temperature_scale_moved(GtkRange *range, gpointer user_data)
{
    gdouble pos = gtk_range_get_value(range);
    s.T = (float)pos;
    set_transition_probs();
}


/* Callback to respond to Gtk scale slide move event */
static void external_field_scale_moved(GtkRange *range, gpointer user_data)
{
    gdouble pos = gtk_range_get_value(range);
    s.B = (float)pos;
    set_transition_probs();
}


/* Callback to respond to Gtk scale slide move event */
static void display_rate_scale_moved(GtkRange *range, gpointer user_data)
{
    gdouble pos = gtk_range_get_value(range);
    s.display_rate = (int)pos;
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
    s.death_rate = (double)DELTA;

    /* Cell differentiation */
    s.differentiation_rate = (double)ALPHA;

    /* Ising Model */
    /* Interaction radius */
    s.Ising_neighborhood = (int)RADIUS;
    set_num_neighbors();
    /* Temperature */
    s.T = (double)TEMPERATURE;
    /* Spin coupling */
    s.J = -1 * (double)COUPLING;
    /* External field */
    s.B = (double)EXTERNAL_FIELD;

    /* Set simulation flags */
    s.running = FALSE;
    s.initialized = FALSE;
    /* Display rate to paint the lattice */
    s.display_rate = (int)SAMPLE_RATE;
    
    /* Initialize double array that holds transition probability values */
    init_transition_probs();
    /* Compute transition probabilities and store them in the array */
    set_transition_probs();
}


/* Activate function */
static void activate(GtkApplication *app, gpointer user_data)
{
    /* Initialize simulation */
    initialize_simulation();

    /* General Gtk widgets for the Window packing */
    GtkWidget *window, *grid, *image_lattice, *label, *frame, *notebook, *box;
    GtkWidget *scale, *radio, *separator;
    GdkPixbuf *pixbuf;

    /* Parameters Section */
    /* Create a Gtk Notebook to hold pages of parameters */
    notebook = gtk_notebook_new();
    /* Pre-Pack all Radio buttons in a box to be framed
    and then packed into the grid */
    grid = gtk_grid_new();

    /* Pre-pack: Coupling Neighborhood & Coupling Strength Radio Buttons */
    
    /* Create a Frame to back the Neighborhood size */
    frame = gtk_frame_new("Neighborhood size");
    gtk_frame_set_label_align(GTK_FRAME(frame), 0, 0);
    /* Create a box to hold stuff together */
    box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    /* Create a radio button for the NN (r=1) choice */
    radio = gtk_radio_button_new_with_label(NULL, "Nearest Neighbors (NN); r = 1");
    g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_NN), (gpointer)"NN interaction selected");
    gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
    /* Create a radio button for the NNN (r=2) choice */
    radio = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(radio), "Next Nearest Neighbors (NNN); r = 2");
    g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_NNN), (gpointer)"NNN interaction selected");
    gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
    /* Create another radio button for the NNNN (r=3) choice */
    radio = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(radio), "Next-Next Nearest Neighbors (NNNN); r = 3");
    g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_NNNN), (gpointer)"NNNN interaction selected");
    gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
    /* Add the packing box to a Neighborhood size Frame */
    gtk_container_add(GTK_CONTAINER(frame), box);
    gtk_grid_attach(GTK_GRID(grid), frame, 0, 0, 1, 1);

    /* Add a vertical separator to the parameter grid for order */
    separator = gtk_separator_new(GTK_ORIENTATION_VERTICAL);
    gtk_grid_attach(GTK_GRID(grid), separator, 1, 0, 1, 2);
    
    /* Create another Frame for the Coupling Strength
    (ferro vs anti-ferro magnetic) */
    frame = gtk_frame_new("Coupling strength");
    gtk_frame_set_label_align(GTK_FRAME(frame), 1, 0);
    /* Create a box to hold stuff together */
    box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    /* Create a radio button for the ferro magnetic case */
    radio = gtk_radio_button_new_with_label(NULL, "  J = - 1 * kB");
    g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_ferro), (gpointer)"Ferro-Magnetic (J < 0) interaction selected");
    gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
    /* Create a radio button for the anti-ferro magnetic case */
    radio = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(radio), "  J = +1 * kB");
    g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_anti_ferro), (gpointer)"Anti-Ferro-Magnetic (J > 0) interaction selected");
    gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
    /* Add the packing box to the Coupling strength Frame */
    gtk_container_add(GTK_CONTAINER(frame), box);
    gtk_grid_attach(GTK_GRID(grid), frame, 2, 0, 1, 2);

    /* Now we are done with pre-packing the radio choices in Frames */

    /* Contact Process (CP) Box */
    /* to control the parameters of the Contact Process */
    box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    /* BIRTH */
    /* Scale bar to set birth rate */
    scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, (gdouble)BETA_MIN, (gdouble)BETA_MAX, (gdouble)BETA_STEP);
    /* Set it to its default value */
    gtk_range_set_value(GTK_RANGE(scale), (gfloat)BETA);
    g_signal_connect(scale, "value-changed", G_CALLBACK(birth_rate_scale_moved), NULL);
    /* Pack it in a Frame */
    frame = gtk_frame_new("Birth rate");
    gtk_container_add(GTK_CONTAINER(frame), scale);
    /* Add that Frame to the CP box */
    gtk_container_add(GTK_CONTAINER(box), frame);
    /* DEATH */
    /* Scale bar to set death rate */
    scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, (gdouble)DELTA_MIN, (gdouble)DELTA_MAX, (gdouble)DELTA_STEP);
    /* Set it to its default value */
    gtk_range_set_value(GTK_RANGE(scale), (gfloat)DELTA);
    g_signal_connect(scale, "value-changed", G_CALLBACK(death_rate_scale_moved), NULL);
    /* Pack it in a Frame */
    frame = gtk_frame_new("Death rate");
    gtk_container_add(GTK_CONTAINER(frame), scale);
    /* Add that Frame to the CP box */
    gtk_container_add(GTK_CONTAINER(box), frame);
    /* ALPHA */
    /* Scale bar to set the cell differentiation rate */
    scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, (gdouble)ALPHA_MIN, (gdouble)ALPHA_MAX, (gdouble)ALPHA_STEP);
    /* Set it to its default value */
    gtk_range_set_value(GTK_RANGE(scale), (gfloat)ALPHA);
    g_signal_connect(scale, "value-changed", G_CALLBACK(differentiation_rate_scale_moved), NULL);
    /* Pack it in a Frame */
    frame = gtk_frame_new("Differentiation rate");
    gtk_container_add(GTK_CONTAINER(frame), scale);
    /* Add that Frame to the CP box */
    gtk_container_add(GTK_CONTAINER(box), frame);

    /* Create a Contact Process label and put it with its box in the Notebook */
    label = gtk_label_new("Contact Process");
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), box, label);

    /* Ising Model (IM) Box */
    /* to control the parameters of the Ising model */
    box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    /* TEMPERATURE */
    /* Create a scale bar to set Temperature */
    scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, (gdouble)TEMPERATURE_MIN, (gdouble)TEMPERATURE_MAX, (gdouble)TEMPERATURE_STEP);
    /* Set it to its default value */
    gtk_range_set_value(GTK_RANGE(scale), (gfloat)TEMPERATURE);
    g_signal_connect(scale, "value-changed", G_CALLBACK(temperature_scale_moved), NULL);
    frame = gtk_frame_new("Temperature");
    gtk_container_add(GTK_CONTAINER(frame), scale);
    gtk_container_add(GTK_CONTAINER(box), frame);
    /* EXTERNAL FIELD */
    scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, (gdouble)EXTERNAL_FIELD_MIN, (gdouble)EXTERNAL_FIELD_MAX, (gdouble)EXTERNAL_FIELD_STEP);
    gtk_range_set_value(GTK_RANGE(scale), (gfloat)EXTERNAL_FIELD);
    g_signal_connect(scale, "value-changed", G_CALLBACK(external_field_scale_moved), NULL);
    frame = gtk_frame_new("External field");
    gtk_container_add(GTK_CONTAINER(frame), scale);
    gtk_container_add(GTK_CONTAINER(box), frame);

    /* Add a vertical separator to the parameter grid for order */
    separator = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
    gtk_container_add(GTK_CONTAINER(box), separator);

    /* Add the grid containing the pre-packed radio buttons to control spin coupling */
    frame = gtk_frame_new("Spin interaction");
    gtk_container_add(GTK_CONTAINER(frame), grid);
    gtk_container_add(GTK_CONTAINER(box), frame);

    /* Create an IM label and put it with its box in the Notebook */
    label = gtk_label_new("Ising Model");
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), box, label);

    /* Initial Conditions (IC) Box */
    /* Create a box to hold a bunch of radio buttons representing different initial configurations */
    box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    /* Create buttons and pack them up */
    /* -1- */
    radio = gtk_radio_button_new_with_label(NULL, "Single differentiated site (spin up or down)");
    g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_initial_condition_1), (gpointer)"option 1 selected");
    gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
    /* -2- */
    radio = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(radio), "Single undifferentiated site");
    g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_initial_condition_2), (gpointer)"option 2 selected");
    gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
    /* -3- */
    radio = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(radio), "Single undifferentiated cluster");
    g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_initial_condition_3), (gpointer)"option 3 selected");
    gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
    /* -4- */
    radio = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(radio), "Single differentiated cluster (spin up or down)");
    g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_initial_condition_4), (gpointer)"option 4 selected");
    gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
    /* -5- */
    radio = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(radio), "Fully occupied undifferentiated lattice");
    g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_initial_condition_5), (gpointer)"option 5 selected");
    gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
    /* -6- */
    radio = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(radio), "Fully occupied differentiated lattice");
    g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_initial_condition_6), (gpointer)"option 6 selected");
    gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
    /* Create IC label for Initial Conditions page and put it in the Notebook */
    label = gtk_label_new("Init Lattice");
    frame = gtk_frame_new("Configuration at genesis (t = 0)");
    gtk_frame_set_label_align(GTK_FRAME(frame), 0, 1);
    /* --- Add the page with the frame and label --- */
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), frame, label);
    /* Add IC box to frame to be put on the frame of the third page of the Notebook */
    gtk_container_add(GTK_CONTAINER(frame), box);

    /* Parameters section final touch */
    /* Create a final frame for the whole Notebook */
    GtkWidget *parameters_frame;
    parameters_frame = gtk_frame_new("Parameters");
    gtk_frame_set_label_align(GTK_FRAME(parameters_frame), 0, 1);

    /* Add the notebook to the recently made frame */
    gtk_container_add(GTK_CONTAINER(parameters_frame), notebook);

    /* MAIN WINDOW */
    /* Create a new window, and set its title */
    window = gtk_application_window_new(app);
    gtk_window_set_title(GTK_WINDOW(window), "Contact Process Ising Model");
    gtk_window_set_resizable(GTK_WINDOW(window), FALSE);

    /* Use Gtk grid to pack our widgets in the Main App Window */
    grid = gtk_grid_new();

    /* Add the frame to the Application window's main grid */
    /* Pack it into the window main grid */
    gtk_grid_attach(GTK_GRID(grid), parameters_frame, 0, 0, 5, 3);

    /* PIX BUFFER */
    /* Pixel buffer @ start up and default canvas display */
    pixbuf = gdk_pixbuf_new(GDK_COLORSPACE_RGB, 0, 8, X_SIZE, Y_SIZE);
    image_lattice = gtk_image_new_from_pixbuf(pixbuf);
    paint_a_background(image_lattice);
    /* Place the image on row 7 of our grid spanning 5 columns */
    gtk_grid_attach(GTK_GRID(grid), image_lattice, 0, 7, 5, 1);
    /* Position (0,3) spanning 5 col and 1 row */

    /* Simulation CONTROLS */
    GtkWidget *button, *ctrl_frame, *button_box;

    button_box = gtk_button_box_new(GTK_ORIENTATION_HORIZONTAL);
    ctrl_frame = gtk_frame_new("Simulation Control");
    /* Initialize Lattice */
    button = gtk_button_new_with_label("Init");
    g_signal_connect(button, "clicked", G_CALLBACK(init_lattice), GTK_IMAGE(image_lattice));
    gtk_container_add(GTK_CONTAINER(button_box), button);
    /* Start */
    button = gtk_button_new_with_label("Start");
    g_signal_connect(button, "clicked", G_CALLBACK(on_button_start_simulation), GTK_IMAGE(image_lattice));
    gtk_container_add(GTK_CONTAINER(button_box), button);
    /* Stop */
    button = gtk_button_new_with_label("Stop");
    g_signal_connect(button, "clicked", G_CALLBACK(on_button_stop_simulation), NULL);
    gtk_container_add(GTK_CONTAINER(button_box), button);
    /* About (to be fair to the user) */
    button = gtk_button_new_with_label("About");
    g_signal_connect(button, "clicked", G_CALLBACK(on_button_show_about), NULL);
    gtk_container_add(GTK_CONTAINER(button_box), button);
    /* Quit */
    button = gtk_button_new_with_label("Quit");
    g_signal_connect(button, "clicked", G_CALLBACK(on_button_stop_simulation), NULL);
    g_signal_connect_swapped(button, "clicked", G_CALLBACK(gtk_widget_destroy), window);
    gtk_container_add(GTK_CONTAINER(button_box), button);

    /* Add all the buttons in their specialized frame */
    gtk_container_add(GTK_CONTAINER(ctrl_frame), button_box);
    /* Place the frame on row 8 of our grid spanning 5 columns */
    gtk_grid_attach(GTK_GRID(grid), ctrl_frame, 0, 8, 5, 1);

    /* FRAME_SAMPLE_RATE */
    /* Scale bar to set display rate */
    scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, 1, 100, 10);
    gtk_range_set_value(GTK_RANGE(scale), (gfloat)SAMPLE_RATE);
    g_signal_connect(scale, "value-changed", G_CALLBACK(display_rate_scale_moved), NULL);
    /* Pack it in a Frame */
    frame = gtk_frame_new("Display rate");
    gtk_container_add(GTK_CONTAINER(frame), scale);
    /* Add that Frame to the CP box */
    gtk_grid_attach(GTK_GRID(grid), frame, 0, 9, 5, 1);

    /* Pack the main grid into the window */
    gtk_container_add(GTK_CONTAINER(window), grid);
    /* Show the window and all widgets in it */
    gtk_widget_show_all(window);
}


/* Implementation of put_pixel function. Code retrieved from:
   https://developer.gnome.org/gdk-pixbuf/stable/gdk-pixbuf-The-GdkPixbuf-Structure.html */
void put_pixel(GdkPixbuf *pixbuf, int x, int y,
               guchar red, guchar green, guchar blue, guchar alpha)
{
    guchar *pixels, *p;
    int rowstride, num_channels;
    num_channels = gdk_pixbuf_get_n_channels(pixbuf);
    rowstride = gdk_pixbuf_get_rowstride(pixbuf);
    pixels = gdk_pixbuf_get_pixels(pixbuf);
    p = pixels + y * rowstride + x * num_channels;
    p[0] = red;
    p[1] = green;
    p[2] = blue;
    p[3] = alpha;
}


/* Creates a pixel buffer and paints an image to display as the default canvas */
static void paint_a_background(gpointer data)
{
    GdkPixbuf *p;
    p = gdk_pixbuf_new(GDK_COLORSPACE_RGB, 0, 8, X_SIZE, Y_SIZE);
    /* Paint a background canvas for the startup image */
    int x, y;
    for (x = 0; x < X_SIZE; x++)
    {
        for (y = 0; y < Y_SIZE; y++)
        {
            put_pixel(p, x, y, (guchar)x, (guchar)y, (guchar)x, 255);
        }
    }
    gtk_image_set_from_pixbuf(GTK_IMAGE(data), GDK_PIXBUF(p));
    g_object_unref(p);
}


/* Function that paints the pixel buffer with the simulation data */
/* The states are:
 0: vacant
-1: spin down
+1: spin up
 2: undifferentiated */
static void paint_lattice(gpointer data)
{
    // Create a Gdk pixbuffer to paint configurations
    GdkPixbuf *p;
    p = gdk_pixbuf_new(GDK_COLORSPACE_RGB, 0, 8, X_SIZE, Y_SIZE);
    /* Paint lattice configuration to a pixel buffer */
    int x, y;
    for (x = 0; x < X_SIZE; x++)
    {
        for (y = 0; y < Y_SIZE; y++)
        {
            switch (s.lattice_configuration[x][y])
            {
                case 0: /* Empty (vacant) site (black) */
                    put_pixel(p, x, y, 0, 0, 0, 255);
                    break;
                case -1: /* Spin down (occupied) site (green) */
                    put_pixel(p, x, y, 0, 255, 0, 255);
                    break;
                case 1: /* Spin up (occupied) site (magenta) */
                    put_pixel(p, x, y, 255, 0, 255, 255);
                    break;
                case 2: /* Undifferentiated (occupied) site (white) */
                    put_pixel(p, x, y, 255, 255, 255, 255);
                    break;
            }
        }
    }
    gtk_image_set_from_pixbuf(GTK_IMAGE(data), GDK_PIXBUF(p));
    g_object_unref(p);
}


/* Main function spanning a Gtk Application object */
int main(int argc, char **argv)
{
    GtkApplication *app;
    int status;
    app = gtk_application_new("keymer.lab.contact_process_ising_model", G_APPLICATION_DEFAULT_FLAGS);
    g_signal_connect(app, "activate", G_CALLBACK(activate), NULL);
    status = g_application_run(G_APPLICATION(app), argc, argv);
    g_object_unref(app);
    free(s.transition_probs);
    return status;
}