
	/***********************************************************************
	 * Run a simple N-body gravitational simulation.
	 *
	 * This code is intended to let me explore facets of different
	 * techniques, not do any serious work.  Not optimized in any
	 * way, not documented well, likely to change on the fly.
	 *
	 * MWR 12/22/2008
	 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#undef DEBUG

#define NBUF                1024
	/* max length of lines in input file(s) */
#define LINELEN             1024
	/* max length of names of objects */
#define NAMELEN               30

	/* any line starting with this in input file is ignored */
#define COMMENT_CHAR         '#'
	/* maximum number of bodies we can handle */
#define MAX_NBODY             20
	/* maximum number of intergration techniques we can handle */
#define MAX_TECHNIQUE         10

	/* physical constants  */
	/* constant of universal gravitation, in MKS units */
#define G     6.673E-11
	/* number of seconds in one year */
#define SEC_IN_YEAR     (86400*365.25)


	/* pointer to a function returning an integer */
typedef int (*PFI)(double, double *);

/* structures */
typedef struct s_VECTOR {
	double val[3];
} VECTOR;


typedef struct s_body {
	int index;
	char name[NAMELEN];
	double mass;
	double pos[3];
	double vel[3];
} BODY;

typedef struct s_technique {
	char name[NAMELEN];
	PFI func;
} TECHNIQUE;



/* these macros can make code easier to read */
	/*
	 * set the given VECTOR elements to have the difference in
	 * positions between bodies i and j, in the sense
	 *     g_body_array[i].pos[0] - g_body_array[j].pos[0]
	 * etc.
	 */
#define SET_POS_DIFF(vec, i, j) \
      vec.val[0] = g_body_array[i].pos[0] - g_body_array[j].pos[0] ; \
      vec.val[1] = g_body_array[i].pos[1] - g_body_array[j].pos[1] ; \
      vec.val[2] = g_body_array[i].pos[2] - g_body_array[j].pos[2] ;

	/* the magnitude of a VECTOR */
#define VECTOR_MAG(vec) \
      ( sqrt(vec.val[0]*vec.val[0] + vec.val[1]*vec.val[1] + \
           vec.val[2]*vec.val[2]) )

	/* the square of the magnitude of a VECTOR */
#define VECTOR_MAG_SQUARED(vec) \
      (vec.val[0]*vec.val[0] + vec.val[1]*vec.val[1] + \
           vec.val[2]*vec.val[2])

	/* the SQUARE of distance between two bodies */
#define DIST_SQUARED(i, j) \
     ( ((g_body_array[i].pos[0]-g_body_array[j].pos[0])*\
        (g_body_array[i].pos[0]-g_body_array[j].pos[0])) + \
       ((g_body_array[i].pos[1]-g_body_array[j].pos[1])*\
        (g_body_array[i].pos[1]-g_body_array[j].pos[1])) + \
       ((g_body_array[i].pos[2]-g_body_array[j].pos[2])*\
        (g_body_array[i].pos[2]-g_body_array[j].pos[2])) )

	/* the difference of two VECTORs */
#define VECTOR_SUBTRACT(a, b, difference) \
       difference.val[0] = a.val[0] - b.val[0]; \
       difference.val[1] = a.val[1] - b.val[1]; \
       difference.val[2] = a.val[2] - b.val[2];

	/* the cross product of two VECTORs */
#define VECTOR_CROSS(a, b, product) \
       product.val[0] = a.val[1]*b.val[2] - a.val[2]*b.val[1] ; \
       product.val[1] = a.val[2]*b.val[0] - a.val[0]*b.val[2] ; \
       product.val[2] = a.val[0]*b.val[1] - a.val[1]*b.val[0] ;



/* global variables */
	/* if set to 1, print messages as we execute */
int g_verbose_flag = 0;
	/* FILE pointer for output */
FILE *g_outfile_fp;
	/*
	 * if set to 1, modify all initial velocities to force
	 *     COM to be motionless?
	 */
int g_recenter_flag = 0;
	/* how frequently should we print positions to output? (seconds) */
double g_print_interval = 0;

int g_body_number = 0;
BODY g_body_array[MAX_NBODY];

	/* these are the diffent methods for performing integration */
int g_technique_number = 0;
TECHNIQUE g_technique_array[MAX_TECHNIQUE];
PFI g_integration_func;

double g_duration_sec = 0;
double g_timestep_sec = 0;

/* functions */
static int read_input(char *filename);
static double calc_total_ke(void);
static double calc_total_gpe(void);
static int tech_euler_1(double suggested_timestep,
                        double *actual_timestep);
static int tech_euler_2(double suggested_timestep,
                        double *actual_timestep);
static int tech_euler_1a(double suggested_timestep,
                        double *actual_timestep);
static int tech_rk4(double suggested_timestep,
                        double *actual_timestep);
static int print_positions(FILE *outfile_fp, double time_sec);
static int recenter_velocities(void);
static int calc_total_angmom(VECTOR *result);
static int calc_grav_forces(int nbody, BODY *body_array, VECTOR *forces);
static void copy_bodies(int num_bodies, BODY *from_array, BODY *to_array);



	/*
	 * usage:
	 *              nbody  input=inputfile [verbose=]
	 *       where
	 *                inputfile        is ASCII text file with
	 *                                    description of each body,
	 *                                    plus some controlling parameters
	 *
	 *                verbose          if present, set verbosity level to 1
	 *
	 *                verbose=N        if present, set verbosity level to N
	 *
	 *       Output is sent to stdout for now, but probably
	 *       will eventually set name of output file(s) in
	 *       the inputfile somewhere.
	 *
	 */

int
main
	(
	int argc,
	char *argv[]
	)
{
	int i;
	char inputfile_name[NBUF];
	double initial_ke, initial_gpe, initial_tot_e;
	double final_ke, final_gpe, final_tot_e;
	double time_sec;
	double printing_clock;
	VECTOR initial_ang_mom, final_ang_mom;

#ifdef TEST
	/* test code, don't use in production */
	{
		VECTOR diff;

		g_body_array[0].pos[0] = 0.0;
		g_body_array[0].pos[1] = 0.0;
		g_body_array[0].pos[2] = 0.0;
		g_body_array[1].pos[0] = 10.0;
		g_body_array[1].pos[1] = -10.0;
		g_body_array[1].pos[2] = 3.0;

		SET_POS_DIFF(diff, 0, 1);
		printf("  diff VECTOR is %lf %lf %lf \n",
				diff.val[0], diff.val[1], diff.val[2]);

		exit(0);
	}
#endif




	/* initialize stuff */
	inputfile_name[0] = '\0';
#ifdef DEBUG
	g_verbose_flag = 1;
#endif
	for (i = 0; i < MAX_NBODY; i++) {
		g_body_array[i].index = -1;
		strcpy(g_body_array[i].name, "");
		g_body_array[i].mass = -1;
	}
	g_outfile_fp = (FILE *) NULL;
	/* initialize the known techniques for integration */
		/* Euler first order method "E1" */
	strcpy(g_technique_array[g_technique_number].name, "E1");
	g_technique_array[g_technique_number++].func = tech_euler_1;
		/*
		 * Euler first order method "E1a" --
		 *   uses new velocities to advance old positions
		 */
	strcpy(g_technique_array[g_technique_number].name, "E1a");
	g_technique_array[g_technique_number++].func = tech_euler_1a;
		/* Euler second order method "E2" (same as Heun's method) */
	strcpy(g_technique_array[g_technique_number].name, "E2");
	g_technique_array[g_technique_number++].func = tech_euler_2;
		/* Runge-Kutta fourth order method "RK4" */
	strcpy(g_technique_array[g_technique_number].name, "RK4");
	g_technique_array[g_technique_number++].func = tech_rk4;


	/* parse arguments */
	if (argc < 2) {
		fprintf(stderr, "usage: nbody input=  [verbose=] \n");
		exit(1);
	}
	for (i = 1; i < argc; i++) {
		if (g_verbose_flag > 0) {
			printf(" arg %d is ..%s.. \n", i, argv[i]);
		}

		/* read input file name */
		if (strncmp(argv[i], "input=", 6) == 0) {
			if (strlen(argv[i]) >= NBUF) {
				fprintf(stderr, "inputfile name ..%s.. is too long, max %d chars. \n",
						argv[i], (int) strlen(argv[i]));
				exit(1);
			}
			if (sscanf(argv[i] + 6, "%s", inputfile_name) != 1) {
				fprintf(stderr, "can't read inputfile_name from ..%s.. \n",
							argv[i] + 6);
				exit(1);
			}
			if (g_verbose_flag > 0) {
				printf(" read input file name ..%s.. \n", inputfile_name);
			}
		}
		if (strcmp(argv[i], "verbose") == 0) {
			g_verbose_flag = 1;
			printf(" set verbose level to ..%d.. \n", g_verbose_flag);
		}
		if (strncmp(argv[i], "verbose=", 8) == 0) {
			if (sscanf(argv[i] + 8, "%d", &g_verbose_flag) != 1) {
				fprintf(stderr, "can't read verbosity level from ..%s.. \n",
							argv[i] + 8);
				exit(1);
			}
			if (g_verbose_flag > 0) {
				printf(" set verbose level to ..%d.. \n", g_verbose_flag);
			}
		}


	}

	/* check to make sure we read all required args */
	if (inputfile_name[0] == '\0') {
		fprintf(stderr, "no input= argument provided ?! \n");
		exit(1);
	}



	/* read all input information */
	if (read_input(inputfile_name) != 0) {
		fprintf(stderr, "read_input fails \n");
		exit(1);
	}

	/*
	 * if desired, add a constant velocity to all bodies
	 *    so that the center of mass of the system
	 *    is motionless
	 */
	if (g_recenter_flag == 1) {
		if (recenter_velocities() != 0) {
			fprintf(stderr, "recenter_velocities returns with error \n");
			exit(1);
		}
	}


	/* compute initial quantities */
	initial_ke = calc_total_ke();
	if (g_verbose_flag > 0) {
		printf(" initial KE  %12.5e \n", initial_ke);
	}
	initial_gpe = calc_total_gpe();
	if (g_verbose_flag > 0) {
		printf(" initial GPE %12.5e \n", initial_gpe);
	}
	initial_tot_e = initial_ke + initial_gpe;
	if (g_verbose_flag > 0) {
		printf(" initial E   %12.5e \n", initial_tot_e);
	}
	if (calc_total_angmom(&initial_ang_mom) != 0) {
		fprintf(stderr, "calc_total_angmom fails for initial \n");
		exit(1);
	}
	if (g_verbose_flag > 0) {
		printf(" initial angmom   %12.5e %12.5e %12.5e \n",
					initial_ang_mom.val[0],
					initial_ang_mom.val[1],
					initial_ang_mom.val[2]);
	}

	/* print initial positions, etc. */
	time_sec = 0;
	if (print_positions(g_outfile_fp, time_sec) != 0) {
		fprintf(stderr, "print_positions fails at time %12.5e \n",
				time_sec);
		exit(1);
	}

	/* enter main loop */
	printing_clock = 0.0;
	while (time_sec < g_duration_sec) {
		int retval;
		double suggested_step, actual_step;

		if (g_verbose_flag > 1) {
			printf("  about to enter integ func   t = %9.4e \n",
						time_sec);
		}

		/*
		 * attempt to move forward by one step.
		 *   This routine is passed a "suggested" timestep,
		 *   and returns the timestep which was actually
		 *   used.  The routine may decide that the given
		 *   timestep is too large, and it may change
		 *   it in either direction.
		 */
		suggested_step = g_timestep_sec;

		/* call update function here */
		retval = (*(g_integration_func))(suggested_step,
					&actual_step);
		if (retval != 0) {
			fprintf(stderr, "integration func fails \n");
			return(1);
		}
		g_timestep_sec = actual_step;

		printing_clock += actual_step;
		if (printing_clock >= g_print_interval) {
			/* print the new positions, etc. */
			if (print_positions(g_outfile_fp, time_sec) != 0) {
				fprintf(stderr, "print_positions fails at time %12.5e \n",
						time_sec);
				exit(1);
			}
			printing_clock = 0.0;
		}

		time_sec += g_timestep_sec;
	}



	/* compute final quantities */
	final_ke = calc_total_ke();
	if (g_verbose_flag > 0) {
		printf(" final   KE  %12.5e \n", final_ke);
	}
	final_gpe = calc_total_gpe();
	if (g_verbose_flag > 0) {
		printf(" final   GPE %12.5e \n", final_gpe);
	}
	final_tot_e = final_ke + final_gpe;
	if (g_verbose_flag > 0) {
		printf(" final   E   %12.5e \n", final_tot_e);
	}
	if (calc_total_angmom(&final_ang_mom) != 0) {
		fprintf(stderr, "calc_total_angmom fails for final \n");
		exit(1);
	}
	if (g_verbose_flag > 0) {
		printf(" final   angmom   %12.5e %12.5e %12.5e \n",
					final_ang_mom.val[0],
					final_ang_mom.val[1],
					final_ang_mom.val[2]);
	}


	/* compute change in quantities */
	{
		double delta_e, fraction_delta_e;
		double start_angmom_mag, end_angmom_mag, delta_angmom_mag;
		double fraction_delta_angmom_mag;

		delta_e = final_tot_e - initial_tot_e;
		fraction_delta_e = delta_e / fabs(initial_tot_e);
		if (g_verbose_flag > 0) {
			printf(" delta_E  %12.5e    frac %12.5e \n",
					delta_e, fraction_delta_e);
		}

		start_angmom_mag = VECTOR_MAG(initial_ang_mom);
		end_angmom_mag = VECTOR_MAG(final_ang_mom);
		delta_angmom_mag = end_angmom_mag - start_angmom_mag;
		fraction_delta_angmom_mag = delta_angmom_mag / fabs(start_angmom_mag);
		if (g_verbose_flag > 0) {
			printf(" delta_angmom  %12.5e    frac %12.5e \n",
					delta_angmom_mag, fraction_delta_angmom_mag);
		}

	}



	/* prepare to exit */
	if (g_verbose_flag > 0) {
		printf(" all done \n");
	}
	fclose(g_outfile_fp);

	exit(0);
}



	/********************************************************************
	 * PROCEDURE: read_input
	 *
	 * DESCRIPTION: Given a file name, we open the file and read
	 *              information from it on the number of objects,
	 *              their initial positions, parameter values, etc.
	 *
	 *              We may make several passes through the file;
	 *              inefficient, but who cares?
	 *
	 * RETURNS:
	 *              0        if all goes well
	 *              1        if an error occurs
	 */

static int
read_input
	(
	char *filename     /* I: name of ASCII text file with input information */
	)
{
	char line[LINELEN];
	int j, nbody;
	FILE *fp;

	/* open file (and see if it exists) */
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "read_input: can't open file ..%s.. \n", filename);
		return(1);
	}

	/* read the number of bodies in the simulation: for example,
	 *    nbody    10
	 */
	nbody = -1;
	rewind(fp);
	while (fgets(line, LINELEN, fp) != NULL) {
		if (g_verbose_flag > 1) {
			printf("next line in input file is:\n..%s..\n", line);
		}
		if (line[0] == COMMENT_CHAR) {
			if (g_verbose_flag > 1) {
				printf(" skipping comment line \n");
			}
			continue;
		}
		if (line[0] == '\n') {
			if (g_verbose_flag > 1) {
				printf(" skipping empty line \n");
			}
			continue;
		}

		if (strncmp(line, "nbody", 5) == 0) {
			if (sscanf(line + 5, "%d", &nbody) != 1) {
				fprintf(stderr, "bad nbody value in line ..%s.. \n", line);
				return(1);
			}
			if (nbody < 1) {
				fprintf(stderr, "nbody must be >= 1 \n");
				return(1);
			}
			if (g_verbose_flag > 0) {
				printf(" nbody is %d \n", nbody);
			}
		}
	}
	g_body_number = nbody;
	if ((g_body_number < 1) || (g_body_number >= MAX_NBODY)) {
		fprintf(stderr, "nbody %d is < 1 or > MAX_NBODY %d \n",
				g_body_number, MAX_NBODY);
		return(1);
	}

	/* read the duration of the simulation, and convert from
	 *    years into seconds
	 *    For example,
	 *
	 *    duration    10.0
	 *
	 * would mean "10 years".  We will set the global variable
	 *    "g_duration_sec" to about 3 x 10^(8)
	 *
	 */
	g_duration_sec = -1;
	rewind(fp);
	while (fgets(line, LINELEN, fp) != NULL) {
		double duration_years;

		if (g_verbose_flag > 1) {
			printf("next line in input file is:\n..%s..\n", line);
		}
		if (line[0] == COMMENT_CHAR) {
			if (g_verbose_flag > 1) {
				printf(" skipping comment line \n");
			}
			continue;
		}
		if (line[0] == '\n') {
			if (g_verbose_flag > 1) {
				printf(" skipping empty line \n");
			}
			continue;
		}

		if (strncmp(line, "duration", 8) == 0) {
			if (sscanf(line + 8, "%lf", &duration_years) != 1) {
				fprintf(stderr, "bad duration value in line ..%s.. \n", line);
				return(1);
			}
			if (duration_years < 0) {
				fprintf(stderr, "duration must be >= 0 \n");
				return(1);
			}
			g_duration_sec = duration_years*SEC_IN_YEAR;
			if (g_verbose_flag > 0) {
				printf(" duration is %9.4e year = %9.4e sec \n",
								duration_years, g_duration_sec);
			}
		}
	}
	if (g_duration_sec < 0) {
		fprintf(stderr, "g_duration_sec %lf is < 0 \n", g_duration_sec);
		return(1);
	}

	/* read the (initial) timestep of the simulation, in seconds.
	 *    For example,
	 *
	 *    timestep    10.0
	 *
	 * would mean "10 seconds".  We will set the global variable
	 *    "g_timestep_sec" to 10.0
	 *
	 */
	g_timestep_sec = -1;
	rewind(fp);
	while (fgets(line, LINELEN, fp) != NULL) {
		double timestep;

		if (g_verbose_flag > 1) {
			printf("next line in input file is:\n..%s..\n", line);
		}
		if (line[0] == COMMENT_CHAR) {
			if (g_verbose_flag > 1) {
				printf(" skipping comment line \n");
			}
			continue;
		}
		if (line[0] == '\n') {
			if (g_verbose_flag > 1) {
				printf(" skipping empty line \n");
			}
			continue;
		}

		if (strncmp(line, "timestep", 8) == 0) {
			if (sscanf(line + 8, "%lf", &timestep) != 1) {
				fprintf(stderr, "bad timestep value in line ..%s.. \n", line);
				return(1);
			}
			if (timestep <= 0) {
				fprintf(stderr, "timestep must be > 0 \n");
				return(1);
			}
			g_timestep_sec = timestep;
			if (g_verbose_flag > 0) {
				printf(" timestep is %9.4e sec \n", g_timestep_sec);
			}
		}
	}
	if (g_timestep_sec <= 0) {
		fprintf(stderr, "g_timestep_sec %lf is < 0 \n", g_timestep_sec);
		return(1);
	}


	/* read the interval we should allow to elapse between printing
	 *    positions to the output file.  The units are seconds.
	 */
	g_print_interval = -1;
	rewind(fp);
	while (fgets(line, LINELEN, fp) != NULL) {
		double interval;

		if (g_verbose_flag > 1) {
			printf("next line in input file is:\n..%s..\n", line);
		}
		if (line[0] == COMMENT_CHAR) {
			if (g_verbose_flag > 1) {
				printf(" skipping comment line \n");
			}
			continue;
		}
		if (line[0] == '\n') {
			if (g_verbose_flag > 1) {
				printf(" skipping empty line \n");
			}
			continue;
		}

		if (strncmp(line, "print_interval", 14) == 0) {
			if (sscanf(line + 14, "%lf", &interval) != 1) {
				fprintf(stderr, "bad print_interval value in line ..%s.. \n", line);
				return(1);
			}
			if (interval <= 0) {
				fprintf(stderr, "print_interval must be > 0 \n");
				return(1);
			}
			g_print_interval = interval;
			if (g_verbose_flag > 0) {
				printf(" print_interval is %9.4e sec \n", g_print_interval);
			}
		}
	}
	if (g_print_interval <= 0) {
		fprintf(stderr, "g_print_interval %lf is < 0 \n", g_print_interval);
		return(1);
	}


	/* look for a line which tells us -- should we leave initial
	 *    velocities as provided by the user, or should we add a
	 *    constant value to force the center-of-mass velocity
	 *    to be zero?
	 */
	g_recenter_flag = -1;
	rewind(fp);
	while (fgets(line, LINELEN, fp) != NULL) {
		char recenter_value[NAMELEN];

		if (g_verbose_flag > 1) {
			printf("next line in input file is:\n..%s..\n", line);
		}
		if (line[0] == COMMENT_CHAR) {
			if (g_verbose_flag > 1) {
				printf(" skipping comment line \n");
			}
			continue;
		}
		if (line[0] == '\n') {
			if (g_verbose_flag > 1) {
				printf(" skipping empty line \n");
			}
			continue;
		}

		if (strncmp(line, "recenter", 8) == 0) {
			if (sscanf(line + 8, "%s", recenter_value) != 1) {
				fprintf(stderr, "bad recenter value in line ..%s.. \n", line);
				return(1);
			}
			if (strcmp(recenter_value, "yes") == 0) {
				g_recenter_flag = 1;
			}
			if (strcmp(recenter_value, "no") == 0) {
				g_recenter_flag = 0;
			}
			if (g_verbose_flag > 0) {
				printf(" recenter_flag is %d \n", g_recenter_flag);
			}
		}
	}
	if (g_recenter_flag < 0) {
		fprintf(stderr, "g_recenter_flag not set?! \n");
		return(1);
	}


	/* read the technique we should use to perform the
	 *    numerical integration during each timestep.
	 *    The input file contains a code, which this
	 *    section of the program must match to one of the
	 *    existing functions.  For example,
	 *
	 *    technique   E1
	 *
	 *    means "use the function associated with
	 *    the code 'E1'."   Near the top of the main()
	 *    routine is a block of code which assigns
	 *    functions to codes.
	 *
	 */
	rewind(fp);
	while (fgets(line, LINELEN, fp) != NULL) {
		char technique[NAMELEN];

		if (g_verbose_flag > 1) {
			printf("next line in input file is:\n..%s..\n", line);
		}
		if (line[0] == COMMENT_CHAR) {
			if (g_verbose_flag > 1) {
				printf(" skipping comment line \n");
			}
			continue;
		}
		if (line[0] == '\n') {
			if (g_verbose_flag > 1) {
				printf(" skipping empty line \n");
			}
			continue;
		}

		if (strncmp(line, "technique", 9) == 0) {
			if (sscanf(line + 9, "%s", technique) != 1) {
				fprintf(stderr, "bad technique value in line ..%s.. \n", line);
				return(1);
			}
			/*
			 * make sure that the technique matches one
			 * of the entries in the
			 */
			g_integration_func = (PFI) NULL;
			for (j = 0; j < g_technique_number; j++) {
				if (strcmp(technique, g_technique_array[j].name) == 0) {
					g_integration_func = g_technique_array[j].func;
					break;
				}
			}
			if (g_integration_func == NULL) {
				fprintf(stderr, "can't find technique %s\n", technique);
				return(1);
			}
			if (g_verbose_flag > 0) {
				printf(" will use the %s technique \n", technique);
			}
		}
	}
	if (g_integration_func == NULL) {
		fprintf(stderr, "no technique provided in %s ?! \n",
					filename);
		return(1);
	}


	/*
	 * Read the name of the file into which to write the
	 *    output.  If the name is "-", then we'll
	 *    write to stdout.
	 */
	rewind(fp);
	while (fgets(line, LINELEN, fp) != NULL) {
		char outfile_name[NAMELEN];

		if (g_verbose_flag > 1) {
			printf("next line in input file is:\n..%s..\n", line);
		}
		if (line[0] == COMMENT_CHAR) {
			if (g_verbose_flag > 1) {
				printf(" skipping comment line \n");
			}
			continue;
		}
		if (line[0] == '\n') {
			if (g_verbose_flag > 1) {
				printf(" skipping empty line \n");
			}
			continue;
		}

		if (strncmp(line, "outfile", 7) == 0) {
			if (sscanf(line + 7, "%s", outfile_name) != 1) {
				fprintf(stderr, "bad output filename in line ..%s.. \n", line);
				return(1);
			}
			/*
			 * open the file for output, set the global variable
			 *    g_outfile_fp to point to the opened file
			 */
			if (strcmp(outfile_name, "-") == 0) {
				g_outfile_fp = stdout;
			}
			else {
				if ((g_outfile_fp = fopen(outfile_name, "w")) == NULL) {
					fprintf(stderr, "can't open file %s for output \n",
									outfile_name);
					return(1);
				}
			}
		}
	}






	/* Now, read the initial information for each body, one at a time.
	 *    Make sure that the number of bodies with information
	 *    is the same as the 'nbody' line indicated.
	 */
	rewind(fp);
	while (fgets(line, LINELEN, fp) != NULL) {
		if (g_verbose_flag > 1) {
			printf("next line in input file is:\n..%s..\n", line);
		}
		if (line[0] == COMMENT_CHAR) {
			if (g_verbose_flag > 1) {
				printf(" skipping comment line \n");
			}
			continue;
		}
		if (line[0] == '\n') {
			if (g_verbose_flag > 1) {
				printf(" skipping empty line \n");
			}
			continue;
		}


		if (strncmp(line, "body ", 5) == 0) {

			int this_index;
			char this_name[LINELEN];
			double this_mass;
			double this_px, this_py, this_pz;
			double this_vx, this_vy, this_vz;

			if (g_verbose_flag > 1) {
				printf(" about to scan for a body's initial info \n");
			}

				if (sscanf(line + 4, " %d %s %lf %lf %lf %lf %lf %lf %lf",
								&this_index, this_name, &this_mass,
								&this_px, &this_py, &this_pz,
								&this_vx, &this_vy, &this_vz) != 9) {
					fprintf(stderr, "bad body in line ..%s.. \n", line);
					return(1);
				}
				if (strlen(this_name) >= NAMELEN) {
					fprintf(stderr, "body with name ..%s.. has name longer than %d chars \n",
								this_name, NAMELEN);
					return(1);
				}
				if ((this_index < 0) || (this_index >= g_body_number)) {
					fprintf(stderr, "body with name ..%s.. has invalid index %d \n",
								this_name, this_index);
					return(1);
				}
				if (this_mass <= 0) {
					fprintf(stderr, "body with name ..%s.. has invalid mass %lf \n",
								this_name, this_mass);
					return(1);
				}
				/* now try to set entries in the appropriate g_body_array[] */
				if (g_body_array[this_index].index != -1) {
					fprintf(stderr, "body with name ..%s.. has repeated index %d \n",
								this_name, this_index);
					return(1);
				}
				else {
					g_body_array[this_index].index = this_index;
					strcpy(g_body_array[this_index].name, this_name);
					g_body_array[this_index].mass = this_mass;
					g_body_array[this_index].pos[0] = this_px;
					g_body_array[this_index].pos[1] = this_py;
					g_body_array[this_index].pos[2] = this_pz;
					g_body_array[this_index].vel[0] = this_vx;
					g_body_array[this_index].vel[1] = this_vy;
					g_body_array[this_index].vel[2] = this_vz;
				}

				if (g_verbose_flag > 0) {
					printf(" nbody %d has \n", this_index);
					printf("      name ..%s.. \n",
									g_body_array[this_index].name);
					printf("      mass ..%12.5le.. \n",
									g_body_array[this_index].mass);
					printf("      pos   %12.5le %12.5le %12.5le  \n",
									g_body_array[this_index].pos[0],
									g_body_array[this_index].pos[1],
									g_body_array[this_index].pos[2]);
					printf("      vel   %12.5le %12.5le %12.5le  \n",
									g_body_array[this_index].vel[0],
									g_body_array[this_index].vel[1],
									g_body_array[this_index].vel[2]);
				}


		}
	}

	/* verify that the input file supplied all required information */
	if (g_outfile_fp == NULL) {
		fprintf(stderr, "read_input: no valid outfile supplied \n");
		return(1);
	}


	fclose(fp);


	return(0);
}



	/********************************************************************
	 * PROCEDURE: calc_total_ke
	 *
	 * DESCRIPTION: Compute the total kinetic energy of all bodies.
	 *              The result will have units of Joules.
	 *
	 * RETURNS:
	 *              KE       in Joules
	 */

static double
calc_total_ke
	(
	void
	)
{
	int i;
	double v;
	double ke;
	double total_ke;

	total_ke = 0;

	for (i = 0; i < g_body_number; i++) {
		v = sqrt( g_body_array[i].vel[0]*g_body_array[i].vel[0]
		        + g_body_array[i].vel[1]*g_body_array[i].vel[1]
		        + g_body_array[i].vel[2]*g_body_array[i].vel[2] );
		ke = 0.5 * g_body_array[i].mass * (v*v);
		total_ke += ke;
		if (g_verbose_flag > 1) {
			printf(" calc_total_ke: body %5d has KE %9.4e  tot %9.4e \n",
						i, ke, total_ke);
		}
	}

	return(total_ke);
}



	/********************************************************************
	 * PROCEDURE: calc_total_gpe
	 *
	 * DESCRIPTION: Compute the total GPE of all bodies.
	 *              The result will have units of Joules.
	 *
	 * RETURNS:
	 *              KE       in Joules
	 */

static double
calc_total_gpe
	(
	void
	)
{
	int i, j;
	double d, d_squared;
	double gpe;
	double total_gpe;

	total_gpe = 0;

	for (i = 0; i < g_body_number - 1; i++) {
		for (j = i + 1; j < g_body_number; j++) {

			d_squared = DIST_SQUARED(i, j);
			d = sqrt(d_squared);
			if (d <= 0) {
				fprintf(stderr, "calc_total_gpe: distance %9.4le is invalid ?! \n",
						d);
				exit(1);
			}

			gpe = 0.0 - ((G*g_body_array[i].mass*g_body_array[j].mass) / d);

			total_gpe += gpe;
			if (g_verbose_flag > 1) {
				printf(" calc_total_gpe: bodies %5d %5d have GPE %9.4e  tot %9.4e \n",
							i, j, gpe, total_gpe);
			}

		}

	}

	return(total_gpe);
}



	/********************************************************************
	 * PROCEDURE: tech_euler_1
	 *
	 * DESCRIPTION: Advance particles one timestep
	 *              using Euler's method to first order.
	 *
	 *              This method does not (yet) attempt to change
	 *              the timestep.
	 *
	 * RETURNS:
	 *              0       if all goes well
	 *              1       if an error occurs
	 */

static int
tech_euler_1
	(
	double suggested_timestep,    /* I: suggested timestep (sec) */
	double *actual_timestep       /* O: actual timestep used (sec) */
	)
{
	int i, j, c;
	double timestep;
	VECTOR new_pos[MAX_NBODY];
	VECTOR new_vel[MAX_NBODY];
	VECTOR cur_force[MAX_NBODY][MAX_NBODY];

	timestep = suggested_timestep;

	/*
	 * compute the forces between all objects
	 *   using the current positions
	 */
	for (i = 0; i < g_body_number; i++) {
		for (j = i; j < g_body_number; j++) {
			VECTOR dist, dist_frac;
			double dist_tot;
			double dist_squared;
			double top, force_mag;

			if (i == j) {
				cur_force[i][j].val[0] = 0.0;
				cur_force[i][j].val[1] = 0.0;
				cur_force[i][j].val[2] = 0.0;
				continue;
			}

			SET_POS_DIFF(dist, i, j);
			dist_squared = VECTOR_MAG_SQUARED(dist);
			dist_tot = sqrt(dist_squared);
			if (dist_tot <= 0.0) {
				fprintf(stderr,
						" tech_euler_1: distance of %le -- quitting \n", dist_tot);
				return(1);
			}

			/* this is where we could add a softening parameter */
			top = G * g_body_array[i].mass * g_body_array[j].mass;
			force_mag = top / dist_squared;

			for (c = 0; c < 3; c++) {
				dist_frac.val[c] = dist.val[c] / dist_tot;
				cur_force[i][j].val[c] = 0.0 - (force_mag * dist_frac.val[c]);
				cur_force[j][i].val[c] = 0.0 - cur_force[i][j].val[c];
			}
		}
	}


	/* do the work here */
	for (i = 0; i < g_body_number; i++) {
		VECTOR tot_force;
		VECTOR tot_accel;

		/*
		 * we will compute a new position and velocity
		 *   for each object in turn, placing the new
		 *   values into the 'new_pos[]' and 'new_vel[]'
		 *   arrays.  Afterwards,
		 *   we'll copy the new values back into the
		 *   original g_body[] array.
		 */


		/* use the current velocity to compute new positions */
		for (c = 0; c < 3; c++) {
			new_pos[i].val[c] = g_body_array[i].pos[c]
							+ g_body_array[i].vel[c]*timestep;
		}

		/* use the current forces to compute new velocities */
		for (c = 0; c < 3; c++) {

			/*
			 * compute the total force on this body
			 *    in each direction
			 */
			tot_force.val[c] = 0.0;
			for (j = 0; j < g_body_number; j++) {
				tot_force.val[c] += cur_force[i][j].val[c];
			}

			/* from total force, compute total acceleration */
			tot_accel.val[c] = tot_force.val[c] / g_body_array[i].mass;

			/* use the acceleration to compute new velocity */
			new_vel[i].val[c] = g_body_array[i].vel[c]
							+ tot_accel.val[c]*timestep;
		}

	}


	/*
	 * having computed the new positions and velocities,
	 *    we now copy those values into the g_body_array[]
	 */
	for (i = 0; i < g_body_number; i++) {

		for (c = 0; c < 3; c++) {
			g_body_array[i].pos[c] = new_pos[i].val[c];
			g_body_array[i].vel[c] = new_vel[i].val[c];
		}

	}



	*actual_timestep = timestep;
	return(0);
}


	/********************************************************************
	 * PROCEDURE: print_positions
	 *
	 * DESCRIPTION: Print out information on each body,
	 *              one body per line.  Each line has format
	 *
	 *                  time  index mass  px py pz  vx vy vz
	 *
	 * RETURNS:
	 *              0       if all goes well
	 *              1       if an error occurs
	 */

static int
print_positions
	(
	FILE *outfile_fp,             /* I: print to this open file */
	double time_sec               /* I: the current time */
	)
{
	int i;

	if (outfile_fp == NULL) {
		fprintf(stderr, "print_positions: given NULL pointer?! \n");
		return(1);
	}

	for (i = 0; i < g_body_number; i++) {
		fprintf(outfile_fp, "%15.8e %5d %12.5e  %12.5e %12.5e %12.5e  %12.5e %12.5e %12.5e \n",
				time_sec,
				g_body_array[i].index,
				g_body_array[i].mass,
				g_body_array[i].pos[0],
				g_body_array[i].pos[1],
				g_body_array[i].pos[2],
				g_body_array[i].vel[0],
				g_body_array[i].vel[1],
				g_body_array[i].vel[2]);
	}

	return(0);
}


	/********************************************************************
	 * PROCEDURE: recenter_velocities
	 *
	 * DESCRIPTION: Compute the center-of-mass velocity
	 *              of the entire system.  Then subtract
	 *              that value from all bodies, so that
	 *              the center of mass of the system should
	 *              remain motionless as time passes.
	 *
	 * RETURNS:
	 *              0       if all goes well
	 *              1       if an error occurs
	 */

static int
recenter_velocities
	(
	void
	)
{
	int i, c;
	double total_mass;
	VECTOR momentum;
	VECTOR com_velocity;

	total_mass = 0.0;
	for (c = 0; c < 3; c++) {
		momentum.val[c] = 0.0;
	}

	for (i = 0; i < g_body_number; i++) {
		total_mass += g_body_array[i].mass;
		for (c = 0; c < 3; c++) {
			momentum.val[c] += g_body_array[i].mass*g_body_array[i].vel[c];
		}
	}

	/* if the total mass is zero, just return now */
	if (total_mass == 0.0) {
		return(0);
	}

	/* compute the velocity of the center of mass */
	for (c = 0; c < 3; c++) {
		com_velocity.val[c] = momentum.val[c] / total_mass;
	}

	/* now subtract this velocity from each body */
	for (i = 0; i < g_body_number; i++) {
		for (c = 0; c < 3; c++) {
			g_body_array[i].vel[c] -= com_velocity.val[c];
		}
	}

	return(0);
}


	/********************************************************************
	 * PROCEDURE: calc_total_angmom
	 *
	 * DESCRIPTION: Compute the total angular momentum
	 *              of all the bodies, as a three-D VECTOR,
	 *              computed around the center of mass of the system.
	 *              The result will have units of kg*m/s.
	 *
	 *              Place the result into the final arg.
	 *
	 * RETURNS:
	 *              0          if all goes well
	 *              1          if an error occurs
	 */

static int
calc_total_angmom
	(
	VECTOR *result        /* O: the total angular momentum of system */
	)
{
	int i, c;
	double total_mass;
	VECTOR com, com_vel;
	VECTOR body_pos, body_vel;
	VECTOR body_ang_mom, total_ang_mom;

	total_mass = 0.0;
	for (c = 0; c < 3; c++) {
		com.val[c] = 0.0;
		com_vel.val[c] = 0.0;
		total_ang_mom.val[c] = 0.0;
	}

	/*
	 * step 1: compute the center of mass of the system
	 *         and the velocity of the center of mass
	 */
	for (i = 0; i < g_body_number; i++) {
		total_mass += g_body_array[i].mass;
		for (c = 0; c < 3; c++) {
			com.val[c] += g_body_array[i].mass*g_body_array[i].pos[c];
			com_vel.val[c] += g_body_array[i].mass*g_body_array[i].vel[c];
		}
	}
	for (c = 0; c < 3; c++) {
		com.val[c] /= total_mass;
		com_vel.val[c] /= total_mass;
	}


	/* step 2: compute the angular momentum of each body around this point */
	for (i = 0; i < g_body_number; i++) {

		/* compute the position VECTOR from center of mass to body i */
		body_pos.val[0] = g_body_array[i].pos[0] - com.val[0];
		body_pos.val[1] = g_body_array[i].pos[1] - com.val[1];
		body_pos.val[2] = g_body_array[i].pos[2] - com.val[2];

		/*
		 * compute the velocity difference between center of mass
		 * and body i, and weight it by the mass of the body
		 */
		body_vel.val[0] = g_body_array[i].mass *
						(g_body_array[i].vel[0] - com_vel.val[0]);
		body_vel.val[1] = g_body_array[i].mass *
						(g_body_array[i].vel[1] - com_vel.val[1]);
		body_vel.val[2] = g_body_array[i].mass *
						(g_body_array[i].vel[2] - com_vel.val[2]);

		/* now calculate the angular momentum of body i around COM */
		VECTOR_CROSS(body_pos, body_vel, body_ang_mom);

		/* and add this body's angular momentum to the total */
		for (c = 0; c < 3; c++) {
			total_ang_mom.val[c] += body_ang_mom.val[c];
		}

	}

	/* all done.  Copy the total angular momentum to output arg */
	for (c = 0; c < 3; c++) {
		result->val[c] = total_ang_mom.val[c];
	}

	return(0);
}



	/********************************************************************
	 * PROCEDURE: tech_euler_2
	 *
	 * DESCRIPTION: Advance particles one timestep
	 *              using Euler's method to second order,
	 *              which is also known as Heun's method.
	 *
	 *              The basic idea:
	 *
	 *                a) compute current accel
	 *                b) use current accel to predict poor future vel
	 *                c) use poor future vel to predict poor future pos
	 *                d) use poor future pos to predict poor future accel
	 *                e) calc average of current and future accel
	 *                f) use average accel to compute better future vel
	 *                g) calc average of current and future vel
	 *                h) use average vel to compute better future pos
	 *
	 *              Thi method does not (yet) attempt to change
	 *              the timestep.
	 *
	 * RETURNS:
	 *              0       if all goes well
	 *              1       if an error occurs
	 */

static int
tech_euler_2
	(
	double suggested_timestep,    /* I: suggested timestep (sec) */
	double *actual_timestep       /* O: actual timestep used (sec) */
	)
{
	int i, j, c;
	double timestep;
	VECTOR poor_new_pos[MAX_NBODY];
	VECTOR better_new_pos[MAX_NBODY];
	VECTOR poor_new_vel[MAX_NBODY];
	VECTOR better_new_vel[MAX_NBODY];
	VECTOR average_vel[MAX_NBODY];
	VECTOR cur_force[MAX_NBODY][MAX_NBODY];
	VECTOR cur_tot_force[MAX_NBODY];
	VECTOR poor_new_force[MAX_NBODY];
	VECTOR average_force[MAX_NBODY];

	timestep = suggested_timestep;

	/*
	 * compute the forces between all objects
	 *   using the current positions
	 */
	for (i = 0; i < g_body_number; i++) {
		for (j = i; j < g_body_number; j++) {
			VECTOR dist, dist_frac;
			double dist_tot;
			double dist_squared;
			double top, force_mag;

			if (i == j) {
				cur_force[i][j].val[0] = 0.0;
				cur_force[i][j].val[1] = 0.0;
				cur_force[i][j].val[2] = 0.0;
				continue;
			}

			SET_POS_DIFF(dist, i, j);
			dist_squared = VECTOR_MAG_SQUARED(dist);
			dist_tot = sqrt(dist_squared);
			if (dist_tot <= 0.0) {
				fprintf(stderr,
						" tech_euler_2: distance of %le -- quitting \n", dist_tot);
				return(1);
			}

			/* this is where we could add a softening parameter */
			top = G * g_body_array[i].mass * g_body_array[j].mass;
			force_mag = top / dist_squared;

			for (c = 0; c < 3; c++) {
				dist_frac.val[c] = dist.val[c] / dist_tot;
				cur_force[i][j].val[c] = 0.0 - (force_mag * dist_frac.val[c]);
				cur_force[j][i].val[c] = 0.0 - cur_force[i][j].val[c];
			}
		}
	}


	/*
	 * first, we go through several steps to predict poor future
	 *    position for all object
	 */
	for (i = 0; i < g_body_number; i++) {

		/*
		 * we will compute a new position and velocity
		 *   for each object in turn, placing the new
		 *   values into the 'new_pos[]' and 'new_vel[]'
		 *   arrays.  Afterwards,
		 *   we'll copy the new values back into the
		 *   original g_body[] array.
		 */


		/* use the current accel to predict poor future vel */
		/*
		 * compute the total force on this body
		 *    in each direction
		 */
		for (c = 0; c < 3; c++) {
			cur_tot_force[i].val[c] = 0.0;
		}
		for (j = 0; j < g_body_number; j++) {
			for (c = 0; c < 3; c++) {
				cur_tot_force[i].val[c] += cur_force[i][j].val[c];
			}
		}
		for (c = 0; c < 3; c++) {
			poor_new_vel[i].val[c] = g_body_array[i].vel[c]
						+ (cur_tot_force[i].val[c]/g_body_array[i].mass)*timestep;
		}

		/* use poor future vel to predict poor future pos */
		for (c = 0; c < 3; c++) {
			poor_new_pos[i].val[c] = g_body_array[i].pos[c]
							+ poor_new_vel[i].val[c]*timestep;
		}
	}

	/*
	 * Now we use poor future pos of all objects to predict
	 *    poor future accel for all objects
	 */
	for (i = 0; i < g_body_number; i++) {
		for (c = 0; c < 3; c++) {
			poor_new_force[i].val[c] = 0.0;
		}
		for (j = 0; j < g_body_number; j++) {
			VECTOR dist, dist_frac;
			double dist_tot;
			double dist_squared;
			double top, force_mag;

			if (i == j) {
				continue;
			}

			for (c = 0; c < 3; c++) {
				dist.val[c] = poor_new_pos[i].val[c] - poor_new_pos[j].val[c];
			}
			dist_squared = VECTOR_MAG_SQUARED(dist);
			dist_tot = sqrt(dist_squared);
			if (dist_tot <= 0.0) {
				fprintf(stderr,
						" tech_euler_2: distance of %le -- quitting \n", dist_tot);
				return(1);
			}

			/* this is where we could add a softening parameter */
			top = G * g_body_array[i].mass * g_body_array[j].mass;
			force_mag = top / dist_squared;

			for (c = 0; c < 3; c++) {
				dist_frac.val[c] = dist.val[c] / dist_tot;
				poor_new_force[i].val[c] += 0.0 - (force_mag * dist_frac.val[c]);
			}
		}
	}

	/*
	 * and now we can walk through the list of all bodies,
	 *    one at a time, and calculate improved future quantities
	 */
	for (i = 0; i < g_body_number; i++) {

		/* calc average of current and future accel */
		for (c = 0; c < 3; c++) {
			average_force[i].val[c] =
					 0.5 * (cur_tot_force[i].val[c] + poor_new_force[i].val[c]);
		}

		/* use average accel to compute better future vel */
		for (c = 0; c < 3; c++) {
			better_new_vel[i].val[c] = g_body_array[i].vel[c] +
						(average_force[i].val[c]/g_body_array[i].mass)*timestep;
		}

		/* calc average of current and future vel */
		for (c = 0; c < 3; c++) {
			average_vel[i].val[c] =
						0.5 * (g_body_array[i].vel[c] + better_new_vel[i].val[c]);
		}

		/* use average vel to compute better future pos */
		for (c = 0; c < 3; c++) {
			better_new_pos[i].val[c] = g_body_array[i].pos[c] +
						average_vel[i].val[c]*timestep;
		}

	}


	/*
	 * having computed the new positions and velocities,
	 *    we now copy those values into the g_body_array[]
	 */
	for (i = 0; i < g_body_number; i++) {

		for (c = 0; c < 3; c++) {
			g_body_array[i].pos[c] = better_new_pos[i].val[c];
			g_body_array[i].vel[c] = better_new_vel[i].val[c];
		}

	}



	*actual_timestep = timestep;
	return(0);
}


	/********************************************************************
	 * PROCEDURE: tech_euler_1a
	 *
	 * DESCRIPTION: Advance particles one timestep
	 *              using Euler's method to first order.
	 *              However, we use the velocity for time N+1
	 *              to advance the position for time N;
	 *              in other words, use the new velocity
	 *              to advance the old position.
	 *
	 *              This method does not (yet) attempt to change
	 *              the timestep.
	 *
	 * RETURNS:
	 *              0       if all goes well
	 *              1       if an error occurs
	 */

static int
tech_euler_1a
	(
	double suggested_timestep,    /* I: suggested timestep (sec) */
	double *actual_timestep       /* O: actual timestep used (sec) */
	)
{
	int i, j, c;
	double timestep;
	VECTOR new_pos[MAX_NBODY];
	VECTOR new_vel[MAX_NBODY];
	VECTOR cur_force[MAX_NBODY][MAX_NBODY];

	timestep = suggested_timestep;

	/*
	 * compute the forces between all objects
	 *   using the current positions
	 */
	for (i = 0; i < g_body_number; i++) {
		for (j = i; j < g_body_number; j++) {
			VECTOR dist, dist_frac;
			double dist_tot;
			double dist_squared;
			double top, force_mag;

			if (i == j) {
				cur_force[i][j].val[0] = 0.0;
				cur_force[i][j].val[1] = 0.0;
				cur_force[i][j].val[2] = 0.0;
				continue;
			}

			SET_POS_DIFF(dist, i, j);
			dist_squared = VECTOR_MAG_SQUARED(dist);
			dist_tot = sqrt(dist_squared);
			if (dist_tot <= 0.0) {
				fprintf(stderr,
						" tech_euler_1: distance of %le -- quitting \n", dist_tot);
				return(1);
			}

			/* this is where we could add a softening parameter */
			top = G * g_body_array[i].mass * g_body_array[j].mass;
			force_mag = top / dist_squared;

			for (c = 0; c < 3; c++) {
				dist_frac.val[c] = dist.val[c] / dist_tot;
				cur_force[i][j].val[c] = 0.0 - (force_mag * dist_frac.val[c]);
				cur_force[j][i].val[c] = 0.0 - cur_force[i][j].val[c];
			}
		}
	}


	/* do the work here */
	for (i = 0; i < g_body_number; i++) {
		VECTOR tot_force;
		VECTOR tot_accel;

		/*
		 * we will compute a new position and velocity
		 *   for each object in turn, placing the new
		 *   values into the 'new_pos[]' and 'new_vel[]'
		 *   arrays.  Afterwards,
		 *   we'll copy the new values back into the
		 *   original g_body[] array.
		 */


		/* use the current forces to compute new velocities */
		for (c = 0; c < 3; c++) {

			/*
			 * compute the total force on this body
			 *    in each direction
			 */
			tot_force.val[c] = 0.0;
			for (j = 0; j < g_body_number; j++) {
				tot_force.val[c] += cur_force[i][j].val[c];
			}

			/* from total force, compute total acceleration */
			tot_accel.val[c] = tot_force.val[c] / g_body_array[i].mass;

			/* use the acceleration to compute new velocity */
			new_vel[i].val[c] = g_body_array[i].vel[c]
							+ tot_accel.val[c]*timestep;
		}

		/* use the NEW velocity to compute new positions */
		for (c = 0; c < 3; c++) {
			new_pos[i].val[c] = g_body_array[i].pos[c]
							+ new_vel[i].val[c]*timestep;
		}


	}


	/*
	 * having computed the new positions and velocities,
	 *    we now copy those values into the g_body_array[]
	 */
	for (i = 0; i < g_body_number; i++) {

		for (c = 0; c < 3; c++) {
			g_body_array[i].pos[c] = new_pos[i].val[c];
			g_body_array[i].vel[c] = new_vel[i].val[c];
		}

	}



	*actual_timestep = timestep;
	return(0);
}



	/********************************************************************
	 * PROCEDURE: compute_body_forces
	 *
	 * DESCRIPTION: Given an array of BODYs, compute the gravitational
	 *              forces between each pair.  We place the forces
	 *              into the array of VECTORs which is given
	 *              as the final argument.
	 *
	 * RETURNS:
	 *              0       if all goes well
	 *              1       if an error occurs
	 */

static int
calc_grav_forces
	(
	int nbody,                    /* I: number of bodies */
	BODY *body_array,             /* I: array with pos and mass of all bodies */
	VECTOR *forces                /* O: we put forces on each body here */
	)
{
	int i, j, c;
	VECTOR cur_force[MAX_NBODY][MAX_NBODY];

	/*
	 * compute the forces between all objects
	 *   using the current positions
	 */
	for (i = 0; i < nbody; i++) {
		for (j = i; j < nbody; j++) {
			VECTOR dist, dist_frac;
			double dist_tot;
			double dist_squared;
			double top, force_mag;

			if (i == j) {
				cur_force[i][j].val[0] = 0.0;
				cur_force[i][j].val[1] = 0.0;
				cur_force[i][j].val[2] = 0.0;
				continue;
			}

			dist.val[0] = body_array[i].pos[0] - body_array[j].pos[0];
			dist.val[1] = body_array[i].pos[1] - body_array[j].pos[1];
			dist.val[2] = body_array[i].pos[2] - body_array[j].pos[2];
			dist_squared = VECTOR_MAG_SQUARED(dist);
			dist_tot = sqrt(dist_squared);
			if (dist_tot <= 0.0) {
				fprintf(stderr,
						" calc_grav_forces: distance of %le -- quitting \n", dist_tot);
				return(1);
			}

			/* this is where we could add a softening parameter */
			top = G * body_array[i].mass * body_array[j].mass;
			force_mag = top / dist_squared;

			for (c = 0; c < 3; c++) {
				dist_frac.val[c] = dist.val[c] / dist_tot;
				cur_force[i][j].val[c] = 0.0 - (force_mag * dist_frac.val[c]);
				cur_force[j][i].val[c] = 0.0 - cur_force[i][j].val[c];
			}
		}
	}

	/*
	 * compute the total force on each body
	 *    in each direction
	 */
	for (i = 0; i < nbody; i++) {
		for (c = 0; c < 3; c++) {
			forces[i].val[c] = 0.0;
			for (j = 0; j < g_body_number; j++) {
				forces[i].val[c] += cur_force[i][j].val[c];
			}
		}
	}

	return(0);
}



	/********************************************************************
	 * PROCEDURE: tech_rk4
	 *
	 * DESCRIPTION: Advance particles one timestep
	 *              using Runge-Kutta fourth-order method.
	 *
	 *              This method does not (yet) attempt to change
	 *              the timestep.
	 *
	 * RETURNS:
	 *              0       if all goes well
	 *              1       if an error occurs
	 */

static int
tech_rk4
	(
	double suggested_timestep,    /* I: suggested timestep (sec) */
	double *actual_timestep       /* O: actual timestep used (sec) */
	)
{
	int i, c;
	double timestep;
	double half_step, sixth_step;
	VECTOR v1[MAX_NBODY];
	VECTOR v2[MAX_NBODY];
	VECTOR v3[MAX_NBODY];
	VECTOR v4[MAX_NBODY];
	VECTOR a1[MAX_NBODY];
	VECTOR a2[MAX_NBODY];
	VECTOR a3[MAX_NBODY];
	VECTOR a4[MAX_NBODY];
	BODY step2_array[MAX_NBODY];
	BODY step3_array[MAX_NBODY];
	BODY step4_array[MAX_NBODY];

	timestep = suggested_timestep;

	/*
	 * step 1: use current positions to compute accelerations.
	 *         Call these "a1".
	 */
	copy_bodies(g_body_number, g_body_array, step2_array);
	if (calc_grav_forces(g_body_number, step2_array, a1) != 0) {
		fprintf(stderr, "tech_rk4: calc_grav_forces fails in step a1\n");
		return(1);
	}
	/* convert those forces to accelerations */
	for (i = 0; i < g_body_number; i++) {
		for (c = 0; c < 3; c++) {
			a1[i].val[c] /= step2_array[i].mass;
		}
	}

	/*
	 * step 2: copy the current velocities into an array for later use.
	 *         Call these "v1".
	 */
	for (i = 0; i < g_body_number; i++) {
		for (c = 0; c < 3; c++) {
			v1[i].val[c] = g_body_array[i].vel[c];
		}
	}

	/*
	 * step 3: use the "v1" velocities to predict positions of objects
	 *         one-half a step into the future.  We'll place the
	 *         positions into the "step2_array[]."
	 */
	half_step = timestep*0.5;
	for (i = 0; i < g_body_number; i++) {
		for (c = 0; c < 3; c++) {
			step2_array[i].pos[c] += half_step*v1[i].val[c];
		}
	}

	/*
	 * step 4: use the "a1" accelerations to predict future velocities
	 *         at half a step into the future.  Call these "v2".
	 */
	for (i = 0; i < g_body_number; i++) {
		for (c = 0; c < 3; c++) {
			v2[i].val[c] = v1[i].val[c] + half_step*a1[i].val[c];
		}
	}

	/*
	 * step 5: use the "step2_array" positions (at half a step in future)
	 *         to compute forces on the bodies in the future.
	 *         Convert these forces to accelerations and place into
	 *         the "a2" array.
	 */
	if (calc_grav_forces(g_body_number, step2_array, a2) != 0) {
		fprintf(stderr, "tech_rk4: calc_grav_forces fails in step a2\n");
		return(1);
	}
	/* convert those forces to accelerations */
	for (i = 0; i < g_body_number; i++) {
		for (c = 0; c < 3; c++) {
			a2[i].val[c] /= step2_array[i].mass;
		}
	}

	/*
	 * step 6: use the "a2" accelerations to predict future velocities
	 *         at half a step into the future.  Call these "v3".
	 */
	for (i = 0; i < g_body_number; i++) {
		for (c = 0; c < 3; c++) {
			v3[i].val[c] = v1[i].val[c] + half_step*a2[i].val[c];
		}
	}

	/*
	 * step 6: use the "v2" velocities to predict positions of objects
	 *         one-half a step into the future.  We'll place the
	 *         positions into the "step3_array[]."
	 */
	copy_bodies(g_body_number, g_body_array, step3_array);
	for (i = 0; i < g_body_number; i++) {
		for (c = 0; c < 3; c++) {
			step3_array[i].pos[c] += half_step*v2[i].val[c];
		}
	}

	/*
	 * step 7: use the "step3_array" positions (at half a step in future)
	 *         to compute forces on the bodies in the future.
	 *         Convert these forces to accelerations and place into
	 *         the "a3" array.
	 */
	if (calc_grav_forces(g_body_number, step3_array, a3) != 0) {
		fprintf(stderr, "tech_rk4: calc_grav_forces fails in step a3\n");
		return(1);
	}
	/* convert those forces to accelerations */
	for (i = 0; i < g_body_number; i++) {
		for (c = 0; c < 3; c++) {
			a3[i].val[c] /= step3_array[i].mass;
		}
	}

	/*
	 * step 8: use the "a3" accelerations to predict future velocities
	 *         at one FULL step into the future.  Call these "v4".
	 */
	for (i = 0; i < g_body_number; i++) {
		for (c = 0; c < 3; c++) {
			v4[i].val[c] = v1[i].val[c] + timestep*a3[i].val[c];
		}
	}

	/*
	 * step 9: use the "v3" velocities to predict positions of objects
	 *         one FULL step into the future.  We'll place the
	 *         positions into the "step4_array[]."
	 */
	copy_bodies(g_body_number, g_body_array, step4_array);
	for (i = 0; i < g_body_number; i++) {
		for (c = 0; c < 3; c++) {
			step4_array[i].pos[c] += timestep*v3[i].val[c];
		}
	}

	/*
	 * step 10: use the "step4_array" positions (at one FULL step in future)
	 *         to compute forces on the bodies in the future.
	 *         Convert these forces to accelerations and place into
	 *         the "a4" array.
	 */
	if (calc_grav_forces(g_body_number, step4_array, a4) != 0) {
		fprintf(stderr, "tech_rk4: calc_grav_forces fails in step a4\n");
		return(1);
	}
	/* convert those forces to accelerations */
	for (i = 0; i < g_body_number; i++) {
		for (c = 0; c < 3; c++) {
			a4[i].val[c] /= step4_array[i].mass;
		}
	}


	/*
	 * Okay, at this point, we should have 4 velocities, and 4 accelerations
	 *   for each object.
	 *
	 *         v1     current velocity
	 *         v2     velocity predicted half a step in future
	 *         v3     improved velocity predicted half a step in future
	 *         v4     velocity predicted one FULL step in future
	 *
	 *         a1     current acceleration
	 *         a2     acceleration predicted half a step in future
	 *         a3     improved acceleration predicted half a step in future
	 *         a4     acceleration predicted one FULL step in future
	 *
	 * We can now combine these 4 measurements to produce one good
	 *   value of velocity (or acceleration) one full step into the future.
	 */
	sixth_step = timestep / 6.0;
	for (i = 0; i < g_body_number; i++) {
		for (c = 0; c < 3; c++) {
			double delta_vel, delta_pos;

			delta_vel = sixth_step *
					(a1[i].val[c] + 2*a2[i].val[c] + 2*a3[i].val[c] + a4[i].val[c]);
			g_body_array[i].vel[c] += delta_vel;

			delta_pos = sixth_step *
					(v1[i].val[c] + 2*v2[i].val[c] + 2*v3[i].val[c] + v4[i].val[c]);
			g_body_array[i].pos[c] += delta_pos;

		}
	}

	*actual_timestep = timestep;

	return(0);
}




	/********************************************************************
	 * PROCEDURE: copy_bodies
	 *
	 * DESCRIPTION: Copy an array of BODY structures into a second array.
	 *              We can then modify the copies and leave the originals
	 *              untouched.
	 *
	 * RETURNS:
	 *              nothing
	 */

static void
copy_bodies
	(
	int num_bodies,               /* I: number of bodies to be copied */
	BODY *from_array,             /* I: copy from this array ... */
	BODY *to_array                /* O: ... into this array  */
	)
{
	int i, c;
	BODY *from_body, *to_body;

	for (i = 0; i < num_bodies; i++) {
		from_body = &(from_array[i]);
		to_body = &(to_array[i]);

		to_body->index = from_body->index;
		strcpy(from_body->name, to_body->name);
		to_body->mass = from_body->mass;
		for (c = 0; c < 3; c++) {
			to_body->pos[c] = from_body->pos[c];
			to_body->vel[c] = from_body->vel[c];
		}
	}

}
