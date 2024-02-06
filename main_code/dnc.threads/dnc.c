#include "dnc.h"
#include "calc_forces.h"
#include "particle_io.h"
#include "move.h"

#include <time.h>
#include <sys/timeb.h>
#include <signal.h>
#include <getopt.h>
#include <string.h>

/* global variables */
float theta, reciptheta;
int maxdepth, bucketSize;
int nThread;
int dump_and_exit;


void print_help()
{
    	fprintf(stderr, "Usage: dnc [options]\n");
    	fprintf(stderr, "Options:\n");
    	fprintf(stderr, "   -i file     Initial data file (binary) (init.dat)\n");
    	fprintf(stderr, "   -o file     Output file pattern (output)\n");
    	fprintf(stderr, "   -t #        Length of timestep (0.1)\n");
    	fprintf(stderr, "   -s #        Time to stop simulation at (10.0)\n");
    	fprintf(stderr, "   -d #        Output time interval (1.0)\n");
    	fprintf(stderr, "   -e #        Softening length (0.1)\n");
    	fprintf(stderr, "   -a #        Opening angle (0.5)\n");
    	fprintf(stderr, "   -b #        Size of tree bucket (1)\n");
    	fprintf(stderr, "   -n #        Number of threads (4)\n");
    	fprintf(stderr, "   -m #        Max mass of BH / Total System Mass (0.01)\n");
}

void set_dump_and_exit(int sig)
{
	fprintf(stderr, "Finishing current step..... ");
	fprintf(stderr, "(Once more to stop now)\n");
	dump_and_exit = 1;
	signal(sig, SIG_DFL);
	return;
}

/*-----------------------------------------------
 * M A I N : 
 *----------------------------------------------*/
int main(int argc, char *argv[])
{
	int N, i;
	Node *particles;
	double t, dt, dt_out, t_stop, t_out;
	float eps;
	int step;
	char outputfile[200], inputfile[200];
	struct timeb tstart2f, tstop2f;
	float tinterval;	
	int option;
	double Mbh;
	double r;
	double Mbh_max = 0.01;

	/* default input parameters */
	sprintf(inputfile, "init.dat");
	sprintf(outputfile, "output");
	t = 0.0;
	dt = 0.1;
	t_out = 0.0;
	dt_out = 1.0;
	t_stop = 10.0;
	theta = 0.5;
	eps = 0.1;
	bucketSize = 1;
	nThread = 4;
	Mbh = 0.0; 

	printf("Welcome to Dave's N-Body Code\n");
	printf("-----------------------------\n");
	printf("This is a C implementation of Dehnen's (2002) ");
	printf("multipole expansion algorithm.\n\n");

	if (argc < 2) {
		print_help();
		exit(9);
	} else {
		/* process flags */
		while ((option = getopt(argc, argv, "hi:o:t:s:a:e:b:n:d:m:")) != EOF){
			switch (option){
				case 'i':
					strncpy(inputfile, optarg, 200);
					break;				
				case 'o':
					strncpy(outputfile, optarg, 200);
					break;				
				case 't':
					dt = atof(optarg);
					break;					
				case 's':
					t_stop = atof(optarg);
					break;					
				case 'a':
					theta = atof(optarg);
					break;					
				case 'e':
					eps = atof(optarg);
					break;					
				case 'b':
					bucketSize = atoi(optarg);
					break;
				case 'n':
					nThread = atoi(optarg);
					break;					
				case 'd':
					dt_out = atof(optarg);
					break;
				case 'm':
					Mbh_max = atof(optarg);
					break;		
				case 'h':
				case '?':
					print_help();
					exit(9);					
			}
		}
		
		particles = read_particles(inputfile, &N, &t);
		t_out = t + dt_out;
	}

	printf("Using single precision\n");
	printf("N:      %d\n", N);
    printf("t:      %f\n", t);
    printf("dt:     %f\n", dt);
    printf("dt_out: %f\n", dt_out);
	printf("Theta:  %f\n", theta);
	printf("eps:    %f\n", eps);	
	printf("bSize:  %d\n", bucketSize);
	printf("Max BH mass:	%f\n",Mbh_max);
	

	dump_and_exit = 0;
	signal(SIGINT, set_dump_and_exit);
	signal(SIGTERM, set_dump_and_exit);

	/* ----------------------------------------------------------
	   M A I N L O O P S T A R T S H E R E
	   ---------------------------------------------------------- */

	step = 0;
	while (t <= t_stop){

		ftime(&tstart2f);

		if (step > 0)
			AdvancePositions(particles, N, dt);

		for (i = 0; i < N; i++){
			VECTOR_CLEAR(ACC(particles + i));
			r = sqrt(POS(particles+i)[0]*POS(particles+i)[0] + POS(particles+i)[1]*POS(particles+i)[1] +POS(particles+i)[2]*POS(particles+i)[2] + eps*eps);
			
			ACC(particles + i)[0] = -Mbh * POS(particles+i)[0]/r/r/r;
			ACC(particles + i)[1] = -Mbh * POS(particles+i)[1]/r/r/r;
			ACC(particles + i)[2] = -Mbh * POS(particles+i)[2]/r/r/r;
		}
		printf("%f %f\n",ACC(particles)[0], POS(particles)[0]);		
		calc_forces(particles, N, eps);

		if (step > 0)
			AdvanceVelocities(particles, N, 0.5 * dt);

		ftime(&tstop2f);
		tinterval = 1.0 * (tstop2f.time - tstart2f.time) + 
		            0.001 * (tstop2f.millitm - tstart2f.millitm);
		printf("Step %04d    --    t = %09.5f    --    CPU time  %06.3f s\n", step, t, tinterval);

		if (step == 0) {
			float amax, a;

			amax = 0;
			for (i = 0; i < N; i++) {
				VECTOR_MAG2(a, ACC(particles + i));
				if (a > amax)
					amax = a;
			}
			amax = sqrt(amax);
			printf("Maximum acceleration is %e\n", amax);
			printf("Time step should be about %e\n", 0.1 * sqrt(eps / amax));

		}
		
		if (t > t_out){
		    save_particles(particles, N, t, outputfile);
		    t_out += dt_out;
		}
		
		if (dump_and_exit)
		    exit(1);
	    
	    step++;
		t += dt; 
	    	
		AdvanceVelocities(particles, N, 0.5 * dt);
		Mbh += Mbh_max * dt/t_stop;
		printf("Mbh = %f\n",Mbh);		
	}
	
	save_particles(particles, N, t, outputfile);

	return 0;
}

