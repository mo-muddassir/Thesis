#include "dnc.h"

Node *read_particles(char *fname, int *N, double *time)
{
	Particle *p;
	Node *n;
	FILE *fp;
	int i, id;

	fp = fopen(fname, "r");
	if (fp == NULL) {
		fprintf(stderr, "ERROR: Unable to open input file %s\n", fname);
		exit(2);
	}
	
    fread(N, sizeof(int), 1, fp);
	fread(time, sizeof(double), 1, fp);

	n = (Node *) calloc(*N, sizeof(Node));
	p = (Particle *) calloc(*N, sizeof(Particle));

	for (i = 0; i < *N; i++) {
		/* link the particle to the node */
		n[i].d.p = p + i;
		
		/* read in the details */
		fread(&(MASS(n + i)), sizeof(float), 1, fp);
		fread(POS(n + i), sizeof(float), 3, fp);
		fread(VEL(n + i), sizeof(float), 3, fp);
		fread(&id, sizeof(int), 1, fp);
		
		TYPE(n + i) = PARTICLE;
		ID(n + i) = id;

		/* assume all forces are to be calculated */
		UPDATE_FORCE(n + i);
	}
	
	fclose(fp);
	return n;
}

/* int cmp_id(const void *_n1, const void *_n2)
{
	Node *n1 = (Node *) _n1;
	Node *n2 = (Node *) _n2;
	if (ID(n1) > ID(n2))
		return 1;
	else
		return -1;

} */

void save_particles(Node * p, int N, double time, char *base)
{
    int i;
	char fname[1000];
	FILE *fp;
	
	sprintf(fname, "%s.%09.4f", base, time);
	printf("Saving to %s\n", fname);
	fflush(stdout);
	
	fp = fopen(fname, "w");
	
	/* make sure original particle order is kept 
	qsort(p, N, sizeof(Node), cmp_id); */

	fwrite(&N, sizeof(int), 1, fp);
	fwrite(&time, sizeof(double), 1, fp);
	for (i = 0; i < N; i++) {
		fwrite(&(MASS(p + i)), sizeof(float), 1, fp);
		fwrite(POS(p + i), sizeof(float), 3, fp);
		fwrite(VEL(p + i), sizeof(float), 3, fp);
		fwrite(&ID(p + i), sizeof(int), 1, fp);
	}
	
	fclose(fp);
}
