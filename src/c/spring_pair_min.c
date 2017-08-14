#include <stdio.h>
#include <errno.h>
#include <math.h>

#include "mol2/benergy.h"
#include "mol2/gbsa.h"
#include "mol2/icharmm.h"
#include "mol2/minimize.h"
#include "mol2/nbenergy.h"
#include "mol2/pdb.h"

#define __TOL__ 5E-4

struct energy_prm {

	struct mol_atom_group *ag;
	struct agsetup  *ag_setup;
	struct acesetup *ace_setup;

    struct springset* sprst;
};

struct spring
{
        struct mol_atom_group *ag;  /**< affected atomgroup */
        int    naspr;      /**< number of affected atoms */
        int   *laspr;      /**< pair of atoms  */
        double fkspr;      /**< force constant */
        double lnspr;      /**< spring length  */
};

struct springset
{
        int nsprings;            /**< number of springs */
        struct spring *springs;  /**< array of springs */
};

void* mymalloc (size_t size)
{
	void* v = (void*) malloc (size);
	if (v == NULL)
	{
		exit (EXIT_FAILURE);
	}
	return v;
}

void read_fix(char *ffile, int *nfix, int **fix);
 
void read_crd(char *crdfile, int *ncrd, double **crd);

void read_springset(struct mol_atom_group* ag, char *sfile, struct springset** sprst);

void free_springset(struct springset *sprst);

void springeng(struct springset *sprst, double* een);

static lbfgsfloatval_t energy_func(
	void* restrict prm,
	const double* restrict array,
	double* restrict gradient,
	const int array_size,
	const lbfgsfloatval_t step);



int main(int argc, char** argv) 
{
	char* pdb   = argv[1];
	char* psf   = argv[2];
	char* prm   = argv[3];
	char* rtf   = argv[4];
	char* ffix  = argv[5];
	char* sfile = argv[6];
	char* out   = argv[7];

	struct mol_atom_group *ag = mol_read_pdb(pdb);

	mol_atom_group_read_geometry(ag, psf, prm, rtf);
	ag->gradients = calloc(ag->natoms, sizeof(struct mol_vector3));
	mol_fixed_init(ag);	

	int  nfix = 0;
	int* fix;
	read_fix(ffix, &nfix, &fix); 
	mol_fixed_update(ag, nfix, fix);
	//mol_fixed_update(ag, 0, NULL);

	struct mol_vector3 *orig_coords = calloc(ag->natoms, sizeof(struct mol_vector3));
	memcpy(orig_coords, ag->coords, ag->natoms*sizeof(struct mol_vector3));
	
	struct agsetup ags;
	init_nblst(ag, &ags);
	update_nblst(ag, &ags);
	
	struct acesetup ace_setup;
	ace_setup.efac = 0.5;
	ace_ini(ag, &ace_setup);
	ace_fixedupdate(ag, &ags, &ace_setup);
	ace_updatenblst(&ags, &ace_setup);
	
	struct springset *sprst;
    read_springset(ag, sfile, &sprst);

	struct energy_prm engpar;
	engpar.ag = ag;
	engpar.ag_setup  = &ags;
	engpar.ace_setup = &ace_setup;
    engpar.sprst = sprst;

	mol_minimize_ag(MOL_LBFGS, 1000, 1E-3, ag, (void *)(&engpar), energy_func);

	mol_write_pdb(out, ag);
	
	return 0;
}


static lbfgsfloatval_t energy_func(
	void* restrict prm,
	const double* restrict array,
	double* restrict gradient,
	const int array_size,
	const lbfgsfloatval_t step)
{
	lbfgsfloatval_t energy = 0.0;
	struct energy_prm* energy_prm = (struct energy_prm*) prm;
    
	if (array != NULL) {
		//ck_assert((size_t) array_size == energy_prm->ag->active_atoms->size * 3);
		mol_atom_group_set_actives(energy_prm->ag, array);
	}
	bool updated = check_clusterupdate(energy_prm->ag, energy_prm->ag_setup);
	if (updated) {
		ace_updatenblst(energy_prm->ag_setup, energy_prm->ace_setup);
	}

	//reset energy
	zero_grads(energy_prm->ag);

	//energy calculations
	aceeng(energy_prm->ag, &energy, energy_prm->ace_setup, energy_prm->ag_setup);
	vdweng(energy_prm->ag, &energy, energy_prm->ag_setup->nblst);
	vdwengs03(1.0, energy_prm->ag_setup->nblst->nbcof, energy_prm->ag, &energy,
		  energy_prm->ag_setup->nf03, energy_prm->ag_setup->listf03);
	beng(energy_prm->ag, &energy);
	aeng(energy_prm->ag, &energy);
	teng(energy_prm->ag, &energy);
	ieng(energy_prm->ag, &energy);

	if(energy_prm->sprst->nsprings <= 0)
    {
       printf("my_en_grad WARNING: no springs\n");
    }

	springeng(energy_prm->sprst, &energy);

	if (gradient != NULL) {
		for (int i = 0; i < array_size / 3; i++ ) {
			int atom_i = energy_prm->ag->active_atoms->members[i];
			gradient[3*i  ] = -energy_prm->ag->gradients[atom_i].X;
			gradient[3*i+1] = -energy_prm->ag->gradients[atom_i].Y;
			gradient[3*i+2] = -energy_prm->ag->gradients[atom_i].Z;
		}
	}

	printf("%f.4\n", energy);
	return energy;
}


void read_fix(char *ffile, int *nfix, int **fix) 
{ 
   int linesz=91; 
   char *buffer=mymalloc(sizeof(char)*linesz); 
   *nfix=0; 
   FILE* fp = fopen (ffile, "r"); 
   while (fgets(buffer, linesz-1, fp)!=NULL) 
   { 
       if(!strncmp(buffer,"ATOM",4))(*nfix)++; 
   } 
   fclose(fp); 
   *fix=mymalloc(*nfix*sizeof(int)); 
   fp = fopen (ffile, "r"); 
   int na = 0; 
   while(fgets(buffer, linesz-1, fp)!=NULL) 
   { 
       if(!strncmp(buffer, "ATOM", 4)) 
       { 
         (*fix)[na]=atoi(buffer+4)-1; 
         na++; 
       } 
   } 
   free(buffer); 
   fclose(fp); 
}

void read_springset(struct mol_atom_group* ag, char *sfile, struct springset** sprst)
{
	int linesz   = 91;
	int nsprings = 0;
	int i, j, nfix;
	int    *fix;
	double f, ln;
	char   *buffer=mymalloc(sizeof(char)*linesz);

	FILE* fp = fopen (sfile, "r");
	if(fgets(buffer, linesz-1, fp) == NULL)
	{
		printf("Error reading springs");
		exit(0);
	}

	nsprings = atoi(buffer);
	if (nsprings <= 0)
	{
		printf("Error reading springs");
		exit(0);
	}
	
	// read name of pdb with springs
	if(fgets(buffer, linesz-1, fp) == NULL)
	{
		printf("Error reading springs");
		exit(0);
	}  

	for (j = linesz-1; j >= 0; j--)
	{
		if (buffer[j] == 'b' || buffer[j] == 'B') 
			break;
			
		buffer[j]='\0';
	}

	// read spring atoms
	read_fix(buffer, &nfix, &fix);
	if (nfix % 2 != 0)
	{
		printf("Spring %i: expecting even number of atoms (%i atoms declared)\n", i, nfix);
		exit(1);
	}


	(*sprst) = mymalloc(sizeof(struct springset));
	(*sprst)->nsprings = nsprings;
	(*sprst)->springs = mymalloc(nsprings*sizeof(struct spring));
	
	for (i = 0; i < nsprings; i++)
	{
		(*sprst)->springs[i].ag = ag;
		(*sprst)->springs[i].naspr = 2;
		(*sprst)->springs[i].laspr = mymalloc(2*sizeof(int));
		
		(*sprst)->springs[i].laspr[0] = fix[i*2  ];
		(*sprst)->springs[i].laspr[1] = fix[i*2+1];

		// read spring length
		if(fgets(buffer, linesz-1, fp) == NULL)
		{
			printf("Error reading springs");
			exit(0);
		}

		ln = atof(buffer);
		(*sprst)->springs[i].lnspr = ln;

		// read spring fk
		if(fgets(buffer, linesz-1, fp) == NULL)
		{
			printf("Error reading springs");
			exit(0);
		}

		f = atof(buffer);
		(*sprst)->springs[i].fkspr = f;
	}    
	    
	free(buffer);
	free(fix);
}

void free_springset(struct springset *sprst)
{
        int i;
        for(i=0; i<sprst->nsprings; i++)
              free(sprst->springs[i].laspr);
        if(sprst->springs != NULL)
               free(sprst->springs);
        if(sprst!= NULL)
               free(sprst);
} 

void springeng(struct springset *sprst, double* een)
{
	int    i, i1, i2, nat;
	double xtot, ytot, ztot, fk, d, d2, ln, coef;
	struct mol_vector3 g;
	struct mol_atom_group *ag;

	for (i = 0; i < sprst->nsprings; i++)
	{
		ag = sprst->springs[i].ag;
		
		nat = sprst->springs[i].naspr;
		ln  = sprst->springs[i].lnspr;
		fk  = sprst->springs[i].fkspr;
		
		i1 = sprst->springs[i].laspr[0];
		i2 = sprst->springs[i].laspr[1];
		
		xtot = ag->coords[i2].X - ag->coords[i1].X;
		ytot = ag->coords[i2].Y - ag->coords[i1].Y;
		ztot = ag->coords[i2].Z - ag->coords[i1].Z;
		
		d2 = xtot*xtot + ytot*ytot + ztot*ztot;
		d  = sqrt(d2);
		(*een) += fk * (d2 - 2 * d * ln + ln * ln);
		
		coef = fk * 2 * (1.0 - ln / d);
		
		g.X = - coef * xtot;
		g.Y = - coef * ytot;
		g.Z = - coef * ztot;
		
		MOL_VEC_SUB(ag->gradients[i1], ag->gradients[i1], g);
		MOL_VEC_ADD(ag->gradients[i2], ag->gradients[i2], g);
	}
}
