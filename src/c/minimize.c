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
	struct mol_atom_group  *ag;
	struct agsetup   *ag_setup;
	struct acesetup *ace_setup;
};

void* mymalloc (size_t size)
{
//	_PRINT_DEPRECATED_
//	fprintf (stderr, "\t(please write your own malloc wrapper)\n");

	void* v = (void*) malloc (size);
	if (v == NULL)
	{
		exit (EXIT_FAILURE);
	}
	return v;
}

/*void read_fix(char *ffile, int *nfix, int **fix) 
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
   int na=0; 
   while(fgets(buffer, linesz-1, fp)!=NULL) 
   { 
       if(!strncmp(buffer,"ATOM",4)) 
       { 
         (*fix)[na]=atoi(buffer+4)-1; 
         na++; 
       } 
   } 
   free(buffer); 
 //  free(fix); 
   fclose(fp); 
}*/

void read_fix(char *ffile, int *nfix, int **fix) 
{ 
   int linesz=10; 
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
   int na=0; 
   while(fgets(buffer, linesz-1, fp)!=NULL) 
   { 
         (*fix)[na]=atoi(buffer)-1; 
         na++;
   } 
   free(buffer); 
   fclose(fp); 
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

	if (gradient != NULL) {
		for (int i = 0; i < array_size / 3; i++ ) {
			int atom_i = energy_prm->ag->active_atoms->members[i];
			gradient[3*i]     = -energy_prm->ag->gradients[atom_i].X;
			gradient[(3*i)+1] = -energy_prm->ag->gradients[atom_i].Y;
			gradient[(3*i)+2] = -energy_prm->ag->gradients[atom_i].Z;
		}
	}

	//printf("%f.4\n", energy);
	return energy;
}

int main(int argc, char** argv) 
{
	char* pdb  = argv[1];
	char* psf  = argv[2];
	char* prm  = argv[3];
	char* rtf  = argv[4];
	char* ffix = argv[5];

	struct mol_atom_group *ag = mol_read_pdb(pdb);
	//struct mol_atom_group *ag = mol_read_pdb("mutated3_nmin.pdb");	
	mol_atom_group_read_geometry(ag, psf, prm, rtf);
	ag->gradients = calloc(ag->natoms, sizeof(struct mol_vector3));
	
	mol_fixed_init(ag);
	
	//int nfix = 0;
	//int* fix;
	//read_fix(ffix, &nfix, &fix); 
	//mol_fixed_update(ag, nfix, fix);
	mol_fixed_update(ag, 0, NULL);
	
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

	struct energy_prm engpar;
	engpar.ag = ag;
	engpar.ag_setup = &ags;
	engpar.ace_setup = &ace_setup;

	mol_minimize_ag(MOL_LBFGS, 1000, 1E-3, ag, (void *)(&engpar), energy_func);
	
	mol_write_pdb("minimized.pdb", ag);
	
	return 0;
}
