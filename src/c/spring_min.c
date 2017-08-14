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
        int   *laspr;      /**< list of atoms */
        double fkspr;      /**< force constant */
        double X0, Y0, Z0; /**< anchor point */
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
   fclose(fp); 
}

void read_crd(char *crdfile, int *ncrd, double **crd)
{
   int linesz=91;
   char *buffer=mymalloc(sizeof(char)*linesz);
   *ncrd=0;
   FILE* fp = fopen (crdfile, "r");
   while (fgets(buffer, linesz-1, fp)!=NULL)
   {
       if(!strncmp(buffer,"ATOM",4))(*ncrd)++;
   }
   fclose(fp);
   *crd=mymalloc(*ncrd*3*sizeof(double));
   fp = fopen (crdfile, "r");
   int na=0;
   while(fgets(buffer, linesz-1, fp)!=NULL)
   {
       if(!strncmp(buffer,"ATOM",4))
       {
         (*crd)[na]=atof(buffer+30);
         na++;
         (*crd)[na]=atof(buffer+38);
         na++;
         (*crd)[na]=atof(buffer+46);
         na++;
       }
   }
   free(buffer);
   fclose(fp);
}

void read_springset(struct mol_atom_group* ag, char *sfile, struct springset** sprst)
{
        int linesz=91;
        int nsprings=0;
        int i, j, nfix, ncrd;
        int    *fix;
        double *crd;
        double  X0, Y0, Z0, f;
        char   *buffer=mymalloc(sizeof(char)*linesz);
        
        FILE* fp = fopen (sfile, "r");
        if(fgets(buffer, linesz-1, fp)==NULL)
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
     
        (*sprst) = mymalloc(sizeof(struct springset));
        (*sprst)->nsprings = nsprings;
        (*sprst)->springs = mymalloc(nsprings*sizeof(struct spring));
        
        for (i = 0; i < nsprings; i++)
        {
           if(fgets(buffer, linesz-1, fp)==NULL)
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
 
           read_fix(buffer, &nfix, &fix);
           (*sprst)->springs[i].ag = ag;
           (*sprst)->springs[i].naspr = nfix;
           (*sprst)->springs[i].laspr = mymalloc(nfix*sizeof(int));
         
           for(j = 0; j < nfix; j++)
           {
           		printf("%i\n", fix[j]);
               (*sprst)->springs[i].laspr[j] = fix[j];
           }
               
           free(fix);
                    
           if(fgets(buffer, linesz-1, fp) == NULL)
           {
                 printf("Error reading springs");
                 exit(0);
           }
           
           for(j = linesz-1; j>=0; j--)
           {
              if(buffer[j] == 'b' || buffer[j] == 'B')
              	break;
              	
              buffer[j]='\0';
           }
           
           read_crd(buffer, &ncrd, &crd);
           X0 = 0.0;
           Y0 = 0.0;
           Z0 = 0.0;
           
           for(j=0; j<ncrd; j++)
           {
               X0 += crd[3*j];
               Y0 += crd[3*j+1];
               Z0 += crd[3*j+2];
           } 
           
           X0 /= ncrd;
           Y0 /= ncrd;
           Z0 /= ncrd;
           free(crd);
           
           (*sprst)->springs[i].X0 = X0;
           (*sprst)->springs[i].Y0 = Y0;
           (*sprst)->springs[i].Z0 = Z0;
           
           if(fgets(buffer, linesz-1, fp) == NULL)
           {
                 printf("Error reading springs");
                 exit(0);
           }
           
           f = atof(buffer);
           (*sprst)->springs[i].fkspr = f;
        }        
        free(buffer);
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
	double xtot, ytot, ztot, fk;
	struct mol_vector3 g;
	struct mol_atom_group *ag;

	for (i = 0; i < sprst->nsprings; i++)
	{
		nat = sprst->springs[i].naspr;
		
		if(nat > 0)
		{
			xtot = 0.0;
			ytot = 0.0;
			ztot = 0.0;
			ag = sprst->springs[i].ag;

			for (i1 = 0; i1 < nat; i1++)
			{
				i2 = sprst->springs[i].laspr[i1];
				xtot += ag->coords[i2].X;
				ytot += ag->coords[i2].Y;
				ztot += ag->coords[i2].Z;
			}

			xtot = xtot / nat - sprst->springs[i].X0;
			ytot = ytot / nat - sprst->springs[i].Y0;
			ztot = ztot / nat - sprst->springs[i].Z0;

			fk = sprst->springs[i].fkspr;
			(*een) += fk * (xtot * xtot + ytot * ytot + ztot * ztot);

			fk = 2 * fk / nat;
			g.X = xtot*fk;
			g.Y = ytot*fk;
			g.Z = ztot*fk;
			
			for (i1 = 0; i1 < nat; i1++)
			{
				i2 = sprst->springs[i].laspr[i1];
				MOL_VEC_SUB(ag->gradients[i2], ag->gradients[i2], g);
			}
		}
	}
}
