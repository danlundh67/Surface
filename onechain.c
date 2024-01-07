#include "structures.h"


struct Structure structure_onechain( struct Structure This_Structure, int chainnumber) {

/************/

  /* Variables */
  struct Structure	New_Structure, surface;

  struct meashure slask;

  /* Counters */
  int		residue , atom, newresid, nrch, sres,satom;
  char chain[2];
  int *chainnr,r,i,j,jpocc,jpocc2,vara;
  char *revalue;
  float x,y,z, mx,my,mz,rad;

/************/

  nrch=1;
  strcpy(chain,This_Structure.Residue[1].chainID);
  for (residue = 1;  residue <= (This_Structure.length) ; residue ++ ) {
    if (strcmp(chain,This_Structure.Residue[residue].chainID)!=0)
      {
	/* printf(" chains %s \n",chain); */
	strcpy(chain,This_Structure.Residue[residue].chainID);
	nrch++;
      }
    strcpy(chain,This_Structure.Residue[residue].chainID);
  }
  /*  printf(" nr of chains %d \n",nrch);  */
  chainnr = (int *)malloc(sizeof(int)*(nrch+1));
  strcpy(chain,This_Structure.Residue[1].chainID);
  chainnr[0]=1;
  r=1;
  for (residue = 1;  residue <= (This_Structure.length) ; residue ++ ) {
    if (strcmp(chain,This_Structure.Residue[residue].chainID)!=0)
      {
	strcpy(chain,This_Structure.Residue[residue].chainID);
	chainnr[r]=residue;
	r++;
      }
    strcpy(chain,This_Structure.Residue[residue].chainID);
  }
  chainnr[r]=This_Structure.length+1;
  


  if( ( New_Structure.Residue = ( struct Amino_Acid * ) malloc ( ( chainnr[chainnumber+1]-chainnr[chainnumber] ) * sizeof_Amino_Acid ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
	}
  strcpy( New_Structure.ident , This_Structure.ident ) ;
  New_Structure.length = chainnr[chainnumber+1]-chainnr[chainnumber];
  newresid=1;
    
  for( residue = chainnr[chainnumber] ; ((residue < (chainnr[chainnumber+1]))) ; residue ++ ) { 
      if ((New_Structure.Residue[newresid].Atom= (struct Atom *)malloc((This_Structure.Residue[residue].size + 1)*sizeof_Atom))==NULL) {
	GENERAL_MEMORY_PROBLEM
	  }
      
      New_Structure.Residue[newresid] = This_Structure.Residue[residue] ;
      for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {
	New_Structure.Residue[newresid].Atom[atom] = This_Structure.Residue[residue].Atom[atom] ;
      }
      newresid++;	
    }
  
  slask=calc_volume(New_Structure,"True",1.4, 1.4);
  revalue=(char *)malloc(sizeof(char)*(strlen(slask.tmppdb)+1));
  strcpy(revalue,slask.tmppdb);

  surface=read_pdb_to_structure(revalue);

  free(revalue);
  for( residue = 1 ; residue <= New_Structure.length ; residue ++ ) {
    for( atom = 1 ; atom <= New_Structure.Residue[residue].size ; atom ++ ) {
      New_Structure.Residue[residue].Atom[atom].surface=0;
      mx= New_Structure.Residue[residue].Atom[atom].coord[1];
      my= New_Structure.Residue[residue].Atom[atom].coord[2];
      mz= New_Structure.Residue[residue].Atom[atom].coord[3];
      rad= New_Structure.Residue[residue].Atom[atom].charge;
      for(sres=1;sres <= surface.length;sres++)
	{
	  for( satom = 1 ; satom <= surface.Residue[sres].size ; satom ++ ) {
	    x=surface.Residue[sres].Atom[satom].coord[1];
	    y=surface.Residue[sres].Atom[satom].coord[2];
	    z=surface.Residue[sres].Atom[satom].coord[3];
	    if (pow(((mx-x)*(mx-x)+(my-y)*(my-y)+(mz-z)*(mz-z)),0.5)<(1.4+2*rad))
	      New_Structure.Residue[residue].Atom[atom].surface=1;
	  }
	}
    }
  }
  /* write_dump_pdb(New_Structure,"surface22.pdb"); */


  return (New_Structure);  
/************/

}



