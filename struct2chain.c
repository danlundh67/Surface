#include "structures.h"


struct the_chain_report structure_chains( struct Structure This_Structure ) {

/************/

  /* Counters */
  int		residue , atom, newresid, nrch;
  char chain[2];
  int *chainnr,r,i,j,jpocc,jpocc2,vara;
  struct the_chain_report myreport;



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
  //  printf(" nr of chains %d \n",nrch);  
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
  
  myreport.nrch=nrch;
  strcpy(myreport.ident,This_Structure.ident);
  myreport.subchains= ( struct chainsrep *)malloc ((nrch*sizeof(struct chainsrep)));


  for( r = 0 ; r < (nrch) ; r++ ) {      
    /* REPORT STUFF */
    myreport.subchains[r].chain_length=chainnr[r+1]-chainnr[r];
    strcpy(myreport.subchains[r].chainID,This_Structure.Residue[chainnr[r]].chainID);
  }

  free(chainnr);
  return (myreport);
  free(myreport.subchains);

/************/

}

