#include "structures.h"


struct the_new_report structure_to_chain( struct Structure This_Structure,int chainnumber ) {

/************/

  /* Variables */
  struct Structure	New_Structure ;

  /* Counters */
  int		residue , atom, newresid, nrch;
  char chain[2];
  int *chainnr,r,i,j,jpocc,jpocc2,vara;
  struct the_new_report myreport;



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


 if( ( New_Structure.Residue = ( struct Amino_Acid * ) malloc ( ( chainnr[chainnumber+1]-chainnr[chainnumber] ) * sizeof_Amino_Acid ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
	}
  strcpy( New_Structure.ident , This_Structure.ident ) ;
  New_Structure.length = chainnr[chainnumber+1]-chainnr[chainnumber];
  newresid=1;

  myreport.subchains.m.secstruct= (char *) malloc ((chainnr[chainnumber+1]-chainnr[chainnumber])*sizeof(char)+1);
    
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
  


    strcpy(myreport.subchains.chainID,This_Structure.Residue[chainnr[chainnumber]].chainID);
    myreport.subchains.volsur=calc_volume(New_Structure,"",0.5,1.4);

    printf("Calc volume/surface done \n");

    myreport.subchains.m=CalcProteinHBonds(New_Structure,NULL );

    printf("Calc hydrogen bonds done \n");

    myreport.subchains.HbondE=mymain(New_Structure);


  return (myreport);

  free(myreport.subchains.m.secstruct);

/************/

}

