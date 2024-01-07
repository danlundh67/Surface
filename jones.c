#include "structures.h"
#include<math.h>

/* struct JONDAT{
  char res_name[4];
  int occur[50];
  int sum;
  }; */

struct JONDAT jp;
struct JONDAT jp2[50];

void read_jondata()
{
  FILE *jf;
  int i,j,x;
  if( ( jf= fopen("m2res2.dat", "r" ) ) == NULL ) {
    printf( "The jones data #1 file could not be opened.\nDying\n\n" ) ;
    exit(  EXIT_FAILURE ) ;
  }
  jp.sum=0;
  for (i=1;i<50;i++)
    {
      fscanf(jf,"%d ",&jp.occur[i-1]);
      jp.sum=jp.sum+jp.occur[i-1];
    }
  /*    for (i=1;i<50;i++)
    printf("%d ",jp.occur[i-1]);
    printf("%d \n",jp.sum);
  */
  fclose(jf);
  if( ( jf= fopen("m2res.dat", "r" ) ) == NULL ) {
    printf( "The jones data #2 file could not be opened.\nDying\n\n" ) ;
    exit(  EXIT_FAILURE ) ;
  }
  x=0;
  while (feof(jf)==0)
    {
      jp2[x].sum=0;
      fscanf(jf,"%s ",jp2[x].res_name);
      for (i=1;i<50;i++)
	{
	  fscanf(jf,"%d ",&jp2[x].occur[i-1]);
	  jp2[x].sum=jp2[x].sum+jp2[x].occur[i-1];
	}
      x++;
    }
  /*   for (i=0;i<20;i++)
    {
      printf("%s ",jp2[i].res_name);
      for (j=0;j<49;j++)
        printf("%d ",jp2[i].occur[j]);
      printf("\n");
      printf("%d\n",jp2[i].sum);
    }
  */
  fflush(jf);
  fclose(jf);

}

double mymain(struct Structure chn1)
{
    int energy, offset;
    struct Structure chn2;
    struct Amino_Acid group1;
    struct Amino_Acid group2;
    struct Atom cb1;
    struct Atom cb2;
    struct Atom pc1;
    struct Atom po1;
    struct Atom n1;
    char myaa[4];
    int pres,nrch;
    int pos1,pos2,hit;
    int dx,dy,dz;
    double dco;
    long dist;
    double sumH;
    int atom,atom2, residue, residue2;
    char chain[2];
    int *chainnr,r,i,j,jpocc,jpocc2,vara;
    double sumJP,pJP;
    hit=0;
    pos1=0;
    sumH=0;
    nrch=1;
    vara=-1;

    /* READ THE JONES OCCURENCES */

    read_jondata();
    sumJP=0.0000;

    /* PREPROCESSING FOR ID OF DIFF CHAINS */

    strcpy(chain,chn1.Residue[1].chainID);
    for (residue = 1;  residue <= (chn1.length) ; residue ++ ) {
      if (strcmp(chain,chn1.Residue[residue].chainID)!=0)
	{
	  /* printf(" chains %s \n",chain); */
	  strcpy(chain,chn1.Residue[residue].chainID);
	  nrch++;
	}
      strcpy(chain,chn1.Residue[residue].chainID);
    }
    /*    printf(" nr of chains %d \n",nrch); */
    chainnr = (int *)malloc(sizeof(int)*(nrch+1));
    strcpy(chain,chn1.Residue[1].chainID);
    chainnr[0]=1;
    r=1;
    

    for (residue = 1;  residue <= (chn1.length) ; residue ++ ) {
      if (strcmp(chain,chn1.Residue[residue].chainID)!=0)
	{
	  strcpy(chain,chn1.Residue[residue].chainID);
	  chainnr[r]=residue;
	  r++;
	}
      strcpy(chain,chn1.Residue[residue].chainID);
    }
    chainnr[r]=chn1.length+1;
    /*  for (r = 0;  r < (nrch) ; r++ ) {
      printf(" chains %d %d-%d \n",r+1,chainnr[r],chainnr[r+1]);
      } */

    chn2=duplicate_structure(chn1);
    for( r = 0 ; r < (nrch) ; r++ ) {      
      hit=0;
      pos1=0;
      sumH=0;
      sumJP=0.0000;

      /* Fix so that DNA chains is skipped */
      /* testing first AA - if nucleotide then skip */

      if ((strcmp(chn1.Residue[chainnr[r]].res_name,"  C")==0)||
	  (strcmp(chn1.Residue[chainnr[r]].res_name,"  A")==0)||
	  (strcmp(chn1.Residue[chainnr[r]].res_name,"  G")==0)||
	  (strcmp(chn1.Residue[chainnr[r]].res_name,"  T")==0))
	{
	  
	  printf("%s-%s DNA  \n",chn1.ident,chn1.Residue[chainnr[r]].chainID);
	  continue;
	}
      for( residue = chainnr[r] ; ((residue < (chainnr[r+1]))) ; residue ++ ) { 
	hit=0;
	pos1=0;
	for( atom = 1 ; atom <= chn1.Residue[residue].size ; atom ++ ) {     
	  if (strcmp(chn1.Residue[residue].Atom[atom].atom_name," CB ")==0)
	    {
	      cb1.coord[1]=chn1.Residue[residue].Atom[atom].coord[1];
	      cb1.coord[2]=chn1.Residue[residue].Atom[atom].coord[2];
	      cb1.coord[3]=chn1.Residue[residue].Atom[atom].coord[3];
	      cb1.serial=chn1.Residue[residue].Atom[atom].serial;
	      hit=1;
	    }
	}
	if (hit==1) {
	  for( residue2 = chainnr[r] ; ((residue2 < (chainnr[r+1]))) ; residue2 ++ ) {
	    if (residue2==residue)
	      {

	      }
	    else
	      {
		for( atom = 1 ; atom <= chn2.Residue[residue2].size ; atom ++ ) {     
		  if (strcmp(chn2.Residue[residue2].Atom[atom].atom_name," CB ")==0)
		    {
		      cb2.coord[1]=chn2.Residue[residue2].Atom[atom].coord[1];
		      cb2.coord[2]=chn2.Residue[residue2].Atom[atom].coord[2];
		      cb2.coord[3]=chn2.Residue[residue2].Atom[atom].coord[3];
		      cb2.serial=chn2.Residue[residue2].Atom[atom].serial;	      
		    }
		}
		if ((pow((cb1.coord[1]-cb2.coord[1]),2)+pow((cb1.coord[2]-cb2.coord[2]),2)+pow((cb1.coord[3]-cb2.coord[3]),2))<100.0)
		  {
		    /* printf("%s-%d ",chn1.Residue[residue2].res_name,residue2); */
		    pos1++;
		  }
	      }
	  }
	  jpocc=jp.occur[pos1-1];
	  for (i=0;i<20;i++)
            {
              if(strcmp(jp2[i].res_name,chn1.Residue[residue].res_name)==0)
                {
                  jpocc2=jp2[i].occur[pos1-1];
		  vara=i;
                }
            }


	  if (jpocc>0)
	    {
	      /*   printf("%d ",r);
	    printf("%s ",chn1.Residue[residue].res_name);
	    printf("%s ",chn1.Residue[residue].res_seq_plus_iCode);
	    printf("%d  %d %d ",pos1,jpocc,jpocc2);
	    printf("%lf ",(float)jpocc2/jp2[vara].sum);
	    printf("%lf ",(float)jpocc/jp.sum);
	    printf("%lf \n",-1.000*0.582*log(((double)((float)jpocc2/jp2[vara].sum)/((float)jpocc/jp.sum)))); */
	      if (vara!=-1)
		sumJP=sumJP+-1.000*0.582*log(((double)((float)jpocc2/jp2[vara].sum)/((float)jpocc/jp.sum)));
	    }
	}

      }
    
       if (nrch==1)
	{
	  /* if (strlen(chn1.Residue[1].chainID)==0)
	   printf("JP %s A %d  %lf \n",chn1.ident,chainnr[r+1]-chainnr[r],sumJP);	 
	  else 
	    printf("JP %s %s %d %lf \n",chn1.ident,chn1.Residue[1].chainID,chn1.length,sumJP);
	  */
	  return sumJP;
	}
      else
	// printf("JP %s %s %d  %lf \n",chn1.ident,chn1.Residue[chainnr[r]].chainID,chainnr[r+1]-chainnr[r],sumJP);
	return sumJP;
    }
  
}
