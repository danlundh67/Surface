#include "structures.h"


struct Structure create_xyzr_structure( struct Structure This_Structure ) {

/************/

  /* Variables */
  struct Structure	New_Structure ;
  char atm[5];
  char resnam[4];

  /* Counters */
  int		residue , atom ;

/************/

  if( ( New_Structure.Residue = ( struct Amino_Acid * ) malloc ( ( This_Structure.length + 1 ) * sizeof_Amino_Acid ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }

  strcpy( New_Structure.ident , This_Structure.ident ) ;
  New_Structure.length = This_Structure.length ;

  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {

    New_Structure.Residue[residue] = This_Structure.Residue[residue] ;

    strcpy(resnam,New_Structure.Residue[residue].res_name);
    
    if( ( New_Structure.Residue[residue].Atom = ( struct Atom * ) malloc ( ( This_Structure.Residue[residue].size + 1 ) * sizeof_Atom ) ) == NULL ) {
      GENERAL_MEMORY_PROBLEM
    }

    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {

      New_Structure.Residue[residue].Atom[atom] = This_Structure.Residue[residue].Atom[atom] ;
      strcpy(atm,New_Structure.Residue[residue].Atom[atom].atom_name);



      if (strcmp(atm," C  ")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.74;
      else if (strcmp(atm," CA ")==0)
	New_Structure.Residue[residue].Atom[atom].charge=2.00;
      else if (strcmp(atm," CB ")==0)
	New_Structure.Residue[residue].Atom[atom].charge=2.00;
      else if ((strcmp(atm," CD ")==0)&&(strcmp(resnam,"GLN")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.74;
      else if ((strcmp(atm," CD ")==0)&&(strcmp(resnam,"GLU")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.74;
      else if (strcmp(atm," CD ")==0)
	New_Structure.Residue[residue].Atom[atom].charge=2.00;
      else if ((strcmp(atm," CD1")==0)&&(strcmp(resnam,"PHE")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.86;
      else if ((strcmp(atm," CD1")==0)&&(strcmp(resnam,"TYR")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.86;
      else if ((strcmp(atm," CD1")==0)&&(strcmp(resnam,"TRP")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.86;
      else if ((strcmp(atm," CD1")==0))
	New_Structure.Residue[residue].Atom[atom].charge=2.00;
      else if ((strcmp(atm," CD2")==0)&&(strcmp(resnam,"TRP")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.74;
      else if ((strcmp(atm," CD2")==0)&&(strcmp(resnam,"PHE")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.86;
      else if ((strcmp(atm," CD2")==0)&&(strcmp(resnam,"HIS")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.86;
      else if ((strcmp(atm," CD2")==0)&&(strcmp(resnam,"TYR")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.86;
      else if ((strcmp(atm," CD2")==0))
	New_Structure.Residue[residue].Atom[atom].charge=2.00;
      else if ((strcmp(atm," CE ")==0))
	New_Structure.Residue[residue].Atom[atom].charge=2.00;
      else if ((strcmp(atm," CE1")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.86;
      else if ((strcmp(atm," CE2")==0)&&(strcmp(resnam,"TRP")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.74;
      else if ((strcmp(atm," CE2")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.86;
      else if ((strcmp(atm," CE3")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.86;
      else if ((strcmp(atm," CG ")==0)&&(strcmp(resnam,"ASN")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.74;
      else if ((strcmp(atm," CG ")==0)&&(strcmp(resnam,"PHE")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.74;
      else if ((strcmp(atm," CG ")==0)&&(strcmp(resnam,"ASP")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.74;
      else if ((strcmp(atm," CG ")==0)&&(strcmp(resnam,"LEU")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.74;
      else if ((strcmp(atm," CG ")==0)&&(strcmp(resnam,"TYR")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.74;
      else if ((strcmp(atm," CG ")==0)&&(strcmp(resnam,"HIS")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.74;
      else if ((strcmp(atm," CG ")==0)&&(strcmp(resnam,"TRP")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.74;
      else if ((strcmp(atm," CG ")==0))
	New_Structure.Residue[residue].Atom[atom].charge=2.00;
      else if ((strcmp(atm," CG1")==0))
	New_Structure.Residue[residue].Atom[atom].charge=2.00;
      else if ((strcmp(atm," CG2")==0))
	New_Structure.Residue[residue].Atom[atom].charge=2.00;
      else if ((strcmp(atm," CH2")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.86;
      else if ((strcmp(atm," CH3")==0))
	New_Structure.Residue[residue].Atom[atom].charge=2.00;
     else if ((strcmp(atm," CZ ")==0)&&(strcmp(resnam,"TYR")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.74;
     else if ((strcmp(atm," CZ ")==0)&&(strcmp(resnam,"ARG")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.74;
     else if ((strcmp(atm," CZ ")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.86;
     else if ((strcmp(atm," CZ2")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.86;
     else if ((strcmp(atm," CZ3")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.86;
     else if (strcmp(atm," N  ")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.70;
     else if (strcmp(atm," ND1")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.54;
     else if (strcmp(atm," ND2")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.80;
     else if (strcmp(atm," NE ")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.70;
     else if (strcmp(atm," NE1")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.70;
     else if ((strcmp(atm," NE2")==0)&&(strcmp(resnam,"HIS")==0))
	New_Structure.Residue[residue].Atom[atom].charge=1.70;
     else if (strcmp(atm," NE2")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.80;
     else if (strcmp(atm," NH1")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.80;
     else if (strcmp(atm," NH2")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.80;
     else if (strcmp(atm," NZ ")==0)
	New_Structure.Residue[residue].Atom[atom].charge=2.00;
     else if (strcmp(atm," O  ")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.40;
     else if (strcmp(atm," OD1")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.40;
     else if (strcmp(atm," OD2")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.40;
     else if (strcmp(atm," OE1")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.40;
     else if (strcmp(atm," OE2")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.40;
     else if (strcmp(atm," OG ")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.60;
     else if (strcmp(atm," OG1")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.60;
     else if (strcmp(atm," OH ")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.60;
     else if (strcmp(atm," SD ")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.80;
     else if (strcmp(atm," SG ")==0)
	New_Structure.Residue[residue].Atom[atom].charge=1.80;
      else 
	New_Structure.Residue[residue].Atom[atom].charge=1.80;

      //      printf ("%lf \n",	New_Structure.Residue[residue].Atom[atom].charge);
    }

  }


  return New_Structure ;

/************/

}

