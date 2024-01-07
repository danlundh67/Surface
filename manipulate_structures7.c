/*
This file is part of ftdock, a program for rigid-body protein-protein docking 
Copyright (C) 1997-2000 Gidon Moont

Biomolecular Modelling Laboratory
Imperial Cancer Research Fund
44 Lincoln's Inn Fields
London WC2A 3PX

+44 (0)20 7269 3348
http://www.bmm.icnet.uk/

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include "structures.h"

struct Structure read_pdb_to_structure( char *pdb_file_name ) {

/************/

  /* Variables */

  /* Counters */
  int	n_residues ;	/* number of residues */
  int	res_size ;	/* number of atoms in single residue */

  /* File stuff */
  FILE	*pdb_file ;
  char	line_buffer[100] ;

  /* What the data is going into */
  struct Structure		This_Structure ;

  /* Variables from the PDB file */
  int	serial ;
  char		atom_name[5] ;
  char		res_name[4] ;
  char		chainID[2] ;
  char		res_seq_plus_iCode[6] ;
  float		coord_x , coord_y , coord_z ;
  float		occupancy, temp_factor ;
  char		olc[2] ;
  int		nc ;

  /* Comparison values */
  char	present_res_seq_plus_iCode[6] ;

/************/

  setvbuf( stdout , (char *)NULL , _IONBF , 0 ) ;

  /* File handling */

  /* Open file */
  /*  printf( " %s", pdb_file_name ) ; */
  if( ( pdb_file = fopen( pdb_file_name, "r" ) ) == NULL ) {
    printf( "This file does not exist here, or is unreadable.\nDying\n\n" ) ;
    exit( EXIT_FAILURE ) ;
  }

/************/

  /* Initialisations */

  /* Counters */
  n_residues = 0 ;
  res_size = 0 ;

  MinZ = MinY = MinX = 90000;
  MaxX= MaxY= MaxZ = -90000;
  

  /* Comparison values */
  strcpy( present_res_seq_plus_iCode , ">" ) ;

  /* Memory allocation */
  if( ( This_Structure.Residue = ( struct Amino_Acid * ) malloc ( sizeof_Amino_Acid ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }

/************/

  /* Read PDB file */

  /* The Atoms */

  while( fgets( line_buffer, 85, pdb_file ) ) {

    if( strncmp( line_buffer, "ATOM", 4 ) == 0 ) {

      /* Have an ATOM */

      /* Get Values */

      /* the following may seem silly, but sscanf convention means that two
         float fields with no white space between them, where the first is
         less than the maximum field width, mucks up everything.
      */

      sscanf( line_buffer +  6 , "%5d" , &serial ) ;
      sscanf( line_buffer + 30 , "%8f" , &coord_x ) ;
      sscanf( line_buffer + 38 , "%8f" , &coord_y ) ;
      sscanf( line_buffer + 46 , "%8f" , &coord_z ) ;
      sscanf( line_buffer + 54 , "%6f" , &occupancy ) ;
      sscanf( line_buffer + 60 , "%6f" , &temp_factor ) ;
      sscanf( line_buffer + 82 , "%2d" , &nc ) ;

      strncpy( atom_name,		line_buffer+12,	4 ) ;
      strncpy( res_name,		line_buffer+17,	3 ) ;
      strncpy( chainID,			line_buffer+21,	1 ) ;
      strncpy( res_seq_plus_iCode,	line_buffer+22,	5 ) ;
      strncpy( olc,			line_buffer+80,	1 ) ;

      strncpy( atom_name + 4,		"\0", 1 ) ;
      strncpy( res_name + 3,		"\0", 1 ) ;
      strncpy( chainID + 1,		"\0", 1 ) ;
      strncpy( res_seq_plus_iCode + 5,	"\0", 1 ) ;
      strncpy( olc + 1,			"\0", 1 ) ;

/************/

      /* New Residue */

      if( strcmp( res_seq_plus_iCode , present_res_seq_plus_iCode ) != 0 ) {

        /* have next residue */

        /* Store old info */
        This_Structure.Residue[n_residues].size = res_size ;

        /* Increment, Reset numbers */
        n_residues ++ ;
        res_size = 0 ;

        /* Memory management */
        if( ( This_Structure.Residue = (struct Amino_Acid * ) realloc ( This_Structure.Residue, ( n_residues + 1 ) * sizeof_Amino_Acid ) ) == NULL ) {
          GENERAL_MEMORY_PROBLEM
        }
        if( ( This_Structure.Residue[n_residues].Atom = ( struct Atom * ) malloc ( sizeof_Atom ) ) == NULL ) {
          GENERAL_MEMORY_PROBLEM
        }

        /* Store new info */
        strcpy( This_Structure.Residue[n_residues].res_seq_plus_iCode , res_seq_plus_iCode );
        strcpy( This_Structure.Residue[n_residues].res_name ,           res_name ) ;
        strcpy( This_Structure.Residue[n_residues].chainID ,            chainID ) ;
        strcpy( This_Structure.Residue[n_residues].olc,                 olc ) ;
        This_Structure.Residue[n_residues].nc = nc ;

      }

      strcpy( present_res_seq_plus_iCode , res_seq_plus_iCode ) ;

/************/

      /* Put Atoms into Structure */

      res_size ++ ;

      if( ( This_Structure.Residue[n_residues].Atom = ( struct Atom * ) realloc ( This_Structure.Residue[n_residues].Atom, ( res_size + 1 ) * sizeof_Atom ) ) == NULL ) {
        GENERAL_MEMORY_PROBLEM
      }

      This_Structure.Residue[n_residues].Atom[res_size].serial = serial ;
      strcpy( This_Structure.Residue[n_residues].Atom[res_size].atom_name, atom_name ) ;
      This_Structure.Residue[n_residues].Atom[res_size].coord[1] = coord_x ;
      This_Structure.Residue[n_residues].Atom[res_size].coord[2] = coord_y ;
      This_Structure.Residue[n_residues].Atom[res_size].coord[3] = coord_z ;
      This_Structure.Residue[n_residues].Atom[res_size].occupancy = occupancy ;
      This_Structure.Residue[n_residues].Atom[res_size].temp_factor = temp_factor ;
      This_Structure.Residue[n_residues].Atom[res_size].surface=1;
      if (coord_x<MinX)
	MinX=coord_x;
      else if (coord_x>MaxX)
	MaxX=coord_x;
      if (coord_y<MinY)
	MinY=coord_y;
      else if (coord_y>MaxY)
	MaxY=coord_y;
      if (coord_z<MinZ)
	MinZ=coord_z;
      else if (coord_z>MaxZ)
	MaxZ=coord_z;
      /*  g_print("%s %s %s %d %d %f %f %f\n",This_Structure.Residue[n_residues].res_name,This_Structure.Residue[n_residues].res_seq_plus_iCode ,This_Structure.Residue[n_residues].Atom[res_size].atom_name,This_Structure.Residue[n_residues].Atom[res_size].serial,This_Structure.Residue[n_residues].Atom[res_size].surface,This_Structure.Residue[n_residues].Atom[res_size].coord[1],This_Structure.Residue[n_residues].Atom[res_size].coord[2],This_Structure.Residue[n_residues].Atom[res_size].coord[3]); 	 */



/************/

    }


  } /* got to end of pdb file */

/************/

  /* Clean up */

  This_Structure.Residue[n_residues].size = res_size ;
  This_Structure.length = n_residues ;
  strcpy( This_Structure.ident , pdb_file_name  );

  /* Finish off */

  fclose( pdb_file ) ;
  /*  printf("x %f %f y %f %f z %f %f\n",MinX, MaxX,MinY,MaxY,MinZ, MaxZ); */
  return This_Structure ;

}



/************************/



void write_structure_to_pdb( struct Structure This_Structure , char *pdb_file_name ) {

/************/

  /* Variables */

  /* Counters */
  int	residue , atom ;

  /* File stuff */
  FILE		*pdb_file ;

/************/

  /* File handling */

  /* Open file */
  printf( "Writing file: %s\n", pdb_file_name ) ;
  if( ( pdb_file = fopen( pdb_file_name, "w" ) ) == NULL ) {
    printf( "This file could not be opened.\nDying\n\n" ) ;
    exit(  EXIT_FAILURE ) ;
  }

/************/

  /* Write PDB file */

  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {

    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {

      fprintf( pdb_file, "ATOM  %5d %4s %3s %1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f              %1s %2d\n", This_Structure.Residue[residue].Atom[atom].serial, This_Structure.Residue[residue].Atom[atom].atom_name, This_Structure.Residue[residue].res_name, This_Structure.Residue[residue].chainID, This_Structure.Residue[residue].res_seq_plus_iCode, This_Structure.Residue[residue].Atom[atom].coord[1], This_Structure.Residue[residue].Atom[atom].coord[2], This_Structure.Residue[residue].Atom[atom].coord[3], This_Structure.Residue[residue].Atom[atom].occupancy, This_Structure.Residue[residue].Atom[atom].temp_factor, This_Structure.Residue[residue].olc, This_Structure.Residue[residue].nc ) ;

    }

  }

/************/

  /* Finish off */

  fclose( pdb_file ) ;

}


void write_dump_pdb( struct Structure This_Structure , char *pdb_file_name ) {

/************/

  /* Variables */

  /* Counters */
  int	residue , atom ;

  /* File stuff */
  FILE		*pdb_file ;

/************/

  /* File handling */

  /* Open file */
  printf( "Writing file: %s\n", pdb_file_name ) ;
  if( ( pdb_file = fopen( pdb_file_name, "w" ) ) == NULL ) {
    printf( "This file could not be opened.\nDying\n\n" ) ;
    exit(  EXIT_FAILURE ) ;
  }

/************/

  /* Write PDB file */

  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {

    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {
      if (This_Structure.Residue[residue].Atom[atom].surface==1)
      fprintf( pdb_file, "ATOM  %5d %4s %3s %1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f              %1s %2d\n", This_Structure.Residue[residue].Atom[atom].serial, This_Structure.Residue[residue].Atom[atom].atom_name, This_Structure.Residue[residue].res_name, This_Structure.Residue[residue].chainID, This_Structure.Residue[residue].res_seq_plus_iCode, This_Structure.Residue[residue].Atom[atom].coord[1], This_Structure.Residue[residue].Atom[atom].coord[2], This_Structure.Residue[residue].Atom[atom].coord[3], This_Structure.Residue[residue].Atom[atom].occupancy, This_Structure.Residue[residue].Atom[atom].temp_factor, This_Structure.Residue[residue].olc, This_Structure.Residue[residue].nc ) ;

    }

  }

/************/

  /* Finish off */

  fclose( pdb_file ) ;

}





/************************/



struct Structure duplicate_structure( struct Structure This_Structure ) {

/************/

  /* Variables */
  struct Structure	New_Structure ;

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

    if( ( New_Structure.Residue[residue].Atom = ( struct Atom * ) malloc ( ( This_Structure.Residue[residue].size + 1 ) * sizeof_Atom ) ) == NULL ) {
      GENERAL_MEMORY_PROBLEM
    }

    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {

      New_Structure.Residue[residue].Atom[atom] = This_Structure.Residue[residue].Atom[atom] ;

    }

  }

  return New_Structure ;

/************/

}



/************************/



struct Structure translate_structure( struct Structure This_Structure , float x_shift , float y_shift , float z_shift ) {

/************/

  /* Variables */
  struct Structure	New_Structure ;

  /* Counters */
  int		residue , atom ;

/************/

  New_Structure = duplicate_structure( This_Structure ) ;

/************/

  for( residue = 1 ; residue <= New_Structure.length ; residue ++ ) {

    for( atom = 1 ; atom <= New_Structure.Residue[residue].size ; atom ++ ) {

      New_Structure.Residue[residue].Atom[atom].coord[1] += x_shift ;
      New_Structure.Residue[residue].Atom[atom].coord[2] += y_shift ;
      New_Structure.Residue[residue].Atom[atom].coord[3] += z_shift ;

    }

  }

  return New_Structure ;

/************/

}



/************************/



struct Structure translate_structure_onto_origin( struct Structure This_Structure ) {

/************/

  /* Variables */
  struct Structure	New_Structure ;

  float			average_x , average_y , average_z ;

  /* Counters */
  int		residue , atom , total_atoms ;

/************/

  New_Structure = duplicate_structure( This_Structure ) ;

/************/

  /* Find current centre */

  total_atoms = 0 ;

  average_x = 0 ;
  average_y = 0 ;
  average_z = 0 ;

  for( residue = 1 ; residue <= New_Structure.length ; residue ++ ) {

    for( atom = 1 ; atom <= New_Structure.Residue[residue].size ; atom ++ ) {

      total_atoms ++ ;

      average_x += New_Structure.Residue[residue].Atom[atom].coord[1] ;
      average_y += New_Structure.Residue[residue].Atom[atom].coord[2] ;
      average_z += New_Structure.Residue[residue].Atom[atom].coord[3] ;

    }

  }

  average_x = average_x / (float)total_atoms ;
  average_y = average_y / (float)total_atoms ;
  average_z = average_z / (float)total_atoms ;

/************/

  /* Translate */

  for( residue = 1 ; residue <= New_Structure.length ; residue ++ ) {

    for( atom = 1 ; atom <= New_Structure.Residue[residue].size ; atom ++ ) {

      New_Structure.Residue[residue].Atom[atom].coord[1] -= average_x ;
      New_Structure.Residue[residue].Atom[atom].coord[2] -= average_y ;
      New_Structure.Residue[residue].Atom[atom].coord[3] -= average_z ;

    }

  }

  return New_Structure ;

/************/

}



/************************/



struct Structure rotate_structure( struct Structure This_Structure , int z_twist , int theta , int phi ) {

/************/

  /* Variables */
  struct Structure	New_Structure ;

  float			post_z_twist_x , post_z_twist_y , post_z_twist_z ;
  float			post_theta_x , post_theta_y , post_theta_z ;

  /* Counters */
  int		residue , atom ;

/************/

  New_Structure = duplicate_structure( This_Structure ) ;

/************/

  for( residue = 1 ; residue <= New_Structure.length ; residue ++ ) {

    for( atom = 1 ; atom <= New_Structure.Residue[residue].size ; atom ++ ) {

      /* Perform Z axis twist */
      post_z_twist_x = New_Structure.Residue[residue].Atom[atom].coord[1] * cos( 0.017453293 * z_twist ) - New_Structure.Residue[residue].Atom[atom].coord[2] * sin( 0.017453293 * z_twist ) ;
      post_z_twist_y = New_Structure.Residue[residue].Atom[atom].coord[1] * sin( 0.017453293 * z_twist ) + New_Structure.Residue[residue].Atom[atom].coord[2] * cos( 0.017453293 * z_twist ) ;
      post_z_twist_z = New_Structure.Residue[residue].Atom[atom].coord[3] ;

      /* Perform theta twist along plane of x-z */
      post_theta_x = post_z_twist_z * sin( 0.017453293 * theta ) + post_z_twist_x * cos( 0.017453293 * theta ) ; 
      post_theta_y = post_z_twist_y ;
      post_theta_z = post_z_twist_z * cos( 0.017453293 * theta ) - post_z_twist_x * sin( 0.017453293 * theta ) ; 

      /* Perform phi twist around z axis */
      New_Structure.Residue[residue].Atom[atom].coord[1] = post_theta_x * cos( 0.017453293 * phi ) - post_theta_y * sin( 0.017453293 * phi ) ;
      New_Structure.Residue[residue].Atom[atom].coord[2] = post_theta_x * sin( 0.017453293 * phi ) + post_theta_y * cos( 0.017453293 * phi ) ;
      New_Structure.Residue[residue].Atom[atom].coord[3] = post_theta_z ;

    }

  }

  return New_Structure ;

/************/

}



/************************/



struct Structure merge_structures( struct Structure Structure_One , struct Structure Structure_Two ) {

/************/

  /* Variables */
  struct Structure	New_Structure ;

  /* Counters */
  int		residue , atom , new_residue ;

/************/

  if( ( New_Structure.Residue = ( struct Amino_Acid * ) malloc ( ( Structure_One.length + Structure_Two.length + 1 ) * sizeof_Amino_Acid ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }

  strcpy( New_Structure.ident , "Complex" ) ;
  New_Structure.length = Structure_One.length + Structure_Two.length ;

  for( residue = 1 ; residue <= Structure_One.length ; residue ++ ) {

    if( ( New_Structure.Residue[residue].Atom = ( struct Atom * ) malloc ( ( Structure_One.Residue[residue].size + 1 ) * sizeof_Atom ) ) == NULL ) {
      GENERAL_MEMORY_PROBLEM
    }
    strcpy( New_Structure.Residue[residue].res_name           , Structure_One.Residue[residue].res_name ) ;
    strcpy( New_Structure.Residue[residue].chainID            , Structure_One.Residue[residue].chainID ) ;
    strcpy( New_Structure.Residue[residue].res_seq_plus_iCode , Structure_One.Residue[residue].res_seq_plus_iCode ) ;
    strcpy( New_Structure.Residue[residue].olc                , Structure_One.Residue[residue].olc ) ;
    New_Structure.Residue[residue].nc                         = Structure_One.Residue[residue].nc   ;
    New_Structure.Residue[residue].size                       = Structure_One.Residue[residue].size ;

    if( ( New_Structure.Residue[residue].Atom = ( struct Atom * ) malloc ( ( Structure_One.Residue[residue].size + 1 ) * sizeof_Atom ) ) == NULL ) {
      GENERAL_MEMORY_PROBLEM
    }

    for( atom = 1 ; atom <= Structure_One.Residue[residue].size ; atom ++ ) {

      New_Structure.Residue[residue].Atom[atom] = Structure_One.Residue[residue].Atom[atom] ;

    }

  }

  for( residue = 1 ; residue <= Structure_Two.length ; residue ++ ) {

    new_residue = residue + Structure_One.length ;

    strcpy( New_Structure.Residue[new_residue].chainID            , Structure_Two.Residue[residue].chainID ) ;
    strcpy( New_Structure.Residue[new_residue].res_seq_plus_iCode , Structure_Two.Residue[residue].res_seq_plus_iCode ) ;
    strcpy( New_Structure.Residue[new_residue].olc                , Structure_Two.Residue[residue].olc ) ;
    New_Structure.Residue[new_residue].nc                         = Structure_Two.Residue[residue].nc   ;
    New_Structure.Residue[new_residue].size                       = Structure_Two.Residue[residue].size ;
    strcpy( New_Structure.Residue[new_residue].res_name           , Structure_Two.Residue[residue].res_name ) ;

    if( ( New_Structure.Residue[new_residue].Atom = ( struct Atom * ) malloc ( ( Structure_Two.Residue[residue].size + 1 ) * sizeof_Atom ) ) == NULL ) {
      GENERAL_MEMORY_PROBLEM
    }
    for( atom = 1 ; atom <= Structure_Two.Residue[residue].size ; atom ++ ) {

      New_Structure.Residue[new_residue].Atom[atom] = Structure_Two.Residue[residue].Atom[atom] ;

    }

  }

  return New_Structure ;

/************/

}



/************************/



float radius_of_structure( struct Structure This_Structure ) {

/************/

  /* Variables */
  float		present , largest ;


  /* Counters */
  int	residue , atom ;

/************/

  largest = 0 ;

  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {

    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {

      present = This_Structure.Residue[residue].Atom[atom].coord[1] * This_Structure.Residue[residue].Atom[atom].coord[1] + This_Structure.Residue[residue].Atom[atom].coord[2] * This_Structure.Residue[residue].Atom[atom].coord[2] + This_Structure.Residue[residue].Atom[atom].coord[3] * This_Structure.Residue[residue].Atom[atom].coord[3] ;

      if( present > largest ) largest = present ;

    }

  }

  return sqrt( largest ) ;

/************/

}


/******************************************************/
/* Derive 2nd structure - Kabsch & Sander             */
/******************************************************/

char* dersec(struct bridge *hbond,int length)
{

  
  int i,j, j1,j2;
  char b;
  char *secvect, *bridgetype;
  int *sheetbridge;
  int ib,ie,jb,je,conb, found;
  struct bridge *mybridge;

  

  mybridge=(struct bridge *) malloc(sizeof(struct bridge)*length+1);
  secvect=(char *) malloc(sizeof(char)*length+1);
  bridgetype=(char *) malloc(sizeof(char)*length+1);
  sheetbridge=(int *) malloc(sizeof(int)*length+1);
  for(i=0;i<length;i++)
    {
      secvect[i]='-';
      sheetbridge[i]=-1;
      bridgetype[i]='n';
      mybridge[i].aprt=-1;
      mybridge[i].bprt=-1;
      
    }
  secvect[length]='\0';
  /*  printf("%s, %d %d\n",secvect,strlen(secvect),length); */
  
 

  /* test for 4, 3, 5 helix */
  for (i=0;i<length-1;i++)
    {
      if ((((i-hbond[i].aprt)==4)&&((i+1)-hbond[i+1].aprt==4))&&(hbond[i].aprt>-1)&&(hbond[i+1].aprt>-1))
	{
	  if (hbond[i].aprt<i)
	    for (j=hbond[i].aprt;j<=(i);j++)
	      secvect[j]='H';
	  else
	    for (j=i;j<=hbond[i].aprt;j++)
	      secvect[j]='H';
	}
      if ((((i-hbond[i].aprt)==3)&&((i+1)-hbond[i+1].aprt==3))&&(hbond[i].aprt>-1)&&(hbond[i+1].aprt>-1))
	{
	  for (j=hbond[i].aprt;j<i;j++)
	    secvect[j]='h';
	}
      if ((((i-hbond[i].aprt)==5)&&((i+1)-hbond[i+1].aprt==5))&&(hbond[i].aprt>-1)&&(hbond[i+1].aprt>-1))
	{
	  for (j=hbond[i].aprt;j<i;j++)
	    secvect[j]='h';
	}
    }
  
  /* test for sheets bridges */

  for (i=1;i<length;i++)
    {
      j1 = 0;
      j2 = 0;
      j = i + 3;
      while (j2 == 0 && j < length) 
	{
	  if ((((hbond[i+1].aprt==j)||(hbond[i+1].bprt==j))&&
	       ((hbond[j].aprt==(i-1))||(hbond[j].bprt==(i-1))))||
	      (((hbond[j+1].aprt==i)||(hbond[j+1].bprt==i))&&
	       ((hbond[i].aprt==(j-1))||(hbond[i].bprt==(j-1)))))
	    {
	      /* b = parallel; */
	      b='p';	      
	    }	 
	  else if ((((hbond[i+1].aprt==(j-1))||(hbond[i+1].bprt==(j-1)))&&
		    ((hbond[j+1].aprt==(i-1))||(hbond[j+1].bprt==(i+1))))||
		   (((hbond[j].aprt==(i))||(hbond[j].bprt==(i)))&&
		    ((hbond[i].aprt==(j))||(hbond[i].bprt==(j)))))
	    {
	      /* b = antiparallel; */
	      b='a'; 
	    }	 	 
	  else
	    b = 'n';

	  if (b != 'n') 
	    {
	      if (j1 == 0) 
		{
		  j1 = j;
		  secvect[j]='E';
		  secvect[i]='E';
		} 
	      else if (j != j1) 
		{
		  j2 = j;
		  secvect[j]='E';
		  secvect[i]='E';
		}
	    }
	
	  j++;
	}
    }

  free(mybridge);
  free(bridgetype);
  free(sheetbridge);

  return(secvect);
}


/******************************************************/
/* Calculating the Bond Energy                        */
/******************************************************/

int CalculateBondEnergy( struct Amino_Acid This_aa )
{
    double dho,dhc;
    double dnc,dno;
    struct Atom cptr;
    long dx,dy,dz;
    long dist;
    int result , atom;

    for( atom = 1 ; atom <= This_aa.size ; atom ++ ) {         
      if (strcmp(This_aa.Atom[atom].atom_name," C  ")==0)
	{
	  cptr.coord[1]=250.0*This_aa.Atom[atom].coord[1];
	  cptr.coord[2]=250.0*This_aa.Atom[atom].coord[2];
	  cptr.coord[3]=250.0*This_aa.Atom[atom].coord[3];
	  cptr.serial=This_aa.Atom[atom].serial;
	}
      if (strcmp(This_aa.Atom[atom].atom_name," O  ")==0)
	{
	  optr.coord[1]=250.0*This_aa.Atom[atom].coord[1];
	  optr.coord[2]=250.0*This_aa.Atom[atom].coord[2];
	  optr.coord[3]=250.0*This_aa.Atom[atom].coord[3];
	  optr.serial=This_aa.Atom[atom].serial;
	}
    }
    /*    if( !(cptr=FindGroupAtom(group,2)) )  return(0);
	  if (!(optr=FindGroupAtom(group,3)) )  return(0); */

    dx = hxorg - optr.coord[1];  
    dy = hyorg - optr.coord[2];  
    dz = hzorg - optr.coord[3];   
    dist = dx*dx+dy*dy+dz*dz;

    if( dist < MinHDist ) 
        return( -9900 );
    dho = sqrt((double)dist);

    dx = hxorg - cptr.coord[1];  
    dy = hyorg - cptr.coord[2];
    dz = hzorg - cptr.coord[3];
    dist = dx*dx+dy*dy+dz*dz;
    if( dist < MinHDist ) 
        return( -9900 );
    dhc = sqrt((double)dist);
    /*    printf("energy dist1 %d %d \n",dist,MinHDist); */

    dx = nxorg - cptr.coord[1]; 
    dy = nyorg - cptr.coord[2]; 
    dz = nzorg - cptr.coord[3];
    dist = dx*dx+dy*dy+dz*dz;
    if( dist < MinHDist ) 
        return( -9900 );
    dnc = sqrt((double)dist);

    /*    printf("energy dist2 %d %d \n",dist,MinHDist); */
    dx = nxorg - optr.coord[1];
    dy = nyorg - optr.coord[2]; 
    dz = nzorg - optr.coord[3];
    dist = dx*dx+dy*dy+dz*dz;
    if( dist < MinHDist ) 
        return( -9900 );
    dno = sqrt((double)dist);

    result = (int)(QConst/dho - QConst/dhc + QConst/dnc - QConst/dno);

    if( result<-9900 ) 
      {   
	return -9900;
      } 
    else if( result>-500 ) 
        return 0;
    return result;
}



/******************************************************/
/* Calculating Protein H Bonds                        */
/******************************************************/


struct hbonder CalcProteinHBonds( struct Structure chn1,char *suplim )
{
    int energy, offset;
    struct Structure chn2;
    struct Amino_Acid group1;
    struct Amino_Acid group2;
    struct Atom ca1;
    struct Atom ca2;
    struct Atom pc1;
    struct Atom po1;
    struct Atom n1;    
    struct hydrobondreport the_report;
    char *secvstr;
    struct bridge *hsecvstr;
    struct hbonder hbon;

    char myaa[4];
    int pres, nrch, recount, recount2, secount,nelement;
    int pos1,pos2,hit,v,secv1;
    int dx,dy,dz,inmed,myco;
    double dco;
    long dist;
    double sumH;
    int atom, residue, residue2;
    char chain[2];
    char CAs1[6],CAs2[6];
    int *chainnr,r;
    hit=0;
    pos1=0;
    sumH=0.00000;
    nrch=1;

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
    /* printf(" nr of chains %d \n",nrch);  */

    /* fix the report structure */

    the_report.nrch=nrch;
    
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
    for( r = 0 ; r < (nrch) ; r++ ) {      
      hit=0;
      pos1=0;
      sumH=0.0000;
      recount=0;
      recount2=0;

      /************************************************************/
      /* allocate sectsr vector and hsecstr vector (hbond vector) */
      /************************************************************/
      secvstr=(char *)malloc((sizeof(char)*(chainnr[r+1]-chainnr[r])+1));
      hsecvstr=(struct bridge *) malloc((sizeof(struct bridge)*(chainnr[r+1]-chainnr[r])+1));
      
      hbon.secstruct=(char *)malloc((sizeof(char)*(chainnr[r+1]-chainnr[r])+1)); 
      for (secv1=0;secv1<(chainnr[r+1]-chainnr[r]);secv1++)
	{
	  secvstr[secv1]='-';
	  hsecvstr[secv1].aprt=-1;
	  hsecvstr[secv1].bprt=-1;
	}
      secvstr[secv1]='\0';

      /* Fix so that DNA chains is skipped */
      /* testing first AA - if nucleotide then skip */

      if ((strcmp(chn1.Residue[chainnr[r]].res_name,"  C")==0)||
	  (strcmp(chn1.Residue[chainnr[r]].res_name,"  A")==0)||
	  (strcmp(chn1.Residue[chainnr[r]].res_name,"  G")==0)||
	  (strcmp(chn1.Residue[chainnr[r]].res_name,"  T")==0))
	{	  
	  continue;
	}
      recount=-1;

      for( residue = chainnr[r] ; ((residue < (chainnr[r+1]))) ; residue ++ ) { 
	pos1++;
	recount++;
	if (hit==1)
	{
	  dx = (int)(pc1.coord[1]-po1.coord[1] );
	  dy = (int)(pc1.coord[2]-po1.coord[2] );
	  dz = (int)(pc1.coord[3]-po1.coord[3] );
	}
      else
	{
	  for( atom = 1 ; atom <= chn1.Residue[residue].size ; atom ++ ) {     
	    if (strcmp(chn1.Residue[residue].Atom[atom].atom_name," C  ")==0)
	      {
		pc1.coord[1]=250.0*chn1.Residue[residue].Atom[atom].coord[1];
		pc1.coord[2]=250.0*chn1.Residue[residue].Atom[atom].coord[2];
		pc1.coord[3]=250.0*chn1.Residue[residue].Atom[atom].coord[3];
		pc1.serial=chn1.Residue[residue].Atom[atom].serial;
	      }
	    if (strcmp(chn1.Residue[residue].Atom[atom].atom_name," O  ")==0)
	      {
		po1.coord[1]=250.0*chn1.Residue[residue].Atom[atom].coord[1];
		po1.coord[2]=250.0*chn1.Residue[residue].Atom[atom].coord[2];
		po1.coord[3]=250.0*chn1.Residue[residue].Atom[atom].coord[3];
		po1.serial=chn1.Residue[residue].Atom[atom].serial;
	      }
	  }
	  hit=1;
	  continue;
	}

      for( atom = 1 ; atom <= chn1.Residue[residue].size ; atom ++ ) {     
	if (strcmp(chn1.Residue[residue].Atom[atom].atom_name," C  ")==0)
	  {
	    pc1.coord[1]=250.0*chn1.Residue[residue].Atom[atom].coord[1];
	    pc1.coord[2]=250.0*chn1.Residue[residue].Atom[atom].coord[2];
	    pc1.coord[3]=250.0*chn1.Residue[residue].Atom[atom].coord[3];
	    pc1.serial=chn1.Residue[residue].Atom[atom].serial;
	  }
	if (strcmp(chn1.Residue[residue].Atom[atom].atom_name," O  ")==0)
	  {
	    po1.coord[1]=250.0*chn1.Residue[residue].Atom[atom].coord[1];
	    po1.coord[2]=250.0*chn1.Residue[residue].Atom[atom].coord[2];
	    po1.coord[3]=250.0*chn1.Residue[residue].Atom[atom].coord[3];
	    po1.serial=chn1.Residue[residue].Atom[atom].serial;
	  }
      }



      if (strstr(chn1.Residue[residue].res_name,"PRO")!=NULL)
	continue;

      for( atom = 1 ; atom <= chn1.Residue[residue].size ; atom ++ ) {     
	if (strcmp(chn1.Residue[residue].Atom[atom].atom_name," CA ")==0)
	  {
	    ca1.coord[1]=250.0*chn1.Residue[residue].Atom[atom].coord[1];
	    ca1.coord[2]=250.0*chn1.Residue[residue].Atom[atom].coord[2];
	    ca1.coord[3]=250.0*chn1.Residue[residue].Atom[atom].coord[3];
	    ca1.serial=chn1.Residue[residue].Atom[atom].serial;
	    /* printf(" Ca %d N %d  \n",ca1.serial,ca1.serial); */
	    pres=residue;
	    strcpy(myaa,chn1.Residue[residue].res_name);
	    strcpy(CAs1,chn1.Residue[residue].res_seq_plus_iCode);
	  }
	if (strcmp(chn1.Residue[residue].Atom[atom].atom_name," N  ")==0)
	  {
	    n1.coord[1]=250.0*chn1.Residue[residue].Atom[atom].coord[1];
	    n1.coord[2]=250.0*chn1.Residue[residue].Atom[atom].coord[2];
	    n1.coord[3]=250.0*chn1.Residue[residue].Atom[atom].coord[3];
	    n1.serial=chn1.Residue[residue].Atom[atom].serial;
	  }
      }
      dist = (long)dx*dx + (long)dy*dy + (long)dz*dz;
      dco = sqrt( (double)dist )/250.0;



      nxorg = (int)(n1.coord[1]);   
      hxorg = nxorg + (int)(dx/dco);
      nyorg = (int)(n1.coord[2]);   
      hyorg = nyorg + (int)(dy/dco);
      nzorg = (int)(n1.coord[3]);   
      hzorg = nzorg + (int)(dz/dco);
      res1 = 0; 
      res2 = 0;

      /* if several chains exist and bonding between the is adressed */
      /* the duplicate structure should be from the next chain       */
      /* i.e. a loop over the different chains are necessary         */
	 

      chn2=duplicate_structure(chn1);

      pos2=0;
      recount2=-1;

      for( residue2 = chainnr[r] ; ((residue2 < (chainnr[r+1]))) ; residue2 ++ ) {
	/* printf(" %s %d \n", chn2.Residue[residue2].res_name,residue2); */
	recount2++;
	pos2++;


	if ((residue==residue2)||(residue==(residue2+1)))
	  continue;
	
	for( atom = 1 ; atom <= chn2.Residue[residue2].size ; atom ++ ) {     
	  if (strcmp(chn2.Residue[residue2].Atom[atom].atom_name," CA ")==0)
	    {
	      ca2.coord[1]=250.0*chn2.Residue[residue2].Atom[atom].coord[1];
	      ca2.coord[2]=250.0*chn2.Residue[residue2].Atom[atom].coord[2];
	      ca2.coord[3]=250.0*chn2.Residue[residue2].Atom[atom].coord[3];
	      ca2.serial=chn2.Residue[residue2].Atom[atom].serial;
	      myco=residue2;
	    }
	}
	dx = (int)(ca1.coord[1]-ca2.coord[1]);
	if( (dist=(long)dx*dx) > MaxHDist )
	  continue;

	dy = (int)(ca1.coord[2]-ca2.coord[2]);
	if( (dist+=(long)dy*dy) > MaxHDist )
	  continue;

	dz = (int)(ca1.coord[3]-ca2.coord[3]);
	if( (dist+=(long)dz*dz) > MaxHDist )
	  continue;
	energy = CalculateBondEnergy(chn2.Residue[residue2]);

	
	if( energy )
	  { 
            offset = pos1 - pos2;
	    if( energy<res1 )
	      {   
		best2CA.serial = best1CA.serial;  
		best1CA.serial = ca2.serial;
		best2.serial = best1.serial;      
		best1.serial = optr.serial;
		res2 = res1;  
		res1 = energy;
		off2 = off1;        
		off1 = offset;
		strcpy(CAs2,chn2.Residue[residue2].res_seq_plus_iCode);
		secount=recount2;

	      } 
	    else if( energy<res2 )
	      {   
		best2CA.serial = ca2.serial;
		best2.serial = optr.serial;
		res2 = energy;
		off2 = offset;
		myco=recount2;
		strcpy(CAs2,chn2.Residue[residue2].res_seq_plus_iCode);
	      }
	  }
      } /* residue2 */

      if( res1<-0.5 ) 
        {   
	  if( res2<-0.5 ) 
	    {
	      sumH=sumH+(double)(float)res2/1000.0;
	      inmed=recount;
	      hsecvstr[inmed].bprt=myco;
	      /* printf("%d %3.1f \n",myco,(float)res2/1000.0); */
	      /* CreateHydrogenBond(ca1,best2CA,n1,best2,res2,off2); */
	    }
	  sumH=sumH+(double)(float)res1/1000.0;
	  inmed=recount;
	  hsecvstr[inmed].aprt=secount;
	}
      /* printf("%d %s %6.2lf %6.2lf\n",residue,chn1.Residue[residue].res_name,(double)(float)res1/1000.0,(double)(float)res2/1000.0); */
	 } /* for all residues in chain */


      
      strcpy(secvstr,dersec(hsecvstr,(chainnr[r+1]-chainnr[r])));
      strcpy(hbon.secstruct,secvstr);
      hbon.hbondv=sumH;
      nelement=0;
      for (secv1=0;secv1<(chainnr[r+1]-chainnr[r]);secv1++)
	{
	  if ((secvstr[secv1]=='H')||(secvstr[secv1]=='h')||(secvstr[secv1]=='E'))
	    nelement++;
	}

      /*   printf("SC %s ",chn1.ident);

      if (nrch==1)
	{
	  if (strlen(chn1.Residue[1].chainID)==0)
	    printf(" A ");
	  else 
	    printf("%s ",chn1.Residue[1].chainID);
	}
      else
	{
	  printf("%s ",chn1.Residue[chainnr[r]].chainID);
	}
      printf ("%d ",chainnr[r+1]-chainnr[r]); 
      printf ("%6.3lf ",sumH); 
      printf("%d ",nelement);
      printf("%s \n",secvstr); */
      
      free(secvstr);
      free(hsecvstr); 
      } 


    free(chainnr);
    /* return the_report; */
    return(hbon);
    free(secvstr);
    free(hbon.secstruct);
    free(hsecvstr);
}      


/******************************************************/
/* General function for calculating angles 3 atm      */
/******************************************************/


float calcangle(struct Atom at1,struct Atom at2,struct Atom at3) {

  double ulen2,vlen2;
  double ux,uy,uz;
  double vx,vy,vz;
  double temp;

  ux= at1.coord[1]-at2.coord[1];
  uy= at1.coord[2]-at2.coord[2];
  uz= at1.coord[3]-at2.coord[3];
  if( !ux && !uy && !uz )
        return 0.0;
  ulen2 = ux*ux + uy*uy + uz*uz;

  vx=at3.coord[1]-at2.coord[1];
  vy=at3.coord[2]-at2.coord[2];
  vz=at3.coord[3]-at2.coord[3];
  if( !vx && !vy && !vz )
    return 0.0;
  vlen2 = vx*vx + vy*vy + vz*vz;
  temp = (ux*vx + uy*vy + uz*vz)/sqrt(ulen2*vlen2);
  return Rad2Deg*acos(temp);
}



void justangle( struct Structure This_Structure ) {

  /* Variables */
  float		present , largest , phi;


  /* Counters */
  int	residue , atom , hit ;
  
  /* AToms */
  struct Atom atomics[7];

/************/



  for( residue = 1 ; residue <= (This_Structure.length) ; residue ++ ) {
    hit=0;
    printf(" %s ", This_Structure.Residue[residue].res_name,residue);
    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {     
      if (strcmp(This_Structure.Residue[residue].Atom[atom].atom_name," N  ")==0)
	{
	  atomics[0].coord[1]=This_Structure.Residue[residue].Atom[atom].coord[1];
	  atomics[0].coord[2]=This_Structure.Residue[residue].Atom[atom].coord[2];
	  atomics[0].coord[3]=This_Structure.Residue[residue].Atom[atom].coord[3];
	}
      if (strcmp(This_Structure.Residue[residue].Atom[atom].atom_name," CA ")==0)
	{
	  atomics[1].coord[1]=This_Structure.Residue[residue].Atom[atom].coord[1];
	  atomics[1].coord[2]=This_Structure.Residue[residue].Atom[atom].coord[2];
	  atomics[1].coord[3]=This_Structure.Residue[residue].Atom[atom].coord[3];
	}
      if (strcmp(This_Structure.Residue[residue].Atom[atom].atom_name," C  ")==0)
	{
	  atomics[2].coord[1]=This_Structure.Residue[residue].Atom[atom].coord[1];
	  atomics[2].coord[2]=This_Structure.Residue[residue].Atom[atom].coord[2];
	  atomics[2].coord[3]=This_Structure.Residue[residue].Atom[atom].coord[3];
	}
      if (strcmp(This_Structure.Residue[residue].Atom[atom].atom_name," O  ")==0)
	{
	  atomics[3].coord[1]=This_Structure.Residue[residue].Atom[atom].coord[1];
	  atomics[3].coord[2]=This_Structure.Residue[residue].Atom[atom].coord[2];
	  atomics[3].coord[3]=This_Structure.Residue[residue].Atom[atom].coord[3];
	}
      if (strcmp(This_Structure.Residue[residue].Atom[atom].atom_name," CB ")==0)
	{
	  atomics[4].coord[1]=This_Structure.Residue[residue].Atom[atom].coord[1];
	  atomics[4].coord[2]=This_Structure.Residue[residue].Atom[atom].coord[2];
	  atomics[4].coord[3]=This_Structure.Residue[residue].Atom[atom].coord[3];
	}
    }
      /* N-Ca-C */
      printf(" %f \n",calcangle(atomics[0],atomics[1],atomics[2]));

      /* N-Ca-Cb */

      /* Ca-C-O */





  }
}

/******************************************************/
/* General function for calculation of torsion angles */
/******************************************************/

float returntors(struct Atom at1,struct Atom at2,struct Atom at3,struct Atom at4) {

  double ax, ay, az;
  double bx, by, bz;
  double cx, cy, cz;
  double px, py, pz;
  double qx, qy, qz;
  double cosom, om;
  double rx, ry, rz;
  double plen,qlen;

  ax= at2.coord[1]-at1.coord[1];
  ay= at2.coord[2]-at1.coord[2];
  az= at2.coord[3]-at1.coord[3];

  /* if any null return 0 */

  /* for psi C - Ca */
  
  bx= at3.coord[1]-at2.coord[1];
  by= at3.coord[2]-at2.coord[2];
  bz= at3.coord[3]-at2.coord[3];

  /* if any null return 0 */
  
  /* for psi N+1 - C */

  cx= at4.coord[1]-at3.coord[1];
  cy= at4.coord[2]-at3.coord[2];
  cz= at4.coord[3]-at3.coord[3];

  /* if any null return 0 */
  
  /* invert if neccessary */
  /*   ay = -ay;  by = -by;  cy = -cy; 
       az = -az;  bz = -bz;  cz = -cz; */

  px = ay*bz - az*by;
  py = az*bx - ax*bz;
  pz = ax*by - ay*bx;

  qx = by*cz - bz*cy;
  qy = bz*cx - bx*cz;
  qz = bx*cy - by*cx;

  plen = px*px + py*py + pz*pz;
  qlen = qx*qx + qy*qy + qz*qz;

  cosom = (px*qx+py*qy+pz*qz)/sqrt(plen*qlen);

  if( cosom > 1.0 )
    { 
      return 0.0; 
    } 
  else if( cosom < -1.0 )
    {
      return 180.0; 
    }
      
  om = -Rad2Deg*acos(cosom);
     
  if ( om < -180. ) om += 360.;
  if ( om > 180. ) om -= 360.;
  
  rx = py*qz - pz*qy;
  ry = pz*qx - px*qz;
  rz = px*qy - py*qx;

  if (ax*rx+ay*ry+az*rz > 0.) 
    {
      return -om; 
    }
  else
    return om; 
}



void calctorsion( struct Structure This_Structure ) {

  /* Variables */
  float		present , largest , phi;


  /* Counters */
  int	residue , atom , hit ;
  
  /* AToms */
  struct Atom atomics[7];

/************/



  for( residue = 1 ; residue <= (This_Structure.length-1) ; residue ++ ) {
    hit=0;
    printf(" %s %d - ", This_Structure.Residue[residue].res_name,residue);
    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {     
      if (strcmp(This_Structure.Residue[residue].Atom[atom].atom_name," N  ")==0)
	{
	  atomics[1].coord[1]=This_Structure.Residue[residue].Atom[atom].coord[1];
	  atomics[1].coord[2]=This_Structure.Residue[residue].Atom[atom].coord[2];
	  atomics[1].coord[3]=This_Structure.Residue[residue].Atom[atom].coord[3];
	  hit=hit+1;
	}
      if (strcmp(This_Structure.Residue[residue].Atom[atom].atom_name," CA ")==0)
	{
	  atomics[2].coord[1]=This_Structure.Residue[residue].Atom[atom].coord[1];
	  atomics[2].coord[2]=This_Structure.Residue[residue].Atom[atom].coord[2];
	  atomics[2].coord[3]=This_Structure.Residue[residue].Atom[atom].coord[3];
	  hit=hit+1;
	}
      if (strcmp(This_Structure.Residue[residue].Atom[atom].atom_name," C  ")==0)
	{
	  atomics[3].coord[1]=This_Structure.Residue[residue].Atom[atom].coord[1];
	  atomics[3].coord[2]=This_Structure.Residue[residue].Atom[atom].coord[2];
	  atomics[3].coord[3]=This_Structure.Residue[residue].Atom[atom].coord[3];
	  hit=hit+1;
	}
    }
    for( atom = 1 ; atom <= This_Structure.Residue[residue+1].size ; atom ++ ) {     
      if (strcmp(This_Structure.Residue[residue+1].Atom[atom].atom_name," N  ")==0)
	{
	  atomics[4].coord[1]=This_Structure.Residue[residue+1].Atom[atom].coord[1];
	  atomics[4].coord[2]=This_Structure.Residue[residue+1].Atom[atom].coord[2];
	  atomics[4].coord[3]=This_Structure.Residue[residue+1].Atom[atom].coord[3];
	  hit=hit+1;
	}
      if (strcmp(This_Structure.Residue[residue+1].Atom[atom].atom_name," C  ")==0)
	{
	  atomics[5].coord[1]=This_Structure.Residue[residue+1].Atom[atom].coord[1];
	  atomics[5].coord[2]=This_Structure.Residue[residue+1].Atom[atom].coord[2];
	  atomics[5].coord[3]=This_Structure.Residue[residue+1].Atom[atom].coord[3];
	  hit=hit+1;
	}
      if (strcmp(This_Structure.Residue[residue+1].Atom[atom].atom_name," CA ")==0)
	{
	  atomics[0].coord[1]=This_Structure.Residue[residue+1].Atom[atom].coord[1];
	  atomics[0].coord[2]=This_Structure.Residue[residue+1].Atom[atom].coord[2];
	  atomics[0].coord[3]=This_Structure.Residue[residue+1].Atom[atom].coord[3];
	  hit=hit+1;
	}
    }
    /*    printf(" %lf %lf %lf \n", atomics[0].coord[1],atomics[0].coord[2],atomics[0].coord[3]); */
    if (residue==1)
      {
	printf(" psi %4.1f \n",returntors(atomics[1],atomics[2],atomics[3],atomics[4]));
	phi=returntors(atomics[3],atomics[4],atomics[0],atomics[5]);
      }
    else
      {
	printf(" phi %4.1f ",phi);
	printf(" psi %4.1f \n",returntors(atomics[1],atomics[2],atomics[3],atomics[4]));	
	phi=returntors(atomics[3],atomics[4],atomics[0],atomics[5]);
      }    
  }
    if (residue==This_Structure.length)
      {
	printf(" %s %d", This_Structure.Residue[residue].res_name,residue);
	printf(" phi %4.1f ",phi);
      }

}




/************************/



float total_span_of_structures( struct Structure Structure_1 , struct Structure Structure_2 ) {

  return  1 + ( ( radius_of_structure( Structure_1 ) + radius_of_structure( Structure_2 ) ) * 2 ) ;

}
