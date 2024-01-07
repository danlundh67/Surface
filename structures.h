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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <malloc.h>
#include <stdlib.h>

/************/

/* These values directly below may be altered, and the programs rebuilt */

#define MAX_ROTATIONS  100000
#define NUMBER_TO_KEEP 10000
#define NUMBER_OF_CONSTRAINTS 50
#define SAVED_HEADER_LINES 1000
#define PI 3.14159265358979323846
#define Rad2Deg (180.0/PI)

/* I do not advise messing with anything below here */

/************/

/* Macros - realy, you don't want to touch these! */


#define GENERAL_MEMORY_PROBLEM printf( "You do not have enough memory ([m|re]alloc failure)\nDying\n\n" ) ; exit( EXIT_FAILURE ) ;

/************/
/* Structure reports and tranformation data */


#define ItemCount 5
#define CMAX 4
#define RMAX 4



typedef struct {
  int type;
  char *name;
  int intvalue;
  char *cvalue;
} ReportData;

static ReportData ReportMap[ItemCount] = {
  {0,"File",0,"None"},
  {1,"Chains",0,""},
  {1,"# Residues",0,""},
  {0,"H-Bond Energy",0,"None"},
  {1,"Hydrophobicity",0,"None"},
};

typedef int  MATRIX[RMAX][CMAX];
typedef double  fMATRIX[RMAX][CMAX];

typedef int  VECTOR[CMAX];
typedef double  fVECTOR[CMAX];

fMATRIX xrot,xrotm, yrot, yrotm,zrot,zrotm;

struct Pval{
  double x;
  double y;
  double z;
  double r;
};

struct Ival{
  int x;
  int y;
  int z;
  int r;
  int visible;
};

struct tripoint{
  int n1,n2,n3;
};


struct Pval *p;
struct Pval *pbak;
struct Ival *h;
struct tripoint *tri;

/* Datastructure to handle the surface dots for each atom */

struct dotpoints{
  struct Ival mydot[10];
};


struct dotpoints *dotatm;






/* The structures comprising a Structure (representation of an organic molecule in 3D) */

struct Atom{
	int		serial ;
	char		atom_name[5] ;
	float		coord[4] ;
	float		occupancy ;
	float		temp_factor ;
        float		charge; 
	int             solvent;
        int             surface;
        int             grida[3];
} ;

struct Amino_Acid{	
	char		res_name[4] ;
	char		chainID[2] ;
	char		res_seq_plus_iCode[6] ;
	char		olc[2] ;
        float           h_val;
	int		nc ;
	int		size ;
	struct Atom	*Atom ;
} ;

struct Structure{
	char			ident[256] ;
	int			length ;
	struct Amino_Acid	*Residue ;	
} ;

/************/


struct hydroreport{
  float solvent;
  float insolv;
  float charge;
  float incharge;
  float seqhydro;
} ;

struct hydrochain {
  double  HbondE;
  int    chain_length;
  char   chainID[2];
};


struct hydrobondreport{
  int                    nrch; 
  struct hydrochain      chainHB[20]; 
  char                   ident[256];
} ;


struct meashure{
  double volume;
  double surface;
  char tmppdb[80];
} ;

struct hbonder{
  char *secstruct;
  double hbondv;
};

struct subreport{
  int                      chain_length;
  char                     chainID[2];
  double                   HbondE;
  struct meashure          volsur;
  struct hbonder           m;
};

struct the_new_report{
  int                    nrch;
  char                   ident[256];
  struct subreport       subchains;
};


struct chainsrep{
  int                      chain_length;
  char                     chainID[2];
};


struct the_chain_report{
  int                    nrch;
  char                   ident[256];
  struct chainsrep       *subchains;
};




/************/
/* fix to store h-bonds */

struct bridge{
    int aprt;
    int bprt;
};
  

/* Angles structure */

struct Angle{
	int	n ;
	int	*z_twist ;
	int	*theta ;
	int	*phi ;
} ;

/************/

/* Score structure */

struct Score{
	int	score ;
	int	coord[4] ;
        int	angle[4] ;
        float	rpscore ;
        int	extra ;
} ;


/************/

/* Matrix Structure */

struct Matrix{
	char	description[100] ;
	float	distance ;
	float	score[21][21] ;
} ;

/************/


struct Hydro{
	char	res_name[4] ;
	float	value ;
        int	natm ;
} ;

/************/

/* Distribution matrix */

struct JONDAT{
  char res_name[4];
  int occur[50];
  int sum;
};

/******************************************/
/*  Stuff relating to plotting structure  */
/******************************************/

int ncords;
double viewx,viewy,viewz,scale_old, plotscale;

float thegrid, theprobe;




/*******************************************/
/* Protein Donor RAtom Coordinates         */
/*******************************************/

#define QConst (-6972000.0)
#define MaxHDist ((long)2250*2250)
#define MinHDist ((long)125*125)

static int hxorg,hyorg,hzorg;
static int nxorg,nyorg,nzorg;
static struct Atom best1CA;
static struct Atom best2CA;
static struct Atom best1;
static struct Atom best2;
static struct Atom optr;
static int res1,res2;
static int off1,off2;
float MinX,MaxX,MinY,MaxY,MinZ,MaxZ;
long SideLen;
float IVoxRatio;


/*******************************************/

/*******************************************/
/* SOLVENT STUFF                           */
/*******************************************/

#define VOXORDER       21
#define VOXORDER2      (VOXORDER*VOXORDER)
#define VOXSIZE        (VOXORDER2*VOXORDER)
#define ProbeRadius  300

extern void *HashTable[VOXSIZE];

static int VoxelCount,InVoxCount;
static int VoxelsDone;



/* Memory allocation sizes */

#define sizeof_Atom		sizeof( struct Atom )
#define sizeof_Amino_Acid	sizeof( struct Amino_Acid )
#define sizeof_Structure	sizeof( struct Structure )

/************/

extern struct Structure read_pdb_to_structure( char *pdb_file_name ) ;
extern void write_structure_to_pdb( struct Structure This_Structure , char *pdb_file_name ) ;
extern struct Structure duplicate_structure( struct Structure This_Structure ) ;
extern struct Structure translate_structure( struct Structure This_Structure , float x_shift , float y_shift , float z_shift ) ;
extern struct Structure translate_structure_onto_origin( struct Structure This_Structure ) ;
extern struct Structure rotate_structure( struct Structure This_Structure , int z_twist , int theta , int phi ) ;
extern struct Structure merge_structures( struct Structure Structure_One , struct Structure Structure_Two ) ;
extern float radius_of_structure( struct Structure This_Structure ) ;
extern float total_span_of_structures( struct Structure Structure_1 , struct Structure Structure_2 ) ;

extern struct Angle generate_global_angles( int angle_step ) ;
extern struct Angle generate_range_of_angles( int angle_step , int angle_range , int z_twist , int theta , int phi ) ;



extern void calctorsion( struct Structure This_Structure ) ;
extern void justangle( struct Structure This_Structure ) ;


extern struct Structure create_xyzr_structure( struct Structure This_Structure );
extern struct the_new_report structure_to_chain( struct Structure This_Structure, int chainnumber);
extern struct meashure calc_volume(struct Structure my,char pdbfile[], float gr, float pr);
extern double mymain(struct Structure chn1);
extern struct Structure  structure_onechain( struct Structure, int);
extern struct the_chain_report structure_chains( struct Structure This_Structure );
struct hbonder CalcProteinHBonds( struct Structure chn1,char *suplim );
