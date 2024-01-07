#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "structures.h"


#define CITATION  fprintf(stderr, "\nCitation: Neil R Voss, et al. J Mol Biol. v360 (4): 2006, pp. 893-906 .\nDOI: http://dx.doi.org/10.1016/j.jmb.2006.05.023 \nE-mail: M Gerstein <mark.gerstein@yale.edu> or NR Voss <vossman77@yahoo.com>.\n\n")

#define COMPILE_INFO fprintf(stderr,"Program:")


#define MAXVDW  2.0

#define MAXBINS 2147483647 

#define MAXLIST 1048576 

struct ind { int i; int j; int k; float b;};
struct real { float x; float y; float z; };
struct sphere { float x; float y; float z; float r; };
struct vector { float x; float y; float z; };


extern float XMIN, YMIN, ZMIN;
extern float XMAX, YMAX, ZMAX;
extern int DX, DY, DZ;
extern int DXY, DXYZ;
extern unsigned int NUMBINS;
extern float MAXPROBE;
extern float GRID;
extern float GRIDVOL;
extern float WATER_RES;
extern float CUTOFF;
extern char XYZRFILE[256];


//typedef unsigned short int gridpt;
#ifndef _GRID_PT
  typedef int gridpt;
  #define _GRID_PT
#endif


//init functions
void finalGridDims (float maxprobe);
float getIdealGrid ();
void assignLimits ();
void testLimits (gridpt grid[]);

//grid util functions
int countGrid (gridpt grid[]);
void zeroGrid (gridpt grid[]);
int copyGridFromTo (gridpt oldgrid[], gridpt newgrid[]);
int copyGrid (gridpt oldgrid[], gridpt newgrid[]);
void inverseGrid (gridpt grid[]);

//file based functions
int read_NumAtoms (struct Structure This_S);
int fill_AccessGrid_fromFile (int numatoms, const float probe, struct Structure This_A, gridpt grid[]);
int get_ExcludeGrid_fromFile (int numatoms, const float probe, struct Structure This_S, gridpt EXCgrid[]);

//generate grids / grid changers
//void expand (gridpt oldgrid[], gridpt newgrid[]);
//void contract (gridpt oldgrid[], gridpt newgrid[]);
void trun_ExcludeGrid (const float probe, gridpt ACCgrid[], gridpt EXCgrid[]); //contract
void grow_ExcludeGrid (const float probe, gridpt ACCgrid[], gridpt EXCgrid[]); //expands
int get_Connected (gridpt grid[], gridpt connect[], const float x, const float y, const float z);
int get_ConnectedRange (gridpt grid[], gridpt connect[], const float x, const float y, const float z);
int get_Connected_Point (gridpt grid[], gridpt connect[], const int gp);
int subt_Grids (gridpt biggrid[], gridpt smgrid[]); //Modifies biggrid; returns difference
int intersect_Grids (gridpt grid1[], gridpt grid2[]); //Modifies grid1; returns final vox num

//point based function
int fill_AccessGrid (const float x, const float y, const float z, const float r, gridpt grid[]);
void empty_ExcludeGrid (const int i, const int j, const int k, const float probe, gridpt grid[]);
void fill_ExcludeGrid (const int i, const int j, const int k, const float probe, gridpt grid[]);
int isEdgePoint (const int i, const int j, const int k, gridpt grid[]);
int ijk2pt(int i, int j, int k);
int isEdgePoint_Fill (const int pt, gridpt grid[]);
int isEdgePoint_Star (const int pt, gridpt grid[]);
//void expand_Point (const int pt, gridpt grid[]);
//void contract_Point (const int pt, gridpt grid[]);
void ijk2pdb (char line[], int i, int j, int k, int n);

//special
int isCloseToVector (const float radius, const int pt);
void limitToTunnelArea(const float radius, gridpt grid[]);
float distFromPt (const float x, const float y, const float z);
float crossSection (struct real p, struct vector v, const gridpt grid[]);
float crossSection2 (const gridpt grid[]);
void real2index (struct real p, struct ind i);

//string functions
void padLeft(char a[], int n);
void padRight(char a[], int n);
void printBar ();
void printVol (int vox);
void printVolCout (int vox);
void basename(char str[], char base[]);

//surface area
float surface_area (gridpt grid[]);
int classifyEdgePoint (const int pt, gridpt grid[]);

//other ideas
//int convex_hull(gridpt grid[], gridpt hull[]);
//int convex_hull(int numatoms, char file[], gridpt hull[]);
int bounding_box(gridpt grid[], gridpt bbox[]);
//int bounding_box(int numatoms, char file[], gridpt bbox[]);
int fill_cavities(gridpt grid[]);

/*************************************************
//output functions (in utils-output.cpp)
**************************************************/
//void write_PDB (gridpt grid[], char outfile[]);
void write_SurfPDB (gridpt grid[], char outfile[]);
struct Structure my_write_SurfPDB (gridpt grid[]);


/* void write_EZD (gridpt grid[], char outfile[]);
void write_BlurEZD (gridpt grid[], char outfile[]);
void write_HalfEZD (gridpt grid[], char outfile[]);
void write_ThirdEZD (gridpt grid[], char outfile[]);
void write_FifthEZD (gridpt grid[], char outfile[]); */


