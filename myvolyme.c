#include "myutils.h" 

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

int fill_AccessGrid_fromFile (int numatoms, const float probe, struct Structure This_A, gridpt grid[]) {

  /* VARIABLES */
  float count;
  float cat;
  float cut;
  int filled=0;
  float x,y,z,r;
  int residue, atom;

  if (grid==NULL) {
    grid = (gridpt*) malloc (NUMBINS*sizeof(int));
    if (grid==NULL) { fprintf(stderr,"GRID IS NULL \n"); exit (1); }
  }
  zeroGrid(grid);

  if(!XYZRFILE[0]) { strcpy(XYZRFILE,This_A.ident); }

  count = 0.0;
  cat = (float)numatoms/60.0;
  cut = cat;

 
  for( residue = 1 ; residue <= This_A.length ; residue ++ ) {
    for( atom = 1 ; atom <= This_A.Residue[residue].size ; atom ++ ) {
      count=count+1.0;
      x = This_A.Residue[residue].Atom[atom].coord[1];
      y = This_A.Residue[residue].Atom[atom].coord[2];
      z = This_A.Residue[residue].Atom[atom].coord[3];
      r = This_A.Residue[residue].Atom[atom].charge;
      filled += fill_AccessGrid(x,y,z,r+probe,grid);
    }
  }

  return filled;
}


int get_ExcludeGrid_fromFile (int numatoms, const float probe,
			      struct Structure This_S, gridpt EXCgrid[]) {
//READ FILE INTO ACCGRID
  gridpt *ACCgrid;
  int voxels;

  ACCgrid = (gridpt*) malloc (NUMBINS*sizeof(int));
  if (ACCgrid==NULL) { fprintf(stderr,"GRID IS NULL \n"); exit (1); }
  fill_AccessGrid_fromFile(numatoms,probe,This_S,ACCgrid);

//TRUNCATE GRID
  trun_ExcludeGrid(probe,ACCgrid,EXCgrid);

//RELEASE ACCGRID
  free (ACCgrid);

  voxels = countGrid(EXCgrid);

  return voxels;
}


int read_NumAtoms (struct Structure This_S) {

  int residue;
  int atom;
  int i;
  float x,y,z,r;
  int count = 0;
  char line[256];
  float minmax[6];
  float FACT;

  minmax[0] = 100;  minmax[1] = 100;  minmax[2] = 100;
  minmax[3] = -100;  minmax[4] = -100;  minmax[5] = -100;
  strcpy( XYZRFILE, This_S.ident ) ;


  for( residue = 1 ; residue <= This_S.length ; residue ++ ) {
    for( atom = 1 ; atom <= This_S.Residue[residue].size ; atom ++ ) {
      count++;

      x = This_S.Residue[residue].Atom[atom].coord[1];
      y = This_S.Residue[residue].Atom[atom].coord[2];
      z = This_S.Residue[residue].Atom[atom].coord[3];
      if (count % 3000 ==0 ) { ; /* printf("."); */}
      if(x < minmax[0]) { minmax[0] = x; }
      if(x > minmax[3]) { minmax[3] = x; }
      if(y < minmax[1]) { minmax[1] = y; }
      if(y > minmax[4]) { minmax[4] = y; }
      if(z < minmax[2]) { minmax[2] = z; }
      if(z > minmax[5]) { minmax[5] = z; }      
    }
  }


  //INCREASE GRID SIZE TO ACCOMODATE SPHERES
  //ALSO ROUND GRID SIZE SO THE XMIN = INTEGER * GRID (FOR BETTER OUTPUT)
  FACT = MAXVDW + MAXPROBE + 2*GRID;
  for(i=0;i<=2;i++) { 
    minmax[i] -= FACT; 
    minmax[i] = (int)(minmax[i]/(4*GRID)-1)*4*GRID;
  }
  for(i=3;i<=5;i++) { 
    minmax[i] += FACT;
    minmax[i] = (int)(minmax[i]/(4*GRID)+1)*4*GRID;
  }
  if(minmax[0] < XMIN) { XMIN = minmax[0]; }
  if(minmax[1] < YMIN) { YMIN = minmax[1]; }
  if(minmax[2] < ZMIN) { ZMIN = minmax[2]; }
  if(minmax[3] > XMAX) { XMAX = minmax[3]; }
  if(minmax[4] > YMAX) { YMAX = minmax[4]; }
  if(minmax[5] > ZMAX) { ZMAX = minmax[5]; }

  return count;
}





struct meashure calc_volume(struct Structure my,char pdbfile[],float gr, float pr) {

  /*  COMPILE_INFO;
      CITATION; */

  // ****************************************************
  // INITIALIZATION
  // ****************************************************

  //HEADER INFO

  double PROBE=pr;
  GRID=gr;
  int numatoms;
  gridpt *EXCgrid;
  int voxels, myvoxels;
  double surf, plotsurf;
  struct meashure mysv;
  double vol;
  char *dumperpdb;
  struct Structure surfacepdb;

  //INITIALIZE GRID
  finalGridDims(PROBE);

  //FIRST PASS, MINMAX
  numatoms = read_NumAtoms(my);

  /*  printf("%s %d \n",my.ident,numatoms); */

  //CHECK LIMITS & SIZE
  assignLimits();

  // ****************************************************
  // STARTING FIRST FILE
  // ****************************************************
  //READ FILE INTO SASGRID



  EXCgrid = (gridpt*) malloc (NUMBINS*sizeof(int));
  if (EXCgrid==NULL) {fprintf(stderr,"GRID IS NULL \n"); exit (1); }
  zeroGrid(EXCgrid);
  if(PROBE > 0.0) { 
    voxels = get_ExcludeGrid_fromFile(numatoms,PROBE,my,EXCgrid);
  } else {
    voxels = fill_AccessGrid_fromFile(numatoms,0.0,my,EXCgrid);
  }


  surf = surface_area(EXCgrid);

  dumperpdb=(char *)malloc(sizeof(char)*80);

  dumperpdb[0]='\0';
  strcat(dumperpdb,my.ident);
  strcat(dumperpdb,"-temp-");
  if (my.Residue[1].chainID[0]!=' ')
    strcat(dumperpdb,my.Residue[1].chainID);
  

  if(pdbfile[0] != '\0') {
    write_SurfPDB(EXCgrid,dumperpdb);  
    
    }
    
  
 
  //RELEASE TEMPGRID
  free (EXCgrid);


  /*  printf("VC %s %s %d vol: ",my.ident,my.Residue[1].chainID,my.length);
  printVolCout(voxels);
  printf("surf: %.0lf  %lf \n",surf,vol); */

  vol= (double)(float)(voxels)*GRIDVOL;
  
  

  mysv.volume= vol;
  mysv.surface=surf;
  strcpy(mysv.tmppdb,dumperpdb);
  free(dumperpdb);
  return mysv;
}

