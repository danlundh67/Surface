#include "myutils.h"


void write_SurfPDB (gridpt grid[], char *outfile) {
  FILE *fp;
  int k, j, i;

  /*  printf("Writing SURFACE PDB to file: %s\n",outfile); */
  if( ( fp = fopen(outfile, "w" ) ) == NULL ) {
    printf( "The temp file %s. could not be opened.\nDying\n\n",outfile ) ;
    exit(  EXIT_FAILURE ) ;
  }

  fprintf(fp,"REMARK (c) Neil Voss, 2005\n");
  fprintf(fp,"REMARK PDB file created from %s\n",XYZRFILE);
  fprintf(fp,"REMARK Grid: %3.2f \tGRIDVOL: %4.2f \tWater_Res: %4.2 \tMaxProbe: %4.2f \tCutoff: %4.2f \n",GRID,GRIDVOL,WATER_RES,MAXPROBE,CUTOFF);

  //  out << "REMARK Date: " << std::ctime(&t) << flush;
  float count = 0;
  int anum=0, pnum=0;
  const float cat = DZ/60.0;
  float cut = cat;
  char line[128];

  /*  printf("Writing the grid to [ %s ]...\n",outfile); */
  /*  printBar(); */

  int sk=-1, sj=-1;
  for(k=0; k<DXYZ; k+=DXY) {
    sk++;
    count++;
    if(count > cut) {
      /* printf("^"); */
      cut += cat;
    }
    sj=-1;
    for(j=0; j<DXY; j+=DX) {
      sj++;
      for(i=0; i<DX; i++) {
        int pt = i+j+k;
        if(grid[pt]) {
          pnum++;
          if(isEdgePoint_Star(pt,grid)) {
            anum++;
            ijk2pdb(line,i,sj,sk,anum);
	    fprintf(fp,"%s \n",line);
          }
        }
  } } }

  fprintf(fp,"\n");
  fflush(fp);
  fclose(fp);
  /* printf("\ndone wrote %d of %d \n\n",anum,pnum); */
  return;
}

struct Structure my_write_SurfPDB (gridpt grid[]) {

  int k, j, i, numb;
  float count = 0;
  int anum=0, pnum=0;
  const float cat = DZ/60.0;
  float cut = cat;
  char line[128];
  float x,y,z;
  struct Structure This_Structure;


    int sk=-1, sj=-1;
  for(k=0; k<DXYZ; k+=DXY) {
    sk++;
    count++;
    if(count > cut) {
      cut += cat;
    }
    sj=-1;
    for(j=0; j<DXY; j+=DX) {
      sj++;
      for(i=0; i<DX; i++) {
        int pt = i+j+k;
        if(grid[pt]) {
          pnum++;
          if(isEdgePoint_Star(pt,grid)) {
            anum++;
          }
        }
      } 
    } 
  }

  printf("\ndone wrote %d of %d \n\n",anum,pnum);

  /* Memory allocation */
  if((This_Structure.Residue =(struct Amino_Acid *)malloc(((anum+2)*sizeof_Amino_Acid)))==NULL) 
    {
      GENERAL_MEMORY_PROBLEM
	;
    }
  for (numb=0;numb<anum;numb++)
    {
      if( ( This_Structure.Residue[numb].Atom =(struct Atom *)malloc(sizeof_Atom*2))== NULL) 
	{
          GENERAL_MEMORY_PROBLEM
	    ;
	}
      This_Structure.Residue[numb].size = 1;
    }
  

  This_Structure.length = anum+1;
  printf(" length %d \n",This_Structure.length);
  strcpy( This_Structure.ident , "Surface");

  count = 0;
  anum=0; 
  pnum=0;
  // cat = DZ/60.0;
  cut = cat;
  numb=1;
  sk=-1, sj=-1;
  for(k=0; k<DXYZ; k+=DXY) {
    sk++;
    count++;
    if(count > cut) {
      cut += cat;
    }
    sj=-1;
    for(j=0; j<DXY; j+=DX) {
      sj++;
      for(i=0; i<DX; i++) {
        int pt = i+j+k;
        if(grid[pt]) {
          pnum++;
          if(isEdgePoint_Star(pt,grid)) {
	    /* insert datapoints into structure data */
	    This_Structure.Residue[numb].Atom[1].serial = numb;
	    strcpy( This_Structure.Residue[numb].Atom[1].atom_name," O  \0");
	    strcpy( This_Structure.Residue[numb].res_name,"HOH") ;
	    strcpy( This_Structure.Residue[numb].chainID,"0") ;
	    strcpy( This_Structure.Residue[numb].res_seq_plus_iCode ,"   1\0");
	    x=(float)(i)*GRID + XMIN;
	    y=(float)(sj)*GRID + YMIN;
	    z=(float)(sk)*GRID + ZMIN;
	    This_Structure.Residue[numb].Atom[1].coord[1] = x;
	    This_Structure.Residue[numb].Atom[1].coord[2] = y;
	    This_Structure.Residue[numb].Atom[1].coord[3] = z;
	    This_Structure.Residue[numb].Atom[1].occupancy = 1.00;
	    This_Structure.Residue[numb].Atom[1].temp_factor = distFromPt(x,y,z);;
	    This_Structure.Residue[numb].Atom[1].surface=1;
	    strcpy(This_Structure.Residue[numb].olc,"0");
	    This_Structure.Residue[numb].nc=numb;
	    // printf("ATOM  %5d %4s %3s %1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f              %1s %2d\n", This_Structure.Residue[numb].Atom[1].serial, This_Structure.Residue[numb].Atom[1].atom_name, This_Structure.Residue[numb].res_name, This_Structure.Residue[numb].chainID, This_Structure.Residue[numb].res_seq_plus_iCode, This_Structure.Residue[numb].Atom[1].coord[1], This_Structure.Residue[numb].Atom[1].coord[2], This_Structure.Residue[numb].Atom[1].coord[3], This_Structure.Residue[numb].Atom[1].occupancy, This_Structure.Residue[numb].Atom[1].temp_factor, This_Structure.Residue[numb].olc, This_Structure.Residue[numb].nc ) ;
	    numb++;
            anum++;
          }
        }
      } 
    } 
  }

  return (This_Structure);
}
