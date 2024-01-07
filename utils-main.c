/*
** utils.cpp
*/
#include "myutils.h"

float XMIN, YMIN, ZMIN;
float XMAX, YMAX, ZMAX;
int DX, DY, DZ;
int DXY, DXYZ;
unsigned int NUMBINS;
float MAXPROBE=15;
float GRID=0.5;
float GRIDVOL;
float WATER_RES;
float CUTOFF=10000;
char XYZRFILE[256]; //XYZR_FILE[0]='\0';

/* float WATER_RES=14.1372/GRIDVOL;


/*********************************************
**********************************************
         INITIALIZE FUNCTIONS
**********************************************
*********************************************/

//init functions
void finalGridDims (float maxprobe) {
  /*  printf("GRID  %f \n",GRID); */
  GRIDVOL=GRID*GRID*GRID;
  /*  printf("GRID  %f GRIDVOL %f\n",GRID,GRIDVOL); */

  WATER_RES=14137.2/GRIDVOL;
  MAXPROBE = maxprobe;
  XYZRFILE[0]='\0';
  XMIN=1000;
  YMIN=1000;
  ZMIN=1000;
  XMAX=-1000;
  YMAX=-1000;
  ZMAX=-1000;
}

float getIdealGrid () {
  int  cont;
  unsigned int numbins;
  int diff;
  unsigned int dx,dy,dz,dxy,dxyz;
  float grid;
  float bng, bpg;
  int bnd, bpd;
  int counts;
  
  cont=1;
  counts=0;
  bnd=1;
  bpd=-1;
  bng=-1.0;
  bpg=-1.0;
  grid=GRID;

  while(cont) {
    counts++;
    if(grid < 0.0001) { grid += 0.01; }
    dx=(int)((XMAX-XMIN)/grid+1);
    dy=(int)((YMAX-YMIN)/grid+1);
    dz=(int)((ZMAX-ZMIN)/grid+1);
    dxy=(dy*dx);
    dxyz=(dz*dxy);
    numbins = dxyz + dxy + dx + 1;

    diff = MAXBINS - numbins;
    //cerr << "grid = " << grid << "\tdiff = " << diff << 
//	"\tnumbins = " << numbins << endl;
    if(diff < 0) {
      if(bng < 0 || diff > bnd) {
        bng = grid;
        bnd = diff;
      }
      grid += 0.0001;
    }
    if(diff > 0) {
      if(bpg < 0 || diff < bpd) {
        bpg = grid;
        bpd = diff;
      }
      grid -= 0.0001;
    }
//    cerr << "\tbpg = " << bpg << "\tbng = " << bng << endl;
//    cerr << "\tbpdiff = " << bpd << "\tbndiff = " << bnd << endl;
    if(fabs(bpg - bng) < 0.0002 || counts > 10000) {
      cont = 0;
    }
  }
  return bpg;
}


void assignLimits () {
  float idealGrid;
  DX=(int)((XMAX-XMIN)/GRID+1);
  DY=(int)((YMAX-YMIN)/GRID+1);
  DZ=(int)((ZMAX-ZMIN)/GRID+1);
  DXY=(DY*DX);
  //unsigned int DYZ=(DY*DZ);
  //unsigned int DXZ=(DZ*DX);
  DXYZ=(DZ*DY*DX);
  //NUMBINS = 3*DXYZ;
  //NUMBINS = DXYZ + DXY + DXZ + DYZ + DX + DY + DZ + 1;
  NUMBINS = DXYZ + DXY + DX + 1;

  /*  printf("Precent filled NUMBINS/2^31: %d %% \n",(int)(NUMBINS*1000.0/MAXBINS)/10.0); */
  idealGrid = getIdealGrid();

  /* printf("Ideal Grid: %f \n",idealGrid); */

  /*if(NUMBINS > MAXBINS) {
    cout << MAXPROBE << "  grid " << GRID << " is too large; use " << 
	idealGrid << " for " << XYZRFILE << endl;
    cerr << "###### grid is too large ######" << endl << endl;
    exit (1);
  }*/
  /* printf("\n"); */

}

void testLimits (gridpt grid[]) {

  unsigned int i;

  /* printf("int(1.2) is %d \n",(int)1.2);
     printf("int(-1.2) is %d \n",(int)-1.2); */

  /*  printf("XMIN: %f\n",XMIN);
  printf("YMIN: %f\n",YMIN);
  printf("ZMIN: %f\n",ZMIN);
  printf(" DX: %d\n",DX);
  printf(" DY: %d\n",DY);
  printf(" DZ: %d\n",DZ);
  printf("DXYZ: %d\n",DXYZ);
  printf("NUMBINS: %d \n",NUMBINS);

  printf("First filled spot: ");
  for(i=0; i<NUMBINS && !grid[i]; i++) { }
  printf(" %d Last filled spot: ",i);
  for(i=NUMBINS-1; i>=0 && !grid[i]; i--) { }
  printf("%d \n\n",i); */
}

/*********************************************
**********************************************
         GRID UTILITY FUNCTIONS
**********************************************
*********************************************/


int countGrid (gridpt grid[]) {
  int voxels;
  voxels=0;
  unsigned int pt;

  /* printf("Counting up Voxels in Grid for Volume...  "); */


  for(pt=0; pt<NUMBINS; pt++) {
    if(grid[pt]) {
      voxels++;
    }
  }
  /*  printf("done [ %d voxels ]\n\n",voxels); */

  return voxels;
}

void zeroGrid (gridpt grid[]) {
  unsigned int pt;

  if (grid==NULL) {
    /* printf("Allocating Grid...\n"); */
    grid = (gridpt*) malloc (NUMBINS*sizeof(int));
    if (grid==NULL) { printf("GRID IS NULL\n"); exit (1); }
  }

  /*  printf("Zero-ing All Voxels in the Grid... %d",NUMBINS); */

  for(pt=0; pt<NUMBINS; pt++) {
    grid[pt] = 0;
  }
  /*  printf("done \n\n"); */

  return;
}

int copyGridFromTo (gridpt oldgrid[], gridpt newgrid[]) {
  return copyGrid(oldgrid,newgrid);
}

int copyGrid (gridpt oldgrid[], gridpt newgrid[]) {
  //Zero Grid Not Required
  int voxels;
  unsigned int pt;  

  voxels=0;
  if (newgrid==NULL) {
    /* printf(" Allocating Grid...\n"); */
    newgrid = (gridpt*) malloc (NUMBINS*sizeof(int));
    if (newgrid==NULL) { printf("GRID IS NULL \n"); exit (1); }
  }

  /* printf("Duplicating Grid and Counting up Voxels...  "); */

  for(pt=0; pt<NUMBINS; pt++) {
    if(oldgrid[pt]) {
      voxels++;
      newgrid[pt] = 1;
    } else {
      newgrid[pt] = 0;
    }
  }

  /* printf("done \n\n"); */

  return voxels;
}

void inverseGrid (gridpt grid[]) {
  unsigned int pt;
  if (grid==NULL) {
    /*  printf("Allocating Grid...\n"); */
    grid = (gridpt*) malloc (NUMBINS*sizeof(int));
    if (grid==NULL) { printf("GRID IS NULL\n"); exit (1); }
  }

  /*  printf("Inversing All Voxels in the Grid... "); */

  for(pt=0; pt<NUMBINS; pt++) {
    if(grid[pt]) {
      grid[pt] = 0;
    } else {
      grid[pt] = 1;
    }
  }
  /*  printf("done \n\n"); */

  return;
}


/*********************************************
**********************************************
        FILE BASED FUNCTIONS
**********************************************
*********************************************/





/*********************************************
**********************************************
        GENERATE GRIDS / GRID CHANGERS
**********************************************
*********************************************/

//void expand (gridpt oldgrid[], gridpt newgrid[]);
//void contract (gridpt oldgrid[], gridpt newgrid[]);

void trun_ExcludeGrid (const float probe, gridpt ACCgrid[], gridpt EXCgrid[]) { 
  /* contract
     limit grid search
     XMIN=(minmax[3] - MAXVDW - PROBE - 2*GRID);
     XMAX=(minmax[3] + MAXVDW + PROBE + 2*GRID); */
  int i,j, k;
  const int imin = 1;
  const int jmin = DX;
  const int kmin = DXY;
  const int imax = DX;
  const int jmax = DXY;
  const int kmax = DXYZ;

  float count = 0;
  const float cat = ((kmax-kmin)/DXY)/60.0;
  float cut = cat;

  if (EXCgrid==NULL) {
    /*    printf("Allocating Grid...\n"); */
    EXCgrid = (gridpt*) malloc (NUMBINS*sizeof(int));
    if (EXCgrid==NULL) { printf("GRID IS NULL\n"); exit (1); }
  }
  copyGridFromTo(ACCgrid,EXCgrid);

  /*  printf("Truncating Excluded Grid from Accessible Grid by Probe %f ...\n",probe);
      printBar(); */

  for(k=kmin; k<kmax; k+=DXY) {
  count++;
  if(count > cut) {
    /*    printf("^"); */
    cut += cat;
  }
  for(j=jmin; j<jmax; j+=DX) {
    for(i=imin; i<imax; i++) {
      if(!ACCgrid[i+j+k]) {
        const int k2 = k/DXY;
        const int j2 = j/DX;
        if(isEdgePoint(i,j,k,ACCgrid)) {
          empty_ExcludeGrid(i,j2,k2,probe,EXCgrid);
        }
      }
    }
  }
  }
  /*  printf("\ndone \n\n"); */
  return;
}

void grow_ExcludeGrid (const float probe, gridpt ACCgrid[], gridpt EXCgrid[]) {

  /* expands limit grid search
     XMIN=(minmax[3] - MAXVDW - PROBE - 2*GRID);
     XMAX=(minmax[3] + MAXVDW + PROBE + 2*GRID); */

  const int imin = 1;
  const int jmin = DX;
  const int kmin = DXY;
  const int imax = DX;
  const int jmax = DXY;
  const int kmax = DXYZ;
  int k, j, i;

  float count = 0;
  const float cat = ((kmax-kmin)/DXY)/60.0;
  float cut = cat;


  if (EXCgrid==NULL) {
    /* printf("Allocating Grid...\n"); */
    EXCgrid = (gridpt*) malloc (NUMBINS*sizeof(int));
    if (EXCgrid==NULL) { printf("GRID IS NULL\n"); exit (1); }
  }
  copyGrid(ACCgrid,EXCgrid);

//MUST USE COPYGRID BEFORE USING GROW_EXC



/*  printf("\nGrowing Excluded Grid from Accessible Grid by Probe %f ....\n",probe);
    printBar(); */

  for(k=kmin; k<kmax; k+=DXY) {
  count++;
  if(count > cut) {
    /* printf("^"); */
    cut += cat;
  }
  for(j=jmin; j<jmax; j+=DX) {
  for(i=imin; i<imax; i++) {
      const int pt = i+j+k;
      if(ACCgrid[pt]) {
//MUST USE COPYGRID BEFORE USING GROW_EXC
        if(isEdgePoint_Star(pt,ACCgrid)) {
          const int k2 = k/DXY;
          const int j2 = j/DX;
//ONLY use of fill_ExcludeGrid()
          fill_ExcludeGrid(i,j2,k2,probe,EXCgrid);
        }
      }
  }}}
  /*  printf("\ndone\n\n"); */
  return;
};

int get_Connected (gridpt grid[], gridpt connect[], const float x, const float y, const float z) {

  const int ip = (int)((x-XMIN)/GRID+0.5);
  const int jp = (int)((y-YMIN)/GRID+0.5);
  const int kp = (int)((z-ZMIN)/GRID+0.5);
  const int gp = ijk2pt(ip,jp,kp);
  int i, j, k, n;

  if (connect==NULL) {
    /* printf("Allocating Grid...\n"); */
    connect = (gridpt*) malloc (NUMBINS*sizeof(int));
    if (connect==NULL) { printf("GRID IS NULL\n"); exit (1); }
    zeroGrid(connect);
  }
  const int max = NUMBINS;
  int steps=0;
  int connected=0;
  if(gp >= 0 && gp <= max && grid[gp]) {
    connect[gp] = 1;
    /* printf("GetConnected..."); */


//    #define MAXLIST 1048576 //2^20
//    #define MAXLIST 32768 //2^15
//    #define MAXLIST 8192  //2^13
    int LIST[MAXLIST];
    int last = 1;
    LIST[0] = gp;
    LIST[1] = 0;

    while(last != 0) {
    //cerr << "." << flush;
      int newlast = 0;
      int NEWLIST[MAXLIST];
      for(n=0; n<last; n++) {
        steps++;
        int p = LIST[n];
        for(i=-1; i<=1; i++) {
        for(j=-DX; j<=DX; j+=DX) {
        for(k=-DXY; k<=DXY; k+=DXY) {
          int pt = p + i + j + k;
          if(grid[pt] && !connect[pt]) {  //isClose2Tunnel???
            connect[pt] = 1;
            connected++;
            if(newlast < MAXLIST-10) {
              NEWLIST[newlast] = pt;
              newlast++;
            }
          }
        }}}
      }
      for(n=0; n<newlast; n++) {
        LIST[n] = NEWLIST[n];
      }
      last=newlast;
      LIST[last]=0;
    }
    //cerr << endl;
    if(steps > 1) {
      /* printf(" performed %d steps \n",steps); */
    } else {
      printf(" done\n");
    }
  } else if(gp > 0 && gp < max) { 
    /* printf("GetConnected: Point OUT OF RANGE \n"); */

  } else {
    /* printf("GetConnected: Point is NOT FILLED\n"); */

  }
  return connected;
};

int get_ConnectedRange (gridpt grid[], gridpt connect[], const float x, const float y, const float z) {
//Get selected point in grid

  const int ip = (int)((x-XMIN)/GRID+0.5);
  const int jp = (int)((y-YMIN)/GRID+0.5);
  const int kp = (int)((z-ZMIN)/GRID+0.5);
  int id, jd, kd, i,j,k, n;
  int gp = ijk2pt(ip,jp,kp);
  if (connect==NULL) {
    /* printf("Allocating Grid...\n"); */
    connect = (gridpt*) malloc (NUMBINS*sizeof(int));
    if (connect==NULL) { printf("GRID IS NULL\n"); exit (1); }
    zeroGrid(connect);
  }

//Oops selected point isn't open! Better get new one
  if(!grid[gp]) {
    const int delta = (int)(1.50/GRID); 
    int stop=0;
    int gd=gp;
    for(id=-delta; !stop && id<=delta; id++) {
    for(jd=-delta; !stop && jd<=delta; jd++) {
    for(kd=-delta; !stop && kd<=delta; kd++) {
      gd = ijk2pt(ip+id,jp+jd,kp+kd);
      if(grid[gd]) {
        stop=1;
        gp=gd;
      }
    }}}
  }

  const int max = NUMBINS;
  int steps=0;
  int connected=0;
  if(gp >= 0 && gp <= max && grid[gp]) {
    connect[gp] = 1;
    /* printf("GetConnected..."); */

    /* #define MAXLIST 1048576 //2^20
       #define MAXLIST 32768 //2^15
       #define MAXLIST 8192  //2^13 */
    int LIST[MAXLIST];
    int last = 1;
    LIST[0] = gp;
    LIST[1] = 0;

    while(last != 0) {
    //cerr << "." << flush;
      int newlast = 0;
      int NEWLIST[MAXLIST];
      for(n=0; n<last; n++) {
        steps++;
        int p = LIST[n];
        for(i=-1; i<=1; i++) {
        for(j=-DX; j<=DX; j+=DX) {
        for(k=-DXY; k<=DXY; k+=DXY) {
          int pt = p + i + j + k;
          if(grid[pt] && !connect[pt]) {  //isClose2Tunnel???
            connect[pt] = 1;
            connected++;
            if(newlast < MAXLIST-10) {
              NEWLIST[newlast] = pt;
              newlast++;
            }
          }
        }}}
      }
      for(n=0; n<newlast; n++) {
        LIST[n] = NEWLIST[n];
      }
      last=newlast;
      LIST[last]=0;
    }
    //cerr << endl;
    if(steps > 1) {
      /*  printf(" performed %d steps \n",steps); */
    } else {
      /* printf(" done\n"); */
    }
  } else if(gp > 0 && gp < max) { 
    /* printf("GetConnected: Point OUT OF RANGE\n"); */
  } else {
    /*  printf("GetConnected: Point is NOT FILLED\n"); */
  }
  return connected;
};


int get_Connected_Point (gridpt grid[], gridpt connect[], const int gp) {
  
  int i,j,k,n;
  if (connect==NULL) {
    /* printf("Allocating Grid...\n"); */
    connect = (gridpt*) malloc (NUMBINS*sizeof(int));
    if (connect==NULL) { printf("GRID IS NULL\n"); exit (1); }
    zeroGrid(connect);
  }
  const int max = NUMBINS;
  int steps=0;
  int connected=0;
  if(gp >= 0 && gp <= max && grid[gp]) {
    connect[gp] = 1;
    /* printf("GetConnected..."); */


    /* #define MAXLIST 1048576 //2^20
       #define MAXLIST 32768 //2^15
       #define MAXLIST 8192  //2^13 */

    int LIST[MAXLIST];
    int last = 1;
    LIST[0] = gp;
    LIST[1] = 0;


    while(last != 0) {
    //cerr << "." << flush;
      int newlast = 0;
      int NEWLIST[MAXLIST];
      for(n=0; n<last; n++) {
        steps++;
        int p = LIST[n];
        for(i=-1; i<=1; i++) {
        for(j=-DX; j<=DX; j+=DX) {
        for(k=-DXY; k<=DXY; k+=DXY) {
          int pt = p + i + j + k;
          if(grid[pt] && !connect[pt]) {  //isClose2Tunnel???
            connect[pt] = 1;
            connected++;
            if(newlast < MAXLIST-10) {
              NEWLIST[newlast] = pt;
              newlast++;
            }
          }
        }}}
      }
      for(n=0; n<newlast; n++) {
        LIST[n] = NEWLIST[n];
      }
      last=newlast;
      LIST[last]=0;
    }
    //cerr << endl;
    if(steps > 1) {
      /* printf(" performed %d steps\n",steps); */
    } else {
      /* printf(" done\n"); */
    }
  }
  return connected;
}

int subt_Grids (gridpt biggrid[], gridpt smgrid[]) {
  int voxels=0;
  int error=0;
  unsigned int pt;
  /* float count = 0;
     const float cat = DZ/60.0;
     float cut = cat; */

  /*  printf("Subtracting Grids (Modifies biggrid)...  "); */
  //printBar();

  for(pt=0; pt<NUMBINS; pt++) {
  /*count++;
  if(count > cut) {
   printf("^"); 
    cut += cat;
  }*/
      if(smgrid[pt]) {
        if(biggrid[pt]) {
          voxels++;
          biggrid[pt] = 0;
        } else {
          error++;
        }
      }
  }
/*  printf("done [  %d  vox changed ]\n\n",voxels); */
  //cerr << "done [ " << error << " errors : " <<
//	int(1000.0*error/(voxels+error))/10.0 << "% ]" << endl;
  return voxels;
}

int intersect_Grids (gridpt grid1[], gridpt grid2[]) {
  /* GRID1 will CHANGE */
  int voxels=0;
  int changed=0;
  unsigned int pt;

  /* float count = 0;
     const float cat = NUMBINS/60.0;
     float cut = cat; */

  /* printf("Intersecting Grids...  "); */
  //printBar();

  for(pt=0; pt<NUMBINS; pt++) {
    /*count++;
    if(count > cut) {
    printf("^"); 
      cut += cat;
    }*/
    if(grid1[pt]) {
      if(!grid2[pt]) {
        changed++;
        grid1[pt] = 0;
      } else {
        voxels++;
      }
    }
  } 

/*  printf("done [ %d vox changed ] ",changed);
    printf("[ %d vox overlap :: ",voxels);
    printf(" %d % ]\n\n",(int)(1000.0*voxels/(voxels+changed))/10.0); */

  return voxels;
};

/*********************************************
**********************************************
        POINT BASED FUNCTIONS
**********************************************
*********************************************/

int fill_AccessGrid (const float x, const float y, const float z, const float R, gridpt grid[]) {


  // could redo to simplify

  const float cutoff = (R / GRID)*(R / GRID);

  int di,dj,dk;
  int filled;
  float distsq;
  int pt;

  const int imin = (int)((x - XMIN - R)/GRID - 1.0);
  const int jmin = (int)((y - YMIN - R)/GRID - 1.0);
  const int kmin = (int)((z - ZMIN - R)/GRID - 1.0);
  const int imax = (int)((x - XMIN + R)/GRID + 1.0);
  const int jmax = (int)((y - YMIN + R)/GRID + 1.0);
  const int kmax = (int)((z - ZMIN + R)/GRID + 1.0);

  const float xk = (x - XMIN)/GRID;
  const float yk = (y - YMIN)/GRID;
  const float zk = (z - ZMIN)/GRID;

  /*  printf(" eq %f %d %f \n",(x - XMIN - R)/GRID - 1.0, imin,cutoff); 
      printf(" xk  %f \n",(x - XMIN)/GRID); */

  filled=0;
  for(di=imin; di<=imax; di++) {
    for(dj=jmin; dj<=jmax; dj++) {
      for(dk=kmin; dk<=kmax; dk++) {
	distsq = (xk-di)*(xk-di) + (yk-dj)*(yk-dj) + (zk-dk)*(zk-dk);
	if(distsq < cutoff) {
	  pt = ijk2pt(di,dj,dk);
	  if(!grid[pt]) {
	    grid[pt] = 1;
	    filled++;
	  }
	}
      }
    }
  }
  /* printf(" filled  %f %f %f %f %d \n",x,y,z,R,filled); */
  return filled;
};

void empty_ExcludeGrid (const int i, const int j, const int k, const float probe, gridpt grid[]) {

  /* provides indexes (i,j,k) of grid where ijk2pt(i,j,k) = gridpt */
  const float R = probe/GRID; //Aug 19: correction for oversize
  const int r = (int)(R+1);
  const float cutoff = R*R;
  int nri,nrj,nrk,pri,prj,prk;

  int di,dj,dk;

//overflow checks (let's not go off the grid)
  if(i < r) { nri = -i; } else { nri = -r;}
  if(j < r) { nrj = -j; } else { nrj = -r;}
  if(k < r) { nrk = -k; } else { nrk = -r;}
  if(i + r >= DX) { pri = DX-i-1; } else { pri = r;}
  if(j + r >= DY) { prj = DY-j-1; } else { prj = r;}
  if(k + r >= DZ) { prk = DZ-k-1; } else { prk = r;}
  float distsq;
  int ind;
  for(di=nri; di<=pri; di++) {
  for(dj=nrj; dj<=prj; dj++) {
  for(dk=nrk; dk<=prk; dk++) {
     ind = ijk2pt(i+di,j+dj,k+dk);
     if(grid[ind]) {
       distsq = di*di + dj*dj + dk*dk;
       if(distsq < cutoff) {
         grid[ind] = 0;
       }
     }
  }}}
  return;
};

void fill_ExcludeGrid (const int i, const int j, const int k, const float probe, gridpt grid[]) {


  /* provides indexes (i,j,k) of grid where ijk2pt(i,j,k) = gridpt */
  const float R = probe/GRID; //Aug 19: correction for oversize
  const int r = (int)(R+1);
  const float cutoff = R*R;
  int nri,nrj,nrk,pri,prj,prk;

  int di,dj,dk;

  /* overflow checks (let's not go off the grid)*/

  if(i < r) { nri = -i; } else { nri = -r;}
  if(j < r) { nrj = -j; } else { nrj = -r;}
  if(k < r) { nrk = -k; } else { nrk = -r;}
  if(i + r >= DX) { pri = DX-i-1; } else { pri = r;}
  if(j + r >= DY) { prj = DY-j-1; } else { prj = r;}
  if(k + r >= DZ) { prk = DZ-k-1; } else { prk = r;}
  float distsq;
  int ind;
  for(di=nri; di<=pri; di++) {
  for(dj=nrj; dj<=prj; dj++) {
  for(dk=nrk; dk<=prk; dk++) {
     ind = ijk2pt(i+di,j+dj,k+dk);
     if(!grid[ind]) {
       distsq = di*di + dj*dj + dk*dk;
       if(distsq < cutoff) {
         grid[ind] = 1;
       }
     }
  }}}
  return;
};

int ijk2pt(int i, int j, int k) {
  return (int)(i+j*DX+k*DXY);
};

/* change bool to int in return type */

int isEdgePoint (const int i, const int j, const int k, gridpt grid[]) {
  //look at neighbors
  short int count=0;
  int dk, dj, di;


  for(dk=k-DXY; dk<=k+DXY; dk+=2*DXY) {
    count++;
    if(grid[i+j+dk]) {
      return 1;
    }
  }
  for(dj=j-DX; dj<=j+DX; dj+=2*DX) {
    count++;
    if(grid[i+dj+k]) {
      return 1;
    }
  }
  for(di=i-1; di<=i+1; di+=2) {
    count++;
    if(grid[di+j+k]) {
      return 1;
    }
  }
  if(count != 6) { printf("EdgePoint count %d != 6 \n"); }
  return 0;
};

int isEdgePoint_Fill (const int pt, gridpt grid[]) {

  /* look at neighbors */


  short int count=0;
  int di, dj, dk;


  for(di=-1; di<=1; di++) {
  for(dj=-DX; dj<=DX; dj+=DX) {
  for(dk=-DXY; dk<=DXY; dk+=DXY) {
    count++;
    if(!grid[pt+di+dj+dk]) {
      return 1;
    }
  }}}
  //cerr << "!" << endl;
  if(count != 27) { printf("EdgePoint count %d != 27\n"); }
  return 0;
};

int isEdgePoint_Star (const int pt, gridpt grid[]) {
  //look at neighbors
  short int count=0;
  int di, dj, dk;

  for(di=-1; di<=1; di+=2) {
    count++;
    if(!grid[pt+di]) {
      return 1;
    }
  }
  for(dj=-DX; dj<=DX; dj+=2*DX) {
    count++;
    if(!grid[pt+dj]) {
      return 1;
    }
  }
  for(dk=-DXY; dk<=DXY; dk+=2*DXY) {
    count++;
    if(!grid[pt+dk]) {
      return 1;
    }
  }
  if(count != 6) { printf("EdgePoint count %d != 6\n"); }
  return 0;
};

//void expand_Point (const int pt, gridpt grid[]);
//void contract_Point (const int pt, gridpt grid[]);

void ijk2pdb (char line[], int i, int j, int k, int n) {


  //char line[128];
  //cerr << "[i = " << i << "] " << flush;
  //cerr << "[j = " << j << "] " << flush;
  //cerr << "[k = " << k << "] " << flush;
  //cerr << "n = " << n << endl;

  line[0] = '\0';

  //LEAD IN
  strcpy(line,"ATOM  ");

  //ATOM NUMBER
  char temp[128];
  temp[0] = '\0';
  sprintf(temp,"%d",n%99999+1);
  padLeft(temp,5);
  strcat(line,temp);

  //ATOM & RESIDUE TYPES
  strcat(line,"  O   HOH  ");

  //RESIDUE NUMBER
  sprintf(temp,"%d",(n/10)%9999+1);
  padLeft(temp,4);
  strcat(line,temp);

  //GAP
  strcat(line,"    ");

  //XYZ COORDINATES 4.3

  float x = (float)(i)*GRID + XMIN;
  sprintf(temp,"%.3f",x);
  padLeft(temp,8);
  strcat(line,temp);
  float y = (float)(j)*GRID + YMIN;
  sprintf(temp,"%.3f",y);
  padLeft(temp,8);
  strcat(line,temp);
  float z = (float)(k)*GRID + ZMIN;
  sprintf(temp,"%.3f",z);
  padLeft(temp,8);
  strcat(line,temp);

  //OCCUPANCY
  sprintf(temp,"  1.00");
  strcat(line,temp);

  //TEMPERATURE
  float dist = distFromPt(x,y,z);
  sprintf(temp,"%.2f",dist);
  padLeft(temp,6);
  strcat(line,temp);

  //PRINT OUT
  //cerr << line << endl;

  return;
}

void limitToTunnelArea(const float radius, gridpt grid[]) {

  int pt;
  /*  printf("Limiting to Cylinder Around Exit Tunnel...  "); */

  for(pt=0; pt<=DXYZ; pt++) {
    if(!isCloseToVector(radius,pt)) {
      grid[pt] = 0;
    }
  }
  /*  printf("done \n\n"); */
  return;
};

int isCloseToVector (const float radius, const int pt) {

//GET QUERY POINT
  const float x = (int)(pt % DX) * GRID + XMIN;
  const float y = (int)((pt % DXY)/ DX) * GRID + YMIN;
  const float z = (int)(pt / DXY) * GRID + ZMIN;

//GET DISTANCE
  float dist = distFromPt(x,y,z);

//RETURN
  if(dist < radius) {
    return 1;
  }
  return 0;
};

float distFromPt (const float x, const float y, const float z) {
//INIT POINT
  //const float xp = 53.652;
  //const float yp = 141.358;
  //const float zp = 66.460;
  const float xp = 58.920;
  const float yp = 140.063;
  const float zp = 80.060;

//VECTOR
  const float xv =  0.58092;
  const float yv = -0.60342;
  const float zv =  0.54627;

//DIFFERENCE VECTOR
  const float dx = x - xp;
  const float dy = y - yp;
  const float dz = z - zp;

  const float lensq = dx*dx + dy*dy + dz*dz;
  const float dot = dx*xv + dy*yv + dz*zv;

//GET CROSS PRODUCT
  const float cross = sqrt(lensq - dot*dot);

  return cross;
}

float crossSection (struct real p, struct vector v, const gridpt grid[])
{
  return crossSection2(grid);
};


float crossSection2 (const gridpt grid[])
{
  float k, i, j;
//INIT POINT
  struct real p;
  p.x =  77.0;
  p.y = 124.0;
  p.z =  99.0;

//TUNNEL VECTOR
  struct vector v;
  v.x = -0.58092;
  v.y =  0.60342;
  v.z = -0.54627;

//GENERATE 2D GRID
//  FIND 2 PERP VECTORS
  struct vector v1, v2;
  v1.x =  0.60342; // =v.y
  v1.y =  0.58092; // =-v.x
  v1.z =  0.00000; // =0
  v2.x = -0.31734; //
  v2.y =  0.32963;
  v2.z =  0.70159;

//  const float x = (int)(pt % DX) * GRID + XMIN;
//  const float y = (int)((pt % DXY)/ DX) * GRID + YMIN;
//  const float z = (int)(pt / DXY) * GRID + ZMIN;
//  return (int)(i+j*DX+k*DXY);
//  const int ip = (int)((x-XMIN)/GRID+0.5);
//  const int jp = (int)((y-YMIN)/GRID+0.5);
//  const int kp = (int)((z-ZMIN)/GRID+0.5);

  struct real r;
  struct ind rp;
  int pt; 
  float count;
  //double mult = GRID*GRID*0.5*0.5*2.0/3.0;
  double mult = GRID*GRID/6.0;
  /*  printf("stepping"); */

  for(k=-5; k<100; k+=0.5) {
    k = (int)(k*4.0)/4.0;
    /*    printf(".");  */
   count = 0.0;
   float total = 0.0;
   for(i=-200; i<=200; i+=GRID*0.5) {
    for(j=-200; j<=200; j+=GRID*0.5) {
     r.x = p.x + v1.x*i + v2.x*j + v.x*k;
     r.y = p.y + v1.y*i + v2.y*j + v.y*k;
     r.z = p.z + v1.z*i + v2.z*j + v.z*k;
     if(r.x >= XMIN && r.x <= XMAX &&
	r.y >= YMIN && r.y <= YMAX &&
	r.z >= ZMIN && r.z <= ZMAX) {
       rp.i = (int)((r.x-XMIN)/GRID+0.5);
       rp.j = (int)((r.y-YMIN)/GRID+0.5);
       rp.k = (int)((r.z-ZMIN)/GRID+0.5);
       pt  = rp.i + rp.j*DX + rp.k*DXY;
       if(pt >= 0 && pt < DXYZ) {
         total++;
         if(grid[pt]) { count++; }
       }
     }
    }
   }

   //cerr << "crossSection:  " << k << " " << count << " of " << total << endl;

   /*   printf(" %f \t %lf \n",k,count*mult); */

  }
  printf("\n");

  return count;
}


/*********************************************
**********************************************
              STRING FUNCTIONS
**********************************************
*********************************************/

void padLeft(char a[], int n) {
  //len = 1 , n = 5
  int i;

  int len = strlen(a);
  if(len < n) {
    for(i=len; i<n; i++) {
      a[i] = ' ';
    }
    for(i=1; i<=len; i++) {
      a[n-i] = a[len-i]; //a[5] = a[1], a[4] = a[0]
      a[len-i] = ' ';
    }
    a[n] = '\0';
  }
}

void padRight(char a[], int n) {
  int len = strlen(a);
  if(len < n) {
    while (len < n) {
      a[len] = ' ';
      len++;
    }
    a[len] = '\0';
  }
}

void printBar () {
  /*  printf("|----+----+----+----+----+---<>---+----+----+----+----+----|\n"); */
  return;
}

void printVol (int vox) {
  //long double vol = vox*GRIDVOL;
  float tenp;
  tenp = 1000000.0;
  if((float)(vox)*GRIDVOL > tenp) {
    int cut = (int)(((float)(vox)/tenp)*GRIDVOL);
    /* printf("%d ,",cut); */
    vox = vox - (int)(cut*tenp/GRIDVOL);
  }
  tenp = 1000.0;
  if((float)(vox)*GRIDVOL > tenp) {
    int cut = (int)(((float)(vox)/tenp)*GRIDVOL);
    /* if(cut >= 100) {
      printf("%d ,", cut);
    } else if(cut >= 10) {
      printf("0%d ,",cut);
    } else if(cut >= 1) {
      printf("00%d",cut);
    } else {
      printf("000");
      } */
    vox = vox - (int)(cut*tenp/GRIDVOL);
  }
  double cut = (float)(vox)*GRIDVOL;
  /* if(cut >= 100) {
     printf("%lf ",cut);
     } else if(cut >= 10) {
     printf("0%lf",cut);
     } else if(cut >= 1) {
     printf("00%lf",cut);
     } else {
     printf("000");
     } */
  return;
};

void printVolCout (int vox) {
  long double vol = (float)(vox)*GRIDVOL;
  long double tenp;
  tenp = 1000000.0; //Millions
  if((float)(vox)*GRIDVOL > tenp) {
    int cut = (int)(((float)(vox)/tenp)*GRIDVOL);
    /* printf("%d",cut); */
    vox = vox - (int)(cut*tenp/GRIDVOL);
  }
  tenp = 1000.0; //Thousands
  if((float)(vox)*GRIDVOL > tenp) {
    int cut = (int)(((float)(vox)/tenp)*GRIDVOL);
    /*
    if(cut >= 100 || vol < 100000) {
       printf("%d",cut);
       } else if(cut >= 10) {
       printf("0%d",cut);
       } else if(cut >= 1) {
       printf("00%d",cut);
       } else {
       printf("000%d",cut);
       }
    */
    vox = vox - (int)(cut*tenp/GRIDVOL);
  }
  double cut = (float)(vox)*GRIDVOL;
  /* if(cut >= 100 || vol < 1000) {
    printf("%.0lf",cut);
  } else if(cut >= 10) {
    printf("0%.0lf",cut);
  } else if(cut >= 1) {
    printf("00%.0lf",cut);
  } else {
    printf("000%.0lf",cut);
  } 
  printf(" "); */ 
  return;
};

void basename (char str[], char base[]) {
  int loc=0;
  int i;
  for(i=0; str[i] != '\0'; i++) {
    if (str[i] == '/') {
      loc = i + 1;
    }
  }
  int max=0;
  for(i=loc; str[i] != '\0'; i++) {
    base[i-loc] = str[i];
    max = i + 1 - loc;
  }
  base[max] = '\0';
  return;
}

/*********************************************
**********************************************
              SURFACE AREA
**********************************************
*********************************************/

float surface_area (gridpt grid[]) {
  //Initialize Variables
  float surf=0.0;
  int i, j, k;
  const float wt[] = { 0.0, 0.894, 1.3409, 1.5879, 4.0, 2.6667, 
		      3.3333, 1.79, 2.68, 4.08, 0}; //weighting factors
/*
  wt[0]=0.0;   wt[1]=0.894; wt[2]=1.3409; wt[3]=1.5879;
  wt[4]=4.0;   wt[5]=2.6667; wt[6]=3.3333;
  wt[7]=1.79;  wt[8]=2.68;   wt[9]=4.08;
*/
  int type; // for return variables
  int edges[10]; //count types
  for(i=0; i<=9; i++) { edges[i] = 0; }
  float count = 0;
  const float cat = DZ/60.0;
  float cut = cat;

  //cerr << "DXY: " << DXY << "\tDX: " << DX << endl;
  /*  printf("Count Surface Voxels for Surface Area...\n");
      printBar(); */

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
          //cerr << pt << "::";
	  type = classifyEdgePoint(pt,grid);
          //surf += wt[type];
          edges[type]++;
        }
  } } }
  /*  printf("\n EDGES: "); */

  float totedge=0;
  for(i=1; i<=9; i++) {
    totedge += (float)(edges[i]);
  }
  for(i=1; i<=9; i++) {
    /* printf("s %d : %3.1f  ",i,((float)(1000*edges[i])/totedge)/10.0); */
    surf += edges[i]*wt[i];
  }
  /* printf("\n\n"); */
  return surf*GRID*GRID;
}

int classifyEdgePoint (const int pt, gridpt grid[]) {
  //look at neighbors
  short int count=0;
  int di,dj,dk;
  short int nb=0; //num of empty neighbors
  for(di=-1; di<=1; di+=2) {
    count++;
    if(!grid[pt+di]) {
      nb++;
    }
  }
  for(dj=-DX; dj<=DX; dj+=2*DX) {
    count++;
    if(!grid[pt+dj]) {
      nb++;
    }
  }
  for(dk=-DXY; dk<=DXY; dk+=2*DXY) {
    count++;
    if(!grid[pt+dk]) {
      nb++;
    }
  }
  //RETURN BASED ON NUMBER OF EMPTY NEIGHBORS
  if(count != 6) {
    printf("classifyEdgePoint count %d != 6\n",count);
  }
  if(pt < DXY) {
    printf("pt < DXY %d < %d \n",pt,DXY);
  }
  //if(pt + DXY > NUMBINS) { 
  //cerr << "pt > NUMBINS " << pt << " > " << NUMBINS << endl;
  //}
  if(nb == 0 || nb == 1) {
    return nb;
  } else if(nb == 2) {
//CHECK FOR BOTH CASES
    //check for cross gaps
    if(!grid[pt+1] && !grid[pt-1]) {
      return 7;
    }
    if(!grid[pt+DX] && !grid[pt-DX]) {
      return 7;
    }
    if(!grid[pt+DXY] && !grid[pt-DXY]) {
      return 7;
    }
    //the normal case
    return 2;
  } else if(nb == 3) {
//CHECK FOR BOTH CASES
    //check for cross gaps
    if(!grid[pt+1] && !grid[pt-1]) {
      return 4;
    }
    if(!grid[pt+DX] && !grid[pt-DX]) {
      return 4;
    }
    if(!grid[pt+DXY] && !grid[pt-DXY]) {
      return 4;
    }
    //the normal case
    return 3;
  } else if(nb == 4) {
//CHECK FOR BOTH CASES
    //check for cross fills
    if(grid[pt+1] && grid[pt-1]) {
      return 8;
    }
    if(grid[pt+DX] && grid[pt-DX]) {
      return 8;
    }
    if(grid[pt+DXY] && grid[pt-DXY]) {
      return 8;
    }
    //the normal case
    return 5;
  } else if(nb == 5) {
    return 6;
  } else if(nb == 6) {
    return 9;
  }
  printf("classifyEdgePoint neighbor count %d is wierd! \n",nb);
  return 0;
};

/*********************************************
**********************************************
              NEW FEATURES
**********************************************
*********************************************/

int fill_cavities(gridpt grid[]) {

  gridpt *cavACC=NULL;
  cavACC = (gridpt*) malloc (NUMBINS*sizeof(int));
  bounding_box(grid,cavACC);
  int pt;

//Create inverse access map
  subt_Grids(cavACC,grid); //modifies cavACC
  //int achanACC_voxels = countGrid(cavACC);

//Get first point
  int stop = 1; int firstpt = 0;
  for(pt=0; pt<NUMBINS && stop; pt++) {
    if(cavACC[pt]) { stop = 0; firstpt = pt;}
  }
  /*  printf("FIRST POINT: %d \n",firstpt); */
//LAST POINT
  stop = 1; int lastpt = 0;
  for(pt=NUMBINS-10; pt>0 && stop; pt--) {
    if(cavACC[pt]) { stop = 0; lastpt = pt;}
  }
  /*  printf("LAST  POINT: %d \n",lastpt); */

//Pull channels out of inverse access map
  gridpt *chanACC=NULL;
  chanACC = (gridpt*) malloc (NUMBINS*sizeof(int));
  zeroGrid(chanACC);
  get_Connected_Point(cavACC,chanACC,firstpt); //modifies chanACC
  get_Connected_Point(cavACC,chanACC,lastpt); //modifies chanACC
  //int chanACC_voxels = countGrid(chanACC);

//Subtract channels from access map leaving cavities
  subt_Grids(cavACC,chanACC); //modifies cavACC
  free (chanACC);
  int cavACC_voxels = countGrid(cavACC);


  int grid_before = countGrid(grid);
//Fill Cavities in grid[];
  for(pt=0; pt<NUMBINS; pt++) {
    if(cavACC[pt]) { grid[pt]=1; }
  }
  int grid_after = countGrid(grid);
  free (cavACC);

  printf("\n CAVITY VOLUME: ");
  printVol(cavACC_voxels);
  printf("\n BEFORE VOLUME: ");
  printVol(grid_before);
  printf("\n AFTER VOLUME:  ");
  printVol(grid_after);
  printf("\n DIFFERENCE:    ");
  printVol(grid_after-grid_before);
  printf("\n\n");

  return cavACC_voxels;
};


int bounding_box(gridpt grid[], gridpt bbox[]) {
  //find min x,y,z and max x,y,z

  int i,j,k;
  zeroGrid(bbox);

//PART I: Determine Extrema
  float count = 0;
  const float cat = DZ/60.0;
  float cut = cat;
  /*  printf("Determining Minima and Maxima...\n");
      printBar(); */
  int xmin=DX, ymin=DXY, zmin=DXYZ;
  int xmax=0, ymax=0, zmax=0;
  for(k=0; k<DXYZ; k+=DXY) {
    count++;
    if(count > cut) {
      /*      printf("^"); */
      cut += cat;
    }
    for(j=0; j<DXY; j+=DX) {
      for(i=0; i<DX; i++) {
        int pt = i+j+k;
        if(grid[pt]) {
          if(i < xmin) { xmin = i; }
          if(j < ymin) { ymin = j; }
          if(k < zmin) { zmin = k; }
          if(i > xmax) { xmax = i; }
          if(j > ymax) { ymax = j; }
          if(k > zmax) { zmax = k; }
        }
  } } }
  /*  printf("\nDONE \n \n"); */


//Grow by one
/*  
  xmin-=1;
  ymin-=DX;
  zmin-=DXY;
  xmax+=1;
  ymax+=DX;
  zmax+=DXY; 
*/

//PART II: FILL BOX
  int vol=0;
  count = 0;
  cut = cat;
  /*  printf("Fill Box... \n");
      printBar(); */
  for(k=zmin; k<=zmax; k+=DXY) {
    count++;
    if(count > cut) {
      /* printf("^"); */
      cut += cat;
    }
    for(j=ymin; j<=ymax; j+=DX) {
      for(i=xmin; i<=xmax; i++) {
        bbox[i+j+k] = 1;
        vol++;
  } } }
  printf("\nBOX VOXELS: ");
  printVol(vol);
  printf("\n\n");

  return vol;
};
