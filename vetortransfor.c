#include "myutils.h"

/* initialize a matrix and set to identity matrix
   (all 0's, 1's on diagonal) */
void identity()
{
   int i,j;
   for(i=0;i<CMAX;i++)
   for(j=0;j<RMAX;j++)
      if(i==j)
	{
         xrotm[j][i] = 1.0;
         xrot[j][i] = 1.0;
	 yrot[j][i] = 1.0;
	 yrotm[j][i] = 1.0;
	 zrot[j][i] = 1.0;
	 zrotm[j][i] = 1.0;
	}
      else
	{
         xrotm[j][i] = 0.0;
         xrot[j][i] = 0.0;
	 yrot[j][i] = 0.0;
	 yrotm[j][i] = 0.0;
	 zrot[j][i] = 0.0;
	 zrotm[j][i] = 0.0;	 	 
	}
}


void mat_mul(fVECTOR v1, fMATRIX mat2, fVECTOR v3)
{
     /* result stored in MATRIX new to avoid problems
        in case parameter mat3 == mat2 or mat 1 */
    fVECTOR new;
    float interm;
     int i,j;
     for(i=0;i<4;i++)
       {
	 interm=0.00000;
	 new[i]=0.00;
       for(j=0;j<4;j++)
	 {
	   interm=interm+v1[j]*mat2[j][i];
	 }
       new[i]=interm;
       }
     memcpy(v3,new,sizeof(new)); 
}


/* rotate about X axis  */
void xrot_init (double theta)
{
   double sintheta,costheta;
   /*   sintheta = sin(theta*0.017453293);
	costheta = cos(theta*0.017453293); */
   identity();   
   xrot[1][1] = cos(theta*0.017453293);
   xrot[1][2] = -sin(theta*0.017453293);
   xrot[2][1] = sin(theta*0.017453293);
   xrot[2][2] = cos(theta*0.017453293);
   xrotm[1][1] = cos(-1.0000*theta*0.017453293);
   xrotm[1][2] = -sin(-1.0000*theta*0.017453293);
   xrotm[2][1] = sin(-1.0000*theta*0.017453293);
   xrotm[2][2] = cos(-1.0000*theta*0.017453293);
}

void yrot_init (double theta)
{
  /* MATRIX rot;
     double sintheta,costheta;
     sintheta = sin(theta);
     costheta = cos(theta); */
   identity();   
   yrot[0][0] = cos(theta*0.017453293);
   yrot[0][2] = sin(theta*0.017453293);
   yrot[2][0] = -sin(theta*0.017453293);
   yrot[2][2] = cos(theta*0.017453293);

   yrotm[0][0] = cos(-1.0000*theta*0.017453293);
   yrotm[0][2] = sin(-1.0000*theta*0.017453293);
   yrotm[2][0] = -sin(-1.0000*theta*0.017453293);
   yrotm[2][2] = cos(-1.0000*theta*0.017453293);
}

void zrot_init (double theta)
{
   identity();   
   zrot[0][0] = cos(theta*0.017453293);
   zrot[0][1] = -sin(theta*0.017453293);
   zrot[1][0] = sin(theta*0.017453293);
   zrot[1][1] = cos(theta*0.017453293);

   zrotm[0][0] = cos(-1.00000*theta*0.017453293);
   zrotm[0][1] = -sin(-1.00000*theta*0.017453293);
   zrotm[1][0] = sin(-1.00000*theta*0.017453293);
   zrotm[1][1] = cos(-1.00000*theta*0.017453293);



}


void print_mat(fMATRIX m)
{
  int i,j;
  for(i=0;i<4;i++)
    {
    for(j=0;j<4;j++)
      printf("%lf ",m[i][j]);
    printf(" \n");
    }
}



/* ROTATE THE COORD ACCORDING TO MATRIX */

void select_viewer(fMATRIX m)
{
  int x,y,z;
  int i,n,j;
  fVECTOR v,v2;

  i=0;
  while (i<ncords)
    {
      v[0]=p[i].x;
      v[1]=p[i].y;
      v[2]=p[i].z;      
      v[3]=1;
      mat_mul(v,m,v2);            
      p[i].x=v2[0];
      p[i].y=v2[1];
      p[i].z=v2[2];
      i++; 
    }
}

void translate_to_origin()
{
  double mx,my,mz;
  int i;

  i=0;
  mx=0.000;
  mz=0.000;
  my=0.000;

  while (i<ncords)
    {
      mx=mx+h[i].x;
      my=my+h[i].y;
      mz=mz+h[i].z;
      i++;
    }
  mz=mz/(ncords);
  mx=mx/(ncords);
  my=my/(ncords);
  i=0;
  while (i<ncords)
    {
      h[i].x=h[i].x-(int)(float)mx;
      h[i].y=h[i].y-(int)(float)my;
      h[i].z=h[i].z-(int)(float)mz;
      i++;
    }  
}


void shift_xy()
{
  int i;
  double m;
  i=0;
   while (i<ncords)
    {
      m=p[i].x;
      p[i].x=p[i].y;
      p[i].y=m;
      i++;
    }
}

void zoom_size(double scale)
{
  int i;
  i=0;
  double rescale;
  rescale=scale/scale_old;

  while (i<ncords)
    {
      /* Todo */
      /* change scale so thta p[i] is always used (use incr scale) */

      p[i].x=p[i].x*scale;
      p[i].y=p[i].y*scale;
      p[i].z=p[i].z*scale;
      i++;
    }
  //  scale_old=scale;
} 



void resize_p(int wi,int hi)
{
  int i;
  int hmax,hmix,vmax,vmix;
  int vsize;
  float scale;
  
  /* if vsize==0 then width is smallest, i.e. determining scalable  */
  /* if vsize==1 then height is smallest i.e. determing the scale */

  vsize=0;
  if (wi>hi)
    vsize=1;
  vmax=-100000;
  hmax=vmax;
  vmix=100000;
  hmix=vmix;
  i=0;
  while (i<ncords)
    {
      if (h[i].x>hmax)
	hmax=h[i].x;      
      if (h[i].x<hmix)
	hmix=h[i].x;
      if (h[i].y>vmax)
	vmax=h[i].y;
      if (h[i].y<vmix)
	vmix=h[i].y;
      i++;
    }
    
  if ((abs(hmax)>=abs(hmix))&&(vsize==0))
    {
      scale=(float)(wi/2)/abs(hmax);
    }
  else if ((abs(hmax)<abs(hmix))&&(vsize==0))
    scale=(float)(wi/2)/hmix;
  else if ((abs(vmax)>=abs(vmix))&&(vsize==1))
    scale=(float)(hi/2)/vmax;
  else if ((abs(vmax)<abs(vmix))&&(vsize==1))
    scale=(float)(hi/2)/vmix;
  else
    scale=1.000;
  i=0;
  printf("%d %d %d %d %d %d %f \n",wi,hi,vmax,vmix,hmax,hmix,scale);
  while (i<ncords)
    {
      h[i].x=(int)(h[i].x*scale);
      h[i].y=(int)(h[i].y*scale);
      i++;
    }
}

void hide2_the_hidden(double scale)
{
  int i,j,k,x,y,z,r;
  int zmax,zmin,ymax,ymin,xmax,xmin;
  float fscale;

  /* simply sort the cordinates on z so that smallet z is ploted first */


  /* fscale=(float)scale; */
  free(h);
  h=(struct Ival *)malloc(sizeof(struct Ival)*ncords+2);
 

  /* convert to int vector since we cannot plot e.g. 0.5 pixels */

  i=0;
  xmax=xmin=(int)(float)p[0].x*scale;
  ymax=ymin=(int)(float)p[0].y*scale;
  zmax=zmin=(int)(float)p[0].z*scale;
  while (i<ncords)
    {      
      h[i].x=((int)((float)(p[i].x*scale)));
      h[i].y=((int)((float)(p[i].y*scale)));
      h[i].z=(int)(float)(p[i].z*scale);  
      h[i].r=(int)(float)(p[i].r*scale);  
      h[i].visible=0;
      i++;
    }
  translate_to_origin();
 i=0;
  
  while (i<(ncords))
    {

      j=i+1;
      while (j<ncords)
	{
	  if (h[i].z>h[j].z)
	    {
	      x=h[i].x;
	      y=h[i].y;
	      z=h[i].z;
	      r=h[i].r;
	      h[i].x=h[j].x;
	      h[i].y=h[j].y;
	      h[i].z=h[j].z;
	      h[i].r=h[j].r;
	      h[j].x=x;
	      h[j].y=y;
	      h[j].z=z;
	      h[j].r=r;
	    }
	  j++;	  
	}
      i++;
    }


  xmax=xmin=h[0].x;
  ymax=ymin=h[0].y;

  while (i<ncords)
    {      
      if (h[i].x>xmax)
	xmax=h[i].x;
      if (h[i].x<xmin)
	xmin=h[i].x;
      if (h[i].y>ymax)
	ymax=h[i].y;
      if (h[i].y<ymin)
	ymin=h[i].y;
    }
  
  for (i=xmin;i<(xmax-xmin);i++)
    {
      for (j=ymin;j<(ymax-ymin);j++)
	{
	  

	}

    }


}


void hide_the_hidden(double scale)
{
  int i,j,k,x,y,z,r;
  int zmax,zmin,ymax,ymin,xmax,xmin;
  float fscale;

  /* simply sort the cordinates on z so that smallet z is ploted first */


  /* fscale=(float)scale; */
  free(h);
  h=(struct Ival *)malloc(sizeof(struct Ival)*ncords+2);
 

  /* convert to int vector since we cannot plot e.g. 0.5 pixels */

  i=0;
  xmax=xmin=(int)(float)p[0].x*scale;
  ymax=ymin=(int)(float)p[0].y*scale;
  zmax=zmin=(int)(float)p[0].z*scale;
  while (i<ncords)
    {      
      h[i].x=((int)((float)(p[i].x*scale)));
      h[i].y=((int)((float)(p[i].y*scale)));
      h[i].z=(int)(float)(p[i].z*scale);  
      h[i].r=(int)(float)(p[i].r*scale);  
      h[i].visible=1;
      i++;
    }

  translate_to_origin();

  i=0;
  
  while (i<(ncords))
    {

      j=i+1;
      while (j<ncords)
	{
	  if (h[i].z>h[j].z)
	    {
	      x=h[i].x;
	      y=h[i].y;
	      z=h[i].z;
	      r=h[i].r;
	      h[i].x=h[j].x;
	      h[i].y=h[j].y;
	      h[i].z=h[j].z;
	      h[i].r=h[j].r;
	      h[j].x=x;
	      h[j].y=y;
	      h[j].z=z;
	      h[j].r=r;
	    }
	  j++;	  
	}
      i++;
      }

  /* The h-vector is now sorted such that the lowest z-value is first 
     in the array, i.e plotting them in this order results in overlaid
     obejects such that they actually hide each other. */
  /* FIX to GET TRIANGLES FOR PLOTTING */

  /* if (tri==NULL)
   {
     tri=(struct tripoint *)malloc(sizeof(struct tripoint)*ncords+2);

     i=0;
     zmin=0;
     while (i<(ncords))
       {
	 tri[i].n1=i;
	 j=i+1;
	 xmax=99999;
	 ymax=99999;
	 while (j<(ncords))
	   {

	     if (((h[i].x-h[j].x)*(h[i].x-h[j].x)+
		 (h[i].y-h[j].y)*(h[i].y-h[j].y)+
		  (h[i].z-h[j].z)*(h[i].z-h[j].z))<xmax)
	       {
		 xmax=((h[i].x-h[j].x)*(h[i].x-h[j].x)+
		       (h[i].y-h[j].y)*(h[i].y-h[j].y)+
		       (h[i].z-h[j].z)*(h[i].z-h[j].z));
		 xmin=j;
	       }
	     j++;
	   }  
	 k=i+1;
	 while (k<ncords)
	   {
	     if ((((h[i].x-h[k].x)*(h[i].x-h[k].x)+
		   (h[i].y-h[k].y)*(h[i].y-h[k].y)+
		   (h[i].z-h[k].z)*(h[i].z-h[k].z))<ymax)&&(xmin!=k))
	       {
		 ymax=((h[i].x-h[k].x)*(h[i].x-h[k].x)+
		       (h[i].y-h[k].y)*(h[i].y-h[k].y)+
		       (h[i].z-h[k].z)*(h[i].z-h[k].z));
		 ymin=k;
	       }
	     k++;
	   }
	 tri[i].n2=xmin;
	 tri[i].n3=ymin;
	 // printf(" %d %d %d \n",tri[i].n1,tri[i].n2,tri[i].n3);
	 i++;


       } 

       }*/ 

  free(dotatm);
  dotatm=NULL;
  if (dotatm==NULL)
    {
      dotatm=(struct dotpoints *)malloc(sizeof(struct dotpoints)*ncords+2);
      i=0;
      while (i<(ncords))
       {
	 dotatm[i].mydot[0].x=h[i].x;
	 dotatm[i].mydot[0].y=h[i].y;	   
	 dotatm[i].mydot[0].z=h[i].z+h[i].r;
	 dotatm[i].mydot[1].x=h[i].x+h[i].r;
	 dotatm[i].mydot[1].y=h[i].y;	   
	 dotatm[i].mydot[1].z=h[i].z;
	 dotatm[i].mydot[2].x=h[i].x-h[i].r;
	 dotatm[i].mydot[2].y=h[i].y;	   
	 dotatm[i].mydot[2].z=h[i].z;
	 dotatm[i].mydot[3].x=h[i].x;
	 dotatm[i].mydot[3].y=h[i].y+h[i].r;	   
	 dotatm[i].mydot[3].z=h[i].z;
	 dotatm[i].mydot[4].x=h[i].x;
	 dotatm[i].mydot[4].y=h[i].y-h[i].r;	   
	 dotatm[i].mydot[4].z=h[i].z;
	 dotatm[i].mydot[5].x=h[i].x+((int)((float)cos(45.0/180*PI)*((float)h[i].r)));
	 dotatm[i].mydot[5].y=h[i].y+((int)((float)sin(45.0/180*PI)*((float)h[i].r)));	   
	 dotatm[i].mydot[5].z=h[i].z;
	 dotatm[i].mydot[6].x=h[i].x+((int)((float)cos(45.0/180*PI)*((float)h[i].r)));
	 dotatm[i].mydot[6].y=h[i].y-((int)((float)sin(45.0/180*PI)*((float)h[i].r)));	   
	 dotatm[i].mydot[6].z=h[i].z;
	 dotatm[i].mydot[7].x=h[i].x-((int)((float)cos(45.0/180*PI)*((float)h[i].r)));
	 dotatm[i].mydot[7].y=h[i].y+((int)((float)sin(45.0/180*PI)*((float)h[i].r)));	   
	 dotatm[i].mydot[7].z=h[i].z;
	 dotatm[i].mydot[8].x=h[i].x-((int)((float)cos(45.0/180*PI)*((float)h[i].r)));
	 dotatm[i].mydot[8].y=h[i].y-((int)((float)sin(45.0/180*PI)*((float)h[i].r)));	   
	 dotatm[i].mydot[8].z=h[i].z;
	 dotatm[i].mydot[9].x=h[i].x+((int)((float)cos(45.0/180*PI)*((float)h[i].r)));
	 dotatm[i].mydot[9].y=h[i].y+((int)((float)sin(45.0/180*PI)*((float)h[i].r)));	   
	 dotatm[i].mydot[9].z=h[i].z+((int)((float)sin(45.0/180*PI)*((float)h[i].r)/2.000));
	 
	 i++;
       }
    }

}







char *tolower(char *in)
{
  int i;
  char c,*out;
  out=(char *)malloc(sizeof(in));
  for (i=0;i<strlen(in);i++)
    {
      c=in[i];
      if (c<91)
	{
	  c=c+32;
	}
      out[i]=c;
      
    }
  out[i]='\0'; 
  return out;

}


