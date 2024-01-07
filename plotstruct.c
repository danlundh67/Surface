#include <windows.h>
#include <wingdi.h>
#include "structures.h"
#include "StrchBlt.h"


void read_to_list(struct Structure chain)
{
  float x,y,z,r;
  int i,n;
  int residue, atom;
  i=0;
  for( residue = 1 ; residue <= chain.length ; residue ++ ) {
    for( atom = 1 ; atom <= chain.Residue[residue].size ; atom ++ ) 
      {	
	if (chain.Residue[residue].Atom[atom].surface==1)
	  i++;	  
      }
  }

  n=i;
  ncords=n;
  p=(struct Pval *)malloc((sizeof(struct Pval)*i)+2);
  pbak=(struct Pval *)malloc((sizeof(struct Pval)*i)+2);
 
  //  printf(" Outer atoms %d \n",ncords);

  i=0;
  for( residue = 1 ; residue <= chain.length ; residue ++ ) {
    for( atom = 1 ; atom <= chain.Residue[residue].size ; atom ++ ) {
      if (chain.Residue[residue].Atom[atom].surface==1)
	{
	  x=chain.Residue[residue].Atom[atom].coord[1];
	  y=chain.Residue[residue].Atom[atom].coord[2];
	  z=chain.Residue[residue].Atom[atom].coord[3];
	  r=chain.Residue[residue].Atom[atom].charge;
	  p[i].x=x/1.00000;
	  p[i].y=y/1.00000;
	  p[i].z=z/1.00000;
	  p[i].r=r/1.00000;
	  pbak[i].x=p[i].x;
	  pbak[i].y=p[i].y;
	  pbak[i].z=p[i].z;
	  pbak[i].r=p[i].r;
	  i++;
	}       
    }
  }
}



void print_data(HWND hwnd, HDC hdc,BOOL bRadio1)
{
  float x,y,z,zscale2;
  int tx,ty,tz, r,k;
  struct Pval s;
  int hi, wi , sz;
  int i,n;
  char *st;
  RECT rcClient; // client area rectangle 
  HPEN hpen, hpenold;
  
  st=(char *)malloc(sizeof(char)*80);
  sprintf(st,"< %4.2lf, %4.2lf ,%4.2lf %d >",viewx,viewy,viewz,ncords);

  
  
  /* sorts the pval vector on z (smallest first) to enable dept drawing */
  /* also fixes so that plot vectors are stored in int vector since only */
  /* entire pixels can be drawn (func: hide_the_hidden)                 */
  
  hide_the_hidden(plotscale); 

  /* moves the coordinates of int vector points so that it centres the plot */

  // translate_to_origin(); 

  GetClientRect(GetDlgItem(hwnd,IDC_PLOT), &rcClient);

  wi=(int)rcClient.right;
  hi=(int)rcClient.bottom;
  /*  resize_p(wi,hi); */
  // sz=(int)widget->allocation.height/50;
  // printf("%d %d \n",wi,hi); 

  i=0;
  
  zscale2=(float)(h[ncords-1].z-h[0].z);

  /* FIX SUCH THAT MOLECULE IS DISPLAYED AT CENTRE OF DRAWING AREA */

  while (i<ncords)
    {

      if (h[i].visible==1)
	{	  
	  h[i].x=(int)(h[i].x+wi/2);
	  h[i].y=(int)(h[i].y+hi/2);
	  h[i].z=h[i].z;
	}
      i++;
    }
  


  i=0;
  if (bRadio1==TRUE)
    {
      while (i<ncords)
	{
	  if (h[i].visible==1)
	    {
	      x=h[i].x;
	      y=h[i].y;
	      z=h[i].z;   
	      tx=h[i].x;
	      ty=h[i].y;
	      r=h[i].r;
	  /* tz=(int)1*plotscale; */
	  /* radius of an oxygen is assumed to be 1.4 A, i.e. diameter appr 3 */
	  /* this radius increase with the zoom (plotscale) factor */

	  // draw a white ellipse, white borders
	      
 
	      MoveToEx(hdc,tx+r,ty,(LPPOINT)NULL);
	      BeginPath(hdc);
	      SelectObject(hdc, GetStockObject(WHITE_BRUSH));
	      SelectObject(hdc, GetStockObject(BLACK_PEN));
	      AngleArc(hdc,tx,ty,r,0,360); 
	      EndPath(hdc);
	      StrokeAndFillPath(hdc);	   
	    }
	  i++;
	}
    }
  i=0;
  /*  if (bRadio1==FALSE)
    {
      while (i<ncords)
	{
	  x=h[tri[i].n1].x;
	  y=h[tri[i].n1].y;
	  z=h[tri[i].n1].z;   
	  r=h[tri[i].n1].r;
	  MoveToEx(hdc,x,y,(LPPOINT)NULL);
	  SelectObject(hdc, GetStockObject(WHITE_BRUSH));
	  SelectObject(hdc, GetStockObject(BLACK_PEN));
	  LineTo(hdc,h[tri[i].n2].x,h[tri[i].n2].y);
	  LineTo(hdc,h[tri[i].n3].x,h[tri[i].n3].y);
	  LineTo(hdc,x,y);
	  i++;
	}

    }
  */

  i=0;
  if (bRadio1==FALSE)
    {
      hpen=CreatePen(PS_DOT,1,RGB(230,230,230));
     
      while (i<ncords)
	{
	  x=h[i].x;
	  y=h[i].y;
	  z=h[i].z;   
	  r=h[i].r;
	  MoveToEx(hdc,x+r,y,(LPPOINT)NULL);
	  BeginPath(hdc);
	  SelectObject(hdc, GetStockObject(WHITE_BRUSH));
	  hpenold=SelectObject(hdc,hpen);
	  AngleArc(hdc,x,y,r,0,360); 
	  EndPath(hdc);
	  StrokeAndFillPath(hdc); 
	  for (k=0;k<10;k++)
	    {
	      x=(int)dotatm[i].mydot[k].x+(int)wi/2;
	      y=(int)dotatm[i].mydot[k].y+(int)hi/2;
	      z=dotatm[i].mydot[k].z;;   
	      /*BeginPath(hdc);
	      MoveToEx(hdc,x+1,y,(LPPOINT)NULL);
	      SelectObject(hdc, GetStockObject(WHITE_BRUSH));
	      SelectObject(hdc, GetStockObject(BLACK_PEN));
	      AngleArc(hdc,x,y,1,0,360);
	      EndPath(hdc);
	      StrokeAndFillPath(hdc);*/
	      
	      SelectObject(hdc, GetStockObject(WHITE_BRUSH));
	      SelectObject(hdc, GetStockObject(BLACK_PEN));
	      MoveToEx(hdc,x,y,(LPPOINT)NULL);
	      LineTo(hdc,x+1,y+1); 

	    } 

	  i++;
	}
      SelectObject(hdc,hpenold);
      DeleteObject(hpen);
    }
  

  // EndPaint(hwnd, &ps); 

}

