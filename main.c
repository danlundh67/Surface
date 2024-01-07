#include <windows.h>  
#include "StrchBlt.h" 
#include <stdio.h>
#include <string.h>
#include "structures.h"

HINSTANCE hInst;   // current instance

LPCTSTR lpszAppName  = "MyApp";
LPCTSTR lpszTitle    = "My Application"; 
char filename[260];
struct Structure chosen_pdb;
struct Structure surface_pdb;
struct Structure org_pdb;
struct the_chain_report the_chains;
struct the_new_report allreport;
UINT chainid;


BOOL bChecked     = TRUE;
BOOL bRadio1      = TRUE;
BOOL bRadio2      = FALSE;
BOOL chainselected = FALSE;

/******************************************************/
/* Opens the standard file selection dialog in windos */
/* it returns the filename of the PDB one whiches to  */
/* Study. This is stored in the global variable file- */                                            
/* name that is used throuout the processing. A max   */
/* 260 chars is allowed including the path            */
/******************************************************/

char*  fileselect(HWND hwnd)
{
        OPENFILENAME ofn;       // common dialog box structure
        char szFile[260];       // buffer for file name
        HANDLE hf;              // file handle

        // Initialize OPENFILENAME
        ZeroMemory(&ofn, sizeof(ofn));
        ofn.lStructSize = sizeof(ofn);
        ofn.hwndOwner = NULL;
        ofn.lpstrFile = szFile;
        //
        // Set lpstrFile[0] to '\0' so that GetOpenFileName does not 
        // use the contents of szFile to initialize itself.
        //
        ofn.lpstrFile[0] = '\0';
        ofn.nMaxFile = sizeof(szFile);
        ofn.lpstrFilter = "PDB\0*.pdb\0ALL\0*.*\0";
        ofn.nFilterIndex = 1;
        ofn.lpstrFileTitle = NULL;
        ofn.nMaxFileTitle = 0;
        ofn.lpstrInitialDir = NULL;
        ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST | OFN_NOCHANGEDIR;

        // Display the Open dialog box. 

        if (GetOpenFileName(&ofn)==TRUE) 
          {
            /* hf = CreateFile(ofn.lpstrFile, GENERIC_READ,
                            0, (LPSECURITY_ATTRIBUTES) NULL,
                            OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL,
                            (HANDLE) NULL);
            */
            printf("%s\n",ofn.lpstrFile);
            /* CloseHandle(hf); */
            return (ofn.lpstrFile);
            
          }
        else
          return NULL;
}


/**********************************************************/
/* Dialog process for handling selection of chain and     */
/* obtaining a surface pdb. Sends message for plotting    */
/* and to write 80 chars of sec structure on selection of */
/* the chain.                                             */
/**********************************************************/ 

LRESULT CALLBACK TestDlgProc( HWND hDlg, UINT uMsg, 
                              WPARAM wParam, LPARAM lParam )
{
  HWND hvalue,scbar,myzoom;

  char dist[80];
  char *surfacefile;
  static SCROLLINFO si;
  int xPos;

  scbar=GetDlgItem(hDlg,IDC_ZOOM);
  myzoom=GetDlgItem(hDlg,IDC_ZOOMVAL);
  switch( uMsg )
    {
    case WM_INITDIALOG :
      {
	si.cbSize=sizeof(SCROLLINFO);
	si.fMask = SIF_RANGE | SIF_PAGE | SIF_POS;
	si.nMin=1;
	si.nMax=100;
	si.nPos=10;
	si.nPage=1;
	SetScrollInfo(scbar,SB_CTL,&si,TRUE);
	plotscale=1.000;
	/* sprintf(dist,"%3.1lf ",plotscale);
	   SetDlgItemText(hDlg,IDC_ZOOMVAL,(LPSTR)dist); */

	CheckRadioButton(hDlg,IDC_RADIO1, IDC_RADIO2,IDC_RADIO1);
	bRadio1=TRUE;
	bRadio2=FALSE;
	printf("INITDIALOG \n");
      }
      break;      
    case WM_COMMAND :
      switch( LOWORD( wParam ) )
	{
	case IDC_CHECKBOX :
	  bChecked = !IsDlgButtonChecked( hDlg, IDC_CHECKBOX );
	  CheckDlgButton( hDlg, IDC_CHECKBOX, 
			  bChecked ? MF_CHECKED : MF_UNCHECKED );
	  break;

	case IDC_CHAIN :
	  if ((HIWORD (wParam)==CBN_SELENDOK))
	    {		
	      hvalue=GetDlgItem(hDlg,IDC_CHAIN);
	      chainid=SendMessage(hvalue,CB_GETCURSEL,0,0);
	      printf("selected chain OK  %d \n",chainid);
	      fflush(stdout);
	      allreport=structure_to_chain(org_pdb,chainid);
	      
	      sprintf(dist,"%4.3lf ",allreport.subchains.HbondE);		    
	      SetDlgItemText(hDlg,IDC_SOLVE,(LPSTR)dist);
	      sprintf(dist,"%.0lf ",allreport.subchains.volsur.volume);
	      SetDlgItemText(hDlg,IDC_VOLUME,(LPSTR)dist);
	      sprintf(dist,"%.0lf ",allreport.subchains.volsur.surface);
	      SetDlgItemText(hDlg,IDC_SURFACE,(LPSTR)dist);

	      sprintf(dist,"%.0lf ",allreport.subchains.m.hbondv);
	      SetDlgItemText(hDlg,IDC_HBOND,(LPSTR)dist);

	      /* bfilename[0]='\0';
		 tmpstruct=create_xyzr_structure(translate_structure_onto_origin(chosen_pdb)); 
	      */

	      chainselected=TRUE;
	      surface_pdb=structure_onechain(org_pdb,chainid); 
	      
	      /* write_structure_to_pdb(tmpstruct ,"orginal.pdb"); 
		 write_structure_to_pdb(surface_pdb ,"surface.pdb");  */
	      /*
		strcpy(bfilename,"surfacefound"); */
	      read_to_list(surface_pdb); 

	      
	      /* Reinitialize the zoom bar */
	      /* and set plotscale to 1.00 */
	      si.cbSize=sizeof(SCROLLINFO);
	      si.fMask = SIF_RANGE | SIF_PAGE | SIF_POS;
	      si.nMin=1;
	      si.nMax=100;
	      si.nPos=10;
	      si.nPage=1;
	      SetScrollInfo(scbar,SB_CTL,&si,TRUE);
	      plotscale=1.000;
	      sprintf(dist,"%3.1lf ",plotscale);
	      SetDlgItemText(hDlg,IDC_ZOOMVAL,(LPSTR)dist);

	      /* SendMessage(swnd,WM_PAINT,0,0); */
	      SendMessage(GetParent(hDlg),WM_USER,(WPARAM)IDM_OPEN,(LPARAM)chainid); 
	      SendMessage(GetParent(hDlg),WM_USER,(WPARAM)IDC_SURFACE,0); 

	    }
	  break;
	case IDC_RADIO1 :
	  bRadio1 = TRUE;
	  bRadio2 = FALSE;
	  CheckRadioButton( hDlg, IDC_RADIO1, IDC_RADIO3,
			    IDC_RADIO1 );
	  /* SendMessage(swnd,WM_PAINT,0,0); */
	  if (chainselected==TRUE)
	    SendMessage(GetParent(hDlg),WM_USER,(WPARAM)IDC_SURFACE,0); 
	  break;
	case IDC_RADIO2 :
	  bRadio1 = FALSE;
	  bRadio2 = FALSE;
	  CheckRadioButton( hDlg, IDC_RADIO1, IDC_RADIO3,
			    IDC_RADIO2 );
	  if (chainselected==TRUE)
	    SendMessage(GetParent(hDlg),WM_USER,(WPARAM)IDC_SURFACE,0); 
	  break;
	case IDC_RADIO3 :
	  bRadio1= FALSE;
	  bRadio2= TRUE;
	  CheckRadioButton( hDlg, IDC_RADIO1, IDC_RADIO3,IDC_RADIO3);
	  if (chainselected==TRUE)
	    SendMessage(GetParent(hDlg),WM_USER,(WPARAM)IDC_SURFACE,0); 
	  break;	  
	case IDC_SOLVE :
	  break;
	  
	case IDCANCEL : 
	  DestroyWindow( hDlg ); 
	  break;
	}
      break;
    case WM_PAINT:
      {
	RECT rcClient;
	HDC hdc2;
	PAINTSTRUCT ps;
	HBITMAP hBitmap1,hBitmap2, hBitmap3, hBitmap4, hBitmap5,hBitmap6;	    

	GetClientRect(hDlg,&rcClient);
	InvalidateRect(hDlg,&rcClient,TRUE);
	hdc2=BeginPaint(hDlg,&ps);
	SetROP2(hdc2,R2_COPYPEN);




	/********************************************************************/
	/* THINGS for creating the different graphs - statistical view      */
	/********************************************************************/


	/*******************************************************/
	/* plot for the STAT1 figure (vol - length distr)      */
	/*******************************************************/
	hBitmap1 = LoadBitmap( hInst, "VOL1BMP" );
	plotgraphs(hDlg,hBitmap1, 1, IDC_ICON1,hdc2);
	DeleteObject( hBitmap1 );

	    
	/*******************************************************/
	/* REPEAT for the next STAT2 figure (vol/length distr) */
	/*******************************************************/
	    
	hBitmap2 = LoadBitmap( hInst, "VOL2BMP" );	    
	plotgraphs(hDlg, hBitmap2, 2, IDC_ICON2,hdc2);
	DeleteObject( hBitmap2 );


	/*******************************************************/
	/* REPEAT for the next STAT3 figure (surf -length )    */
	/*******************************************************/
	    
	hBitmap3 = LoadBitmap( hInst, "SURF1BMP" );	    
	plotgraphs(hDlg, hBitmap3, 3, IDC_ICON3,hdc2);
	DeleteObject( hBitmap3 );


	/*******************************************************/
	/* REPEAT for the next STAT3 figure (surf/length dist) */
	/*******************************************************/
	    
	hBitmap4 = LoadBitmap( hInst, "SURF2BMP" );	    
	plotgraphs(hDlg, hBitmap4, 4, IDC_ICON4,hdc2);
	DeleteObject( hBitmap4 );
	     

	/*****************************************************/
	/* REPEAT for the next STAT5 figure (hbond - length) */
	/*****************************************************/

	hBitmap5 = LoadBitmap( hInst, "HB1BMP" );
	plotgraphs(hDlg, hBitmap5, 5, IDC_ICON5,hdc2);
	DeleteObject( hBitmap5 );

	/*******************************************************/
	/* REPEAT for the next STAT6 figure (hbond/length distr) */
	/*******************************************************/
	      
	hBitmap6 = LoadBitmap( hInst, "HB2BMP" );	    	    
	plotgraphs(hDlg, hBitmap6, 6, IDC_ICON6,hdc2);
	DeleteObject( hBitmap6 );

	EndPaint(hDlg,&ps);
	ReleaseDC(hDlg,hdc2);
      }
      return 0;

    case WM_HSCROLL:
      {
      si.cbSize=sizeof(si);
      si.fMask= SIF_ALL;
      GetScrollInfo(scbar,SB_CTL,&si);
      xPos=si.nPos;
      switch (LOWORD(wParam))
	{
	case SB_LINELEFT:
	  si.nPos -=1;
	  break;
	case SB_LINERIGHT:
	  si.nPos +=1;
	  break;
	case SB_PAGERIGHT:
	  si.nPos +=si.nPage;
	  break;
	case SB_PAGELEFT:
	  si.nPos -=si.nPage;
	  break;
	case SB_THUMBTRACK:
	  si.nPos=si.nTrackPos;
	  break;
	default:
	  break;
	}
      // Set position and then retrieve it. 
      // May not be the same as the value set - windiows adjustment

      si.fMask=SIF_POS;
      SetScrollInfo(scbar,SB_CTL,&si,TRUE);
      GetScrollInfo(scbar,SB_CTL,&si);
      
      // if position change, rotate structure
      // since this is HSCROLL the rotation occurs 
      // on z-axix


      if (si.nPos != xPos)
	{

	  plotscale=(double)(float)(si.nPos/10.000);
	  sprintf(dist,"%3.1lf ",plotscale);
	  SetDlgItemText(hDlg,IDC_ZOOMVAL,(LPSTR)dist);
	  if (chainselected==TRUE)
	    SendMessage(GetParent(hDlg),WM_USER,(WPARAM)IDC_SURFACE,0); 
	}     
      }
      break;

    case WM_DESTROY :
      /* hDlgModeless = NULL; */
      DestroyWindow( hDlg );
      break;
    default :
      return( FALSE );
    }

  return( TRUE );
}


/*****************************************/
/* Functions for plotting the graphs     */
/*****************************************/

void plotgraphs(HWND hwnd, HBITMAP hBitmap1, int val, int where, HDC hDC)
{
  BITMAP bm1;
  HDC     hMemDC;
  RECT    rect;
  POINT   p;
  int     x1,y1,x2,y2, i ,j,z;
  char dista[80];

  //  hDC=GetDC(hwnd);
  GetObject( hBitmap1, sizeof( BITMAP ), &bm1 );
  hMemDC = CreateCompatibleDC( hDC );	    
  SelectObject( hMemDC, hBitmap1 );	  
  
  if (GetStretchBltMode(hDC) != HALFTONE)
    SetStretchBltMode( hDC, HALFTONE );
  GetClientRect(GetDlgItem(hwnd,where),&rect);
  p.x=rect.left;
  p.y=rect.top;
  ClientToScreen(GetDlgItem(hwnd,where),&p);
  ScreenToClient(hwnd,&p);   
  z=0;
  if (chainselected==TRUE)
    {      
      if (val==1)
	{
	  /* setting cross in volume - length plot */
	  x1=82+(int)the_chains.subchains[chainid].chain_length*31.0/200.0;  
	  // GetDlgItemText(pwnd,IDC_VOLUME,(LPSTR)dista,80);
	  sprintf(dista,"%.0lf ",allreport.subchains.volsur.volume);
	  y1=160-(int)(atoi(dista))*29.5/50000.0;
	  z=1;
	}
      else if (val==2)
	{
	  /* setting cross in volume/length distr plot */
	  x2=(int)(the_chains.subchains[chainid].chain_length);  
	  // GetDlgItemText(pwnd,IDC_VOLUME,(LPSTR)dista,80);
	  sprintf(dista,"%.0lf ",allreport.subchains.volsur.volume);
	  y2=(int)(atoi(dista));
	  x1=(int)(((float)y2/x2-126.500)*159.0/33)+20; 
	  // y1=160;
	  y1=160;
	  for (i=0;i<bm1.bmHeight;i++)
	    {
	    if (GetPixel(hMemDC,x1,i)==RGB(255,0,0))
	      {
		y1=i;
	      }
	    }
	   if (y1==160)
	    {
	      if ((x1>20)&&(x1<159))
		for (i=0;i<bm1.bmHeight;i++)
		  if (GetPixel(hMemDC,x1+1,i)==RGB(255,0,0))
		    y1=i;
	    }

	  z=2;
	}
      else if (val==3)
	{
	  /* setting cross in surface length plot */
	  sprintf(dista,"%.0lf ",allreport.subchains.volsur.surface);
	  y2=(int)(atoi(dista));

	  x1=75;
	  y1=160;
	  // y1=13;
	  x1=330;
	  x1=75+(int)((float)the_chains.subchains[chainid].chain_length/1600.0*(329-75));
	  y1= 160-(int)((float)y2/70000.0*(160-13));
	}
      else if (val==4)
	{
	  /* setting cross in surface/length dist plot */
	  sprintf(dista,"%.0lf ",allreport.subchains.volsur.surface);
	  y2=(int)(atoi(dista));
	  x2=(int)(the_chains.subchains[chainid].chain_length);
	  x1=20+(int)(((float)y2/x2-29)*159/(130-29));
	  y1=160;
	  for (i=0;i<bm1.bmHeight;i++)
	    {
	    if (GetPixel(hMemDC,x1,i)==RGB(255,0,0))
	      {
		y1=i;
	      }
	    }
	   if (y1==160)
	    {
	      if ((x1>20)&&(x1<159))
		for (i=0;i<bm1.bmHeight;i++)
		  if (GetPixel(hMemDC,x1+1,i)==RGB(255,0,0))
		    y1=i;
	    }
	  
	}
      else if (val==5)
	{
	  /* setting cross in hbond energy length plot */
	  x1=68+(int)the_chains.subchains[chainid].chain_length*(33.0/200.0-0.2/200);
	  sprintf(dista,"%.0lf ",allreport.subchains.m.hbondv);		    
	  //GetDlgItemText(hchng,IDC_HBOND,(LPSTR)dista,80);
	  y1=13+(int)(-1*atoi(dista))*(15.0/200.0-0.2/200);
	  // x1=68;
	  // y1=13;

	}
      else if (val==6)
	{
	  /* setting cross in hbond energy-length dist plot */
	  sprintf(dista,"%.0lf ",allreport.subchains.m.hbondv);		    
	  x1=20+(int)(((float)atoi(dista)/the_chains.subchains[chainid].chain_length+1.8)*159/1.1);

	  y1=160;
	  for (i=0;i<bm1.bmHeight;i++)
	    if (GetPixel(hMemDC,x1,i)==RGB(255,0,0))
	      y1=i;
	  if (y1==160)
	    {
	      if ((x1>20)&&(x1<159))
		for (i=0;i<bm1.bmHeight;i++)
		  if (GetPixel(hMemDC,x1+1,i)==RGB(255,0,0))
		    y1=i;
	    }
	}
      /* testa så att dessa x och y ligger inom bilden,    */
      /* annars så skriv ut felmeddelande alt. inget kryss */

      for (i=(x1-5);i<(x1+5);i++)
	SetPixel(hMemDC,i,y1,RGB(0,0,255));
      for (j=(y1-5);j<(y1+5);j++)
	SetPixel(hMemDC,x1,j,RGB(0,0,255));
    }

  StretchBlt( hDC,
	      p.x+3,
	      p.y+3,
	      (rect.right-rect.left-5),
	      (rect.bottom-rect.top-5),
	      hMemDC, 0, 0, bm1.bmWidth, bm1.bmHeight,
	      SRCCOPY );
	   
  DeleteDC( hMemDC );
  // ReleaseDC( hwnd, hDC );
}







/*******************************************************/
/* The main function, i.e. register a window and set   */
/* appropriate window procedure. The window procedure  */
/* then calls the different dialogs such that input    */
/* and display information is carried out              */
/*******************************************************/

int APIENTRY WinMain( HINSTANCE hInstance, HINSTANCE hPrevInstance,
                      LPTSTR lpCmdLine, int nCmdShow)
{
   MSG      msg;
   HWND     hWnd; 
   WNDCLASSEX wc;

   // Register the main application window class.
   //............................................
   wc.style         = CS_HREDRAW | CS_VREDRAW;
   wc.lpfnWndProc   = (WNDPROC)WndProc;       
   wc.cbClsExtra    = 0;                      
   wc.cbWndExtra    = 0;                      
   wc.hInstance     = hInstance;              
   wc.hIcon         = LoadIcon( hInstance, lpszAppName ); 
   wc.hCursor       = LoadCursor(NULL, IDC_ARROW);
   wc.hbrBackground = (HBRUSH)(COLOR_WINDOW+1);
   wc.lpszMenuName  = lpszAppName;              
   wc.lpszClassName = lpszAppName;              
   wc.cbSize        = sizeof( WNDCLASSEX );
   wc.hIconSm       = LoadImage( hInstance, lpszAppName,
                                 IMAGE_ICON, 16, 16,
                                 LR_DEFAULTCOLOR );

   if ( !RegisterClassEx( &wc ) )
      return( FALSE );

   hInst = hInstance; 

   // Create the main application window.
   //....................................
   hWnd = CreateWindow( lpszAppName, 
                        lpszTitle,    
                        WS_OVERLAPPEDWINDOW, 
                        CW_USEDEFAULT, 0, 
                        CW_USEDEFAULT, 0,  
                        NULL,              
                        NULL,              
                        hInstance,         
                        NULL               
                      );

   if ( !hWnd ) 
      return( FALSE );

   ShowWindow( hWnd, nCmdShow ); 
   UpdateWindow( hWnd );         

   while( GetMessage( &msg, NULL, 0, 0) )   
   {
      TranslateMessage( &msg ); 
      DispatchMessage( &msg );  
   }

   return( msg.wParam ); 
}



/*****************************************/
/* Windows procedure for handling the    */
/* main window events.                   */
/*****************************************/


LRESULT CALLBACK WndProc( HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam )
{

  int i,j,x,v,yPos;
  static HWND hchng=NULL;
  HWND hlist=NULL;
  static HWND zoominfo=NULL;
  HWND hvalue,scbar,myzoom;
  char zoomf[60];
  char dist[80];
  int maxstring;
  char c;
  static HWND secdial=NULL;
  char dista[81], dista2[81];
  static int chainid;
  char *outbuffer=NULL;
  static HWND swnd=NULL;
  static SCROLLINFO si;
  RECT myrc;
  


  /* Variables for plotting surface  */
  
  static RECT rcClient; // client area rectangle
  HDC hdc,hdc2;          // handle to device context (DC) 
  PAINTSTRUCT ps;   // paint data for Begin/EndPaint 

  /* rotation variable */
  double theta,distance;




  switch( uMsg )
    {
    case WM_CREATE: 
      {
      /***********************************************/
      /* inserts a dialogbox in the parameter window */
      /* the layout is specified in test4.rc         */
      /* it is also resized to fit the window        */
      /***********************************************/
      //(DLGPROC)TestDlgProc 
      hchng= CreateDialog( hInst, "TestDialog", hWnd,(DLGPROC)TestDlgProc);
      GetClientRect(hWnd, &myrc);
      printf("client area  %ld  %ld \n",myrc.right,myrc.bottom);
      SetWindowPos(hchng,HWND_TOP,0,0,myrc.right-1,myrc.bottom-1,SWP_DRAWFRAME);

      //(DLGPROC)ChildSWndProc

      /***********************************************/
      /* inserts a dialogbox in the parameter window */
      /* the layout is specified in StrchBlt.rc      */
      /* it is also resized to fit the window        */
      /***********************************************/

      }
      break;
    case WM_SIZE:
      {
	GetClientRect(hWnd, &myrc);
	SetWindowPos(hchng,HWND_TOP,0,0,myrc.right-1,myrc.bottom-1,SWP_DRAWFRAME);
      }
      break;
    case WM_COMMAND :
      switch( LOWORD( wParam ) )
	{
	case IDM_TEST :
	  {
	  }
	  break;
	case IDM_FILE :
	  {
	    strcpy(filename,fileselect(hWnd));
	    if (filename!=NULL)
	      {
		SetDlgItemText(hchng,IDC_URL1,(LPSTR)filename);
		fflush(stdout);
		free(dotatm);
		dotatm=NULL;
		chosen_pdb=read_pdb_to_structure(filename);

		/* Calculate number of chains and allocate reportstructure 
		   call struct_to_chain which fills reportstructure */
		org_pdb = create_xyzr_structure(translate_structure_onto_origin(chosen_pdb));

		the_chains=structure_chains(org_pdb);

		hlist=GetDlgItem(hchng,IDC_CHAIN);
		SendMessage(hlist,CB_RESETCONTENT,0,0);
	  

		for(i=0;i<the_chains.nrch;i++)
		  {
		    if ((the_chains.nrch==1)&&(strcmp(the_chains.subchains[i].chainID,"A")!=0))
		      sprintf(dist,"A: %d",the_chains.subchains[i].chain_length);
		    else 
		      sprintf(dist,"%s: %d",the_chains.subchains[i].chainID,the_chains.subchains[i].chain_length);
		    fflush(stdout);
		    SendMessage(hlist,CB_ADDSTRING,0,(LPARAM)dist);
		    SendMessage(hlist,CB_SETITEMDATA,i,(LPARAM)i);

		  }

		//SendMessage(hWnd,WM_USER,(WPARAM)IDM_TEST,0);  

		/* SendMessage(swnd,WM_PAINT,0,0);  */
	      }

	  }
	  return 0;
	case IDM_ABOUT :
	  DialogBox( hInst, "AboutBox", hWnd, (DLGPROC)About );
	  break;

	case IDM_EXIT :
	  DestroyWindow( hWnd );
	  break;
	} /* WM_COMMAND switch Loword wParam */
      break;
    case WM_USER:
      switch(wParam) 
	{ 	  
	case IDM_OPEN:
	  {
	    chainid=(int)(lParam);
	    printf("feature window WM_USER %d\n",(int)(lParam));
	    // sprintf(dist,"chain %d",lParam);
	    if (strlen(allreport.subchains.m.secstruct)>80)
	      {
		free(outbuffer);
		outbuffer=(char *)malloc(sizeof(char)*(strlen(allreport.subchains.m.secstruct)+1));
		sprintf(outbuffer,"%s",allreport.subchains.m.secstruct);
		for (i=0;i<80;i++)
		  {
		    c=(char)outbuffer[i];
		    dista[i]=c;
		  }
		dista[80]='\0';
	      }
	    else
	      {
		sprintf(dista,"%s",allreport.subchains.m.secstruct);
	      } 
	    free(outbuffer);
	    SetDlgItemText(secdial,IDC_SECSTR,(LPSTR)dista); 		  
	    maxstring=0;

	  } /* IDM_OPEN */
	  break;
	case IDC_SECSCROLL:
	  {
	    j=(int)lParam-1;
	    /* printf("secscroll  %d \n",(int)lParam-1);
	       fflush(stdout); */
	    i=0;
	    free(outbuffer);
	    outbuffer=(char *)malloc(sizeof(char)*(strlen(allreport.subchains.m.secstruct)+1));
	    
	    /********************************************************************/
	    /* this may seem unneccessary but since the we are dealing with 
	       events we need to copy the sec str again since in this event all
	       previous data is lost besides the chainid */
	    /********************************************************************/
	  
	    sprintf(outbuffer,"%s",allreport.subchains.m.secstruct);
	    // printf("scroll pos %d %d %d %d \n",i,j,strlen(outbuffer),strlen(outbuffer)-j);
	    
	    if (strlen(outbuffer)>80)
	      {
		x=strlen(outbuffer)-j;
		if (x>80)
		  {
		    x=80;
		    for (i=0;i<x;i++)
		      {
			if (i<(strlen(outbuffer)-j))
			  {
			    c=(char)outbuffer[j+i];
			    dista[i]=c;
			  }
		      }
		    // printf("%d %c \n",(x),dist[x-1]);
		    dista[x-1]='\0';
		  }
		else
		  {
		    j=0;
		    v=strlen(outbuffer);
		    x=(strlen(outbuffer)-80);
		    // printf("%d %d %d \n",x,strlen(outbuffer),v);
		    for (i=x;i<=v;i++)
		      {
			c=(char)outbuffer[i];
			dista[j]=c;
			j++;
		      }
		    dista[j]='\0';
		    //  printf("1 %s \n",dista);
		  }
	      }
	    else
	      {
		sprintf(dista,"%s",outbuffer);
	      }
	    
	    // printf("2 %s \n",dista);
	    SetDlgItemText(secdial,IDC_SECSTR,(LPSTR)dista);
	  } /* IDE_SECSCROLL */
	  break;
	case IDC_SURFACE :
	  {
	    /* SendMessage(swnd,WM_PAINT,0,0); */
	    SendMessage(secdial,WM_PAINT,0,0);

	  }
	  break;
	} /* WM_USER wParam*/
      break;
    case WM_VSCROLL:    
      return 0;
    case WM_DESTROY :
      PostQuitMessage(0);
      break;

    default :
      return( DefWindowProc( hWnd, uMsg, wParam, lParam ) );
    }

  return( 0L );
}


LRESULT CALLBACK About( HWND hDlg,           
                        UINT message,        
                        WPARAM wParam,       
                        LPARAM lParam)
{
   switch (message) 
   {
       case WM_INITDIALOG: 
               return (TRUE);

       case WM_COMMAND:                              
               if (   LOWORD(wParam) == IDOK         
                   || LOWORD(wParam) == IDCANCEL)    
               {
                       EndDialog(hDlg, TRUE);        
                       return (TRUE);
               }
               break;
   }

   return (FALSE); 
}
