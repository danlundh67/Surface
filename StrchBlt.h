#define IDM_EXIT           100
#define IDM_TEST           200
#define IDM_ABOUT          301
#define IDM_FILE           500
#define IDM_OPEN           501


#define ID_EXIT                         3000
#define IDC_EXIT                        3000
#define IDC_Download                    3001
#define IDC_URL1                        3002
#define IDC_Header1                     3003
#define IDC_Resource1                   3004
#define IDC_LIST1                       3005
#define IDC_CallbackList                3005
#define IDC_URL2                        3011
#define IDC_Header2                     3012
#define IDC_Resource2                   3013

#define IDC_CHECKBOX                    3014
#define IDC_CHAIN                       3015

#define IDC_CANCEL                      3017
#define IDC_SOLVE                       3018
#define IDC_VOLUME                      3019
#define IDC_SURFACE                     3020
#define IDC_ZOOM                        3021
#define IDC_ZOOMVAL                     3022
#define IDC_HBOND                       3023
#define IDC_RADIO1                      400
#define IDC_RADIO2                      401
#define IDC_RADIO3                      402

/* included defs for feature dialog */

#define IDC_SECSTR                      2501
#define IDC_SECSCROLL                   2502
#define IDC_ICON1                       2503
#define IDC_ICON2                       2504
#define IDC_ICON3                       2505
#define IDC_ICON4                       2506
#define IDC_ICON5                       2507
#define IDC_ICON6                       2508
#define IDC_ICON7                       2507
#define IDC_ICON8                       2508


/* included defs for struct window */

#define IDC_HSCROLL                      2701
#define IDC_VSCROLL                      2702
#define IDC_PLOT                         2703


LRESULT CALLBACK WndProc  (HWND, UINT, WPARAM, LPARAM);
LRESULT CALLBACK About    (HWND, UINT, WPARAM, LPARAM);
LRESULT CALLBACK TestDlgProc( HWND, UINT, WPARAM, LPARAM);
LRESULT CALLBACK FeatDlgProc( HWND, UINT, WPARAM, LPARAM);
LRESULT CALLBACK ChildSWndProc( HWND, UINT, WPARAM, LPARAM);
void plotgraphs(HWND, HBITMAP, int, int, HDC);
