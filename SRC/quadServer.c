#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>		// usleep

#include "wsserver.h"

#ifdef WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#ifndef snprintf
#define snprintf _snprintf
#endif
#endif
#include <winsock2.h>
#endif


int main(int argc, char *argv[])
{
  int       round, i, j, f, stat, *nvert = NULL , *nquad = NULL, quad[4],x,
      vt = 0, qt = 0, *segs = NULL, *tris = NULL, *nf = NULL, *piv = NULL, min;
  float     focus[4], colorLine[3], colorBox[3];
  double    size, box[6], bigBox[6],  *xyzs = NULL;
  char      gpname[34], *startapp, tmp[256], str[256];
  wvContext *cntxt;
  wvData    items[5];
  float     eye[3]    = {0.0, 0.0, 7.0};
  float     center[3] = {0.0, 0.0, 0.0};
  float     up[3]     = {0.0, 1.0, 0.0};

  FILE      *fp;
  /* get our starting application line
   *
   * for example on a Mac:
   * setenv WV_START "open -a /Applications/Firefox.app ../client/wv.html"
   */
  startapp = getenv("WV_START");

  if (argc < 2 ) {
      printf("\n Usage: quadServer ndatafile(s) RGB colorbox (3 floats) \n\n");
      return 1;
  }
  printf(" Color Boxes: 1 = green 2 = acdl\n\n");
  scanf ("%d", &i );
  switch(i) {
    case 1:
       printf(" YOU HAVE CHOSEN GREENS \n ");
       colorBox[0]  = 0.9;
       colorBox[1]  = 1.0;
       colorBox[2]  = 0.4;
       colorLine[0] = 0.2;
       colorLine[1] = 0.0;
       colorLine[2] = 0.2;
       break;
    case 2:
    default:
       printf(" YOU HAVE CHOSEN ACDL\n ");
       colorBox[0]  = 255.0/255.0;
       colorBox[1]  = 255.0/255.0;
       colorBox[2]  = 255.0/255.0;
       colorLine[0] = 138.0/255.0;
       colorLine[1] = 23.0/255.0;
       colorLine[2] = 50.0/255.0;
       break;

  }


  nvert = (int *)malloc(argc * sizeof (int));
  nquad = (int *)malloc(argc * sizeof (int));
  nf    = (int *)malloc(argc * sizeof (int));
  piv   = (int *)malloc(argc * sizeof (int));
  if ( nvert == NULL || nquad == NULL || nf == NULL || piv == NULL) {
      printf(" Error allocating memo \n ");
      return 1;
  }
  /* create the WebViewer context */
  cntxt = wv_createContext(1, 30.0, 1.0, 10.0, eye, center, up);
  if (cntxt == NULL) {
      printf(" failed to create wvContext!\n");
      return 1;
  }
  for ( f = 1 ; f < argc; f++ ) {
      tmp[0]='\0';
      snprintf(str, 256,"%s", argv[f]);
      while (sscanf(str,"%[^0123456789]%s",tmp, str) > 1
	  ||sscanf(str,"%d%s",&x,str))
	{
	  if (tmp[0]=='\0') nf[f] = x;
	  tmp[0]='\0';
	}
      piv[f] = f;
  }
  for ( f = 1 ; f < argc; f++ ) {
      for ( i = 1; i < argc - f ; i++ ) {
	  if ( nf[piv[i]] > nf[piv[i + 1]]) {
	      min      = piv[i    ];
	      piv[i ]  = piv[i + 1];
	      piv[i+1] = min;
	  }
      }
  }
  for ( round = 0 ; round < 2; round++) {
      for ( f = 1 ; f < argc; f++ ) {
	  fp = fopen(argv[piv[f]], "r");
	  printf("openfile %s\n ", argv[piv[f]]);
	  if (fp == NULL) {
	      printf("\n ERROR: Opening %s!\n\n", argv[piv[f]]);
	      return 1;
	  }
	  j = fscanf(fp, "%d %d", &nvert[f], &nquad[f]);
	  if (j != 2) {
	      printf("\n ERROR: reading header\n\n");
	      fclose(fp);
	      return 1;
	  }
	  vt  += nvert[f];
	  qt  += nquad[f];
	  printf(" FACE %d ---> QUADS %d\n ", piv[f], nquad[f]);
	  xyzs = (double *) malloc(3*nvert[f]*sizeof(double));
	  tris = (int    *) malloc(6*nquad[f]*sizeof(int));
	  segs = (int    *) malloc(8*nquad[f]*sizeof(int));
	  if ((xyzs == NULL) || (tris == NULL) || (segs == NULL)) {
	      printf("\n ERROR: malloc!\n\n");
	      if (xyzs != NULL) free(xyzs);
	      if (tris != NULL) free(tris);
	      if (segs != NULL) free(segs);
	      fclose(fp);
	      return 1;
	  }
	  /* get the verts and the bounding box */
	  box[0] = box[1] = box[2] =  1.e100;
	  box[3] = box[4] = box[5] = -1.e100;
	  for (i = 0; i < nvert[f]; i++) {
	      j = fscanf(fp, "%lf %lf %lf", &xyzs[3*i  ], &xyzs[3*i+1], &xyzs[3*i+2]);
	      if (j != 3) {
		  printf("\n ERROR: reading vert #%d!\n\n", i+1);
		  free(xyzs);
		  free(tris);
		  free(segs);
		  fclose(fp);
		  return 1;
	      }
	      if (xyzs[3*i  ] < box[0]) box[0] = xyzs[3*i  ];
	      if (xyzs[3*i  ] > box[3]) box[3] = xyzs[3*i  ];
	      if (xyzs[3*i+1] < box[1]) box[1] = xyzs[3*i+1];
	      if (xyzs[3*i+1] > box[4]) box[4] = xyzs[3*i+1];
	      if (xyzs[3*i+2] < box[2]) box[2] = xyzs[3*i+2];
	      if (xyzs[3*i+2] > box[5]) box[5] = xyzs[3*i+2];
	  }
	  if ( round == 0 ) {
	      if ( f == 1 )
		for ( j = 0 ; j < 6; j++ ) bigBox[j] = box[j];
	      else {
		  if ( box[0] < bigBox[0] ) bigBox[0] = box[0];
		  if ( box[1] < bigBox[1] ) bigBox[1] = box[1];
		  if ( box[2] < bigBox[2] ) bigBox[2] = box[2];
		  if ( box[3] > bigBox[3] ) bigBox[3] = box[3];
		  if ( box[4] > bigBox[4] ) bigBox[4] = box[4];
		  if ( box[5] > bigBox[5] ) bigBox[5] = box[5];
	      }
	  }
	  if ( round == 0 ) continue;
	  for ( j = 0 ; j < 6; j++ ) box[j] = bigBox[j];
	  size = box[3]-box[0];
	  if (size < box[4]-box[1]) size = box[4]-box[1];
	  if (size < box[5]-box[2]) size = box[5]-box[2];
	  focus[0] = 0.5*(box[0]+box[3]);
	  focus[1] = 0.5*(box[1]+box[4]);
	  focus[2] = 0.5*(box[2]+box[5]);
	  focus[3] = size;
	  /* get the quads; make triangles and line segments */
	  for (i = 0; i < nquad[f]; i++) {
	      j = fscanf(fp, "%d %d %d %d", &quad[0], &quad[1], &quad[2], &quad[3]);
	      if (j != 4) {
		  printf("\n ERROR: reading quad #%d!\n\n", i+1);
		  free(xyzs);
		  free(tris);
		  free(segs);
		  fclose(fp);
		  return 1;
	      }
	      tris[6*i  ] = quad[0];
	      tris[6*i+1] = quad[1];
	      tris[6*i+2] = quad[2];
	      tris[6*i+3] = quad[0];
	      tris[6*i+4] = quad[2];
	      tris[6*i+5] = quad[3];
	      segs[8*i  ] = quad[0];
	      segs[8*i+1] = quad[1];
	      segs[8*i+2] = quad[1];
	      segs[8*i+3] = quad[2];
	      segs[8*i+4] = quad[2];
	      segs[8*i+5] = quad[3];
	      segs[8*i+6] = quad[3];
	      segs[8*i+7] = quad[0];
	  }
	  fclose(fp);
	  /* make the scene */
          if ( f > 1 && nf[piv[f]] == nf[piv[f-1]])
	  snprintf(gpname, 34, "Body %d Face %d%d", 1, nf[piv[f]], f);
          else
            snprintf(gpname, 34, "Body %d Face %d", 1, nf[piv[f]]);
	  stat = wv_setData(WV_REAL64, nvert[f], (void *) xyzs,  WV_VERTICES, &items[0]);
	  if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", stat, gpname);
	  wv_adjustVerts(&items[0], focus);
	  stat = wv_setData(WV_INT32, 6*nquad[f], (void *) tris, WV_INDICES, &items[1]);
	  if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", stat, gpname);
	  /* set the foreground color -- red */
	  stat = wv_setData(WV_REAL32, 1, (void *) colorBox,  WV_COLORS, &items[2]);
	  if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", stat, gpname);
	  stat = wv_setData(WV_INT32, 8*nquad[f], (void *) segs, WV_LINDICES, &items[3]);
	  if (stat < 0) printf(" wv_setData = %d for %s/item 3!\n", stat, gpname);
	  /* line color -- black */
	  stat = wv_setData(WV_REAL32, 1, (void *) colorLine,  WV_LCOLOR, &items[4]);
	  if (stat < 0) printf(" wv_setData = %d for %s/item 4!\n", stat, gpname);
	  stat = wv_addGPrim(cntxt, gpname, WV_TRIANGLE,
			     WV_ON|WV_ORIENTATION, 5, items);
	  if (stat < 0)
	    printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
	  free (xyzs);
	  free (segs);
	  free (tris);
	  xyzs = NULL;
	  segs = NULL;
	  tris = NULL;
      }
  }
  /* start the server code */

  stat = 0;
  if (wv_startServer(7681, NULL, NULL, NULL, 0, cntxt) == 0) {

      /* we have a single valid server -- stay alive a long as we have a client */
      while (wv_statusServer(0)) {
	  usleep(500000);
	  if (stat == 0) {
	      if (startapp != NULL) system(startapp);
	      stat++;
	  }
      }
  }
  wv_cleanupServers();

  /* finish up */
  free(xyzs);
  free(tris);
  free(segs);
  return 0;
}


/* call-back invoked when a message arrives from the browser */

void browserMessage(/*@unused@*/ void *wsi, /*@unused@*/ char *text,
		    /*@unused@*/ int  lena)
{

}
