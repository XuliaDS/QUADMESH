/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             EGADS Tessellation using wv with Quad Tessellation
 *
 *      Copyright 2011-2018, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <math.h>
#include <string.h>
#include <unistd.h>		// usleep
#include <nlopt.h>
#include <time.h>
#ifdef  WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <winsock2.h>
#endif
#include "egads.h"
#include "wsserver.h"
int TOTSWAPS     = 0;
int OVERALLSWAPS = 0;
#define EPS  1.E-03
#define DEPS 1.E-03
//#define DEBUG_AREA
int NOUT = 10000;
#define NLOPTMAXEVAL 100000

#define DTOL 1.e-9
#define TRIPLEOPTI 1
#define PI              3.1415926535897931159979635
#define DEBUG
int GLO = 0 ;
int QUADTESS  = 0;
int OPTCOUNT  = 0;
int TESSCOUNT = 0;
#define ANGPASS 0.3
#define ERRCTT  3.554147

int SCOUNT = 0;
#define FORCETRIANGULATION
#define MAX_LINKS 100
#define FORCEQUADS
/* structure to hold on to the EGADS triangulation per Body */
typedef struct {
  ego    *faces;
  ego    *edges;
  ego    body;
  ego    tess;
  int    mtype;
  int    nfaces;
  int    nedges;
  int    plen;
  double *uvs;      /* face #1 updated uvs */
} bodyData;


/* globals used in these functions */
static wvContext *cntxt;
static bodyData  *bodydata;
static int       nbody;
static float     focus[4];



int main(int argc, char *argv[])
{
  int          i, j, k, ibody, stat, oclass, mtype, len, ntri, sum;
  int          nseg, ngp, atype, alen, quad, *segs, *senses;
  const int    *tris, *tric, *ptype, *pindex, *ints;
  float        arg, colorBox[3], colorLine[3];
  double       box[6], size, params[3];
  const double *xyzs, *uvs, *ts, *reals;
  char         gpname[66], *startapp;
  const char   *OCCrev, *string;
  ego          context, tess, model, geom, *bodies, *dum;
  wvData       items[5];
  float        eye[3]      = {0.0, 0.0, 7.0};
  float        center[3]   = {0.0, 0.0, 0.0};
  float        up[3]       = {0.0, 1.0, 0.0};
  static int   sides[3][2] = {{1,2}, {2,0}, {0,1}};
  static int   sideq[4][2] = {{1,2}, {2,5}, {5,0}, {0,1}};
  static int   neigq[4]    = {  0,     3,     4,     2};
  FILE *fil = NULL;
  char fn[100];
  /* get our starting application line
   *
   * for example on a Mac:
   * setenv WV_START "open -a /Applications/Firefox.app ../client/wv.html"
   */
  startapp = getenv("WV_START");

  if ((argc != 2) && (argc != 5)) {
    printf("\n Usage: vAttr filename [angle maxlen sag]\n\n");
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



  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %2d.%02d with %s\n\n", i, j, OCCrev);

  /* initialize */
  printf(" EG_open           = %d\n", EG_open(&context));
  stat = 0;
  if(argv[1])  stat = EG_loadModel(context, 0, argv[1], &model);
  printf(" EG_loadModel      = %d\n", stat);
  printf(" EG_getBoundingBox = %d\n", EG_getBoundingBox(model, box));
  printf("       BoundingBox = %lf %lf %lf\n", box[0], box[1], box[2]);
  printf("                     %lf %lf %lf\n", box[3], box[4], box[5]);
  printf(" \n");

  size = box[3]-box[0];
  if (size < box[4]-box[1]) size = box[4]-box[1];
  if (size < box[5]-box[2]) size = box[5]-box[2];

  focus[0] = 0.5*(box[0]+box[3]);
  focus[1] = 0.5*(box[1]+box[4]);
  focus[2] = 0.5*(box[2]+box[5]);
  focus[3] = size;

  /* get all bodies */
  stat = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody,
      &bodies, &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getTopology = %d\n", stat);
    return 1;
  }
  printf(" EG_getTopology:     nBodies = %d\n", nbody);
  bodydata = (bodyData *) malloc(nbody*sizeof(bodyData));
  if (bodydata == NULL) {
    printf(" MALLOC Error on Body storage!\n");
    return 1;
  }

  params[0] =  0.1 *size;
  params[1] =  0.01*size;
  params[2] =  20.0;
  if (argc == 5) {
    sscanf(argv[2], "%f", &arg);
    params[2] = arg;
    sscanf(argv[3], "%f", &arg);
    params[0] = arg;
    sscanf(argv[4], "%f", &arg);
    params[1] = arg;
    params[0] *= size;
    params[1] *= size;
  }

  /* fill our structure a body at at time */
  for (ibody = 0; ibody < nbody; ibody++) {
#ifdef FORCETRIANGULATION
    stat = EG_attributeAdd(bodies[ibody], ".qParams", ATTRSTRING, 4, NULL, NULL,
        "off");
    if (stat != EGADS_SUCCESS)
      printf(" Body %d: attributeAdd = %d\n", ibody, stat);
#endif


    mtype = 0;
    EG_getTopology(bodies[ibody], &geom, &oclass,
        &mtype, NULL, &j, &dum, &senses);
    bodydata[ibody].body  = bodies[ibody];
    bodydata[ibody].mtype = mtype;
    bodydata[ibody].uvs   = NULL;
    if (mtype == WIREBODY) {
      printf(" Body %d: Type = WireBody\n", ibody+1);
    } else if (mtype == FACEBODY) {
      printf(" Body %d: Type = FaceBody\n", ibody+1);
    } else if (mtype == SHEETBODY) {
      printf(" Body %d: Type = SheetBody\n", ibody+1);
    } else {
      printf(" Body %d: Type = SolidBody\n", ibody+1);
    }

    stat = EG_getBodyTopos(bodies[ibody], NULL, FACE,
        &bodydata[ibody].nfaces, &bodydata[ibody].faces);
    i    = EG_getBodyTopos(bodies[ibody], NULL, EDGE,
        &bodydata[ibody].nedges, &bodydata[ibody].edges);
    if ((stat != EGADS_SUCCESS) || (i != EGADS_SUCCESS)) {
      printf(" EG_getBodyTopos Face = %d\n", stat);
      printf(" EG_getBodyTopos Edge = %d\n", i);
      return 1;
    }

    stat = EG_makeTessBody(bodies[ibody], params, &bodydata[ibody].tess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_makeTessBody %d = %d\n", ibody, stat);
      continue;
    }
    stat = EG_makeTessBody(bodies[ibody], params, &bodydata[ibody].tess);

#ifdef FORCEQUADS
    tess = bodydata[ibody].tess;
    stat = EG_quadTess(tess, &bodydata[ibody].tess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_quadTess %d = %d  -- reverting...\n", ibody, stat);
      bodydata[ibody].tess = tess;
      continue;
    }
    EG_deleteObject(tess);
    QUADTESS = 1;
#endif
  }
  printf(" Using angle = %lf,  relSide = %lf,  relSag = %lf\n",
          params[2], params[0], params[1]);
  printf(" \n");

  /* create the WebViewer context */
  cntxt = wv_createContext(1, 30.0, 1.0, 10.0, eye, center, up);
  if (cntxt == NULL) {
    printf(" failed to create wvContext!\n");
    for (ibody = 0; ibody < nbody; ibody++) {
      EG_deleteObject(bodydata[ibody].tess);
      EG_free(bodydata[ibody].edges);
      EG_free(bodydata[ibody].faces);
    }
    free(bodydata);

    printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
    printf(" EG_close          = %d\n", EG_close(context));
    return 1;
  }

  /* make the scene */
  for (ngp = sum = stat = ibody = 0; ibody < nbody; ibody++) {
    quad = 0;
    stat = EG_attributeRet(bodydata[ibody].tess, ".tessType", &atype,
        &alen, &ints, &reals, &string);
    if (stat == EGADS_SUCCESS)
      if (atype == ATTRSTRING)
        if (strcmp(string, "Quad") == 0) quad = 1;

    /* get faces */
    for (i = 0; i < bodydata[ibody].nfaces; i++) {
	stat = EG_getTessFace(bodydata[ibody].tess, i+1, &len,
			      &xyzs, &uvs, &ptype, &pindex, &ntri,
			      &tris, &tric);
	bodydata[ibody].plen = len;
	bodydata[ibody].uvs  = (double *) EG_alloc(2*len*sizeof(double));
	if (bodydata[ibody].uvs != NULL)
	  for (j = 0; j < 2*len; j++) bodydata[ibody].uvs[j] = uvs[j];
	if (stat != EGADS_SUCCESS) continue;
	sprintf(gpname, "Body %d Face %d", ibody+1, i+1);
	stat = wv_setData(WV_REAL64, len, (void *) xyzs,  WV_VERTICES, &items[0]);
	if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
	wv_adjustVerts(&items[0], focus);
	stat = wv_setData(WV_INT32, 3*ntri, (void *) tris, WV_INDICES, &items[1]);
	if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
	stat = wv_setData(WV_REAL32, 1, (void *) colorBox,  WV_COLORS, &items[2]);
	if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
	if (quad == 0) {
	    for (nseg = j = 0; j < ntri; j++) {
		for (k = 0; k < 3; k++) {
		    if (tric[3*j+k] < j+1) {
			nseg++;
		    }
		}
	    }
	} else {
	    for (nseg = j = 0; j < ntri/2; j++) {
		for (k = 0; k < 4; k++) {
		    if (tric[6*j+neigq[k]] < 2*j+1) {
			nseg++;
		    }
		}
	    }
	}
	segs = (int *) malloc(2*nseg*sizeof(int));
	if (segs == NULL) {
	    printf(" Can not allocate %d Sides!\n", nseg);
	    continue;
	}
	if (quad == 0) {
	    for (nseg = j = 0; j < ntri; j++)
	      for (k = 0; k < 3; k++)
		if (tric[3*j+k] < j+1) {
		    segs[2*nseg  ] = tris[3*j+sides[k][0]];
		    segs[2*nseg+1] = tris[3*j+sides[k][1]];
		    nseg++;
		}
	} else {
	    for (nseg = j = 0; j < ntri/2; j++)
	      for (k = 0; k < 4; k++)
		if (tric[6*j+neigq[k]] < 2*j+1) {
		    segs[2*nseg  ] = tris[6*j+sideq[k][0]];
		    segs[2*nseg+1] = tris[6*j+sideq[k][1]];
		    nseg++;
		}
	}
	stat = wv_setData(WV_INT32, 2*nseg, (void *) segs, WV_LINDICES, &items[3]);
	if (stat < 0) printf(" wv_setData = %d for %s/item 3!\n", i, gpname);
	stat = wv_setData(WV_REAL32, 1, (void *) colorLine,  WV_LCOLOR, &items[4]);
	if (stat < 0) printf(" wv_setData = %d for %s/item 4!\n", i, gpname);
	stat = wv_addGPrim(cntxt, gpname, WV_TRIANGLE,
			   WV_ON|WV_ORIENTATION, 5, items);
	if (stat < 0)
	  printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
	if (stat > 0) ngp = stat+1;
	sum += ntri;
	snprintf(fn, 100, "fqFace_%d", i + 1 );
	fil = fopen(fn, "w");
	if ( fil == NULL ) continue;
	for (k = 0; k < nseg; k++) {
	    fprintf(fil, "%lf %lf %lf\n", xyzs[ 3 * ( segs[2 * k] - 1 )],
		    xyzs[ 3 * ( segs[2 * k] - 1 ) + 1], xyzs[ 3 * ( segs[2 * k] - 1 ) + 2]);
	    fprintf(fil, "%lf %lf %lf\n\n\n", xyzs[ 3 * ( segs[2 * k + 1] - 1 )],
		    xyzs[ 3 * ( segs[2 * k + 1] - 1 ) + 1], xyzs[ 3 * ( segs[2 * k + 1] - 1 ) + 2]);
	}
	free(segs);
	fclose (fil);
    }
    // get edges //
    for (i = 0; i < bodydata[ibody].nedges; i++) {
      stat  = EG_getTessEdge(bodydata[ibody].tess, i+1, &len, &xyzs, &ts);
      if (stat != EGADS_SUCCESS) continue;
      if (len == 0) continue;
      nseg = len-1;
      segs = (int *) malloc(2*nseg*sizeof(int));
      if (segs == NULL) {
        printf(" Can not allocate %d segments for Body %d Edge %d\n",
            nseg, ibody, i+1);
        continue;
      }
      for (j = 0; j < len-1; j++) {
        segs[2*j  ] = j + 1;
        segs[2*j+1] = j + 2;
      }

      sprintf(gpname, "Body %d Edge %d", ibody+1, i+1);
      stat = wv_setData(WV_REAL64, len, (void *) xyzs, WV_VERTICES, &items[0]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
      wv_adjustVerts(&items[0], focus);
      stat = wv_setData(WV_REAL32, 1, (void *) colorLine,  WV_COLORS,   &items[1]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
      stat = wv_setData(WV_INT32, 2*nseg, (void *) segs, WV_INDICES,  &items[2]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
      free(segs);
      stat = wv_addGPrim(cntxt, gpname, WV_LINE, WV_ON, 3, items);
      if (stat < 0) {
        printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
      } else {
        if (cntxt != NULL)
          if (cntxt->gPrims != NULL) {
            cntxt->gPrims[stat].lWidth = 1.5;
            if (wv_addArrowHeads(cntxt, stat, 0.05, 1, &nseg) != 0)
              printf(" wv_addArrowHeads Error\n");
          }
        ngp = stat+1;
      }
    }
  }

  printf(" ** %d gPrims with %d triangles **\n", ngp, sum);
  // start the server code //

  stat = 0;
  if (wv_startServer(7681, NULL, NULL, NULL, 0, cntxt) == 0) {

    // we have a single valid server -- stay alive a long as we have a client //
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
  for (ibody = 0; ibody < nbody; ibody++) {
    EG_deleteObject(bodydata[ibody].tess);
    EG_free(bodydata[ibody].edges);
    EG_free(bodydata[ibody].faces);
    EG_free(bodydata[ibody].uvs);
  }
  free(bodydata);

  printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
  printf(" EG_close          = %d\n", EG_close(context));


  return 0;
}




/* call-back invoked when a message arrives from the browser */

void browserMessage(/*@unused@*/ void *wsi, char *text, /*@unused@*/ int lena)
{
}
