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
#include <time.h>
#ifdef  WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <winsock2.h>
#endif
#include "wsserver.h"
#include "egads.h"

#define MAX(a,b)    (((a) > (b)) ? (a) : (b))
#define MIN(a,b)    (((a) < (b)) ? (a) : (b))


#define CROSS(a,b,c)      c[0] = (a[1]*b[2]) - (a[2]*b[1]);\
    c[1] = (a[2]*b[0]) - (a[0]*b[2]);\
    c[2] = (a[0]*b[1]) - (a[1]*b[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

#define DEG160 2.79252680319093
#define DEG10  0.1745329251994329

#define PI     3.1415926535897931159979635
static int MESHPLOTS = 0;

#define EPS11   1.E-11
#define DEBUG

#define SWAP     0
#define COLLAPSE 1
#define SPLIT    2

#define SWAPCOLLAPSE        0
#define DOUBLECOLLAPSE      1
#define SWAPDOUBLECOLLAPSE  2
#define SWAPSPLIT           3
#define DOUBLESPLIT         4
#define SWAPDOUBLESPLIT     5
#define DOUBLESWAP          6
static int FACECHOICE = -1;

typedef struct{
  int verts[4], qadj[4], id ;
} Quad ;


typedef struct{
  int  *verts, *quads;
  int   nV, nQ; // nV = origin(1) + peaks (n)
  int  *idxV, *idxQ;
} vStar;

typedef struct {
  int     *qIdx, *qAdj, **valence, *vType, *remQuads, *remVerts;
  int     totQuads, totVerts, sizeVerts, sizeQuads ;
  double  *xyzs, *uvs ;
} meshData;


typedef struct {
  int      fID, oriQ, oriV, minQ, extraQuads, *vFix;
  ego      face;
  double   minsize, range[4], minAngle, maxAngle ;
  meshData *mesh, *bestMesh;
} meshMap;


typedef struct {
  ego      *faces, *edges, body, tess;
  int       mtype, nfaces, nedges, plen;
  meshMap  **qm;
} bodyData;

typedef struct{
  int  verts[6], vals[6];
  int  q[2];
} quadGroup  ;

/* globals used in these functions */

/**************************************/
static void swapInt ( int *a, int *b ) {
  int c;
  c  = *a;
  *a = *b;
  *b =  c;
}
static void unitVector(double *v ) {
  double n;
  n     = sqrt (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] );
  if( n > EPS11 ) {
      v[0] /=n; v[1] /=n; v[2] /=n;
  }
  else {
      v[0] = 0.0; v[1] = 0.0; v[2] = 0.0;
  }
  return;
}


static void
cross_product(double *A, double *B, double *cross) {
  cross[0] = A[1] * B[2] - A[2] * B[1];
  cross[1] = A[2] * B[0] - A[0] * B[2];
  cross[2] = A[0] * B[1] - A[1] * B[0];
}




static int
inList ( int n, int *list, int p) {
  int i ;
  for ( i = 0 ; i < n ; i++ )
    if ( list [i] == p ) return i;
  return -1;
}



/* IO FUNCTIONS */
//#ifdef DEBUG
static void printVertexCoords (meshMap *qm, int v ) {
  v--;
  printf(" #vertex %d  ==============================================\n" , v+1);
  printf("%lf %lf %lf %lf %lf %d\n",
	 qm -> mesh -> xyzs [3 * v ],qm -> mesh -> xyzs [3 * v +1], qm -> mesh -> xyzs [3 * v +2 ],
	 qm -> mesh -> uvs [2 * v + 1],  qm -> mesh -> uvs [2 * v ], v + 1);
  printf(" #==============================================\n" );
}
static void printQuadCoords (meshMap *qm, int qID ) {
  int i = 0, v;
  printf(" #=================== quad %d ==================\n", qID );
  for ( i = 0 ; i <= 4; i++ ) {
      v = qm -> mesh -> qIdx [ 4 * ( qID - 1) +  i%4] - 1;
      printf("%lf %lf %lf %d %lf %lf \n",
	     qm -> mesh -> xyzs [3 * v ],qm -> mesh -> xyzs [3 * v +1], qm -> mesh -> xyzs [3 * v + 2 ], v + 1,
	     qm -> mesh -> uvs [2 * v + 1],  qm -> mesh -> uvs [2 * v ]);
  }
  printf(" #==============================================\n" );
}
//#endif

static void   printMeshStats           (meshMap *qm, int sweep);
static void   printMesh                (meshMap *qm, char *);
static int    checkMesh                (meshMap *qm) ;
static int  EG_wvsData                 (meshData *qm, char *buffer);
static void printQuadSpecs             (meshData *qm, int id) ;
static void printQuadGroup             (meshData *qm, quadGroup qg ) ;
//static void printVertSpecs             (meshData *qm, int id);

/* MEMO FUNCTIONS */
static void   EG_freeStar              (vStar **star );
static void   EG_freeMeshData          (meshData **mesh );
static int    EG_allocMeshData         (meshData **mesh, int nQ, int nV );
static int    EG_restoreMeshMap        (meshMap *qm, meshData *m1, meshData *m2 ) ;
static int    EG_createMeshMap         (bodyData *bodydata);

/* GEOMETRY FUNCTIONS */
static int    EG_projectToTangentPlane (double *normal, double *O, double *p, double *proj) ;
static double EG_angleAtVnormalPlane   (meshMap *qm, int vC, int v1, int v2 );
static int    EG_angleAtBoundaryVertex (meshMap *qm, int v, int *links, double *size ) ;
static int    EG_quadAngleOrientation  (meshMap *qm, double minT, double maxT,  int qID,  int *ori, int *order, double *theta );
static int    EG_makePositiveAngles    (meshMap *qm, double minAngle, double maxAngle, int qID);
static void   EG_3Drotation            (meshMap *qm, double minT, double maxT, double dtheta, int iA, int iB, int iC);
static void   rotationUV               (meshMap *qm, double theta, int counterClockWise, double *t0, double *t1);
static double EG_segment               (ego face,    double *uv1, double *uv2 );
static int    EG_moveDistAway          (meshMap *qm, double alpha, int iA, int iB);
static int    EG_makeValidMesh         (meshMap *qm, int nP, /*@unused@*/ /*@null@*/int *pList, int fullRegularization);
static void   EG_averageCoords         (meshMap *qm, int vC, double min, double max );
/***********************************************************/
/* MESH MAP FUNCTIONS*/
static int  EG_cleanNeighborhood      (meshMap *qm, int qID, int transfer, int *activity );
static int  EG_cleanQuad              (meshMap *qm, int qID,  int adj, int tansfer, int forcing, int *activity ) ;
static int  setValence                (meshData *qm, int vID );
static int  getValence                (meshData *qm, int v );
static int  EG_buildStar              (meshData *qm, vStar **star, int vID );
static int    vertexLinksToBounds     (meshMap *qm, int vID, int *nb);
static void   updateVertex            (meshMap *qm, int vID, double *uv );
static int  EG_swappingOperation      (meshMap *qm, quadGroup qg, int swap, int *activity);
static int  EG_splittingOperation     (meshMap *qm, int vC, int vL, int vR, int *activity);
static int  EG_mergeVertices          (meshMap *qm, int qC, int centre, int *activity);
static int  EG_nValenceCount          (meshData *qm, int q, int n );
static int  checkQuad                 (meshData *qm, int q );
static int  checkVertex               (meshData *qm, int v );
static int  EG_quadIsBoundary         (meshData *qm, int qID );
static int  EG_adjQtoPair             (meshData *qm, int qID, int v1, int v2, int *adj) ;
static void EG_commonVerts            (meshData *qm, int q1, int q2, int *v );
static int  EG_createQuadGroup        (meshData *qm, quadGroup *qg, int q0, int q1 ) ;
/**************/
static int  EG_swap                   (meshMap *qm, int qID, int *activity);
static int  EG_doubleSwap             (meshMap *qm, quadGroup qg, int forcing, int *activity );
static int  EG_collapse               (meshMap *qm, int  vID, int *activity, int forcing);
static int  EG_split                  (meshMap *qm, int  qID, int *activity);
static int  EG_swapSplit              (meshMap *qm, quadGroup qg, int forcing, int *activity );
static int  EG_swapCollapse           (meshMap *qm, quadGroup qg, int forcing, int *activity );
static int  EG_doubleCollapse         (meshMap *qm, quadGroup qg, int forcing, int *activity );
static int  EG_doubleSplit            (meshMap *qm, quadGroup qg, int forcing, int *activity );
static int  EG_swapDoubleSplit        (meshMap *qm, quadGroup qg, int *activity );
static int  EG_swapDoubleCollapse     (meshMap *qm, quadGroup qg, int *activity );
static int  EG_transferValences       (meshMap *qm, int *qID, int try5533, int *transfering, int *activity  );
/* MESH DATA FUNCTIONS */
static int  validSwap                 (meshData *qm, int v1, int v2 );
static int  validCollapse             (meshMap *qm, int qID, int v );
/*********************/

static void EG_freeStar ( vStar **star ) {
  if ( (*star) == NULL ) return;
  EG_free ( (*star) -> verts);
  EG_free ( (*star) -> quads);
  EG_free ( (*star) -> idxQ);
  EG_free ( (*star) -> idxV);
  EG_free ( *star );
  *star = NULL;
}
static void EG_freeMeshData (meshData **mesh ) {
  int j;
  if ( *mesh  == NULL ) return ;
  for ( j = 0 ; j < (*mesh) -> sizeVerts; ++j) EG_free(( *mesh) -> valence[j]) ;
  EG_free(( *mesh) -> valence) ;
  EG_free(( *mesh) -> qIdx );
  EG_free(( *mesh) -> qAdj );
  EG_free(( *mesh) -> uvs     );
  EG_free(( *mesh) -> xyzs    );
  EG_free(( *mesh) -> remQuads);
  EG_free(( *mesh) -> remVerts);
  EG_free(( *mesh) -> vType   );
  EG_free(( *mesh)            );
  *mesh = NULL;
}
static int  EG_allocMeshData (meshData **mesh, int nQ, int nV ) {
  int j;
  if ( *mesh ) EG_freeMeshData( mesh );
  (*mesh)               = (meshData*) EG_alloc (sizeof (meshData) );
  if ( (*mesh) == NULL ) return EGADS_MALLOC;
  (*mesh) -> sizeVerts  = nV;
  (*mesh) -> sizeQuads  = nQ;
  (*mesh) -> xyzs       = (double *) EG_alloc(3*nV *sizeof(double ));
  (*mesh) -> uvs        = (double *) EG_alloc(2*nV *sizeof(double ));
  (*mesh) -> vType      = (int    *) EG_alloc(  nV *sizeof(   int ));
  (*mesh) -> qIdx    = (int    *) EG_alloc(4*nQ *sizeof(   int ));
  (*mesh) -> qAdj    = (int    *) EG_alloc(4*nQ *sizeof(   int ));
  (*mesh) -> remQuads   = (int    *) EG_alloc(  nQ *sizeof(   int ));
  (*mesh) -> remVerts   = (int    *) EG_alloc(  nV *sizeof(   int ));
  (*mesh) -> valence    = (int   **) EG_alloc(  nV *sizeof(   int*));
  if (    (*mesh)->qIdx  == NULL || (*mesh)->qAdj  == NULL ||
      (*mesh)->xyzs     == NULL || (*mesh)->uvs      == NULL ||
      (*mesh)->vType    == NULL || (*mesh)->remQuads == NULL ||
      (*mesh)->remVerts == NULL || (*mesh)->valence  == NULL )
    return EGADS_MALLOC;
  (*mesh) -> totVerts    = nV;
  (*mesh) -> totQuads    = nQ;
  (*mesh) -> remQuads[0] = 0;
  (*mesh) -> remVerts[0] = 0;
  for (j = 0; j < nV; j++) {
      (*mesh) -> valence[j] = (int    *) EG_alloc ((2 + nV )*sizeof(int   ));
      if ((*mesh) -> valence[j] == NULL )
	return EGADS_MALLOC;
  }
  return EGADS_SUCCESS;
}

static int EG_backupQuads (meshMap *qm, int *nq, int *qlist, Quad **quad) {
  int i, j, q, v, qcount, *qaux = NULL, k, add = 0, stat ;
  vStar *star = NULL;
  qaux = EG_alloc ( qm -> mesh ->  totQuads * sizeof ( int ));
  if ( qaux == NULL ) return EGADS_MALLOC;
  for ( qcount = q = 0 ; q < *nq; q++ ) {
      if ( qlist[q] == -1 ) continue;
      for ( i = 0 ; i < 4; i++) {
	  v    = qm -> mesh -> qIdx [ 4 * ( qlist[q] - 1) + i];
	  stat = EG_buildStar(qm -> mesh, &star, v);
	  if ( star == NULL || stat != EGADS_SUCCESS ) return stat;
	  for ( j = 0 ; j < star -> nQ; j++ ) {
	      if ( star-> quads[j] == -1 ) continue;
	      for ( add = k = 0 ; k < qcount; k++ ) {
		  if ( qaux[k] == star -> quads[j] ) {
		      add = 1;
		      break;
		  }
	      }
	      if ( add == 0 || qcount == 0 ) qaux[qcount++] = star ->quads[j];
	  }
      }
  }
  EG_freeStar (&star);
  (*quad)       = EG_alloc ( qcount * sizeof (Quad ));
  if ( (*quad) == NULL ) return EGADS_MALLOC;
  for ( q = 0 ; q < qcount; q++ ) {
      (*quad)[q].id = qaux[q];
      for ( j = 0 ; j < 4; j++ ) {
	  (*quad)[q].qadj [j] = qm -> mesh -> qAdj [ 4 * ( qaux[q] - 1) + j ];
	  (*quad)[q].verts[j] = qm -> mesh -> qIdx [ 4 * ( qaux[q] - 1) + j ];
      }
  }
  EG_free (qaux);
  *nq = qcount;
  return EGADS_SUCCESS;
}



static int EG_restoreQuads (meshMap *qm, Quad *quad, int nq ) {
  int i, j, *vid = NULL, k, stat;
  vid            = EG_alloc ( qm -> mesh -> totVerts * sizeof ( int) );
  if ( vid == NULL ) return EGADS_MALLOC;
  for ( k = i = 0 ; i < nq; i++ ) {
      if ( quad[i].id == -1 ) continue;
      for ( j = 0 ; j < 4; j++ ) {
	  qm -> mesh -> qAdj   [ 4 * ( quad[i].id - 1) + j] = quad[i].qadj [j];
	  qm -> mesh -> qIdx   [ 4 * ( quad[i].id - 1) + j] = quad[i].verts[j];
	  if ( k == qm -> mesh -> totVerts ) continue;
	  if ( k == 0 || ( k > 0 && inList (k, vid, quad[i].verts[j] ) == -1 )) {
	      qm -> mesh -> valence[ quad[i].verts[j] - 1][1]   = quad[i].id;
	      vid[k++] = quad[i].verts[j];
	  }
      }
  }
  for ( i  = 0 ; i < k; i++ ) {
      stat = setValence ( qm -> mesh, vid[i]);
      if ( stat != EGADS_SUCCESS ) {
	  EG_free ( vid );
	  return stat;
      }
  }
  EG_free ( vid );
  return EGADS_SUCCESS;
}

static int  EG_restoreMeshMap (meshMap *qm, meshData *m1, meshData *m2 ) {
  int i, j, k;
  m1 -> totQuads = m2 -> totQuads;
  m1 -> totVerts = m2 -> totVerts;
  for ( j = 0 ; j < m2 -> totQuads; ++j) {
      for ( k = 0 ; k < 4; ++k) {
	  m1 -> qIdx[4*j + k] = m2 -> qIdx[4*j+ k];
	  m1 -> qAdj[4*j + k] = m2 -> qAdj[4*j+ k];
      }
  }
  for ( j = 0 ; j < m2 -> totVerts; ++j) {
      m1 -> vType[  j    ] = m2 -> vType[  j    ];
      m1 -> uvs  [2*j    ] = m2 -> uvs  [2*j    ];
      m1 -> uvs  [2*j + 1] = m2 -> uvs  [2*j + 1];
      m1 -> xyzs [3*j    ] = m2 -> xyzs [3*j    ];
      m1 -> xyzs [3*j + 1] = m2 -> xyzs [3*j + 1];
      m1 -> xyzs [3*j + 2] = m2 -> xyzs [3*j + 2];
      if (m1 -> vType[j] == -2 ) {
	  m1 -> valence[j][0] = -1;
      } else {
	  for ( i = 0 ; i < m2->valence[j][0] + 2; ++i)
	    m1 -> valence[j][i] = m2 -> valence[j][i];
      }
  }
  for ( i = 0 ; i <= m2 -> remQuads[0]; i++ ) {
      m1 -> remVerts[i] = m2 -> remVerts[i];
      m1 -> remQuads[i] = m2 -> remQuads[i];
  }
#ifdef DEBUG
  return checkMesh (qm );
#endif
  return EGADS_SUCCESS;
}
static int EG_createMeshMap(bodyData *bodydata)
{
  int         f, iA, iB, stat = 0, stat2 = 0, j, q, i, auxID, k, kk, kOK, len,  ntri, nquad, e4[4], *faceEdges = NULL;
  const int   *tris, *tric, *ptype, *pindex;
  double      angle ;
  const double *xyzs, *uvs;
  char buffer[100];
  int    qV[6]    = { 0, 1, 2, 5, 0, 1};
  int    qLoop[5] = { 0, 1, 2, 3, 0   };
  vStar *star     = NULL;
  FILE *fill = NULL;
  bodydata->qm = (meshMap**) EG_alloc (bodydata->nfaces * sizeof(meshMap*));
  if (bodydata->qm == NULL ) return  EGADS_MALLOC;
  faceEdges = (int *)EG_alloc ( bodydata -> nedges * sizeof ( int));
  if ( faceEdges == NULL ) return EGADS_MALLOC;
  for ( f = 0 ; f < bodydata->nfaces; ++f) {
      bodydata ->qm[f] = (meshMap*) EG_alloc(sizeof(meshMap));
      if ( bodydata ->qm[f] == NULL ) {
	  printf("Create Quad Map: MALLOC error!! \n ");
	  EG_free (faceEdges );
	  return EGADS_MALLOC;
      }
      bodydata ->qm[f] ->mesh       = NULL;
      bodydata ->qm[f] ->bestMesh   = NULL;
      bodydata ->qm[f] ->fID        = f + 1;
      // Edges associated to face //
      stat = EG_getTessFace (bodydata->tess, f + 1, &len,
			     &xyzs, &uvs, &ptype, &pindex, &ntri,
			     &tris, &tric);
      if ( stat != EGADS_SUCCESS ) {
	  printf("EG_createMeshMap :: EG_getTessFace %d !!\n", stat );
	  stat2 = stat;
      }
      for ( i = 0 ; i < bodydata -> nedges; i++ ) faceEdges[i] = 0;
      for ( i = 0 ; i < len; i++ ) {
	  if ( pindex[i] == -1 ) continue;
	  faceEdges[pindex[i] -1]++;
      }
      nquad = (int)ntri/2;
      stat  = EG_allocMeshData ( &bodydata ->qm[f] ->mesh      ,2 * nquad,  2 * len);
      stat += EG_allocMeshData ( &bodydata ->qm[f] ->bestMesh  ,2 * nquad,  2 * len);
      bodydata ->qm[f] -> vFix = EG_alloc ( 2 * len * sizeof ( int ));
      if ( stat != EGADS_SUCCESS ||
	  bodydata ->qm[f] ->mesh       == NULL ||
	  bodydata ->qm[f] ->bestMesh   == NULL ||
	  bodydata ->qm[f] ->vFix       == NULL  ) {
	  fprintf(stderr,"In createMeshMap  EG_allocMeshData = %d\n ", stat );
	  stat2 = stat;
	  continue;
      }
      bodydata->qm[f] ->vFix[0]    = 0;
      bodydata->qm[f] ->extraQuads = 0;
      bodydata->qm[f] ->face       = bodydata->faces[f];
      stat                         = EG_getRange ( bodydata ->qm[f] ->face, bodydata ->qm[f] ->range, &i);
      if ( stat != EGADS_SUCCESS ) {
	  printf(" Face %d EG_getRange %d !!\n ", f + 1, stat );
	  stat2  = stat ;
	  continue;
      }
      bodydata->qm[f] ->mesh ->totVerts = len;
      bodydata->qm[f] ->oriV            = len;
      bodydata->qm[f] ->mesh ->totQuads = nquad;
      bodydata->qm[f] ->oriQ            = nquad;
      for (j = 0; j < len; j++) {
	  bodydata->qm[f] ->mesh ->valence[j][0]        = 0;
	  bodydata->qm[f] ->mesh ->vType  [j]           = ptype [j];
	  bodydata->qm[f] ->mesh ->uvs    [ 2 * j     ] = uvs[ 2 * j    ];
	  bodydata->qm[f] ->mesh ->uvs    [ 2 * j + 1 ] = uvs[ 2 * j + 1];
	  updateVertex (bodydata ->qm[f], j + 1, &bodydata->qm[f] ->mesh ->uvs[ 2 * j ] );
      }
      for (j = 0; j < nquad; j++)
	for ( k = 0; k < 4; ++k) {
	    bodydata ->qm[f] ->mesh ->qIdx[4*j + k ] = tris[6*j + qV[k+1]];
	}
      for ( j = 0; j < nquad; j++) {
	  kk  = 0;
	  kOK = 0;
	  q   = 0;
	  while ( q < nquad )
	    {
	      if ( q == j ) {
		  if ( j == nquad - 1 && kk < 4 ) {
		      bodydata ->qm[f] ->mesh ->qAdj[4*j + kk++] = -1;
		      if (kk == 4) break;
		      q  = 0;
		  }
		  else q++;
	      }
	      if  ( q == nquad ) break;
	      for ( k = 0 ; k < 4; ++k ) {
		  if(( bodydata ->qm[f] ->mesh ->qIdx[4*j + qLoop[kk]]     == bodydata ->qm[f] ->mesh ->qIdx[4*q + qLoop[k    ]]  ||
		      bodydata ->qm[f] ->mesh ->qIdx[4*j + qLoop[kk]]     == bodydata ->qm[f] ->mesh ->qIdx[4*q + qLoop[k + 1]]) &&
		      (bodydata ->qm[f] ->mesh ->qIdx[4*j + qLoop[kk + 1]] == bodydata ->qm[f] ->mesh ->qIdx[4*q + qLoop[k    ]] ||
			  bodydata ->qm[f] ->mesh ->qIdx[4*j + qLoop[kk + 1]] == bodydata ->qm[f] ->mesh ->qIdx[4*q + qLoop[k + 1]]) )
		    {
		      bodydata->qm[f] ->mesh ->qAdj[4*j + kk++] = q + 1;
		      q   =  -1;
		      kOK =   1;
		      k   =   4;
		      if (kk == 4)  q = nquad;
		    }
	      }
	      if ( (kOK == 0) && (q >= nquad -1) ){
		  bodydata ->qm[f] ->mesh  ->qAdj[4*j + kk++] = -1;
		  q                       = -1;
		  if (kk == 4) break;
	      }
	      else  kOK = 0 ;
	      q++;
	    }
	  if (kOK == 0 && kk < 4) {
	      while (kk < 4)
		bodydata ->qm[f] ->mesh ->qAdj[4*j + kk++] = -1;
	  }
      }
      for ( j = 0 ; j < nquad; ++j) {
	  for ( k = 0 ; k < 4; ++k) {
	      if ( bodydata ->qm[f] ->mesh ->qAdj[4*j +k] > j + 1  || bodydata ->qm[f] ->mesh ->qAdj[4*j +k] == -1) {
		  auxID = bodydata ->qm[f] ->mesh ->qIdx [4*j + qLoop[k]] -1;
		  bodydata ->qm[f] ->mesh ->valence[auxID][1] = j + 1;
		  bodydata ->qm[f] ->mesh ->valence[auxID][2 + bodydata ->qm[f] ->mesh ->valence[auxID][0]++] = bodydata ->qm[f] ->mesh ->qIdx[4*j+ qLoop[k+1]];
		  auxID = bodydata->qm[f] ->mesh  -> qIdx [4*j + qLoop[k+1]] -1;
		  bodydata ->qm[f] ->mesh ->valence[auxID][1] = j + 1;
		  bodydata ->qm[f] ->mesh ->valence[auxID][2 + bodydata ->qm[f] ->mesh ->valence[auxID][0]++] = bodydata->qm[f] -> mesh  -> qIdx[4*j+ qLoop[k]];
	      }
	  }
      }
      sprintf(buffer,"gnuInit_%i",f+1);
      printMesh(bodydata->qm[f] , buffer);
      sprintf(buffer,"wvsInit_%i.txt",f+1);
      EG_wvsData(bodydata->qm[f] -> mesh , buffer);
      sprintf(buffer, "face_%d_edges", f + 1);
      fill = fopen ( buffer, "w");
      if ( fill == NULL ) {
	  EG_free(faceEdges);
	  return EGADS_MALLOC;
      }
      bodydata ->qm[f] ->minsize  = 1000000000.0;
      bodydata ->qm[f] ->minAngle = DEG10 * 0.1;
      bodydata ->qm[f] ->maxAngle =          PI;
      for ( j = 0 ; j < len; j++ ) {
	  if ( bodydata ->qm[f] ->mesh ->vType[j] == -1 ) continue;
	  k = bodydata ->qm[f] ->mesh ->valence[j][0];
	  for ( i   = 0 ; i < k; i++ ) {
	      auxID = bodydata ->qm[f] ->mesh ->valence[j][2 + i] - 1;
	      if ( bodydata ->qm[f] ->mesh ->vType[auxID] == -1 ) continue;
	      angle = EG_segment (bodydata ->qm[f] ->face, &bodydata ->qm[f] ->mesh ->uvs [ 2 *j ], &bodydata ->qm[f] ->mesh ->uvs [ 2 *auxID ]);
	      bodydata ->qm[f] ->minsize = MIN (angle, bodydata ->qm[f] ->minsize);
	  }
	  fprintf(fill, "%lf %lf %lf %d\n", xyzs[3 * j ], xyzs[3 * j +1], xyzs[3 * j +2], j + 1);
	  stat = EG_angleAtBoundaryVertex (bodydata ->qm[f], j + 1, e4, &angle);
	  if ( stat != EGADS_SUCCESS || angle < EPS11 ) {
	      stat2 = EGADS_GEOMERR;
	      printf(" Stat in EG_angleAtBoundaryVertex %d angle %f\n ", stat, angle );
	      printf(" Vertices: %d %d %d \n ", j+ 1, e4[0], e4[1] );
	      break;
	  }
	  else if ( angle < 0.85 * PI ) {
	      if ( bodydata ->qm[f] ->mesh ->vType[j] > 0 ) {
		  printf(" Node %d is now consider edge so that its regular valence is 3 \n", j + 1);
		  bodydata ->qm[f] ->mesh ->vType[j] = 0;
	      }
	  }
	  else if ( angle < 1.25 * PI  ) {
	      if ( bodydata ->qm[f] ->mesh ->valence[j][0] >= 3 )
		bodydata ->qm[f] ->mesh ->vType  [j] = 3;
	  }
	  else if ( angle < 1.85 * PI ){
	      if ( bodydata ->qm[f] ->mesh ->valence[j][0] >= 4 )
		bodydata ->qm[f] ->mesh ->vType  [j] = 4 ;
	  } else {
	      if ( bodydata->qm[f] ->mesh ->valence[j][0] >= 5 )
		bodydata->qm[f] ->mesh ->vType  [j] = 5 ;
	  }
	  stat = EG_buildStar ( bodydata->qm[f] ->mesh, &star, j + 1 );
	  if ( stat != EGADS_SUCCESS || star == NULL ) {
	      stat2 = stat;
	      continue;
	  }
	  for ( k = 0; k < star -> nQ; k++ ) {
	      if ( star -> quads[k] == -1 ) continue;
	      iA = star -> verts[  2 * k + 1] ;
	      iB = star -> verts[ star -> idxV [ 2 * k + 3] ] ;
	      angle = EG_angleAtVnormalPlane (bodydata -> qm[f], j + 1, iA, iB);
	      bodydata ->qm[f] ->minAngle = MIN (angle, bodydata ->qm[f] ->minAngle);
	      bodydata ->qm[f] ->maxAngle = MAX (angle, bodydata ->qm[f] ->maxAngle);
	  }
	  EG_freeStar ( & star ) ;
      }
      fclose (fill);
      stat = checkMesh (bodydata->qm[f]);
      printf(" Min max angle at bounds %lf  %lf   ( 10 DEG %lf  160 DEG %lf ) \n ", bodydata ->qm[f] ->minAngle, bodydata ->qm[f] ->maxAngle, DEG10, DEG160 );
      if ( stat != EGADS_SUCCESS ) {
	  printf("In EG_createMeshMap :: checkMesh at face %d --> %d!!\n", f + 1, stat );
	  stat2 = stat;
      }
  }
  EG_free(faceEdges);
  return stat2;
}

static int EG_projectToTangentPlane(double normal[], double *nO, double *p, double *proj) {
  double c, dotNN = 0.0, dotNP = 0.0, dist, lambda;
  c       = DOT (normal, nO); // Equation plane: a*x + b*y + c*z = C
  dotNN   = DOT (normal, normal);
  dotNP   = DOT (normal, p);
  dist    = (p[0] - nO[0]) * (p[0] - nO[0]) + (p[1] - nO[1]) * (p[1] - nO[1]) + (p[2] - nO[2]) * (p[2] - nO[2]);
  if ( fabs(dotNP - c) <= EPS11 || sqrt (dist) < EPS11) {
      proj[0] = p[0]; proj[1] = p[1]; proj[2] = p[2];
      return EGADS_SUCCESS;
  }
  lambda  = (c - dotNP) / dotNN;
  proj[0] = p[0] + lambda * normal[0];
  proj[1] = p[1] + lambda * normal[1];
  proj[2] = p[2] + lambda * normal[2];
  dist    = DOT ( normal, proj);
  if( fabs(dist - c) < EPS11)
    return EGADS_SUCCESS;
  else{
      printf(" ORIGIN %lf %lf  %lf  NORMAL %lf %lf  %lf  TARGET %lf %lf %lf\n", nO[0], nO[1], nO[2], normal[0], normal[1], normal[2], p[0], p[1], p[2]);
      printf(" POINT SHOULD BELONG TO PLANE!!!!! %.16e ~= %.16e\n",dist,c);
      printf(" DOT NN %lf PN %lf LAMBDA %lf  POINT %lf %lf %lf\n", dotNN, dotNP, lambda, proj[0], proj[1], proj[2]);
      return EGADS_GEOMERR;
  }
}

static double EG_angleAtVnormalPlane ( meshMap *qm, int vC, int v1, int v2 ) {
  int i, stat;
  double dot1, dot2, xyz [18], normal[3], cross[3], proj1[3], proj2[3], centre[18], v01[3], v02[3];
  stat       = EG_evaluate(qm -> face, &qm -> mesh -> uvs [ 2 * ( vC - 1 ) ], centre);
  if ( stat != EGADS_SUCCESS) return stat;
  if ( qm -> face -> mtype == SREVERSE )
    cross_product(&centre[6], &centre[3], normal);
  else
    cross_product(&centre[3], &centre[6], normal);
  unitVector(normal );
  stat  = EG_evaluate(qm -> face, &qm -> mesh -> uvs [ 2 * ( v1 - 1 ) ], xyz);
  stat += EG_projectToTangentPlane(normal, centre, xyz, proj1);
  stat += EG_evaluate(qm -> face, &qm -> mesh -> uvs [ 2 * ( v2 - 1 ) ], xyz);
  stat += EG_projectToTangentPlane(normal, centre, xyz, proj2);
  if ( stat != EGADS_SUCCESS ) return stat;
  for ( i = 0 ; i < 3 ; i++ ) {
      v01[i] = proj1[i] - centre[i];
      v02[i] = proj2[i] - centre[i];
  }
  unitVector (v01);
  unitVector (v02);
  CROSS      (v01, v02, cross );
  dot1 = DOT (v01, v02 );
  dot2 = DOT (normal, cross);
  if      ( dot1 >=  1.0 ) return 0.0;
  else if ( dot1 <= -1.0 ) return PI;
  if ( dot2 > 0 ) return  acos ( dot1 );
  else  return ( 2.0*PI - acos ( dot1 ));
}

static int EG_angleAtBoundaryVertex(meshMap *qm, int v, int *links, double *size ) {
  int i, j, k;
  vStar *star = NULL;
  *size = 0.0;
  if ( qm -> mesh -> vType[v - 1 ] < 0 ) return EGADS_INDEXERR;
  i = EG_buildStar ( qm -> mesh, &star, v);
  if ( i != EGADS_SUCCESS || star == NULL ) {
      printf(" Looking at corners: buildstar %d is %d \n ", v, i );
      return i;
  }
  for ( links[0] = links[1] = k = i = 0 ;i < star -> nQ; i++ ) {
      j = star -> verts[ 2 * i + 1 ] - 1 ;
      if ( qm -> mesh -> vType[j] != -1 ) k++;
      if ( star -> quads[i] == -1  ) {
	  links[1] = star -> verts [ 2 * i + 1];
	  links[0] = star -> verts [star -> idxV[ 2 * i + 3]];
	  continue;
      }
  }
  EG_freeStar ( &star);
  if ( k >=3 ) *size = PI * 0.5; // boundary vertex is connected to more than two bounds. Angle is fine
  else if ( k != 2 ) {
      printf(" EG_angleAtBoundaryVertex:: vertex %d is at surface bounds and connected only to another boundary vertex !!\n ", v);
      return EGADS_GEOMERR;
  }
  *size = EG_angleAtVnormalPlane (qm, v , links[0], links[1] );
  return EGADS_SUCCESS;
}


/* Returns
 *  1 if all angles are between the range
 *  0 invalid quad || too big / small
 *  < 0 something terrible happened
 */

static int EG_quadAngleOrientation(meshMap *qm, double minT, double maxT,  int qID, int *ori, int *order, double *theta ) {
  int     i, qV[4], k, k1, k2, sign, stat, vA, vB, vC, area = 0, count = 0, nb[10];
  double cross[3], qNormal[3], vAB[3], vAC[3], dot, quv[2], xyz[18], projABCD[12], signedAngle ;
  quv[0] = 0.0; quv[1] = 0.0;
  for (i = 0; i < 4; i++ ) {
      qV[i]   = qm -> mesh ->qIdx[4*(qID - 1) + i] - 1;
      if ( vertexLinksToBounds (qm, qV[i] + 1 , nb) > 0 && maxT > PI  ) maxT = PI;
      quv[0] += 0.25 * qm -> mesh -> uvs [2 * qV[i]    ];
      quv[1] += 0.25 * qm -> mesh -> uvs [2 * qV[i] + 1];
  }
  i        = EG_evaluate(qm->face, quv, xyz);
  if  ( i != EGADS_SUCCESS ) {
      printf(" In EG_quadAngleOrientation:: EG_evaluate for quad %d  is %d !!\n", qID, i );
      return i;
  }
  unitVector (&xyz[3]);
  unitVector (&xyz[6]);
  if ( qm -> face -> mtype == SREVERSE )
    cross_product(&xyz[6], &xyz[3], qNormal);
  else
    cross_product(&xyz[3], &xyz[6], qNormal);
  for ( i = 0 ; i < 4; i++ ) {
      stat       = EG_projectToTangentPlane(qNormal, xyz, &qm -> mesh -> xyzs [ 3 * qV[i] ], &projABCD[3 * i]);
      if ( stat != EGADS_SUCCESS ) {
	  printf("EG_quadAngleOrientation :: EG_projectToTangentPlane quad %d vert %d --> %d !!\n", qID, qV[i], stat );
	  return stat;
      }
  }
  for ( sign = k = 0 ; k < 4; ++k) {
      vA   =   k        ;
      vB   = ( k + 1 )%4;
      vC   = ( k + 3 )%4;
      for ( i = 0 ; i < 3; ++i) {
	  vAB[i] = projABCD[3 * vB + i] - projABCD[3 * vA + i];
	  vAC[i] = projABCD[3 * vC + i] - projABCD[3 * vA + i] ;
      }
      unitVector (vAB);
      unitVector (vAC);
      CROSS      (vAB, vAC, cross);
      dot         =  DOT(vAB, vAC);
      if      ( fabs(dot - 1.0) < EPS11 ) theta[k] = 0.0;
      else if ( fabs(dot + 1.0) < EPS11 ) theta[k] = PI;
      else                                theta[k] = acos(dot);
      if ( DOT (qNormal, cross ) < 0 ) {
	  signedAngle = 2.0 * PI - theta[k];
	  ori[k]      = -1;
      }
      else {
	  ori[k] =  1;
	  signedAngle = theta[k];
      }
      order[k] = k;
      if ( ori[k] == -1 ) sign = 1;
      if ( qm -> mesh -> vType [qV[vA]] == -1 && qm -> mesh -> valence [ qV[vA] ][0] == -2 ) {
	  sign = 1;
	  count++;
      } else if ( signedAngle <= maxT && signedAngle >= minT ) count++;
      else {
#ifdef DEBUG
	  printf(" angle %f max %f min %f\n ", signedAngle, maxT, minT );
#endif
      }
  }
  if   ( count == 4 ) area = 1;
  for  ( k = 0; k < 3; k++) {
      for  ( i = k + 1; i < 4; i++) {
	  k1   = 0;
	  if ( sign == 0 ) {
	      if ( theta[ order [i] ] > theta[ order[k] ] ) k1 = 1;
	  } else {
	      if (theta[ order[i] ] * (double)ori[order[i]] <
		  theta[ order[k] ] * (double)ori[order[k]] )
		k1 = 1;
	  }
	  if ( k1 == 0 ) continue;
	  k2       = order[i];
	  order[i] = order[k];
	  order[k] = k2;
      }
  }
#ifdef DEBUG
  if ( area != 1 ) {
      printf(" ************   ATENCION:: INVALID QUAD  %lf  %lf  ***************\n", minT, maxT );
      printf("------------ AREA QUAD %d  IS %d Internal angles ordered by size --------------\n", qID, area );
      for ( sign = k = 0 ; k < 4; ++k)
	printf("Vertex %d has angle %f and orientation %d \n ", qV[order[k]] + 1, theta[order[k]], ori[order[k]]);
  }
#endif
  return area;
}

/* clockwise = -1 , counterclockwise = 1
 * reverse = 0 assumes invalid quads e.g., ( AB x AC ) * normal < 0 so rotated B towards C (clockwise ) and
 * C towards B ( counterclockwise).if <ABC is small, it will increase the angle.
 * reverse = 1 assumes <ABC is big and will move B counter-clockwise and C clockwise
 */
static void
EG_3Drotation ( meshMap *qm, double minT, double maxT, double theta, int iA, int iB, int iC ) {
  int  i, a, v,  area, counterClockWise = 1, stat, fix = 0, reverse = 0, piv[4], ori[4], links[2];
  double ang0, ang1, angles[4], uvf[2], angle;
  vStar *star = NULL;
  angle = EG_angleAtVnormalPlane( qm, iA, iB, iC ) ;
#ifdef DEBUG
  printf(" LOOKING AT ANGLE %d %d  %d  --> %lf   MIN MAX %lf  %lf  \n ", iA, iB, iC, angle, minT, maxT );
#endif
  if ( angle > maxT ) {
      if ( qm -> mesh -> vType[iA - 1] >= 0 ) {
	  i    = EG_angleAtBoundaryVertex (qm, iA, links,  &ang0);
	  ang0 = EG_angleAtVnormalPlane   (qm, iA, links[0], iB );
	  ang1 = EG_angleAtVnormalPlane   (qm, iA, links[0], iC );
	  if ( ang0 < ang1 || links[0] == iC ) reverse  = 1;
      } else if ( angle < 1.5 * PI ) reverse = 1;
  }
  else if ( angle >= minT ) return;
  ang0 = EG_segment (qm -> face , &qm -> mesh -> uvs [ 2 * ( iA - 1) ],
		     &qm -> mesh -> uvs [ 2 * ( iB - 1) ]);
  ang1 = EG_segment (qm -> face, &qm -> mesh -> uvs [ 2 * ( iA - 1) ],
		     &qm -> mesh -> uvs [ 2 * ( iC - 1) ]);
  if (      ang0 > ang1 && ang1 / ang0 < 0.95 ) EG_moveDistAway (qm, 0.95, iA, iB);
  else if ( ang0 < ang1 && ang0 / ang1 < 0.95 ) EG_moveDistAway (qm, 0.95, iA, iC);
  if ( reverse == 1 ) counterClockWise *= -1;
  for ( a = 0 ; a < 2; a++) {
      counterClockWise *= -1;
      if ( a == 0 ) v = iB - 1;
      else          v = iC - 1;
      fix = 0;
      if (qm -> mesh -> vType [v] != -1 ) continue;
      if (qm -> vFix[0] > 0 && inList (qm -> vFix[0], &qm -> vFix[1], v + 1) >= 0 ) {
	  stat    = EG_buildStar ( qm -> mesh, &star, v + 1 ) ;
	  if ( stat != EGADS_SUCCESS || star == NULL ) return ;
	  fix     = 1;
	  uvf[0]  = qm -> mesh -> uvs[ 2 * v    ];
	  uvf[1]  = qm -> mesh -> uvs[ 2 * v  +1];
	  for ( i = 0 ; i < star -> nQ; i++ ) {
	      star -> verts[i] = EG_quadAngleOrientation(qm, minT, maxT, star -> quads[i], ori, piv, angles);
#ifdef DEBUG
	      printf(" i %d -> %d AReA %d \n ", i, star -> quads[i], star -> verts[i] );
#endif
	  }
      }
#ifdef DEBUG
      printf(" ROTATE %d  %d   %d \n ", iA, v + 1, counterClockWise);
#endif
      rotationUV (qm, theta, counterClockWise, &qm -> mesh -> uvs[ 2 * ( iA -1 ) ],
		  &qm -> mesh -> uvs[ 2 * v ] );
      updateVertex (qm, v + 1, &qm -> mesh -> uvs[ 2 * v ]);
      if ( fix == 1 && star != NULL ) {
	  for ( i = 0 ; i < star -> nQ; i++ ) {
	      area = EG_quadAngleOrientation(qm, minT, maxT, star -> quads[i], ori, piv, angles);
	      if ( area == 0 && star -> verts[i] == 1 ) {
#ifdef DEBUG
		  printf (" ROTATED COOORDS MAKE INVALID QUAD \n ");
		  printVertexCoords (qm, v + 1 );
#endif
		  updateVertex ( qm, v + 1, uvf ) ;
		  break;
	      }
	  }
	  EG_freeStar ( & star ) ;
      }
  }
  return;
}


static void
rotationUV ( meshMap *qm, double theta, int counterClockWise, double *t0, double *t1) {
  int d, stat;
  double normal[4], proj[3], cross[4], xyz1[18], dt, duv[2], xyz0[18], prevt[2], v1[3], v2[3];
  stat      = EG_evaluate ( qm -> face, t0, xyz0);
  if ( stat != EGADS_SUCCESS ) {
#ifdef DEBUG
      printf("EG_rotationUV :: EG_normalToSurface ->%d!!\n",  stat);
#endif
      return;
  }
  unitVector (&xyz0[3]);
  unitVector (&xyz0[6]);
  if ( qm -> face -> mtype == SREVERSE )
    cross_product(&xyz0[6], &xyz0[3], normal);
  else
    cross_product(&xyz0[3], &xyz0[6], normal);
  prevt[0]  = t1[0]; prevt[1]  = t1[1];
  stat      = EG_evaluate ( qm -> face, t0, xyz0);
  stat     += EG_evaluate ( qm -> face, t1, xyz1);
  if ( stat != EGADS_SUCCESS ){
#ifdef DEBUG
      printf(" RotationUV rotation vector has EG_evaluate %d!!\n", stat);
#endif
      return;
  }
  stat   = EG_projectToTangentPlane(normal, xyz0, xyz1, proj);
  if ( stat != EGADS_SUCCESS ) {
#ifdef DEBUG
      printf("EG_rotateUV :: EG_projectToTangentPlane %d !!\n", stat);
#endif
      return ;
  }
  v1[0]  = proj[0] - xyz0[0];
  v1[1]  = proj[1] - xyz0[1];
  v1[2]  = proj[2] - xyz0[2];
  unitVector (v1);
  for ( d = 0 ; d < 2; d++ ) {
      dt     = theta;
      if ( d == 1 ) dt *= -1.0;
      duv[0] = ( prevt[0] - t0[0] ) * cos ( dt ) - ( prevt[1] - t0[1] ) * sin ( dt );
      duv[1] = ( prevt[0] - t0[0] ) * sin ( dt ) + ( prevt[1] - t0[1] ) * cos ( dt );
      t1[0]  = t0[0] + duv[0];
      t1[1]  = t0[1] + duv[1];
      if (    t1[0] < qm ->range[0] || t1[0] > qm ->range[1] ||
	  t1[1] < qm ->range[2] || t1[1] > qm ->range[3] ) {
	  if ( d == 0 ) continue;
	  t1[0] = prevt[0]; t1[1] = prevt[1];
	  return;
      }
      stat       = EG_evaluate ( qm -> face, t1, xyz1);
      if ( stat != EGADS_SUCCESS ) {
	  if ( d == 0 ) continue;
	  t1[0] = prevt[0]; t1[1] = prevt[1];
	  return;
      }
      stat = EG_projectToTangentPlane(normal, xyz0, xyz1, proj);
      if ( stat != EGADS_SUCCESS ) {
#ifdef DEBUG
	  printf("EG_rotateUV :: EG_projectToTangentPlane %d !!\n", stat);
#endif
	  t1[0] = prevt[0]; t1[1] = prevt[1];
	  return ;
      }
      v2[0] = proj[0] - xyz0[0];
      v2[1] = proj[1] - xyz0[1];
      v2[2] = proj[2] - xyz0[2];
      unitVector (v2 );
      CROSS (v1, v2, cross);
      if ( DOT (normal, cross ) * (double)counterClockWise > 0 ) return;
      else {
	  if ( d == 1 ) {
	      t1[0] = prevt[0]; t1[1] = prevt[1];
	      return;
	  }
      }
  }
  return;
}

static double EG_segment (ego face, double *uv1, double *uv2 ) {
  int i, n = 10, stat ;
  double p1[18], p2[18],  dist[2], uvEps[2], seg = 0.0, totArc = 0.0, dt;
  stat       = EG_evaluate (face, uv1, p1);
  stat      += EG_evaluate (face, uv2, p2);
  if ( stat != EGADS_SUCCESS ) return -1.0;
  dist[0]    = uv2[0] - uv1[0];
  dist[1]    = uv2[1] - uv1[1];
  for ( i = 0; i < n; i++ ) {
      dt         = (double) ( i + 1 ) / (double) n;
      uvEps[0]   = uv1[0] + dt * dist[0];
      uvEps[1]   = uv1[1] + dt * dist[1];
      stat       = EG_evaluate (face, uvEps, p2 );
      if ( stat != EGADS_SUCCESS ) return -1.0;
      seg        = ( p1[0] - p2[0] ) * ( p1[0] - p2[0] );
      seg       += ( p1[1] - p2[1] ) * ( p1[1] - p2[1] );
      seg       += ( p1[2] - p2[2] ) * ( p1[2] - p2[2] );
      p1[0]      = p2[0]; p1[1] = p2[1]; p1[2] = p2[2];
      totArc    += sqrt (seg);
  }
  return totArc;
}


static int EG_moveDistAway ( meshMap *qm, double alpha, int iA, int iB) {
  int    i, it;
  double seg, uvm[2], uvcopy[4], norm = 0.0, uv0[2], uv1[2], size, dt, segprev;
  if ( qm -> mesh -> vType[iB -1] != -1 ||
      (qm -> vFix[0] > 0 && inList (qm -> vFix[0], &qm -> vFix[1], iB) >= 0 )) return EGADS_SUCCESS;
  uv0[0]     = qm -> mesh -> uvs[ 2 * ( iA - 1 )    ];
  uv0[1]     = qm -> mesh -> uvs[ 2 * ( iA - 1 ) + 1];
  uv1[0]     = qm -> mesh -> uvs[ 2 * ( iB - 1 )    ];
  uv1[1]     = qm -> mesh -> uvs[ 2 * ( iB - 1 ) + 1];
  seg        = EG_segment (qm -> face, uv0, uv1 );
  if ( seg < 0.0 ) return EGADS_GEOMERR;
  size       = seg * alpha;
  if ( alpha * seg  < qm -> minsize * 0.25 ) size = 0.25 *  qm -> minsize;
#ifdef DEBUG
  printf(" Homotecia AB %d %d %lf --> %lf \n", iA, iB, seg, size );
#endif
  uvm   [0] = uv1[0]; uvm   [1] = uv1[1];
  uvcopy[0] = uv0[0]; uvcopy[1] = uv0[1];
  uvcopy[2] = uv1[0]; uvcopy[3] = uv1[1];
  it        = 0;
  dt        = 0.25;
  segprev   = seg;
  while ( seg < size ) {
      uvm[0] += ( uvcopy[2] - uvcopy[0] ) * dt;
      uvm[1] += ( uvcopy[3] - uvcopy[1] ) * dt;
      if (uvm[0] < qm -> range[0] || uvm[0] > qm -> range[1] ||
	  uvm[1] < qm -> range[2] || uvm[1] > qm -> range[3] ) return EGADS_SUCCESS;
      uv1[0] = uvm[0] ; uv1[1] = uvm[1];
      it++;
      seg    = EG_segment (qm -> face, uv0, uvm );
      if ( it > 100) {
#ifdef DEBUG
	  printf(" Seg < size Stuck. Break !\n");
#endif
	  break;
      }
      if ( fabs ( segprev - seg ) < size / 10.0 ) dt += 0.05;
      segprev = seg;
  }
  uvcopy[2] = uv1[0]; uvcopy[3] = uv1[1];
  it        = 0 ;
  while ( fabs( seg - size ) / size > 0.005 ) {
      uvm[0] = ( uvcopy[0] + uvcopy[2] ) * 0.5;
      uvm[1] = ( uvcopy[1] + uvcopy[3] ) * 0.5;
      seg    = EG_segment (qm -> face, uv0, uvm );
      i      = 0;
      if ( seg > size ) i = 1;
      uvcopy[ 2 * i    ] = uvm[0];
      uvcopy[ 2 * i + 1] = uvm[1];
      it++;
      norm  = ( uvcopy[0] - uvcopy[2] )* ( uvcopy[0] - uvcopy[2] );
      norm += ( uvcopy[1] - uvcopy[3] )* ( uvcopy[1] - uvcopy[3] );
      if ( sqrt (norm ) < EPS11 || it > 100 ) {
#ifdef DEBUG
	  printf(" POINTS ARE TOO CLOSE TO EACH OTHER: %f it %d  BISECTION DIDN't CONVERGE\n", sqrt (norm), it );
#endif
	  break;
      }
  }
  updateVertex (qm, iB, uvm);
  return EGADS_SUCCESS;
}

static int EG_makeValidMesh(meshMap *qm, int nP, /*@unused@*/ /*@null@*/int *pList, int fullRegularization)
{
  int  i, j, k, q, stat = EGADS_SUCCESS, it = 0, itMax,  ori[4], piv[4], *quads = NULL;
  double dthetam, dthetaM, minT, maxT, angles[4], *uvxyz = NULL;
#ifdef DEBUG
  char buffer[100];
#endif
  vStar *star   = NULL;
  qm -> vFix[0] = 0;
  if ( fullRegularization == 0 ) { // move around only affected vertices
      if ( nP == 0 || pList == NULL ) return EGADS_SUCCESS;
      quads     = EG_alloc ( qm -> mesh -> totQuads * sizeof ( int ) ) ;
      if (quads == NULL ) return EGADS_MALLOC;
#ifdef DEBUG
      printf(" Make valid mesh for %d points \n ", nP );
      printMesh (qm, NULL);
#endif
      for ( k = j = 0 ; j < nP; j++) {
	  if ( qm -> mesh -> vType [ pList[j] -1 ] == -2 ) continue;
	  printf(" BUILD STAR %d \n ", pList[j] );
	  stat  = EG_buildStar (qm -> mesh, &star, pList[j]);
	  if ( stat != EGADS_SUCCESS || star == NULL ) {
	      EG_free(quads );
	      return stat;
	  }
	  for ( i = 0 ; i < star -> nQ; i++ ) {
	      if ( star -> quads[i] == -1 ) continue;
	      if ( k == 0 || ( k > 0 && inList (k, quads, star -> quads[i]) == -1 ) ) quads[k++] = star -> quads[i];
	  }
      }
      EG_freeStar(&star);
      qm -> vFix[0] = 0;
      uvxyz = EG_alloc ( 4 * 5 * k * sizeof ( double ));
      if ( uvxyz == NULL ) return EGADS_MALLOC;
      for ( it = i = 0 ; i < k; i++ ) {
	  for ( j = 0 ; j < 4; j++ ) {
	      q = qm -> mesh -> qIdx [ 4 * ( quads[i] - 1) + j ] - 1;
	      uvxyz[it++] = qm -> mesh -> uvs  [ 2 * q    ];
	      uvxyz[it++] = qm -> mesh -> uvs  [ 2 * q + 1];
	      uvxyz[it++] = qm -> mesh -> xyzs [ 3 * q    ];
	      uvxyz[it++] = qm -> mesh -> xyzs [ 3 * q + 1];
	      uvxyz[it++] = qm -> mesh -> xyzs [ 3 * q + 2];
	  }
      }
      for ( i = 0 ; i < k; i++ ) {
	  stat = EG_makePositiveAngles(qm,  qm -> minAngle,  qm -> maxAngle, quads[i] );
	  if ( stat < 0 ) {
	      qm -> vFix[0] = 0;
#ifdef DEBUG
	      printf(" Mesh is invalid for MIN MAX TOTALS [0, 200 deg]\n");
	      sprintf(buffer, "InvalidMesh_%d",MESHPLOTS++);
	      printMesh(qm , buffer);
#endif
	      for ( it = i = 0 ; i < k; i++ ) {
		  for ( j = 0 ; j < 4; j++ ) {
		      q = qm -> mesh -> qIdx [ 4 * ( quads[i] - 1) + j ] - 1;
		      qm -> mesh -> uvs  [ 2 * q    ] = uvxyz[it++];
		      qm -> mesh -> uvs  [ 2 * q + 1] = uvxyz[it++];
		      qm -> mesh -> xyzs [ 3 * q    ] = uvxyz[it++];
		      qm -> mesh -> xyzs [ 3 * q + 1] = uvxyz[it++];
		      qm -> mesh -> xyzs [ 3 * q + 2] = uvxyz[it++];
		  }
	      }
	      EG_free (uvxyz);
	      EG_free (quads );
	      return stat;
	  }
      }
      for ( i = 0; i < k; i++ ) {
	  for ( j = 0 ; j <= 4; j++ ) {
	      if ( j == 0 ) q = quads[i];
	      else q = qm -> mesh -> qAdj [ 4 *( quads[i] - 1) + j - 1 ];
	      if ( q == -1 ) continue;
	      if ( EG_quadAngleOrientation (qm,qm -> minAngle, qm -> maxAngle, q, ori, piv, angles) <= 0 )  {
#ifdef DEBUG
		  printf(" Mesh is invalid for MIN MAX TOTALS [0, 200 deg]\n");
		  sprintf(buffer,"InvalidMesh_%d",MESHPLOTS++);
		  printMesh(qm , buffer);
#endif
		  for ( it = i = 0 ; i < k; i++ ) {
		      for (  j = 0 ; j < 4; j++ ) {
			  q = qm -> mesh -> qIdx [ 4 * ( quads[i] - 1) + j ] - 1;
			  qm -> mesh -> uvs  [ 2 * q    ] = uvxyz[it++];
			  qm -> mesh -> uvs  [ 2 * q + 1] = uvxyz[it++];
			  qm -> mesh -> xyzs [ 3 * q    ] = uvxyz[it++];
			  qm -> mesh -> xyzs [ 3 * q + 1] = uvxyz[it++];
			  qm -> mesh -> xyzs [ 3 * q + 2] = uvxyz[it++];
		      }
		  }
		  EG_free (quads );
		  EG_free (uvxyz);
		  return EGADS_GEOMERR;
	      }
	  }
      }
      EG_free(uvxyz);
      EG_free(quads);
  } else {
      itMax   = 100;
      minT    = 0.0;
      maxT    = PI;
      dthetam = (DEG10 - minT) / (double)itMax;
      dthetaM = (PI - DEG160 ) / (double)itMax;
      stat    = EG_restoreMeshMap (qm, qm -> bestMesh, qm -> mesh);
      while ( it <= itMax ) {
	  qm -> vFix[0] = 0;
#ifdef DEBUG
	  printf("SMOOTH ROUND %d max %d \n ", it, itMax );
#endif
	  for ( i = 0; i < qm -> mesh -> totQuads; i++ ) {
	      if ( qm -> mesh -> qIdx [4 * i ] == -2 ) continue;
	      stat = EG_makePositiveAngles( qm, minT, maxT, i + 1 ) ;
#ifdef DEBUG
	      printMesh (qm, NULL);
#endif
	      if ( stat < 0 ) return stat;
	  }
	  for ( j = i = 0; i < qm -> mesh -> totQuads; i++ ) {
	      if ( qm -> mesh -> qIdx [4 * i ] == -2 ) continue;
	      //if ( EG_quadAngleOrientation (qm,qm -> minAngle, qm -> maxAngle, 1 + i, ori, piv, angles) == 0 ) j = 1;
	      if ( EG_quadAngleOrientation (qm,qm -> minAngle, PI, 1 + i, ori, piv, angles) == 0 ) {
		  printf(" ************ makeValid  ATENCION:: INVALID QUAD  %lf  %lf  ***************\n", minT, maxT );
		  printf("------------ AREA QUAD %d INVALID Internal angles ordered by size --------------\n", i + 1 );
		  for ( k = 0 ; k < 4; ++k)
		    printf("Vertex %d has angle %f and orientation %d \n ", qm -> mesh -> qIdx [ 4 * i + k],
			   angles[piv[k]], ori[piv[k]]);
		  j = 1;
	      }
	  }
	  if ( j     == 0 ) {
	      stat       = EG_restoreMeshMap (qm, qm -> bestMesh, qm -> mesh );
	      printf(" update mesh \n ");
	      printMesh (qm, NULL );
	      EG_wvsData(qm ->  mesh, NULL );
	      if ( stat != EGADS_SUCCESS ) {
		  printf( " EG_makeValidMesh ::: EG_restoreMeshMap stat %d !!\n ", stat );
		  return stat;
	      }
	  }
	  it++;
	  minT += dthetam;
	  maxT -= dthetaM;
      }
      for ( i = 0; i < qm -> mesh -> totQuads; i++ ) {
	  if ( qm -> mesh -> qIdx [4 * i ] == -2 ) continue;
	  if ( EG_makePositiveAngles (qm, qm -> minAngle, PI, 1 + i ) <= 0 )  {
	      EG_free(quads);
#ifdef DEBUG
	      printf(" Mesh is invalid for MIN MAX TOTALS %lf %lf \n", minT, maxT);
	      sprintf(buffer,"InvalidMesh_%d",MESHPLOTS++);
	      printMesh(qm , buffer);
#endif
	      printf(" Mesh is invalid for MIN MAX TOTALS %lf %lf \n", minT, maxT);
	      printMesh(qm , NULL);
	      printf(" Stop at angle interval %f  %f  \n", minT, maxT );
	      EG_wvsData(qm ->  mesh, NULL );
	      return EG_restoreMeshMap (qm, qm -> mesh, qm -> bestMesh );
	  }
      }
  }
  printf(" LEAVE VALID MESH SUCCESSFULLY \n ");
  return EGADS_SUCCESS;
}


static int EG_makePositiveAngles (meshMap *qm, double minAngle, double maxAngle, int qID) {
  int i, j, k, kk, iO, iA, iB, it , itMAX = 50, i0, area = 0,  *sa = NULL, fullNeg = 0 ;
  int  ori[8], piv[8], quad[4], type[4], vfix[4], nb[10];
  static int pivTri[6] = { 1, 2, 2, 3, 1, 3 };
  double ang, angles[8], dt = 0.01, ds = 0.95, xyzA[3], xyzB[3], dista, distb, uv[2], totDist = 0.0;
  vStar *star = NULL;
  for ( j = 0; j < 4; j++ ) {
      quad[j] = qm -> mesh -> qIdx [ 4 * ( qID - 1 ) + j];
      type[j] = qm -> mesh -> vType[    quad[j] -1 ] ;
      vfix[j] = 0;
      //if ( vertexLinksToBounds (qm, quad[j], nb) > 0 && maxAngle > PI  ) maxAngle = PI;
      if (maxAngle > PI  ) maxAngle = PI;
  }
#ifdef DEBUG
  printf("makePositiveAngle for quad %d \n ", qID );
#endif
  for ( it = 0 ; it < itMAX; it++ ) {
      area      = EG_quadAngleOrientation (qm, minAngle, maxAngle, qID, ori, piv, angles);
      for ( j   = 0 ; j < 4; j++ ) {
	  if ( type [ piv[j] ] != -1 || vfix[piv[j]]  == 1 ) continue;
	  k     = quad[ piv[j]] - 1;
	  uv[0] = qm -> mesh -> uvs [ 2 * k   ];
	  uv[1] = qm -> mesh -> uvs [ 2 * k +1];
	  i     = EG_buildStar (qm -> mesh, &star, k + 1);
	  if(i != EGADS_SUCCESS || star == NULL ) return EGADS_MALLOC;
	  sa    = EG_alloc ( star -> nQ * sizeof ( int ) );
	  if ( sa == NULL ) {
	      EG_freeStar ( &star );
	      return EGADS_MALLOC;
	  }
	  for ( i = 0 ; i < star -> nQ; i++ ) {
	      if ( star -> quads[i] == -1 || star -> quads[i] == qID ) continue;
	      sa[i] = EG_quadAngleOrientation (qm, minAngle, maxAngle, star -> quads[i], &ori[4], &piv[4], &angles[4]);
#ifdef DEBUG
	      printf(" i %d Q %d ---> area %d \n", i, star -> quads[i], sa[i]);
#endif
	  }
	  EG_averageCoords (qm, k + 1, minAngle, maxAngle);
	  for ( i = 0 ; i < star -> nQ; i++ ) {
	      if ( star -> quads[i] == -1 || star -> quads[i] == qID ) continue;
	      area = EG_quadAngleOrientation (qm, minAngle, maxAngle, star -> quads[i], &ori[4], &piv[4], &angles[4]);
	      if ( area <= 0 && sa[i] == 1 ) {
#ifdef DEBUG
		  printf(" At vertex %d !! invalid area quad %d -> %d\n ", k + 1, star -> quads[i], area );
		  printf(" RESET COORDS %d \n ", k + 1 );
		  printQuadCoords (qm, star -> quads[i] );
#endif
		  updateVertex ( qm, k + 1, uv );
#ifdef DEBUG
		  printQuadCoords (qm, star -> quads[i] );
#endif
		  break;
	      }
	  }
	  EG_freeStar ( & star );
	  EG_free ( sa );
      }
      if ( EG_quadAngleOrientation (qm, minAngle, maxAngle, qID, ori, piv,angles) != 0 ) break;
#ifdef DEBUG
      printf(" it %d area %d min %f max %f\n", it, area,  minAngle, maxAngle);
#endif
      for ( fullNeg = i0 = j = 0 ; j < 4; j++ ) {
	  fullNeg += ori[j];
	  if ( type[piv[j]] >= 0 && ( type[piv[i0]] == -1 || angles[piv[j]] > angles [piv[i0]] ) ) i0 = j;
      }
#ifdef DEBUG
      printQuadCoords (qm, qID);
#endif
      for ( j = 0; j < 4; j++ ) piv[j] = (piv[i0] +  j ) %4 ;
#ifdef DEBUG
      if ( it %5 == 0 ) printMesh(qm, NULL);
#endif
      totDist   = 0.0;
      for ( j   = 0; j < 4; j++ ) {
	  iO    = piv[j];
	  if ( type [ iO ] == -1  && ori[iO] == -1 && fullNeg != 4 ) {
	      nb[0] = vertexLinksToBounds (qm, quad[iO], nb  );
	      if ( nb[0] == 0 || ( nb[0] > 0 && qm -> mesh -> vType [ nb[1] -1 ] < 4 )) continue;
	  }
#ifdef DEBUG
	  printf ("=========== O ( %d ) = %d ::: MAX ANGLE PASS %f %f   ==============================\n", iO, quad[iO],  minAngle, maxAngle );
#endif
	  for (i  = 0; i < 3; i++ ) {
	      iA  = (iO + pivTri[2 * i    ]) %4;
	      iB  = (iO + pivTri[2 * i + 1]) %4;
	      for ( kk = 0; kk < 3; kk ++ ) {
		  xyzA[kk] = qm -> mesh -> xyzs [ 3 * ( quad[iA] -1 ) + kk];
		  xyzB[kk] = qm -> mesh -> xyzs [ 3 * ( quad[iB] -1 ) + kk];
	      }
#ifdef DEBUG
	      printQuadCoords (qm, qID );
#endif
	      EG_3Drotation (qm, minAngle, maxAngle, dt,  quad[iO], quad[iA], quad[iB]);
	      ang   = EG_angleAtVnormalPlane  (qm, quad[iO], quad[iA], quad[iB]);
	      dista = 0.0; distb = 0.0;
	      for ( kk = 0 ; kk < 3; kk++ ) {
		  dista += (xyzA[kk] - qm -> mesh -> xyzs [ 3 * ( quad[iA] -1 ) + kk] ) *
		      (xyzA[kk] - qm -> mesh -> xyzs [ 3 * ( quad[iA] -1 ) + kk] );
		  distb += (xyzB[kk] - qm -> mesh -> xyzs [ 3 * ( quad[iB] -1 ) + kk] ) *
		      (xyzB[kk] - qm -> mesh -> xyzs [ 3 * ( quad[iB] -1 ) + kk] );
	      }
	      dista = sqrt ( dista);
	      distb = sqrt ( distb);
	      if ( dista + distb > 0.0 && dista  < 0.1 * qm -> minsize && distb < 0.1 * qm -> minsize) {
		  dt += 0.01;
		  ds -= 0.01;
	      }
#ifdef DEBUG
	      printf( " it %d ----> dt %f ds %f distances %f  %f  :\n Angle i %d -> %d %d %d is now %f  \n ",
		      it, dt, ds, dista, distb, i, quad[iO], quad[iA], quad[iB], ang );
	      printQuadCoords (qm, qID );
#endif
	      totDist += dista + distb;
	      if ( (ang > maxAngle || ang < minAngle ) && dista + distb > 0.0  ) {
		  if ( it == 0 ) {
		      if ( dista > 0.0 ) vfix[iA] = 1;
		      if ( distb > 0.0 ) vfix[iB] = 1;
#ifdef DEBUG
		      printf(" --------- ATTENTION IN QUAD %d fixing angles %d %d \n ", qID, quad[iA], quad[iB] );
#endif
		  }
		  break;
	      }
	  }
	  if ( (ang > maxAngle || ang < minAngle ) && dista + distb > 0.0 ) break;
      }
      if ( totDist < EPS11 ) break;
  }
  area = EG_quadAngleOrientation (qm, minAngle, maxAngle,  qID, ori, piv, angles);
#ifdef DEBUG
  printQuadCoords (qm, qID );
  printf("\n\n=============================================================\n");
  printf(" LEAVE MAKEPOSITIVE QUAD %d IS NOW AREA %d \n ", qID, area );
  printf("\n\n=============================================================\n");
#endif
  if ( area == 1 )
    for ( j = 0 ; j < 4; j++ ) {
	if ( vfix[j] == 0 ) continue;
	if (qm -> vFix[0] == 0 || inList ( qm -> vFix[0], &qm -> vFix[1], quad[j] ) == -1 ) {
#ifdef DEBUG
	    printf(" ATENTION: FIXING vert %d \n ", quad[j] );
#endif
	    qm -> vFix[++qm -> vFix[0] ] = quad[j];
	}
    }
  return area;
}


static void EG_adjustPercentage (meshMap *qm, double totarc, double maxPerc, int n, double *list ) {
  int k, i;
  double extra, x;
  for ( k = 0 ; k < n; k++ ) {
      x = list[k] * 100.0 / totarc ;
      if ( x > maxPerc ) {
	  extra   = list[k] - totarc * maxPerc / 100.0;
	  list[k] = totarc * maxPerc / 100.0;
	  for ( i = 0 ; i < n; i++ ) {
	      if ( i == k ) continue;
	      list[i] += extra / (double)( n - 1);
	  }
      }
  }
}

static void EG_averageCoords (meshMap *qm, int vID, double minAngle, double maxAngle ) {
  int  vA, vB, i, j, k, auxID = 0, stat, i0, piv[4], ori[4], nb[10], n;
  double uvc[2], angle, angles[4], *length = NULL, totarc, xyz[18],  extra, arc, uva[2] ;
  vStar *star = NULL;
  if (qm -> mesh -> vType[ vID -1] != -1  ||
      (qm -> vFix[0] > 0 && inList ( qm -> vFix[0], &qm -> vFix[1], vID ) >= 0  )) return;
  stat   = EG_buildStar ( qm -> mesh, &star, vID);
  if ( star == NULL || stat != EGADS_SUCCESS) return;
  nb[0]      = vertexLinksToBounds (qm, vID, nb );
  length     = EG_alloc ( star -> nQ * sizeof ( double ) ) ;
  if (length == NULL ) {
      EG_free (star );
      return;
  }
  totarc = 0.0;
  uva[0] = qm -> mesh -> uvs [ 2 * ( vID -1 )    ];
  uva[1] = qm -> mesh -> uvs [ 2 * ( vID -1 ) + 1];
  for (k = 0; k < star -> nQ; k++) {
      vA =  star -> verts [ 2 * k + 1] - 1;
      if ( qm -> mesh -> vType [ vA ] != -1 ) {
	  uvc[0] = qm -> mesh -> uvs [  2 * vA     ] +
	      0.1 * (qm -> mesh -> uvs [ 2 * ( vID - 1 )    ] - qm -> mesh -> uvs [  2 * vA    ]);
	  uvc[1] = qm -> mesh -> uvs [  2 * vA + 1 ] +
	      0.1 * (qm -> mesh -> uvs [ 2 * ( vID - 1 ) + 1] - qm -> mesh -> uvs [  2 * vA + 1]);
	  i      = EG_curvature ( qm  -> face, uvc, xyz );
      }
      else i     = EG_curvature ( qm  -> face, &qm -> mesh -> uvs [ 2 * vA ], xyz );
      length[k]  = MAX (fabs ( xyz[0]), fabs (xyz[4]) );
      arc        = EG_segment   ( qm -> face, &qm -> mesh -> uvs [ 2 * ( vID -1 ) ],
				  &qm -> mesh -> uvs [ 2 * vA ] );
      printf(" V %d curvature %f arc %f VAL %f \n ", vA + 1, length[k], arc , arc * length[k] );
      if ( length[k] < EPS11 ) length[k]  = arc;
      else                     length[k] *= arc;
      totarc += length[k];
      uva[0] += qm -> mesh -> uvs [ 2 * vA     ] ;
      uva[1] += qm -> mesh -> uvs [ 2 * vA + 1 ] ;
  }
  uva[0] /= (double) ( star -> nQ + 1 );
  uva[1] /= (double) ( star -> nQ + 1 );
  EG_adjustPercentage (qm, totarc, 50.0, star -> nQ, length);
  extra   = 0.0;
  for ( k = 0 ; k < star -> nQ; k++ ) {
      extra += length[k] * 100.0 / totarc ;
      printf(" V %d = %d has length %f perc %f accc %f\n ", k, star -> verts [ 2 * k + 1] , length[k], length[k] * 100.0 / totarc, extra );
  }
  uvc[0] = 0.0; uvc[1] = 0.0;
#ifdef DEBUG
  printVertexCoords (qm, vID );
#endif
  for  ( k = i = 0; i < star -> nQ; i++ ) {
      if ( EG_quadAngleOrientation (qm, minAngle, maxAngle,  star -> quads[i], ori, piv, angles) == 0 ) continue;
      vA    = star -> verts [ 2 * i + 1];
      vB    = star -> verts [ star -> idxV [ 2 * i + 3]];
      angle = EG_angleAtVnormalPlane (qm, vID, vA, vB);
      if      ( angle > maxAngle ) k =  1;
      else if ( angle < minAngle ) k = -1;
      else continue;
      totarc   = 0.0;
      for ( n = j = 0; j < star -> nQ; j++ ) {
	  auxID = star -> verts [ 2 * j + 1 ] - 1;
	  if ((k ==  1 && ( auxID + 1 == vA || auxID + 1 == vB )) ||
	      (k == -1 &&   auxID + 1 != vA && auxID + 1 != vB )) continue;
	  stat = EG_curvature ( qm  -> face, &qm -> mesh -> uvs [ 2 * ( star -> verts [ 2 * j + 1] - 1 ) ], xyz );
	  printf(" Curvature %f  %f \n ", xyz[0], xyz[4] );
	  length[n] = MAX (fabs ( xyz[0]), fabs (xyz[4]) );
	  if (length[n] < EPS11 )
	    length[n] = EG_segment ( qm -> face, &qm -> mesh -> uvs [ 2 * ( vID -1 ) ],
				     &qm -> mesh -> uvs [ 2 * ( star -> verts [ 2 * j + 1] - 1 ) ] );
	  printf(" Weight %d => %lf \n ", j, length[n]);
	  totarc += length[n];
	  n++;
      }
      EG_adjustPercentage (qm, totarc, 50.0, n, length);
      extra = 0.0;
      for ( n = j = 0; j < star -> nQ; j++ ) {
	  auxID = star -> verts [ 2 * j + 1 ] - 1;
	  if ((k ==  1 && ( auxID + 1 == vA || auxID + 1 == vB )) ||
	      (k == -1 &&   auxID + 1 != vA && auxID + 1 != vB )) continue;
	  extra  += length[n] * 100.0 / totarc;
	  printf(" PULLING V %d = %d has length %f perc %f accc %f\n ", k, star -> verts [ 2 * j + 1] , length[n], length[n] * 100.0 / totarc, extra );
	  uvc[0] += length[n  ] * qm -> mesh -> uvs [2 * auxID    ];
	  uvc[1] += length[n++] * qm -> mesh -> uvs [2 * auxID + 1];
      }
      //uvc[0] /= totarc;
      //uvc[1] /= totarc;
      uvc[0] = 0.5 * (uvc[0] / totarc + uva[0] );
      uvc[1] = 0.5 * (uvc[1] / totarc + uva[1] );
      updateVertex (qm, vID, uvc);
      break;
  }
  if ( k == 0 && nb[0] > 0 ) {
      if ( nb[0] != 2 || ( nb[0] == 2 &&
	  ( qm -> mesh -> vType [ nb[1] -1 ] < 4 && qm -> mesh -> vType [ nb[2] -1] < 4 ) ) ) {
	  for  ( j  = i0 = k = 0; k < star -> nQ; k++ ) {
	      auxID = star -> verts [ 2 * k + 1] - 1;
	      if ( qm -> mesh -> vType [ auxID] >= 0 ) {
		  if ( i0 == 0 ) i0 = 2 * k + 1;
		  else {
		      j = 2 * k + 1;
		      if ( j - i0 > 2 && ( i0 + ( star -> nV -1 ) - j != 2 ) ) {
			  i0 = star -> nV;
			  break;
		      }
		  }
	      }
	  } if ( i0 == star -> nV ) {
	      uvc[0]      = 0.0;
	      uvc[1]      = 0.0;
	      for  ( k    = 0; k < nb[0]; k++ ) {
		  uvc[0] += qm -> mesh -> uvs [2 * ( nb[1 + k] -1  )    ];
		  uvc[1] += qm -> mesh -> uvs [2 * ( nb[1 + k] -1  ) + 1];
	      }
	      uvc[0] /= (double)nb[0];
	      uvc[1] /= (double)nb[0];
	      updateVertex (qm, vID, uvc);
	      k = 1;
	  } else k = 0;
      }
  }
  if ( k == 0  ) {
      totarc = 0.0;
      for (k = 0; k < star -> nQ; k++) {
	  auxID   = star -> verts [ 2 * k + 1] - 1;
	  totarc += length[k];
	  uvc[0] += length[k] * qm ->mesh ->uvs [2 * auxID    ];
	  uvc[1] += length[k] * qm ->mesh ->uvs [2 * auxID + 1];
      }
      //uvc[0] /= totarc;
      //uvc[1] /= totarc;
      uvc[0] = 0.5 * (uvc[0] / totarc + uva[0] );
      uvc[1] = 0.5 * (uvc[1] / totarc + uva[1] );
      printVertexCoords (qm, vID );
      updateVertex (qm, vID, uvc);
      printVertexCoords (qm, vID );
  }
  EG_free ( length );
  EG_freeStar (&star );
#ifdef DEBUG
  printf(" LEAVE WITH COORDS\n" );
  printVertexCoords (qm, vID );
#endif
  return;
}

/* v should be dimension 3: v[0] = n common verts: v[1] (v[2]) idx*/
static void EG_commonVerts ( meshData *mesh, int q1, int q2, int *v ) {
  int i, j, k;
  v[0] = 0;
#ifdef DEBUG
  if (checkQuad ( mesh, q1 ) != EGADS_SUCCESS || checkQuad ( mesh, q2 ) != EGADS_SUCCESS ) {
      printf(" EG_commonVerts wrong quad pair %d %d\n", q1, q2);
      printQuadSpecs (mesh, q1);
      printQuadSpecs (mesh, q2);
      v[0] = -1;
      return;
  }
#endif
  if ( q1 == q2 ) {
      for ( k = j = 0 ; j < 4 ; j++ )
	v[++k] = mesh -> qIdx [ 4 * ( q2 - 1 ) + j ];
      return;
  }
  for (k = i = 0 ; i < 4 ; i++ ) {
      for ( j = 0 ; j < 4 ; j++ ) {
	  if (mesh -> qIdx [ 4 * ( q1 - 1 ) + i ] ==
	      mesh -> qIdx [ 4 * ( q2 - 1 ) + j ] )
	    v[++k] = mesh -> qIdx [ 4 * ( q2 - 1 ) + j ];
      }
  }
  v[0] = k;
}

static int EG_quadVertIdx ( meshData *mesh, int q, int v ) {
  int i = 0 ;
  for ( i = 0 ; i < 4; i++ )
    if ( mesh -> qIdx [ 4 * ( q - 1 ) + i ] == v ) return i;
  return -1;
}

static int getValence ( meshData *mesh, int v )  {
  int val;
#ifdef DEBUG
  val = checkVertex (mesh, v );
  if ( val != EGADS_SUCCESS ) return val;
#endif
  val =     mesh -> valence [ v - 1 ][0];
  if      ( mesh -> vType   [ v - 1 ] == 0 ) val += 2;
  else if ( mesh -> vType   [ v - 1 ] == 3 ) val++;
  else if ( mesh -> vType   [ v - 1 ] == 5 ) val--;
  return val;
}

static int setValence ( meshData *mesh, int vID ) {
  int i;
  vStar *star = NULL;
#ifdef DEBUG
  i = checkVertex ( mesh, vID );
  if ( i != EGADS_SUCCESS ) {
      printf(" set valence for vertex %d is %d\n ", vID, i );
      return i;
  }
#endif
  i       = EG_buildStar ( mesh, &star, vID );
  if ( i != EGADS_SUCCESS || star == NULL ) {
      printf(" In setValence for vertex %d stat from EG_buildStar = %d\n", vID, i );
      return i;
  }
  mesh -> valence[vID -1 ][0] = star -> nQ;
  for ( i = 0 ; i < star -> nQ; i++ ) {
      if ( star -> verts[ 2 * i + 1] == -1 ) continue;
      mesh -> valence[vID -1 ][2 + i] = star -> verts [ 2 * i + 1] ;
  }
  EG_freeStar ( &star );
  return EGADS_SUCCESS;
}

static int EG_nValenceCount ( meshData *mesh, int q, int n ) {
  int i, count, val;
  for ( count = i = 0 ; i < 4; i++ ) {
      val = getValence ( mesh, mesh ->qIdx [ 4 * ( q - 1) + i ] );
      if ( val < 0 ) return val;
      if ( n > 5 ) {
	  if ( val > 5 ) count++;
      }
      else if ( n == val ) count++;
  }
  return count;
}

static int checkVertex ( meshData *mesh, int v ) {
  int q, i;
  if ( v == -2  ) return EGADS_EMPTY;
  if ( v <= 0  || v > mesh ->totVerts )           {
      printf(" Vertex for v %d is out of bounds !! \n ", v);
      return EGADS_INDEXERR;
  }
  if (  mesh ->vType [ v - 1 ] == -2 ) {
      printf(" V %d is empty !\n ", v );
      return EGADS_EMPTY;
  }
  if (( mesh ->vType [ v - 1 ] == 0 && mesh ->valence [ v - 1][0] < 2 ) ||
      ( mesh ->vType [ v - 1 ]  > 0 && mesh ->valence [ v - 1][0] < 3 ) ) {
      printf(" checkVertex:: Vertex %d Type %d is a doublet!!! We have big problems!!!\n", v , mesh ->vType[v - 1]);
      return EGADS_GEOMERR;
  }
  q = mesh ->valence [ v - 1][1] - 1;
  if ( q < 0 || q >= mesh ->totQuads ) {
      printf(" checkVertex:: Vertex for v %d has associated a quad %d out of bounds ( NMAX %d ) !! \n ", v, q + 1, mesh ->totQuads -1 );
      return EGADS_INDEXERR;
  }
  for ( i = 0 ; i < 4; i++) {
      if ( mesh ->qIdx [ 4 * q  + i] == -2 ) {
	  printf("checkVertex::  Vertex for v %d has associated quad %d but  quad is empty !! \n ", v, q + 1);
	  return EGADS_EMPTY;
      }
      if ( mesh ->qIdx [ 4 * q  + i] ==  v )   return EGADS_SUCCESS;
  }
  printf("checkVertex::  Vertex for v %d has associated quad %d but !! \n ", v, q + 1);
  printQuadSpecs ( mesh, q + 1);
  return EGADS_INDEXERR;
}

static int checkQuad ( meshData *mesh, int q ) {
  int i, j;
  if ( q > mesh -> totQuads || q <= 0  || mesh ->qIdx [ 4 * ( q - 1)  ] == -2 ) {
      printf(" q %d -> max %d ass %d \n ", q, mesh -> totQuads, mesh ->qIdx [ 4 * ( q - 1)  ]);
      return EGADS_EMPTY;
  }
  for ( i = 0 ; i < 4; i++) {
      j = mesh ->qIdx [ 4 * ( q - 1) + i];
      if ( j == -2 ) {
	  printf("j = %d \n ", j );
	  return EGADS_EMPTY;
      }
#ifdef DEBUG
      j = checkVertex (mesh, mesh ->qIdx [ 4 * ( q - 1) + i] );
      if (j != EGADS_SUCCESS ) {
	  printf(" VERTEX %d is %d \n ", mesh ->qIdx [ 4 * ( q - 1) + i] , j );
	  return j;
      }
#endif
  }
  return EGADS_SUCCESS;
}

static int checkMesh(meshMap *qm) {
  int stat, i, j, k,  val1, val2, v1, v2;
  for ( i = 0 ; i < qm -> mesh ->totVerts; i++ ) {
      if (qm -> mesh -> vType[i] == -2 )  continue;
      val1 = qm -> mesh  -> valence[i][0];
      val2 = qm -> mesh  -> valence[i][1];
      if ( val1 <= 0 || val1 > qm -> mesh  -> totVerts ) {
	  stat = setValence ( qm -> mesh, i + 1 );
	  if ( stat != EGADS_SUCCESS ) {
	      fprintf(stderr," checkMesh vertex %d setValence --> %d!!\n", i + 1, stat );
	      return stat;
	  }
      }
      if ( val2 <= 0 || val2 > qm -> mesh -> totQuads ) {
	  printf("In checkMesh quad for %d is out of bounds!!! %d > %d \n ", i + 1, val2, qm -> mesh -> totQuads);
	  return EGADS_INDEXERR;
      }
      stat = checkQuad (qm -> mesh, val2 );
      if ( stat != EGADS_SUCCESS ) {
	  fprintf(stderr," checkMesh vertex %d has associated quad %d --> %d!!\n", i + 1, val2, stat);
	  return stat;
      }
      v1  = i + 1;
      for ( j = 0 ; j < val1; ++j ) {
	  v2  = qm -> mesh  -> valence[i][2 + j] - 1;
	  for ( k = 0 ; k < qm -> mesh ->valence [v2][0]; k++ )
	    if ( qm -> mesh ->valence [v2][2 + k] == v1 ) break;
	  if ( qm -> mesh ->valence [ v2 ][2 + k] != v1 ) {
	      fprintf(stderr," checkMesh Vertex %d has assigned %d as link but %d doesn't point at %d\n", v1, v2, v2, v1 );
	      return EGADS_INDEXERR;
	  }
      }
  }
  for ( i = 0 ; i < qm -> mesh  -> totQuads; i++ ) {
      if ( checkQuad (qm -> mesh, i+1) == EGADS_EMPTY ) continue;
      for ( j = 0; j < 4; ++j ) {
	  v1 = qm -> mesh  -> qAdj [ 4 * i  + j] - 1;
	  if ( v1 < 0 ) continue;
	  val1 = -1;
	  for ( k = 0 ; k < 4; ++k ) {
	      if ( qm -> mesh  -> qAdj [ 4 * v1  + k] == i + 1 ) {
		  val1 = 1;
		  break;
	      }
	  }
	  if ( val1 == -1 ) {
	      fprintf (stderr, " checkMesh quads %d and %d don't point at each other\n", i +1, v1 + 1);
	      printQuadSpecs ( qm -> mesh , i + 1 ) ;
	      printQuadSpecs ( qm -> mesh , v1 + 1 ) ;
	      return EGADS_INDEXERR;
	  }
      }
  }
  return EGADS_SUCCESS;
}

static void meshCount( meshData *mesh, int *nI, int *nV, int *nQ ) {
  int i, qSum, vSum = 0, vSum2 = 0;
  for ( vSum = i = qSum = 0; i< mesh -> totQuads;i++ )
    if ( mesh -> qIdx [4 * i ] != -2 ) qSum++;
  for ( vSum =  i = 0 ; i < mesh -> totVerts; i++ ) {
      if ( mesh -> vType[i] == -2 ) continue;
      vSum2++;
      if ( getValence ( mesh, i + 1 ) != 4) vSum++;
  }
  *nQ = qSum;
  *nI = vSum;
  *nV = vSum2;
}

/* Assuming qID is collapsing through v */
static int validCollapse ( meshMap *qm, int qID, int v ) {
  int j, k, kk, id, link, aux, aux2, type, val ;
  id  = EG_quadVertIdx ( qm -> mesh, qID, v);
  aux = qm -> mesh ->qIdx [ 4 * ( qID - 1 ) + (id + 2)%4];
  for ( j = 0 ; j < 2; j++ ) {
      link = qm -> mesh ->qIdx [ 4 * ( qID - 1 ) + (id + 2 * j + 1)%4];
      type = qm -> mesh ->vType  [link - 1];
      val  = qm -> mesh ->valence[link - 1][0];
      if      ((type >   0 && val <= 3 ) || (type ==  0 && val <= 2 )) return 0;
      else if  (type == -1 && val == 3 ) {
	  for ( k = 0 ; k < 3; k++ ) {
	      aux = qm -> mesh ->valence [ link - 1 ][2 + k] - 1;
	      if ( EG_quadVertIdx ( qm -> mesh, qID,  aux + 1 ) >= 0 || qm -> mesh ->valence[ aux ][0] != 3) continue;
	      if (      qm -> mesh ->vType[aux] > 0 ) return 0;
	      else if ( qm -> mesh ->vType[aux] != -1 ) continue;
	      for ( kk = 0 ; kk < 3; kk++ ) {
		  aux2 = qm -> mesh ->valence [ aux ][2 + kk] - 1;
		  if ( aux2 + 1 == link ) continue;
		  if ( qm -> mesh ->valence[ aux2 ][0] == 3 && qm -> mesh ->vType[aux] != 0 ) return 0;
	      }
	  }
      }
  }
  return 1;
}
/* Assuming we will break link v1-v2 */
static int validSwap ( meshData *mesh, int v1, int v2 ) {
  int i, vs[2];
  vs[0] = v1 - 1; vs[1] = v2 - 1;
  for ( i = 0; i < 2; i++ ) {
#ifdef DEBUG
      if      (checkVertex ( mesh, vs[i] + 1 ) != EGADS_SUCCESS )                return 0;
#endif
      if      (mesh ->vType[ vs[i] ] >=  4 && mesh ->valence[ vs[i] ][0] == 4 )  return 0;
      if      (mesh ->vType[ vs[i] ] !=  0 && mesh ->valence[ vs[i] ][0] == 3 )  return 0;
      else if (mesh ->vType[ vs[i] ] ==  0 ) {
	  if      ( mesh ->valence[ vs[i] ][0] == 2  ) return 0;
	  else if ( mesh ->valence[ vs[i] ][0] == 3 &&
	      mesh ->valence[ vs[( i + 1) %2] ][0] == 3 ) return 0;
      }
  }
  return 1;
}

/*  0  = interior
 *  1  = edge
 *  2  = corner // irregular vertex at boundary
 * < 0 = quad doesnt exist
 * NOTE: We consider a quad to be boundary if it has at least one vertex type >=0
 * even though its adjacents might all live inside the domain. This happens for high
 * valences
 */
static int EG_quadIsBoundary ( meshData *mesh, int qID ) {
  int i, q, v, bV = 0 , bQ = 0 ;
#ifdef DEBUG
  i = checkQuad ( mesh, qID );
  if ( i != EGADS_SUCCESS)  return i;
#endif
  for (i  = 0 ; i < 4; i++ ) {
      v   = mesh ->qIdx [ 4 * ( qID - 1) + i ];
      q   = mesh ->qAdj [ 4 * ( qID - 1) + i ];
      if ( mesh ->vType[v - 1] >=  0 ) bV = 1;
      if ( q                   == -1 ) bQ = 1;
  }
  if ( bQ == 1 ) return 1;
  if ( bV == 1 ) return 2;
  else return 0;
}


static int EG_createQuadGroup ( meshData *mesh, quadGroup *qg, int q0, int q1 ) {
  int i,  ids[2], piv = 0, aux, vaux[6], common[3], centre, c1 = 0, c2 = 0;
  qg -> q[0]     = q0;
  qg -> q[1]     = q1;
  for ( i = 0 ; i < 6; i++) {
      qg -> verts[i] = -1;
      qg -> vals [i] = -1;
  }
  if ( q0 < 0 || q0 > mesh -> totQuads || q1 < 0 || q1 > mesh -> totQuads ) {
      printf(" Q0 %d Q1 %d beyond limits !!!! \n ", q0, q1 );
      return EGADS_INDEXERR;
  }
  EG_commonVerts (mesh, q0, q1, common );
  if ( common[0] != 2 ) {
      printf(" You have tried to create a quad group where quads are not adjacent\n");
      printQuadSpecs(mesh, q0);
      printQuadSpecs(mesh, q1);
      return EGADS_INDEXERR;
  }
  centre  = common[1];
  if ( getValence (mesh, common[2] ) > getValence (mesh, common[1] ) ) centre = common[2];
  ids[0] = EG_quadVertIdx (mesh, qg -> q[0], centre);
  ids[1] = EG_quadVertIdx (mesh, qg -> q[1], centre);
  piv    = 0;
  if ( mesh -> qAdj [ 4 * ( qg ->q[0] - 1 ) + ids[0] ] == qg ->q[1] ) piv = 1;
  for ( i = 0 ; i < 4; i++ )
    qg -> verts[i] = mesh -> qIdx [ 4 * (qg ->q[piv] - 1) + (ids[piv] + i ) % 4 ];
  aux = ( piv + 1 ) % 2;
  qg -> verts[4] = mesh -> qIdx [ 4 * (qg ->q[aux] - 1)  + ( ids[ aux ] + 2 ) % 4 ];
  qg -> verts[5] = mesh -> qIdx [ 4 * (qg ->q[aux] - 1)  + ( ids[ aux ] + 3 ) % 4 ];
  if ( piv == 1 ) swapInt ( &qg ->q[0], &qg ->q[1] );
  for ( i = 0 ; i < 6; i++ ) qg -> vals  [i] = getValence ( mesh, qg->verts[i]);
  if ( qg -> vals[0] == qg -> vals[3] ) {
      if ( qg -> vals[1] != 4 ) c1++;
      if ( qg -> vals[5] != 4 ) c1++;
      if ( qg -> vals[2] != 4 ) c2++;
      if ( qg -> vals[4] != 4 ) c2++;
      if ( c2 > c1 ) {
	  for ( i = 0 ; i < 6; i++ ) vaux[i] = qg -> verts[i];
	  for ( i = 0 ; i < 6; i++ ) qg -> verts[i] = vaux[ (i + 3) % 6] ;
	  for ( i = 0 ; i < 6; i++ ) qg -> vals  [i] = getValence ( mesh, qg->verts[i]);
	  swapInt ( &qg ->q[0], &qg ->q[1] );
      }
  }
  return EGADS_SUCCESS;
}

static int EG_doubleCollapse ( meshMap *qm, quadGroup qg, int forcing, int *activity ) {
  int i, stat, i3, i5, q3, f = 0, remq = 0;
  *activity  = 0;
  if (qg.vals[0] * qg.vals[3] == 16 ) {
      for ( i = i3 = 0 ; i3 < 6; i3++ ) {
	  if ( qg.vals[i3] == 3 ) {
	      i = 1;
	      break;
	  }
      }
      if ( i == 0 ) return EGADS_SUCCESS;
      i5 = ( i3 + 1 ) % 6; if ( i5 % 3 == 0 ) i5 = ( i3 + 5 ) % 6;
      if (                    forcing == 0 && (qg.vals[i5]!= 5 ||
	  (qg.vals[(i5 + 3) %6]   != 5 && qg.vals [ (i3 + 3) % 6 ] != 3 ) ) ) return EGADS_SUCCESS;
      else if ( forcing == 1) {
	  if ( qm -> extraQuads < 0 ) return EGADS_SUCCESS;
	  if ( qg.vals[(i3 + 3)%6] == 3 )      i5 = (i3 + 3)%6;
	  else if ( qg.vals[(i3 + 2)%6] == 5 ) i5 = (i3 + 2)%6;
	  else if ( qg.vals[(i3 + 4)%6] == 5 ) i5 = (i3 + 4)%6;
	  else return EGADS_SUCCESS;
      }
  }
  else if (qg.vals[0] * qg.vals[3] == 12 ) {
      i3 = 0 ;            if ( qg.vals[3 ] == 3 ) i3 = 3;
      i5 = ( i3 + 1 )% 6; if ( qg.vals[i5] != 5 ) i5 = ( i3 + 5 )% 6;
      if ( qg.vals[i5] != 5 || qg.vals [( i5 + 3 ) % 5] < 5 ) return EGADS_SUCCESS;
  }
  else return EGADS_SUCCESS;
  q3 = 0;
  if ( forcing == 1 && qg.vals[0] * qg.vals[3] == 16 ) {
      f = 1;
      if ( i3 >= 3 ) q3 = 1;
  } else {
      if ( i5 >= 3 ) q3 = 1;
  }
  for ( i = 0 ; i < 2; i++ ) {
      if ( qm -> mesh -> qIdx[ 4 * ( qg.q[( q3 + i ) %2] -1 )  ] == -2 ) continue;
      stat       = EG_collapse (qm, qg.q[( q3 + i ) %2], &i3, 1 );
      if ( stat == EGADS_SUCCESS ) {
	  if ( i3 == 0 ) break;
	  else if ( f == 1 ) remq++;
	  (*activity)++;
      }
      if ( stat != EGADS_SUCCESS ) {
	  printf("EG_doubleCollapse after EG_forceColapse quad %d -> %d !!\n ",qg.q[( q3 + i ) %2], stat  );
	  return stat;
      }
  }
  qm -> extraQuads -= remq;
  return EGADS_SUCCESS;
}

static int
EG_swapDoubleCollapse (meshMap *qm, quadGroup qg, int *activity ) {
  int  k, swap = 1, id, j, stat, v0;
  *activity  = 0;
  if (    qg.vals[0] * qg.vals[3] != 20 ||
      qg.vals[2] * qg.vals[4] != 9  ||
      validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0) return EGADS_SUCCESS;
  if ( qg.vals[swap] * qg.vals[(swap+3)%6] != 12 ) swap = 2;
  if ( qg.vals[swap] * qg.vals[(swap+3)%6] != 12 ) return EGADS_SUCCESS;
  v0         = qg.verts[0];
  stat       = EG_swappingOperation(qm, qg, swap, &(*activity) );
  if ( stat != EGADS_SUCCESS || *activity == 0 ) {
#ifdef DEBUG
      printf( "In EG_swapDoubleCollapse: Swapping went %d activity %d \n", stat, *activity );
#endif
      return stat;
  }
  *activity = 1;
  for ( k   = 0; k < 2; k++ ) {
      id    = EG_quadVertIdx ( qm -> mesh, qg.q[k], qg.verts[3] );
      if ( id >= 0 ) break;
  }
  if ( id < 0 ) {
      printf(" I can't find vertex %d in quads!!! \n", qg.verts[3]);
      printQuadGroup (qm -> mesh, qg);
      return EGADS_INDEXERR;
  }
  stat        = EG_collapse ( qm, qg.q[k], &( *activity), 1 );
  if ( stat == EGADS_SUCCESS && *activity > 0 ) return stat;
  stat    = EG_createQuadGroup (qm -> mesh, &qg, qg.q[0], qg.q[1]);
  if ( stat != EGADS_SUCCESS ) {
      printf("EG_swapDoubleCollapse after invalid force collapse :: EG_createQuadGroup is %d \n ",stat);
      printQuadGroup (qm -> mesh, qg);
  }
  for ( swap = 0 ; swap < 6; swap++ ) if ( qg.verts[swap] == v0 ) break;
  stat    = EG_swappingOperation (qm, qg, swap, &j );
  if ( stat == EGADS_SUCCESS && j == 1 ) *activity = 0;
  else {
#ifdef DEBUG
      printf(" In EG swapDoubleCollapse: I failed to collapse after swapping! s = %d act = %d \n ", stat, j );
#endif
  }
  return stat;
}


static int
EG_swapDoubleSplit (meshMap *qm, quadGroup qg, int *activity ) {
  int  i5, q, i55, i3, val3, i0, v30[2], i,  stat, adj[2], q0[2], v0;
  quadGroup sqg;
  *activity  = 0;
  for ( q = i3 = 0 ; i3 < 6; i3++)
    if ( qg.vals[i3] == 3 && i3 %3 != 0 ) {
	q = 1;
	break;
    }
  if ( q == 0 ) return EGADS_SUCCESS;
  i5  = (i3 + 3)% 6;
  if ( qg.vals[i5] < 5 ) return EGADS_SUCCESS;
  i55 = (i5 + 1)%6;
  if ( qg.vals[i55] < 5 ) {
      i55 = (i5 + 5)%6 ;
      if ( qg.vals[i55] < 5 ) return EGADS_SUCCESS;
  }
  q = 0;
  if ( EG_quadVertIdx  (qm -> mesh, qg.q[q], qg.verts[i5]) < 0 ) q = 1;
  stat = EG_adjQtoPair (qm -> mesh, qg.q[q], qg.verts[i5], qg.verts[i55], adj );
  if ( stat != EGADS_SUCCESS || adj[1] == -1 ) {
      if ( stat != EGADS_SUCCESS ) printf("In EG_swapDoubleSplit: adjacent to pair %d --> !!\n ", stat);
      return stat;
  }
  if ( (i3 + 1) % 6 == 0 || (i3 + 5) % 6 == 0 ) i0 = 0;
  else                                          i0 = 3;
  v30[0]     = qg.verts[ i0];
  v30[1]     = qg.verts[(i0 + 3)%6];
  q0[0]      = qg.q[0];
  q0[1]      = qg.q[1];
  val3       = qg.verts[i3];
  stat       = EG_createQuadGroup (qm -> mesh, &qg, qg.q[q], adj[1] );
  if ( stat != EGADS_SUCCESS ) {
      printf("Inside EG_swapDoubleSplit: before swapping EG_createQuadGroup stat %d\n ", stat );
      printQuadGroup (qm -> mesh, qg );
      return stat;
  }
  if ( validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0 ) return EGADS_SUCCESS;
  for ( i0 = 0 ; i0 < 6; i0++) if ( qg.verts[i0] == v30[0] ) break;
  v0         = qg.verts[0];
  stat       = EG_swappingOperation (qm, qg, i0, &(*activity) );
  if ( stat != EGADS_SUCCESS || *activity == 0 ) {
#ifdef DEBUG
      printf(" EG_swapDoubleSplit error at swap: %d activity %d!!\n ", stat, *activity );
#endif
      return stat;
  }
  i = 0 ; if ( EG_quadVertIdx (qm -> mesh, q0[0], val3) < 0 ) i = 1;
  stat       = EG_adjQtoPair (qm -> mesh, q0[i], v30[0], v30[1], adj );
  if ( stat != EGADS_SUCCESS || adj[1] == -1 ) {
      if ( stat != EGADS_SUCCESS ) printf(" EG_swapDoubleSplit after swapping adjacent to pair %d !!\n ", stat);
      return stat;
  }
  stat = EG_createQuadGroup (qm -> mesh, &sqg,  q0[i], adj[1]);
  if ( stat != EGADS_SUCCESS ) {
      printf("Inside EG_swapDoubleSplit: before splitting EG_createQuadGroup stat %d\n ", stat );
      printQuadGroup (qm -> mesh, qg );
      return stat;
  }
  stat = EG_doubleSplit (qm, sqg, 0, &i);
  if ( i == 0 ) {
      stat = EG_createQuadGroup (qm -> mesh, &qg, qg.q[0], qg.q[1]);
      for ( i = 0 ; i < 6; i++ ) if ( qg.verts[i] == v0 ) break;
      stat = EG_swappingOperation (qm, qg, i, &q );
      if ( stat == EGADS_SUCCESS && q == 1 ) *activity = 0;
  }
  return stat;
}


static int EG_transferValences ( meshMap *qm, int *qID, int try5533, int *transfering, int *activity  ) {
  int i, j, swap = 0, stat, min;
  int links[5], qAux[2];
  quadGroup qg;
  *activity   = 0;
  if ( qID[0] < 0 || qm -> mesh -> qIdx [ 4 * ( qID[0] - 1)]  == -2 ) {
#ifdef DEBUG
      printf(" EG_transferValences for quad %d is bad quad %d \n ", qID[0], checkQuad ( qm -> mesh, qID[0] ) );
#endif
      return EGADS_SUCCESS;
  }
  if (*transfering == 0 ) {
      stat          = EG_cleanQuad ( qm, qID[0], 1, 0, try5533 , &(*activity) );
      if ( stat    != EGADS_SUCCESS || *activity > 0 ) {
	  qID[0]    = -1;
	  if ( stat != EGADS_SUCCESS)
	    printf("EG_transferValences: EG_cleanQuad %d --> %d!!\n", qID[0], stat);
	  return stat ;
      }
  }
  if ( qID[0] == -1 || qm -> mesh -> qIdx [ 4 * ( qID[0]- 1) ] == -2 ||
      EG_nValenceCount (qm -> mesh, qID[0], 4) == 4 ) return EGADS_SUCCESS;
  for ( j = 0 ; j < 4; j++ ) {
      i       = qm -> mesh -> qAdj[ 4 * ( qID[0] - 1 ) + j];
      if ( i == -1 || i == qID[1] ) continue;
      if ( *transfering == 0 ) {
	  stat           = EG_cleanQuad ( qm, i, 1, 0,  try5533 , &(*activity) );
	  if ( stat     != EGADS_SUCCESS || *activity > 0 ) {
	      if ( stat != EGADS_SUCCESS) printf(" EG_TransferValence EG_cleanQuad --> %d !!\n ", stat );
	      qID[0]     = -1;
	      return stat;
	  }
      }
      if ( qm -> mesh -> qAdj[ 4 * ( qID[0] - 1 ) + j] != i ) continue;
      stat       = EG_createQuadGroup  (qm -> mesh, &qg, qID[0], i );
      if ( stat != EGADS_SUCCESS  ) {
	  printf(" Inside EG_transferValences EG_createQuadGroup %d !!\n", stat );
	  printQuadGroup (qm -> mesh, qg );
	  return stat;
      }
      if ( *transfering == 0 && qg.vals[0] * qg.vals[3] == 15 && qID[1] != -1  ) {
	  j = 0;
	  if ( qg.q[j] == qID[0] ) j = 1;
	  qID[0] = qg.q[j];
	  qID[1] = qg.q[(j + 1)%2];
	  if ( EG_quadIsBoundary(qm -> mesh, qg.q[1] ) == 1 && qm -> extraQuads >= 0 ) {
	      stat       = EG_collapse (qm, qID[0], &( *activity ), 1 );
	      if ( stat != EGADS_SUCCESS )
		printf("EG_transferValences: forceCollapse gave %d !!\n", stat);
	      qID[0] = qID[1];
	      qID[1] = -1;
	      return stat;
	  }
	  return EG_transferValences (qm, qID, 0, &(*transfering),  &(*activity) );
      }
      if (validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0 ) continue;
      if      (qg.vals[0] * qg.vals[3] >= 20 ) min = 12;
      else if ( try5533 == 1 && *transfering == 0  ) {
	  if      (qg.vals[0] * qg.vals[3] >= 25) min = 16;
	  else if (qg.vals[0] * qg.vals[3] == 16) min =  9;
      } else continue;
      if (      qg.vals[1] * qg.vals[4] <= min ) swap = 1;
      else if ( qg.vals[2] * qg.vals[5] <= min ) swap = 2;
      else continue;
      stat = EG_swappingOperation (qm, qg, swap, &(*activity) );
      if (  stat != EGADS_SUCCESS || *activity == 0 ) {
	  if ( stat != EGADS_SUCCESS)
	    printf(" EG_transferValences EG_swappingOperation stat %d activity %d !!\n", stat, *activity);
#ifdef DEBUG
	  printf(" EG_transferValences EG_swappingOperation stat %d activity %d !!\n", stat, *activity);
#endif
	  return stat;
      }
      i         = 0;
      if ( min == 9 || min == 16) *transfering = 1;
      else                        *transfering = 0;
      if (qID[1] > 0 && qm -> mesh -> qIdx [ 4 * ( qID[1] - 1)] != -1 ) {
	  if      ( qg.q[0] == qID[1] ) i = 1;
	  else if ( qg.q[1] == qID[1] ) i = 0;
	  else {
	      EG_commonVerts ( qm -> mesh, qg.q[0], qID[1], links);
	      if ( links[0] > 0 ) {
		  if ( links[0] == 1 ) {
		      EG_commonVerts ( qm -> mesh, qg.q[1], qID[1], links);
		      if ( links[0] == 0 ) i = 1;
		  }
		  else i = 1;
	      }
	  }
      }
      else if(EG_nValenceCount (qm -> mesh,qg.q[0],5) == 2 &&
	  EG_nValenceCount (qm -> mesh,qg.q[1],5) >= 1 &&
	  EG_nValenceCount (qm -> mesh,qg.q[1],3) == 1) i = 1;
      else if(EG_nValenceCount(qm -> mesh,qg.q[1],4) <
	  EG_nValenceCount(qm -> mesh,qg.q[0],4) ) i = 1;
      qID[0] = qg.q[  i          ];
      qID[1] = qg.q[( i + 1 ) %2 ];
      if (*transfering == 1 ) {
	  for ( j = 0 ; j < 2; j++ ) {
	      i          = 1;
	      qAux[0]    = qID[ j       ];
	      qAux[1]    = qID[(j + 1)%2];
	      stat       = EG_transferValences ( qm, qAux, 0, &i, &min);
	      if ( stat != EGADS_SUCCESS ) {
		  printf(" EG_transferValences: separating valences after forcing %d!!\n", stat);
		  return stat;
	      }
	  }
	  *transfering = 0;
	  qID[0]       = qAux[0];
	  qID[1]       = qAux[1];
      }
      if ( *activity > 0 ) break;
  }
  return stat;
}


static int EG_basicOperation (meshMap *qm, int qID, int type, int *activity ) {
  int stat = EGADS_SUCCESS;
  //*activity = 0;
  switch (type ) {
    case SWAP:
      stat  = EG_swap ( qm, qID, & (*activity));
      break;
    case COLLAPSE:
      stat  = EG_collapse ( qm, qID, & (*activity), 0);
      break;
    case SPLIT:
      stat  = EG_split ( qm, qID, & (*activity));
      break;
  }
  if ( stat != EGADS_SUCCESS && *activity == 0 ) return EGADS_SUCCESS;
  return stat;
}

static int EG_composeOperation (meshMap *qm, quadGroup qg, int type, int forcing, int *activity ) {
  int stat = EGADS_SUCCESS;
  switch (type ) {
    case DOUBLESWAP:
      stat  = EG_doubleSwap ( qm, qg, forcing, & (*activity));
      break;
    case SWAPCOLLAPSE:
      stat     = EG_swapCollapse ( qm, qg,forcing, & (*activity));
      break;
    case DOUBLECOLLAPSE:
      stat  = EG_doubleCollapse ( qm, qg, forcing,  & (*activity));
      break;
    case SWAPDOUBLECOLLAPSE:
      stat  = EG_swapDoubleCollapse ( qm, qg, & (*activity));
      break;
    case SWAPSPLIT:
      stat     = EG_swapSplit ( qm, qg, forcing, & (*activity));
      break;
    case DOUBLESPLIT:
      stat  = EG_doubleSplit ( qm, qg, forcing, & (*activity));
      break;
    case SWAPDOUBLESPLIT:
      stat  = EG_swapDoubleSplit ( qm, qg, & (*activity));
      break;
  }
  if ( stat != EGADS_SUCCESS )
    printf(" OPERATION %d gave %d and activity %d \n ", type, stat, *activity );
  return stat;
}

static int EG_cleanQuad (meshMap *qm, int qID, int useAdj, int transfer, int forcing, int *activity ) {
  int stat, i, q, qadj, act = 0;
  int opBasic[3] = {COLLAPSE, SWAP, SPLIT};
  int opComp [7] = {SWAPCOLLAPSE, DOUBLECOLLAPSE, SWAPDOUBLECOLLAPSE, DOUBLESWAP,
      SWAPSPLIT   , DOUBLESPLIT,     SWAPDOUBLESPLIT};
  quadGroup qg;
  *activity      = 0;
  if ( qm -> mesh -> qIdx [ 4 * (qID - 1 ) ] == -2  ) {
#ifdef DEBUG
      printf(" EG_cleanQuad %d stat %d\n", qID, stat );
#endif
      return EGADS_SUCCESS;
  }
  if ( EG_nValenceCount (qm -> mesh, qID, 4) == 4 ) return EGADS_SUCCESS;
#ifdef DEBUG
  printf(" EG CLEAN QUAD %d  transfer %d forcing %d\n ", qID , transfer, forcing );
#endif
  if ( transfer == 0 ) {
      if ( EG_nValenceCount (qm -> mesh, qID, 6) > 0 || qm -> extraQuads < 0 ) swapInt( &opBasic[2], &opBasic[0]);
      for ( i = 0 ; i < 3; i++ ) {
	  stat       = EG_basicOperation(qm, qID, opBasic[i], &act );
	  *activity += act;
	  if ( stat != EGADS_SUCCESS ) return stat;
	  if ( act > 0 || qID < 0 || qm -> mesh -> qIdx [ 4 * (qID - 1 ) ] == -2 ) return EGADS_SUCCESS;
      }
  }
  if ( useAdj == 0 )  return EGADS_SUCCESS;
  if ( qm -> extraQuads < 0 )
    for ( i = 0 ; i < 3; i++) swapInt( &opComp[i], &opComp[i+4]);
  if ( qID < 0 || qm -> mesh -> qIdx [ 4 * (qID - 1 ) ] == -2) {
#ifdef DEBUG
      printf(" EG_cleanQuad %d checkQuad is stat %d!!\n", qID, stat );
#endif
      return EGADS_SUCCESS;
  }
  if ( forcing == 1 && EG_nValenceCount (qm -> mesh, qID, 3) > 1 ) {
#ifdef DEBUG
      printf(" FORCE COLLAPSE AROUND QUAD %d \n ", qID  );
#endif
      stat = EG_collapse (qm, qID, &(*activity), 1 );
      if ( stat != EGADS_SUCCESS ) return stat;
      if ( *activity > 0 ) {
	  --qm -> extraQuads;
	  return EGADS_SUCCESS;
      }
  }
  if ( transfer == 0 ) {
      for ( q        = 0 ; q < 4; q++) {
	  qadj       = qm -> mesh -> qAdj [ 4 * ( qID -1 ) + q];
	  if ( qadj == -1 ) continue;
	  stat       = EG_cleanQuad (qm, qadj, 0, 0, 0, &act);
	  *activity += act;
	  if ( stat != EGADS_SUCCESS ) {
	      printMesh(qm, NULL);
	      printf(" EG_cleanQuad restore mesh after adjacent stat --> %d !!\n", stat );
	      return stat;
	  }
	  if ( act > 0 || qID < 0 || qm -> mesh -> qIdx [ 4 * (qID - 1 ) ] == -2) return EGADS_SUCCESS;
      }
  }
  for ( i  = 0 ; i < 7; i++ ) {
      if  ( transfer == 1 && opComp[i] == DOUBLESWAP ) continue;
      for ( q  = 0 ; q < 4; q++) {
	  if ( qID < 0 || qm -> mesh -> qIdx [ 4 * (qID - 1 ) ] == -2 ) return EGADS_SUCCESS;
	  qadj = qm -> mesh -> qAdj [ 4 * ( qID -1 ) + q];
	  if ( qadj == -1 ) continue;
	  stat = EG_createQuadGroup (qm -> mesh, &qg, qID, qadj);
	  if ( stat != EGADS_SUCCESS ) {
	      printf("Inside EG_cleanQuad: EG_createQuadGroup stat --> %d !!\n ", stat );
	      printQuadGroup (qm -> mesh, qg );
	      return stat;
	  }
	  stat       = EG_composeOperation(qm, qg, opComp[i], forcing, &act );
	  *activity += act;
	  if ( stat != EGADS_SUCCESS ) {
	      printMesh(qm, NULL);
	      printf("In EG_cleanQuads restor mesh after EG_composeOperation stat -->  %d !!\n ", stat );
	  }
	  if ( act > 0 || qID < 0 || qm -> mesh -> qIdx [ 4 * (qID - 1 ) ] == -2 ) return EGADS_SUCCESS;
      }
  }
  return EGADS_SUCCESS;
}

static void EG_destroymeshMap(meshMap ***qm, int nfaces ) {
  int i ;
  if ( *qm  == NULL ) return ;
  for ( i = 0 ; i < nfaces; ++i) {
      if ( ( *qm)[i] ) {
	  EG_freeMeshData (&((*qm)[i] -> mesh      ));
	  EG_freeMeshData (&((*qm)[i] -> bestMesh  ));
	  EG_free((*qm)[i] -> vFix );
	  EG_free((*qm)[i]);
      }
  }
  EG_free (*qm );
}

static int
EG_adjQtoPair(meshData *mesh, int qID, int v1, int v2, int *adj) {
  int i, aux = -1;
  adj[0] = -1; adj[1] = -1;
  for ( i = 0 ; i < 4; ++i) {
      if ( mesh -> qIdx[4*(qID - 1) + i ] == v1 ) adj[0] = i;
      if ( mesh -> qIdx[4*(qID - 1) + i ] == v2 ) aux    = i;
      if ( aux != -1 && adj[0] != -1 ) break;
  }
  if ( aux == -1 || adj[0] == -1 ) return EGADS_SUCCESS;
  if      ( abs (adj[0] - aux ) == 3 ) adj[0] = 3;
  else if ( aux < adj[0]             ) adj[0] = aux;
  adj[1] = mesh -> qAdj[4*(qID - 1) + adj[ 0 ] ];
  return EGADS_SUCCESS;
}

static int
resizeQm( meshMap *qm) {
  int    stat, nV, vRem, nQ, qRem,  i, j, k, *vpiv = NULL, *qpiv = NULL;
  double  eval[18];
  nV = 0 ; vRem = 0 , nQ = 0 , qRem = 0;
  for ( i = 0 ; i < qm -> mesh ->totVerts; i++) {
      if(qm -> mesh -> vType[i] != -2 ) nV++;
      else vRem++;
  }
  for ( i = 0 ; i < qm -> mesh ->totQuads; i++) {
      if(qm -> mesh ->qIdx[4*i] != -2 ) nQ++;
      else qRem++;
  }
  if ( vRem != qRem ) {
      if ( qRem > 0 ) {
	  printf(" In resizeQm: I have %d removed vertices and %d quads!! they should match!!!!!\n ", vRem, qRem);
	  return EGADS_INDEXERR;
      }
  }
  if ( vRem == 0 )return EGADS_SUCCESS;
  stat = EG_restoreMeshMap (qm, qm -> bestMesh, qm -> mesh );
  if ( stat != EGADS_SUCCESS) return stat;
  vpiv = (int*) EG_alloc(qm -> mesh -> totVerts * sizeof(int));
  qpiv = (int*) EG_alloc(qm -> mesh -> totQuads * sizeof(int));
  if (vpiv == NULL || qpiv == NULL) return EGADS_MALLOC;
  for ( j = i = 0 ; i < qm -> bestMesh -> totQuads; i++) {
      if(qm -> mesh -> qIdx[4*i] != -2 ) {
	  qpiv[i] = j;
	  ++j;
      } else qpiv[i] = -2;
  }
  for ( j = i = 0 ; i < qm -> bestMesh -> totVerts; i++) {
      if(qm -> mesh -> vType[i] != -2 ) {
	  vpiv[i] = j;
	  ++j;
      } else vpiv[i] = -2;
  }
  nV   = qm -> mesh -> totVerts - vRem;
  nQ   = qm -> mesh -> totQuads - qRem;
  stat = EG_allocMeshData ( &qm -> mesh, nQ, nV );
  if ( stat != EGADS_SUCCESS ) {
      EG_free (vpiv);
      EG_free (qpiv);
      printf("Inside resizeQM: EG_allocMeshData %d !!\n", stat );
      return stat;
  }
  for ( i = 0 ; i < qm -> bestMesh -> totQuads; i++){
      if ( qpiv[i] == -2) continue;
      for ( k = 0 ; k < 4; ++k) {
	  qm -> mesh -> qIdx[4*qpiv[i] + k ] = vpiv[ qm -> bestMesh -> qIdx[4*i + k] - 1] + 1;
	  if ( qm -> bestMesh -> qAdj[4*i + k] == -1 ) qm -> mesh -> qAdj[4*qpiv[i] + k ] = - 1;
	  else qm -> mesh       -> qAdj[4*qpiv[i] + k ] = qpiv[ qm -> bestMesh -> qAdj[4*i + k] - 1] + 1;
      }
  }
  for ( i = 0 ; i < qm -> bestMesh -> totVerts; i++) {
      if ( vpiv[i] == -2 ) continue;
      j   = vpiv[i];
      qm -> mesh -> vType[j      ] = qm -> bestMesh -> vType[  i    ];
      qm -> mesh -> uvs  [2*j    ] = qm -> bestMesh -> uvs  [2*i    ];
      qm -> mesh -> uvs  [2*j + 1] = qm -> bestMesh -> uvs  [2*i + 1];
      stat                         = EG_evaluate(qm -> face, &qm -> mesh -> uvs[2*j], eval);
      qm -> mesh -> xyzs[3*j    ]  = eval[0];
      qm -> mesh -> xyzs[3*j + 1]  = eval[1];
      qm -> mesh -> xyzs[3*j + 2]  = eval[2];
      qm -> mesh -> valence[j][0]  = qm -> bestMesh -> valence[i][0];
      qm -> mesh -> valence[j][1]  = qpiv[qm -> bestMesh -> valence[i][1] - 1] + 1;
      for ( k = 0 ; k < qm -> bestMesh -> valence[i][0]; ++k)
	qm -> mesh -> valence[j][2 + k] = vpiv[qm -> bestMesh -> valence[i][2 + k] -1] + 1;
  }
  EG_free(vpiv);
  EG_free(qpiv);
  return checkMesh (qm);
}


static void updateVertex ( meshMap *qm, int vID, double *uv ) {
  int i;
  double eval[18];
  i  = EG_evaluate ( qm -> face, uv, eval);
  if ( i != EGADS_SUCCESS ) return ;
  qm -> mesh -> uvs [ 2 * ( vID -1 )    ] = uv  [0];
  qm -> mesh -> uvs [ 2 * ( vID -1 ) + 1] = uv  [1];
  qm -> mesh -> xyzs[ 3 * ( vID -1 )    ] = eval[0];
  qm -> mesh -> xyzs[ 3 * ( vID -1 ) + 1] = eval[1];
  qm -> mesh -> xyzs[ 3 * ( vID -1 ) + 2] = eval[2];
}



static int vertexLinksToBounds (meshMap *qm, int vID, int *nb ) {
  int i, j;
  for ( nb[0] = i = 0 ; i < qm -> mesh -> valence[ vID -1 ][0]; i++ ) {
      j  =  qm -> mesh -> valence[ vID -1 ][2 + i] - 1;
      if ( qm -> mesh -> vType[ j ] >= 4 ) nb[ ++nb[0] ] = j + 1;
  }
  for ( i = 0 ; i < qm -> mesh -> valence[ vID -1 ][0]; i++ ) {
      j   =  qm -> mesh -> valence[ vID -1 ][2 + i] - 1;
      if ( qm -> mesh -> vType[ j ] >= 0 && qm -> mesh -> vType [ j ] < 4 ) nb[ ++nb[0] ] = j + 1;
  }
  return nb[0];
}

static int
EG_buildStar(meshData *mesh, vStar **star, int vID ) {
  int  i = 0 , id0 = -1, q = 0, auxV, auxQ,  v = 0, quadID, prevQuad, it = 0, it2 = 0;
  int adj[2], *vertex = NULL, *quads = NULL;
  int qLoop[8] = {0, 1, 2, 3, 0, 1, 2, 3};
#ifdef DEBUG
  i       = checkVertex ( mesh , vID ) ;
  if ( i != EGADS_SUCCESS ) {
      printf(" Vertex %d is %d \n ", vID, i);
      return i;
  }
#endif
  vertex       = (int * ) EG_alloc ( mesh -> totVerts * sizeof ( int ) );
  quads        = (int * ) EG_alloc ( mesh -> totQuads * sizeof ( int ) );
  if ( vertex == NULL || quads == NULL ) {
      printf("EG_buildStar MALLOC at quads & verts!!\n ");
      return EGADS_MALLOC;
  }
  // quads are -1 bias
  quadID  = mesh -> valence[ vID - 1 ][1] - 1;
#ifdef DEBUG
  i       = checkQuad ( mesh, quadID + 1);
  if ( i != EGADS_SUCCESS){
      printf(" EG_buildStar at vertex %d has associated a bad quad: %d --> %d!!\n", vID, quadID + 1, i );
      printQuadSpecs ( mesh, quadID + 1);
      EG_free ( vertex );
      EG_free ( quads  );
      return i;
  }
#endif
  vertex [v++] = vID;
  it           = 0;
  do {
#ifdef DEBUG
      i = checkQuad ( mesh, quadID + 1);
      if ( i != EGADS_SUCCESS ) {
	  printf(" In EG_buildStar quad %d is bad quad --> %d!!\n",  quadID + 1, i);
	  printQuadSpecs( mesh, quadID + 1);
	  EG_free (quads );
	  EG_free (vertex);
	  return i;
      }
#endif
      id0 = EG_quadVertIdx ( mesh, quadID + 1, vID );
      if ( id0 < 0 ) {
	  printf(" In EG_buildStar id for Vert %d in Quad %d is %d !!\n", vID, quadID + 1, id0);
#ifdef DEBUG
	  printQuadSpecs( mesh, quadID + 1);
#endif
	  EG_free ( vertex ) ;
	  EG_free ( quads ) ;
	  return EGADS_INDEXERR;
      }
      for ( i = 1 ; i <= 2; ++i) vertex[v++] = mesh -> qIdx[ 4 * quadID   + qLoop[id0 + i ] ] ;
      quads[q++] = quadID + 1;
      prevQuad   = quadID;
      quadID     = mesh -> qAdj[ 4 * prevQuad + qLoop[id0 + 3 ] ] - 1;
      if ( quadID  < 0 ) { //make a "ghost" quad
	  auxQ        = prevQuad;
	  vertex[v++] = mesh -> qIdx[4*auxQ + qLoop[id0 + 3]] ;
	  auxV        = mesh -> qIdx[4*auxQ + qLoop[id0 + 1]] ;
	  it2 = 0;
	  do
	    {
	      i = EG_adjQtoPair(mesh, auxQ + 1, vID, auxV, adj );
	      if ( adj[1] == -1 || i != EGADS_SUCCESS ) break;
	      auxQ = adj[1] - 1;
	      i    = EG_quadVertIdx ( mesh, auxQ + 1, vID ) ;
	      if ( i < 0 ){
		  printf(" In buildStar vertex Id %d in quad %d is %d\n", vID, quadID + 1, i);
#ifdef DEBUG
		  printQuadSpecs ( mesh, quadID + 1);
#endif
		  EG_free ( vertex ) ;
		  EG_free ( quads ) ;
		  return EGADS_INDEXERR;
	      }
	      auxV = mesh -> qIdx [ 4 * auxQ + qLoop [ i + 1 ] ];
	      it2++;
	      if ( it2 > 200 ) {
		  printf(" stuck in interior loop of build star!!!!!!!!\n");
		  EG_free (quads);
		  EG_free (vertex);
		  return EGADS_RANGERR;
	      }
	    }
	  while ( adj[1] != - 1 );
	  quads  [q++]    = -1 ;
	  vertex [v++]    = -1 ;
	  quadID          = auxQ;
      }
      if ( quadID < 0 ) {
	  printf(" I am stuck in build star. Pointing a NULL quad \n");
	  EG_free ( vertex ) ;
	  EG_free ( quads ) ;
	  return EGADS_INDEXERR;
      }
      it++;
      if ( it > 200 ) {
	  printf(" EG_buildStar:: stuck in outer loop of build star!!!!!!!!\n");
	  EG_free ( vertex ) ;
	  EG_free ( quads ) ;
	  return EGADS_RANGERR;
      }
  }
  while ( quadID + 1 != quads[0] );
  if ( *star != NULL ) EG_freeStar (& (*star ) );
  *star               = (vStar *)EG_alloc ( sizeof(vStar));
  if ((*star) == NULL ) {
      EG_free ( vertex ) ;
      EG_free ( quads ) ;
      return EGADS_MALLOC;
  }
  (*star) -> nQ    = q;
  (*star) -> nV    = v;
  (*star) -> verts = (int*)  EG_alloc (    v * sizeof(int  ));
  (*star) -> quads = (int*)  EG_alloc (    q * sizeof(int  ));
  (*star) -> idxV  = (int*)  EG_alloc (2 * v * sizeof(int  ));
  (*star) -> idxQ  = (int*)  EG_alloc (2 * q * sizeof(int  ));
  if ((*star) -> verts == NULL || (*star) -> quads == NULL ||
      (*star) -> idxV  == NULL || (*star) -> idxQ  == NULL )
    {
      EG_free((*star  ));
      EG_free ( vertex) ;
      EG_free ( quads ) ;
      return EGADS_MALLOC;
    }
  for ( i = 0 ; i < q; ++i) {
      (*star) -> quads[i    ] = quads [i];
      (*star) -> idxQ [i    ] = i;
      (*star) -> idxQ [q + i] = i;
  }
  for ( i = 0 ; i < v; ++i) {
      (*star) -> verts[i    ] = vertex [i];
      (*star) -> idxV [i    ] = i    ;
      (*star) -> idxV [v + i] = i + 1;
  }
  EG_free ( vertex ) ;
  EG_free ( quads ) ;
  return EGADS_SUCCESS;
}

static int
EG_collapse (meshMap *qm, int qID, int *activity, int forcing )  {
  int i, vC, v5[5], v3[5], qb[4], stat, i3, j, links[2], val[2];
  *activity     = 0 ;
  if ( qm -> mesh -> qIdx [ 4 * (qID - 1 ) ] == -2 ) {
#ifdef DEBUG
      printf(" EG_collapse for quad %d is %d \n ", qID, checkQuad (qm -> mesh, qID) );
#endif
      return EGADS_SUCCESS;
  }
  for ( v3[0] = v5[0] = i = 0 ; i < 4; i++ ) {
      vC =  qm -> mesh -> qIdx [ 4 *( qID - 1) + i];
      if      ( getValence (qm ->mesh, vC) == 3) v3[++v3[0]] = i;
      else if ( getValence (qm ->mesh, vC) == 5) v5[++v5[0]] = i;
      qb[i] = vC ;
  }
  if ( v3[0] + v5[0] <= 2 || v3[0] * v5[0] == 0 ) return EGADS_SUCCESS;
  if ( v5[0] == 1 || ( v5[0] == 2 && abs ( v5[1] - v5[2] ) % 2 == 0 ) ) {
      vC = qm -> mesh -> qIdx [ 4 * ( qID -1 ) + (v5[1] + 1 ) %4 ];
      if ( validCollapse ( qm, qID, vC ) == 1 ) {
	  stat = EG_mergeVertices(qm, qID, vC, &(*activity));
	  if ( stat != EGADS_SUCCESS ) {
	      printf(" EG_collapse :: EG_mergeVertices %d activity %d \n", stat, *activity );
	      return stat ;
	  }
      }
  }
  if ( *activity == 0 && forcing == 1 && EG_nValenceCount ( qm -> mesh, qID, 3 ) > 0 ) {
      j = 0; i3 = 0;
      while ( j < 2 ) {
	  vC = qm -> mesh -> qIdx [ 4 * ( qID - 1 ) + i3 ];
	  if ( qm -> mesh -> vType [vC -1 ] == -1 && qm -> mesh -> valence[vC -1][0] == 3 ) {
	      links[0] = qm -> mesh -> qIdx [ 4 * ( qID - 1 ) + (i3 + 1 ) % 4 ];
	      links[1] = qm -> mesh -> qIdx [ 4 * ( qID - 1 ) + (i3 + 3 ) % 4 ];
	      val  [0] = getValence ( qm -> mesh, links[0] );
	      val  [1] = getValence ( qm -> mesh, links[1] );
	      if ( (val[0] >= 5 || val[1] >= 5 || j == 1) && validCollapse (qm, qID, vC ) == 1) {
		  stat = EG_mergeVertices ( qm, qID, vC, &(*activity) );
		  if ( stat != EGADS_SUCCESS ) {
		      printf(" EG_collapse forcing :: EG_mergeVertices %d activity %d \n", stat, *activity );
		      return stat ;
		  }
		  if ( *activity > 0 ) break;
	      }
	  }
	  i3++;
	  if ( i3 == 4 ) {
	      j++;
	      i3 = 0;
	  }
      }

  }
  for ( i = 0 ; i < 4; i++ ) {
      vC  = qm -> mesh -> vType [ qb[i] -1 ];
      if ( vC == -1 && qm -> mesh -> valence [ qb[i] - 1][0] == 2 ) {
	  stat = EG_mergeVertices(qm, qm -> mesh -> valence [ qb[i] - 1][1], qb[i], &(*activity));
	  if ( stat != EGADS_SUCCESS ) {
	      printf("EG_collapse doublet EG_mergeVertices %d activity %d \n ", stat, *activity );
	      return stat ;
	  }
      }
  }
  return EGADS_SUCCESS;
}





static int
EG_mergeVertices (meshMap *qm, int qC, int centre, int *activity  ) {
  int  stat, i, j, q,  adjq, adjPair[2], auxQ, oldQ[8], nb[10], nq;
  int piv[4]  = {1, 0, 3, 2} ;
  double uv[2], uvxyz[10];
  Quad *quad  = NULL;
  vStar *star = NULL;
  nq          = 1;
  stat        = EG_backupQuads (qm, &nq, &qC, &quad );
  printQuadCoords (qm, qC );
  if ( stat != EGADS_SUCCESS ) {
      printf(" In EG_mergeVertices stat in EG_backpQuads %d !!\n ", stat );
      return stat;
  }
  *activity   = 0;
  i           = EG_quadVertIdx  ( qm -> mesh, qC, centre );
  if ( qm -> mesh -> vType [ centre - 1] != -1 ) {
      centre = qm -> mesh -> qIdx [ 4 * ( qC - 1) + ( i + 2 ) %4];
      i      = EG_quadVertIdx  ( qm -> mesh, qC, centre );
  }
  if ( qm -> mesh -> vType [ centre - 1] != -1 || i < 0) {
#ifdef DEBUG
      printf(" EG_mergeVertices collapse through %d with index at %d = %d !!\n ", centre, qC, i );
      printQuadSpecs ( qm -> mesh, qC );
#endif
      return EGADS_INDEXERR;
  }
  for ( q = 0; q < 4; ++q ) {
      oldQ[q    ] = qm -> mesh -> qIdx[4 * (qC - 1) + (q + i)%4 ];
      oldQ[q + 4] = qm -> mesh -> qAdj[4 * (qC - 1) + (q + i)%4 ];
  }
  for ( j  = 0 ; j < 3; j++ ) {
      if ( j < 2 )  {
	  uvxyz[j    ] = qm -> mesh -> uvs  [ 2 * (oldQ[0] - 1) + j];
	  uvxyz[5 + j] = qm -> mesh -> uvs  [ 2 * (oldQ[2] - 1) + j];
      }
      uvxyz[2 + j]     = qm -> mesh -> xyzs [ 3 * (oldQ[0] - 1) + j];
      uvxyz[7 + j]     = qm -> mesh -> xyzs [ 3 * (oldQ[2] - 1) + j];
  }
  // save old quads vertices and adjacents to update map correctly
  nb[0] = vertexLinksToBounds (qm , centre , nb );
  if (      qm -> mesh -> valence [oldQ[0] - 1][0] !=  2 &&
      (     qm -> mesh -> vType   [oldQ[2] - 1]    == -1 &&
	  ( nb[0] == 0 || ( nb[0] > 0 && qm -> mesh -> vType [ nb[1] - 1] < 4 ) ) ) ) {
      uv[0]  = 0.0;
      uv[1]  = 0.0;
      for (i = 0; i < 4; i++ ) {
	  j  = qm -> mesh -> qIdx [4 * ( qC - 1 ) + i ] - 1;
	  uv[0] += 0.25 * qm -> mesh -> uvs [2 * j    ];
	  uv[1] += 0.25 * qm -> mesh -> uvs [2 * j + 1];
      }
      updateVertex (qm, oldQ[2], uv);
  }
  stat       = EG_buildStar ( qm -> mesh, &star, centre );
  if ( stat != EGADS_SUCCESS || star == NULL ) {
      printf(" Attention:: building star at %d --> %d !!\n ", centre, stat );
      return stat;
  }
  for ( i    = 0 ; i < 4; i++ ) {
      q      = oldQ [ 4 + i ];
      if ( q < 0 ) continue;
      adjq = oldQ [ 4 + piv[ i ] ]; // - 1 bias
      stat = EG_adjQtoPair ( qm -> mesh, q , oldQ[ i ], oldQ[ ( i + 1 ) % 4  ] , adjPair );
      if ( stat != EGADS_SUCCESS || adjPair[1] != qC ) {
	  EG_freeStar (&star);
	  printf(" STAT %d  adj %d to pair %d %d is %d \n ", stat, q ,oldQ[ i ], oldQ[ ( i + 1 ) % 4  ] , adjPair[1]);
	  return  EGADS_INDEXERR;
      }
      qm -> mesh -> qAdj [ 4 * (q - 1 ) + adjPair[0] ] = adjq;
  }
  // Eliminate vertex p[0] from all the quads and its valences
  for ( i = 0 ; i < star -> nQ; ++i) {
      q = star -> quads[i];
      if      ( q == -1 ) continue; // ghost quad
      else if ( q == qC ) {  // eliminate quad qC
	  for ( j = 0 ; j < 4; ++ j) {
	      qm -> mesh -> qIdx[ 4 * (qC - 1) + j ] = -2;
	      qm -> mesh -> qAdj[ 4 * (qC - 1) + j ] = -2;
	  }
      } else {
	  for ( j = 0 ; j < 4; ++ j)
	    if ( qm -> mesh -> qIdx[ 4 * (q - 1)  + j] == oldQ[0])
	      qm -> mesh -> qIdx[4 * (q - 1) + j] = oldQ[2];
      }
  }
  // Point all the collapsed quad vertices to a valid quad
  for ( i = 1 ; i < 4; i++) {
      if ( qm -> mesh -> valence[ oldQ[i] -1 ][1] != qC ) continue;
      for ( q = 0 ; q < 4; q++ ) {
	  auxQ = oldQ[ 4 + q ];
	  if ( auxQ < 0 ) continue;
	  if ( EG_quadVertIdx ( qm -> mesh, auxQ, oldQ[i]) >= 0 ) {
	      qm -> mesh -> valence[ oldQ[i] - 1 ][1] = auxQ;
	      break;
	  }
      }
      if ( qm -> mesh -> valence[ oldQ[i] -1 ][1] == qC ) {
	  printf("VOOPP EG_mergeVertices vertex still points at collapsed quad %d !!\n", EGADS_INDEXERR );
	  EG_freeStar (&star);
	  return EGADS_INDEXERR ;
      }
  }
  // Set valences for merged stuff
  for ( i = 1 ; i < 4; i++) {
      stat = setValence ( qm -> mesh, oldQ[i] );
      if ( stat != EGADS_SUCCESS) {
	  printf( "LINKS Inside EG_mergeVertices stat in setValence %d = %d\n", i, stat );
	  EG_freeStar(&star);
	  return stat;
      }
  }
  // set valences to links
  for ( i = 0 ; i < star -> nQ; ++i) {
      j   = star -> verts [ 2 * i + 1 ];
      if ( j < 0 || qm -> mesh -> vType [ j - 1] == -2 ) continue; // ghost vertex (-1) or collapsed vertex (-2)
      stat       = setValence ( qm -> mesh, j );
      if ( stat != EGADS_SUCCESS) {
	  printf( "Inside EG_mergeVertices stat in setValence %d = %d\n", i, stat );
	  EG_freeStar(&star);
	  return stat;
      }
  }
  EG_freeStar (&star);
  // delete vertex vC2
  qm -> mesh -> vType   [  oldQ[0] - 1]               = -2; // -2 = removed
  qm -> mesh -> remQuads[++qm -> mesh -> remQuads[0]] = qC;
  qm -> mesh -> remVerts[++qm -> mesh -> remVerts[0]] = oldQ[0];
  // remove possible doublets from collapsing: look at links valences
  /*for  ( i = 0 ; i < 2; i++ ) {
      stat       = EG_removeDoublet(qm, oldQ[1 + 2 * i ], &j);
   *activity += j;
      if ( stat != EGADS_SUCCESS ) {
	  printf("In EG_mergeVertices I failed at removing doublet from %d -> %d !!\n ", oldQ[2 * i + 1] , stat);
	  return stat;
      }
  }*/
  *activity  = 1;
  if ( EG_makeValidMesh ( qm, 3, &oldQ[1], 0) == EGADS_SUCCESS ) {
      EG_free (quad );
      return EGADS_SUCCESS;
  }
#ifdef DEBUG
  printf(" COLLAPSE MESH FAILED %d (quad %d vertex %d )!\n", stat, qC, oldQ[0]);
  printMesh (qm , NULL );
#endif
  qm -> mesh -> vType   [  oldQ[0] -1]  = -1;
  stat = EG_restoreQuads ( qm, quad, nq );
  EG_free (quad);
  if ( stat != EGADS_SUCCESS ) {
      printf(" EG_restoreQuads stat %d \n", stat );
      return stat;
  }
  i = inList (qm -> mesh -> remQuads[0], &qm -> mesh -> remQuads[1], qC );
  for ( j = i; j < qm -> mesh ->remQuads[0] - 1; j++ ) {
      qm -> mesh -> remQuads[i] = qm -> mesh -> remQuads[i + 1];
      qm -> mesh -> remVerts[i] = qm -> mesh -> remVerts[i + 1];
  }
  --qm -> mesh -> remQuads[0];
  --qm -> mesh -> remVerts[0];
  for ( i = 0 ; i < 3; i++ ) {
      if ( i < 2 )  {
	  qm -> mesh -> uvs  [ 2 * (oldQ[0] -1 ) + i] = uvxyz[    i];
	  qm -> mesh -> uvs  [ 2 * (oldQ[2] -1 ) + i] = uvxyz[5 + i];
      }
      qm -> mesh -> xyzs [ 3 * (oldQ[0] - 1) + i] = uvxyz[2 + i];
      qm -> mesh -> xyzs [ 3 * (oldQ[2] - 1) + i] = uvxyz[7 + i];
  }
#ifdef DEBUG
  printMesh (qm , NULL );
#endif
  *activity = 0;
  return EGADS_SUCCESS;
}


static int EG_swap (meshMap *qm, int qIn, int *activity) {
  int   stat, q, swap = 0;
  quadGroup qg;
  *activity  = 0;
  if ( qm -> mesh -> qIdx [ 4 * (qIn - 1 ) ] == -2 ) {
#ifdef DEBUG
      printf(" EG_swap for quad %d is %d \n ", qIn, stat );
#endif
      return EGADS_SUCCESS;
  }
  for ( swap = q = 0 ; q < 4; q++ ) {
      qg.q[0]       = qIn;
      qg.q[1]       = qm -> mesh -> qAdj [ 4 * ( qIn - 1 ) + q ];
      if ( qg.q[1] == -1 ) continue;
      stat          = EG_createQuadGroup  ( qm -> mesh, &qg, qg.q[0], qg.q[1] );
      if ( stat    != EGADS_SUCCESS ) {
	  printf(" IN EG_swap -> EG_createQuadGroup  %d !! \n", stat );
	  printQuadGroup (qm ->mesh, qg);
	  return stat;
      }
      if      (qg.vals[0] <= 4 || validSwap ( qm -> mesh, qg.verts[0], qg.verts [3] ) != 1 ) continue;
      if      (qg.vals[1] * qg.vals[4] == 9 ) swap = 1;
      else if (qg.vals[2] * qg.vals[5] == 9 ) swap = 2;
      else if (qg.vals[3] >= 5 ||
	  ( qm -> mesh -> vType [ qg.verts[0] - 1] == 2 * qm -> mesh -> sizeQuads && qg.vals[0] > 5 ) ||
	  ( qm -> mesh -> vType [ qg.verts[3] - 1] == 2 * qm -> mesh -> sizeQuads && qg.vals[3] > 5 ) ){
	  if  (     qg.vals[1] * qg.vals[4] == 12 ) swap = 1;
	  else if ( qg.vals[2] * qg.vals[5] == 12 ) swap = 2;
      }
      if (swap != 0 ) break;
  }
  if (swap  == 0 ) return EGADS_SUCCESS;
  return EG_swappingOperation (qm, qg, swap, &(*activity) );
}


static int
EG_swappingOperation (meshMap *qm, quadGroup qg, int swap, int *activity ) {
  int stat, i0, i1, i, j, q, adj, idx[4], qID[2], adjQmap[6], piv[4], ori[4], area, nq;
  double angles[4];
  Quad *quad = NULL;
  swap       = swap %3;
  *activity  = 0;
  printf(" swapping opration thru %d \n", swap);
  if ( swap == 0 ) {
#ifdef DEBUG
      printf(" swapping throu 0-3 will result in the same pair!! \n ");
#endif
      return EGADS_INDEXERR;
  }
  qID[0]     = qg.q[0]; qID[1] = qg.q[1];
  nq         = 2;
  stat       = EG_backupQuads ( qm, &nq, qID, &quad);
  if ( stat != EGADS_SUCCESS ) {
      EG_free(quad );
      printf(" In EG_swappingOperation :: EG_backupQuad %d !!\n ", stat );
      return stat;
  }
  i0 = EG_quadVertIdx ( qm -> mesh, qID[0], qg.verts[0]); // centre
  i1 = EG_quadVertIdx ( qm -> mesh, qID[1], qg.verts[3]); // opposite
  // Modify Quads and get adj map
  qm -> mesh -> qIdx [ 4 * ( qID[0] - 1 )    ] = qg.verts[swap ];
  qm -> mesh -> qIdx [ 4 * ( qID[1] - 1 )    ] = qg.verts[swap ];
  for ( i = 0 ; i < 3; i++ ) {
      adjQmap[i    ] = qm -> mesh -> qAdj [ 4 * ( qID[0] - 1 ) + ( i + i0 )%4 ];
      adjQmap[i + 3] = qm -> mesh -> qAdj [ 4 * ( qID[1] - 1 ) + ( i + i1 )%4 ];
      qm -> mesh -> qIdx [ 4 * ( qID[0] - 1 ) + i + 1 ] = qg.verts[( swap + 3 + i )%6 ];
      qm -> mesh -> qIdx [ 4 * ( qID[1] - 1 ) + i + 1 ] = qg.verts[( swap + 1 + i )   ];
  }
  qm -> mesh -> qAdj [ 4 * ( qID[0] - 1 )    ] = qID[1];
  qm -> mesh -> qAdj [ 4 * ( qID[1] - 1 ) + 3] = qID[0];
  for ( i = 0 ; i < 3; i++ ) {
      adj = adjQmap[ (3 + i + swap) % 6 ] - 1;
      qm -> mesh -> qAdj [ 4 * ( qID[0] - 1 ) + i + 1] = adj + 1;
      if ( adj >= 0 ) {
	  for ( j = 0; j < 4; j++ ) if ( qm -> mesh -> qAdj [ 4 * adj + j ] == qID[1] )
	    qm -> mesh -> qAdj [ 4 * adj + j ] = qID[0];
      }
      adj = adjQmap[i + swap] -1;
      qm -> mesh -> qAdj [ 4 * ( qID[1] - 1 ) + i    ] = adj + 1;
      if ( adj >= 0 ) {
	  for ( j = 0; j < 4; j++ ) if ( qm -> mesh -> qAdj [ 4 * adj + j ] == qID[0] )
	    qm -> mesh -> qAdj [ 4 * adj + j ] = qID[1];
      }
  }
  for ( i = 0 ; i < 4; i++ ) {
      j = qm -> mesh -> qIdx [ 4 * (qID[0] - 1) + i ] - 1;
      qm -> mesh -> valence[ j ][1] = qID[0];
      j = qm -> mesh -> qIdx [ 4 * (qID[1] - 1) + i ] - 1;
      qm -> mesh -> valence[ j ][1] = qID[1];
  }
  idx[0] = 0 ; idx[1] = 3; idx[2] = swap; idx[3] = swap + 3;
  for ( i = 0 ; i < 4; i++ ) {
      stat = setValence ( qm -> mesh , qg.verts[idx[i]] );
      if ( stat != EGADS_SUCCESS ) {
	  printf(" Inside swapping operation set valence for %d --> %d !!!\n ", qg.verts[idx[i]], stat );
	  EG_free(quad );
	  return stat;
      }
  }
  for ( j  = 0 ; j < 2; j++ ) {
      i = qg.verts[ j * 3 ];
      if ( qm -> mesh -> vType [ i - 1] != -1 ) continue;
      q    = 0 ; if ( EG_quadVertIdx ( qm ->mesh, qID[q], i ) < 0 ) q = 1;
      area = EG_quadAngleOrientation(qm, qm -> minAngle, qm -> maxAngle, qID[q], ori, piv, angles);
      if ( area < 0  ) return area;
      if ( area == 0 ) EG_averageCoords (qm, i, qm -> minAngle, qm -> maxAngle);
  }
  *activity = 1;
  if ( EG_makeValidMesh(qm, 6, qg.verts,0) == EGADS_SUCCESS ) {
      EG_free(quad );
      printf(" SWAP WENT WELL \n ");
      return EGADS_SUCCESS;
  }
#ifdef DEBUG
  printf(" SWAP DIDNT WORK MESH IS INVALID \n ");
  printMesh (qm , NULL );
#endif
  stat = EG_restoreQuads (qm, quad, nq);
  EG_free(quad );
  if ( stat != EGADS_SUCCESS ) return stat;
#ifdef DEBUG
  printMesh (qm , NULL );
#endif
  *activity = 0;
  return EGADS_SUCCESS;
}

static int EG_split (meshMap *qm, int qID, int *activity) {
  int poly[3], val[3], v,  q, id0 = 0, i, stat, dist = 0, validSplit = 0, id6[2], links[4];
  vStar *star = NULL;
  *activity   = 0;
  if ( qm -> mesh -> qIdx [ 4 * (qID - 1 ) ] == -2 ) {
#ifdef DEBUG
      printf(" EG_split for quad %d is %d \n ", qID, stat );
#endif
      return EGADS_SUCCESS;
  }
  for ( v = 0; v < 4; v++ ) {
      poly[0]     = qm -> mesh -> qIdx [ 4 * ( qID - 1 ) + v ];
      val[0]      = getValence ( qm -> mesh, poly[0] );
      if  ( val[0] < 5 || qm -> mesh -> vType [ poly[0] -1 ] == 4 ) continue;
      stat        = EG_buildStar(qm -> mesh, &star, poly[0] );
      if  ( stat != EGADS_SUCCESS || star == NULL ) return stat;
      id6 [0]     = -1;
      id6 [1]     = -1;
      for ( validSplit = q = 0 ; q < star -> nQ; ++q) {
	  if ( star -> quads[q] == - 1 ) continue;
	  id0     = 2 * q + 1;
	  poly[1] = star -> verts[id0];
	  val[1]  = getValence ( qm -> mesh, poly[1] ) ;
	  for ( i = 0 ; i <= qm -> mesh -> valence [ poly[0] - 1 ][0] - 4; ++i) {
	      dist    = 4 + 2*i;
	      poly[2] = star -> verts  [ star -> idxV [id0 + dist] ];
	      if ( poly[2] < 0 ) continue; // ghost vertex
	      val[2]  = getValence ( qm -> mesh, poly[2] ) ;
	      if ( val[1] == 3 && val[2] == 3 && qm -> mesh -> vType [ poly[0] -1 ] == -1 ) {
		  if ( ( star -> nQ == 6 && dist == 6 ) || star -> nQ != 6 ) {
		      validSplit = 1;
		      break;
		  }
	      }
	      else if ( val[1] <= 3 && val[2] <= 4 && id6[0] == -1 ) id6[0] = id0;
	      else if ( val[2] <= 3 && val[1] <= 4 && id6[0] == -1 ) id6[0] = star -> idxV[id0 + dist];
	      else if ( val[1] <= 4 && val[2] <= 4 && id6[0] == -1 ) id6[1] = id0;
	  }
	  if ( validSplit == 1) break;
      }
      if ( validSplit == 0 && val[0] >= 6 && ( id6[0] != -1 || id6[1] != -1 ) ) {
	  validSplit  = 1;
	  if ( qm -> mesh -> vType [ poly[0] - 1] == 3 ) { // boundary vertex: Special split since regular = 3 ( not val 4)
	      for ( q = 0 ; q < star -> nQ; q++ )
		if ( star -> quads[ q ] == -1 ) break;
	      id6[0]  = - 1; id6[1] = -1;
	      for ( i = 0 ; i < 2 ; i++ ) {
		  id0 = star -> idxQ [ q + i ] ;
		  if ( i == 1 ) dist = (star -> nV - 1) - 4;
		  else          dist = 4 + 2 * i;
		  links[2 * i     ]  = 2 * id0 + 1;
		  links[2 * i + 1 ]  = star -> idxV [ 2 * id0 + 1 + dist ];
		  if (         getValence (qm -> mesh, star -> verts [ links[2*i    ] ]) == 4 ) {
		      if      (getValence (qm -> mesh, star -> verts [ links[2*i + 1] ]) == 3 && id6[0] == -1 ) id6[0] = i;
		      else if (getValence (qm -> mesh, star -> verts [ links[2*i + 1] ]) == 4 && id6[1] == -1 ) id6[1] = i;
		  }
	      }
	      dist    = 4;
	      if      ( id6[0] != -1 ) id0 = links[ 3 * id6[0]];
	      else if ( id6[1] != -1 ) id0 = links[ 3 * id6[1]];
	      else              validSplit = 0;
	  }
	  else {
	      dist = 6;
	      id0  = id6[0]; if ( id0 < 0 ) id0 = id6[1];
	  }
      }
      if ( validSplit == 1)  {
	  stat      = EG_splittingOperation (qm, star -> verts[0],
					     star -> verts[id0], star -> verts[ star -> idxV[id0 + dist]],
					     &(*activity ));
	  EG_freeStar ( &star);
	  return stat;
      }
      EG_freeStar ( &star);
  }
  return EGADS_SUCCESS;
}


static int
EG_splittingOperation(meshMap *qm, int vC, int vL, int vR, int *activity ) {
  int qIdx[4], modQ[4], verts[4], adj[2], poly[4], q, newQ, i, j, stat, id0 = -1, id1 = -1, dist, links[4], vals[4], addedV = 0, nq;
  vStar  *star = NULL;
  Quad   *quad = NULL;
  *activity    = 0;
  stat         = EG_buildStar(qm -> mesh, &star, vC );
#ifdef DEBUG
  printf(" SPLITTING %d thru %d %d \n ", vC, vL, vR );
#endif
  if ( stat   != EGADS_SUCCESS || star == NULL ) {
      printf(" In splittingOperation build star %d is NULL %d !!\n", vC, stat );
      return stat;
  }
  nq         = star -> nQ;
  stat       = EG_backupQuads (qm, &nq, star -> quads, & quad );
  if ( stat != EGADS_SUCCESS ) {
      EG_free (quad );
      printf("In EG_splittingOperation :: EG_backupQuads %d !!\n", stat );
      return stat;
  }
  if ( qm -> mesh -> remQuads[0] > 0 ) {
      poly[3] = qm -> mesh -> remVerts[qm -> mesh -> remVerts[0]--];
      newQ    = qm -> mesh -> remQuads[qm -> mesh -> remQuads[0]--];
  }
  else {
      poly[3] = qm -> mesh -> totVerts + 1;
      newQ    = qm -> mesh -> totQuads + 1;
      if ( qm -> mesh -> totVerts > 2 * qm -> oriV ) {
	  printf(" We have duplicated the number of initial vertices. This is too much. \n");
	  EG_free (quad );
	  return EGADS_INDEXERR;
      }
      ++qm -> mesh -> totVerts; ++qm -> mesh -> totQuads;
      addedV = 1;
  }
  qm -> mesh -> vType [poly[3] - 1] = -1;
  if ( qm -> mesh -> vType [ vC - 1] != -1 ) {
      for ( q = 0 ; q < star -> nQ; q++ )
	if ( star -> quads[ q ] == -1 ) break;
      for ( i = 0 ; i < 2 ; i++ ) {
	  if ( i == 1 ) dist = (star -> nV - 1) - 4;
	  else          dist = 4;
	  id0 = star -> idxQ [ q + i ] ;
	  links[2 * i     ]  = star -> verts [ star -> idxV[2 * id0 + 1       ]];
	  links[2 * i + 1 ]  = star -> verts [ star -> idxV[2 * id0 + 1 + dist]];
	  vals [2 * i     ]  = getValence (qm -> mesh, links[2*i    ]);
	  vals [2 * i + 1 ]  = getValence (qm -> mesh, links[2*i + 1]);
      }
      dist = 4;
      i = 0;
      if ( vals[0] * vals[1] > vals[2] * vals[3] ) {
	  vL = links[3];
	  vR = links[2];
      }
      else {
	  vL = links[0];
	  vR = links[1];
      }
  }
  id0 = - 1; id1 = -1;
  for ( j = 0 ; j < star -> nQ; j++ ) {
      if ( star -> verts [ 2 * j + 1] == vL ) id0 = j;
      if ( star -> verts [ 2 * j + 1] == vR ) id1 = j;
  }
  if ( id0 == -1 || id1 == -1 ) {
      EG_freeStar (&star ) ;
#ifdef DEBUG
      printf(" Problem finding vertices in EG_splittingOperation %d !!\n ", EGADS_INDEXERR );
#endif
      EG_free ( quad );
      return EGADS_INDEXERR;
  }
  poly [0] = star -> verts[   0];
  poly [1] = star -> verts[ 2 * id0 + 1];
  poly [2] = star -> verts[ 2 * id1 + 1];
  qIdx [0] = id0;
  qIdx [1] = star -> idxQ[ id1 + star -> nQ -1];
  qIdx [2] = id1;
  qIdx [3] = star -> idxQ[ id0 + star -> nQ -1];
  verts[0] = poly[1];
  verts[1] = poly[2];
  verts[2] = poly[2];
  verts[3] = poly[1];
  qm -> mesh -> qIdx[ 4 * ( newQ - 1)    ] = poly[1];
  qm -> mesh -> qIdx[ 4 * ( newQ - 1) + 1] = poly[0];
  qm -> mesh -> qIdx[ 4 * ( newQ - 1) + 2] = poly[2];
  qm -> mesh -> qIdx[ 4 * ( newQ - 1) + 3] = poly[3];
  for ( i = 0 ; i < 4; ++i) {
      modQ[i]  = star -> quads [ qIdx[i] ] ;
      qm -> mesh -> qAdj[ 4 * ( newQ - 1) + i ] = modQ[i];
      qm -> mesh -> valence[ qm -> mesh -> qIdx [ 4 * ( newQ - 1 )  + i ] -1 ] [1]  = newQ;
      if ( modQ[i] == -1 ) continue;
      stat = EG_adjQtoPair ( qm -> mesh, modQ[i], poly[0], verts[i], adj );
      qm -> mesh -> qAdj [ 4 * ( modQ[i] - 1) + adj[0] ]  = newQ;
  }
  j = qIdx[2];
  q = star -> quads[j++];
  while ( q != star -> quads [ qIdx[0]] ) {
      for ( i = 0 ; i < 4; ++i)
	if (qm -> mesh -> qIdx[ 4 * (q - 1) + i] == poly[0] )
	  qm -> mesh -> qIdx[ 4 * (q - 1) + i] = poly[3];
      q = star -> quads [ star -> idxQ[j++] ];
      if ( q == -1 ) {
	  if ( star -> quads [ qIdx[0] ]== -1 ) break;
	  else q = star -> quads [ star -> idxQ[j++] ];
      }
  }
  // Add valences to splitting vertices
  for ( i = 0 ; i < 4; i++ ) {
      stat = setValence( qm -> mesh, poly[i] ) ;
      if ( stat != EGADS_SUCCESS ) {
	  EG_freeStar(&star);
	  EG_free ( quad );
	  return stat;
      }
  }
  // Update valences at links
  for ( i = 0 ; i < 4; i++ ) {
      for ( j = 0 ; j < qm -> mesh -> valence [ poly[i] - 1][0]; j++ ) {
	  stat = setValence( qm -> mesh, qm -> mesh -> valence [ poly[i] - 1 ][ 2 + j ] );
	  if ( stat != EGADS_SUCCESS ) {
	      EG_freeStar(&star);
	      EG_free ( quad );
	      return stat;
	  }
      }
  }
  EG_freeStar (&star );
  updateVertex (qm, poly[3], &qm -> mesh -> uvs[ 2 * ( poly[0] - 1 )]);
  EG_averageCoords ( qm, poly[0], qm -> minAngle, qm -> maxAngle);
  EG_averageCoords ( qm, poly[3], qm -> minAngle, qm -> maxAngle);
  *activity = 1;
  if ( EG_makeValidMesh(qm, 4, poly, 0) == EGADS_SUCCESS ) {
      EG_free (quad );
      return EGADS_SUCCESS;
  }
#ifdef DEBUG
  printf(" MESH SPLIT FAILED \n ");
  printMesh (qm , NULL );
#endif
  stat = EG_restoreQuads (qm, quad, nq);
  EG_free(quad );
  if ( stat != EGADS_SUCCESS ) return stat;
  if ( addedV == 1 ) {
      --qm -> mesh -> totVerts;
      --qm -> mesh -> totQuads;
  }
  else {
      qm -> mesh -> remVerts[0]++;
      qm -> mesh -> remQuads[0]++;
  }
  for ( i = 0 ; i < 4; i++ )
    qm -> mesh -> qIdx[ 4 * ( newQ - 1) + i ] = -2;
  qm -> mesh -> vType [ poly[3] -1 ] = -2;
#ifdef DEBUG
  printMesh (qm , NULL );
#endif
  *activity = 0;
  // New Vertex location: take mid points of quads to which they belong
  return stat;
}


static int
EG_doubleSwap( meshMap *qm, quadGroup qg, int forcing, int *activity ) {
  int  piv5 = -1, piv3 = -1, q5, i, adjPiv5, stat, adj[2], swap = 0, vopp3;
  *activity = 0;
  if      (validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0 ) return EGADS_SUCCESS;
  else if (qg. vals[0] == 4 && forcing  == 0 ) return EGADS_SUCCESS;
  else if (qg. vals[0] == 5 ) {
      if  (qg. vals[2] * qg. vals[4] == 15 ) {
	  piv3 = 4;
	  if ( qg. vals[2] == 3 ) piv3 = 2;
	  swap = 0;
      }
      else if ( qg.vals[1] == 3 || qg.vals[5] == 3 ) {
	  piv3 = 1; if ( qg.vals[1] != 3 ) piv3 = 5;
	  swap = 1;
      }
  }
  else if ( forcing == 1 ) {
      if  ( qg.vals[1] != 3 && qg.vals[5] != 3 ) return EGADS_SUCCESS;
      piv3 = 1; if ( qg.vals[1] != 3 ) piv3 = 5;
      if ( qg.vals[1] != 5 && qg.vals[5] != 5 ) swap = 1;
  }
  if ( piv3 == -1 ) return EGADS_SUCCESS;
  vopp3 = (piv3 + 3)%6;
  if ( getValence(qm -> mesh, qg.verts[vopp3] ) != 4 ) return EGADS_SUCCESS;
  if ( swap == 0 ) {
      piv5 = (vopp3 + 1 )%6;
      if ( piv5 %3 == 0 ) piv5 = (vopp3 + 5 )%6;
  }
  else {
      piv5 = (vopp3 + 1 )%6;
      if ( piv5 %3 != 0 ) piv5 = (vopp3 + 5 )%6;
  }
  q5 = 0;
  if ( vopp3 > 3 ) q5 = 1;
  stat       = EG_adjQtoPair ( qm -> mesh, qg.q[q5], qg.verts[vopp3], qg.verts[piv5], adj);
  if ( stat != EGADS_SUCCESS || adj[1] == -1 ) {
#ifdef DEBUG
      if ( stat != EGADS_SUCCESS ) printf(" EG_doubleSwap adjToPair %d --> %d !!\n", adj[1], stat );
#endif
      return stat;
  }
  i       = EG_quadVertIdx ( qm -> mesh, adj[1], qg.verts[vopp3]);
  adjPiv5 = qm ->mesh -> qIdx [ 4 * ( adj[1] -1 ) + ( i + 1)%4];
  if (adjPiv5 == qg.verts[piv5] )
    adjPiv5 = qm ->mesh -> qIdx [ 4 * ( adj[1] -1 ) + ( i + 3)%4];
  if (     swap == 0 &&  getValence (qm ->mesh, adjPiv5) > 4 ) return EGADS_SUCCESS;
  else if (swap == 1 && (getValence (qm ->mesh, adjPiv5) <  5 ||
      (qm -> mesh -> vType[adjPiv5 -1 ] != -1 &&
	  qm -> mesh -> vType[qg.verts[vopp3] -1 ] !=-1 ) )) return EGADS_SUCCESS;
  piv5       = qg.verts[0];
  stat       = EG_swappingOperation (qm, qg, piv3, &(*activity) );
  if ( stat != EGADS_SUCCESS || *activity == 0 ) {
#ifdef DEBUG
      printf(" EG_doubleSwap: at first swap activity %d, stat --> %d !!\n ", *activity, stat );
#endif
      return stat;
  }
  stat = EG_swap (qm, adj[1], &i );
  if ( i == 0 && stat == EGADS_SUCCESS ) {
      stat = EG_createQuadGroup (qm -> mesh, &qg, qg.q[0], qg.q[1]);
      for ( swap = 0 ; swap < 6; swap++ ) if ( qg.verts[swap] == piv5) break;
      stat       = EG_swappingOperation (qm, qg, swap, &i);
      if ( stat == EGADS_SUCCESS && i == 1 ) *activity = 0;
      else {
#ifdef DEBUG
	  printf("In EG_doubleSwap undoing first swap activity %d stat %d \n!!", i, stat );
#endif
	  return stat;
      }
  }
  return stat;
}

static int EG_doubleSplit(meshMap *qm, quadGroup qg, int forcing, int *activity ) {
  int i, j, stat, f = 0 ;
  int piv[2] = {1, 5} ;
  *activity  = 0;
  if ( qg.vals[1] != 3 ) {
      piv[0] = 5; piv[1] = 1;
  }
  if ( qg.vals[piv[0] ] != 3 ) return EGADS_SUCCESS;
  if ((qm -> mesh -> vType[ qg.verts[0]      - 1] >= 0 &&
      qm -> mesh -> vType [ qg.verts[piv[1]] - 1] >= 0 ) ||
      qm -> mesh -> vType [ qg.verts[0]      - 1] == 4) return EGADS_SUCCESS;
  if ((forcing == 0 && (qg.vals[0] != 5 || qg.vals[1] * qg.vals[5]  != 15)) ||
      (forcing == 1 && (qm ->extraQuads > 0 || qg.vals[0] * qg.vals[piv[1]] <= 16))) return EGADS_SUCCESS;
  if ( qg.vals[1] * qg.vals[5] != 15 || qg.vals[0] != 5 ) f = 1;
  stat       = EG_splittingOperation (qm, qg.verts[0], qg.verts[piv[0]], qg.verts[piv[1]], &(*activity));
  if ( stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
      printf("In EG_doubleSplit: force 1st split through %d - %d --> %d activity %d !!\n ",
	     qg.verts[0], qg.verts[piv[0]], stat, *activity );
#endif
      return stat;
  }
  for ( j    = 0 ; j < 2; j++)
    if ( EG_quadVertIdx (qm -> mesh, qg.q[j], qg.verts[piv[1]] ) >= 0 ) break;
  stat = EG_split (qm, qg.q[j], &i);
  if ( stat == EGADS_SUCCESS && i == 0 ) {
      j     = qm -> mesh -> valence[ qg.verts[0] -1][1];
      stat  = EG_collapse (qm, j, &i, 1 );
      if ( stat == EGADS_SUCCESS && i == 1 ) *activity = 0;
      return stat;
  }
  if ( f == 1 && *activity > 0 ) {
      qm -> extraQuads += 2;
  }
  return stat;
}

static int EG_swapSplit(meshMap *qm,quadGroup qg, int forcing, int *activity  ) {
  int  stat, i, j, i3 = -1, i5 = -1, v3opp = -1, q5, vL5, vL5adj, swap = 0, adj[2], f = 0, v0;
  vStar *star = NULL;
  *activity   = 0;
  if  ( qg.vals[0] * qg.vals[3] != 20 ||
      validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0 ) return EGADS_SUCCESS;
  if      (qg.vals[1] == 3 ) i3 = 1;
  else if (qg.vals[5] == 3 ) i3 = 5;
  for ( i = 1; i < 6; i++)
    if (qg.vals[i] == 5 )  i5 = i;
  if (i3   != -1 && ( i5 == -1 ||( i5 != -1 && i5 == (i3 + 3)%6 ) )) {
      v3opp = ( i3 +3)%6;
      i5    =  v3opp;;
  }
  else if ( i3 == -1 && forcing == 1 && ( i5 == -1 || i5 == 2 || i5 == 4 )) {
      if ( qm -> extraQuads > 0 ) return EGADS_SUCCESS;
      if ( i5 == -1 ) v3opp = -1;
      else            v3opp = i5;
      f = 1;
  }
  else return EGADS_SUCCESS;
  if ( v3opp == -1 ) {
      for ( i  = 0 ; i < 2; i++ ) {
	  j    = 2 + 2 * i;
	  if ( i == 0 ) vL5 = 1;
	  else          vL5 = 5;
	  if (  getValence (qm ->mesh, qg.verts[j]) == 3 ) continue;
	  stat              = EG_adjQtoPair (qm -> mesh, qg.q[i], qg.verts[j], qg.verts[vL5], adj );
	  if ( stat        != EGADS_SUCCESS || adj[1] == -1 ) continue;
	  q5                = EG_quadVertIdx (qm -> mesh, adj[1], qg.verts[j]);
	  vL5adj            = qm -> mesh -> qIdx [ 4 * ( adj[1] - 1 ) + ( q5 + 1 ) %4 ];
	  if ( vL5adj == qg.verts[vL5] )
	    vL5adj = qm -> mesh -> qIdx [ 4 * ( adj[1] - 1 ) + ( q5 + 3) %4 ];
	  if ( getValence (qm ->mesh, vL5adj) == 3 ) {
	      i5   = j;
	      swap = j;
	      break;
	  }
      }
  } else {
      vL5 = ( v3opp + 1 ) % 6;
      if ( vL5 %3 == 0 ) vL5 = (v3opp + 5) % 6;
      q5 = 0;
      if ( EG_quadVertIdx  (qm -> mesh, qg.q[q5], qg.verts[vL5] ) < 0 ) q5 = 1;
      stat = EG_adjQtoPair (qm -> mesh, qg.q[q5], qg.verts[v3opp], qg.verts[vL5], adj );
      if ( stat     != EGADS_SUCCESS || adj[1] == -1 ) {
	  if ( stat != EGADS_SUCCESS) printf("EG_swapSplit: EG_adjQtoPair from quad %d is %d\n!!", qg.q[q5],stat );
	  return stat;
      }
      i           = EG_quadVertIdx (qm -> mesh, adj[1], qg.verts[v3opp]);
      vL5adj      = qm -> mesh -> qIdx [ 4 * ( adj[1] - 1 ) + ( i + 1 ) %4 ];
      if ( vL5adj == qg.verts[vL5] )
	vL5adj  = qm -> mesh -> qIdx [ 4 * ( adj[1] - 1 ) + ( i + 3) %4 ];
      if ( i3 != -1 && ( qg.vals[v3opp] == 5 || getValence (qm ->mesh, vL5adj) == 3 ) ) swap = i3;
      else if ( forcing == 1 && qm -> extraQuads <= 0 && ( qg.vals[v3opp] == 5 || getValence (qm ->mesh, vL5adj) == 3 ) ) swap = v3opp;
  }
  if ( swap %3 == 0 || (qm -> mesh -> vType[ qg.verts[swap] - 1] >= 0 &&
      qm -> mesh -> vType[ qg.verts[(swap + 3 ) % 6] - 1] >= 0 ) ) return EGADS_SUCCESS;
  stat       = EG_swappingOperation( qm, qg, swap, &(*activity) );
  v0         = qg. verts[0];
  if ( stat != EGADS_SUCCESS || *activity == 0 ) {
#ifdef DEBUG
      printf(" In swapSplit thru %d : EG_swappingOperation went %d activity %d !!\n ", swap, stat, *activity );
      printQuadGroup(qm ->mesh, qg);
#endif
      return stat;
  }
  stat       = EG_buildStar (qm -> mesh, &star, qg.verts[i5] );
  if ( stat != EGADS_SUCCESS || star == NULL ) {
      printf("In swapSplit build star for %d --> %d!!\n", qg.verts[i5], stat);
      return stat;
  }
  for ( i = 0 ; i < star -> nQ; i++ )
    if ( star -> verts [2 * i + 1] == qg.verts[3] ) break;
  if ( star -> verts [2 * i + 1] != qg.verts[3] ) {
      EG_freeStar (&star );
      return EGADS_INDEXERR;
  }
  adj[0] = star -> verts[ star -> idxV[2 * i + 1 + 6] ];
  if ( star -> nQ == 5 ) {
      adj[1] = star -> verts[ star -> idxV[2 * i + 1 + 4] ];
      if ( getValence ( qm -> mesh, adj[1] ) < getValence ( qm -> mesh, adj[0] ) ) adj[0] = adj[1];
  }
  EG_freeStar (&star );
  if ( qm -> mesh -> vType [ qg.verts[i5] - 1] == 4 ) return EGADS_GEOMERR;
  stat = EG_splittingOperation (qm, qg.verts[i5], qg.verts[3], adj[0], &i);
  if ( stat != EGADS_SUCCESS || i == 0 ) {
#ifdef DEBUG
      printf("In swapSplit splittingOperation --> %d !!\n", stat);
#endif
      stat = EG_createQuadGroup (qm -> mesh, &qg, qg.q[0], qg.q[1]);
      for ( swap = 0 ; swap < 6; swap++ ) if ( qg.verts[swap] == v0 ) break;
      stat       = EG_swappingOperation (qm, qg,swap, &i );
      if ( stat == EGADS_SUCCESS && i == 1 ) *activity = 0;
  }
  if ( f == 1 && *activity > 0 ) qm -> extraQuads++;
  return stat;
}

static int EG_swapCollapse (meshMap *qm,quadGroup qg, int forcing, int *activity  ) {
  int  stat, i, i3 = -1, q5, qC, vL5, vL5adj, vOpp3, swap = 0, adj[2], f = 0, v0;
  *activity   = 0;
  if (validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0 ) return EGADS_SUCCESS;
  if (qg.vals[1] * qg.vals[2] == 15 ) {
      i3 = 1;
      if ( qg.vals[i3] == 5 )i3 = 2;
  }
  else if (qg.vals[4] * qg.vals[5] == 15 ) {
      i3 = 4;
      if ( qg.vals[i3] == 5 ) i3 = 5;
  }
  else if(qg.vals[0] * qg.vals[3] == 20 ) {
      if      (qg.vals[1] == 3 ) i3 = 1;
      else if (qg.vals[5] == 3 ) i3 = 5;
      else return EGADS_SUCCESS;
      if ( i3 == 1 ) {
	  q5   = 0;
	  stat = EG_adjQtoPair (qm -> mesh, qg.q[0], qg.verts[2], qg.verts[3], adj );
      }
      else {
	  q5   = 1;
	  stat = EG_adjQtoPair (qm -> mesh, qg.q[1], qg.verts[4], qg.verts[3], adj );
      }
      if ( stat != EGADS_SUCCESS || adj[1] == -1 ) {
	  if ( stat != EGADS_SUCCESS ) printf("EG_swapCollapse centre 4-5 adj to pair %d !!\n", stat );
	  return stat;
      }
      stat = EG_createQuadGroup (qm ->mesh, &qg, qg.q[q5], adj[1] );
      if (stat != EGADS_SUCCESS ) {
	  printf("Inside EG_swapCollapse :: EG_createGroup %d !!\n ", stat );
	  return stat;
      }
      return EG_swapCollapse (qm, qg, forcing, &(*activity) );
  }
  else if ( forcing == 1 ) {
      f = 1;
      if ( qm  -> extraQuads < 0 ) return EGADS_SUCCESS;
      if( qg.vals[1] * qg.vals[5] == 25 || qg.vals[2] * qg.vals[4] == 25 ) {
	  i3 = 2; if ( qg.vals[1] != 5 ) i3 = 1;
      }
      else if( qg.vals[1] * qg.vals[4] == 15 ) {
	  i3 = 1; if ( qg.vals[1] != 3 ) i3 = 4;
      }
      else if ( qg.vals[2] * qg.vals[5] == 15 ) {
	  i3 = 2; if ( qg.vals[2] != 3 ) i3 = 5;
      }
      else if(qg.vals[1] * qg.vals[5] == 9 ||
	  qg.vals[2] * qg.vals[4] == 9 ) {
	  v0         = qg.verts[0];
	  stat       = EG_swappingOperation (qm, qg, 1, &(*activity) );
	  if ( stat != EGADS_SUCCESS || *activity == 0 ) {
#ifdef DEBUG
	      printf("forcing swapcollapse:: after swapping %d activity %d !!\\n ", stat, *activity);
#endif
	      return stat;
	  }
	  qC   = qg.q[0];
	  if ( EG_nValenceCount ( qm -> mesh, qC, 3 ) < 2 ) qC = qg.q[1];
	  stat = EG_collapse (qm, qC, &(*activity), 1 );
	  if ( stat  != EGADS_SUCCESS || *activity == 0 ) {
	      stat    = EG_createQuadGroup (qm -> mesh, &qg, qg.q[0], qg.q[1]);
	      for ( swap = 0 ; swap < 6; swap++ ) if ( qg.verts[swap] == v0 ) break;
	      stat       = EG_swappingOperation (qm, qg, swap, &i );
	      if ( stat == EGADS_SUCCESS && i == 1 ) *activity = 0;
	  }
	  else  qm ->extraQuads--;
	  return stat;
      }
      else return EGADS_SUCCESS;
  } else return EGADS_SUCCESS;
  q5     = 0; if ( i3 > 3 ) q5 = 1;
  if (validCollapse (qm, qg.q[q5], qg.verts[i3]) == 0) return EGADS_SUCCESS;
  vOpp3  = (i3 + 3 ) % 6;
  q5     = 0; if ( vOpp3 > 3 ) q5 = 1;
  qC     = qg.q[ ( q5 +1)%2];
  vL5    = ( vOpp3 + 1) %6;
  if ( vL5 %3 != 0 ) vL5 = ( vOpp3 + 5) %6;
  if ( qg.vals[vOpp3] == 3 ) return EGADS_SUCCESS;
  stat       = EG_adjQtoPair (qm -> mesh, qg.q[q5], qg.verts[vOpp3], qg.verts[vL5], adj );
  if ( stat != EGADS_SUCCESS || adj[1] == -1 ) return stat;
  i      = EG_quadVertIdx ( qm -> mesh, adj[1], qg.verts[vOpp3]);
  vL5adj = qm -> mesh -> qIdx [ 4 * ( adj[1] - 1 ) + ( i + 1 ) %4 ];
  if (  vL5adj == qg.verts[vL5] ) vL5adj = qm -> mesh -> qIdx [ 4 * ( adj[1] - 1 ) + ( i + 3 ) %4 ];
  if ((     forcing == 0 && qg.vals[vL5] == 4 && getValence (qm -> mesh, vL5adj) >= 4 ) ||
      ( forcing == 1 && qm -> extraQuads >= 0 && getValence (qm -> mesh, vL5adj) > 4 ) ) return EGADS_SUCCESS;
  stat  = EG_createQuadGroup (qm -> mesh, &qg, qg.q[q5], adj[1]);
  if ( stat != EGADS_SUCCESS ) {
      printf("EG_swapCollapse before swap: EG_createQuadGroup --> %d !!\n", stat );
      printQuadGroup (qm -> mesh, qg );
      return stat;
  }
  if (validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0 ) return EGADS_SUCCESS;
  for ( swap = 0 ; swap < 6; swap++ )
    if ( qg.verts[swap] == vL5adj ) break;
  v0         = qg.verts[0];
  stat       = EG_swappingOperation (qm, qg, swap, &(*activity) );
  if ( stat != EGADS_SUCCESS || *activity == 0 ) {
#ifdef DEBUG
      printf("EG_swapCollapse after swapping %d activity %d !!\n", swap, *activity );
#endif
      return stat;
  }
  vL5        = qg.verts[0];
  stat       = EG_collapse (qm, qC, &(*activity), 1 );
  if ( stat != EGADS_SUCCESS || *activity == 0 ) {
#ifdef DEBUG
      printf("EG_swapCollapse after collapsing %d !!\n", stat );
#endif
      stat       = EG_createQuadGroup (qm -> mesh, &qg, qg.q[0], qg.q[1]);
      for ( swap = 0 ; swap < 6; swap++ ) if ( qg.verts[swap] == v0 ) break;
      stat       = EG_swappingOperation (qm, qg, swap, &i );
      if ( stat == EGADS_SUCCESS && i == 1 ) *activity = 0;
  }
  if ( f == 1 && stat == EGADS_SUCCESS ) qm -> extraQuads--;
  return stat;
}

static int EG_cleanNeighborhood (meshMap *qm, int qID,  int transfer, int *activity ) {
  int i, act, stat = 0, count = 0, j, v[4];
  vStar *star = NULL;
  *activity   = 0;
  if ( qm -> mesh -> qIdx [ 4 * (qID - 1 ) ] == -2 ) return EGADS_SUCCESS;
  for ( i = 0 ; i < 4; i++ ) v[i]  = qm -> mesh -> qIdx[ 4 * ( qID -1 ) + i ];
  stat        = EG_cleanQuad ( qm, qID, 1, transfer, 0, &act );
  if ( stat  != EGADS_SUCCESS ) return stat;
  *activity  += act;
  for ( i = 0 ; i < 4; i++ ) {
      if ( qm -> mesh -> vType[ v[i] -1 ] == -2 ) continue;
      stat       = EG_buildStar (qm -> mesh, &star, v[i]);
      if ( star == NULL || stat != EGADS_SUCCESS ) {
	  printf("EG_cleanNeighborhood : EG_buildStar  = %d \n", stat );
	  return stat;
      }
      for ( j = 0 ; j < star -> nQ; j++ ) {
	  if ( star -> quads[j] == -1  ) continue;
	  if ( EG_cleanQuad (qm, star -> quads[j], 1, transfer, 0, &act ) != EGADS_SUCCESS ) {
	      EG_freeStar(&star);
	      return stat;
	  }
	  count += act;
      }
  }
  EG_freeStar(&star);
  *activity += count;
#ifdef DEBUG
  return checkMesh (qm );
#endif
  return EGADS_SUCCESS;
}

static int
EG_wvsData(meshData *mesh, char *buffer) {
  int i;
  FILE *fil = NULL;
  char aux[100];
  if ( buffer == NULL ) {
      snprintf ( aux, 100, "wvs_%d", MESHPLOTS++ ) ;
      fil = fopen(aux, "w" );
      printf(" Writing in File %s  \n", aux);
  } else {
      fil = fopen (buffer, "w" );
      printf(" Writing in File %s  \n", buffer);
  }
  if ( fil == NULL ) return EGADS_MALLOC;
  fprintf(fil,"%d %d\n", mesh -> totVerts, mesh -> totQuads);
  for ( i = 0 ; i < mesh -> totVerts; i++ ) {
      fprintf(fil, "%lf %lf %lf \n",
	      mesh -> xyzs[ 3 * i],
	      mesh -> xyzs[ 3 * i + 1],
	      mesh -> xyzs[ 3 * i + 2] );
  }
  fprintf(fil,"\n");
  for ( i = 0 ; i < mesh -> totQuads; i++ ) {

      fprintf(fil, "%d %d %d %d\n",
	      mesh -> qIdx[4 * i    ],
	      mesh -> qIdx[4 * i + 1],
	      mesh -> qIdx[4 * i + 2],
	      mesh -> qIdx[4 * i + 3]);
  }
  fclose (fil);
  printf(" \n Written in file %s\n ", buffer);
  return EGADS_SUCCESS;
}


static int EG_fullMeshRegularization(meshMap *qm )
{
  int   ITMAX, it = 0, stat, activity = 0, totActivity = 0, loopact, i, j, k, kk, q, transfer = 0;
  int   best_iV, iV, quadPair[2], prevPair[2], totV, vQ;
  char buffer[100];
  // GET RANGE FOR EACH POINT
  meshCount (qm -> mesh, &iV, &totV, &vQ);
  fprintf(stderr,"Original mesh has %d QUADS, %d / %d irregular vertices (%f percent) ============\n", vQ, iV, totV, (double) iV * 100.0 / (double)totV );
  ITMAX = 20;
  it    = 0;
  do {
      totActivity   = 0;
      for (i = 0 ; i < qm -> mesh -> totQuads; i++ ) {
	  if ( qm ->mesh -> qIdx[ 4 * i ] == -2 ) continue; //can be a deleted quad
	  stat         = EG_cleanNeighborhood ( qm, i + 1, 0, &activity );
	  totActivity += activity;
	  if ( stat  != EGADS_SUCCESS ) {
	      printf("In EG_cleanMesh: EG_CleanNeighborhood for quad %d --> %d!!\n ",i + 1, stat );
	      return stat;
	  }
      }
      it++;
      meshCount(qm -> mesh, &iV, &totV, &vQ);
  } while (totActivity > 0 && it < ITMAX && iV > 2);
  stat       = EG_restoreMeshMap (qm, qm -> bestMesh, qm -> mesh);
  if ( stat != EGADS_SUCCESS ) {
      printf(" EG_restoreMeshMap for qb -> bestMesh after basic %d !!\n ", stat );
      return stat;
  }
  best_iV    = iV;
  if ( iV > 2 ) {
      for ( kk  = 0 ; kk <= 1; kk++ ) {
	  k     = 1;
	  qm -> extraQuads = 0;
	  it  = 0 ;
	  if ( kk == 0 ) fprintf(stderr,"Translating valences around mesh round 1\n ");
	  else           fprintf(stderr,"Translating valences around mesh round 2 :: forcing \n ");
	  do {
	      it++;
	      totActivity  =  0;
	      prevPair[0]  = -1;
	      prevPair[1]  = -1;
	      for (q       = 0 ; q < qm -> mesh-> totQuads; q++) {
		  transfer = 0 ;
		  if ( qm ->mesh -> qIdx[ 4 * q ] == -2 || q + 1 == prevPair[0] || q + 1 == prevPair[1] ) continue;
		  quadPair[0] = q + 1;
		  quadPair[1] = -1;
		  if ( EG_quadIsBoundary ( qm ->mesh, q + 1 ) != 0 && k == 0 ) continue;
		  stat          = EG_transferValences ( qm, quadPair, kk, &transfer, &activity );
		  if (stat     != EGADS_SUCCESS ) return stat;
		  if (activity == 0 || quadPair[0] < 0 || qm ->mesh -> qIdx[ 4 * ( quadPair[0] -1 )] == -2 ) continue;
		  totActivity += activity;
		  meshCount(qm -> mesh, &iV, &totV, &vQ );
		  if ( iV < best_iV ) {
		      stat = EG_restoreMeshMap ( qm, qm -> bestMesh, qm -> mesh);
		      if ( stat != EGADS_SUCCESS )  return stat;
		      best_iV = iV;
		  }
		  for ( j = 0 ; j < 20 ; j++ ) {
		      activity = 0;
		      for ( i  = 0 ; i < 2; i++ ) {
			  if (quadPair[i] < 0 || qm ->mesh -> qIdx[ 4 * ( quadPair[i] -1 )] == -2 ) continue;
			  stat       = EG_cleanNeighborhood (qm, quadPair[i],  transfer, &loopact);
			  if ( stat != EGADS_SUCCESS ) {
			      return stat;
			  }
			  activity  +=loopact;
		      }
		      if ( activity > 0 ) break;
		      stat       = EG_transferValences ( qm, quadPair, kk, &transfer, &activity );
		      if ( stat != EGADS_SUCCESS ) {
			  return stat;
		      }
		      if (activity == 0 || quadPair[0] < 0 || qm ->mesh -> qIdx[ 4 * ( quadPair[0] -1 )] == -2 ) break;
		  }
		  meshCount(qm -> mesh, &iV, &totV, &vQ );
		  if ( iV < best_iV ) {
		      stat = EG_restoreMeshMap ( qm, qm -> bestMesh, qm -> mesh);
		      if ( stat != EGADS_SUCCESS ) {
			  return stat;
		      }
		      best_iV = iV;
		  }
		  if ( iV <=2 ) break;
		  prevPair[0] = quadPair[0];
		  prevPair[1] = quadPair[1];
	      }
	      if ( iV <=2 ) break;
	  } while ( totActivity > 0 && it < ITMAX && iV > 2);
#ifdef DEBUG
	  stat = checkMesh (qm );
#endif
	  if ( stat != EGADS_SUCCESS ) {
	      return stat;
	  }
	  meshCount(qm -> mesh, &iV, &totV, &vQ );
	  if ( iV  < best_iV ) {
	      EG_restoreMeshMap ( qm, qm -> bestMesh, qm -> mesh);
	      best_iV = iV;
	  }
	  if ( iV <=2 ) break;
      }
      if ( iV  > best_iV ) {
	  EG_restoreMeshMap ( qm, qm -> mesh, qm -> bestMesh);
	  best_iV = iV;
      }
  }
  meshCount(qm -> mesh, &iV, &totV, &vQ );
  printf("\n\n*******************************************************************\n");
  printf(" AFTER FULL REGULARIZATION: mesh has %d QUADS, %d / %d irregular vertices (%f percent) ============\n", vQ, iV, totV, (double) iV * 100.0 / (double)totV );
  printf("*******************************************************************\n");
  fprintf(stderr,"\n\n*******************************************************************\n");
  fprintf(stderr," AFTER FULL REGULARIZATION: mesh has %d QUADS, %d / %d irregular vertices (%f percent) ============\n", vQ, iV, totV, (double) iV * 100.0 / (double)totV );
  fprintf(stderr,"*******************************************************************\n");
  stat = resizeQm (qm) ;
  if ( stat != EGADS_SUCCESS ) {
      printf(" After resizing mesh:: %d \n", stat );
      return stat;
  }
  printMeshStats(qm, 2);
  sprintf(buffer,"wvsRegular_%d.txt", qm -> fID);
  stat = EG_wvsData(qm -> mesh, buffer);
  if ( stat != EGADS_SUCCESS ) {
      printf(" writing in wvs file %d !! \n ", stat );
      return stat;
  }
  sprintf(buffer,"gnuRegular_%d.txt", qm -> fID);
  printMesh(qm, buffer);
  stat  += EG_makeValidMesh(qm, 0, NULL, 1);
  sprintf(buffer,"gnuFinalMesh_%d.txt", qm ->fID);
  printMesh (qm, buffer );
  sprintf(buffer,"wvsFinalMesh_%d.txt", qm ->fID);
  stat += EG_wvsData(qm -> mesh, buffer);
  return stat;
}





//#ifdef STANDALONE
int main(int argc, char *argv[])
{
  clock_t      start_t, end_t, total_t;
  int          stat = 0,  f , i, j, iBody, oclass, mtype, nbody;
  int          atype, alen, *senses, min;
  const int    *ints;
  float        arg;
  double       box[6], size, params[3], time, frac ;
  const double *reals ;
  const char   *OCCrev, *string;
  ego          context, tess, model, geom, *bodies, *dum;
  bodyData     *bodydata;
  start_t = clock();
  if (argc < 5 ) {
      printf("\n Usage: %s = (1) filename, Tesselation:  (2) angle (3) maxlen and chord (4) \n\n", argv[0]);
      printf(" You have entered : \n");
      for ( i = 0 ; i < argc; i++ ) printf(" argv[%d] = %s\n ", i, argv[i] );
      return 1;
  }
  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %2d.%02d with %s\n\n", i, j, OCCrev);

  /* initialize */
  stat = EG_open(&context);
  if ( stat != EGADS_SUCCESS) return 1;
  stat = EG_loadModel(context, 0, argv[1], &model);
  printf(" EG_loadModel      = %d\n", stat);
  printf(" EG_getBoundingBox = %d\n", EG_getBoundingBox(model, box));
  printf("       BoundingBox = %lf %lf %lf\n", box[0], box[1], box[2]);
  printf("                     %lf %lf %lf\n", box[3], box[4], box[5]);
  printf(" \n");
  size = box[3]-box[0];
  if (size < box[4]-box[1]) size = box[4]-box[1];
  if (size < box[5]-box[2]) size = box[5]-box[2];

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
  // default granularity
  sscanf(argv[2], "%f", &arg);
  params[2] = arg;
  sscanf(argv[3], "%f", &arg);
  params[0] = arg;
  sscanf(argv[4], "%f", &arg);
  params[1] = arg;
  printf(" Using angle = %lf,  relSide = %lf,  relSag = %lf\n",
	 params[2], params[0], params[1]);
  params[0] *= size;
  params[1] *= size;
  if ( argc == 6 ) sscanf (argv[5], "%d", &FACECHOICE);
  /* fill our structure a body at at time */
  for (iBody = 0; iBody < nbody; iBody++) {
      stat = EG_attributeAdd(bodies[iBody], ".qParams",
			     ATTRSTRING, 4, NULL, NULL, "off");
      if (stat != EGADS_SUCCESS)
	printf(" Body %d: attributeAdd = %d\n", iBody, stat);
      EG_getTopology(bodies[iBody], &geom, &oclass,
		     &mtype, NULL, &j, &dum, &senses);
      bodydata[iBody].body  = bodies[iBody];
      bodydata[iBody].mtype = mtype;
      if (mtype == WIREBODY) {
	  printf(" Body %d: Type = WireBody\n", iBody+1);
      } else if (mtype == FACEBODY) {
	  printf(" Body %d: Type = FaceBody\n", iBody+1);
      } else if (mtype == SHEETBODY) {
	  printf(" Body %d: Type = SheetBody\n", iBody+1);
      } else {
	  printf(" Body %d: Type = SolidBody\n", iBody+1);
      }
      stat = EG_getBodyTopos(bodies[iBody], NULL, FACE,
			     &bodydata[iBody].nfaces, &bodydata[iBody].faces);
      i    = EG_getBodyTopos(bodies[iBody], NULL, EDGE,
			     &bodydata[iBody].nedges, &bodydata[iBody].edges);
      if ((stat != EGADS_SUCCESS) || (i != EGADS_SUCCESS)) {
	  printf(" EG_getBodyTopos Face = %d\n", stat);
	  printf(" EG_getBodyTopos Edge = %d\n", i);
	  return 1;
      }
      stat = EG_makeTessBody(bodies[iBody], params, &bodydata[iBody].tess);
      if (stat != EGADS_SUCCESS) {
	  printf(" EG_makeTessBody %d = %d\n", iBody, stat);
	  continue;
      }
      stat = EG_makeTessBody(bodies[iBody], params, &bodydata[iBody].tess);
      tess = bodydata[iBody].tess;
      stat = EG_quadTess(tess, &bodydata[iBody].tess);
      if (stat != EGADS_SUCCESS) {
	  printf(" EG_quadTess %d = %d  -- reverting...\n", iBody, stat);
	  bodydata[iBody].tess = tess;
	  continue;
      }
      EG_deleteObject(tess);
  }
  for (  iBody = 0; iBody < nbody; iBody++) {
      stat = EG_attributeRet(bodydata[iBody].tess, ".tessType", &atype,
			     &alen, &ints, &reals, &string);
      if (stat == EGADS_SUCCESS)
	/* get faces */
	stat = EG_createMeshMap(&bodydata[iBody]);
      if ( stat != EGADS_SUCCESS) goto cleanup;
      for (iBody = 0; iBody < nbody; iBody++) {
	  for (f = 0; f < bodydata[iBody].nfaces; ++f) {
	      if ( FACECHOICE >= 0 )  f = FACECHOICE;
	      printMeshStats ( bodydata[iBody].qm[f] , -1 );
	      printf("===================================================================\n");
	      printf("===================================================================\n");
	      printf(" FULL MESH REGULARIZATION FACE %d \n ", f + 1);
	      printf("===================================================================\n");
	      printf("===================================================================\n");
	      stat = EG_fullMeshRegularization(bodydata[iBody].qm[f] );
	      fprintf(stderr, " EG_fullMeshRegularization face %d / %d = %d \n ", f + 1, bodydata[iBody].nfaces,  stat );
	      if ( stat != EGADS_SUCCESS ) return stat;
	      printMeshStats (  bodydata[iBody].qm[f] , -1 );
	      if ( FACECHOICE >= 0  ) break;
	  }
      }
  }
  cleanup:
  printf(" STAT IN CLEANUP  %d  ", stat );
  for (iBody = 0; iBody < nbody; iBody++) {
      EG_destroymeshMap( &bodydata[iBody].qm, bodydata[iBody].nfaces);
      EG_free ( bodydata[iBody].faces );
      EG_free ( bodydata[iBody].edges );
  }
  EG_free(bodydata);
  end_t   = clock();
  total_t =(end_t - start_t);
  time    = (double) total_t / CLOCKS_PER_SEC;
  min     = floor   (time) / 60 ;
  frac    = time - min * 60;
  fprintf(stderr, "Total time taken by CPU: %d minutes and %f seconds\n", min, frac);
  fprintf(stderr, "Total M FILES %d\n", MESHPLOTS - 1);

  EG_close ( context );
  return 0;
}
//#endif



/********************   IO FUNCTIONS **********************/

static void
printQuadSpecs(meshData *mesh, int id) {
  int i, v = 0 , val = 0 ;
  --id;
  printf(" QUAD %d HAS VERTICES ",id +1);
  for ( i = 0 ; i < 4; ++i ) {
      v   = mesh -> qIdx[4 * id + i];
      if ( v <= 0 ) printf ( " DEL %d %d \n ", v, 0);
      else {
	  val = getValence ( mesh, v );
	  printf(" %d ( val = %d )  ",v, val);
      }
  }
  printf("\t AND ADJACENT QUADS ");
  for ( i = 0 ; i < 4; ++i ) printf(" %d ",mesh -> qAdj[4*id + i]);
  printf("\n");
}
/*
static void
printVertSpecs (meshData *mesh, int id) {
	if ( checkVertex (mesh, id) ) {
		printf(" !!!!!!    checkVertex = %d !! \n ", checkVertex (mesh, id) );
		return;
	}
	printf("-------------------------------------------------------\n");
	printf(" Vertex %d is type %d and has valence %d (actual %d quad %d )\n ",
			id, mesh -> vType[ id - 1], getValence (mesh, id ), mesh -> valence[id -1][0], mesh -> valence[id -1][1]);
	printf("-------------------------------------------------------\n");
	return;
}
 */


static void
printMeshStats(meshMap *qm, int sweep) {
  int i,len, val ;
  int intVal[100], boundVal[100];
  char buffer[100];
  FILE *fout;
  sprintf(buffer, "MESH_STATS_%d.txt", sweep);
  printf(" WRITING ON FILE %s\n",buffer);
  fout = fopen(buffer, "w");
  if ( fout == NULL )
    return;
  len = qm -> mesh ->totVerts;
  for ( i = 0; i < 100; ++i) {
      intVal  [i] = 0;
      boundVal[i] = 0;
  }
  for ( i = 0; i < len; ++i) {
      if ( qm -> mesh -> vType[i] !=  -2){
	  if ( qm -> mesh -> vType[i] == -1){
	      val = qm -> mesh -> valence[i][0];
	      ++intVal[val];
	  } else {
	      ++boundVal[qm -> mesh ->valence[i][0]];
	  }
      }
  }
  fprintf(fout,"--------------------- MESH STATS AFTER %d SWEEP --------------\n",sweep);
  fprintf(fout,"---- TOTAL VERTICES %d TOTAL QUADS %d --------------\n",
	  qm -> mesh -> totVerts - qm -> mesh -> remVerts[0],qm -> mesh -> totQuads - qm -> mesh -> remQuads[0] );
  fprintf(fout," INTERIOR VERTICES\n");
  for ( i = 0 ; i < 100; ++i) {
      if ( intVal[i]  > 0 ) fprintf(fout," VALENCE %d = %d VERTICES\n", i, intVal[i]);
  }
  fprintf(fout," BOUNDARY VERTICES\n");
  for ( i = 0 ; i < 100; ++i) {
      if ( boundVal[i]  > 0 ) fprintf(fout," VALENCE %d = %d VERTICES\n", i, boundVal[i]);
  }
  fclose(fout);
  return ;
}


static void
printMesh(meshMap *qm, char *name  ) {
  int i, j, k, v, d;
  double eval[18], average[2], dist;
  char name2[500];
  FILE *fout = NULL;
  if ( name == NULL ) {
      snprintf ( name2, 500, "M_%d",MESHPLOTS++ ) ;
      printf(" Writing in File %s  \n", name2);
      fout = fopen(name2, "w" );
  }
  else {
      printf(" Writing in File %s  \n", name);
      fout = fopen(name, "w" );
  }
  if ( fout == NULL ) return;
  for ( i = 0 ; i < qm -> mesh -> totQuads; ++i) {
      if ( qm -> mesh -> qIdx [ 4 * i] == -2 || qm -> mesh -> vType [ qm -> mesh -> qIdx [ 4 * i - 1] ] == -2 ) continue;
      for ( k = 0; k < 4; ++k) {
	  v  =   qm -> mesh -> qIdx[4*i + k] - 1;
	  fprintf(fout, "%lf %lf %lf %d %lf %lf \n",qm -> mesh->xyzs[3*v], qm -> mesh->xyzs[3*v +1], qm -> mesh->xyzs[3*v + 2], v + 1,
		  qm -> mesh -> uvs[2*v] , qm -> mesh -> uvs[2*v + 1]);
	  dist = 0.0;
	  EG_evaluate(qm -> face, &qm -> mesh -> uvs[2*v ], eval);
	  for ( d = 0 ; d < 3; ++d) dist += ( eval[d] - qm -> mesh->xyzs[3*v + d]) * ( eval[d] - qm -> mesh->xyzs[3*v + d]);
	  dist = sqrt  (dist);
	  if ( dist > EPS11 ) {
	      printf(" DIST = %11.2e  IN QUAD %d  VERTEX %d. UVs and xyzs are mismatched.  UV %lf  %lf \n",
		     dist,i+1, v+1, qm -> mesh -> uvs[2*v], qm -> mesh -> uvs[2*v + 1]);
	      for ( d = 0 ; d < 3; ++d) printf( "%lf  != %lf \t", eval[d], qm -> mesh->xyzs[3*v + d]);
	      fclose(fout);
	      exit (1);
	  }
      }
      v  =   qm -> mesh -> qIdx[4*i ] - 1;
      fprintf(fout, "%lf %lf %lf %d %lf %lf\n",qm -> mesh->xyzs[3*v], qm -> mesh->xyzs[3*v +1], qm -> mesh->xyzs[3*v + 2],
	      v + 1, qm -> mesh -> uvs[2*v] , qm -> mesh -> uvs[2*v + 1]);
      fprintf(fout,"\n\n");
      average[0] = 0.0; average[1] = 0.0;
      for (j = 0; j < 4; j++ ) {
	  k = qm -> mesh -> qIdx [4 * i + j ] - 1;
	  average[0] += 0.25 * qm -> mesh -> uvs [2 * k    ];
	  average[1] += 0.25 * qm -> mesh -> uvs [2 * k + 1];
      }
      EG_evaluate(qm -> face, average, eval);
      fprintf(fout, "%lf %lf %lf %d %lf %lf\n", eval[0], eval[1], eval[2], i + 1, average[0], average[1]) ;
      fprintf(fout,"\n\n");
  }
  fclose(fout);
}

static void printQuadGroup ( meshData *mesh, quadGroup qg ) {
  int i, stat ;
  for ( i = 0 ; i < 2; i++ ){
      stat   = checkQuad (mesh, qg.q[i]);
      if ( stat != EGADS_SUCCESS ) {
	  printf(" Quad %d is %d !!! \n", qg.q[i], stat );
	  return;
      }
      printQuadSpecs ( mesh, qg.q[i] );
  }
  for ( i = 0 ; i < 6; i++) {
      stat = checkVertex ( mesh, qg.verts[i]);
      if ( stat != EGADS_SUCCESS ) {
	  printf(" Vertex %d is %d !!! \n", qg.verts[i], stat );
	  return;
      }
      printf(" QV(%d) = %d type %d valence %d \n",
	     i, qg.verts[i], mesh -> vType [qg.verts[i] - 1], qg.vals[i]);
  }
}



