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

#define PI            3.1415926535897931159979635
#define MINVALIDANGLE 0.0
#define MAXVALIDANGLE 3.1415926535897931159979635
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
	int  *verts, *quads;
	int   nV, nQ; // nV = origin(1) + peaks (n)
	int  *idxV, *idxQ;
} vStar;

typedef struct {
	int     *quadIdx, *quadAdj, **valence, *vType, *remQuads, *remVerts;
	int     totQuads, totVerts, sizeVerts, sizeQuads ;
	double  *xyzs, *uvs ;
} meshData;


typedef struct {
	int      fID, actualQuads, oriQ, oriV, minQ, extraQuads,  boundaryVerts;
	ego      face;
	double   MINANGLE, MAXANGLE, minsize, maxsize;
	meshData *mesh, *backupMesh, *bestMesh;
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
static bodyData  *bodydata;
static int        nbody;
static float      focus[4];
static char       buffer[500];

static void swapInt ( int *a, int *b ) {
	int c;
	c  = *a;
	*a = *b;
	*b =  c;
}
static void unitVector(double *v, double *norm) {
	double n;
	n     = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	n     = sqrt(n);
	if(n > EPS11 ) {
		v[0] /=n; v[1] /=n; v[2] /=n;
	}
	else {
		v[0] = 0.0; v[1] = 0.0; v[2] = 0.0;
	}
	*norm = n;
	return;
}

static void
cross_product(double *A, double *B, double *cross) {
	cross[0] = A[1] * B[2] - A[2] * B[1];
	cross[1] = A[2] * B[0] - A[0] * B[2];
	cross[2] = A[0] * B[1] - A[1] * B[0];
}


static double
dotProduct (double *v1, double *v2 ) {
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] ;
}


static void EG_freeStar ( /*@null@*/ vStar **star ) {
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
	EG_free(( *mesh) -> quadIdx );
	EG_free(( *mesh) -> quadAdj );
	EG_free(( *mesh) -> uvs     );
	EG_free(( *mesh) -> xyzs    );
	EG_free(( *mesh) -> remQuads);
	EG_free(( *mesh) -> remVerts);
	EG_free(( *mesh) -> vType   );
	EG_free(( *mesh)            );
	*mesh = NULL;
}
static int  EG_allocMeshData ( meshData **mesh, int nQ, int nV ) {
	int j;
	if ( *mesh ) EG_freeMeshData( mesh );
	(*mesh)               = (meshData*) EG_alloc (sizeof (meshData) );
	if ( (*mesh) == NULL ) return EGADS_MALLOC;
	(*mesh) -> sizeVerts  = nV;
	(*mesh) -> sizeQuads  = nQ;
	(*mesh) -> xyzs       = (double *) EG_alloc(3*nV *sizeof(double ));
	(*mesh) -> uvs        = (double *) EG_alloc(2*nV *sizeof(double ));
	(*mesh) -> vType      = (int    *) EG_alloc(  nV *sizeof(   int ));
	(*mesh) -> quadIdx    = (int    *) EG_alloc(4*nQ *sizeof(   int ));
	(*mesh) -> quadAdj    = (int    *) EG_alloc(4*nQ *sizeof(   int ));
	(*mesh) -> remQuads   = (int    *) EG_alloc(  nQ *sizeof(   int ));
	(*mesh) -> remVerts   = (int    *) EG_alloc(  nV *sizeof(   int ));
	(*mesh) -> valence    = (int   **) EG_alloc(  nV *sizeof(   int*));
	if (    (*mesh)->quadIdx  == NULL || (*mesh)->quadAdj  == NULL ||
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


#define AREA2D(a,b,c)   ((a[0]-c[0])*(b[1]-c[1]) -  (a[1]-c[1])*(b[0]-c[0]))
#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
		a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
		a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define MAX(a,b)	(((a) > (b)) ? (a) : (b))
#define MIN(a,b)        (((a) < (b)) ? (a) : (b))


typedef struct {
	double uv[2];         /* (u,v) for node */
	double duv[2];        /* delta for coordinate update */
	double area;          /* accumulated area; -1 is boundary node */
	double xyz[3];	/* xyz for the node */
} Node;


/* perform the laplacian smoothing on the grid vertices */













/* IO FUNCTIONS */
static void printStar               (vStar *star) ;
//static void printStarFile(vStar *star, meshMap *qm) ;
static void printMeshStats(meshMap *qm, int sweep);
static void printMesh     (meshMap *qm,char[], int  );
static int  checkMesh     (meshMap *qm) ;

static int EG_normalToSurface       (meshMap *qm, double *uv, double *normal )  ;
static int EG_projectToTangentPlane (double normal[], double *O, double *p, double *proj) ;
/***********************************************************/
/* MESH MAP FUNCTIONS*/
static int    EG_createMeshMap        (bodyData *bodydata);
static double EG_angleAtVnormalPlane  (meshMap *qm, int vC, int v1, int v2 );
static int    EG_angleAtBoundaryVertex(meshMap *qm, int v, int *links, double *size ) ;
static int    vertexLinksToBounds (meshMap *qm, int vID );
static int    quadAngleOrientation    (meshMap *qm, int qID, int *validArea, int *ori, int *order, double *theta );
static int    quadAverageCoords       (meshMap *qm, int q, double *uv, double *p) ;
static void   averageCoordsMinusLinks (meshMap *qm,  int vc, int l1, int l2 ) ;
static void   averageCoordsUsingLinks ( meshMap *qm, int vc, int l1, int l2 );
static void   weightedAverage         (meshMap *qm, int vC );
static int    optimize_angles         (meshMap *qm, int nP, /*@unused@*/ /*@null@*/int *pList, int fullRegularization);
static int    moveDistAway            (meshMap *qm, double size, double *uv0, double *uv1 ) ;
/* REGULARIZATION FUNCTIONS */
/* FUNDAMENTAL */
static int  EG_swappingOperation      (meshMap *qm, quadGroup qg, int swap);
static int  EG_splittingOperation     (meshMap *qm, int vC, int vL, int vR);
static int  EG_mergeVertices          (meshMap *qm, int qC, int centre);

/**************/
static int  EG_cleanNeighborhood      (meshMap *qm, int qID, int transfer, int *activity );
static int  EG_swap                   (meshMap *qm, int qID, int *activity);
static int  EG_doubleSwap             (meshMap *qm, quadGroup qg, int forcing, int *activity );
static int  EG_collapse               (meshMap *qm, int  vID, int *activity);
static int  EG_forceCollapse          (meshMap *qm, int  vID, int *activity);
static int  EG_split                  (meshMap *qm, int  qID, int *activity);
static int  EG_removeDoublet          (meshMap *qm, int  vID ) ;
static int  EG_swapSplit              (meshMap *qm, quadGroup qg, int forcing, int *activity );
static int  EG_swapCollapse           (meshMap *qm, quadGroup qg, int forcing, int *activity );
static int  EG_doubleCollapse         (meshMap *qm, quadGroup qg, int forcing, int *activity );
static int  EG_doubleSplit            (meshMap *qm, quadGroup qg, int forcing, int *activity );
static int  EG_swapDoubleSplit        (meshMap *qm, quadGroup qg, int *activity );
static int  EG_swapDoubleCollapse     (meshMap *qm, quadGroup qg, int *activity );
static int  EG_transferValences       (meshMap *qm, int *qID, int try5533, int *transfering, int *activity  );
static int  EG_cleanQuad              (meshMap *qm, int qID,  int adj, int tansfer, int forcing, int *activity ) ;
static void EG_smoothQuads(meshMap *qm, int npass );
/* MESH DATA FUNCTIONS */
static int  EG_wvsData                (meshData *qm, char *buffer);
static void printQuadSpecs            (meshData *qm, int id) ;
static void printQuadGroup            (meshData *qm, quadGroup qg ) ;
static void printVertSpecs            (meshData *qm, int id);
static int    setValence              (meshData *qm, int vID );
static int  EG_buildStar              (meshData *qm, vStar **star, int vID );
static int  EG_nValenceCount          (meshData *qm, int q, int n );
static int  checkQuad                 (meshData *qm, int q );
static int  checkVertex               (meshData *qm, int v );
static int  EG_quadIsBoundary         (meshData *qm, int qID );
static int  EG_adjQtoPair             (meshData *qm, int qID, int v1, int v2, int *adj) ;
static void EG_commonVerts            (meshData *qm, int q1, int q2, int *v );
static int  EG_createQuadGroup        (meshData *qm, quadGroup *qg, int q0, int q1 ) ;
static int  getValence                (meshData *qm, int v );
static int  validSwap                 (meshData *qm, int v1, int v2 );
static int  validCollapse             (meshData *qm, int qID, int v );
/*********************/



static void
EG_smoothQuads(meshMap *qm, int npass )
{
	int           i, j, i0, i1, i2, i3, status, pass, nv, nq;
	double        qarea, sum, big, delta1 = 0, sums[2], x1[3], x2[3], xn[3];
	double        tAreaUV, tAreaXYZ, holdArea, results[18];
	static double wXYZ = 0.75;
	Node          *verts = NULL;
	nv    = qm -> mesh -> totVerts;
	nq    = qm -> mesh -> totQuads;
	verts = (Node *)EG_alloc ( nv * sizeof ( Node ));
	printMesh (qm, buffer, 0);
	if (verts == NULL ) return;
	for (i = 0; i < nv; i++) {
		if ( qm -> mesh -> vType[i] == -2 ) continue;
		verts[i].uv[0]  = qm -> mesh ->uvs  [ 2 * i    ];
		verts[i].uv[1]  = qm -> mesh ->uvs  [ 2 * i + 1];
		verts[i].xyz[0] = qm -> mesh ->xyzs [ 3 * i    ];
		verts[i].xyz[1] = qm -> mesh ->xyzs [ 3 * i + 1];
		verts[i].xyz[2] = qm -> mesh ->xyzs [ 3 * i + 2];
	}
	/* outer iteration -- pass 1 (uv only) */
	npass = 0;
	for (i = 0; i < nv; i++) {

		/* initialize deltas */
		for (j = 0; j < nv; j++) {
			verts[j].duv[0] = verts[j].duv[1] = 0.0;
			if ( qm -> mesh -> vType[j] >= 0 ) verts[j].area = -1.0;
			else verts[j].area = 0.0 ;
		}
		/* calculate and distribute change */
		for (j = 0; j < nq; j++) {
			if ( checkQuad ( qm -> mesh, j + 1 ) != EGADS_SUCCESS ) continue;
			i0    = qm -> mesh-> quadIdx [ 4 * j    ] - 1;
			i1    = qm -> mesh-> quadIdx [ 4 * j + 1] - 1;
			i2    = qm -> mesh-> quadIdx [ 4 * j + 2] - 1;
			i3    = qm -> mesh-> quadIdx [ 4 * j + 3] - 1;
			printf(" QUAD %d is %d %d %d %d\n ", j + 1, i0, i1, i2, i3 );
			qarea = -(AREA2D(verts[i0].uv, verts[i1].uv, verts[i2].uv) +
					AREA2D(verts[i0].uv, verts[i2].uv, verts[i3].uv));
			if (qarea <= 0.0) {
#ifdef DEBUG
				printf(" Quad %d: NonPositive Area = %le  %d %d %d %d\n",
						j, qarea, i0, i1, i2, i3);
				printf("          %lf %lf  %lf %lf  %lf %lf  %lf %lf\n",
						verts[i0].uv[0],verts[i0].uv[1], verts[i1].uv[0],verts[i1].uv[1],
						verts[i2].uv[0],verts[i2].uv[1], verts[i3].uv[0],verts[i3].uv[1]);
#endif
				qarea = -qarea;
			}

			sum = qarea*(verts[i0].uv[0] + verts[i1].uv[0] + verts[i2].uv[0] +
					verts[i3].uv[0])/4.0;
			verts[i0].duv[0] += sum;
			verts[i1].duv[0] += sum;
			verts[i2].duv[0] += sum;
			verts[i3].duv[0] += sum;
			sum = qarea*(verts[i0].uv[1] + verts[i1].uv[1] + verts[i2].uv[1] +
					verts[i3].uv[1])/4.0;
			verts[i0].duv[1] += sum;
			verts[i1].duv[1] += sum;
			verts[i2].duv[1] += sum;
			verts[i3].duv[1] += sum;

			if (verts[i0].area >= 0.0) verts[i0].area += qarea;
			if (verts[i1].area >= 0.0) verts[i1].area += qarea;
			if (verts[i2].area >= 0.0) verts[i2].area += qarea;
			if (verts[i3].area >= 0.0) verts[i3].area += qarea;
		}
		/* update distributions */
		big = 0.0;
		for (j = 0; j < nv; j++) {
			if (qm -> mesh -> vType[j] != -1 || verts[j].area <= 0.0) continue;
			sums[0] = verts[j].duv[0]/verts[j].area;
			sums[1] = verts[j].duv[1]/verts[j].area;
			sum     = fabs(sums[0] - verts[j].uv[0]);
			if (big < sum) big = sum;
			sum     = fabs(sums[1] - verts[j].uv[1]);
			if (big < sum) big = sum;
			verts[j].uv[0] = sums[0];
			verts[j].uv[1] = sums[1];
		}
		if (i == 0) {
			delta1 = big;
			if (delta1 == 0.0) break;
		} else {
			if (big/delta1 < 1.e-3) break;
		}
	}

	for (j = 0; j < nv; j++) {
		if (qm -> mesh -> vType[j] != -1 ) continue;
		status = EG_evaluate(qm -> face, verts[j].uv, results);
		if (status != EGADS_SUCCESS)
			return;
		verts[j].xyz[0] = results[0];
		verts[j].xyz[1] = results[1];
		verts[j].xyz[2] = results[2];
	}
	for (j = 0; j < nv; j++) {
		if (qm -> mesh -> vType[j] == -2 )continue;
		qm -> mesh -> uvs  [ 2 * j    ] = verts[j].uv[0];
		qm -> mesh -> uvs  [ 2 * j + 1] = verts[j].uv[1];
		qm -> mesh -> xyzs [ 3 * j    ] = verts[j].xyz[0];
		qm -> mesh -> xyzs [ 3 * j + 1] = verts[j].xyz[1];
		qm -> mesh -> xyzs [ 3 * j + 2] = verts[j].xyz[2];
	}

	/* pseudo non-linear loop */
	for (pass = 0; pass < npass; pass++) {
		/* get xyz */
		for (j = 0; j < nv; j++) {
			if (qm -> mesh -> vType[j] != -1 ) continue;
			status = EG_evaluate(qm -> face, verts[j].uv, results);
			if (status != EGADS_SUCCESS)
				return;
			verts[j].xyz[0] = results[0];
			verts[j].xyz[1] = results[1];
			verts[j].xyz[2] = results[2];
		}
		for (j = 0; j < nv; j++) {
			if (qm -> mesh -> vType[j] == -2 )continue;
			qm -> mesh -> uvs  [ 2 * j    ] = verts[j].uv[0];
			qm -> mesh -> uvs  [ 2 * j + 1] = verts[j].uv[1];
			qm -> mesh -> xyzs [ 3 * j    ] = verts[j].xyz[0];
			qm -> mesh -> xyzs [ 3 * j + 1] = verts[j].xyz[1];
			qm -> mesh -> xyzs [ 3 * j + 2] = verts[j].xyz[2];
		}
		printMesh (qm, buffer, 0);

		tAreaUV = tAreaXYZ = 0.0;
		for (j  = 0; j < nq; j++) {
			if ( checkQuad ( qm -> mesh, j + 1 ) != EGADS_SUCCESS ) continue;
			i0    = qm -> mesh-> quadIdx [ 4 * j    ] - 1;
			i1    = qm -> mesh-> quadIdx [ 4 * j + 1] - 1;
			i2    = qm -> mesh-> quadIdx [ 4 * j + 2] - 1;
			i3    = qm -> mesh-> quadIdx [ 4 * j + 3] - 1;
			holdArea  = -(AREA2D(verts[i0].uv, verts[i1].uv, verts[i2].uv) +
					AREA2D(verts[i0].uv, verts[i2].uv, verts[i3].uv));
			if (holdArea < 0.0) holdArea = -holdArea;
			tAreaUV  += holdArea;
			x1[0]     = verts[i1].xyz[0] - verts[i0].xyz[0];
			x2[0]     = verts[i2].xyz[0] - verts[i0].xyz[0];
			x1[1]     = verts[i1].xyz[1] - verts[i0].xyz[1];
			x2[1]     = verts[i2].xyz[1] - verts[i0].xyz[1];
			x1[2]     = verts[i1].xyz[2] - verts[i0].xyz[2];
			x2[2]     = verts[i2].xyz[2] - verts[i0].xyz[2];
			CROSS(xn, x1, x2);
			holdArea  = DOT(xn, xn);
			if (holdArea <  0.0) holdArea = -holdArea;
			tAreaXYZ += holdArea;
			x1[0]     = verts[i3].xyz[0] - verts[i0].xyz[0];
			x1[1]     = verts[i3].xyz[1] - verts[i0].xyz[1];
			x1[2]     = verts[i3].xyz[2] - verts[i0].xyz[2];
			CROSS(xn, x2, x1);
			holdArea  = DOT(xn, xn);
			if (holdArea <  0.0) holdArea = -holdArea;
			tAreaXYZ += holdArea;
		}
#ifdef DEBUG
		printf(" ** %d   Areas = %le  %le **\n", pass, tAreaUV, tAreaXYZ);
#endif
		printf(" OUTER ITERATION MIX PASS  it %d \n ", pass);
		/* outer iteration -- pass 2 (mix) */
		for (i = 0; i < nv; i++) {

			/* initialize deltas */
			for (j = 0; j < nv; j++) {
				verts[j].duv[0] = verts[j].duv[1] = 0.0;
				if ( qm -> mesh -> vType[j] >= 0 ) verts[j].area = -1.0;
				else verts[j].area = 0.0 ;
			}
			/* calculate and distribute change */
			for (j = 0; j < nq; j++) {
				if ( checkQuad ( qm -> mesh, j + 1 ) != EGADS_SUCCESS ) continue;
				i0    = qm -> mesh-> quadIdx [ 4 * j    ] - 1;
				i1    = qm -> mesh-> quadIdx [ 4 * j + 1] - 1;
				i2    = qm -> mesh-> quadIdx [ 4 * j + 2] - 1;
				i3    = qm -> mesh-> quadIdx [ 4 * j + 3] - 1;
				qarea = -(AREA2D(verts[i0].uv, verts[i1].uv, verts[i2].uv) +
						AREA2D(verts[i0].uv, verts[i2].uv, verts[i3].uv));
				if (qarea <= 0.0) {
					qarea = -qarea;
#ifdef DEBUG
					printf(" Quad %d: Neg Area = %le\n", j, qarea);
#endif
				}
				qarea *= (1.0-wXYZ)/tAreaUV;

				x1[0]     = verts[i1].xyz[0] - verts[i0].xyz[0];
				x2[0]     = verts[i2].xyz[0] - verts[i0].xyz[0];
				x1[1]     = verts[i1].xyz[1] - verts[i0].xyz[1];
				x2[1]     = verts[i2].xyz[1] - verts[i0].xyz[1];
				x1[2]     = verts[i1].xyz[2] - verts[i0].xyz[2];
				x2[2]     = verts[i2].xyz[2] - verts[i0].xyz[2];
				CROSS(xn, x1, x2);
				holdArea  = DOT(xn, xn);
				if (holdArea <  0.0) holdArea = -holdArea;
				qarea    += holdArea*wXYZ/tAreaXYZ;
				x1[0]     = verts[i3].xyz[0] - verts[i0].xyz[0];
				x1[1]     = verts[i3].xyz[1] - verts[i0].xyz[1];
				x1[2]     = verts[i3].xyz[2] - verts[i0].xyz[2];
				CROSS(xn, x2, x1);
				holdArea  = DOT(xn, xn);
				if (holdArea <  0.0) holdArea = -holdArea;
				qarea    += holdArea*wXYZ/tAreaXYZ;

				sum = qarea*(verts[i0].uv[0] + verts[i1].uv[0] + verts[i2].uv[0] +
						verts[i3].uv[0])/4.0;
				verts[i0].duv[0] += sum;
				verts[i1].duv[0] += sum;
				verts[i2].duv[0] += sum;
				verts[i3].duv[0] += sum;
				sum = qarea*(verts[i0].uv[1] + verts[i1].uv[1] + verts[i2].uv[1] +
						verts[i3].uv[1])/4.0;
				verts[i0].duv[1] += sum;
				verts[i1].duv[1] += sum;
				verts[i2].duv[1] += sum;
				verts[i3].duv[1] += sum;

				if (verts[i0].area >= 0.0) verts[i0].area += qarea;
				if (verts[i1].area >= 0.0) verts[i1].area += qarea;
				if (verts[i2].area >= 0.0) verts[i2].area += qarea;
				if (verts[i3].area >= 0.0) verts[i3].area += qarea;
			}
			/* update distributions */
			big = 0.0;
			for (j = 0; j < nv; j++) {
				if (qm -> mesh -> vType[j] != -1 || verts[j].area <= 0.0) continue;
				sums[0] = verts[j].duv[0]/verts[j].area;
				sums[1] = verts[j].duv[1]/verts[j].area;
				sum     = fabs(sums[0] - verts[j].uv[0]);
				if (big < sum) big = sum;
				sum     = fabs(sums[1] - verts[j].uv[1]);
				if (big < sum) big = sum;
				verts[j].uv[0] = sums[0];
				verts[j].uv[1] = sums[1];
			}
			if (i == 0) {
				delta1 = big;
				if (delta1 == 0.0) break;
			} else {
				if (big/delta1 < 1.e-3) break;
			}
		}
	}
	printf(" LAST WRITE \n ");
	EG_free ( verts );
	printMesh (qm, buffer, 0);
	return;
}












static double EG_segment ( meshMap *qm, double *uv1, double *uv2 ) {
	int i, j, n = 100, stat ;
	double p1[18], p2[18], dist[2], uvEps[2], seg = 0.0, totArc = 0.0, dt ;
	stat    = EG_evaluate (qm -> face, uv1, p1);
	if ( stat != EGADS_SUCCESS ) return -1.0;
	dist[0] = uv2[0] - uv1[0];
	dist[1] = uv2[1] - uv1[1];
	for ( i = 0; i < n; i++ ) {
		dt  = (double ) ( i + 1 ) / (double ) n;
		uvEps[0] = uv1[0] + dt * dist[0];
		uvEps[1] = uv1[1] + dt * dist[1];
		stat     = EG_evaluate (qm -> face, uvEps, p2 );
		if ( stat != EGADS_SUCCESS ) return -1.0;
		seg      = 0.0;
		for ( j = 0 ; j < 3; j++ ) {
			seg  += (p1[j] - p2[j] ) * (p1[j] - p2[j] );
			p1[j] = p2[j];
		}
		totArc += sqrt ( seg );
	}
	return totArc;
}

/* v should be dimension 3: v[0] = n common verts: v[1] (v[2]) idx*/
static void EG_commonVerts ( meshData *mesh, int q1, int q2, int *v ) {
	int i, j, k;
	v[0] = 0;
	if (checkQuad ( mesh, q1 ) != EGADS_SUCCESS ){
		printf(" EG_commonVerts you are trying to operate on a bad quad: %d --> %d!!\n", q1, checkQuad ( mesh, q1 ) );
		printQuadSpecs (mesh, q1 );
		v[0] = -1;
	}
	if (  checkQuad ( mesh, q2 ) != EGADS_SUCCESS ) {
		printf(" EG_commonVerts you are trying to operate on a bad quad: %d --> %d!!\n", q2, checkQuad ( mesh, q2 ) );
		printQuadSpecs (mesh, q2 );
		v[0] = -1;
	}
	if ( q1 == q2 ) {
		for ( k = j = 0 ; j < 4 ; j++ )
			v[++k] = mesh -> quadIdx [ 4 * ( q2 - 1 ) + j ];
		return;
	}
	for ( v[0] = k = i = 0 ; i < 4 ; i++ ) {
		for ( j = 0 ; j < 4 ; j++ ) {
			if (mesh -> quadIdx [ 4 * ( q1 - 1 ) + i ] ==
					mesh -> quadIdx [ 4 * ( q2 - 1 ) + j ] )
				v[++k] = mesh -> quadIdx [ 4 * ( q2 - 1 ) + j ];
		}
	}
	v[0] = k;
}

static int EG_quadVertIdx ( meshData *mesh, int q, int v ) {
	int i = 0 ;
	for ( i = 0 ; i < 4; i++ )
		if ( mesh -> quadIdx [ 4 * ( q - 1 ) + i ] == v ) return i;
	return -1;
}

static int getValence ( meshData *mesh, int v )  {
	int val;
	val = checkVertex (mesh, v );
	if ( val != EGADS_SUCCESS ) return val;
	val =     mesh -> valence [ v - 1 ][0];
	if      ( mesh -> vType   [ v - 1 ] == 0 ) val += 2;
	else if ( mesh -> vType   [ v - 1 ] == 3 ) val++;
	else if ( mesh -> vType   [ v - 1 ] == 5 ) val--;
	return val;
}

static int setValence ( meshData *mesh, int vID ) {
	int i;
	vStar *star = NULL;
	i = checkVertex ( mesh, vID );
	if ( i != EGADS_SUCCESS ) return i;
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
	return checkVertex ( mesh, vID );
}

static int EG_nValenceCount ( meshData *mesh, int q, int n ) {
	int i, count, val;
	for ( count = i = 0 ; i < 4; i++ ) {
		val = getValence ( mesh, mesh ->quadIdx [ 4 * ( q - 1) + i ] );
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
	if (  mesh ->vType [ v - 1 ] == -2 ) return EGADS_EMPTY;
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
		if ( mesh ->quadIdx [ 4 * q  + i] == -2 ) {
			printf("checkVertex::  Vertex for v %d has associated quad %d but  quad is empty !! \n ", v, q + 1);
			return EGADS_EMPTY;
		}
		if ( mesh ->quadIdx [ 4 * q  + i] ==  v )   return EGADS_SUCCESS;
	}
	printf("checkVertex::  Vertex for v %d has associated quad %d but !! \n ", v, q + 1);
	printQuadSpecs ( mesh, q + 1);
	return EGADS_INDEXERR;
}

static int checkQuad ( meshData *mesh, int q ) {
	int i, j;
	if ( q <= 0 || q > mesh ->totQuads ) {
		printf("checkQuad q = %d but MAX = %d !! \n", q, mesh ->totQuads );
		return EGADS_INDEXERR;
	}
	for ( i = 0 ; i < 4; i++) {
		j = checkVertex (mesh, mesh ->quadIdx [ 4 * ( q - 1) + i] );
		if (j != EGADS_SUCCESS ) return j;
	}
	return EGADS_SUCCESS;
}

static int checkMesh(meshMap *qm) {
	int stat, i, j, k,  val1, val2, v1, v2;
	for ( i = 0 ; i < qm -> mesh ->totVerts; i++ ) {
		if (qm -> mesh -> vType[i] == -2 )  continue;
		if (qm -> mesh -> vType[i] != 0 && qm -> mesh -> valence[i][0] == 2 ) {
			j = EG_removeDoublet (qm, i + 1 );
			if (qm -> mesh -> vType[i] == -2 )  continue;
			printQuadSpecs ( qm -> mesh , qm -> mesh -> valence [i][1] );
			fprintf(stderr," checkMesh is invalid: we have a  doublet at vertex %d type %d\n", i + 1, qm -> mesh -> vType[i] );
			snprintf(buffer,500,"InvalidMesh_%d", qm -> fID);
			printMesh(qm, buffer, 1);
			return EGADS_GEOMERR;
		}
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
			v1 = qm -> mesh  -> quadAdj [ 4 * i  + j] - 1;
			if ( v1 < 0 ) continue;
			val1 = -1;
			for ( k = 0 ; k < 4; ++k ) {
				if ( qm -> mesh  -> quadAdj [ 4 * v1  + k] == i + 1 ) {
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
		if ( checkQuad ( mesh, i + 1 ) == EGADS_SUCCESS ) qSum++;
	for ( vSum =  i = 0 ; i < mesh -> totVerts; i++ ) {
		if ( mesh -> vType[i] == -2 ) continue;
		vSum2++;
		if ( getValence ( mesh, i + 1 ) != 4) vSum++;
	}
	printf( " Current mesh has total quads %d and %d are irregular ( %d TOTAL )\n",
			qSum, vSum, vSum2);
	*nQ = qSum;
	*nI = vSum;
	*nV = vSum2;
}

/* Assuming qID is collapsing through v */
static int validCollapse ( meshData *mesh, int qID, int v ) {
	int j, k, kk, id, link, aux, aux2 ;
	id = EG_quadVertIdx ( mesh, qID, v);
	for ( j = 0 ; j < 2; j++ ) {
		link = mesh ->quadIdx [ 4 * ( qID - 1 ) + (id + 2 * j + 1)%4];
		printVertSpecs(mesh, link);
		if      (mesh ->vType [link - 1]  >  0 && mesh ->valence[link - 1][0] <= 3 ) return 0;
		else if (mesh ->vType [link - 1] ==  0 && mesh ->valence[link - 1][0] <= 2 ) return 0;
		else if (mesh ->vType [link - 1] == -1 && mesh ->valence[link - 1][0] == 3 ) {
			for ( k = 0 ; k < 3; k++ ) {
				aux = mesh ->valence [ link - 1 ][2 + k] - 1;
				if ( EG_quadVertIdx ( mesh, qID,  aux + 1 ) >= 0 ) continue;
				if ( mesh ->valence[ aux ][0] == 3 ) {
					if (  mesh ->vType[aux] > 0 ) return 0;
					else if ( mesh ->vType[aux] == -1 ) {
						for ( kk = 0 ; kk < 3; kk++ ) {
							aux2 = mesh ->valence [ aux ][2 + kk] - 1;
							if ( aux2 + 1 == link ) continue;
							if ( mesh ->valence[ aux2 ][0] == 3 && mesh ->vType[aux] != 0 ) return 0;
						}
					}
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
		if      (checkVertex ( mesh, vs[i] + 1 ) != EGADS_SUCCESS )                return 0;
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
	i = checkQuad ( mesh, qID );
	if ( i != EGADS_SUCCESS)  return i;
	for (i  = 0 ; i < 4; i++ ) {
		v   = mesh ->quadIdx [ 4 * ( qID - 1) + i ];
		q   = mesh ->quadAdj [ 4 * ( qID - 1) + i ];
		if ( mesh ->vType[v - 1] >=  0 ) bV = 1;
		if ( q                   == -1 ) bQ = 1;
	}
	if ( bQ == 1 ) return 1;
	if ( bV == 1 ) return 2;
	else return 0;

}

static int quadAverageCoords(meshMap *qm, int q, double *uv, double *xyz ) {
	int i, v;
	double eval[18], centre[2];
	centre[0] = 0.0; centre[1] = 0.0;
	i = checkQuad ( qm -> mesh, q );
	if ( i != EGADS_SUCCESS ) {
		printf(" In quadAverageCoords  stat for quad %d  is %d !!\n", q, i );
		return i;
	}
	for (i = 0; i < 4; i++ ) {
		v      = qm -> mesh -> quadIdx[4*( q - 1) + i] - 1;
		centre[0] += 0.25 * qm -> mesh -> uvs [2 * v    ];
		centre[1] += 0.25 * qm -> mesh -> uvs [2 * v + 1];
	}
	i     = EG_evaluate(qm->face, centre, eval);
	if  ( i != EGADS_SUCCESS ) {
		printf(" In quadAverageCoords  stat for quad %d  is %d !!\n", q, i );
		return i;
	}
	for ( i = 0 ; i < 3; i++) xyz[i] = eval[i];
	uv[0] = centre[0]; uv[1] = centre[1];
	return EGADS_SUCCESS;
}

static int EG_createQuadGroup ( meshData *mesh, quadGroup *qg, int q0, int q1 ) {
	int i,  ids[2], piv = 0, aux, vaux[6], common[3], centre, c1 = 0, c2 = 0;
	qg -> q[0]     = q0;
	qg -> q[1]     = q1;
	for ( i = 0 ; i < 6; i++) {
		qg -> verts[i] = -1;
		qg -> vals [i] = -1;
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
	if ( mesh -> quadAdj [ 4 * ( qg ->q[0] - 1 ) + ids[0] ] == qg ->q[1] ) piv = 1;
	for ( i = 0 ; i < 4; i++ )
		qg -> verts[i] = mesh -> quadIdx [ 4 * (qg ->q[piv] - 1) + (ids[piv] + i ) % 4 ];
	aux = ( piv + 1 ) % 2;
	qg -> verts[4] = mesh -> quadIdx [ 4 * (qg ->q[aux] - 1)  + ( ids[ aux ] + 2 ) % 4 ];
	qg -> verts[5] = mesh -> quadIdx [ 4 * (qg ->q[aux] - 1)  + ( ids[ aux ] + 3 ) % 4 ];
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
	int i, act = 0, stat, i3, i5, q3, f = 0;
	*activity  = 0;
	if (qg.vals[0] * qg.vals[3] == 16 ) {
		for ( i3 = 0 ; i3 < 6; i3++ ) if ( qg.vals[i3] == 3 ) break;
		if  (qg.vals[i3] != 3 ) return EGADS_SUCCESS;
		i5 = ( i3 + 1 ) % 6; if ( i5 % 3 == 0 ) i5 = ( i3 + 5 ) % 6;
		if (                    forcing == 0 && (qg.vals[i5]!= 5 ||
				(qg.vals[(i5 + 3) %6]   != 5 && qg.vals [ (i3 + 3) % 6 ] != 3 ) ) ) return EGADS_SUCCESS;
		else if ( forcing == 1) {
			if ( qm -> extraQuads < 0 ) return EGADS_SUCCESS;
			if ( qg.vals[(i3 + 3)%6] == 3 )      i5 = (i3 + 3)%6;
			else if ( qg.vals[(i3 + 2)%6] == 5 ) i5 = (i3 + 2)%6;
			else if ( qg.vals[(i3 + 4)%6] == 5 ) i5 = (i3 + 4)%6;
			else return EGADS_SUCCESS;
			printf(" found something:: i5 %d \n ", i5 );
		}
	}
	else if (qg.vals[0] * qg.vals[3] == 12 ) {
		i3 = 0 ;            if ( qg.vals[3 ] == 3 ) i3 = 3;
		i5 = ( i3 + 1 )% 6; if ( qg.vals[i5] != 5 ) i5 = ( i3 + 5 )% 6;
		if ( qg.vals[i5] != 5 || qg.vals [( i5 + 3 ) % 5] < 5 ) return EGADS_SUCCESS;
	}
	else return EGADS_SUCCESS;
	if ( forcing == 1 && qg.vals[0] * qg.vals[3] == 16 ) {
		f = 1;
		q3 = 0; if ( i3 >= 3 ) q3 = 1;
	}else {
		q3 = 0; if ( i5 >= 3 ) q3 = 1;
	}
	printf(" DOUBLE COLLAPSE %d %d \n ", qg.q[0], qg.q[1] );
	for ( i = 0 ; i < 2; i++ ) {
		if (checkQuad ( qm -> mesh, qg.q[( q3 + i ) %2] ) != EGADS_SUCCESS )continue;
		stat       = EG_forceCollapse (qm, qg.q[( q3 + i ) %2], &act);
		if ( act > 0 && f == 1 ) qm -> extraQuads--;
		*activity += act;
		if ( stat != EGADS_SUCCESS ) {
			printf(" EG_doubleCollapse: EG_forceCollapse went %d !! \n", stat);
			return stat;
		}
	}
	return stat;
}


static int
EG_swapDoubleCollapse (meshMap *qm, quadGroup qg, int *activity ) {
	int  k, swap = 1, id, j, stat;
	*activity  = 0;
	if (    qg.vals[0] * qg.vals[3] != 20 ||
			qg.vals[2] * qg.vals[4] != 9  ||
			validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0) return EGADS_SUCCESS;
	if ( qg.vals[swap] * qg.vals[(swap+3)%6] != 12 ) swap = 2;
	if ( qg.vals[swap] * qg.vals[(swap+3)%6] != 12 ) return EGADS_SUCCESS;
	stat       = EG_swappingOperation(qm, qg, swap);
	*activity  = 1;
	if ( stat != EGADS_SUCCESS ) {
		printf( "In EG_swapDoubleCollapse: Swapping went %d \n", stat );
		return stat;
	}
	for ( k = 0; k < 2; k++ ) {
		id  = EG_quadVertIdx ( qm -> mesh, qg.q[k], qg.verts[3] );
		if ( id >= 0 ) break;
	}
	if ( id < 0 ) {
		printf(" I can't find vertex %d in quads!!! \n", qg.verts[3]);
		printQuadGroup (qm -> mesh, qg);
		return EGADS_INDEXERR;
	}
	stat       = EG_collapse ( qm, qg.q[k], &j );
	if ( stat != EGADS_SUCCESS )
		printf(" In EG swapDoubleCollapse: I failed to collapse after swapping! s = %d act = %d \n ", stat, j );
	return stat;
}


static int
EG_swapDoubleSplit (meshMap *qm, quadGroup qg, int *activity ) {
	int  i5, q, i55, i3, val3, i0, v30[2], i,  stat, adj[2], q0[2];
	*activity  = 0;
	for ( i3 = 0 ; i3 < 6; i3++)
		if ( qg.vals[i3] == 3 ) break;
	if ( qg.vals[i3] != 3 || i3 % 3 == 0 ) return EGADS_SUCCESS;
	i5  =(i3 + 3)% 6;
	if ( qg.vals[i5] < 5 ) return EGADS_SUCCESS;
	i55 = (i5 + 1)%6;
	if ( qg.vals[i55] < 5 ) {
		i55 = (i5 + 5)%6 ;
		if ( qg.vals[i55] != 5 ) return EGADS_SUCCESS;
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
	v30[0]     = qg.verts[i0];
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
	printQuadGroup (qm -> mesh, qg );
	for ( i0 = 0 ; i0 < 6; i0++) if ( qg.verts[i0] == v30[0] ) break;
	*activity  = 1;
	stat       = EG_swappingOperation (qm, qg, i0 );
	if ( stat != EGADS_SUCCESS ) {
		printf(" EG_swapDoubleSplit error at swap: %d !!\n ", stat );
		return stat;
	}
	i = 0 ; if ( EG_quadVertIdx (qm -> mesh, q0[0], val3) < 0 ) i = 1;
	stat = EG_adjQtoPair (qm -> mesh, q0[i], v30[0], v30[1], adj );
	if ( stat != EGADS_SUCCESS || adj[1] == -1 ) {
		if ( stat != EGADS_SUCCESS ) printf(" EG_swapDoubleSplit after swapping adjacent to pair %d !!\n ", stat);
		return stat;
	}
	stat = EG_createQuadGroup (qm -> mesh, &qg,  q0[i], adj[1]);
	if ( stat != EGADS_SUCCESS ) {
		printf("Inside EG_swapDoubleSplit: before splitting EG_createQuadGroup stat %d\n ", stat );
		printQuadGroup (qm -> mesh, qg );
		return stat;
	}
	return EG_doubleSplit (qm, qg, 0, &i);
}



static int EG_transferValences ( meshMap *qm, int *qID, int try5533, int *transfering, int *activity  ) {
	int i, j, swap = 0, stat, min;
	int links[5], qAux[2];
	quadGroup qg;
	*activity   = 0;
	stat        = checkQuad ( qm -> mesh, qID[0] ) ;
	if ( stat  != EGADS_SUCCESS ) {
		printf(" EG_transferValences for quad %d is bad quad %d \n ", qID[0], stat );
		return stat;
	}
	if (*transfering == 0 ) {
		stat       = EG_cleanQuad ( qm, qID[0], 1, 0, try5533 , &(*activity) );
		if ( stat != EGADS_SUCCESS || *activity > 0 ) {
			if ( stat != EGADS_SUCCESS)
				printf("EG_transferValences: EG_cleanQuad %d --> %d!!\n", qID[0], stat);
			if ( checkQuad (qm -> mesh, qID[0] ) != EGADS_SUCCESS)  qID[0] = -1;
			return stat;
		}
	}
	if ( EG_nValenceCount (qm -> mesh, qID[0], 4) == 4 ) return EGADS_SUCCESS;
	for ( swap  = j = 0 ; j < 4; j++ ) {
		qg.q[0] = qID[0];
		qg.q[1] = qm -> mesh -> quadAdj[ 4 * ( qID[0] - 1 ) + j];
		if ( qg.q[1] < 0 || qg.q[1] == qID[1] ) continue;
		if (*transfering == 0 ) {
			stat           = EG_cleanQuad ( qm, qg.q[1], 1, 0,  try5533 , &(*activity) );
			if ( stat     != EGADS_SUCCESS || *activity > 0 ) {
				if ( stat != EGADS_SUCCESS) printf(" EG_TransferValence clean quad %d !!\n ", stat );
				qID[0]     = -1;
				return stat;
			}
		}
		qg.q[0]    = qID[0];
		qg.q[1]    = qm -> mesh -> quadAdj[ 4 * ( qID[0] - 1 ) + j];
		stat       = EG_createQuadGroup  (qm -> mesh, &qg, qg.q[0], qg.q[1] );
		if ( stat != EGADS_SUCCESS  ) {
			printf(" Inside EG_transferValences EG_createQuadGroup %d !!\n", stat );
			printQuadGroup (qm -> mesh, qg );
			return stat;
		}
		if ( *transfering == 0 && qg.vals[0] * qg.vals[3] == 15 ) {
			j = 0;
			if ( qg.q[j] == qID[0] ) j = 1;
			qID[0] = qg.q[j]; qID[1] = qg.q[(j + 1)%2];
			if ( EG_quadIsBoundary(qm -> mesh, qg.q[1] ) == 1 && qm -> extraQuads >= 0 ) {
				stat = EG_forceCollapse (qm, qID[0], &(*activity));
				printMesh (qm, buffer,0);
				qID[0] = qID[1];
				qID[1] = -1;
				if ( stat != EGADS_SUCCESS ) printf("EG_transferValences: forceCollapse gave %d !!\n", stat);
				return stat;
			}
			else if  ( EG_quadIsBoundary   (qm -> mesh, qg.q[1] ) == 1 )
				return EG_transferValences (qm, qID, 0, &(*transfering),  &(*activity) );
		}
		if (validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0 ) continue;
		min     = 0;
		if      (qg.vals[0] * qg.vals[3] >= 20 ) min = 12;
		else if ( try5533 == 1 && *transfering == 0  ) {
			if      (qg.vals[0] * qg.vals[3] >= 25) min = 16;
			else if (qg.vals[0] * qg.vals[3] == 16) min =  9;
		}
		if ( min == 0 ) continue;
		for( i = 1 ; i <= 2; i++ )
			if (qg.vals[i] * qg.vals [3 + i] <= min ) swap = i;
		if ( swap == 0 ) continue;
		stat       = EG_swappingOperation (qm, qg, swap);
		*activity  = 1;
		printMesh(qm, buffer,0);
		if ( stat != EGADS_SUCCESS ) return stat;
		i = 0;
		if ( min == 9 || min == 16) *transfering = 1;
		else                        *transfering = 0;
		if (qID[1] != -1 &&
				checkQuad ( qm -> mesh, qID[1]) == EGADS_SUCCESS ) {
			if      ( qg.q[0] == qID[1] ) i = 1;
			else if ( qg.q[1] == qID[1] ) i = 0;
			else {
				if ( checkQuad ( qm -> mesh, qg.q[0]) != EGADS_SUCCESS ||
						checkQuad ( qm -> mesh, qID[1]) != EGADS_SUCCESS )
				{
					printf(" Quads are invalid %d =%d !!\n ",qg.q[0], checkQuad ( qm -> mesh, qg.q[0]));
					printf(" Quads are invalid %d =%d !!\n ",qID[1], checkQuad ( qm -> mesh, qID[1]));
					exit(1);
				}
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
		qID[0] = qg.q[i];
		qID[1] = qg.q[ ( i + 1 ) %2 ];
		if (*transfering == 1 ) {
			for ( j = 0 ; j < 2; j++ ) {
				i       = 1;
				qAux[0] = qID[ j       ];
				qAux[1] = qID[(j + 1)%2];
				stat    = EG_transferValences ( qm, qAux, 0, &i, &min);
				if ( stat != EGADS_SUCCESS ) {
					printf(" EG_transferValences: separating valences after forcing %d!!\n", stat);
					return stat;
				}
			}
			*transfering = 0;
		}
		qID[0] = qAux[0];
		qID[1] = qAux[1];
		if ( *activity > 0 ) break;
	}
	return stat;
}


static int EG_basicOperation (meshMap *qm, int qID, int type, int *activity ) {
	int stat = EGADS_SUCCESS;
	switch (type ) {
	case SWAP:
		stat     = EG_swap ( qm, qID, & (*activity));
		break;
	case COLLAPSE:
		stat  = EG_collapse ( qm, qID, & (*activity));
		break;
	case SPLIT:
		stat  = EG_split ( qm, qID, & (*activity));
		break;
	}
	if ( stat != EGADS_SUCCESS) printf(" EG_basicOperation %d  around %d --> %d!!\n", type, qID, stat);
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
	if ( stat != EGADS_SUCCESS) {
		printf(" EG_composeOperationOperation %d --> %d!!\n", type, stat);
		printQuadGroup (qm -> mesh, qg);
	}
	return stat;
}



static int EG_cleanQuad (meshMap *qm, int qID, int useAdj, int transfer, int forcing, int *activity ) {
	int stat, i, q, qadj, act = 0;
	int opBasic[3] = {COLLAPSE, SWAP, SPLIT};
	int opComp [7] = {SWAPCOLLAPSE, DOUBLECOLLAPSE, SWAPDOUBLECOLLAPSE, DOUBLESWAP,
			SWAPSPLIT,    DOUBLESPLIT,    SWAPDOUBLESPLIT};
	quadGroup qg;
	*activity      = 0;
	stat           = checkQuad ( qm -> mesh, qID  );
	if ( stat     != EGADS_SUCCESS ) {
		printf(" EG_cleanQuad %d stat %d\n", qID, stat );
		return stat;
	}
	if ( EG_nValenceCount (qm -> mesh, qID, 4) == 4 ) return EGADS_SUCCESS;
	if ( transfer == 0 ) {
		if ( EG_nValenceCount (qm -> mesh, qID, 6) > 0 || qm -> extraQuads < 0 ) swapInt( &opBasic[2], &opBasic[0]);
		for ( i = 0 ; i < 3; i++ ) {
			stat      = EG_basicOperation(qm, qID, opBasic[i], &act );
			*activity = act;
			if ( stat != EGADS_SUCCESS || act > 0 ) {
				if ( stat != EGADS_SUCCESS ) printf("EG_cleanQuad %d  stat in basic %d!! \n ",qID, stat);
				return stat;
			}
			if ( checkQuad ( qm -> mesh, qID  )   != EGADS_SUCCESS ) return EGADS_SUCCESS;
		}
	}
	if ( useAdj == 0  ) return EGADS_SUCCESS;
	if ( qm -> extraQuads < 0 )
		for ( i = 0 ; i < 3; i++) swapInt( &opComp[i], &opComp[i+4]);
	stat           = checkQuad ( qm -> mesh, qID  );
	if ( stat     != EGADS_SUCCESS ) {
		printf(" EG_cleanQuad %d checkQuad is stat %d!!\n", qID, stat );
		return stat;
	}
	if ( forcing == 1 && EG_nValenceCount (qm -> mesh, qID, 3) > 1 ) {
		//	fprintf(stderr,"cleanQuad forcing extra quads %d\n ", qm -> extraQuads );
		stat      = EG_forceCollapse (qm, qID, &(*activity ) );
		if ( *activity > 0 || stat != EGADS_SUCCESS ) {
			--qm -> extraQuads;
			if ( stat != EGADS_SUCCESS ) printf("EG_cleanQuad %d forceCollapse %d!! \n ",qID, stat);
			return stat;
		}
	}
	if ( transfer == 0 ) {
		for ( q  = 0 ; q < 4; q++) {
			qadj = qm -> mesh -> quadAdj [ 4 * ( qID -1 ) + q];
			if ( qadj == -1 ) continue;
			stat = EG_cleanQuad (qm, qadj, 0, 0, 0, &act);
			if ( act > 0 || stat != EGADS_SUCCESS ) {
				*activity = act;
				if ( stat != EGADS_SUCCESS )printf("Inside EG_cleanQuad: in cleaning neighbor quad:  stat %d\n ", stat );
				return stat;
			}
		}
	}
	for ( i  = 0 ; i < 7; i++ ) {
		if ( transfer == 1 && opComp[i] == DOUBLESWAP ) continue;
		for ( q  = 0 ; q < 4; q++) {
			if ( checkQuad ( qm -> mesh, qID  )   != EGADS_SUCCESS ) return EGADS_SUCCESS;
			qadj = qm -> mesh -> quadAdj [ 4 * ( qID -1 ) + q];
			if ( qadj == -1 ) continue;
			stat = EG_createQuadGroup (qm -> mesh, &qg, qID, qadj);
			if ( stat != EGADS_SUCCESS ) {
				printf("Inside EG_cleanQuad: EG_createQuadGroup stat %d\n ", stat );
				printQuadGroup (qm -> mesh, qg );
				return stat;
			}
			stat      = EG_composeOperation(qm, qg, opComp[i], forcing, &act );
			*activity = act;
			if ( stat != EGADS_SUCCESS || act > 0  ) return stat;
		}
	}
	return stat;
}

static int EG_createMeshMap(bodyData *bodydata)
{
	int         f, stat = 0, j, q, i, auxID, k, kk, kOK, len,  ntri, nquad, e4[4];
	int         *faceEdges = NULL ;
	const int   *tris, *tric, *ptype, *pindex;
	double      eval[18], angle;
	const double *xyzs, *uvs;
	int    qV[6]    = { 0, 1, 2, 5, 0, 1};
	int    qLoop[5] = { 0, 1, 2, 3, 0};
	FILE *fill = NULL;
	bodydata->qm = (meshMap**) EG_alloc(bodydata->nfaces*sizeof(meshMap*));
	if (bodydata->qm == NULL ) return  EGADS_MALLOC;
	faceEdges = (int *)EG_alloc ( bodydata -> nedges * sizeof ( int));
	if ( faceEdges == NULL ) return EGADS_MALLOC;
	for ( f = 0 ; f < bodydata->nfaces; ++f) {
		bodydata -> qm[f] = (meshMap*) EG_alloc(sizeof(meshMap));
		if ( bodydata -> qm[f] == NULL ) {
			printf("Create Quad Map: MALLOC error!! \n ");
			EG_free (faceEdges );
			return EGADS_MALLOC;
		}
		bodydata -> qm[f] -> mesh       = NULL;
		bodydata -> qm[f] -> backupMesh = NULL;
		bodydata -> qm[f] -> bestMesh   = NULL;
		bodydata -> qm[f] -> fID        = f + 1;
		// Edges associated to face //
		stat = EG_getTessFace(bodydata->tess, f + 1, &len,
				&xyzs, &uvs, &ptype, &pindex, &ntri,
				&tris, &tric);
		for ( i = 0 ; i < bodydata -> nedges; i++ ) faceEdges[i] = 0;
		for ( j = i = 0 ; i < len; i++ ) {
			if ( pindex[i] == -1 ) continue;
			j++;
			faceEdges[pindex[i] -1]++;
		}
		bodydata ->qm[f] -> boundaryVerts = j;
		for ( j = i = 0 ; i < bodydata -> nedges; i++ )
			if ( faceEdges[i] > 0 ) {
				if ( j < 4 ) e4[j] = faceEdges[i];
				j++;
			}
		printf(" Face %d is bounded by the following edges (TOTAL %d ) \n ", f + 1, j);
		nquad = (int)ntri/2;
		stat  = EG_allocMeshData ( &bodydata ->qm[f] -> mesh      ,2 * nquad,  2 * len);
		stat += EG_allocMeshData ( &bodydata ->qm[f] -> backupMesh,2 * nquad,  2 * len);
		stat += EG_allocMeshData ( &bodydata ->qm[f] -> bestMesh  ,2 * nquad,  2 * len);
		if ( stat != EGADS_SUCCESS ||
				bodydata ->qm[f] -> mesh       == NULL ||
				bodydata ->qm[f] -> backupMesh == NULL ||
				bodydata ->qm[f] -> bestMesh   == NULL ) {
			fprintf(stderr,"In createMeshMap  EG_allocMeshData = %d\n ", stat );
			EG_free(faceEdges );
			return stat;
		}
		bodydata->qm[f] -> extraQuads       = 0;
		bodydata->qm[f] -> face             = bodydata->faces[f];
		bodydata->qm[f] -> mesh -> totVerts = len;
		bodydata->qm[f] -> oriV             = len;
		bodydata->qm[f] -> mesh -> totQuads = nquad;
		bodydata->qm[f] -> oriQ             = nquad;
		for (j = 0; j < len; j++) {
			bodydata->qm[f] -> mesh  ->  valence[j][0]   = 0;
			bodydata->qm[f] -> mesh  ->  uvs   [2*j    ] = uvs   [2*j   ];
			bodydata->qm[f] -> mesh  ->  uvs   [2*j + 1] = uvs   [2*j +1];
			stat = EG_evaluate(bodydata->qm[f] -> face, &bodydata->qm[f] -> mesh  ->uvs[2*j], eval);
			bodydata->qm[f] -> mesh  ->  xyzs  [3*j    ] = eval[0];
			bodydata->qm[f] -> mesh  ->  xyzs  [3*j + 1] = eval[1];
			bodydata->qm[f] -> mesh  ->  xyzs  [3*j + 2] = eval[2];
			bodydata->qm[f] -> mesh  ->  vType [j      ] = ptype [j];
		}
		for (j = 0; j < nquad; j++)
			for ( k = 0; k < 4; ++k) {
				bodydata->qm[f] -> mesh  ->  quadIdx[4*j + k ] = tris[6*j + qV[k+1]];
			}
		for ( j = 0; j < nquad; j++) {
			kk  = 0;
			kOK = 0;
			q   = 0;
			while ( q < nquad )
			{
				if ( q == j ) {
					if ( j == nquad - 1 && kk < 4 ) {
						bodydata->qm[f] -> mesh  ->  quadAdj[4*j + kk++] = -1;
						if (kk == 4) break;
						q  = 0;
					}
					else q++;
				}
				if ( q == nquad ) break;
				for ( k = 0 ; k < 4; ++k ) {
					if((bodydata->qm[f] -> mesh-> quadIdx[4*j + qLoop[kk    ]] == bodydata->qm[f]-> mesh-> quadIdx[4*q + qLoop[k    ]] ||
							bodydata->qm[f] -> mesh-> quadIdx[4*j + qLoop[kk    ]] == bodydata->qm[f]-> mesh-> quadIdx[4*q + qLoop[k + 1]]
					)&&(bodydata->qm[f] -> mesh-> quadIdx[4*j + qLoop[kk + 1]] == bodydata->qm[f]-> mesh-> quadIdx[4*q + qLoop[k    ]] ||
							bodydata->qm[f] -> mesh-> quadIdx[4*j + qLoop[kk + 1]] == bodydata->qm[f]-> mesh-> quadIdx[4*q + qLoop[k + 1]]) )
					{
						bodydata->qm[f] -> mesh  ->  quadAdj[4*j + kk++] = q + 1;
						q   =  -1;
						kOK =   1;
						k   =   4;
						if (kk == 4)  q = nquad;
					}
				}
				if ( (kOK == 0) && (q >= nquad -1) ){
					bodydata->qm[f] -> mesh  ->  quadAdj[4*j + kk++] = -1;
					q                       = -1;
					if (kk == 4) break;
				}
				else  kOK = 0 ;
				q++;
			}
			if (kOK == 0 && kk < 4) {
				while (kk < 4)
					bodydata->qm[f] -> mesh  ->  quadAdj[4*j + kk++] = -1;
			}
		}
		for ( j = 0 ; j < nquad; ++j) {
			for ( k = 0 ; k < 4; ++k) {
				if ( bodydata->qm[f] -> mesh-> quadAdj[4*j +k] > j + 1  || bodydata->qm[f] -> mesh  -> quadAdj[4*j +k] == -1) {
					auxID = bodydata->qm[f] -> mesh  -> quadIdx [4*j + qLoop[k]] -1;
					bodydata->qm[f] -> mesh    -> valence[auxID][1] = j + 1;
					bodydata->qm[f] -> mesh    -> valence[auxID][2 + bodydata->qm[f] -> mesh  -> valence[auxID][0] ] = bodydata->qm[f] -> mesh  -> quadIdx[4*j+ qLoop[k+1]];
					++bodydata->qm[f] -> mesh  -> valence[auxID][0];
					auxID = bodydata->qm[f] -> mesh  -> quadIdx [4*j + qLoop[k+1]] -1;
					bodydata->qm[f] -> mesh    -> valence[auxID][1] = j + 1;
					bodydata->qm[f] -> mesh    -> valence[auxID][2 + bodydata->qm[f] -> mesh  -> valence[auxID][0] ] = bodydata->qm[f] -> mesh  -> quadIdx[4*j+ qLoop[k]];
					++bodydata->qm[f] -> mesh  -> valence[auxID][0];
				}
			}
		}
		snprintf(buffer,500,"gnuInit_%i",f+1);
		printMesh(bodydata->qm[f] , buffer, 1);
		snprintf(buffer,500,"wvsInit_%i.txt",f+1);
		EG_wvsData(bodydata->qm[f] -> mesh , buffer);
		snprintf(buffer, 500, "face_%d_edges", f + 1);
		fill = fopen ( buffer, "w");
		bodydata->qm[f] -> minsize = 1000;
		bodydata->qm[f] -> maxsize = 0.0;
		for ( j = 0 ; j < len; j++ ) {
			if ( bodydata->qm[f] -> mesh ->vType[j] == -1 ) continue;
			k = bodydata->qm[f] -> mesh -> valence[j][0];
			for ( i = 0 ; i < k; i++ ) {
				printf(" i %d / %d \n ", i + 1, k);
				auxID = bodydata->qm[f] -> mesh -> valence[j][2 + i] - 1;
				if ( bodydata->qm[f] -> mesh -> vType[auxID] == -1 ) continue;
				angle = EG_segment (bodydata->qm[f], &bodydata->qm[f] -> mesh -> uvs [ 2 *j ], &bodydata->qm[f] -> mesh -> uvs [ 2 *auxID ]);
				if ( angle > bodydata->qm[f] -> maxsize )bodydata->qm[f] -> maxsize = angle;
				if ( angle < bodydata->qm[f] -> minsize )bodydata->qm[f] -> minsize = angle;
			}
			fprintf(fill, "%lf %lf %lf %d\n", xyzs[3 * j ], xyzs[3 * j +1], xyzs[3 * j +2], j + 1);
			stat = EG_angleAtBoundaryVertex (bodydata -> qm[f], j + 1, e4, &angle);
			if ( stat != EGADS_SUCCESS || angle < EPS11 ) {
				EG_free(faceEdges);
				return EGADS_GEOMERR;
			}
			else if ( angle < 0.85 * PI ) {
				if ( bodydata->qm[f] -> mesh ->vType[j] > 0 ) {
					printf(" Node %d is now consider edge so that its regular valence is 3 \n", j + 1);
					bodydata->qm[f] -> mesh ->vType[j] = 0;
				}
			}
			else if ( angle < 1.25 * PI  ) {
				if ( bodydata->qm[f] -> mesh -> valence[j][0] >= 3 )
					bodydata->qm[f] -> mesh -> vType  [j] = 3;
			}
			else if ( angle < 1.85 * PI ){
				if ( bodydata->qm[f] -> mesh -> valence[j][0] >= 4 )
					bodydata->qm[f] -> mesh -> vType  [j] = 4 ;
			} else {
				if ( bodydata->qm[f] -> mesh -> valence[j][0] >= 5 )
					bodydata->qm[f] -> mesh -> vType  [j] = 5 ;
			}
			if ( ptype [j] != bodydata ->  qm[f] -> mesh -> vType  [j])
				printf(" ANGLE %lf Vertex %d in Face %d is now edge TYPE %d was %d \n ", angle, j + 1, f + 1,   bodydata->qm[f] -> mesh  -> vType  [j], ptype[j]);
		}
		fclose (fill);
		if ( stat != EGADS_SUCCESS ) stat = EGADS_SUCCESS;
		printf(" Mesh min max sizes are %f  %f  \n ",  bodydata->qm[f] -> minsize , bodydata->qm[f] -> maxsize );
		stat = checkMesh (bodydata->qm[f]);
	}
	EG_free(faceEdges);
	return stat;
}


static void EG_destroymeshMap(meshMap ***qm, int nfaces ) {
	int i ;
	if ( *qm  == NULL ) return ;
	for ( i = 0 ; i < nfaces; ++i) {
		if ( ( *qm)[i] ) {
			EG_freeMeshData (&((*qm)[i] -> mesh      ));
			EG_freeMeshData (&((*qm)[i] -> backupMesh));
			EG_freeMeshData (&((*qm)[i] -> bestMesh  ));
		}
		EG_free(( *qm)[i]);
	}
	EG_free (*qm );
}

static int
EG_adjQtoPair(meshData *mesh, int qID, int v1, int v2, int *adj) {
	int i, aux = -1;
	adj[0] = -1; adj[1] = -1;
	for ( i = 0 ; i < 4; ++i) {
		if ( mesh -> quadIdx[4*(qID - 1) + i ] == v1 ) adj[0] = i;
		if ( mesh -> quadIdx[4*(qID - 1) + i ] == v2 ) aux    = i;
		if ( aux != -1 && adj[0] != -1 ) break;
	}
	if ( aux == -1 || adj[0] == -1 ) return EGADS_SUCCESS;
	if      ( abs (adj[0] - aux ) == 3 ) adj[0] = 3;
	else if ( aux < adj[0]             ) adj[0] = aux;
	adj[1] = mesh -> quadAdj[4*(qID - 1) + adj[ 0 ] ];
	return EGADS_SUCCESS;
}

static int restoreMeshData(meshMap *qm, meshData *m1, meshData *m2 ) {
	int i, j, k;
	m1 -> totQuads = m2 -> totQuads;
	m1 -> totVerts = m2 -> totVerts;
	for ( j = 0 ; j < m2 -> totQuads; ++j) {
		for ( k = 0 ; k < 4; ++k) {
			m1 -> quadIdx[4*j + k] = m2 -> quadIdx[4*j+ k];
			m1 -> quadAdj[4*j + k] = m2 -> quadAdj[4*j+ k];
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
	return checkMesh (qm );
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
		if(qm -> mesh ->quadIdx[4*i] != -2 ) nQ++;
		else qRem++;
	}
	if ( vRem != qRem ) {
		if ( qRem > 0 ) {
			printf(" In resizeQm: I have %d removed vertices and %d quads!! they should match!!!!!\n ", vRem, qRem);
			return EGADS_INDEXERR;
		}
	}
	if ( vRem == 0 )return EGADS_SUCCESS;
	stat = restoreMeshData (qm, qm -> backupMesh, qm -> mesh );
	if ( stat != EGADS_SUCCESS) return stat;
	vpiv = (int*) EG_alloc(qm -> mesh -> totVerts * sizeof(int));
	qpiv = (int*) EG_alloc(qm -> mesh -> totQuads * sizeof(int));
	if (vpiv == NULL || qpiv == NULL) return EGADS_MALLOC;
	for ( j = i = 0 ; i < qm -> backupMesh -> totQuads; i++) {
		if(qm -> mesh -> quadIdx[4*i] != -2 ) {
			qpiv[i] = j;
			++j;
		} else qpiv[i] = -2;
	}
	for ( j = i = 0 ; i < qm -> backupMesh -> totVerts; i++) {
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
	for ( i = 0 ; i < qm -> backupMesh -> totQuads; i++){
		if ( qpiv[i] == -2) continue;
		for ( k = 0 ; k < 4; ++k) {
			qm -> mesh -> quadIdx[4*qpiv[i] + k ] = vpiv[ qm -> backupMesh -> quadIdx[4*i + k] - 1] + 1;
			if ( qm -> backupMesh -> quadAdj[4*i + k] == -1 ) qm -> mesh -> quadAdj[4*qpiv[i] + k ] = - 1;
			else qm -> mesh       -> quadAdj[4*qpiv[i] + k ] = qpiv[ qm -> backupMesh -> quadAdj[4*i + k] - 1] + 1;
		}
	}
	for ( i = 0 ; i < qm -> backupMesh -> totVerts; i++) {
		if ( vpiv[i] == -2 ) continue;
		j   = vpiv[i];
		qm -> mesh -> vType[j      ] = qm -> backupMesh -> vType[  i    ];
		qm -> mesh -> uvs  [2*j    ] = qm -> backupMesh -> uvs  [2*i    ];
		qm -> mesh -> uvs  [2*j + 1] = qm -> backupMesh -> uvs  [2*i + 1];
		stat                         = EG_evaluate(qm -> face, &qm -> mesh -> uvs[2*j], eval);
		qm -> mesh -> xyzs[3*j    ]  = eval[0];
		qm -> mesh -> xyzs[3*j + 1]  = eval[1];
		qm -> mesh -> xyzs[3*j + 2]  = eval[2];
		qm -> mesh -> valence[j][0]  = qm -> backupMesh -> valence[i][0];
		qm -> mesh -> valence[j][1]  = qpiv[qm -> backupMesh -> valence[i][1] - 1] + 1;
		for ( k = 0 ; k < qm -> backupMesh -> valence[i][0]; ++k)
			qm -> mesh -> valence[j][2 + k] = vpiv[qm -> backupMesh -> valence[i][2 + k] -1] + 1;
	}
	EG_free(vpiv);
	EG_free(qpiv);
	return checkMesh (qm);
}


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
	printf(" =================== quad %d ==================\n", qID );
	for ( i = 0 ; i <= 4; i++ ) {
		v = qm -> mesh -> quadIdx [ 4 * ( qID - 1) +  i%4] - 1;
		printf("%lf %lf %lf %lf %lf %d\n",
				qm -> mesh -> xyzs [3 * v ],qm -> mesh -> xyzs [3 * v +1], qm -> mesh -> xyzs [3 * v +2 ],
				qm -> mesh -> uvs [2 * v + 1],  qm -> mesh -> uvs [2 * v ], v + 1);
	}
	printf(" ==============================================\n" );
}


static int EG_normalToSurface (meshMap *qm, double *uv, double *normal ) {
	int stat;
	double eval[18],norm;
	stat = EG_evaluate(qm -> face, uv, eval);
	if ( stat != EGADS_SUCCESS) {
		printf(" EG_normalToSurface EG_evaluate at %lf %lf --> %d!!\n", uv[0], uv[1], stat);
		return stat;
	}
	//
	unitVector (&eval[3], &norm );
	unitVector (&eval[6], &norm );
	if ( qm -> face -> mtype == SREVERSE )
		cross_product(&eval[6], &eval[3], normal);
	else
		cross_product(&eval[3], &eval[6], normal);
	return EGADS_SUCCESS;
}

static void weightedAverage (meshMap *qm, int vID ) {
	int i, auxID, stat;
	double uvc[2], eval[18], *length = NULL, totLength = 0.0;
	if ( qm -> mesh -> vType[ vID -1] != -1 ) return;
	length = (double * ) EG_alloc ( qm -> mesh -> valence [ vID - 1][0] * sizeof (double ) ) ;
	if ( length == NULL ) return;
	printf(" Average coordinates arclength scaled %d \n ", vID ) ;
	for ( i   = 0; i < qm -> mesh -> valence [ vID - 1][0]; i++) {
		auxID     = qm -> mesh -> valence [ vID - 1][2 + i];
		length[i] = EG_segment ( qm, & qm -> mesh -> uvs [ 2 * ( vID - 1 ) ] ,
				&qm -> mesh -> uvs [ 2 * ( auxID - 1 ) ]);
		totLength += length[i];
	}
	uvc[0] = 0.0; uvc[1]  = 0.0;
	for ( i   = 0; i < qm -> mesh -> valence [ vID - 1][0]; i++) {
		auxID = qm -> mesh -> valence [ vID - 1][2 + i] - 1;
		uvc[0] += (length[i] / totLength ) * qm -> mesh -> uvs    [2*auxID    ];
		uvc[1] += (length[i] / totLength ) * qm -> mesh -> uvs    [2*auxID + 1];

	}
	EG_free (length);
	stat    = EG_evaluate ( qm -> face, uvc, eval);
	if ( stat != EGADS_SUCCESS ) return;
	for ( i = 0 ; i < 3; i++ ) {
		if ( i < 2 )
			qm -> mesh -> uvs [ 2 * ( vID -1 ) + i] = uvc[i];
		qm -> mesh   -> xyzs[ 3 * ( vID -1 ) + i] = eval[i];
	}
	return;
}


static void
averageCoordsUsingLinks ( meshMap *qm, int vc, int l1, int l2 ) {
	int i ;
	double uv[2], eval[18] ;
	printf(" AVERAGE COORDS USING LINKS %d %d %d\n ", vc, l1, l2 );
	if ( qm ->mesh -> vType[ vc ] != -1 ) return;
	uv[0]  = qm -> mesh -> uvs[2 * (vc - 1)  ];
	uv[1]  = qm -> mesh -> uvs[2 * (vc - 1) + 1];
	uv[0] += qm -> mesh -> uvs[2 * (l1 - 1)   ];
	uv[1] += qm -> mesh -> uvs[2 * (l1 - 1) + 1];
	uv[0] += qm -> mesh -> uvs[2 * (l2 - 1)   ];
	uv[1] += qm -> mesh -> uvs[2 * (l2 - 1) + 1];
	uv[0] /= 3; uv[1] /= 3;
	i = EG_evaluate (qm -> face, uv, eval);
	vc--;
	for ( i = 0 ; i < 3; i++ ) {
		if ( i < 2 ) qm -> mesh -> uvs[2 * vc + i] = uv[i];
		qm -> mesh -> xyzs [3 * vc + i] = eval[i];
	}
}




static void
averageCoordsMinusLinks ( meshMap *qm, int vc, int l1, int l2 ) {
	int i, j, k, n, stat, *links = NULL ;
	double uv[2], eval[18] ;
	vc--;
	printf(" AVERAGE COORDS MINUS LINKS %d %d %d\n ", vc + 1, l1, l2 );
	if ( qm ->mesh -> vType[ vc ] != -1 ) return;
	uv[0]  = qm -> mesh -> uvs[2 *vc    ];
	uv[1]  = qm -> mesh -> uvs[2* vc + 1];
	n      = 1;
	links  = (int * ) EG_alloc ( qm -> mesh -> valence[vc][0] * sizeof ( int ) );
	printVertexCoords (qm, vc + 1 ) ;
	printVertexCoords (qm , l1 ) ;
	printVertexCoords (qm , l2 ) ;
	if ( links == NULL ) return;
	for ( i = k = 0; i < qm -> mesh -> valence[vc][0]; i++ ) {
		j = qm -> mesh -> valence[vc][2 + i];
		if ( j == l1 || j == l2 ) continue;
		uv[0] += qm -> mesh -> uvs[2 * ( j - 1)    ];
		uv[1] += qm -> mesh -> uvs[2 * ( j - 1) + 1];
		links[k++] = j;
		n++;
	}
	uv[0] /= n; uv[1] /= n;
	stat = EG_evaluate (qm -> face, uv, eval);
	if ( stat != EGADS_SUCCESS ) {
		EG_free (links );
		return;
	}
	for ( i = 0 ; i < 3; i++ ) {
		if ( i < 2 ) qm -> mesh -> uvs[2 * vc + i] = uv[i];
		qm -> mesh -> xyzs [3 * vc + i] = eval[i];
	}
	EG_free(links);
}


static int vertexLinksToBounds (meshMap *qm, int vID ) {
	int i, j, k = 0;
	if ( qm -> mesh -> vType[ vID -1 ] != -1 ) k = 1;
	for (i = 0 ; i < qm -> mesh -> valence[ vID -1 ][0]; i++ ) {
		j  =  qm -> mesh -> valence[ vID -1 ][2 + i] - 1;
		if ( qm -> mesh -> vType[ j ] >= 4 )  return 4;
		else if ( qm -> mesh -> vType[ j ] != -1 ) k = 1;
	}
	return k;
}



static int validInterior(meshMap *qm, int qID, double minAngle, double maxAngle ) {
	int    vj, k, i, i0, sum , stat, area, links[2], orient[8], piv[8], it = 0, itMAX = 20, aux ;
	double angles[8], daux ;
	printf(" \n\n-----------   CheckInvalid quad %d min max angles %lf %lf------------\n", qID, minAngle, maxAngle );
	printMesh (qm, buffer,0);
	stat        = quadAngleOrientation(qm, qID, &area, orient, piv, angles);
	if ( stat != EGADS_SUCCESS ) {
		printf("In checkInvalidElement: quadAlgebraic area at %d is %d !!\n", qID, stat);
		return stat;
	}
	if ( area == 2 && maxAngle < MAXVALIDANGLE) area = 0;
	else if ( area == 1 ) {
		for ( k = 0 ; k < 4; k++ ) {
			aux = qm -> mesh -> quadIdx [ 4 * (  qID -1 ) + piv[k] ] - 1;
			if ( qm -> mesh -> vType[aux] != -1 ||
					( maxAngle < PI * 0.7 && qm -> mesh -> valence[aux][0] != 4 )) continue;
			if (angles [ piv[k] ] > maxAngle || angles[ piv[k] ] < minAngle ) area = 0 ;
		}
	}
	printf (" area %d angle 0 %f angle 3 %f \n ", area, angles[piv[0]], angles[piv[3]]);
	while ( area == 0 && it < itMAX ) {
		for  ( k = 0; k < 4; k++) {
			vj       = qm -> mesh -> quadIdx [ 4 * ( qID - 1) +   piv[k]        ];
			links[0] = qm -> mesh -> quadIdx [ 4 * ( qID - 1) + ( piv[k] + 1)%4 ];
			links[1] = qm -> mesh -> quadIdx [ 4 * ( qID - 1) + ( piv[k] + 3)%4 ];
			if ( qm -> mesh -> vType [ vj -1 ] != -1 ) continue;
			if ( vertexLinksToBounds (qm, vj ) == 4 ) {
				continue;
				if (orient[piv[k]] == 1 ) continue;
				for ( sum = i = 0 ; i < qm -> mesh -> valence[vj][0]; i++ ) {
					aux = qm -> mesh -> valence[vj][ 2 + i ] - 1;
					if ( qm -> mesh -> vType[aux] == 4 ) {
						i0 = aux;
						sum++;
					}
				}
				if ( sum != 1 ) continue;
				daux = EG_segment ( qm, &qm -> mesh -> uvs [ 2 * vj ], & qm -> mesh -> uvs [ 2 * i0 ] );
				if ( daux < qm -> minsize * 0.5 ) continue;
			}
			if (  vertexLinksToBounds (qm, vj ) == 4 && orient[piv[k]] == 1 ) continue;
			if ( orient [piv[k]] == -1  ) {
				printQuadCoords (qm, qID );
				averageCoordsMinusLinks (qm, vj, links[0], links[1]);
				printQuadCoords (qm, qID );
			}
			else if ( angles[ piv[k]] < minAngle && minAngle < PI * 0.25) {
				printQuadCoords (qm, qID );
				averageCoordsUsingLinks (qm, vj, links[0], links[1]);
				printQuadCoords (qm, qID );
			}
			else weightedAverage ( qm, vj ) ;
			stat = quadAngleOrientation(qm, qID, &area, &orient[4], &piv[4], &angles[4]);
			if ( area == 1 ) break;
		}
		printQuadCoords (qm, qID );
		stat = quadAngleOrientation(qm, qID, &area, orient, piv, angles);
		if ( area != 1 ) area = 0;
		else {
			for ( k = 0 ; k < 4; k++ ) {
				aux = qm -> mesh -> quadIdx [ 4 * (  qID -1 ) + piv[k] ] - 1;
				if ( qm -> mesh -> vType[aux] != -1 || qm -> mesh -> valence[aux][0] != 4 ) continue;
				if (angles [ piv[k] ] > maxAngle || angles[ piv[k] ] < minAngle ) area = 0 ;
			}
		}
		it++;
	}
	if ( area == 0 ) return EGADS_GEOMERR;
	return EGADS_SUCCESS;
}



static int validBoundary ( meshMap *qm, double minAngle, double maxAngle ) {
	int v, i, j, k, k1, k2, aux, aux2, q, stat, links[2], i0, nv, dir, area, orient[4], piv[4], it, itMax = 10;
	double range[4], eval[18], angleAtV, thetaRot, cross[4], normal[4], angles[5], tAB[4], tn[4], proj[6], vec[8], segAv,
	angle0, angle1, angle2, angle3, angle4, dt, theta;
	vStar *star = NULL, *star2 = NULL;
	stat = EG_getRange ( qm -> face, range, &i);
	printf(" VALID BOUNDARY CHECK \n ");
	printMesh(qm, buffer,0);
	for ( v = 0 ; v < qm -> mesh -> totVerts; v++ ) {
		if ( qm -> mesh -> vType[v] < 0 ) continue;
		stat      = EG_buildStar (qm ->mesh, &star, v + 1);
		if (stat != EGADS_SUCCESS || star == NULL ) return stat;
		stat      = EG_angleAtBoundaryVertex (qm, v + 1, links, &angleAtV);
		if (stat != EGADS_SUCCESS ) {
			printf("validBoundary: EG_angleAtBoundaryVertex %d --> %d!!\n ", v +1, stat );
			EG_freeStar ( &star );
			return stat;
		}
		stat = EG_normalToSurface (qm, &qm -> mesh -> uvs [ 2 * v ], normal);
		if ( stat != EGADS_SUCCESS ) {
			printf(" validBoundary: EG_normalToSurface at %d --> %d!!\n", v + 1, stat);
			return stat;
		}
		tAB[0] = qm ->mesh -> uvs [2 * v    ];
		tAB[1] = qm ->mesh -> uvs [2 * v + 1];
		tAB[2] = qm ->mesh -> uvs [2 * ( links[0] - 1)    ];
		tAB[3] = qm ->mesh -> uvs [2 * ( links[0] - 1) + 1];
		stat   = EG_projectToTangentPlane(normal,
				&qm -> mesh -> xyzs [ 3 * v ],
				&qm -> mesh -> xyzs [ 3 * (links[0] -1)], proj);
		for ( i   = 0 ; i < 3; i++ )
			vec[i] = proj[i] - qm -> mesh -> xyzs [ 3 * v + i];
		unitVector (vec, &vec[3] );
		printf(" CENTRE %d LINKS %d %d angle %f \n ", v + 1, links[0], links[1], angleAtV );
		for ( nv = j = q = 0; q < star -> nQ; q++ ) {
			if ( star -> quads[q] == -1 ) continue;
			stat = quadAngleOrientation(qm, star -> quads[q], &area, orient, piv, angles);
			if ( stat != EGADS_SUCCESS) {
				printf("Inside validBoundary: quadANgleOrientation  %d = %d!!\n", star -> quads[q], stat);
				EG_free(star);
				return stat;
			}
			if ( area != 1 ) j = 1;
			else {
				for ( k = 0 ; k < 4; k++ ) {
					aux = qm -> mesh -> quadIdx [ 4 * (  star -> quads[q] -1 ) + piv[k] ] - 1;
					if ( qm -> mesh -> vType[aux] != -1 || qm -> mesh -> valence[aux][0] != 4 ) continue;
					if      (qm -> mesh -> vType[v] == 4 && (angles [ piv[k] ] > PI * 0.8 || angles[ piv[k] ] < PI * 0.1 )) j = 1;
					else if (qm -> mesh -> vType[v] != 4 && (angles [ piv[k] ] > maxAngle || angles[ piv[k] ] < minAngle )) j = 1;

				}
			}
			if ( star -> verts [ 2 * q + 1 ] == links[0] ) i0 = q;
			if ( qm -> mesh -> vType [star -> verts [ 2 * q + 1 ] - 1] == -1 ) nv++;
		}
		if ( j == 0 ) continue;
		printf(" V %d in boundary is going to move\n ", v + 1 );
		for ( dir = 0 ; dir < 2; dir++ ) {
			thetaRot  = (1.0) / (double)(nv + 1) * angleAtV ;
			printf(" TOT LEGS %d  ANGLE TEST %f \n ", nv + 1, thetaRot );
			if ( dir  ==  1 ) thetaRot *= -1.0;
			tn[0]      = ( tAB[2] - tAB[0] ) * cos ( thetaRot ) - ( tAB[3] - tAB[1] ) * sin ( thetaRot );
			tn[1]      = ( tAB[2] - tAB[0] ) * sin ( thetaRot ) + ( tAB[3] - tAB[1] ) * cos ( thetaRot );
			tn[2]      = tAB[0] + tn[0];
			tn[3]      = tAB[1] + tn[1];
			printf(" Move from %d rotated %d dir %d\n ", v + 1, links[0], dir );
			stat       = moveDistAway ( qm, qm -> minsize , tAB, &tn[2] );
			if ( stat != EGADS_SUCCESS ) {
				if ( dir == 1 ) {
					printf("we have gone both directions  and still out of domain !!!! \n ");
					EG_freeStar ( &star);
					return EGADS_GEOMERR;
				}
				dir = ( dir + 1 ) %2 ;
				break;
			}
			stat       = EG_evaluate ( qm -> face, &tn[2], eval);
			printf(" Evaluate %lf %lf --> %lf %lf %lf \n ", tn[2], tn[3], eval[0], eval[1], eval[2]);
			printf(" uA %f %f --- uB %f  %f  \n ", tAB[0], tAB[1], tAB[2], tAB[3] );
			if ( stat != EGADS_SUCCESS ) {
				printf(" We can't evaluate at EPS distance from vertex !\n ");
				if ( dir == 1 ) {
					printf("we have gone both directions  and still out of domain !!!! \n ");
					EG_freeStar ( &star);
					return EGADS_GEOMERR;
				}
				continue;
			}
			stat = EG_projectToTangentPlane(normal, &qm -> mesh -> xyzs [ 3 * v], eval, &proj[3]);
			if ( stat != EGADS_SUCCESS ) {
				printf("EG_checkInvalidElement: EG_projectToTangentPlane %d, %d --> %d!!\n", v, links[0], stat);
				return stat;
			}
			for ( i = 0 ; i < 3; i++ )
				vec[4 + i] = proj[3 + i] - qm -> mesh -> xyzs [ 3 * v + i];
			unitVector    (&vec[4], &vec[7] );
			cross_product (vec, &vec[4], cross);
			if ( dotProduct (normal, cross ) > 0 ) break;
			else {
				if ( dir == 1 ) {
					printf(" both directions give negative dot product !!\n ");
					EG_freeStar ( &star);
					return EGADS_GEOMERR;
				}
				dir = ( dir + 1 ) %2 ;
				break;
			}
		}
		printf("\n\n d = %d -- NOW COMPUTE ALL COORDINATES \n ", dir ) ;
		angle0 = 0.0;
		for (j = 0 ; j < nv ; j++ ) {
			aux        = star -> verts [ star -> idxV [ 2 * (i0 + j + 1) + 1 ] ] - 1;
			printf(" j %d VERTEX %d type %d \n ", j, aux, qm -> mesh -> vType [ aux ] );
			thetaRot = (double ) ( j + 1) / (double)(nv + 1) * angleAtV ;
			if ( dir ==  1 ) thetaRot *= -1.0;
			k1      = star -> verts [ star -> idxV [ 2 * (i0 + j    ) + 1 ] ] - 1;
			k2      = star -> verts [ star -> idxV [ 2 * (i0 + j + 2) + 1 ] ] - 1;
			angle1  = EG_angleAtVnormalPlane ( qm, star -> verts[0], links[0], k1 + 1);
			if ( dir == 1 ) angle1 *= -1.0;
			for (it = 0 ; it < itMax; it++ ) {
				dt      = angle1 + (double)(it + 1)/ (double)itMax * thetaRot;
				printf(" IT %d rotate %f \n ", it, dt );
				angle2  = EG_angleAtVnormalPlane ( qm, star -> verts[0], k1 + 1  , aux + 1);
				angle3  = EG_angleAtVnormalPlane ( qm, star -> verts[0], links[0], aux + 1);
				angle4  = EG_angleAtVnormalPlane ( qm, star -> verts[0], aux + 1 ,  k2 + 1);
				segAv   = EG_segment ( qm, &qm -> mesh -> uvs [ 2 * v ], & qm -> mesh -> uvs [ 2 * aux ]);
				if ( angle2 < PI && angle3 > angle0 && angle4 < PI && segAv >= qm -> minsize * 0.5) {
					printf(" ROTATAION WORKED. break loop \n");
					break;
				}
				printf(" Vertex has angle wrp prev %f and wrp ori %f need to move \n ", angle2, angle3);
				tn[0] = ( tAB[2] - tAB[0] ) * cos ( dt ) - ( tAB[3] - tAB[1] ) * sin ( dt );
				tn[1] = ( tAB[2] - tAB[0] ) * sin ( dt ) + ( tAB[3] - tAB[1] ) * cos ( dt );
				printf("ROTATE %f ORI %d %d move vertex %d\n ", dt, star ->verts[0], links[0], aux + 1);
				tn[2] =   tAB[0] + tn[0];
				tn[3] =   tAB[1] + tn[1];
				if ( angle2 < PI && angle3 > angle0 && angle4 < PI && segAv < qm -> minsize * 0.5) {
					printf(" just the sizing is wrong. push point\n");
					tn[2] =	qm -> mesh -> uvs[ 2 * aux    ];
					tn[3] = qm -> mesh -> uvs[ 2 * aux + 1];
				}
				stat   = moveDistAway ( qm, qm -> minsize * 0.5, tAB, &tn[2] );
				stat   = EG_evaluate ( qm -> face, &tn[2], eval);
				printf(" Evaluate %lf %lf --> %lf %lf %lf \n ", tn[2], tn[3], eval[0], eval[1], eval[2]);
				if ( stat != EGADS_SUCCESS ) {
					printf(" We have messed up with recomputing  !!!! \n ");
					EG_freeStar ( &star);
					return EGADS_GEOMERR;
				}
				for ( i = 0; i < 3; i++ ) {
					if ( i < 2 ) qm -> mesh -> uvs[ 2 * aux + i] = tn[2 + i];
					qm -> mesh -> xyzs[ 3 * aux + i]             = eval[i];
				}
				stat = EG_buildStar ( qm -> mesh, &star2, aux + 1 );
				if (stat != EGADS_SUCCESS || star2 == NULL ) {
					EG_freeStar ( &star2 );
					return stat;
				}
				printf(" BUILD STAR AROUND %d \n ", aux + 1 );
				for ( i    = 0; i < star2 -> nQ; i++ ) {
					stat   = quadAngleOrientation(qm, star2 -> quads[i], &area, orient, piv, angles);
					angle4 = EG_angleAtVnormalPlane ( qm, aux + 1, star2 -> verts [ 2 * i + 1], star2 -> verts [ star2 -> idxV[ 2 * i + 3] ] );
					if ( area != 1 || angle4 > PI )
						stat = validInterior (qm, star2 -> quads[i], 0.0, PI  );
				}
				EG_freeStar (&star2);
			}
			angle0 = EG_angleAtVnormalPlane ( qm, star -> verts[0], links[0], aux + 1);
			printMesh (qm, buffer,0);
		}
		EG_freeStar (&star);
		printf(" ------------- NEXT VERTEX -----------\n");
	}
	printMesh(qm, buffer,0);
	return EGADS_SUCCESS;
}







static int
EG_projectToTangentPlane(double normal[], double *O, double *p, double *proj) {
	double c, dotNN = 0.0, dotNP = 0.0, dist, lambda;
	dist    = (p[0]- O[0]) * (p[0]- O[0]) + (p[1]- O[1])*(p[1]- O[1]) + (p[2]- O[2])*(p[2]- O[2]);
	dist    = sqrt(dist);
	if (dist < EPS11 ) {
		proj[0] = p[0]; proj[1] = p[1]; proj[2] = p[2];
		return EGADS_SUCCESS;
	}
	c       = normal[0] *      O[0] + normal[1] *      O[1] + normal[2] *      O[2]; // Equation plane: a*x + b*y + c*z = C
	dotNN   = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
	dotNP   = normal[0] *      p[0] + normal[1] *      p[1] + normal[2] *      p[2];
	if ( fabs(dotNP - c) <= EPS11) {
		proj[0] = p[0]; proj[1] = p[1]; proj[2] = p[2];
		return EGADS_SUCCESS;
	}
	lambda  =   (c - dotNP) / dotNN;
	proj[0] = p[0] + lambda * normal[0];
	proj[1] = p[1] + lambda * normal[1];
	proj[2] = p[2] + lambda * normal[2];
	// check that point belongs to plane.
	dist  = normal[0] * proj[0] + normal[1] * proj[1] +  normal[2] * proj[2];
	if( fabs(dist - c) < EPS11) {
		return EGADS_SUCCESS;
	}
	else{
		printf(" ORIGIN %lf %lf  %lf  NORMAL %lf %lf  %lf  TARGET %lf %lf %lf\n", O[0], O[1], O[2], normal[0], normal[1], normal[2], p[0], p[1], p[2]);
		printf(" POINT SHOULD BELONG TO PLANE!!!!! %.16e ~= %.16e\n",dist,c);
		printf(" DOT NN %lf PN %lf LAMBDA %lf  POINT %lf %lf %lf\n", dotNN, dotNP, lambda, proj[0], proj[1], proj[2]);
		return EGADS_GEOMERR;
	}
}


static double
EG_angleAtVnormalPlane ( meshMap *qm, int vC, int v1, int v2 ) {
	int i, stat;
	double dot1, dot2, evalVert[18], normal[4], cross[3], proj0[3], proj1[3], proj2[3], v01[4], v02[4];
	stat   = EG_evaluate(qm -> face, &qm -> mesh -> uvs[ 2 * ( vC - 1 ) ], evalVert);
	if ( stat != EGADS_SUCCESS) return stat;
	if ( qm -> face -> mtype == SREVERSE )
		cross_product(&evalVert[6], &evalVert[3], normal);
	else
		cross_product(&evalVert[3], &evalVert[6], normal);
	unitVector(normal, &normal[3]);
	stat = EG_projectToTangentPlane(normal, evalVert, &qm -> mesh -> xyzs [ 3 * ( vC - 1)], proj0);
	stat = EG_projectToTangentPlane(normal, evalVert, &qm -> mesh -> xyzs [ 3 * ( v1 - 1)], proj1);
	stat = EG_projectToTangentPlane(normal, evalVert, &qm -> mesh -> xyzs [ 3 * ( v2 - 1)], proj2);
	//printVertexCoords ( qm, vC );
	//printVertexCoords ( qm, v1 );
	//printVertexCoords ( qm, v2 );
	for ( i = 0 ; i < 3 ; i++ ) {
		v01[i] = proj1[i] - proj0[i];
		v02[i] = proj2[i] - proj0[i];
	}
	unitVector (v01, &v01[3] );
	unitVector (v02, &v02[3] );
	cross_product (v01, v02, cross );
	dot1 = dotProduct ( v01, v02 );
	dot2 = dotProduct ( normal, cross);
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
		printf("EG_angleAtBoundaryVertex:: vertex %d is at surface bounds and connected only to another boundary vertex !!\n ", v);
		return EGADS_GEOMERR;
	}
	*size = EG_angleAtVnormalPlane (qm, v , links[0], links[1] );
	return EGADS_SUCCESS;
}

static int moveDistAway ( meshMap *qm, double size, double *uv0, double *uv1 ) {
	int    i, stat;
	double dt, eval[18], eval2[18], range[4], vec[2], seg, uvPrev[2], norm;
	vec[0]     = uv1[0] - uv0[0];
	vec[1]     = uv1[1] - uv0[1];
	stat       = EG_getRange ( qm -> face, range, &i);
	if ( stat != EGADS_SUCCESS ) {
		printf("In moveEpsAway: getRange at face %d -> %d!!\n ", qm -> fID, stat );
		return stat;
	}
	if ( size  < qm -> minsize * 0.5) size = 0.5 * qm -> minsize;
	if (uv1[0] < range[0] || uv1[0] > range[1] ||
			uv1[1] < range[2] || uv1[1] > range[3] )
	{
		uv1[0] = uv0[0];
		uv1[1] = uv0[1];
		return EGADS_EMPTY;
	}
	stat      = EG_evaluate ( qm -> face, uv0, eval );
	stat      = EG_evaluate ( qm -> face, uv1, eval2);
	//printf(" GO FROM\n %f %f %f  \n %f %f %f\n ", eval[0], eval[1], eval[2], eval2[0], eval2[1], eval2[2]);
	seg       = EG_segment ( qm, uv0, uv1 );
	norm      = sqrt (vec[0] * vec[0] + vec[1] * vec[1] );
	dt        = 0.1 * norm;
	while ( seg > size ) {
		dt     *= 0.5;
		uv1[0] = uv0[0] + vec[0] * dt;
		uv1[1] = uv0[1] + vec[1] * dt;
		seg    = EG_segment ( qm, uv0, uv1 );
	}
	if ( fabs ( seg - size ) < 0.05 * size ) return EGADS_SUCCESS;
	uv1[0]     = uv0[0];
	uv1[1]      = uv0[1];
	uvPrev[0]   = uv1[0];
	uvPrev[1]   = uv1[1];
	while ( size > seg ) {
		uv1[0] += vec[0] * dt;
		uv1[1] += vec[1] * dt;
		norm    = vec[0] * ( uv1[0] - uvPrev[0]) + vec[1] * ( uv1[1] - uvPrev[1]);
		stat    = EG_evaluate ( qm -> face, uv1, eval);
		//    printf("%lf %lf %lf %lf %lf IN MOVE DIST \n ", eval[0], eval[1], eval[2], uv1[0], uv1[1] );
		seg    = EG_segment ( qm, uv0, uv1 );
		//  printf("#  %f %f < %f\n ", dt, seg, size);
		if (norm < 0 || uv1[0] < range[0] || uv1[0] > range[1] ||
				uv1[1] < range[2] || uv1[1] > range[3]) {
			//  printf(" OUT OF RANGE || norm < 0 %f \n", norm);
			uv1[0] = uvPrev[0];
			uv1[1] = uvPrev[1];
			return EGADS_SUCCESS;
		}
		uvPrev[0] = uv1[0];
		uvPrev[1] = uv1[1];
	}
	//printf(" OUB BECAUSE size < seg ::: %f   %f\n ", size, seg );
	return EGADS_SUCCESS;
}



// Will detect if quad intersects but allows obtuse angles during regularization
static int
quadAngleOrientation(meshMap *qm, int qID, int *validArea, int *ori, int *order, double *theta ) {
	int     i, qV[4], k, k1, k2, sign, stat, vA, vB, vC;
	double cross[4],qNormal[4], vABCD[12], dot, quv[5], projABCD[12];
	*validArea = 0;
	stat = checkQuad ( qm -> mesh, qID );
	if ( stat != EGADS_SUCCESS ) {
		printf(" In quadAngleOrientation:: checking area of quad %d  !!! %d \n", qID, stat );
		return stat;
	}
	stat  = quadAverageCoords (qm, qID, quv, &quv[2]);
	if ( stat != EGADS_SUCCESS) {
		printf("vectorAtVertexNplane :: Average coordinates for quad %d --> %d!!\n ", qID, stat);
		return stat;
	}
	stat = EG_normalToSurface ( qm, quv, qNormal);
	if ( stat != EGADS_SUCCESS) {
		printf("quadAngleOrientation at %d: EG_normalToSurface %d !!\n ", qID, stat);
		return stat;
	}
	for ( i = 0 ; i < 4; i++ ) {
		qV[i]      = qm -> mesh ->quadIdx[4*(qID - 1) + i] - 1;
		stat       = EG_projectToTangentPlane(qNormal, &quv[2], &qm -> mesh -> xyzs [ 3 * qV[i] ], &projABCD[3 * i]);
		if ( stat != EGADS_SUCCESS ) {
			printf("quadAngleOrientation :: EG_projectToTangentPlane quad %d vert %d --> %d !!\n", qID, qV[i], stat );
			return stat;
		}
	}
	for ( k = 0 ; k < 4; ++k) {
		ori[k] = 1;
		vA     = k;
		for ( k1 = 1; k1 < 3; k1++ ) {
			vB   = ( k + k1     )%4;
			for ( k2 = k1 + 1; k2 < 4; k2++ ) {
				vC = ( k + k2 )%4;
				for ( i = 0 ; i < 3; ++i) {
					vABCD[    i] = projABCD[3 * vB + i] - projABCD[3 * vA + i];
					vABCD[4 + i] = projABCD[3 * vC + i] - projABCD[3 * vA + i] ;
				}
				unitVector     (&vABCD[0], &vABCD[3]   );
				unitVector     (&vABCD[4], &vABCD[7]   );
				cross_product  (vABCD, &vABCD[4], cross);
				if (dotProduct (qNormal, cross ) < 0 ) ori[k] = -1;
			}
		}
		if      ( ori[k] ==  1 && *validArea == 0 ) *validArea = 1;
		else if ( ori[k] == -1 && *validArea == 1 ) *validArea = 2;

	}
	for ( sign = k = 0 ; k < 4; ++k) {
		vA   =   k        ;
		vB   = ( k + 1 )%4;
		vC   = ( k + 3 )%4;
		for ( i = 0 ; i < 3; ++i) {
			vABCD[    i] = projABCD[3 * vB + i] - projABCD[3 * vA + i];
			vABCD[4 + i] = projABCD[3 * vC + i] - projABCD[3 * vA + i] ;
		}
		unitVector    (&vABCD[0], &vABCD[3] );
		unitVector    (&vABCD[4], &vABCD[7] );
		cross_product (vABCD, &vABCD[4], cross);
		dot         = dotProduct (vABCD, &vABCD[4]);
		if ( dotProduct (qNormal, cross ) < 0 ) ori[k] = -1;
		else                                    ori[k] =  1;
		if      ( fabs(dot - 1.0) < EPS11 ) theta[k] = 0.0;
		else if ( fabs(dot + 1.0) < EPS11 ) theta[k] = PI;
		else                                theta[k] = acos(dot);
		if (EG_quadIsBoundary ( qm -> mesh, qID) != 0 &&
				ori[k] == -1 && ( 2.0 * PI - theta[k] ) > PI * 1.1 ) *validArea = 0;
		order[k] = k;
		if ( ori[k] == -1 ) sign = 1;
	}
	if ( sign == 0 ) *validArea = 1;
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
	printf("------------ Internal angles at quad %d ordered by size --------------\n", qID);
	for ( sign = k = 0 ; k < 4; ++k)
		printf("Vertex %d has angle %f and orientation %d \n ", qV[order[k]] + 1, theta[order[k]], ori[order[k]]);
	return EGADS_SUCCESS;;
}

static int
EG_buildStar(meshData *mesh, vStar **star, int vID ) {
	int  stat, i = 0 , id0 = -1, q = 0, auxV, auxQ,  v = 0, quadID, prevQuad, it = 0, it2 = 0;
	int adj[2], *vertex = NULL, *quads = NULL;
	int qLoop[8] = {0, 1, 2, 3, 0, 1, 2, 3};
	stat       = checkVertex ( mesh , vID ) ;
	if ( stat != EGADS_SUCCESS ) return stat;
	vertex     = (int * ) EG_alloc ( mesh -> totVerts * sizeof ( int ) );
	quads      = (int * ) EG_alloc ( mesh -> totQuads * sizeof ( int ) );
	if ( vertex == NULL || quads == NULL ) {
		printf("EG_buildStar MALLOC at quads & verts!!\n ");
		return EGADS_MALLOC;
	}
	// quads are -1 bias
	quadID  = mesh -> valence[ vID - 1 ][1] - 1;
	i       = checkQuad ( mesh, quadID + 1);
	if ( i != EGADS_SUCCESS){
		printf(" EG_buildStar at vertex %d has associated a bad quad: %d --> %d!!\n", vID, quadID + 1, i );
		printQuadSpecs ( mesh, quadID + 1);
		EG_free ( vertex );
		EG_free ( quads  );
		return i;
	}
	vertex [v++] = vID;
	it           = 0;
	do {
		stat = checkQuad ( mesh, quadID + 1);
		if ( stat != EGADS_SUCCESS ) {
			printf(" In EG_buildStar quad %d is bad quad --> %d!!\n",  quadID + 1, stat);
			printQuadSpecs( mesh, quadID + 1);
			EG_free (quads );
			EG_free (vertex);
			return i;
		}
		id0 = EG_quadVertIdx ( mesh, quadID + 1, vID );
		if ( id0 < 0 ) {
			printf(" In EG_buildStar id for Vert %d in Quad %d is %d !!\n", vID, quadID + 1, id0);
			printQuadSpecs( mesh, quadID + 1);
			EG_free ( vertex ) ;
			EG_free ( quads ) ;
			return EGADS_INDEXERR;
		}
		for ( i = 1 ; i <= 2; ++i) vertex[v++] = mesh -> quadIdx[ 4 * quadID   + qLoop[id0 + i ] ] ;
		quads[q++] = quadID + 1;
		prevQuad   = quadID;
		quadID     = mesh -> quadAdj[ 4 * prevQuad + qLoop[id0 + 3 ] ] - 1;
		if ( quadID  < 0 ) { //make a "ghost" quad
			auxQ        = prevQuad;
			vertex[v++] = mesh -> quadIdx[4*auxQ + qLoop[id0 + 3]] ;
			auxV        = mesh -> quadIdx[4*auxQ + qLoop[id0 + 1]] ;
			it2 = 0;
			do
			{
				stat = EG_adjQtoPair(mesh, auxQ + 1, vID, auxV, adj );
				if ( adj[1] == -1 ) break;
				auxQ = adj[1] - 1;
				i    = EG_quadVertIdx ( mesh, auxQ + 1, vID ) ;
				if ( i < 0 ){
					printf(" In buildStar vertex Id %d in quad %d is %d\n", vID, quadID + 1, i);
					printQuadSpecs ( mesh, quadID + 1);
					EG_free ( vertex ) ;
					EG_free ( quads ) ;
					return EGADS_INDEXERR;
				}
				auxV = mesh -> quadIdx [ 4 * auxQ + qLoop [ i + 1 ] ];
				it2++;
				if ( it2 > 200 ) {
					printf(" stuck in interior loop of build star!!!!!!!!\n");
					EG_free (quads);
					EG_free (vertex);
					return EGADS_RANGERR;
				}
			}
			while ( adj[1] != - 1 );
			quads  [q++] = -1 ;
			vertex [v++] = -1 ;
			quadID       = auxQ;
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
EG_removeDoublet(meshMap *qm, int vID) {
	int  qC, stat, link[2], i;
	if ( qm -> mesh -> vType[ vID - 1] != -1 || qm -> mesh -> valence[vID - 1][0]  != 2 ) return EGADS_SUCCESS;
	qC      = qm -> mesh -> valence[vID - 1][1];
	link[0] = qm -> mesh -> valence[vID - 1][2];
	link[1] = qm -> mesh -> valence[vID - 1][3];
	stat    = EG_mergeVertices(qm, qC, vID);
	if ( stat != EGADS_SUCCESS ) {
		printf("EG_removeDoublet I failed at removing doublet %d !! \n", vID);
		return EGADS_GEOMERR;
	}
	// Check that it has left valid elements at each old link
	for ( i = 0 ; i < 2; i++ ) {
		if ( checkVertex ( qm -> mesh, link[i] ) == EGADS_SUCCESS ) {
			stat  = EG_removeDoublet (qm, link[i] ) ;
			if ( stat != EGADS_SUCCESS ) {
				printf("EG_removeDoublet I failed at removing doublet inside doublet %d !! \n", link[i]);
				printMesh(qm, buffer,0);
				return stat;
			}
		}
	}
	return EGADS_SUCCESS;
}

static int
EG_forceCollapse ( meshMap *qm, int qID, int *activity ) {
	int j = 0, i3 = 0, centre, links[2], val[2];
	*activity     = 0;
	if ( EG_nValenceCount ( qm -> mesh, qID, 3 ) == 0 ) return EGADS_SUCCESS;
	while ( j < 2 ) {
		centre = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + i3 ];
		if ( qm -> mesh -> vType [centre -1 ] == -1 && qm -> mesh -> valence[centre -1][0] == 3 ) {
			links[0] = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + (i3 + 1 ) % 4 ];
			links[1] = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + (i3 + 3 ) % 4 ];
			val  [0] = getValence ( qm -> mesh, links[0] );
			val  [1] = getValence ( qm -> mesh, links[1] );
			if ( (val[0] >= 5 || val[1] >= 5 || j == 1) && validCollapse (qm -> mesh, qID, centre) == 1) {
				*activity = 1;
				printf(" FORCE collapse mesh\n");
				printMesh (qm, buffer,0);
				j = EG_mergeVertices ( qm, qID, centre );
				printMesh (qm, buffer,0);
				return j;
			}
		}
		i3++;
		if ( i3 == 4 ) {
			j++;
			i3 = 0;
		}
	}
	return EGADS_SUCCESS;
}

static int
EG_collapse (meshMap *qm, int qID, int *activity )  {
	int i, vC, stat, v5[5], v3[5] ;
	*activity     = 0 ;
	stat          = checkQuad ( qm -> mesh , qID );
	if ( stat    != EGADS_SUCCESS ) {
		printf(" EG_collapse for quad %d is %d \n ", qID, stat );
		return stat;
	}
	for ( v3[0] = v5[0] = i = 0 ; i < 4; i++ ) {
		if      ( getValence (qm ->mesh, qm -> mesh -> quadIdx [ 4 *( qID - 1) + i]) == 3) v3[++v3[0]] = i;
		else if ( getValence (qm ->mesh, qm -> mesh -> quadIdx [ 4 *( qID - 1) + i]) == 5) v5[++v5[0]] = i;
	}
	if ( v3[0] + v5[0] <= 2 || v3[0] * v5[0] == 0 ) return EGADS_SUCCESS;
	if ( v5[0] == 1 || ( v5[0] == 2 && abs ( v5[1] - v5[2] ) % 2 == 0 ) ) {
		vC = qm -> mesh -> quadIdx [ 4 * ( qID -1 ) + (v5[1] + 1 ) %4 ];
		if ( validCollapse ( qm -> mesh, qID, vC ) == 0 ) return EGADS_SUCCESS;
		stat      = EG_mergeVertices(qm, qID, vC);
		*activity = 1;
		return  stat;
	}
	return EGADS_SUCCESS;
}


static int
EG_mergeVertices (meshMap *qm, int qC, int centre ) {
	int  stat, i, j, q,  adjq, adjPair[2], auxQ, oldQ[8], links[2];
	int piv[4]  = {1, 0, 3, 2} ;
	vStar *star = NULL;
	// assuming ALWAYS that we collapse vertex 0 to 2
	i = EG_quadVertIdx  ( qm -> mesh, qC, centre );
	if ( qm -> mesh -> vType [ centre - 1] != -1 ) centre = qm -> mesh -> quadIdx [ 4 * ( qC - 1) + ( i + 2 ) %4];
	i = EG_quadVertIdx  ( qm -> mesh, qC, centre );
	if ( qm -> mesh -> vType [ centre - 1] != -1 || i < 0) {
		printf(" EG_mergeVertices collapse through %d with index at %d = %d !!\n ", centre, qC, i );
		printQuadSpecs ( qm -> mesh, qC );
		return EGADS_INDEXERR;
	}
	stat       = EG_buildStar ( qm -> mesh, &star, centre );
	if ( stat != EGADS_SUCCESS || star == NULL ) return stat;
	links[0]   = qm -> mesh -> quadIdx [ 4 * ( qC - 1) + ( i + 1 ) %4];
	links[1]   = qm -> mesh -> quadIdx [ 4 * ( qC - 1) + ( i + 3 ) %4];
	// save old quads vertices and adjacents to update map correctly
	for ( q = 0; q < 4; ++q ) {
		oldQ[q    ] = qm -> mesh -> quadIdx[4 * (qC - 1) + (q + i)%4 ] - 1;
		oldQ[q + 4] = qm -> mesh -> quadAdj[4 * (qC - 1) + (q + i)%4 ] - 1;
	}
#ifdef DEBUGG
	printQuadSpecs ( qm -> mesh, qC );
	printf(" In merge vertices: QUAD %d collapse vertex  centre %d at %d  v(%d) = %d (type %d ) to v(%d) = %d ( type %d) "
			"and remove one valence from vertices v(%d) = %d (type %d ) v(%d) = %d (type %d )\n",qC, centre, i,
			oldQ[0] + 1, qm -> mesh -> valence[ oldQ[0 ]][0], qm -> mesh -> vType[ oldQ[0] ] ,
			oldQ[2] + 1, qm -> mesh -> valence[ oldQ[2 ]][0], qm -> mesh -> vType[ oldQ[2] ] ,
			oldQ[1] + 1, qm -> mesh -> valence[ oldQ[1 ]][0], qm -> mesh -> vType[ oldQ[1] ] ,
			oldQ[3] + 1, qm -> mesh -> valence[ oldQ[3 ]][0], qm -> mesh -> vType[ oldQ[3] ] );
#endif
	j = qm -> mesh -> quadIdx [ 4 * (qC - 1 ) + (i + 2 )%4] - 1;
	if ( qm -> mesh -> vType [j] == -1 ) {
		stat = quadAverageCoords(qm, qC, &qm -> mesh -> uvs[ 2 * oldQ[2] ], &qm -> mesh -> xyzs[ 3 * oldQ[2] ] );
		if ( stat != EGADS_SUCCESS ) {
			EG_freeStar(&star);
			return stat;
		}
	}
	for ( i = 0 ; i < 4; i++ ) {
		q   = oldQ [ 4 + i ]; // - 1
		if ( q < 0 ) continue;
		adjq = oldQ [ 4 + piv[ i ] ]; // - 1 bias
		stat = EG_adjQtoPair ( qm -> mesh, q + 1, oldQ[ i ] + 1, oldQ[ ( i + 1 ) % 4  ] + 1, adjPair );
		if ( stat != EGADS_SUCCESS || adjPair[1] != qC ) {
			EG_freeStar (&star);
			return  EGADS_INDEXERR;
		}
		qm -> mesh -> quadAdj [ 4 * q + adjPair[0] ] = adjq + 1;
	}
	// Eliminate vertex p[0] from all the quads and its valences
	for ( i = 0 ; i < star -> nQ; ++i) {
		q = star -> quads[i] - 1;
		if ( q <  0 ) continue; // ghost quad
		if ( q + 1 == qC ) {  // eliminate quad qC
			for ( j = 0 ; j < 4; ++ j) {
				qm -> mesh -> quadIdx[ 4 * (qC - 1) + j ] = -2;
				qm -> mesh -> quadAdj[ 4 * (qC - 1) + j ] = -2;
			}
		}
		else {
			for ( j = 0 ; j < 4; ++ j)
				if ( qm -> mesh -> quadIdx[ 4 * q  + j] == oldQ[0] + 1) qm -> mesh -> quadIdx[4 * q + j] = oldQ[2] +  1;
		}
	}
	// Point all the collapsed quad vertices to a valid quad
	for ( i = 1 ; i < 4; i++) {
		if ( qm -> mesh -> valence[ oldQ[i] ][1] == qC ) {
			for ( q = 0 ; q < 4; q++ ) {
				auxQ = oldQ[ 4 + q ];
				if ( auxQ < 0 ) continue;
				if ( EG_quadVertIdx ( qm -> mesh, auxQ + 1, oldQ[i] + 1 ) >= 0 ) {
					qm -> mesh -> valence[ oldQ[i] ][1] = auxQ + 1;
					break;
				}
			}
			if ( qm -> mesh -> valence[ oldQ[i] ][1] == qC ) {
				stat = EGADS_INDEXERR;
				EG_freeStar (&star);
				return stat ;
			}
		}
	}
	// Set valences for merged stuff
	for ( i = 1 ; i < 4; i++) {
		stat = setValence ( qm -> mesh, oldQ[i] + 1 );
		if ( stat != EGADS_SUCCESS) {
			printf( "Inside EG_mergeVertices stat in setValence %d = %d\n", i, stat );
			EG_freeStar(&star);
			return stat;
		}
	}
	// set valences to links
	for ( i = 0 ; i < star -> nQ; ++i) {
		j = star -> verts [ 2 * i + 1 ];
		if ( j < 0 || qm -> mesh -> vType [ j - 1] == -2 ) continue;
		stat = setValence ( qm -> mesh, j );
		if ( stat != EGADS_SUCCESS) {
			printf( "Inside EG_mergeVertices stat in setValence %d = %d\n", i, stat );
			EG_freeStar(&star);
			return stat;
		}
	}
	EG_freeStar (&star);
	// delete vertex vC2
	qm -> mesh -> valence [  oldQ[0]    ][0] = -1;
	qm -> mesh -> vType   [  oldQ[0]    ]    = -2; // -2 = removed
	for ( i = 0 ; i < 3; i++ ) {
		if ( i < 2 ) qm -> mesh-> uvs [2*oldQ[0] + i] = qm -> mesh->uvs [2*oldQ[2] + i ];
		qm -> mesh-> xyzs[3*oldQ[0] + i]              = qm -> mesh->xyzs[3*oldQ[2] + i ];
	}
	qm -> mesh -> remQuads[++qm -> mesh -> remQuads[0]] = qC;
	qm -> mesh -> remVerts[++qm -> mesh -> remVerts[0]] = oldQ[0] + 1;
	// remove possible doublets from collapsing: look at links valences
	for  ( i = 0 ; i < 2; i++ ) {
		stat = EG_removeDoublet(qm, links[i]);
		if ( stat != EGADS_SUCCESS) {
			printf("In EG_mergeVertices I failed at removing doublet from %d -> %d !!\n ", links[i], stat);
			printMesh(qm, buffer, 0);
			return stat;
		}
	}
	stat = optimize_angles ( qm, 2, links, 0);
	if ( stat != EGADS_SUCCESS ) {
		printf(" In merge vertices:: failed during optimization %d \n ", stat);
		return stat;
	}
	for ( i = j = 0 ; i < qm -> mesh -> totQuads; i++ ) {
		if ( checkQuad (qm -> mesh, i + 1) == EGADS_SUCCESS ) j++;
	}
	return EGADS_SUCCESS;
}


static int
EG_swap (meshMap *qm, int qIn, int *activity) {
	int   stat, q, swap = 0;
	quadGroup qg;
	*activity  = 0;
	stat       = checkQuad (qm -> mesh, qIn);
	if ( stat != EGADS_SUCCESS ) {
		printf(" EG_swap for quad %d is %d \n ", qIn, stat );
		return stat;
	}
	for ( swap = q = 0 ; q < 4; q++ ) {
		qg.q[0]    = qIn;
		qg.q[1]    = qm -> mesh -> quadAdj [ 4 * ( qIn - 1 ) + q ];
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
	if (swap == 0 ) return EGADS_SUCCESS;
	stat      = EG_swappingOperation (qm, qg, swap);
	*activity = 1;
	return EGADS_SUCCESS;
}


static int
EG_swappingOperation(meshMap *qm, quadGroup qg, int swap ) {
	int stat, i0, i1, vc, i, j, q, adj, idx[4], qID[2], adjQmap[6], area, ori[4], order[4];
	double  angles[4];
	swap   = swap %3;
	qID[0] = qg.q[0]; qID[1] = qg.q[1];
	if ( swap == 0 ) {
		printf(" swapping throu 0-3 will result in the same pair!! \n ");
		return EGADS_INDEXERR;
	}
#ifdef DEBUGG
	printQuadGroup (qm ->mesh, qg );
	printf(" Inside EG_swappingOperation: quads %d  %d swapping %d-%d to %d-%d\n",
			qID[0], qID[1], qg.verts[0], qg.verts[3], qg.verts[swap], qg.verts [ swap + 3]);
#endif
	i0 = EG_quadVertIdx ( qm -> mesh, qID[0], qg.verts[0]); // centre
	i1 = EG_quadVertIdx ( qm -> mesh, qID[1], qg.verts[3]); // opposite
	// Modify Quads and get adj map
	qm -> mesh -> quadIdx [ 4 * ( qID[0] - 1 )    ] = qg.verts[swap ];
	qm -> mesh -> quadIdx [ 4 * ( qID[1] - 1 )    ] = qg.verts[swap ];
	for ( i = 0 ; i < 3; i++ ) {
		adjQmap[i    ] = qm -> mesh -> quadAdj [ 4 * ( qID[0] - 1 ) + ( i + i0 )%4 ];
		adjQmap[i + 3] = qm -> mesh -> quadAdj [ 4 * ( qID[1] - 1 ) + ( i + i1 )%4 ];
		qm -> mesh -> quadIdx [ 4 * ( qID[0] - 1 ) + i + 1 ] = qg.verts[( swap + 3 + i )%6 ];
		qm -> mesh -> quadIdx [ 4 * ( qID[1] - 1 ) + i + 1 ] = qg.verts[( swap + 1 + i )   ];
	}
	//for ( i = 0 ; i < 6; i++ ) printf(" ADJ %d = %d\n", i, adjQmap[i]);
	qm -> mesh -> quadAdj [ 4 * ( qID[0] - 1 )    ] = qID[1];
	qm -> mesh -> quadAdj [ 4 * ( qID[1] - 1 ) + 3] = qID[0];
	for ( i = 0 ; i < 3; i++ ) {
		adj = adjQmap[ (3 + i + swap) % 6 ] - 1;
		qm -> mesh -> quadAdj [ 4 * ( qID[0] - 1 ) + i + 1] = adj + 1;
		if ( adj >= 0 ) {
			for ( j = 0; j < 4; j++ ) if ( qm -> mesh -> quadAdj [ 4 * adj + j ] == qID[1] )
				qm -> mesh -> quadAdj [ 4 * adj + j ] = qID[0];
		}
		adj = adjQmap[i + swap] -1;
		qm -> mesh -> quadAdj [ 4 * ( qID[1] - 1 ) + i    ] = adj + 1;
		if ( adj >= 0 ) {
			for ( j = 0; j < 4; j++ ) if ( qm -> mesh -> quadAdj [ 4 * adj + j ] == qID[0] )
				qm -> mesh -> quadAdj [ 4 * adj + j ] = qID[1];
		}
	}
	for ( i = 0 ; i < 4; i++ ) {
		j = qm -> mesh -> quadIdx [ 4 * (qID[0] - 1) + i ] - 1;
		qm -> mesh -> valence[ j ][1] = qID[0];
		j = qm -> mesh -> quadIdx [ 4 * (qID[1] - 1) + i ] - 1;
		qm -> mesh -> valence[ j ][1] = qID[1];
	}
	idx[0] = 0 ; idx[1] = 3; idx[2] = swap; idx[3] = swap + 3;
	for ( i = 0 ; i < 4; i++ ) {
		stat = setValence ( qm -> mesh , qg.verts[idx[i]] );
		if ( stat != EGADS_SUCCESS ) {
			printf(" Inside swapping operation set valence for %d --> %d !!!\n ", qg.verts[idx[i]], stat );
			return stat;
		}
	}
	//printMesh(qm, buffer,0);
	for ( j  = 0 ; j < 2; j++ ) {
		vc   = qg.verts[ j * 3 ];
		q    = 0 ; if ( EG_quadVertIdx ( qm ->mesh, qID[q], vc ) < 0 ) q = 1;
		stat = quadAngleOrientation(qm, qID[q], &area, ori, order, angles);
		if ( area == 1 && ori[0] == 1 && angles[0] <= PI ) continue;
		i0   = qg.verts[(3 * j + 1)%6];
		i1   = qg.verts[(3 * j + 5)%6];
		averageCoordsMinusLinks ( qm, vc, i0, i1 );
	}
	//printMesh(qm, buffer,0);
	stat = optimize_angles(qm, 6, qg.verts,0);
	return stat;
}

static int
EG_split (meshMap *qm, int qID, int *activity) {
	int poly[3], val[3], v,  q, id0 = 0, i, stat, dist = 0, validSplit = 0, id6[2], links[4];
	vStar *star = NULL;
	*activity   = 0;
	stat        = checkQuad (qm -> mesh, qID );
	if ( stat != EGADS_SUCCESS ) {
		printf(" EG_split for quad %d is %d \n ", qID, stat );
		return stat;
	}
	for ( v = 0; v < 4; v++ ) {
		poly[0]     = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + v ];
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
#ifdef DEBUG
			printf("WE FOUND A CANDIDATE FOR SPLITTING: SPLIT VERTEX V %d FROM LINKS %d  and %d dist %d \n",
					star -> verts[0], star -> verts[id0], star -> verts[star->idxV[id0 + dist]], dist);
			printVertSpecs ( qm -> mesh, star -> verts[0]);
			printVertSpecs ( qm -> mesh, star -> verts[id0]);
			printVertSpecs ( qm -> mesh, star -> verts[star->idxV[id0 + dist]]);
#endif
			//  printMesh(qm, buffer, 0);
			stat      = EG_splittingOperation(qm, star -> verts[0], star -> verts[id0], star -> verts[ star -> idxV[id0 + dist]]);
			*activity = 1;
			return stat;
		}
		EG_freeStar ( &star);
	}
	return EGADS_SUCCESS;
}


static int
EG_splittingOperation(meshMap *qm, int vC, int vL, int vR ) {
	int qIdx[4], modQ[4], verts[4], adj[2], poly[4], q, newQ, i, j, stat, id0 = -1, id1 = -1, dist, links[4], vals[4];
	vStar  *star  = NULL;
	double eval[18], uvAv[2];
	if ( qm -> mesh -> remQuads[0] > 0 ) {
		poly[3] = qm -> mesh -> remVerts[qm -> mesh -> remVerts[0]--];
		newQ    = qm -> mesh -> remQuads[qm -> mesh -> remQuads[0]--];
	}
	else {
		poly[3] = qm -> mesh -> totVerts + 1;
		newQ    = qm -> mesh -> totQuads + 1;
		if ( qm -> mesh -> totVerts > 2 * qm -> oriV ) {
			printf(" We have duplicated the number of initial vertices. This is too much. \n");
			return EGADS_SUCCESS;
		}
		++qm -> mesh -> totVerts; ++qm -> mesh -> totQuads;
	}
	qm -> mesh -> vType [poly[3] - 1] = -1;
	stat = EG_buildStar(qm -> mesh, &star, vC );
	if ( stat != EGADS_SUCCESS || star == NULL ) {
		printf(" In splittingOperation build star %d is NULL %d !!\n", vC, stat );
		return stat;
	}
	if ( stat != EGADS_SUCCESS || star == NULL ) return stat;
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
	if ( id0 == -1 || id1 == -1 ) return EGADS_INDEXERR;
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
	qm -> mesh -> quadIdx[ 4 * ( newQ - 1)    ] = poly[1];
	qm -> mesh -> quadIdx[ 4 * ( newQ - 1) + 1] = poly[0];
	qm -> mesh -> quadIdx[ 4 * ( newQ - 1) + 2] = poly[2];
	qm -> mesh -> quadIdx[ 4 * ( newQ - 1) + 3] = poly[3];
	for ( i = 0 ; i < 4; i++ ) printf(" v %d %d \n", i, qm -> mesh -> quadIdx[ 4 * ( newQ - 1)  +i  ]);
	for ( i = 0 ; i < 4; ++i) {
		modQ[i]  = star -> quads [ qIdx[i] ] ;
		qm -> mesh -> quadAdj[ 4 * ( newQ - 1) + i ] = modQ[i];
		qm -> mesh -> valence[ qm -> mesh -> quadIdx [ 4 * ( newQ - 1 )  + i ] -1 ] [1]  = newQ;
		if ( modQ[i] == -1 ) continue;
		stat = EG_adjQtoPair ( qm -> mesh, modQ[i], poly[0], verts[i], adj );
		qm -> mesh -> quadAdj [ 4 * ( modQ[i] - 1) + adj[0] ]  = newQ;
	}
	j = qIdx[2];
	q = star -> quads[j++];
	while ( q != star -> quads [ qIdx[0]] ) {
		for ( i = 0 ; i < 4; ++i)
			if (qm -> mesh -> quadIdx[ 4 * (q - 1) + i] == poly[0] )
				qm -> mesh -> quadIdx[ 4 * (q - 1) + i] = poly[3];
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
			return stat;
		}
	}
	// Update valences at links
	for ( i = 0 ; i < 4; i++ ) {
		for ( j = 0 ; j < qm -> mesh -> valence [ poly[i] - 1][0]; j++ ) {
			stat = setValence( qm -> mesh, qm -> mesh -> valence [ poly[i] - 1 ][ 2 + j ] );
			if ( stat != EGADS_SUCCESS ) {
				return stat;
			}
		}
	}
	// New Vertex location: take mid points of quads to which they belong
	for ( j = 0 ; j < 3; j++ ) {
		if ( j < 2 ) qm -> mesh -> uvs[ 2 * ( poly[3] - 1 ) + j ] = qm -> mesh -> uvs[ 2 * ( poly[0] - 1 ) +j];
		qm -> mesh -> xyzs[ 3 * ( poly[3] - 1 ) + j ] = qm -> mesh -> xyzs[ 3 * ( poly[0] - 1 ) +j];
	}
	for ( j = 0 ; j < 4; j = j + 3) {
		if ( qm -> mesh -> vType [ poly[j] -1 ] != -1 ) continue;
		stat = EG_buildStar (qm -> mesh, &star, poly[j] );
		if ( stat != EGADS_SUCCESS || star == NULL ) return stat;
		uvAv[0] = qm -> mesh -> uvs[ 2 * ( poly[j] - 1 )    ];
		uvAv[1] = qm -> mesh -> uvs[ 2 * ( poly[j] - 1 ) + 1];
		q       = 1;
		for ( i = 0 ; i < star -> nQ; i++ ) {
			id0 = star -> verts [ 2 * i + 1] ;
			if ( id0 == poly[1] || id0 == poly[2] ) continue;
			uvAv[0] += qm -> mesh -> uvs[ 2 * ( id0 - 1 )    ];
			uvAv[1] += qm -> mesh -> uvs[ 2 * ( id0 - 1 ) + 1];
			q++;
		}
		uvAv[0] /= (double) q;
		uvAv[1] /= (double) q;
		stat = EG_evaluate ( qm -> face, uvAv, eval );
		for ( i = 0 ; i < 3; i++ ) {
			if ( i < 2 ) qm -> mesh -> uvs [ 2 * ( poly[j] - 1) + i] = uvAv[i];
			qm -> mesh -> xyzs [ 3 * ( poly[j] - 1) + i] = eval[i];
		}
	}
	EG_freeStar (&star);
	stat = optimize_angles(qm, 4, poly, 0);
	for ( i = j = 0 ; i < qm -> mesh -> totQuads; i++ ) {
		if ( checkQuad (qm -> mesh, i + 1) == EGADS_SUCCESS ) j++;
	}
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
		if ( stat != EGADS_SUCCESS ) printf(" EG_doubleSwap adjToPair %d --> %d !!\n", adj[1], stat );
		return stat;
	}
	i       = EG_quadVertIdx ( qm -> mesh, adj[1], qg.verts[vopp3]);
	adjPiv5 = qm ->mesh -> quadIdx [ 4 * ( adj[1] -1 ) + ( i + 1)%4];
	if (adjPiv5 == qg.verts[piv5] )
		adjPiv5 = qm ->mesh -> quadIdx [ 4 * ( adj[1] -1 ) + ( i + 3)%4];
	if (     swap == 0 &&  getValence (qm ->mesh, adjPiv5) > 4 ) return EGADS_SUCCESS;
	else if (swap == 1 && (getValence (qm ->mesh, adjPiv5) <  5 ||
			(qm -> mesh -> vType[adjPiv5 -1 ] != -1 &&
					qm -> mesh -> vType[qg.verts[vopp3] -1 ] !=-1 ) )) return EGADS_SUCCESS;
	piv5       = qg.verts[0];
	stat       = EG_swappingOperation (qm, qg, piv3 );
	*activity  = 1;
	if ( stat != EGADS_SUCCESS ) {
		printf(" EG_doubleSwap: at first swap --> %d !!\n ", stat );
		return stat;
	}
	stat = EG_swap (qm, adj[1], &i );
	if ( stat != EGADS_SUCCESS || i == 0 ) {
		// undo swap and leave
		stat = EG_createQuadGroup (qm -> mesh, &qg, qg.q[0], qg.q[1]);
		for ( i = 0 ; i < 6; i++ ) if ( qg.verts[i] == piv5) break;
		stat       = EG_swappingOperation (qm, qg, i );
		if ( stat == EGADS_SUCCESS ) *activity = 0;
	}
	return stat;
}

static int EG_doubleSplit(meshMap *qm, quadGroup qg, int forcing, int *activity ) {
	int i, j, stat, f = 0 ;
	*activity  = 0;
	int piv[2] = {1, 5} ;
	if ( qg.vals[1] != 3 ) {
		piv[0] = 5; piv[1] = 1;
	}
	if ( qg.vals[piv[0] ] != 3 ) return EGADS_SUCCESS;
	if ((qm -> mesh -> vType [ qg.verts[0]      - 1 ] >= 0 &&
			qm -> mesh -> vType [ qg.verts[piv[1]] - 1  ] >= 0 ) ||
			qm -> mesh -> vType [ qg.verts[0]      - 1  ] == 4) return EGADS_SUCCESS;
	if ((forcing == 0 && (qg.vals[0] != 5     || qg.vals[1] * qg.vals[5]      != 15)) ||
			(forcing == 1 && (qm ->extraQuads > 0 || qg.vals[0] * qg.vals[piv[1]] <= 16))) return EGADS_SUCCESS;
	*activity  = 1;
	if ( qg.vals[1] * qg.vals[5] != 15 || qg.vals[0] != 5 ) f = 1;
	printf(" DOUBLE SPLIT %d \n ", qg.verts[0] );
	stat       = EG_splittingOperation (qm, qg.verts[0], qg.verts[piv[0]], qg.verts[piv[1]]);
	if ( stat != EGADS_SUCCESS ) {
		printf("In EG_doubleSplit: force 1st split through %d - %d --> %d !!\n ",
				qg.verts[0], qg.verts[piv[0]], stat );
		return stat;
	}
	for ( j = 0 ; j < 2; j++)
		if ( EG_quadVertIdx (qm -> mesh, qg.q[j], qg.verts[piv[1]] ) >= 0 ) break;
	stat = EG_split (qm, qg.q[j], &i);
	if ( i == 0 || stat != EGADS_SUCCESS ) {
		if ( i == 0 ) {
			// undo split by collapsing
			stat = EG_forceCollapse (qm, qm ->mesh -> valence[ qg.verts[0] - 1][1], &i );
			if ( i > 0 && stat == EGADS_SUCCESS ) {
				*activity = 0;
			}
		}
	}
	if ( f == 1 && *activity > 0 ) {
		qm -> extraQuads += 2;
		//	fprintf(stderr,"doubleSplit forcing extra quads %d\n ", qm -> extraQuads );
	}
	return stat;
}

static int EG_swapSplit(meshMap *qm,quadGroup qg, int forcing, int *activity  ) {
	int  stat, i, j, i3 = -1, i5 = -1, v3opp = -1, q5, vL5, vL5adj, swap = 0, adj[2], f = 0;
	vStar * star = NULL;
	*activity = 0;
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
			stat              = EG_adjQtoPair (qm -> mesh, qg.q[i], qg.verts[j], qg.verts[vL5], adj );
			if ( stat        != EGADS_SUCCESS || adj[1] == -1 ) continue;
			q5                = EG_quadVertIdx (qm -> mesh, adj[1], qg.verts[j]);
			vL5adj            = qm -> mesh -> quadIdx [ 4 * ( adj[1] - 1 ) + ( q5 + 1 ) %4 ];
			if ( vL5adj == qg.verts[vL5] )
				vL5adj = qm -> mesh -> quadIdx [ 4 * ( adj[1] - 1 ) + ( q5 + 3) %4 ];
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
		if ( stat != EGADS_SUCCESS || adj[1] == -1 ) {
			if ( stat != EGADS_SUCCESS) printf("EG_swapSplit: EG_adjQtoPair from quad %d is %d\n!!", qg.q[q5],stat );
			return stat;
		}
		i        = EG_quadVertIdx (qm -> mesh, adj[1], qg.verts[v3opp]);
		vL5adj   = qm -> mesh -> quadIdx [ 4 * ( adj[1] - 1 ) + ( i + 1 ) %4 ];
		if ( vL5adj == qg.verts[vL5] )
			vL5adj = qm -> mesh -> quadIdx [ 4 * ( adj[1] - 1 ) + ( i + 3) %4 ];
		if ( i3 != -1 && ( qg.vals[v3opp] == 5 || getValence (qm ->mesh, vL5adj) == 3 ) ) swap = i3;
		else if ( forcing == 1 && qm -> extraQuads <= 0 && ( qg.vals[v3opp] == 5 || getValence (qm ->mesh, vL5adj) == 3 ) ) swap = v3opp;
	}
	if ( swap %3 == 0 ) return EGADS_SUCCESS;
	stat       = EG_swappingOperation( qm, qg, swap );
	*activity  = 1;
	if ( stat != EGADS_SUCCESS) {
		printf(" In swapSplit thru %d : EG_swappingOperation went %d !!\n ", swap, stat );
		printQuadGroup(qm ->mesh, qg);
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
		printf("In swapSplit vertex %d is not linked to %d !!\n", qg.verts[3], qg.verts[i5]);
		return EGADS_INDEXERR;
	}
	adj[0] = star -> verts[ star -> idxV[2 * i + 1 + 6] ];
	if ( star -> nQ == 5 ) {
		adj[1] = star -> verts[ star -> idxV[2 * i + 1 + 4] ];
		if ( getValence ( qm -> mesh, adj[1] ) < getValence ( qm -> mesh, adj[0] ) ) adj[0] = adj[1];
	}
	printf(" SwapSplit Call splitting operation from %d  \n",qg.verts[i5]);
	if ( qm -> mesh -> vType [ qg.verts[i5] - 1] == 4 ) {
		// undo swapping: this is a sharp vertex. don't split thrugh here
		j = qg.verts[0];
		stat = EG_createQuadGroup (qm -> mesh, &qg, qg.q[0], qg.q[1] );
		for ( i = 0 ; i < 6; i++ )
			if ( qg.verts[i] == j ) break;
		stat       = EG_swappingOperation( qm, qg, i );
		if ( stat == EGADS_SUCCESS )*activity  = 0;
		return stat;
	}
	stat = EG_splittingOperation (qm, qg.verts[i5], qg.verts[3], adj[0]);
	if ( stat != EGADS_SUCCESS) {
		printf("In swapSplit splittingOperation --> %d !!\n", stat);
		return EGADS_INDEXERR;
	}
	EG_freeStar(&star);
	if ( f == 1 && *activity > 0 ) {
		qm -> extraQuads++;
		//	fprintf(stderr,"swapSplit forcing extra quads %d\n ", qm -> extraQuads );
		printf("swapSplit forcing extra quads %d\n ", qm -> extraQuads );
	}
	printf(" LEAVE SWAP SPLIT \n ");
	//printMesh(qm, buffer,0);
	return stat;
}


static int EG_swapCollapse (meshMap *qm,quadGroup qg, int forcing, int *activity  ) {
	int  stat, i, i3 = -1, q5, qC, vL5, vL5adj, vOpp3, swap = 0, adj[2], f = 0;
	*activity   = 0;
	if (validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0 ) return EGADS_SUCCESS;
	if      (qg.vals[1] * qg.vals[2] == 15 ) {
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
			stat       = EG_swappingOperation (qm, qg, 1 );
			*activity  = 1;
			if ( stat != EGADS_SUCCESS ) {
				printf("forcing swapcollapse:: after swapping %d \n ", stat);
				return stat;
			}
			qC = qg.q[0];
			if ( EG_nValenceCount ( qm -> mesh, qC, 3 ) < 2 ) qC = qg.q[1];
			stat = EG_forceCollapse (qm, qC, &i);
			if ( i > 0 ) {
				qm ->extraQuads--;
				//		fprintf(stderr,"swapCollapse forcing extra quads %d\n ", qm -> extraQuads );
			}
			if ( stat != EGADS_SUCCESS )printf("forcing swapcollapse:: after swapping %d \n ", stat);
			return stat;
		}
		else return EGADS_SUCCESS;
	}
	else return EGADS_SUCCESS;
	q5     = 0; if ( i3 > 3 ) q5 = 1;
	if (validCollapse (qm -> mesh, qg.q[q5], qg.verts[i3]) == 0) return EGADS_SUCCESS;
	vOpp3  = (i3 + 3 ) % 6;
	q5     = 0; if ( vOpp3 > 3 ) q5 = 1;
	qC     = qg.q[ ( q5 +1)%2];
	vL5    = ( vOpp3 + 1) %6;
	if ( vL5 %3 != 0 ) vL5 = ( vOpp3 + 5) %6;
	if ( qg.vals[vOpp3] == 3 ) return EGADS_SUCCESS;
	stat       = EG_adjQtoPair (qm -> mesh, qg.q[q5], qg.verts[vOpp3], qg.verts[vL5], adj );
	if ( stat != EGADS_SUCCESS || adj[1] == -1 ) return stat;
	i      = EG_quadVertIdx ( qm -> mesh, adj[1], qg.verts[vOpp3]);
	vL5adj = qm -> mesh -> quadIdx [ 4 * ( adj[1] - 1 ) + ( i + 1 ) %4 ];
	if (  vL5adj == qg.verts[vL5] ) vL5adj = qm -> mesh -> quadIdx [ 4 * ( adj[1] - 1 ) + ( i + 3 ) %4 ];
	if (( forcing == 0 && qg.vals[vL5] == 4 && getValence (qm -> mesh, vL5adj) >= 4 ) ||
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
	printf(" SWAP COLLAPSE \n ");
	printQuadGroup (qm ->mesh, qg);
	stat       = EG_swappingOperation (qm, qg, swap );
	*activity  = 1;
	if ( stat != EGADS_SUCCESS ) {
		printf("EG_swapCollapse after swapping %d !!\n", swap );
		return stat;
	}
	vL5  = qg.verts[0];
	stat = EG_forceCollapse (qm, qC, &i);
	if ( stat != EGADS_SUCCESS || i > 0 ) {
		if ( i > 0 && f == 1 ) {
			qm -> extraQuads--;
			//	fprintf(stderr,"swapCollapse forcing extra quads %d\n ", qm -> extraQuads );
			printf("swapCollapse forcing extra quads %d\n ", qm -> extraQuads );
		}
		return stat;
	}
	// undo swap and leave
	stat    = EG_createQuadGroup (qm -> mesh, &qg, qg.q[0], qg.q[1]);
	if ( stat != EGADS_SUCCESS ) {
		printf("EG_swapCollapse undo swap: EG_createQuadGroup --> %d !!\n", stat );
		printQuadGroup (qm -> mesh, qg );
		return stat;
	}
	for ( i = 0 ; i < 6; i++ ) if ( qg.verts[i] == vL5 ) break;
	stat      = EG_swappingOperation (qm, qg, i );
	*activity = 0;
	return stat;
}


static int EG_cleanNeighborhood (meshMap *qm, int qID,  int transfer, int *activity ) {
	int i, act, stat = 0, count = 0, j, v[4];
	vStar *star = NULL;
	*activity   = 0;
	stat        = checkQuad (qm -> mesh, qID );
	if ( stat  != EGADS_SUCCESS ) {
		printf("EG_cleanNeighborhood : Quad %d is invalid = %d \n", qID, stat );
		return checkMesh (qm);
	}
	stat       = restoreMeshData(qm, qm -> backupMesh, qm -> mesh);
	if ( stat != EGADS_SUCCESS ) return stat;
	for ( i = 0 ; i < 4; i++ ) v[i]  = qm -> mesh -> quadIdx[ 4 * ( qID -1 ) + i ];
	stat       = EG_cleanQuad ( qm, qID, 1, transfer, 0, &act );
	if ( stat != EGADS_SUCCESS ) {
		printf("\n\nInside EG_cleanNeighborhood: EG_cleanQuad %d --> %d!!\n", qID, stat );
		stat = restoreMeshData(qm, qm -> mesh, qm -> backupMesh );
		if ( stat != EGADS_SUCCESS ) return stat;
		act = 0;
	}
	if ( act > 0 ) stat = restoreMeshData(qm, qm -> backupMesh, qm -> mesh);
	if ( stat != EGADS_SUCCESS ) return stat;
	*activity += act;
	for ( i = 0 ; i < 4; i++ ) {
		if ( checkVertex ( qm -> mesh, v[i] ) != EGADS_SUCCESS ) continue;
		stat       = EG_buildStar (qm -> mesh, &star, v[i]);
		if ( star == NULL || stat != EGADS_SUCCESS ) {
			printf("EG_cleanNeighborhood : EG_buildStar  = %d \n", stat );
			return stat;
		}
		for ( j = 0 ; j < star -> nQ; j++ ) {
			if ( star -> quads[j] == -1 || checkQuad ( qm -> mesh, star -> quads[j]) != EGADS_SUCCESS ) continue;
			if ( EG_cleanQuad (qm, star -> quads[j], 1, transfer, 0, &act ) != EGADS_SUCCESS ) {
				stat = restoreMeshData (qm, qm -> mesh, qm -> backupMesh);
				if ( stat != EGADS_SUCCESS ) {
					EG_freeStar(&star);
					printf(" In cleanNeighborhood:: restoreMeshData inside star loop %d \n ", stat);
					return stat;
				}
				act  = 0;
			}
			if ( act > 0 ) {
				stat = restoreMeshData(qm, qm -> backupMesh, qm -> mesh);
				if ( stat != EGADS_SUCCESS ) {
					printf(" In cleanNeighborhood:: restoreMeshData inside star loop:: update %d \n ", stat);
					EG_freeStar ( &star );
					return stat;
				}
				count += act;
			}
		}
		EG_freeStar ( &star );
	}
	*activity += count;
	return checkMesh (qm );
}

static int
EG_wvsData(meshData *mesh, char *buffer) {
	int i;
	FILE *fil = NULL;
	printf(" Write on WVS DATA FILE %s \n ", buffer);
	fil = fopen (buffer, "w" );
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
				mesh -> quadIdx[4 * i    ],
				mesh -> quadIdx[4 * i + 1],
				mesh -> quadIdx[4 * i + 2],
				mesh -> quadIdx[4 * i + 3]);
	}
	fclose (fil);
	printf(" All good for %s\n ", buffer);
	return EGADS_SUCCESS;
}



static int EG_fullMeshRegularization(meshMap *qm )
{
	int     ITMAX, it = 0, stat, activity = 0, totActivity = 0, loopact, i, j, k, kk, q, transfer = 0;
	int     best_iV, iV, quadPair[2], prevPair[2], totV, vQ;
	// GET RANGE FOR EACH POINT
	stat       = checkMesh (qm);
	if ( stat != EGADS_SUCCESS) {
		fprintf(stderr,"At meshRegularization: starting mesh is invalid --> %d !! I'm leaving program \n ", stat);
		return stat;
	}
	meshCount (qm -> mesh, &iV, &totV, &vQ);
	fprintf(stderr,"Original mesh has %d QUADS, %d / %d irregular vertices (%f percent) ============\n", vQ, iV, totV, (double) iV * 100.0 / (double)totV );
	ITMAX = 100;
	stat  = restoreMeshData (qm, qm -> backupMesh, qm -> mesh);
	it    = 0;
	do {
		totActivity   = 0;
		for (i = 0 ; i < qm -> mesh -> totQuads; i++ ) {
			if ( checkQuad (qm ->mesh, i + 1) != EGADS_SUCCESS ) continue; //can be a deleted quad
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
	fprintf(stderr,"*******************************************************************\n");
	fprintf(stderr," AFTER BASIC: mesh has %d QUADS, %d / %d irregular vertices (%f percent) ============\n", vQ, iV, totV, (double) iV * 100.0 / (double)totV );
	fprintf(stderr, "*******************************************************************\n");
	snprintf(buffer,500,"BASIC_%d", qm -> fID );
	printMesh ( qm, buffer, 1 );
	stat = restoreMeshData (qm, qm -> bestMesh, qm -> mesh);
	if ( stat != EGADS_SUCCESS ) {
		printf(" restoreMeshData for qb -> bestMesh after basic %d !!\n ", stat );
		return stat;
	}
	best_iV = iV;
	if ( iV > 2 ) {
		stat       = restoreMeshData (qm, qm -> backupMesh, qm -> mesh);
		stat      += optimize_angles (qm, 0, NULL, 2);
		if ( stat != EGADS_SUCCESS ) {
			printf("Before applying transfer valences: restore + optimize %d !!\n", stat);
			return stat;
		}
		for ( k = 0 ; k <= 2; k++ ) {
			kk = 0;
			if ( k >= 2 ) {
				kk = 1;
				qm -> extraQuads = 0;
			}
			it  = 0 ;
			do {
				totActivity = 0;
				stat = checkMesh (qm);
				stat = restoreMeshData (qm, qm -> backupMesh, qm -> mesh);
				printf("********************  ITERATION %d : TOTACT %d : iV  %d  MESH %d ************************\n",
						it, totActivity, iV, MESHPLOTS);
				prevPair[0] = -1;
				prevPair[1] = -1;
				for (q = 0 ; q < qm -> mesh-> totQuads; q++) {
					transfer = 0 ;
					if ( checkQuad( qm -> mesh, q + 1 ) != EGADS_SUCCESS ) continue;
					if ( q + 1 != prevPair[0] && q + 1 != prevPair[1] ) {
						quadPair[0] = q + 1;
						quadPair[1] = -1;
					}
					else continue;
					printf(" Quad %d %d previous %d  %d \n ", quadPair[0], quadPair[1], prevPair[0], prevPair[1] );
					if ( EG_quadIsBoundary ( qm ->mesh, q+ 1 ) != 0 && k == 0  ) continue;
					stat           = EG_transferValences ( qm, quadPair, kk, &transfer, &activity );
					if (stat      != EGADS_SUCCESS ) {
						stat       = restoreMeshData (qm, qm -> mesh, qm -> backupMesh);
						if ( stat != EGADS_SUCCESS ) return stat;
						activity  = 0;
					}
					if ( activity > 0 ) printMesh(qm, buffer,0);
					if ( activity == 0 || quadPair[0] == -1 || checkQuad( qm -> mesh, quadPair[0] ) != EGADS_SUCCESS  ) continue;
					totActivity += activity;
					j = 0;
					while ( activity > 0 && iV > 2 && j < 20 ) {
						activity  = 0;
						for ( i = 0 ; i < 2; i++ ) {
							if ( checkQuad (qm -> mesh, quadPair[i] ) != EGADS_SUCCESS ) continue;
							stat       = EG_cleanNeighborhood (qm, quadPair[i],  transfer, &loopact);
							if ( stat != EGADS_SUCCESS ) return stat;
							activity +=loopact;
						}
						if ( activity > 0 ) break;
						stat       = EG_transferValences ( qm, quadPair, kk, &transfer, &activity );
						if ( stat != EGADS_SUCCESS || checkQuad (qm -> mesh, quadPair[0] ) != EGADS_SUCCESS  ) {
							stat       = restoreMeshData (qm, qm -> mesh, qm -> backupMesh);
							if ( stat != EGADS_SUCCESS ) return stat;
							break;
						}
						j++;
					}
					if ( iV < best_iV ) {
						stat = restoreMeshData ( qm, qm -> bestMesh, qm -> mesh);
						if ( stat != EGADS_SUCCESS )  return stat;
						best_iV = iV;
					}
					meshCount(qm -> mesh, &iV, &totV, &vQ );
					if ( iV <=2 ) break;
					prevPair[0] = quadPair[0];
					prevPair[1] = quadPair[1];
				}
				prevPair[0] = -1;
				prevPair[1] = -1;
				stat        = checkMesh (qm );
				if ( stat != EGADS_SUCCESS ) return stat;
				meshCount(qm -> mesh, &iV, &totV, &vQ );
				if ( iV  < best_iV ) {
					restoreMeshData ( qm, qm -> bestMesh, qm -> mesh);
					best_iV = iV;
				}
				if ( iV <=2 ) break;
				++it;
			} while ( totActivity > 0 && it < ITMAX && iV > 2);
			stat = checkMesh (qm );
			if ( stat != EGADS_SUCCESS )  return stat;
			meshCount(qm -> mesh, &iV, &totV, &vQ );
			if ( iV  < best_iV ) {
				restoreMeshData ( qm, qm -> bestMesh, qm -> mesh);
				best_iV = iV;
			}
			if ( iV <=2 ) break;
		}
		if ( iV > best_iV ) {
			stat = restoreMeshData (qm, qm -> mesh, qm -> bestMesh);
			printf(" Best round had %d quads compared to current %d \n", iV, best_iV);
			fprintf(stderr, " Best round had %d quads compared to current %d \n", iV, best_iV);
			if ( stat != EGADS_SUCCESS ) return stat;
		}
		meshCount(qm -> mesh, &iV, &totV, &vQ );
	}
	fprintf(stderr,"*******************************************************************\n");
	fprintf(stderr," FULL REGULARIZATION: mesh has %d QUADS, %d / %d irregular vertices (%f percent) ============\n", vQ, iV, totV, (double) iV * 100.0 / (double)totV );
	fprintf(stderr, "*******************************************************************\n");
	fprintf(stderr," FINAL MESH HAS %d QUADS AND %d are IRREGULAR ( %d total) EXTRAQUADS %d \n",
			vQ, iV, totV,  qm -> extraQuads);
	printf(" FINAL MESH HAS %d QUADS AND %d are IRREGULAR ( %d total) EXTRAQUADS %d \n",
			vQ, iV, totV,  qm -> extraQuads);
	stat = resizeQm (qm) ;
	if ( stat != EGADS_SUCCESS ) {
		printf(" After resizing mesh:: %d \n", stat );
		return stat;
	}
	printMeshStats(qm, 2);
	snprintf(buffer,500,"wvsRegular_%d.txt", qm -> fID);
	stat = EG_wvsData(qm -> mesh, buffer);
	if ( stat != EGADS_SUCCESS ) {
		printf(" writing in wvs file %d !! \n ", stat );
		return stat;
	}
	snprintf(buffer,500,"BEFORE_GLOBAL%d.txt", qm -> fID);
	printMesh(qm, buffer, 1);
	stat  += optimize_angles(qm, 0, NULL, 1);
	snprintf(buffer,500,"finalMesh_%d.txt", qm ->fID);
	printMesh (qm, buffer, 1 );
	snprintf(buffer,500,"wvsFinalMesh_%d.txt", qm ->fID);
	stat += EG_wvsData(qm -> mesh, buffer);
	return stat;
}


static int optimize_angles(meshMap *qm, int nP, /*@unused@*/ /*@null@*/int *pList, int fullRegularization)
{
	int  i, v, q, stat = EGADS_SUCCESS, nV = 0, it = 0, itMax = 100, *vID = NULL;
	double dtheta, minT, maxT;
	vStar *star = NULL;
	if ( fullRegularization == 0 ) { // move around only affected vertices
#ifdef DEBUG
		vID = (int *) EG_alloc ( qm -> mesh  -> totVerts * sizeof ( int ) );
		if ( vID == NULL ) return EGADS_MALLOC;
		printf(" POINT LIST %d \n", nP );
		for ( i = 0 ; i < nP; i++ ) {
			printf(" P %d = %d \n", i, pList[i]);
		}
		printMesh(qm, buffer, 0);
#endif
		if ( nP == 0 || pList == NULL ) {
			EG_free (vID );
			return EGADS_SUCCESS;
		}
		for ( nV = i = 0 ; i < nP; ++i) {
			if (  qm -> mesh  -> vType[ pList[i] - 1 ] == - 1 ) vID[nV++] = pList[i];
		}
		for ( v = 0 ; v < nV; v++) {
			if ( vertexLinksToBounds (qm, vID[v] ) == 1 ) {
				stat = validBoundary (qm, MINVALIDANGLE, MAXVALIDANGLE);
				if ( stat != EGADS_SUCCESS ) {
					EG_free (vID );
					return stat;
				}
			}
		}
		for ( v = 0 ; v < nV; v++) {
			if ( vertexLinksToBounds (qm, vID[v] ) == 1 ) continue;
			stat = EG_buildStar (qm -> mesh, &star, vID[v] );
			if ( stat != EGADS_SUCCESS || star == NULL ) return stat;
			for ( i = 0 ; i < star -> nQ; i++ ) {
				stat = validInterior ( qm, star -> quads[i], MINVALIDANGLE, MAXVALIDANGLE);
				if ( stat != EGADS_SUCCESS ) {
					printf(" Mesh is invalid for MIN MAX TOTALS [0, 200 deg]\n");
					printMesh(qm , buffer, 1);
					EG_free(vID);
					return stat;
				}
			}
			EG_freeStar(&star);
		}
		EG_free (vID);
		return EGADS_SUCCESS;
	}
	else {
		//EG_smoothQuads (qm, 50 );
		itMax  =  50;
		minT   = 0.0;
		maxT   = PI;
		dtheta = 0.5 * PI / (double)itMax;
		printMesh(qm, buffer,0);
		printf("FULL  OPTIMIZE ANGLES -----------------------------   \n ");
		while ( it <= itMax ) {
			for ( q = v = 0 ; v < qm -> mesh -> totQuads; v++) {
				if (checkQuad ( qm -> mesh, v + 1) != EGADS_SUCCESS ||
						EG_quadIsBoundary ( qm -> mesh, v + 1) == 0 ) continue;
				printf(" Valid Interior for quad %d -> %f %f\n ", v+1, minT, maxT);
				stat = validInterior ( qm, v + 1, minT, maxT);
			}
			stat = validBoundary (qm, minT, maxT);
			if ( stat != EGADS_SUCCESS ) return stat;
			printMesh (qm, buffer,0);
			printf("\n\n INTERIOR  OPTIMIZE ANGLES -----------------------------   \n ");
			for ( q = v = 0 ; v < qm -> mesh -> totQuads; v++) {
				if (checkQuad ( qm -> mesh, v + 1) != EGADS_SUCCESS ||
						EG_quadIsBoundary ( qm -> mesh, v+ 1) != 0 ) continue;
				printf(" Valid Interior for quad %d -> %f %f\n ", v+1, minT, maxT);
				stat = validInterior ( qm, v + 1, minT, maxT);
				if ( stat == EGADS_SUCCESS ) q++;
				else {
					stat = validInterior ( qm, v + 1,MINVALIDANGLE, MAXVALIDANGLE);
					if ( stat != EGADS_SUCCESS ) {
						printf(" Quad %d is invalid for interval [ %f ,%f ]\n", v+1, MINVALIDANGLE, MAXVALIDANGLE);
						snprintf(buffer,500, "badmesh_%d.txt", qm -> fID);
						printMesh(qm , buffer, 1);
					}
				}
			}
			if ( stat == EGADS_SUCCESS ) stat = validBoundary (qm, MINVALIDANGLE, MAXVALIDANGLE);
			if ( stat != EGADS_SUCCESS || q == 0 ) break;

			minT += dtheta;
			maxT -= dtheta;
			if ( maxT < minT ) break;
			it++;
		}
		if ( stat != EGADS_SUCCESS) {
			snprintf(buffer,500, "badmesh_%d.txt", qm -> fID);
			printMesh(qm , buffer, 1);
			snprintf(buffer,500, "wvsBad_%i.txt", qm -> fID );
			i = EG_wvsData(qm -> mesh , buffer);
		}
		else if ( fullRegularization == 1 ) {
			snprintf(buffer,500, "centroid_%d.txt", qm -> fID);
			printMesh(qm , buffer, 1);
			snprintf(buffer,500, "wvscentroid_%i.txt", qm -> fID );
			i = EG_wvsData(qm -> mesh , buffer);
		}
	}
	if (stat != EGADS_SUCCESS ) printf(" MESH VALIDIDTY %d \n ", stat );
	printMesh(qm, buffer, 0);
	return stat;
}





//#ifdef STANDALONE
int main(int argc, char *argv[])
{
	clock_t      start_t, end_t, total_t;
	int          stat = 0,  f , i, j, iBody, oclass, mtype;
	int          atype, alen, *senses;
	const int    *ints;
	float        arg;
	double       box[6], size, params[3];
	const double *reals ;
	const char   *OCCrev, *string;
	ego          context, tess, model, geom, *bodies, *dum;
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
				fprintf(stderr,"===================================================================\n");
				fprintf(stderr,"===================================================================\n");
				fprintf(stderr," FULL MESH REGULARIZATION FACE %d \n ", f + 1);
				fprintf(stderr,"===================================================================\n");
				fprintf(stderr,"===================================================================\n");
				printf("===================================================================\n");
				printf("===================================================================\n");
				printf(" FULL MESH REGULARIZATION FACE %d \n ", f + 1);
				printf("===================================================================\n");
				printf("===================================================================\n");
				stat = EG_fullMeshRegularization(bodydata[iBody].qm[f] );
				fprintf(stderr, " EG_fullMeshRegularization face %d / %d = %d \n ", f + 1, bodydata[iBody].nfaces,  stat );
				if ( stat != EGADS_SUCCESS ) return stat;
				snprintf(buffer,500,"FINALMESH_%d", f + 1);
				printMesh(bodydata[iBody].qm[f],buffer, 1);
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
	fprintf(stderr, "Total time taken by CPU: %f\n", (double ) total_t / ( 60.0 * CLOCKS_PER_SEC) );
	fprintf(stderr, "Total M FILES %d\n", MESHPLOTS);
	EG_close ( context );
	return 0;
}
//#endif



/********************   IO FUNCTIONS **********************/

static void
printQuadSpecs(meshData *mesh, int id) {
	int i, v = 0 , val = 0 ;
	i = checkQuad ( mesh, id );
	if ( i != EGADS_SUCCESS ) return ;
	--id;
	printf(" QUAD %d HAS VERTICES ",id +1);
	for ( i = 0 ; i < 4; ++i ) {
		v   = mesh -> quadIdx[4 * id + i];
		val = checkVertex ( mesh , v ) ;
		if ( val != EGADS_SUCCESS) {
			printf(" HOLD ON: QUAD %d is fine but vertex %d = %d\n", id+ 1, v, val);
		}
		val = getValence ( mesh, v );
		printf(" %d ( val = %d )  ",v, val);
	}
	printf("\t AND ADJACENT QUADS ");
	for ( i = 0 ; i < 4; ++i ) printf(" %d ",mesh -> quadAdj[4*id + i]);
	printf("\n");
}

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



static void
printMeshStats(meshMap *qm, int sweep) {
	int i,len, val ;
	int intVal[100], boundVal[100];
	FILE *fout;
	snprintf(buffer, sizeof(char) * 32, "MESH_STATS_%d.txt", sweep);
	printf(" WRITING ON FILE %s\n",buffer);
	fout = fopen(buffer, "w");
	if ( fout == NULL ) return;
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
printMesh(meshMap *qm, char *name, int usename ) {
	int i,k, v, d;
	double eval[18], average[5], dist;
	FILE *fout = NULL;
	if ( usename == 0 ) snprintf ( name, 100,"M_%d",MESHPLOTS++ ) ;
	printf(" Writing in File %s  ", name);
	fout = fopen(name, "w" );
	if ( fout == NULL ) return;
	for ( i = 0 ; i < qm -> mesh -> totQuads; ++i) {
		if ( checkQuad ( qm -> mesh, i + 1) != EGADS_SUCCESS )
			continue;
		for ( k = 0; k < 4; ++k) {
			v  =   qm -> mesh -> quadIdx[4*i + k] - 1;
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
				return;
			}
		}
		v  =   qm -> mesh -> quadIdx[4*i ] - 1;
		fprintf(fout, "%lf %lf %lf %d %lf %lf\n",qm -> mesh->xyzs[3*v], qm -> mesh->xyzs[3*v +1], qm -> mesh->xyzs[3*v + 2],
				v + 1, qm -> mesh -> uvs[2*v] , qm -> mesh -> uvs[2*v + 1]);
		fprintf(fout,"\n\n");
		quadAverageCoords(qm, i + 1, average, &average[2]);
		fprintf(fout, "%lf %lf %lf %d %lf %lf\n", average[2], average[3], average[4], i + 1, average[0], average[1]) ;
		fprintf(fout,"\n\n");
	}
	fclose(fout);
}



static void
printStar(vStar *star) {
	int i;
	printf(" ===== STAR CENTRED AT VERTEX %d \n", star->verts[0]);
	printf(" ROSA  DE %d PICOS \n",star->nV);
	for ( i = 1 ; i < star->nV; ++ i) {
		if ( i%2 == 0)
			printf(" OPP  PICO %d\n",star->verts[i]);
		else
			printf(" LINK  PICO %d\n",star->verts[i]);
	}
	for ( i = 0 ; i < star ->nQ; ++i) printf(" QUAD %d = %d\n",i,star->quads[i]);
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



void sampleNormalPlane(meshMap *qm, char *name, int vID) {
	double min[3], max[3], c, p[3], dt[3], aux, point[18], normal[4];
	int  i, j, k, nP;
	FILE *f;
	printf(" WRITIN IN %s \n", name);
	f = fopen(name,"w");
	if ( f == NULL) return ;
	nP = 20;
	c = 0.0;
	vID--;
	i = EG_evaluate(qm -> face, &qm -> mesh -> uvs[2 * vID], point);
	if ( i != EGADS_SUCCESS) return ;
	//
	if ( qm -> face -> mtype == SREVERSE )
		cross_product(&point[6], &point[3], normal);
	else
		cross_product(&point[3], &point[6], normal);
	for (i = 0 ; i < 3; ++i) {
		min  [i] = qm -> mesh ->  xyzs[3*vID + i] -1; //qm -> xyzs[3*vID + i] - 100 *qm -> xyzs[3*vID + i] ;
		max  [i] = qm -> mesh -> xyzs[3*vID + i]  +1;//0.5;// qm -> xyzs[3*vID + i] + 100 *qm -> xyzs[3*vID + i] ;
		if ( min[i] > max[i]) {
			dt[0] = min[i];
			min[i] = max[i];
			max[i] = dt[0];
		}
		c       += normal[i] * point[i];
	}
	for (i = 0 ; i < 3; ++i)
		dt[i] = (max[i] - min[i] )/(double)(nP -1);
	fprintf(f,"%lf %lf %lf\n",point[0], point[1], point[2]);
	for (i = 0 ; i < nP; ++i)
	{
		for ( j = 0 ; j < nP; ++j)
		{
			for ( k = 0 ; k < nP; ++k)
			{
				p[0] = min[0] + (double)i*dt[0];
				p[1] = min[1] + (double)j*dt[1];
				aux  = c - p[0] * normal[0] - p[1] * normal[1];
				if ( fabs(normal[2]) > EPS11) p[2] = aux / normal[2];
				else {
					p[2] = min[2] + (double)k*dt[2];
					if (fabs(normal[1]) < EPS11) {
						p[1] = min[1] + (double)j*dt[1];
						p[0] = point[0];
					}
					else {
						if (fabs(normal[0]) < EPS11) p[1] = point[1];
						else {
							aux  = c - p[0] * normal[0];
							p[1] = aux / normal[1];
						}
					}
				}
				fprintf(f,"%lf %lf %lf\n",p[0], p[1], p[2]);
			}
		}
	}
	fclose(f);
}







