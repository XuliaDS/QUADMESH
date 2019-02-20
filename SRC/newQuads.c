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
#define DEG160 2.79252680319093
#define DEG10  0.1745329251994329
#define PI            3.1415926535897931159979635
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
	int      fID, oriQ, oriV, minQ, extraQuads, smooth, *vFix;
	ego      face;
	double   minsize, maxsize, range[4];
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
static int        nb[10];
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

static int
inList ( int n, int *list, int p) {
	int i ;
	for ( i = 0 ; i < n ; i++ )
		if ( list [i] == p ) return i;
	return -1;
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








/* IO FUNCTIONS */
//static void printStar               (vStar *star) ;
//static void printStarFile(vStar *star, meshMap *qm) ;
static void printVertexCoords (meshMap *qm, int v );
static void printMeshStats(meshMap *qm, int sweep);
static void printMesh     (meshMap *qm,char[], int  );
static int  checkMesh     (meshMap *qm) ;
void sampleNormalPlane(meshMap *qm, char *name, double *uv ) ;
static int EG_normalToSurface       (meshMap *qm, double *uv, double *normal )  ;
static int EG_projectToTangentPlane (double normal[], double *O, double *p, double *proj) ;
/***********************************************************/
/* MESH MAP FUNCTIONS*/
static int    EG_createMeshMap        (bodyData *bodydata);
static double EG_angleAtPlaneUsingUVs ( meshMap *qm, double *uvC, double *uv1, double *uv2);
static double EG_angleAtVnormalPlane  (meshMap *qm, int vC, int v1, int v2 );
static int    EG_angleAtBoundaryVertex(meshMap *qm, int v, int *links, double *size ) ;
static int    vertexLinksToBounds     (meshData *mesh, int vID, int *nb );
static int    quadAngleOrientation    (meshMap *qm, int qID,  int *ori, int *order, double *theta );
static int    quadAverageCoords       (meshMap *qm, int q, double *uv, double *p) ;
static void   averageCoordsMinusLinks (meshMap *qm,  int vc, int l1, int l2 ) ;
static void   weightedAverage         (meshMap *qm, int vC, double max, int *fix );
static void   updateVertex            (meshMap *qm, int vID, double *uv );
static int    optimize_angles         (meshMap *qm, int nP, /*@unused@*/ /*@null@*/int *pList, int fullRegularization);
static int    moveDistAway            (meshMap *qm, double alpha, int iA, int iB);
/* REGULARIZATION FUNCTIONS */
/* FUNDAMENTAL */
static int  EG_swappingOperation      (meshMap *qm, quadGroup qg, int swap);
static int  EG_splittingOperation     (meshMap *qm, int vC, int vL, int vR);
static int  EG_mergeVertices          (meshMap *qm, int qC, int centre);
static void  EG_reduceAngle            (meshMap *qm,  double dtheta, int iA, int iB, int iC, int reverse);
static double EG_segment              (meshMap *qm, double *uv1, double *uv2 );
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
rotationUV ( meshMap *qm, double theta, int counterClockWise, double *t0, double *t1) {
	int d, stat;
	double normal[4], proj[3], cross[4], xyz1[18], dt, duv[2], xyz0[18], prevt[2], vec[8];
	stat = EG_normalToSurface (qm, t0, normal);
	if ( stat != EGADS_SUCCESS ) {
		printf("EG_rotationUV :: EG_normalToSurface ->%d!!\n",  stat);
		return;
	}
	prevt[0]  = t1[0]; prevt[1]  = t1[1];
	stat      = EG_evaluate ( qm -> face, t0, xyz0);
	stat     += EG_evaluate ( qm -> face, t1, xyz1);
	if ( stat != EGADS_SUCCESS ){
		printf(" RotationUV rotation vector has EG_evaluate %d!!\n", stat);
		return;
	}
	stat   = EG_projectToTangentPlane(normal, xyz0, xyz1, proj);
	if ( stat != EGADS_SUCCESS ) {
		printf("EG_rotateUV :: EG_projectToTangentPlane %d !!\n", stat);
		return ;
	}
	vec[0]  = proj[0] - xyz0[0];
	vec[1]  = proj[1] - xyz0[1];
	vec[2]  = proj[2] - xyz0[2];
	unitVector (vec, &vec[3]);
	for ( d = 0 ; d < 2; d++ ) {
		dt = theta;
		if ( d == 1 ) dt *= -1.0;
		duv[0] = ( prevt[0] - t0[0] ) * cos ( dt ) - ( prevt[1] - t0[1] ) * sin ( dt );
		duv[1] = ( prevt[0] - t0[0] ) * sin ( dt ) + ( prevt[1] - t0[1] ) * cos ( dt );
		t1[0]  = t0[0] + duv[0];
		t1[1]  = t0[1] + duv[1];
		if (    t1[0] < qm ->range[0] || t1[0] > qm ->range[1] ||
				t1[1] < qm ->range[2] || t1[1] > qm ->range[3] ) {
			if ( d == 0 ) continue;
			printf(" Both directions seem invalid! BUG!\n");
			t1[0] = prevt[0]; t1[1] = prevt[1];
			return;
		}
		stat       = EG_evaluate ( qm -> face, t1, xyz1);
		if ( stat != EGADS_SUCCESS ) {
			if ( d == 0 ) continue;
			printf(" Both directions seem invalid! BUG!\n");
			t1[0] = prevt[0]; t1[1] = prevt[1];
			return;
		}
		stat = EG_projectToTangentPlane(normal, xyz0, xyz1, proj);
		if ( stat != EGADS_SUCCESS ) {
			printf("EG_rotateUV :: EG_projectToTangentPlane %d !!\n", stat);
			t1[0] = prevt[0]; t1[1] = prevt[1];
			return ;
		}
		vec[4]  = proj[0] - xyz0[0];
		vec[5]  = proj[1] - xyz0[1];
		vec[6]  = proj[2] - xyz0[2];
		unitVector (&vec[4], &vec[7] );
		cross_product (vec, &vec[4], cross);
		if ( dotProduct (normal, cross ) * (double)counterClockWise > 0 ) return;
		else {
			if ( d == 1 ) {
				printf(" both directions give negative dot product !!\n ");
				t1[0] = prevt[0]; t1[1] = prevt[1];
				return;
			}
		}
	}
	return;
}





/* clockwise = -1 , counterclockwise = 1
 * reverse = 0 assumes invalid quads e.g., ( AB x AC ) * normal < 0 so rotated B towards C (clockwise ) and
 * C towards B ( counterclockwise).if <ABC is small, it will increase the angle.
 * reverse = 1 assumes <ABC is big and will move B counter-clockwise and C clockwise
 */
static void
EG_reduceAngle ( meshMap *qm, double theta, int iA, int iB, int iC, int reverse) {
	int d, i, a, v,ii, counterClockWise = 1, links[2];
	double uvback[2], dist0, distm, ang0, ang1, seg0, seg1;
	if ( reverse == 1 ) counterClockWise *= -1;
	printf("\n\n REDUCE ANGLE %d %d %d reverse %d  \n ", iA, iB, iC, reverse);
	for ( a = 0 ; a < 2; a++) {
		counterClockWise *= -1;
		if ( a == 0 ) v = iB - 1;
		else          v = iC - 1;
		printf(" VERTEX %d direction %d  A TYPE %d  \n ", v + 1, counterClockWise, qm -> mesh -> vType [ iA -1 ] );
		if ( qm -> mesh -> vType [v] != -1 ) continue;
		if ( qm -> vFix[0] > 0 && inList (qm -> vFix[0], &qm -> vFix[1], v + 1) >= 0 ) continue;
		/*if ( qm -> mesh -> vType [ iA -1 ] == -1 && vertexLinksToBounds (qm -> mesh, v + 1, nb) > 0 ) {
			distm = 0.0;
			for ( d = ii = i = 0 ; i < nb[0]; i++ ) {
				if ( qm -> mesh -> vType [ nb[ 1 + i] -1 ] >=4 ) d++;
				dist0 = EG_segment (qm, &qm -> mesh -> uvs[ 2 * v ], &qm -> mesh -> uvs[ 2 * ( nb[1 + i] - 1)  ] );
				if ( dist0 > distm ) {
					ii    = i + 1;
					distm = dist0;
				}
			}
			if ( d >= 0 ) {
				moveDistAway (qm, 0.95, nb[ii], v + 1 );
				continue;
			}
		}*/
		if ( a == 0 )
			printf("# iB %d --> %f %f %f \n ", iB, qm -> mesh -> xyzs [ 3 * v],
					qm -> mesh -> xyzs [ 3 * v + 1], qm -> mesh -> xyzs [ 3 * v + 2] );
		else
			printf("# iC %d --> %f %f %f \n ", iC, qm -> mesh -> xyzs [ 3 * v],
					qm -> mesh -> xyzs [ 3 * v + 1], qm -> mesh -> xyzs [ 3 * v + 2] );
		rotationUV (qm, theta, counterClockWise, &qm -> mesh -> uvs[ 2 * ( iA -1 ) ],
				&qm -> mesh -> uvs[ 2 * v ]);
		updateVertex (qm, v + 1, &qm -> mesh -> uvs[ 2 * v ]);
		printVertexCoords (qm, v+1 );
		nb[0]   = vertexLinksToBounds (qm ->mesh, v + 1, nb );
		for ( i = 0 ; i < nb[0]; i++ ) {
			EG_angleAtBoundaryVertex (qm, nb[i + 1], links, &distm );
			ang0 = EG_angleAtPlaneUsingUVs(qm, &qm -> mesh -> uvs[ 2 * ( nb[i + 1] - 1 )] ,
					&qm -> mesh -> uvs[ 2 * ( links[0] - 1 )],
					&qm -> mesh -> uvs[ 2 * v]);
			printf("l0 ANGLE WITH RESPECT %d %d --> %f\n ", nb[i + 1], links[0], ang0);
			ang1 = EG_angleAtPlaneUsingUVs(qm, &qm -> mesh -> uvs[ 2 * ( nb[i + 1] - 1 )] ,
					&qm -> mesh -> uvs[ 2 * v],
					&qm -> mesh -> uvs[ 2 * ( links[1] - 1 )]);
			printf("l1 ANGLE WITH RESPECT %d %d --> %f\n ", nb[i + 1], links[1], ang1);
			if ( ang0 > distm ) {
				seg0 = EG_segment(qm, &qm -> mesh -> uvs[ 2 * ( nb[1 + i] - 1 ) ], &qm -> mesh -> uvs[ 2 * v ]);
				printf(" CORRECT ROTATION BY FLIPPING Lo\n");
				uvback[0] =  qm -> mesh -> uvs[ 2 * (links[1] -1 )    ];
				uvback[1] =  qm -> mesh -> uvs[ 2 * (links[1] -1 ) + 1];
				rotationUV (qm, theta, -1, &qm -> mesh -> uvs[ 2 * ( nb[ i + 1] -1 ) ],
						uvback);
				updateVertex (qm, v + 1, uvback);
				printVertexCoords (qm, v+1 );
				seg1 = EG_segment(qm, &qm -> mesh -> uvs[ 2 * ( nb[i + 1] -1 ) ], &qm -> mesh -> uvs[ 2 * v ]);
				printf(" SEGMENT 0 %f  segment 1 %f ratio %f \n ", seg0, seg1, seg0 / seg1);
				moveDistAway (qm, seg0 / seg1, nb[i + 1], v + 1 );
				printVertexCoords (qm, v+1 );
			}
			else if ( ang1 >= distm ) {
				seg0 = EG_segment(qm, &qm -> mesh -> uvs[ 2 * ( nb[1 + i] - 1 ) ], &qm -> mesh -> uvs[ 2 * v ]);
				printf(" CORRECT ROTATION BY FLIPPING L1\n");
				uvback[0] =  qm -> mesh -> uvs[ 2 * (links[0] -1 )    ];
				uvback[1] =  qm -> mesh -> uvs[ 2 * (links[0] -1 ) + 1];
				rotationUV (qm, theta, 1, &qm -> mesh -> uvs[ 2 * ( nb[ i + 1] -1 ) ],
						uvback);
				updateVertex (qm, v + 1, uvback);
				printVertexCoords (qm, v+1 );
				seg1 = EG_segment(qm, &qm -> mesh -> uvs[ 2 * ( nb[i + 1] -1 ) ], &qm -> mesh -> uvs[ 2 * v ]);
				printf(" SEGMENT 0 %f  segment 1 %f ratio %f \n ", seg0, seg1, seg0 / seg1);
				moveDistAway (qm, seg0 / seg1, nb[i + 1], v + 1 );
				printVertexCoords (qm, v+1 );
			}
			printVertexCoords (qm, v+1 );
		}
	}
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
	if ( q1 <= 0 || q2 <= 0 || q1 > mesh -> totQuads || q2 > mesh -> totQuads ) {
		fprintf(stderr, "ERROR %d %d !!!!!\n", q1, q2);
		exit(1);
	}
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
		j = mesh ->quadIdx [ 4 * ( q - 1) + i];
		if ( j == -2 ) return j;
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
	if ( qID < 0 || qID > mesh -> totQuads ) {
		fprintf(stderr,"VALID COLLAPSE AT %d !!!\n ", qID );
		return EGADS_INDEXERR;
	}
	id  = EG_quadVertIdx ( mesh, qID, v);
	aux = mesh ->quadIdx [ 4 * ( qID - 1 ) + (id + 2)%4];
	if ( mesh -> vType [ aux -1 ] != -1 && vertexLinksToBounds ( mesh, v, nb ) != 0 ) return 0;
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
	printf(" DOUBLE COLLAPSE %d %d FORCING %d f %d \n ", qg.q[0], qg.q[1], forcing, f );
	for ( i = 0 ; i < 2; i++ ) {
		printf(" DOUBLE COLLAPSE QUAD %d \n ", qg.q[( q3 + i ) %2] );
		if (checkQuad ( qm -> mesh, qg.q[( q3 + i ) %2] ) != EGADS_SUCCESS )continue;
		stat       = EG_forceCollapse (qm, qg.q[( q3 + i ) %2], &act);
		printf(" FORCE COLLAPSE %d is %d \n ",qg.q[( q3 + i ) %2], act);
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
		printf(" EG transfer valences :: clean quad %d forcing %d \n ",qID[0], try5533 );
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
		printf(" TRANSFER VALENCES::: QUAD GROUP \n ");
		printQuadGroup (qm -> mesh, qg );
		if ( stat != EGADS_SUCCESS  ) {
			printf(" Inside EG_transferValences EG_createQuadGroup %d !!\n", stat );
			printQuadGroup (qm -> mesh, qg );
			return stat;
		}
		if ( *transfering == 0 && qg.vals[0] * qg.vals[3] == 15 && qID[1] != -1  ) {
			printf(" COMING FROM %d %d ---> ", qID[0], qID[1] );
			j = 0;
			if ( qg.q[j] == qID[0] ) j = 1;
			qID[0] = qg.q[j]; qID[1] = qg.q[(j + 1)%2];
			printf(" NOW  %d %d ---> \n", qID[0], qID[1] );
			if ( EG_quadIsBoundary(qm -> mesh, qg.q[1] ) == 1 && qm -> extraQuads >= 0 ) {
				printf(" EXTRAQUADS %d ---> FORCE COLLAPSE TRANSFEr AT %d\n ", qm -> extraQuads, qID[0]);
				stat = EG_forceCollapse (qm, qID[0], &(*activity));
				printMesh (qm, buffer,0);
				qID[0] = qID[1];
				qID[1] = -1;
				if ( stat != EGADS_SUCCESS ) printf("EG_transferValences: forceCollapse gave %d !!\n", stat);
				return stat;
			}
			else if  ( EG_quadIsBoundary   (qm -> mesh, qg.q[1] ) != 1 )
				if ( qID[0] <= 0 || qID[0] > qm -> mesh -> totQuads ) {
					fprintf(stderr, " TRANSFER POR %d!!!\n ", qID[0] );
					exit(1);
				}
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
				if ( qAux[0] <= 0 || qAux[0] > qm -> mesh -> totQuads ) {
					fprintf(stderr, " TRANSFER POR %d!!!\n ", qID[0] );
					exit(1);
				}

				stat    = EG_transferValences ( qm, qAux, 0, &i, &min);
				if ( stat != EGADS_SUCCESS ) {
					printf(" EG_transferValences: separating valences after forcing %d!!\n", stat);
					return stat;
				}
			}
			*transfering = 0;
			qID[0] = qAux[0];
			qID[1] = qAux[1];
		}
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
			if ( checkQuad ( qm -> mesh, qID  ) != EGADS_SUCCESS ) return EGADS_SUCCESS;
		}
	}
	printf(" CLEAN QUAD %d --> focing %d adj %d\n ", qID, forcing, useAdj );
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
		printf(" OP TYPE %d \n ", opComp[i]);
		if ( transfer == 1 && opComp[i] == DOUBLESWAP ) continue;
		for ( q  = 0 ; q < 4; q++) {
			if ( checkQuad ( qm -> mesh, qID  )   != EGADS_SUCCESS ) return EGADS_SUCCESS;
			qadj = qm -> mesh -> quadAdj [ 4 * ( qID -1 ) + q];
			if ( qadj == -1 ) continue;
			stat = EG_createQuadGroup (qm -> mesh, &qg, qID, qadj);
			printf(" CALLING COMPOSE OPERATION USING %d %d forcing %d\n ", qg.q[0], qg.q[1], forcing);
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
	int         f, stat = 0, stat2 = 0, j, q, i, auxID, k, kk, kOK, len,  ntri, nquad, e4[4];
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
		bodydata -> qm[f] -> smooth     = 0;
		// Edges associated to face //
		stat = EG_getTessFace(bodydata->tess, f + 1, &len,
				&xyzs, &uvs, &ptype, &pindex, &ntri,
				&tris, &tric);
		if ( stat != EGADS_SUCCESS ) {
			printf("EG_createMeshMap :: EG_getTessFace %d !!\n", stat );
			stat2 = stat;
		}
		for ( i = 0 ; i < bodydata -> nedges; i++ ) faceEdges[i] = 0;
		for ( j = i = 0 ; i < len; i++ ) {
			if ( pindex[i] == -1 ) continue;
			j++;
			faceEdges[pindex[i] -1]++;
		}
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
		bodydata -> qm[f] -> vFix    = EG_alloc( 2 * len * sizeof ( int ));
		if (bodydata->qm[f] -> vFix == NULL ) return EGADS_MALLOC;
		bodydata -> qm[f] -> vFix[0] = 0 ;
		if ( stat != EGADS_SUCCESS ||
				bodydata ->qm[f] -> mesh       == NULL ||
				bodydata ->qm[f] -> backupMesh == NULL ||
				bodydata ->qm[f] -> bestMesh   == NULL ) {
			fprintf(stderr,"In createMeshMap  EG_allocMeshData = %d\n ", stat );
			stat2 = stat;
			continue;
		}
		bodydata->qm[f] -> extraQuads       = 0;
		bodydata->qm[f] -> face             = bodydata->faces[f];
		stat                                = EG_getRange ( bodydata -> qm[f] -> face, bodydata -> qm[f] -> range, &i);
		if ( stat != EGADS_SUCCESS ) {
			printf(" Face %d EG_getRange %d !!\n ", f + 1, stat );
			stat2  = stat ;
			continue;
		}
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
				stat2 = EGADS_GEOMERR;
				printf(" Stat in EG_angleAtBOundaryVertex %d angle %f\n ", stat, angle );
				break;
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
		printf(" Mesh min max sizes are %f  %f  \n ",  bodydata->qm[f] -> minsize , bodydata->qm[f] -> maxsize );
		stat = checkMesh (bodydata->qm[f]);
		if ( stat != EGADS_SUCCESS ) {
			printf("In EG_createMeshMap :: checkMesh at face %d --> %d!!\n", f + 1, stat );
			stat2 = stat;
		}
	}
	EG_free(faceEdges);
	return stat2;
}


static void EG_destroymeshMap(meshMap ***qm, int nfaces ) {
	int i ;
	if ( *qm  == NULL ) return ;
	printf(" total faces %d \n ", nfaces );
	for ( i = 0 ; i < nfaces; ++i) {
		if ( ( *qm)[i] ) {
			EG_freeMeshData (&((*qm)[i] -> mesh      ));
			EG_freeMeshData (&((*qm)[i] -> backupMesh));
			EG_freeMeshData (&((*qm)[i] -> bestMesh  ));
			EG_free ( (*qm)[i] -> vFix );
			EG_free(( *qm)[i]);
		}
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
	printf(" #=================== quad %d ==================\n", qID );
	for ( i = 0 ; i <= 4; i++ ) {
		v = qm -> mesh -> quadIdx [ 4 * ( qID - 1) +  i%4] - 1;
		printf("%lf %lf %lf %lf %lf %d\n",
				qm -> mesh -> xyzs [3 * v ],qm -> mesh -> xyzs [3 * v +1], qm -> mesh -> xyzs [3 * v +2 ],
				qm -> mesh -> uvs [2 * v + 1],  qm -> mesh -> uvs [2 * v ], v + 1);
	}
	printf(" #==============================================\n" );
}


static int EG_normalToSurface (meshMap *qm, double *uv, double *normal ) {
	int   stat;
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

static void weightedAverage (meshMap *qm, int vID, double maxAngle, int *fix ) {
	int i, j, auxID = 0, stat, links[2], nI, i0, l1, l2, it, itMAX = 100;
	double uvc[2], angle, optangle, angL0, ang1, ang2, xyz0[3], seg, dt, dist ;
	vStar *star = NULL;
	*fix = 0;
	if ( qm -> mesh -> vType[ vID -1] != -1 ) return;
	printVertexCoords (qm, vID );
	printf("========== weighted Average coordinates %d fix %d smooth %d \n ", vID, j, qm -> smooth ) ;
	uvc[0]    = qm -> mesh -> uvs [ 2 * ( vID - 1 )     ];
	uvc[1]    = qm -> mesh -> uvs [ 2 * ( vID - 1 ) + 1 ];
	j      = 1;
	for (i = 0; i < qm -> mesh -> valence [ vID - 1][0]; i++) {
		auxID   = qm -> mesh -> valence [ vID - 1][2 + i] - 1;
		uvc[0] += qm -> mesh -> uvs     [2 * auxID    ];
		uvc[1] += qm -> mesh -> uvs     [2 * auxID + 1];
		j++;
	}
	uvc[0] /= (double)j;
	uvc[1] /= (double)j;
	printf(" centred coordinates using point as weight too \n ");
	updateVertex (qm, vID, uvc);
	printVertexCoords (qm, vID );
	if ( qm -> smooth > 0  ) {
		if ( maxAngle < DEG160 ) maxAngle = DEG160;
		if ( vertexLinksToBounds (qm -> mesh, vID, nb ) == 0 ) return;
		else if ( nb[0] >= 2 ) {
			stat      = EG_buildStar ( qm -> mesh, &star, vID);
			if ( star == NULL || stat != EGADS_SUCCESS) return;
			for  ( j = i0 = i = 0; i < star -> nQ; i++ ) {
				auxID = star -> verts [ 2* i + 1] - 1;
				if ( qm -> mesh -> vType [ auxID] >= 0 ) {
					if ( i0 == 0 ) i0 = 2 * i + 1;
					else {
						j = 2 * i + 1;
						if ( j - i0 > 2 && ( i0 + ( star -> nV -1 ) - j != 2 ) ) {
							i0 = star -> nV;
							break;
						}
					}
				}
			}
			if ( i0 != star -> nV ) {
				EG_freeStar (&star);
				return;
			}
			*fix = 1;
			printf(" FIXING VERTEX TO BOUNDS !!!!\n ");
			uvc[0]   = 0.0;
			uvc[1]   = 0.0;
			for  ( j = i = 0; i < star -> nQ; i++ ) {
				auxID = star -> verts [ 2 * i + 1] - 1;
				if ( qm -> mesh -> vType [ auxID] >= 0 ) {
					uvc[0] += qm -> mesh -> uvs     [2 * auxID    ];
					uvc[1] += qm -> mesh -> uvs     [2 * auxID + 1];
					j++;
				}
			}
			uvc[0] /= (double)j;
			uvc[1] /= (double)j;
			updateVertex (qm, vID, uvc);
			EG_freeStar (&star);
			return;
		} else {
			auxID = nb[1];
			stat      = EG_buildStar ( qm -> mesh, &star, auxID );
			if ( star == NULL || stat != EGADS_SUCCESS) return;
			stat      = EG_angleAtBoundaryVertex (qm, auxID, links, &angle);
			for  ( i0 = nI = i = 0; i < star -> nQ; i++ ) {
				if ( star -> quads [i] == -1 ) i0 = i + 1;
				else if (qm ->mesh -> vType[ star -> verts [ 2 * i + 1 ] - 1] == -1 ) nI++;
			}
			for ( i = 0 ; i < star -> nQ; i++ )
				if (  star -> verts [ star -> idxV [ 2 * ( i0 + i + 1) + 1 ]] == vID ) break;
			optangle    = angle * ( i + 1 ) / (double)( nI + 1 );
			l1          = star -> verts [ star -> idxV [ 2 * ( i0 + i     ) + 1 ]];
			l2          = star -> verts [ star -> idxV [ 2 * ( i0 + i + 2 ) + 1 ]];
			printf("j = %d  BOUND LINK %d --> VERTEX NEIGHBORS %d %d\n", j, auxID, l1, l2 );
			angL0    = EG_angleAtVnormalPlane ( qm, auxID, links[0], vID);
			ang1     = EG_angleAtVnormalPlane ( qm, auxID, l1, vID);
			ang2     = EG_angleAtVnormalPlane ( qm, auxID, vID, l2);
			printf(" Angles %f  %f o opt %f ( OPTI %f ) MAXANGLE %f\n", ang1, ang2, angL0, optangle, maxAngle);
			it  = 0;
			dt  = 0.01;
			seg = EG_segment ( qm, &qm -> mesh -> uvs[2 * (auxID - 1) ], &qm -> mesh -> uvs [2 * (vID - 1) ] );
			if ( seg < 0.25 * qm -> minsize )
				stat = moveDistAway (qm, 1.0, auxID, vID );
			xyz0[0] = qm -> mesh -> xyzs [ 3 * ( vID -1 )     ];
			xyz0[1] = qm -> mesh -> xyzs [ 3 * ( vID -1 ) + 1 ];
			xyz0[2] = qm -> mesh -> xyzs [ 3 * ( vID -1 ) + 2 ];
			while  ( ang1 > maxAngle || ang2 > maxAngle  ||
					fabs ( angL0 - optangle ) > maxAngle      ) {
				*fix = 1;
				if ( ang1 > optangle ) EG_reduceAngle(qm, dt, auxID, links[0], vID, 1);
				else EG_reduceAngle(qm, dt, auxID, links[0], vID, 0);
				printf(" WE ARE USING CENTRE %d PRV NXt %d %d links %d %d\n ",
						auxID, l1, l2, links[0], links[1]);
				angL0     = EG_angleAtVnormalPlane ( qm, auxID, links[0], vID);
				ang1      = EG_angleAtVnormalPlane ( qm, auxID, l1, vID);
				ang2      = EG_angleAtVnormalPlane ( qm, auxID, vID, l2);
				++it;
				printf(" FIRST ANGLE %f > %f \n SECOND ANGLE %f > %f\n DISTANCE TO OPTI %f %f\n ",
						ang1, 0.8 * PI, ang2, 0.8 * PI, fabs (angL0 - optangle ), 0.8 * PI );
				printVertexCoords (qm, vID );
				dist  = ( xyz0[0] - qm -> mesh -> xyzs [ 3 * ( vID -1 )   ] ) * ( xyz0[0] - qm -> mesh -> xyzs [ 3 * ( vID -1 )    ] );
				dist += ( xyz0[1] - qm -> mesh -> xyzs [ 3 * ( vID -1 ) +1] ) * ( xyz0[1] - qm -> mesh -> xyzs [ 3 * ( vID -1 ) + 1] );
				dist += ( xyz0[2] - qm -> mesh -> xyzs [ 3 * ( vID -1 ) +2] ) * ( xyz0[2] - qm -> mesh -> xyzs [ 3 * ( vID -1 ) + 2] );
				if ( sqrt ( dist ) < 0.1 * qm -> minsize ) dt += 0.01;
				printf( " it %d ----> dt %f distance %f  vs %f \n ", it, dt, sqrt (dist), 0.1 * qm -> minsize );
				if ( it > itMAX || (
						ang1 < PI && ang2 < PI && fabs ( angL0 - optangle ) < PI && sqrt ( dist ) >= 0.1 * qm -> minsize ))
					break;
			}
		}
		EG_freeStar ( &star );
	}
	printf(" LEAVE WITH COORDINATES\n");
	printVertexCoords (qm, vID );
	return ;
}





static void
averageCoordsMinusLinks ( meshMap *qm, int vc, int l1, int l2 ) {
	int i, j, k, n, *links = NULL ;
	double uv[2];
	vc--;
	printf(" AVERAGE COORDS MINUS LINKS %d %d %d\n ", vc + 1, l1, l2 );
	if ( qm ->mesh -> vType[ vc ] != -1 ) return;
	if ( vertexLinksToBounds (qm -> mesh, vc + 1, nb) > 0 &&
			qm ->mesh -> vType[ nb[1] - 1 ] >= 4 && ( nb[1] != l1 && nb[1] != l2 )) {
		moveDistAway (qm, 0.5, nb[1], vc + 1) ;
	}
	else {
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
		updateVertex (qm, vc + 1, uv );
		EG_free(links);
	}
}


static int vertexLinksToBounds (meshData *mesh, int vID, int *nb ) {
	int i, j;
	for ( nb[0] = i = 0 ; i < mesh -> valence[ vID -1 ][0]; i++ ) {
		j  =  mesh -> valence[ vID -1 ][2 + i] - 1;
		if ( mesh -> vType[ j ] >= 4 ) nb[ ++nb[0] ] = j + 1;
	}
	for ( i = 0 ; i < mesh -> valence[ vID -1 ][0]; i++ ) {
		j  =  mesh -> valence[ vID -1 ][2 + i] - 1;
		if ( mesh -> vType[ j ] >= 0 && mesh -> vType [ j ] < 4 ) nb[ ++nb[0] ] = j + 1;
	}
	return nb[0];
}


static int makePositiveAngles (meshMap *qm, int qID, double minAngle, double maxAngle ) {
	int reverse, i, j, k, kk, iO, iA, iB, iC, tO, tA, tB, fullNegative = 0, it , itMAX = 20, ori[4], piv[4], boundary = -1, area = 0, areaOpt = 0, bl[2];
	static int pivTri[6] = { 1, 2, 2, 3, 1, 3 };
	double ang, angles[4], dt = 0.01, ds = 0.95, xyzab[6], dista, distb;
	printMesh (qm, buffer, 0);
	for ( k = j = 0; j < 4; j++ ) {
		iO      = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + j] - 1;
		if ( qm -> mesh -> vType [ iO ] >= 0 &&
				( boundary == -1 || ( boundary >= 0 && qm -> mesh -> vType [ k ] < 4 ) )) {
			boundary = j;
			k = iO;
		}
	}
	bl[0] = 0; bl[1] = 0 ;
	if ( maxAngle < DEG160 ) maxAngle = DEG160;
	if ( minAngle > DEG10  ) minAngle = DEG10;
	printf(" Inside make Positive function for %d \n ", qID );
	for ( it = 0 ; it < itMAX; it++ ) {
		printQuadCoords (qm, qID);
		area = quadAngleOrientation (qm, qID, ori, piv, angles);
		if ( area < 0 ) {
			printf("makePositiveAngles: area in %d --> %d!!\n ", qID, area);
			printMesh(qm, buffer,0);
			break;
		}
		if ( area == 1 ) {
			for ( j = 0 ; j < 4; j++ )
				if (angles [j] > maxAngle || angles[j] < minAngle ) {
					area    = 0 ;
					areaOpt = 1;
				}
			if ( area == 0 )
				printf(" makePositive %d --> min max angle %f %f \n ", area, minAngle, maxAngle);
		}
		if ( area == 1 || ( area == 2 && maxAngle > PI )) break;
		fullNegative = ori[0] + ori[1] + ori[2] + ori[3];
		if ( boundary != -1 )
			for ( j = 0; j < 4; j++ ) piv[j] = ( boundary + j ) %4 ;
		if ( it %5 == 0 ) printMesh(qm, buffer, 0);
		for ( j   = 0; j < 4; j++ ) {
			k     = piv[j];
			iO    = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + k];
			tO    = qm -> mesh -> vType [ iO - 1] ;
			i     = vertexLinksToBounds (qm -> mesh, iO, nb );
			if ( area == 0 && ori[k] == -1 && qm -> mesh -> vType [ iO - 1 ] == -1 && fullNegative != -4 &&
					( i == 0 || ( i > 0 && qm -> mesh -> vType [ nb[1] -1 ] < 4 ))) {
				printf(" SKIP VERTEX \n ");
				continue;
			}
			printf ("=========== O ( %d ) = %d ::: MAX ANGLE PASS %f %f   ==============================\n",k, iO,  minAngle, maxAngle );
			for (i  = 0; i < 3; i++ ) {
				iA  = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + ( k + pivTri[2 * i]    ) %4];
				iB  = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + ( k + pivTri[2 * i + 1]) %4];
				tA  = qm -> mesh -> vType [ iA - 1] ;
				tB  = qm -> mesh -> vType [ iB - 1] ;
				ang = EG_angleAtVnormalPlane( qm, iO, iA, iB);
				if ( ang >= maxAngle ) {
					printf("\n\n i = %d --> Group %d %d %d --> ANGLE %f \n",i, iO, iA, iB , ang );
					printQuadCoords (qm, qID );
					if ( tO >= 0 || area == 0 ) {
						dista = EG_segment (qm , &qm -> mesh -> uvs [ 2 * ( iO - 1) ],
								&qm -> mesh -> uvs [ 2 * ( iA - 1) ]);
						distb = EG_segment (qm , &qm -> mesh -> uvs [ 2 * ( iO - 1) ],
								&qm -> mesh -> uvs [ 2 * ( iB - 1) ]);
						printf(" DIST %f  %f  RATIO %f  %f \n ", dista, distb, distb / dista, dista / distb);
						i       = vertexLinksToBounds (qm -> mesh, iA, nb );
						if ( i == 0 || ( i == 1 && qm -> mesh -> vType [ nb[1] -1] < 4 ) ) i = 1;
						if ( i == 1 && dista > distb && tA == -1 && ( distb / dista) < 0.95 )
							moveDistAway (qm, ds, iO, iA);
						i       = vertexLinksToBounds (qm -> mesh, iB, nb );
						if ( i == 0 || ( i == 1 && qm -> mesh -> vType [ nb[1] -1] < 4 ) ) i = 1;
						if ( i == 1 && dista < distb && tB == -1 && ( dista / distb) < 0.95 )
							moveDistAway (qm, ds, iO, iB);
					}
					for ( kk = 0; kk < 3; kk ++ ) {
						xyzab[    kk] = qm -> mesh -> xyzs [ 3 * ( iA -1 ) + kk];
						xyzab[3 + kk] = qm -> mesh -> xyzs [ 3 * ( iB -1 ) + kk];
					}
					reverse = 0;
					if ( tO >= 0 ) {
						i     = EG_angleAtBoundaryVertex ( qm, iO, bl, &ang );
						dista = EG_angleAtVnormalPlane (qm, iO, bl[0], iA );
						distb = EG_angleAtVnormalPlane (qm, iO, bl[0], iB );
						printf(" ORIGIN IS BOUNDARY VERTEX ::  ANGLE A %d %d %d %f \n", iO, bl[0], iA, dista);
						printf(" ORIGIN IS BOUNDARY VERTEX ::  ANGLE B %d %d %d %f \n", iO, bl[0], iB, distb);
						if  ( dista < distb || bl[0] == iB ) reverse = 1;
					}
					else if ( ang < 1.5 * PI ) reverse = 1;
					EG_reduceAngle (qm, dt, iO, iA, iB , reverse);
					ang   = EG_angleAtVnormalPlane  (qm, iO, iA, iB);
					dista = 0.0;
					distb = 0.0;
					for ( kk = 0; kk < 3; kk ++ ) {
						dista +=( xyzab[    kk] - qm -> mesh -> xyzs [ 3 * ( iA -1 ) + kk ] ) *
								( xyzab[    kk] - qm -> mesh -> xyzs [ 3 * ( iA -1 ) + kk ] ) ;
						distb +=( xyzab[3 + kk] - qm -> mesh -> xyzs [ 3 * ( iB -1 ) + kk ] ) *
								( xyzab[3 + kk] - qm -> mesh -> xyzs [ 3 * ( iB -1 ) + kk ] ) ;
					}
					if ( sqrt ( dista ) < 0.1 * qm -> minsize && sqrt ( distb ) < 0.1 * qm -> minsize) {
						dt += 0.01;
						ds -= 0.01;
					}
					printf( " it %d ----> dt %f ds %f distances %f  %f  : Angle i %d -> %d %d %d is now %f  \n ", it, dt, ds, sqrt ( dista ) , sqrt ( distb ),
							i, iO, iA, iB , ang );
					if ( ang >= maxAngle ) break;
				}
			}
			if ( ang >= maxAngle ) break;
			iA  = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + ( k + 1) %4];
			iB  = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + ( k + 3) %4];
			ang = EG_angleAtVnormalPlane( qm, iO, iA, iB);
			for ( kk = 0; kk < 3; kk ++ ) {
				xyzab[    kk] = qm -> mesh -> xyzs [ 3 * ( iA -1 )];
				xyzab[3 + kk] = qm -> mesh -> xyzs [ 3 * ( iB -1 )];
			}
			if (areaOpt == 1 && ang > maxAngle ) {
				printf("OPTIANGLE ang %f maxangle %f   reduce %d %d  %d\n ", ang, maxAngle, iO, iA, iB );
				EG_reduceAngle(qm, 0.01, iO, iA, iB , 1);
				printf("\n\n AFTER i = %d --> angle centre,%d  %d %d --> %f \n ",i, iO, iA, iB , ang );
			} else if (areaOpt == 1 &&  ang < minAngle ) {
				printf("OPTIANGLE  ang %f  mainagle %f  reduce %d %d  %d\n ", ang, minAngle, iO, iA, iB );
				EG_reduceAngle(qm, 0.01, iO, iA, iB , 0);
				printf("\n\n AFTER i = %d --> angle centre,%d  %d %d --> %f \n ",i, iO, iA, iB , ang );
			} else continue;
			dista = 0.0;
			distb = 0.0;
			for ( kk = 0; kk < 3; kk ++ ) {
				dista += ( xyzab[    kk] - qm -> mesh -> xyzs [ 3 * ( iA -1 ) + kk ] ) *
						(  xyzab[    kk] - qm -> mesh -> xyzs [ 3 * ( iA -1 ) + kk ] ) ;
				distb += ( xyzab[3 + kk] - qm -> mesh -> xyzs [ 3 * ( iB -1 ) + kk ] ) *
						(  xyzab[3 + kk] - qm -> mesh -> xyzs [ 3 * ( iB -1 ) + kk ] ) ;
			}
			if ( sqrt ( dista ) < 0.1 * qm -> minsize && sqrt ( distb ) < 0.1 * qm -> minsize)
				dt += 0.01;
		}
		if ( boundary != -1  ) {
			iO  = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) +  piv[0]];
			iC  = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + (piv[0] + 2)%4];
			iA  = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + (piv[0] + 3)%4];
			iB  = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + (piv[0] + 5)%4];
			ang = EG_angleAtVnormalPlane( qm, iC, iA, iB);
			if ( ang < minAngle )
				moveDistAway (qm, ds, iO, iC );
			printf(" BOUNDARY QUAD:::: SHRINK AND STRETCH \n ");
			if (areaOpt == 0 ){
				for ( j = 0; j < 4; j++ ) {
					k   = piv[j];
					iO  = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + k];
					if ( qm -> mesh -> vType [ iO -1 ] != -1 ) continue;
					i   = vertexLinksToBounds (qm -> mesh, iO, nb );
					if ( i > 0 )
						moveDistAway (qm, 0.95, nb[1], iO);
				}
			}
		}
	}
	area = quadAngleOrientation (qm, qID, ori, piv, angles);
	printf("\n\n=============================================================\n");
	printf(" LEAVE MAKEPOSITIVE QUAD %d IS NOW AREA %d \n ", qID, area );
	printf("\n\n=============================================================\n");
	printMesh (qm, buffer, 0);
	return area;

}

static int validBoundary ( meshMap *qm, double minAngle, double maxAngle ) {
	int  jj, v, i, k, j, q, qq, stat = 0 , aux, i0 = -1, area, area2, orient[4], piv[4];
	double  angles[4],  uvback[2], qback[8];
	vStar *star = NULL, *qStar = NULL;
	printf("\n\n=============================================================\n");
	printf(" VALID BOUNDARY CHECK  MIN MAX ANGLE %f  %f  \n ", minAngle, maxAngle );
	printf("=============================================================\n");
	printMesh(qm, buffer,0);
	for ( k = 0 ; k < qm -> vFix[0]; k++ ) printf(" VERTEX %d is FIXED \n ", qm -> vFix [ 1 + k] );
	for ( v = 0 ; v < qm -> mesh -> totVerts; v++ ) {
		if ( qm -> mesh -> vType[v] < 0 ) continue;
		printf(" VALID BOUNDARY FOR VERTEX %d SMOOTHING %d \n ", v + 1, qm -> smooth );
		stat      = EG_buildStar (qm ->mesh, &star, v + 1);
		if (stat != EGADS_SUCCESS || star == NULL ) goto cleanup;
		for  ( q = 0; q < star -> nQ; q++ ) {
			if ( star -> quads [q] == -1 ) i0 = q + 1;
		}
		for ( qq  = 0; qq < star -> nQ; qq++ )
		{
			q     = star -> idxQ [ qq + i0 ] ;
			if ( star -> quads[q] == -1 ) break;
			printf("VERTEX %d q = %d  Quad %d check orientation \n ", v + 1, q, star -> quads[q] );
			area = quadAngleOrientation(qm, star -> quads[q], orient, piv, angles);
			if ( area < 0 ) {
				printf("Inside validBoundary: quadANgleOrientation  %d = %d!!\n", star -> quads[q], stat);
				goto cleanup;
			}
			printMesh (qm, buffer, 0);
			if ( qm -> smooth > 0 && area >= 1 )  {
				printf(" Laplacian on boundary quad %d  \n", star -> quads[q] );
				for ( i = 0 ; i < 3; i++ ) {
					aux = star -> verts [ star -> idxV [ 2 * q + i + 1] ] -1 ;
					if ( qm -> mesh -> vType [aux] != -1 ) continue;
					if ( qm -> vFix[0] > 0 && inList (qm -> vFix[0], &qm -> vFix[1], aux + 1) >= 0 ) continue;
					printf(" WEIGHT VERTEX %d \n ", aux + 1 );
					uvback [0] = qm -> mesh -> uvs [ 2 * aux    ];
					uvback [1] = qm -> mesh -> uvs [ 2 * aux + 1];
					weightedAverage (qm, aux + 1, maxAngle, &k);
					area2 = quadAngleOrientation(qm, star -> quads[q], orient, piv, angles);
					if ( area2 > area ) {
						updateVertex (qm, aux + 1, uvback);
						qm -> vFix[++qm -> vFix[0]] = aux + 1;
						continue;
					}
					printMesh (qm, buffer, 0);
					if ( k == 1 )  qm -> vFix[++qm -> vFix[0]] = aux + 1;
					stat = EG_buildStar (qm -> mesh, &qStar, aux + 1 ) ;
					if ( stat != EGADS_SUCCESS || qStar == NULL ) {
						printf(" validBoundary In Laplacian boundary :: EG_buildStar at %d --> %d!!\n", aux + 1, stat );
						goto cleanup;
					}
					printf(" CHEcK THAT ALL GOOD !!!! ");
					for ( k = 0 ; k < qStar -> nQ; k++ ) {
						area = quadAngleOrientation(qm, qStar -> quads[k], orient, piv, angles);
						if ( area == 0  ) {
							for ( j  = 0 ; j < 3; j++ ){
								jj = qStar -> verts[ qStar -> idxV  [2 * k + j+ 1] ] -1 ;
								qback[2 * j    ] = qm -> mesh -> uvs [ 2 * jj    ];
								qback[2 * j + 1] = qm -> mesh -> uvs [ 2 * jj + 1];
							}
							printf(" VALID BOUNDARY ::::::  %d make positive \n", qStar -> quads[k] );
							area = makePositiveAngles(qm, qStar -> quads[k], 0.0, PI);//minAngle, maxAngle);
							if ( area == 0 ) {
								updateVertex (qm, aux + 1, uvback);
								updateVertex (qm, qStar -> verts[ qStar -> idxV [2 * k + 1] ], &qback[0]);
								updateVertex (qm, qStar -> verts[ qStar -> idxV [2 * k + 2] ], &qback[2]);
								updateVertex (qm, qStar -> verts[ qStar -> idxV [2 * k + 3] ], &qback[4]);
								printf(" Bad idea: Reset vertex %d  and make it fix \n", aux + 1);
								qm -> vFix[++qm -> vFix[0]] = aux + 1;
								break;
							}
						}
					}
					printf(" DONE EVERYTHING FOR VERTEX %d \n", aux  + 1);
					EG_freeStar ( &qStar );
				}
			}
			area = makePositiveAngles(qm, star -> quads[q], minAngle, maxAngle);
			if ( area <= 0 ) {
				fprintf(stderr, " Negative area at boundary quads %d. PROBLEM !!!\n ", star -> quads[q]);
				printf ( " Negative area at boundary quads %d. PROBLEM !!!\n ", star -> quads[q]);
				printMesh (qm, buffer, 0);
				if ( area < 0 ) return area;
			}
		}
	}
	if ( qm -> smooth > 0 ) {
		for ( q = 0 ; q < qm -> mesh -> totQuads; q++ ) {
			if (checkQuad ( qm -> mesh, q + 1)          != EGADS_SUCCESS ) continue;
			printf(" After boundary: Laplacian interior quad %d \n ", q + 1);
			area = quadAngleOrientation(qm, q + 1, orient, piv, angles);
			if (area == 1 ) {
				for ( i = 0 ; i < 4; i++ ) {
					j = qm -> mesh -> quadIdx [ 4 * q + i];
					if ( qm -> vFix[0] > 0 && inList (qm -> vFix[0], &qm -> vFix[1], j) >= 0 ) continue;
					if ( vertexLinksToBounds ( qm -> mesh, j, nb ) == 0 ) weightedAverage (qm, j, maxAngle, &k);
				}
			} else {
				area = makePositiveAngles(qm, q + 1, minAngle, maxAngle );
				if ( area == 1 ) {
					for ( i = 0 ; i < 4; i++ ) {
						j = qm -> mesh -> quadIdx [ 4 * q + i];
						if ( qm -> mesh -> vType [ j - 1] != -1 ) continue;
						if ( qm -> vFix[0] == 0 ||
								(qm -> vFix[0] > 0 && inList (qm -> vFix[0], &qm -> vFix[1], j ) == -1 ))
							qm -> vFix[++qm -> vFix[0]] = j;
					}
				}
				printQuadCoords (qm, q+1 );
			}
		}
	}
	cleanup:
	EG_freeStar(&star);
	return stat;
}






static int
EG_projectToTangentPlane(double normal[], double *O, double *p, double *proj) {
	double c, dotNN = 0.0, dotNP = 0.0, dist, lambda;
	dist    = (p[0]- O[0]) * (p[0]- O[0]) + (p[1]- O[1])*(p[1]- O[1]) + (p[2]- O[2])*(p[2]- O[2]);
	dist    = sqrt(dist);
	//printf(" PROJECT %f  %f  %f onto PLANE MADE BY %f %f %f normal %f %f %f\n ",
	//	p[0], p[1], p[2], O[0], O[1], O[2], normal[0], normal[1], normal[2]);
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
	//printf(" PROJCETION %f %f  %f\n ", proj[0], proj[1], proj[2]);
	dist  = normal[0] * proj[0] + normal[1] * proj[1] +  normal[2] * proj[2];
	//printf(" Distance point to plane %f \n ", dist - c);
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
EG_angleAtPlaneUsingUVs ( meshMap *qm, double *uvC, double *uv1, double *uv2) {
	int i, stat;
	double dot1, dot2, evalC[18], normal[4], cross[3], proj1[3], proj2[3], v01[4], v02[4], eval1[18], eval2[18];
	stat       = EG_evaluate(qm -> face, uvC, evalC);
	stat      += EG_evaluate(qm -> face, uv1, eval1);
	stat      += EG_evaluate(qm -> face, uv2, eval2);
	if ( stat != EGADS_SUCCESS) return stat;
	if ( qm -> face -> mtype == SREVERSE )
		cross_product(&evalC[6], &evalC[3], normal);
	else
		cross_product(&evalC[3], &evalC[6], normal);
	unitVector(normal, &normal[3]);
	//printf(" %f  %f  %f O\n ", evalC[0], evalC[1], evalC[2] );
	//printf(" %f  %f  %f A\n ", eval1[0], eval1[1], eval1[2] );
	//printf(" %f  %f  %f B\n ", eval2[0], eval2[1], eval2[2] );
	stat = EG_projectToTangentPlane(normal, evalC, eval1, proj1);
	stat = EG_projectToTangentPlane(normal, evalC, eval2, proj2);
	//printf(" %f  %f  %f PA\n ", proj1[0], proj1[1], proj1[2] );
	//printf(" %f  %f  %f PB\n ", proj2[0], proj2[1], proj2[2] );
	//printf(" NORMAL TO SURFACE %f  %f  %f \n ", normal[0], normal[1], normal[2] );
	for ( i = 0 ; i < 3 ; i++ ) {
		v01[i] = proj1[i] - evalC[i];
		v02[i] = proj2[i] - evalC[i];
	}
	//printf(" V1 %f %f %f  \n ", v01[0], v01[1], v01[2] );
	//printf(" V1 %f %f %f  \n ", v02[0], v02[1], v02[2] );
	unitVector (v01, &v01[3] );
	unitVector (v02, &v02[3] );
	//printf(" V1 %f %f %f  \n ", v01[0], v01[1], v01[2] );
	//printf(" V1 %f %f %f  \n ", v02[0], v02[1], v02[2] );
	cross_product (v01, v02, cross );
	dot1 = dotProduct ( v01, v02 );
	dot2 = dotProduct ( normal, cross);
	//printf(" DOTPRODUCTS ANGLE %f ORIENTATION %f\n ", dot1, dot2);
	if      ( dot1 >=  1.0 ) return 0.0;
	else if ( dot1 <= -1.0 ) return PI;
	if ( dot2 > 0 ) return  acos ( dot1 );
	else  return ( 2.0*PI - acos ( dot1 ));
}


static double
EG_angleAtVnormalPlane ( meshMap *qm, int vC, int v1, int v2 ) {
	return EG_angleAtPlaneUsingUVs(qm, &qm -> mesh -> uvs[ 2 * ( vC - 1 )] , &qm -> mesh -> uvs[ 2 * ( v1 - 1 )], &qm -> mesh -> uvs[ 2 * ( v2 - 1 ) ]);
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

static int moveDistAway ( meshMap *qm, double alpha, int iA, int iB) {
	int    i, stat, it;
	double seg, uvm[2], uvcopy[4], norm = 0.0, uv0[2], uv1[2], eval[18], size, dt, segprev;
	if ( qm -> mesh -> vType[iB -1] != -1 ||
			(qm -> vFix[0] > 0 && inList (qm -> vFix[0], &qm -> vFix[1], iB) >= 0 )) return EGADS_SUCCESS;
	uv0[0]     = qm -> mesh -> uvs[ 2 * ( iA - 1 )    ];
	uv0[1]     = qm -> mesh -> uvs[ 2 * ( iA - 1 ) + 1];
	uv1[0]     = qm -> mesh -> uvs[ 2 * ( iB - 1 )    ];
	uv1[1]     = qm -> mesh -> uvs[ 2 * ( iB - 1 ) + 1];
	seg        = EG_segment ( qm, uv0, uv1 );
	size = seg * alpha;
	if ( alpha * seg  < qm -> minsize * 0.25) size = 0.25 *  qm -> minsize;
	printf(" HOMOTECIA %d -> %d  SEG %f target %f \n", iA, iB, seg, size);
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
		seg    = EG_segment ( qm, uv0, uvm );
		if ( it > 100) {
			printf(" Seg < size Stuck. Break !\n");
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
		seg    = EG_segment ( qm, uv0, uvm );
		i      = 0;
		if ( seg > size ) i = 1;
		uvcopy[ 2 * i    ] = uvm[0];
		uvcopy[ 2 * i + 1] = uvm[1];
		it++;
		norm  = ( uvcopy[0] - uvcopy[2] )* ( uvcopy[0] - uvcopy[2] );
		norm += ( uvcopy[1] - uvcopy[3] )* ( uvcopy[1] - uvcopy[3] );
		if ( sqrt (norm ) < 1.e-08 || it > 100 ) {
			printf(" POINTS ARE TOO CLOSE TO EACH OTHER: %f it %d  BISECTION DIDN't CONVERGE\n", sqrt (norm), it );
			break;
		}
	}
	stat = EG_evaluate ( qm -> face, uvm, eval);
	if ( stat != EGADS_SUCCESS ) return stat;
	updateVertex (qm, iB, uvm);
	printf(" IT TOOK %d iTERTATIONS ::: NOW SIZE %f error %f\n", it, seg, fabs( seg - size ) / size );
	return EGADS_SUCCESS;
}



// Will detect if quad intersects but allows obtuse angles during regularization
static int
quadAngleOrientation(meshMap *qm, int qID, int *ori, int *order, double *theta ) {
	int     i, qV[4], k, k1, k2, sign, stat, vA, vB, vC, orip, orin, area = 0;
	double cross[4],qNormal[4], vABCD[12], dot, quv[5], projABCD[12];
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
	for ( orip = orin = k = 0 ; k < 4; ++k) {
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
		if      ( ori[k] ==  1 && orip == 0) orip = 1;
		else if ( ori[k] == -1 && orin == 0) orin = 1;
	}
	if      ( orip == 1 && orin == 0 ) area = 1;
	else if ( orip == 1 && orin == 1 ) area = 2;
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
				ori[k] == -1 && ( 2.0 * PI - theta[k] ) > PI * 1.1 ) area = 0;
		order[k] = k;
		if ( ori[k] == -1 ) sign = 1;
	}
	if ( sign == 0 ) area = 1;
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
	if ( area != 1 )
		printf(" ************   ATENCION:: INVALID QUAD ***************\n");
	printf("------------ AREA QUAD %d  IS %d Internal angles ordered by size --------------\n", qID, area );
	for ( sign = k = 0 ; k < 4; ++k)
		printf("Vertex %d has angle %f and orientation %d \n ", qV[order[k]] + 1, theta[order[k]], ori[order[k]]);
	return area;
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
	if ( qID <= 0 || qID > qm -> mesh -> totQuads ) {
		fprintf(stderr," ATTENTION!!!! %d", qID );
		exit(1);
	}
	if ( EG_nValenceCount ( qm -> mesh, qID, 3 ) == 0 ) return EGADS_SUCCESS;
	while ( j < 2 ) {
		centre = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + i3 ];
		if ( qm -> mesh -> vType [centre -1 ] == -1 && qm -> mesh -> valence[centre -1][0] == 3 ) {
			links[0] = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + (i3 + 1 ) % 4 ];
			links[1] = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + (i3 + 3 ) % 4 ];
			val  [0] = getValence ( qm -> mesh, links[0] );
			val  [1] = getValence ( qm -> mesh, links[1] );
			printf(" VALID COLLAPSE AT %d through %d -->%d \n", qID, centre, validCollapse (qm -> mesh, qID, centre) );
			if ( (val[0] >= 5 || val[1] >= 5 || j == 1) && validCollapse (qm -> mesh, qID, centre) == 1) {
				*activity = 1;
				printf(" FORCE collapse mesh quad %d \n", qID);
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
	int  stat, i, j, q,  adjq, adjPair[2], auxQ, oldQ[8], links[2], mergeToMidPoint = 0;
	int piv[4]  = {1, 0, 3, 2} ;
	vStar *star = NULL;
	// assuming ALWAYS that we collapse vertex 0 to 2
	i = EG_quadVertIdx  ( qm -> mesh, qC, centre );
	j = qm -> mesh -> quadIdx [ 4 * ( qC - 1) + ( i + 2 ) %4] - 1;
	if ( qm -> mesh -> vType [ centre - 1] == -1 &&  qm -> mesh -> vType [ j - 1 ] == -1 ) mergeToMidPoint = 1;
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
	i = vertexLinksToBounds (qm -> mesh, oldQ[2] + 1, nb );
	if ( mergeToMidPoint == 1 && (i == 0 || ( i > 0 && qm -> mesh -> vType [ nb[1] -1 ] < 4 ))) {
		qm -> mesh->uvs [2*oldQ[2]     ] = 0.5 * ( qm -> mesh->uvs [2*oldQ[2]     ] + qm -> mesh->uvs [2*oldQ[0]     ]);
		qm -> mesh->uvs [2*oldQ[2] + 1 ] = 0.5 * ( qm -> mesh->uvs [2*oldQ[2] + 1 ] + qm -> mesh->uvs [2*oldQ[0] + 1 ]);
		updateVertex (qm, oldQ[2] + 1, &qm -> mesh->uvs [2*oldQ[2]] );
	}
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
		if ( qg.q[0] <= 0 || qg.q[1] <= 0 ||
				qg.q[0] > qm -> mesh ->totQuads || qg.q[1] > qm -> mesh ->totQuads ) {
			printf(" quads %d %d \n", qg.q[0], qg.q[1] ); exit (1);
		}
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
	printMesh(qm, buffer,0);
	printf(" SWAPPING OPERATION ADAPT \n ");
	for ( j  = 0 ; j < 2; j++ ) {
		vc   = qg.verts[ j * 3 ];
		printf(" VERTEX %d type  %d linke %d \n", vc,
				qm -> mesh -> vType [ vc - 1], vertexLinksToBounds (qm -> mesh, vc, nb ) );
		if ( qm -> mesh -> vType [ vc - 1] != -1 ) continue;
		q    = 0 ; if ( EG_quadVertIdx ( qm ->mesh, qID[q], vc ) < 0 ) q = 1;
		printf(" Orientation area %d \n ", qID[q] );
		area = quadAngleOrientation(qm, qID[q], ori, order, angles);
		if ( area < 0  ) return area;
		if ( area == 1 ) continue;
		i0   = qg.verts[(3 * j + 1)%6];
		i1   = qg.verts[(3 * j + 5)%6];
		averageCoordsMinusLinks ( qm, vc, i0, i1 );
	}
	printMesh(qm, buffer,0);
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
			EG_freeStar ( &star);
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
	//weightedAverage ( qm, vC, &i );
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
				EG_freeStar(&star);
				return stat;
			}
		}
	}
	EG_freeStar(&star);
	// New Vertex location: take mid points of quads to which they belong
	updateVertex (qm, poly[3], &qm -> mesh -> uvs[ 2 * ( poly[0] - 1 )]);
	printMesh(qm, buffer,0);
	printf(" SPLITTINGOPERATION ADAPT\n");
	averageCoordsMinusLinks ( qm, poly[0], poly[1], poly[2] );
	averageCoordsMinusLinks ( qm, poly[3], poly[1], poly[2] );
	printf(" AFTER MINUS LINKS \n");
	printMesh (qm, buffer, 0 );
	stat = optimize_angles(qm, 4, poly, 0);
	printMesh(qm, buffer,0);
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
	printf(" DOUBLE SPLIT %d  FIRST SPLIT THRU %d %d \n ", qg.verts[0], qg.verts[piv[0]], qg.verts[piv[1]]);
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
	printf(" SWAP ANG SPLIT GROUP \n ");
	printQuadGroup (qm ->mesh, qg );
	printMesh(qm, buffer, 0);
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
		printf(" forcing:: v3oppp%d \n ", v3opp);
	}
	else return EGADS_SUCCESS;
	if ( v3opp == -1 ) {
		for ( i  = 0 ; i < 2; i++ ) {
			j    = 2 + 2 * i;
			if ( i == 0 ) vL5 = 1;
			else          vL5 = 5;
			if (  getValence (qm ->mesh, qg.verts[j]) == 3 ) continue;
			printf(" i %d : j %d LOOK AT ADJ  PAIR %d %d \n ", i, j, qg.verts[j], qg.verts[vL5] );
			stat              = EG_adjQtoPair (qm -> mesh, qg.q[i], qg.verts[j], qg.verts[vL5], adj );
			printf(" ADJACENT QUAD %d \n ", adj[1]);
			if ( stat        != EGADS_SUCCESS || adj[1] == -1 ) continue;
			q5                = EG_quadVertIdx (qm -> mesh, adj[1], qg.verts[j]);
			vL5adj            = qm -> mesh -> quadIdx [ 4 * ( adj[1] - 1 ) + ( q5 + 1 ) %4 ];

			if ( vL5adj == qg.verts[vL5] )
				vL5adj = qm -> mesh -> quadIdx [ 4 * ( adj[1] - 1 ) + ( q5 + 3) %4 ];
			printf(" VERTEX %d adj to %d not %d \n ", vL5adj, qg.verts[j], qg.verts[vL5] );
			if ( getValence (qm ->mesh, vL5adj) == 3 ) {
				i5   = j;
				swap = j;
				printf(" vertex %d valence 3 swap through %d \n ", vL5adj, swap );
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
	if ( swap %3 == 0 ||
			(qm -> mesh -> vType[ qg.verts[swap] - 1] >= 0 &&
					qm -> mesh -> vType[ qg.verts[(swap + 3 ) % 6] - 1] >= 0 ) ) return EGADS_SUCCESS;
	printf(" SWAPSPLIT SWAP \n ");
	stat       = EG_swappingOperation( qm, qg, swap );
	printMesh(qm, buffer, 0);
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
	printf(" Split %d thru %d %d\n ", qg.verts[i5], qg.verts[3], adj[0]);
	stat = EG_splittingOperation (qm, qg.verts[i5], qg.verts[3], adj[0]);
	if ( stat != EGADS_SUCCESS) {
		printf("In swapSplit splittingOperation --> %d !!\n", stat);
		return EGADS_INDEXERR;
	}
	printMesh(qm, buffer, 0);
	EG_freeStar(&star);
	if ( f == 1 && *activity > 0 ) {
		qm -> extraQuads++;
		//	fprintf(stderr,"swapSplit forcing extra quads %d\n ", qm -> extraQuads );
		printf("swapSplit forcing extra quads %d\n ", qm -> extraQuads );
	}
	printf(" LEAVE SWAP SPLIT \n ");
	printMesh(qm, buffer,0);
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
	if ( qID < 0 || qID > qm -> mesh -> totQuads ) {
		if ( qID > qm -> mesh -> totQuads ) fprintf(stderr, " quad is %d !!!! \n ", qID );
		return EGADS_SUCCESS;
	}
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
	ITMAX = 20;
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
		//stat      += optimize_angles (qm, 0, NULL, 2);
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
			printf(" 888888888888   TRANSFER VALENCES TYPE %d kk %d  \n ", k, kk );
			do {
				it++;
				totActivity = 0;
				stat = checkMesh (qm);
				stat = restoreMeshData (qm, qm -> backupMesh, qm -> mesh);
				fprintf(stderr, "********************  ITERATION %d : TOTACT %d : iV  %d  MESH %d ************************\n",
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
					if ( EG_quadIsBoundary ( qm ->mesh, q+ 1 ) != 0 && k == 0  ) continue;
					stat           = EG_transferValences ( qm, quadPair, kk, &transfer, &activity );
					if (stat      != EGADS_SUCCESS ) {
						stat       = restoreMeshData (qm, qm -> mesh, qm -> backupMesh);
						if ( stat != EGADS_SUCCESS ) return stat;
						activity  = 0;
					}
					if ( activity > 0 ) printMesh(qm, buffer,0);
					if ( activity == 0 || checkQuad( qm -> mesh, quadPair[0] ) != EGADS_SUCCESS  ) continue;
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
				printf(" MADE IT TILL HERE !!! \n ");
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
	int  i, j, stat = EGADS_SUCCESS,  it = 0, itMax, area;
	double dtheta, minT, maxT;
	vStar *star = NULL;
	qm -> vFix[0] = 0;
	if ( fullRegularization == 0 ) { // move around only affected vertices
		if ( nP == 0 || pList == NULL ) return EGADS_SUCCESS;
		printMesh(qm, buffer, 0);
		printf( "START OPTIMIZATION REGULARIZATION\n ");
		for ( j = 0 ; j < nP; j++) {
			printf(" Optimize angles around vertex %d \n ", pList [j] );
			if ( qm -> mesh -> vType [ pList[j] -1 ] == -2 ) continue;
			stat = EG_buildStar (qm -> mesh, &star, pList[j]);
			if ( stat != EGADS_SUCCESS || star == NULL ) return stat;
			for ( i = 0 ; i < star -> nQ; i++ ) {
				if ( star -> quads[i] == -1 ) continue;
				area = makePositiveAngles(qm, star -> quads[i], 0.0, 1.1 * PI );
				if ( area < 0 ) {
					printf(" Mesh is invalid for MIN MAX TOTALS [0, 200 deg]\n");
					printMesh(qm , buffer, 1);
					return EGADS_GEOMERR;
				}
			}
		}
		printf( "FINISH OPTIMIZATION REGULARIZATION\n ");
		EG_freeStar(&star);
	} else {
		itMax  = 20;
		minT   = 0.0;
		maxT   = PI;
		dtheta = 0.5 * PI / (double)itMax;
		printMesh(qm, buffer,0);
		printf("FULL  OPTIMIZE ANGLES -----------------------------   \n ");
		qm -> smooth  = 0;
		stat          = restoreMeshData (qm, qm -> backupMesh, qm -> mesh);
		while ( it <= itMax ) {
			if ( it > 5          ) qm -> smooth = 1;
			if ( it >= itMax - 2 ) qm -> smooth = 2;
			qm -> vFix[0] = 0;
			printf(" it %d max %d smooth %d\n ", it, itMax, qm -> smooth);
			stat = validBoundary (qm, minT, maxT);
			if ( stat != EGADS_SUCCESS ) {
				printf(" */*/*/*/*/ INVALID BOUNDARY DURING ITERATION %d LEAVE WITH BEST \n ", it );
				EG_free ( qm -> vFix );
				stat   = restoreMeshData (qm, qm -> mesh, qm -> backupMesh);
				if ( stat != EGADS_SUCCESS )
					return stat;
				return validBoundary (qm, minT, maxT);
			}
			minT += dtheta;
			maxT -= dtheta;
			if ( maxT < minT ) break;
			it++;
		}
		if ( stat != EGADS_SUCCESS) {
			snprintf(buffer,500, "thisMesh_%d.txt", qm -> fID);
			printMesh(qm , buffer, 1);
			stat   = restoreMeshData (qm, qm -> mesh, qm -> backupMesh );
			snprintf(buffer,500, "bestMesh_%d.txt", qm -> fID);
			printMesh(qm , buffer, 1);
			snprintf(buffer,500, "wvsBad_%i.txt", qm -> fID );
			stat  = EG_wvsData(qm -> mesh , buffer);
		}
	}
	if (stat != EGADS_SUCCESS ) printf(" MESH INVALID %d \n ", stat );
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
	printf(" Writing in File %s  \n", name);
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



/*static void
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
 */
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



void sampleNormalPlane(meshMap *qm, char *name, double *uv ) {
	double min[3], max[3], c, p[3], dt[3], aux, point[18], normal[4];
	int  i, j, k, nP;
	FILE *f;
	printf(" WRITIN IN %s \n", name);
	f = fopen(name,"w");
	if ( f == NULL) return ;
	nP = 20;
	c = 0.0;
	i = EG_evaluate(qm -> face, uv, point);
	if ( i != EGADS_SUCCESS) return ;
	//
	if ( qm -> face -> mtype == SREVERSE )
		cross_product(&point[6], &point[3], normal);
	else
		cross_product(&point[3], &point[6], normal);
	for (i = 0 ; i < 3; ++i) {
		min  [i] = point[i] -1; //qm -> xyzs[3*vID + i] - 100 *qm -> xyzs[3*vID + i] ;
		max  [i] = point[i]  +1;//0.5;// qm -> xyzs[3*vID + i] + 100 *qm -> xyzs[3*vID + i] ;
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







