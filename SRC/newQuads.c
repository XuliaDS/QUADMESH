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
#include <nlopt.h>
#include <time.h>
#ifdef  WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <winsock2.h>
#endif
#include "wsserver.h"
#include "egads.h"

#define FATAL_ERROR -100
#define PI            3.1415926535897931159979635
#define MINVALIDANGLE 0.0//0.0872664625997165   // 5 degrees
#define MAXVALIDANGLE 3.1415926535897931159979635//3.05432619099008     // 175 degrees
static int PRINTAREA, MESHPLOTS = 0;
#define ISVERTEX 0
#define ISQUAD   1
#define DIAGONAL 1
#define ADJACENT 0


#define INTERIOR 0
#define FULL     1
#define SINGLE   0
#define DOUBLE   1


#define EPSAREA 1.E-01
#define EPSANGLE 0.0349065850398866  // +- 2 degrees
#define EPS08   1.E-08
#define EPS11   1.E-11
#define EPS14   1.E-14
#define EPS04   1.E-04
#define EPS02   1.E-02
#define EPS06   1.E-06

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
	int      fID, actualQuads, oriQ, oriV, minQ, extraQuads, optNquads;
	ego      face;
	double   MINANGLE, MAXANGLE;
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
static float     focus[4];
static char buffer[500];

static void
swapInt ( int *a, int *b ) {
	int c;
	c  = *a;
	*a = *b;
	*b =  c;
}

static void
EG_unitVector(double *v, double *norm) {
	double n;
	n     = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	n     = sqrt(n);
	if(n > EPS14 ) {
		v[0] /=n; v[1] /=n; v[2] /=n;
	}
	else {
		v[0] = 0.0; v[1] = 0.0; v[2] = 0.0;
	}
	*norm = n;
	return;
}


static void
EG_freeStar ( /*@null@*/ vStar **star ) {
	if ( *star == NULL ) return;
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

static int
EG_allocMeshData ( meshData **mesh, int nQ, int nV ) {
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
	if ((*mesh)->quadIdx  == NULL || (*mesh)->quadAdj  == NULL ||
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

/* IO FUNCTIONS */
static void printStar               (vStar *star) ;
//static void printStarFile(vStar *star, meshMap *qm) ;
static void printMeshStats(meshMap *qm, int sweep);
static void printMesh     (meshMap *qm,char[], int  );
static int  checkMesh     (meshMap *qm) ;

static void cross_product            (double *A, double *B, double *cross);
static double dotProduct             (double *v1, double *v2 );
static int  EG_projectToTangentPlane (double normal[], double *O, double *p, double *proj) ;
/***********************************************************/
/* MESH MAP FUNCTIONS*/
static int    EG_createMeshMap        (bodyData *bodydata);
static double EG_angleAtVnormalPlane  (meshMap *qm, int vC, int v1, int v2 );
static double EG_angleAtBoundaryVertex(meshMap *qm, int v ) ;
static int    checkInvalidElement     (meshMap *qm, int qID, double minAngle, double maxAngle ) ;
static int    quadAlgebraicArea       (meshMap *qm, double minTheta, double maxTheta, int qID, int vID, int *validArea, double *angles,double *area) ;
static int    vectorAtVertexNplane    (meshMap *qm, int v0, int v1, double *uvEPS );
static int    quadAverageCoords       (meshMap *qm, int q, double *uv, double *p) ;
static int    EG_segment              (meshMap *qm, int v1, int v2, int qv, double *seg );
static int    optimize_angles         (meshMap *qm, int nP, /*@unused@*/ /*@null@*/int *pList, int fullRegularization);
/* REGULARIZATION FUNCTIONS */
/* FUNDAMENTAL */
static int  EG_swappingOperation      (meshMap *qm, quadGroup qg, int swap);
static int  EG_splittingOperation     (meshMap *qm, vStar *star, int id0, int distanceToId0 ) ;
static int  EG_mergeVertices          (meshMap *qm, int qC, int centre);

/**************/
static int  EG_swap                   (meshMap *qm, int qID, int *activity);
static int  EG_cleanMesh              (meshMap *qm, int interior, int adj, int *activity );
static int  EG_cleanNeighborhood      (meshMap *qm, int qID, int adjacent, int transfer, int *activity );
static int  EG_forceCollapse          (meshMap *qm, int  vID, int *activity);
static int  EG_forceSplit             (meshMap *qm, int v0, int vL, int *activity );
static int  EG_collapse               (meshMap *qm, int  vID, int *activity);
static int  EG_split                  (meshMap *qm, int  qID, int *activity);
static int  EG_removeDoublet          (meshMap *qm, int  vID ) ;
static int  EG_doubleSwap             (meshMap *qm, quadGroup qg, int forcing, int *activity );
//static int  EG_doubleSwapDiag         (meshMap *qm, quadGroup qg, int *activity );
//static int  EG_collapseSwapDiag       (meshMap *qm, quadGroup qg, int *activity );
static int  EG_swapSplit              (meshMap *qm, quadGroup qg, int forcing, int *activity );
static int  EG_swapCollapse           (meshMap *qm, quadGroup qg, int forcing, int *activity );
static int  EG_doubleCollapse         (meshMap *qm, quadGroup qg, int forcing, int *activity );
static int  EG_doubleSplit            (meshMap *qm, quadGroup qg, int *activity );
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
static int  quadValenceDistribution   (meshData *qm, int qID, int *idx, int *v ) ;
static int  checkQuad                 (meshData *qm, int q );
static int  checkVertex               (meshData *qm, int v );
static int  EG_quadIsBoundary         (meshData *qm, int qID );
static int  EG_adjQtoPair             (meshData *qm, int qID, int v1, int v2, int *adj) ;
static void EG_commonVerts            (meshData *qm, int q1, int q2, int *v );
static int  EG_createQuadGroup        (meshData *qm, quadGroup *qg, int q0, int q1 ) ;
static int  getValence                (meshData *qm, int v );

static int linkedVertices             (meshData *qm, int v1, int v2 );
static int  validSwap                 (meshData *qm, int v1, int v2 );
static int  validCollapse             (meshData *qm, int qID, int v );
/*********************/

static int
EG_segment ( meshMap *qm, int id1, int id2, int datatype, double *segment  ) {
	int i, j, stat , n = 100;
	double p1[18], p2[18], dist[2], uv1[2], uv2[2], uvEps[2], seg = 0.0, totArc = 0.0, dt ;
	if ( datatype == ISVERTEX ) {
		for ( i = 0 ; i < 2; i++ ) {
			uv1[i] = qm -> mesh -> uvs [ 2 * ( id1 - 1 ) + i];
			uv2[i] = qm -> mesh -> uvs [ 2 * ( id2 - 1 ) + i];
		}
	}
	else {
		stat  = quadAverageCoords(qm, id1, uv1, p1);
		stat += quadAverageCoords(qm, id2, uv2, p2);
		if ( stat != EGADS_SUCCESS ) {
			printf(" In EG_segment quad option for %d %d average coords %d \n ", id1, id2, stat );
			return stat;
		}
	}
	stat = EG_evaluate (qm -> face, uv1, p1);
	dist[0] = uv2[0] - uv1[0];
	dist[1] = uv2[1] - uv1[1];
	if ( stat != EGADS_SUCCESS ) {
		printf(" In EG_segment evaluate at %lf %lf -> %d \n ", uv1[0], uv1[1], stat );
		return stat;
	}
	for ( i = 0; i < n; i++ ) {
		dt = (double ) ( i + 1 ) / (double ) n;
		uvEps[0] = uv1[0] + dt * dist[0];
		uvEps[1] = uv1[1] + dt * dist[1];
		stat     = EG_evaluate (qm -> face, uvEps, p2 );
		seg      = 0.0;
		for ( j = 0 ; j < 3;j ++ ) {
			seg  += (p1[j] - p2[j] ) * (p1[j] - p2[j] );
			p1[j] = p2[j];
		}
		totArc += sqrt ( seg );
	}
	seg = 0.0;
	*segment = totArc;
	return EGADS_SUCCESS;
}

static int linkedVertices ( meshData *mesh, int v1, int v2 ) {
	int i;
	for ( i = 0 ; i < mesh ->valence [v1 -1 ][0]; i++ )
		if ( mesh ->valence [ v1 -1 ][2 + i] == v2 ) return 1;
	return 0;
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
		printf(" EG_commonVerts you are trying to operate on a bad quad: %d --> %d!!\n", q1, checkQuad ( mesh, q1 ) );
		printQuadSpecs (mesh, q1 );
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

static int
EG_quadVertIdx ( meshData *mesh, int q, int v ) {
	int i = 0 ;
	for ( i = 0 ; i < 4; i++ )
		if ( mesh -> quadIdx [ 4 * ( q - 1 ) + i ] == v ) return i;
	return -1;
}

static int
quadValenceDistribution ( meshData *mesh, int qID, int *idx, int *v ) {
	int stat, i, vPair[2], k, i3, i5;
	stat = checkQuad ( mesh, qID );
	if ( stat != EGADS_SUCCESS ) return stat;
	for ( i = 0 ; i < 15; i++ ) {
		if ( i % 5 == 0 ) {
			idx[i] = 0;
			v  [i] = 0;
		}
		else {
			idx[i] = -1;
			v  [i] = -1;
		}
	}
	/* v[0] = total vertices with valence 3, v[1], .., v[4] ids of valence 3
	 * v[5] = total vertices with valence 4, etc
	 */
	for ( i = 0 ; i < 4; i++ ) {
		vPair[0] = mesh -> quadIdx [ 4 * ( qID - 1 ) + i ];
		vPair[1] = getValence ( mesh, vPair[0] );
		if      ( vPair[1] == 3 ) k = 0;
		else if ( vPair[1] == 4 ) k = 1;
		else                      k = 2;
		idx[ 5 * k ]++;
		v  [ 5 * k ]++;
		idx[ 5 * k + idx[5 * k ] ] = i;
		v  [ 5 * k +   v[5 * k ] ] = vPair[0];
	}
	if ( idx[10] > 0 && idx[0] > 0 ) {
		for (     i5 = 0 ; i5 < idx[10]; i5++ ) {
			for ( i3 = 0 ; i3 < idx [0]; i3++ ) {
				if ( abs ( idx[1 + i3] - idx[11 + i5 ] % 2 == 1 ) ) {
					swapInt (&idx[1 ], &idx[ 1  + i3]);
					swapInt (&idx[11], &idx[ 11 + i5]);
					swapInt (&v  [1 ], &v  [ 1  + i3]);
					swapInt (&v  [11], &v  [ 11 + i5]);
					return EGADS_SUCCESS;
				}
			}
		}
	}
	return EGADS_SUCCESS;
}

static int getValence ( meshData *mesh, int v )  {
	int val;
	val = checkVertex (mesh, v );
	if ( val != EGADS_SUCCESS ) return val;
	val =     mesh -> valence [ v - 1 ][0];
	if      ( mesh -> vType   [ v - 1 ] == 0 ) val += 2;
	else if ( mesh -> vType   [ v - 1 ] >  0 && mesh -> vType [ v -1 ] != 2 * mesh -> sizeQuads ) val++;
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
	if ( mesh ->vType [ v - 1 ] == -2 ) return EGADS_EMPTY;
	if (( mesh ->vType [ v -1  ] == 0 && mesh ->valence [ v - 1][0] < 2 ) ||
			( mesh ->vType [ v -1  ]  > 0 && mesh ->valence [ v - 1][0] < 3 ) ) {
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


static int
checkQuad ( meshData *mesh, int q ) {
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


static int
checkMesh(meshMap *qm) {
	int stat, i, j, k,  val1, val2, v1, v2;
	for ( i = 0 ; i < qm -> mesh ->totVerts; i++ ) {
		if (qm -> mesh -> vType[i] == -2 )  continue;
		if (qm -> mesh -> vType[i] != 0 && qm -> mesh -> valence[i][0] == 2 ) {
			printQuadSpecs ( qm -> mesh , qm -> mesh -> valence [i][1] );
			fprintf(stderr," checkMesh is invalid: we have a  doublet at vertex %d type %d\n", i + 1, qm -> mesh -> vType[i] );
			snprintf(buffer,500,"WRONGMESH_%d", qm -> fID);
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
			v2  = qm -> mesh  -> valence[i ][2 + j];
			if ( linkedVertices (qm -> mesh , v2, v1 ) == 0  ) {
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




static void
meshCount( meshData *mesh, int *nI, int *nV, int *nQ ) {
	int i, qSum, vSum = 0, vSum2 = 0;
	for ( vSum = i = qSum = 0; i< mesh -> totQuads;i++ )
		if ( checkQuad ( mesh, i + 1 ) == EGADS_SUCCESS ) qSum++;
	for ( vSum =  i = 0 ; i < mesh -> totVerts; i++ ) {
		if ( mesh -> vType[i] == -2 ) continue;
		vSum2++;
		if ( getValence ( mesh, i + 1 ) != 4) vSum++;
	}
	// fprintf(stderr, " In this round we have total quads %d (ori %d) and %d are irregular\n",
	//	  qSum,mesh -> q0,vSum);
	printf( " Current mesh has total quads %d and %d are irregular ( %d TOTAL )\n",
			qSum, vSum, vSum2);
	*nQ = qSum;
	*nI = vSum;
	*nV = vSum2;
}


/*
 * Assuming qID = v1, vC1, v2, vC2 and collapsing vC1 -> vC2
 */
static int
validCollapse ( meshData *mesh, int qID, int v ) {
	int j, k, id, link, aux ;
	id = EG_quadVertIdx ( mesh, qID, v);
	for ( j = 0 ; j < 2; j++ ) {
		link = mesh ->quadIdx [ 4 * ( qID - 1 ) + (id + 2 * j + 1)%4];
		printVertSpecs(mesh, link);
		if      (mesh ->vType [link - 1]  >  0 && mesh ->valence[link - 1][0] <= 3 ) return 0;
		else if (mesh ->vType [link - 1] ==  0 && mesh ->valence[link - 1][0] <= 2 ) return 0;
		else if (mesh ->vType [link - 1] == -1 && mesh ->valence[link - 1][0] == 3 ) {
			for ( k = 0 ; k < mesh ->valence[ link - 1 ][0]; k++ ) {
				aux = mesh ->valence [ link - 1 ][2 + k] - 1;
				if ( EG_quadVertIdx ( mesh, qID,  aux + 1 ) >= 0 ) continue;
				if ( mesh ->valence[ aux ][0] == 3 && mesh ->vType[aux] != 0 ) return 0;
			}
		}
	}
	return 1;
}

static int
validSwap ( meshData *mesh, int v1, int v2 ) {
	int i, vs[2];
	vs[0] = v1 - 1; vs[1] = v2 - 1;
	for ( i = 0; i < 2; i++ ) {
		if      (checkVertex ( mesh, vs[i] + 1 ) != EGADS_SUCCESS )                return 0;
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
static int
EG_quadIsBoundary ( meshData *mesh, int qID ) {
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



static int
quadAverageCoords(meshMap *qm, int q, double *uv, double *xyz ) {
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




static int
EG_createQuadGroup ( meshData *mesh, quadGroup *qg, int q0, int q1 ) {
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
	if ( common[0] == 2 ) {
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
	}
	return EGADS_SUCCESS;
}

static int
EG_doubleCollapse ( meshMap *qm, quadGroup qg, int forcing, int *activity ) {
	int i, act = 0, stat, i3, i5, q3;
	*activity = 0;
	printf("DOUBLE COLLAPSE  \n");
	printQuadGroup (qm -> mesh, qg);
	if (qg.vals[0] * qg.vals[3] == 16 ) {
		for ( i3 = 0 ; i3 < 6; i3++ ) if ( qg.vals[i3] == 3 ) break;
		if  (qg.vals[i3] != 3 ) return EGADS_SUCCESS;
		i5 = ( i3 + 1 ) % 6; if ( i5 % 3 == 0 ) i5 = ( i3 + 5 ) % 6;
		if ( forcing == 0 && ( qg.vals[i5]!= 5 ||
				(qg.vals[( i5 + 3) % 6]  != 5 && qg.vals [ ( i3 + 3 ) % 6 ] != 3 ) ) ) return EGADS_SUCCESS;
		else if ( forcing == 1) {
			//if ( qm -> extraQuads < 0 ) return EGADS_SUCCESS;
			printf(" FORCING \n ");
			if ( qg.vals[(i3 + 3)%6] == 3 ) i5 = (i3 + 3)%6;
			else if ( qg.vals[(i3 + 2)%6] == 5 ) i5 = (i3 + 2)%6;
			else if ( qg.vals[(i3 + 4)%6] == 5 ) i5 = (i3 + 4)%6;
			else return EGADS_SUCCESS;
			printf(" found something:: i5 %d \n ", i5 );
		}
	}
	else if (qg.vals[0] * qg.vals[3] == 12 ) {
		i3 = 0 ; if ( qg.vals[3] == 3 ) i3 = 3;
		i5 = ( i3 + 1 )% 6; if ( qg.vals[i5] != 5 ) i5 = ( i3 + 5 )% 6;
		if ( qg.vals[i5] != 5 ) return EGADS_SUCCESS;
		if ( qg.vals [ ( i5 + 3 ) % 5] < 5 ) return EGADS_SUCCESS;
	}
	else return EGADS_SUCCESS;
	printf(" DOUBLE COLLAPSE IS HAPPENING \n");
	printMesh(qm, buffer,0);
	if ( forcing == 1 && qg.vals[0] * qg.vals[3] == 16 ) {
		q3 = 0; if ( i3 >= 3 ) q3 = 1;
	}else {
		q3 = 0; if ( i5 >= 3 ) q3 = 1;
	}
	for ( i = 0 ; i < 2; i++ ) {
		printf(" Try collapsing at quad %d i = %d\n", qg.q[(q3+i)%2], i);
		stat       = EG_forceCollapse (qm, qg.q[( q3 + i ) %2], &act);
		*activity += act;
		if ( stat != EGADS_SUCCESS || act == 0 ) {
			printf(" EG_doubleCollapse: Force collapse went %d and act %d !! \n", stat, *activity );
			return stat;
		}
		printMesh(qm, buffer,0);
	}
	return stat;
}


static int
EG_swapDoubleCollapse (meshMap *qm, quadGroup qg, int *activity ) {
	int  k, swap, id, j, stat;
	*activity  = 0;
	if ( qg.vals[0] * qg.vals[3] != 20 || qg.vals[2] * qg.vals[4] != 9 ) return EGADS_SUCCESS;
	swap = 1;
	if ( qg.vals[swap] * qg.vals[(swap+3)%6] != 12 ) swap = 2;
	if ( qg.vals[swap] * qg.vals[(swap+3)%6] != 12 ) return EGADS_SUCCESS;
	if ( validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0 ) return EGADS_SUCCESS;
	stat       = EG_swappingOperation(qm, qg, swap);
	*activity  = 1;
	printMesh(qm, buffer,0);
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
	stat    = EG_collapse ( qm, qg.q[k], &j );
	printf(" Swap-Double Collapse result   %d   ACTIVITY %d \n ", stat, j );
	printMesh(qm, buffer,0);
	if ( stat != EGADS_SUCCESS || j != 1 ) {
		if  ( stat != EGADS_SUCCESS)
			printf(" In EG swapDoubleCollapse: I failed to collapse after swapping! s = %d act = %d \n ", stat, j );
		else printf(" LEAVE SWAPPED AND DO NOTHING\n");
	}
	return stat;
}





static int
EG_swapDoubleSplit (meshMap *qm, quadGroup qg, int *activity ) {
	int  i5, q, i55, i3, val3, i0, v30[2], i,  stat, adj[2], q0[2];
	*activity  = 0;
	printf(" IN SWAP DOUBLE SPLIT FUNCTION QUADS %d  %d \n", qg.q[0], qg.q[1]);
	printQuadGroup (qm->mesh, qg);
	for ( i3 = 0 ; i3 < 6; i3++)
		if ( qg.vals[i3] == 3 ) break;
	if ( qg.vals[i3] != 3 || i3 % 3 == 0 ) {
		printf(" RETURN \n ");
		return EGADS_SUCCESS;
	}
	i5  =(i3 + 3)% 6;
	if ( qg.vals[i5] < 5 ) return EGADS_SUCCESS;
	i55 = (i5 + 1)%6;
	printQuadGroup (qm->mesh, qg);
	printf(" i3 %d i5 %d i55 %d\n ", i3, i5, i55 );
	if ( qg.vals[i55] < 5 ) {
		i55 = (i5 + 5)%6 ;
		printf(" i3 %d i5 %d i55 %d\n ", i3, i5, i55 );
		if ( qg.vals[i55] != 5 ) return EGADS_SUCCESS;
	}
	printQuadGroup (qm->mesh, qg);
	q = 0;
	if ( EG_quadVertIdx  (qm -> mesh, qg.q[q], qg.verts[i5]) < 0 ) q = 1;
	stat = EG_adjQtoPair (qm -> mesh, qg.q[q], qg.verts[i5], qg.verts[i55], adj );
	if ( stat != EGADS_SUCCESS || adj[1] == -1 ) return stat;
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
	printMesh (qm, buffer, 0);
	stat = EG_swappingOperation (qm, qg, i0 );
	printMesh (qm, buffer, 0);
	if ( stat != EGADS_SUCCESS ) {
		printf(" EG_swapDoubleSplit error at swap: %d !!\n ", stat );
		return stat;
	}
	i = 0 ; if ( EG_quadVertIdx (qm -> mesh, q0[0], val3) < 0 ) i = 1;
	stat = EG_adjQtoPair (qm -> mesh, q0[i], v30[0], v30[1], adj );
	if ( stat != EGADS_SUCCESS || adj[1] == -1 ) {
		printf(" EG_swapDoubleSplit after swapping adjacent to pair %d %d is %d --> %d!!\n ", v30[0], v30[1], adj[1], stat);
		return stat;
	}
	stat = EG_createQuadGroup (qm -> mesh, &qg,  q0[i], adj[1]);
	if ( stat != EGADS_SUCCESS ) {
		printf("Inside EG_swapDoubleSplit: before splitting EG_createQuadGroup stat %d\n ", stat );
		printQuadGroup (qm -> mesh, qg );
		return stat;
	}
	printf(" GOTO DOUBLE SPLIT\n ");
	return EG_doubleSplit (qm, qg, &i);
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
			if ( stat != EGADS_SUCCESS) printf("EG_transferValences: EG_cleanQuad %d --> %d!!\n", qID[0], stat);
			if ( checkQuad (qm -> mesh, qID[0] ) != EGADS_SUCCESS)  qID[0] = -1;
			return stat;
		}
	}
	if (EG_nValenceCount (qm -> mesh, qID[0], 4) == 4 ) return EGADS_SUCCESS;
	printf("START TRANSFER VAL TRANSFER VALENCES AT Q %d coming from %d TRY 5533 %d  transfering %d  \n",
			qID[0], qID[1], try5533, *transfering );
	printQuadSpecs   ( qm -> mesh, qID[0] );
	for ( swap = j = 0 ; j < 4; j++ ) {
		qg.q[0] = qID[0];
		qg.q[1] = qm -> mesh -> quadAdj[ 4 * ( qID[0] - 1 ) + j];
		if ( qg.q[1] < 0 || qg.q[1] == qID[1] ) continue;
		if (*transfering == 0 ) {
			stat = EG_cleanQuad ( qm, qg.q[1], 1, 0,  try5533 , &(*activity) );
			printf(" EG_CLEAN QUAD PRINT %d activity %d \n ", stat, *activity );
			if ( stat != EGADS_SUCCESS || *activity > 0 ) {
				if ( stat != EGADS_SUCCESS) printf(" EG_TransferValence clean quad %d !!\n ", stat );
				qID[0]    = -1;
				return stat;
			}
		}
		qg.q[0] = qID[0];
		qg.q[1] = qm -> mesh -> quadAdj[ 4 * ( qID[0] - 1 ) + j];
		printf(" Create Quad group qith %d %d \n ", qg.q[0] , qg.q[1]);
		stat       = EG_createQuadGroup  (qm -> mesh, &qg, qg.q[0], qg.q[1] );
		if ( stat != EGADS_SUCCESS  ) {
			printf(" Inside EG_transferValences EG_createQuadGroup %d !!\n", stat );
			printQuadGroup (qm -> mesh, qg );
			return stat;
		}
		printf(" We are at %d coming from %d \n\n", qID[0], qID[1] );
		if ( *transfering == 0 && qg.vals[0] * qg.vals[3] == 15 ) {
			j = 0;
			if ( qg.q[j] == qID[0] ) j = 1;
			qID[0] = qg.q[j]; qID[1] = qg.q[(j + 1)%2];
			printf("-> rather use the following FROM TRANSFER 3,5 pair call now CALL TRANSFER WITH %d  %d  TRANSFER %d \n", qID[0], qID[1], *transfering );
			if ( EG_quadIsBoundary(qm -> mesh, qg.q[1] ) == 1 && qm -> extraQuads >= 0 ) {
				printf(" Extra quads %d -> force a collapse at %d\n ", qm -> extraQuads, qID[0] );
				stat = EG_forceCollapse (qm, qID[0], &(*activity));
				printMesh (qm, buffer,0);
				qID[0] = qID[1];
				qID[1] = -1;
				if ( stat != EGADS_SUCCESS ) printf("EG_transferValences: forceCollapse gave %d !!\n", stat);
				return stat;
			}
			else if ( EG_quadIsBoundary(qm -> mesh, qg.q[1] ) == 1 )
				return EG_transferValences ( qm, qID, 0, &(*transfering),  &(*activity) );
		}
		if (validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0 ) continue;
		min     = 0;
		if (qg.vals[0] * qg.vals[3] >= 20 ) {
			if (qm -> mesh -> vType [ qg.verts[3] - 1] * qm -> mesh -> vType [ qg.verts[0] - 1] == 0 ) min = 16;
			else                                                                                       min = 12;
		}
		else if ( try5533 == 1 && *transfering == 0  ) {
			if      (qg.vals[0] * qg.vals[3] >= 25) min = 16;
			else if (qg.vals[0] * qg.vals[3] == 16) min =  9;
		}
		if ( min == 0 ) continue;
		for( i = 1 ; i <= 2; i++ )
			if (qg.vals[i] * qg.vals [3 + i] <= min ) swap = i;
		if ( swap == 0 ) continue;
		printf(" Min swap was  %d swap through %d\n ", min, swap );
		printMesh(qm, buffer,0);
		printf(" QUAD GROUP %d  %d  ( PREV %d ) \n ", qg.q[0], qg.q[1], qID[1] );
		stat       = EG_swappingOperation (qm, qg, swap);
		*activity  = 1;
		printMesh(qm, buffer,0);
		if ( stat != EGADS_SUCCESS ) return stat;
		i = 0;
		if ( min == 9 || min == 16) *transfering = 1;
		else                        *transfering = 0;
		printf(" NOW Qs ARE %d %d and transfer  %d\n ", qg.q[0], qg.q[1], *transfering  );
		if ( qg.q[0] != qID[0] && qg.q[1] != qID[0] ) {
			fprintf(stderr, " THIS DOES HAPPEN EVENTUALLY....\n ");
			qID[1] = qID[0];
			EG_commonVerts ( qm -> mesh, qg.q[0], qID[1], links);
			if ( links[0] != 2 ) {
				EG_commonVerts ( qm -> mesh, qg.q[1], qID[1], links);
				if ( links[0] == 2 )
					i = 1;
			}
			qID[0] = qg.q[i];
		}
		else {
			if ( qID[1] != -1 ) {
				printf(" Common verts to quads %d  %d \n ", qg.q[0], qID[1] );
				if      ( qg.q[0] == qID[1] ) i = 1;
				else if ( qg.q[1] == qID[1] ) i = 0;
				else {
					EG_commonVerts ( qm -> mesh, qg.q[0], qID[1], links);
					printf(" Common to %d %d = %d \n ", qg.q[0], qID[1], links[0] );
					if ( links[0] > 0 ) {
						if ( links[0] == 1 ) {
							EG_commonVerts ( qm -> mesh, qg.q[1], qID[1], links);
							if ( links[0] == 0 )
								i = 1;
						}
						else i = 1;
					}
				}
			}
			else if (EG_nValenceCount(qm -> mesh,qg.q[0],5) == 2 &&
					EG_nValenceCount (qm -> mesh,qg.q[1],5) == 1 && EG_nValenceCount(qm -> mesh,qg.q[1],3) == 1) i = 1;
			else if ( EG_nValenceCount(qm -> mesh,qg.q[1], 4) < EG_nValenceCount(qm -> mesh,qg.q[0],4) ) i = 1;
			printf(" i  = %d \n ", i );
			qID[0] = qg.q[i];
			qID[1] = qg.q[ ( i + 1 ) %2 ];
		}
		if (*transfering == 1 ) {
			qAux[0] = qID[1];
			qAux[1] = qID[0];
			printf(" TRANSFERRING::: SEPARATE VALENCES %d %d\n", qAux[0] , qAux[1] );
			stat = EG_transferValences ( qm, qAux, 0, &i, &j);
			if ( stat != EGADS_SUCCESS ) return stat;
		}
		if ( *activity > 0 ) break;
	}
	printf(" LEAVE TRANSFER WITH  %d  %d  TRANSFERRING %d ACTIVITY %d!\n", qID[0], qID[1], *transfering, *activity );
	return stat;
}













static int EG_basicOperation (meshMap *qm, int qID, int type, int *activity ) {
	int stat = EGADS_SUCCESS;
	switch (type ) {
	case SWAP:
		//printf(" Try Swapping \n");
		stat     = EG_swap ( qm, qID, & (*activity));
		break;
	case COLLAPSE:
		//printf(" Try Collapsing \n");
		stat  = EG_collapse ( qm, qID, & (*activity));
		break;
	case SPLIT:
		//printf(" Try Splitting \n");
		stat  = EG_split ( qm, qID, & (*activity));
		break;
	}
	if ( stat != EGADS_SUCCESS) printf(" EG_basicOperation %d  around %d --> %d!!\n", type, qID, stat);
	if ( *activity >0  ) printf(" basicOperation for quad %d type %d DONE\n ", qID, type );
	return stat;
}

static int EG_composeOperation (meshMap *qm, quadGroup qg, int type, int forcing, int *activity ) {
	int stat = EGADS_SUCCESS;
	switch (type ) {
	case DOUBLESWAP:
		//printf(" Try doubleSwap \n");
		stat  = EG_doubleSwap ( qm, qg, forcing, & (*activity));
		break;
	case SWAPCOLLAPSE:
		//printf(" Try swapCollapse \n");
		stat     = EG_swapCollapse ( qm, qg,forcing, & (*activity));
		break;
	case DOUBLECOLLAPSE:
		// printf(" Try doubleCollapse \n");
		stat  = EG_doubleCollapse ( qm, qg, forcing,  & (*activity));
		break;
	case SWAPDOUBLECOLLAPSE:
		// printf(" Try swapDoubleCollapse \n");
		stat  = EG_swapDoubleCollapse ( qm, qg, & (*activity));
		break;
	case SWAPSPLIT:
		// printf(" Try swapSplit \n");
		stat     = EG_swapSplit ( qm, qg, forcing, & (*activity));
		break;
	case DOUBLESPLIT:
		//   printf(" Try doubleSplit\n");
		stat  = EG_doubleSplit ( qm, qg, & (*activity));
		break;
	case SWAPDOUBLESPLIT:
		//     printf(" Try swapDoubleSplit\n");
		stat  = EG_swapDoubleSplit ( qm, qg, & (*activity));
		break;
	}
	if ( stat != EGADS_SUCCESS) {
		printf(" EG_composeOperationOperation %d --> %d!!\n", type, stat);
		printQuadGroup (qm -> mesh, qg);
	}
	if ( *activity > 0 ) printf(" composeOperation for quad pair %d %d type %d DONE SOMETHING\n ", qg.q[0], qg.q[1], type );
	return stat;
}



static int EG_cleanQuad (meshMap *qm, int qID, int useAdj, int transfer, int forcing, int *activity ) {
	int stat, i, q, qadj, act = 0, imax;
	int opBasic[3] = {COLLAPSE, SWAP, SPLIT};
	int opComp [7] = {SWAPCOLLAPSE, DOUBLECOLLAPSE, SWAPDOUBLECOLLAPSE, DOUBLESWAP, SWAPSPLIT,  DOUBLESPLIT,  SWAPDOUBLESPLIT};
	quadGroup qg;
	*activity      = 0;
	stat           = checkQuad ( qm -> mesh, qID  );
	if ( stat     != EGADS_SUCCESS ) {
		printf(" EG_cleanQuad %d stat %d\n", qID, stat );
		return stat;
	}
	if ( EG_nValenceCount (qm -> mesh, qID, 4) == 4 ) {
		// printf(" Regular Quad. Leave\n");
		return EGADS_SUCCESS;
	}
	printf(" Inside clean quad %d using adj %d BASIC forcing functions %d \n ", qID, useAdj, forcing );
	if ( transfer == 0 ) {
		if ( EG_nValenceCount (qm -> mesh, qID, 6) > 0 || qm -> extraQuads <= 0 ) swapInt( &opBasic[2], &opBasic[0]);
		for ( i = 0 ; i < 3; i++ ) {
			stat      = EG_basicOperation(qm, qID, opBasic[i], &act );
			*activity = act;
			if ( stat != EGADS_SUCCESS || act > 0 ) {
				printf(" stat in basic %d act %d \n ", stat, act);
				return stat;
			}
		}
	}
	printf(" Done with basic\n ");
	if ( useAdj == 0  ) return EGADS_SUCCESS;
	printf(" EXTRAQUADS %d \n ", qm -> extraQuads );
	if ( qm -> extraQuads <= 0 ) {
		printf(" EXTRAQUADS %d SO REVERSEORDER\n ", qm -> extraQuads );
		for ( i = 0 ; i < 3; i++) swapInt( &opComp[i], &opComp[i+4]);
	}
	//if ( qm -> extraQuads == 0 ) imax = 7;
	//else                         imax = 4;
	imax = 7;
	printf(" Inside clean quad %d using adj %d COMPOSE \n ", qID, useAdj );
	stat           = checkQuad ( qm -> mesh, qID  );
	if ( stat     != EGADS_SUCCESS ) {
		printf(" EG_cleanQuad %d stat %d\n", qID, stat );
		return stat;
	}
	if ( forcing == 1 && qm -> extraQuads >= 0 && EG_nValenceCount (qm -> mesh, qID, 3) > 1 ) {
		printf(" FORCE COLLAPSE FOR QUAD %d valencecount %d\n ", qID, EG_nValenceCount (qm -> mesh, qID, 3));
		stat = EG_forceCollapse (qm, qID, &(*activity ) );
		if ( *activity > 0 || stat != EGADS_SUCCESS ) return stat;
	}
	for ( i  = 0 ; i < imax; i++ ) printf(" GROUP %d\n", opComp[i]);
	for ( i  = 0 ; i < imax; i++ ) {
		for ( q  = 0 ; q < 4; q++) {
			printf(" q = %d i %d \n ", q, i );
			qadj = qm -> mesh -> quadAdj [ 4 * ( qID -1 ) + q];
			if ( qadj == -1 ) continue;
			printf(" USING GROUP %d  %d  OP TYPE %d \n ", qID, qadj , opComp[i]);
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
	printf("LEAVE CC\n ");
	return stat;
}

static int EG_createMeshMap(bodyData *bodydata)
{
	int         f, stat = 0, j, q, i, auxID, k, kk, kOK, len,  ntri, nquad, leave = 0, e4[4];
	int         *faceEdges = NULL, n1, n2;
	const int   *tris, *tric, *ptype, *pindex;
	double      angles[4], eval[18];
	const double *xyzs, *uvs;
	int    qV[6]    = { 0, 1, 2, 5, 0, 1};
	int    qLoop[5] = { 0, 1, 2, 3, 0};
	bodydata->qm = (meshMap**) EG_alloc(bodydata->nfaces*sizeof(meshMap*));
	if (bodydata->qm == NULL ) {
		return  EGADS_MALLOC;
	}
	faceEdges = (int *)EG_alloc ( bodydata -> nedges * sizeof ( int));
	if ( faceEdges == NULL ) return EGADS_MALLOC;
	for ( f = 0 ; f < bodydata->nfaces; ++f) {
		bodydata -> qm[f] = (meshMap*) EG_alloc(sizeof(meshMap));
		bodydata -> qm[f] -> mesh       = NULL;
		bodydata -> qm[f] -> backupMesh = NULL;
		bodydata -> qm[f] -> bestMesh   = NULL;
		bodydata->qm[f] -> fID = f + 1;
		if (bodydata->qm[f] == NULL ) return EGADS_MALLOC;
		// Edges associated to face //
		stat = EG_getTessFace(bodydata->tess, f + 1, &len,
				&xyzs, &uvs, &ptype, &pindex, &ntri,
				&tris, &tric);
		for ( i = 0 ; i < bodydata -> nedges; i++ ) faceEdges[i] = 0;
		for ( i = 0 ; i < len; i++ ) {
			if ( pindex[i] == -1 ) continue;
			faceEdges[pindex[i] -1]++;
		}
		for ( j = i = 0 ; i < bodydata -> nedges; i++ )
			if ( faceEdges[i] > 0 ) {
				if ( j < 4 ) e4[j] = faceEdges[i];
				j++;
			}
		printf(" Face %d is bounded by the following edges (TOTAL %d ) \n ", f + 1, j);
		nquad                        = (int)ntri/2;
		bodydata->qm[f] -> optNquads = nquad - nquad/3;
		printf(" Based on mesh %d optimal \n ", bodydata->qm[f] -> optNquads);
		if ( j == 4 ) {
			n1 = e4[0];
			n2 = e4[1];
			if ( e4[2] > n1 ) n1 = e4[2];
			if ( e4[3] > n2 ) n2 = e4[3];
			fprintf(stderr, " N0 = %d  %d and N1 %d %d \n", e4[0], e4[1], e4[2], e4[3] );
			bodydata->qm[f] -> optNquads = n1 * n2;
			printf(" 4 bounds take optimal %d x %d \n ", n1, n2);
		}
		printf(" Initial mesh has %d quads and %d vertices OPTI EST %d  \n ", nquad, len, bodydata->qm[f] -> optNquads);
		bodydata->qm[f] -> extraQuads = nquad -  bodydata->qm[f] -> optNquads;
		printf(" ESTIMATED OPTIMAL NUMBER OF QUADS IS %d OVERHEAD %d \n ",
				bodydata->qm[f] -> optNquads,  nquad - bodydata->qm[f] -> optNquads);
		stat  = EG_allocMeshData ( &bodydata ->qm[f] -> mesh      ,2 * nquad,  2 * len);
		stat += EG_allocMeshData ( &bodydata ->qm[f] -> backupMesh,2 * nquad,  2 * len);
		stat += EG_allocMeshData ( &bodydata ->qm[f] -> bestMesh  ,2 * nquad,  2 * len);
		if ( stat != EGADS_SUCCESS) {
			fprintf(stderr,"In createMeshMap  EG_allocMeshData = %d\n ", stat );
			return stat;
		}
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
					if((bodydata->qm[f] -> mesh  -> quadIdx[4*j + qLoop[kk    ]] == bodydata->qm[f] -> mesh  -> quadIdx[4*q + qLoop[k    ]] ||
							bodydata->qm[f] -> mesh  -> quadIdx[4*j + qLoop[kk    ]] == bodydata->qm[f] -> mesh  -> quadIdx[4*q + qLoop[k + 1]]
					)&&(bodydata->qm[f] -> mesh  -> quadIdx[4*j + qLoop[kk + 1]] == bodydata->qm[f] -> mesh  -> quadIdx[4*q + qLoop[k    ]] ||
							bodydata->qm[f] -> mesh  -> quadIdx[4*j + qLoop[kk + 1]] == bodydata->qm[f] -> mesh  -> quadIdx[4*q + qLoop[k + 1]]) )
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
		bodydata->qm[f]-> MINANGLE = MINVALIDANGLE;
		bodydata->qm[f]-> MAXANGLE = MAXVALIDANGLE;
		PRINTAREA = 1;
		for ( j = 0 ; j < nquad; ++j) {
			for ( k = 0 ; k < 4; ++k) {
				if ( bodydata->qm[f] -> mesh  -> quadAdj[4*j +k] > j + 1  || bodydata->qm[f] -> mesh  -> quadAdj[4*j +k] == -1) {
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
		for ( j = 0 ; j < len; j++ ) {
			printf(" looking at vertex: %d type %d ---> Might become edge \n ", j + 1, bodydata->qm[f] -> mesh ->vType[j] );
			if ( bodydata->qm[f] -> mesh ->vType[j] == -1 ) continue;
			angles[0] = EG_angleAtBoundaryVertex (bodydata->qm[f] , j + 1);
			//angles[1] = EG_angleAtVnormalPlane (bodydata->qm[f] -> mesh , j + 1, links[1], links[0] );
			//if (angles[0] > 2.0 * M_PI ) angles[0] = angles[1];
			printf(" TOTAL ANGLE SUM %lf  \n", angles[0]);
			if ( angles[0]  >  0.65 * PI  ) {
				if ( bodydata->qm[f] -> mesh  -> vType  [j] == 0 ) {
					if (angles[0] <= PI * 1.25  )
						bodydata->qm[f] -> mesh  -> vType  [j] = bodydata->qm[f] -> mesh-> sizeQuads;
					else
						bodydata->qm[f] -> mesh  -> vType  [j] = 2 * bodydata->qm[f] -> mesh-> sizeQuads;
					fprintf(stderr," Vertex %d in Face %d is now edge TYPE %d ( MAX %d )\n ",
							j + 1, f + 1,   bodydata->qm[f] -> mesh  -> vType  [j], bodydata->qm[f] -> mesh-> sizeQuads);
					printf("\n\n Vertex %d in Face %d is now edge TYPE %d ( MAX %d )\n ",
							j + 1, f + 1,   bodydata->qm[f] -> mesh  -> vType  [j], bodydata->qm[f] -> mesh-> sizeQuads);
				}
				else if ( bodydata->qm[f] -> mesh  -> vType  [j] > 0 && bodydata->qm[f] -> mesh  -> valence[j][0] == 2 )  {
					fprintf(stderr," Unsuitable for regularization. refine mesh and try again \n");
					EG_free(faceEdges);
					return FATAL_ERROR;
				}
			}
		}
		if ( stat != EGADS_SUCCESS ) stat = EGADS_SUCCESS;
		stat = checkMesh (bodydata->qm[f]);
		printf(" MESH OK %d \n",stat );
		snprintf(buffer,500,"ORIFACE_%i",f+1);
		printMesh(bodydata->qm[f] , buffer, 1);
		snprintf(buffer,500,"wvsInit_%i.txt",f+1);
		EG_wvsData(bodydata->qm[f] -> mesh , buffer);
		if (leave == 1 ) return FATAL_ERROR;
	}
	PRINTAREA = 0;
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

static int EG_averageCoords(meshMap *qm, int vID ) {
	int i, ii, count, auxID, stat, bds, vs[3], links[2], it = 0, ib;
	double uvc[2], uvl[2], eval[18], eval2[18],  *length = NULL, totLength = 0.0, norm[3], dt, boundangle = 0.0, dot[2], speed;
	vStar *star = NULL;
	uvc[0]      = 0.0;
	uvc[1]      = 0.0;
	if ( qm -> mesh -> vType[ vID -1] != -1 ) return EGADS_SUCCESS;
	stat   = EG_buildStar ( qm -> mesh, &star,vID);
	length = (double * ) EG_alloc ( star -> nQ * sizeof (double ) ) ;
	if ( stat != EGADS_SUCCESS || star == NULL || length == NULL ) return stat;
	for ( bds  = i = 0; i < star -> nQ; i++) {
		auxID    = star -> verts[2 * i + 1] - 1;
		if ( qm -> mesh -> vType[ auxID] != -1  ) bds++;
		stat = EG_segment ( qm, vID, auxID + 1, ISVERTEX, &length[i]);
		if ( stat != EGADS_SUCCESS ) {
			printf(" EG_averageCoords:: EG_segment at %d %d is %d\n ", vID, auxID + 1, stat);
			EG_free ( length);
			EG_freeStar(&star);
			return stat;
		}
		totLength += length[i];
		//printf(" POINT %d is type %d val %d Distance to vert %lf \n ", auxID +1, qm -> mesh -> vType[auxID], qm -> mesh -> valence[auxID][0] , length[i]);
		//printf("%lf %lf %lf %d %lf %lf\n ",
		//     qm -> mesh -> xyzs[ 3 * auxID], qm -> mesh -> xyzs[ 3 * auxID + 1],qm -> mesh -> xyzs[ 3 * auxID + 2],
		//   auxID + 1, qm -> mesh -> uvs [ 2 * auxID], qm -> mesh -> uvs [ 2 * auxID]);
	}
	uvc[0] = 0.0; uvc[1]  = 0.0;
	for ( i = 0; i < star -> nQ; i++) {
		auxID   = star -> verts[2 * i + 1] - 1;
		uvc[0] += (length[i] / totLength ) * qm -> mesh -> uvs    [2*auxID    ];
		uvc[1] += (length[i] / totLength ) * qm -> mesh -> uvs    [2*auxID + 1];
		//printf("%lf %lf %d %lf %lf %lf  ADDED WEIGHT %lf  FROM \n ", qm -> mesh -> uvs    [2*auxID    ], qm -> mesh -> uvs    [2*auxID  +1  ], auxID + 1,
		//     qm -> mesh -> xyzs[ 3 * auxID],qm -> mesh -> xyzs[ 3 * auxID + 1],qm -> mesh -> xyzs[ 3 * auxID + 2], length[i] / totLength);
	}
	EG_freeStar ( &star ) ;
	EG_free (length);
	stat    = EG_evaluate ( qm -> face, uvc, eval);
	if ( stat != EGADS_SUCCESS ) return stat;
	for ( i = 0 ; i < 3; i++ ) {
		if ( i < 2 )
			qm -> mesh -> uvs [ 2 * ( vID -1 ) + i] = uvc[i];
		qm -> mesh   -> xyzs[ 3 * ( vID -1 ) + i] = eval[i];
	}
	//printf("%lf %lf %d %lf %lf %lf NEW \n ",uvc[0], uvc[1], vID , eval[0], eval[1], eval[2] );
	return EGADS_SUCCESS;
	if ( bds != 2 ) {
		for ( ii = i = 0; i < star -> nQ; i++) {
			auxID    = star -> verts[2 * i + 1];
			if ( qm -> mesh -> vType[ auxID -1 ] != -1 )
				boundangle =  EG_angleAtBoundaryVertex(qm, auxID );
			if ( boundangle > PI ) {
				ib = i;
				ii = auxID;
				break;
			}
		}
		if ( ii > 0 && qm -> mesh -> valence[ii -1][0] == 3 ) {
			links[0] = -1; links[1] = -1;
			for ( count = i = 0 ;i < qm -> mesh -> valence[ii -1][0]; i++ ) {
				auxID = qm -> mesh -> valence[ii -1][2 + i];
				if ( qm -> mesh -> vType[auxID] == -1 ) continue;
				links[count++] = auxID;
				if ( count > 2 ) break;
			}
			//	  printf(" BOUND ANGLE IS %lf \n ", boundangle );
			uvl[0]  = qm -> mesh -> uvs [ 2 * ( vID - 1)    ] -  qm -> mesh -> uvs [ 2 * ( ii - 1)    ];
			uvl[1]  = qm -> mesh -> uvs [ 2 * ( vID - 1) + 1] -  qm -> mesh -> uvs [ 2 * ( ii - 1) +1 ];
			uvc[0]  = qm -> mesh -> uvs [ 2 * ( ii - 1)    ];//  + uvl[0] * 0.1;
			uvc[1]  = qm -> mesh -> uvs [ 2 * ( ii - 1) + 1];//  + uvl[1] * 0.1;
			norm[0] =  sqrt ( uvl[0] * uvl[0] + uvl[1] * uvl[1]);
			//	  printf(" NORM0 %lf \n", norm[0] );
			norm[1] = 0.0;
			if ( norm[0] < EPS04 ) {
				if ( boundangle - PI > PI * 5 / 180.0 ) {
					uvl[0]  = qm -> mesh -> uvs [ 2 * ( links[0] - 1)     ] - qm -> mesh -> uvs [ 2 * ( ii - 1)    ];
					uvl[1]  = qm -> mesh -> uvs [ 2 * ( links[0] - 1)  + 1] - qm -> mesh -> uvs [ 2 * ( ii - 1) + 1];
					uvl[2]  = qm -> mesh -> uvs [ 2 * ( links[1] - 1)     ] - qm -> mesh -> uvs [ 2 * ( ii - 1)    ];
					uvl[3]  = qm -> mesh -> uvs [ 2 * ( links[1] - 1)  + 1] - qm -> mesh -> uvs [ 2 * ( ii - 1) + 1];
					uvl[0]  = -uvl[0] - uvl[2];
					uvl[1]  = -uvl[1] - uvl[3];
					norm[1] = sqrt (uvl[0] * uvl[0] + uvl[1] * uvl[1]);
					printf(" TAKE PERPENDICULAR DIRECTION USING VERTEX %d %d %lf %lf \n ", links[0], ii,uvl[0] , uvl[1] );
					if ( norm[1] < EPS04 ) boundangle = PI;
					printf(" NORM %lf SO %f\n", norm[1], boundangle );
				}
				if ( norm[1] < EPS04 ) {
					uvl[2] = qm -> mesh -> uvs [ 2 * ( links[0] - 1)     ] - qm -> mesh -> uvs [ 2 * ( ii - 1)    ];
					uvl[3] = qm -> mesh -> uvs [ 2 * ( links[0] - 1)  + 1] - qm -> mesh -> uvs [ 2 * ( ii - 1) + 1];
					dot[0] = uvc[0] * uvl[2] + uvc[1] + uvl[3];
					printf(" DIR %d %d -> %lf %lf \n ", ii, links[0], uvl[2], uvl[3] );
					if ( dot[0] < 0.0 ) {
						printf(" DIR %d %d -> %lf %lf \n ", ii, links[1], uvl[2], uvl[3] );
						uvl[2] = qm -> mesh -> uvs [ 2 * ( links[1] - 1)     ] - qm -> mesh -> uvs [ 2 * ( ii - 1)    ];
						uvl[3] = qm -> mesh -> uvs [ 2 * ( links[1] - 1)  + 1] - qm -> mesh -> uvs [ 2 * ( ii - 1) + 1];
						dot[0] = uvc[0] * uvl[2] + uvc[1] * uvl[3];
						if (  dot[0] < 0.0 ) {
							printf(" IMPOSSIBLE!!!\n ");
							exit (1);
						}
					}
					uvl[0] = -uvl[3];
					uvl[1] =  uvl[2];
				}
				dot[1] =  qm -> mesh -> uvs [ 2 * ( vID -1 )] * uvl[0] + qm -> mesh -> uvs [ 2 * ( vID -1 ) + 1] * uvl[1];
				if ( dot[1] < 0.0 ) {
					uvl[0] *= -1.0;
					uvl[1] *= -1.0;
				}
			}
			printf(" LINK %d %d -> vector director %lf %lf \n ", links[0], ii, uvl[0], uvl[1] );
			stat      = EG_evaluate (qm -> face, uvc, eval );
			totLength = 0.0;
			for ( i = 0 ; i < star -> nQ; i++) {
				if ( i == ib ) continue;
				totLength += length[i];
			}
			totLength /= (double)( star -> nQ - 1);
			speed      = 0.01 * totLength;
			printf(" AT VERTEX %d (ii %d) -> MOVING WITH %lf %lf \n ", vID, ii, uvl[0], uvl[1] );
			printf("%lf %lf %lf %d \n ", eval[0], eval[1], eval[2], vID );
			printf(" Targt LENGTH %lf local %lf MOVING WITH SPEED %lf \n ", totLength, norm[0], speed );
			norm[1] = 0.0;
			while ( fabs ( totLength - norm[1] ) > totLength * 0.1 && it < 100  ) {
				printf("\n\n NORM %lf length %lf diff %lf target %lf\n ", norm[1], totLength,  norm[1]- totLength, totLength * 0.1);
				it++;
				uvc[0]      += speed * uvl[0];
				uvc[1]      += speed * uvl[1];
				stat         = EG_evaluate ( qm -> face, uvc, eval2 );
				norm[2]      = 0.0;
				for ( i      = 0; i < 3; i++) {
					norm[2] += (eval2[i] - eval[i]) * (eval2[i] - eval[i]) ;
					eval[i]  =  eval2[i];
				}
				norm[2]  = sqrt ( norm[2] );
				norm[1] += norm[2];
				if ( norm[1] > totLength  ) break;
				printf(" #POINT NOW HAS %lf  TARGET %lf LOCAL %lf dt %lf \n ", norm[1], totLength, norm[2], dt );
				//speed = 0.5 * ( fabs (norm[1] - totLength) );
				printf("%lf %lf %lf %d %lf %lf \n ", eval[0], eval[1], eval[2], vID, uvc[0], uvc[1] );
			}
		}
		else {
			uvc[0] = 0.0; uvc[1]  = 0.0;
			for ( i = 0; i < star -> nQ; i++) {
				auxID   = star -> verts[2 * i + 1] - 1;
				uvc[0] += (length[i] / totLength ) * qm -> mesh -> uvs    [2*auxID    ];
				uvc[1] += (length[i] / totLength ) * qm -> mesh -> uvs    [2*auxID + 1];
				//    printf("%lf %lf %d %lf %lf %lf  ADDED WEIGHT %lf  FROM \n ", qm -> mesh -> uvs    [2*auxID    ], qm -> mesh -> uvs    [2*auxID  +1  ], auxID + 1,
				//     qm -> mesh -> xyzs[ 3 * auxID],qm -> mesh -> xyzs[ 3 * auxID + 1],qm -> mesh -> xyzs[ 3 * auxID + 2], length[i] / totLength);
			}
		}
	} else {
		vs[0] = -1; vs[1] = -1; ii = 0;
		for ( count = i = 0 ; i < star -> nQ; i++ ) {
			vs[0] = star -> verts [2 * i + 1 ] - 1;
			if ( qm -> mesh -> vType[ vs[0] ] == -1 ) continue;
			for ( ii = i + 1 ; ii < star -> nQ; ii++ ) {
				vs[1] = star -> verts [2 * i + 1 ] - 1;
				if ( qm -> mesh -> vType[ vs[1] ] != -1 ) {
					count = 1;
					break;
				}
			}
			if ( count == 1 ) break;
		}
		if ( ii - i == 1 ) {
			for ( count = i =  0 ; i < qm -> mesh ->valence[vID - 1][0]; i++) {
				auxID   = qm -> mesh -> valence[vID - 1][2 + i] - 1;
				count++;
				uvc[0] += qm -> mesh -> uvs    [2*auxID    ];
				uvc[1] += qm -> mesh -> uvs    [2*auxID + 1];
			}
		}
		else {
			uvc[0] = qm -> mesh -> uvs  [2*vs[0]     ] + qm -> mesh -> uvs  [2*vs[1]    ];
			uvc[1] = qm -> mesh -> uvs  [2*vs[0] + 1 ] + qm -> mesh -> uvs  [2*vs[1] + 1];
			count  = 2;
		}
		uvc[0] /= (double)count;
		uvc[1] /= (double)count;
	}
	EG_freeStar ( &star ) ;
	EG_free (length);
	stat    = EG_evaluate ( qm -> face, uvc, eval);
	if ( stat != EGADS_SUCCESS ) return stat;
	for ( i = 0 ; i < 3; i++ ) {
		if ( i < 2 )
			qm -> mesh -> uvs [ 2 * ( vID -1 ) + i] = uvc[i];
		qm -> mesh   -> xyzs[ 3 * ( vID -1 ) + i] = eval[i];
	}
	printf("%lf %lf %d %lf %lf %lf NEW \n ",uvc[0], uvc[1], vID , eval[0], eval[1], eval[2] );
	return EGADS_SUCCESS;
}

static int
EG_adjQtoPair(meshData *mesh, int qID, int v1, int v2, int *adj) {
	int i, aux = -1;
	adj[0] = -1; adj[1] = -1;
#ifdef DEBUGG
	printf(" Looking for adjacent to vertex pair %d %d in quad %d\n", v1,v2, qID );
#endif
	for ( i = 0 ; i < 4; ++i) {
		if ( mesh -> quadIdx[4*(qID - 1) + i ] == v1 ) adj[0] = i;
		if ( mesh -> quadIdx[4*(qID - 1) + i ] == v2 ) aux    = i;
		if ( aux != -1 && adj[0] != -1 ) break;
	}
	if ( aux == -1 || adj[0] == -1 ) return EGADS_SUCCESS;
	if      ( abs (adj[0] - aux ) == 3 ) adj[0] = 3;
	else if ( aux < adj[0]             ) adj[0] = aux;
	adj[1] = mesh -> quadAdj[4*(qID - 1) + adj[ 0 ] ];
#ifdef DEBUGG
	printf(" Adjacent is at position %d and is quad %d \n", adj[0], adj[1] );
#endif
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
			printf(" I HAVE %d REMOVED VERTICES BUT ACTUALLY %d REMOVED QUADS!!!!!\n ", vRem, qRem);
			return FATAL_ERROR;
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

static int checkInvalidElement(meshMap *qm, int qID, double minAngle, double maxAngle ) {
	int    j, k, vaux, val, stat, area, it = 0, itMAX = 20, v[4];
	double angles[4], quadArea;
	k = checkQuad ( qm -> mesh, qID );
	if ( k != EGADS_SUCCESS ) return k;
	PRINTAREA = 0;
	for ( k = 0; k < 4; ++k) v[k] = qm -> mesh -> quadIdx[4 * ( qID - 1 ) + k];
	stat = quadAlgebraicArea(qm, minAngle, maxAngle, qID, v[0], &area, angles, &quadArea);
	if ( stat != EGADS_SUCCESS || area == 1 )  {
		if ( stat != EGADS_SUCCESS) printf(" Stat in area %d --> = %d\n ", area, stat );
		return stat;
	}
	while ( it < itMAX )
	{
		for  ( k = 0; k < 4; k++) {
			stat = EG_averageCoords(qm, v[k]);
			if ( stat != EGADS_SUCCESS) {
				printf(" Stat EG_averageCoords %d\n ", stat );
				return stat;
			}
			val = qm -> mesh -> valence [v[k] - 1] [0];
			for (j = 0; j < val; j++) {
				vaux = qm -> mesh -> valence [ v[k] - 1] [ 2 + j];
				if ( vaux == -1 || qm -> mesh -> vType[vaux -1] != -1 ) continue;
				stat = EG_averageCoords(qm , vaux );
			}
			stat = EG_averageCoords(qm , v[k] );
		}
		stat = quadAlgebraicArea(qm, minAngle, maxAngle, qID, v[0], &area, angles, &quadArea ); // detecting invalid/distorted elements for optimizer
		if ( area == 1 ) break;
		++it;
	}
	return stat;
}














/*static int checkInvalidElement(meshMap *qm, int qID, double minAngle, double maxAngle ) {
  int    i, j, k, kk, v[4], piv[4], stat, area, it = 0, itMax = 20, vaux, vaux2, r;
  double angles[4], quadArea;
  i = checkQuad ( qm, qID );
  if ( i != EGADS_SUCCESS ) return i;
  for ( i = 0; i < 4; ++i) {
      v[i]   = qm -> quadIdx[4 * ( qID - 1 ) + i];
      piv[i] = i;
  }
  if ( EG_nValenceCount ( qm , qID, 0, 4 ) != 4 ) {
      minAngle = MINVALIDANGLE;
      maxAngle = MAXVALIDANGLE;
  }
  stat = quadAlgebraicArea(qm, minAngle, maxAngle, qID, v[0], &area, angles, &quadArea);
  if ( stat != EGADS_SUCCESS)  {
      printf(" Stat in area = %d\n ", stat );
      return stat;
  }
  if ( area == 1 )               return EGADS_SUCCESS;
  printf("\n\n ================== INITI CHECKINVALID Q %d =================== \n\n ", qID );
  printMesh(qm, buffer,0);
  do
    {
      printf(" CHECKINVALID QUAD R %d IT %d -> MIN MAX %lf %lf\n", r, it, minAngle, maxAngle );
      for  ( k = 0; k < 4; k++)
	while ( angles[k] > 2.0 * PI ) angles[k] -= 2.0*PI;
      for  ( k = 0; k < 3; k++) {
	  for  ( j = k + 1; j < 4; j++) {
	      if ( angles[ piv [j] ] > angles[ piv[k ] ]) {
		  i      = piv[j];
		  piv[j] = piv[k];
		  piv[k] = i;
	      }
	  }
      }
      if ( angles[ piv[0]] < maxAngle) return EGADS_SUCCESS;
      for  ( i = 0; i < 4; i++) {
	  if ( angles [ piv[i] ] > minAngle && angles[piv[i] ] < maxAngle ) continue;
	  printf("\n %lf %lf = %d before IN FIRST LOOP \n", qm -> uvs[ 2 * ( v[piv[i]] -1 ) ], qm -> uvs[ 2 * ( v[piv[i]] -1 ) + 1],  v[piv[i]] );
	  stat = EG_averageCoords(qm, v[ piv[i]]);
	  printf("%lf %lf %d after IN FIRST LOOP \n\n",  qm -> uvs[ 2 * ( v[piv[i]] -1 ) ], qm -> uvs[ 2 * ( v[piv[i]] -1 ) + 1],  v[piv[i]]);
	  stat = quadAlgebraicArea(qm, minAngle, maxAngle, qID, v[piv[i]], &area, angles, &quadArea ); // detecting invalid/distorted elements for optimizer
	  if ( area == 1 ) return stat;
      }
      for  ( i = 0; i < 4; i++) {
	  for (j = 0; j < qm -> valence[ v[ piv[ i] ] - 1 ][0]; ++j) {
	      vaux = qm -> valence[ v[ piv[ i] ] - 1 ][2 + j];
	      printf("\n\n Moving Linkt %d from vertex %d (i - %d\n", vaux, v[piv[i]], i );
	      printf(" %lf %lf = %d before IN second LOOP \n", qm -> uvs[ 2 * ( vaux -1 ) ], qm -> uvs[ 2 * ( vaux -1 ) + 1], vaux);
	      EG_averageCoords(qm, vaux);
	      printf(" %lf %lf %d after IN second LOOP  \n\n\n",qm -> uvs[ 2 * ( vaux -1 ) ], qm -> uvs[ 2 * ( vaux -1 ) + 1], vaux);
	      stat = quadAlgebraicArea(qm, minAngle, maxAngle, qID, v[ piv[ i] ], &area, angles, &quadArea ); // detecting invalid/distorted elements for optimizer
	      if ( area == 1 ) return stat;
	  }
      }
      printMesh(qm,buffer,0);
      ++it;
    }
  while ( it < itMax) ;
  it = 0 ;
  printMesh(qm, buffer,0);
  PRINTAREA =1;
  stat = quadAlgebraicArea(qm, minAngle, maxAngle, qID, v[0], &area, angles, &quadArea );
  PRINTAREA =0;
  return EGADS_GEOMERR;
}
 */
static int
EG_projectToTangentPlane(double normal[], double *O, double *p, double *proj) {
	double c, dotNN = 0.0, dotNP = 0.0, dist, lambda;
	dist    = (p[0]- O[0]) * (p[0]- O[0]) + (p[1]- O[1])*(p[1]- O[1]) + (p[2]- O[2])*(p[2]- O[2]);
	dist    = sqrt(dist);
	if (dist < EPS14 ) {
		proj[0] = p[0]; proj[1] = p[1]; proj[2] = p[2];
		return EGADS_SUCCESS;
	}
	c       = normal[0] *      O[0] + normal[1] *      O[1] + normal[2] *      O[2]; // Equation plane: a*x + b*y + c*z = C
	dotNN   = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
	dotNP   = normal[0] *      p[0] + normal[1] *      p[1] + normal[2] *      p[2];
	if ( fabs(dotNP - c) <= EPS14) {
		proj[0] = p[0]; proj[1] = p[1]; proj[2] = p[2];
		return EGADS_SUCCESS;
	}
	lambda  =   (c - dotNP) / dotNN;
	proj[0] = p[0] + lambda * normal[0];
	proj[1] = p[1] + lambda * normal[1];
	proj[2] = p[2] + lambda * normal[2];
	// check that point belongs to plane.
	dist  = normal[0] * proj[0] + normal[1] * proj[1] +  normal[2] * proj[2];
	if( fabs(dist - c) < EPS14) {
		return EGADS_SUCCESS;
	}
	else{
		printf(" ORIGIN %lf %lf  %lf  NORMAL %lf %lf  %lf  TARGET %lf %lf %lf\n", O[0], O[1], O[2], normal[0], normal[1], normal[2], p[0], p[1], p[2]);
		printf(" POINT SHOULD BELONG TO PLANE!!!!! %lf ~= %lf\n",dist,c);
		printf(" DOT NN %lf PN %lf LAMBDA %lf  POINT %lf %lf %lf\n", dotNN, dotNP, lambda, proj[0], proj[1], proj[2]);
		return FATAL_ERROR;
	}
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





static double
EG_angleBetweenVectors ( double *normal, double *v1, double *v2) {
	double cross[4], dotV1V2, dotCN, norm, theta;
	EG_unitVector (v1, &norm );
	EG_unitVector (v2, &norm );
	EG_unitVector (normal, &norm);
	cross_product (v1, v2, cross );
	dotV1V2 = dotProduct ( v1, v2 );
	dotCN   = dotProduct ( normal, cross);
	if      ( dotV1V2 >=  1.0 ) theta = 0.0;
	else if ( dotV1V2 <= -1.0 ) theta = PI;
	else                        theta = acos ( dotV1V2 );
	if ( dotCN < 0 ) theta = 2.0*PI - theta;
	// printf(" Angle Is %lf ( wrp normal %lf ) DOT P %lf \n ", theta, dotCN, dotV1V2 );
	//printf(" V1 %lf %lf %lf \n ", v1[0], v1[1], v1[2]  );
	//printf(" V2 %lf %lf %lf \n ", v2[0], v2[1], v2[2]  );
	return theta;
}


static double
EG_angleAtVnormalPlane ( meshMap *qm, int vC, int v1, int v2 ) {
	int i, stat;
	double normal[4], v1EPS[5], v2EPS[5], eval[18], vec1[4], vec2[4], a1, a2 ;
	//printf(" ANGLE %d,%d - %d,%d\n ", v1,vC, v2,vC );
	stat  = vectorAtVertexNplane (qm , vC, v1, v1EPS);
	stat += vectorAtVertexNplane (qm , vC, v2, v2EPS);
	stat += EG_evaluate ( qm -> face, &qm -> mesh -> uvs [ 2 * ( vC - 1) ], eval );
	if ( stat != EGADS_SUCCESS ) return stat;
	//
	if ( qm -> face -> mtype == SREVERSE )
		cross_product(&eval[6], &eval[3], normal);
	else
		cross_product(&eval[3], &eval[6], normal);
	for ( i = 0 ; i < 3; i++ ) {
		vec1[i] = v1EPS [ 2 + i ] - qm -> mesh -> xyzs[ 3 * ( vC - 1 ) + i ];
		vec2[i] = v2EPS [ 2 + i ] - qm -> mesh -> xyzs[ 3 * ( vC - 1 ) + i ];
	}
	EG_unitVector (vec1, &vec1[3] );
	EG_unitVector (vec2, &vec2[3] );
	/*printf("%lf %lf %lf %d \n ",
	 qm -> xyzs[ 3 * ( vC - 1 )], qm -> xyzs[ 3 * ( vC - 1 ) + 1], qm -> xyzs[ 3 * ( vC - 1 ) + 2], vC);
  printf("%lf %lf %lf %d \n\n\n ",
	 qm -> xyzs[ 3 * ( vC - 1 )    ] + vec1[0],
	 qm -> xyzs[ 3 * ( vC - 1 ) + 1] + vec1[1],
	 qm -> xyzs[ 3 * ( vC - 1 ) + 2] + vec1[2], v1);
  printf("%lf %lf %lf %d \n ",
	 qm -> xyzs[ 3 * ( vC - 1 )], qm -> xyzs[ 3 * ( vC - 1 ) + 1], qm -> xyzs[ 3 * ( vC - 1 ) + 2], vC);
  printf("%lf %lf %lf %d \n\n\n ",
	 qm -> xyzs[ 3 * ( vC - 1 )    ] + vec2[0],
	 qm -> xyzs[ 3 * ( vC - 1 ) + 1] + vec2[1],
	 qm -> xyzs[ 3 * ( vC - 1 ) + 2] + vec2[2], v2);
	 */
	a1 =  EG_angleBetweenVectors (normal, vec1, vec2 );
	for ( i = 0 ; i < 3; i++ ) {
		vec1[i] = qm -> mesh -> xyzs[ 3 * ( v1 - 1 ) + i ] - qm -> mesh -> xyzs[ 3 * ( vC - 1 ) + i ];
		vec2[i] = qm -> mesh -> xyzs[ 3 * ( v2 - 1 ) + i ] - qm -> mesh -> xyzs[ 3 * ( vC - 1 ) + i ];
	}
	EG_unitVector (vec1, &vec1[3] );
	EG_unitVector (vec2, &vec2[3] );
	/*printf("%lf %lf %lf %d \n ",
	 qm -> xyzs[ 3 * ( vC - 1 )], qm -> xyzs[ 3 * ( vC - 1 ) + 1], qm -> xyzs[ 3 * ( vC - 1 ) + 2], vC);
  printf("%lf %lf %lf %d \n\n\n ",
	 qm -> xyzs[ 3 * ( vC - 1 )    ] + vec1[0],
	 qm -> xyzs[ 3 * ( vC - 1 ) + 1] + vec1[1],
	 qm -> xyzs[ 3 * ( vC - 1 ) + 2] + vec1[2], v1);
  printf("%lf %lf %lf %d \n ",
	 qm -> xyzs[ 3 * ( vC - 1 )], qm -> xyzs[ 3 * ( vC - 1 ) + 1], qm -> xyzs[ 3 * ( vC - 1 ) + 2], vC);
  printf("%lf %lf %lf %d \n\n\n ",
	 qm -> xyzs[ 3 * ( vC - 1 )    ] + vec2[0],
	 qm -> xyzs[ 3 * ( vC - 1 ) + 1] + vec2[1],
	 qm -> xyzs[ 3 * ( vC - 1 ) + 2] + vec2[2], v2);*/
	a2 =  EG_angleBetweenVectors (normal, vec1, vec2 );
	//printf(" ANGLE AT NODES %lf AT PLANE %lf \n ", a2, a1 );
	return a2;
}



static double EG_angleAtBoundaryVertex(meshMap *qm, int v ) {
	int stat, i, j, k, links[2];
	double angle;
	vStar *star = NULL;
	if ( qm -> mesh -> vType[v -1 ] < 0 ) return (double)FATAL_ERROR;
	stat = EG_buildStar ( qm -> mesh, &star, v);
	if ( stat != EGADS_SUCCESS || star == NULL ) {
		printf(" Looking at corners: buildstar %d is %d \n ", v, stat );
		return stat;
	}
	links[0] = -1; links[1] = -1;
	for ( i = j = 0 ;j < star -> nQ; j++ ) {
		k = star -> verts[ 2 * j + 1 ] - 1 ;
		if ( star -> quads[j] == -1  ) {
			links[1] = star -> verts [ 2 * j + 1];
			links[0] = star -> verts [star -> idxV[ 2 * j + 3]];
			//  printf(" star %d -> links %d %d \n ", v, links[0], links[1] );
		}
		if ( qm -> mesh -> vType[k] == -2 ) fprintf(stderr,"WHAT THE FUCK!! \n ");
		if ( qm -> mesh -> vType[k] != -1 )
			i++;
	}
	if ( i > 2) {
		EG_freeStar(&star);
		return (double)PI * 0.5;
	}
	angle = EG_angleAtVnormalPlane (qm, v , links[0], links[1] );
	EG_freeStar(&star);
	return angle;
}


static double
triang_cross_area(double *A, double *B, double *C, double *cross){
	int i;
	double AB[3], AC[3], norm;
	for ( i = 0 ; i < 3; ++i){
		AB[i] = B[i] - A[i];
		AC[i] = C[i] - A[i];
	}
	cross_product(AB, AC, cross);
	norm = cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2];
	norm = 0.5*sqrt(norm);
	return (norm);
}

static int
vectorAtVertexNplane (meshMap *qm, int v0, int v1, double *uvEPS ) {
	int stat, i;
	double evalVert[18], v0_normal[4], proj[3], nEPS, dir[2];
	stat = EG_evaluate(qm -> face, &qm -> mesh -> uvs[2 * ( v0 - 1 ) ], evalVert);
	if ( stat != EGADS_SUCCESS) return stat;
	//
	if ( qm -> face -> mtype == SREVERSE )
		cross_product(&evalVert[6], &evalVert[3], v0_normal);
	else
		cross_product(&evalVert[3], &evalVert[6], v0_normal);
	//sampleNormalPlane(qm, v0_normal, &qm -> xyzs[3 * ( v0 - 1)], v0);
	EG_unitVector(v0_normal, &v0_normal[3]);
	dir[0] = qm -> mesh -> uvs[ 2 * ( v1 - 1 )      ] - qm -> mesh -> uvs[ 2 * ( v0 - 1 )      ];
	dir[1] = qm -> mesh -> uvs[ 2 * ( v1 - 1 )  + 1 ] - qm -> mesh -> uvs[ 2 * ( v0 - 1 )  + 1 ];
	nEPS   = sqrt ( dir[0] * dir[0] + dir[1] * dir[1] ) ;
	if ( nEPS < EPS08 ) {
		printf(" VErtices %d and %d are %.15e distance apart in parametric space!!! \n ", v0, v1, nEPS);
		return FATAL_ERROR;
	}
	for ( i = 0 ; i < 2; i++ ) {
		//dir[i]  /= nEPS;
		uvEPS[i] = qm -> mesh -> uvs[ 2 * ( v0 - 1 ) + i];
	}
	do {
		uvEPS[0]  += dir[0]/10.0;// nEPS/10.0 * dir[0];
		uvEPS[1]  += dir[1]/10.0;//nEPS/10.0 * dir[1];
		stat       = EG_evaluate(qm -> face, uvEPS, evalVert);
		if( stat  != EGADS_SUCCESS ) return stat;
		stat       = EG_projectToTangentPlane(v0_normal, &qm -> mesh -> xyzs [ 3 * ( v0 - 1 ) ], evalVert, proj);
		if( stat  != EGADS_SUCCESS ) return stat;
		nEPS       = 0.0;
		for ( i = 0 ; i < 3; ++i)
			nEPS += (proj[i] - qm -> mesh -> xyzs[ 3 * ( v0 - 1) + i ] ) * ( proj[i] - qm -> mesh -> xyzs[ 3 * ( v0 - 1) + i ]);
	} while(sqrt ( nEPS ) < 1.E-02);
	for ( i = 0 ; i < 3 ; i++ )  uvEPS [ 2 + i ] = proj [ i ];
	return EGADS_SUCCESS;
}

static int
quadAlgebraicArea(meshMap *qm, double minTheta, double maxTheta, int qID, int vID,  int *validArea, double *quadAngles, double *quadArea ) {
	int i, id, vPos, qV[4], k, tri, piv[4], stat, degenerated[4];
	double evalVert[18], cross[4], dotCross[2*4], area[2*4], theta[2*4], vNormal[4], uv_eps[10], vEps[12], theta2;
	double dot;
	vPos = -1;
	stat = checkQuad ( qm -> mesh, qID );
	if ( stat != EGADS_SUCCESS ) {
		printf(" In quadAlgebraicArea:: checking area of quad %d  !!! %d \n", qID, stat );
		return stat;
	}
	if ( PRINTAREA == 1 ) printQuadSpecs (qm -> mesh, qID );
	for ( i = 0 ; i < 4; ++i) {
		qV[i] = qm -> mesh ->quadIdx[4*(qID - 1) + i] - 1;
		if (qV[i] == (vID -1) )
			vPos = i;
	}
	if ( vPos == -1 ) {
		printf(" VERTEX %d is STAR CENTRE AND SHOULD BELONG TO QUAD %d\n",vID, qID);
		printQuadSpecs ( qm -> mesh, qID );
		return EGADS_INDEXERR;
	}
	for ( k = 0 ; k < 4; ++k) {
		degenerated[k       ] = 0;
		// Compute area starting from all vertices to capture invalid shapes
		for ( tri = 0 ; tri < 2; ++ tri) { // quad area as a sum of 2 triangles
			dotCross   [2 *k      ] = 0.0;
			dotCross   [2 *k + tri] = 0.0;
			piv [0] = qV[(k  + vPos          )%4];
			piv [1] = qV[(k  + vPos + 1 + tri)%4];
			piv [2] = qV[(k  + vPos + 2 + tri)%4];
			//  printf(" CENTRE %d -> %d %d ->%d \n", piv[0] + 1, piv[1] + 1, piv[0] + 1, piv[2] + 1 );
			for ( id = 0 ; id < 2; ++id) { // move away from vertex towards B and C respectively: we are computing the angle very close to the centre
				stat = vectorAtVertexNplane (qm , piv[0] + 1, piv[id + 1] + 1, &uv_eps[ 5 * id ] );
				if ( stat != EGADS_SUCCESS) { // vector OB or OC is too small => O ~B  or O ~C INVALID
					printf(" stat in vector_at Normal plane = %d \n", stat );
					degenerated[  k    ] =    1;
					theta      [2*k    ] =  0.0;
					theta      [2*k + 1] =  0.0;
					area       [2*k    ] =  0.0;
					area       [2*k + 1] =  0.0;
					break;
				}
			}
			if ( degenerated[k] == 1 ) continue;
			for ( i = 0 ; i < 3; ++i) {
				vEps[    i] = uv_eps[    2 + i] - qm -> mesh -> xyzs [ 3 * piv[0] + i ] ;
				vEps[4 + i] = uv_eps[5 + 2 + i] - qm -> mesh -> xyzs [ 3 * piv[0] + i ] ;
			}
			// vEps live in normal plane at piv[0]
			EG_unitVector (&vEps[0], &vEps[3] );
			EG_unitVector (&vEps[4], &vEps[7] );
			//	  printf(" DISTANCE %lf  %lf  \n ", vEps[3], vEps[7] );
			if ( vEps[3] < EPS11 || vEps[7] < EPS11 ) {
				theta[2*k + tri] = 0.0;
				dot              = 1.0;
			}
			else {
				dot = 0.0;
				for     ( i = 0 ; i < 3; ++i)       dot             += vEps[i] * vEps[i + 4];
				if      ( fabs(dot - 1.0) < EPS11 ) theta[2*k + tri] = 0.0;
				else if ( fabs(dot + 1.0) < EPS11 ) theta[2*k + tri] = PI;
				else                                theta[2*k + tri] = acos(dot);
			}
			area [ 2 * k + tri ] = triang_cross_area(&qm -> mesh -> xyzs[ 3 * piv[0] ] ,
					&qm -> mesh -> xyzs[ 3 * piv[1] ] ,
					&qm -> mesh -> xyzs[ 3 * piv[2] ] , cross);
			EG_unitVector ( cross, &cross[3] );
			stat = EG_evaluate(qm -> face, &qm -> mesh -> uvs[ 2 * piv[0] ], evalVert);
			if ( stat != EGADS_SUCCESS) return stat;
			//
			if ( qm -> face -> mtype == SREVERSE )
				cross_product(&evalVert[6], &evalVert[3], vNormal);
			else
				cross_product(&evalVert[3], &evalVert[6], vNormal);
			EG_unitVector(vNormal, &vNormal[3]);
			dotCross[ 2 * k + tri ] = cross[0] * vNormal[0] + cross[1] * vNormal[1] + cross[2] * vNormal[2];
			//  printf(" NORMAL %lf %lf %lf cross %lf %lf  %lf  DOT %lf \n ",
			//	 vNormal[0], vNormal[1], vNormal[2], cross[0], cross[1], cross[2], dotCross[ 2 * k + tri]);
			if ( dotCross[ 2 * k + tri ] < 0) {
				theta[ 2 * k + tri ]  =  2.0 * PI - theta[ 2 * k + tri];  // we have anti-clockwise orientation. Postive area implies reversed.
				area [ 2 * k + tri ] *= -1.0;
			}
		}
	}
	*validArea = 1;
	//if ( fabs( (area[0] + area[1]) - (area[2] + area[3]) )  > EPSAREA || ( area[0] * area[1] < 0 ) ) *validArea = -1;
	if ( ( area[0] * area[1] < 0 ) ) *validArea = -1;
	if ( (PRINTAREA == 1 ) && *validArea == -1)
		printf(" INVALID AREA BECAUSE AREA SUM  %lf != %lf || area * area < 0  %lf \n", area[0] + area[1] , area[2] + area[3], area[0] * area[1] );
	*quadArea = area[0] + area[1];
	for ( k = 0 ; k < 4; ++k) {
		quadAngles[k] = theta[2*k] + theta[2*k+1];
		if (( qm -> mesh -> vType[qV[(k  + vPos )%4]] == -1 ) &&
				( degenerated[k] == 1 || quadAngles[k] > maxTheta || quadAngles[k] < minTheta ) ) *validArea = - 1;
	}
	if ( (PRINTAREA == 1 ) && *validArea == -1) {
		printf(" ========== INVALID QUAD   %d  ==========  \n", qID ) ;
		dot    = 0.0;
		theta2 = EG_angleAtVnormalPlane ( qm, qV[0] + 1, qV[1] + 1, qV[3] + 1 );
		printf(" CENTRAL ANGLE %d %d %d --- > %lf \n ", qV[1] + 1, qV[0] + 1, qV[3] + 1, theta2 );
		for ( k = 0 ; k < 4; ++k) {
			dot += quadAngles[k];
			printf(" VERT %d  QUADANGLE %.15lf = %.15lf (DEG) -> ANGLE BOUNDS [%lf , %lf]  DEGENERATED %d \n",
					qV[(k  + vPos )%4] + 1, quadAngles[k], quadAngles[k] * 180.0/ PI, minTheta, maxTheta, degenerated[k]);
			printf(" angle %d  dot %lf  %lf  Theta %lf + %lf = %lf, area = %lf  +   %lf = %lf \n",k,
					dotCross[ 2 * k], dotCross[ 2 * k + 1],
					theta[ 2 * k ], theta[ 2 * k + 1], theta[ 2 * k ] + theta[ 2 * k + 1 ] ,
					area[ 2 * k ], area[ 2 * k + 1], area[ 2 * k ] + area[ 2 * k + 1 ] );

		}
		printf(" QUAD INTERNAL ANGLES SUM %.16f ~ %.16f (?)  VALIDAREA %d  \n", dot, 2*PI, *validArea);
	}
	return EGADS_SUCCESS;
}

static int
EG_buildStar(meshData *mesh, vStar **star, int vID ) {
	int  stat, i = 0 , id0 = -1, q = 0, auxV, auxQ,  v = 0, quadID, prevQuad, it = 0, it2 = 0;
	int adj[2], *vertex = NULL, *quads = NULL;
	int qLoop[8] = {0, 1, 2, 3, 0, 1, 2, 3};
	stat   = checkVertex ( mesh , vID ) ;
	if ( stat != EGADS_SUCCESS ) return stat;
	vertex = (int * ) EG_alloc ( mesh -> totVerts * sizeof ( int ) );
	quads  = (int * ) EG_alloc ( mesh -> totQuads * sizeof ( int ) );
	if ( vertex == NULL || quads == NULL ) return EGADS_MALLOC;
	// quads are -1 bias
	quadID       = mesh -> valence[ vID - 1 ][1] - 1;
	i            = checkQuad ( mesh, quadID + 1);
	if ( i != EGADS_SUCCESS){
		printf(" EG_buildStar at vertex %d has associated a bad quad: %d --> %d!!\n", vID, quadID + 1, i );
		printQuadSpecs ( mesh, quadID + 1);
		EG_free ( vertex );
		EG_free ( quads  );
		return i;
	}
	vertex [v++] = vID;
	it = 0;
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
		//   printf(" NEXT QUAD IN BUILD STAR %d is %d \n", vID, quadID + 1);
		//  printQuadSpecs ( mesh, quadID + 1 );
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
					printf(" In buildStar:: vertex Id %d in quad %d is %d\n", vID, quadID + 1, i);
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
		printf(" I failed at removing doublet \n");
		return EGADS_GEOMERR;
	}
	// Check that it has left valid elements at each old link
	for ( i = 0 ; i < 2; i++ ) {
		if ( checkVertex ( qm -> mesh, link[i] ) == EGADS_SUCCESS ) {
			stat  = EG_removeDoublet (qm, link[i] ) ;
			if ( stat != EGADS_SUCCESS ) {
				printf(" I failed at removing doublet inside doublet \n");
				printMesh(qm, buffer,0);
				return stat;
			}
		}
	}
	return EGADS_SUCCESS;
}

static int
EG_forceCollapse ( meshMap *qm, int qID, int *activity ) {
	int i, i3, centre;
	*activity  = 0;
	printQuadSpecs (qm -> mesh, qID );
	if ( EG_nValenceCount ( qm -> mesh, qID, 3 ) == 0 ) return EGADS_SUCCESS;
	for ( i = i3 = 0 ; i3 < 4; i3++ ) {
		centre = qm -> mesh -> quadIdx [ 4 * ( qID - 1 ) + i3 ];
		if ( qm -> mesh -> vType [centre -1 ] == -1 && qm -> mesh -> valence[centre -1][0] == 3 ) {
			i = 1;
			break;
		}
	}
	if ( i == 0 ) return EGADS_SUCCESS;
	printf(" collapse through %d\n ", centre);
	if (validCollapse (qm -> mesh, qID, centre) == 0) {
		printf("invalid collapse %d\n", qID );
		return EGADS_SUCCESS;
	}
	*activity = 1;
	return EG_mergeVertices ( qm, qID, centre );
}

static int
EG_collapse (meshMap *qm, int qID, int *activity )  {
	int vC, stat, ids[15], vs[15] ;
	*activity     = 0 ;
	stat          = checkQuad ( qm -> mesh , qID );
	if ( stat != EGADS_SUCCESS ) {
		printf(" EG_collapse for quad %d is %d \n ", qID, stat );
		return stat;
	}
	stat          = quadValenceDistribution (qm -> mesh, qID, ids, vs );
	if ( stat    != EGADS_SUCCESS || ids[5] > 1 || ids[10] * ids[0] == 0    ) return stat; // 3 or more vertices have valence 4
	if ( ids[10] == 1 || (ids[10] == 2 && abs ( ids[11] - ids[12] ) % 2 == 0 ) ) {
		vC = qm -> mesh -> quadIdx [ 4 * ( qID -1 ) + (ids[11] + 1 ) %4 ];
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
	links[0] = qm -> mesh -> quadIdx [ 4 * ( qC - 1) + ( i + 1 ) %4];
	links[1] = qm -> mesh -> quadIdx [ 4 * ( qC - 1) + ( i + 3 ) %4];
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
	qm ->extraQuads = j - qm -> optNquads;
	//fprintf(stderr, " MERGE REMOVED QUADS / VERTICES :: %d  %d   extraquads %d \n ", qm -> mesh -> remQuads[0], qm -> mesh -> remVerts[0], qm ->extraQuads );
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
	int stat, i0, i1, i, j, adj, idx[4], qID[2];
	int adjQmap[6];
	swap   = swap %3;
	qID[0] = qg.q[0]; qID[1] = qg.q[1];
	if ( swap == 0 ) {
		printf(" swapping throu 0-3 will result in the same pair!! \n ");
		return EGADS_INDEXERR;
	}
#ifdef DEBUG
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
	stat = optimize_angles(qm, 6, qg.verts,0);
	return stat;
}


static int
EG_forceSplit ( meshMap *qm, int v0, int vL, int *activity ) {
	int q, id0, i, i3, dist[2], vOpp, stat, links[4], id6[2];
	vStar *star = NULL;
	*activity   = 0;
	stat        = EG_buildStar (qm -> mesh , &star, v0 ) ;
	if ( stat  != EGADS_SUCCESS || star == NULL ) return stat;
	for ( i3 = i = 0 ; i < star -> nQ; i++ ) {
		if ( star -> quads[i] == -1 ) continue;
		printf(" Looking for %d i %d -> %d \n ", vL, 2 * i + 1, star -> verts[2 * i + 1]);
		if ( star -> verts[ 2 * i + 1] == vL ) {
			i3 = 2 * i + 1;
			break;
		}
	}
	if ( i3 == 0 ) {
		printf(" In forceSplit: I can't find vertex %d as link of star %d \n ", vL, v0);
		printStar(star);
		return EGADS_INDEXERR;
	}
	dist[0] = -1; dist[1] = -1;
	if ( qm -> mesh -> vType [ v0 - 1] >= 0 ) { // boundary vertex: Special split since regular = 3 ( not val 4)
		for ( q = 0 ; q < star -> nQ; q++ )
			if ( star -> quads[ q ] == -1 ) break;
		id6[0]  = - 1; id6[1] = -1;
		for ( i = 0 ; i < 2 ; i++ ) {
			id0 = star -> idxQ [ q + i ] ;
			if ( i == 1 ) dist[0] = (star -> nV - 1) - 4;
			else          dist[0] = 4 + 2 * i;
			links[2 * i     ]  = 2 * id0 + 1;
			links[2 * i + 1 ]  = star -> idxV [ 2 * id0 + 1 + dist[0] ];
			if (         getValence (qm -> mesh, star -> verts [ links[2*i    ] ]) == 4 ) {
				if      (getValence (qm -> mesh, star -> verts [ links[2*i + 1] ]) == 3 && id6[0] == -1 ) id6[0] = i;
				else if (getValence (qm -> mesh, star -> verts [ links[2*i + 1] ]) == 4 && id6[1] == -1 ) id6[1] = i;
			}
		}
		dist[0]    = 4;
		if      ( id6[0] != -1 ) id0 = links[ 3 * id6[0]];
		else if ( id6[1] != -1 ) id0 = links[ 3 * id6[1]];
		else {
			EG_freeStar(&star);
			return EGADS_SUCCESS;
		}
		stat      = EG_splittingOperation (qm, star, id0, dist[0] );
		printMesh(qm, buffer,0);
	}
	else {
		if      ( star -> nQ <= 5 ) dist[0] = 4;
		if      ( star -> nQ == 5 ) dist[1] = 6;
		else if ( star -> nQ == 6 ) dist[0] = 6;
		id6[0] = star -> verts [ star -> idxV [ i3 + dist[0]] ];
		id6[1] = star -> verts [ star -> idxV [ i3 + dist[1]] ];
		if ( getValence ( qm -> mesh, id6[0] ) > getValence ( qm -> mesh, id6[1] ) )
			swapInt ( &dist[0], &dist[1]);
		for ( i = 0 ; i < 2; i++ ) {
			if ( dist[i] == -1 ) break;
			vOpp      = star -> verts [ star -> idxV [ i3 + dist[0]] ];
			//if ( getValence ( qm -> mesh, vOpp ) > 4 ) continue;
			printf(" Forced Split at centre %d and links %d  %d \n",
					v0, vL, vOpp );
			printVertSpecs (qm -> mesh, v0);
			printVertSpecs (qm -> mesh, vL);
			printVertSpecs (qm -> mesh, vOpp);
			stat      = EG_splittingOperation (qm, star, i3, dist[0] );
			*activity = 1;
			break;
		}
	}
	EG_freeStar ( &star );
	if ( stat != EGADS_SUCCESS )
		printf(" EG_forceSplit %d !!!! \n ", stat );
	printMesh(qm, buffer,0);
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
		if  ( val[0] < 5 || qm -> mesh -> vType [ poly[0] -1 ] == qm -> mesh -> sizeQuads * 2 ) continue;
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
			if ( qm -> mesh -> vType [ poly[0] - 1] >= 0 ) { // boundary vertex: Special split since regular = 3 ( not val 4)
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
				printf(" I6 %d %d \n ", id6[0], id6[1] );
				dist = 6;
				id0  = id6[0]; if ( id0 < 0 ) id0 = id6[1];
			}
		}
		if ( validSplit == 1)  {
			printf("WE FOUND A CANDIDATE FOR SPLITTING: SPLIT VERTEX V %d FROM LINKS %d  and %d dist %d \n",
					star -> verts[0], star -> verts[id0], star -> verts[star->idxV[id0 + dist]], dist);
			stat      = EG_splittingOperation(qm, star, id0, dist);
#ifdef DEBUG
			printVertSpecs ( qm -> mesh, star -> verts[0]);
			printVertSpecs ( qm -> mesh, star -> verts[id0]);
			printVertSpecs ( qm -> mesh, star -> verts[star->idxV[id0 + dist]]);
			printMesh (qm, buffer, 0);
#endif
			*activity = 1;
			EG_freeStar (&star);
			return stat;
		}
		EG_freeStar ( &star);
	}
	return EGADS_SUCCESS;
}


static int
EG_splittingOperation(meshMap *qm, vStar *star, int id0, int distanceToId0) {
	int qIdx[4], modQ[4], verts[4], adj[2], poly[4], q, newQ, i, j, stat;
	vStar *inStar = NULL;
	double eval[18], uvAv[2];
	if ( star == NULL ) {
		printf(" I entered with a NULL STAR \n ");
		return EGADS_MALLOC;
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
			return EGADS_SUCCESS;
		}
		++qm -> mesh -> totVerts; ++qm -> mesh -> totQuads;
	}
	qm -> mesh -> vType [poly[3] - 1] = -1;
	poly [0] = star -> verts[   0];
	poly [1] = star -> verts[ id0];
	poly [2] = star -> verts[ star -> idxV [ id0 + distanceToId0]];
	qIdx [0] = (int)(id0 - 1)/2;
	qIdx [1] = star -> idxQ  [ qIdx[0] + distanceToId0/2  - 1];
	qIdx [2] = star -> idxQ  [ qIdx[0] + distanceToId0/2     ];
	qIdx [3] = qIdx[0] - 1; if ( qIdx[3] == -1 ) qIdx[3] =  star -> nQ - 1;
	verts[0] = poly[1];
	verts[1] = poly[2];
	verts[2] = poly[2];
	verts[3] = poly[1];
	qm -> mesh -> quadIdx[ 4 * ( newQ - 1)    ] = poly[1];
	qm -> mesh -> quadIdx[ 4 * ( newQ - 1) + 1] = poly[0];
	qm -> mesh -> quadIdx[ 4 * ( newQ - 1) + 2] = poly[2];
	qm -> mesh -> quadIdx[ 4 * ( newQ - 1) + 3] = poly[3];
	//printf(" New quad %d \n", newQ);
	//for ( i = 0 ; i < 4; i++ ) printf(" v %d %d \n", i, qm -> mesh -> quadIdx[ 4 * ( newQ - 1)  +i  ]);
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
	qm -> mesh -> uvs[ 2 * ( poly[3] - 1 )    ] = qm -> mesh -> uvs[ 2 * ( poly[0] - 1 )    ];
	qm -> mesh -> uvs[ 2 * ( poly[3] - 1 ) + 1] = qm -> mesh -> uvs[ 2 * ( poly[0] - 1 ) + 1];
	for ( j = 0 ; j < 4; j = j + 3) {
		if ( qm -> mesh -> vType [ poly[j] -1 ] != -1 ) continue;
		stat = EG_buildStar (qm -> mesh, &inStar, poly[j] );
		if ( stat != EGADS_SUCCESS || inStar == NULL ) return stat;
		uvAv[0] = 0.0;
		uvAv[1] = 0.0;
		for ( i = 0 ; i < inStar -> nQ; i++ ) {
			if ( inStar -> quads[i] == newQ ) continue;
			stat = quadAverageCoords(qm, inStar -> quads[i], eval, &eval[2]);
			uvAv[0] += eval[0];
			uvAv[1] += eval[1];
		}
		uvAv[0] /= (double) ( inStar -> nQ -1 );
		uvAv[1] /= (double) ( inStar -> nQ -1 );
		stat = EG_evaluate ( qm -> face, uvAv, eval );
		for ( i = 0 ; i < 3; i++ ) {
			if ( i < 2 ) qm -> mesh -> uvs [ 2 * ( poly[j] - 1) + i] = uvAv[i];
			qm -> mesh -> xyzs [ 3 * ( poly[j] - 1) + i] = eval[i];
		}
	}
	EG_freeStar (&inStar);
	stat = optimize_angles(qm, 4, poly, 0);
	for ( i = j = 0 ; i < qm -> mesh -> totQuads; i++ ) {
		if ( checkQuad (qm -> mesh, i + 1) == EGADS_SUCCESS ) j++;
	}
	qm ->extraQuads = j - qm -> optNquads;
	//fprintf(stderr, " REMOVED QUADS / VERTICES :: %d  %d   extraquads %d \n ", qm -> mesh -> remQuads[0], qm -> mesh -> remVerts[0], qm ->extraQuads );
	return stat;
}






static int
EG_doubleSwap( meshMap *qm, quadGroup qg, int forcing, int *activity ) {
	int  piv5 = -1, piv3 = -1, q5, i, link5, stat, adj[2], swap = 0, v30;
	*activity = 0;
	quadGroup auxQg;
	if ( forcing ==1 ) printf("DOUBLE SWAP FORCING SWAP\n");
	if      (validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0 ) return EGADS_SUCCESS;
	if      (qg. vals[0] == 3 || qg. vals[3] == 3 ) return EGADS_SUCCESS;
	else if (qg. vals[0] == 4 &&    forcing  == 0 ) return EGADS_SUCCESS;
	else if (qg. vals[0] == 5 ) {
		if  (qg. vals[2] * qg. vals[4] == 15 ) {
			piv5 = 2; swap = 4;
			if ( qg. vals[2] == 3 ) {
				piv5 = 4; swap = 2;
			}
		}
		else if ( qg.vals[1] == 3 || qg.vals[5] == 3 ) {
			piv3 = 1; if ( qg.vals[1] != 3 ) piv3 = 5;
			piv5 = (piv3 + 3)%6;
			swap = piv3;
			if ( forcing == 0 && getValence(qm -> mesh, qg.verts[piv5] ) != 4 ) return EGADS_SUCCESS;
			printf("DIAGONAL SWAPPING PIV 5 %d = %d  piv 3 = %d\n", piv5, qg.verts[piv5], piv3 );
		}
	}
	else if ( forcing == 1 ) {
		if  ( qg. vals[1] * qg. vals[5] != 15 )     return EGADS_SUCCESS;
		piv5 = 1; swap = 5;
		if ( qg. vals[1] == 3 ) {
			piv5 = 5; swap = 1;
		}
	}
	if ( piv5 == -1 || qm -> mesh -> vType[ qg.verts[piv5] -1 ] != -1 ) return EGADS_SUCCESS;
	link5 = ( piv5 + 1 )%6;
	v30   = qg. verts[ ( piv5 + 5 )%6];
	if ( link5 % 3 == 0 ) {
		v30   = qg. verts[ link5 ];
		link5 = (piv5 + 5 )%6;
	}
	q5 = 0;
	if ( piv5 > 3 ) q5 = 1;
	printQuadGroup (qm -> mesh, qg);
	if (piv3 != -1 )   // diagonal quad swapping
		stat     = EG_adjQtoPair ( qm -> mesh, qg.q[q5], qg.verts[piv5], v30, adj);
	else
		stat     = EG_adjQtoPair ( qm -> mesh, qg.q[q5], qg.verts[piv5], qg. verts[link5], adj);
	if ( stat != EGADS_SUCCESS || adj[1] == -1 ) {
		printf(" EG_doubleSwap adjToPair %d --> %d !!\n", adj[1], stat );
		return stat;
	}
	printf(" AUX GROUP\n");
	stat       = EG_createQuadGroup (qm -> mesh, &auxQg, qg. q[q5], adj[1] );
	printQuadGroup(qm -> mesh, auxQg);
	if ( stat != EGADS_SUCCESS ) {
		printf("EG_doubleSwap aux swap: EG_createQuadGroup --> %d !!\n", stat );
		printQuadGroup (qm -> mesh, auxQg);
		return stat;
	}
	/*if ( piv3 != -1 ) {
		for ( i = 0 ; i < 6; i++)
			if (auxQg.verts[i] == qg.verts[link5] ) break;
		printf("found link5 %d = %d in i %d \n ", link5, qg.verts[link5], i );
		if (  auxQg.vals [i] <= 4 ) return EGADS_SUCCESS;

	}
	else {*/
	for ( i = 0 ; i < 6; i++)
		if (auxQg.verts[i] == v30 ) {
			printf(" OPPO TO v30 %d -> %d\n ", v30, auxQg.verts [ ( i + 3 ) %6 ] );
			if (auxQg.vals [ ( i + 3 ) %6 ] < 5 || qm -> mesh -> vType [ auxQg.verts[ ( i + 3 ) %6 ] -1 ] == 0 )  return EGADS_SUCCESS;
		}
	//}
	printf(" ----- > DOUBLE SWAP THRU %d\n ", swap);
	printMesh(qm, buffer, 0);
	if ( piv3!= -1 )			fprintf(stderr,"DIAGONAL SWAPPING \n" );
	stat       = EG_swappingOperation (qm, qg, swap );
	printMesh(qm, buffer, 0);
	*activity  = 1;
	if ( stat != EGADS_SUCCESS ) {
		printf(" EG_doubleSwap: at first swap --> %d !!\n ", stat );
		return stat;
	}
	if ( piv3 != -1 ) {
		i = EG_quadVertIdx ( qm -> mesh, adj[1], qg. verts[piv5] );
		link5 = qm ->mesh -> quadIdx [ 4 * ( adj[1] -1 ) + ( i+ 1)%4];
		if ( link5 == v30) link5 = qm ->mesh -> quadIdx [ 4 * ( adj[1] -1 ) + ( i+ 3)%4];
		q5    = adj[1];
		stat  = EG_adjQtoPair ( qm -> mesh, q5, qg.verts[piv5], link5, adj);
		if ( stat != EGADS_SUCCESS || adj[1] == -1 ) {
			printf("EG_doubleSwap: afterswapping using diagonal swap adjacent for group is %d stat %d\n ", adj[1], stat );
			return stat;
		}
		qg.q[q5] = q5;
	}
	else {
		for ( q5 = 0; q5 < 2; q5++)
			if ( EG_quadVertIdx ( qm -> mesh, qg. q[q5], qg. verts[piv5]) >= 0 ) break;
	}
	stat       = EG_createQuadGroup (qm -> mesh, &qg, qg. q[q5], adj[1] );
	if ( stat != EGADS_SUCCESS ) {
		printf("EG_doubleSwap after first swap: EG_createQuadGroup --> %d !!\n", stat );
		printQuadGroup (qm -> mesh , qg );
		return stat;
	}
	for ( i = 0 ; i < 6; i++)
		if ( qg. verts[i] == v30 ) break;
	if (   qg. verts[i] != v30 ) {
		printf(" EG_doubleSwap:: I should have found vertex %d in new quad group \n ", v30 );
		printQuadGroup (qm -> mesh, qg );
		return EGADS_INDEXERR;
	}
	printf(" TRYING SECOND OPERATION \n");
	printMesh(qm, buffer, 0);
	if ( validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0 ) return EGADS_SUCCESS;
	return EG_swappingOperation (qm, qg, i );
}

static int EG_doubleSplit(meshMap *qm, quadGroup qg, int *activity ) {
	int i, j, stat;
	*activity  = 0;
	int piv[2] = {1, 5} ;
	if ( qg.vals[0] != 5 || qg.vals[1] * qg.vals[5] != 15 ) return EGADS_SUCCESS;
	if ( qg.vals[1] >= 5 ) {
		piv[0] = 5; piv[1] = 1;
	}
	if (qm -> mesh -> vType [ qg.verts[0]      - 1 ] >= 0 ||
			qm -> mesh -> vType [ qg.verts[piv[1]] - 1 ] >= 0 ) return EGADS_SUCCESS;
	*activity  = 1;
	printf(" FORCE SPLIT 1 \n ");
	printMesh(qm, buffer,0);
	stat       = EG_forceSplit (qm, qg. verts[0], qg.verts[piv[0]], &i);
	if ( stat != EGADS_SUCCESS || i == 0 ) {
		printf("In EG_doubleSplit: force 1st split through %d - %d --> %d !!\n ",
				stat, qg.verts[0], qg.verts[piv[0]] );
		return stat;
	}
	printMesh(qm, buffer,0);
	for ( j = 0 ; j < 2; j++)
		if ( EG_quadVertIdx (qm -> mesh, qg.q[j], qg.verts[piv[1]] ) >= 0 ) break;
	printf(" CALL SPLIT FROM QUAD %d \n ", qg.q[j] );
	stat = EG_split (qm, qg.q[j], &i);
	if ( i == 0 || stat != EGADS_SUCCESS ) {
		printf("In EG_doubleSplit: 2nd split around %d didn't work  --> act = %d stat = %d !!\n ",
				qg.verts[piv[1]], i, stat );
		return stat;
	}
	return stat;
}

static int EG_swapSplit(meshMap *qm,quadGroup qg, int forcing, int *activity  ) {
	int  stat, i, j, i3 = -1, i5 = -1, v3opp = -1, q5, vL5, vL5adj, swap = 0, adj[2];
	*activity = 0;
	printQuadGroup (qm -> mesh, qg );

	if  ( qg.vals[0] * qg.vals[3] != 20 ||
			validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0 )return EGADS_SUCCESS;
	if      (qg.vals[1] == 3 ) i3 = 1;
	else if (qg.vals[5] == 3 ) i3 = 5;
	for ( i = 1; i < 6; i++)
		if (qg.vals[i] == 5 )  i5 = i;
	if (i3   != -1 && ( i5 == -1 ||( i5 != -1 && i5 == (i3 + 3)%6 ) )) {
		i5    = ( v3opp +3)%6;
		v3opp = i3;
	}
	else if ( i3 == -1 && forcing == 1 && ( i5 == -1 || i5 == 2 || i5 == 4 )) {
		//if ( qm -> extraQuads >= 0 ) return EGADS_SUCCESS;
		if ( i5 == -1 ) v3opp = -1;
		else            v3opp = i5;
		printf(" VOPP %d look for nearby things\n ", v3opp );
	}
	else return EGADS_SUCCESS;
	if ( forcing ==1 ) printf("SWAP SPLIT FORCING SWAP  i3 %d i5 %d\n", i3, i5);
	if ( v3opp == -1 ) {
		for ( i  = 0 ; i < 2; i++ ) {
			j    = 2 + 2 * i;
			if ( i == 0 ) vL5 = 1;
			else          vL5 = 5;
			stat              = EG_adjQtoPair (qm -> mesh, qg.q[i], qg.verts[j], qg.verts[vL5], adj );
			printf(" Adj to %d from edge %d %d is %d \n ", qg.q[i], qg.verts[j], qg.verts[vL5],adj[1]);
			if ( stat        != EGADS_SUCCESS || adj[1] == -1 ) continue;
			q5                = EG_quadVertIdx (qm -> mesh, adj[1], qg.verts[j]);
			vL5adj            = qm -> mesh -> quadIdx [ 4 * ( adj[1] - 1 ) + ( q5 + 1 ) %4 ];
			if ( vL5adj == qg.verts[vL5] )
				vL5adj = qm -> mesh -> quadIdx [ 4 * ( adj[1] - 1 ) + ( q5 + 3) %4 ];
			if      ( getValence (qm ->mesh, vL5adj) == 3 ) {
				printf(" forcing thru %d \n ", j );
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
		printf(" v3 != -1  Adj to %d edge %d %d is %d \n ", qg.q[q5], qg.verts[v3opp], qg.verts[vL5], adj[1]);
		if ( stat != EGADS_SUCCESS || adj[1] == -1 ) {
			printf("EG_swapSplit: EG_adjQtoPair from quad %d is %d --> %d\n!!", qg.q[q5], adj[1], stat );
			return stat;
		}
		i        = EG_quadVertIdx (qm -> mesh, adj[1], qg.verts[v3opp]);
		vL5adj   = qm -> mesh -> quadIdx [ 4 * ( adj[1] - 1 ) + ( i + 1 ) %4 ];
		if ( vL5adj == qg.verts[vL5] )
			vL5adj = qm -> mesh -> quadIdx [ 4 * ( adj[1] - 1 ) + ( i + 3) %4 ];
		printf(" vL5adj is %d \n ", vL5adj);
		if ( i3 != -1 && ( qg.vals[v3opp] == 5 || getValence (qm ->mesh, vL5adj) == 3 ) ) swap = i3;
		else if ( forcing == 1 /*&& qm -> extraQuads <= 0*/ && ( qg.vals[v3opp] == 5 || getValence (qm ->mesh, vL5adj) == 3 ) ) swap = v3opp;
	}
	if ( swap == 0 ) return EGADS_SUCCESS;
	printMesh(qm, buffer,0);
	printf("Swap split::: swapping thru %d \n ", swap);
	printQuadGroup (qm -> mesh, qg );

	stat       = EG_swappingOperation( qm, qg, swap );
	*activity  = 1;
	if ( stat != EGADS_SUCCESS) {
		printf(" In swapSplit thru %d : EG_swappingOperation went %d !!\n ", swap, stat );
		printQuadGroup(qm ->mesh, qg);
		return stat;
	}
	printMesh(qm, buffer,0);
	printf(" FORCE SPLIT THRU %d \n ", qg.verts[3]);
	if ( i5 == -1 ) {
		printf(" WWWWWW %d \n ", i5 );
		exit (1);
	}
	stat = EG_forceSplit (qm, qg.verts[i5], qg.verts[3], &i);
	printf(" Stat force split %d !\n ", stat);
	return stat;
}


static int EG_swapCollapse (meshMap *qm,quadGroup qg, int forcing, int *activity  ) {
	int  stat, i, i3 = -1, q5, qC, vL5, vL5adj, vOpp3, swap = 0, adj[2];
	*activity   = 0;
	printf(" Swap collapse group %d %d \n ", qg.q[0], qg.q[1] );
	printf(" VALID SWAP %d %d  ?? %d \n ",
			qg.verts[0], qg.verts[3],validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ));
	if (qg.vals[0] * qg.vals[3] != 16 ||
			validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0 ) return EGADS_SUCCESS;
	if (qg.vals[1] * qg.vals[2] == 15 ) {
		i3 = 1;
		if ( qg.vals[i3] == 5 ) {
			i3 = 2;
		}
	}
	else if (qg.vals[4] * qg.vals[5] == 15 ) {
		i3 = 4;
		if ( qg.vals[i3] == 5 ) {
			i3 = 5;
		}
	}
	else if ( forcing == 1 ) {
		//if ( qm  -> extraQuads < 0 ) return EGADS_SUCCESS;
		if( qg.vals[1] * qg.vals[5] == 25 ||
				qg.vals[2] * qg.vals[4] == 25 ) {
			i3 = 2; if ( qg.vals[1] != 5 ) i3 = 1;
		}
		else if( qg.vals[1] * qg.vals[4] == 15 ) {
			i3 = 1; if ( qg.vals[1] != 3 ) i3 = 4;
			printf(" 3,5 pair vertex 3 %d and 5 %d \n ", qg.verts[i3], qg.verts[(i3 + 3)%6] );
		}
		else if ( qg.vals[2] * qg.vals[5] == 15 ) {
			i3 = 2; if ( qg.vals[2] != 3 ) i3 = 5;
			printf(" 3,5 pair vertex 3 %d and 5 %d \n ", qg.verts[i3], qg.verts[(i3 + 3)%6] );
		}
		else if(qg.vals[1] * qg.vals[5] == 9 ||
				qg.vals[2] * qg.vals[4] == 9 ) {
			printMesh(qm, buffer,0);
			printf(" Forcing In swapCollapse -> using 33 pairs \n ");
			stat      = EG_swappingOperation (qm, qg, 1 );
			*activity = 1;
			if ( stat != EGADS_SUCCESS ) {
				printf("forcing swapcollapse:: after swapping %d \n ", stat);
				return stat;
			}
			qC = qg.q[0];
			if ( EG_nValenceCount ( qm -> mesh, qC, 3 ) < 2 ) qC = qg.q[1];
			printMesh(qm, buffer,0);
			stat = EG_forceCollapse (qm, qC, &i);
			if ( stat != EGADS_SUCCESS )printf("forcing swapcollapse:: after swapping %d \n ", stat);
			printMesh(qm, buffer,0);
			return stat;
		}
		else return EGADS_SUCCESS;
		printf(" forcing found something\n ");
		printQuadGroup (qm -> mesh, qg);
	}
	else return EGADS_SUCCESS;
	if ( forcing ==1 ) printf("SWAP COLLAPSE FORCING SWAP\n");
	vOpp3  = (i3 + 3 ) % 6;
	q5     = 0; if ( vOpp3 > 3 ) q5 = 1;
	qC     = qg.q[ ( q5 +1)%2];
	vL5    = ( vOpp3 + 1) %6;
	if ( vL5 %3 != 0 ) vL5 = ( vOpp3 + 5) %6;
	if ( qg.vals[vOpp3] == 3 ) return EGADS_SUCCESS;
	stat   = EG_adjQtoPair (qm -> mesh, qg.q[q5], qg.verts[vOpp3], qg.verts[vL5], adj );
	if ( stat != EGADS_SUCCESS || adj[1] == -1 ) return stat;
	i      = EG_quadVertIdx ( qm -> mesh, adj[1], qg.verts[vOpp3]);
	vL5adj = qm -> mesh -> quadIdx [ 4 * ( adj[1] - 1 ) + ( i + 1 ) %4 ];
	if ( vL5adj == qg.verts[vL5] ) vL5adj = qm -> mesh -> quadIdx [ 4 * ( adj[1] - 1 ) + ( i + 3 ) %4 ];
	if (( forcing == 0 && qg.vals[vL5] == 4 && getValence (qm -> mesh, vL5adj) >= 4 ) ||
			( forcing == 1 /*&& qm -> extraQuads >= 0*/ && getValence (qm -> mesh, vL5adj) > 4 ) ) return EGADS_SUCCESS;
	stat  = EG_createQuadGroup (qm -> mesh, &qg, qg.q[q5], adj[1]);
	if ( stat != EGADS_SUCCESS ) {
		printf("EG_swapCollapse before swap: EG_createQuadGroup --> %d !!\n", stat );
		printQuadGroup (qm -> mesh, qg );
		return stat;
	}
	if (validSwap (qm -> mesh, qg.verts[0], qg.verts[3] ) == 0 ) return EGADS_SUCCESS;
	for ( swap = 0 ; swap < 6; swap++ )
		if ( qg.verts[swap] == vL5adj ) break;
	printf(" SWAP COLLAPSE: SWAPPING  through %d !\n ", swap);
	if ( forcing == 1 ) printf(" SWAP COLLAPSE forcing: SWAPPING !\n ");
	printMesh(qm, buffer,0);
	stat      += EG_swappingOperation (qm, qg, swap );
	printMesh(qm, buffer,0);
	*activity  = 1;
	if ( stat != EGADS_SUCCESS ) {
		printf("EG_swapCollapse after swapping %d !!\n", swap );
		return stat;
	}
	vL5  = qg.verts[0];
	stat = EG_forceCollapse (qm, qC, &i);
	printMesh(qm, buffer,0);
	if ( stat != EGADS_SUCCESS || i > 0 ) return stat;
	printf(" NOW CREATE QUAD GROUP %d %d \n ", qg.q[q5], adj[1] );
	stat  = EG_createQuadGroup (qm -> mesh, &qg, qg.q[0], qg.q[1]);
	for ( i = 0 ; i < 6; i++ ) if ( qg.verts[i] == vL5 ) break;
	stat      += EG_swappingOperation (qm, qg, i );
	*activity  = 0;
	return stat;
}

static int EG_cleanMesh (meshMap *qm, int domain, int operType, int *activity ) {
	int i, act, stat, count = 0;
	*activity  = 0;
	stat       = restoreMeshData(qm, qm -> backupMesh, qm -> mesh);
	if ( stat != EGADS_SUCCESS ) return stat;
	for ( count = i = 0 ; i < qm -> mesh -> totQuads; i++ ) {
		if ( checkQuad ( qm -> mesh, i + 1) != EGADS_SUCCESS ||
				( domain == SINGLE && EG_quadIsBoundary (qm -> mesh, i + 1) != 0 ) ) continue;
		printf(" in clean mesh call clean quad \n ");
		if ( EG_cleanQuad (qm, i + 1, operType, 0, 0, &act ) != EGADS_SUCCESS ) {
			printf("\n\nInside EG_cleanMesh: EG_cleanQuad %d --> %d!!\n", i + 1, stat );
			stat = restoreMeshData(qm, qm -> mesh, qm -> backupMesh );
			if ( stat != EGADS_SUCCESS ) return stat;
			act = 0 ;
		}
		if ( act > 0 ) {
			stat = restoreMeshData(qm, qm -> backupMesh, qm -> mesh);
			if ( stat != EGADS_SUCCESS ) return stat;
		}
		count += act;
		if ( act == 0 && checkQuad ( qm -> mesh, i + 1) == EGADS_SUCCESS ) {
			if ( EG_cleanNeighborhood ( qm, i + 1, operType, 0, &act ) != EGADS_SUCCESS ) return stat;
		}
		count += act;
	}
	*activity = count;
	return EGADS_SUCCESS;
}



static int EG_cleanNeighborhood (meshMap *qm, int qID, int adjacent, int transfer, int *activity ) {
	int i, act, stat = 0, count = 0, j, v[4];
	vStar *star = NULL;
	*activity   = 0;
	stat        = checkQuad (qm -> mesh, qID );
	printf(" IN CLEAN NEIGHBORHOOD %d -> %d\n ", qID, stat);
	if ( stat  != EGADS_SUCCESS ) {
		printf("EG_cleanNeighborhood : Quad %d is invalid = %d \n", qID, stat );
		return checkMesh (qm);
	}
	stat       = restoreMeshData(qm, qm -> backupMesh, qm -> mesh);
	if ( stat != EGADS_SUCCESS ) return stat;
	for ( i = 0 ; i < 4; i++ ) v[i]  = qm -> mesh -> quadIdx[ 4 * ( qID -1 ) + i ];
	printf(" in clean neighbor call clean quad \n ");
	stat = EG_cleanQuad ( qm, qID, adjacent, transfer, 0, &act );
	if ( stat != EGADS_SUCCESS ) {
		printf("\n\nInside EG_cleanNeighborhood: EG_cleanQuad %d --> %d!!\n", qID, stat );
		stat = restoreMeshData(qm, qm -> mesh, qm -> backupMesh );
		if ( stat != EGADS_SUCCESS ) return stat;
		act = 0;
	}
	if ( act > 0 ) stat = restoreMeshData(qm, qm -> backupMesh, qm -> mesh);
	if ( stat != EGADS_SUCCESS ) return stat;
	*activity += act;
	printf(" BEFORE NEIGHBORS: ACTIVITY %d\n ", *activity);
	for ( i = 0 ; i < 4; i++ ) {
		printf(" v[%d] = %d\n ", i, v[i] );
		if ( checkVertex ( qm -> mesh, v[i] ) != EGADS_SUCCESS ) continue;
		stat = EG_buildStar (qm -> mesh, &star, v[i]);
		if ( star == NULL || stat != EGADS_SUCCESS ) {
			printf("EG_cleanNeighborhood : EG_buildStar  = %d \n", stat );
			return stat;
		}
		for ( j = 0 ; j < star -> nQ; j++ ) {
			if ( star -> quads[j] == -1 ) continue;
			printf(" in clean neighbor in loop call clean quad j %d = %d activity %d\n ", j, star -> quads[j], count );
			if (  checkQuad ( qm -> mesh, star -> quads[j]) != EGADS_SUCCESS ) {
				printf(" QUAD %d is %d \n ", star -> quads[j], checkQuad ( qm -> mesh, star -> quads[j]));
				continue;
			}
			//printf(" Inside star loop centred at %d quad %d neighbor\n", star -> verts[0], star -> quads[j]);
			if ( EG_cleanQuad (qm, star -> quads[j], adjacent, transfer, 0, &act ) != EGADS_SUCCESS ) {
				stat = restoreMeshData (qm, qm -> mesh, qm -> backupMesh);
				if ( stat != EGADS_SUCCESS ) {
					EG_freeStar(&star);
					return stat;
				}
				act  = 0;
			}
			if ( act > 0 ) {
				stat = restoreMeshData(qm, qm -> backupMesh, qm -> mesh);
				if ( stat != EGADS_SUCCESS ) {
					EG_freeStar ( &star );
					return stat;
				}
				break;
			}
			count += act;
			printf(" j %d / %d\n", j, star -> nQ);
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
	int     ITMAX, it = 0, stat, activity = 0, totActivity = 0, aux, i, j, k, addQuad = 0, q, transfer = 0;
	int     best_iV, iV, quadPair[2], totV, vQ, iV0, kt, kb;
	// GET RANGE FOR EACH POINT
	stat       = checkMesh (qm);
	if ( stat != EGADS_SUCCESS) {
		fprintf(stderr,"At meshRegularization: starting mesh is invalid --> %d !! I'm leaving program \n ", stat);
		return stat;
	}
	meshCount (qm -> mesh, &iV, &totV, &vQ);
	fprintf(stderr,"============= Original mesh has %d QUADS, %d / %d irregular vertices (%f percent) ============\n", vQ, iV, totV, (double) iV * 100.0 / (double)totV );
	ITMAX = 100;
	restoreMeshData (qm, qm -> backupMesh, qm -> mesh);
	for ( k = 0 ; k < 2; k++ ) {
		kb                   = INTERIOR;
		kt                   = SINGLE;
		if ( k % 2 == 1 ) kb = FULL;
		if ( k     >= 2 ) kt = DOUBLE;
		kb = FULL;
		it = 0;
		do {
			stat       = EG_cleanMesh (qm, kb, kt, &totActivity );
			if ( stat != EGADS_SUCCESS ) {
				restoreMeshData (qm, qm -> mesh, qm -> backupMesh );
				stat        = EGADS_SUCCESS;
				totActivity = 0;
			}
			it++;
			meshCount(qm -> mesh, &iV, &totV, &vQ);
			if ( totActivity > 0 ) restoreMeshData (qm, qm -> backupMesh, qm -> mesh);
		}
		while (totActivity > 0 && it < ITMAX && iV > 2);
		printf(" Using %d domain ( 0 = interior, 1 = full) and operations type %d ( 0 = basic, 1 = compo) it %d act %d\n", kb, kt, it, totActivity );
		snprintf(buffer,500,"DONE%d_%d_%d",kb, kt, qm -> fID );
		printMesh ( qm, buffer, 1 );
		fprintf(stderr," In First round: k = %d OUT BECAUSE TOTACT %d it %d iv %d \n", k, totActivity, it, iV);
		if ( iV <= 2 ) break;
	}
	stat = restoreMeshData (qm, qm -> bestMesh, qm -> mesh);
	if ( stat != EGADS_SUCCESS ) goto cleanup;
	best_iV = iV;
	fprintf(stderr,"*******************************************************************\n");
	fprintf(stderr," AFTER BASIC: mesh has %d QUADS, %d / %d irregular vertices (%f percent) ============\n", vQ, iV, totV, (double) iV * 100.0 / (double)totV );
	fprintf(stderr, "*******************************************************************\n");
	snprintf(buffer,500,"BASIC_%d", qm -> fID );
	printMesh ( qm, buffer, 1 );
	if ( iV > 2 ) {
		totActivity = 0 ;
		stat        = restoreMeshData (qm, qm -> backupMesh, qm -> mesh);
		iV0         = iV;
		it          = 0 ;
		stat        = optimize_angles(qm, 0, NULL, 2);
		do {
			totActivity = 0;
			printf("*******************************************************************\n");
			for ( k = 0 ; k < 2; k++ ) {
				stat = checkMesh (qm);
				stat = restoreMeshData (qm, qm -> backupMesh, qm -> mesh);
				fprintf(stderr,"********************  ITERATION %d : TOTACT %d : iV  %d  MESH %d ************************\n",
						it, totActivity, iV, MESHPLOTS);
				printf("********************  ITERATION %d : TOTACT %d : iV  %d  MESH %d ************************\n",
						it, totActivity, iV, MESHPLOTS);
				if ( k == 0 ) fprintf(stderr," TRANSFER VALENCES USING 35 PAIRS \n");
				else          fprintf(stderr," TRANSFER VALENCES ALSO USING 55, 33 PAIRS AND FORCING FUNCTIONS \n");
				for ( transfer = q = 0 ; q < qm -> mesh-> totQuads; q++) {
					printf("============== TRANSFER VALENCES AROUND QUAD  %d -> %d \n", q + 1, checkQuad( qm -> mesh, q + 1 ) );
					if ( checkQuad( qm -> mesh, q + 1 ) != EGADS_SUCCESS ) continue;
					quadPair[0] = q + 1;
					quadPair[1] = -1;
					printf("============== TRANSFER VALENCES AROUND QUAD  %d FORCING %d  \n", q + 1, k);
					stat       = EG_transferValences ( qm, quadPair, k, &transfer, &activity );
					stat      += checkMesh(qm);
					if ( stat != EGADS_SUCCESS ) {
						printf(" ERRTRANSFER VALENCE  %d  ", quadPair[0]);
						fprintf(stderr," ERRTRANSFER VALENCE  %d  \n", quadPair[0]);
						stat = restoreMeshData (qm, qm -> mesh, qm -> backupMesh);
						if ( stat != EGADS_SUCCESS ) goto cleanup;
						activity  = 0;
					}
					if ( activity == 0) continue;
					totActivity += activity;
					printf(" transfer valences has done something : now quad pairs %d %d stat = %d activity %d  TRANSFER %d \n",
							quadPair[0], quadPair[1], stat, activity, transfer );
					j         = 0;
					if ( quadPair[0] == -1 ) continue;
					while ( activity > 0 && iV > 2 && j < 20 ) {
						activity  = 0;
						printf(" Clean Neighborhood in quad %d TRANSFER %d \n ", quadPair[0], transfer);
						for ( i = 0 ; i < 2; i++ ) {
							printf(" i = %d -> quad %d transfer %d \n ", i, quadPair[i], transfer );
							if ( quadPair[i] == -1 || checkQuad (qm -> mesh, quadPair[i] ) != EGADS_SUCCESS ) continue;
							stat      = EG_cleanNeighborhood (qm, quadPair[i], 1, transfer, &activity );
							printf(" IN i loop stat %d neighbor activity %d\n ", stat, activity);
							if ( stat != EGADS_SUCCESS) goto cleanup;
							if ( activity > 0 ) break;
						}
						totActivity += activity ;
						if ( activity > 0 ) break;
						printf(" Inside while loop for transfer keep on moving  STILL AT  %d %d stat = %d  transfer %d addq  %d\n",
								quadPair[0], quadPair[1], stat, transfer, addQuad);
						stat      = EG_transferValences ( qm, quadPair, k, &transfer, &activity );
						printf(" AFTER TRANSFER Inside while loop for transfer :  WE ARE STILL AT  %d %d stat = %d ACT %d \n",
								quadPair[0], quadPair[1], stat, aux );
						if ( stat != EGADS_SUCCESS ) {
							restoreMeshData (qm, qm -> mesh, qm -> backupMesh);
							if ( stat != EGADS_SUCCESS) goto cleanup;
						}
						if ( quadPair[0] == -1 ) break;
					}
					meshCount(qm -> mesh, &iV, &totV, &vQ );
					if ( iV < iV0 ) {
						stat = restoreMeshData (qm, qm -> backupMesh, qm -> mesh);
						if ( stat != EGADS_SUCCESS ) goto cleanup;
						iV0 = iV;
						fprintf(stderr, " In transfer valences :: mesh has improved to  %d quads \n ", iV );
					}
					if ( iV < best_iV ) {
						stat = restoreMeshData ( qm, qm -> bestMesh, qm -> mesh);
						if ( stat != EGADS_SUCCESS ) goto cleanup;
						fprintf(stderr, " In transfer valences :: best mesh has improved to  %d quads \n ", iV );
						best_iV = iV;
					}
					printMesh ( qm, buffer , 0);
					if ( iV <=2 ) break;
				}
				meshCount(qm -> mesh, &iV, &totV, &vQ );
				if ( iV  < best_iV ) {
					restoreMeshData ( qm, qm -> bestMesh, qm -> mesh);
					best_iV = iV;
				}
				if ( iV <=2 ) break;
			}
			++it;
			stat  = checkMesh (qm );
			if ( stat != EGADS_SUCCESS ) goto cleanup;
			fprintf(stderr, "Tot activity %d it %d iv %d \n", totActivity, it, iV);
			printf("THIS ROUND END Tot activity %d it %d iv %d \n", totActivity, it, iV);
		}
		while ( totActivity > 0 && it < ITMAX && iV > 2);// && iV < iV0);
		if ( iV > best_iV ) {
			stat = restoreMeshData (qm, qm -> mesh, qm -> bestMesh);
			if ( stat != EGADS_SUCCESS ) goto cleanup;
			printf(" Best round had %d quads compared to current %d \n", iV0, iV);
			fprintf(stderr, " Best round had %d quads compared to current %d \n", iV0, iV);
		}
		meshCount(qm -> mesh, &iV, &totV, &vQ );
		fprintf(stderr," NOW MESH HAS %d QUADS AND %d are IRREGULAR ( %d total ). OPTI QUADS %d OVERHEAD %d \n",
				vQ, iV, totV, qm -> optNquads, qm -> extraQuads);
	}
	fprintf(stderr,"*******************************************************************\n");
	fprintf(stderr," FULL REGULARIZATION: mesh has %d QUADS, %d / %d irregular vertices (%f percent) ============\n", vQ, iV, totV, (double) iV * 100.0 / (double)totV );
	fprintf(stderr, "*******************************************************************\n");

	stat = checkMesh(qm);
	if ( stat != EGADS_SUCCESS ) goto cleanup;
	fprintf(stderr," FINAL MESH HAS %d QUADS AND %d are IRREGULAR ( %d total). OPTI QUADS %d OVERHEAD %d \n",
			vQ, iV, totV, qm -> optNquads, qm -> extraQuads);
	printf("\n\n PERFORM GLOBAL OPTIMIZATION\n");
	//goto cleanup;
	stat = resizeQm (qm) ;
	if ( stat != EGADS_SUCCESS ) goto cleanup;
	meshCount(qm -> mesh, &iV, &totV, &vQ );
	printf(" FINAL MESH HAS %d QUADS AND %d are IRREGULAR ( %d total). OPTI QUADS %d OVERHEAD %d \n",
			vQ, iV, totV, qm -> optNquads, qm -> extraQuads);
	printMeshStats(qm, 2);
	snprintf(buffer,500,"wvsRegular_%d.txt", qm -> fID);
	stat = EG_wvsData(qm -> mesh, buffer);
	if ( stat != EGADS_SUCCESS ) goto cleanup;
	snprintf(buffer,500,"BEFORE_GLOBAL%d.txt", qm -> fID);
	printMesh(qm, buffer, 1);
	for ( j = i = 0 ; i < qm -> mesh  -> totVerts; i++ )
		if ( qm -> mesh  -> vType[i] == -1 ) j++;
	stat      = optimize_angles(qm, 0, NULL, 1);
	snprintf(buffer,500,"finalMesh_%d.txt", qm ->fID);
	printMesh (qm, buffer, 1 );
	snprintf(buffer,500,"wvsFinalMesh_%d.txt", qm ->fID);
	stat = EG_wvsData(qm -> mesh, buffer);
	cleanup:
	return stat;
}


static int optimize_angles(meshMap *qm, int nP, /*@unused@*/ /*@null@*/int *pList, int fullRegularization)
{
	int  i, k, kk, v, q, stat = EGADS_SUCCESS, nV = 0, it = 0, itMax = 100;
	int *vID = NULL;
	vStar *star = NULL;
	vID = (int *) EG_alloc ( qm -> mesh  -> totVerts * sizeof ( int ) );
	if ( vID == NULL ) return EGADS_MALLOC;
	if ( fullRegularization == 0 ) { // move around only affected vertices
#ifdef DEBUGG
		printf(" POINT LIST %d \n", nP );
		for ( i = 0 ; i < nP; i++ ) {
			printf(" P %d = %d \n", i, pList[i]);
		}
#endif
		if ( nP == 0 ) return EGADS_SUCCESS;
		for ( nV = i = 0 ; i < nP; ++i) {
			if (  qm -> mesh  -> vType[ pList[i] - 1 ] == - 1 ) vID[nV++] = pList[i];
		}
	}
	else {  // full mesh smoothing
		for ( nV = i = 0 ; i < qm -> mesh  -> totVerts; i++ )
			if ( qm -> mesh  -> vType[i] == -1 ) vID[nV++] = i + 1;
	}
	if ( nV == 0 ) {
		EG_free (vID );
		return EGADS_SUCCESS;
	}
	qm -> MINANGLE = MINVALIDANGLE;
	qm -> MAXANGLE = MAXVALIDANGLE;
	while ( it < itMax ) {
		it++;
		for ( q = v = 0 ; v < nV; v++) {
			//	  printf(" BUILD STAR %d / %d \n", v, vID[v] );
			stat = EG_buildStar (qm -> mesh , &star, vID[v] );
			for ( k = 0 ; k < star ->nQ; k++) {
				if ( star -> quads[k] == -1 ) continue;
				stat = checkInvalidElement(qm, star -> quads[k], qm -> MINANGLE, qm -> MAXANGLE );
				//	      printf(" -+-+-+-+-+-+   INVALID CHECK %d = %d -->  %d  ( %lf %lf )\n ",
				//		     k, star -> quads[k], stat,  qm -> MINANGLE, qm -> MAXANGLE);
				if ( stat == EGADS_SUCCESS ) {
					for ( kk = 0 ; kk < k; kk++ ) {
						stat = checkInvalidElement(qm, star -> quads[kk], MINVALIDANGLE, MAXVALIDANGLE);
						//    printf(" -+-+-+-+-+-+   INVALID CHECK %d = %d -->  %d  ( %lf %lf )\n ",
						//     kk, star -> quads[kk], stat, MINVALIDANGLE, MAXVALIDANGLE);
						if ( stat != EGADS_SUCCESS ) break;
					}
					if ( stat == EGADS_SUCCESS )  q++;
				}
			}
		}
		if ( stat != EGADS_SUCCESS || q == 0 ) break;
		snprintf(buffer,500, "wvsCentroid%d_%i.txt", it, qm -> fID );
		stat = EG_wvsData(qm -> mesh , buffer);
		qm -> MINANGLE += 2.0 * PI/ 180.0;
		qm -> MAXANGLE -= 2.0 * PI/ 180.0;
		//	fprintf(stderr, " AT AVERAGING MIN MAX ANGLES %lf  %lf  \n ", qm ->MINANGLE, qm ->MAXANGLE);
		if ( fullRegularization == 0  ) break;
		if ( qm -> MAXANGLE < qm -> MINANGLE  ) break;
		if ( fullRegularization == 2 && it == 10 ) break;
		printMesh(qm, buffer,0);
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
	EG_free(vID);
	EG_freeStar (&star);
	printf(" LEAVE OPTIMIZE \n ");
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
	if (argc != 5 ) {
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
				f = 3;
				printMeshStats ( bodydata[iBody].qm[f] , -1 );
				stat = EG_fullMeshRegularization(bodydata[iBody].qm[f] );
				fprintf(stderr, " EG_fullMeshRegularization face %d / %d = %d \n ", f + 1, bodydata[iBody].nfaces,  stat );
				if ( stat == FATAL_ERROR ) return stat;
				snprintf(buffer,500,"FINALMESH_%d", f + 1);
				printMesh(bodydata[iBody].qm[f],buffer, 1);
				printMeshStats (  bodydata[iBody].qm[f] , -1 );
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
	printf(" Vertex %d is type %d and has valence %d \n ",
			id, mesh -> vType[ id - 1], getValence (mesh, id ));
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


static  void
printMesh(meshMap *qm, char *name, int usename ) {
	int i,k, v, d;
	double eval[18], average[5], dist;
	FILE *fout = NULL;
	if ( usename == 1 ) fout = fopen(name, "w");
	else    snprintf ( name, 100,"M_%d",MESHPLOTS++ ) ;
	printf(" Writing in File %s  ", name);
	fout = fopen(name, "w" );
	for ( i = 0 ; i < qm -> mesh -> totQuads; ++i) {
		if ( checkQuad ( qm -> mesh, i + 1) != EGADS_SUCCESS ) continue;
		for ( k = 0; k < 4; ++k) {
			v  =   qm -> mesh -> quadIdx[4*i + k] - 1;
			fprintf(fout, "%lf %lf %lf %d %lf %lf \n",qm -> mesh->xyzs[3*v], qm -> mesh->xyzs[3*v +1], qm -> mesh->xyzs[3*v + 2], v + 1,
					qm -> mesh -> uvs[2*v] , qm -> mesh -> uvs[2*v + 1]);
			dist = 0.0;
			EG_evaluate(qm -> face, &qm -> mesh -> uvs[2*v ], eval);
			for ( d = 0 ; d < 3; ++d) dist += ( eval[d] - qm -> mesh->xyzs[3*v + d]) * ( eval[d] - qm -> mesh->xyzs[3*v + d]);
			dist = sqrt  (dist);
			if ( dist > EPS06 ) {
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











