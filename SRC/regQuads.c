#ifdef STANDALONE
#include <time.h>
#define IOFILES
#endif

#include "egads.h"
#include "regQuads.h"


//#define DEBUG
//#define DEBUGG


/* BASE-LEVEL FUNCTIONS */
static void swapInt(int *a, int *b)
{
  int c;
  c  = *a;
  *a = *b;
  *b =  c;
}


static void unitVector(double *v)
{
  double n;
  
  n = sqrt (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  if (n > EPS11) {
    v[0] /=n; v[1] /=n; v[2] /=n;
  } else {
    v[0] = 0.0; v[1] = 0.0; v[2] = 0.0;
  }
}


static int inList(int n, int *list, int p)
{
  int i;

  if (n == 0) return -1;
  for (i = 0; i < n; i++)
    if (list[i] == p) return i;
  return -1;
}


/* IO FUNCTIONS */
#ifdef DEBUG
static void printVertexCoords(meshMap *qm, int v)
{
  v--;
  printf(" #vertex %d  ==============================================\n" , v+1);
  printf("%lf %lf %lf %lf %lf %d\n",
	 qm -> xyzs[3 * v    ], qm -> xyzs[3 * v + 1],
         qm -> xyzs[3 * v + 2], qm -> uvs [2 * v + 1],
         qm -> uvs [2 * v    ], v + 1);
  printf(" #==============================================\n");
}


static void printQuadCoords(meshMap *qm, int qID)
{
  int i = 0, v;
  double uv[3], uva[2], xyz[18];
  uv[0] = uv[1] = uv[2] = 0.0;
  printf(" #=================== quad %d ==================\n", qID);
  for (i = 0; i <= 4; i++) {
      v = qm -> qIdx[4 * (qID - 1) +  i%4] - 1;
      printf("%lf %lf %lf %d %lf %lf \n",
	     qm -> xyzs[3 * v    ], qm -> xyzs[3 * v + 1],
             qm -> xyzs[3 * v + 2], v + 1,
	     qm -> uvs [2 * v + 1], qm -> uvs [2 * v    ]);
      if (i < 4) {
	  uv[0] += 0.25 * qm -> xyzs[3 * v    ];
	  uv[1] += 0.25 * qm -> xyzs[3 * v + 1];
	  uv[2] += 0.25 * qm -> xyzs[3 * v + 2];
      }
  }
  i = EG_invEvaluate(qm -> face, uv, uva, xyz);
  printf("\n\n\n %lf %lf %lf %d %lf %lf \n",
         xyz[0], xyz[1], xyz[2], qID, uv[0], uv[1]);
  printf("\n\n #==============================================\n");
}
#endif


static int  checkMesh                  (meshMap *qm);
#ifdef IOFILES
static int  EG_wvsData                 (meshMap *qm, /*@null@*/ char *buffer);
static void printMesh                  (meshMap *qm, /*@null@*/ char *);
#endif
static void printStarFile( meshMap *qm, vStar *star, char *name);
#ifdef DEBUG

static void printStarFile              (meshMap *m, vStar *star, char *name);
#endif
static void printQuadSpecs             (meshMap *qm, int id);
static void printQuadGroup             (meshMap *qm, quadGroup qg);

/* GEOMETRY FUNCTIONS */
static int    EG_quadAngleOrientation  (meshMap *qm, double minT, double maxT,
                                        int qID, int *ori, int *order, double *theta,
                                        double *quadSize);
static int    EG_angleAtVnormalPlane   (meshMap *qm, int vC, int v1, int v2,
                                        double *angle);
static int    EG_angleAtBoundaryVertex (meshMap *qm, int v, int *links,
                                        double *size);
static int    EG_makePositiveAngles    (meshMap *qm, double minAngle,
                                        double maxAngle, int qID, int);
static void   rotationUV               (meshMap *qm, double theta,
                                        int counterClockWise, double *t0,
                                        double *t1);
static double EG_segment               (ego face,    double *uv1, double *uv2);
static int    EG_moveDistAway          (meshMap *qm, double alpha, int iA,
                                        int iB);
static void   EG_computeCoords         (meshMap *qm, int vC, double min,
                                        double max);
static int    EG_centroid              (meshMap *qm, int qID, double *quv,
                                        int type);
/***********************************************************/
/* MESH MAP FUNCTIONS*/
static int  EG_cleanQuad              (meshMap *qm, int qID,  int adj,
                                       int tansfer, int forcing, int *activity);
static int  setValence                (meshMap *qm, int vID);

static int  EG_buildStar              (meshMap *qm, vStar **star, int vID);
static void updateVertex              (meshMap *qm, int vID, double *uv);
static int  EG_swappingOperation      (meshMap *qm, quadGroup qg, int swap,
                                       int *activity);
static int  EG_splittingOperation     (meshMap *qm, int vC, int vL, int vR,
                                       int *activity);
static int  EG_mergeVertices          (meshMap *qm, int qC, int centre,
                                       int *activity);
static int  checkQuad                 (meshMap *qm, int q);
static int  checkVertex               (meshMap *qm, int v);
static int  EG_quadIsBoundary         (meshMap *qm, int qID);
static int  EG_adjQtoPair             (meshMap *qm, int qID, int v1, int v2,
                                       int *adj);
/**************/
static int  EG_swap                   (meshMap *qm, int qID, int *activity);
static int  EG_doubleSwap             (meshMap *qm, quadGroup qg, int forcing,
                                       int *activity);
static int  EG_collapse               (meshMap *qm, int  vID, int *activity,
                                       int forcing);
static int  EG_split                  (meshMap *qm, int  qID, int *activity);
static int  EG_swapSplit              (meshMap *qm, quadGroup qg, int forcing,
                                       int *activity);
static int  EG_swapCollapse           (meshMap *qm, quadGroup qg, int forcing,
                                       int *activity);
static int  EG_doubleCollapse         (meshMap *qm, quadGroup qg, int forcing,
                                       int *activity);
static int  EG_doubleSplit            (meshMap *qm, quadGroup qg, int forcing,
                                       int *activity);
static int  EG_swapDoubleSplit        (meshMap *qm, quadGroup qg, int *activity);
static int  EG_swapDoubleCollapse     (meshMap *qm, quadGroup qg, int *activity);
/**************/


static void meshCount(meshMap *qm, int *nI, int *nV, int *nQ)
{
  int i, qSum, vSum = 0, vSum2 = 0;

  for (vSum = i = qSum = 0; i < qm -> totQ; i++)
    if (qm -> qIdx[4 * i] != -2) qSum++;
  for (vSum = i = 0; i < qm -> totV; i++) {
      if (qm -> vType[i] == -2) continue;
      vSum2++;
      if (qm -> valence[i][1] != 4) vSum++;
  }
  *nQ = qSum;
  *nI = vSum;
  *nV = vSum2;
}


#ifdef IOFILES
static void printMeshStats(meshMap *qm, int sweep)
{
  int i,len, val, intVal[100], boundVal[100];
  char buffer[100];
  FILE *fout;
  
  if (sweep == 0)
    sprintf(buffer, "OriginalMeshStats_face_%d_.txt", qm -> fID);
  else if (sweep == 10)
    sprintf(buffer, "PreprocessMeshStats_face_%d_.txt", qm -> fID);
  else
    sprintf(buffer, "RegularlMeshStats_face_%d_.txt", qm -> fID);
  fout = fopen(buffer, "w");
  if (fout == NULL) return;
  len = qm ->totV;
  for (i = 0; i < 100; ++i) {
      intVal  [i] = 0;
      boundVal[i] = 0;
  }
  for (i = 0; i < len; ++i) {
      if (qm -> vType[i] !=  -2) {
	  if (qm -> vType[i] == -1) {
	      val = qm -> valence[i][1];
	      ++intVal[val];
	  } else {
	      ++boundVal[qm ->valence[i][1]];
	  }
      }
  }
  fprintf(fout,"---- TOTAL VERTICES %d TOTAL QUADS %d --------------\n",
	  qm -> totV - qm -> remV[0],
          qm -> totQ - qm -> remQ[0]);
  fprintf(fout," INTERIOR VERTICES\n");
  for (i = 0; i < 100; ++i) {
      if (intVal[i]  > 0)
        fprintf(fout," VALENCE %d = %d VERTICES\n", i, intVal[i]);
  }
  fprintf(fout," BOUNDARY VERTICES\n");
  for (i = 0; i < 100; ++i) {
      if (boundVal[i]  > 0)
        fprintf(fout," VALENCE %d = %d VERTICES\n", i, boundVal[i]);
  }
  meshCount(qm, &i, &len, &val);
  fprintf(fout,"\n\nMesh has %d QUADS, %d / %d irregular vertices (%.2f %% ) ============\n",
          val, i, len, (double) i * 100.0 / (double) len);
  fclose(fout);
}
#endif


static void EG_freeStar(vStar **star)
{
  if ((*star) == NULL) return;
  EG_free((*star) -> verts);
  EG_free((*star) -> quads);
  EG_free((*star) -> idxQ);
  EG_free((*star) -> idxV);
  EG_free(*star);
  *star = NULL;
}


static int EG_backupQuads(meshMap *qm, int *nq, int *qlist, Quad **quad)
{
  int   i, j, q, v, qcount, *qaux = NULL, k, add = 0, stat;
  vStar *star = NULL;
  
  qaux       = EG_alloc(qm -> totQ * sizeof(int));
  if  (qaux == NULL) return EGADS_MALLOC;
  for (qcount = q = 0; q < *nq; q++) {
      if (qlist[q] == -1) continue;
      for (i = 0; i < 4; i++) {
	  v    = qm -> qIdx[4 * (qlist[q] - 1) + i];
	  stat = EG_buildStar(qm, &star, v);
	  if (star == NULL || stat != EGADS_SUCCESS) {
	      EG_free (qaux);
	      return stat;
	  }
	  for (j = 0; j < star -> nQ; j++) {
	      if (star-> quads[j] == -1) continue;
	      for (add = k = 0; k < qcount; k++) {
		  if (qaux[k] == star -> quads[j]) {
		      add = 1;
		      break;
		  }
	      }
	      if (add == 0 || qcount == 0) qaux[qcount++] = star ->quads[j];
	  }
      }
  }
  EG_freeStar(&star);
  (*quad)      = EG_alloc(qcount * sizeof(Quad));
  if ((*quad) == NULL) {
      EG_free(qaux);
      return EGADS_MALLOC;
  }
  for (q = 0; q < qcount; q++) {
      (*quad)[q].id = qaux[q];
      for (j = 0; j < 4; j++) {
	  (*quad)[q].qadj [j] = qm -> qAdj[4 * (qaux[q] - 1) + j];
	  (*quad)[q].verts[j] = qm -> qIdx[4 * (qaux[q] - 1) + j];
      }
  }
  EG_free (qaux);
  *nq = qcount;
  return EGADS_SUCCESS;
}


static int EG_restoreQuads(meshMap *qm, Quad *quad, int nq)
{
  int i, j, *vid = NULL, k, stat;
  
  vid      = EG_alloc(qm -> totV * sizeof(int));
  if (vid == NULL) return EGADS_MALLOC;
  for (k  = i = 0; i < nq; i++) {
      if (quad[i].id == -1) continue;
      for (j = 0; j < 4; j++) {
	  qm -> qAdj[4 * (quad[i].id - 1) + j] = quad[i].qadj [j];
	  qm -> qIdx[4 * (quad[i].id - 1) + j] = quad[i].verts[j];
	  if (k == qm -> totV) continue;
	  if (inList(k, vid, quad[i].verts[j]) == -1) {
	      qm -> valence[quad[i].verts[j] - 1][0] = quad[i].id;
	      vid[k++] = quad[i].verts[j];
	  }
      }
  }
  for (i   = 0; i < k; i++) {
      stat = setValence(qm, vid[i]);
      if (stat != EGADS_SUCCESS) {
	  EG_free (vid);
	  return stat;
      }
  }
  EG_free (vid);
  
  return EGADS_SUCCESS;
}


int EG_createMeshMap(bodyQuad *bodydata, int uvtype)
{
  int          f, iA, iB, stat = 0, stat2 = 0, j, q, i, auxID, k, kk, kOK, len, ne;
  int          ntri, nquad, e4[4], ori[8];
  const int    *tris, *tric, *ptype, *pindex;
  double       angle,  minQarea, averagesize, surfarea, angles[4], qsize;
  const double *xyzs, *uvs;
  int          qV[6]    = { 0, 1, 2, 5, 0, 1};
  int          qLoop[5] = { 0, 1, 2, 3, 0   };
  vStar        *star = NULL;

  bodydata->qm = (meshMap **) EG_alloc(bodydata->nfaces * sizeof(meshMap *));
  if (bodydata->qm == NULL) return  EGADS_MALLOC;
  for (f = 0; f < bodydata->nfaces; ++f) {
      bodydata ->qm[f] = (meshMap *) EG_alloc(sizeof(meshMap));
      if (bodydata ->qm[f] == NULL) {
	  printf("Create Quad Map: MALLOC error!! \n ");
	  stat2 = EGADS_MALLOC;
	  continue;
      }
      bodydata -> qm[f] -> fID       = f + 1;
      bodydata -> qm[f] -> pp        = 0;
      bodydata -> qm[f] -> plotcount = 0;
      bodydata -> qm[f] -> uvtype    = uvtype;
      /* Edges associated to face */
      stat = EG_getTessFace(bodydata->tess, f + 1, &len,
                            &xyzs, &uvs, &ptype, &pindex, &ntri,
                            &tris, &tric);
      if (stat != EGADS_SUCCESS) {
	  printf("EG_createMeshMap :: EG_getTessFace %d !!\n", stat);
	  stat2 = stat;
	  continue;
      }
      nquad = (int) ntri/2;
      bodydata -> qm[f] -> oriV    = len;
      bodydata -> qm[f] -> oriQ    = nquad;
      bodydata -> qm[f] -> sizeV   = 2 * len;
      bodydata -> qm[f] -> sizeQ   = 2 * nquad;
      bodydata -> qm[f] -> totV    = len;
      bodydata -> qm[f] -> totQ    = nquad;
      bodydata -> qm[f] -> xyzs    = (double *) EG_alloc(3*(2 * len  )*sizeof(double ));
      bodydata -> qm[f] -> uvs     = (double *) EG_alloc(2*(2 * len  )*sizeof(double ));
      bodydata -> qm[f] -> vType   = (int    *) EG_alloc(  (2 * len  )*sizeof(   int ));
      bodydata -> qm[f] -> qIdx    = (int    *) EG_alloc(4*(2 * nquad)*sizeof(   int ));
      bodydata -> qm[f] -> qAdj    = (int    *) EG_alloc(4*(2 * nquad)*sizeof(   int ));
      bodydata -> qm[f] -> remQ    = (int    *) EG_alloc(  (2 * nquad)*sizeof(   int ));
      bodydata -> qm[f] -> remV    = (int    *) EG_alloc(  (2 * len  )*sizeof(   int ));
      bodydata -> qm[f] -> valence = (int   **) EG_alloc(  (2 * len  )*sizeof(   int *));
      bodydata -> qm[f] -> vFix    = (int    *) EG_alloc(  (2 * len  )*sizeof(   int ));
      if (bodydata -> qm[f]->qIdx  == NULL || bodydata -> qm[f]->qAdj    == NULL ||
	  bodydata -> qm[f]->xyzs  == NULL || bodydata -> qm[f]->uvs     == NULL ||
	  bodydata -> qm[f]->vType == NULL || bodydata -> qm[f]->remQ    == NULL ||
	  bodydata -> qm[f]->remV  == NULL || bodydata -> qm[f]->valence == NULL ||
	  bodydata -> qm[f]->vFix  == NULL) {
	  stat2 = EGADS_MALLOC;
	  continue;
      }
      bodydata -> qm[f] -> remQ[0] = 0;
      bodydata -> qm[f] -> remV[0] = 0;
      bodydata -> qm[f] -> invsteps = 0;
      for (j = 0; j < 2 * len; j++) {
	  bodydata -> qm[f] -> valence[j] = (int *) EG_alloc(len * sizeof(int));
	  if (bodydata -> qm[f] -> valence[j] == NULL) {
	      stat2 = EGADS_MALLOC;
	      continue;
	  }
      }
      bodydata->qm[f] ->vFix[0] = 0;
      bodydata->qm[f] ->face    = bodydata->faces[f];
      stat = EG_getRange (bodydata->qm[f] ->face, bodydata->qm[f]->range, &j);
      for (j = 0; j < len; j++) {
	  bodydata->qm[f] ->valence[j][0]      = 0; //quad associated with vertex
	  bodydata->qm[f] ->valence[j][1]      = 0; // adjusted valence
	  bodydata->qm[f] ->valence[j][2]      = 0; // real valence ( [3], [4],.. are the vertices linked)
	  bodydata->qm[f] ->vType  [j]         = ptype[j];
	  bodydata->qm[f] ->uvs    [2 * j    ] = uvs[2 * j    ];
	  bodydata->qm[f] ->uvs    [2 * j + 1] = uvs[2 * j + 1];
	  if (uvs[2 * j    ] < bodydata ->qm[f] -> range[0] - 1.e-6 ||
	      uvs[2 * j    ] > bodydata ->qm[f] -> range[1] + 1.e-6||
	      uvs[2 * j + 1] < bodydata ->qm[f] -> range[2] - 1.e-6||
	      uvs[2 * j + 1] > bodydata ->qm[f] -> range[3] + 1.e-6) {
	      printf(" VERTEX %d = %f  %f  is out of range !! \n ",
		     j + 1, uvs[2 * j], uvs[2 * j + 1]);
	      printf(" RANGE FACE %d --- > %f  %f  %f  %f\n ", bodydata ->qm[f] -> fID,
		     bodydata ->qm[f] -> range[0],bodydata ->qm[f] -> range[1],
		     bodydata ->qm[f] -> range[2],bodydata ->qm[f] -> range[3]);
	      stat2 = EGADS_INDEXERR;
	      break;
	  }
	  bodydata->qm[f] ->xyzs[3 * j    ] = xyzs[3 * j    ];
	  bodydata->qm[f] ->xyzs[3 * j + 1] = xyzs[3 * j + 1];
	  bodydata->qm[f] ->xyzs[3 * j + 2] = xyzs[3 * j + 2];
      }
      if (stat2 != EGADS_SUCCESS) continue;
      for (j = 0; j < nquad; j++)
	for (k = 0; k < 4; ++k) {
	    bodydata ->qm[f] ->qIdx[4*j + k] = tris[6*j + qV[k+1]];
	}
      for (j = 0; j < nquad; j++)
	{
	  kk  = 0;
	  kOK = 0;
	  q   = 0;
	  while (q < nquad) {
	      if (q == j) {
		  if (j == nquad - 1 && kk < 4) {
		      bodydata ->qm[f] ->qAdj[4*j + kk++] = -1;
		      if (kk == 4) break;
		      q = 0;
		  }
		  else q++;
	      }
	      if  (q == nquad) break;
	      for (k = 0; k < 4; ++k) {
		  if ((bodydata ->qm[f] ->qIdx[4*j + qLoop[kk    ]] ==
                       bodydata ->qm[f] ->qIdx[4*q + qLoop[k     ]] ||
		       bodydata ->qm[f] ->qIdx[4*j + qLoop[kk    ]] ==
                       bodydata ->qm[f] ->qIdx[4*q + qLoop[k  + 1]]) &&
		      (bodydata ->qm[f] ->qIdx[4*j + qLoop[kk + 1]] ==
                       bodydata ->qm[f] ->qIdx[4*q + qLoop[k     ]] ||
                       bodydata ->qm[f] ->qIdx[4*j + qLoop[kk + 1]] ==
                       bodydata ->qm[f] ->qIdx[4*q + qLoop[k  + 1]])) {
		       bodydata ->qm[f] ->qAdj[4*j + kk++] = q + 1;
		      q   = -1;
		      kOK =  1;
		      k   =  4;
		      if (kk == 4) q = nquad;
		    }
	      }
	      if ((kOK == 0) && (q >= nquad -1)) {
		  bodydata ->qm[f] -> qAdj[4*j + kk++] = -1;
		  q = -1;
		  if (kk == 4) break;
	      }
	      else kOK = 0;
	      q++;
	    }
	  if (kOK == 0 && kk < 4) {
	      while (kk < 4)
		bodydata ->qm[f] ->qAdj[4*j + kk++] = -1;
	  }
      }
      for (j = 0; j < nquad; j++)
	    for (q = 0; q < 4; ++q)
	    	bodydata->qm[f] ->valence[bodydata->qm[f] -> qIdx[4 * j + q] - 1][0] = j + 1;
      bodydata ->qm[f] ->minsize  = 1000000000.0;
      bodydata ->qm[f] ->minAngle = DEG10 * 0.1;
      bodydata ->qm[f] ->maxAngle =          PI;
      averagesize                 = 0.0;
      for (j = 0; j < len; j++) {
	  stat = setValence (bodydata ->qm[f], j + 1);
	  if (stat != EGADS_SUCCESS) {
	      printf("In EG_createMeshMap :: set valence at %d is %d!!\n ", j + 1, stat);
	      stat2 = stat;
	      break;
	  }
      }
      for (ne = j = 0; j < len; j++) {
	  if (bodydata ->qm[f] ->vType[j] == -1) continue;
	  k = bodydata ->qm[f] ->valence[j][2];
	  for (i    = 0; i < k; i++) {
	      auxID = bodydata ->qm[f] ->valence[j][3 + i] - 1;
	      if (bodydata ->qm[f] ->vType[auxID] == -1) continue;
	      angle = EG_segment( bodydata ->qm[f] ->face,
                                 &bodydata ->qm[f] ->uvs[2*j],
                                 &bodydata ->qm[f] ->uvs[2*auxID]);
	      bodydata ->qm[f] ->minsize = MIN(angle, bodydata ->qm[f] ->minsize);
	      averagesize += angle;
	      ne++;
	  }
	  stat = EG_angleAtBoundaryVertex(bodydata ->qm[f], j + 1, e4, &angle);
	  if (stat != EGADS_SUCCESS || angle < EPS11) {
	      stat2 = EGADS_GEOMERR;
	      printf(" Stat in EG_angleAtBoundaryVertex %d angle %f\n ", stat, angle);
	      printf(" Vertices: %d %d %d \n ", j+ 1, e4[0], e4[1]);
	      break;
	  }
	  else if (angle < 0.85 * PI) bodydata ->qm[f] ->vType[j] = 2;
	  else if (angle < 1.25 * PI) bodydata ->qm[f] ->vType[j] = 3;
	  else if (angle < 1.85 * PI) bodydata ->qm[f] ->vType[j] = 4;
	  else                        bodydata ->qm[f] ->vType[j] = 5;
	  if (bodydata -> qm[f] -> valence[j][2] < bodydata ->qm[f] ->vType[j])
	    bodydata ->qm[f] ->vType[j] = bodydata -> qm[f] -> valence[j][2];
	  if      (bodydata ->qm[f]-> vType[j] == 2)
	    bodydata ->qm[f] -> valence[j][1] = bodydata ->qm[f] -> valence[j][2] + 2; //adjusted based on vertex type, eg, a corner is regular when val = 2
	  else if (bodydata ->qm[f]-> vType[j] == 3)
	    bodydata ->qm[f] -> valence[j][1] = bodydata ->qm[f] -> valence[j][2] + 1;
	  else if (bodydata ->qm[f]-> vType[j] >= 5)
	    bodydata ->qm[f] -> valence[j][1] = bodydata ->qm[f] -> valence[j][2] - 1;
	  stat = EG_buildStar(bodydata->qm[f], &star, j + 1);
	  if (stat != EGADS_SUCCESS || star == NULL) {
	      stat2 = stat;
	      continue;
	  }
	  for (k = 0; k < star -> nQ; k++) {
	      if (star -> quads[k] == -1) continue;
	      iA   = star -> verts[2 * k + 1];
	      iB   = star -> verts[star -> idxV[2 * k + 3]];
	      stat =  EG_angleAtVnormalPlane(bodydata -> qm[f], j + 1, iA, iB, &angle);
	      if (stat != EGADS_SUCCESS) continue;
	      bodydata ->qm[f] ->minAngle = MIN(angle, bodydata ->qm[f] ->minAngle);
	      bodydata ->qm[f] ->maxAngle = MAX(angle, bodydata ->qm[f] ->maxAngle);
	  }
	  EG_freeStar(&star);
      }
      if (stat2 != EGADS_SUCCESS) continue;
      // ESTIMATE AREA AND TOTAL QUADS
      surfarea = 0.0;
      for (q = 0; q < nquad; q++) {
    	  stat      = EG_quadAngleOrientation(bodydata -> qm[f], 0.0, PI, q + 1,
                                              ori, &ori[4], angles, &qsize);
    	  surfarea += qsize;
      }
      averagesize /= (double)ne;
      bodydata ->qm[f] ->avsize = averagesize;
      minQarea     = bodydata ->qm[f] ->minsize * bodydata ->qm[f] ->minsize;
      minQarea     = averagesize * averagesize;
      bodydata -> qm[f] -> estQ = ceil(surfarea / minQarea);
      stat = checkMesh(bodydata->qm[f]);

#ifdef DEBUG
      printf(" Min max angle at bounds %lf  %lf ( 10 DEG %lf  160 DEG %lf ) \n ",
             bodydata ->qm[f] ->minAngle, bodydata ->qm[f] ->maxAngle, DEG10, DEG160);
#endif
      if (stat != EGADS_SUCCESS) {
	  printf("In EG_createMeshMap :: checkMesh at face %d --> %d!!\n",
                 f + 1, stat);
	  stat2 = stat;
      }
  }
  return stat2;
}


static int EG_projectToTangentPlane(double normal[], double *nO, double *p,
                                    double *proj)
{
  double c, dotNN = 0.0, dotNP = 0.0, dist, lambda;
  
  c       = DOT(normal, nO); // Equation plane: a*x + b*y + c*z = C
  dotNN   = DOT(normal, normal);
  dotNP   = DOT(normal, p);
  dist    = (p[0] - nO[0]) * (p[0] - nO[0]) + (p[1] - nO[1]) * (p[1] - nO[1]) +
            (p[2] - nO[2]) * (p[2] - nO[2]);
  if (fabs(dotNP - c) <= EPS11 || sqrt(dist) < EPS11) {
      proj[0] = p[0]; proj[1] = p[1]; proj[2] = p[2];
      return EGADS_SUCCESS;
  }
  lambda  =   (c - dotNP) / dotNN;
  proj[0] = p[0] + lambda * normal[0];
  proj[1] = p[1] + lambda * normal[1];
  proj[2] = p[2] + lambda * normal[2];
  dist    = DOT (normal, proj);
  if (fabs(dist - c) < EPS11)
    return EGADS_SUCCESS;
  else{
      printf(" ORIGIN %lf %lf  %lf  NORMAL %lf %lf  %lf  TARGET %lf %lf %lf\n",
             nO[0], nO[1], nO[2], normal[0], normal[1], normal[2], p[0], p[1], p[2]);
      printf(" POINT SHOULD BELONG TO PLANE!!!!! %.16e ~= %.16e\n", dist, c);
      printf(" DOT NN %lf PN %lf LAMBDA %lf  POINT %lf %lf %lf\n",
             dotNN, dotNP, lambda, proj[0], proj[1], proj[2]);
      return EGADS_GEOMERR;
  }
}


static int EG_angleAtVnormalPlane(meshMap *qm, int vC, int v1, int v2,
                                  double *angle)
{
  int    stat;
  double dot1, dot2, xyz[18], normal[3], cross[3], proj1[3], proj2[3];
  double centre[18], v01[3], v02[3], quv[2], norm1, norm2;
  
  stat = EG_evaluate(qm -> face, &qm -> uvs[2 * (vC - 1)], centre);
  if (stat != EGADS_SUCCESS) return stat;
  dot1 = centre[3]*centre[3] + centre[4]*centre[4] + centre[5]*centre[5];
  dot2 = centre[6]*centre[6] + centre[7]*centre[7] + centre[8]*centre[8];
  if (dot1 < EPS11 || dot2 < EPS11) {
      stat = EG_centroid(qm, vC, quv, 1);
      if (stat != EGADS_SUCCESS) return stat;
      stat = EG_evaluate(qm -> face, quv, centre);
      if (stat != EGADS_SUCCESS) return stat;
      dot1 = centre[3]*centre[3] + centre[4]*centre[4] + centre[5]*centre[5];
      dot2 = centre[6]*centre[6] + centre[7]*centre[7] + centre[8]*centre[8];
      if (dot1 < EPS11 || dot2 < EPS11) return EGADS_GEOMERR;
  }
  v01[0] = centre[3]; v01[1] = centre[4]; v01[2] = centre[5];
  v02[0] = centre[6]; v02[1] = centre[7]; v02[2] = centre[8];
  if (qm -> face -> mtype == SREVERSE) {
      CROSS(v02, v01, normal);
  } else {
      CROSS(v01, v02, normal);
  }
  stat  = EG_evaluate(qm -> face, &qm -> uvs[2 * (v1 - 1)], xyz);
  stat += EG_projectToTangentPlane(normal, centre, xyz, proj1);
  stat += EG_evaluate(qm -> face, &qm -> uvs[2 * (v2 - 1)], xyz);
  stat += EG_projectToTangentPlane(normal, centre, xyz, proj2);
  if (stat != EGADS_SUCCESS) return stat;
  v01[0] = proj1[0] - centre[0];
  v01[1] = proj1[1] - centre[1];
  v01[2] = proj1[2] - centre[2];
  v02[0] = proj2[0] - centre[0];
  v02[1] = proj2[1] - centre[1];
  v02[2] = proj2[2] - centre[2];
  norm1  = DOT (v01, v01);
  norm2  = DOT (v02, v02);
  CROSS     (v01, v02, cross);
  dot1 = DOT(v01, v02) / sqrt ( norm1 * norm2 );
  dot2 = DOT(normal, cross);
  if      (dot1 >=  1.0) *angle = 0.0;
  else if (dot1 <= -1.0) *angle = PI;
  else                   *angle = acos(dot1);
  if (dot2 < 0)          *angle =  2.0 * PI - *angle;
  return EGADS_SUCCESS;
}


static int EG_angleAtBoundaryVertex(meshMap *qm, int v, int *links, double *size)
{
  int   i, j, k;
  vStar *star = NULL;
  
  *size = 0.0;
  if (qm -> vType[v - 1] < 0) return EGADS_INDEXERR;
  i = EG_buildStar(qm, &star, v);
  if (i != EGADS_SUCCESS || star == NULL) {
      printf(" Looking at corners: buildstar %d is %d \n ", v, i);
      return i;
  }
  for (links[0] = links[1] = k = i = 0; i < star -> nQ; i++) {
      j = star -> verts[2 * i + 1] - 1;
      if (qm -> vType[j] != -1) k++;
      if (star -> quads[i] == -1) {
	  links[1] = star -> verts[2 * i + 1];
	  links[0] = star -> verts[star -> idxV[2 * i + 3]];
	  continue;
      }
  }
  EG_freeStar(&star);
  if (k >= 3) {
      *size = PI * 0.5; // boundary vertex is connected to more than two bounds. Angle is fine
      return EGADS_SUCCESS;
  }
  else if (k != 2) {
      printf(" EG_angleAtBoundaryVertex:: vertex %d is at surface bounds and connected only to another boundary vertex !!\n ", v);
      return EGADS_GEOMERR;
  }
  else return EG_angleAtVnormalPlane(qm, v, links[0], links[1], &(*size));
}


/*
 * Returns <0 EGADS_ERROR
 * = 0 self intersects
 * = 1 Good area
 * = 2 contains obtuse angles ( if quad is boundary obtuse angles => area = 0 )
 */



static int EG_quadAngleOrientation(meshMap *qm, double minT, double maxT,
                                   int qID, int *ori, int *order, double *theta, double *quadSize)
{
  int    i, qV[4], k, k1, k2, sign, stat, vA, vB, vC, area = 0, count = 0, doublet = 0, pos;
  int    qb, selfint = 1;
  double projABCD[12], cross[3], qNormal[3], vAB[3], vAC[3], quv[2], xyz[18];
  double signedAngle, norm1, norm2, norm3, angle, dot, size, c, dotNN, dotNP;
  qV[0]  = qm ->qIdx[4*(qID - 1)    ] - 1;
  qV[1]  = qm ->qIdx[4*(qID - 1) + 1] - 1;
  qV[2]  = qm ->qIdx[4*(qID - 1) + 2] - 1;
  qV[3]  = qm ->qIdx[4*(qID - 1) + 3] - 1;
  quv[0] = 0.25 * (qm -> uvs [ 2 * qV[0]    ] + qm -> uvs [ 2 * qV[1]] +
                   qm -> uvs [ 2 * qV[2]    ] + qm -> uvs [ 2 * qV[3]]);
  quv[1] = 0.25 * (qm -> uvs [ 2 * qV[0] + 1] + qm -> uvs [ 2 * qV[1] + 1] +
                   qm -> uvs [ 2 * qV[2] + 1] + qm -> uvs [ 2 * qV[3] + 1]);
  stat      = EG_evaluate(qm->face, quv, xyz);
  if (stat != EGADS_SUCCESS)
    {
 	  printf(" In EG_quadAngleOrientation:: EG_evaluate for quad %d  is %d !!\n",
                qID, stat);
 	  return stat;
    }
#ifdef DEBUG
  for ( k = 0; k < 4; k++ ) {
      c = 0.0;
      for ( k1 = 0; k1 < 4; k1++) {
	  if ( k1 == k ) continue;
	  if ( qV[k] == qV[k1] ) {
	      printf(" quad has repeated vertices !!\n ");
	      printQuadSpecs (qm , qID ) ;
	      exit (10);
	  }
	  c =  ( qm -> xyzs [ 3 * qV[k]]     - qm -> xyzs [ 3 * qV[k1]] ) *
	      ( qm -> xyzs [ 3 * qV[k]]     - qm -> xyzs [ 3 * qV[k1]] ) ;
	  c += ( qm -> xyzs [ 3 * qV[k] + 1] - qm -> xyzs [ 3 * qV[k1] + 1] ) *
	      ( qm -> xyzs [ 3 * qV[k] + 1] - qm -> xyzs [ 3 * qV[k1] + 1] ) ;
	  c += ( qm -> xyzs [ 3 * qV[k] + 2] - qm -> xyzs [ 3 * qV[k1] + 2] ) *
	      ( qm -> xyzs [ 3 * qV[k] + 2] - qm -> xyzs [ 3 * qV[k1] + 2] ) ;
	  if ( c < 1.e-08 ) {
	      printf(" quad %d as vertices %d = %d !!!!!\n ", qID,
		     qV[k] + 1, qV[k1] + 1);
	      return EGADS_GEOMERR;
	  }

      }
  }
#endif
  norm1 = xyz[3] * xyz[3] + xyz[4] * xyz[4] + xyz[5] * xyz[5];
  norm2 = xyz[6] * xyz[6] + xyz[7] * xyz[7] + xyz[8] * xyz[8];
  if (norm1 < EPS11 || norm2 < EPS11) {
	  stat       = EG_centroid(qm, qID, quv, 0);
	  stat      += EG_evaluate(qm->face, quv, xyz);
	  if (stat  != EGADS_SUCCESS) {
		  printf(" In EG_quadAngleOrientation:: EG_evaluate for quad %d  is %d !!\n",
				  qID, stat);
		  return stat;
	  }
	  norm1 = xyz[3] * xyz[3] + xyz[4] * xyz[4] + xyz[5] * xyz[5];
	  norm2 = xyz[6] * xyz[6] + xyz[7] * xyz[7] + xyz[8] * xyz[8];
	  if (norm1 < EPS11 || norm2 < EPS11) return EGADS_GEOMERR;
  }
  vAB[0] = xyz[3]; vAB[1] = xyz[4]; vAB[2] = xyz[5];
  vAC[0] = xyz[6]; vAC[1] = xyz[7]; vAC[2] = xyz[8];
  if (qm -> face -> mtype == SREVERSE) {
      CROSS(vAC, vAB, qNormal);
  } else {
      CROSS(vAB, vAC, qNormal);
  }
  norm1       = sqrt (qNormal[0] * qNormal[0] + qNormal[1] * qNormal[1] +
		      qNormal[2] * qNormal[2] );
  qNormal[0] /= norm1;
  qNormal[1] /= norm1;
  qNormal[2] /= norm1;
  c           = DOT(qNormal, xyz);
  dotNN       = DOT(qNormal, qNormal);
  for (i   = 0; i < 4; i++) {
	 dotNP               = qNormal[0] * qm -> xyzs[3 * qV[i]    ] +
                               qNormal[1] * qm -> xyzs[3 * qV[i] + 1] +
                               qNormal[2] * qm -> xyzs[3 * qV[i] + 2];
	 angle               = (c - dotNP) / dotNN;
	 projABCD[3 * i    ] = qm -> xyzs[3 * qV[i]    ] + angle * qNormal[0];
	 projABCD[3 * i + 1] = qm -> xyzs[3 * qV[i] + 1] + angle * qNormal[1];
	 projABCD[3 * i + 2] = qm -> xyzs[3 * qV[i] + 2] + angle * qNormal[2];
  }
  qb         = EG_quadIsBoundary(qm, qID);
  if (qb != 0 && maxT > PI) maxT = PI;
  for (pos   = doublet = count = sign = k = 0; k < 4; ++k) {
      ori[k] = 1;
      vA     = k;
      size   = 0.0;
      if (qm -> valence[qV[vA]][2] == 2 &&
          qm -> vType  [qV[vA]]    == -1) {
	  ori[k]      = 1;
	  theta[k]    = 0.5 * PI;
	  signedAngle = theta[k];
	  doublet     = 1;
      } else {
	  theta[k]    = 0.0;
	  signedAngle = 0.0;
	  for (k1 = 1; k1 < 3; k1++) {
	      vB  = (k + k1    )%4;
	      vC  = (k + k1 + 1)%4;
	      vAB[0] = projABCD[3 * vB    ] - projABCD[3 * vA    ];
	      vAB[1] = projABCD[3 * vB + 1] - projABCD[3 * vA + 1];
	      vAB[2] = projABCD[3 * vB + 2] - projABCD[3 * vA + 2];
	      vAC[0] = projABCD[3 * vC    ] - projABCD[3 * vA    ];
	      vAC[1] = projABCD[3 * vC + 1] - projABCD[3 * vA + 1];
	      vAC[2] = projABCD[3 * vC + 2] - projABCD[3 * vA + 2];
	      CROSS(vAB, vAC, cross);
	      norm1  = vAB[0] * vAB[0] + vAB[1] * vAB[1] + vAB[2] * vAB[2];
	      norm2  = vAC[0] * vAC[0] + vAC[1] * vAC[1] + vAC[2] * vAC[2];
	      norm3  = sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
	      size  += 0.5 * norm3;
	      dot    = DOT(vAB, vAC) / sqrt(norm1 * norm2);
	      if      (fabs(dot - 1.0) < EPS11) angle = 0.0;
	      else if (fabs(dot + 1.0) < EPS11) angle = PI;
	      else                              angle = acos(dot);
	      if (DOT(qNormal, cross) < 0) {
		  signedAngle += (2.0 * PI - angle);
		  ori[k] = -1;
	      }
	      else signedAngle += angle;
#ifdef DEBUGG
	      printf(" Angle %d %d %d ---> %lf (SIGNED %lf ) size %lf ORI %d \n ",
		     qV[vB]+1, qV[vA]+1, qV[vC]+1, angle, signedAngle, size, ori[k]);
#endif
	      theta[k] += angle;
	  }
      }
      order[k]  = k;
      if      (ori[k] == -1                  ) sign = 1;
      else if (ori[k] ==  1                  ) {
	  selfint = 0;
	  if (pos == 0 && fabs(size) > 1.e-08) {
	    *quadSize = fabs(size);
	    pos       = 1;
	  }
	  else if (fabs(size            ) > 1.e-08 &&
	           fabs(size - *quadSize) > 1.e-08) {
	      printf("---------------- Areas :: %.16e  %.16e  !!\n ",
                     size, *quadSize);
	  }
      }
      if (signedAngle < maxT + 1.e-08 && signedAngle > minT - 1.e-08) count++;
      else {
#ifdef DEBUGG
	  printf(" angle %f max %f min %f\n ", signedAngle, maxT, minT);
#endif
      }
  }
  if      (count                                == 4) area = 1;
  else if ((doublet == 1 || qb == 0 || qm -> pp == 1) && selfint == 0) area = 2;
  else area = 0;
  for  (k = 0; k < 3; k++) {
      for  (i = k + 1; i < 4; i++) {
	  k1   = 0;
	  if (sign == 0) {
	      if (theta[order[i]] > theta[order[k]]) k1 = 1;
	  } else {
	      if (theta[order[i]] * (double)ori[order[i]] <
                  theta[order[k]] * (double)ori[order[k]])
		k1 = 1;
	  }
	  if (k1 == 0) continue;
	  k2       = order[i];
	  order[i] = order[k];
	  order[k] = k2;
      }
  }
#ifdef DEBUGG
  if (area != 1) {
      printf(" ************   ATENCION:: INVALID QUAD  %lf  %lf  ( AREA %d )  ***************\n",
             minT, maxT, area);
      printf("------------ AREA QUAD %d  IS %d Internal angles ordered by size --------------\n",
             qID, area);
      for (sign = k = 0; k < 4; ++k)
	printf("Vertex %d has angle %f and orientation %d \n ",
               qV[order[k]] + 1, theta[order[k]], ori[order[k]]);
      //printQuadCoords(qm , qID);
  }
#endif
  return area;
}


/* clockwise = -1 , counterclockwise = 1
 * reverse = 0 assumes invalid quads e.g., ( AB x AC ) * normal < 0 so rotated B towards C (clockwise ) and
 * C towards B ( counterclockwise).if <ABC is small, it will increase the angle.
 * reverse = 1 assumes <ABC is big and will move B counter-clockwise and C clockwise
 */
static void EG_3Drotation(meshMap *qm, double minT, double maxT, double theta,
                          int iA, int iB, int iC)
{
  int    i, a, v,  area, counterClockWise = 1, stat, fix = 0, reverse = 0;
  int    piv[4], ori[4], links[2];
  double ang0, ang1, angles[4], uvf[2], angle, qsize;
  vStar  *star = NULL;
  
  stat = EG_angleAtVnormalPlane(qm, iA, iB, iC, &angle);
#ifdef DEBUGG
  printf(" LOOKING AT ANGLE %d %d  %d  --> %lf   MIN MAX %lf  %lf  \n ",
         iA, iB, iC, angle, minT, maxT);
#endif
  if (stat != EGADS_SUCCESS) return;
  if (angle > maxT) {
      if (qm -> vType[iA - 1] >= 0) {
	  i     = EG_angleAtBoundaryVertex(qm, iA, links,  &ang0);
	  stat  = EG_angleAtVnormalPlane  (qm, iA, links[0], iB, &ang0);
	  stat += EG_angleAtVnormalPlane  (qm, iA, links[0], iC, &ang1);
	  if (stat != EGADS_SUCCESS) return;
	  if (ang0 < ang1 || links[0] == iC) reverse = 1;
      } else if (angle < 1.5 * PI) reverse = 1;
  }
  else if (angle >= minT) return;
  ang0 = EG_segment( qm -> face, &qm -> uvs[2 * ( iA - 1)],
                    &qm -> uvs[2 * (iB - 1)]);
  ang1 = EG_segment( qm -> face, &qm -> uvs[2 * ( iA - 1)],
                    &qm -> uvs[2 * (iC - 1)]);
  if (     ang0 > ang1 && ang1 / ang0 < 0.95) EG_moveDistAway(qm, 0.95, iA, iB);
  else if (ang0 < ang1 && ang0 / ang1 < 0.95) EG_moveDistAway(qm, 0.95, iA, iC);
  if (reverse == 1) counterClockWise *= -1;
  for (a = 0; a < 2; a++) {
      counterClockWise *= -1;
      if (a == 0) v = iB - 1;
      else        v = iC - 1;
      fix = 0;
      if (qm -> vType[v] != -1) continue;
      if (inList(qm -> vFix[0], &qm -> vFix[1], v + 1) >= 0) {
	  stat      = EG_buildStar(qm, &star, v + 1);
	  if (stat != EGADS_SUCCESS || star == NULL) return;
	  fix     = 1;
	  uvf[0]  = qm -> uvs[2 * v    ];
	  uvf[1]  = qm -> uvs[2 * v + 1];
	  for (i  = 0; i < star -> nQ; i++) {
	      star->verts[i] = EG_quadAngleOrientation(qm, minT, maxT, star->quads[i],
                                                       ori, piv, angles, &qsize);
#ifdef DEBUGG
	      printf(" i %d -> %d AReA %d \n ",
                     i, star -> quads[i], star -> verts[i]);
#endif
	  }
      }
#ifdef DEBUGG
      printf(" ROTATE %d  %d   %d \n ", iA, v + 1, counterClockWise);
#endif
      rotationUV(qm, theta, counterClockWise, &qm -> uvs[2 * (iA - 1)],
                 &qm -> uvs[2 * v]);
      updateVertex(qm, v + 1, &qm -> uvs[2 * v]);
      if (fix == 1 && star != NULL) {
	  for (i = 0; i < star -> nQ; i++) {
	      area = EG_quadAngleOrientation(qm, minT, maxT, star -> quads[i],
                                             ori, piv, angles, &qsize);
	      if (area < 0 || (area == 0 && star -> verts[i] > 0)) {
#ifdef DEBUGG
		  printf (" ROTATED COOORDS MAKE INVALID QUAD \n ");
		  printVertexCoords(qm, v + 1);
#endif
		  updateVertex(qm, v + 1, uvf);
		  break;
	      }
	  }
	  EG_freeStar(&star);
      }
  }
  return;
}


static void rotationUV(meshMap *qm, double theta, int counterClockWise,
                       double *t0, double *t1)
{
  int    d, stat;
  double normal[4], proj[3], cross[4], xyz1[18], dt, duv[2], xyz0[18];
  double prevt[2], v1[3], v2[3], norm1, norm2;
  
  stat      = EG_evaluate(qm -> face, t0, xyz0);
  if (stat != EGADS_SUCCESS) {
#ifdef DEBUG
      printf("EG_rotationUV :: EG_normalToSurface ->%d!!\n",  stat);
#endif
      return;
  }
  norm1 = sqrt(xyz0[3] * xyz0[3] + xyz0[4] * xyz0[4] + xyz0[5] * xyz0[5]);
  norm2 = sqrt(xyz0[6] * xyz0[6] + xyz0[7] * xyz0[7] + xyz0[8] * xyz0[8]);
  if (norm1 < EPS11 || norm2 < EPS11) return;
  for (d = 0; d < 3; d++) {
      v1[d] = xyz0[3 + d] / norm1;
      v2[d] = xyz0[6 + d] / norm2;
  }
  if (qm -> face -> mtype == SREVERSE) {
      CROSS(v2, v1, normal);
  } else {
      CROSS(v1, v2, normal);
  }
  prevt[0]  = t1[0]; prevt[1]  = t1[1];
  stat      = EG_evaluate(qm -> face, t0, xyz0);
  stat     += EG_evaluate(qm -> face, t1, xyz1);
  if (stat != EGADS_SUCCESS) {
#ifdef DEBUG
      printf(" RotationUV rotation vector has EG_evaluate %d!!\n", stat);
#endif
      return;
  }
  stat      = EG_projectToTangentPlane(normal, xyz0, xyz1, proj);
  if (stat != EGADS_SUCCESS) {
#ifdef DEBUG
      printf("EG_rotateUV :: EG_projectToTangentPlane %d !!\n", stat);
#endif
      return;
  }
  v1[0] = proj[0] - xyz0[0];
  v1[1] = proj[1] - xyz0[1];
  v1[2] = proj[2] - xyz0[2];
  unitVector(v1);
  for (d = 0; d < 2; d++) {
      dt     = theta;
      if (d == 1) dt *= -1.0;
      duv[0] = (prevt[0] - t0[0]) * cos(dt) - (prevt[1] - t0[1]) * sin(dt);
      duv[1] = (prevt[0] - t0[0]) * sin(dt) + (prevt[1] - t0[1]) * cos(dt);
      t1[0]  = t0[0] + duv[0];
      t1[1]  = t0[1] + duv[1];
      if (t1[0] < qm ->range[0] || t1[0] > qm ->range[1] ||
	  t1[1] < qm ->range[2] || t1[1] > qm ->range[3]) {
	  if (d == 0) continue;
	  t1[0] = prevt[0]; t1[1] = prevt[1];
	  return;
      }
      stat      = EG_evaluate(qm -> face, t1, xyz1);
      if (stat != EGADS_SUCCESS) {
	  if (d == 0) continue;
	  t1[0] = prevt[0]; t1[1] = prevt[1];
	  return;
      }
      stat      = EG_projectToTangentPlane(normal, xyz0, xyz1, proj);
      if (stat != EGADS_SUCCESS) {
#ifdef DEBUG
	  printf("EG_rotateUV :: EG_projectToTangentPlane %d !!\n", stat);
#endif
	  t1[0] = prevt[0]; t1[1] = prevt[1];
	  return;
      }
      v2[0] = proj[0] - xyz0[0];
      v2[1] = proj[1] - xyz0[1];
      v2[2] = proj[2] - xyz0[2];
      unitVector(v2);
      CROSS(v1, v2, cross);
      if (DOT(normal, cross) * (double) counterClockWise > 0) return;
      else {
	  if (d == 1) {
	      t1[0] = prevt[0]; t1[1] = prevt[1];
	      return;
	  }
      }
  }
}


static double EG_segment(ego face, double *uv1, double *uv2)
{
  int    i, n = 10, stat;
  double p1[18], p2[18], dist[2], uvEps[2], seg = 0.0, totArc = 0.0, dt, range[4];

  i = EG_getRange(face, range,  &i);
  if (uv1[0] < range[0] - 1.e-08 || uv1[0] > range[1] + 1.e-08 ||
      uv1[1] < range[2] - 1.e-08 || uv1[1] > range[3] + 1.e-08) {
      printf("EG_segment 1 %lf %lf is out of bounds U ( %lf %lf ) V ( %lf %lf )!! \n ",
	     uv1[0], uv1[1], range[0], range[1], range[2], range[3]);
      return EGADS_RANGERR;
  }
  if (uv2[0] < range[0] - 1.e-08 || uv2[0] > range[1] + 1.e-08 ||
      uv2[1] < range[2] - 1.e-08 || uv2[1] > range[3] + 1.e-08) {
      printf("EG_segment 2 %lf %lf is out of bounds U ( %lf %lf ) V ( %lf %lf ) !! \n ",
	     uv2[0], uv2[1], range[0], range[1], range[2], range[3]);
      return EGADS_RANGERR;
  }
  stat       = EG_evaluate (face, uv1, p1);
  stat      += EG_evaluate (face, uv2, p2);
  if (stat != EGADS_SUCCESS) return -1.0;
  dist[0]    = uv2[0] - uv1[0];
  dist[1]    = uv2[1] - uv1[1];
  seg        = (p1[0] - p2[0]) * (p1[0] - p2[0]);
  seg       += (p1[1] - p2[1]) * (p1[1] - p2[1]);
  seg       += (p1[2] - p2[2]) * (p1[2] - p2[2]);
  for (i = 0; i < n; i++) {
      dt         = (double) (i + 1) / (double) n;
      uvEps[0]   = uv1[0] + dt * dist[0];
      uvEps[1]   = uv1[1] + dt * dist[1];
      stat       = EG_evaluate(face, uvEps, p2);
      if (stat != EGADS_SUCCESS) return -1.0;
      seg        = (p1[0] - p2[0]) * (p1[0] - p2[0]);
      seg       += (p1[1] - p2[1]) * (p1[1] - p2[1]);
      seg       += (p1[2] - p2[2]) * (p1[2] - p2[2]);
      p1[0]      = p2[0]; p1[1] = p2[1]; p1[2] = p2[2];
      totArc    += sqrt (seg);
  }
  return totArc;
}


static int EG_moveDistAway(meshMap *qm, double alpha, int iA, int iB)
{
  int    i, it;
  double seg, uvm[2], uvcopy[4], norm = 0.0, uv0[2], uv1[2], size, dt, segprev;

  if ( qm -> vType[iB -1] != -1 ||
      inList (qm -> vFix[0], &qm -> vFix[1], iB) >= 0) return EGADS_SUCCESS;
  uv0[0] = qm -> uvs[2 * (iA - 1)    ];
  uv0[1] = qm -> uvs[2 * (iA - 1) + 1];
  uv1[0] = qm -> uvs[2 * (iB - 1)    ];
  uv1[1] = qm -> uvs[2 * (iB - 1) + 1];
  seg    = EG_segment(qm -> face, uv0, uv1);
  if (seg < 0.0) return EGADS_GEOMERR;
  size   = seg * alpha;
  if (alpha * seg  < qm -> minsize * 0.25) size = 0.25 * qm -> minsize;
#ifdef DEBUGG
  printf(" Homotecia AB %d %d %lf --> %lf \n", iA, iB, seg, size);
#endif
  uvm   [0] = uv1[0]; uvm   [1] = uv1[1];
  uvcopy[0] = uv0[0]; uvcopy[1] = uv0[1];
  uvcopy[2] = uv1[0]; uvcopy[3] = uv1[1];
  it        = 0;
  dt        = 0.25;
  segprev   = seg;
  while (seg < size) {
      uvm[0] += (uvcopy[2] - uvcopy[0]) * dt;
      uvm[1] += (uvcopy[3] - uvcopy[1]) * dt;
      if (uvm[0] < qm -> range[0] || uvm[0] > qm -> range[1] ||
	  uvm[1] < qm -> range[2] || uvm[1] > qm -> range[3]) return EGADS_SUCCESS;
      uv1[0] = uvm[0]; uv1[1] = uvm[1];
      it++;
      seg    = EG_segment(qm -> face, uv0, uvm);
      if (it > 100) {
#ifdef DEBUG
	  printf(" Seg < size Stuck. Break !\n");
#endif
	  break;
      }
      if (fabs(segprev - seg) < size / 10.0) dt += 0.05;
      segprev = seg;
  }
  uvcopy[2] = uv1[0]; uvcopy[3] = uv1[1];
  it        = 0;
  while (fabs(seg - size) / size > 0.005) {
      uvm[0] = (uvcopy[0] + uvcopy[2]) * 0.5;
      uvm[1] = (uvcopy[1] + uvcopy[3]) * 0.5;
      seg    = EG_segment(qm -> face, uv0, uvm);
      i      = 0;
      if (seg > size) i = 1;
      uvcopy[2 * i    ] = uvm[0];
      uvcopy[2 * i + 1] = uvm[1];
      it++;
      norm  = (uvcopy[0] - uvcopy[2]) * (uvcopy[0] - uvcopy[2]);
      norm += (uvcopy[1] - uvcopy[3]) * (uvcopy[1] - uvcopy[3]);
      if (sqrt(norm) < EPS11 || it > 100) {
#ifdef DEBUG
	  printf(" POINTS ARE TOO CLOSE TO EACH OTHER: %f it %d  BISECTION DIDN't CONVERGE\n",
		 sqrt (norm), it);
#endif
	  break;
      }
  }
  updateVertex(qm, iB, uvm);
  return EGADS_SUCCESS;
}


static int EG_makeValidMesh(meshMap *qm, int nP, /*@null@*/ int *pList,
                            int fullReg)
{
  int    i, j, kq, kv, q, stat = EGADS_SUCCESS, it = 0, itMax,  ori[4], piv[4];
  int    *quads = NULL, *verts = NULL, recover = 0;
  double dthetam, dthetaM, minT, maxT, angles[4], *uvxyz = NULL, qsize;
#ifdef DEBUG
  char buffer[100];
#endif
  vStar *star   = NULL;

  uvxyz   = EG_alloc(5 * qm -> totV * sizeof(double));
  if (uvxyz == NULL) return EGADS_MALLOC;
  for (j = 0; j < qm -> totV; j++) {
      if (qm -> vType[j] != -1) continue;
      for (i = 0; i < 3; i++) {
	  if (i < 2) uvxyz[5 * j + i] = qm -> uvs[2 * j + i];
	  uvxyz[5 * j + 2 + i] = qm -> xyzs[3 * j + i];
      }
  }
  if (fullReg == 0) { // move around only affected vertices
      if (nP == 0 || pList == NULL) {
	  EG_free (uvxyz);
	  return EGADS_SUCCESS;
      }
      quads   = (int *) EG_alloc(qm -> totQ * sizeof(int));
      verts   = (int *) EG_alloc(qm -> totV * sizeof(int));
      if (quads == NULL || verts == NULL) {
          if (quads != NULL) EG_free(quads);
          if (verts != NULL) EG_free(verts);
	  EG_free(uvxyz);
	  return EGADS_MALLOC;
      }
#ifdef DEBUG
      printf(" Make valid mesh for %d points \n ", nP);
      printMesh(qm, NULL);
#endif
      for (kv = kq = j = 0; j < nP; j++) {
	  if (qm -> vType[pList[j] - 1] == -2) continue;
	  stat = EG_buildStar(qm, &star, pList[j]);
	  if (stat != EGADS_SUCCESS || star == NULL) {
	      EG_free(uvxyz);
	      EG_free(quads);
	      EG_free(verts);
	      return stat;
	  }
	  if (qm -> vType[star->verts[0] - 1] == -1 &&
	      inList(kv, verts, star->verts[0]) == -1) verts[kv++] = star -> verts[0];
	  for (i = 0; i < star -> nQ; i++) {
	      if (star -> quads[i] == -1) continue;
	      if (inList(kq, quads, star->quads[i]) == -1)
		quads[kq++] = star -> quads[i];
	      if (qm -> vType[star->verts[2 * i + 1] - 1] == -1 &&
		  inList(kv, verts, star->verts[2 * i + 1]) == -1)
		verts[kv++] = star -> verts[2 * i + 1];
	      if (qm -> vType[star->verts[2 * i + 2] - 1] == -1 &&
		  inList(kv, verts, star->verts[2 * i + 2]) == -1)
		verts[kv++] = star -> verts[2 * i + 2];
	  }
      }
      qm -> vFix[0] = 0;
      for (i = 0; i < kq; i++) {
	  if (EG_makePositiveAngles(qm, qm->minAngle, qm->maxAngle, quads[i], 0) <= 0) {
	      qm -> vFix[0] = 0;
	      recover       = 1;
	      break;
	  }
      }
      EG_free(quads);
      if (recover == 0) {
	  for (i = 0; i < kv; i++) {
	      stat = EG_buildStar(qm, &star, verts[i]);
	      if (stat != EGADS_SUCCESS || star == NULL) {
		  EG_free (uvxyz);
		  EG_free (verts);
		  return stat;
	      }
	      for (j = 0; j < star -> nQ; j++) {
		  if (star -> quads[j] == -1) continue;
		  if (EG_quadAngleOrientation(qm, qm->minAngle, qm->maxAngle,
					      star -> quads[j], ori, piv, angles, &qsize) <= 0)
		    {
		      recover = 1;
		      break;
		    }
	      }
	      if (recover == 1) break;
	  }
      }
      EG_freeStar(&star);
      EG_free(verts);
      if (recover == 1) {
#ifdef DEBUG
	  printf(" Mesh is invalid for MIN MAX TOTALS [0, 200 deg]\n");
	  snprintf(buffer, 100,"face_%d_InvalidMesh_%d",qm -> fID, qm -> plotcount);
	  printMesh(qm , buffer);
	  for (j = 0; j < nP; j++) {
	      if (qm -> vType[pList[j] - 1] == -2) continue;
	      stat = EG_buildStar(qm, &star, pList[j]);
	      if (stat != EGADS_SUCCESS || star == NULL) continue;
	      snprintf(buffer, 100,"face_%d_star_%d_%d",qm -> fID, pList[j], qm -> plotcount++);
	      printStarFile(qm, star, buffer);
	  }
	  EG_freeStar(&star);
#endif
	  stat = EGADS_GEOMERR;
	  for (j = 0; j < qm -> totV; j++) {
	      if (qm -> vType[j] != -1) continue;
	      for (i = 0; i < 3; i++) {
		  if (i < 2) qm -> uvs[2 * j + i] = uvxyz[5 * j + i];
		  qm -> xyzs[3 * j + i] = uvxyz[5 * j + 2 + i];
              }
	  }
	  qm -> invsteps++;
      }
      EG_free(uvxyz);
  } else {
      itMax   = 100;
      minT    = 0.0;
      maxT    = PI;
      dthetam = (DEG10 - minT) / (double)itMax;
      dthetaM = (PI - DEG160) / (double)itMax;
      qm -> vFix[0] = 0;
      for (it = 0; it < itMax; it++) {
#ifdef DEBUG
	  printf("SMOOTH ROUND %d max %d \n ", it, itMax);
#endif
	  for (q = 0; q < qm -> totQ; q++) {
	      if (qm -> qIdx[4 * q] == -2) continue;
	      stat = EG_makePositiveAngles(qm, minT, maxT, q + 1, 1);
#ifdef DEBUG
	      printMesh(qm, NULL);
#endif
	      if (stat < 0) {
		  EG_free(uvxyz);
		  return stat;
	      }
	  }
	  for (q = 0; q < qm -> totV; q++) {
	      if (qm -> vType[q] != -1) continue;
	      stat      = EG_buildStar(qm, &star, q + 1);
	      if (stat != EGADS_SUCCESS || star == NULL) continue;
	      for (j = i = 0; i < star -> nQ; i++)
		if (EG_quadAngleOrientation(qm, minT, maxT, star -> quads[i], ori,
					    piv, angles, &qsize) == 1) j++;
	      if (j == star -> nQ) {
		  for (i = 0; i < 3; i++) {
		      if (i < 2)
			uvxyz[5 * q     + i] = qm -> uvs [2 * q + i];
		      uvxyz  [5 * q + 2 + i] = qm -> xyzs[3 * q + i];
		  }
		  kv = inList(qm->vFix[0], &qm->vFix[1], q + 1);
		  if (kv != -1) {
		      qm -> vFix[0]--;
		      for (i = kv + 1; i < qm -> vFix[0]; i++) qm -> vFix[i] = qm -> vFix[i + 1];
		  }
	      } else {
		  for (i = 0; i < 3; i++) {
		      if (i < 2)
			qm -> uvs[2 * q + i] = uvxyz[5 * q     + i];
		      qm -> xyzs [3 * q + i] = uvxyz[5 * q + 2 + i];
		  }
		  for (j = i = 0; i < star -> nQ; i++)
		    if (EG_quadAngleOrientation(qm, minT, maxT, star -> quads[i], ori,
						piv, angles, &qsize) == 1) j++;
		  if (j == star -> nQ && inList(qm->vFix[0], &qm->vFix[1], q + 1) == -1) {
#ifdef DEBUG
		      printf(" makeValid :: fixing vertex %d \n ", q + 1);
#endif
		      qm -> vFix[++qm -> vFix[0]] = q + 1;
		  }
	      }
	  }
	  minT += dthetam;
	  maxT -= dthetaM;
	  if (fullReg != 1) break;
      }
      EG_free(uvxyz);
      EG_freeStar(&star);
      stat = EGADS_SUCCESS;
  }
#ifdef DEBUG
  printf(" LEAVING PEACEFULLY FROM OPTIMIZER ??? %d \n ", stat);
#endif
  qm -> vFix[0] = 0;
  return stat;
}


static int EG_makePositiveAngles(meshMap *qm, double minAngle, double maxAngle,
                                 int qID, int smooth)
{
  int  i, j, k, kk, iO, iA, iB, it , itMAX = 50, i0, area = 0, doublet = 0;
  int  ori[8], piv[8], quad[4], type[4], vfix[4], * sa = NULL;
  static int pivTri[6] = { 1, 2, 2, 3, 1, 3 };
  double ang, angles[8], dt = 0.01, ds = 0.95, xyzA[3], xyzB[3], dista, distb, maxT, minT;
  double uv[2], totDist = 0.0, qsize;
  vStar *star = NULL;
  
  for (j = 0; j < 4; j++) {
      quad[j] = qm -> qIdx[4 * (qID - 1) + j];
      type[j] = qm -> vType[quad[j] - 1];
      vfix[j] = 0;
      if (qm -> valence[quad[j] - 1][2] == 2 && type[j] == -1) doublet = 1;
  }
  if (maxAngle > PI) maxAngle = PI;
#ifdef DEBUGG
  printf("makePositiveAngle for quad %d \n ", qID);
#endif
  for (it  = 0; it < itMAX; it++) {
      area = EG_quadAngleOrientation(qm, minAngle, maxAngle, qID, ori, piv, angles, &qsize);
      if ((qm -> pp == 1 || doublet == 1) && area >= 1)  break;
      if (doublet == 0) {
	  for (j = 0; j < 4; j++) {
	      if (type[piv[j]] != -1 || vfix[piv[j]]  == 1) {
#ifdef DEBUG
		  printf(" V %d type %d vix %d \n ", j, type[piv[j]], vfix[piv[j]]);
#endif
		  continue;
	      }
	      k      = quad[piv[j]] - 1;
	      uv[0]  = qm -> uvs[2 * k    ];
	      uv[1]  = qm -> uvs[2 * k + 1];
	      i      = EG_buildStar(qm, &star, k + 1);
	      if (i != EGADS_SUCCESS || star == NULL) return EGADS_MALLOC;
	      sa     = EG_alloc(star -> nQ * sizeof(int));
	      if (sa == NULL) {
		  EG_freeStar(&star);
		  return EGADS_MALLOC;
	      }
	      for (i = 0; i < star -> nQ; i++) {
		  if (star -> quads[i] == -1 || star -> quads[i] == qID) continue;
		  sa[i] = EG_quadAngleOrientation(qm, minAngle, maxAngle, star->quads[i],
						  &ori[4], &piv[4], &angles[4], &qsize);
	      }
	      EG_computeCoords(qm, k + 1, minAngle, maxAngle);
	      for (i = 0; i < star -> nQ; i++) {
		  if (star -> quads[i] == -1 || star -> quads[i] == qID) continue;
		  area = EG_quadAngleOrientation(qm, minAngle, maxAngle, star->quads[i],
						 &ori[4], &piv[4], &angles[4], &qsize);
		  if ((smooth == 0 && area == 0 && sa[i] >= 1) ||
		      (smooth == 1 && area < sa[i])) {
#ifdef DEBUGG
		      printf(" At vertex %d !! invalid area quad %d -> %d\n ",
			     k + 1, star -> quads[i], area);
		      printf(" RESET COORDS %d \n ", k + 1);
		      printQuadCoords(qm, star -> quads[i]);
#endif
		      updateVertex(qm, k + 1, uv);
#ifdef DEBUGG
		      printQuadCoords(qm, star -> quads[i]);
#endif
		      break;
		  }
	      }
	      EG_freeStar(&star);
	      EG_free(sa);
	      sa = NULL;
	  }
	  area = EG_quadAngleOrientation(qm, minAngle, maxAngle, qID, ori, piv, angles, &qsize);
      }
      if (smooth == 1) {
	  for (j = i = 0; i < 4; i++) {
	      if      (angles[i] >= minAngle && angles[i] <= maxAngle) j++;
	      else if (type[i] == -1 && qm -> valence [quad[i]-1 ][2] != 4 &&
		       angles[i] > 0 && angles[i] < PI)  j++;
	  }
	  if (j == 4) area = 1;
      }
      if (area == 1 || (smooth == 0 && area == 2)) break;
#ifdef DEBUGG
      printf(" it %d area %d min %f max %f\n", it, area,  minAngle, maxAngle);
#endif
      for (i0 = j = 0; j < 4; j++) {
	  if (type[piv[j]] >= 0 && (type[piv[i0]] == -1 ||
	      angles[piv[j]] > angles[piv[i0]])) i0 = j;
      }
#ifdef DEBUGG
      printQuadCoords(qm, qID);
#endif
      for (j = 0; j < 4; j++) piv[j] = (piv[i0] + j)%4;
#ifdef DEBUG
      if (it%5 == 0) printMesh(qm, NULL);
#endif
      totDist = 0.0;
      for (j  = 0; j < 4; j++) {
	  iO  = piv[j];
	  if (type[iO] == -1  && ori[iO] == -1 &&
	      qm -> vType[qm -> valence[quad[iO] -1][3] -1 ] == -1) continue;
#ifdef DEBUGG
	  printf ("=========== O ( %d ) = %d ::: MAX ANGLE PASS %f %f   ==============================\n",
		  iO, quad[iO],  minAngle, maxAngle);
#endif
	  if (type[i0] >= 4 ||
	      (type[iO] == -1 && qm -> valence[quad[iO] - 1][2] != 4)) {
	      maxT = PI;
	      minT = 0.0;
	  }
	  else {
	      maxT = maxAngle;
	      minT = minAngle;
	  }
	  for (i = 0; i < 3; i++) {
	      iA = (iO + pivTri[2 * i    ])%4;
	      iB = (iO + pivTri[2 * i + 1])%4;
	      for (kk = 0; kk < 3; kk ++) {
		  xyzA[kk] = qm -> xyzs[3 * (quad[iA] - 1) + kk];
		  xyzB[kk] = qm -> xyzs[3 * (quad[iB] - 1) + kk];
	      }
#ifdef DEBUGG
	      printQuadCoords(qm, qID);
#endif
	      EG_3Drotation(qm, minT, maxT, dt,
                            quad[iO], quad[iA], quad[iB]);
	      dista = 0.0; distb = 0.0;
	      if (EG_angleAtVnormalPlane(qm, quad[iO], quad[iA], quad[iB],
                                         &ang) != EGADS_SUCCESS) continue;
	      for (kk = 0; kk < 3; kk++) {
		  dista += (xyzA[kk] - qm -> xyzs[3 * (quad[iA] - 1) + kk]) *
		           (xyzA[kk] - qm -> xyzs[3 * (quad[iA] - 1) + kk]);
		  distb += (xyzB[kk] - qm -> xyzs[3 * (quad[iB] - 1) + kk]) *
		           (xyzB[kk] - qm -> xyzs[3 * (quad[iB] - 1) + kk]);
	      }
	      dista = sqrt(dista);
	      distb = sqrt(distb);
	      if (dista + distb > 0.0 && dista < 0.1 * qm -> minsize &&
		  distb < 0.1 * qm -> minsize) {
		  dt += 0.01;
		  ds -= 0.01;
	      }
#ifdef DEBUGG
	      printf("it %d ----> dt %f ds %f distances %f  %f  :\n Angle i %d -> %d %d %d is now %f  \n ",
		      it, dt, ds, dista, distb, i, quad[iO], quad[iA], quad[iB], ang);
	      printQuadCoords(qm, qID);
#endif
	      totDist += dista + distb;
	      if ((ang > maxT || ang < minT) && dista + distb > 0.0) {
		  if (it == 0) {
		      if (dista > 0.0) vfix[iA] = 1;
		      if (distb > 0.0) vfix[iB] = 1;
#ifdef DEBUGG
		      printf(" --------- ATTENTION IN QUAD %d fixing angles %d %d \n ",
                             qID, quad[iA], quad[iB]);
#endif
		  }
		  break;
	      }
	  }
	  if ((ang > maxT || ang < minT) && dista + distb > 0.0) break;
      }
      if (totDist < EPS11) break;
  }
#ifdef DEBUGG
  printQuadCoords(qm, qID);
  printf("\n\n=============================================================\n");
  printf(" LEAVE MAKEPOSITIVE QUAD %d IS NOW AREA %d \n ", qID, area);
  printf("\n\n=============================================================\n");
#endif
  if (area == 1 || (area == 2 && smooth == 0)) {
      for (j = 0; j < 4; j++) {
	  if (vfix[j] == 0 || inList(qm->vFix[0], &qm->vFix[1], quad[j]) != -1) continue;
#ifdef DEBUGG
	  printf(" ATENTION: FIXING vert %d \n ", quad[j]);
#endif
	  qm -> vFix[++qm -> vFix[0] ] = quad[j];
      }
  }
  return area;
}


/*
 * type:  0 = quad , 1 = vertex
 */
static int EG_centroid(meshMap *qm, int id, double *quv, int type)
{
  int i, v, j, n;
  double xyz0[3], xyz1[3];
  
  quv[0] = quv[1] = xyz0[0] = xyz0[1] = xyz0[2] = 0.0;
  if (type == 0) n = 4;
  else           n = qm -> valence[id - 1][2] + 1;
  for (j = i = 0; i < n; i++) {
      if (type == 0)  v = qm ->qIdx[4*(id - 1) + i] - 1;
      else {
	  if (i == 0) v = id - 1;
	  else        v = qm ->valence[id - 1][3 + i - 1] - 1;
      }
      xyz0[0] += qm -> xyzs[3 * v    ];
      xyz0[1] += qm -> xyzs[3 * v + 1];
      xyz0[2] += qm -> xyzs[3 * v + 2];
      quv[0]  += qm -> uvs [2 * v    ];
      quv[1]  += qm -> uvs [2 * v + 1];
      j++;
  }
  xyz0[0] /= (double)j;
  xyz0[1] /= (double)j;
  xyz0[2] /= (double)j;
  quv [0] /= (double)j;
  quv [1] /= (double)j;
  i = EG_invEvaluateGuess (qm -> face, xyz0, quv, xyz1);
  if (i != EGADS_SUCCESS || quv[0] < qm -> range[0] || quv[0] > qm -> range[1] ||
      quv[1] < qm -> range[2] || quv[1] > qm -> range[3])
      i  = EG_invEvaluate (qm -> face, xyz0, quv, xyz1);
  if (i != EGADS_SUCCESS) {
      printf("EG_centroid :: EG_invEvaluate at %lf %lf %lf is %d!!\n ",
	     xyz0[0], xyz0[1], xyz0[2], i);
  }
  if      (quv[0] < qm -> range[0]) quv[0] = qm -> range[0];
  else if (quv[0] > qm -> range[1]) quv[0] = qm -> range[1];
  if      (quv[1] < qm -> range[2]) quv[1] = qm -> range[2];
  else if (quv[1] > qm -> range[3]) quv[1] = qm -> range[3];
  
  return EGADS_SUCCESS;
}


static void EG_computeCoords(meshMap *qm, int vID, double minAngle, double maxAngle)
{
  int    vA, vB, vC, i, j, k, auxID = 0, stat, i0, piv[4], ori[4];
  int    n = 0, *verts = NULL, bdry = 0 ;
  double uvc[2], angle, angles[4], totarc, arc, uva[2], qsize;
  vStar  *star = NULL;

  if (qm->vType[ vID -1] != -1 ||
      inList(qm->vFix[0], &qm->vFix[1], vID) >= 0) {
#ifdef DEBUG
      printf(" VERTEX %d is fix !\n ", vID);
#endif
      return;
  }
#ifdef DEBUG
  printf(" AVERAGE COORDS FOR VERTEX %d \n ", vID);
  printVertexCoords(qm, vID);
#endif
  if (qm->valence[vID - 1][2] == 2) {
        vA     = qm -> valence[vID -1][3] - 1;
        vB     = qm -> valence[vID -1][4] - 1;
        uva[0] = 0.5 * (qm -> uvs[2 * vA    ] + qm -> uvs[2 * vB    ]);
        uva[1] = 0.5 * (qm -> uvs[2 * vA + 1] + qm -> uvs[2 * vB + 1]);
        updateVertex(qm, vID, uva);
  #ifdef DEBUGG
        printf(" LEAVE WITH COORDS\n");
        printVertexCoords(qm, vID);
  #endif
        return;
    }
  stat  = EG_buildStar(qm, &star, vID);
  if (star == NULL || stat != EGADS_SUCCESS) return;
  verts = EG_alloc(star -> nQ * sizeof(int));
  if (verts == NULL) return;
  for (k = i = 0; i < star -> nQ; i++) {
      if (EG_quadAngleOrientation(qm, minAngle, maxAngle, star -> quads[i],
                                  ori, piv, angles, &qsize) <= 0) continue;
      vA         = star -> verts[             2 * i + 1];
      vB         = star -> verts[star -> idxV[2 * i + 3]];
      stat       = EG_angleAtVnormalPlane(qm, vID, vA, vB, &angle);
      if ( stat != EGADS_SUCCESS) {
	  printf(" Angle At normal plane %f -> %d \n ",angle, stat );
	  continue;
      }
#ifdef DEBUG
      printf(" Quad %d angle at normal plane %f ( %d %d %d ) \n ",
	     star -> quads[i], angle, vID, vA, vB);
#endif
      if      (angle > PI)  k =  1;
      else if (angle < 0.0) k = -1;
      else continue;
      for (n = j = 0; j < star -> nQ; j++) {
	  auxID  = star -> verts [ 2 * j + 1 ] - 1;
	  if ((k ==  1 && (auxID + 1 == vA || auxID + 1 == vB)) ||
	      (k == -1 &&  auxID + 1 != vA && auxID + 1 != vB)) continue;
	  verts[n++] = auxID;
	  break;
      }
  }
  vA = qm -> valence[vID - 1][3] -1;
  vB = qm -> valence[vID - 1][4] -1;
  vC = qm -> valence[vID - 1][5] -1;
  printf(" FIRST LINKS %d %d %d\n ", vA + 1, vB +1, vC + 1 ) ;
  if (k == 0               && qm -> vType[vC] == -1 &&
      (qm -> vType[vA] > 0 && qm -> vType[vA]  < 4) &&
      (qm -> vType[vB] > 0 && qm -> vType[vB]  < 4))
    {
      for (j    = i0 = k = 0; k < star -> nQ; k++) {
	  auxID = star -> verts [2 * k + 1] - 1;
	  if (qm -> vType[auxID] >= 0) {
	      if (i0 == 0) i0 = 2 * k + 1;
	      else {
		  j = 2 * k + 1;
		  if (j - i0 > 2 && (i0 + (star -> nV - 1) - j != 2)) {
		      i0 = star -> nV;
		      break;
		  }
	      }
	  }
      }
      k = 0;
      if (i0 == star -> nV) {
	  printf(" CENTRE POINT \n ");
	  verts[0] = vA + 1;
	  verts[1] = vB + 1;
	  n        = 2;
	  uva[0]   = 0.5 * ( qm -> uvs[2 * vA    ] + qm -> uvs[2 * vB    ] );
	  uva[1]   = 0.5 * ( qm -> uvs[2 * vA + 1] + qm -> uvs[2 * vB + 1] );
	  bdry     = 1;
      }
    }
  printf(" k = %d \n ", k ) ;
  if (k == 0) {
      for (j = 0; j < star -> nQ; j++) verts[j] = star -> verts[2 * j + 1] - 1;
      n = star -> nQ;
  }
  EG_freeStar(&star);
  if ( bdry != 1 ) {
      if (qm -> uvtype == INVEVAL && EG_centroid(qm, vID, uva, 1) != EGADS_SUCCESS) {
	  EG_free(verts);
	  return;
      } else {
	  if (n == 0) {
	      EG_free(verts);
	      return;
	  }
	  totarc = 0.0;
	  uva[0] = qm -> uvs[2 * (vID - 1)    ];
	  uva[1] = qm -> uvs[2 * (vID - 1) + 1];
	  uvc[0] = uvc[1] = 0.0;
	  for (j = 0; j < n; j++) {
	      arc  = EG_segment(qm->face, &qm ->uvs[2 * (vID - 1)],
				&qm ->uvs[2 * verts[j]]);
	      totarc += arc;
	      uvc[0] += arc * qm -> uvs[2 * verts[j]    ];
	      uvc[1] += arc * qm -> uvs[2 * verts[j] + 1];
	      uva[0] +=       qm -> uvs[2 * verts[j]    ];
	      uva[1] +=       qm -> uvs[2 * verts[j] + 1];
	  }
	  uva[0] /= (double)(n + 1);
	  uva[1] /= (double)(n + 1);
	  if (qm -> uvtype == ARCLENGTH && n > 1) {
	      uva[0] = uvc[0] / totarc;
	      uva[1] = uvc[1] / totarc;
	  }
      }
  }
  EG_free(verts);
  if (uva[0] >= qm -> range[0] && uva[0] <= qm -> range[1] &&
      uva[1] >= qm -> range[2] && uva[1] <= qm -> range[3]) updateVertex(qm, vID, uva);
#ifdef DEBUGG
  printf(" LEAVE WITH COORDS\n");
  printVertexCoords(qm, vID);
#endif
return;
}


/* v should be dimension 3: v[0] = n common verts: v[1] (v[2]) */
static void EG_commonVerts(meshMap *qm, int q1, int q2, int *v)
{
  int i, j, k;

  v[0] = 0;
#ifdef DEBUG
  if (checkQuad(qm, q1) != EGADS_SUCCESS || checkQuad(qm, q2) != EGADS_SUCCESS) {
      printf(" EG_commonVerts wrong quad pair %d %d\n", q1, q2);
      printQuadSpecs(qm, q1);
      printQuadSpecs(qm, q2);
      v[0] = -1;
      return;
  }
#endif
  if (q1 == q2) {
      for (k = j = 0; j < 4; j++)
	v[++k] = qm -> qIdx[4 * (q2 - 1) + j];
      return;
  }
  for (k = i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
	  if (qm -> qIdx[4 * (q1 - 1) + i] ==
	      qm -> qIdx[4 * (q2 - 1) + j])
	    v[++k] = qm -> qIdx[4 * (q2 - 1) + j];
      }
  }
  v[0] = k;
}


static int EG_quadVertIdx(meshMap *qm, int q, int v)
{
  int i = 0;

  for (i = 0; i < 4; i++)
    if (qm -> qIdx[4 * (q - 1) + i] == v) return i;

  return -1;
}


static int setValence(meshMap *qm, int vID)
{
  int   i, i4[10], i3[10], i1[100], j;
  vStar *star = NULL;

  if (vID <= 0 || vID > qm -> totV) return EGADS_INDEXERR;
  i      = EG_buildStar(qm, &star, vID);
  if (i != EGADS_SUCCESS || star == NULL) {
      printf(" In setValence for vertex %d stat from EG_buildStar = %d\n", vID, i);
      return i;
  }

  qm -> valence[vID - 1][2] = star -> nQ; // real valence
  qm -> valence[vID - 1][1] = qm -> valence[vID - 1][2];
  if (qm -> vType[vID - 1] == 2 ||
      qm -> vType[vID - 1] == 0) qm -> valence[vID - 1][1] = star -> nQ + 2; //adjusted based on vertex type, eg, a corner is regular when val = 2
  else if (qm -> vType[vID - 1] == 3) qm -> valence[vID -1][1] = star -> nQ + 1;
  else if (qm -> vType[vID - 1] >= 5) qm -> valence[vID -1][1] = star -> nQ  -1;
  for (i3[0] = i4[0] = i1[0] = i = 0; i < star -> nQ; i++) {
      if (star -> verts[2 * i + 1] == -1) continue;
      if (qm -> vType[star -> verts[2 * i + 1] -1 ] >= 4)
	i4[++i4[0]] = star -> verts[2 * i + 1];
      else if (qm -> vType[star -> verts[2 * i + 1] - 1] >= 2)
	i3[++i3[0]] = star -> verts[2 * i + 1];
      else
	i1[++i1[0]] = star -> verts[2 * i + 1];
  }
  for (j = i = 0; i < i4[0]; i++)
    qm -> valence[vID - 1][3 + (j++)] = i4[1 + i];
  for (i = 0; i < i3[0]; i++)
    qm -> valence[vID - 1][3 + (j++)] = i3[1 + i];
  for (i = 0; i < i1[0]; i++)
    qm -> valence[vID - 1][3 + (j++)] = i1[1 + i];
  EG_freeStar(&star);
  
  return EGADS_SUCCESS;
}


static int EG_nValenceCount(meshMap *qm, int q, int n)
{
  int i, count, val;

  for (count = i = 0; i < 4; i++) {
      val    = qm -> valence[qm ->qIdx[4 * (q - 1) + i] - 1][1];
      if (val < 0) return val;
      if (n > 5) {
	  if (val > 5) count++;
      }
      else if (n == val) count++;
  }
  return count;
}


static int checkVertex(meshMap *qm, int v)
{
  int q, i;

  if (v == -2) return EGADS_EMPTY;
  if (v <= 0  || v > qm -> totV) {
      printf(" Vertex for v %d is out of bounds !! \n ", v);
      return EGADS_INDEXERR;
  }
  if (qm -> vType[v-1] == -2) {
      printf(" V %d is empty !\n ", v);
      return EGADS_EMPTY;
  }
  if (qm -> vType[v-1] > 0 && qm -> valence[v-1][2] < qm -> vType[v-1]) {
      printf(" checkVertex:: Vertex %d Type %d is a doublet!!! We have big problems!!!\n",
	     v , qm ->vType[v - 1]);
      return EGADS_GEOMERR;
  }
  q = qm ->valence[v-1][0] - 1;
  if (q < 0 || q >= qm->totQ) {
      printf(" checkVertex:: Vertex for v %d has associated a quad %d out of bounds ( NMAX %d ) !! \n ",
	     v, q + 1, qm -> totQ-1);
      return EGADS_INDEXERR;
  }
  for (i = 0; i < 4; i++) {
      if (qm -> qIdx[4*q +i] == -2) {
	  printf("checkVertex::  Vertex for v %d has associated quad %d but  quad is empty !! \n ",
		 v, q + 1);
	  return EGADS_EMPTY;
      }
      if (qm -> qIdx[4*q+i] ==  v) return EGADS_SUCCESS;
  }
  printf("checkVertex::  Vertex for v %d has associated quad %d but !! \n ",
	 v, q + 1);
  printQuadSpecs(qm, q+1);
  return EGADS_INDEXERR;
}


static int checkQuad(meshMap *qm, int q)
{
  int i, j;

  if (q > qm -> totQ || q <= 0  || qm -> qIdx[4*(q-1)] == -2) {
      printf(" q %d -> max %d ass %d \n ",
	     q, qm -> totQ, qm -> qIdx[4*(q-1)]);
      return EGADS_EMPTY;
  }
  for (i = 0; i < 4; i++) {
      j = qm ->qIdx[4*(q-1)+i];
      if (j == -2) {
	  printf("j = %d \n ", j);
	  return EGADS_EMPTY;
      }
#ifdef DEBUG
      j = checkVertex(qm, qm ->qIdx[4*(q-1)+i]);
      if (j != EGADS_SUCCESS) {
	  printf(" VERTEX %d is %d \n ", qm ->qIdx[4 * (q - 1) + i] , j);
	  return j;
      }
#endif
  }
  return EGADS_SUCCESS;
}


static int checkMesh(meshMap *qm)
{
  int stat, i, j, k, val1, val2, v1, v2;

  for (i = 0; i < qm->totV; i++) {
      if (qm->vType[i] == -2)  continue;
      val1 = qm->valence[i][2];
      val2 = qm->valence[i][0];
      if (val1 <= 0 || val1 > qm -> totV) {
	  stat = setValence(qm, i + 1);
	  if (stat != EGADS_SUCCESS) {
	      printf(" checkMesh vertex %d setValence --> %d!!\n", i + 1, stat);
	      return stat;
	  }
      }
      if (val2 <= 0 || val2 > qm -> totQ) {
	  printf("In checkMesh quad for %d is out of bounds!!! %d > %d \n ",
		 i + 1, val2, qm -> totQ);
	  return EGADS_INDEXERR;
      }
      stat = checkQuad(qm, val2);
      if (stat != EGADS_SUCCESS) {
	  printf(" checkMesh vertex %d has associated quad %d --> %d!!\n",
		 i + 1, val2, stat);
	  return stat;
      }
      v1 = i + 1;
      for (j = 0; j < val1; ++j) {
	  v2 = qm -> valence[i][3 + j];
	  for (k = 0; k < qm ->valence[v2 - 1][2]; k++)
	    if (qm ->valence[v2 - 1][3 + k] == v1) {
		k = -1;
		break;
	    }
	  if (k != -1) {
	      printf(" checkMesh Vertex %d has assigned %d as link but %d doesn't point at %d\n",
		     v1, v2, v2, v1);
	      for (k = 0; k < qm ->valence[v1 - 1][2]; k++)
		printf(" V1 %d ==> %d \n ", k, qm -> valence[v1-1][3+k]);
	      printf("\n\n");
	      for (k = 0; k < qm ->valence[v2 - 1][2]; k++)
		printf(" V1 %d ==> %d \n ", k, qm -> valence[v2-1][3+k]);
	      return EGADS_INDEXERR;
	  }
      }
  }
  for (i = 0; i < qm -> totQ; i++) {
      if (checkQuad(qm, i+1) == EGADS_EMPTY) continue;
      for (j = 0; j < 4; ++j) {
	  v1 = qm -> qAdj[4 * i  + j] - 1;
	  if (v1 < 0) continue;
	  val1 = -1;
	  for (k = 0; k < 4; ++k) {
	      if (qm->qAdj[4*v1+k] == i + 1) {
		  val1 = 1;
		  break;
	      }
	  }
	  if (val1 == -1) {
	      fprintf (stderr, " checkMesh quads %d and %d don't point at each other\n",
		       i +1, v1 + 1);
	      printQuadSpecs(qm, i + 1);
	      printQuadSpecs(qm, v1 + 1);
	      return EGADS_INDEXERR;
	  }
      }
  }
  return EGADS_SUCCESS;
}



/* Assuming qID is collapsing through v */
static int validCollapse(meshMap *qm, int qID, int v)
{
  int j, k, kk, id, link, aux, aux2, type, val;

  if (v <= 0 || v > qm -> totV) return EGADS_INDEXERR;
  id     = EG_quadVertIdx(qm, qID, v);
  aux    = qm ->qIdx[4 * (qID - 1) + (id + 2)%4] - 1;
  if ( qm -> vType[aux] != -1 ) {
      for (j = 0; j < qm -> valence[v - 1][2]; j++) {
	  k  = qm -> valence[v - 1][3 + j] -1;
	  if (qm -> vType[k] == -1) continue;
	  if (EG_quadVertIdx(qm, qID, k + 1) < 0) {
#ifdef DEBUG
	      printf(" Vertex %d is bounds and not in quad %d\n ", k + 1, qID );
#endif
	      return 0;
	  }
      }
  }
  for (j   = 0; j < 2; j++) {
      link = qm ->qIdx[4 * (qID - 1) + (id + 2 * j + 1)%4];
      type = qm ->vType  [link - 1];
      val  = qm ->valence[link - 1][2];
      if      (type > 0 && val <= type) {
#ifdef DEBUG
	  printf(" link %d has valence %d type %d !\n ", link, val, type);
#endif
	  return 0;
      }
      else if (val == 3) {
	  for (k = 0; k < 3; k++) {
	      aux = qm ->valence[link - 1][3 + k] - 1;
	      if (EG_quadVertIdx(qm, qID,  aux + 1) >= 0) continue;
	      if (qm ->valence[aux][2] == qm -> vType[aux]) {
#ifdef DEBUG
	  printf(" Chain Effect Doublets %d --> %d %d \n",
		 aux + 1, qm ->valence[aux][2] , qm -> vType[aux]) ;
#endif
		  return 0;
	      }
	      if (qm ->valence[aux][2] != 3) continue;
	      for (kk  = 0; kk < 3; kk++) {
		  aux2 = qm ->valence[aux][3 + kk] - 1;
		  if (aux2 + 1 == link) continue;
		  if ( qm -> valence[aux2][2] == qm ->vType[aux] ||
		      (qm -> valence[aux2][2] == 3 && qm ->vType[aux] == -1)) {
#ifdef DEBUG
	  printf(" KK Chain Effect Doublets %d --> %d %d \n",
		 aux2 + 1, qm ->valence[aux][2] , qm -> vType[aux]) ;
#endif
		      return 0;
		  }
	      }
	  }
      }
  }
  return 1;
}


/* Assuming we will break link v1-v2 */
static int validSwap(meshMap *qm, int v1, int v2)
{
  int i, vs[2], type, val;

  vs[0] = v1 - 1; vs[1] = v2 - 1;
  for (i = 0; i < 2; i++) {
#ifdef DEBUG
      if (checkVertex(qm, vs[i] + 1) != EGADS_SUCCESS) return 0;
#endif
      type = qm ->vType  [vs[i]];
      val  = qm ->valence[vs[i]][2];
      if ((type > 0 && val <= type) || (type == -1 && val <= 3)) return 0;
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
static int EG_quadIsBoundary (meshMap *qm, int qID)
{
  int i, q, v, bV = 0 , bQ = 0;
  
#ifdef DEBUG
  i = checkQuad(qm, qID);
  if (i != EGADS_SUCCESS)  return i;
#endif
  for (i = 0; i < 4; i++) {
      v  = qm ->qIdx[4 * (qID - 1) + i];
      q  = qm ->qAdj[4 * (qID - 1) + i];
      if (qm ->vType[v - 1] >=  0) bV = 1;
      if (q                 == -1) bQ = 1;
  }
  if (bQ == 1) return 1;
  if (bV == 1) return 2;
  else return 0;
}


static int EG_createQuadGroup(meshMap *qm, quadGroup *qg, int q0, int q1)
{
  int i,  ids[2], piv = 0, aux, vaux[6], common[3], centre, c1 = 0, c2 = 0;

  qg -> q[0] = q0;
  qg -> q[1] = q1;
  for (i = 0; i < 6; i++) {
      qg -> verts[i] = -1;
      qg -> vals [i] = -1;
  }
  if (q0 < 0 || q0 > qm -> totQ || q1 < 0 || q1 > qm -> totQ) {
      printf(" Q0 %d Q1 %d beyond limits !!!! \n ", q0, q1);
      return EGADS_INDEXERR;
  }
  EG_commonVerts(qm, q0, q1, common);
  if (common[0] != 2) {
      printf(" You have tried to create a quad group where quads are not adjacent\n");
      printQuadSpecs(qm, q0);
      printQuadSpecs(qm, q1);
      return EGADS_INDEXERR;
  }
  centre = common[1];
  if (qm -> valence[common[2] -1][1] > qm -> valence[common[1] -1][1]) centre = common[2];
  ids[0] = EG_quadVertIdx(qm, qg -> q[0], centre);
  ids[1] = EG_quadVertIdx(qm, qg -> q[1], centre);
  piv    = 0;
  if (qm -> qAdj[4 * (qg ->q[0] - 1) + ids[0]] == qg ->q[1]) piv = 1;
  for (i = 0; i < 4; i++)
    qg -> verts[i] = qm -> qIdx[4 * (qg ->q[piv] - 1) + (ids[piv] + i)%4];
  aux = (piv + 1)%2;
  qg -> verts[4] = qm -> qIdx[4 * (qg ->q[aux] - 1) + (ids[aux] + 2)%4];
  qg -> verts[5] = qm -> qIdx[4 * (qg ->q[aux] - 1) + (ids[aux] + 3)%4];
  if (piv == 1) swapInt(&qg ->q[0], &qg ->q[1]);
  for (i = 0; i < 6; i++) qg -> vals[i] = qm -> valence[qg -> verts[i] -1][1];
  if (qg -> vals[0] == qg -> vals[3]) {
      if (qg -> vals[1] != 4) c1++;
      if (qg -> vals[5] != 4) c1++;
      if (qg -> vals[2] != 4) c2++;
      if (qg -> vals[4] != 4) c2++;
      if (c2 > c1) {
	  for (i = 0; i < 6; i++) vaux[i] = qg -> verts[i];
	  for (i = 0; i < 6; i++) qg -> verts[i] = vaux[(i + 3)%6];
	  for (i = 0; i < 6; i++) qg -> vals[i] = qm -> valence[qg -> verts[i] - 1][1];
	  swapInt(&qg ->q[0], &qg ->q[1]);
      }
  }

  return EGADS_SUCCESS;
}


static int EG_doubleCollapse(meshMap *qm, quadGroup qg, int forcing, int *activity)
{
  int i, stat, i3, i5, q3;

  *activity  = 0;
  if (qg.vals[0] * qg.vals[3] == 16) {
      for (i5 = 0; i5 < 2; i5++) {
	  for (i = i3 = 0; i3 < 6; i3++) {
	      if (qg.vals[i3] != 3) continue;
	      if (i5 == 1 || qg.vals[(i3+1)%6] >= 5 || qg.vals[(i3+5)%6] >= 5) {
		  i  =  1;
		  i5 = -1;
		  break;
	      }
	  }
	  if (i5 == -1) break;
      }
      if (i == 0) return EGADS_SUCCESS;

      i5 = (i3 + 1)%6;
      if (i5%3 == 0) i5 = (i3 + 5)%6;
      if (forcing == 0 && (qg.vals[i5] < 5 ||
	  ((qg.vals[(i3 + 3)%6] >= 5 && qg.vals[(i5 + 3)%6] != 3) ||
	      qg.vals[(i3 + 3)%6] * qg.vals[(i5 + 3)%6] == 16)))
	return EGADS_SUCCESS;
      if (forcing == 0 && (qg.vals[i5] < 5 && qg.vals[(i5 + 3)%6] < 5 &&
	  qg.vals[(i3 + 3)%6] == 4)) return EGADS_SUCCESS;
      else if (forcing == 1) {
	  if (     qg.vals[(i3 + 3)%6] == 3) i5 = (i3 + 3)%6;
	  else if (qg.vals[(i3 + 2)%6] >= 5) i5 = (i3 + 2)%6;
	  else if (qg.vals[(i3 + 4)%6] >= 5) i5 = (i3 + 4)%6;
	  else return EGADS_SUCCESS;
      }
  }
  else if (qg.vals[0] * qg.vals[3] == 12) {
      i3 = 0;
      if (qg.vals[3] == 3) i3 = 3;
      i5 = (i3 + 1)%6;
      if (qg.vals[i5] != 5) i5 = (i3 + 5)%6;
      if (qg.vals[i5] != 5 || qg.vals[(i5 + 3)%6] == 4) return EGADS_SUCCESS;
  }
  else return EGADS_SUCCESS;

  q3 = 0;
  if (forcing == 1 && qg.vals[0] * qg.vals[3] == 16) {
      if (i3 >= 3) q3 = 1;
  } else {
      if (i5 >= 3) q3 = 1;
  }
  for (i = 0; i < 2; i++) {
      if (qm -> qIdx[4 * (qg.q[(q3 + i)%2] -1)] == -2) continue;
      stat        = EG_collapse(qm, qg.q[(q3 + i)%2], &i3, 1);
      if (stat   == EGADS_SUCCESS) {
	  if (i3 == 0) break;
	  (*activity)++;
      }
      if (stat != EGADS_SUCCESS) {
	  printf("EG_doubleCollapse after EG_forceColapse quad %d -> %d !!\n ",
		 qg.q[(q3 + i)%2], stat);
	  return stat;
      }
  }
  return EGADS_SUCCESS;
}


static int EG_swapDoubleCollapse(meshMap *qm, quadGroup qg, int *activity)
{
  int  k, swap = 1, id, j, stat, v0;

  *activity  = 0;
  if (qg.vals[0] * qg.vals[3] != 20 ||
      qg.vals[2] * qg.vals[4] != 9  ||
      validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
  if (qg.vals[swap] * qg.vals[(swap+3)%6] != 12) swap = 2;
  if (qg.vals[swap] * qg.vals[(swap+3)%6] != 12) return EGADS_SUCCESS;
  v0        = qg.verts[0];
  stat      = EG_swappingOperation(qm, qg, swap, &(*activity));
  if (stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
      printf("In EG_swapDoubleCollapse: Swapping went %d activity %d \n",
	      stat, *activity);
#endif
      return stat;
  }
  *activity = 1;
  for (k = 0; k < 2; k++) {
      id = EG_quadVertIdx(qm, qg.q[k], qg.verts[3]);
      if (id >= 0) break;
  }
  if (id < 0) {
      printf(" I can't find vertex %d in quads!!! \n", qg.verts[3]);
      printQuadGroup(qm, qg);
      return EGADS_INDEXERR;
  }
  stat      = EG_collapse(qm, qg.q[k], &(*activity), 1);
  if (stat == EGADS_SUCCESS && *activity > 0) return stat;
  stat      = EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
  if (stat != EGADS_SUCCESS) {
      printf("EG_swapDoubleCollapse after invalid force collapse :: EG_createQuadGroup is %d \n ",
	     stat);
      printQuadGroup(qm, qg);
      return stat;
  }
  for (swap = 0; swap < 6; swap++) if (qg.verts[swap] == v0) break;
  stat = EG_swappingOperation(qm, qg, swap, &j);
  if (stat == EGADS_SUCCESS && j == 1) *activity = 0;
  else {
#ifdef DEBUG
      printf(" In EG swapDoubleCollapse: I failed to collapse after swapping! s = %d act = %d \n ",
	     stat, j);
#endif
  }
  return stat;
}


static int EG_swapDoubleSplit(meshMap *qm, quadGroup qg, int *activity)
{
  int       i5, q, i55, i3, val3, i0, v30[2], i,  stat, adj[2], q0[2], v0;
  quadGroup sqg;

  *activity = 0;
  for (q = i3 = 0; i3 < 6; i3++)
    if (qg.vals[i3] == 3 && i3%3 != 0) {
	q = 1;
	break;
    }
  if (q == 0) return EGADS_SUCCESS;
  i5  = (i3 + 3)%6;
  if (qg.vals[i5] < 5) return EGADS_SUCCESS;
  i55 = (i5 + 1)%6;
  if (qg.vals[i55] < 5) {
      i55 = (i5 + 5)%6;
      if (qg.vals[i55] < 5) return EGADS_SUCCESS;
  }
  q = 0;
  if (EG_quadVertIdx(qm, qg.q[q], qg.verts[i5]) < 0) q = 1;
  stat = EG_adjQtoPair(qm, qg.q[q], qg.verts[i5], qg.verts[i55], adj);
  if (stat != EGADS_SUCCESS || adj[1] == -1) {
      if (stat != EGADS_SUCCESS)
	printf("In EG_swapDoubleSplit: adjacent to pair %d --> !!\n ", stat);
      return stat;
  }
  if ((i3 + 1)%6 == 0 || (i3 + 5)%6 == 0) i0 = 0;
  else                                    i0 = 3;
  v30[0] = qg.verts[i0];
  v30[1] = qg.verts[(i0 + 3)%6];
  q0[0]  = qg.q[0];
  q0[1]  = qg.q[1];
  val3   = qg.verts[i3];
  stat   = EG_createQuadGroup(qm, &qg, qg.q[q], adj[1]);
  if (stat != EGADS_SUCCESS) {
      printf("Inside EG_swapDoubleSplit: before swapping EG_createQuadGroup stat %d\n ",
	     stat);
      printQuadGroup(qm, qg);
      return stat;
  }
  if (validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
  for (i0 = 0; i0 < 6; i0++) if (qg.verts[i0] == v30[0]) break;
  v0            = qg.verts[0];
  stat          = EG_swappingOperation(qm, qg, i0, &(*activity));
  if (stat     != EGADS_SUCCESS || *activity == 0) {
      if (stat != EGADS_SUCCESS)
	  printf(" EG_swapDoubleSplit error at swap: %d activity %d!!\n ",
	  	     stat, *activity);
#ifdef DEBUG
      printf(" EG_swapDoubleSplit error at swap: %d activity %d!!\n ",
	     stat, *activity);
#endif
      return stat;
  }
  i = 0;
  if (EG_quadVertIdx(qm, q0[0], val3) < 0) i = 1;
  stat = EG_adjQtoPair(qm, q0[i], v30[0], v30[1], adj);
  if (stat != EGADS_SUCCESS || adj[1] == -1) {
      if (stat != EGADS_SUCCESS)
	printf(" EG_swapDoubleSplit after swapping adjacent to pair %d !!\n ",
	       stat);
      return stat;
  }
  stat = EG_createQuadGroup(qm, &sqg, q0[i], adj[1]);
  if (stat != EGADS_SUCCESS) {
      printf("Inside EG_swapDoubleSplit: before splitting EG_createQuadGroup stat %d\n ",
	     stat);
      printQuadGroup(qm, qg);
      return stat;
  }
  stat   = EG_doubleSplit(qm, sqg, 0, &i);
  if (i == 0 && stat == EGADS_SUCCESS) {
      stat      = EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
      if (stat != EGADS_SUCCESS) {
	  printf(" In EG_swapDoubleSplit: recover swap create group %d\n ", stat);
	  printQuadGroup(qm, qg);
	  return EGADS_SUCCESS;
      }
      for (i = 0; i < 6; i++) if (qg.verts[i] == v0) break;
      stat = EG_swappingOperation(qm, qg, i, &q);
      if (stat == EGADS_SUCCESS && q == 1) *activity = 0;
  }
  return stat;
}


static int EG_transferValences(meshMap *qm, int *qID, int try5533,
			       int *transfering, int *activity)
{
  int       i, j, swap = 0, stat, min, links[5], qAux[2];
  quadGroup qg;

  *activity = 0;
  if (qID[0] < 0 || qm -> qIdx[4 * (qID[0] - 1)] == -2) {
#ifdef DEBUG
      printf(" EG_transferValences for quad %d is bad quad %d \n ",
	     qID[0], checkQuad(qm, qID[0]));
#endif
      return EGADS_SUCCESS;
  }

  if (*transfering == 0) {
      stat      = EG_cleanQuad(qm, qID[0], 1, 0, try5533, &(*activity));
      if (stat != EGADS_SUCCESS || *activity > 0) {
	  qID[0] = -1;
	  if (stat != EGADS_SUCCESS)
	      printf("EG_transferValences: EG_cleanQuad %d --> %d!!\n",
		     qID[0], stat);
	  return stat;
      }
  }
  if (qID[0] == -1 || qm -> qIdx[4 * (qID[0]- 1)] == -2 ||
      EG_nValenceCount(qm, qID[0], 4) == 4) return EGADS_SUCCESS;
  for (j = 0; j < 4; j++) {
      i = qm -> qAdj[4 * (qID[0] - 1) + j];
      if (i == -1 || i == qID[1]) continue;
      if (*transfering == 0) {
	  stat      = EG_cleanQuad(qm, i, 1, 0, try5533, &(*activity));
	  if (stat != EGADS_SUCCESS || *activity > 0) {
	      if (stat != EGADS_SUCCESS)
		  printf(" EG_TransferValence EG_cleanQuad --> %d !!\n ", stat);
	      qID[0] = -1;
	      return stat;
	  }
      }
      if (qm -> qAdj[4 * (qID[0] - 1) + j] != i) continue;
      stat = EG_createQuadGroup(qm, &qg, qID[0], i);
      if (stat != EGADS_SUCCESS) {
	  printf(" Inside EG_transferValences EG_createQuadGroup %d !!\n", stat);
	  printQuadGroup(qm, qg);
	  return stat;
      }
      if (*transfering == 0 && qg.vals[0] * qg.vals[3] == 15 && qID[1] != -1) {
	  j = 0;
	  if (qg.q[j] == qID[0]) j = 1;
	  qID[0] = qg.q[j];
	  qID[1] = qg.q[(j + 1)%2];
	  if (EG_quadIsBoundary(qm, qg.q[1]) == 1) {
	      stat      = EG_collapse(qm, qID[0], &(*activity), 1);
	      if (stat != EGADS_SUCCESS)
		printf("EG_transferValences: forceCollapse gave %d !!\n", stat);
	      qID[0]    = qID[1];
	      qID[1]    = -1;
	      return stat;
	  }
	  return EG_transferValences(qm, qID, 0, &(*transfering), &(*activity));
      }
      if (validSwap(qm, qg.verts[0], qg.verts[3]) == 0) continue;
      if          (qg.vals[0] * qg.vals[3] >= 20) min = 12;
      else if (try5533 == 1 && *transfering == 0) {
	  if      (qg.vals[0] * qg.vals[3] >= 25) min = 16;
	  else if (qg.vals[0] * qg.vals[3] == 16) min =  9;
      } else continue;
      if (     qg.vals[1] * qg.vals[4] <= min) swap = 1;
      else if (qg.vals[2] * qg.vals[5] <= min) swap = 2;
      else continue;
      stat          = EG_swappingOperation(qm, qg, swap, &(*activity));
      if (stat     != EGADS_SUCCESS || *activity == 0) {
	  if (stat != EGADS_SUCCESS)
	    printf(" EG_transferValences EG_swappingOperation stat %d activity %d !!\n",
		   stat, *activity);
#ifdef DEBUG
	  printf(" EG_transferValences EG_swappingOperation stat %d activity %d !!\n",
		 stat, *activity);
#endif
	  return stat;
      }
      i        = 0;
      if (min == 9 || min == 16) *transfering = 1;
      else                       *transfering = 0;
      if (qID[1] > 0 && qm -> qIdx[4 * (qID[1] - 1)] != -1) {
	  if      (qg.q[0] == qID[1]) i = 1;
	  else if (qg.q[1] == qID[1]) i = 0;
	  else {
	      EG_commonVerts(qm, qg.q[0], qID[1], links);
	      if (links[0] > 0) {
		  if (links[0] == 1) {
		      EG_commonVerts(qm, qg.q[1], qID[1], links);
		      if (links[0] == 0) i = 1;
		  }
		  else i = 1;
	      }
	  }
      }
      else if (EG_nValenceCount(qm,qg.q[0], 5) == 2 &&
	       EG_nValenceCount(qm,qg.q[1], 5) >= 1 &&
	       EG_nValenceCount(qm,qg.q[1], 3) == 1) i = 1;
      else if (EG_nValenceCount(qm,qg.q[1], 4) <
	       EG_nValenceCount(qm,qg.q[0], 4)) i = 1;
      qID[0] = qg.q[ i       ];
      qID[1] = qg.q[(i + 1)%2];
      if (*transfering == 1) {
	  for (j = 0; j < 2; j++) {
	      i         = 1;
	      qAux[0]   = qID[ j       ];
	      qAux[1]   = qID[(j + 1)%2];
	      stat      = EG_transferValences(qm, qAux, 0, &i, &min);
	      if (stat != EGADS_SUCCESS) {
		  printf(" EG_transferValences: separating valences after forcing %d!!\n",
			 stat);
		  return stat;
	      }
	  }
	  *transfering = 0;
	  qID[0]       = qAux[0];
	  qID[1]       = qAux[1];
      }
      if (*activity > 0) break;
  }
  return stat;
}


static int EG_basicOperation (meshMap *qm, int qID, int type, int *activity)
{
  int stat = EGADS_SUCCESS;

  switch (type) {
    case SWAP:
      stat = EG_swap(qm, qID, &(*activity));
      break;
    case COLLAPSE:
      stat = EG_collapse(qm, qID, &(*activity), 0);
      break;
    case SPLIT:
      stat = EG_split(qm, qID, &(*activity));
      break;
  }
  if (stat != EGADS_SUCCESS && *activity == 0) return EGADS_SUCCESS;

  return stat;
}


static int EG_composeOperation(meshMap *qm, quadGroup qg, int type, int forcing,
			       int *activity)
{
  int stat = EGADS_SUCCESS;

  switch (type) {
    case DOUBLESWAP:
      stat = EG_doubleSwap(qm, qg, forcing, &(*activity));
      break;
    case SWAPCOLLAPSE:
      stat = EG_swapCollapse(qm, qg,forcing, &(*activity));
      break;
    case DOUBLECOLLAPSE:
      stat = EG_doubleCollapse(qm, qg, forcing,  &(*activity));
      break;
    case SWAPDOUBLECOLLAPSE:
      stat = EG_swapDoubleCollapse(qm, qg, &(*activity));
      break;
    case SWAPSPLIT:
      stat = EG_swapSplit(qm, qg, forcing, &(*activity));
      break;
    case DOUBLESPLIT:
      stat = EG_doubleSplit(qm, qg, forcing, &(*activity));
      break;
    case SWAPDOUBLESPLIT:
      stat = EG_swapDoubleSplit(qm, qg, &(*activity));
      break;
  }
  if (stat != EGADS_SUCCESS)
    printf(" OPERATION %d gave %d and activity %d \n ", type, stat, *activity);

  return stat;
}


static int EG_cleanQuad(meshMap *qm, int qID, int useAdj, int transfer,
			int forcing, int *activity)
{
  int stat, i, q, qadj, act = 0;
  int opBasic[3] = {COLLAPSE, SWAP, SPLIT};
  int opComp [7] = {SWAPCOLLAPSE, DOUBLECOLLAPSE, SWAPDOUBLECOLLAPSE, DOUBLESWAP,
                    SWAPSPLIT,    DOUBLESPLIT,    SWAPDOUBLESPLIT};
  quadGroup qg;

  *activity = 0;
  if (qm -> qIdx[4 * (qID - 1)] == -2) {
#ifdef DEBUG
      printf(" EG_cleanQuad %d stat %d\n", qID, stat);
#endif
      return EGADS_SUCCESS;
  }
  if (EG_nValenceCount(qm, qID, 4) == 4) return EGADS_SUCCESS;
#ifdef DEBUG
  printf(" EG CLEAN QUAD %d  transfer %d forcing %d\n ", qID, transfer, forcing);
  printQuadSpecs (qm, qID ) ;
#endif
  if (transfer == 0) {
      for (i = 0; i < 3; i++) {
	  stat       = EG_basicOperation(qm, qID, opBasic[i], &act);
	  *activity += act;
	  if (stat != EGADS_SUCCESS) {
	      printf(" EG_cleanQuad Basic operation %d !!\n ", stat);
	      return stat;
	  }
	  if (act > 0 || qID < 0 || qm -> qIdx[4 * (qID - 1)] == -2)
            return EGADS_SUCCESS;
      }
  }
  if (useAdj == 0)  return EGADS_SUCCESS;
  if (qID < 0 || qm -> qIdx[4 * (qID - 1)] == -2) {
#ifdef DEBUG
      printf(" EG_cleanQuad %d checkQuad is stat %d!!\n", qID, stat);
#endif
      return EGADS_SUCCESS;
  }
  if (forcing == 1 && EG_nValenceCount(qm, qID, 3) > 1) {
#ifdef DEBUG
      printf(" FORCE COLLAPSE AROUND QUAD %d \n ", qID);
#endif
      stat = EG_collapse(qm, qID, &(*activity), 1);
      if (stat != EGADS_SUCCESS || *activity > 0) {
	  if (stat != EGADS_SUCCESS)
              printf("EG_cleanQuad force Collapse went %d!1\n ", stat);
	  return stat;
      }
  }
  if (transfer == 0) {
      for (q = 0; q < 4; q++) {
	  qadj      = qm -> qAdj[4 * (qID - 1) + q];
	  if (qadj == -1) continue;
	  stat       = EG_cleanQuad(qm, qadj, 0, 0, 0, &act);
	  *activity += act;
	  if (stat != EGADS_SUCCESS) {
	      printf(" EG_cleanQuad restore mesh after adjacent stat --> %d !!\n",
		     stat);
	      return stat;
	  }
	  if (act > 0 || qID < 0 || qm -> qIdx[4 * (qID - 1)] == -2)
	    return EGADS_SUCCESS;
      }
  }
  for (i = 0; i < 7; i++) {
      if (transfer == 1 && opComp[i] == DOUBLESWAP) continue;
      for (q  = 0; q < 4; q++) {
	  if (qID < 0 || qm -> qIdx[4 * (qID - 1)] == -2)
	    return EGADS_SUCCESS;
	  qadj = qm -> qAdj[4 * (qID - 1) + q];
	  if (qadj == -1) continue;
	  stat = EG_createQuadGroup(qm, &qg, qID, qadj);
	  if (stat != EGADS_SUCCESS) {
	      printf("Inside EG_cleanQuad: EG_createQuadGroup stat --> %d !!\n ",
		     stat);
	      printQuadGroup(qm, qg);
	      return stat;
	  }
	  stat       = EG_composeOperation(qm, qg, opComp[i], forcing, &act);
	  *activity += act;
	  if (stat != EGADS_SUCCESS) {
	      printf("In EG_cleanQuads restor mesh after EG_composeOperation stat -->  %d !!\n ",
		     stat);
	      return stat;
	  }
	  if (act > 0 || qID < 0 || qm -> qIdx[4 * (qID - 1)] == -2)
	    return EGADS_SUCCESS;
      }
  }

  return EGADS_SUCCESS;
}


void EG_destroymeshMap(bodyQuad *bodydata)
{
  int i, j;

  if (bodydata->qm == NULL) return;
  for (i = 0; i < bodydata->nfaces; ++i) {
      if (bodydata->qm[i]) {
	  for (j = 0; j < bodydata -> qm[i] -> sizeV; ++j)
            EG_free(bodydata -> qm[i] -> valence[j]);
	  EG_free(bodydata -> qm[i] -> valence);
	  EG_free(bodydata -> qm[i] -> qIdx);
	  EG_free(bodydata -> qm[i] -> qAdj);
	  EG_free(bodydata -> qm[i] -> uvs);
	  EG_free(bodydata -> qm[i] -> xyzs);
	  EG_free(bodydata -> qm[i] -> remQ);
	  EG_free(bodydata -> qm[i] -> remV);
	  EG_free(bodydata -> qm[i] -> vType);
	  EG_free(bodydata -> qm[i] -> vFix);
	  EG_free(bodydata -> qm[i]);
      }
  }
  EG_free(bodydata->qm);
  bodydata->qm = NULL;
}


static int EG_adjQtoPair(meshMap *qm, int qID, int v1, int v2, int *adj)
{
  int i, aux = -1;

  adj[0] = -1; adj[1] = -1;
  for (i = 0; i < 4; ++i) {
      if (qm -> qIdx[4*(qID - 1) + i] == v1) adj[0] = i;
      if (qm -> qIdx[4*(qID - 1) + i] == v2) aux    = i;
      if (aux != -1 && adj[0] != -1) break;
  }
  if (aux == -1 || adj[0] == -1) return EGADS_SUCCESS;
  if      (abs(adj[0] - aux) == 3) adj[0] = 3;
  else if (aux < adj[0]          ) adj[0] = aux;
  adj[1] = qm -> qAdj[4*(qID - 1) + adj[0]];

  return EGADS_SUCCESS;
}


static int resizeQm(meshMap *qm)
{
  int    stat = EGADS_SUCCESS, nV, vRem, nQ, qRem, i, j, k;
  int    *vpiv = NULL, *qpiv = NULL;
  int    *verts = NULL, *adjs = NULL, *type = NULL, bV, bQ;
  double *uvxyz = NULL;

  vpiv = (int *) EG_alloc(qm -> totV * sizeof(int));
  qpiv = (int *) EG_alloc(qm -> totQ * sizeof(int));
  if (vpiv == NULL || qpiv == NULL) return EGADS_MALLOC;
  for (nV = vRem = i = 0; i < qm ->totV; i++) {
      if (qm -> vType[i] != -2) vpiv[i] = nV++;
      else {
	  vpiv[i] = -2;
	  vRem++;
      }
  }
  for (nQ = qRem = i = 0; i < qm ->totQ; i++) {
      if (qm ->qIdx[4*i] != -2) qpiv[i] = nQ++;
      else {
	  qpiv[i] = -2;
	  qRem++;
      }
  }
  if (vRem != qRem) {
      printf(" In resizeQm: I have %d removed vertices and %d quads!! they should match!!!!!\n ",
	     vRem, qRem);
      stat = EGADS_INDEXERR;
      goto cleanup;
  }
  if (vRem == 0) {
      EG_free(vpiv);
      EG_free(qpiv);
      return EGADS_SUCCESS;
  }
  bV    = qm -> totV;
  bQ    = qm -> totQ;
  uvxyz = (double *) EG_alloc(5 * bV *sizeof(double));
  type  = (int    *) EG_alloc(    bV *sizeof(int));
  verts = (int    *) EG_alloc(4 * bQ *sizeof(int));
  adjs  = (int    *) EG_alloc(4 * bQ *sizeof(int));
  if (uvxyz == NULL || type == NULL ||
      verts == NULL || adjs == NULL ) {
      stat = EGADS_MALLOC;
      goto cleanup;
  }
 for (j = 0; j < bQ; ++j) {
      for (k = 0; k < 4; ++k) {
	  verts[4*j + k] = qm -> qIdx[4*j + k];
	  adjs [4*j + k] = qm -> qAdj[4*j + k];
      }
  }
  for (j = 0; j < bV; ++j) {
      type [  j    ] = qm -> vType[  j    ];
      uvxyz[5*j    ] = qm -> uvs  [2*j    ];
      uvxyz[5*j + 1] = qm -> uvs  [2*j + 1];
      uvxyz[5*j + 2] = qm -> xyzs [3*j    ];
      uvxyz[5*j + 3] = qm -> xyzs [3*j + 1];
      uvxyz[5*j + 4] = qm -> xyzs [3*j + 2];
  }
  qm -> totV    = nV;
  qm -> totQ    = nQ;
  qm -> remQ[0] = 0;
  qm -> remV[0] = 0;
  for (i = 0; i < bQ; i++) {
      if (qpiv[i] == -2) continue;
      for (k = 0; k < 4; ++k) {
	  qm -> qIdx[4*qpiv[i] + k ] = vpiv[verts[4*i + k] - 1] + 1;
	  if (adjs[4*i + k] == -1) qm -> qAdj[4*qpiv[i] + k] = - 1;
	  else                     qm -> qAdj[4*qpiv[i] + k] = qpiv[adjs[4*i + k] - 1] + 1;
	  qm -> valence[qm -> qIdx[4*qpiv[i] + k ] - 1][0] = qpiv[i] + 1;
      }
  }
  for (i = 0; i < bV; i++) {
      if (vpiv[i] == -2) continue;
      j = vpiv[i];
      qm -> vType  [j      ] = type [  i    ];
      qm -> uvs    [2*j    ] = uvxyz[5*i    ];
      qm -> uvs    [2*j + 1] = uvxyz[5*i + 1];
      qm -> xyzs   [3*j    ] = uvxyz[5*i + 2];
      qm -> xyzs   [3*j + 1] = uvxyz[5*i + 3];
      qm -> xyzs   [3*j + 2] = uvxyz[5*i + 4];
      setValence (qm, j + 1);
  }
cleanup:
  EG_free(vpiv);
  EG_free(qpiv);
  EG_free(uvxyz);
  EG_free(type);
  EG_free(verts);
  EG_free(adjs);
  return stat;
}


static void updateVertex(meshMap *qm, int vID, double *uv)
{
  int    i;
  double eval[18];

  i = EG_evaluate(qm -> face, uv, eval);
  if (i != EGADS_SUCCESS || qm -> vType [vID -1] != -1) return;
  qm -> uvs [2 * (vID - 1)    ] = uv  [0];
  qm -> uvs [2 * (vID - 1) + 1] = uv  [1];
  qm -> xyzs[3 * (vID - 1)    ] = eval[0];
  qm -> xyzs[3 * (vID - 1) + 1] = eval[1];
  qm -> xyzs[3 * (vID - 1) + 2] = eval[2];
}


static int EG_buildStar(meshMap *qm, vStar **star, int vID)
{
  int i = 0, id0 = -1, q = 0, auxV, auxQ, v = 0, quadID, prevQuad;
  int it = 0, it2 = 0, adj[2], *vertex = NULL, *quads = NULL;
  int qLoop[8] = {0, 1, 2, 3, 0, 1, 2, 3};

#ifdef DEBUGGG
  i = checkVertex(qm , vID);
  if (i != EGADS_SUCCESS) {
      printf(" Vertex %d is %d \n ", vID, i);
      return i;
  }
#endif
  vertex = (int *) EG_alloc(qm -> totV * sizeof(int));
  quads  = (int *) EG_alloc(qm -> totQ * sizeof(int));
  if (vertex == NULL || quads == NULL) {
      printf("EG_buildStar MALLOC at quads & verts!!\n ");
      if (vertex != NULL) EG_free(vertex);
      if (quads  != NULL) EG_free(quads);
      return EGADS_MALLOC;
  }
  // quads are -1 bias
  quadID = qm -> valence[vID - 1][0] - 1;
#ifdef DEBUGGG
  i = checkQuad(qm, quadID + 1);
  if (i != EGADS_SUCCESS) {
      printf(" EG_buildStar at vertex %d has associated a bad quad: %d --> %d!!\n",
	     vID, quadID + 1, i);
      printQuadSpecs(qm, quadID + 1);
      EG_free(vertex);
      EG_free(quads);
      return i;
  }
#endif
  vertex[v++] = vID;
  it          = 0;
  do {
#ifdef DEBUGGG
      i = checkQuad(qm, quadID + 1);
      if (i != EGADS_SUCCESS) {
	  printf(" In EG_buildStar quad %d is bad quad --> %d!!\n",
		 quadID + 1, i);
	  printQuadSpecs(qm, quadID + 1);
	  EG_free(quads);
	  EG_free(vertex);
	  return i;
      }
#endif
      id0 = EG_quadVertIdx(qm, quadID + 1, vID);
      if (id0 < 0) {
	  printf(" In EG_buildStar id for Vert %d in Quad %d is %d !!\n",
		 vID, quadID + 1, id0);
#ifdef DEBUG
	  printQuadSpecs(qm, quadID + 1);
#endif
	  EG_free(vertex);
	  EG_free(quads);
	  return EGADS_INDEXERR;
      }
      for (i = 1; i <= 2; ++i)
	vertex[v++] = qm -> qIdx[4 * quadID + qLoop[id0 + i]];
      quads[q++] = quadID + 1;
      prevQuad   = quadID;
      quadID     = qm -> qAdj[4 * prevQuad + qLoop[id0 + 3 ]] - 1;
      if (quadID  < 0) { //make a "ghost" quad
	  auxQ        = prevQuad;
	  vertex[v++] = qm -> qIdx[4*auxQ + qLoop[id0 + 3]];
	  auxV        = qm -> qIdx[4*auxQ + qLoop[id0 + 1]];
	  it2 = 0;
	  do {
	      i = EG_adjQtoPair(qm, auxQ + 1, vID, auxV, adj);
	      if (adj[1] == -1 || i != EGADS_SUCCESS) break;
	      auxQ = adj[1] - 1;
	      i    = EG_quadVertIdx(qm, auxQ + 1, vID);
	      if (i < 0) {
		  printf(" In buildStar vertex Id %d in quad %d is %d\n",
			 vID, quadID + 1, i);
#ifdef DEBUG
		  printQuadSpecs(qm, quadID + 1);
#endif
		  EG_free(vertex);
		  EG_free( quads);
		  return EGADS_INDEXERR;
	      }
	      auxV = qm -> qIdx[4 * auxQ + qLoop[i + 1]];
	      it2++;
	      if (it2 > 200) {
		  printf(" stuck in interior loop of build star!!!!!!!!\n");
		  EG_free(quads);
		  EG_free(vertex);
		  return EGADS_RANGERR;
	      }
	  } while (adj[1] != - 1);
	  quads[q++]  = -1;
	  vertex[v++] = -1;
	  quadID      = auxQ;
      }
      if (quadID < 0) {
	  printf(" I am stuck in build star. Pointing a NULL quad \n");
	  EG_free(vertex);
	  EG_free(quads);
	  return EGADS_INDEXERR;
      }
      it++;
      if (it > 200) {
	  printf(" EG_buildStar:: stuck in outer loop of build star!!!!!!!!\n");
	  EG_free(vertex);
	  EG_free(quads);
	  return EGADS_RANGERR;
      }
  } while (quadID + 1 != quads[0]);
  if (*star != NULL) EG_freeStar(&(*star));
  *star = (vStar *) EG_alloc(sizeof(vStar));
  if ((*star) == NULL) {
      EG_free(vertex);
      EG_free(quads);
      return EGADS_MALLOC;
  }
  (*star) -> nQ    = q;
  (*star) -> nV    = v;
  (*star) -> verts = (int *) EG_alloc (    v * sizeof(int));
  (*star) -> quads = (int *) EG_alloc (    q * sizeof(int));
  (*star) -> idxV  = (int *) EG_alloc (2 * v * sizeof(int));
  (*star) -> idxQ  = (int *) EG_alloc (2 * q * sizeof(int));
  if ((*star) -> verts == NULL || (*star) -> quads == NULL ||
      (*star) -> idxV  == NULL || (*star) -> idxQ  == NULL) {
      if ((*star) -> verts != NULL) EG_free((*star) -> verts);
      if ((*star) -> quads != NULL) EG_free((*star) -> quads);
      if ((*star) -> idxV  != NULL) EG_free((*star) -> idxV);
      if ((*star) -> idxQ  != NULL) EG_free((*star) -> idxQ);
      EG_free((*star));
      EG_free(vertex);
      EG_free( quads);
      return EGADS_MALLOC;
  }
  for (i = 0; i < q; ++i) {
      (*star) -> quads[i    ] = quads[i];
      (*star) -> idxQ [i    ] = i;
      (*star) -> idxQ [q + i] = i;
  }
  for (i = 0; i < v; ++i) {
      (*star) -> verts[i    ] = vertex[i];
      (*star) -> idxV [i    ] = i;
      (*star) -> idxV [v + i] = i + 1;
  }
  EG_free (vertex);
  EG_free (quads);
  return EGADS_SUCCESS;
}


static int EG_collapse(meshMap *qm, int qID, int *activity, int forcing)
{
  int i, vC, v5[5], v3[5], qb[4], stat, i3 = 0, val[2], links[2], doublet;

  *activity = 0;
  if (qm -> qIdx[4 * (qID - 1)] == -2) {
#ifdef DEBUG
      printf(" EG_collapse for quad %d is %d \n ",
	     qID, checkQuad(qm, qID));
#endif
      return EGADS_SUCCESS;
  }
  doublet    = -1;
  for (v3[0] = v5[0] = i = 0; i < 4; i++) {
      qb[i]  = qm -> qIdx[4 * (qID - 1) + i];
      val[0] = qm -> valence[qb[i] -1][1];
      if      ( val[0] == 3) v3[++v3[0]] = i;
      else if ( val[0] >= 5) v5[++v5[0]] = i;
      else if ( val[0] == 2) doublet     = i;
  }
  printMesh(qm, NULL);

  if (doublet != -1) {
#ifdef DEBUG
      printMesh(qm, NULL);

#endif
      stat = EG_mergeVertices(qm, qID, qb[doublet], &(*activity));
      if (stat != EGADS_SUCCESS || *activity == 0) {
	  printf(" I FAILED WITH DOUBLETS !!!!!!!! \n ");
	  return EGADS_GEOMERR;
      }
  }
  if (forcing == 0 ) {
      if (v3[0] + v5[0] <= 2 || v3[0] * v5[0] == 0 ) return EGADS_SUCCESS;
      if (v5[0] == 1 || (v5[0] == 2 && abs(v5[1] - v5[2])%2 == 0)) {
	  vC   = qm -> qIdx[4 * (qID - 1) + (v5[1] + 1)%4];
	  if (validCollapse(qm, qID, vC) == 1) {
	      stat         = EG_mergeVertices(qm, qID, vC, &(*activity));
	      if (stat    != EGADS_SUCCESS) {
		  printf(" EG_collapse :: EG_mergeVertices %d activity %d \n",
			 stat, *activity);
		  return stat;
	      }
	  }
      }
  }
  if (*activity == 0 && forcing == 1 ) {
      if ( EG_nValenceCount(qm, qID, 3) == 0) return EGADS_SUCCESS;
      i = 0; i3 = 0;
      while (i < 2) {
	  vC = qm -> qIdx[4 * (qID - 1) + i3];
	  if (qm->vType[vC - 1] == -1 && qm->valence[vC - 1][2] == 3) {
	      links[0] = qm -> qIdx[4 * (qID - 1) + (i3 + 1)%4];
	      links[1] = qm -> qIdx[4 * (qID - 1) + (i3 + 3)%4];
	      val  [0] = qm -> valence[ links[0] - 1][1];
	      val  [1] = qm -> valence[ links[1] - 1][1];
	      if ((qm -> valence[ links[0] - 1][1] >= 5 ||
		   qm -> valence[ links[1] - 1][1] >= 5 || i == 1) &&
		  validCollapse(qm, qID, vC) == 1) {
		  stat         = EG_mergeVertices(qm, qID, vC, &(*activity));
		  if (stat    != EGADS_SUCCESS) {
		      printf(" EG_collapse :: forcing EG_mergeVertices at quad %d thru %d -- %d activity %d \n",
			     qID, vC, stat, *activity);
		      printQuadSpecs(qm, qID ) ;
		      vStar * star = NULL;
		      char buff[100];
		      for ( i = 0; i < 4; i++ ) {
			  stat = EG_buildStar(qm, &star,  qm -> qIdx[4 * (qID - 1) + i]);

			  snprintf(buff, 100, "star_%d_%d",star -> verts[0], qm -> plotcount++ );
			  printStarFile( qm, star, buff);
		      }
		      exit(1);
		  }
		  if (*activity > 0) break;
	      }
	  }
	  i3++;
	  if (i3 == 4) {
	      i++;
	      i3 = 0;
	  }
      }
  }
  if (stat == EGADS_SUCCESS && *activity == 0) return stat;
  for (i = 0; i < 4; i++) {
      vC = qm -> vType[qb[i] - 1];
      if ( vC == -2 ) continue;
      printf(" i = %d LEAVING QUAD %d VALENCE %d  \n ",i,  qb[i], qm -> valence[qb[i] - 1][2] );
      if (vC == -1 && qm -> valence[qb[i] - 1][2] == 2) {
	  stat = EG_collapse(qm, qm -> valence[qb[i] - 1][0],  &(*activity), 0);
	  if (stat != EGADS_SUCCESS || *activity == 0) {
	      printf(" I FAILED WITH DOUBLETS !!!!!!!! \n ");
	      return EGADS_GEOMERR;
	  }
	  i = 0 ;
      }
  }
  return EGADS_SUCCESS;
}


static int EG_mergeVertices(meshMap *qm, int qC, int centre, int *activity)
{
  int    stat, i, j, q,  adjq, adjPair[2], auxQ, oldQ[8], nq, doublet = 0;
  int    piv[4] = {1, 0, 3, 2};
  double uv[2], uvxyz[10];
  Quad   *quad = NULL;
  vStar  *star = NULL;

  printf(" COLLAPSING QUAD %d THRU %d \n ", qC, centre);
  printQuadSpecs(qm, qC ) ;
  nq        = 1;
  stat      = EG_backupQuads(qm, &nq, &qC, &quad);
  if (stat != EGADS_SUCCESS || quad == NULL) {
      printf(" In EG_mergeVertices stat in EG_backpQuads %d !!\n ", stat);
      return stat;
  }
  *activity  = 0;
  i          = EG_quadVertIdx(qm, qC, centre);
  if (qm -> vType[centre - 1] != -1) {
      centre = qm -> qIdx[4 * (qC - 1) + (i + 2)%4];
      i      = EG_quadVertIdx(qm, qC, centre);
  }
#ifdef DEBUG
      printf(" EG_mergeVertices collapsing through %d with index at %d = %d \n ",
	     centre, qC, i);
      printQuadSpecs(qm, qC);
      printMesh(qm, NULL);
#endif

  if (qm -> vType[centre - 1] != -1 || i < 0) {
      printf(" EG_mergeVertices collapse through %d is type %d and id in quad %d -> %d !!\n",
	     centre, qm -> vType[centre - 1], qC, i);
      printQuadSpecs(qm, qC);
      return EGADS_INDEXERR;
  }
  if (qm -> valence[centre -1][2] == 2) doublet = 1;
  for (q = 0; q < 4; ++q) {
      oldQ[q    ] = qm -> qIdx[4 * (qC - 1) + (q + i)%4];
      oldQ[q + 4] = qm -> qAdj[4 * (qC - 1) + (q + i)%4];
  }
  for (j  = 0; j < 3; j++) {
      if (j < 2)  {
	  uvxyz[j    ] = qm -> uvs[2 * (oldQ[0] - 1) + j];
	  uvxyz[5 + j] = qm -> uvs[2 * (oldQ[2] - 1) + j];
      }
      uvxyz[2 + j]     = qm -> xyzs[3 * (oldQ[0] - 1) + j];
      uvxyz[7 + j]     = qm -> xyzs[3 * (oldQ[2] - 1) + j];
  }
  // save old quads vertices and adjacents to update map correctly
  q = qm -> valence[ centre - 1][3] - 1;
  if ( doublet == 0 && qm -> vType[oldQ[2] - 1] == -1 &&
      (qm -> vType[q] == -1 || qm -> vType[q] < 4)) {
       stat = EG_centroid(qm, oldQ[2] , uv, 1);
       updateVertex(qm, oldQ[2], uv);
  }
  stat      = EG_buildStar(qm, &star, centre);
  if (stat != EGADS_SUCCESS || star == NULL) {
      printf(" Attention:: building star at %d --> %d !!\n ", centre, stat);
      return stat;
  }
  for (i = 0; i < 4; i++) {
      q = oldQ[4 + i];
      if (q < 0) continue;
      adjq = oldQ[4 + piv[i]]; // - 1 bias
      stat = EG_adjQtoPair(qm, q, oldQ[i], oldQ[(i + 1)%4], adjPair);
      if (stat != EGADS_SUCCESS || adjPair[1] != qC) {
	  EG_freeStar(&star);
	  printf(" STAT %d  adj %d to pair %d %d is %d \n ",
		 stat, q ,oldQ[i], oldQ[(i + 1)%4], adjPair[1]);
	  return  EGADS_INDEXERR;
      }
      qm -> qAdj[4 * (q - 1) + adjPair[0]] = adjq;
  }
  // Eliminate vertex p[0] from all the quads and its valences
  for (i = 0; i < star -> nQ; ++i) {
      q = star -> quads[i];
      if      (q == -1) continue; // ghost quad
      else if (q == qC) {  // eliminate quad qC
	  for (j = 0; j < 4; ++ j) {
	      qm -> qIdx[4 * (qC - 1) + j] = -2;
	      qm -> qAdj[4 * (qC - 1) + j] = -2;
	  }
      } else {
	  for (j = 0; j < 4; ++ j)
	    if (qm -> qIdx[4 * (q - 1) + j] == oldQ[0])
	        qm -> qIdx[4 * (q - 1) + j] = oldQ[2];
      }
  }
  // Point all the collapsed quad vertices to a valid quad
  for (i = 1; i < 4; i++) {
      if (qm -> valence[oldQ[i] - 1][0] != qC) continue;
      for (q = 0; q < 4; q++) {
	  auxQ = oldQ[4 + q];
	  if (auxQ < 0) continue;
	  if (EG_quadVertIdx(qm, auxQ, oldQ[i]) >= 0) {
	      qm -> valence[oldQ[i] - 1][0] = auxQ;
	      break;
	  }
      }
      if (qm -> valence[oldQ[i] - 1][0] == qC) {
	  printf("EG_mergeVertices opposite vertex still points at collapsed quad %d !!\n",
		 EGADS_INDEXERR);
	  EG_freeStar(&star);
	  return EGADS_INDEXERR;
      }
  }
  // Set valences for merged stuff
  for (i = 1; i < 4; i++) {
      stat = setValence(qm, oldQ[i]);
      if (stat != EGADS_SUCCESS) {
	  printf( "EG_mergeVertices merged setValence %d = %d\n", i, stat);
	  EG_freeStar(&star);
	  return stat;
      }
  }
  // set valences to links
  for (i = 0; i < star -> nQ; ++i) {
      j = star -> verts[2 * i + 1];
      if (j < 0 || qm -> vType[j - 1] == -2) continue; // ghost vertex (-1) or collapsed vertex (-2)
      stat = setValence(qm, j);
      if (stat != EGADS_SUCCESS) {
	  printf( "Inside EG_mergeVertices stat in setValence %d = %d\n",
                 i, stat);
	  EG_freeStar(&star);
	  return stat;
      }
  }
  EG_freeStar(&star);
  // delete vertex vC2
  qm -> vType[oldQ[0] - 1]                    = -2; // -2 = removed
  qm -> remQ[++qm -> remQ[0]] = qC;
  qm -> remV[++qm -> remV[0]] = oldQ[0];
  *activity  = 1;
  if (qm ->vType  [oldQ[1] - 1]    == -1 &&
      qm ->valence[oldQ[1] - 1][2] ==  2)
    EG_computeCoords(qm, oldQ[1], 0.0, PI);
  if (qm ->vType  [oldQ[3] - 1]    == -1 &&
      qm ->valence[oldQ[3] - 1][2] ==  2)
    EG_computeCoords(qm, oldQ[3], 0.0, PI);
  if (doublet == 1 || (EG_makeValidMesh(qm, 3, &oldQ[1], 0) == EGADS_SUCCESS)) {
      EG_free(quad);
      return EGADS_SUCCESS;
  }
#ifdef DEBUG
  printf(" COLLAPSE MESH FAILED %d (quad %d vertex %d )!\n", stat, qC, oldQ[0]);
  printMesh(qm, NULL);
#endif
  qm -> vType[oldQ[0] - 1]  = -1;
  stat = EG_restoreQuads(qm, quad, nq);
  EG_free(quad);
  if (stat != EGADS_SUCCESS) {
      printf("EG_mergeVertices :: EG_restoreQuads stat %d \n", stat);
      return stat;
  }
  i = -1;
  if (qm -> remQ[0] > 0)
    i = inList(qm -> remQ[0], &qm -> remQ[1], qC);
  if (i != -1) {
	  for (j = i + 1; j < qm ->remQ[0] - 1; j++) {
		  qm -> remQ[j] = qm -> remQ[j + 1];
		  qm -> remV[j] = qm -> remV[j + 1];
	  }
  }
  --qm -> remQ[0];
  --qm -> remV[0];
  for (i = 0; i < 3; i++) {
      if (i < 2)  {
	  qm -> uvs[2 * (oldQ[0] - 1) + i] = uvxyz[    i];
	  qm -> uvs[2 * (oldQ[2] - 1) + i] = uvxyz[5 + i];
      }
      qm -> xyzs[3 * (oldQ[0] - 1) + i] = uvxyz[2 + i];
      qm -> xyzs[3 * (oldQ[2] - 1) + i] = uvxyz[7 + i];
  }
#ifdef DEBUG
  printMesh(qm, NULL);
#endif
  *activity = 0;
  return EGADS_SUCCESS;
}


static int EG_swap(meshMap *qm, int qIn, int *activity)
{
  int       stat, q, swap = 0;
  quadGroup qg;

  *activity  = 0;
  if (qm -> qIdx[4 * (qIn - 1)] == -2) {
#ifdef DEBUG
      printf("EG_swap for quad %d is %d \n ", qIn, stat);
#endif
      return EGADS_SUCCESS;
  }
  for (swap = q = 0; q < 4; q++) {
      qg.q[0]      = qIn;
      qg.q[1]      = qm -> qAdj[4 * (qIn - 1) + q];
      if (qg.q[1] == -1) continue;
      stat         = EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
      if (stat    != EGADS_SUCCESS) {
	  printf("IN EG_swap -> EG_createQuadGroup  %d !! \n", stat);
	  printQuadGroup(qm, qg);
	  return stat;
      }
      if      (qg.vals[0] <= 4 ||
	  validSwap(qm, qg.verts[0], qg.verts[3]) != 1) continue;
      if      (qg.vals[1] * qg.vals[4] == 9) swap = 1;
      else if (qg.vals[2] * qg.vals[5] == 9) swap = 2;
      else if (qg.vals[3] >= 5 ||
	  (qm -> vType[qg.verts[0] - 1] ==
	      2 * qm -> sizeQ && qg.vals[0] > 5) ||
	      (qm -> vType[qg.verts[3] - 1] ==
		  2 * qm -> sizeQ && qg.vals[3] > 5)) {
	  if  (    qg.vals[1] * qg.vals[4] == 12) swap = 1;
	  else if (qg.vals[2] * qg.vals[5] == 12) swap = 2;
      }
      if (swap != 0) break;
  }
  if (swap  == 0) return EGADS_SUCCESS;
  return EG_swappingOperation(qm, qg, swap, &(*activity));
}


static int
EG_swappingOperation(meshMap *qm, quadGroup qg, int swap, int *activity)
{
  int    area, nq, stat, i0, i1, i, j, q, adj, idx[4], qID[2], adjQmap[6], piv[4], ori[4];
  double angles[4], qsize;
  Quad   *quad = NULL;

  swap      = swap%3;
  *activity = 0;
  if (swap == 0) {
#ifdef DEBUG
      printf(" swapping throu 0-3 will result in the same pair!! \n ");
#endif
      return EGADS_INDEXERR;
  }
  qID[0]     = qg.q[0]; qID[1] = qg.q[1];
  nq         = 2;
  stat       = EG_backupQuads(qm, &nq, qID, &quad);
  if (stat != EGADS_SUCCESS || quad == NULL) {
      EG_free(quad);
      printf(" In EG_swappingOperation :: EG_backupQuad %d !!\n ", stat);
      return stat;
  }
  i0 = EG_quadVertIdx(qm, qID[0], qg.verts[0]); // centre
  i1 = EG_quadVertIdx(qm, qID[1], qg.verts[3]); // opposite
  // Modify Quads and get adj map
  qm -> qIdx[4 * (qID[0] - 1)] = qg.verts[swap];
  qm -> qIdx[4 * (qID[1] - 1)] = qg.verts[swap];
  for (i = 0; i < 3; i++) {
      adjQmap[i    ] = qm -> qAdj[4 * (qID[0] - 1) + (i + i0)%4];
      adjQmap[i + 3] = qm -> qAdj[4 * (qID[1] - 1) + (i + i1)%4];
      qm -> qIdx[4 * (qID[0] - 1) + i + 1] = qg.verts[(swap + 3 + i)%6];
      qm -> qIdx[4 * (qID[1] - 1) + i + 1] = qg.verts[(swap + 1 + i)  ];
  }
  qm -> qAdj[4 * (qID[0] - 1)    ] = qID[1];
  qm -> qAdj[4 * (qID[1] - 1) + 3] = qID[0];
  for (i = 0; i < 3; i++) {
      adj = adjQmap[(3 + i + swap)%6] - 1;
      qm -> qAdj[4 * (qID[0] - 1) + i + 1] = adj + 1;
      if (adj >= 0) {
	  for (j = 0; j < 4; j++) if (qm -> qAdj[4 * adj + j] == qID[1])
	    qm -> qAdj[4 * adj + j] = qID[0];
      }
      adj = adjQmap[i + swap] -1;
      qm -> qAdj[4 * (qID[1] - 1) + i] = adj + 1;
      if (adj >= 0) {
	  for (j = 0; j < 4; j++) if (qm -> qAdj[4 * adj + j] == qID[0])
	    qm -> qAdj[4 * adj + j] = qID[1];
      }
  }
  for (i = 0; i < 4; i++) {
      j = qm -> qIdx[4 * (qID[0] - 1) + i] - 1;
      qm -> valence[j][0] = qID[0];
      j = qm -> qIdx[4 * (qID[1] - 1) + i] - 1;
      qm -> valence[j][0] = qID[1];
  }
  idx[0] = 0; idx[1] = 3; idx[2] = swap; idx[3] = swap + 3;
  for (i = 0; i < 4; i++) {
      stat = setValence(qm, qg.verts[idx[i]]);
      if (stat != EGADS_SUCCESS) {
	  printf(" Inside swapping operation set valence for %d --> %d !!!\n ",
		 qg.verts[idx[i]], stat);
	  EG_free(quad);
	  return stat;
      }
  }
  for (j  = 0; j < 2; j++) {
      i    = qg.verts[j * 3];
      if (qm -> vType[i - 1] != -1) continue;
      q    = 0;
      if (EG_quadVertIdx(qm, qID[q], i) < 0) q = 1;
      area = EG_quadAngleOrientation(qm, qm -> minAngle, qm -> maxAngle, qID[q],
                                     ori, piv, angles, &qsize);
      if (area <  0) return area;
      if (area == 0) EG_computeCoords(qm, i, qm -> minAngle, qm -> maxAngle);
  }
  *activity = 1;
  if (EG_makeValidMesh(qm, 6, qg.verts, 0) == EGADS_SUCCESS) {
      EG_free(quad);
      return EGADS_SUCCESS;
  }
#ifdef DEBUG
  printf(" EG_swapping operation didn't work well: restore quads\n ");
  printMesh(qm, NULL);
#endif
  stat = EG_restoreQuads(qm, quad, nq);
  EG_free(quad);
  if (stat != EGADS_SUCCESS) {
      printf(" EG_swappingOperation :: restore quads %d !!\n ", stat);
      return stat;
  }
#ifdef DEBUG
  printMesh(qm, NULL);
#endif
  *activity = 0;
  return EGADS_SUCCESS;
}


static int EG_split(meshMap *qm, int qID, int *activity)
{
  int   poly[3], val[3], v, q, id0 = 0, i, stat, dist = 0, validSplit = 0;
  int   id6[2], links[4];
  vStar *star = NULL;

  *activity   = 0;
  if (qm -> qIdx[4 * (qID - 1)] == -2) {
#ifdef DEBUG
      printf(" EG_split for quad %d is %d \n ", qID, stat);
#endif
      return EGADS_SUCCESS;
  }
  for (v = 0; v < 4; v++) {
      poly[0]   = qm -> qIdx[4 * (qID - 1) + v];
      val[0]    = qm -> valence[ poly[0] - 1][1];
      if (val[0] < 5 || qm -> vType[poly[0] - 1] == 4) continue;
      stat      = EG_buildStar(qm, &star, poly[0]);
      if (stat != EGADS_SUCCESS || star == NULL) return stat;
      id6[0]    = -1;
      id6[1]    = -1;
      for (validSplit = q = 0; q < star -> nQ; ++q) {
	  if (star -> quads[q] == - 1) continue;
	  id0     = 2 * q + 1;
	  poly[1] = star -> verts[id0];
	  val[1]  = qm -> valence[ poly[1] -1][1];
	  for (i = 0; i <= qm -> valence[poly[0] - 1][2] - 4; ++i) {
	      dist    = 4 + 2*i;
	      poly[2] = star -> verts[star -> idxV[id0 + dist]];
	      if (poly[2] < 0) continue; // ghost vertex
	      val[2]  = qm -> valence[ poly[2] -1][1];
	      if (val[1] == 3 && val[2] == 3 &&
		  qm -> vType[poly[0] - 1] == -1) {
		  if ((star -> nQ == 6 && dist == 6) || star -> nQ != 6) {
		      validSplit = 1;
		      break;
		  }
	      }
	      else if (val[1] <= 3 && val[2] <= 4 && id6[0] == -1) id6[0] = id0;
	      else if (val[2] <= 3 && val[1] <= 4 && id6[0] == -1)
		id6[0] = star -> idxV[id0 + dist];
	      else if (val[1] <= 4 && val[2] <= 4 && id6[0] == -1) id6[1] = id0;
	  }
	  if (validSplit == 1) break;
      }
      if (validSplit == 0 && val[0] >= 6 && (id6[0] != -1 || id6[1] != -1)) {
	  validSplit  = 1;
          /* boundary vertex: Special split since regular = 3 ( not val 4)  */
	  if (qm -> vType[poly[0] - 1] == 3) {
	      for (q = 0; q < star -> nQ; q++)
		if (star -> quads[q] == -1) break;
	      id6[0] = - 1; id6[1] = -1;
	      for (i = 0; i < 2; i++) {
		  id0 = star -> idxQ[q + i];
		  if (i == 1) dist = (star -> nV - 1) - 4;
		  else        dist = 4 + 2 * i;
		  links[2 * i    ] = 2 * id0 + 1;
		  links[2 * i + 1] = star -> idxV[2 * id0 + 1 + dist];
		  if (  qm -> valence[ star->verts[ links[2*i     ]] - 1][1] == 4) {
		    if (qm -> valence[ star->verts[ links[2*i + 1 ]] - 1][1] == 3 &&
		             id6[0] == -1) id6[0] = i;
		    else if (qm -> valence[ star->verts[ links[2*i + 1  ]] -1][1] == 4  &&
		  	     id6[1] == -1) id6[1] = i;
		  }
	      }
	      dist    = 4;
	      if      (id6[0] != -1) id0 = links[3 * id6[0]];
	      else if (id6[1] != -1) id0 = links[3 * id6[1]];
	      else  validSplit = 0;
	  } else {
	      dist = 6;
	      id0  = id6[0]; if (id0 < 0) id0 = id6[1];
	  }
      }
      if (validSplit == 1)  {
	  stat = EG_splittingOperation(qm, star -> verts[0], star -> verts[id0],
				       star -> verts[star -> idxV[id0 + dist]],
				       &(*activity));
	  EG_freeStar(&star);
	  return stat;
      }
      EG_freeStar(&star);
  }
  return EGADS_SUCCESS;
}


static int EG_splittingOperation(meshMap *qm, int vC, int vL, int vR,
                                 int *activity)
{
  int   qIdx[4], modQ[4], verts[4], adj[2], poly[4], q, newQ, i, j, stat;
  int   id0 = -1, id1 = -1, dist, links[4], vals[4], addedV = 0, nq;
  double uv[2];
  vStar *star = NULL;
  Quad  *quad = NULL;

  *activity   = 0;
  stat        = EG_buildStar(qm, &star, vC);
  if (stat   != EGADS_SUCCESS || star == NULL) {
      printf(" In splittingOperation build star %d is NULL %d !!\n", vC, stat);
      return stat;
  }
  nq        = star -> nQ;
  stat      = EG_backupQuads(qm, &nq, star -> quads, &quad);
  if (stat != EGADS_SUCCESS || quad == NULL) {
      EG_free(quad);
      printf("In EG_splittingOperation :: EG_backupQuads %d !!\n", stat);
      return stat;
  }
  if (qm -> remQ[0] > 0) {
      poly[3] = qm -> remV[qm -> remV[0]--];
      newQ    = qm -> remQ[qm -> remQ[0]--];
  }
  else {
      poly[3] = qm -> totV + 1;
      newQ    = qm -> totQ + 1;
      if (poly[3] >  qm -> sizeV) {
	  printf(" We have duplicated the number of initial vertices. This is too much. \n");
	  EG_free(quad);
	  return EGADS_INDEXERR;
      }
      ++qm -> totV; ++qm -> totQ;
      addedV = 1;
  }
  qm -> vType[poly[3] - 1] = -1;
  if (qm -> vType[vC - 1] != -1) {
      for (q = 0; q < star -> nQ; q++)
	if (star -> quads[q] == -1) break;
      for (i = 0; i < 2; i++) {
	  if (i == 1) dist = (star -> nV - 1) - 4;
	  else        dist = 4;
	  id0 = star -> idxQ[q + i];
	  links[2 * i    ] = star -> verts[star -> idxV[2 * id0 + 1       ]];
	  links[2 * i + 1] = star -> verts[star -> idxV[2 * id0 + 1 + dist]];
	  vals [2 * i    ] = qm -> valence [ links[2*i    ] - 1][1];
	  vals [2 * i + 1] = qm -> valence [ links[2*i + 1] - 1][1];
      }
      dist = 4;
      i    = 0;
      if (vals[0] * vals[1] > vals[2] * vals[3]) {
	  vL = links[3];
	  vR = links[2];
      } else {
	  vL = links[0];
	  vR = links[1];
      }
  }
  id0 = - 1; id1 = -1;
  for (j = 0; j < star -> nQ; j++) {
      if (star -> verts[2 * j + 1] == vL) id0 = j;
      if (star -> verts[2 * j + 1] == vR) id1 = j;
  }
  if (id0 == -1 || id1 == -1) {
      EG_freeStar(&star);
#ifdef DEBUG
      printf("EG_splittingOperation  Problem finding vertices %d !!\n ",
             EGADS_INDEXERR);
#endif
      EG_free(quad);
      return EGADS_INDEXERR;
  }
  poly [0] = star -> verts[0];
  poly [1] = star -> verts[2 * id0 + 1];
  poly [2] = star -> verts[2 * id1 + 1];
  qIdx [0] = id0;
  qIdx [1] = star -> idxQ[id1 + star -> nQ - 1];
  qIdx [2] = id1;
  qIdx [3] = star -> idxQ[id0 + star -> nQ - 1];
  verts[0] = poly[1];
  verts[1] = poly[2];
  verts[2] = poly[2];
  verts[3] = poly[1];
  qm -> qIdx[4 * (newQ - 1)    ] = poly[1];
  qm -> qIdx[4 * (newQ - 1) + 1] = poly[0];
  qm -> qIdx[4 * (newQ - 1) + 2] = poly[2];
  qm -> qIdx[4 * (newQ - 1) + 3] = poly[3];
  for (i = 0; i < 4; ++i) {
      modQ[i] = star -> quads[qIdx[i]];
      qm -> qAdj[4 * (newQ - 1) + i] = modQ[i];
      qm -> valence[qm -> qIdx[4*(newQ - 1)+i]-1][0] = newQ;
      if (modQ[i] == -1) continue;
      stat = EG_adjQtoPair(qm, modQ[i], poly[0], verts[i], adj);
      qm -> qAdj[4 * (modQ[i] - 1) + adj[0]] = newQ;
  }
  j = qIdx[2];
  q = star -> quads[j++];
  while (q != star -> quads[qIdx[0]]) {
      for (i = 0; i < 4; ++i)
	if (qm -> qIdx[4 * (q - 1) + i] == poly[0])
	    qm -> qIdx[4 * (q - 1) + i] = poly[3];
      q = star -> quads[star -> idxQ[j++]];
      if (q == -1) {
	  if (star -> quads[qIdx[0]] == -1) break;
	  else q = star -> quads[star -> idxQ[j++]];
      }
  }
  // Add valences to splitting vertices
  for (i = 0; i < 4; i++) {
      stat = setValence(qm, poly[i]);
      if (stat != EGADS_SUCCESS) {
	  EG_freeStar(&star);
	  EG_free (quad);
	  return stat;
      }
  }
  // Update valences at links
  for (i = 0; i < 4; i++) {
      for (j = 0; j < qm -> valence[poly[i] - 1][2]; j++) {
	  stat = setValence(qm, qm -> valence[poly[i] - 1][3 + j]);
	  if (stat != EGADS_SUCCESS) {
	      EG_freeStar(&star);
	      EG_free(quad);
	      return stat;
	  }
      }
  }
  EG_freeStar(&star);
  updateVertex(qm, poly[3], &qm -> uvs[2 * (poly[0] - 1)]);
  if ( qm -> vType [ poly[0] - 1] == -1 )
    {
      EG_centroid (qm, poly[0], uv, 1);
      updateVertex(qm, poly[0], uv);
    }
  if ( qm -> vType [ poly[3] - 1] == -1 )
    {
      EG_centroid (qm, poly[3], uv, 1);
      updateVertex(qm, poly[3], uv);
    }
  EG_computeCoords(qm, poly[0], qm -> minAngle, qm -> maxAngle);
  EG_computeCoords(qm, poly[3], qm -> minAngle, qm -> maxAngle);
  *activity = 1;
  if (EG_makeValidMesh(qm, 4, poly, 0) == EGADS_SUCCESS) {
      EG_free(quad);
      return EGADS_SUCCESS;
  }
#ifdef DEBUG
  printf("EG_splittingOperation mesh valid failed: restore quads\n");
  printMesh(qm, NULL);
#endif
  stat = EG_restoreQuads(qm, quad, nq);
  EG_free(quad);
  if (stat != EGADS_SUCCESS) {
      printf(" EG_splittingOperation :: EG_restoreQuads %d !!\n ", stat);
      return stat;
  }
  if (addedV == 1) {
      --qm -> totV;
      --qm -> totQ;
  } else {
      qm -> remV[++qm -> remV[0]] = poly[3];
      qm -> remQ[++qm -> remQ[0]] = newQ;
  }
  qm -> vType[poly[3] - 1] = -2;
  for (i = 0; i < 4; i++)
    qm -> qIdx[4 * (newQ - 1) + i] = -2;
#ifdef DEBUG
  printMesh(qm, NULL);
#endif
  *activity = 0;
  return stat;
}


static int EG_doubleSwap(meshMap *qm, quadGroup qg, int forcing, int *activity)
{
  int piv5 = -1, piv3 = -1, q5, i, adjPiv5, stat, adj[2], swap = 0, vopp3;

  *activity = 0;
  if      (validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
  else if (qg.vals[0] == 4 && forcing == 0) return EGADS_SUCCESS;
  else if (qg.vals[0] == 5) {
      if  (qg.vals[2] * qg. vals[4] == 15) {
	  piv3 = 4;
	  if (qg.vals[2] == 3) piv3 = 2;
	  swap = 0;
      }
      else if (qg.vals[1] == 3 || qg.vals[5] == 3) {
	  piv3 = 1; if (qg.vals[1] != 3) piv3 = 5;
	  swap = 1;
      }
  }
  else if (forcing == 1) {
      if  (qg.vals[1] != 3 && qg.vals[5] != 3) return EGADS_SUCCESS;
      piv3 = 1; if (qg.vals[1] != 3) piv3 = 5;
      if (qg.vals[1] != 5 && qg.vals[5] != 5) swap = 1;
  }
  if (piv3 == -1) return EGADS_SUCCESS;
  vopp3 = (piv3 + 3)%6;
  if (qm -> valence[qg.verts[vopp3] -1][1] != 4) return EGADS_SUCCESS;
  if (swap == 0) {
      piv5 = (vopp3 + 1)%6;
      if (piv5%3 == 0) piv5 = (vopp3 + 5)%6;
  }
  else {
      piv5 = (vopp3 + 1)%6;
      if (piv5%3 != 0) piv5 = (vopp3 + 5)%6;
  }
  q5   = 0;
  if (vopp3 > 3) q5 = 1;
  stat = EG_adjQtoPair(qm, qg.q[q5], qg.verts[vopp3], qg.verts[piv5], adj);
  if (stat != EGADS_SUCCESS || adj[1] == -1) {
#ifdef DEBUG
      if (stat != EGADS_SUCCESS)
	printf(" EG_doubleSwap adjToPair %d --> %d !!\n", adj[1], stat);
#endif
      return stat;
  }
  i       = EG_quadVertIdx(qm, adj[1], qg.verts[vopp3]);
  adjPiv5 = qm -> qIdx[4 * (adj[1] - 1) + (i + 1)%4];
  if (adjPiv5 == qg.verts[piv5])
    adjPiv5 = qm -> qIdx[4 * (adj[1] - 1) + (i + 3)%4];
  if (      swap == 0 && qm -> valence[adjPiv5 -1][1] > 4) return EGADS_SUCCESS;
  else if ((swap == 1 && qm -> valence[adjPiv5 -1][1] < 5) ||
           (qm -> vType[adjPiv5 - 1] != -1 &&
	    qm -> vType[qg.verts[vopp3] - 1] != -1)) return EGADS_SUCCESS;
  piv5      = qg.verts[0];
  stat      = EG_swappingOperation(qm, qg, piv3, &(*activity));
  if (stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
      printf(" EG_doubleSwap: at first swap activity %d, stat --> %d !!\n ",
             *activity, stat);
#endif
      return stat;
  }
  stat = EG_swap(qm, adj[1], &i);
  if (i == 0 && stat == EGADS_SUCCESS) {
      stat      = EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
      for (swap = 0; swap < 6; swap++) if (qg.verts[swap] == piv5) break;
      stat      = EG_swappingOperation(qm, qg, swap, &i);
      if (stat == EGADS_SUCCESS && i == 1) *activity = 0;
      else {
#ifdef DEBUG
	  printf("In EG_doubleSwap undoing first swap activity %d stat %d \n!!",
		 i, stat);
#endif
	  return stat;
      }
  }
  return stat;
}


static int EG_doubleSplit(meshMap *qm, quadGroup qg, int forcing, int *activity)
{
  int i, j, stat, piv[2] = {1, 5};

  *activity = 0;
  if (qg.vals[1] != 3) {
      piv[0] = 5; piv[1] = 1;
  }
  if (qg.vals[piv[0]] != 3) return EGADS_SUCCESS;
  if ((qm -> vType[qg.verts[0]      - 1] >= 0 &&
       qm -> vType[qg.verts[piv[1]] - 1] >= 0) ||
       qm -> vType[qg.verts[0]      - 1] == 4) return EGADS_SUCCESS;
  if ((forcing == 0 && (qg.vals[0] < 5 || qg.vals[1] * qg.vals[5]  != 15)) ||
      (forcing == 1 &&  qg.vals[0] * qg.vals[piv[1]] <= 16)) return EGADS_SUCCESS;
  stat         = EG_splittingOperation(qm, qg.verts[0], qg.verts[piv[0]],
			       qg.verts[piv[1]], &(*activity));
  if (stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
      printf("In EG_doubleSplit: force 1st split through %d - %d --> %d activity %d !!\n ",
	     qg.verts[0], qg.verts[piv[0]], stat, *activity);
#endif
      return stat;
  }
  for (j    = 0; j < 2; j++)
    if (EG_quadVertIdx(qm, qg.q[j], qg.verts[piv[1]]) >= 0) break;
  stat      = EG_split(qm, qg.q[j], &i);
  if (stat == EGADS_SUCCESS && i == 0) {
      j     = qm -> valence[qg.verts[0] - 1][0];
      stat  = EG_collapse(qm, j, &i, 1);
      if (stat == EGADS_SUCCESS && i == 1) *activity = 0;
      return stat;
  }
  return stat;
}


static int EG_swapSplit(meshMap *qm, quadGroup qg, int forcing, int *activity)
{
  int  stat, i, j, i3 = -1, i5 = -1, v3opp = -1, q5, vL5, vL5adj, swap = 0, adj[2], v0;
  vStar *star = NULL;

  *activity = 0;
  if  (qg.vals[0] * qg.vals[3] != 20 ||
      validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
  if      (qg.vals[1] == 3) i3 = 1;
  else if (qg.vals[5] == 3) i3 = 5;
  for (i = 1; i < 6; i++)
    if (qg.vals[i] == 5) i5 = i;
  if (i3 != -1 && (i5 == -1 || (i5 != -1 && i5 == (i3 + 3)%6))) {
      v3opp = (i3 +3)%6;
      i5    = v3opp;
  }
  else if (i3 == -1 && forcing == 1 && (i5 == -1 || i5 == 2 || i5 == 4)) {
      if ( i5 == -1) v3opp = -1;
      else           v3opp = i5;
  }
  else return EGADS_SUCCESS;
  if (v3opp == -1) {
      for (i  = 0; i < 2; i++) {
	  j = 2 + 2 * i;
	  if (i == 0) vL5 = 1;
	  else        vL5 = 5;
	  if (qm -> valence[qg.verts[j] -1][1]  == 3) continue;
	  stat = EG_adjQtoPair(qm, qg.q[i], qg.verts[j],
			       qg.verts[vL5], adj);
	  if (stat != EGADS_SUCCESS || adj[1] == -1) continue;
	  q5 = EG_quadVertIdx(qm, adj[1], qg.verts[j]);
	  vL5adj = qm -> qIdx[4 * (adj[1] - 1) + (q5 + 1)%4];
	  if (vL5adj == qg.verts[vL5])
	    vL5adj = qm -> qIdx[4 * (adj[1] - 1) + (q5 + 3)%4];
	  if (qm -> valence[vL5adj - 1][1] == 3) {
	      i5   = j;
	      swap = j;
	      break;
	  }
      }
  } else {
      vL5 = (v3opp + 1)%6;
      if (vL5%3 == 0) vL5 = (v3opp + 5)%6;
      q5 = 0;
      if (EG_quadVertIdx(qm, qg.q[q5], qg.verts[vL5]) < 0) q5 = 1;
      stat = EG_adjQtoPair(qm, qg.q[q5], qg.verts[v3opp],
			   qg.verts[vL5], adj);
      if (stat     != EGADS_SUCCESS || adj[1] == -1) {
	  if (stat != EGADS_SUCCESS)
	    printf("EG_swapSplit: EG_adjQtoPair from quad %d is %d\n!!",
		   qg.q[q5], stat);
	  return stat;
      }
      i           = EG_quadVertIdx(qm, adj[1], qg.verts[v3opp]);
      vL5adj      = qm -> qIdx[4 * (adj[1] - 1) + (i + 1)%4];
      if (vL5adj == qg.verts[vL5])
	vL5adj  = qm -> qIdx[4 * (adj[1] - 1) + (i + 3)%4];
      if (i3 != -1 && (qg.vals[v3opp] == 5 || qm -> valence[vL5adj - 1][1] == 3))
	swap = i3;
      else if (forcing == 1 && (qg.vals[v3opp] == 5 || qm -> valence[vL5adj - 1][1] == 3))
	swap = v3opp;
  }
  if (swap%3 == 0 || (qm -> vType[qg.verts[swap] - 1] >= 0 &&
      qm -> vType[qg.verts[(swap + 3)%6] - 1] >= 0) ||
      qm -> vType[qg.verts[i5] - 1] == 4) return EGADS_SUCCESS;
  stat         = EG_swappingOperation(qm, qg, swap, &(*activity));
  v0           = qg. verts[0];
  if (stat    != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
      printf(" In swapSplit thru %d : EG_swappingOperation went %d activity %d !!\n ",
	     swap, stat, *activity);
      printQuadGroup(qm, qg);
#endif
      return stat;
  }
  stat      = EG_buildStar(qm, &star, qg.verts[i5]);
  if (stat != EGADS_SUCCESS || star == NULL) {
      printf("In swapSplit build star for %d --> %d!!\n", qg.verts[i5], stat);
      return stat;
  }
  for (i = 0; i < star -> nQ; i++)
    if (star -> verts[2 * i + 1] == qg.verts[3]) break;
  if (star -> verts[2 * i + 1] != qg.verts[3]) {
      EG_freeStar(&star);
      return EGADS_INDEXERR;
  }
  adj[0] = star -> verts[star -> idxV[2 * i + 1 + 6]];
  if (star -> nQ == 5) {
      adj[1] = star -> verts[star -> idxV[2 * i + 1 + 4]];
      if (qm -> valence[adj[1] - 1][1] < qm -> valence[adj[0] - 1][1])
	adj[0] = adj[1];
  }
  EG_freeStar(&star);
  stat         = EG_splittingOperation(qm, qg.verts[i5], qg.verts[3], adj[0], &i);
  if (stat    != EGADS_SUCCESS || i == 0) {
#ifdef DEBUG
      printf("In swapSplit splittingOperation --> %d !!\n", stat);
#endif
      stat      = EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
      for (swap = 0; swap < 6; swap++) if (qg.verts[swap] == v0) break;
      stat      = EG_swappingOperation(qm, qg,swap, &i);
      if (stat == EGADS_SUCCESS && i == 1) *activity = 0;
  }
  return stat;
}


static int EG_swapCollapse(meshMap *qm, quadGroup qg, int forcing, int *activity)
{
  int  stat, i, i3 = -1, q5, qC, vL5, vL5adj, vOpp3, swap = 0, adj[2], v0;

  *activity = 0;
  if (validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
  if (qg.vals[1] * qg.vals[2] == 15) {
      i3 = 1;
      if (qg.vals[i3] == 5) i3 = 2;
  }
  else if (qg.vals[4] * qg.vals[5] == 15) {
      i3 = 4;
      if (qg.vals[i3] == 5) i3 = 5;
  }
  else if (qg.vals[0] * qg.vals[3] == 20) {
      if      (qg.vals[1] == 3) i3 = 1;
      else if (qg.vals[5] == 3) i3 = 5;
      else return EGADS_SUCCESS;
      if (i3 == 1) {
	  q5   = 0;
	  stat = EG_adjQtoPair(qm, qg.q[0], qg.verts[2], qg.verts[3], adj);
      }
      else {
	  q5   = 1;
	  stat = EG_adjQtoPair(qm, qg.q[1], qg.verts[4], qg.verts[3], adj);
      }
      if (stat != EGADS_SUCCESS || adj[1] == -1) {
	  if (stat != EGADS_SUCCESS)
	    printf("EG_swapCollapse centre 4-5 adj to pair %d !!\n", stat);
	  return stat;
      }
      stat = EG_createQuadGroup(qm, &qg, qg.q[q5], adj[1]);
      if (stat != EGADS_SUCCESS) {
	  printf("Inside EG_swapCollapse :: EG_createGroup %d !!\n ", stat);
	  return stat;
      }
      return EG_swapCollapse(qm, qg, forcing, &(*activity));
  }
  else if (forcing == 1) {
      if (qg.vals[1] * qg.vals[5] == 25 || qg.vals[2] * qg.vals[4] == 25) {
	  i3 = 2; if (qg.vals[1] != 5) i3 = 1;
      }
      else if (qg.vals[1] * qg.vals[4] == 15) {
	  i3 = 1; if (qg.vals[1] != 3) i3 = 4;
      }
      else if (qg.vals[2] * qg.vals[5] == 15) {
	  i3 = 2; if (qg.vals[2] != 3) i3 = 5;
      }
      else if (qg.vals[1] * qg.vals[5] == 9 || qg.vals[2] * qg.vals[4] == 9) {
	  v0        = qg.verts[0];
	  stat      = EG_swappingOperation(qm, qg, 1, &(*activity));
	  if (stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
	      printf("forcing swapcollapse:: after swapping %d activity %d !!\\n ",
		     stat, *activity);
#endif
	      return stat;
	  }
	  qC = qg.q[0];
	  if (EG_nValenceCount(qm, qC, 3) < 2) qC = qg.q[1];
	  stat = EG_collapse(qm, qC, &(*activity), 1);
	  if (stat  != EGADS_SUCCESS || *activity == 0) {
	      stat   = EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
	      for (swap = 0; swap < 6; swap++) if (qg.verts[swap] == v0) break;
	      stat         = EG_swappingOperation(qm, qg, swap, &i);
	      if (stat    == EGADS_SUCCESS && i == 1) *activity = 0;
	  }
	  return stat;
      }
      else return EGADS_SUCCESS;
  } else return EGADS_SUCCESS;
  q5    = 0; if (i3 > 3) q5 = 1;
  if (validCollapse(qm, qg.q[q5], qg.verts[i3]) == 0) return EGADS_SUCCESS;
  vOpp3 = (i3 + 3)%6;
  qC    = qg.q[q5];
  vOpp3 = (i3 + 3)%6;
  vL5   = (vOpp3 + 1)%6;
  if (vL5%3 != 0) vL5 = (vOpp3 + 5)%6;
  if (qg.vals[vOpp3] == 3) return EGADS_SUCCESS;
  stat = EG_adjQtoPair(qm, qg.q[q5], qg.verts[vOpp3], qg.verts[vL5], adj);
  if (stat != EGADS_SUCCESS || adj[1] == -1) return stat;
  i      = EG_quadVertIdx(qm, adj[1], qg.verts[vOpp3]);
  vL5adj = qm -> qIdx[4 * (adj[1] - 1) + (i + 1)%4];
  if (vL5adj == qg.verts[vL5]) vL5adj = qm->qIdx[4 * (adj[1] - 1) + (i + 3)%4];
  if ((forcing == 0 && qg.vals[vL5] == 4 && qm -> valence[vL5adj - 1][1] >= 4) ||
      (forcing == 1 && qm -> valence[vL5adj - 1][1] > 4)) return EGADS_SUCCESS;
  stat = EG_createQuadGroup(qm, &qg, qg.q[q5], adj[1]);
  if (stat != EGADS_SUCCESS) {
      printf("EG_swapCollapse before swap: EG_createQuadGroup --> %d !!\n", stat);
      printQuadGroup(qm, qg);
      return stat;
  }
  if (validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
  for (swap = 0; swap < 6; swap++)
    if (qg.verts[swap] == vL5adj) break;
  v0        = qg.verts[0];
  stat      = EG_swappingOperation(qm, qg, swap, &(*activity));
  if (stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
      printf("EG_swapCollapse after swapping %d activity %d !!\n",
	     swap, *activity);
#endif
      return stat;
  }
  vL5       = qg.verts[0];
  stat      = EG_collapse(qm, qC, &(*activity), 1);
  if (stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
      printf("EG_swapCollapse after collapsing %d !!\n", stat);
#endif
      stat      = EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
      for (swap = 0; swap < 6; swap++) if (qg.verts[swap] == v0) break;
      stat      = EG_swappingOperation(qm, qg, swap, &i);
      if (stat == EGADS_SUCCESS && i == 1) *activity = 0;
  }
  return stat;
}


static int EG_cleanNeighborhood(meshMap *qm, int qID,  int transfer,
                                int *activity)
{
  int   i, act, stat = 0, count = 0, j, v[4];
  vStar *star = NULL;

  *activity = 0;
  if (qm -> qIdx[4 * (qID - 1)] == -2) return EGADS_SUCCESS;
  for (i = 0; i < 4; i++) v[i]  = qm -> qIdx[4 * (qID - 1) + i];
  stat      = EG_cleanQuad(qm, qID, 1, transfer, 0, &act);
  if (stat != EGADS_SUCCESS) return stat;
  *activity += act;
  for (i = 0; i < 4; i++) {
      if (qm -> vType[v[i] - 1] == -2) continue;
      stat      = EG_buildStar(qm, &star, v[i]);
      if (star == NULL || stat != EGADS_SUCCESS) {
	  printf("EG_cleanNeighborhood : EG_buildStar  = %d \n", stat);
	  return stat;
      }
      for (j = 0; j < star -> nQ; j++) {
	  if (star -> quads[j] == -1) continue;
	  if (EG_cleanQuad(qm, star->quads[j], 1, transfer, 0, &act) != EGADS_SUCCESS) {
	      EG_freeStar(&star);
	      return stat;
	  }
	  count += act;
      }
  }
  EG_freeStar(&star);
  *activity += count;
#ifdef DEBUG
  return checkMesh(qm);
#else
  return EGADS_SUCCESS;
#endif
}


int EG_fullMeshRegularization(meshMap *qm)
{
  int    i, j, k, q, stat;
  int    ITMAX = 20 , it = 0, activity = 0, totActivity = 0, loopact;
  int    iV, quadPair[2], prevPair[2], totV, vQ, transfer = 0, ori[4], piv[4];
  double angle[4], seg, qsize;
  vStar  *star = NULL;
#ifdef IOFILES
  char   buffer[100];
#endif

  // GET RANGE FOR EACH POINT
  meshCount(qm, &iV, &totV, &vQ);
  printf("\n\n Face %d --> original mesh has %d QUADS, %d / %d irregular vertices (%.2f %%) \n",
	 qm -> fID, vQ, iV, totV, (double) iV * 100.0 / (double)totV);
  // PREPROCESS: Reduce the number of quads by forcing collapses everywhere
  qm -> pp = 1;
  for (it = 0; it < ITMAX; it++) {
      totActivity = 0;
#ifdef DEBUG
      printf(" ROUND %d  BAR %lf \n ", it + 1,  qm -> avsize * qm -> avsize * (double) (it + 1)  / (double) ( ITMAX )) ;
      printMesh (qm, NULL );
      resizeQm(qm ) ;
      EG_wvsData (qm, NULL);
#endif
      for (i = 0; i < qm -> totQ; i++) {
	  if ( qm -> qIdx [ 4 * i ] == -2 ) continue;

	  q  = EG_quadAngleOrientation (qm, DEG20, PI, i + 1,
					ori, piv, angle, &qsize);
	  if ( q != 1 || qsize < qm -> avsize * qm -> avsize * (double) (it + 1)  / (double) ( ITMAX ) ) {
	      printf(" quad %d is below average %lf < %lf ! \n ", i + 1,
		     qsize, qm -> avsize * qm -> avsize * (double) (it + 1)  / (double) ( ITMAX ));
	      stat = EG_collapse (qm, i + 1, &activity, 1);
	  }
	  for ( q = 0 ; q < qm -> totQ; q++ ) {
	      if ( qm -> qIdx [ 4 * q ] == -2 ) continue;
	      for ( k = 0 ; k < 4; k++ ) {
		  if ( qm -> valence [ qm -> qIdx [ 4 * q + k] - 1][2] == 1 ){
		      printf("FOUND IT!!!!!!!!!!!!!!\n ");
		      printQuadSpecs (qm, q + 1 ) ;
		      exit (1);
		  }
	      }

	  }
	  totActivity += activity;
	  if (stat  != EGADS_SUCCESS) {
	      printf("In EG_cleanMesh: EG_collapse for quad %d --> %d!!\n ",
		     i + 1, stat);
	      return stat;
	  }
	  meshCount(qm, &iV, &totV, &vQ);
	  if (vQ <= qm -> estQ) break;
      }
      meshCount(qm, &iV, &totV, &vQ);
      printf(" ACTIVITY %d \n ", totActivity ) ;
      if ( iV <= 2 || vQ <= qm -> estQ) break;
  }
  qm -> pp = 0;
  stat      = resizeQm(qm );
  if (stat != EGADS_SUCCESS) {
        printf(" After resizing mesh:: %d \n", stat);
        return stat;
    }
  EG_makeValidMesh(qm, 0, NULL, 2);
#ifdef IOFILES
  snprintf(buffer,100, "gnuPreprocess_%d.txt", qm -> fID);
  printMesh (qm, buffer);
  snprintf(buffer,100, "wvsPreprocess_%d.txt", qm -> fID);
  EG_wvsData (qm, buffer);
  snprintf(buffer,100, "preprocessMeshStats_%d.txt", qm -> fID);
  printMeshStats(qm, 10);
  printf("\n\n After Preprocessing Face %d --> has %d QUADS  (EST %d ), %d / %d irregular vertices (%.2f %%) \n\n",
	 qm -> fID, vQ, qm -> estQ, iV, totV, (double) iV * 100.0 / (double)totV);
#endif
  for (it = 0; it < ITMAX; it++) {
      totActivity = 0;
      for (i = 0; i < qm -> totQ; i++) {
	  if (qm -> qIdx[4 * i] == -2) continue; //can be a deleted quad
	  stat         = EG_cleanNeighborhood(qm, i + 1, 0, &activity);
	  totActivity += activity;
	  if (stat  != EGADS_SUCCESS) {
	      printf("In EG_cleanMesh: EG_CleanNeighborhood for quad %d --> %d!!\n ",
		     i + 1, stat);
	      return stat;
	  }
      }
      meshCount(qm, &iV, &totV, &vQ);
      if ( totActivity == 0 || iV <= 2) break;
  }
#ifdef DEBUG
  printf(" TRANSFER VALENCES \n ");
#endif
  if (iV > 2) {
      for (k  = 0; k <= 1; k++) {
	  it = 0;
#ifdef DEBUG
	  if (k == 0) printf("Translating valences around mesh round 1\n ");
	  else        printf("Translating valences around mesh round 2 :: forcing \n ");
#endif
	  do {
	      it++;
	      totActivity  =  0;
	      prevPair[0]  = -1;
	      prevPair[1]  = -1;
	      for (q       = 0; q < qm ->  totQ; q++) {
		  transfer = 0;
		  if (qm -> qIdx[4 * q] == -2 || q + 1 == prevPair[0] ||
		      q + 1 == prevPair[1]) continue;
		  quadPair[0]   = q + 1;
		  quadPair[1]   = -1;
		  stat          = EG_transferValences(qm, quadPair, k,
						      &transfer, &activity);
		  if (stat     != EGADS_SUCCESS) return stat;
		  if (activity == 0 || quadPair[0] < 0 ||
		      qm -> qIdx[4*(quadPair[0]-1)] == -2) continue;
		  totActivity += activity;
		  for (j = 0; j < 20; j++) {
		      activity = 0;
		      for (i  = 0; i < 2; i++) {
			  if (quadPair[i] < 0 ||
			      qm -> qIdx[4*(quadPair[i]-1)] == -2) continue;
			  stat      = EG_cleanNeighborhood(qm, quadPair[i],
                                                           transfer, &loopact);
			  if (stat != EGADS_SUCCESS) return stat;
			  activity +=loopact;
		      }
		      if (activity > 0) break;
		      stat      = EG_transferValences(qm, quadPair, k,
                                                      &transfer, &activity);
		      if (stat != EGADS_SUCCESS) return stat;
		      if (activity == 0 || quadPair[0] < 0 ||
			  qm -> qIdx[4*(quadPair[0]-1)] == -2) break;
		  }
		  if (iV <= 2) break;
		  prevPair[0] = quadPair[0];
		  prevPair[1] = quadPair[1];
	      }
	      if (iV <= 2) break;
	      for (q = 0; q < qm ->  totQ; q++) {
		  if (qm -> qIdx[4 * q] == -2) continue;
		  stat = EG_cleanQuad(qm, q + 1, 1, 0, 0, &i);
		  if (stat != EGADS_SUCCESS) {
		      printf(" In EG_fullMeshRegularization clean quad %d !!\n ",
			     stat);
		      return stat;
		  }
		  totActivity += i;
	      }
	      meshCount(qm, &iV, &totV, &vQ);
	      if (iV <= 2) break;
	  } while (totActivity > 0 && it < ITMAX && iV > 2);
#ifdef DEBUG
	  stat = checkMesh(qm);
	  if (stat != EGADS_SUCCESS) return stat;
#endif
      }
  }
  meshCount(qm, &iV, &totV, &vQ);
  printf("\n\n*******************************************************************\n");
  printf(" Final mesh has %d QUADS, %d / %d irregular vertices (%.2f %%) ============\n",
	 vQ, iV, totV, (double) iV * 100.0 / (double) totV);
  printf(" Total invalid steps  %d\n ", qm -> invsteps);
  printf("*******************************************************************\n");
  stat      = resizeQm(qm );
  if (stat != EGADS_SUCCESS) {
      printf(" After resizing mesh:: %d \n", stat);
      return stat;
  }
#ifdef IOFILES
  snprintf(buffer,100, "wvsRegular_%d.txt", qm -> fID);
  stat = EG_wvsData(qm, buffer);
  if (stat != EGADS_SUCCESS) {
      printf(" writing in wvs file %d !! \n ", stat);
      return stat;
  }
  snprintf(buffer,100,"gnuRegular_%d.txt", qm -> fID);
  printMesh(qm, buffer);
#endif
  return EG_makeValidMesh(qm, 0, NULL, 1);
}


/********************   IO FUNCTIONS **********************/

static void printQuadSpecs(meshMap *qm, int id)
{
  int i, v = 0 , val = 0;

  --id;
  printf(" QUAD %d HAS VERTICES ", id+1);
  for (i = 0; i < 4; ++i) {
      v   = qm -> qIdx[4 * id + i];
      if (v <= 0) printf( " DEL %d %d \n ", v, 0);
      else {
	  val = qm -> valence[v -1][1];
	  printf(" %d ( val = %d )  ", v, val);
      }
  }
  printf("\t AND ADJACENT QUADS ");
  for (i = 0; i < 4; ++i) printf(" %d ", qm -> qAdj[4*id + i]);
  printf("\n");
}


//#ifdef DEBUG

static void printStarFile( meshMap *qm, vStar *star, char *name)
{
  int   i, v0, v, j, stat, p, id;
  vStar *st = NULL;
  FILE  *fil= NULL;
  
  fil = fopen (name, "w");
  if (fil == NULL) return;
  printf(" Writing in %s\n ", name);
  for (i = 0; i < star -> nQ; i++) {
      if (star -> quads[i] == -1) continue;
      id = EG_quadVertIdx(qm, star -> quads[i], star -> verts[0]);
      for (p = 0; p <= 4; p++)  {
	  v = qm -> qIdx[4 * (star -> quads[i] - 1) + (id + p)%4] - 1;
	  fprintf(fil,"%lf %lf %lf %d\n ",
                  qm -> xyzs[3 * v    ], qm -> xyzs[3 * v + 1],
                  qm -> xyzs[3 * v + 2], v + 1);
      }
      fprintf(fil, "\n\n\n ");
  }
  fprintf(fil, "\n\n\n ");
  for (j = 1; j < star -> nV; j++) {
      v0 = star -> verts[j];
      if (v0 < 0) continue;
      stat = EG_buildStar(qm, &st, v0);
      if (stat != EGADS_SUCCESS || st == NULL) continue;
      for (i = 0; i < st -> nQ; i++) {
	  if (st -> quads[i] == -1) continue;
	  id = EG_quadVertIdx(qm, st -> quads[i], v0);
	  for (p = 0; p <= 4; p++)  {
	      v = qm -> qIdx[4 * (st -> quads[i] - 1) + (id + p)%4] - 1;
	      if (v < 0)  {
		  printf(" VERTEX IS %d !!!! \n ", v);
		  continue;
	      }
	      fprintf(fil,"%lf %lf %lf %d\n ", qm-> xyzs[3 * v    ],
                      qm-> xyzs[3 * v + 1], qm-> xyzs[3 * v + 2], v + 1);
	  }
	  fprintf(fil, "\n\n\n ");
      }
      fprintf(fil, "\n\n\n ");
  }
  EG_freeStar(&st);
  fclose(fil);
  return;
}
//#endif


#ifdef IOFILES
static void printMesh(meshMap *qm, /*@null@*/ char *name)
{
  int    i, j, k, v, d;
  double eval[18], average[2], dist;
  char   name2[500];
  FILE   *fout = NULL;

  if (name == NULL) {
      snprintf (name2, 500, "FACE_%d_M_%d",qm -> fID, qm -> plotcount++);
      printf(" Writing in File %s  \n", name2);
      fout = fopen(name2, "w");
  } else {
      printf(" Writing in File %s  \n", name);
      fout = fopen(name, "w");
  }
  if (fout == NULL) return;

  for (i = 0; i < qm -> totQ; ++i) {
      if (qm -> qIdx[4 * i] == -2 ||
	  qm -> vType[qm -> qIdx[4 * i] - 1] == -2) continue;
      for (k = 0; k < 4; ++k) {
	  v = qm -> qIdx[4*i + k] - 1;
	  fprintf(fout, "%lf %lf %lf %d %lf %lf \n", qm -> xyzs[3*v  ],
		  qm -> xyzs[3*v+1], qm -> xyzs[3*v+2], v + 1,
		  qm ->uvs[2*v  ], qm ->uvs[2*v+1]);
	  if ( qm -> vType[v] != -1) continue;
	  dist = 0.0;
	  EG_evaluate(qm -> face, &qm -> uvs[2*v], eval);
	  for (d = 0; d < 3; ++d)
	    dist += (eval[d] - qm -> xyzs[3*v + d]) *
	            (eval[d] - qm -> xyzs[3*v + d]);
	  dist = sqrt(dist);
	  if (dist > EPS11) {
	      printf(" DIST = %11.2e  IN QUAD %d  VERTEX %d. UVs and xyzs are mismatched.  UV %lf  %lf \n",
		     dist,i+1, v+1, qm -> uvs[2*v], qm -> uvs[2*v + 1]);
	      for (d = 0; d < 3; ++d)
		printf( "%lf  != %lf \t", eval[d], qm -> xyzs[3*v + d]);
	      fclose(fout);
	      return;
	  }
      }
      v = qm -> qIdx[4*i] - 1;
      fprintf(fout, "%lf %lf %lf %d %lf %lf\n",  qm -> xyzs[3*v  ],
	      qm -> xyzs[3*v+1], qm -> xyzs[3*v+2], v + 1,
	      qm ->uvs[2*v  ], qm ->uvs[2*v+1]);
      fprintf(fout,"\n\n");
      average[0] = 0.0; average[1] = 0.0;
      for (j = 0; j < 4; j++) {
	  k = qm -> qIdx[4 * i + j] - 1;
	  average[0] += 0.25 * qm -> uvs[2 * k    ];
	  average[1] += 0.25 * qm -> uvs[2 * k + 1];
      }
      EG_evaluate(qm -> face, average, eval);
      fprintf(fout, "%lf %lf %lf %d %lf %lf\n", eval[0], eval[1], eval[2],
	      i + 1, average[0], average[1]);
      fprintf(fout,"\n\n");
  }
  fclose(fout);
}


static int EG_wvsData(meshMap *qm, /*@null@*/ char *buffer)
{
  int  i;
  FILE *fil = NULL;
  char aux[100];

  if (buffer == NULL) {
      snprintf(aux, 100, "wvsFace_%d_M_%d", qm -> fID, qm -> plotcount++);
      fil = fopen(aux, "w");
  } else {
      fil = fopen (buffer, "w");
      printf(" Writing in File %s  \n", buffer);
  }
  if (fil == NULL) return EGADS_MALLOC;
  fprintf(fil, "%d %d\n", qm -> totV, qm -> totQ);
  for (i = 0; i < qm -> totV; i++) {
      fprintf(fil, "%lf %lf %lf \n",
	      qm -> xyzs[3 * i    ],
	      qm -> xyzs[3 * i + 1],
	      qm -> xyzs[3 * i + 2]);
  }
  fprintf(fil,"\n");
  for (i = 0; i < qm -> totQ; i++) {

      fprintf(fil, "%d %d %d %d\n",
	      qm -> qIdx[4 * i    ],
	      qm -> qIdx[4 * i + 1],
	      qm -> qIdx[4 * i + 2],
	      qm -> qIdx[4 * i + 3]);
  }
  fclose (fil);
  return EGADS_SUCCESS;
}
#endif


static void printQuadGroup(meshMap *qm, quadGroup qg)
{
  int i, stat;

  for (i = 0; i < 2; i++) {
      stat = checkQuad(qm, qg.q[i]);
      if ( stat != EGADS_SUCCESS) {
	  printf(" Quad %d is %d !!! \n", qg.q[i], stat);
	  return;
      }
      printQuadSpecs(qm, qg.q[i]);
  }
  for (i = 0; i < 6; i++) {
      stat = checkVertex (qm, qg.verts[i]);
      if (stat != EGADS_SUCCESS) {
	  printf(" Vertex %d is %d !!! \n", qg.verts[i], stat);
	  return;
      }
      printf(" QV(%d) = %d type %d valence %d \n",
	     i, qg.verts[i], qm -> vType[qg.verts[i] - 1], qg.vals[i]);
  }
}


int EG_makeQuadTess(bodyQuad bodydata, ego *quadTess)
{
  int          i, j, npts, nt, stat, outLevel, nedges, mQ, *tris;
  const int    *ptype, *pindex, *trs, *trc;
  const double *xyzs, *ts, *uvs;
  ego          tess, newTess, obj, *edges;
  egTessel     *btess;
  
  *quadTess = NULL;
  tess      = bodydata.tess;
  btess     = (egTessel *) tess->blind;
  outLevel  = EG_outLevel(tess);
  
  /* initialize the new tessellation object */
  stat = EG_initTessBody(btess->src, &newTess);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_initTessBody = %d (EG_makeQuadTess)!\n", stat);
    return stat;
  }
  
  stat = EG_getBodyTopos(btess->src, NULL, EDGE, &nedges, &edges);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getBodyTopos E = %d (EG_makeQuadTess)!\n", stat);
    EG_deleteObject(newTess);
    return stat;
  }
  
  /* copy the Edges */
  for (i = 0; i < nedges; i++) {
    if (edges[i]->mtype == DEGENERATE) continue;
    stat = EG_getTessEdge(tess, i+1, &npts, &xyzs, &ts);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_getTessEdge %d = %d (EG_makeQuadTess)!\n",
               i+1, stat);
      EG_free(edges);
      EG_deleteObject(newTess);
      return stat;
    }
    if (npts == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_getTessEdge %d -- no points (EG_makeQuadTess)!\n",
               i+1);
      EG_free(edges);
      EG_deleteObject(newTess);
      return EGADS_INDEXERR;
    }
    
    stat = EG_setTessEdge(newTess, i+1, npts, xyzs, ts);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_setTessEdge %d = %d (EG_makeQuadTess)!\n",
               i+1, stat);
      EG_free(edges);
      EG_deleteObject(newTess);
      return stat;
    }
  }
  EG_free(edges);
  
  /* make the Face tessellations */
  
  for (mQ = i = 0; i < bodydata.nfaces; i++) {
    if (bodydata.qm[i]->totQ > mQ) mQ = bodydata.qm[i]->totQ;
    if (btess->tess2d[i].tfi != 1) continue;
    stat = EG_getTessFace(tess, i+1, &npts, &xyzs, &uvs, &ptype, &pindex,
                          &nt, &trs, &trc);
    if ((stat != EGADS_SUCCESS) || (nt == 0)) {
      if (outLevel > 0)
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: EG_getTessFace %d = %d (EG_quadTess)!\n",
                 i+1, stat);
        } else {
          printf(" EGADS Error: Face %d has no tessellation (EG_quadTess)!\n",
                 i+1);
        }
      EG_deleteObject(newTess);
      return stat;
    }
    
    stat = EG_setTessFace(newTess, i+1, npts, xyzs, uvs, nt, trs);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_setTessFace %d = %d (EG_makeQuadTess)!\n",
               i+1, stat);
      EG_deleteObject(newTess);
      return stat;
    }
  }
#ifdef DEBUG
  printf(" Max Quads = %d\n", mQ);
#endif
  tris = (int *) EG_alloc(6*mQ*sizeof(int));
  if (tris == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Quad Allocation = %d (EG_makeQuadTess)!\n", mQ);
    EG_deleteObject(newTess);
    return EGADS_MALLOC;
  }
  
  /* fill in quads as triangle pairs */
  for (i = 0; i < bodydata.nfaces; i++) {
    if (btess->tess2d[i].tfi == 1) continue;
    for (j = 0; j < bodydata.qm[i]->totQ; j++) {
      tris[6*j  ] = bodydata.qm[i]->qIdx[4*j  ];
      tris[6*j+1] = bodydata.qm[i]->qIdx[4*j+1];
      tris[6*j+2] = bodydata.qm[i]->qIdx[4*j+2];
      tris[6*j+3] = bodydata.qm[i]->qIdx[4*j  ];
      tris[6*j+4] = bodydata.qm[i]->qIdx[4*j+2];
      tris[6*j+5] = bodydata.qm[i]->qIdx[4*j+3];
    }
    stat = EG_setTessFace(newTess, i+1, bodydata.qm[i]->totV,
                          bodydata.qm[i]->xyzs,
                          bodydata.qm[i]->uvs,
                          bodydata.qm[i]->totQ*2, tris);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_setTessFace %d = %d (EG_makeQuadTess)!\n",
               i+1, stat);
      EG_free(tris);
      EG_deleteObject(newTess);
      return stat;
    }
  }
  EG_free(tris);
  
  /* close up the open tessellation */
  stat = EG_statusTessBody(newTess, &obj, &i, &npts);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_statusTessBody = %d (EG_makeQuadTess)!\n", stat);
    EG_deleteObject(newTess);
    return stat;
  }
  if (i != 1) {
    if (outLevel > 0)
      printf(" EGADS Error: New Tessellation Object is Open (EG_makeQuadTess)!\n");
    EG_deleteObject(newTess);
    return EGADS_TESSTATE;
  }
#ifndef LITE
  stat = EG_attributeAdd(newTess, ".tessType", ATTRSTRING, 4,
                         NULL, NULL, "Quad");
  if (stat != EGADS_SUCCESS)
    if (outLevel > 0)
      printf(" EGADS Warning: EG_attributeAdd = %d (EG_makeQuadTess)!\n", stat);
#endif
  
  *quadTess = newTess;
  return EGADS_SUCCESS;
}


#ifdef STANDALONE
int main(int argc, char *argv[])
{
  clock_t      start_t, end_t, total_t;
  int          stat = 0,  f , i, j, iBody, oclass, mtype, nbody;
  int          atype, alen, *senses, min, uvt, FACECHOICE = -1;
  const int    *ints;
  float        arg;
  double       box[6], size, params[3], time, frac;
  const double *reals;
  const char   *OCCrev, *string;
#ifdef IOFILES
  char buffer[100];
#endif
  ego          context, tess, model, geom, *bodies, *dum;
  bodyQuad     *bodydata;

  start_t = clock();
  if (argc < 6) {
      printf("\n Usage: %s = (1) filename, Tesselation:  (2) angle (3) maxlen and chord (4) laplacian  type: 0 = average 1 = arc length 2 = invEval \n\n", argv[0]);
      printf(" You have entered : \n");
      for (i = 0; i < argc; i++) printf(" argv[%d] = %s\n ", i, argv[i]);
      return 1;
  }
  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %2d.%02d with %s\n\n", i, j, OCCrev);

  /* initialize */
  stat = EG_open(&context);
  if (stat != EGADS_SUCCESS) return 1;
  stat = EG_loadModel(context, 0, argv[1], &model);
  printf(" EG_loadModel      = %d\n", stat);
  printf(" EG_getBoundingBox = %d\n", EG_getBoundingBox(model, box));
  printf("       BoundingBox = %lf %lf %lf\n", box[0], box[1], box[2]);
  printf("                     %lf %lf %lf\n", box[3], box[4], box[5]);
  printf(" \n");
                            size = box[3]-box[0];
  if (size < box[4]-box[1]) size = box[4]-box[1];
  if (size < box[5]-box[2]) size = box[5]-box[2];

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
  uvt        = atoi(argv[5]);
  if (uvt == AVERAGE)
    printf("******** USING UV AVERAGES ************\n");
  else if (uvt== ARCLENGTH)
    printf("******** USING UV ARCLENGTH WEIGHTS ************\n");
  else if (uvt == INVEVAL)
    printf("******** USING POLYGON CENTROID ************\n");
  else {
      printf("******** YOU HAVE ENTERED THE WRONG CONDITION FOR LAPLACIAN %d !!! ************\n", uvt);
      return EGADS_INDEXERR;
  }
  if (argc == 7) sscanf(argv[6], "%d", &FACECHOICE);
  
  /* get all bodies */
  stat = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody,
			&bodies, &senses);
  if (stat != EGADS_SUCCESS) {
      printf(" EG_getTopology = %d\n", stat);
      return 1;
  }
  printf(" EG_getTopology:     nBodies = %d\n", nbody);
  bodydata = (bodyQuad *) malloc(nbody*sizeof(bodyQuad));
  if (bodydata == NULL) {
      printf(" MALLOC Error on Body storage!\n");
      return 1;
  }
 
  /* fill our structure a body at at time */
  for (iBody = 0; iBody < nbody; iBody++) {
      stat = EG_attributeAdd(bodies[iBody], ".qParams",
			     ATTRSTRING, 4, NULL, NULL, "off");
      if (stat != EGADS_SUCCESS)
	printf(" Body %d: attributeAdd = %d\n", iBody, stat);
      EG_getTopology(bodies[iBody], &geom, &oclass,
		     &mtype, NULL, &j, &dum, &senses);
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
			     &bodydata[iBody].nedges, NULL);
      if ((stat != EGADS_SUCCESS) || (i != EGADS_SUCCESS)) {
	  printf(" EG_getBodyTopos Face = %d\n", stat);
	  printf(" EG_getBodyTopos Edge = %d\n", i);
	  continue;
      }
      stat = EG_makeTessBody(bodies[iBody], params, &bodydata[iBody].tess);
      if (stat != EGADS_SUCCESS) {
	  printf(" EG_makeTessBody %d = %d\n", iBody, stat);
	  continue;
      }
      tess = bodydata[iBody].tess;
      /* disable regularization in EGADS */
      stat = EG_attributeAdd(tess, ".qRegular", ATTRSTRING, 3, NULL, NULL, "Off");
      if (stat != EGADS_SUCCESS)
          printf(" EG_attributeAdd qRegular %d = %d\n", iBody, stat);
      stat = EG_quadTess(tess, &bodydata[iBody].tess);
      if (stat != EGADS_SUCCESS) {
	  printf(" EG_quadTess %d = %d  -- reverting...\n", iBody, stat);
	  bodydata[iBody].tess = tess;
	  continue;
      }
      EG_deleteObject(tess);
  }
  for (iBody = 0; iBody < nbody; iBody++) {
      stat = EG_attributeRet(bodydata[iBody].tess, ".tessType", &atype,
			     &alen, &ints, &reals, &string);
      if (stat == EGADS_SUCCESS)
	stat = EG_createMeshMap(&bodydata[iBody], uvt);
      if (stat != EGADS_SUCCESS) goto cleanup;
      for (iBody = 0; iBody < nbody; iBody++) {
	  for (f = 0; f < bodydata[iBody].nfaces; ++f) {
	      if (FACECHOICE >= 0) f = FACECHOICE;
#ifdef IOFILES
	      printMeshStats(bodydata[iBody].qm[f], 0);
	      sprintf(buffer,"gnuInit_%i.txt",f+1);
	      printMesh(bodydata->qm[f] , buffer);
	      sprintf(buffer,"wvsInit_%i.txt",f+1);
	      EG_wvsData(bodydata->qm[f] , buffer);
#endif
#ifdef DEBUG
	      printf("===================================================================\n");
	      printf("===================================================================\n");
	      printf(" FULL MESH REGULARIZATION FACE %d \n ", f + 1);
	      printf("===================================================================\n");
	      printf("===================================================================\n");
#endif
	      stat      = EG_fullMeshRegularization(bodydata[iBody].qm[f]);
	      if (stat != EGADS_SUCCESS) {
		  fprintf(stderr, " EG_fullMeshRegularization face %d / %d = %d !!!\n ", f + 1, bodydata[iBody].nfaces,  stat);
		  break;
	      }
#ifdef IOFILES
	      printMeshStats(bodydata[iBody].qm[f], 1);
	      snprintf(buffer,100,"gnuFinal_%i.txt", f+1);
	      printMesh(bodydata->qm[f] , buffer);
	      snprintf(buffer,100,"wvsFinal_%i.txt",f+1);
	      EG_wvsData(bodydata->qm[f], buffer);
#endif
	      if (FACECHOICE >= 0) break;
	  }
          if (FACECHOICE == -1) {
              ego newTess;
            
              stat = EG_makeQuadTess(bodydata[iBody], &newTess);
              if (stat != EGADS_SUCCESS) {
                  fprintf(stderr, " EG_makeQuadTess = %d!\n ", stat);
                  goto cleanup;
              }
              stat = EG_deleteObject(newTess);
              if (stat != EGADS_SUCCESS) {
                  fprintf(stderr, " EG_deleteObject = %d!\n ", stat);
                  goto cleanup;
              }
          }
      }
  }

cleanup:
  printf(" STAT IN CLEANUP  %d\n", stat);
  for (iBody = 0; iBody < nbody; iBody++) {
      EG_destroymeshMap(&bodydata[iBody]);
      EG_free(bodydata[iBody].faces);
      EG_deleteObject(bodydata[iBody].tess);
  }
  EG_free(bodydata);

  end_t   = clock();
  total_t = end_t - start_t;
  time    = (double) total_t / CLOCKS_PER_SEC;
  min     = floor(time) / 60;
  frac    = time -  min * 60;
  fprintf(stderr, "Total time taken by CPU: %d minutes and %f seconds\n",
          min, frac);
  printf("Total time taken by CPU: %d minutes and %f seconds\n",
	  min, frac);

  EG_deleteObject(model);
  EG_close(context);
  return 0;
}
#endif
