/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Quad Regularization Functions
 *
 *      Copyright 2018-2019, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#ifdef STANDALONE
#include <time.h>
#define REPORT
#endif

#include "egads.h"
#include "xdQuads.h"

#define DEBUG
#ifdef DEBUG
#define REPORT
#endif


/* BASE-LEVEL FUNCTION */
static int inList(int n, int *list, int p)
{
  int i;

  if (n == 0) return -1;
  for (i = 0; i < n; i++)
    if (list[i] == p) return i;
  return -1;
}


/* IO FUNCTIONS */
static void meshCount(meshMap *qm, int *nI, int *nV, int *nQ)
{
  int i, qSum, vSum = 0, vSum2 = 0;

  for (i = qSum = 0; i < qm->totQ; i++)
    if (qm->qIdx[4 * i] != -2) qSum++;
  for (i = vSum = 0; i < qm->totV; i++) {
      if (qm->vType[i] == -2) continue;
      vSum2++;
      if (qm->valence[i][1] != 4) vSum++;
  }
  *nQ = qSum;
  *nI = vSum;
  *nV = vSum2;
}


#ifdef REPORT
static void meshStats(meshMap *qm, int sweep)
{
  int  i,len, val, intVal[100], boundVal[100];
  char buffer[100];
  FILE *fout;

  if (sweep == 0)
    sprintf(buffer, "OriginalMeshStats_face_%d_.txt", qm->fID);
  else if (sweep == 10)
    sprintf(buffer, "PreprocessMeshStats_face_%d_.txt", qm->fID);
  else
    sprintf(buffer, "RegularlMeshStats_face_%d_.txt", qm->fID);
  fout = fopen(buffer, "w");
  if (fout == NULL) return;
  len = qm->totV;
  for (i = 0; i < 100; ++i) {
      intVal  [i] = 0;
      boundVal[i] = 0;
  }
  for (i = 0; i < len; ++i) {
      if (qm->vType[i] !=  -2) {
	  if (qm->vType[i] == -1) {
	      val = qm->valence[i][1];
	      ++intVal[val];
	  } else {
	      ++boundVal[qm->valence[i][1]];
	  }
      }
  }
  fprintf(fout,"---- TOTAL VERTICES %d TOTAL QUADS %d --------------\n",
	  qm->totV - qm->remV[0], qm->totQ - qm->remQ[0]);
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


static void gnuData(meshMap *qm, /*@null@*/ char *name)
{
  int    i, j, k, v, d;
  double eval[18], average[2], dist;
  char   name2[100];
  FILE   *fout = NULL;

  if (name == NULL) {
      snprintf(name2, 100, "FACE_%d_M_%d",qm->fID, qm->plotcount++);
      printf(" Writing in File %s  \n", name2);
      fout = fopen(name2, "w");
  } else {
      printf(" Writing in File %s  \n", name);
      fout = fopen(name, "w");
  }
  if (fout == NULL) return;

  for (i = 0; i < qm->totQ; ++i) {
      if (qm->qIdx[4 * i] == -2 ||
	  qm->vType[qm->qIdx[4 * i] - 1] == -2) continue;
      for (k = 0; k < 4; ++k) {
	  v  = qm->qIdx[4*i + k] - 1;
	  fprintf(fout, "%lf %lf %lf %d %lf %lf \n", qm->xyzs[3*v  ],
		  qm->xyzs[3*v+1], qm->xyzs[3*v+2], v + 1,
		  qm->uvs [2*v  ], qm->uvs [2*v+1]);
	  if (qm->vType[v] != -1) continue;
	  dist = 0.0;
	  EG_evaluate(qm->face, &qm->uvs[2*v], eval);
	  for (d = 0; d < 3; ++d)
	    dist += (eval[d] - qm->xyzs[3*v + d]) *
	    (eval[d] - qm->xyzs[3*v + d]);
	  dist = sqrt(dist);
	  if (dist > EPS08) {
	      printf(" DIST = %11.2e  IN QUAD %d  VERTEX %d. UVs and xyzs are mismatched.  UV %lf  %lf \n",
		     dist,i+1, v+1, qm->uvs[2*v], qm->uvs[2*v + 1]);
	      for (d = 0; d < 3; ++d)
		printf(" %lf  != %lf \t", eval[0], qm->xyzs[3*v + d]);
	      fclose(fout);
	      return;
	  }
      }
      v = qm->qIdx[4*i] - 1;
      fprintf(fout, "%lf %lf %lf %d %lf %lf\n",  qm->xyzs[3*v  ],
	      qm->xyzs[3*v + 1], qm->xyzs[3*v + 2], v + 1,
	      qm->uvs [2*v    ], qm->uvs [2*v + 1]);
      fprintf(fout,"\n\n");
      average[0] = 0.0; average[1] = 0.0;
      for (j = 0; j < 4; j++) {
	  k = qm->qIdx[4 * i + j] - 1;
	  average[0] += 0.25 * qm->uvs[2 * k    ];
	  average[1] += 0.25 * qm->uvs[2 * k + 1];
      }
      EG_evaluate(qm->face, average, eval);
      fprintf(fout, "%lf %lf %lf %d %lf %lf\n", eval[0], eval[1], eval[2],
	      i + 1, average[0], average[1]);
      fprintf(fout,"\n\n");
  }
  fclose(fout);
}


static void wvsData(meshMap *qm, /*@null@*/ char *buffer)
{
  int  i;
  FILE *fil = NULL;
  char aux[100];

  if (buffer == NULL) {
      snprintf(aux, 100, "wvsFace_%d_M_%d", qm->fID, qm->plotcount++);
      fil = fopen(aux, "w");
  } else {
      fil = fopen (buffer, "w");
      printf(" Writing in File %s  \n", buffer);
  }
  if (fil == NULL) {
      printf(" I couldn't write in wvsData file %s!!!!\n ", buffer);
      return ;
  }
  fprintf(fil, "%d %d\n", qm->totV, qm->totQ);
  for (i = 0; i < qm->totV; i++) {
      fprintf(fil, "%lf %lf %lf \n", qm->xyzs[3 * i    ],
	      qm->xyzs[3 * i + 1], qm->xyzs[3 * i + 2]);
  }
  fprintf(fil,"\n");
  for (i = 0; i < qm->totQ; i++) {
      fprintf(fil, "%d %d %d %d\n",  qm->qIdx[4 * i    ],
	      qm->qIdx[4 * i + 1], qm->qIdx[4 * i + 2],
	      qm->qIdx[4 * i + 3]);
  }
  fclose (fil);
  return ;
}
#endif


/* DEBUG FUNCTIONS */
#ifdef DEBUG
static void printVertex(meshMap *qm, int v)
{
  v--;
  printf(" #vertex %d  =================================== \n" , v+1);
  printf("%lf %lf %lf %lf %lf %d\n",
	 qm->xyzs[3 * v    ], qm->xyzs[3 * v + 1],
	 qm->xyzs[3 * v + 2], qm->uvs [2 * v + 1],
	 qm->uvs [2 * v    ], v + 1);
  printf(" #=============================================\n");
}


static void printQuad (meshMap *qm, int id)
{
  int    i, v = 0 , val = 0;
  double uv[3], uva[2], xyz[18];

  uv[0] = uv[1] = uv[2] = 0.0;
  --id;
  printf("---------------------------------------\n");
  printf(" QUAD %d HAS VERTICES ", id+1);
  for (i = 0; i < 4; ++i) {
      v   = qm->qIdx[4 * id + i];
      if (v <= 0) printf(" DEL %d %d \n ", v, 0);
      else {
	  val = qm->valence[v - 1][1];
	  printf(" %d ( val = %d )  ", v, val);
      }
  }
  printf("\t AND ADJACENT QUADS ");
  for (i = 0; i < 4; ++i) printf(" %d ", qm->qAdj[4*id + i]);
  printf("\n\n Coordinates \n");
  for (i = 0; i <= 4; i++) {
      v = qm->qIdx[4 * id + i%4] - 1;
      printf("%lf %lf %lf %d %lf %lf \n",
	     qm->xyzs[3 * v    ], qm->xyzs[3 * v + 1],
	     qm->xyzs[3 * v + 2], v + 1,
	     qm->uvs [2 * v + 1], qm->uvs [2 * v    ]);
      if (i < 4) {
	  uv[0] += 0.25 * qm->xyzs[3 * v    ];
	  uv[1] += 0.25 * qm->xyzs[3 * v + 1];
	  uv[2] += 0.25 * qm->xyzs[3 * v + 2];
      }
  }
  i = EG_invEvaluate(qm->face, uv, uva, xyz);
  printf("\n\n\n %lf %lf %lf %d %lf %lf \n",
	 xyz[0], xyz[1], xyz[2], id + 1, uv[0], uv[1]);
  printf(" \n\n---------------------------------------\n");
}


static void printQuadGroup(meshMap *qm, quadGroup qg)
{
  int i;

  for (i = 0; i < 2; i++) {
      if (qg.q[i] < 0 || qg.q[i] > qm->totQ ) {
	  printf(" Quad %d is out of bounds !!! \n", qg.q[i] );
	  return;
      }
      printQuad(qm, qg.q[i]);
  }
  for (i = 0; i < 6; i++) {
      printf(" QV (%d) = %d ", i, qg.verts[i] );
      if (qg.verts[i] < 0) printf(" NO VAL !! NEGATIVE! \n");
      if (qg.verts[i] ==0 || qg.verts[i] > qm->totV) {
	  printf(" QUAD HAS VERTEX %d !!\n ", qg.verts[i] );
	  exit(1);
      }
      else printf("type %d valence %d \n", qm->vType[qg.verts[i] - 1],
		  qg.vals[i]);
  }
}


static int checkMesh(meshMap *qm)
{
  int i, j, k, quad, val, v1, v2;

  for (i = 0; i < qm->totV; i++) {
      if (qm->vType[i] == -2) continue;
      val = qm->valence[i][2];
      if (val <= 0 || val > qm->totV) {
	  printf(" checkMesh vertex %d valence out of bounds %d (max %d ) !!\n",
		 i + 1, qm->valence[i][2], qm->totV);
	  return EGADS_INDEXERR;
      }
      quad = qm->valence[i][0];
      if (quad <= 0 || quad > qm->totQ) {
	  printf("In checkMesh quad for %d is out of bounds!!! %d > %d \n ",
		 i + 1, quad, qm->totQ);
	  return EGADS_INDEXERR;
      }
      for (j = 0; j < val; ++j) {
	  v2 = qm->valence[i][3 + j] - 1;
	  for (k = 0; k < qm->valence[v2][2]; k++)
	    if (qm->valence[v2][3 + k] == i + 1) {
		k = -1;
		break;
	    }
	  if (k != -1) {
	      printf(" checkMesh Vertex %d has assigned %d as link but %d doesn't point at %d\n",
		     i + 1, v2 + 1, v2 + 1,  i + 1);
	      for (k = 0; k < qm->valence[i][2]; k++)
		printf(" V1 %d ==> %d \n ", k, qm->valence[i][3+k]);
	      printf("\n\n");
	      for (k = 0; k < qm->valence[v2][2]; k++)
		printf(" V1 %d ==> %d \n ", k, qm->valence[v2][3+k]);
	      return EGADS_INDEXERR;
	  }
      }
  }
  for (i = 0; i < qm->totQ; i++) {
      if (qm->qIdx[4 * i] == -2) continue;
      for (j = 0; j < 4; ++j) {
	  v1 = qm->qAdj[4 * i  + j] - 1;
	  if (v1 < 0) continue;
	  val = -1;
	  for (k = 0; k < 4; ++k) {
	      if (qm->qAdj[4 * v1 + k] == i + 1) {
		  val = 1;
		  break;
	      }
	  }
	  if (val == -1) {
	      fprintf (stderr, " checkMesh quads %d and %d don't point at each other\n",
		       i +1, v1 + 1);
	      printQuad(qm, i + 1);
	      printQuad(qm,v1 + 1);
	      return EGADS_INDEXERR;
	  }
      }
  }
  return EGADS_SUCCESS;
}
#endif


static void updateVertex(meshMap *qm, int vID, double *uv)
{
  double eval[18];

  if (EG_evaluate(qm->face, uv, eval) != EGADS_SUCCESS ||
      qm->vType[vID -1]   != -1) return;
  qm->uvs [2 * (vID - 1)    ] = uv  [0];
  qm->uvs [2 * (vID - 1) + 1] = uv  [1];
  qm->xyzs[3 * (vID - 1)    ] = eval[0];
  qm->xyzs[3 * (vID - 1) + 1] = eval[1];
  qm->xyzs[3 * (vID - 1) + 2] = eval[2];
}

static int EG_adjQtoPair(meshMap *qm, int qID, int v1, int v2, int *adj)
{
  int i, aux = -1;

  adj[0] = -1; adj[1] = -1;
  for (i = 0; i < 4; ++i) {
      if (qm->qIdx[4*(qID - 1) + i] == v1) adj[0] = i;
      if (qm->qIdx[4*(qID - 1) + i] == v2) aux    = i;
      if (aux != -1 && adj[0] != -1) break;
  }
  if (aux == -1 || adj[0] == -1) return EGADS_SUCCESS;
  if      (abs(adj[0] - aux) == 3) adj[0] = 3;
  else if (aux < adj[0]          ) adj[0] = aux;
  adj[1] = qm->qAdj[4*(qID - 1) + adj[0]];

  return EGADS_SUCCESS;
}

static int EG_quadVertIdx(meshMap *qm, int q, int v)
{
  int i = 0;

  for (i = 0; i < 4; i++)
    if (qm->qIdx[4 * (q - 1) + i] == v) return i;

  return -1;
}


static void EG_freeStar(vStar **star)
{
  if ((*star) == NULL) return;
  EG_free((*star)->verts);
  EG_free((*star)->quads);
  EG_free((*star)->idxQ);
  EG_free((*star)->idxV);
  EG_free(*star);
  *star = NULL;
}


static int EG_buildStar(meshMap *qm, vStar **star, int vID)
{
  int i = 0, id0 = -1, q = 0, auxV, auxQ, v = 0, quadID, prevQuad;
  int it = 0, it2 = 0, adj[2], *vertex = NULL, *quads = NULL;
  int qLoop[8] = {0, 1, 2, 3, 0, 1, 2, 3};

  vertex = (int *) EG_alloc(qm->totV * sizeof(int));
  quads  = (int *) EG_alloc(qm->totQ * sizeof(int));
  if (vertex == NULL || quads == NULL) {
      printf(" EG_buildStar MALLOC at quads & verts!!\n ");
      if (vertex != NULL) EG_free(vertex);
      if (quads  != NULL) EG_free(quads);
      return EGADS_MALLOC;
  }
  // quads are -1 bias
  quadID      = qm->valence[vID - 1][0] - 1;
  vertex[v++] = vID;
  it          = 0;
  do {
      id0 = EG_quadVertIdx(qm, quadID + 1, vID);
      if (id0 < 0) {
	  printf(" In EG_buildStar id for Vert %d in Quad %d is %d !!\n",
		 vID, quadID + 1, id0);
#ifdef DEBUG
	  printQuad(qm, quadID + 1);
#endif
	  EG_free(vertex);
	  EG_free(quads);
	  return EGADS_INDEXERR;
      }
      for (i = 1; i <= 2; ++i)
	vertex[v++] = qm->qIdx[4 * quadID + qLoop[id0 + i]];
      quads[q++] = quadID + 1;
      prevQuad   = quadID;
      quadID     = qm->qAdj[4 * prevQuad + qLoop[id0 + 3 ]] - 1;
      if (quadID  < 0) { //make a "ghost" quad
	  auxQ        = prevQuad;
	  vertex[v++] = qm->qIdx[4*auxQ + qLoop[id0 + 3]];
	  auxV        = qm->qIdx[4*auxQ + qLoop[id0 + 1]];
	  it2         = 0;
	  do {
	      i    = EG_adjQtoPair(qm, auxQ + 1, vID, auxV, adj);
	      if (adj[1] == -1 || i != EGADS_SUCCESS) break;
	      auxQ = adj[1] - 1;
	      i    = EG_quadVertIdx(qm, auxQ + 1, vID);
	      if (i < 0) {
		  printf(" In buildStar vertex Id %d in quad %d is %d\n",
			 vID, quadID + 1, i);
#ifdef DEBUG
		  printQuad(qm, quadID + 1);
#endif
		  EG_free(vertex);
		  EG_free(quads);
		  return EGADS_INDEXERR;
	      }
	      auxV = qm->qIdx[4 * auxQ + qLoop[i + 1]];
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
  (*star)->nQ    = q;
  (*star)->nV    = v;
  (*star)->verts = (int *) EG_alloc(    v * sizeof(int));
  (*star)->quads = (int *) EG_alloc(    q * sizeof(int));
  (*star)->idxV  = (int *) EG_alloc(2 * v * sizeof(int));
  (*star)->idxQ  = (int *) EG_alloc(2 * q * sizeof(int));
  if ((*star)->verts == NULL || (*star)->quads == NULL ||
      (*star)->idxV  == NULL || (*star)->idxQ  == NULL) {
      if ((*star)->verts != NULL) EG_free((*star)->verts);
      if ((*star)->quads != NULL) EG_free((*star)->quads);
      if ((*star)->idxV  != NULL) EG_free((*star)->idxV);
      if ((*star)->idxQ  != NULL) EG_free((*star)->idxQ);
      EG_free((*star));
      EG_free(vertex);
      EG_free( quads);
      return EGADS_MALLOC;
  }
  for (i = 0; i < q; ++i) {
      (*star)->quads[i    ] = quads[i];
      (*star)->idxQ [i    ] = i;
      (*star)->idxQ [q + i] = i;
  }
  for (i = 0; i < v; ++i) {
      (*star)->verts[i    ] = vertex[i];
      (*star)->idxV [i    ] = i;
      (*star)->idxV [v + i] = i + 1;
  }
  EG_free (vertex);
  EG_free (quads);
  return EGADS_SUCCESS;
}


static int setValence(meshMap *qm, int vID)
{
  int   i, i4[10], i3[10], i1[100], j, si;

  if (vID <= 0 || vID > qm->totV) {
      printf(" VERTEX %d is out of bounds !\n ", vID);
      return EGADS_INDEXERR;
  }
  si = vID - 1;
  i  = EG_buildStar(qm, &qm->star[si], vID);
  if (i != EGADS_SUCCESS || qm->star[si] == NULL) {
      printf(" In setValence for vertex %d stat from EG_buildStar = %d\n", vID, i);
      return i;
  }
  if (qm->valence[si][2] < qm->star[si]->nQ ) {
      qm->valence[si] = EG_reall (qm->valence[si], (3 + qm->star[si]->nQ ) * sizeof(int) );
      if (qm->valence[si] == NULL ) {
	  printf(" VALENCE is NULL %d \n ", si + 1 );
	  return EGADS_MALLOC;
      }
  }
  qm->valence[si][2] = qm->star[si]->nQ; // real valence
  qm->valence[si][1] = qm->star[si]->nQ;
  if (qm->vType[si] == 2 ||
      qm->vType[si] == 0) qm->valence[si][1] = qm->star[si]->nQ + 2; //adjusted based on vertex type, eg, a corner is regular when val = 2
  else if (qm->vType[si] == 3) qm->valence[si][1] = qm->star[si]->nQ + 1;
  else if (qm->vType[si] >= 5) qm->valence[si][1] = qm->star[si]->nQ - 1;
  for (i3[0] = i4[0] = i1[0] = i = 0; i < qm->star[si]->nQ; i++) {
      if (qm->star[si]->verts[2 * i + 1] == -1) continue;
      if      (qm->vType[qm->star[si]->verts[2 * i + 1] - 1] >= 4)
	i4[++i4[0]] = qm->star[si]->verts[2 * i + 1];
      else if (qm->vType[qm->star[si]->verts[2 * i + 1] - 1] >= 2)
	i3[++i3[0]] = qm->star[si]->verts[2 * i + 1];
      else
	i1[++i1[0]] = qm->star[si]->verts[2 * i + 1];
  }
  for (j = i = 0; i < i4[0]; i++)
    qm->valence[si][3 + (j++)] = i4[1 + i];
  for (i = 0; i < i3[0]; i++)
    qm->valence[si][3 + (j++)] = i3[1 + i];
  for (i = 0; i < i1[0]; i++)
    qm->valence[si][3 + (j++)] = i1[1 + i];
  return EGADS_SUCCESS;
}


static int EG_backupQuads(meshMap *qm, int *nq, int *qlist, Quad **quad)
{
  int   i, j, q, v, qcount, *qaux = NULL, k, add = 0;

  qaux      = EG_alloc(qm->totQ * sizeof(int));
  if (qaux == NULL) return EGADS_MALLOC;
  for (qcount = q = 0; q < *nq; q++) {
      if (qlist[q] == -1) continue;
      for (i = 0; i < 4; i++) {
	  v    = qm->qIdx[4 * (qlist[q] - 1) + i] - 1;
	  if (qm->star[v] == NULL) {
	      EG_free (qaux);
	      printf("Star for vertex %d is NULL !!\n ", v+ 1);
	      return EGADS_MALLOC;
	  }
	  for (j = 0; j < qm->star[v]->nQ; j++) {
	      if (qm->star[v]->quads[j] == -1) continue;
	      for (add = k = 0; k < qcount; k++) {
		  if (qaux[k] == qm->star[v]->quads[j]) {
		      add = 1;
		      break;
		  }
	      }
	      if (add == 0 || qcount == 0) qaux[qcount++] = qm->star[v]->quads[j];
	  }
      }
  }
  (*quad)      = EG_alloc(qcount * sizeof(Quad));
  if ((*quad) == NULL) {
      EG_free(qaux);
      return EGADS_MALLOC;
  }
  for (q = 0; q < qcount; q++) {
      (*quad)[q].id = qaux[q];
      for (j = 0; j < 4; j++) {
	  (*quad)[q].qadj [j] = qm->qAdj[4 * (qaux[q] - 1) + j];
	  (*quad)[q].verts[j] = qm->qIdx[4 * (qaux[q] - 1) + j];
      }
  }
  EG_free(qaux);
  *nq = qcount;
  return EGADS_SUCCESS;
}


static int EG_restoreQuads(meshMap *qm, Quad *quad, int nq)
{
  int i, j, *vid = NULL, k, stat;

  vid      = EG_alloc(qm->totV * sizeof(int));
  if (vid == NULL) return EGADS_MALLOC;
  for (k  = i = 0; i < nq; i++) {
      if (quad[i].id == -1) continue;
      for (j = 0; j < 4; j++) {
	  qm->qAdj[4 * (quad[i].id - 1) + j] = quad[i].qadj [j];
	  qm->qIdx[4 * (quad[i].id - 1) + j] = quad[i].verts[j];
	  if (k == qm->totV) continue;
	  if (inList(k, vid, quad[i].verts[j]) == -1) {
	      qm->valence[quad[i].verts[j] - 1][0] = quad[i].id;
	      vid[k++] = quad[i].verts[j];
	  }
      }
  }
  for (i   = 0; i < k; i++) {
      stat = setValence(qm, vid[i]);
      if (stat != EGADS_SUCCESS) {
	  printf(" FAILED AT SETTING VALENCE IN %d \n ", vid[i] );
	  EG_free (vid);
	  return stat;
      }
  }
  EG_free (vid);

  return EGADS_SUCCESS;
}


static int EG_centroid(meshMap *qm, int n, int *list, double *quv, int usequv, int noV)
{
  int i, sum = 0;
  double  xyz0[18], xyz1[18];

  if (usequv == 1) {
      i      = EG_evaluate(qm->face, quv, xyz0);
      if (i != EGADS_SUCCESS ) return i;
      sum    = 1;
  } else {
      quv [0] = quv [1] = 0.0;
      xyz0[0] = xyz0[1] = xyz0[2] = 0.0;
  }
  for (i = 0; i < n; i++) {
      if ((qm->vType  [list[i] -1]    == -1 &&
	   qm->valence[list[i] -1][2] ==  2) ||
	                      list[i] == noV) continue;
      xyz0[0] += qm->xyzs[3 * (list[i] - 1)    ];
      xyz0[1] += qm->xyzs[3 * (list[i] - 1) + 1];
      xyz0[2] += qm->xyzs[3 * (list[i] - 1) + 2];
      quv[0]  += qm->uvs [2 * (list[i] - 1)    ];
      quv[1]  += qm->uvs [2 * (list[i] - 1) + 1];
      sum++;
  }
  xyz0[0] /= (double) sum;
  xyz0[1] /= (double) sum;
  xyz0[2] /= (double) sum;
  quv [0] /= (double) sum;
  quv [1] /= (double) sum;
  i        = EG_invEvaluateGuess(qm->face, xyz0, quv, xyz1);
  if (i   != EGADS_SUCCESS  || quv[0] < qm->range[0] ||
      quv[0] > qm->range[1] || quv[1] < qm->range[2] ||
      quv[1] > qm->range[3]) i = EG_invEvaluate(qm->face, xyz0, quv, xyz1);
  if (i != EGADS_SUCCESS) {
      printf("EG_centroid :: EG_invEvaluate at %lf %lf %lf is %d!!\n ",
	     xyz0[0], xyz0[1], xyz0[2], i);
      return i;
  }
  if      (quv[0] < qm->range[0]) quv[0] = qm->range[0];
  else if (quv[0] > qm->range[1]) quv[0] = qm->range[1];
  if      (quv[1] < qm->range[2]) quv[1] = qm->range[2];
  else if (quv[1] > qm->range[3]) quv[1] = qm->range[3];

  return EGADS_SUCCESS;
}

static int EG_projectToTangentPlane(double *normal, double *nO, double *p,
				    double *proj)
{
  double c, dotNN = 0.0, dotNP = 0.0, dist, lambda;

  dotNN      = sqrt (DOT(normal, normal));
  normal[0] /= dotNN;
  normal[1] /= dotNN;
  normal[2] /= dotNN;
  c          = DOT(normal, nO); // Equation plane: a*x + b*y + c*z = C
  dotNP      = DOT(normal, p);
  dist       = (p[0] - nO[0]) * (p[0] - nO[0]) + (p[1] - nO[1]) * (p[1] - nO[1]) +
      (p[2] - nO[2]) * (p[2] - nO[2]);
  if (fabs(dotNP - c) < EPS08 || sqrt(dist) < EPS08) {
      proj[0] = p[0]; proj[1] = p[1]; proj[2] = p[2];
      return EGADS_SUCCESS;
  }
  lambda  = (c - dotNP);
  proj[0] = p[0] + lambda * normal[0];
  proj[1] = p[1] + lambda * normal[1];
  proj[2] = p[2] + lambda * normal[2];
  dist    = DOT(normal, proj);
  if (fabs(dist - c) < EPS08) return EGADS_SUCCESS;
  printf("EG_projectToTangentPlane:  Plane O %lf %lf  %lf normal %lf %lf  %lf Poit %lf %lf %lf\n",
	 nO[0], nO[1], nO[2], normal[0], normal[1], normal[2], p[0], p[1], p[2]);
  printf("Projection is not in plane : distance %.16e \n", fabs( dist - c ) );
  return EGADS_GEOMERR;
}


static int EG_angleAtVnormalPlane(meshMap *qm, int vC, int v1, int v2,
				  double *angle)
{
  int    stat;
  double dot1, dot2, xyz[18], normal[3], cross[3], proj1[3], proj2[3];
  double centre[18], v01[3], v02[3], quv[2], norm1, norm2;

  stat = EG_evaluate(qm->face, &qm->uvs[2 * (vC - 1)], centre);
  if (stat != EGADS_SUCCESS) return stat;
  dot1 = centre[3]*centre[3] + centre[4]*centre[4] + centre[5]*centre[5];
  dot2 = centre[6]*centre[6] + centre[7]*centre[7] + centre[8]*centre[8];
  if (dot1 < EPS08 || dot2 < EPS08) {
      quv[0] = qm->uvs[ 2 * (vC - 1)    ];
      quv[1] = qm->uvs[ 2 * (vC - 1) + 1];
      stat = EG_centroid(qm, qm->valence[vC -1][2], &qm->valence[vC -1][3], quv, 1, 0);
      if (stat != EGADS_SUCCESS) return stat;
      stat = EG_evaluate(qm->face, quv, centre);
      if (stat != EGADS_SUCCESS) return stat;
      dot1 = centre[3]*centre[3] + centre[4]*centre[4] + centre[5]*centre[5];
      dot2 = centre[6]*centre[6] + centre[7]*centre[7] + centre[8]*centre[8];
      if (dot1 < EPS08 || dot2 < EPS08) return EGADS_GEOMERR;
  }
  v01[0] = centre[3]; v01[1] = centre[4]; v01[2] = centre[5];
  v02[0] = centre[6]; v02[1] = centre[7]; v02[2] = centre[8];
  if (qm->face->mtype == SREVERSE) {
      CROSS(v02, v01, normal);
  } else {
      CROSS(v01, v02, normal);
  }
  stat  = EG_evaluate(qm->face, &qm->uvs[2 * (v1 - 1)], xyz);
  stat += EG_projectToTangentPlane(normal, centre, xyz, proj1);
  stat += EG_evaluate(qm->face, &qm->uvs[2 * (v2 - 1)], xyz);
  stat += EG_projectToTangentPlane(normal, centre, xyz, proj2);
  if (stat != EGADS_SUCCESS) return stat;
  v01[0] = proj1[0] - centre[0];
  v01[1] = proj1[1] - centre[1];
  v01[2] = proj1[2] - centre[2];
  v02[0] = proj2[0] - centre[0];
  v02[1] = proj2[1] - centre[1];
  v02[2] = proj2[2] - centre[2];
  norm1  = DOT(v01, v01);
  norm2  = DOT(v02, v02);
  CROSS       (v01, v02, cross);
  dot1 = DOT(v01, v02) / sqrt ( norm1 * norm2 );
  dot2 = DOT(normal, cross);
  if      (dot1 >=  1.0) *angle = 0.0;
  else if (dot1 <= -1.0) *angle = PI;
  else                   *angle = acos(dot1);
  if (dot2 < 0)          *angle =  2.0 * PI - *angle;
  return EGADS_SUCCESS;
}


static int EG_angAtBdVert(meshMap *qm, int v, int *links, double *size)
{
  int i, j, k;

  *size = 0.0;
  if (qm->vType[v - 1] < 0) return EGADS_INDEXERR;
  if ( qm->star[v -1] == NULL) {
      printf(" Looking at corners: buildstar %d is NULL \n ", v);
      return EGADS_MALLOC;
  }
  for (links[0] = links[1] = k = i = 0; i < qm->star[v-1]->nQ; i++) {
      j = qm->star[v-1]->verts[2 * i + 1] - 1;
      if (qm->vType[j] != -1) k++;
      if (qm->star[v-1]->quads[i] == -1) {
	  links[1] = qm->star[v-1]->verts[2 * i + 1];
	  links[0] = qm->star[v-1]->verts[qm->star[v-1]->idxV[2 * i + 3]];
      }
  }
  if (k >= 3) {
      *size = PI * 0.5; // boundary vertex is connected to more than two bounds. Angle is fine
      return EGADS_SUCCESS;
  }
  else if (k != 2) {
      printf(" EG_angAtBdVert:: vertex %d is at surface bounds and connected only to another boundary vertex !!\n ", v);
      return EGADS_GEOMERR;
  }
  else return EG_angleAtVnormalPlane(qm, v, links[0], links[1], &(*size));
}


static int EG_nValenceCount(meshMap *qm, int q, int n)
{
  int i, count, val;

  for (count = i = 0; i < 4; i++) {
      val    = qm->valence[qm->qIdx[4 * (q - 1) + i] - 1][1];
      if (val < 0) return val;
      if (n > 5) {
	  if (val > 5) count++;
      }
      else if (n == val) count++;
  }
  return count;
}


static int EG_normalAtVertex(meshMap *qm, int v, double *normal, double *xyz) {
  int    i, j, stat;
  double norm1, norm2, vAB[3], vAC[3];

  if ( qm->star[v -1] == NULL ) {
      printf("EG_normalAtPoint star %d is NULL \n ", v );
      return EGADS_MALLOC;
  }
  j = qm->star[v - 1]->verts[0] - 1;
  stat      = EG_evaluate(qm->face, &qm->uvs[2 * j], xyz);
  if (stat != EGADS_SUCCESS) {
      printf(" EG_normalAtVertex :: EG_evaluate %d -->%d !!\n",
	     v, stat);
      return stat;
  }
  norm1 = xyz[3] * xyz[3] + xyz[4] * xyz[4] + xyz[5] * xyz[5];
  norm2 = xyz[6] * xyz[6] + xyz[7] * xyz[7] + xyz[8] * xyz[8];
  if (norm1 < EPS08 || norm2 < EPS08)
    {
      for (i = 1 ; i < qm->star[v-1]->nV; i++ ) {
	  j  = qm->star[v - 1]->verts[i] - 1;
	  if (j < 0 || qm->vType[j] != -1 ) continue;
	  stat      = EG_evaluate(qm->face, &qm->uvs[2 * j], xyz);
	  if (stat != EGADS_SUCCESS) {
	      printf(" EG_normalAtVertex :: EG_evaluate %d -->%d !!\n",
		     v, stat);
	      return stat;
	  }
	  norm1 = xyz[3] * xyz[3] + xyz[4] * xyz[4] + xyz[5] * xyz[5];
	  norm2 = xyz[6] * xyz[6] + xyz[7] * xyz[7] + xyz[8] * xyz[8];
	  if (norm1 > EPS08 && norm2 > EPS08) {
	      i = 0;
	      break;
	  }
      }
      if ( i != 0 ) {
	  printf("EG_normalAtVertex :: I couldn't find any non degenerate point around %d\n!!", v);
	  return EGADS_DEGEN;
      }
    }
  vAB[0] = xyz[3]; vAB[1] = xyz[4]; vAB[2] = xyz[5];
  vAC[0] = xyz[6]; vAC[1] = xyz[7]; vAC[2] = xyz[8];
  if (qm->face->mtype == SREVERSE) {
      CROSS(vAC, vAB, normal);
  } else {
      CROSS(vAB, vAC, normal);
  }
  norm1       = sqrt (normal[0] * normal[0] + normal[1] * normal[1] +
		      normal[2] * normal[2]);
  if ( norm1 < EPS08) {
      printf(" ATENTION !!!! NORM OF NORMAL IS %lf!! \n ", norm1 );
      printf(" du %lf %lf %lf \n ", xyz[3], xyz[4], xyz[5] );
      printf(" dv %lf %lf %lf \n ", xyz[6], xyz[7], xyz[8] );
      printf(" normal %lf %lf %lf \n ", normal[0], normal[1], normal[2]);
      return EGADS_GEOMERR;
  }
  normal[0] /= norm1;
  normal[1] /= norm1;
  normal[2] /= norm1;
  return EGADS_SUCCESS;
}


/*
 * Returns <0 EGADS_ERROR
 * = 0         (QA0) area is good
 * = 100       (QA1) area has obtuse angles
 * = 10000     (QA2) corner has obtuse angle or quad is very DELTA like
 * = 1000000   (QA3) selfintersects
 * = 100000000 (QACB) crosses the domain boundary!
 * ratio = ABC / ACD (ideal =1 so triangle split is equal forming parallelogram)
 */
static int EG_quadArea(meshMap *qm, /*@null@*/ double *normal, /*@null@*/ double *xyzO,
		       int qID, int vID, double *tr)
{
  int i, k, k1, stat, vA, vB, vC, vD, count, doublet, qV[4], ori[4], bv[4], lr[2], iv = -1, bvp = 0;
  int area = QA0, selfint = 1, cw, s1, s2, bound = 0;
//  int piv[6] = {1, 2, 2, 3, 1, 3};
  double pABCD[12], cross[3], qNormal[3], vAB[3], vAC[3], vAD[3], xyz[18], tria[8], vr[4];
  double c, dotNP, lambda, norm1, norm2, ang[4], dot, qa[2], pen = 1.0, ma = 0.0;

  qV[0] = qm->qIdx[4*(qID - 1)    ] - 1;
  qV[1] = qm->qIdx[4*(qID - 1) + 1] - 1;
  qV[2] = qm->qIdx[4*(qID - 1) + 2] - 1;
  qV[3] = qm->qIdx[4*(qID - 1) + 3] - 1;
  *tr   = 0.0;
  if (qNormal == NULL ) {
      qNormal[0] = normal[0];
      qNormal[1] = normal[1];
      qNormal[2] = normal[2];
      xyz    [0] = xyzO[0];
      xyz    [1] = xyzO[1];
      xyz    [2] = xyzO[2];
  } else {
      stat = EG_normalAtVertex(qm, qV[0] + 1, qNormal, xyz);
      if ( stat != EGADS_SUCCESS ) {
          printf("EG_quadArea :: EG_normalAtVertex %d !!\n ", stat);
          return stat;
      }
  }
  if ( vID <= 0 ){
      vA  = EG_quadArea(qm, qNormal, xyz, qID, qV[0] + 1, &vr[0]);
      i   = vA;
      *tr = vr[0];
      vB = EG_quadArea(qm, qNormal, xyz, qID, qV[1] + 1, &vr[1]);
      if (i < vB) {
          i = vB;
          *tr = vr[1];
      } else if (i == vB) *tr = MAX(*tr, vr[1]);
      vC = EG_quadArea(qm, qNormal, xyz, qID, qV[2] + 1, &vr[2]);
      if (i < vC) {
          i = vC;
          *tr = vr[2];
      } else if (i == vC) *tr = MAX(*tr, vr[2]);
      vD = EG_quadArea(qm, qNormal, xyz, qID, qV[3] + 1, &vr[3]);
      if (i < vD) {
          i = vD;
          *tr = vr[1];
      } else if (i == vD) *tr = MAX(*tr, vr[1]);
#ifdef DEBUG
      printf("\n=================================\n");
      printf(" AREAS %d %d %d %d \n ", vA, vB, vC, vD);
      printf(" ratios %lf %lf %lf %lf \n ", vr[0], vr[1], vr[2], vr[3]);
      printf(" RETURN %d RATIO %lf\n", i, *tr);
      printf("\n=================================\n");
#endif
      return i;
   }
  c      = DOT(qNormal, xyz);
  for (i = 0; i < 4; i++) {
      dotNP            = qNormal[0] * qm->xyzs[3 * qV[i]    ] +
                         qNormal[1] * qm->xyzs[3 * qV[i] + 1] +
                         qNormal[2] * qm->xyzs[3 * qV[i] + 2];
      lambda           = (c - dotNP);
      pABCD[3 * i    ] = qm->xyzs[3 * qV[i]    ] + lambda * qNormal[0];
      pABCD[3 * i + 1] = qm->xyzs[3 * qV[i] + 1] + lambda * qNormal[1];
      pABCD[3 * i + 2] = qm->xyzs[3 * qV[i] + 2] + lambda * qNormal[2];
  }
  qa[0]   = qa[1] = 0.0;
  tria[0] = tria[1] = tria[2] = tria[3] = 0.0;
  bv[0]   = bv[1] = bv[2] = bv[3] = -1;
  for (doublet = count = k = 0; k < 4; ++k) {
      ori[k]   = 1;
      vA       = k;
      if (qV[vA] + 1 == vID) {
          iv = k;
          if (qm->vType[qm -> valence[qV[vA]][3] -1] > 0) bv[k] = k;
      }
      if (qm -> vType[qV[vA]] != -1) {
          bound = 1;
          bv[k] = k;
      }
      if (qm->valence[qV[vA]][2] * qm -> vType[qV[vA]] == -2) {
          ori[k]     = 1;
          doublet    = 1;
      } else {
          vB     = (k + 1)%4;
          vC     = (k + 2)%4;
          vD     = (k + 3)%4;
          vAB[0] = pABCD[3 * vB    ] - pABCD[3 * vA    ];
          vAB[1] = pABCD[3 * vB + 1] - pABCD[3 * vA + 1];
          vAB[2] = pABCD[3 * vB + 2] - pABCD[3 * vA + 2];
          vAC[0] = pABCD[3 * vC    ] - pABCD[3 * vA    ];
          vAC[1] = pABCD[3 * vC + 1] - pABCD[3 * vA + 1];
          vAC[2] = pABCD[3 * vC + 2] - pABCD[3 * vA + 2];
          vAD[0] = pABCD[3 * vD    ] - pABCD[3 * vA    ];
          vAD[1] = pABCD[3 * vD + 1] - pABCD[3 * vA + 1];
          vAD[2] = pABCD[3 * vD + 2] - pABCD[3 * vA + 2];
          if (DOT(vAB,vAB) < EPS08 ||
              DOT(vAC,vAC) < EPS08 ||
              DOT(vAD,vAD) < EPS08 ) {
              selfint = 2;
              break;
          }
          CROSS(vAB, vAC, cross);
          tria[k] =  0.5 * sqrt ( DOT (cross, cross) ) ;
          if (DOT(qNormal, cross) < 0.0) ori[k] = -1;
          else {
              CROSS(vAC, vAD, cross);
              if (DOT(qNormal, cross) < 0.0) ori[k] = -1;
          }
          CROSS(vAB, vAD, cross);
          norm1    = DOT(vAB, vAB);
          norm2    = DOT(vAD, vAD);
          if (norm1 < EPS08 || norm2 < EPS08) dot = 1.0;
          else dot = DOT(vAB, vAD) / sqrt(norm1 * norm2);
          if      (fabs(dot - 1.0) < EPS08) ang[k] = 0.0;
          else if (fabs(dot + 1.0) < EPS08) ang[k] = PI;
          else                              ang[k] = acos(dot);
          if (DOT(qNormal, cross) < 0)      ang[k] = (2.0 * PI - ang[k]);
          ma = MAX(ma, ang[k]);
      }
      if (selfint == 2) {
          selfint = 1;
          break;
      }
      if (ori[k] == -1) continue;
      selfint = 0;
      count++;
  }
  qa[0] = tria[0] + tria[2];
  qa[1] = tria[1] + tria[3];
  if      (selfint == 1) area = QA3;
  else if (count   != 4) area = QA1;
  // now look for special cases: sharp interior corners, boundary quads...
  vr[0] = 0.0;
  vr[1] = 1.0;
  if ( doublet == 0 && qa[0] > EPS08 && qa[1] > EPS08
      && area <= QA1) {
      k  =  iv    % 2; // id of corresponding diagonal split
      k1 = (iv + 3)%4; // id of triangle with centre iv
      if (area == QA0) {
          for (s1 = i = 0; i < qm->star[vID-1]->nQ; i++)
            if(qm->vType[qm->star[vID-1]->verts[2 *i + 1] -1] != -1)
              s1++;
          if ( s1 >= 2 ) dot = 0.00;
          else dot = MAX(MAX(ang[iv], ang[(iv+1)%4]),ang[(iv+3)%4]);
          pen      = 3.554147 * (DEG175 - dot) / DEG5;
          pen      = (1.0 - erfc(pen)*0.5);
          if(tria[k] > EPS08 && tria[k + 2] > EPS08)
            vr[0] = MIN(tria[k] / tria[k + 2], tria[k + 2] / tria[k]);
          if( tria[(k1 + 2)%4] > EPS08 &&
              tria[ k1       ] < tria[(k1 + 2)%4])
            vr[1] = tria[k1] / tria[(k1+2)%4];
      } else {
          vA = 0;
          if (qa[0] > qa[1] ) vA = 1;
          vB = (vA + 1)%2;
          if (tria[k] + tria [k + 2] > qa[vA]) { // MIN AREA = qa[vA]
              if(tria[k] > EPS08 && tria[k + 2] > EPS08)
                vr[0] = -MIN(tria[k    ] / tria[k + 2],
                             tria[k + 2] / tria[k    ]);
          } else if   (tria[k1] < tria[(k1 + 2)%4]) {
              vr[0] = -tria[k1] / tria[(k1 + 2)%4];
          } else {
              vr[0] = -MIN(tria[k    ] / tria[k + 2],
                           tria[k + 2] / tria[k    ]);
          }
      }
      *tr = MIN(vr[0], vr[1]);
      if ( *tr > 0.0 && *tr < 1.e-05 ) *tr = 0.0;
  }
  for (k = 0 ; k < 4; k++ ) {
      if (bv[k] == -1) continue;
      if (bvp == 0 ) bvp = 1;
      else bvp = MAX(qm->vType[qV[bv[k]]], bvp);
      if (qm->vType[qV[bv[k]]]  == -1 ) continue;
      // check if we cross the domain boundary. This can happen even for valid quads!!
      for (lr[0] = lr[1] = i = 0; i < qm->star[qV[bv[k]]]->nQ; i++) {
          k1 = qm->star[qV[bv[k]]]->verts[2 * i + 1] - 1;
          if (qm->star[qV[bv[k]]]->quads[i] == -1) {
              lr[1] = qm->star[qV[bv[k]]]->verts[2 * i + 1];
              lr[0] = qm->star[qV[bv[k]]]->verts[qm->star[qV[bv[k]]]->idxV[2 * i + 3]];
              break;
          }
      }
      dotNP  = qNormal[0] * qm->xyzs[3 * (lr[0] -1)    ] +
               qNormal[1] * qm->xyzs[3 * (lr[0] -1) + 1] +
               qNormal[2] * qm->xyzs[3 * (lr[0] -1) + 2];
      lambda = (c - dotNP);
      xyz[0] = qm->xyzs[3 * (lr[0] -1)    ] + lambda * qNormal[0];
      xyz[1] = qm->xyzs[3 * (lr[0] -1) + 1] + lambda * qNormal[1];
      xyz[2] = qm->xyzs[3 * (lr[0] -1) + 2] + lambda * qNormal[2];
      if (lr[1] == qV[(bv[k] +1)%4] + 1 || lr[0] == qV[(bv[k] +3)%4] + 1) continue;
      for ( k1 = 1; k1 <= 2; k1++) {
          vB     = (bv[k] + k1    )%4;
          vC     = (bv[k] + k1 + 1)%4;
          if (qV[vB] + 1 == lr[0] || qV[vB] + 1 == lr[1] ||
              qV[vC] + 1 == lr[0] || qV[vC] + 1 == lr[1] )continue;
          vAB[0] = pABCD[3*vB    ] - pABCD[3 * bv[k]    ];
          vAB[1] = pABCD[3*vB + 1] - pABCD[3 * bv[k] + 1];
          vAB[2] = pABCD[3*vB + 2] - pABCD[3 * bv[k] + 2];
          vAC[0] = pABCD[3*vC    ] - pABCD[3 * bv[k]    ];
          vAC[1] = pABCD[3*vC + 1] - pABCD[3 * bv[k] + 1];
          vAC[2] = pABCD[3*vC + 2] - pABCD[3 * bv[k] + 2];
          vAD[0] = xyz[0]          - pABCD[3 * bv[k]    ];
          vAD[1] = xyz[1]          - pABCD[3 * bv[k] + 1];
          vAD[2] = xyz[2]          - pABCD[3 * bv[k] + 2];
          CROSS(vAB, vAC, cross);
          cw = 0;
          if (DOT(qNormal, cross) < 0 ) cw = 1;
          // check that BC and CD dont cross the boundary
          CROSS(vAB, vAD, cross);
          s1       = signbit(DOT(qNormal, cross));
          CROSS(vAD, vAC, cross);
          s2     = signbit(DOT(qNormal, cross));
          norm1    = DOT(vAC, vAC);
          if ((cw == 0 && s1 == 0 && s2 == 0) ||
              (cw == 1 && s1 == 1 && s2 == 1)) {
#ifdef DEBUG
              printf(" OJO! QUAD %d LINK %d %d CROSSES THE DOMAIN BOUNDS ! REJECT\n ", qID, qV[vB] + 1, qV[vC] + 1);
              printf(" LINK AB AL %d %d  %d %d = %d\n ", qV[bv[k]] + 1, qV[vB] + 1, qV[bv[k]] + 1, lr[0], s1);
              printf(" LINK AL AC %d %d  %d %d = %d\n ", qV[bv[k]] + 1, lr[0], qV[bv[k]] + 1, qV[vC] + 1, s2);
              printQuad(qm, qID);
#endif
              selfint = 2;
              break;
          }
      }
      if (selfint == 2) break;
  }
  if ( selfint == 2 ||
      (selfint == 1 && bvp == 4  )) area = QACB;
  if (area >= QA3) *tr = 0.0;
  if (doublet == 1 ) {
      if (area < QA3) {
          area = QA0;
          *tr  = 1.0;
      }
  } else {
      if ( area < QA3 &&
          ((bvp == 4 && ma      > PI) ||
           (bvp == 1 && ang[iv] > PI) ||
                            *tr < -0.5)) area = QA2;
      if (area == QA0 && *tr < EPS08) area = QA1;
      if (area == QA1 && bound == 1 ) area *= 2;
  }
#ifdef DEBUG
  if (area != QA0) {
      printf("\n\n QUAD %d SIGN AREA %d VAL %lf %lf --> RATIO %lf\n",
             qID, area, qa[0],qa[1], *tr);
      if (area == QA2 )
        printf(" ************ SHARP CORNER %d = %lf ***************\n", bvp, ma);
      else if (area == QA3)
        printf(" ************ SELFINTERSECTS        ***************\n");
      else if (area == QACB)
        printf(" ************ CROSSES DOMAIN BOUNDS ***************\n");
      else if ( area % QA1 == 0)
        printf(" ************ OBTUSE QUAD           ***************\n");
      printQuad(qm, qID);
      /*for (k = 0; k < 4; ++k)
	printf("Vertex %d (T %d val %d) and orientation %d \n ",
	       qV[k] + 1,qm->vType[qV[k]], qm->valence[qV[k]][2], ori[k]);
      for (k = 0; k < 4; ++k)
        printf("%lf %lf %lf %d\n ", pABCD[3 * k], pABCD[3 * k + 1], pABCD[3 * k + 2], qV[k] + 1);
      k = 0;
      printf("%lf %lf %lf %d\n ", pABCD[3 * k], pABCD[3 * k + 1], pABCD[3 * k + 2], qV[k] + 1);*/
  }
#endif
  /*if ( iv % 2 == 0) {
      printf(" TRIANGLE SPLIT\n ");
      printf("%lf %lf %lf %d\n",pABCD[0],  pABCD[1], pABCD[2], qV[0] + 1);
      printf("%lf %lf %lf %d\n",pABCD[3],  pABCD[4], pABCD[5], qV[1] + 1);
      printf("%lf %lf %lf %d\n",pABCD[6],  pABCD[7], pABCD[8], qV[2] + 1);
      printf("%lf %lf %lf %d\n\n\n",pABCD[0],  pABCD[1], pABCD[2], qV[0] + 1);
      printf("%lf %lf %lf %d\n",pABCD[0],  pABCD[1], pABCD[2], qV[0] + 1);
      printf("%lf %lf %lf %d\n",pABCD[6],  pABCD[7], pABCD[8], qV[2] + 1);
      printf("%lf %lf %lf %d\n",pABCD[9],  pABCD[10], pABCD[11], qV[3] + 1);
      printf("%lf %lf %lf %d\n\n\n",pABCD[0],  pABCD[1], pABCD[2], qV[0] + 1);
  }else {
      printf("%lf %lf %lf %d\n",pABCD[3],  pABCD[4], pABCD[5], qV[1] + 1);
      printf("%lf %lf %lf %d\n",pABCD[6],  pABCD[7], pABCD[8], qV[2] + 1);
      printf("%lf %lf %lf %d\n",pABCD[9],  pABCD[10], pABCD[11], qV[3] + 1);
      printf("%lf %lf %lf %d\n\n\n",pABCD[3],  pABCD[4], pABCD[5], qV[1] + 1);
      printf("%lf %lf %lf %d\n",pABCD[3],  pABCD[4], pABCD[5], qV[1] + 1);
      printf("%lf %lf %lf %d\n",pABCD[9],  pABCD[10], pABCD[11], qV[3] + 1);
      printf("%lf %lf %lf %d\n",pABCD[0],  pABCD[1], pABCD[2], qV[0] + 1);
      printf("%lf %lf %lf %d\n\n\n",pABCD[3],  pABCD[4], pABCD[5], qV[1] + 1);
  }*/
#ifdef DEBUG
      printf("\n\nQ %d  FROM %d %lf %lf ---> %.16f (PENALTY %lf)\n ", qID, vID, vr[0], vr[1], *tr, pen);
      printf(" VERTEX %d bvs %d %d %d %d\n", vID, bv[0], bv[1], bv[2], bv[3]);
      printf(" ANGLES %lf %lf %lf %lf\n\n", ang[0], ang[1], ang[2], ang[3]);
#endif
  if (area == QA0) *tr *= pen;
  return area;
}


static void EG_placeVertex(meshMap *qm, int vID, double pass ) {
  int q, v, j, i, k, i0, i1, ja, jb, stat, qi, update = 0,
      la = -1, lb, nt, block = 0, round, doublet = 0, bt = -1, links[2];
  double uv[6], normal[3], xyz[18], sra[3], srb[3], qr;
  int *aArea = NULL, *bArea = NULL;

#ifdef DEBUG
  char buffer[100];
  int d1, d2;
  double p[2], pos[18];
  FILE *fout;
#endif
  v          = vID - 1;
  if (qm->vType[v] != -1 || qm ->star[v] == NULL) {
#ifdef DEBUG
      printf(" VERTEX %d is boundary star %p!\n ", vID, qm-> star[v]);
#endif
      return;
  }
#ifdef DEBUG
  printf(" Place Coords for vertex %d\n",vID);
  snprintf(buffer, 100,"PLACEStart_%d_%d", qm->plotcount, vID);
  printf("Writing in %s\n ", buffer);
  fout = fopen(buffer,"w");
  if (fout != NULL ) {
      for (d1 = 0; d1 < qm->star[v]->nQ; d1++) {
          p[0] = p[1] = 0.0;
          for (d2 = 0; d2 <= 4; d2++) {
              k = qm->qIdx[ 4 * (qm->star[v]->quads[d1] -1) + d2%4 ] - 1;
              fprintf(fout, "%lf %lf %lf %d\n",  qm->xyzs[3*k  ],
                      qm->xyzs[3*k + 1], qm->xyzs[3*k + 2], k + 1 );
              if ( d2 == 4 ) break;
              p[0] += 0.25 * qm->uvs[2 * k    ];
              p[1] += 0.25 * qm->uvs[2 * k + 1];
          }
          fprintf(fout,"\n\n");
          EG_evaluate(qm->face, p, pos);
          fprintf(fout, "%lf %lf %lf %d\n", pos[0], pos[1], pos[2],
                  qm->star[v]->quads[d1]);
          fprintf(fout,"\n\n");
      }
      fclose (fout);
  }
#endif
  if (qm->valence[v][2] == 2) {
      i0    = qm->valence[v][3] - 1;
      i1    = qm->valence[v][4] - 1;
      uv[0] = 0.5 * (qm->uvs[2 * i0    ] + qm->uvs[2 * i1    ]);
      uv[1] = 0.5 * (qm->uvs[2 * i0 + 1] + qm->uvs[2 * i1 + 1]);
      updateVertex(qm, vID, uv);
#ifdef DEBUG
      printf("PLACE COORD DOUBLET -->LEAVE WITH COORDS\n");
      printVertex(qm, vID);
#endif
      return;
  }
  i0 = i1 = 0;
  for (i = 1 ; i < qm->star[v]->nV; i++) {
      k = qm->vType[qm->star[v]->verts[i] - 1];
      bt = MAX(bt, k);
      if (qm -> valence[qm->star[v]->verts[i] - 1][2] == 2) doublet = 1;
      if (i%2 == 0) continue;
      if (k >= 4 && la == -1) la = i;
  }
  uv [0] = qm->uvs[2 * v    ];
  uv [1] = qm->uvs[2 * v + 1];
  uv [4] = uv[0];
  uv [5] = uv[1];
  sra[0] = sra[1] = sra[2] = 1.0;
  stat  = EG_normalAtVertex(qm, vID, normal, xyz);
  aArea = EG_alloc(qm->star[v]->nQ * sizeof(int));
  bArea = EG_alloc(qm->star[v]->nQ * sizeof (int   ));
  if (stat  != EGADS_SUCCESS || aArea == NULL ||
      bArea == NULL) return ;
  for (ja = j = 0 ;j < qm->star[v]->nQ; j++ ) {
      qi  = qm->star[v]->quads[j] -1;
      aArea[j] = EG_quadArea(qm, normal, xyz, qi + 1, vID, &qr);
      ja  += aArea[j];
      if (aArea[j] >= QA3) continue;
      if (qr > 0.0) {
          if (qm->vType[qm -> qIdx[4 * qi    ]-1] != -1 ||
              qm->vType[qm -> qIdx[4 * qi + 1]-1] != -1 ||
              qm->vType[qm -> qIdx[4 * qi + 2]-1] != -1 ||
              qm->vType[qm -> qIdx[4 * qi + 3]-1] != -1 )
            sra[2] = MIN(sra[2], qr);
          else sra[0] = MIN(sra[0], qr);
      } else sra[1] = MIN(sra[1], qr);
  }
  if ( pass > 0.0 ) {
      sra[0] = MIN(sra[0], 1.0 - MIN(pass, 0.25));
      sra[2] = MIN(sra[2], 1.0 - MIN(pass, 0.25));
  }
  srb[0] = srb[1] = srb[2] = 1.0;
  for (round = 0; round < 2; round++) {
      if (round == 0) {
          nt = qm->star[v]->nQ;
          if (nt %2 == 0) nt /= 2;
      }
      else nt = qm->star[v]->nV;
#ifdef DEBUG
      printf (" ReSET ROUND %d -------------- la %d \n ", round, la);
#endif
      for (i1 = i0 = q = 0; q <= nt; q++) {
          uv[2] = uv[0];
          uv[3] = uv[1];
          updateVertex(qm, vID, uv);
#ifdef DEBUG
          printf(" q %d / %d\n ", q, nt);
#endif
          block = 0;
          if (q == 0 ) {
              EG_centroid( qm, qm->star[v]->nQ,
                           &qm->valence[v][3],&uv[2], (round + 1)%2, 0);
              if (round == 0) {
                  stat      = EG_normalAtVertex(qm, vID, normal, xyz);
                  if (stat != EGADS_SUCCESS) break;
              }
          } else {
              if ( round == 0) {
                  if (qm->star[v]->nQ == 3 ) break;
                  i0 = qm->star[v]->verts[2 * (q - 1) + 1] - 1;
                  i1 = qm->star[v]->verts[qm->star[v]->idxV[2 * (q - 1) + 5]] - 1;
                  j  = qm->star[v]->verts[qm->star[v]->idxV[2 * (q - 1) + 7]] - 1;
                  if (j == i0 || qm->vType[i0] * qm->valence[i0][2] == -2 ||
                                 qm->vType[i1] * qm->valence[i1][2] == -2 ) continue;
                  if ((qm -> vType[i0]  > 0 && qm -> vType[i1]  > 0 ) &&
                      (qm -> vType[i0] <= 4 && qm -> vType[i1] <= 4)) block = 1;
                  uv[2] = 0.5 *(qm->uvs[2 * i0    ] + qm->uvs[2 * i1    ]);
                  uv[3] = 0.5 *(qm->uvs[2 * i0 + 1] + qm->uvs[2 * i1 + 1]);
              } else {
                  if ( la != -1 ) {
                      lb = qm->star[v]->idxV [la + q];
#ifdef DEBUG
                      printf(" LA %d --> LOOK FOR VERTICES %d ( LB %d)\n ", la, q,
                             qm->star[v]->verts[lb]);
#endif
                      if((abs(lb - la) < 3 ||
                          abs(lb - la) > nt - 1 - 3 ))  continue;
                      i0 = qm->star[v]->verts[la] - 1;
                      i1 = qm->star[v]->verts[lb] - 1;
                      if (qm->vType[i1] * qm->valence[i1][2] == -2 ) continue;
                      uv[2] = 0.5 *(qm->uvs[2 * i0    ] + qm->uvs[2 * i1    ]);
                      uv[3] = 0.5 *(qm->uvs[2 * i0 + 1] + qm->uvs[2 * i1 + 1]);
                  } else {
                      if (q - 1 >= qm->star[v]-> nQ ||
                          (ja == QA0 && sra[0] >= 0.1 &&
                              sra[2] >= 0.1)) break;
                      if (qm->star[v]-> nQ  == 3) {
                          i0 = v;
                          i1 = qm->star[v]->verts[2 * (q - 1) + 1] - 1;
                          if (ja == QA0 && qm->vType[i1] != -1) continue;
                          uv[2] = 0.5 *(qm->uvs[2 * i0    ] + qm->uvs[2 * i1    ]);
                          uv[3] = 0.5 *(qm->uvs[2 * i0 + 1] + qm->uvs[2 * i1 + 1]);
                      } else {
#ifdef DEBUG
                          printf(" CENTROID EXCLUDING %d \n",
                                 qm->star[v]->verts[2 *(q - 1) + 1]);
#endif
                          EG_centroid( qm, qm->star[v]->nQ,
                                       &qm->valence[v][3],&uv[2], 1,
                                       qm->star[v]->verts[2 *(q - 1) + 1]);
                      }
                  }
              }
          }
          updateVertex(qm, vID, &uv[2]);
          srb[0] = srb[1] = srb[2] = 1.0;
          for (jb = j = 0; j < qm->star[v]->nQ; j++) {
              qi        = qm->star[v]->quads[j] -1;
              bArea[j] = EG_quadArea(qm, normal, xyz, qi + 1,vID,
                                      &qr);
              jb += bArea[j];
              if (bArea[j] >= QA3) continue;
              if (qr > 0.0) {
                  if (qm->vType[qm -> qIdx[4 * qi    ]-1] != -1 ||
                      qm->vType[qm -> qIdx[4 * qi + 1]-1] != -1 ||
                      qm->vType[qm -> qIdx[4 * qi + 2]-1] != -1 ||
                      qm->vType[qm -> qIdx[4 * qi + 3]-1] != -1 )
                    srb[2]    = MIN(srb[2], qr);
                  else srb[0] = MIN(srb[0], qr);
              }
              else srb[1] = MIN(srb[1], qr);
          }
          if ( pass > 0.0) {
              srb[0] = MIN(srb[0], 1.0 - MIN(pass, 0.25));
              srb[2] = MIN(srb[2], 1.0 - MIN(pass, 0.25));
          }
#ifdef DEBUG
          if (q == 0 ) printf (" CENTROID\n");
          else printf("\n\n q - %d / %d MID POINTS BETWEEN %d %d \n ",q ,nt,  i0 + 1, i1 +1);
          printVertex(qm, vID);
          printf("NOW  AREA %d RATIOS POS %lf NEG %lf QB %lf\n "
                 "BEST AREA %d RATIOS POS %lf NEG %lf QB %lf\n ",
                 jb, srb[0], srb[1], srb[2],
                 ja, sra[0], sra[1], sra[2]);
#endif
          k = 0;
#ifdef DEBUG
              printf(" PASS ANGLE %lf QINT %d!!! \n ", pass, bt);
#endif
          if (block == 1 && (jb != QA0 ||
              srb[0] < 0.1 || srb[2] < 0.1)) block = 0;
          if (pass < 0.0 && jb < QA2 && bt < 4) block = 2;
          if (jb == QA0 && srb[2] > 0.75 && srb[0] > 0.75 ) block = 1;
          if      (block > 0) k = 1;
          else if (jb <  ja) k = 1;
          else if (ja == jb) {
              if  (ja == QA0) {
#ifdef DEBUG
                  printf(" AREA 0 COMPARE GROUPS %lf %lf VS %lf %lf\n",
                         srb[0], sra[0], srb[2], sra[2]);
#endif
                  if (q == 0 && bt < 4 && pass > 0.0 &&
                         srb[0] > pass && srb[2] > pass) k = 1;
                  else if (srb[2] > sra[2] &&
                           srb[0] > sra[0]) k = 1;
                  if (k == 0 && fabs (srb[2] - sra[2]) < 1.e-05) {
                                  if (srb[0] > sra[0]) k = 1;
                  }
                  if (k == 0 && fabs (srb[0] - sra[0]) < 1.e-05) {
                                  if (srb[2] > sra[2]) k = 1;
                  }
                  if (k == 0 && (sra[2] < 0.1 || sra[0] < 0.1)) {
                      if (srb[0] >= 0.1 && srb[2] >= 0.1) k = 1;
                  }
                  if (k == 0 && srb[2] < 0.1) {
                      if (sra[2] > 0.1 ) continue;
                      if (srb[2] > sra[2] ||
                          srb[0] > sra[0] + 0.1) k = 1;
                  }
                  if (k == 0 && srb[0] < 0.1) {
                      if (sra[0] > 0.1 ) continue;
                      if (srb[0] > sra[0] ||
                          srb[2] > sra[2] + 0.1) k = 1;
                  }
                  if (k == 0) {
                      if (fabs(srb[2] - sra[2]) < 0.1) {
                          if (srb[0] >= 0.1 + sra[0] ) k = 1;
                          else if (srb[0] - sra[0] > sra[2] - srb[2]) k = 1;
                      } else if (fabs(srb[0] - sra[0]) < 0.1) {
                          if (srb[2] >= 0.1 + sra[2] ) k = 1;
                          else if (srb[2] - sra[2] > sra[0] - srb[0]) k = 1;
                      }else if (MIN (srb[0], srb[2] ) >
                      MIN (sra[0], sra[2])) k = 1;
                  }
                      /*
                  } else if (sra[2] + srb[2] > sra[0] + srb[0] ) {
                          if (srb[0] > sra[0]) k = 1;
                          else if (srb[2] > sra[2] &&
                                  (srb[2] - sra[2]) > (sra[0] - srb[0])) k = 1;
                  } else {
                      if (srb[2] > sra[2]) k = 1;
                      else if (srb[0] > sra[0] &&
                          (srb[0] - sra[0]) > (sra[2] - srb[2])) k = 1;
                  }*/
              } else if (ja < QA3) {
                  if (srb[1] < 0.0) {
                      if (sra[1] < 0.0 && srb[1] > sra[1]) k = 1;
                  } else {
                      if      (sra[1] < 0.0) k = 1;
                      else if (
                          (srb[0] > sra[0] && fabs(srb[0] - sra[0] ) > EPS08) ||
                          (srb[2] > sra[2] && fabs(srb[2] - sra[2] ) > EPS08))
                           k = 1;
                  }
              }
              else if (srb[2] + srb[0] >  sra[2] + sra[0]) k = 1;
          }
          if (k == 0) continue;
#ifdef DEBUG
          printf("UPDATE BEST AREA ja %d jb %d \n", ja,jb);
          printf("A    %lf %lf %lf\n ", sra[0], sra[1], sra[2]);
          printf("BEST %lf %lf %lf\n ", srb[0], srb[1], srb[2]);
#endif
          update = 1;
          uv[4]  = uv[2];
          uv[5]  = uv[3];
          sra[0] = srb[0];
          sra[1] = srb[1];
          sra[2] = srb[2];
          for (j = 0; j < qm->star[v]->nQ; j++) aArea [j] = bArea [j];
          ja = jb;
          if (ja == QA0     &&
              (doublet == 1 || block == 1 ||
              (sra[0] > pass && sra[2] > pass)))  block = 1;
          if (block > 0){
#ifdef DEBUG
              printf("BLOCK POSITION QINT %d block %d (pass %lf) src %lf %lf %lf\n",
                     bt, block,pass, srb[0], srb[1], srb[2]);
#endif
              break;
          }
      }
      if (block > 0) break;
  }
  if (update == 1) updateVertex(qm, vID, &uv[4]);
  else updateVertex(qm, vID, uv);
#ifdef DEBUG
  snprintf(buffer, 100,"PLACE_%d_%d", qm->plotcount, vID);
  printf("Writing in %s\n ", buffer);
  fout = fopen(buffer,"w");
  if (fout != NULL ) {
      for (d1  = 0; d1 < qm->star[v]->nQ; d1++ ) {
          p[0] = p[1] = 0.0;
          for (d2 = 0; d2 <= 4; d2++) {
              k   = qm->qIdx[ 4 * (qm->star[v]->quads[d1] -1) + d2%4 ] - 1;
              fprintf(fout, "%lf %lf %lf %d\n",  qm->xyzs[3*k  ],
                      qm->xyzs[3*k + 1], qm->xyzs[3*k + 2], k + 1 );
              if (d2 == 4) break;
              p[0] += 0.25 * qm->uvs[2 * k    ];
              p[1] += 0.25 * qm->uvs[2 * k + 1];
          }
          fprintf(fout,"\n\n");
          EG_evaluate(qm->face, p, pos);
          fprintf(fout, "%lf %lf %lf %d\n", pos[0], pos[1], pos[2],
                  qm->star[v]->quads[d1]);
          fprintf(fout,"\n\n");
      }
      fclose (fout);
  }
#endif
  EG_free(bArea);
  EG_free(aArea);
  return;
}


int EG_createMeshMap(bodyQuad *bodydata, int uvtype)
{
  int          f, stat = 0, stat2 = 0, j, q, k, kk, kOK, len;
  int          ntri, nquad, e4[4] ;
  const int    *tris, *tric, *ptype, *pindex;
  double       angle;
  const double *xyzs, *uvs;
  int          qV[6]    = { 0, 1, 2, 5, 0, 1};
  int          qLoop[5] = { 0, 1, 2, 3, 0   };

  bodydata->qm = (meshMap **) EG_alloc(bodydata->nfaces * sizeof(meshMap *));
  if (bodydata->qm == NULL) return  EGADS_MALLOC;
  for (f = 0; f < bodydata->nfaces; ++f) {
      bodydata->qm[f] = (meshMap *) EG_alloc(sizeof(meshMap));
      if (bodydata->qm[f] == NULL) {
	  printf("Create Quad Map: MALLOC error!! \n ");
	  stat2 = EGADS_MALLOC;
	  continue;
      }
      bodydata->qm[f]->fID       = f + 1;
      bodydata->qm[f]->plotcount = 0;
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
      bodydata->qm[f]->oriV    = len;
      bodydata->qm[f]->oriQ    = nquad;
      bodydata->qm[f]->sizeV   = 2 * len;
      bodydata->qm[f]->sizeQ   = 2 * nquad;
      bodydata->qm[f]->totV    = len;
      bodydata->qm[f]->totQ    = nquad;
      bodydata->qm[f]->qInv    = NULL;
      bodydata->qm[f]->xyzs    = (double *) EG_alloc(3*(2 * len  )*sizeof(double));
      bodydata->qm[f]->uvs     = (double *) EG_alloc(2*(2 * len  )*sizeof(double));
      bodydata->qm[f]->vType   = (int    *) EG_alloc(  (2 * len  )*sizeof(  int ));
      bodydata->qm[f]->qIdx    = (int    *) EG_alloc(4*(2 * nquad)*sizeof(  int ));
      bodydata->qm[f]->qAdj    = (int    *) EG_alloc(4*(2 * nquad)*sizeof(  int ));
      bodydata->qm[f]->remQ    = (int    *) EG_alloc(  (2 * nquad)*sizeof(  int ));
      bodydata->qm[f]->remV    = (int    *) EG_alloc(  (2 * len  )*sizeof(  int ));
      bodydata->qm[f]->valence = (int   **) EG_alloc(  (2 * len  )*sizeof(  int*));
      bodydata->qm[f]->star    = (vStar **) EG_alloc(  (2 * len  )*sizeof(vStar*));
      if (bodydata->qm[f]->qIdx  == NULL || bodydata->qm[f]->qAdj    == NULL ||
	  bodydata->qm[f]->xyzs  == NULL || bodydata->qm[f]->uvs     == NULL ||
	  bodydata->qm[f]->vType == NULL || bodydata->qm[f]->remQ    == NULL ||
	  bodydata->qm[f]->remV  == NULL || bodydata->qm[f]->valence == NULL ||
	  bodydata->qm[f]->star  == NULL) {
	  stat2 = EGADS_MALLOC;
	  continue;
      }
      bodydata->qm[f]->remQ[0]  = 0;
      bodydata->qm[f]->remV[0]  = 0;
      bodydata->qm[f]->invsteps = 0;
      for (j = 0; j < 2 * len; j++) {
	  bodydata->qm[f]->star[j]    = NULL;
	  bodydata->qm[f]->valence[j] = (int *) EG_alloc(3 * sizeof(int));
	  if (bodydata->qm[f]->valence[j] == NULL) {
	      stat2 = EGADS_MALLOC;
	      continue;
	  }
	  bodydata->qm[f]->valence[j][2] = 0;
      }
      bodydata->qm[f]->face    = bodydata->faces[f];
      stat   = EG_getRange(bodydata->qm[f]->face, bodydata->qm[f]->range, &j);
      for (j = 0; j < len; j++) {
	  bodydata->qm[f]->vType  [j]         = ptype[j];
	  bodydata->qm[f]->uvs    [2 * j    ] = uvs[2 * j    ];
	  bodydata->qm[f]->uvs    [2 * j + 1] = uvs[2 * j + 1];
	  if (uvs[2 * j    ] < bodydata->qm[f]->range[0] - 1.e-6 ||
	      uvs[2 * j    ] > bodydata->qm[f]->range[1] + 1.e-6 ||
	      uvs[2 * j + 1] < bodydata->qm[f]->range[2] - 1.e-6 ||
	      uvs[2 * j + 1] > bodydata->qm[f]->range[3] + 1.e-6) {
	      printf(" VERTEX %d = %f  %f  is out of range !! \n ",
		     j + 1, uvs[2 * j], uvs[2 * j + 1]);
	      printf(" RANGE FACE %d --- > %f  %f  %f  %f\n ", bodydata->qm[f]->fID,
		     bodydata->qm[f]->range[0],bodydata->qm[f]->range[1],
		     bodydata->qm[f]->range[2],bodydata->qm[f]->range[3]);
	      stat2 = EGADS_INDEXERR;
	      break;
	  }
	  bodydata->qm[f]->xyzs[3 * j    ] = xyzs[3 * j    ];
	  bodydata->qm[f]->xyzs[3 * j + 1] = xyzs[3 * j + 1];
	  bodydata->qm[f]->xyzs[3 * j + 2] = xyzs[3 * j + 2];
      }
      if (stat2 != EGADS_SUCCESS) continue;
      for (j = 0; j < nquad; j++)
	for (k = 0; k < 4; ++k) {
	    bodydata->qm[f]->qIdx[4*j + k] = tris[6*j + qV[k+1]];
	}
      for (j = 0; j < nquad; j++)
	{
	  kk  = 0;
	  kOK = 0;
	  q   = 0;
	  while (q < nquad) {
	      if (q == j) {
		  if (j == nquad - 1 && kk < 4) {
		      bodydata->qm[f]->qAdj[4*j + kk++] = -1;
		      if (kk == 4) break;
		      q = 0;
		  }
		  else q++;
	      }
	      if  (q == nquad) break;
	      for (k = 0; k < 4; ++k) {
		  if ((bodydata->qm[f]->qIdx[4*j + qLoop[kk    ]] ==
		      bodydata->qm[f]->qIdx[4*q + qLoop[k     ]] ||
		      bodydata->qm[f]->qIdx[4*j + qLoop[kk    ]] ==
			  bodydata->qm[f]->qIdx[4*q + qLoop[k  + 1]]) &&
		      (bodydata->qm[f]->qIdx[4*j + qLoop[kk + 1]] ==
			  bodydata->qm[f]->qIdx[4*q + qLoop[k     ]] ||
			  bodydata->qm[f]->qIdx[4*j + qLoop[kk + 1]] ==
			      bodydata->qm[f]->qIdx[4*q + qLoop[k  + 1]])) {
		      bodydata->qm[f]->qAdj[4*j + kk++] = q + 1;
		      q   = -1;
		      kOK =  1;
		      k   =  4;
		      if (kk == 4) q = nquad;
		  }
	      }
	      if ((kOK == 0) && (q >= nquad -1)) {
		  bodydata->qm[f]->qAdj[4*j + kk++] = -1;
		  q = -1;
		  if (kk == 4) break;
	      }
	      else kOK = 0;
	      q++;
	  }
	  if (kOK == 0 && kk < 4) {
	      while (kk < 4)
		bodydata->qm[f]->qAdj[4*j + kk++] = -1;
	  }
	}
      for (j = 0; j < nquad; j++) {
	  for (q = 0; q < 4; ++q)
	    bodydata->qm[f]->valence[bodydata->qm[f]->qIdx[4 * j + q] - 1][0] = j + 1;
      }
      for (j = 0; j < len; j++) {
	  stat = setValence (bodydata->qm[f], j + 1);
	  if (stat != EGADS_SUCCESS) {
	      printf("In EG_createMeshMap :: set valence at %d is %d!!\n ",
		     j + 1, stat);
	      stat2 = stat;
	      break;
	  }
      }
      for (e4[0]= e4[1] = j = 0; j < len; j++) {
	  if (bodydata->qm[f]->vType[j] == -1) continue;
	  stat = EG_angAtBdVert(bodydata->qm[f], j + 1, e4, &angle);
	  if (stat != EGADS_SUCCESS || angle < EPS08) {
	      stat2   = EGADS_GEOMERR;
	      printf(" FACE %d EG_angAtBdVert %d angle %f\n ",f + 1, stat, angle);
	      printf(" Vertices: %d %d %d \n ", j+ 1, e4[0], e4[1]);
	      break;
	  }
	  else if (angle < 0.85 * PI) bodydata->qm[f]->vType[j] = 2;
	  else if (angle < 1.25 * PI) bodydata->qm[f]->vType[j] = 3;
	  else if (angle < 1.85 * PI) bodydata->qm[f]->vType[j] = 4;
	  else                        bodydata->qm[f]->vType[j] = 5;
	  if (bodydata->qm[f]->valence[j][2] < bodydata->qm[f]->vType[j])
	    bodydata->qm[f]->vType[j] = bodydata->qm[f]->valence[j][2];
	  if (bodydata->qm[f]->vType[j] == 2)
	    bodydata->qm[f]->valence[j][1] = bodydata->qm[f]->valence[j][2] + 2; //adjusted based on vertex type, eg, a corner is regular when val = 2
	  else if (bodydata->qm[f]->vType[j] == 3)
	    bodydata->qm[f]->valence[j][1] = bodydata->qm[f]->valence[j][2] + 1;
	  else if (bodydata->qm[f]->vType[j] >= 5)
	    bodydata->qm[f]->valence[j][1] = bodydata->qm[f]->valence[j][2] - 1;
      }
      if (stat2 != EGADS_SUCCESS) continue;
      // ESTIMATE AREA AND TOTAL QUADS
#ifdef DEBUG
      stat = checkMesh(bodydata->qm[f]);
      if (stat != EGADS_SUCCESS) {
	  printf("In EG_createMeshMap :: checkMesh at face %d -->%d!!\n",
		 f + 1, stat);
	  stat2 = stat;
      }
#endif
  }
  return stat2;
}

static int EG_makeValidMesh(meshMap *qm, int nP, /*@null@*/ int *pList,
			    int fullReg)
{
  int    si, v, q, i, j, k, kv, it = 0, itMax, sum = 0;
  int    recover = 1, stat = EGADS_SUCCESS, *mv = NULL;
  double *uvxyz = NULL, pass, qr;
  int *area = NULL;
#ifdef DEBUG
  char   buffer[100];
  FILE *fout = NULL;
  double pos[18], uv[2];
#endif

  uvxyz = (double*)EG_alloc(5 * qm->totV * sizeof(double));
  mv    = (int   *)EG_alloc(    qm->totV * sizeof(int));
  if (uvxyz == NULL || mv == NULL ) return EGADS_MALLOC;
  area  = (int   *)EG_alloc(qm->totQ * sizeof(int   ));
  for (k = j = 0; j < qm->totV; j++) {
      if (qm->vType[j] != -1) continue;
      uvxyz[5 * j    ] = qm->uvs [2 * j    ];
      uvxyz[5 * j + 1] = qm->uvs [2 * j + 1];
      uvxyz[5 * j + 2] = qm->xyzs[3 * j    ];
      uvxyz[5 * j + 3] = qm->xyzs[3 * j + 1];
      uvxyz[5 * j + 4] = qm->xyzs[3 * j + 2];
  }
  for (j = 0; j < qm->totQ; j++) area[j] = -1;
  if (fullReg == 0) { // move around only affected vertices
      if (nP  == 0 || pList == NULL) goto cleanup;
      for (kv = j = 0; j < nP; j++) {
          si = pList[j] - 1;
          if (qm->vType[pList[j] - 1] == -2) continue;
          if (qm->vType[pList[j] - 1] == -1 &&
              inList (kv, mv, pList[j] - 1) == -1 )
            mv[kv++] = pList[j] - 1;
          if( qm->star[si] == NULL) {
              printf(" STAR %d is NULL !!!\n ", si + 1);
              stat = EGADS_MALLOC;
              goto cleanup;
          }
          for (i = 0; i < qm->star[si]->nQ; i++) {
              if (qm->star[si]->quads[i] == -1) continue;
              v = qm->star[si]->verts[2 * i + 1] - 1;
              if (qm->vType[v] == -1 && inList (kv, mv, v) == -1 ) mv[kv++] = v;
              v = qm->star[si]->verts[2 * i + 2] - 1;
              if (qm->vType[v] == -1 && inList (kv, mv, v) == -1 ) mv[kv++] = v;
          }
      }
      itMax = 5;
  } else {
      for (kv = i = 0 ; i < qm->totV; i++ ) {
          if ( qm->vType[i] != -1 ) continue;
          mv[kv++] = i;
      }
      itMax = 50;
  }
  for (i = 0; i < kv; i++) {
      if (qm->star[mv[i]] == NULL ) {
          printf(" Star in v %d is nULL \n ", i + 1 );
          goto cleanup;
      }
      for (j = 0; j < qm->star[mv[i]]->nQ; j++) {
          k = qm->star[mv[i]]->quads[j] - 1;
          if (qm->qIdx[4 * k] == -2) {
              printf(" Attention!! star %d has assigned a void quad!!\n ", i + 1 );
#ifdef DEBUG
              printQuad(qm, k + 1);
#endif
          }
          if (k < 0) continue;
          area[k] = 0;
      }
  }
#ifdef DEBUG
  printf(" EG_makeValid mesh for %d regular  ---- \n ", fullReg);
  snprintf(buffer, 100,"MKSA_%d", qm->plotcount++);
  printf("Writing in %s\n ", buffer);
  fout = fopen(buffer,"w");
  if (fout != NULL ) {
      for (i = 0 ; i < qm->totQ; i++ ) {
          if (area[i] == -1) continue;
          uv[0] = uv[1] = 0.0;
          for (j = 0; j <= 4; j++) {
              v = qm->qIdx[ 4 * i + j%4 ] - 1;
              fprintf(fout, "%lf %lf %lf %d\n",  qm->xyzs[3*v  ],
                      qm->xyzs[3*v + 1], qm->xyzs[3*v + 2], v + 1 );
              if ( j == 4 ) break;
              uv[0] += 0.25 * qm->uvs[2 * v    ];
              uv[1] += 0.25 * qm->uvs[2 * v + 1];
          }
          fprintf(fout,"\n\n");
          EG_evaluate(qm->face, uv, pos);
          fprintf(fout, "%lf %lf %lf %d\n", pos[0], pos[1], pos[2],
                  i + 1);
          fprintf(fout,"\n\n");
      }
      fclose (fout);
  }
#endif
  for (it = 0 ; it <= itMax; it++) {
      if (it == 0 ) pass = -1.0;
      else {
          if (fullReg < 2) pass = 0.1;
          else             pass = 0.5 - (double)it/(double)itMax * 0.4;
      }
      for (k = 0 ; k < kv; k++)
        EG_placeVertex(qm,  mv[k] + 1, pass);
#ifdef DEBUG
      printf(" MAKEVALID IT %d ===================  \n ", it );
      snprintf(buffer, 100,"MKLA_%d", qm->plotcount++);
      printf("Writing in %s\n ", buffer);
      fout       = fopen(buffer,"w");
      if (fout  != NULL ) {
          for (i = 0 ; i < qm->totQ; i++ ) {
              if (area[i] == -1) continue;
              uv[0]  = uv[1] = 0.0;
              for (j = 0; j <= 4; j++) {
                  v  = qm->qIdx[4 * i + j%4] - 1;
                  fprintf(fout, "%lf %lf %lf %d\n", qm->xyzs[3*v  ],
                          qm->xyzs[3*v + 1], qm->xyzs[3*v + 2], v + 1 );
                  if (j == 4) break;
                  uv[0] += 0.25 * qm->uvs[2 * v    ];
                  uv[1] += 0.25 * qm->uvs[2 * v + 1];
              }
              fprintf(fout,"\n\n");
              EG_evaluate(qm->face, uv, pos);
              fprintf(fout, "%lf %lf %lf %d\n", pos[0], pos[1], pos[2],
                      i + 1);
              fprintf(fout,"\n\n");
          }
          fclose (fout);
      }
#endif
      if (it == 0) continue;
      for (sum = q = 0; q < qm->totQ; q++) {
          if (area[q] == -1) continue;
          area[q] = EG_quadArea(qm, NULL, NULL, q + 1, -1, &qr);
          if (area[q] == QA0 && fullReg == 2 && qr < 0.1) area[q] = QA1;
#ifdef DEBUG
          printf(" QUAD %d AREA %d RATIO %lf\n ", q + 1, area[q], qr);
#endif
          if ((fullReg  >=  2 && area[q] != QA0) ||
                area[q] > QA1) sum = 1;
      }
      if (qm->qInv && fullReg < 2) {
          for (k = 0 ; k < qm->qInv[0]; k++)  {
              if (qm->qInv[1 + k] == -2) continue;
              if (area[qm->qInv[1 + k] -1] <= QA1 &&
                  area[qm->qInv[1 + k] -1] >= QA0) qm->qInv[1 + k] = -2;
          }
      }
      if (sum == 0) {
          recover = 0 ;
          if (fullReg < 2) break;
          for (j = 0; j < qm->totV; j++) {
              if (qm->vType[j] != -1) continue;
              uvxyz[5 * j    ]  = qm->uvs [2 * j    ];
              uvxyz[5 * j + 1]  = qm->uvs [2 * j + 1];
              uvxyz[5 * j + 2]  = qm->xyzs[3 * j    ];
              uvxyz[5 * j + 3]  = qm->xyzs[3 * j + 1];
              uvxyz[5 * j + 4]  = qm->xyzs[3 * j + 2];
          }
      }
  }
  if ((recover == 1 && fullReg <  2) ||
      (recover == 0 && fullReg >= 2  && sum != 0)) {
      if (recover == 1 ) {
          stat = EGADS_GEOMERR;
#ifdef DEBUG
          if (sum < QA2)
            printf(" Mesh has obtuse angles\n");
          else
            printf(" Mesh has invalid quads\n");
          snprintf(buffer, 100,"face_%d_InvalidMesh_%d",
                   qm->fID, qm->invsteps);
          gnuData(qm , buffer);
#endif
          qm->invsteps++;
      }
      for (j = 0; j < qm->totV; j++) {
          if (qm->vType[j] != -1) continue;
          qm->uvs  [2 * j    ] = uvxyz[ 5 * j    ];
          qm->uvs  [2 * j + 1] = uvxyz[ 5 * j + 1];
          qm->xyzs [3 * j    ] = uvxyz[ 5 * j + 2];
          qm->xyzs [3 * j + 1] = uvxyz[ 5 * j + 3];
          qm->xyzs [3 * j + 2] = uvxyz[ 5 * j + 4];
      }
  } else if (sum != 0 && fullReg >= 2) {
      stat = EGADS_GEOMERR;
      gnuData(qm , NULL);
      if(qm->qInv) EG_free(qm->qInv);
      for (sum = q = 0; q < qm->totQ; q++) {
          if (area[q] == -1 || area[q] == QA0)continue;
          sum++;
      }
      qm->qInv = EG_alloc((1 + sum) * sizeof(int));
      if(qm->qInv == NULL) {
          stat = EGADS_MALLOC;
          goto cleanup;
      }
      for (qm->qInv[0] = q = 0; q < qm->totQ; q++) {
          if (area[q] == -1 || area[q] == QA0)continue;
           qm->qInv[++qm->qInv[0]] = q + 1;
      }
#ifdef DEBUG
      for (j = 0 ; j < qm->qInv[0]; j++) printf(" INV %d / %d = %d\n", j + 1, qm->qInv[0], qm->qInv[1+ j]);
#endif
  }
#ifdef DEBUG
  printf(" EG_makeValidMesh stat %d\n", stat);
  gnuData(qm, NULL);
#endif

  cleanup:
  EG_free(uvxyz);
  EG_free(mv);
  EG_free(area);
  return stat;
}


static int EG_createQuadGroup(meshMap *qm, quadGroup *qg, int q0, int q1)
{
  int i, j, k, ids[8], piv = 0, aux, vaux[6], vc[4], c1 = 0, c2 = 0;

  if (q0 < 0 || q1 < 0 || q0 > qm->totQ || q1 > qm->totQ ||
      qm->qIdx[4 * ( q0 - 1 )] == -2 ||
      qm->qIdx[4 * ( q1 - 1 )] == -2 ) {
#ifdef DEBUG
      printf(" EG_commonVerts wrong quad pair %d %d\n", q0, q1);
      printQuad(qm, q0);
      printQuad(qm, q1);
#endif
      return EGADS_INDEXERR;
  }
  qg->q[0] = q0;
  qg->q[1] = q1;
  if (q0 == q1) {
      printf("EG_createQuadGroup :: Quads %d are the same!!!\n", q0);
      return EGADS_INDEXERR;
  }
  for (k = i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
	  if (qm->qIdx[4 * (q0 - 1) + i] == qm->qIdx[4 * (q1 - 1) + j] &&
	      qm->qIdx[4 * (q0 - 1) + i] != -1) {
	      ids[2 *k    ] = i;
	      ids[2 *k + 1] = j;
	      vc [   k++  ] = qm->qIdx[4 * (q1 - 1) + j];
	  }
      }
  }
  if (k != 2) {
      printf("EG_createQuadGroup: Quad %d %d are not adjacent!!\n", q0, q1);
#ifdef DEBUG
      printQuad(qm, q0);
      printQuad(qm, q1);
#endif
      return EGADS_INDEXERR;
  }
  if (qm->valence[vc[1] -1][1] >
  qm->valence[vc[0] -1][1]) {
      ids[0] = ids[2];
      ids[1] = ids[3];
  }
  piv    = 0;
  if (qm->qAdj[4 * (qg->q[0] - 1) + ids[0]] == qg->q[1]) piv = 1;
  qg->verts[0] = qm->qIdx[4 * (qg->q[piv] - 1) +  ids[piv]];
  qg->verts[1] = qm->qIdx[4 * (qg->q[piv] - 1) + (ids[piv] + 1)%4];
  qg->verts[2] = qm->qIdx[4 * (qg->q[piv] - 1) + (ids[piv] + 2)%4];
  qg->verts[3] = qm->qIdx[4 * (qg->q[piv] - 1) + (ids[piv] + 3)%4];
  aux = (piv + 1)%2;
  qg->verts[4] = qm->qIdx[4 * (qg->q[aux] - 1) + (ids[aux] + 2)%4];
  qg->verts[5] = qm->qIdx[4 * (qg->q[aux] - 1) + (ids[aux] + 3)%4];
  if (piv == 1) {
      aux      = qg->q[0];
      qg->q[0] = qg->q[1];
      qg->q[1] = aux;
  }
  if (qm->valence[qg->verts[0] -1][1] ==
      qm->valence[qg->verts[3] -1][1]) {
      if (qm->valence[qg->verts[1] -1][1] != 4) c1++;
      if (qm->valence[qg->verts[5] -1][1] != 4) c1++;
      if (qm->valence[qg->verts[2] -1][1] != 4) c2++;
      if (qm->valence[qg->verts[4] -1][1] != 4) c2++;
      if (c2 > c1) {
	  vaux[0]      = qg->verts[0];
	  vaux[1]      = qg->verts[1];
	  vaux[2]      = qg->verts[2];
	  qg->verts[0] = qg->verts[3];
	  qg->verts[1] = qg->verts[4];
	  qg->verts[2] = qg->verts[5];
	  qg->verts[3] = vaux[0];
	  qg->verts[4] = vaux[1];
	  qg->verts[5] = vaux[2];
	  aux          = qg->q[0];
	  qg->q[0]    = qg->q[1];
	  qg->q[1]    = aux;
      }
  }
  qg->vals[0] = qm->valence[qg->verts[0] -1][1];
  qg->vals[1] = qm->valence[qg->verts[1] -1][1];
  qg->vals[2] = qm->valence[qg->verts[2] -1][1];
  qg->vals[3] = qm->valence[qg->verts[3] -1][1];
  qg->vals[4] = qm->valence[qg->verts[4] -1][1];
  qg->vals[5] = qm->valence[qg->verts[5] -1][1];
  return EGADS_SUCCESS;
}

/* Assuming qID is collapsing through v */
static int EG_validCollapse(meshMap *qm, int qID, int v)
{
  int j, k, kk, id, link[2], aux, aux2, type, val;

  if (v <= 0 || v > qm->totV) return EGADS_INDEXERR;
  id     = EG_quadVertIdx(qm, qID, v);
  aux    = qm->qIdx[4 * (qID - 1) + (id + 2)%4] - 1;
  if (qm->vType[aux] != -1) {
      for (j = 0; j < qm->valence[v - 1][2]; j++) {
	  k  = qm->valence[v - 1][3 + j] -1;
	  if (qm->vType[k] == -1) continue;
	  if (EG_quadVertIdx(qm, qID, k + 1) < 0) {
#ifdef DEBUG
	      printf(" Vertex %d is bounds and not in quad %d. "
		  "This will produce a new artificial edge, ie boundary vertices from "
		  "different edges will be linked. Not at all desirable \n ", k + 1, qID );
#endif
	      return 0;
	  }
      }
  }
  link[0]  = qm->qIdx[4 * (qID - 1) + (id + 1)%4];
  link[1]  = qm->qIdx[4 * (qID - 1) + (id + 3)%4];
  for (j   = 0; j < 2; j++) {
      type = qm->vType  [link[j] - 1];
      val  = qm->valence[link[j] - 1][2];
      if (type > 0 ) {
	  if ( val <= type) {
#ifdef DEBUG
	      printf(" link %d has valence %d type %d !\n ", link[j], val, type);
#endif
	      return 0;
	  }else continue;
      }
      else if (val != 3) continue;
      if (qm->vType[link[(j+1)%2] -1] >=4 ) return 0;
      for (k  = 0; k < 3; k++) {
	  aux = qm->valence[link[j] - 1][3 + k] - 1;
	  if (EG_quadVertIdx(qm, qID,  aux + 1) >= 0) continue;
	  if (qm->valence[aux][2] == qm->vType[aux]) {
#ifdef DEBUG
	      printf(" Chain Effect Doublets %d -->%d %d \n",
		     aux + 1, qm->valence[aux][2], qm->vType[aux]);
#endif
	      return 0;
	  }
	  if (qm->valence[aux][2] != 3) continue;
	  for (kk = 0; kk < 3; kk++) {
	      aux2 = qm->valence[aux][3 + kk] - 1;
	      if (aux2 + 1 == link[j]) continue;
	      if (qm->vType[aux] >= 4 ||
		  qm->vType[aux] == 2 ||
		   qm->valence[aux2][2] == qm->vType[aux] ||
		  (qm->valence[aux2][2] == 3 && qm->vType[aux] == -1)) {
#ifdef DEBUG
		  printf(" KK Chain Effect Doublets %d -->%d %d \n",
			 aux2 + 1, qm->valence[aux][2], qm->vType[aux]);
#endif
		  return 0;
	      }
	  }
      }
  }
  return 1;
}


/* Assuming we will break link v1-v2 */
static int EG_validSwap(meshMap *qm, int v1, int v2)
{
  int i, vs[2], type, val;

  vs[0]    = v1 - 1; vs[1] = v2 - 1;
  for (i   = 0; i < 2; i++) {
      type = qm->vType  [vs[i]];
      val  = qm->valence[vs[i]][2];
      if ((type > 0 && val <= type) || (type == -1 && val <= 3)) return 0;
  }
  return 1;
}


static int EG_swappingOperation(meshMap *qm, quadGroup qg, int swap,
				int *activity)
{
  int   nq, stat, i0, i1, i, j, k, adj;
  int   *list, qID[2], adjQmap[6];
  double ratio;
  Quad *quad = NULL;

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
  qm->qIdx[4 * (qID[0] - 1)] = qg.verts[swap];
  qm->qIdx[4 * (qID[1] - 1)] = qg.verts[swap];
  for (i = 0; i < 3; i++) {
      adjQmap[i    ] = qm->qAdj[4 * (qID[0] - 1) + (i + i0)%4];
      adjQmap[i + 3] = qm->qAdj[4 * (qID[1] - 1) + (i + i1)%4];
      qm->qIdx[4 * (qID[0] - 1) + i + 1] = qg.verts[(swap + 3 + i)%6];
      qm->qIdx[4 * (qID[1] - 1) + i + 1] = qg.verts[(swap + 1 + i)  ];
  }
  qm->qAdj[4 * (qID[0] - 1)    ] = qID[1];
  qm->qAdj[4 * (qID[1] - 1) + 3] = qID[0];
  for (i = 0; i < 3; i++) {
      adj = adjQmap[(3 + i + swap)%6] - 1;
      qm->qAdj[4 * (qID[0] - 1) + i + 1] = adj + 1;
      if (adj >= 0) {
	  for (j = 0; j < 4; j++) if (qm->qAdj[4 * adj + j] == qID[1])
	    qm->qAdj[4 * adj + j] = qID[0];
      }
      adj = adjQmap[i + swap] -1;
      qm->qAdj[4 * (qID[1] - 1) + i] = adj + 1;
      if (adj >= 0) {
	  for (j = 0; j < 4; j++) if (qm->qAdj[4 * adj + j] == qID[0])
	    qm->qAdj[4 * adj + j] = qID[1];
      }
  }
  for (i = 0; i < 4; i++) {
      j = qm->qIdx[4 * (qID[0] - 1) + i] - 1;
      qm->valence[j][0] = qID[0];
      j = qm->qIdx[4 * (qID[1] - 1) + i] - 1;
      qm->valence[j][0] = qID[1];
  }

  list = EG_alloc (qm->totV * sizeof (int ));
  for (i0 = i = 0; i < 4; i++) {
      if ( inList(i0, list, qg.verts[i] ) == -1 ) list[i0++] = qg.verts[i];
      for (j = 1; j < qm->star[qg.verts[i] - 1]->nV; j++) {
	  k  = qm->star[qg.verts[i] - 1]->verts[j];
	  if ( k == -1) continue;
	  if (inList(i0, list, k ) == -1 ) list[i0++] = k;
      }
  }
  for ( i = 0 ; i < i0; i++ ) {
      stat = setValence(qm, list[i]);
      if (stat != EGADS_SUCCESS) {
	  fprintf(stderr," ERROR SETTING VALENCE %d \n ", list[i]);
	  EG_free(list);
	  return stat;
      }
  }
  EG_free(list);
  *activity = 1;
  i = EG_quadArea(qm, NULL, NULL, qg.q[0],-1, &ratio);
  j = EG_quadArea(qm, NULL, NULL, qg.q[1],-1, &ratio);
  if (i != QACB && j != QACB) {
      if (EG_makeValidMesh(qm, 6, qg.verts, 0) == EGADS_SUCCESS) {
	  EG_free(quad);
	  return EGADS_SUCCESS;
      }
  }
  else {
#ifdef DEBUG
      printf(" SWAP CROSSES BOUND %d  %d  \n ", i, j);
      gnuData(qm, NULL);
#endif
  }
#ifdef DEBUG
  printf(" EG_swapping operation didn't work well: restore quads\n ");
#endif
  stat = EG_restoreQuads(qm, quad, nq);
  EG_free(quad);
  if (stat != EGADS_SUCCESS) {
      printf(" EG_swappingOperation :: restore quads %d !!\n ", stat);
      return stat;
  }
  *activity = 0;
  return EGADS_SUCCESS;
}


static int EG_splittingOperation(meshMap *qm, int vC, int vL, int vR,
				 int *activity)
{
  int   qIdx[4], modQ[4], verts[4], adj[2], poly[4], q, newQ, i, j, stat;
  int   id0 = -1, id1 = -1, dist, links[4], vals[4], addedV = 0, nq, si, *list = NULL, n;
  double uv[2];

  Quad  *quad = NULL;

  si          = vC - 1;
  *activity   = 0;
  if (qm->star[si] == NULL ) {
      printf(" In splittingOperation build star %d is NULL !!\n", vC);
      return EGADS_MALLOC;
  }
  nq        = qm->star[si]->nQ;
  stat      = EG_backupQuads(qm, &nq, qm->star[si]->quads, &quad);
  if (stat != EGADS_SUCCESS || quad == NULL) {
      EG_free(quad);
      printf("In EG_splittingOperation :: EG_backupQuads %d !!\n", stat);
      return stat;
  }
  if (qm->remQ[0] > 0) {
      poly[3] = qm->remV[qm->remV[0]--];
      newQ    = qm->remQ[qm->remQ[0]--];
  } else {
      poly[3] = qm->totV + 1;
      newQ    = qm->totQ + 1;
      if (poly[3] >  qm->sizeV) {
	  printf(" We have duplicated the number of initial vertices. This is too much. \n");
	  EG_free(quad);
	  return EGADS_INDEXERR;
      }
      ++qm->totV; ++qm->totQ;
      addedV = 1;
  }
  qm->vType[poly[3] - 1] = -1;
  if (qm->vType[vC - 1] != -1) {
      for (q = 0; q < qm->star[si]->nQ; q++)
	if (qm->star[si]->quads[q] == -1) break;
      for (i = 0; i < 2; i++) {
	  if (i == 1) dist = (qm->star[si]->nV - 1) - 4;
	  else        dist = 4;
	  id0 = qm->star[si]->idxQ[q + i];
	  links[2 * i    ] = qm->star[si]->verts[qm->star[si]->idxV[2 * id0 + 1       ]];
	  links[2 * i + 1] = qm->star[si]->verts[qm->star[si]->idxV[2 * id0 + 1 + dist]];
	  vals [2 * i    ] = qm->valence [links[2*i    ] - 1][1];
	  vals [2 * i + 1] = qm->valence [links[2*i + 1] - 1][1];
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
  for (j = 0; j < qm->star[si]->nQ; j++) {
      if (qm->star[si]->verts[2 * j + 1] == vL) id0 = j;
      if (qm->star[si]->verts[2 * j + 1] == vR) id1 = j;
  }
  if (id0 == -1 || id1 == -1) {
#ifdef DEBUG
      printf("EG_splittingOperation  Problem finding vertices %d !!\n ",
	     EGADS_INDEXERR);
      printf(" vL %d vR %d \n ", vL, vR);
      for (j = 0; j < qm->star[si]->nV; j++)
	printf(" j %d =%d\n ", j, qm->star[si]->verts[j]);
#endif
      EG_free(quad);
      return EGADS_INDEXERR;
  }
  poly [0] = qm->star[si]->verts[0];
  poly [1] = qm->star[si]->verts[2 * id0 + 1];
  poly [2] = qm->star[si]->verts[2 * id1 + 1];
  qIdx [0] = id0;
  qIdx [1] = qm->star[si]->idxQ[id1 + qm->star[si]->nQ - 1];
  qIdx [2] = id1;
  qIdx [3] = qm->star[si]->idxQ[id0 + qm->star[si]->nQ - 1];
  verts[0] = poly[1];
  verts[1] = poly[2];
  verts[2] = poly[2];
  verts[3] = poly[1];
  qm->qIdx[4 * (newQ - 1)    ] = poly[1];
  qm->qIdx[4 * (newQ - 1) + 1] = poly[0];
  qm->qIdx[4 * (newQ - 1) + 2] = poly[2];
  qm->qIdx[4 * (newQ - 1) + 3] = poly[3];
  for (i = 0; i < 4; ++i) {
      modQ[i] = qm->star[si]->quads[qIdx[i]];
      qm->qAdj[4 * (newQ - 1) + i] = modQ[i];
      qm->valence[qm->qIdx[4*(newQ - 1)+i]-1][0] = newQ;
      if (modQ[i] == -1) continue;
      stat = EG_adjQtoPair(qm, modQ[i], poly[0], verts[i], adj);
      qm->qAdj[4 * (modQ[i] - 1) + adj[0]] = newQ;
  }
  j = qIdx[2];
  q = qm->star[si]->quads[j++];
  while (q != qm->star[si]->quads[qIdx[0]]) {
      for (i = 0; i < 4; ++i)
	if (qm->qIdx[4 * (q - 1) + i] == poly[0])
	  qm->qIdx[4 * (q - 1) + i] = poly[3];
      q = qm->star[si]->quads[qm->star[si]->idxQ[j++]];
      if (q == -1) {
	  if (qm->star[si]->quads[qIdx[0]] == -1) break;
	  else q = qm->star[si]->quads[qm->star[si]->idxQ[j++]];
      }
  }
  // update stars
  list = EG_alloc(qm->totV * sizeof ( int ));
  if (list == NULL) {
      EG_free(quad);
      return EGADS_MALLOC;
  }
  setValence(qm, poly[3]);
  for (n = i = 0; i < 4; i++) {
      if (inList(n, list, poly[i]) == -1) list[n++] = poly[i];
      for (j = 1; j < qm->star[poly[i] - 1]->nV; j++) {
	  q  = qm->star[poly[i] - 1]->verts[j];
	  if ( q == -1 ) continue;
	  if (inList(n, list, q ) == -1 ) list[n++] = q;
      }
  }
  for ( i = 0 ; i < n; i++ ) {
      stat = setValence(qm, list[i]);
      if (stat != EGADS_SUCCESS) {
	  fprintf(stderr," ERROR SETTING VALENCE %d \n ", poly[i]);
	  EG_free(quad);
	  EG_free(list);
	  return stat;
      }
  }
  EG_free(list);
  si = poly[3] - 1;
  if (qm->star[si] == NULL) {
      fprintf(stderr," star %d ois NULL\n", poly[3]);
      return EGADS_MALLOC;
  }
  id0    = -1;
  stat = EG_centroid(qm, qm->star[si]->nQ,
		       &qm->valence[poly[3] -1][3], uv, 0, 0);
  updateVertex  (qm, poly[3], uv);
  EG_placeVertex(qm, poly[3], 0.1);
  *activity = 1;
#ifdef DEBUG
  printf(" SPLIT THRU %d\n ", poly[0] );
  gnuData(qm, NULL);
#endif

  if (EG_makeValidMesh(qm, 4, poly, 0) == EGADS_SUCCESS) {
      EG_free(quad);
      return EGADS_SUCCESS;
  }
#ifdef DEBUG
  printf("EG_splittingOperation mesh valid failed: restore quads\n");
#endif
  stat = EG_restoreQuads(qm, quad, nq);
  EG_free(quad);
  if (stat != EGADS_SUCCESS) {
      printf(" EG_splittingOperation :: EG_restoreQuads %d !!\n ", stat);
      return stat;
  }
  if (addedV == 1) {
      --qm->totV;
      --qm->totQ;
  } else {
      qm->remV[++qm->remV[0]] = poly[3];
      qm->remQ[++qm->remQ[0]] = newQ;
  }
  qm->vType[poly[3] - 1] = -2;
  for (i = 0; i < 4; i++)
    qm->qIdx[4 * (newQ - 1) + i] = -2;
  *activity = 0;
#ifdef DEBUG
  gnuData(qm, NULL);
#endif
  return stat;
}


static int EG_mergeVertices(meshMap *qm, int qC, int centre, int *activity)
{
  int    stat, i, j, q, adjq, adjPair[2], auxQ, oldQ[8], nq, doublet = 0;
  int    piv[4] = {1, 0, 3, 2}, n, *list = NULL;
  double uv[2], uvxyz[10];
  Quad   *quad = NULL;

  nq        = 1;
  stat      = EG_backupQuads(qm, &nq, &qC, &quad);
  stat     += EG_centroid(qm, 4, &qm->qIdx[4 * (qC - 1)], uv, 0, 0);
  if (stat != EGADS_SUCCESS || quad == NULL) {
      printf(" In EG_mergeVertices stat in EG_backpQuads %d !!\n ", stat);
      return stat;
  }
  *activity  = 0;
  i          = EG_quadVertIdx(qm, qC, centre);
  if (qm->vType[centre - 1] != -1) {
      i = (i+ 2)%4;
      centre = qm->qIdx[4 * (qC - 1) + i];
  }
#ifdef DEBUG
  printf(" EG_mergeVertices collapsing %d through %d = %d \n ",
	 qC, centre, i);
  gnuData(qm, NULL);
#endif

  if (qm->vType[centre - 1] != -1 || i < 0) {
      printf(" EG_mergeVertices collapse through %d is type %d and id in quad %d->%d !!\n",
	     centre, qm->vType[centre - 1], qC, i);
#ifdef DEBUG
      printQuad(qm, qC);
#endif
      return EGADS_INDEXERR;
  }
  if (qm->valence[centre -1][2] == 2) doublet = 1;
  for (q = 0; q < 4; ++q) {
      oldQ[q    ] = qm->qIdx[4 * (qC - 1) + (q + i)%4];
      oldQ[q + 4] = qm->qAdj[4 * (qC - 1) + (q + i)%4];
  }
  uvxyz[0] = qm->uvs [2 * (oldQ[0] - 1)    ];
  uvxyz[1] = qm->uvs [2 * (oldQ[0] - 1) + 1];
  uvxyz[2] = qm->xyzs[3 * (oldQ[0] - 1)    ];
  uvxyz[3] = qm->xyzs[3 * (oldQ[0] - 1) + 1];
  uvxyz[4] = qm->xyzs[3 * (oldQ[0] - 1) + 2];
  uvxyz[5] = qm->uvs [2 * (oldQ[2] - 1)    ];
  uvxyz[6] = qm->uvs [2 * (oldQ[2] - 1) + 1];
  uvxyz[7] = qm->xyzs[3 * (oldQ[2] - 1)    ];
  uvxyz[8] = qm->xyzs[3 * (oldQ[2] - 1) + 1];
  uvxyz[9] = qm->xyzs[3 * (oldQ[2] - 1) + 2];
  // save old quads vertices and adjacents to update map correctly
  for (i = 0; i < 4; i++) {
      q = oldQ[4 + i];
      if (q < 0) continue;
      adjq = oldQ[4 + piv[i]]; // - 1 bias
      stat = EG_adjQtoPair(qm, q, oldQ[i], oldQ[(i + 1)%4], adjPair);
      if (stat != EGADS_SUCCESS || adjPair[1] != qC) {
	  printf("EG_mergeVertices stat %d  adj %d to pair %d %d is %d \n ",
		 stat, q ,oldQ[i], oldQ[(i + 1)%4], adjPair[1]);
	  return  EGADS_INDEXERR;
      }
      qm->qAdj[4 * (q - 1) + adjPair[0]] = adjq;
  }
  // Eliminate vertex oldQ[0] from all the quads and its valences
  for (i = 0; i < qm->star[oldQ[0] - 1]->nQ; ++i) {
      q = qm->star[oldQ[0] - 1]->quads[i];
      if      (q == -1) continue; // ghost quad
      else if (q == qC) {  // eliminate quad qC
	  for (j = 0; j < 4; ++ j) {
	      qm->qIdx[4 * (qC - 1) + j] = -2;
	      qm->qAdj[4 * (qC - 1) + j] = -2;
	  }
      } else {
	  for (j = 0; j < 4; ++ j)
	    if (qm->qIdx[4 * (q - 1) + j] == oldQ[0])
	      qm->qIdx[4 * (q - 1) + j] = oldQ[2];
      }
  }
  // Point all the collapsed quad vertices to a valid quad
  for (i = 1; i < 4; i++) {
      if (qm->valence[oldQ[i] - 1][0] != qC) continue;
      for (q = 0; q < 4; q++) {
	  auxQ = oldQ[4 + q];
	  if (auxQ < 0) continue;
	  if (EG_quadVertIdx(qm, auxQ, oldQ[i]) >= 0) {
	      qm->valence[oldQ[i] - 1][0] = auxQ;
	      break;
	  }
      }
      if (qm->valence[oldQ[i] - 1][0] == qC) {
	  printf(" EG_mergeVertices opposite vertex still points at collapsed quad %d !!\n",
		 EGADS_INDEXERR);
	  EG_free(list);
	  return EGADS_INDEXERR;
      }
  }
  setValence(qm, oldQ[2]);
  setValence(qm, oldQ[1]);
  setValence(qm, oldQ[3]);
  list = EG_alloc(qm->totV * sizeof ( int ));
  if (list == NULL) {
      EG_free(quad);
      return EGADS_MALLOC;
  }
  n = 0;
  for (i = 1; i < 4; i++) {
      for (j = 1; j < qm->star[oldQ[i] - 1]->nV; j++) {
	  q  = qm->star[oldQ[i] - 1]->verts[j];
	  if ( q == -1 ) continue;
	  if (inList(n, list, q ) == -1 ) list[n++] = q;
      }
  }
  for ( i = 0 ; i < n; i++ ) {
      stat = setValence(qm, list[i]);
      if (stat != EGADS_SUCCESS) {
	  printf(" EG_mergeVertices error setting valence %d \n ", list[i]);
	  EG_free(quad);
	  EG_free(list);
	  return stat;
      }
  }
  EG_free(list);
  if (qm->vType[oldQ[3] - 1] >= 4 ||
      qm->vType[oldQ[1] - 1] >= 4 ) doublet = 0;
  if (qm->vType[oldQ[2] - 1] == -1 && doublet == 0 ) {
      updateVertex(qm, oldQ[2], uv);
      EG_placeVertex(qm, oldQ[2], 0.1);
  }
  qm->vType[oldQ[0] - 1]  = -2; // -2 = removed
  qm->remQ[++qm->remQ[0]] = qC;
  qm->remV[++qm->remV[0]] = oldQ[0];
  *activity  = 1;
#ifdef DEBUG
  printf(" BEFORE CALLING MAKEVALID \n ");
  gnuData(qm, NULL);
#endif

  if ( doublet == 1 ) {
#ifdef DEBUG
      printf(" COLLAPSING DOUBLET\n ");
#endif
      EG_free(quad);
      return EGADS_SUCCESS;
  }
  if (EG_makeValidMesh(qm, 3, &oldQ[1], 0) == EGADS_SUCCESS) {
      EG_free(quad);
      return EGADS_SUCCESS;
  }
#ifdef DEBUG
  printf("EG_mergeVertices collapse quad %d thru vertex %d -->%d !!\n",
	 qC, oldQ[0], stat);
#endif
  qm->vType[oldQ[0] - 1]  = -1;
  stat = EG_restoreQuads(qm, quad, nq);
  EG_free(quad);
  if (stat != EGADS_SUCCESS) {
      printf("EG_mergeVertices :: EG_restoreQuads stat %d \n", stat);
      return stat;
  }
  i = -1;
  if (qm->remQ[0] > 0)
    i = inList(qm->remQ[0], &qm->remQ[1], qC);
  if (i != -1) {
      for (j = i + 1; j < qm->remQ[0] - 1; j++) {
	  qm->remQ[j] = qm->remQ[j + 1];
	  qm->remV[j] = qm->remV[j + 1];
      }
  }
  --qm->remQ[0];
  --qm->remV[0];
  for (i = 0; i < 3; i++) {
      if (i < 2)  {
	  qm->uvs[2 * (oldQ[0] - 1) + i] = uvxyz[    i];
	  qm->uvs[2 * (oldQ[2] - 1) + i] = uvxyz[5 + i];
      }
      qm->xyzs[3 * (oldQ[0] - 1) + i] = uvxyz[2 + i];
      qm->xyzs[3 * (oldQ[2] - 1) + i] = uvxyz[7 + i];
  }
  *activity = 0;
  return EGADS_SUCCESS;
}


static int EG_swap(meshMap *qm, int qIn, int *activity)
{
  int       stat, q, swap = 0;
  quadGroup qg;

  *activity  = 0;
  if (qm->qIdx[4 * (qIn - 1)] == -2) {
#ifdef DEBUG
      printf(" EG_swap for quad %d is %d \n ", qIn, stat);
#endif
      return EGADS_SUCCESS;
  }
  for (swap = q = 0; q < 4; q++) {
      qg.q[0]      = qIn;
      qg.q[1]      = qm->qAdj[4 * (qIn - 1) + q];
      if (qg.q[1] == -1) continue;
      stat         = EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
      if (stat    != EGADS_SUCCESS) {
#ifdef DEBUG
	  printf(" EG_swap->EG_createQuadGroup  %d !! \n", stat);
	  printQuadGroup(qm, qg);
#endif
	  return stat;
      }
      if      (qg.vals[0] <= 4 ||
	  EG_validSwap(qm, qg.verts[0], qg.verts[3]) != 1) continue;
      if      (qg.vals[1] * qg.vals[4] == 9) swap = 1;
      else if (qg.vals[2] * qg.vals[5] == 9) swap = 2;
      else if (qg.vals[3] >= 5 ||
	  (qm->vType[qg.verts[0] - 1] ==
	      2 * qm->sizeQ && qg.vals[0] > 5) ||
	      (qm->vType[qg.verts[3] - 1] ==
		  2 * qm->sizeQ && qg.vals[3] > 5)) {
	  if  (    qg.vals[1] * qg.vals[4] == 12) swap = 1;
	  else if (qg.vals[2] * qg.vals[5] == 12) swap = 2;
      }
      if (swap != 0) break;
  }
  if (swap  == 0) return EGADS_SUCCESS;
  return EG_swappingOperation(qm, qg, swap, &(*activity));
}


static int EG_doubleSwap(meshMap *qm, quadGroup qg, int forcing, int *activity)
{
  int piv5 = -1, piv3 = -1, q5, i, adjPiv5, stat, adj[2], swap = 0, vopp3;

  *activity = 0;
  if      (EG_validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
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
  if (qm->valence[qg.verts[vopp3] -1][1] != 4) return EGADS_SUCCESS;
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
	printf(" EG_doubleSwap adjToPair %d -->%d !!\n", adj[1], stat);
#endif
      return stat;
  }
  i       = EG_quadVertIdx(qm, adj[1], qg.verts[vopp3]);
  adjPiv5 = qm->qIdx[4 * (adj[1] - 1) + (i + 1)%4];
  if (adjPiv5 == qg.verts[piv5])
    adjPiv5 = qm->qIdx[4 * (adj[1] - 1) + (i + 3)%4];
  if (      swap == 0 && qm->valence[adjPiv5 -1][1] > 4) return EGADS_SUCCESS;
  else if ((swap == 1 && qm->valence[adjPiv5 -1][1] < 5) ||
      (qm->vType[adjPiv5 - 1] != -1 &&
	  qm->vType[qg.verts[vopp3] - 1] != -1)) return EGADS_SUCCESS;
  piv5      = qg.verts[0];
  stat      = EG_swappingOperation(qm, qg, piv3, &(*activity));
  if (stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
      printf(" EG_doubleSwap: at first swap activity %d, stat -->%d !!\n ",
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


static int EG_collapse(meshMap *qm, int qID, int *activity, int forcing)
{
  int i, vC, vO, v5[5], v3[5], qb[4], stat, i3 = 0, val[2], links[3], act = 0;

  *activity = 0;
  if (qm->qIdx[4 * (qID - 1)] == -2) {
#ifdef DEBUG
      printf(" EG_collapse for quad %d is empty \n ", qID);
#endif
      return EGADS_SUCCESS;
  }
  for (i = 0; i < 4; i++) qb[i] = qm->qIdx[4 * (qID - 1) + i];
  for (v3[0] = v5[0] = i = 0; i < 4; i++) {
      val[0]   = qm->valence[qb[i] -1][1];
      if      (val[0] == 3) v3[++v3[0]] = i;
      else if (val[0] >= 5) v5[++v5[0]] = i;
      else if (val[0] == 2) {
	  stat = EG_mergeVertices(qm, qID, qb[i], &(*activity));
	  if (stat != EGADS_SUCCESS || *activity == 0) {
	      printf("EG_collapse :: I can't fail whilst removing doublets !!!!!!!! \n ");
#ifdef DEBUG
	      gnuData(qm, NULL);
	      printQuad (qm, qID );
#endif
	      return EGADS_GEOMERR;
	  }
	  break;
      }
  }
  if (forcing == 0 && *activity == 0) {
      if (v3[0] + v5[0] <= 2 || v3[0] * v5[0] == 0) return EGADS_SUCCESS;
      if (v5[0] == 1 || (v5[0] == 2 && abs(v5[1] - v5[2])%2 == 0)) {
	  vC   = qb[ (v5[1] + 1)%4];
	  if (EG_validCollapse(qm, qID, vC) == 1) {
	      stat      = EG_mergeVertices(qm, qID, vC, &(*activity));
	      if (stat != EGADS_SUCCESS) {
		  printf(" EG_collapse :: EG_mergeVertices %d activity %d \n",
			 stat, *activity);
		  return stat;
	      }
	  }
      }
  }
  if (*activity == 0 && forcing > 0) {
      if (EG_nValenceCount(qm, qID, 3) == 0 && forcing == 1) return EGADS_SUCCESS;
      i = 0; i3 = 0;
      while (i < 2) {
	  vC = qm->qIdx[4 * (qID - 1) +  i3       ];
	  vO = qm->qIdx[4 * (qID - 1) + (i3 + 2)%4];
#ifdef DEBUG
	  printQuad(qm, qID ) ;
	  printf(" PAIRS %d -->OPP %d VAL %d \n ", vC, vO, qm->valence[vO - 1][2]);
#endif
	  if (qm->vType  [vC - 1]  == -1 && ( forcing == 3 ||
	      (qm->valence[vO - 1][2] <  10 && qm->valence[vC - 1][2] ==  3)) )
	    {
	      links[0] = qm->qIdx[4 * (qID - 1) + (i3 + 1)%4];
	      links[1] = qm->qIdx[4 * (qID - 1) + (i3 + 3)%4];
	      val  [0] = qm->valence[links[0] - 1][1];
	      val  [1] = qm->valence[links[1] - 1][1];
	      if ((qm->valence[links[0] - 1][1] >= 5 ||
		  qm->valence[links[1] - 1][1] >= 5 || i == 1) &&
		  EG_validCollapse(qm, qID, vC) == 1) {
		  stat      = EG_mergeVertices(qm, qID, vC, &(*activity));
		  // printf(" MERGE QUAD %d stat %d activity %d \n ", qID, stat, *activity );
		  if (stat != EGADS_SUCCESS) {
		      printf("EG_collapse :: forcing EG_mergeVertices quad %d !! \n",stat);
#ifdef DEBUG
		      printQuad(qm, qID);
#endif
		      return stat;
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
  if (stat  == EGADS_SUCCESS && *activity == 0) return stat;
  for (i = 0; i < 4; i++) {
      links[0] = qb[i] - 1;
      if (qm->vType[links[0]] == -2) continue;
      links[1] = qm->valence[links[0]][3] - 1;
      links[2] = qm->valence[links[0]][4] - 1;
      for (i3  = 0; i3 < 3; i3++) {
	  if (qm->vType[links[i3]] != -1 ||
	      qm->valence[links[i3]][2] != 2) continue;
	  stat = EG_collapse(qm, qm->valence[links[i3]][0], &act, 0);
	  if (stat != EGADS_SUCCESS || act == 0) {
	      printf(" EG_collapse removing doublet in doublet loop : V %d stat %d act %d !!!!\n ",
		     links[i3] + 1, stat, act );
	      return EGADS_GEOMERR;
	  }
      }
  }
  return EGADS_SUCCESS;
}


static int EG_swapCollapse(meshMap *qm, quadGroup qg, int forcing, int *activity)
{
  int  stat, i, i3 = -1, q5, qC, vL5, vL5adj, vOpp3, swap = 0, adj[2], v0;

  *activity = 0;
  if (EG_validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
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
  if (EG_validCollapse(qm, qg.q[q5], qg.verts[i3]) == 0) return EGADS_SUCCESS;
  vOpp3 = (i3 + 3)%6;
  qC    = qg.q[q5];
  vOpp3 = (i3 + 3)%6;
  vL5   = (vOpp3 + 1)%6;
  if (vL5%3 != 0) vL5 = (vOpp3 + 5)%6;
  if (qg.vals[vOpp3] == 3) return EGADS_SUCCESS;
  stat = EG_adjQtoPair(qm, qg.q[q5], qg.verts[vOpp3], qg.verts[vL5], adj);
  if (stat != EGADS_SUCCESS || adj[1] == -1) return stat;
  i      = EG_quadVertIdx(qm, adj[1], qg.verts[vOpp3]);
  vL5adj = qm->qIdx[4 * (adj[1] - 1) + (i + 1)%4];
  if (vL5adj == qg.verts[vL5]) vL5adj = qm->qIdx[4 * (adj[1] - 1) + (i + 3)%4];
  if ((forcing == 0 && qg.vals[vL5] == 4 && qm->valence[vL5adj - 1][1] >= 4) ||
      (forcing == 1 && qm->valence[vL5adj - 1][1] > 4)) return EGADS_SUCCESS;
  stat      = EG_createQuadGroup(qm, &qg, qg.q[q5], adj[1]);
  if (stat != EGADS_SUCCESS) {
#ifdef DEBUG
      printf("EG_swapCollapse before swap: EG_createQuadGroup -->%d !!\n", stat);
      printQuadGroup(qm, qg);
#endif
      return stat;
  }
  if (EG_validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
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
      if (qg.vals[i5] < 5) i5 = (i3 + 5)%6;
      if (qg.vals[i5] < 5 || qg.vals[(i5 + 3)%6] == 4) return EGADS_SUCCESS;
  }
  else return EGADS_SUCCESS;

  q3 = 0;
  if (forcing == 1 && qg.vals[0] * qg.vals[3] == 16) {
      if (i3 >= 3) q3 = 1;
  } else {
      if (i5 >= 3) q3 = 1;
  }
  for (i = 0; i < 2; i++) {
      if (qm->qIdx[4 * (qg.q[(q3 + i)%2] -1)] == -2) continue;
      stat        = EG_collapse(qm, qg.q[(q3 + i)%2], &i3, 1);
      if (stat   == EGADS_SUCCESS) {
	  if (i3 == 0) break;
	  (*activity)++;
      }
      if (stat != EGADS_SUCCESS) {
	  printf(" EG_doubleCollapse after EG_forceColapse quad %d->%d !!\n ",
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
      EG_validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
  if (qg.vals[swap] * qg.vals[(swap+3)%6] != 12) swap = 2;
  if (qg.vals[swap] * qg.vals[(swap+3)%6] != 12) return EGADS_SUCCESS;
  v0        = qg.verts[0];
  stat      = EG_swappingOperation(qm, qg, swap, &(*activity));
  if (stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
      printf(" In EG_swapDoubleCollapse: Swapping went %d activity %d \n",
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
#ifdef DEBUG
      printQuadGroup(qm, qg);
#endif
      return EGADS_INDEXERR;
  }
  stat      = EG_collapse(qm, qg.q[k], &(*activity), 1);
  if (stat == EGADS_SUCCESS && *activity > 0) return stat;
  if (*activity == 0 ) {
      j = qg.q[k] - 1;
      if ( qm->qIdx[4 * j] == -2 ) {
	  printf(" QUAD IS VOID BUT ACTIVITY 0 \n ");
	  return EGADS_INDEXERR;
      }
  }
  stat      = EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
  if (stat != EGADS_SUCCESS) {
      printf("EG_swapDoubleCollapse after invalid force collapse :: EG_createQuadGroup is %d \n ",
	     stat);
#ifdef DEBUG
      printQuadGroup(qm, qg);
#endif
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


static int EG_split(meshMap *qm, int qID, int *activity)
{
  int si, v, q, id0 = 0, i, stat, dist = 0, validSplit = 0;
  int poly[3], val[3], id6[2], links[4];

  *activity   = 0;
  if (qm->qIdx[4 * (qID - 1)] == -2) {
#ifdef DEBUG
      printf(" EG_split for quad %d is %d \n ", qID, stat);
#endif
      return EGADS_SUCCESS;
  }
  for (v = 0; v < 4; v++) {
      poly[0]   = qm->qIdx[4 * (qID - 1) + v];
      val [0]   = qm->valence[poly[0] - 1][1];
      if (val[0] < 5 || qm->vType[poly[0] - 1] == 4) continue;
      si = poly[0] - 1;
      if (qm->star[si] == NULL ) {
	  printf("EG_split star at vertex %d is NULL !!\n ", si + 1);
	  return EGADS_MALLOC;
      }
      id6[0]    = -1;
      id6[1]    = -1;
      for (validSplit = q = 0; q < qm->star[si]->nQ; ++q) {
	  if (qm->star[si]->quads[q] == - 1) continue;
	  id0     = 2 * q + 1;
	  poly[1] = qm->star[si]->verts[id0];
	  val [1] = qm->valence[poly[1] -1][1];
	  for (i = 0; i <= qm->valence[poly[0] - 1][2] - 4; ++i) {
	      dist    = 4 + 2*i;
	      poly[2] = qm->star[si]->verts[qm->star[si]->idxV[id0 + dist]];
	      if (poly[2] < 0) continue; // ghost vertex
	      val[2]  = qm->valence[poly[2] -1][1];
	      if (val[1] == 3 && val[2] == 3 &&
		  qm->vType[poly[0] - 1] == -1) {
		  if ((qm->star[si]->nQ == 6 && dist == 6) || qm->star[si]->nQ != 6) {
		      validSplit = 1;
		      break;
		  }
	      }
	      else if (val[1] <= 3 && val[2] <= 4 && id6[0] == -1) id6[0] = id0;
	      else if (val[2] <= 3 && val[1] <= 4 && id6[0] == -1)
		id6[0] = qm->star[si]->idxV[id0 + dist];
	      else if (val[1] <= 4 && val[2] <= 4 && id6[0] == -1) id6[1] = id0;
	  }
	  if (validSplit == 1) break;
      }
      if (validSplit == 0 && val[0] >= 6 && (id6[0] != -1 || id6[1] != -1)) {
	  validSplit  = 1;
	  /* boundary vertex: Special split since regular = 3 ( not val 4)  */
	  if (qm->vType[poly[0] - 1] == 3) {
	      for (q = 0; q < qm->star[si]->nQ; q++)
		if (qm->star[si]->quads[q] == -1) break;
	      id6[0] = - 1; id6[1] = -1;
	      for (i = 0; i < 2; i++) {
		  id0 = qm->star[si]->idxQ[q + i];
		  if (i == 1) dist = (qm->star[si]->nV - 1) - 4;
		  else        dist = 4 + 2 * i;
		  links[2 * i    ] = 2 * id0 + 1;
		  links[2 * i + 1] = qm->star[si]->idxV[2 * id0 + 1 + dist];
		  if (  qm->valence[qm->star[si]->verts[links[2*i     ]] - 1][1] == 4) {
		      if (qm->valence[qm->star[si]->verts[links[2*i + 1 ]] - 1][1] == 3 &&
			  id6[0] == -1) id6[0] = i;
		      else if (qm->valence[qm->star[si]->verts[links[2*i + 1]] -1][1] == 4  &&
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
	  stat = EG_splittingOperation(qm, qm->star[si]->verts[0], qm->star[si]->verts[id0],
				       qm->star[si]->verts[qm->star[si]->idxV[id0 + dist]],
				       &(*activity));
	  return stat;
      }
  }
  return EGADS_SUCCESS;
}


static int EG_doubleSplit(meshMap *qm, quadGroup qg, int forcing, int *activity)
{
  int i, j, stat, piv[2] = {1, 5};

  *activity = 0;
  if (qg.vals[1] != 3) {
      piv[0] = 5; piv[1] = 1;
  }
  if (qg.vals[piv[0]] != 3) return EGADS_SUCCESS;
  if ((qm->vType[qg.verts[0]      - 1] >= 0 &&
      qm->vType[qg.verts[piv[1]] - 1] >= 0) ||
      qm->vType[qg.verts[0]      - 1] == 4) return EGADS_SUCCESS;
  if ((forcing == 0 && (qg.vals[0] < 5 || qg.vals[1] * qg.vals[5]  != 15)) ||
      (forcing == 1 &&  qg.vals[0] * qg.vals[piv[1]] <= 16)) return EGADS_SUCCESS;
#ifdef DEBUG
  printf(" DOUBLE SPLIT FUNCTION \n ");
#endif
  stat         = EG_splittingOperation(qm, qg.verts[0], qg.verts[piv[0]],
				       qg.verts[piv[1]], &(*activity));
  if (stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
      printf("In EG_doubleSplit: force 1st split through %d - %d -->%d activity %d !!\n ",
	     qg.verts[0], qg.verts[piv[0]], stat, *activity);
#endif
      return stat;
  }
  for (j    = 0; j < 2; j++)
    if (EG_quadVertIdx(qm, qg.q[j], qg.verts[piv[1]]) >= 0) break;
  stat      = EG_split(qm, qg.q[j], &i);
  if (stat == EGADS_SUCCESS && i == 0) {
      j     = qm->valence[qg.verts[0] - 1][0];
      stat  = EG_collapse(qm, j, &i, 1);
      if (stat == EGADS_SUCCESS && i == 1) *activity = 0;
      return stat;
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
	printf("In EG_swapDoubleSplit: adjacent to pair %d -->!!\n ", stat);
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
#ifdef DEBUG
      printQuadGroup(qm, qg);
#endif
      return stat;
  }
  if (EG_validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
  for (i0 = 0; i0 < 6; i0++) if (qg.verts[i0] == v30[0]) break;
  v0            = qg.verts[0];
  stat          = EG_swappingOperation(qm, qg, i0, &(*activity));
  if (stat     != EGADS_SUCCESS || *activity == 0) {
      if (stat != EGADS_SUCCESS)
	printf(" EG_swapDoubleSplit error at swap: %d activity %d!!\n ",
	       stat, *activity);
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
#ifdef DEBUG
      printQuadGroup(qm, qg);
#endif
      return stat;
  }
  stat   = EG_doubleSplit(qm, sqg, 0, &i);
  if (i == 0 && stat == EGADS_SUCCESS) {
      stat      = EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
      if (stat != EGADS_SUCCESS) {
	  printf(" In EG_swapDoubleSplit: recover swap create group %d\n ", stat);
#ifdef DEBUG
	  printQuadGroup(qm, qg);
#endif
	  return EGADS_SUCCESS;
      }
      for (i = 0; i < 6; i++) if (qg.verts[i] == v0) break;
      stat = EG_swappingOperation(qm, qg, i, &q);
      if (stat == EGADS_SUCCESS && q == 1) *activity = 0;
  }
  return stat;
}


static int EG_swapSplit(meshMap *qm, quadGroup qg, int forcing, int *activity)
{
  int   stat, i, j, i3 = -1, i5 = -1, v3opp = -1, q5, vL5, vL5adj, swap = 0;
  int   adj[2], v0, si;

  *activity = 0;
  if  (qg.vals[0] * qg.vals[3] != 20 ||
      EG_validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
  if      (qg.vals[1] == 3) i3 = 1;
  else if (qg.vals[5] == 3) i3 = 5;
  for (i = 1; i < 6; i++)
    if (qg.vals[i] == 5) i5 = i;
  if (i3 != -1 && (i5 == -1 || (i5 != -1 && i5 == (i3 + 3)%6))) {
      v3opp   = (i3 +3)%6;
      i5      = v3opp;
      forcing = 0;
  }
  else if (i3 == -1 && forcing == 1 && (i5 == -1 || i5 == 2 || i5 == 4)) {
      if ( i5 == -1) v3opp = -1;
      else           v3opp = i5;
  }
  else return EGADS_SUCCESS;

  if (v3opp == -1) {
      for (i = 0; i < 2; i++) {
	  j  = 2 + 2 * i;
	  if (i == 0) vL5 = 1;
	  else        vL5 = 5;
	  if (qm->valence[qg.verts[j] -1][1]  == 3) continue;
	  stat = EG_adjQtoPair(qm, qg.q[i], qg.verts[j],
			       qg.verts[vL5], adj);
	  if (stat != EGADS_SUCCESS || adj[1] == -1) continue;
	  q5 = EG_quadVertIdx(qm, adj[1], qg.verts[j]);
	  vL5adj = qm->qIdx[4 * (adj[1] - 1) + (q5 + 1)%4];
	  if (vL5adj == qg.verts[vL5])
	    vL5adj = qm->qIdx[4 * (adj[1] - 1) + (q5 + 3)%4];
	  if (qm->valence[vL5adj - 1][1] == 3) {
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
      vL5adj      = qm->qIdx[4 * (adj[1] - 1) + (i + 1)%4];
      if (vL5adj == qg.verts[vL5])
	vL5adj  = qm->qIdx[4 * (adj[1] - 1) + (i + 3)%4];
      if (i3 != -1 && (qg.vals[v3opp] == 5 || qm->valence[vL5adj - 1][1] == 3))
	swap = i3;
      else if (forcing == 1 && (qg.vals[v3opp] == 5 || qm->valence[vL5adj - 1][1] == 3))
	swap = v3opp;
  }
  if (swap%3 == 0 || (qm->vType[qg.verts[swap] - 1] > 0 &&
      qm->vType[qg.verts[(swap + 3)%6] - 1] > 0)) return EGADS_SUCCESS;
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

  si = qg.verts[i5] - 1;
  if ( qm->star[si] == NULL ) {
      printf("In swapSplit star at %d is NULL !!\n", qg.verts[i5]);
      return EGADS_MALLOC;
  }
  for (i = 0; i < qm->star[si]->nQ; i++)
    if (qm->star[si]->verts[2 * i + 1] == qg.verts[3]) break;
  if (qm->star[si]->verts[2 * i + 1] != qg.verts[3]) {
      return EGADS_INDEXERR;
  }
  adj[0] = qm->star[si]->verts[qm->star[si]->idxV[2 * i + 1 + 6]];
  if (qm->star[si]->nQ == 5) {
      adj[1] = qm->star[si]->verts[qm->star[si]->idxV[2 * i + 1 + 4]];
      if (qm->valence[adj[1] - 1][1] < qm->valence[adj[0] - 1][1])
	adj[0] = adj[1];
  }
#ifdef DEBUG
  printf("In swapSplit call splitting op using %d %d %d !!\n",
	 qg.verts[i5], qg.verts[3], adj[0]);
#endif

  stat         = EG_splittingOperation(qm, qg.verts[i5], qg.verts[3], adj[0], &i);
  if (stat    != EGADS_SUCCESS || i == 0) {
#ifdef DEBUG
      printf("In swapSplit splittingOperation -->%d !!\n", stat);
#endif
      stat      = EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
      for (swap = 0; swap < 6; swap++) if (qg.verts[swap] == v0) break;
      stat      = EG_swappingOperation(qm, qg,swap, &i);
      if (stat == EGADS_SUCCESS && i == 1) *activity = 0;
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
  if (qm->qIdx[4 * (qID - 1)] == -2) {
#ifdef DEBUG
      printf(" EG_cleanQuad %d stat %d\n", qID, stat);
#endif
      return EGADS_SUCCESS;
  }
  if (EG_nValenceCount(qm, qID, 4) == 4) return EGADS_SUCCESS;
  if (transfer == 0) {
      for (i = 0; i < 3; i++) {
	  stat       = EG_basicOperation(qm, qID, opBasic[i], &act);
	  *activity += act;
	  if (stat != EGADS_SUCCESS) {
	      printf(" EG_cleanQuad Basic operation %d !!\n ", stat);
	      return stat;
	  }
	  if (act > 0 || qID < 0 || qm->qIdx[4 * (qID - 1)] == -2)
	    return EGADS_SUCCESS;
      }
  }
  if (useAdj == 0)  return EGADS_SUCCESS;
  if (qID < 0 || qm->qIdx[4 * (qID - 1)] == -2) {
#ifdef DEBUG
      printf(" EG_cleanQuad %d is empty !!\n", qID);
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
	  qadj      = qm->qAdj[4 * (qID - 1) + q];
	  if (qadj == -1) continue;
	  stat       = EG_cleanQuad(qm, qadj, 0, 0, 0, &act);
	  *activity += act;
	  if (stat != EGADS_SUCCESS) {
	      printf(" EG_cleanQuad restore mesh after adjacent stat -->%d !!\n",
		     stat);
	      return stat;
	  }
	  if (act > 0 || qID < 0 || qm->qIdx[4 * (qID - 1)] == -2)
	    return EGADS_SUCCESS;
      }
  }
  for (i = 0; i < 7; i++) {
      if (transfer == 1 && opComp[i] == DOUBLESWAP) continue;
      for (q  = 0; q < 4; q++) {
	  if (qID < 0 || qm->qIdx[4 * (qID - 1)] == -2)
	    return EGADS_SUCCESS;
	  qadj = qm->qAdj[4 * (qID - 1) + q];
	  if (qadj == -1) continue;
	  stat = EG_createQuadGroup(qm, &qg, qID, qadj);
	  if (stat != EGADS_SUCCESS) {
	      printf("Inside EG_cleanQuad: EG_createQuadGroup stat -->%d !!\n ",
		     stat);
#ifdef DEBUG
	      printQuadGroup(qm, qg);
#endif
	      return stat;
	  }
	  stat       = EG_composeOperation(qm, qg, opComp[i], forcing, &act);
	  *activity += act;
	  if (stat != EGADS_SUCCESS) {
	      printf("In EG_cleanQuads restor mesh after EG_composeOperation stat -->%d !!\n ",
		     stat);
	      return stat;
	  }
	  if (act > 0 || qID < 0 || qm->qIdx[4 * (qID - 1)] == -2)
	    return EGADS_SUCCESS;
      }
  }

  return EGADS_SUCCESS;
}


static int EG_transferValences(meshMap *qm, int *qID, int try5533,
			       int *transfering, int *activity)
{
  int       i, j, swap = 0, stat, min, qAux[2];
  quadGroup qg;

  *activity = 0;
  if (qID[0] < 0 || qm->qIdx[4 * (qID[0] - 1)] == -2) {
#ifdef DEBUG
      printf(" EG_transferValences for quad %d is out of bounds !!\n ", qID[0]);
#endif
      return EGADS_SUCCESS;
  }

  if (*transfering == 0) {
      stat      = EG_cleanQuad(qm, qID[0], 1, 0, try5533, &(*activity));
      if (stat != EGADS_SUCCESS || *activity > 0) {
	  qID[0] = -1;
	  if (stat != EGADS_SUCCESS)
	    printf("EG_transferValences: EG_cleanQuad %d -->%d!!\n",
		   qID[0], stat);
	  return stat;
      }
  }
  if (qID[0] == -1 || qm->qIdx[4 * (qID[0] - 1)] == -2 ||
      EG_nValenceCount(qm, qID[0], 4) == 4) return EGADS_SUCCESS;
  for (j = 0; j < 4; j++) {
      i = qm->qAdj[4 * (qID[0] - 1) + j];
      if (i == -1 || i == qID[1]) continue;
      if (*transfering == 0) {
	  stat      = EG_cleanQuad(qm, i, 1, 0, try5533, &(*activity));
	  if (stat != EGADS_SUCCESS || *activity > 0) {
	      if (stat != EGADS_SUCCESS)
		printf(" EG_TransferValence EG_cleanQuad -->%d !!\n ", stat);
	      qID[0] = -1;
	      return stat;
	  }
      }
      if (qm->qAdj[4 * (qID[0] - 1) + j] != i) continue;
      stat = EG_createQuadGroup(qm, &qg, qID[0], i);
      if (stat != EGADS_SUCCESS) {
	  printf(" Inside EG_transferValences EG_createQuadGroup %d !!\n", stat);
#ifdef DEBUG
	  printQuadGroup(qm, qg);
#endif
	  return stat;
      }
      if (*transfering == 0 && qg.vals[0] * qg.vals[3] == 15 && qID[1] != -1) {
	  j = 0;
	  if (qg.q[j] == qID[0]) j = 1;
	  qID[0] = qg.q[j];
	  qID[1] = qg.q[(j + 1)%2];
	  for ( i = 0 ; i < 4; i++ ) {
	      j = qm->qIdx [4 * (qID[0] - 1 ) + i ] - 1;
	      if (qm->vType[qm->valence[j][3] - 1] > 0 ) {
		  stat      = EG_collapse(qm, qID[0], &(*activity), 1);
		  if (stat != EGADS_SUCCESS)
		    printf("EG_transferValences: forceCollapse gave %d !!\n", stat);
		  qID[0]    = qID[1];
		  qID[1]    = -1;
		  return stat;
	      }
	  }
	  return EG_transferValences(qm, qID, 0, &(*transfering), &(*activity));
      }
      if (EG_validSwap(qm, qg.verts[0], qg.verts[3]) == 0) continue;
      if          (qg.vals[0] * qg.vals[3] >= 20) min = 12;
      else if (try5533 == 1 && *transfering == 0) {
	  if      (qg.vals[0] * qg.vals[3] >= 25) min = 16;
	  else if (qg.vals[0] * qg.vals[3] == 16) min =  9;
      } else continue;
      if (     qg.vals[1] * qg.vals[4] <= min) swap = 1;
      else if (qg.vals[2] * qg.vals[5] <= min) swap = 2;
      else continue;
#ifdef DEBUG
      printf("Transfer valences thru %d MIN %d\n ", swap, min);
      printQuadGroup(qm, qg);
#endif
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
#ifdef DEBUG
      printf(" TF CURRENT Q %d %d  ( q0 %d q1 %d) \n",
             qg.q[0], qg.q[1], qID[0], qID[1]);
#endif
      if (qID[1] > 0 && qm->qIdx[4 * (qID[1] - 1)] != -1 &&
          qg.q[0] == qID[1]) i = 1;
      else if (EG_nValenceCount(qm,qg.q[0], 5) == 2 &&
               EG_nValenceCount(qm,qg.q[1], 5) >= 1 &&
               EG_nValenceCount(qm,qg.q[1], 3) == 1) i = 1;
      else if (EG_nValenceCount(qm,qg.q[1], 4) <
               EG_nValenceCount(qm,qg.q[0], 4)) i = 1;
      qID[0] = qg.q[ i       ];
      qID[1] = qg.q[(i + 1)%2];
#ifdef DEBUG
      printf(" ASSIGNED QUAS %d %d TRANSFER %d\n", qID[0], qID[1], *transfering);
#endif
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


void EG_destroyMeshMap(bodyQuad *bodydata)
{
  int i, j;

  if (bodydata->qm == NULL) return;
  for (i = 0; i < bodydata->nfaces; ++i) {
      if (bodydata->qm[i]) {
	  for (j = 0; j < bodydata->qm[i]->sizeV; ++j) {
	      EG_free(bodydata->qm[i]->valence[j]);
	      EG_freeStar(&bodydata->qm[i]->star[j]);
	  }
	  EG_free(bodydata->qm[i]->valence);
	  EG_free(bodydata->qm[i]->qIdx);
	  EG_free(bodydata->qm[i]->qAdj);
	  EG_free(bodydata->qm[i]->uvs);
	  EG_free(bodydata->qm[i]->xyzs);
	  EG_free(bodydata->qm[i]->remQ);
	  EG_free(bodydata->qm[i]->remV);
	  EG_free(bodydata->qm[i]->vType);
	  EG_free(bodydata->qm[i]->star);
	  EG_free(bodydata->qm[i]);
      }
  }
  EG_free(bodydata->qm);
  bodydata->qm = NULL;
}


static int resizeQm(meshMap *qm)
{
  int    stat = EGADS_SUCCESS, nV, vRem, nQ, qRem, i, j, k;
  int    *vpiv = NULL, *qpiv = NULL;
  int    *verts = NULL, *adjs = NULL, *type = NULL, bV, bQ;
  double *uvxyz = NULL;

  vpiv = (int *) EG_alloc(qm->totV * sizeof(int));
  qpiv = (int *) EG_alloc(qm->totQ * sizeof(int));
  if (vpiv == NULL || qpiv == NULL) return EGADS_MALLOC;
  for (nV = vRem = i = 0; i < qm->totV; i++) {
      if (qm->vType[i] != -2) vpiv[i] = nV++;
      else {
	  vpiv[i] = -2;
	  vRem++;
      }
  }
  for (nQ = qRem = i = 0; i < qm->totQ; i++) {
      if (qm->qIdx[4*i] != -2) qpiv[i] = nQ++;
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
  bV    = qm->totV;
  bQ    = qm->totQ;
  uvxyz = (double *) EG_alloc(5 * bV *sizeof(double));
  type  = (int    *) EG_alloc(    bV *sizeof(int));
  verts = (int    *) EG_alloc(4 * bQ *sizeof(int));
  adjs  = (int    *) EG_alloc(4 * bQ *sizeof(int));
  if (uvxyz == NULL || type == NULL || verts == NULL || adjs == NULL) {
      stat = EGADS_MALLOC;
      goto cleanup;
  }
  for (j = 0; j < bQ; ++j) {
      for (k = 0; k < 4; ++k) {
	  verts[4*j + k] = qm->qIdx[4*j + k];
	  adjs [4*j + k] = qm->qAdj[4*j + k];
      }
  }
  for (j = 0; j < bV; ++j) {
      type [  j    ] = qm->vType[  j    ];
      uvxyz[5*j    ] = qm->uvs  [2*j    ];
      uvxyz[5*j + 1] = qm->uvs  [2*j + 1];
      uvxyz[5*j + 2] = qm->xyzs [3*j    ];
      uvxyz[5*j + 3] = qm->xyzs [3*j + 1];
      uvxyz[5*j + 4] = qm->xyzs [3*j + 2];
  }
  qm->totV    = nV;
  qm->totQ    = nQ;
  qm->remQ[0] = 0;
  qm->remV[0] = 0;
  for (i = 0; i < bQ; i++) {
      if (qpiv[i] == -2) continue;
      for (k = 0; k < 4; ++k) {
	  qm->qIdx[4*qpiv[i] + k ] = vpiv[verts[4*i + k] - 1] + 1;
	  if (adjs[4*i + k] == -1) qm->qAdj[4*qpiv[i] + k] = - 1;
	  else                     qm->qAdj[4*qpiv[i] + k] = qpiv[adjs[4*i + k] - 1] + 1;
	  qm->valence[qm->qIdx[4*qpiv[i] + k] - 1][0] = qpiv[i] + 1;
      }
  }
  for (i = 0; i < bV; i++) {
      if (vpiv[i] == -2) continue;
      j = vpiv[i];
      qm->vType  [j      ] = type [  i    ];
      qm->uvs    [2*j    ] = uvxyz[5*i    ];
      qm->uvs    [2*j + 1] = uvxyz[5*i + 1];
      qm->xyzs   [3*j    ] = uvxyz[5*i + 2];
      qm->xyzs   [3*j + 1] = uvxyz[5*i + 3];
      qm->xyzs   [3*j + 2] = uvxyz[5*i + 4];
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


static int EG_cleanNeighborhood(meshMap *qm, int qID,  int transfer,
				int *activity)
{
  int si, i, act, stat = 0, count = 0, j, v[4];

  *activity = 0;
  if (qm->qIdx[4 * (qID - 1)] == -2) return EGADS_SUCCESS;
  for (i = 0; i < 4; i++) v[i]  = qm->qIdx[4 * (qID - 1) + i];
  stat      = EG_cleanQuad(qm, qID, 1, transfer, 0, &act);
  if (stat != EGADS_SUCCESS) return stat;
  *activity += act;
  for (i = 0; i < 4; i++) {
      if (qm->vType[v[i] - 1] == -2) continue;
      si = v[i] - 1;
      if ( qm->star[si] == NULL ) {
	  printf(" EG_cleanNeighborhood star at %d is NULL !!\n", v[i]);
	  return EGADS_MALLOC;
      }
      for (j = 0; j < qm->star[si]->nQ; j++) {
	  if (qm->star[si]->quads[j] == -1) continue;
	  if (EG_cleanQuad(qm, qm->star[si]->quads[j], 1, transfer, 0, &act) != EGADS_SUCCESS) {
	      return stat;
	  }
	  count += act;
      }
  }
  *activity += count;
#ifdef DEBUG
  return checkMesh(qm);
#else
  return EGADS_SUCCESS;
#endif
}


static void EG_distToSurfCentre (meshMap *qm, double *l2, double *linf) {
  int i, j, k;
  double quv[2], pos[3], pInv[3], aux;

  *l2 = 0.0; *linf = 0.0;
  for (i = 0; i < qm->totQ; i++) {
      if ( qm->qIdx [4 * i] == -2 ) continue;
      pos[0] = pos[1] = pos[2] = 0.0;
      quv[0] = quv[1] = 0.0;
      for ( k = 0; k < 4; k++ ) {
	  j       = qm->qIdx [ 4 * i + k] - 1;
	  pos[0] += 0.25 * qm->xyzs[ 3 * j    ];
	  pos[1] += 0.25 * qm->xyzs[ 3 * j + 1];
	  pos[2] += 0.25 * qm->xyzs[ 3 * j + 2];
	  quv[0] += 0.25 * qm->uvs [ 2 * j    ];
	  quv[1] += 0.25 * qm->uvs [ 2 * j + 1];
      }
      if (quv[0] < qm->range[0] || quv[0] > qm->range[1] ||
	  quv[1] < qm->range[2] || quv[1] > qm->range[3] ) {
	  quv[0] = qm->uvs [2 * j    ];
	  quv[1] = qm->uvs [2 * j + 1];
      }
      j = EG_invEvaluateGuess(qm->face, pos, quv, pInv);
      if (j != EGADS_SUCCESS ) {
	  j  = EG_invEvaluate(qm->face, pos, quv, pInv);
	  if (j != EGADS_SUCCESS ) continue;
      }
      aux = (pos[0] - pInv[0]) * ( pos[0] - pInv[0] ) +
	  (pos[1] - pInv[1]) * ( pos[1] - pInv[1] ) +
	  (pos[2] - pInv[2]) * ( pos[2] - pInv[2] ) ;
      *linf = MAX (*linf, aux );
      *l2  += aux;
  }
  return;
}



int EG_meshRegularization(meshMap *qm)
{
  int    i, j, k, q, stat = EGADS_SUCCESS, round, s1, s2, vq[4], si, ni = 0, n0 = 0;
  int    ITMAX, it = 0, activity = 0, totActivity = 0, loopact;
  int    iV0, iV, qPair[4], prevPair[2], totV, totV0, vQ0, vQ, transfer = 0, *skipQuad = NULL, sq;
  double l2Init, linfInit, l2Actual, linfActual, stol, norm1, norm2, err;
  double pos[18*2], quv[2], cross[3], v1[3], v2[3], qN[3], *qArea = NULL;
  quadGroup qg;
#ifdef REPORT
  char   buffer[100];
  clock_t time0, time1, ttime;
  int min;
  float secs;
#endif

  // GET RANGE FOR EACH POINT
  meshCount(qm, &iV0, &totV0, &vQ0);
  // PREPROCESS: Reduce the number of quads by forcing collapses everywhere
  EG_distToSurfCentre(qm, &l2Init, &linfInit);
#ifdef REPORT
  time0 = clock();
#endif
  if ( l2Init > 1.e-14 && linfInit > 1.e-14 ) {
      n0 = ni = it = sq = 0;
      skipQuad   = EG_alloc (qm->totQ * sizeof(int    ));
      qArea      = EG_alloc (qm->totQ * sizeof(double ));
      if (skipQuad == NULL || qArea == NULL) return EGADS_MALLOC;
      stol  = 1.1;
      ITMAX = qm->totQ;
      round = 0;
      while ( it < ITMAX ) {
	  EG_distToSurfCentre(qm, &l2Actual, &linfActual);
	  k = -1;
	  for (si = i = 0; i < qm->totQ; i++) {
	      if (qm->qIdx[ 4 * i ]           == -2 ||
		  inList (sq, skipQuad, i + 1 ) != -1 ) continue;
	      pos[0] = pos[1] = pos[2] = 0.0;
	      quv[0] = quv[1] = 0.0;
	      for ( q = 0; q < 4; q++ ) {
		  vq[q]   = qm->qIdx [4 * i + q] - 1;
		  quv[0] += 0.25 * qm->uvs [ 2 * vq[q]    ];
		  quv[1] += 0.25 * qm->uvs [ 2 * vq[q] + 1];
		  pos[0] += 0.25 * qm->xyzs[ 3 * vq[q]    ];
		  pos[1] += 0.25 * qm->xyzs[ 3 * vq[q] + 1];
		  pos[2] += 0.25 * qm->xyzs[ 3 * vq[q] + 2];
	      }
	      if (quv[0] < qm->range[0] || quv[0] > qm->range[1] ||
		  quv[1] < qm->range[2] || quv[1] > qm->range[3] ) {
		  quv[0] = qm->uvs [2 * vq[0]    ];
		  quv[1] = qm->uvs [2 * vq[0] + 1];
	      }
	      stat  = EG_invEvaluateGuess(qm->face, pos, quv, &pos[18] );
	      err   = (pos[0] - pos[18]) * (pos[0] - pos[18]) +
		  (pos[1] - pos[19]) * (pos[1] - pos[19]) +
		  (pos[2] - pos[20]) * (pos[2] - pos[20]);
	      stat  = EG_evaluate(qm->face, quv, pos );
	      norm1 = pos[3] * pos[3] + pos[4] * pos[4] + pos[5] * pos[5];
	      norm2 = pos[6] * pos[6] + pos[7] * pos[7] + pos[8] * pos[8];
	      if (norm1 < EPS08 || norm2 < EPS08 ) continue;
	      v1[0] = pos[3]; v1[1] = pos[4]; v1[2] = pos[5];
	      v2[0] = pos[6]; v2[1] = pos[7]; v2[2] = pos[8];
	      if (qm->face->mtype == SREVERSE) {
		  CROSS(v2, v1, qN);
	      } else {
		  CROSS(v1, v2, qN);
	      }
	      v1[0] = qm->xyzs[3 * vq[1]    ] - qm->xyzs[3 * vq[0]    ];
	      v1[1] = qm->xyzs[3 * vq[1] + 1] - qm->xyzs[3 * vq[0] + 1];
	      v1[2] = qm->xyzs[3 * vq[1] + 2] - qm->xyzs[3 * vq[0] + 2];
	      v2[0] = qm->xyzs[3 * vq[2]    ] - qm->xyzs[3 * vq[0]    ];
	      v2[1] = qm->xyzs[3 * vq[2] + 1] - qm->xyzs[3 * vq[0] + 1];
	      v2[2] = qm->xyzs[3 * vq[2] + 2] - qm->xyzs[3 * vq[0] + 2];
	      s1    = s2  = 1;
	      CROSS(v1, v2, cross);
	      if (DOT(qN, cross) < 0) s1 = -1;
	      qArea[i] = 0.5 * sqrt ( DOT (cross, cross) ) ;
	      v1[0]    = qm->xyzs[3 * vq[3]    ] - qm->xyzs[3 * vq[0]    ];
	      v1[1]    = qm->xyzs[3 * vq[3] + 1] - qm->xyzs[3 * vq[0] + 1];
	      v1[2]    = qm->xyzs[3 * vq[3] + 2] - qm->xyzs[3 * vq[0] + 2];
	      CROSS(v2, v1, cross);
	      if (DOT(qN, cross) < 0) s2 = -1;
	      qArea[i] += 0.5 * sqrt ( DOT (cross, cross) ) ;
	      if ( s1 * s2 == -1 ) {
#ifdef DEBUG
		  printf(" Quad %d self-intersects: try destroy\n ", i + 1 );
		  gnuData(qm, NULL);
		  printQuad(qm, i + 1 ) ;
#endif
		  qArea[i] = -1.0; //invalid quad: self-intersects.
		  k   = i;
		  si  = 1;
		  if ( n0 == 0 ) ni++;
	      }
	      if (err >= l2Actual) continue;
	      if (si == 0 && ((k == -1 && EG_nValenceCount(qm, i + 1, 3) > 0) ||
	          (k != -1 && fabs(qArea[i]) < fabs(qArea[k]) ))) k = i;
	  }
	  if ( n0 == 0 && ni > 0 ) {
	      qm->qInv = EG_alloc ((ni + 1) * sizeof (int));
	      if (qm->qInv == NULL) {
	          EG_free (skipQuad);
	          EG_free (qArea   );
	          return EGADS_MALLOC;
	      }
	      qm->qInv[0] = ni;
	      for (j = i = 0 ; i < qm->totQ;i++ )
	        if (qArea[i] < 0.0) qm->qInv[++j] = i + 1;
	  }
	  n0 = 1; // only first round is allowed to have invalid quads (self int)
	  if (qm->qInv && qm->qInv[0] > 0) {
	      for (j = 0 ; j < qm->qInv[0]; j++) {
		  if (qm->qInv[1 + j] == -2 ||
		      qm->qIdx[4 * (qm->qInv[ 1 + j ] - 1)] == -2 )  continue;
		  stat  = EG_collapse(qm, qm->qInv[1 + j], &activity, 3);
		  if (stat != EGADS_SUCCESS) {
		      printf(" EG_meshRegularization :: preprocess EG_collapse for quad %d -->%d!!\n ",
			     qm->qInv[1 + j], stat);
		      EG_free (skipQuad);
		      EG_free (qArea   );
		      EG_free(qm->qInv);
		      return stat;
		  }
		  if (activity > 0) {
		      it++;
		      qm->qInv[1 + j] = -2;
		  }
	      }
	  }
#ifdef DEBUG
	  printf("IT %d (ROUND %d) L2 err %1.2e  vs %1.2e Linf %1.2e VS %1.2e -->%.2f %%\n",
		 it, round, l2Init, l2Actual, linfInit, linfActual, ( l2Actual - l2Init ) / l2Init * 100.0 );
#endif
	  if (l2Actual > l2Init * stol ) break;
	  if (k == -1 || qm->qIdx[ 4 * k ] == -2 ) {
	      for (k = j = i = 0; i < qm->totQ; i++) {
		  if (qm->qIdx [ 4 * i ]          == -2 ||
		      inList (sq, skipQuad, i + 1 ) != -1 ) continue;
		  if (fabs(qArea[i]) < fabs(qArea[k]) || j == 0  ) {
		      k = i;
		      j = 1;
		  }
	      }
	      if ( j == 0 ) {
		  if ( round < 20 ) sq = 0 ;
		  else break;
		  round++;
		  continue;
	      }
	  }
	  stat  = EG_collapse(qm, k + 1, &activity, 2);
	  if (stat != EGADS_SUCCESS) {
	      printf(" EG_meshRegularization :: preprocess EG_collapse for quad %d -->%d!!\n ",
		     k + 1, stat);
	      EG_free(skipQuad);
	      EG_free(qArea   );
	      EG_free(qm->qInv);
	      return stat;
	  }
	  if (activity > 0) it++;
	  else {
	      if ( sq > qm->totQ -1 ) {
		  printf(" beyond bounds :: size sq %d and current %d\n ", qm->totQ, sq );
		  break;
	      } skipQuad[sq++] = k + 1;
	  }
      }
      EG_free (skipQuad);
      EG_free (qArea   );
      if (qm->qInv && qm->qInv[0] > 0 ) {
	  for (j = 0 ; j < qm->qInv[0]; j++ ) {
	      if ( qm->qInv[ 1 + j ] == -2 ) continue;
	      fprintf(stderr,"After Preprocessing there are still invalid quads from the original tessellation "
		      "that couldn't be eliminated. Mesh is in valid. Change parameters! \n ");
	      EG_free(qm->qInv);
	      qm->qInv = NULL;
	      return EGADS_GEOMERR;
	  }
      }
      EG_free(qm->qInv);
      qm->qInv = NULL;
      stat       = resizeQm(qm );
      if (stat != EGADS_SUCCESS) {
	  printf(" EG_meshRegularization resizing after preprocess %d !!\n", stat);
	  return stat;
      }

  }
  EG_makeValidMesh(qm, 0, NULL, 1);
#ifdef REPORT
  time1 = clock();
  ttime = (double)(time1 - time0) / CLOCKS_PER_SEC;
  min   = floor(ttime) / 60;
  secs  = ttime -  min * 60;
  printf("\n\n");
  printf("****************** PREPROCESSING: %d collapses ******************\n", it);
  EG_distToSurfCentre(qm, &l2Actual, &linfActual);
  printf(" Centroid error L2 %1.2e vs %1.2e %1.2f %%\n", l2Actual, l2Init, (l2Actual - l2Init) / l2Init * 100.0);
  meshCount(qm, &iV, &totV, &vQ);
  printf(" Original mesh: %d quads and %d / %d irregular vertices (%.2f %%)\n",
	 vQ0, iV0, totV0, (double) iV0 * 100.0 / (double)totV0);
  printf(" Current  mesh: %d quads and %d / %d irregular vertices (%.2f %%)\n",
	 vQ, iV, totV, (double) iV * 100.0 / (double) totV);
  printf(" Invalid steps  %d\t CPU time: %d mins and %f sec\n",
	 qm->invsteps, min, secs);
  printf("__________________________IO FILES ______________________________\n");
  snprintf(buffer,100, "gnuPreprocess_%d.txt", qm->fID);
  gnuData(qm, buffer);
  snprintf(buffer,100, "wvsPreprocess_%d.txt", qm->fID);
  wvsData(qm, buffer);
  snprintf(buffer,100, "preprocessMeshStats_%d.txt", qm->fID);
  meshStats(qm, 10);
  qm->invsteps = 0;
  printf("*****************************************************************\n");
  time0 = clock();
#endif
  ITMAX = 100;
  for (it = 0; it < ITMAX; it++) {
      totActivity = 0;
      for (i = 0; i < qm->totQ; i++) {
	  if (qm->qIdx[4 * i] == -2) continue; //can be a deleted quad
	  stat         = EG_cleanNeighborhood(qm, i + 1, 0, &activity);
	  totActivity += activity;
	  if (stat    != EGADS_SUCCESS) {
	      printf(" In EG_cleanMesh: EG_CleanNeighborhood for quad %d -->%d!!\n ",
		     i + 1, stat);
	      return stat;
	  }
      }
      meshCount(qm, &iV, &totV, &vQ);
      if (totActivity == 0 || iV <= 2) break;
  }
#ifdef DEBUG
  printf(" TRANSFER VALENCES \n ");
  snprintf(buffer,100, "notransfer_%d.txt", qm->fID);
  gnuData(qm, buffer);
#endif
  if (iV > 2) {
      ITMAX = 100;
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
	      for (q       = 0; q < qm->totQ; q++) {
		  transfer = 0;
		  if (qm->qIdx[4 * q] == -2 || q + 1 == prevPair[0] ||
		      q + 1 == prevPair[1]) continue;
		  qPair[0]   = q + 1;
		  qPair[1]   = -1;
		  stat          = EG_transferValences(qm, qPair, k,
						      &transfer, &activity);
		  if (stat     != EGADS_SUCCESS) return stat;
		  if (activity == 0 || qPair[0] < 0 ||
		      qm->qIdx[4*(qPair[0]-1)] == -2) continue;
		  totActivity += activity;
		  for (j = 0; j < 20; j++) {
		      activity = 0;
		      for (i  = 0; i < 2; i++) {
			  if (qPair[i] < 0 ||
			      qm->qIdx[4*(qPair[i]-1)] == -2) continue;
			  stat      = EG_cleanNeighborhood(qm, qPair[i],
							   transfer, &loopact);
			  if (stat != EGADS_SUCCESS) return stat;
			  activity +=loopact;
		      }
		      if (activity > 0) break;
		      stat      = EG_transferValences(qm, qPair, k,
						      &transfer, &activity);
		      if (stat != EGADS_SUCCESS) return stat;
		      if (activity == 0 || qPair[0] < 0 ||
			  qm->qIdx[4*(qPair[0]-1)] == -2) break;
		  }
		  if (iV <= 2) break;
		  prevPair[0] = qPair[0];
		  prevPair[1] = qPair[1];
	      }
	      if (iV <= 2) break;
	      for (q = 0; q < qm->totQ; q++) {
		  if (qm->qIdx[4 * q] == -2) continue;
		  stat = EG_cleanQuad(qm, q + 1, 1, 0, 0, &i);
		  if (stat != EGADS_SUCCESS) {
		      printf(" In EG_meshRegularization clean quad %d !!\n ",
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
  stat      = resizeQm(qm);
  if (stat != EGADS_SUCCESS) {
      printf("EG_meshRegularization final resize %d !!\n", stat);
      return stat;
  }
#ifdef REPORT
  time1 = clock();
  ttime = (double)(time1 - time0) / CLOCKS_PER_SEC;
  min   = floor(ttime) / 60;
  secs  = ttime -  min * 60;
  meshCount(qm, &iV, &totV, &vQ);

  printf("\n\n");
  printf("********************* REGULARIZATION RESULTS *********************\n");
  EG_distToSurfCentre(qm, &l2Actual, &linfActual);
  printf(" Centroid errors L2 %1.2e vs %1.2e Linf %1.2e vs %1.2e\n",
	 l2Actual,l2Init,  linfActual, linfInit );
  meshCount(qm, &iV, &totV, &vQ);
  printf(" Original mesh: %d quads and %d / %d irregular vertices (%.2f %%) \n",
	 vQ0, iV0, totV0, (double) iV0 * 100.0 / (double)totV0);
  printf(" Final    mesh: %d quads and %d / %d irregular vertices (%.2f %%)\n",
	 vQ, iV, totV, (double) iV * 100.0 / (double) totV);
  printf(" Invalid steps  %d\t CPU time: %d mins and %f sec\n",
	 qm->invsteps, min, secs);
  printf("__________________________ IO FILES _______________________________\n");
  snprintf(buffer,100, "regularizationMeshStats_%d.txt", qm->fID);
  meshStats(qm, 10);
  snprintf(buffer,100, "wvsRegular_%d.txt", qm->fID);
  wvsData(qm, buffer);
  snprintf(buffer,100,"gnuRegular_%d.txt", qm->fID);
  gnuData(qm, buffer);
  printf("*******************************************************************\n");
#endif
  stat = EG_makeValidMesh(qm, 0, NULL, 2);
  if (qm->qInv == NULL) return stat;
  skipQuad = EG_alloc(qm ->qInv[0] * sizeof(int));
  if (skipQuad == NULL) return EGADS_MALLOC;
  for (i = q = 0; q < qm->qInv[0]; q++) {
      if (qm->vType[qm->qIdx[4 * (qm->qInv[1+ q] - 1)    ] -1] >= 4 ||
          qm->vType[qm->qIdx[4 * (qm->qInv[1+ q] - 1) + 1] -1] >= 4 ||
          qm->vType[qm->qIdx[4 * (qm->qInv[1+ q] - 1) + 2] -1] >= 4 ||
          qm->vType[qm->qIdx[4 * (qm->qInv[1+ q] - 1) + 3] -1] >= 4 )
        skipQuad[i++] = qm->qInv[1+ q];
  }
  for (q = 0; q < qm->qInv[0]; q++) {
      if (qm->vType[qm->qIdx[4 * (qm->qInv[1+ q] - 1)    ] -1] >= 4 ||
          qm->vType[qm->qIdx[4 * (qm->qInv[1+ q] - 1) + 1] -1] >= 4 ||
          qm->vType[qm->qIdx[4 * (qm->qInv[1+ q] - 1) + 2] -1] >= 4 ||
          qm->vType[qm->qIdx[4 * (qm->qInv[1+ q] - 1) + 3] -1] >= 4 ) continue;
      skipQuad[i++] = qm->qInv[1+ q];
  }
  ni = qm->qInv[0];
  EG_free(qm->qInv);
  qm->qInv = NULL;
#ifdef DEBUG
  printf(" check invalid quad list: Try to instert / collapse quads to make mesh valid\n");
#endif
  for ( q = 0; q < ni; q++) {
#ifdef DEBUG
      printf(" OPERATE AROUND QUAD %d \n ",skipQuad[q] );
#endif
      if (skipQuad[q] <= 0 ||
          skipQuad[q] > qm->totQ ||
          qm->qIdx[4 * (skipQuad[q] - 1)] == -2 ) {
          fprintf(stderr," WHAT ARE YOU TRYING TO DO ??? Q %d \n", skipQuad[q]);
          printf(" WHAT ARE YOU TRYING TO DO ??? Q %d \n", skipQuad[q]);
          continue;
      }
      for (j = i = 0 ; i < 4; i++) {
          k = qm->qIdx[4 * (skipQuad[q] - 1) + i] - 1;
          if (qm->vType[k] == -1) EG_placeVertex(qm, k + 1, 0.5);
          if (qm->vType[k] == -1 && qm->valence[k][2] > 3 &&
              qm->vType[qm->valence[k][3] -1] >= 4) {
#ifdef DEBUG
              printf(" IT IS A CORNER CASE Q %d VERT %d!!!\n",
                     qm->valence[k][3], k + 1);
#endif
              j = k + 1;
          }
      }
#ifdef DEBUG
      printf(" FINISHED DOING STUFF AROUND %d\n", skipQuad[q]);
#endif
      i = EG_quadArea(qm, NULL, NULL, skipQuad[q], -1, pos);
      if (i == QA0 && pos[0] >= 0.1) {
#ifdef DEBUG
          printf(" QUAD %d IS FIXED !!!\n", skipQuad[q]);
          gnuData(qm, NULL);
#endif
          continue;
      }
#ifdef DEBUG
      printf("LOOK FOR SWAPS FOR QUAD %d j %d\n ", skipQuad[q], j);
#endif
      if (j == 0) {
          s1     = skipQuad[q] - 1;
          vq[0]  = qm->qIdx[4*(skipQuad[q] -1)    ];
          vq[1]  = qm->qIdx[4*(skipQuad[q] -1) + 1];
          vq[2]  = qm->qIdx[4*(skipQuad[q] -1) + 2];
          vq[3]  = qm->qIdx[4*(skipQuad[q] -1) + 3];
          stat   = EG_angleAtVnormalPlane(qm, vq[0], vq[1], vq[3], &pos[0]);
          stat  += EG_angleAtVnormalPlane(qm, vq[1], vq[2], vq[0], &pos[1]);
          stat  += EG_angleAtVnormalPlane(qm, vq[2], vq[3], vq[1], &pos[2]);
          stat  += EG_angleAtVnormalPlane(qm, vq[3], vq[0], vq[2], &pos[3]);
          printf(" ANGLE %d = %lf %d = %lf %d = %lf %d = %lf \n",
                 vq[0], pos[0], vq[1], pos[1], vq[2], pos[2], vq[3], pos[3]);
          if (stat != EGADS_SUCCESS ) continue;
          i = 0;
          if (pos[1] > pos[i] ) i = 1;
          if (pos[2] > pos[i] ) i = 2;
          if (pos[3] > pos[i] ) i = 3;
          j  = vq[(i + 2)%4];
          s1 = vq[(i + 3)%4];
          s2 = vq[(i + 3)%4];
          stat = EG_adjQtoPair(qm, skipQuad[q], j, s1,  qPair);
          stat = EG_adjQtoPair(qm, skipQuad[q], j, s2, &qPair[2]);
          k    = 1;
          if (qPair[3] != -1 ||
              qm -> qAdj[4 *(qPair[3] - 1)    ] == -1 ||
              qm -> qAdj[4 *(qPair[3] - 1) + 1] == -1 ||
              qm -> qAdj[4 *(qPair[3] - 1) + 2] == -1 ||
              qm -> qAdj[4 *(qPair[3] - 1) + 3] == -1) k = 3;
          stat     = EG_createQuadGroup(qm, &qg, skipQuad[q], qPair[k]);
          if (EG_validSwap(qm, qg.verts[0], qg.verts[3])== 0 ) {
              stat = EG_createQuadGroup(qm, &qg, skipQuad[q], qPair[(k + 2)%4]);
              stat = EG_createQuadGroup(qm, &qg, skipQuad[q], qPair[k]);
              if (EG_validSwap(qm, qg.verts[0], qg.verts[3]) == 0) continue;
          }
          if(qg.verts[1] == vq[i] || qg.verts[4] == vq[i] ) j = 1;
          else j = 2;
#ifdef DEBUG
          printf(" ATTENTION: WE ARE Swaping quads %d %d thru corner SW %d\n ", qg.q[0], qg.q[1], j);
          printf(" FIXING MESH : SWAP QUAD NOT CORNERN TYPE!\n");
          printQuadGroup(qm, qg);
          gnuData(qm, NULL);
#endif
          stat = EG_swappingOperation(qm, qg, j, &s1);
          if ( stat != EGADS_SUCCESS) return stat;
      } else {
          k  = j-1;
          n0 = qm->valence[k][3];
          i  = EG_quadVertIdx(qm, skipQuad[q], n0);
#ifdef DEBUG
          printf(" k %d ---> n0 %d  QIDX %d\n", k, n0, i);
#endif
          if (i >= 0) {
              s1   = k + 1;
              s2   = qm->qIdx [ 4 * (skipQuad[q] -1) + (i + 2) %4];
              stat = EG_adjQtoPair(qm, skipQuad[q], s1, s2, qPair);
#ifdef DEBUG
              printf(" Adjacent to pair %d %d from %d is %d\n",
                     s1, s2, skipQuad[q], qPair[1]);
#endif
              if (qPair[1] == -1 ||
                  qm -> qAdj[4 *(qPair[1] - 1)    ] == -1 ||
                  qm -> qAdj[4 *(qPair[1] - 1) + 1] == -1 ||
                  qm -> qAdj[4 *(qPair[1] - 1) + 2] == -1 ||
                  qm -> qAdj[4 *(qPair[1] - 1) + 3] == -1) {
                  s1 = EG_quadVertIdx(qm, skipQuad[q], j);
                  s1 = qm -> qIdx[4 * (skipQuad[q] -1) + (s1 + 2)%4];
	          stat = EG_adjQtoPair(qm, skipQuad[q], s1, s2, qPair);
#ifdef DEBUG
              printf(" Adjacent to pair %d %d from %d is %d\n",
                     s1, s2, skipQuad[q], qPair[1]);
#endif
	      }
	      if (qPair[1] == -1 ||
	          qm -> qAdj[4 *(qPair[1] - 1)    ] == -1 ||
	          qm -> qAdj[4 *(qPair[1] - 1) + 1] == -1 ||
	          qm -> qAdj[4 *(qPair[1] - 1) + 2] == -1 ||
	          qm -> qAdj[4 *(qPair[1] - 1) + 3] == -1) continue;
	      stat = EG_createQuadGroup(qm, &qg, skipQuad[q], qPair[1]);
	      if (EG_validSwap(qm, qg.verts[0], qg.verts[3])== 0 ) continue;
	      if (qg.verts[1] == n0 || qg.verts[4] == n0 ) i = 1;
	      else i = 2;
	  } else {
	      if (qm->star[k] == NULL || qm->vType[k] != -1) {
#ifdef DEBUG
	          printf("Star %d for insterting quad is NULL\n ", k + 1);
#endif
	          continue;
	      }
	      for ( i = 0 ; i < qm->star[k]->nQ; i++)
	        if (qm->star[k]->quads[i] == skipQuad[q]) break;
	      s2 = qm->star[k]->quads[qm->star[k]->idxQ[i + 1]];
	      if(s2 != -1) {
	          s1 = EG_quadVertIdx(qm, s2, n0);
	          if (s1 < 0 )
	            s2 = qm->star[k]->quads[qm->star[k]->idxQ[i + qm->star[k]->nQ - 1]];
	      } else s2 = qm->star[k]->quads[qm->star[k]->idxQ[i + qm->star[k]->nQ - 1]];
	      if(s2 == -1) continue;
	      s1 = EG_quadVertIdx(qm, s2, qm->valence[k][3]);
#ifdef DEBUG
	      printf(" ADJ Q %d containing Vert %d\n", s2, n0);
	      printQuad (qm, s2);
#endif
	      if ( s1 < 0 ) continue;
	      stat = EG_createQuadGroup(qm, &qg, skipQuad[q], s2);
	      if((qg.verts[1] == n0 ||
	          qg.verts[4] == n0 ) &&
	          (qm->vType[qg.verts[1] - 1] != -1 ||
	              qm->vType[qg.verts[4] - 1] != -1) ) i = 1;
	      else i = 2;
	      if ((qm->vType[qg.verts[i       ] - 1] != -1 &&
	          qm->vType[qg.verts[(i + 3)%6] - 1] != -1 ) ||
	          EG_validSwap(qm, qg.verts[0], qg.verts[3])== 0 ) continue;
	  }
#ifdef DEBUG
	  printf(" ATTENTION: CENTRE %d WE ARE Swaping quads %d %d thru corner SW %d\n ",k + 1, skipQuad[q], s2, i);
	  printQuadGroup(qm, qg);
	  gnuData(qm, NULL);
#endif
	  stat = EG_swappingOperation(qm, qg, i, &j);
	  if ( stat != EGADS_SUCCESS) {
	      EG_free(skipQuad);
	      return stat;
	  }
      }
  }
  EG_free(skipQuad);
  return EG_makeValidMesh(qm, 0, NULL, 3);
}


int EG_makeQuadTess(bodyQuad bodydata, ego *quadTess)
{
  int          i, j, npts, nt, stat, outLevel, nedges, mQ, *tris;
  const int    *ptype, *pindex, *trs, *trc;
  const double *xyzs, *ts, *uvs;
  ego          tess, newTess, obj, *edges;
  egTessel     *btess, *ntess;

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
  ntess = (egTessel *) newTess->blind;

  stat  = EG_getBodyTopos(btess->src, NULL, EDGE, &nedges, &edges);
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
      ntess->tess2d[i].tfi = 1;
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
  int          atype, alen, *senses, min, FACECHOICE = -1;
  const int    *ints;
  float        arg;
  double       box[6], size, params[3], time, frac;
  const double *reals;
  const char   *OCCrev, *string;
#ifdef REPORT
  char buffer[100];
#endif
  ego          context, tess, model, geom, *bodies, *dum;
  bodyQuad     *bodydata;

  start_t = clock();
  if (argc != 3 && argc != 2 && argc != 5 & argc != 6) {
      printf("\n Usage: %s = (1) filename, [Tesselation:  (2) angle (3) maxlen and chord (4)] [facechoice] \n\n", argv[0]);
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


  params[0] = 0.025;
  params[1] = 0.001;
  params[2] = 15.0;
  if (argc > 3) {
      sscanf(argv[2], "%f", &arg);
      params[2] = arg;
      sscanf(argv[3], "%f", &arg);
      params[0] = arg;
      sscanf(argv[4], "%f", &arg);
      params[1] = arg;
      printf(" Using angle = %lf,  relSide = %lf,  relSag = %lf\n",
	     params[2], params[0], params[1]);

  }
  params[0] *= size;
  params[1] *= size;
  if (argc == 3) sscanf(argv[2], "%d", &FACECHOICE);
  if (argc == 6) sscanf(argv[5], "%d", &FACECHOICE);

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
      if (stat != EGADS_SUCCESS) {
	  printf(" Tessellation is NOT Quadded!\n");
	  continue;
      }
      if (atype != ATTRSTRING) {
	  printf(" Tessellation Flag is the Wrong Type!\n");
	  continue;
      }
      if (strcmp(string, "Quad") != 0) {
	  printf(" Tessellation Flag is Not Quad = %s!\n", string);
	  continue;
      }
      stat = EG_createMeshMap(&bodydata[iBody], 0);
      if (stat != EGADS_SUCCESS) {
#ifdef REPORT
	  for (f = 0; f < bodydata[iBody].nfaces; ++f) {
	      if (FACECHOICE >= 0) f = FACECHOICE - 1;
	      meshStats(bodydata->qm[f], 0);
	      sprintf(buffer,"gnuInit_%i.txt",f + 1);
	      gnuData(bodydata->qm[f] , buffer);
	      sprintf(buffer,"wvsInit_%i.txt",f + 1);
	      wvsData(bodydata->qm[f] , buffer);
	  }
#endif
	  goto cleanup;
      }
      for (f = 0; f < bodydata[iBody].nfaces; ++f) {
	  if (FACECHOICE >= 0) f = FACECHOICE - 1;
#ifdef REPORT
	  meshStats(bodydata->qm[f], 0);
	  sprintf(buffer,"gnuInit_%i.txt",f + 1);
	  gnuData(bodydata->qm[f] , buffer);
	  sprintf(buffer,"wvsInit_%i.txt",f + 1);
	  wvsData(bodydata->qm[f] , buffer);
	  printf("\n\n======================================\n");
	  printf(" Mesh Regularization FACE %d \n ", f + 1);
	  printf("======================================\n");
#endif
	  stat      = EG_meshRegularization(bodydata[iBody].qm[f]);
	  if (stat != EGADS_SUCCESS) {
	      fprintf(stderr, " EG_meshRegularization face %d / %d = %d !!!\n ",
		      f + 1, bodydata[iBody].nfaces,  stat);
	      break;
	  }
#ifdef REPORT
	  meshStats(bodydata[iBody].qm[f], 1);
	  snprintf(buffer,100,"gnuFinal_%i.txt", f+1);
	  gnuData(bodydata->qm[f] , buffer);
	  snprintf(buffer,100,"wvsFinal_%i.txt",f+1);
	  wvsData(bodydata->qm[f], buffer);
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

  cleanup:
  if (stat != EGADS_SUCCESS ) printf(" EG_main stat in cleanup %d !!\n", stat);
  for (iBody = 0; iBody < nbody; iBody++) {
      EG_destroyMeshMap(&bodydata[iBody]);
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
