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
#include "espQuads.h"

//#define DEBUG
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
      if (dist > EPS11) {
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
    int  i, j, k, n;
    FILE *fil = NULL;
    char aux[100];

    if (qm->fID == 0 ) return;
    if (buffer == NULL) {
        snprintf(aux, 100, "wvsFace_%d_M_%d", qm->fID, qm->plotcount++);
        fil = fopen(aux, "w");
    } else {
        fil = fopen (buffer, "w");
        printf(" Writing in File %s  \n", buffer);
    }
    if (fil == NULL) {
        printf(" I couldn't write in wvsData file !!!!\n ");
        return ;
    }
    fprintf(fil, "%d %d\n", qm->totV, qm->totQ);
    for (i = 0; i < qm->totV; i++) {
        fprintf(fil, "%lf %lf %lf \n", qm->xyzs[3 * i    ],
                qm->xyzs[3 * i + 1], qm->xyzs[3 * i + 2]);
    }
    fprintf(fil,"\n");
    for (n = i = 0; i < qm->totQ; i++) {
        fprintf(fil, "%d %d %d %d\n",  qm->qIdx[4 * i    ],
                qm->qIdx[4 * i + 1], qm->qIdx[4 * i + 2],
                qm->qIdx[4 * i + 3]);
        if (qm->qAdj[4 * i    ] == -1 ) n++;
        if (qm->qAdj[4 * i + 1] == -1 ) n++;
        if (qm->qAdj[4 * i + 2] == -1 ) n++;
        if (qm->qAdj[4 * i + 3] == -1 ) n++;
    }
    fprintf(fil,"%d\n", n );
    for (k = i = 0; i < qm->totQ; i++) {
        for (j = 0; j < 4; j++) {
            if (qm->qAdj[4 * i + j] == -1 ) {
                fprintf(fil, "%lf %lf %lf\n",
                        qm->xyzs[3 * (qm->qIdx[4 * i + j]- 1)    ],
                        qm->xyzs[3 * (qm->qIdx[4 * i + j]- 1) + 1],
                        qm->xyzs[3 * (qm->qIdx[4 * i + j]- 1) + 2]);
                fprintf(fil, "%lf %lf %lf\n",
                        qm->xyzs[3 * (qm->qIdx[4 * i + (j + 1)%4]- 1)    ],
                        qm->xyzs[3 * (qm->qIdx[4 * i + (j + 1)%4]- 1) + 1],
                        qm->xyzs[3 * (qm->qIdx[4 * i + (j + 1)%4]- 1) + 2]);
                k++;
            }
        }
    }
    if (k != n)
      printf(" wvsData Error filling the edges segments!! %d != %d\n", k, n);
    fclose(fil);
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
  int   i, i4[10], i3[10], i1[100], j;
  vStar *star = NULL;

  if (vID <= 0 || vID > qm->totV) return EGADS_INDEXERR;
  i      = EG_buildStar(qm, &star, vID);
  if (i != EGADS_SUCCESS || star == NULL) {
      printf(" In setValence for vertex %d stat from EG_buildStar = %d\n", vID, i);
      return i;
  }
  if (qm->valSize[vID -1] < star->nQ ) {
      qm->valSize[vID -1] = star->nQ;
      qm->valence[vID -1] = EG_reall (qm->valence[vID -1], (3 + star->nQ ) * sizeof(int) );
      if (qm->valence[vID -1] == NULL ) return EGADS_MALLOC;
  }
  qm->valence[vID - 1][2] = star->nQ; // real valence
  qm->valence[vID - 1][1] = qm->valence[vID - 1][2];
  if (qm->vType[vID - 1] == 2 ||
      qm->vType[vID - 1] == 0) qm->valence[vID - 1][1] = star->nQ + 2; //adjusted based on vertex type, eg, a corner is regular when val = 2
  else if (qm->vType[vID - 1] == 3) qm->valence[vID -1][1] = star->nQ + 1;
  else if (qm->vType[vID - 1] >= 5) qm->valence[vID -1][1] = star->nQ - 1;
  for (i3[0] = i4[0] = i1[0] = i = 0; i < star->nQ; i++) {
      if (star->verts[2 * i + 1] == -1) continue;
      if      (qm->vType[star->verts[2 * i + 1] - 1] >= 4)
    i4[++i4[0]] = star->verts[2 * i + 1];
      else if (qm->vType[star->verts[2 * i + 1] - 1] >= 2)
    i3[++i3[0]] = star->verts[2 * i + 1];
      else
    i1[++i1[0]] = star->verts[2 * i + 1];
  }
  for (j = i = 0; i < i4[0]; i++)
    qm->valence[vID - 1][3 + (j++)] = i4[1 + i];
  for (i = 0; i < i3[0]; i++)
    qm->valence[vID - 1][3 + (j++)] = i3[1 + i];
  for (i = 0; i < i1[0]; i++)
    qm->valence[vID - 1][3 + (j++)] = i1[1 + i];
  EG_freeStar(&star);

  return EGADS_SUCCESS;
}


static int EG_backupQuads(meshMap *qm, int *nq, int *qlist, Quad **quad)
{
  int   i, j, q, v, qcount, *qaux = NULL, k, add = 0, stat;
  vStar *star = NULL;
  
  qaux      = EG_alloc(qm->totQ * sizeof(int));
  if (qaux == NULL) return EGADS_MALLOC;
  for (qcount = q = 0; q < *nq; q++) {
    if (qlist[q] == -1) continue;
    for (i = 0; i < 4; i++) {
      v    = qm->qIdx[4 * (qlist[q] - 1) + i];
      stat = EG_buildStar(qm, &star, v);
      if (star == NULL || stat != EGADS_SUCCESS) {
        EG_free (qaux);
        return stat;
      }
      for (j = 0; j < star->nQ; j++) {
        if (star->quads[j] == -1) continue;
        for (add = k = 0; k < qcount; k++) {
          if (qaux[k] == star->quads[j]) {
            add = 1;
            break;
          }
        }
        if (add == 0 || qcount == 0) qaux[qcount++] = star->quads[j];
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
      EG_free (vid);
      return stat;
      }
  }
  EG_free (vid);
  
  return EGADS_SUCCESS;
}

static int EG_centroid(meshMap *qm, int n, int *list, double *quv, int usequv)
{
  int i, j = 0;
  double xyz0[3], xyz1[3], pos[18];

  if (usequv == 1) {
      i = EG_evaluate(qm->face, quv, pos);
      if ( i != EGADS_SUCCESS ) return i;
      xyz0[0] = pos[0]; xyz0[1] = pos[1]; xyz0[2] = pos[2];
      j = 1;
  }
  else {
      quv [0] = quv [1] = 0.0;
      xyz0[0] = xyz0[1] = xyz0[2] = 0.0;
  }
  for (i  = 0; i < n; i++) {
      xyz0[0] += qm->xyzs[3 * (list[i] - 1)    ];
      xyz0[1] += qm->xyzs[3 * (list[i] - 1) + 1];
      xyz0[2] += qm->xyzs[3 * (list[i] - 1) + 2];
      quv[0]  += qm->uvs [2 * (list[i] - 1)    ];
      quv[1]  += qm->uvs [2 * (list[i] - 1) + 1];
      j++;
  }
  xyz0[0] /= (double) j;
  xyz0[1] /= (double) j;
  xyz0[2] /= (double) j;
  quv [0] /= (double) j;
  quv [1] /= (double) j;
  i        = EG_invEvaluateGuess(qm->face, xyz0, quv, xyz1);
  if (i   != EGADS_SUCCESS    || quv[0] < qm->range[0] ||
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


static int EG_projectToTangentPlane(double normal[], double *nO, double *p,
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
  if (fabs(dotNP - c) < EPS11 || sqrt(dist) < EPS11) {
      proj[0] = p[0]; proj[1] = p[1]; proj[2] = p[2];
      return EGADS_SUCCESS;
  }
  lambda  = (c - dotNP);
  proj[0] = p[0] + lambda * normal[0];
  proj[1] = p[1] + lambda * normal[1];
  proj[2] = p[2] + lambda * normal[2];
  dist    = DOT(normal, proj);
  if (fabs(dist - c) < EPS11) return EGADS_SUCCESS;
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
  if (dot1 < EPS11 || dot2 < EPS11) {
      quv[0] = qm->uvs[ 2 * (vC - 1)    ];
      quv[1] = qm->uvs[ 2 * (vC - 1) + 1];
      stat = EG_centroid(qm, qm->valence[vC -1][2], &qm->valence[vC -1][3], quv, 1);
      if (stat != EGADS_SUCCESS) return stat;
      stat = EG_evaluate(qm->face, quv, centre);
      if (stat != EGADS_SUCCESS) return stat;
      dot1 = centre[3]*centre[3] + centre[4]*centre[4] + centre[5]*centre[5];
      dot2 = centre[6]*centre[6] + centre[7]*centre[7] + centre[8]*centre[8];
      if (dot1 < EPS11 || dot2 < EPS11) return EGADS_GEOMERR;
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


static int EG_angleAtBoundaryVertex(meshMap *qm, int v, int *links, double *size)
{
  int   i, j, k;
  vStar *star = NULL;

  *size = 0.0;
  if (qm->vType[v - 1] < 0) return EGADS_INDEXERR;
  i = EG_buildStar(qm, &star, v);
  if (i != EGADS_SUCCESS || star == NULL) {
      printf(" Looking at corners: buildstar %d is %d \n ", v, i);
      return i;
  }
  for (links[0] = links[1] = k = i = 0; i < star->nQ; i++) {
      j = star->verts[2 * i + 1] - 1;
      if (qm->vType[j] != -1) k++;
      if (star->quads[i] == -1) {
      links[1] = star->verts[2 * i + 1];
      links[0] = star->verts[star->idxV[2 * i + 3]];
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


/*
 * Returns <0 EGADS_ERROR
 * = 0 self intersects
 * = 1 Good area
 * = 2 contains obtuse angles ( if quad is boundary obtuse angles => area = 0 )
 */
static int EG_quadAngleOrientation(meshMap *qm, int qID, int *ori, double *theta )
{
    int i, k, k1, stat, vA, vB, vC, count, doublet,   qV[4] ;
    int bID = -1, area = 1, selfint = 1, sharpBound = 0 ;
    int piv[6] = {1, 2, 2, 3, 1, 3};
    double projABCD[12], cross[3], qNormal[3], vAB[3], vAC[3], quv[2], xyz[18];
    double norm1, norm2, angle, dot, c, dotNP;

    //printf(" QUAD %d -----\n ", qID);
    qV[0]  = qm->qIdx[4*(qID - 1)    ] - 1;
    qV[1]  = qm->qIdx[4*(qID - 1) + 1] - 1;
    qV[2]  = qm->qIdx[4*(qID - 1) + 2] - 1;
    qV[3]  = qm->qIdx[4*(qID - 1) + 3] - 1;
    quv[0] = 0.25 * (qm->uvs[2 * qV[0]    ] + qm->uvs[2 * qV[1]    ] +
                     qm->uvs[2 * qV[2]    ] + qm->uvs[2 * qV[3]    ]);
    quv[1] = 0.25 * (qm->uvs[2 * qV[0] + 1] + qm->uvs[2 * qV[1] + 1] +
                     qm->uvs[2 * qV[2] + 1] + qm->uvs[2 * qV[3] + 1]);
    stat      = EG_evaluate(qm->face, quv, xyz);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_quadAngleOrientation:: EG_evaluate for quad %d  is %d !!\n",
             qID, stat);
      return stat;
    }
    norm1 = xyz[3] * xyz[3] + xyz[4] * xyz[4] + xyz[5] * xyz[5];
    norm2 = xyz[6] * xyz[6] + xyz[7] * xyz[7] + xyz[8] * xyz[8];
    if (norm1 < EPS11 || norm2 < EPS11) {
      stat      = EG_centroid(qm, 4, &qm->qIdx[ 4 * (qID - 1)], quv, 0);
      stat     += EG_evaluate(qm->face, quv, xyz);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_quadAngleOrientation:: EG_evaluate for quad %d  is %d !!\n",
               qID, stat);
        return stat;
      }
      norm1 = xyz[3] * xyz[3] + xyz[4] * xyz[4] + xyz[5] * xyz[5];
      norm2 = xyz[6] * xyz[6] + xyz[7] * xyz[7] + xyz[8] * xyz[8];
      if (norm1 < EPS11 || norm2 < EPS11) return EGADS_DEGEN;
    }
    vAB[0] = xyz[3]; vAB[1] = xyz[4]; vAB[2] = xyz[5];
    vAC[0] = xyz[6]; vAC[1] = xyz[7]; vAC[2] = xyz[8];
    if (qm->face->mtype == SREVERSE) {
      CROSS(vAC, vAB, qNormal);
    } else {
      CROSS(vAB, vAC, qNormal);
    }
    norm1       = sqrt (qNormal[0] * qNormal[0] + qNormal[1] * qNormal[1] +
                        qNormal[2] * qNormal[2]);
    qNormal[0] /= norm1;
    qNormal[1] /= norm1;
    qNormal[2] /= norm1;
    c           = DOT(qNormal, xyz);
    for (i      = 0; i < 4; i++) {
      dotNP               = qNormal[0] * qm->xyzs[3 * qV[i]    ] +
                            qNormal[1] * qm->xyzs[3 * qV[i] + 1] +
                            qNormal[2] * qm->xyzs[3 * qV[i] + 2];
      angle               = (c - dotNP);
      projABCD[3 * i    ] = qm->xyzs[3 * qV[i]    ] + angle * qNormal[0];
      projABCD[3 * i + 1] = qm->xyzs[3 * qV[i] + 1] + angle * qNormal[1];
      projABCD[3 * i + 2] = qm->xyzs[3 * qV[i] + 2] + angle * qNormal[2];
    }
    for (doublet = count = k = 0; k < 4; ++k) {
      ori[k] = 1;
      vA     = k;
      if (qm->valence[qV[vA]][2] ==  2 &&
          qm->vType  [qV[vA]]    == -1) {
        ori[k]   = 1;
        theta[k] = 0.5 * PI;
        doublet  = 1;
      } else {
        if (qm->vType[qV[vA]] >= 4 ) {
          bID        = k ;
          sharpBound = 1;
        }
        for (k1  = 0; k1 <= 2; k1++) {
          vB     = (k + piv[2 * k1    ])%4;
          vC     = (k + piv[2 * k1 + 1])%4;
          vAB[0] = projABCD[3 * vB    ] - projABCD[3 * vA    ];
          vAB[1] = projABCD[3 * vB + 1] - projABCD[3 * vA + 1];
          vAB[2] = projABCD[3 * vB + 2] - projABCD[3 * vA + 2];
          vAC[0] = projABCD[3 * vC    ] - projABCD[3 * vA    ];
          vAC[1] = projABCD[3 * vC + 1] - projABCD[3 * vA + 1];
          vAC[2] = projABCD[3 * vC + 2] - projABCD[3 * vA + 2];
          CROSS(vAB, vAC, cross);
          if ( k1 != 2 ) {
            if (DOT(qNormal, cross) < 0) ori[k] = -1;
            continue;
          }
          norm1  = vAB[0] * vAB[0] + vAB[1] * vAB[1] + vAB[2] * vAB[2];
          norm2  = vAC[0] * vAC[0] + vAC[1] * vAC[1] + vAC[2] * vAC[2];
          if (norm1 < EPS11 || norm2 < EPS11) {
              ori[k] = -1;
              dot    = 1.0;
          } else dot = DOT(vAB, vAC) / sqrt(norm1 * norm2);
          if      (fabs(dot - 1.0) < EPS11) theta[k] = 0.0;
          else if (fabs(dot + 1.0) < EPS11) theta[k] = PI;
          else                              theta[k] = acos(dot);
          if (DOT(qNormal, cross) < 0) theta[k]  = (2.0 * PI - theta[k]);
        }
      }
      if (ori[k] == -1) continue;
      selfint = 0;
      count++;
    }
    if      (selfint == 1) area = 0;
    else if (count   != 4) area = 2;
    if (sharpBound == 1 && area > 0 ){
    if       (theta[bID] > PI    ) area = 0;
    else if ( theta[bID] > DEG170) area = 2;
    }

    else if (doublet == 1 && area != 0 ) area = 1;
#ifdef DEBUG
    if (area != 1 ) {
        if ( area == 2 )
      printf(" ************ ATENCION QUAD %d DEGENERATED QUAD ( AREA 2 )       ***************\n", qID);
        else if (sharpBound == 2 )
      printf(" ************ ATENCION QUAD %d HAS A SHARP INVALID CORNER ( AREA 0 )   ***************\n", qID);
        else
      printf(" ************ ATENCION QUAD %d SELFINTERSECTING QUAD ( AREA 0 )   ***************\n", qID);
        printf("------------    Internal angles --------------\n");
        for (k = 0; k < 4; ++k)
      printf("Vertex %d (T %d val %d) has angle %f and orientation %d \n ",
             qV[k] + 1,qm->vType[qV[k]], qm->valence[qV[k]][2], theta[k], ori[k]);
        printQuad(qm , qID);
        for (k = 0; k < 4; ++k)
      printf("%lf %lf %lf %d\n ", projABCD[3 * k], projABCD[3 * k + 1], projABCD[3 * k + 2], qV[k] + 1);
        k = 0;
        printf("%lf %lf %lf %d\n ", projABCD[3 * k], projABCD[3 * k + 1], projABCD[3 * k + 2], qV[k] + 1);
    }
#endif
    return area;
}

static void EG_computeCoords(meshMap *qm, vStar *star, int weight )
{
  int    i, j, bdry, k, kk, i0, i1, vID, ori[4];
  double uvc[2], angles[4], uva[2];
  vID = star->verts[0] - 1;
  if (qm->vType[vID] != -1 ) {
#ifdef DEBUG
      printf(" VERTEX %d is boundary !\n ", vID + 1);
#endif
      return;
  }

#ifdef DEBUG
  if ( weight != -1 )
  printf(" AVERAGE COORDS FOR VERTEX %d WEIGHT %d\n ", vID + 1, star->verts[weight]);
  else
      printf(" AVERAGE COORDS FOR VERTEX %d CENTROID NO BIAS\n ", vID + 1);
  printVertex(qm, vID + 1);
#endif
  if (qm->valence[vID][2] == 2) {
      i0     = qm->valence[vID][3] - 1;
      i1     = qm->valence[vID][4] - 1;
      uva[0] = 0.5 * (qm->uvs[2 * i0    ] + qm->uvs[2 * i1    ]);
      uva[1] = 0.5 * (qm->uvs[2 * i0 + 1] + qm->uvs[2 * i1 + 1]);
      updateVertex(qm, vID + 1, uva);
#ifdef DEBUG
      printf(" LEAVE WITH COORDS\n");
      printVertex(qm, vID + 1);
#endif
      return;
  }
  for (bdry = i = 0 ; i < star->nQ; i++) {
      if ( qm->vType[star->verts[2 * i + 1] - 1] != -1 )
      bdry++;
  }
  if (bdry >= 2 && star->nQ > 3 ) {
      for (i0 = i1 = k = 0; k < star->nQ; k++) {
      i0 = star->verts[2 * k + 1] - 1;
      if (qm->vType[i0] == -1 ) continue;
      i1 = star->verts[star->idxV[2 * k + 5]] - 1;
      j  = star->verts[star->idxV[2 * k + 7]] - 1;
      if (qm->vType[i1] == -1 && j != i0 ) continue;
      uva[0] = 0.5 * ( qm->uvs[2 * i0    ] + qm->uvs[2 * i1    ] );
      uva[1] = 0.5 * ( qm->uvs[2 * i0 + 1] + qm->uvs[2 * i1 + 1] );
      uvc[0] = qm->uvs[2 * vID    ];
      uvc[1] = qm->uvs[2 * vID + 1];
      updateVertex (qm, vID + 1, uva);
#ifdef DEBUG
      printf(" MID POINTS BETWEENM BOUNDS %d  %d \n ", i0 + 1, i1 + 1);
      for (i = 0; i < star->nQ; i++)
        printQuad(qm, star->quads[i]);
#endif
      for (j = i = 0; i < star->nQ; i++) {
          EG_quadAngleOrientation(qm, star->quads[i], ori, angles);
          for ( kk = 0 ; kk < 4; kk++ ) {
          if ( angles[kk] > PI || ori[kk] == -1 ) {
#ifdef DEBUG
              printf(" ANGLE %d is %lf \n ", qm->qIdx[ 4 * ( star->quads[i] - 1 ) + kk], angles[kk] );
#endif
              j = 1;
              break;
          }
          }
          if ( j != 0 ) break;
      }
      if (j == 0 ) {
#ifdef DEBUG
          printf(" VERTEX IS GOOD \n ");
#endif
          return;
      }
      updateVertex (qm, vID + 1, uvc);
      }
  }
  //uva[0] = uva[1] = 0.0;
  //i = EG_centroid(qm, star->nV, star->verts, uva, 0 );
  uva[0] = qm->uvs[2 * vID    ];
  uva[1] = qm->uvs[2 * vID + 1];
  i      = EG_centroid(qm, qm->valence[vID][2],
                          &qm->valence[vID][3], uva, 1) ;
  if ( i != EGADS_SUCCESS ) return;
  if ( star->nQ > 5 && weight == -1 ) weight = 0;
  if ( weight != -1 ) {
      i = star->verts[weight] - 1;
      uva[0] = 0.5 * (qm->uvs[2 * i    ] + uva[0]);
      uva[1] = 0.5 * (qm->uvs[2 * i + 1] + uva[1]);
  }
  if (uva[0] >= qm->range[0] && uva[0] <= qm->range[1] &&
      uva[1] >= qm->range[2] && uva[1] <= qm->range[3]) {
      updateVertex(qm, vID + 1, uva);
  }
#ifdef DEBUG
  printf(" LEAVE WITH COORDS\n");
  printVertex(qm, vID + 1);
#endif
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
  vStar        *star = NULL;

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
    bodydata->qm[f]->xyzs    = (double *) EG_alloc(3*(2 * len  )*sizeof(double ));
    bodydata->qm[f]->uvs     = (double *) EG_alloc(2*(2 * len  )*sizeof(double ));
    bodydata->qm[f]->vType   = (int    *) EG_alloc(  (2 * len  )*sizeof(   int ));
    bodydata->qm[f]->qIdx    = (int    *) EG_alloc(4*(2 * nquad)*sizeof(   int ));
    bodydata->qm[f]->qAdj    = (int    *) EG_alloc(4*(2 * nquad)*sizeof(   int ));
    bodydata->qm[f]->remQ    = (int    *) EG_alloc(  (2 * nquad)*sizeof(   int ));
    bodydata->qm[f]->remV    = (int    *) EG_alloc(  (2 * len  )*sizeof(   int ));
    bodydata->qm[f]->valence = (int   **) EG_alloc(  (2 * len  )*sizeof(   int *));
    bodydata->qm[f]->valSize = (int    *) EG_alloc(  (2 * len  )*sizeof(   int ));
    if (bodydata->qm[f]->qIdx    == NULL || bodydata->qm[f]->qAdj    == NULL ||
        bodydata->qm[f]->xyzs    == NULL || bodydata->qm[f]->uvs     == NULL ||
        bodydata->qm[f]->vType   == NULL || bodydata->qm[f]->remQ    == NULL ||
        bodydata->qm[f]->remV    == NULL || bodydata->qm[f]->valence == NULL ||
    bodydata->qm[f]->valSize == NULL ) {
      stat2 = EGADS_MALLOC;
      continue;
    }
    bodydata->qm[f]->remQ[0]  = 0;
    bodydata->qm[f]->remV[0]  = 0;
    bodydata->qm[f]->invsteps = 0;
    for (j = 0; j < 2 * len; j++) {
        bodydata->qm[f]->valence[j] = (int *) EG_alloc(sizeof(int));
        bodydata->qm[f]->valSize[j] = 1;
        if (bodydata->qm[f]->valence[j] == NULL) {
            stat2 = EGADS_MALLOC;
            continue;
        }
    }
    bodydata->qm[f]->face = bodydata->faces[f];
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
      stat = EG_angleAtBoundaryVertex(bodydata->qm[f], j + 1, e4, &angle);
      if (stat != EGADS_SUCCESS || angle < EPS11) {
        stat2   = EGADS_GEOMERR;
        printf(" Stat in EG_angleAtBoundaryVertex %d angle %f\n ", stat, angle);
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
      stat = EG_buildStar(bodydata->qm[f], &star, j + 1);
      if (stat != EGADS_SUCCESS || star == NULL) {
        stat2 = stat;
        continue;
      }
      EG_freeStar(&star);
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
  int    v, q,  i, j, k, it = 0, itMax,  sum = 0, ori[4], *qa = NULL, ta[50],taux[50],  m, a1, m2;
  int     *verts = NULL, recover = 0, stat = EGADS_SUCCESS, wi;
  double  angles[4], *uvxyz = NULL, uv[8];
#ifdef DEBUG
  char   buffer[100];
  int qi;
#endif
  vStar  **star = NULL;

  uvxyz = (double*)EG_alloc(5 * qm->totV * sizeof(double));
  star  = (vStar**)EG_alloc(    qm->totV * sizeof(vStar*));
  qa    = (int   *)EG_alloc(    qm->totQ * sizeof(int));
  if (uvxyz == NULL || star  == NULL || qa == NULL) return EGADS_MALLOC;
  for (j = 0; j < qm->totV; j++) {
      star[j] = NULL;
      if (qm->vType[j] != -1) continue;
      uvxyz[5 * j    ] = qm->uvs [2 * j    ];
      uvxyz[5 * j + 1] = qm->uvs [2 * j + 1];
      uvxyz[5 * j + 2] = qm->xyzs[3 * j    ];
      uvxyz[5 * j + 3] = qm->xyzs[3 * j + 1];
      uvxyz[5 * j + 4] = qm->xyzs[3 * j + 2];
  }
  for ( j = 0 ; j < qm->totQ; j++ ) qa[j] = -1;
#ifdef DEBUG
  printf(" Make valid mesh for %d points \n ", nP);
  if ( nP > 0 )
    for ( i = 0 ; i < nP; i++ ) printf(" list %d = %d \n ", i, pList[i]);
  gnuData(qm, NULL);
#endif
  recover = 1;
  if (fullReg == 0) { // move around only affected vertices
      if (nP  == 0 || pList == NULL) {
      stat = EGADS_SUCCESS;
      goto cleanup;
      }
      for (j = 0; j < nP; j++) {
      if (qm->vType[pList[j] - 1] == -2) continue;
      stat      = EG_buildStar(qm, &star[pList[j] - 1], pList[j]);
      if (stat != EGADS_SUCCESS ||  star[pList[j] - 1] == NULL) goto cleanup;
      for (i    = 0; i < star[pList[j] -1]->nQ; i++) {
          if (star[pList[j] -1]->quads[i] == -1) continue;
          v = star[pList[j] -1]->verts[2 * i + 1] - 1;
          if (qm->vType[v] == -1 && star [v] == NULL ) {
          stat      = EG_buildStar(qm, &star[v], v + 1);
          if (stat != EGADS_SUCCESS ||  star[v] == NULL) goto cleanup;
          }
          v = star[pList[j] -1]->verts[2 * i + 2] - 1;
          if (qm->vType[v] == -1 && star [v] == NULL ) {
          stat      = EG_buildStar(qm, &star[v], v + 1);
          if (stat != EGADS_SUCCESS ||  star[v] == NULL) goto cleanup;
          }
      }
      }
      itMax = 5;
  } else {
      for (i = 0 ; i < qm->totV; i++ ) {
      if ( qm->vType[i] != -1 ) continue;
      stat        = EG_buildStar(qm, &star[i], i + 1);
      if (stat   != EGADS_SUCCESS || star[i] == NULL) goto cleanup;
      }
      itMax = 50;
  }
  for (i = 0; i < qm->totV; i++) {
      if ( qm->vType[i] != -1 || star[i] == NULL ) continue;
      for (j = 0; j < star[i]->nQ; j++) {
      if (   star[i]->quads[j] == -1) continue;
      if (qa[star[i]->quads[j]    -1] == -1) {
          qa[star[i]->quads[j]    -1] = EG_quadAngleOrientation(qm, star[i]->quads[j], ori, angles);
      }
      }
  }
  for (it = 0 ; it <  itMax; it++) {
#ifdef DEBUG
      printf(" MAKEVALID IT %d ===================  \n ", it );
      gnuData(qm, NULL);
#endif
      for (q = 0 ; q < qm->totQ; q++ ) {
        if (qa[q] == -1 || (it == 0 && qa[q] == 1)) continue;
      for (i = 0; i < 4; i++) {
          v  = qm->qIdx [4 * q + i ] - 1;
          if ( star[v] == NULL || qm->vType[v] != -1) continue;
          uv[0]   = qm->uvs [ 2 * v    ];
          uv[1]   = qm->uvs [ 2 * v + 1];
          for (m2 = k = 0; k <= star[v]->nQ; k++ ) {
              wi  = -1 + 2 * k;
              if (k > 0 && qm->vType[star[v]->verts[wi] -1] != -1 )  continue;
              EG_computeCoords (qm, star[v], wi);
              for (a1 = m = j = 0; j < star[v]->nQ; j++ ) {
                  ta[j] = EG_quadAngleOrientation(qm, star[v]->quads[j], ori, angles);
                  if ( ta[j] == 1 ) {
                      a1++;
                      continue;
                  }
                  else if ((qa[star[v]->quads[j] - 1] == 0 && ta[j] > 0 ) ||
                            qa[star[v]->quads[j] - 1] ==      ta[j] ) continue;
#ifdef DEBUG
                  printf(" RESET !!! AREA %d was %d now %d \n ", star[v]->quads[j],
                          qa[star[v]->quads[j] - 1], ta[j] );
                  for ( qi = 0 ; qi < star[v]->nQ; qi++ ) {
                      printQuad(qm, star[v]->quads[qi]);
                  }
#endif
                  updateVertex(qm, v + 1, uv) ;
                  m    = 1;
                  break;
              }
              if (a1 == star[v]->nQ ) {
#ifdef DEBUG
                  printf(" GOOD MOVE V %d -. ALL AREAS ARE 1 \n ", v+1);
#endif
                  m2 = 1;
                  break;
              }
              else if (m == 0 && m2 == 0) {
                  uv[2] = qm->uvs[2 * v    ];
                  uv[3] = qm->uvs[2 * v + 1];
                  for ( j = 0 ; j < star[v]->nQ; j++ )
                      taux[j] = ta[j];
                 m2    = 2;
              }
          }
          if (m2 == 1 )
              for ( j = 0 ; j < star[v]->nQ; j++ )
                  qa[star[v]->quads[j] -1 ] = ta[j];
          else if ( m2 == 2 ) {
              if (m2 == 2) updateVertex(qm, v + 1, &uv[2]) ;
              for ( j = 0 ; j < star[v]->nQ; j++ )
                  qa[star[v]->quads[j] -1 ] = taux[j];
          }
      }
      }
      for (sum = i = 0 ; i < qm->totQ; i++ ) {
          if ( qa[i] == -1) continue;
          if ( qa[i] == 0                 ) sum = 1;
          if ( qa[i] == 2 && fullReg == 1 ) sum = 2;
#ifdef DEBUG
          printf(" QUAD %d has area %d \n ", i +1, qa[i] );
          EG_quadAngleOrientation(qm, i + 1, ori, angles);
          printQuad(qm, i + 1 );
#endif
          if ( sum != 0 ) break;
      }
#ifdef DEBUG
      printf(" ROUND %d -->sum %d \n ", it, sum );
#endif
      if (sum == 0 ) {
          recover = 0 ;
          if (fullReg != 1 ) break;
          for (j = 0; j < qm->totV; j++) {
              if (qm->vType[j] != -1) continue;
              uvxyz[5 * j    ] = qm->uvs [2 * j    ];
              uvxyz[5 * j + 1] = qm->uvs [2 * j + 1];
              uvxyz[5 * j + 2] = qm->xyzs[3 * j    ];
              uvxyz[5 * j + 3] = qm->xyzs[3 * j + 1];
              uvxyz[5 * j + 4] = qm->xyzs[3 * j + 2];
          }
      }
  }
  stat = EGADS_SUCCESS;
  if ( recover == 1 && sum == 1) {
#ifdef DEBUG
      printf(" Mesh is invalid for MIN MAX TOTALS [0, 200 deg]\n");
      snprintf(buffer, 100,"face_%d_InvalidMesh_%d",
           qm->fID, qm->plotcount);
      gnuData(qm , buffer);
#endif
      stat = EGADS_GEOMERR;
      for (j = 0; j < qm->totV; j++) {
      if (qm->vType[j] != -1) continue;
      qm->uvs  [ 2 * j    ] = uvxyz[ 5 * j    ];
      qm->uvs  [ 2 * j + 1] = uvxyz[ 5 * j + 1];
      qm->xyzs [ 3 * j    ] = uvxyz[ 5 * j + 2];
      qm->xyzs [ 3 * j + 1] = uvxyz[ 5 * j + 3];
      qm->xyzs [ 3 * j + 2] = uvxyz[ 5 * j + 4];
      }
      qm->invsteps++;
      if ( fullReg != 0 ) {
      for (j = i = 0; i < qm->totQ; i++) {
          if (qm->vType[qm->qIdx[4 * i]] == -2) continue;
          if (EG_quadAngleOrientation(qm, i + 1, ori, angles) != 1 ) {
          j = 1;
          break;
          }
          if ( j == 0 ) stat = EGADS_SUCCESS;
      }
      }
  }

#ifdef DEBUG
  printf(" LEAVING PEACEFULLY FROM OPTIMIZER ??? %d \n ", stat);
  gnuData(qm, NULL);
#endif

  cleanup:
  EG_free(uvxyz);
  EG_free(verts);
  EG_free(qa);
  for ( i = 0 ; i < qm->totV; i++ )
    if ( star[i] ) EG_freeStar(&star[i]);
  EG_free(star);
  return stat;
}


static int EG_createQuadGroup(meshMap *qm, quadGroup *qg, int q0, int q1)
{
  int i, j, k, ids[2], piv = 0, aux, vaux[6], vc[3], centre, c1 = 0, c2 = 0;

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
  for (i = 0; i < 6; i++) {
    qg->verts[i] = -1;
    qg->vals [i] = -1;
  }
  if (q0 == q1) {
      printf("EG_createQuadGroup :: Quads %d are the same!!!\n", q0);
      return EGADS_INDEXERR;
  }
  for (k = i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
        if (qm->qIdx[4 * (q0 - 1) + i] == qm->qIdx[4 * (q1 - 1) + j] &&
            qm->qIdx[4 * (q0 - 1) + i] != -1)
          vc[++k] = qm->qIdx[4 * (q1 - 1) + j];
        }
    }
  vc[0] = k;
  if (vc[0] != 2) {
    printf("EG_createQuadGroup: Quad %d %d are not adjacent!!\n", q0, q1);
#ifdef DEBUG
    printQuad(qm, q0);
    printQuad(qm, q1);
#endif
    return EGADS_INDEXERR;
  }
  centre = vc[1];
  if (qm->valence[vc[2] -1][1] >
      qm->valence[vc[1] -1][1]) centre = vc[2];
  ids[0] = EG_quadVertIdx(qm, qg->q[0], centre);
  ids[1] = EG_quadVertIdx(qm, qg->q[1], centre);
  piv    = 0;
  if (qm->qAdj[4 * (qg->q[0] - 1) + ids[0]] == qg->q[1]) piv = 1;
  for (i = 0; i < 4; i++)
    qg->verts[i] = qm->qIdx[4 * (qg->q[piv] - 1) + (ids[piv] + i)%4];
  aux = (piv + 1)%2;
  qg->verts[4] = qm->qIdx[4 * (qg->q[aux] - 1) + (ids[aux] + 2)%4];
  qg->verts[5] = qm->qIdx[4 * (qg->q[aux] - 1) + (ids[aux] + 3)%4];
  if (piv == 1) {
    aux       = qg->q[0];
    qg->q[0] = qg->q[1];
    qg->q[1] = aux;
  }
  for (i = 0; i < 6; i++) qg->vals[i] = qm->valence[qg->verts[i] -1][1];
  if (qg->vals[0] == qg->vals[3]) {
    if (qg->vals[1] != 4) c1++;
    if (qg->vals[5] != 4) c1++;
    if (qg->vals[2] != 4) c2++;
    if (qg->vals[4] != 4) c2++;
    if (c2 > c1) {
      for (i = 0; i < 6; i++) vaux[i]      = qg->verts[i];
      for (i = 0; i < 6; i++) qg->verts[i] = vaux[(i + 3)%6];
      for (i = 0; i < 6; i++) qg->vals[i]  = qm->valence[qg->verts[i] - 1][1];
      aux       = qg->q[0];
      qg->q[0] = qg->q[1];
      qg->q[1] = aux;
    }
  }
  
  return EGADS_SUCCESS;
}

/* Assuming qID is collapsing through v */
static int validCollapse(meshMap *qm, int qID, int v)
{
  int j, k, kk, id, link, aux, aux2, type, val;
  
  if (v <= 0 || v > qm->totV) return EGADS_INDEXERR;
  id     = EG_quadVertIdx(qm, qID, v);
  aux    = qm->qIdx[4 * (qID - 1) + (id + 2)%4] - 1;
  if (qm->vType[aux] != -1) {
    for (j = 0; j < qm->valence[v - 1][2]; j++) {
      k  = qm->valence[v - 1][3 + j] -1;
      if (qm->vType[k] == -1) continue;
      if (EG_quadVertIdx(qm, qID, k + 1) < 0) {
#ifdef DEBUG
        printf(" Vertex %d is bounds and not in quad %d\n ", k + 1, qID );
#endif
        return 0;
      }
    }
  }
  for (j   = 0; j < 2; j++) {
    link = qm->qIdx[4 * (qID - 1) + (id + 2 * j + 1)%4];
    type = qm->vType  [link - 1];
    val  = qm->valence[link - 1][2];
    if (type > 0 ) {
    if ( val <= type) {
#ifdef DEBUG
        printf(" link %d has valence %d type %d !\n ", link, val, type);
#endif
        return 0;
    }else continue;
    }
    else if (val != 3) continue;
    for (k  = 0; k < 3; k++) {
      aux = qm->valence[link - 1][3 + k] - 1;
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
        if (aux2 + 1 == link) continue;
        if (qm->valence[aux2][2] == qm->vType[aux] ||
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
static int validSwap(meshMap *qm, int v1, int v2)
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
  int   nq, stat, i0, i1, i, j, adj, idx[4], qID[2], adjQmap[6];
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
  idx[0] = 0; idx[1] = 3; idx[2] = swap; idx[3] = swap + 3;
  for (i = 0; i < 4; i++) {
      stat = setValence(qm, qg.verts[idx[i]]);
      if (stat != EGADS_SUCCESS) {
      printf(" Inside swapping operation set valence for %d -->%d !!!\n ",
         qg.verts[idx[i]], stat);
      EG_free(quad);
      return stat;
      }
  }
  *activity = 1;
  if (EG_makeValidMesh(qm, 6, qg.verts, 0) == EGADS_SUCCESS) {
      EG_free(quad);
      return EGADS_SUCCESS;
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
  int   id0 = -1, id1 = -1, dist, links[4], vals[4], addedV = 0, nq;
  double uvall[2], uvback[2];
  vStar *star = NULL;
  Quad  *quad = NULL;

  *activity   = 0;
  stat        = EG_buildStar(qm, &star, vC);
  if (stat   != EGADS_SUCCESS || star == NULL) {
      printf(" In splittingOperation build star %d is NULL %d !!\n", vC, stat);
      return stat;
  }
  nq        = star->nQ;
  stat      = EG_backupQuads(qm, &nq, star->quads, &quad);
  if (stat != EGADS_SUCCESS || quad == NULL) {
      EG_free(quad);
      printf("In EG_splittingOperation :: EG_backupQuads %d !!\n", stat);
      return stat;
  }
  uvall[0] = qm->uvs [ 2 * ( vC - 1)    ];
  uvall[1] = qm->uvs [ 2 * ( vC - 1) + 1];
  EG_centroid (qm, qm->valence[vC -1][2], &qm->valence[vC -1][3], uvall, 1);
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
      for (q = 0; q < star->nQ; q++)
    if (star->quads[q] == -1) break;
      for (i = 0; i < 2; i++) {
      if (i == 1) dist = (star->nV - 1) - 4;
      else        dist = 4;
      id0 = star->idxQ[q + i];
      links[2 * i    ] = star->verts[star->idxV[2 * id0 + 1       ]];
      links[2 * i + 1] = star->verts[star->idxV[2 * id0 + 1 + dist]];
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
  for (j = 0; j < star->nQ; j++) {
      if (star->verts[2 * j + 1] == vL) id0 = j;
      if (star->verts[2 * j + 1] == vR) id1 = j;
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
  poly [0] = star->verts[0];
  poly [1] = star->verts[2 * id0 + 1];
  poly [2] = star->verts[2 * id1 + 1];
  qIdx [0] = id0;
  qIdx [1] = star->idxQ[id1 + star->nQ - 1];
  qIdx [2] = id1;
  qIdx [3] = star->idxQ[id0 + star->nQ - 1];
  verts[0] = poly[1];
  verts[1] = poly[2];
  verts[2] = poly[2];
  verts[3] = poly[1];
  qm->qIdx[4 * (newQ - 1)    ] = poly[1];
  qm->qIdx[4 * (newQ - 1) + 1] = poly[0];
  qm->qIdx[4 * (newQ - 1) + 2] = poly[2];
  qm->qIdx[4 * (newQ - 1) + 3] = poly[3];
  for (i = 0; i < 4; ++i) {
      modQ[i] = star->quads[qIdx[i]];
      qm->qAdj[4 * (newQ - 1) + i] = modQ[i];
      qm->valence[qm->qIdx[4*(newQ - 1)+i]-1][0] = newQ;
      if (modQ[i] == -1) continue;
      stat = EG_adjQtoPair(qm, modQ[i], poly[0], verts[i], adj);
      qm->qAdj[4 * (modQ[i] - 1) + adj[0]] = newQ;
  }
  j = qIdx[2];
  q = star->quads[j++];
  while (q != star->quads[qIdx[0]]) {
      for (i = 0; i < 4; ++i)
    if (qm->qIdx[4 * (q - 1) + i] == poly[0])
        qm->qIdx[4 * (q - 1) + i] = poly[3];
      q = star->quads[star->idxQ[j++]];
      if (q == -1) {
      if (star->quads[qIdx[0]] == -1) break;
      else q = star->quads[star->idxQ[j++]];
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
      for (j = 0; j < qm->valence[poly[i] - 1][2]; j++) {
      stat = setValence(qm, qm->valence[poly[i] - 1][3 + j]);
      if (stat != EGADS_SUCCESS) {
          EG_freeStar(&star);
          EG_free(quad);
          return stat;
      }
      }
  }
  EG_freeStar(&star);
  uvback[0] = qm->uvs[2 * (poly[0] - 1)    ];
  uvback[1] = qm->uvs[2 * (poly[0] - 1) + 1];
  updateVertex(qm, poly[3], uvback);
  /*if (qm->vType[poly[0] - 1] == -1 ) {
      EG_centroid (qm, qm->valence[ poly[0] -1][2],
              &qm->valence[ poly[0] -1][3],
              &qm->uvs[2 * (poly[0] -1)], 1);
      updateVertex(qm, poly[0], &qm->uvs[ 2 * (poly[0] - 1)]);
#ifdef DEBUG
      printf(" SPLIT CENTRE 0 %d \n ", poly[0]);
      printVertex(qm, poly[0]);
#endif
  }*/
  EG_centroid (qm, qm->valence[poly[3] -1][2],
              &qm->valence[poly[3] -1][3],
              &qm->uvs[2 *(poly[3] -1)], 0);
  updateVertex(qm, poly[3], &qm->uvs[2 * (poly[3] - 1)]);
#ifdef DEBUG
  printf(" SPLIT CENTRE 0 %d \n ", poly[3]);
  printVertex(qm, poly[3]);
#endif
  *activity = 1;
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
  updateVertex (qm, poly[0], uvback);
  qm->vType[poly[3] - 1] = -2;
  for (i = 0; i < 4; i++)
    qm->qIdx[4 * (newQ - 1) + i] = -2;
#ifdef DEBUG
  gnuData(qm, NULL);
#endif
  *activity = 0;
  return stat;
}


static int EG_mergeVertices(meshMap *qm, int qC, int centre, int *activity)
{
  int    stat, i, j, q,  adjq, adjPair[2], auxQ, oldQ[8], nq, doublet = 0;
  int    piv[4] = {1, 0, 3, 2};
  double uv[2], uvxyz[10];
  Quad   *quad = NULL;
  vStar  *star = NULL;

#ifdef DEBUG
  printf(" COLLAPSING QUAD %d THRU %d \n ", qC, centre);
  printQuad(qm, qC ) ;
#endif
  nq        = 1;
  stat      = EG_backupQuads(qm, &nq, &qC, &quad);
  if (stat != EGADS_SUCCESS || quad == NULL) {
      printf(" In EG_mergeVertices stat in EG_backpQuads %d !!\n ", stat);
      return stat;
  }
  *activity  = 0;
  i          = EG_quadVertIdx(qm, qC, centre);
  if (qm->vType[centre - 1] != -1) {
      centre = qm->qIdx[4 * (qC - 1) + (i + 2)%4];
      i      = EG_quadVertIdx(qm, qC, centre);
  }
#ifdef DEBUG
      printf(" EG_mergeVertices collapsing through %d with index at %d = %d \n ",
         centre, qC, i);
      printQuad(qm, qC);
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
  for (j  = 0; j < 3; j++) {
      if (j < 2)  {
      uvxyz[j    ] = qm->uvs[2 * (oldQ[0] - 1) + j];
      uvxyz[5 + j] = qm->uvs[2 * (oldQ[2] - 1) + j];
      }
      uvxyz[2 + j]     = qm->xyzs[3 * (oldQ[0] - 1) + j];
      uvxyz[7 + j]     = qm->xyzs[3 * (oldQ[2] - 1) + j];
  }
  // save old quads vertices and adjacents to update map correctly
  q = qm->valence[ centre - 1][3] - 1;
  if (doublet == 0 && qm->vType[oldQ[2] - 1] == -1 &&
      (qm->vType[q] == -1 || qm->vType[q] < 4)) {
      uv[0] = qm->uvs[ 2 * ( oldQ[0] -1 )    ];
      uv[1] = qm->uvs[ 2 * ( oldQ[0] -1 ) + 1];
      stat  = EG_centroid(qm, qm->valence[oldQ[2] -1][2],
              &qm->valence[oldQ[2] -1][3], uv, 1);
       updateVertex(qm, oldQ[2], uv);
  }
  stat      = EG_buildStar(qm, &star, centre);
  if (stat != EGADS_SUCCESS || star == NULL) {
      printf(" EG_mergeVertices :: building star at %d -->%d !!\n ",
             centre, stat);
      return stat;
  }
  for (i = 0; i < 4; i++) {
      q = oldQ[4 + i];
      if (q < 0) continue;
      adjq = oldQ[4 + piv[i]]; // - 1 bias
      stat = EG_adjQtoPair(qm, q, oldQ[i], oldQ[(i + 1)%4], adjPair);
      if (stat != EGADS_SUCCESS || adjPair[1] != qC) {
      EG_freeStar(&star);
      printf("EG_mergeVertices stat %d  adj %d to pair %d %d is %d \n ",
         stat, q ,oldQ[i], oldQ[(i + 1)%4], adjPair[1]);
      return  EGADS_INDEXERR;
      }
      qm->qAdj[4 * (q - 1) + adjPair[0]] = adjq;
  }
  // Eliminate vertex p[0] from all the quads and its valences
  for (i = 0; i < star->nQ; ++i) {
      q = star->quads[i];
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
      EG_freeStar(&star);
      return EGADS_INDEXERR;
      }
  }
  // Set valences for merged stuff
  for (i = 1; i < 4; i++) {
      stat = setValence(qm, oldQ[i]);
      if (stat != EGADS_SUCCESS) {
      printf(" EG_mergeVertices merged setValence %d = %d !!\n", i, stat);
      EG_freeStar(&star);
      return stat;
      }
  }
  // set valences to links
  for (i = 0; i < star->nQ; ++i) {
      j = star->verts[2 * i + 1];
      if (j < 0 || qm->vType[j - 1] == -2) continue; // ghost vertex (-1) or collapsed vertex (-2)
      stat = setValence(qm, j);
      if (stat != EGADS_SUCCESS) {
      printf(" EG_mergeVertices stat in setValence %d = %d !!\n", i, stat);
      EG_freeStar(&star);
      return stat;
      }
  }
  EG_freeStar(&star);
  // delete vertex vC2
  qm->vType[oldQ[0] - 1]    = -2; // -2 = removed
  qm->remQ[++qm->remQ[0]] = qC;
  qm->remV[++qm->remV[0]] = oldQ[0];
  *activity  = 1;
  if (qm->vType  [oldQ[1] - 1]    == -1 &&
      qm->valence[oldQ[1] - 1][2] ==  2) {
      i = qm->valence[oldQ[1] -1][3] - 1;
      j = qm->valence[oldQ[1] -1][4] - 1;
      qm->uvs[ 2 * (oldQ[1] - 1)     ] = 0.5 * (qm->uvs[2 * i    ] +
                                              qm->uvs[2 * j    ]);
      qm->uvs[ 2 * (oldQ[1] - 1) + 1 ] = 0.5 * (qm->uvs[2 * i + 1] +
                                              qm->uvs[2 * j + 1]);
      updateVertex(qm, oldQ[1], &qm->uvs[ 2 * (oldQ[1] - 1)]);
  }
  if (qm->vType  [oldQ[3] - 1]    == -1 &&
      qm->valence[oldQ[3] - 1][2] ==  2) {
      i = qm->valence[oldQ[3] -1][3] - 1;
      j = qm->valence[oldQ[3] -1][4] - 1;
      qm->uvs[ 2 * (oldQ[3] - 1)     ] = 0.5 * (qm->uvs[2 * i    ] +
                                              qm->uvs[2 * j    ]);
      qm->uvs[ 2 * (oldQ[3] - 1) + 1 ] = 0.5 * (qm->uvs[2 * i + 1] +
                                              qm->uvs[2 * j + 1]);
      updateVertex(qm, oldQ[3], &qm->uvs[ 2 * (oldQ[3] - 1)]);
  }
  if (doublet == 1 || EG_makeValidMesh(qm, 3, &oldQ[1], 0) == EGADS_SUCCESS) {
      EG_free(quad);
      return EGADS_SUCCESS;
  }
#ifdef DEBUG
  printf("EG_mergeVertices collapse quad %d thru vertex %d -->%d !!\n",
         qC, oldQ[0], stat);
  gnuData(qm, NULL);
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
#ifdef DEBUG
  gnuData(qm, NULL);
#endif
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
      validSwap(qm, qg.verts[0], qg.verts[3]) != 1) continue;
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
  int i, vC, v5[5], v3[5], qb[4], stat, i3 = 0, val[2], links[3], act = 0;

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
      if (validCollapse(qm, qID, vC) == 1) {
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
      vC = qm->qIdx[4 * (qID - 1) + i3];
      if (qm->vType[vC - 1] == -1 && qm->valence[vC - 1][2] == 3) {
        links[0] = qm->qIdx[4 * (qID - 1) + (i3 + 1)%4];
        links[1] = qm->qIdx[4 * (qID - 1) + (i3 + 3)%4];
        val  [0] = qm->valence[links[0] - 1][1];
        val  [1] = qm->valence[links[1] - 1][1];
        if ((qm->valence[links[0] - 1][1] >= 5 ||
             qm->valence[links[1] - 1][1] >= 5 || i == 1) &&
             validCollapse(qm, qID, vC) == 1) {
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
      validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
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
  int   poly[3], val[3], v, q, id0 = 0, i, stat, dist = 0, validSplit = 0;
  int   id6[2], links[4];
  vStar *star = NULL;

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
    stat      = EG_buildStar(qm, &star, poly[0]);
    if (stat != EGADS_SUCCESS || star == NULL) return stat;
    id6[0]    = -1;
    id6[1]    = -1;
    for (validSplit = q = 0; q < star->nQ; ++q) {
      if (star->quads[q] == - 1) continue;
      id0     = 2 * q + 1;
      poly[1] = star->verts[id0];
      val [1] = qm->valence[poly[1] -1][1];
      for (i = 0; i <= qm->valence[poly[0] - 1][2] - 4; ++i) {
        dist    = 4 + 2*i;
        poly[2] = star->verts[star->idxV[id0 + dist]];
        if (poly[2] < 0) continue; // ghost vertex
        val[2]  = qm->valence[poly[2] -1][1];
        if (val[1] == 3 && val[2] == 3 &&
            qm->vType[poly[0] - 1] == -1) {
          if ((star->nQ == 6 && dist == 6) || star->nQ != 6) {
            validSplit = 1;
            break;
          }
        }
        else if (val[1] <= 3 && val[2] <= 4 && id6[0] == -1) id6[0] = id0;
        else if (val[2] <= 3 && val[1] <= 4 && id6[0] == -1)
          id6[0] = star->idxV[id0 + dist];
        else if (val[1] <= 4 && val[2] <= 4 && id6[0] == -1) id6[1] = id0;
      }
      if (validSplit == 1) break;
    }
    if (validSplit == 0 && val[0] >= 6 && (id6[0] != -1 || id6[1] != -1)) {
      validSplit  = 1;
      /* boundary vertex: Special split since regular = 3 ( not val 4)  */
      if (qm->vType[poly[0] - 1] == 3) {
        for (q = 0; q < star->nQ; q++)
          if (star->quads[q] == -1) break;
        id6[0] = - 1; id6[1] = -1;
        for (i = 0; i < 2; i++) {
          id0 = star->idxQ[q + i];
          if (i == 1) dist = (star->nV - 1) - 4;
          else        dist = 4 + 2 * i;
          links[2 * i    ] = 2 * id0 + 1;
          links[2 * i + 1] = star->idxV[2 * id0 + 1 + dist];
          if (  qm->valence[star->verts[links[2*i     ]] - 1][1] == 4) {
            if (qm->valence[star->verts[links[2*i + 1 ]] - 1][1] == 3 &&
                id6[0] == -1) id6[0] = i;
            else if (qm->valence[star->verts[links[2*i + 1]] -1][1] == 4  &&
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
      stat = EG_splittingOperation(qm, star->verts[0], star->verts[id0],
                                   star->verts[star->idxV[id0 + dist]],
                                   &(*activity));
      EG_freeStar(&star);
      return stat;
    }
    EG_freeStar(&star);
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
  if (validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
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
  int   adj[2], v0;
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
  
  if (swap%3 == 0 || (qm->vType[qg.verts[swap] - 1] >= 0 &&
                      qm->vType[qg.verts[(swap + 3)%6] - 1] >= 0) ||
      qm->vType[qg.verts[i5] - 1] == 4) return EGADS_SUCCESS;
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
    printf("In swapSplit build star for %d -->%d!!\n", qg.verts[i5], stat);
    return stat;
  }
  for (i = 0; i < star->nQ; i++)
    if (star->verts[2 * i + 1] == qg.verts[3]) break;
  if (star->verts[2 * i + 1] != qg.verts[3]) {
    EG_freeStar(&star);
    return EGADS_INDEXERR;
  }
  adj[0] = star->verts[star->idxV[2 * i + 1 + 6]];
  if (star->nQ == 5) {
    adj[1] = star->verts[star->idxV[2 * i + 1 + 4]];
    if (qm->valence[adj[1] - 1][1] < qm->valence[adj[0] - 1][1])
      adj[0] = adj[1];
  }
  EG_freeStar(&star);
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
#ifdef DEBUG
  printf(" EG CLEAN QUAD %d  transfer %d forcing %d\n ", qID, transfer, forcing);
  printQuad (qm, qID ) ;
#endif
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
      if (qID[1] > 0 && qm->qIdx[4 * (qID[1] - 1)] != -1 &&
      qg.q[0] == qID[1]) i = 1;
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


void EG_destroyMeshMap(bodyQuad *bodydata)
{
  int i, j;

  if (bodydata->qm == NULL) return;
  for (i = 0; i < bodydata->nfaces; ++i) {
      if (bodydata->qm[i]) {
      for (j = 0; j < bodydata->qm[i]->sizeV; ++j)
            EG_free(bodydata->qm[i]->valence[j]);
      EG_free(bodydata->qm[i]->valSize);
      EG_free(bodydata->qm[i]->valence);
      EG_free(bodydata->qm[i]->qIdx);
      EG_free(bodydata->qm[i]->qAdj);
      EG_free(bodydata->qm[i]->uvs);
      EG_free(bodydata->qm[i]->xyzs);
      EG_free(bodydata->qm[i]->remQ);
      EG_free(bodydata->qm[i]->remV);
      EG_free(bodydata->qm[i]->vType);
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
  int   i, act, stat = 0, count = 0, j, v[4];
  vStar *star = NULL;

  *activity = 0;
  if (qm->qIdx[4 * (qID - 1)] == -2) return EGADS_SUCCESS;
  for (i = 0; i < 4; i++) v[i]  = qm->qIdx[4 * (qID - 1) + i];
  stat      = EG_cleanQuad(qm, qID, 1, transfer, 0, &act);
  if (stat != EGADS_SUCCESS) return stat;
  *activity += act;
  for (i = 0; i < 4; i++) {
    if (qm->vType[v[i] - 1] == -2) continue;
    stat      = EG_buildStar(qm, &star, v[i]);
    if (star == NULL || stat != EGADS_SUCCESS) {
      printf("EG_cleanNeighborhood : EG_buildStar  = %d \n", stat);
      return stat;
    }
    for (j = 0; j < star->nQ; j++) {
      if (star->quads[j] == -1) continue;
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
        quv[0] = qm->uvs [ 2 * j    ];
        quv[1] = qm->uvs [ 2 * j + 1];
        }
        j = EG_invEvaluateGuess(qm->face, pos, quv, pInv);
        if ( j != EGADS_SUCCESS ) {
            j = EG_invEvaluate(qm->face, pos, quv, pInv);
            if ( j != EGADS_SUCCESS ) continue;
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
  int    i, j, k, q, stat = EGADS_SUCCESS, round, s1, s2, vq[4], si;
  int    ITMAX, it = 0, activity = 0, totActivity = 0, loopact;
  int    iV0, iV, quadPair[2], prevPair[2], totV, totV0, vQ0, vQ, transfer = 0, *skipQuad = NULL, sq;
  double l2Init, linfInit, l2Actual, linfActual, stol, norm1, norm2, err;
  double pos[18*2], quv[2],cross[3], v1[3], v2[3], qN[3], *qArea = NULL;
#ifdef REPORT
  char   buffer[100];
#endif

  // GET RANGE FOR EACH POINT
  meshCount(qm, &iV0, &totV0, &vQ0);
  // PREPROCESS: Reduce the number of quads by forcing collapses everywhere
  EG_distToSurfCentre(qm, &l2Init, &linfInit);
#ifdef REPORT
  printf(" Initial surface %d Quads ::  L2 error approx %1.2e  Linf %1.2e \n", qm->totQ, l2Init, linfInit );
#endif
  if ( l2Init > 1.e-14 && linfInit > 1.e-14 ) {
      it = 0; sq = 0;
      skipQuad   = EG_alloc (qm->totQ * sizeof(int    ));
      qArea      = EG_alloc (qm->totQ * sizeof(double ));
      if ( skipQuad == NULL || qArea == NULL ) return EGADS_MALLOC;
      stol  = 1.1;
      ITMAX = qm->totQ;
      round = 0;
      while ( it < ITMAX ) {
      EG_distToSurfCentre(qm, &l2Actual, &linfActual);
      k = -1;
      for (si = j = i = 0; i < qm->totQ; i++) {
          if (qm->qIdx[ 4 * i ] == -2 ) continue;
          pos[0]  = pos[1] = pos[2] = 0.0;
          quv[0]  = quv[1] = 0.0;
          for ( q = 0; q < 4; q++ ) {
          vq[q]   = qm->qIdx [ 4 * i + q] - 1;
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
          stat   = EG_invEvaluateGuess(qm->face, pos, quv, &pos[18] );
          err    = (pos[0] - pos[18]) * (pos[0] - pos[18]) +
               (pos[1] - pos[19]) * (pos[1] - pos[19]) +
               (pos[2] - pos[20]) * (pos[2] - pos[20]);
          stat   = EG_evaluate(qm->face, quv, pos );
          norm1  = pos[3] * pos[3] + pos[4] * pos[4] + pos[5] * pos[5];
          norm2  = pos[6] * pos[6] + pos[7] * pos[7] + pos[8] * pos[8];
          if (norm1 < EPS11 || norm2 < EPS11 ||
          inList (sq, skipQuad, i + 1 ) != -1 ) continue;
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
          CROSS(v1, v2, cross);
          qArea[i] = 0.5 * sqrt ( DOT (cross, cross) ) ;
          s1       = 1;
          if (DOT(qN, cross) < 0) s1 = -1;
          v1[0]    = qm->xyzs[3 * vq[3]    ] - qm->xyzs[3 * vq[0]    ];
          v1[1]    = qm->xyzs[3 * vq[3] + 1] - qm->xyzs[3 * vq[0] + 1];
          v1[2]    = qm->xyzs[3 * vq[3] + 2] - qm->xyzs[3 * vq[0] + 2];
          CROSS(v2, v1, cross);
          qArea[i] += 0.5 * sqrt ( DOT (cross, cross) ) ;
          s2        = 1;
          if (DOT(qN, cross) < 0) s2 = -1;
          if ( s1 * s2 == -1 && si == 0 ) {
#ifdef DEBUG
          printf(" Quad %d self-intersects: try destroy\n ", i + 1 );
          gnuData(qm, NULL);
          printQuad(qm, i + 1 ) ;
#endif
          qArea[i] = 0.0; //invalid quad: self-intersects.
          k  = i;
          si = 1;
          }
          if (err >= linfActual ) continue;
          if ( si ==  0 &&
          ((k == -1 && EG_nValenceCount(qm, i + 1, 3) > 0) ||
           (k != -1 && fabs(qArea[i]) < fabs(qArea[k])))) k = i;
      }
#ifdef DEBUG
      printf("IT %d (ROUND %d) L2 err %1.2e  vs %1.2e Linf %1.2e VS %1.2e -->%.2f %%\n",
         it, round, l2Init, l2Actual, linfInit, linfActual, l2Actual / linfActual * 100.0 );
#endif
      if (l2Actual / l2Init > stol ) {
#ifdef DEBUG
          printf(" LEAVE %.12f VS %.12f !\n ", l2Actual / l2Init, stol);
#endif
          break;
      }
      if (k == -1 ) {
          for (k = j = i = 0; i < qm->totQ; i++) {
          if (qm->qIdx [ 4 * i ] == -2 ||
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
          EG_free (skipQuad);
          EG_free (qArea   );
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
      stat      = resizeQm(qm );
      if (stat != EGADS_SUCCESS) {
      printf(" EG_meshRegularization resizing after preprocess %d !!\n", stat);
      return stat;
      }
  }
  EG_makeValidMesh(qm, 0, NULL, 2);
#ifdef REPORT
      EG_distToSurfCentre(qm, &l2Actual, &linfActual);
      printf(" Preprocessing step :: total collapses %d\n "
      "surface errors L2 %1.2e Linf %1.2e\n",it, l2Actual, l2Init);
#endif
#ifdef REPORT
      snprintf(buffer,100, "gnuPreprocess_%d.txt", qm->fID);
      gnuData(qm, buffer);
      snprintf(buffer,100, "wvsPreprocess_%d.txt", qm->fID);
      wvsData(qm, buffer);
      snprintf(buffer,100, "preprocessMeshStats_%d.txt", qm->fID);
      meshStats(qm, 10);
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
          quadPair[0]   = q + 1;
          quadPair[1]   = -1;
          stat          = EG_transferValences(qm, quadPair, k,
                                              &transfer, &activity);
          if (stat     != EGADS_SUCCESS) return stat;
          if (activity == 0 || quadPair[0] < 0 ||
              qm->qIdx[4*(quadPair[0]-1)] == -2) continue;
          totActivity += activity;
          for (j = 0; j < 20; j++) {
            activity = 0;
            for (i  = 0; i < 2; i++) {
              if (quadPair[i] < 0 ||
                  qm->qIdx[4*(quadPair[i]-1)] == -2) continue;
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
                qm->qIdx[4*(quadPair[0]-1)] == -2) break;
          }
          if (iV <= 2) break;
          prevPair[0] = quadPair[0];
          prevPair[1] = quadPair[1];
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
#ifdef REPORT
  meshCount(qm, &iV, &totV, &vQ);
  printf("\n\n Original mesh had %d QUADS, %d / %d irregular vertices (%.2f %%) \n",
         vQ0, iV0, totV0, (double) iV0 * 100.0 / (double)totV0);
  printf("\n\n*******************************************************************\n");
  printf(" Final mesh has %d QUADS, %d / %d irregular vertices (%.2f %%) ============\n",
         vQ, iV, totV, (double) iV * 100.0 / (double) totV);
  printf(" Total invalid steps  %d\n ", qm->invsteps);
  printf("*******************************************************************\n");
#endif
  stat      = resizeQm(qm);
  if (stat != EGADS_SUCCESS) {
    printf("EG_meshRegularization final resize %d !!\n", stat);
    return stat;
  }
#ifdef REPORT
  snprintf(buffer,100, "wvsRegular_%d.txt", qm->fID);
  wvsData(qm, buffer);
  snprintf(buffer,100,"gnuRegular_%d.txt", qm->fID);
  gnuData(qm, buffer);
#endif
  return EG_makeValidMesh(qm, 0, NULL, 1);
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
    if (stat != EGADS_SUCCESS) goto cleanup;
    for (f = 0; f < bodydata[iBody].nfaces; ++f) {
      if (FACECHOICE >= 0) f = FACECHOICE -1;
#ifdef REPORT
          meshStats(bodydata[iBody].qm[f], 0);
          sprintf(buffer,"gnuInit_%i.txt",f+1);
          gnuData(bodydata->qm[f] , buffer);
          sprintf(buffer,"wvsInit_%i.txt",f+1);
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
