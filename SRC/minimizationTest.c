#include "egads.h"
#include "regQuads.h"


#define DEBUG


#ifdef DEBUG
#define REPORT
#endif

#define REPORT

#ifdef REPORT
  #include <time.h>
#endif


/* BASE-LEVEL FUNCTION */

// Solves the square linear system A x = b using simple LU decomposition
// Returns CAPS_BADVALUE for a singular matrix
int solveLU(int n, double A[], double b[], double x[] )
{
     int i,j,k;
     double y;

     // Compute the LU decomposition in place
     for(k = 0; k < n-1; k++) {
         if (A[k*n+k] == 0) return EGADS_DEGEN;
         for(j = k+1; j < n; j++) {
             y = A[j*n+k]/A[k*n+k];
             for(i = k; i < n; i++) {
                 A[j*n+i] = A[j*n+i]-y*A[k*n+i];
             }
             A[j*n+k] = y;
         }
     }

     // Forward solve
     for(i = 0; i < n; i++) {
         y=0.0;
         for(j = 0 ;j < i;j++) {
             y += A[i*n+j]*x[j];
         }
         x[i]=(b[i]-y);
     }

     // Back substitution
     for(i = n-1; i >=0; i--) {
         y = 0.0;
         for(j = i+1; j < n; j++) {
             y += A[i*n+j]*x[j];
         }
         x[i] = (x[i]-y)/A[i*n+i];
     }

     return EGADS_SUCCESS;
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
  EG_free((*star)->area);
  EG_free((*star)->angle);
  EG_free((*star)->ratio);
  EG_free (*star);
  *star = NULL;
}


static int EG_buildStar(meshMap *qm, vStar **star, int vID)
{
  int i = 0, id0 = -1, q = 0, auxV, auxQ, v = 0, quadID, prevQuad,
      it = 0, it2 = 0, adj[2], *vertex = NULL, *quads = NULL;
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
  (*star)->type  = -1;
  (*star)->verts = (int    *) EG_alloc(    v * sizeof(int));
  (*star)->quads = (int    *) EG_alloc(    q * sizeof(int));
  (*star)->idxV  = (int    *) EG_alloc(2 * v * sizeof(int));
  (*star)->idxQ  = (int    *) EG_alloc(2 * q * sizeof(int));
  (*star)->area  = (int    *) EG_alloc(    q * sizeof(int));
  (*star)->ratio = (double *) EG_alloc(    q * sizeof(double));
  (*star)->angle = (double *) EG_alloc(    q * sizeof(double));
  if ((*star)->verts == NULL || (*star)->quads == NULL ||
      (*star)->idxV  == NULL || (*star)->idxQ  == NULL ||
      (*star)->area  == NULL || (*star)->ratio == NULL ||
      (*star)->angle == NULL ) {
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




#define EPS 1.e-10
// Find the point P = S(u,v) midpoint of uv0 and uv1 that minimizes the distance between
// min 1/2 ( l_0 ^ 2 + l_1 ^ 2) + lambda * (l_1 - l_0) = L(u, v, lambda)
// grad (L) = (L_1, L_2, L_3) = (0, 0, 0) --> solution using Newton x_n+1 = x_n + delta_n (3x3 system)
static void EG_minArc(ego face, double fact, double *uv0, double *uv1, double *uv)
{
  int    i, it, nT = 100;
  double tol, l[2],  p0[18], p1[18], pn[18], J[3][3], ATJ[3][3], A[3 * 3],
         r0[3], r1[3],b, delta[3],  dlu0, dlu1, dlv0, dlv1, F[3],
         ddl0, ddl1, L[3], detJ, x0, x1, x2,  e1, e2, uvn[3], ctt;

  //Initial guess uv = 0.5 (uv0 + uv1)
  uvn[0] =  0.5 * (uv0[0] + uv1[0]);
  uvn[1] =  0.5 * (uv0[1] + uv1[1]);
  uvn[2] = -0.0;  //lambda

  i      = EG_evaluate(face, uv0, p0);
  i     += EG_evaluate(face, uv1, p1);
  i     += EG_evaluate(face, uvn, pn);
  if (i != EGADS_SUCCESS) {
   printf(" EG_minArc :: EG_evaluate %d !!\n ", i);
   return;
  }

#ifdef DEBUG
  r0[0] = pn[0] - p0[0];
  r0[1] = pn[1] - p0[1];
  r0[2] = pn[2] - p0[2];

  r1[0] = pn[0] - p1[0];
  r1[1] = pn[1] - p1[1];
  r1[2] = pn[2] - p1[2];

  l [0] = DOT(r0, r0);
  l [1] = DOT(r1, r1);
  printf(" INITIAL GUESS %lf  %lf  %lf  %lf  %lf\n",
  pn[0], pn[1], pn[2], uvn[0], uvn[1]);
  b    = sqrt (l[0]) + sqrt(l[1]);
  printf(" --------------------------------------------------- \n");
  printf(" TOTAL ARC %lf final chords ||PM - P0||^2 = %lf ||PM - P1||^2 = %lf\n ratio %lf %lf (EXPECTED %lf %lf)\n",
  b, l[0], l[1], l[0] / (b * b), l[1]/(b * b), fact * fact, (1.0 - fact) * (1.0 - fact) );
  l[0] = sqrt(l[0]);
  l[1] = sqrt(l[1]);
  b    = l[0] + l[1];
  printf(" As distances ||PM - P0|| = %lf ||PM - P1|| = %lf\n ratio %lf %lf (EXPECTED %lf %lf)\n",
  l[0], l[1], l[0] / b, l[1]/b, fact, (1.0 - fact));
  printf(" --------------------------------------------------- \n");
  printf("%lf %lf %lf P0\n", p0[0], p0[1], p0[2]);
  printf("%lf %lf %lf P1\n", p1[0], p1[1], p1[2]);
  printf("%lf %lf %lf Pn\n", pn[0], pn[1], pn[2]);
#endif


  tol  = 1.e-10;
  ctt  = (1 - fact) * (1.0 - fact) / (fact * fact);
  for (it = 0; it < nT; it++) {

#ifdef DEBUG
  printf("NEW POINT \n %lf %lf %lf %lf %lf\n",
  pn[0], pn[1], pn[2], uvn[0], uvn[1]);
  printf(" du  %lf %lf %lf\n dv  %lf %lf %lf\n"
  " duu %lf %lf %lf\n duv %lf %lf %lf\n"
  " dvv %lf %lf %lf\n",   pn [3], pn [4], pn[5],
  pn [6], pn [7], pn [8], pn [9], pn[10], pn[11],
  pn[12], pn[13], pn[14], pn[15], pn[16], pn[17]);
#endif

  r0[0] = pn[0] - p0[0];
  r0[1] = pn[1] - p0[1];
  r0[2] = pn[2] - p0[2];

  r1[0] = pn[0] - p1[0];
  r1[1] = pn[1] - p1[1];
  r1[2] = pn[2] - p1[2];

  l [0] = DOT(r0, r0);
  l [1] = DOT(r1, r1);

#ifdef DEBUG
        printf(" DISTANCE FUNCTIONS \n r0 %lf %lf %lf \n r1 %lf %lf %lf\n",
        r0[0], r0[1], r0[2], r1[0], r1[1], r1[2]);
#endif
  dlu0    = 2.0 * (r0[0] * pn[3] + r0[1] * pn[4] + r0[2] * pn[5]);
  dlu1    = 2.0 * (r1[0] * pn[3] + r1[1] * pn[4] + r1[2] * pn[5]);
  dlv0    = 2.0 * (r0[0] * pn[6] + r0[1] * pn[7] + r0[2] * pn[8]);
  dlv1    = 2.0 * (r1[0] * pn[6] + r1[1] * pn[7] + r1[2] * pn[8]);

  L[0]    = dlu0 * (l[0] - ctt * uvn[2]) + dlu1 * (l[1] + uvn[2]);
  L[1]    = dlv0 * (l[0] - ctt * uvn[2]) + dlv1 * (l[1] + uvn[2]);
  L[2]    = l[1] - ctt * l[0];
#ifdef DEBUG
  printf(" DLU %lf %lf  DLV %lf %lf \n", dlu0, dlu1, dlv0, dlv1);
  printf(" LU %lf %lf %lf ", L[0], L[1], L[2]);
#endif
//duu
  ddl0    = 2.0 * (pn[9 ] * r0[0] + pn[3] * pn[3] +
                   pn[10] * r0[1] + pn[4] * pn[4] +
                   pn[11] * r0[2] + pn[5] * pn[5]);
  ddl1    = 2.0 * (pn[9 ] * r1[0] + pn[3] * pn[3] +
                   pn[10] * r1[1] + pn[4] * pn[4] +
                   pn[11] * r1[2] + pn[5] * pn[5]);
  b       = ddl0 * (l[0] - ctt * uvn[2]) + ddl1* (l[1] + uvn[2]);
  J[0][0] = dlu0 * dlu0 + dlu1 * dlu1 + b;

  //duv
  ddl0    = 2.0 * (pn[12] * r0[0] + pn[3] * pn[6] +
                   pn[13] * r0[1] + pn[4] * pn[7] +
                   pn[14] * r0[2] + pn[5] * pn[8]);

  ddl1    = 2.0 * (pn[12] * r1[0] + pn[3] * pn[6] +
                   pn[13] * r1[1] + pn[4] * pn[7] +
                   pn[14] * r1[2] + pn[5] * pn[8]);
  b       =  ddl0 * (l[0] - ctt * uvn[2]) + ddl1 * (l[1] + uvn[2]);
  J[0][1] =  dlu0 * dlv0 + dlu1 * dlv1 + b;

  J[0][2] = dlu1 - ctt * dlu0;

  //dvv
  ddl0    = 2.0 * (pn[15] * r0[0] + pn[6] * pn[6] +
                   pn[16] * r0[1] + pn[7] * pn[7] +
                   pn[17] * r0[2] + pn[8] * pn[8]);

  ddl1    = 2.0 * (pn[15] * r1[0] + pn[6] * pn[6] +
                   pn[16] * r1[1] + pn[7] * pn[7] +
                   pn[17] * r1[2] + pn[8] * pn[8]);
  b       = ddl0 * (l[0] - ctt * uvn[2]) + ddl1 * (l[1] + uvn[2]);
  J[1][0] = J[0][1];
  J[1][1] = dlv0 * dlv0 + dlv1 * dlv1 + b;
  J[1][2] = dlv1 - ctt * dlv0;

  J[2][0] = J[0][2];
  J[2][1] = J[1][2];
  J[2][2] = 0.0;

  // Solve Linear System: J * delta = - f_n
  // For now: Invert Jacobian directly --> delta = J^-1 f_n

  detJ = J[0][0] * J[1][1] * J[2][2] + J[0][1] * J[1][2] * J[2][0] +
         J[0][2] * J[2][1] * J[1][0] - J[0][2] * J[1][1] * J[2][0] -
         J[1][2] * J[2][1] * J[0][0] - J[2][2] * J[1][0] * J[0][1];
#ifdef DEBUG
  printf(" JACOBIAN MATRIX 3 x 3  DET %1.2e\n", detJ);
  printf(" %lf %lf %lf \n",J[0][0], J[0][1], J[0][2]);
  printf(" %lf %lf %lf \n",J[1][0], J[1][1], J[1][2]);
  printf(" %lf %lf %lf \n",J[2][0], J[2][1], J[2][2]);
  printf(" ----------------------\n");
#endif
  if (fabs(detJ) < 1.e-08) {
      printf(" DETERMINANT SIZE %lf  \n", detJ);
      break;
  }
  ATJ[0][0] =    J[1][1] * J[2][2] - J[2][1] * J[1][2];
  ATJ[0][1] = - (J[0][1] * J[2][2] - J[2][1] * J[0][2]);
  ATJ[0][2] =    J[0][1] * J[1][2] - J[1][1] * J[0][2];
  ATJ[1][0] = - (J[1][0] * J[2][2] - J[2][0] * J[1][2]);
  ATJ[1][1] =    J[0][0] * J[2][2] - J[2][0] * J[0][2];
  ATJ[1][2] = - (J[0][0] * J[1][2] - J[1][0] * J[0][2]);
  ATJ[2][0] =    J[1][0] * J[2][1] - J[2][0] * J[1][1];
  ATJ[2][1] = - (J[0][0] * J[2][1] - J[2][0] * J[0][1]);
  ATJ[2][2] =    J[0][0] * J[1][1] - J[1][0] * J[0][1];

#ifdef DEBUG
  printf(" INVERSE MATRIX 3 x 3  WITHOUT 1/DET\n");
  printf(" %lf %lf %lf \n",ATJ[0][0], ATJ[0][1], ATJ[0][2]);
  printf(" %lf %lf %lf \n",ATJ[1][0], ATJ[1][1], ATJ[1][2]);
  printf(" %lf %lf %lf \n",ATJ[2][0], ATJ[2][1], ATJ[2][2]);
  printf(" ----------------------\n");
#endif
  printf(" LO %lf  L1 %lf  L2 %lf\n", L[0], L[1], L[2]);
  A[0] = J[0][0]; A[1] = J[0][1]; A[2] = J[0][2];
  A[3] = J[1][0]; A[4] = J[1][1]; A[5] = J[1][2];
  A[6] = J[2][0]; A[7] = J[2][1]; A[8] = J[2][2];
  F[0] = -L[0]; F[1] = -L[1]; F[2] = -L[2];
  solveLU(3, A, F, delta);
  printf(" DELTA %lf %lf %lf\n", delta[0], delta[1], delta[2]);
  detJ      = 1.0 / detJ;
  delta[0]  = -detJ * (ATJ[0][0] * L[0] + ATJ[0][1] * L[1] + ATJ[0][2] * L[2]);
  delta[1]  = -detJ * (ATJ[1][0] * L[0] + ATJ[1][1] * L[1] + ATJ[1][2] * L[2]);
  delta[2]  = -detJ * (ATJ[2][0] * L[0] + ATJ[2][1] * L[1] + ATJ[2][2] * L[2]);


  // Using LU
  printf(" DELTA %lf %lf %lf\n", delta[0], delta[1], delta[2]);
  uvn[0]   += delta[0];
  uvn[1]   += delta[1];
  uvn[2]   += delta[2];
  x2        = sqrt(DOT(delta, delta));
  if      (it == 0) x0 = x2;
  else if (it == 1) x1 = x2;
  else {
      e1 = fabs(x1 / x0);
      e2 = fabs(x2 / x1);
      printf(" ERR IT %d e1 %lf  e2  %lf alpha %lf x0 %lf x1 %lf x2 %lf \n",
      it, e1, e2, log(e2) / log(e1), x0, x1, x2);
      x0 = x1;
      x1 = x2;
  }
  i = EG_evaluate(face, uvn, pn);
  if (i != EGADS_SUCCESS || x2 < tol ) {
      printf("EG_evaluate %d  DELTA SIZE %1.2e < %1.2e \n",
      EG_evaluate(face, uvn, pn), x2, tol);
      break;
  }

#ifdef DEBUG
   printf(" it %d --xn %lf %lf %lf DELTA %lf %lf %lf   SIZE  %1.2e < %1.2e\n",
   it, uvn[0], uvn[1], uvn[2], delta[0], delta[1], delta[2], x2, tol );
#endif
  }
  uv[0] = uvn[0];
  uv[1] = uvn[1];
#ifdef DEBUG
  printf(" LEAVING WITH FINAL POINT %lf  %lf  %lf  %lf  %lf\n",
  pn[0], pn[1], pn[2], uv[0], uv[1]);
  b    = sqrt (l[0]) + sqrt(l[1]);
  printf(" --------------------------------------------------- \n");
  printf(" Leaving frac solver TOTAL ARC %lf with final chords ||PM - P0|| = %lf ||PM - P1|| = %lf\n ratio %lf %lf (EXPECTED %lf %lf)\n",
  b, l[0], l[1], l[0] / (b * b), l[1]/(b * b), fact * fact, (1.0 - fact) * (1.0 - fact) );
  l[0] = sqrt(l[0]);
  l[1] = sqrt(l[1]);
  b    = l[0] + l[1];
  printf(" Leaving frac solver with final chords ||PM - P0|| = %lf ||PM - P1|| = %lf\n ratio %lf %lf (EXPECTED %lf %lf)\n",
  l[0], l[1], l[0] / b, l[1]/b, fact, (1.0 - fact));
  printf(" --------------------------------------------------- \n");
  printf("%lf %lf %lf P0\n", p0[0], p0[1], p0[2]);
  printf("%lf %lf %lf P1\n", p1[0], p1[1], p1[2]);
  printf("%lf %lf %lf Pn\n", pn[0], pn[1], pn[2]);
#endif
  return;
}


int EG_createMeshMap(bodyQuad *bodydata)
{
  int          f, stat = 0, stat2 = 0, j, q, k, kk, kOK, len;
  int          ntri, nquad;
  const int    *tris, *tric, *ptype, *pindex;
  const double *xyzs, *uvs;
  int          qV[6]    = { 0, 1, 2, 5, 0, 1};
  int          qLoop[5] = { 0, 1, 2, 3, 0   };
  ego          tess;
  egTessel     *btess;

  tess  = bodydata->tess;
  btess = (egTessel *) tess->blind;

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
      bodydata->qm[f]->oriV   = len;
      bodydata->qm[f]->oriQ   = nquad;
      bodydata->qm[f]->sizeV  = 2 * len;
      bodydata->qm[f]->sizeQ  = 2 * nquad;
      bodydata->qm[f]->totV   = len;
      bodydata->qm[f]->totQ   = nquad;
      bodydata->qm[f]->vInv   = NULL;
      bodydata->qm[f]->regBd0 = 1;
      if( btess->tess2d[f].tfi == 1) {
          bodydata->qm[f]->xyzs    = NULL;
          bodydata->qm[f]->uvs     = NULL;
          bodydata->qm[f]->vType   = NULL;
          bodydata->qm[f]->qIdx    = NULL;
          bodydata->qm[f]->qAdj    = NULL;
          bodydata->qm[f]->remQ    = NULL;
          bodydata->qm[f]->remV    = NULL;
          bodydata->qm[f]->valence = NULL;
          bodydata->qm[f]->star    = NULL;
          bodydata->qm[f]->bdAng   = NULL;
          bodydata->qm[f]->fID     = 0;
          continue;
      }
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
          bodydata->qm[f]->fID = 0;
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
              bodydata->qm[f]->fID = 0;
              continue;
          }
          bodydata->qm[f]->valence[j][2] = 0;
      }
      bodydata->qm[f]->face = bodydata->faces[f];
      stat   = EG_getRange(bodydata->qm[f]->face, bodydata->qm[f]->range, &j);
      for (k = j = 0; j < len; j++) {
          bodydata->qm[f]->vType  [j]         = ptype[j];
          if (ptype[j] >= 0 ) k++;
          bodydata->qm[f]->uvs    [2 * j    ] = uvs[2 * j    ];
          bodydata->qm[f]->uvs    [2 * j + 1] = uvs[2 * j + 1];
          if (uvs[2 * j    ] < bodydata->qm[f]->range[0] - 1.e-6 ||
              uvs[2 * j    ] > bodydata->qm[f]->range[1] + 1.e-6 ||
              uvs[2 * j + 1] < bodydata->qm[f]->range[2] - 1.e-6 ||
              uvs[2 * j + 1] > bodydata->qm[f]->range[3] + 1.e-6) {
              printf(" EG_createMeshMap :: vertex %d = %f  %f  is out of range !! \n ",
                     j + 1, uvs[2 * j], uvs[2 * j + 1]);
              printf(" Range FACE %d --- > %f  %f  %f  %f\n ", bodydata->qm[f]->fID,
                     bodydata->qm[f]->range[0],bodydata->qm[f]->range[1],
                     bodydata->qm[f]->range[2],bodydata->qm[f]->range[3]);
              stat2 = EGADS_INDEXERR;
              bodydata->qm[f]->fID = 0;
              break;
          }
          bodydata->qm[f]->xyzs[3 * j    ] = xyzs[3 * j    ];
          bodydata->qm[f]->xyzs[3 * j + 1] = xyzs[3 * j + 1];
          bodydata->qm[f]->xyzs[3 * j + 2] = xyzs[3 * j + 2];
      }
      bodydata->qm[f]->bdAng = (double *) EG_alloc(k * sizeof(double));
      if (bodydata->qm[f]->bdAng == NULL ) {
          bodydata->qm[f]->fID = 0;
          stat2 = EGADS_MALLOC;
          continue;
      }
      if (stat2 != EGADS_SUCCESS) continue;
      for (j = 0; j < nquad; j++) {
          for (k = 0; k < 4; ++k)
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
              bodydata->qm[f]->fID = 0;
              break;
          }
      }
  }
  return stat2;
}


#ifdef STANDALONE
int main (int argc, char *argv[])
{
  clock_t      start_t, end_t, total_t;
  int          stat = 0,  f , i, j, iBody, oclass, mtype, nbody;
  int          atype, alen, *senses, min, FACECHOICE = -1;
  const int    *ints;
  float        arg;
  double       box[6], size, params[3], time, frac;
  const double *reals;
  const char   *OCCrev, *string;
  ego          context, tess, model, geom, *bodies, *dum;
  bodyQuad     *bodydata;

  start_t = clock();
  if (argc != 3 && argc != 2 && argc != 5 && argc != 6 && argc != 7) {
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


  params[0] = 0.1;
  params[1] = 0.09;
  params[2] = 19.0;
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
  if (argc >= 6) sscanf(argv[5], "%d", &FACECHOICE);
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
          printf(" Body %d: Type = WireBody\n",  iBody+1);
      } else if (mtype == FACEBODY) {
          printf(" Body %d: Type = FaceBody\n",  iBody+1);
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
      stat = EG_createMeshMap(&bodydata[iBody] );
      if (stat != EGADS_SUCCESS) {
          printf(" EG_createMeshMap %d !!!\n", stat);
          goto cleanup;
      }
      if (FACECHOICE > bodydata[iBody].nfaces) FACECHOICE = -1;
      for (f = 0; f < bodydata[iBody].nfaces; ++f) {
          if (FACECHOICE >= 0) f = FACECHOICE - 1;

          int iA, iB;
          double uv0[2], uv1[2], uv[2];
          for (i=0; i<bodydata[iBody].qm[f]->totQ; i++) {
              if (bodydata[iBody].qm[f]->vType[i] != -1 ||
                  bodydata[iBody].qm[f]->star[i] == NULL )continue;
              iA     = bodydata[iBody].qm[f]->star[i]->verts[1] - 1;
              iB     = bodydata[iBody].qm[f]->star[i]->verts[3] - 1;
              uv0[0] = bodydata[iBody].qm[f]->uvs[2 * iA    ];
              uv0[1] = bodydata[iBody].qm[f]->uvs[2 * iA + 1];
              uv1[0] = bodydata[iBody].qm[f]->uvs[2 * iB    ];
              uv1[1] = bodydata[iBody].qm[f]->uvs[2 * iB + 1];
              EG_minArc(bodydata[iBody].qm[f]->face, 0.25, uv0, uv1, uv);
          }
          if (stat != EGADS_SUCCESS) {
              fprintf(stderr, " EG_meshRegularization face %d / %d = %d DEACTIVATE !!!\n ",
                      f + 1, bodydata[iBody].nfaces,  stat);
              bodydata[iBody].qm[f]->fID = 0;
              continue;
          }
          if (FACECHOICE >= 0) break;
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
#ifdef REPORT
  end_t   = clock();
  total_t = end_t - start_t;
  time    = (double) total_t / CLOCKS_PER_SEC;
  min     = floor(time) / 60;
  frac    = time -  min * 60;
  fprintf(stderr, "Total time taken by CPU: %d minutes and %f seconds\n",
          min, frac);
  printf("Total time taken by CPU: %d minutes and %f seconds\n",
         min, frac);
#endif
  EG_deleteObject(model);
  EG_close(context);
  return 0;
}
#endif
