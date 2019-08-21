#include "egads.h"
#include <math.h>
#include <string.h>
#include "Surreal/SurrealS.h"
#define CROSS(a,b,c)      c[0] = (a[1]*b[2]) - (a[2]*b[1]);\
                          c[1] = (a[2]*b[0]) - (a[0]*b[2]);\
                          c[2] = (a[0]*b[1]) - (a[1]*b[0])
#define DOT(a,b)          (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define DEBUG


#ifdef DEBUG
#define REPORT
#endif

#define REPORT

#ifdef REPORT
  #include <time.h>
#endif

typedef struct {
  ego      tess,   *faces;
  int      nedges, nfaces;
} bodyQuad;

/* BASE-LEVEL FUNCTION */

// Solves the square linear system A x = b using simple LU decomposition
// Returns CAPS_BADVALUE for a singular matrix
int solveLU(int n, double A[], double b[], double x[] )
{
     int i,j,k;
     double y;


     printf(" MATRIX IN \n");
     printf(" %lf %lf %lf %lf\n", A[0], A[1], A[2], A[3]);
     printf(" %lf %lf %lf %lf\n", A[4], A[5], A[6], A[7]);
     printf(" %lf %lf %lf %lf\n", A[8], A[9], A[10], A[11]);
     printf(" %lf %lf %lf %lf\n", A[12], A[13], A[14], A[15]);
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

template <class T>
static void EG_minL (double w1, double w2, double w3, T lambda,
                     T mu, T *pO, double *pA, double *pB, double *pC, T *L)
{
    T da, db, dc, crA[3], crB[3], crC[3], duA[3], duB[3], duC[3],dvA[3], dvB[3],
      dvC[3], a, b, c, v[3], BA[3], CB[3], AC[3], OA[3], OB[3], OC[3];

    da = 1.0 + lambda;
    dc = 1.0 + mu;
    db = 1.0 - (w1 * w1) / (w2 * w2) * lambda - (w3 * w3) / (w2 * w2) * mu;
    OA[0] = pA[0] - pO[0]; OA[1] = pA[1] - pO[1]; OA[2] = pA[2] - pO[2];
    OB[0] = pB[0] - pO[0]; OB[1] = pB[1] - pO[1]; OB[2] = pB[2] - pO[2];
    OC[0] = pC[0] - pO[0]; OC[1] = pC[1] - pO[1]; OC[2] = pC[2] - pO[2];
    CROSS (OA, OB, crC);
    c     = DOT(crC,crC);
    CROSS (OB, OC, crA);
    a     = DOT(crA,crA);
    CROSS (OC, OA, crB);
    b     = DOT(crB,crB);
    BA[0] = pA[0] - pB[0]; BA[1] = pA[1] - pB[1]; BA[2] = pA[2] - pB[2];
    CB[0] = pB[0] - pC[0]; CB[1] = pB[1] - pC[1]; CB[2] = pB[2] - pC[2];
    AC[0] = pC[0] - pA[0]; AC[1] = pC[1] - pA[1]; AC[2] = pC[2] - pA[2];
    v [0] = pO[3]; v[1] = pO[4]; v[2] = pO[5];
    CROSS (v, CB, duA);
    CROSS (v, AC, duB);
    CROSS (v, BA, duC);
    v [0] = pO[6]; v[1] = pO[7]; v[2] = pO[8];
    CROSS (v, CB, dvA);
    CROSS (v, AC, dvB);
    CROSS (v, BA, dvC);

    L[0] = da * DOT(crA,duA) + db * DOT(crB,duB) + dc * DOT(crC,duC);
    L[1] = da * DOT(crA,dvA) + db * DOT(crB,dvB) + dc * DOT(crC,dvC);
    L[2] = 0.5 * a - 0.5 * (w1 * w1) / (w2 * w2) * b;
    L[3] = 0.5 * c - 0.5 * (w3 * w3) / (w2 * w2) * b;
    return;
}


static void EG_baricentricSplit(ego face, const double w1, const double w2,const  double w3,
                                const double *uvA, const double *uvB,const  double *uvC, double *uvOUT)
{
  int    i, it, nT = 100, stat;
  double tol, pA[18], pB[18], pC[18], OA[3], OB[3], OC[3], pIT[18], J[16], JI[16],
         crA[3], crB[3], crC[3], grad[3], a, b, c, d, delta[4], L[4], x0, ka, kc,
         x1, x2, e1 = 0.0, e2 = 0.0, uvIT[8], duA[3], duB[3], duC[3], CB[3],
         dvA[3], dvB[3], dvC[3], AC[3], BA[3], da, db, dc, det;
  uvIT[0]   = w1 * uvA[0] + w2 * uvB[0];// + w3 * uvC[0];
  uvIT[1]   = w1 * uvA[1] + w2 * uvB[1];// + w3 * uvC[1];
  uvIT[0]   = 0.25 * uvA[0] + 0.5 * uvB[0] + 0.25 * uvC[0];
  uvIT[1]   = 0.25 * uvA[1] + 0.5 * uvB[1] + 0.25 * uvC[1];
  //Initial guess uv = 0.3 (uvA + uvB + uvC)
  uvIT[2]   = uvIT[3] = 1.0;
  ka        = (w1 * w1) / (w2 * w2);
  kc        = (w3 * w3) / (w2 * w2);
  stat      = EG_evaluate(face, uvA , pA);
  stat     += EG_evaluate(face, uvB , pB);
  stat     += EG_evaluate(face, uvC , pC);
  stat     += EG_evaluate(face, uvIT, pIT);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_baricentricSplit EG_evaluate triangle vertices stat %d !!\n", stat);
    return;
  }
#ifdef DEBUG
  OA[0] = pA[0] - pIT[0];
  OA[1] = pA[1] - pIT[1];
  OA[2] = pA[2] - pIT[2];

  OB[0] = pB[0] - pIT[0];
  OB[1] = pB[1] - pIT[1];
  OB[2] = pB[2] - pIT[2];

  OC[0] = pC[0] - pIT[0];
  OC[1] = pC[1] - pIT[1];
  OC[2] = pC[2] - pIT[2];
  /* get areas */
  CROSS(OB, OC, crA);
  a = sqrt(DOT(crA, crA));
  CROSS(OC, OA, crB);
  b = sqrt(DOT(crB, crB));
  CROSS(OA, OB, crC);
  c = sqrt(DOT(crC, crC));
  d = a + b + c;
  printf(" Initial areas TOT %lf  a %lf  b  %lf  c %lf \n ratios %lf %lf %lf (expected %lf %lf %lf\n",
           d, a, b, c, a / d, b/d, c/d, w1, w2, w3);
#endif
  tol   = 1.e-12;
  CB[0] = pB[0] - pC[0]; CB[1] = pB[1] - pC[1]; CB[2] = pB[2] - pC[2];
  AC[0] = pC[0] - pA[0]; AC[1] = pC[1] - pA[1]; AC[2] = pC[2] - pA[2];
  BA[0] = pA[0] - pB[0]; BA[1] = pA[1] - pB[1]; BA[2] = pA[2] - pB[2];
  for (it = 0; it < nT; it++) {
     /* Get vectors */
     OA[0] = pA[0] - pIT[0]; OA[1] = pA[1] - pIT[1]; OA[2] = pA[2] - pIT[2];
     OB[0] = pB[0] - pIT[0]; OB[1] = pB[1] - pIT[1]; OB[2] = pB[2] - pIT[2];
     OC[0] = pC[0] - pIT[0]; OC[1] = pC[1] - pIT[1]; OC[2] = pC[2] - pIT[2];
     printf(" ---------------------------------------\n");
     printf("%lf %lf %lf O\n", pIT[0], pIT[1], pIT[2]);
     printf("%lf %lf %lf A\n\n\n",  pA[0],  pA[1], pA[2]);
     printf("%lf %lf %lf O\n", pIT[0], pIT[1], pIT[2]);
     printf("%lf %lf %lf B\n\n\n",  pB[0],  pB[1], pB[2]);
     printf("%lf %lf %lf O\n", pIT[0], pIT[1], pIT[2]);
     printf("%lf %lf %lf C\n\n\n",  pC[0],  pC[1], pC[2]);
     printf("%lf %lf %lf A\n",  pA[0],  pA[1], pA[2]);
     printf("%lf %lf %lf B\n",  pB[0],  pB[1], pB[2]);
     printf("%lf %lf %lf C\n",  pC[0],  pC[1], pC[2]);
     printf("%lf %lf %lf A\n",  pA[0],  pA[1], pA[2]);
     printf(" ---------------------------------------\n");
     /* get areas */
    CROSS(OB, OC, crA);
    CROSS(OC, OA, crB);
    CROSS(OA, OB, crC);
    da = 1.0 + uvIT[2];
    db = 1.0 - ka * uvIT[2] - kc * uvIT[3];
    dc = 1.0 + uvIT[3];
    a  = DOT(crA, crA);
    b  = DOT(crB, crB);
    c  = DOT(crC, crC);
#ifdef DEBUG
    d = a + b + c;
    printf("IT %d TOT %lf A B C %lf %lf %lf RATIOS %lf %lf %lf %lf \n",
    it, d, a,b,c, a/b, ka,  c/b, kc);
#endif
    /* Evaluate fn */
    /* a = w1 (a + b + c)
     * b = w2 (a + b + c)
     * c = w3 (a + b + c)
     * a + b + c = b / w2
     * a = w1/w2 * b; c = w3 / w2 *  b */
    grad[0] = pIT[3]; grad[1] = pIT[4]; grad[2] = pIT[5];
    CROSS(grad, CB, duA);
    CROSS(grad, AC, duB);
    CROSS(grad, BA, duC);
    grad[0] = pIT[6]; grad[1] = pIT[7]; grad[2] = pIT[8];
    CROSS(grad, CB, dvA);
    CROSS(grad, AC, dvB);
    CROSS(grad, BA, dvC);
    printf(" CROSS A %lf %lf %lf \n", crA[0], crA[1], crA[2]);
    printf(" du    A %lf %lf %lf \n", duA[0], duA[1], duA[2]);
    printf(" CROSS B %lf %lf %lf \n", crB[0], crB[1], crB[2]);
    printf(" du    A %lf %lf %lf \n", duB[0], duB[1], duB[2]);
    printf(" CROSS C %lf %lf %lf \n", crC[0], crC[1], crC[2]);
    printf(" du    A %lf %lf %lf \n", duC[0], duC[1], duC[2]);

    L[0] = DOT(duA, crA) * da + DOT(duB, crB) * db +
           DOT(duC, crC) * dc;
    L[1] = DOT(dvA, crA) * da + DOT(dvB, crB) * db +
           DOT(dvC, crC) * dc;
    L[2] = 0.5 * (a - ka * b);
    L[3] = 0.5 * (c - kc * b);
    printf("it %d ----- L %lf %lf %lf %lf  \n",it,  L[0], L[1], L[2], L[3]);
   /*
    SurrealS <4> p0[9], lambda, mu, LS[4];
    lambda = uvIT[2]; mu = uvIT[3];
    p0[0]  = pIT[0]; p0[1] = pIT[1];p0[2] = pIT[2];
    p0[3]  = pIT[3]; p0[4] = pIT[4];p0[5] = pIT[5];
    p0[6]  = pIT[6]; p0[7] = pIT[7];p0[8] = pIT[8];

    p0[0].deriv(0) = pIT[3]; //x u
    p0[1].deriv(0) = pIT[4]; //y u
    p0[2].deriv(0) = pIT[5]; //z u

    p0[0].deriv(1) = pIT[6]; // x v
    p0[1].deriv(1) = pIT[7];// y v
    p0[2].deriv(1) = pIT[8];// z v

    p0[3].deriv(0) = pIT[9]; //x uu
    p0[4].deriv(0) = pIT[10]; // y uu
    p0[5].deriv(0) = pIT[11];// z uu

    p0[3].deriv(1) = pIT[12]; // x uv
    p0[4].deriv(1) = pIT[13]; // y uv
    p0[5].deriv(1) = pIT[14]; // z uv

    p0[6].deriv(0) = pIT[12]; // x vu
    p0[7].deriv(0) = pIT[13]; // y vu
    p0[8].deriv(0) = pIT[14]; // z vu

    p0[6].deriv(1) = pIT[15]; //x vv
    p0[7].deriv(1) = pIT[16]; // y vv
    p0[8].deriv(1) = pIT[17];// z vv

    lambda.deriv(2) = 1.0;
    mu.deriv(3)     = 1.0;

    EG_minL(w1, w2, w3, lambda, mu, p0, pA, pB, pC, LS);

    L[0] = LS[0].value();
    L[1] = LS[1].value();
    L[2] = LS[2].value();
    L[3] = LS[3].value();
    double JS[16];
    JS[0] = LS[0].deriv(0);
    JS[1] = LS[0].deriv(1);
    JS[2] = LS[0].deriv(2);
    JS[3] = LS[0].deriv(3);

    JS[4] = LS[1].deriv(0);
    JS[5] = LS[1].deriv(1);
    JS[6] = LS[1].deriv(2);
    JS[7] = LS[1].deriv(3);

    JS[8]  = LS[2].deriv(0);
    JS[9]  = LS[2].deriv(1);
    JS[10] = LS[2].deriv(2);
    JS[11] = LS[2].deriv(3);

    JS[12] = LS[3].deriv(0);
    JS[13] = LS[3].deriv(1);
    JS[14] = LS[3].deriv(2);
    JS[15] = LS[3].deriv(3);

*/


    printf("it %d ----- L %lf %lf %lf %lf  \n",it,  L[0], L[1], L[2], L[3]);
    J[2 ] = DOT(duA, crA) - DOT(duB, crB) * ka;
    J[3 ] = DOT(duC, crC) - DOT(duB, crB) * kc;
    J[6 ] = DOT(dvA, crA) - DOT(dvB, crB) * ka;
    J[7 ] = DOT(dvC, crC) - DOT(dvB, crB) * kc;
    J[8 ] = J[2];
    J[12] = J[3];
    J[9 ] = J[6];
    J[13] = J[7];
    J[10] = J[11] = 0.0;
    J[14] = J[15] = 0.0;
    J[0]  = DOT(duA, duA) * da + DOT(duB, duB) * db +
            DOT(duC, duC) * dc;
    J[1]  = DOT(duA, dvA) * da + DOT(duB, dvB) * db +
            DOT(duC, dvC) * dc;
    J[5]  = DOT(dvA, dvA) * da + DOT(dvB, dvB) * db +
            DOT(dvC, dvC) * dc;
    grad[0] = pIT[9]; grad[1] = pIT[10]; grad[2] = pIT[11];
    CROSS(grad, CB, duA);
    CROSS(grad, AC, duB);
    CROSS(grad, BA, duC);
    J[0]   += DOT(duA, crA) * da + DOT(duB, crB) * db +
              DOT(duC, crC) * dc;
    grad[0] = pIT[12]; grad[1] = pIT[13]; grad[2] = pIT[14];
    CROSS(grad, CB, duA);
    CROSS(grad, AC, duB);
    CROSS(grad, BA, duC);
    J[1]    += DOT(duA, crA) * da + DOT(duB, crB) * db +
               DOT(duC, crC) * dc;
    grad[0]  = pIT[15]; grad[1] = pIT[16]; grad[2] = pIT[17];
    CROSS(grad, CB, duA);
    CROSS(grad, AC, duB);
    CROSS(grad, BA, duC);
    J[5] += DOT(duA, crA) * da + DOT(duB, crB) * db +
            DOT(duC, crC) * dc;
    J[4]  = J[1];
#ifdef DEBUG
    printf("\n JACOBIAN MATRIX 4 x 4 \n");
    printf(" %lf %lf %lf %lf\n",J[0], J[1], J[2], J[3]);
    printf(" %lf %lf %lf %lf\n",J[4], J[5], J[6], J[7]);
    printf(" %lf %lf %lf %lf\n",J[8], J[9], J[10], J[11]);
    printf(" %lf %lf %lf %lf\n",J[12], J[13], J[14], J[15]);
    printf(" ----------------------\n");
#endif
/*#ifdef DEBUG
    printf("\n SURREAL JACOBIAN MATRIX 4 x 4 \n");
    printf(" %lf %lf %lf %lf\n",JS[0], JS[1], JS[2], JS[3]);
    printf(" %lf %lf %lf %lf\n",JS[4], JS[5], JS[6], JS[7]);
    printf(" %lf %lf %lf %lf\n",JS[8], JS[9], JS[10], JS[11]);
    printf(" %lf %lf %lf %lf\n",JS[12], JS[13], JS[14], JS[15]);
    printf(" ----------------------\n");
#endif*/
    /* Directly invert 4x4 matrix */
    JI[0 ] = J[6] * J[11] * J[13] - J[ 7] * J[10] * J[13] + J[ 7] * J[ 9] * J[14] -
             J[5] * J[11] * J[14] - J[ 6] * J[ 9] * J[15] + J[ 5] * J[10] * J[15];
    JI[1 ] = J[3] * J[10] * J[13] - J[ 2] * J[11] * J[13] - J[ 3] * J[ 9] * J[14] +
             J[1] * J[11] * J[14] + J[ 2] * J[ 9] * J[15] - J[ 1] * J[10] * J[15];
    JI[2 ] = J[2] * J[ 7] * J[13] - J[ 3] * J[ 6] * J[13] + J[ 3] * J[ 5] * J[14] -
             J[1] * J[ 7] * J[14] - J[ 2] * J[ 5] * J[15] + J[ 1] * J[ 6] * J[15];
    JI[3 ] = J[3] * J[ 6] * J[ 9] - J[ 2] * J[ 7] * J[ 9] - J[ 3] * J[ 5] * J[10] +
             J[1] * J[ 7] * J[10] + J[ 2] * J[ 5] * J[11] - J[ 1] * J[ 6] * J[11];
    JI[4 ] = J[7] * J[10] * J[12] - J[ 6] * J[11] * J[12] - J[ 7] * J[ 8] * J[14] +
             J[4] * J[11] * J[14] + J[ 6] * J[ 8] * J[15] - J[ 4] * J[10] * J[15];
    JI[5 ] = J[2] * J[11] * J[12] - J[ 3] * J[10] * J[12] + J[ 3] * J[ 8] * J[14] -
             J[0] * J[11] * J[14] - J[ 2] * J[ 8] * J[15] + J[ 0] * J[10] * J[15];
    JI[6 ] = J[3] * J[ 6] * J[12] - J[ 2] * J[ 7] * J[12] - J[ 3] * J[4] * J[14] +
             J[0] * J[ 7] * J[14] + J[ 2] * J[ 4] * J[15] - J[ 0] * J[6] * J[15];
    JI[7 ] = J[2] * J[ 7] * J[ 8] - J[ 3] * J[ 6] * J[ 8] + J[ 3] * J[4] * J[10] -
             J[0] * J[ 7] * J[10] - J[ 2] * J[ 4] * J[11] + J[ 0] * J[6] * J[11];
    JI[8 ] = J[5] * J[11] * J[12] - J[ 7] * J[ 9] * J[12] + J[ 7] * J[8] * J[13] -
             J[4] * J[11] * J[13] - J[ 5] * J[ 8] * J[15] + J[ 4] * J[9] * J[15];
    JI[9 ] = J[3] * J[ 9] * J[12] - J[ 1] * J[11] * J[12] - J[ 3] * J[8] * J[13] +
             J[0] * J[11] * J[13] + J[ 1] * J[ 8] * J[15] - J[ 0] * J[9] * J[15];
    JI[10] = J[1] * J[ 7] * J[12] - J[ 3] * J[ 5] * J[12] + J[ 3] * J[4] * J[13] -
             J[0] * J[ 7] * J[13] - J[ 1] * J[ 4] * J[15] + J[ 0] * J[5] * J[15];
    JI[11] = J[3] * J[ 5] * J[ 8] - J[ 1] * J[ 7] * J[ 8] - J[ 3] * J[4] * J[ 9] +
             J[0] * J[ 7] * J[ 9] + J[ 1] * J[ 4] * J[11] - J[ 0] * J[5] * J[11];
    JI[12] = J[6] * J[ 9] * J[12] - J[ 5] * J[10] * J[12] - J[ 6] * J[8] * J[13] +
             J[4] * J[10] * J[13] + J[ 5] * J[ 8] * J[14] - J[ 4] * J[9] * J[14];
    JI[13] = J[1] * J[10] * J[12] - J[ 2] * J[ 9] * J[12] + J[ 2] * J[8] * J[13] -
             J[0] * J[10] * J[13] - J[ 1] * J[ 8] * J[14] + J[ 0] * J[9] * J[14];
    JI[14] = J[2] * J[ 5] * J[12] - J[ 1] * J[ 6] * J[12] - J[ 2] * J[4] * J[13] +
             J[0] * J[ 6] * J[13] + J[ 1] * J[ 4] * J[14] - J[ 0] * J[5] * J[14];
    JI[15] = J[1] * J[ 6] * J[ 8] - J[ 2] * J[ 5] * J[ 8] + J[ 2] * J[4] * J[ 9] -
             J[0] * J[ 6] * J[ 9] - J[ 1] * J[ 4] * J[10] + J[ 0] * J[5] * J[10];
    det    = J[3] * J[ 6] * J[ 9] * J[12] - J[ 2] * J[ 7] * J[ 9] * J[12] -
             J[3] * J[ 5] * J[10] * J[12] + J[ 1] * J[ 7] * J[10] * J[12] +
             J[2] * J[ 5] * J[11] * J[12] - J[ 1] * J[ 6] * J[11] * J[12] -
             J[3] * J[ 6] * J[ 8] * J[13] + J[ 2] * J[ 7] * J[ 8] * J[13] +
             J[3] * J[ 4] * J[10] * J[13] - J[ 0] * J[ 7] * J[10] * J[13] -
             J[2] * J[ 4] * J[11] * J[13] + J[ 0] * J[ 6] * J[11] * J[13] +
             J[3] * J[ 5] * J[ 8] * J[14] - J[ 1] * J[ 7] * J[ 8] * J[14] -
             J[3] * J[ 4] * J[ 9] * J[14] + J[ 0] * J[ 7] * J[ 9] * J[14] +
             J[1] * J[ 4] * J[11] * J[14] - J[ 0] * J[ 5] * J[11] * J[14] -
             J[2] * J[ 5] * J[ 8] * J[15] + J[ 1] * J[ 6] * J[ 8] * J[15] +
             J[2] * J[ 4] * J[ 9] * J[15] - J[ 0] * J[ 6] * J[ 9] * J[15] -
             J[1] * J[ 4] * J[10] * J[15] + J[ 0] * J[ 5] * J[10] * J[15];
    det    = 1.0 / det;
    #ifdef DEBUG
        printf("\n INVERSE JACOBIAN MATRIX 4 x 4 DET %lf ORI DET %1.12e \n", det,  1.0 / det);
        printf(" %lf %lf %lf %lf\n",JI[0]*det, JI[1]*det, JI[2]*det, JI[3]*det);
        printf(" %lf %lf %lf %lf\n",JI[4]*det, JI[5]*det, JI[6]*det, JI[7]*det);
        printf(" %lf %lf %lf %lf\n",JI[8]*det, JI[9]*det, JI[10]*det, JI[11]*det);
        printf(" %lf %lf %lf %lf\n",JI[12]*det, JI[13]*det, JI[14]*det, JI[15]*det);
        printf(" ----------------------\n");
    #endif
    solveLU (4, J, L, delta);
    printf(" START DELTA %lf %lf %lf %lf\n", delta[0], delta[1], delta[2], delta[3]);
    delta[0]  = det * (JI[0 ] * L[0] + JI[1 ] * L[1] + JI[2 ] * L[2] + JI[ 3] * L[3]);
    delta[1]  = det * (JI[4 ] * L[0] + JI[5 ] * L[1] + JI[6 ] * L[2] + JI[ 7] * L[3]);
    delta[2]  = det * (JI[8 ] * L[0] + JI[9 ] * L[1] + JI[10] * L[2] + JI[11] * L[3]);
    delta[3]  = det * (JI[12] * L[0] + JI[13] * L[1] + JI[14] * L[2] + JI[15] * L[3]);
    printf(" DELTA %lf %lf %lf %lf\n", delta[0], delta[1], delta[2], delta[3]);
    uvIT [4]  =  uvIT[0];
    uvIT [5]  =  uvIT[1];
    uvIT [6]  =  uvIT[2];
    uvIT [7]  =  uvIT[3];
    uvIT [0] -= delta[0];
    uvIT [1] -= delta[1];
    uvIT [2] -= delta[2];
    uvIT [3] -= delta[3];
    x2        = sqrt(delta[0] * delta[0] + delta[1] * delta[1] +
                     delta[2] * delta[2] + delta[3] * delta[3]);
    if      (it == 0) x0 = x2;
    else if (it == 1) x1 = x2;
    else {
        e1 = fabs(x1 / x0);
        e2 = fabs(x2 / x1);
        x0 = x1;
        x1 = x2;
      }
#ifdef DEBUG
    printf(" Xn = (%lf %lf %lf %lf )\n L = %lf %lf %lf %lf DELTA %lf %lf %lf %lf SIZE  %1.8e < %1.8e\n",
    uvIT[0], uvIT[1], uvIT[2], uvIT[3], L[0], L[1], L[2], L[3],
   delta[0], delta[1], delta[2], delta[3],x2, tol );
#endif
    i = EG_evaluate(face, uvIT, pIT);
    if (i != EGADS_SUCCESS || x2 < tol ) {
        printf("EG_evaluate %d  DELTA SIZE %1.2e < %1.2e \n",
        i, x2, tol);
        break;
    }
  }
  #ifdef DEBUG
  printf(" \n\n Found point %lf %lf %lf %lf %lf\n",
           pIT[0], pIT[1], pIT[2], uvIT[0], uvIT[1]);
  if (e1 > 0.0 && e2 > 0.0)
       printf(" Newton iterations %d Rate %lf \n", it, log(e2) / log(e1));
  OA[0] = pA[0] - pIT[0];
  OA[1] = pA[1] - pIT[1];
  OA[2] = pA[2] - pIT[2];
  OB[0] = pB[0] - pIT[0];
  OB[1] = pB[1] - pIT[1];
  OC[0] = pC[0] - pIT[0];
  OC[1] = pC[1] - pIT[1];
  OC[2] = pC[2] - pIT[2];

  /* get areas */
  CROSS(OB, OC, crA);
  a = 0.5 * sqrt(DOT(crA, crA));
  CROSS(OC, OA, crB);
  b = 0.5 * sqrt(DOT(crB, crB));
  CROSS(OA, OB, crC);
  c = 0.5 * sqrt(DOT(crC, crC));
  d = a + b + c;
  printf(" Final areas a %lf  b  %lf  c %lf \n ratios %lf %lf %lf (expected %lf %lf %lf\n",
           a, b, c, a / d, b/d, c/d, w1, w2, w3);
#endif
  uvOUT[0] = uvIT[0]; uvOUT[1] = uvIT[1];
  if ( it > 15) exit(1);
  fprintf(stderr, "IT %d RATE %lf\n", it,log(e2) / log(e1) );
  return;
}







static void nodesWeights(double *w, double *s, double *t)
{
  // degree 4 polynomial; 6 points
  // ref: PXQuadTriangle (PXE_QuadRule3), quadrature_15oct14.nb
    double x, wi;
    x = 4./9. - sqrt(10)/18. + sqrt(950 - 220*sqrt(10))/90.;
    wi = 2.*(1./12. + 3*sqrt(95 - 22*sqrt(10))/496. - sqrt(950 - 220*sqrt(10))/7440.);
    t[0] = x;
    s[0] = x;
    w[0] = wi;

    t[1] = x;
    s[1] = 1 - 2*x;
    w[1] = wi;

    t[2] = 1 - 2*x;
    s[2] = x;
    w[2] = wi;

    x    = 4./9. - sqrt(10)/18. - sqrt(950 - 220*sqrt(10))/90.;
    wi   = 2.*(1./12. - 3*sqrt(95 - 22*sqrt(10))/496. + sqrt(950 - 220*sqrt(10))/7440.);
    t[3] = x;
    s[3] = x;
    w[3] = wi;

    t[4] = x;
    s[4] = 1 - 2*x;
    w[4] = wi;

    t[5] = 1 - 2*x;
    s[5] = x;
    w[5] = wi;
    printf("%.16f %.16f %.16f\n", w[0], s[0], t[0]);
    printf("%.16f %.16f %.16f\n", w[1], s[1], t[1]);
    printf("%.16f %.16f %.16f\n", w[2], s[2], t[2]);
    printf("%.16f %.16f %.16f\n", w[3], s[3], t[3]);
    printf("%.16f %.16f %.16f\n", w[4], s[4], t[4]);
    printf("%.16f %.16f %.16f\n", w[5], s[5], t[5]);
}



void phiEvaluate(double s, double t, double phi[]) {
  // phi
  phi[0] =  s*(-1.1E1/2.0)-t*(1.1E1/2.0)+s*t*1.8E1-s*(t*t)*(2.7E1/2.0)-(s*s)*t*(2.7E1/2.0)+(s*s)*9.0-
              (s*s*s)*(9.0/2.0)+(t*t)*9.0-(t*t*t)*(9.0/2.0)+1.0;

  phi[1] =  s-(s*s)*(9.0/2.0)+(s*s*s)*(9.0/2.0);

  phi[2] =  t-(t*t)*(9.0/2.0)+(t*t*t)*(9.0/2.0);
// Want this ?
  phi[3] =  s*t*(-9.0/2.0)+(s*s)*t*(2.7E1/2.0);
// or this
  phi[4] =  s*t*(-9.0/2.0)+s*(t*t)*(2.7E1/2.0);
////
  phi[5] =  t*(-9.0/2.0)+s*t*(9.0/2.0)-s*(t*t)*(2.7E1/2.0)+(t*t)*1.8E1-(t*t*t)*(2.7E1/2.0);

  phi[6] =  t*9.0-s*t*(4.5E1/2.0)+s*(t*t)*2.7E1+(s*s)*t*(2.7E1/2.0)-(t*t)*(4.5E1/2.0)+(t*t*t)*(2.7E1/2.0);

  phi[7] =  s*9.0-s*t*(4.5E1/2.0)+s*(t*t)*(2.7E1/2.0)+(s*s)*t*2.7E1-(s*s)*(4.5E1/2.0)+(s*s*s)*(2.7E1/2.0);

  phi[8] =  s*(-9.0/2.0)+s*t*(9.0/2.0)-(s*s)*t*(2.7E1/2.0)+(s*s)*1.8E1-(s*s*s)*(2.7E1/2.0);

  phi[9] =  s*t*2.7E1-s*(t*t)*2.7E1-(s*s)*t*2.7E1;
}

static void EG_minL2insert(ego face, double *tuv) { // assuming tuv is length 9 x 2
  /* code */
  int i, j, Q = 6, nv = 10, kx = 18 * 9, ku = 18, it, nT = 100;
  double  txyz[18 * 10], s[6], t[6], w[6], L[2], J[4], delta[2], nrm1,
          nrm2, uvj[2], Puvj[18], v[3], duu[3], duv[3], dvv[3], du[3],
          dv[3], xp[3], phi[10], x0 = 0.0, x1 = 0.0, x2 = 0.0, e1 = 0.0,
          e2 = 0.0, EPS10 = 1.e-10;

  // TWEAK TRIANGLE //
  xp[0] = 0.0; xp[1] =0.0; xp[2] = 0.0;
  EG_invEvaluate(face, xp, tuv, txyz);
  xp[0] = 1.0; xp[1] =0.0; xp[2] = 0.0;
  EG_invEvaluate(face, xp,   &tuv[2], txyz);
  xp[0] = 0.5; xp[1] = 0.5; xp[2] = 0.0;
  EG_invEvaluate(face, xp,   &tuv[4], txyz);
  xp[0] = 0.8; xp[1] = 0.2; // 3
  EG_invEvaluate(face, xp,   &tuv[6], txyz);
  xp[0] = 0.6; xp[1] = 0.4; // 3
  EG_invEvaluate(face, xp,   &tuv[8], txyz);
  xp[0] = 0.4; xp[1] = 0.4; // 3
  EG_invEvaluate(face, xp,   &tuv[10], txyz);
  xp[0]   = 0.2; xp[1] = 0.2; // 3
  EG_invEvaluate(face, xp,   &tuv[12], txyz);
  xp[0]   = 0.3; xp[1] = 0.0; // 3
  EG_invEvaluate(face, xp,   &tuv[14], txyz);
  xp[0]   = 0.6; xp[1] = 0.0; // 3
  EG_invEvaluate(face, xp,   &tuv[16], txyz);
  tuv[ku] = tuv[ku + 1] = 0.0;
  for (j = i = 0; i < 9; i++) {
    if (EG_evaluate(face, &tuv[2 *i], &txyz[18 * i]) != EGADS_SUCCESS) {
      printf(" EG_evaluate vertex %d stat %d !!\n", i + 1,
               EG_evaluate(face, &tuv[2 *i], &txyz[18 * i]));
               return;
    }
    printf("%lf %lf %lf %d %lf %lf\n", txyz[18*i], txyz[18*i + 1], txyz[18*i + 2], i + 1, tuv[2 * i], tuv[2 * i + 1] );
    nrm1 = txyz[18 * i + 3] * txyz[18 * i + 3] + txyz[18 * i + 4] * txyz[18 * i + 4] +
           txyz[18 * i + 5] * txyz[18 * i + 5];
    nrm2 = txyz[18 * i + 6] * txyz[18 * i + 6] + txyz[18 * i + 7] * txyz[18 * i + 7] +
           txyz[18 * i + 8] * txyz[18 * i + 8];
    if (sqrt(nrm1) < 1.e-09 || sqrt(nrm2) < 1.e-09) continue;
    tuv[ku] += tuv[2 * i]; tuv[ku + 1] += tuv[2 * i + 1];
    j++;
  }
  tuv[ku] /= (double)j; tuv[ku + 1] /= (double)j;
  if (EG_evaluate(face, &tuv[ku], &txyz[kx]) != EGADS_SUCCESS) {
    printf(" EG_evaluate centroid stat %d !!\n", EG_evaluate(face, &tuv[ku], &txyz[kx]));
    return;
  }
  printf("%lf %lf %lf %d\n", txyz[kx], txyz[kx + 1], txyz[kx + 2], 10);
  //nodesWeights(w, s, t);
  w[0] = 0.2233815896780115; s[0] = 0.4459484909159649; t[0] = 0.4459484909159649;
  w[1] = 0.2233815896780115; s[1] = 0.1081030181680702; t[1] = 0.4459484909159649;
  w[2] = 0.2233815896780115; s[2] = 0.4459484909159649; t[2] = 0.1081030181680702;
  w[3] = 0.1099517436553219; s[3] = 0.0915762135097707; t[3] = 0.0915762135097707;
  w[4] = 0.1099517436553219; s[4] = 0.8168475729804585; t[4] = 0.0915762135097707;
  w[5] = 0.1099517436553219; s[5] = 0.0915762135097707; t[5] = 0.8168475729804585;

  for (int it = 0; it < nT; it++) {
     L[0] = L[1] = 0.0; // 0 = u 1 = v
     J[0] = J[1] = J[3] = 0.0;
     for (j = 0; j < Q; j++) {
       phiEvaluate(s[j], t[j], phi);
       uvj[0] = uvj[1] = 0.0; // uvj will be basis interpolated uvs for S(u,v)
       xp[0] = xp[1] = xp[2] = 0.0; // moving point = sum_i^N phi(i) * tyxz[i]
        for (i = 0 ; i < nv; i++) {
         uvj[0] += phi[i] * tuv [2 * i    ];
         uvj[1] += phi[i] * tuv [2 * i + 1];
         xp [0] += phi[i] * txyz[3 * i    ];
         xp [1] += phi[i] * txyz[3 * i + 1];
         xp [2] += phi[i] * txyz[3 * i + 2];
       }
       printf(" UV %lf  %lf  \n", uvj[0], uvj[1] );
       if (EG_evaluate(face, uvj, Puvj) != EGADS_SUCCESS) {
         printf(" EG_evaluate %d !!\n", EG_evaluate(face, uvj, Puvj));
         return;
       }
       printf("Du INTER %lf %lf %lf du T %lf %lf %lf\n",
        Puvj[3],  Puvj[4], Puvj[5], txyz[kx + 3],txyz[kx + 4],txyz[kx + 5] );
       du [0]  = Puvj[3 ] - txyz[kx +  3]; du [1] = Puvj[ 4] - txyz[kx +  4]; du [2] = Puvj[ 5] - txyz[kx +  5]; // du (S - X)
       dv [0]  = Puvj[6 ] - txyz[kx +  6]; dv [1] = Puvj[ 7] - txyz[kx +  7]; dv [2] = Puvj[ 8] - txyz[kx +  8]; // dv (S - X)
       duu[0]  = Puvj[9 ] - txyz[kx +  9]; duu[1] = Puvj[10] - txyz[kx + 10]; duu[2] = Puvj[11] - txyz[kx + 11]; // duu (S - X)
       duv[0]  = Puvj[12] - txyz[kx + 12]; duv[1] = Puvj[13] - txyz[kx + 13]; duv[2] = Puvj[14] - txyz[kx + 17]; // duv (S - X)
       dvv[0]  = Puvj[15] - txyz[kx + 15]; dvv[1] = Puvj[16] - txyz[kx + 16]; dvv[2] = Puvj[17] - txyz[kx + 14]; // dvv (S - X)
       v  [0]  = Puvj[ 0] -   xp[0];   v[1] = Puvj[ 1] -  xp[ 1];   v[2] = Puvj[ 2] -  xp[ 2]; //     (S - X)
       L  [0] += w[j] * phi[9] * DOT(du, v);
       L  [1] += w[j] * phi[9] * DOT(dv, v);
       J  [0] += w[j] * phi[9] * phi[9] * (DOT(duu, v) + DOT(du, du) );
       J  [1] += w[j] * phi[9] * phi[9] * (DOT(duv, v) + DOT(du, dv) );
       J  [3] += w[j] * phi[9] * phi[9] * (DOT(dvv, v) + DOT(dv, dv) );
     }
     J[2] = J[1];
     // need duu, duv, dvv
     nrm1 = J[0] * J[3] - J[1] * J[2];
#ifdef DEBUG
     printf(" JACOBIAN MATRIX 2 x 2 DET %1.2e\n", nrm1);
     printf(" %lf %lf \n",J[0], J[1]);
     printf(" %lf %lf \n",J[2], J[3]);
     printf(" ----------------------\n");
#endif
     if (fabs(nrm1) < 1.e-10) {
         printf(" IT %d DETERMINANT SIZE %lf  \n", it, nrm1);
         break;
     }
     nrm1         =   1.0 / nrm1;
     delta[0]     = -nrm1 * ( J[3] * L[0] - J[1] * L[1]);
     delta[1]     = -nrm1 * (-J[2] * L[0] + J[0] * L[1]);
     tuv[ku    ] += delta[0];
     tuv[ku + 1] += delta[1];
     x2           = sqrt(delta[0] * delta[0] + delta[1] * delta[1]);
      if      (it == 0) x0 = x2;
      else if (it == 1) x1 = x2;
      else {
          e1 = fabs(x1 / x0);
          e2 = fabs(x2 / x1);
          x0 = x1;
          x1 = x2;
      }
      i = EG_evaluate(face, &tuv[ku], &txyz[kx]);
      if (i != EGADS_SUCCESS || x2 < EPS10 ) break;
  }
#ifdef DEBUG
  printf(" NEW POINT %lf %lf %lf\n", txyz[kx + 0], txyz[kx + 1], txyz[kx + 2]);
  if (i != EGADS_SUCCESS) printf("EG_evaluate %d !!\n", i);
  printf("IT %d DELTA SIZE %1.2e < %1.2e e1 %lf %lf\n",
  it, x2, EPS10, e1, e2);
  if (e1 > EPS10 && e2 > EPS10)
      printf("CONVERGENCE RATE %lf \n", log(e2) / log(e1));
      #endif
 return;
}





int main (int argc, char *argv[])
{
  clock_t      start_t, end_t, total_t;
  int          stat = 0,  f , i, j, iBody, oclass, mtype, nbody, len, iA, iB, iC;
  int          atype, alen, *senses, min, FACECHOICE = -1, ntri;
  const int    *ints,  *tris, *tric, *ptype, *pindex;
  float        arg;
  double       box[6], size, params[3], time, frac, w1, w2, w3, uvOUT[2];
  const double *reals, *uvs, *xyzs;
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
      for ( f = 0 ; f < bodydata->nfaces; f++) {
      stat = EG_getTessFace(bodydata->tess, f + 1, &len,
                            &xyzs, &uvs, &ptype, &pindex, &ntri,
                            &tris, &tric);
       for (i = 0 ; i < ntri; i++) {
         iA   = tris[3 * i     ] - 1;
         iB   = tris[3 * i  + 1] - 1;
         iC   = tris[3 * i  + 2] - 1;
         w1   = 1.0 / 6.0; w2 = 2.0 / 6.0; w3 = 3.0 / 6.0;
         double uvIN[18 * 10];
         EG_minL2insert(bodydata->faces[f], uvIN);
         //EG_baricentricSplit(bodydata->faces[f], w1, w2, w3, &uvs[2 * iA],
          //                   &uvs[2 * iB], &uvs[2 * iC], uvOUT);
       }
     }
  }
  if (stat != EGADS_SUCCESS ) printf(" EG_main stat in cleanup %d !!\n", stat);
  for (iBody = 0; iBody < nbody; iBody++) {
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
