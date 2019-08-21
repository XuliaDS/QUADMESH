#include "egads.h"
#include <math.h>
#include <string.h>

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
static void EG_baricentricSplit(ego face, const double w1, const double w2,const  double w3,
                                const double *uvA, const double *uvB,const  double *uvC, double *uvOUT)
{
  int    i, it, nT = 100, stat;
  double tol, pA[18], pB[18], pC[18], OA[3], OB[3], OC[3], pIT[18], J[16], JI[16],
         crA[3], crB[3], crC[3], grad[3], a, b, c, d, delta[4], L[4], x0 = 0.0, ka, kc,
         x1 = 0.0, x2 = 0.0, e1 = 0.0, e2 = 0.0, uvIT[4], duA[3], duB[3], duC[3], CB[3],
         dvA[3], dvB[3], dvC[3], AC[3], BA[3], da, db, dc, det;
  //uvIT[0]   = w1 * uvA[0] + w2 * uvB[0] + w3 * uvC[0];
  //uvIT[1]   = w1 * uvA[1] + w2 * uvB[1] + w3 * uvC[1];
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
    L[0] = DOT(duA, crA) * da + DOT(duB, crB) * db +
           DOT(duC, crC) * dc;
    L[1] = DOT(dvA, crA) * da + DOT(dvB, crB) * db +
           DOT(dvC, crC) * dc;
    L[2] = 0.5 * (a - ka * b);
    L[3] = 0.5 * (c - kc * b);
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
#ifdef DEBUGG
    printf("\n JACOBIAN MATRIX 4 x 4 \n");
    printf(" %lf %lf %lf %lf\n",J[0], J[1], J[2], J[3]);
    printf(" %lf %lf %lf %lf\n",J[4], J[5], J[6], J[7]);
    printf(" %lf %lf %lf %lf\n",J[8], J[9], J[10], J[11]);
    printf(" %lf %lf %lf %lf\n",J[12], J[13], J[14], J[15]);
    printf(" ----------------------\n");
#endif
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
    #ifdef DEBUGG
        printf("\n INVERSE JACOBIAN MATRIX 4 x 4 DET %lf ORI DET %1.12e \n", det,  1.0 / det);
        printf(" %lf %lf %lf %lf\n",JI[0]*det, JI[1]*det, JI[2]*det, JI[3]*det);
        printf(" %lf %lf %lf %lf\n",JI[4]*det, JI[5]*det, JI[6]*det, JI[7]*det);
        printf(" %lf %lf %lf %lf\n",JI[8]*det, JI[9]*det, JI[10]*det, JI[11]*det);
        printf(" %lf %lf %lf %lf\n",JI[12]*det, JI[13]*det, JI[14]*det, JI[15]*det);
        printf(" ----------------------\n");
    #endif
    delta[0]  = det * (JI[0 ] * L[0] + JI[1 ] * L[1] + JI[2 ] * L[2] + JI[ 3] * L[3]);
    delta[1]  = det * (JI[4 ] * L[0] + JI[5 ] * L[1] + JI[6 ] * L[2] + JI[ 7] * L[3]);
    delta[2]  = det * (JI[8 ] * L[0] + JI[9 ] * L[1] + JI[10] * L[2] + JI[11] * L[3]);
    delta[3]  = det * (JI[12] * L[0] + JI[13] * L[1] + JI[14] * L[2] + JI[15] * L[3]);
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
#ifdef DEBUGG
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
       printf(" Final areas TOT %lf  a %lf  b  %lf  c %lf \n ratios %lf %lf %lf (expected %lf %lf %lf\n",
                d, a, b, c, a / d, b/d, c/d, w1, w2, w3);
#endif
  uvOUT[0] = uvIT[0]; uvOUT[1] = uvIT[1];
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
         EG_baricentricSplit(bodydata->faces[f], w1, w2, w3, &uvs[2 * iA],
                             &uvs[2 * iB], &uvs[2 * iC], uvOUT);
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
