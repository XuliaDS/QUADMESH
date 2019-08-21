#include "egads.h"
#include <math.h>
#include <string.h>

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
  ego      tess,   *faces, *edges;
  int      nedges, nfaces;
} bodyQuad;

static void EG_getEdgePoint(ego edge, const double w, const double tm,
                                      const double tp, double *tOUT)
{
   int nT = 50, stat, i, it, nS = 20;
   double pIT[18], pm[18], pp[18], dt = 0.0, vt0[3], vt1[3], r0[3], r1[3], k, xyz[3],
          res = 1.0, t, tIT, l0, l1, lt, f, ff, e1 = -1.0, e2 = -1.0, x0 = 0.0,
          s, x1 = 0.0, x2 = 0.0, nEPS = 1.e-10;

   k     = (w * w) / ((1.0 - w ) * (1.0 - w));
   tIT   = tp * w + tm * ( 1.0 - w );
   *tOUT = tIT;
   stat  = EG_evaluate(edge, &tm, pm);
   stat += EG_evaluate(edge, &tp, pp);
   stat += EG_evaluate(edge, &tIT, pIT);
   if ( stat != EGADS_SUCCESS) {
     printf("EG_getEdgePoint before Newton: EG_evaluate %d\n ", stat);
     *tOUT = t;
     return;
   }
#ifdef DEBUG
   r0[0] = pIT[0] - pm[0]; r1[0] = pIT[0] - pp[0];
   r0[1] = pIT[1] - pm[1]; r1[1] = pIT[1] - pp[1];
   r0[2] = pIT[2] - pm[2]; r1[2] = pIT[2] - pp[2];
   l0    = sqrt(DOT(r0, r0));
   l1    = sqrt(DOT(r1, r1));
   lt    = l0 + l1;
   printf("\n Initial l0 = %lf l1 = %lf split edge in fracs %lf %lf (expected %lf %lf)\n",
   l0, l1, l0 / lt, l1 / lt, w, 1.0 - w);
#endif
   for ( it = 0 ; it < nT; it++) {
     /* perform a line search such that residual decreases from the previous iteration */
     s = 1.0; f = 0.0;
     for (i = 0; i < nS; i++) {
       t    = tIT + s * dt;
       stat = EG_evaluate(edge, &t, pIT);
       if (stat != EGADS_SUCCESS) {
         printf(" EG_getSidepoint: EG_evaluate = %d\n", stat);
         s /= 2.0;
         continue;
       }
       if (t < tm || t > tp ) {
         if (it > 0) {
            s /= 2.0;
            continue;
          }
          else if (it == 0) {
 #ifdef DEBUG
         printf(" t %lf OUT OF RANGE: [%lf %lf]\n",
         t, tm, tp);
 #endif
          if (i == 0) {
            pIT[0]    = w * pp[0] + (1.0 - w) * pm[0];
            pIT[1]    = w * pp[1] + (1.0 - w) * pm[1];
            pIT[2]    = w * pp[2] + (1.0 - w) * pm[2];
            stat      = EG_invEvaluateGuess(edge, pIT, &tIT, xyz);
            if (stat != EGADS_SUCCESS ||
                         t < tm || t > tp ) {
              printf(" EG_getSidepoint: EG_invEvaluateGuess out-of-range!\n");
              /* we are hopelessly out of range... */
              return;
            }
            continue; /* try again */
          } else {
             /* already tried better guess with inverse evaluate... */
             printf(" EG_getSidepoint: Initial out-of-range!\n");
             return;
           }
         }
       }
       // distance vector
       r0 [0] = pIT[0] - pm[0]; r1[0]  = pIT[0] - pp[0];
       r0 [1] = pIT[1] - pm[1]; r1[1]  = pIT[1] - pp[1];
       r0 [2] = pIT[2] - pm[2]; r1[2]  = pIT[2] - pp[2];
       l0     = DOT(r0, r0);
       l1     = DOT(r1, r1);
       f      = 0.5 * (l0 - k * l1);
       if (it == 0 || fabs(f) < res) {
         res = fabs(f); /* save off the last residual and exit */
         break;
       }  else {
         s /= 2.0; /* try again if the residual grew */
         continue;
       }
     }
     tIT    = t;
     if ( i == nS || res < nEPS ) break;
     // l0 = DOT(p - pm, p - pm) l1 = DOT(p - pp, p - pp)
     vt0[0] = pIT[3]; vt1[0] = pIT[3];
     vt0[1] = pIT[4]; vt1[1] = pIT[4];
     vt0[2] = pIT[5]; vt1[2] = pIT[5];
     ff     = DOT(vt0, r0) - k * DOT(vt1, r1);
     // dl0 = DOT(dt, p - pm) dl1 = DOT(dt - pp, p - pp)
     dt     = - (f / ff);
     x2     = fabs (dt);
     if      (it == 0) x0 = x2;
     else if (it == 1) x1 = x2;
     else {
         e1 = fabs(x1 / x0);
         e2 = fabs(x2 / x1);
         x0 = x1;
         x1 = x2;
     }
     if (fabs (dt) < nEPS) break;
   }
#ifdef DEBUG
   if (e1 >= 0.0 && e2 >= 0.0)
   printf("Newton iterations %d convergence %lf (e1 %lf e2 %lf) \n", it , log(e2) / log(e1), e1, e2);
   else printf("Newton iterations %d convergence NA \n", it);
      r0[0] = pIT[0] - pm[0]; r0[1] = pIT[1] - pm[1]; r0[2] = pIT[2] - pm[2];
      r1[0] = pIT[0] - pp[0]; r1[1] = pIT[1] - pp[1]; r1[2] = pIT[2] - pp[2];
      l0    = sqrt(DOT(r0, r0));
      l1    = sqrt(DOT(r1, r1));
      lt    = l0 + l1;
      printf(" FINAL l0 = %lf l1 = %lf split edge in fracs %lf %lf (expected %lf %lf)\n",
      l0, l1, l0 / lt, l1 / lt, w, 1.0 - w);
#endif

   if (res >= nEPS && fabs(dt) >= nEPS)
     printf(" EG_getSidepoint: not converged -- residual %1.2le delta %1.2le (%1.2le)!\n",
            res, fabs(dt), nEPS);

   /* found a solution out of range -- report it (for now) */
   if (tIT < tm || tIT > tp )
     printf(" EG_getSidepoint: solution out-of-range!\n");
   *tOUT = tIT;
}

#ifdef STANDALONE
int main (int argc, char *argv[])
{
  clock_t      start_t, end_t, total_t;
  int          stat = 0,  f , i, j, iBody, oclass, mtype, nbody, len, iA, iB, iC;
  int          atype, alen, *senses, min, FACECHOICE = -1, ntri;
  const int    *ints,  *tris, *tric, *ptype, *pindex;
  float        arg;
  double       box[6], size, params[3], time, frac, w1, w2, w3, uvOUT[2];
  const double *reals, *uvs, *xyzs, *ts;
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
                             &bodydata[iBody].nedges,&bodydata[iBody].edges);
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
      for (f = 0; f < bodydata[iBody].nedges; f++) {
        stat  = EG_getTessEdge(bodydata->tess, f+1, &len, &xyzs, &ts);
        if (stat != EGADS_SUCCESS) continue;
        if (len == 0) continue;
        for (i = 0; i < len - 1; i++) {
          w1   = 0.15;
          EG_getEdgePoint(bodydata->edges[f], w1, ts[i], ts[i+1],uvOUT);
       }
     }
  }
  if (stat != EGADS_SUCCESS ) printf(" EG_main stat in cleanup %d !!\n", stat);
  for (iBody = 0; iBody < nbody; iBody++) {
      EG_free(bodydata[iBody].faces);
      EG_free(bodydata[iBody].edges);
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
