static void EG_minArcN (const ego face, int n, const double *uvlist, double *uv)
{
  int    stat, i, j, i1,  it, nT = 100, m, iC, iR;
  double x0 = 0.0, x1 = 0.0, x2 = 0.0, e1 = 0.0, e2 = 0.0, pIT[18], lt;
  double *l = NULL, *plist = NULL, *dlu = NULL, *dlv = NULL, *ALU = NULL, *L = NULL, *uvIT = NULL,
          *delta = NULL, *r = NULL, *dluu = NULL, *dluv = NULL, *dlvv = NULL;

  //Initial guess uv = 0.5 (uv0 + uv1)
  printf(" CHECK WITH MIN ARC\n");
  EG_minArc(face, 0.5, uvlist, &uvlist[2], uv);
  m       = n + 1;
  uvIT    = EG_alloc(     m * sizeof(double));
  delta   = EG_alloc(     m * sizeof(double));
  L       = EG_alloc(     m * sizeof(double));
  plist   = EG_alloc(18 * n * sizeof(double));
  l       = EG_alloc(     n * sizeof(double));
  dlu     = EG_alloc(     n * sizeof(double));
  dlv     = EG_alloc(     n * sizeof(double));
  dluu    = EG_alloc(     n * sizeof(double));
  dluv    = EG_alloc(     n * sizeof(double));
  dlvv    = EG_alloc(     n * sizeof(double));
  r       = EG_alloc( 3 * n * sizeof(double));
  ALU     = EG_alloc( m * m * sizeof(double));
  if (uvIT == NULL || delta == NULL || L   == NULL || plist == NULL ||
         l == NULL ||  dlu  == NULL || dlv == NULL || dluv  == NULL ||
      dlvv == NULL ||     r == NULL || ALU == NULL ) return;
  uvIT[0] = uvIT[1] = 0.0;
  printf(" TOTAL POINTS %d INCOGINTAS 2 + %d\n", n, n - 1);
  // get all point coordinates
  for (i = 0; i < n; i++) {
      uvIT[0]  += uvlist[2 * i    ];
      uvIT[1]  += uvlist[2 * i + 1];
      stat      = EG_evaluate(face, &uvlist[2 * i], &plist[18 * i]);
      if (stat != EGADS_SUCCESS) {
          printf(" EG_minArc :: EG_evaluate %d !!\n ", stat);
          goto cleanup;
      }
      if (2 + i < m) uvIT[2 + i] = 1.0;
      printf("%lf %lf %lf %d\n", plist[18 * i], plist[18 * i + 1], plist[18 * i + 2], i);
  }
  uvIT[0]  /= (double)(n + 1);
  uvIT[1]  /= (double)(n + 1);
  for (i = 0; i < m; i++) printf(" U %lf \n", uvIT[i]);
  stat      = EG_evaluate(face, uvIT, pIT);
  if (stat != EGADS_SUCCESS) {
      printf(" EG_minArc :: EG_evaluate %d !!\n ", stat);
      goto cleanup;
  }
#ifdef DEBUG
  lt = 0.0;
  for (i = 0; i < n; i++) {
      r[3 * i    ] = pIT[0] - plist[18 * i     ];
      r[3 * i + 1] = pIT[1] - plist[18 * i  + 1];
      r[3 * i + 2] = pIT[2] - plist[18 * i  + 2];
      l[i]         = r[3 * i    ] * r[3 * i     ] + r[3 * i + 1] * r[3 * i + 1] +
                     r[3 * i + 2] * r[3 * i + 2];
      lt          += sqrt (l[i]);
  }
  printf(" \n\n --------------------------------------------------- \n");
  printf(" INITIAL GUESS ARCS DISTRIBUTION \n");
  for ( i = 0 ; i < n; i++)
      printf(" ARC l%d = %lf (%lf) RATIO %lf (%lf) EXP %lf \n",
             i, l[i], sqrt(l[i]), l[i] / (lt * lt), sqrt(l[i]) / lt, 1.0 / (double)n);
#endif
  for (it  = 0; it < nT; it++) {
      for (i = 0 ; i < n; i++) {
          j        = 3 * i;
          r[j    ] = pIT[0] - plist[18 * i      ];
          r[j + 1] = pIT[1] - plist[18 * i   + 1];
          r[j + 2] = pIT[2] - plist[18 * i   + 2];
          l    [i] = r[j] * r[j] + r[j + 1] * r[j + 1] + r[j + 2] * r[j + 2];
          dlu  [i] = 2.0 * (r[j] * pIT[3] + r[j + 1] * pIT[ 4] + r[j + 2] * pIT[5]);
          dlv  [i] = 2.0 * (r[j] * pIT[6] + r[j + 1] * pIT[ 7] + r[j + 2] * pIT[8]);

          dluu [i] = 2.0 * (r[j] * pIT[9] + r[j + 1] * pIT[10] + r[j + 2] * pIT[11] +
                     pIT[3] * pIT[3] + pIT[4] * pIT[4] + pIT[5] * pIT[5]);

          dluv [i] = 2.0 * (r[j] * pIT[12] + r[j + 1] * pIT[13] + r[j + 2] * pIT[14] +
                     pIT[3] * pIT[6] + pIT[4] * pIT[7] + pIT[5] * pIT[8]);

          dlvv [i] = 2.0 * (r[j] * pIT[15] + r[j + 1] * pIT[16] + r[j + 2] * pIT[17] +
                     pIT[6] * pIT[6] + pIT[7] * pIT[7] + pIT[8] * pIT[8]);
          for (iC = 0; iC < m; iC++) ALU[m * i + iC] = 0.0;
#ifdef DEBUG
          printf("\n\n r%d %lf %lf %lf NORM %lf dlu %lf dlv %lf\n",i,
                 r[j],r[j + 1],r[j + 2], l[i], dlu[i], dlv[i]);
#endif
      }
      for (iC = 0; iC < m; iC++) ALU[(m - 1) * m + iC] = 0.0;
      printf(" A IS %d x %d MATRIX \n", m, m);
      for (iR = 0; iR < m; iR++) {
        for (iC = 0; iC < m; iC++) {
            if (iC == m - 1) printf("A(%d) = %lf\n",iR * m + iC,ALU[iR * m + iC]);
            else             printf("A(%d) = %lf " ,iR * m + iC,ALU[iR * m + iC]);
        }
      }
      L[0] = L[1] = 0.0;
      for (i = 0; i < n; i++) {
          L  [0]     +=  dlu[i] * l[i];
          L  [1]     +=  dlv[i] * l[i];
          ALU[0]     += dluu[i] * l[i] + dlu[i] * dlu[i];
          ALU[1]     += dluv[i] * l[i] + dlu[i] * dlv[i];
          ALU[m + 1] += dlvv[i] * l[i] + dlv[i] * dlv[i];
          if (i == n - 1) break;
          i1        = i + 1;
          L[0    ] += uvIT[2 + i] * (dlu[i] - dlu[i1]);
          L[1    ] += uvIT[2 + i] * (dlv[i] - dlv[i1]);
          printf(" POS %d TOTAL %d\n", 2 + i, m - 1);
          L  [2 + i]     = l[i] - l[i1];
          ALU[0]        += uvIT[2 + i] * (dluu[i] - dluu[i1]);
          ALU[1]        += uvIT[2 + i] * (dluv[i] - dluv[i1]);
          ALU[m + 1]    += uvIT[2 + i] * (dlvv[i] - dlvv[i1]);
          ALU[    2 + i] =  dlu[i] - dlu[i1];
          ALU[m + 2 + i] =  dlv[i] - dlv[i1];
          printf(" DIAGONAL %d \n", (2 + i) * m + (2 + i));
          //ALU[(2 + i) * m + (2 + i)] = l[i] - l[i1];
      }
#ifdef DEBUG
      for (i = 0; i < m; i++)
          printf("L[%d] = %lf\n", i, L[i]);
#endif
      for (iR = 0; iR < m; iR++) {
        for (iC = 0; iC < iR; iC++) {
            printf(" ASSIGNING %d <-- %d\n",m * iR + iC, m*iC + iR );
            ALU[m * iR + iC] = ALU[m*iC + iR];
          }
      }
      printf(" A IS %d x %d MATRIX \n", m, m);
      for (iR = 0; iR < m; iR++) {
        for (iC = 0; iC < m; iC++) {
            if (iC == m - 1) printf("A(%d) = %lf\n",iR * m + iC,ALU[iR * m + iC]);
            else             printf("A(%d) = %lf " ,iR * m + iC,ALU[iR * m + iC]);
        }
      }
      solveLU(m, ALU, L, delta);
      printf(" DELTA  ");
      x2 = 0.0;
      for (i = 0; i < m; i++) {
          uvIT[i]  -= delta[i];
          printf("%lf\t", uvIT[i]);
          x2 += delta[i] * delta[i];
      }
      x2     = sqrt(x2);
      i = EG_evaluate(face, uvIT, pIT);
      if (i != EGADS_SUCCESS || x2 < EPS10 ) break;
      if      (it == 0) x0 = x2;
      else if (it == 1) x1 = x2;
      else {
          e1 = fabs(x1 / x0);
          e2 = fabs(x2 / x1);
          x0 = x1;
          x1 = x2;
      }

#ifdef DEBUGG
      printf(" NEW POINT %lf %lf %lf\n", pIT[0], pIT[1], pIT[2]);
      printf("du  %lf %lf %lf dv %lf %lf %lf\n", pIT[3], pIT[4], pIT[5], pIT[6], pIT[7], pIT[8]);
      printf("duu %lf %lf %lf duv %lf %lf %lf dvv %lf %lf %lf\n", pIT[9], pIT[10], pIT[11], pIT[12], pIT[13], pIT[14], pIT[15], pIT[16], pIT[17]);
#endif

  }
#ifdef DEBUG
  printf("\n\n--------------- REPORT ------------------------------------ \n");
  printf(" NEW POINT %lf %lf %lf\n", pIT[0], pIT[1], pIT[2]);
  if (i != EGADS_SUCCESS) printf("EG_evaluate %d !!\n", i);
  printf("IT %d DELTA SIZE %1.2e < %1.2e e1 %lf %lf\n",
  it, x2, EPS10, e1, e2);
  if (e1 > EPS10 && e2 > EPS10)
      printf("CONVERGENCE RATE %lf \n", log(e2) / log(e1));
  if (e1 > EPS10 && e2 > EPS10)
      printf("CONVERGENCE RATE %lf \n", log(e2) / log(e1));
  lt = 0.0;
  for ( i = 0 ; i < n; i++) {
      r[3 * i    ] = pIT[0] - plist[18 * i    ];
      r[3 * i + 1] = pIT[1] - plist[18 * i + 1];
      r[3 * i + 2] = pIT[2] - plist[18 * i + 2];
      l[i]         = r[3 * i    ] * r[3* i     ] + r[3 * i + 1] * r[3 * i + 1] +
                     r[3 * i + 2] * r[3 * i + 2];
      lt          += sqrt (l[i]);
  }
  printf(" \n\n --------------------------------------------------- \n");
  printf(" INITIAL GUESS ARCS DISTRIBUTION \n");
  for ( i = 0 ; i < n; i++)
      printf(" ARC l%d = %lf (%lf) RATIO %lf (%lf) EXP %lf \n",
             i, l[i], sqrt(l[i]), l[i] / (lt * lt), sqrt(l[i]) / lt, 1.0 / (double)n);
  printf(" --------------------------------------------------- \n");
#endif
  uv[0] = uvIT[0];
  uv[1] = uvIT[1];
cleanup:
  EG_free(uvIT);
  EG_free(delta);
  EG_free(L);
  EG_free(plist);
  EG_free(l);
  EG_free(dlu);
  EG_free(dlv);
  EG_free(dluu);
  EG_free(dluv);
  EG_free(dlvv);
  EG_free(r);
  EG_free(ALU);
  exit(1);
  return;
}



