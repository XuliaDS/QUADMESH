#define EPS 1.e-10
// Find the point P = S(u,v) midpoint of uv0 and uv1 that minimizes the distance between
// min 1/2 ( l_0 ^ 2 + l_1 ^ 2) + lambda * (l_1 - l_0) = L(u, v, lambda)
// grad (L) = (L_1, L_2, L_3) = (0, 0, 0) --> solution using Newton x_n+1 = x_n + delta_n (3x3 system)
static void EG_minArcFact(ego face, double *uv0, double *uv1, double *uv)
{
  int    i, it, nT = 100;
  double tol, l[2], p0[18], p1[18], pn[18], J[2][2], ATJ[2][2], fact,
         r[2], b, delta[2], L[2], detJ, x0, x1, x2,  e1, e2, uvn[2];

  //Initial guess uv = 0.5 (uv0 + uv1)
  uvn[0] =  0.5 * (uv0[0] + uv1[0]);
  uvn[1] =  0.5 * (uv0[1] + uv1[1]);
  printf(" INITIAL GUESS %lf %lf %lf\n", uvn[0], uvn[1], uvn[2]);
  i      = EG_evaluate(face, uv0, p0);
  i     += EG_evaluate(face, uv1, p1);
  i     += EG_evaluate(face, uvn, pn);
  if (i != EGADS_SUCCESS) {
   printf(" EG_minArc :: EG_evaluate %d !!\n ", i);
   return;
  }


  #ifdef DEBUG
  printf(" MIN ARC FACT  \n Looking for midpoint between %lf %lf  and %lf  %lf  \n",
  uv0[0], uv0[1], uv1[0], uv1[1]);
  printf(" %lf %lf %lf \n %lf %lf %lf\n",
  p0[0], p0[1], p0[2], p1[0], p1[1], p1[2]);
  #endif // DEBUG

  tol  = 1.e-14;
  fact = 0.25;
  for (it = 0; it < nT; it++) {

#ifdef DEBUG
  printf(" \n\n NEW POINT \n %lf %lf %lf %lf %lf\n",
  pn[0], pn[1], pn[2], uvn[0], uvn[1]);
/*  printf(" du  %lf %lf %lf\n dv  %lf %lf %lf\n"
  " duu %lf %lf %lf\n duv %lf %lf %lf\n"
  " dvv %lf %lf %lf\n",   pn [3], pn [4], pn[5],
  pn [6], pn [7], pn [8], pn [9], pn[10], pn[11],
  pn[12], pn[13], pn[14], pn[15], pn[16], pn[17]);*/
#endif

  r[0] = pn[0] - p0[0];
  r[1] = pn[1] - p0[1];
  r[2] = pn[2] - p0[2];
  l[0] = sqrt(DOT(r, r));

  r[0] = pn[0] - p1[0];
  r[1] = pn[1] - p1[1];
  r[2] = pn[2] - p1[2];

  l[1] = sqrt(DOT(r, r));
#ifdef DEBUG
printf("-------------   ARC LENGTHS L %lf %lf  \n",
        l[0], l[1]);
#endif

  r[0]    = pn[0] - fact * p0[0] + (fact - 1.0) * p1[0];
  r[1]    = pn[1] - fact * p0[1] + (fact - 1.0) * p1[1];
  r[2]    = pn[2] - fact * p0[2] + (fact - 1.0) * p1[2];

  printf(" R1 %lf %lf %lf\n", r[0], r[1], r[2]);
  b       =       r[0] * pn[9 ] +  r[1] * pn[10] +  r[2] * pn[11];
  J[0][0] = b + (pn[3] * pn[3 ] + pn[4] * pn[4 ] + pn[5] * pn[5 ]);

  b       =       r[0] * pn[12] +  r[1] * pn[13] +  r[2] * pn[14];
  J[0][1] = b + (pn[3] * pn[6 ] + pn[4] * pn[7 ] + pn[5] * pn[8 ]);

  J[1][0] = J[0][1];

  b       =       r[0] * pn[15] +  r[1] * pn[16] +  r[2] * pn[17];
  J[1][1] = b + (pn[6] * pn[6 ] + pn[7] * pn[7 ] + pn[8] * pn[8 ]);

  // Solve Linear System: J * delta = - f_n
  // For now: Invert Jacobian directly --> delta = J^-1 f_n
  detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
#ifdef DEBUGG
  printf(" JACOBIAN MATRIX 2 x 2 DET %1.2e\n", detJ);
  printf(" %lf %lf \n",J[0][0], J[0][1]);
  printf(" %lf %lf \n",J[1][0], J[1][1]);
  printf(" ----------------------\n");
#endif
  if (fabs(detJ) < 1.e-08) {
      printf(" DETERMINANT SIZE %lf  \n", detJ);
      break;
  }
  ATJ[0][0] =   J[1][1];
  ATJ[0][1] = - J[0][1];
  ATJ[1][0] = - J[1][0];
  ATJ[1][1] =   J[0][0];

#ifdef DEBUGG
  printf(" INVERSE MATRIX 2 x 2  WITHOUT 1/DET\n");
  printf(" %lf %lf \n",ATJ[0][0], ATJ[0][1]);
  printf(" %lf %lf \n",ATJ[1][0], ATJ[1][1]);
  printf(" ----------------------\n");
#endif

  L[0]     =  r[0] * pn[3] + r[1] * pn[4] + r[2] * pn[5];
  L[1]     =  r[0] * pn[6] + r[1] * pn[7] + r[2] * pn[8];
  detJ     = 1.0 / detJ;
  delta[0] = -detJ * (ATJ[0][0] * L[0] + ATJ[0][1] * L[1]);
  delta[1] = -detJ * (ATJ[1][0] * L[0] + ATJ[1][1] * L[1]);
  uvn[0]  += delta[0];
  uvn[1]  += delta[1];
  x2        = sqrt(delta[0] * delta[0] + delta[1] * delta[1]);
  if      (it == 0) x0 = x2;
  else if (it == 1) x1 = x2;
  else {
      e1 = fabs(x1 / x0);
      e2 = fabs(x2 / x1);
      printf(" IT %d  L0 %lf %lf  \n (xn - xn-1)/(xn-1 - xn-2) %lf (xn+1 - xn)/(xn - xn-1) %lf CONVERGENCE %lf \n",
      it, L[0], L[1], e1, e2, log(e2) / log(e1));
      x0 = x1;
      x1 = x2;
  }
  i = EG_evaluate(face, uvn, pn);
  if (i != EGADS_SUCCESS || x2 < tol ) {
      printf("EG_evaluate %d  DELTA SIZE %1.2e < %1.2e \n",
      i, x2, tol);
      break;
  }

#ifdef DEBUG
   printf(" Xn = (%lf %lf) DELTA %lf %lf SIZE  %1.2e < %1.2e\n",
   uvn[0], uvn[1], delta[0], delta[1], x2, tol );
#endif
  }
  uv[0] = uvn[0];
  uv[1] = uvn[1];
  r[0] = pn[0] - p0[0];
  r[1] = pn[1] - p0[1];
  r[2] = pn[2] - p0[2];
  l[0] = sqrt(DOT(r, r));
  r[0] = pn[0] - p1[0];
  r[1] = pn[1] - p1[1];
  r[2] = pn[2] - p1[2];
  l[1] = sqrt(DOT(r, r));
  r[0] = p1[0] - p0[0];
  r[1] = p1[1] - p0[1];
  r[2] = p1[2] - p0[2];
  b    = sqrt(DOT(r, r));
  #ifdef DEBUG
  printf(" LEAVING WITH FINAL ARCS ||PM - P0|| = %lf ||PM - P1|| = %lf RATIOS %lf %lf (EXPECTED %lf %lf)\n",
           l[0], l[1], l[0] / b, l[1]/b, fact, 1.0 - fact);
   #endif

   //exit(1);

  return;
}




