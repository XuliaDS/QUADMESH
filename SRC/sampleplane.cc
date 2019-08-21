/*
void sampleNormalPlane(double normal[], double point[], int vID, meshMap *qm) {
  double min[3], max[3], c, p[3], dt[3], r = 1.0;
  int  i, j, k, nP;
  FILE *f;
  char name[32];
  snprintf(name, sizeof(char) * 32, "PLANE%i_%i.txt", vID, NC++);
  printf("QM PLOT COUNT %d PLANE %d WRITING IN %s\n",qm->plotcount, NC, name);
  f = fopen(name,"w");
  if ( f == NULL) return ;
  nP = 40;
  c = 0.0;
  vID--;
  for (i = 0 ; i < 3; ++i) {
    min  [i] = qm -> xyzs[3*vID + i]-0.25; //qm -> xyzs[3*vID + i] - 100 *qm -> xyzs[3*vID + i] ;
    max  [i] = qm -> xyzs[3*vID + i]+0.25;//0.5;// qm -> xyzs[3*vID + i] + 100 *qm -> xyzs[3*vID + i] ;
    if ( min[i] > max[i]) {
      dt[0] = min[i];
      min[i] = max[i];
      max[i] = dt[0];
    }
    c       += normal[i] * point[i];
  }
  for (i = 0 ; i < 3; ++i)
    dt[i] = (max[i] - min[i] )/(double)(nP -1);
  fprintf(f,"%lf %lf %lf\n",point[0], point[1], point[2]);
  for (i = 0 ; i < nP; ++i)
  {
    for ( j = 0 ; j < nP; ++j)
    {
      for ( k = 0 ; k < nP; ++k)
      {
        p[0] = min[0] + (double)i*dt[0];
        p[1] = min[1] + (double)j*dt[1];
        double aux = c - p[0] * normal[0] - p[1] * normal[1];
        if ( fabs(normal[2]) > qEPS) p[2] = aux / normal[2];
        else {
          p[2] = min[2] + (double)k*dt[2];
          if (fabs(normal[1]) < qEPS) {
            p[1] = min[1] + (double)j*dt[1];
            p[0] = point[0];
          }
          else {
            if (fabs(normal[0]) < qEPS) p[1] = point[1];
            else {
              aux  = c - p[0] * normal[0];
              p[1] = aux / normal[1];
            }
          }
        }
        fprintf(f,"%lf %lf %lf\n",p[0], p[1], p[2]);
      }
    }
  }
  fclose(f);
}
*/
