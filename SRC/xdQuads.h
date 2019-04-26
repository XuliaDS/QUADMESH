#include <math.h>
#include <string.h>


#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))


#define AVERAGE   0
#define ARCLENGTH 1
#define INVEVAL   2
#define CROSS(a,b,c)  c[0] = (a[1]*b[2]) - (a[2]*b[1]);\
                      c[1] = (a[2]*b[0]) - (a[0]*b[2]);\
                      c[2] = (a[0]*b[1]) - (a[1]*b[0])
#define DOT(a,b)     (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

#define PI     3.1415926535897931159979635

//#define EPS11  1.E-11
#define EPS08  1.E-08

#define SWAP                0
#define COLLAPSE            1
#define SPLIT               2
#define SWAPCOLLAPSE        0
#define DOUBLECOLLAPSE      1
#define SWAPDOUBLECOLLAPSE  2
#define SWAPSPLIT           3
#define DOUBLESPLIT         4
#define SWAPDOUBLESPLIT     5
#define DOUBLESWAP          6

#define QA0  0
#define QA1  100
#define QA2  10000
#define QA3  1000000
#define QACB 100000000   // crosses domain bounds: SUPER INVALID


#define ANGPASS 3.1
typedef struct{
  int verts[4], qadj[4], id;
} Quad;

typedef struct{
  int    *area;
  double *ratio;
}mArea;

typedef struct{
  int  *verts, *quads;
  int   nV, nQ; // nV = n + 1 =  origin(1) + peaks (n)
  int  *idxV, *idxQ;
} vStar;


typedef struct {
  int      fID, oriQ, oriV, plotcount, totQ, totV,
  sizeV, sizeQ, *qIdx, *qAdj, **valence,
  *qInv, *vType, *remQ, *remV, invsteps;
  ego      face;
  double   range[4],  *xyzs, *uvs;
  vStar **star;
} meshMap;


typedef struct {
  ego      tess,   *faces;
  int      nedges, nfaces;
  meshMap  **qm;
} bodyQuad;


typedef struct{
  int  verts[6], vals[6];
  int  q[2];
} quadGroup;


extern int  EG_outLevel(const egObject *object);

extern int  EG_createMeshMap(bodyQuad *bodydata, int uvtype);
extern int  EG_meshRegularization(meshMap *qm);
extern int  EG_makeQuadTess(bodyQuad bodydata, ego *quadTess);
extern void EG_destroyMeshMap(bodyQuad *bodydata);
