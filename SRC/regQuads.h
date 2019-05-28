#include <math.h>
#include <string.h>


#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define ERFC(a,b,x) (0.5 * erfc(3.554147 * ((x -(0.5 * (a + b))) / (0.5 * (a - b)))))
// NORMALIZED erfc = 0 at x = 1 and 1 at x = -1. min is better than max (eg angles) a < b
// max is better than min (sizes) b < a

#define CROSS(a,b,c)  c[0] = (a[1]*b[2]) - (a[2]*b[1]);\
                      c[1] = (a[2]*b[0]) - (a[0]*b[2]);\
                      c[2] = (a[0]*b[1]) - (a[1]*b[0])
#define DOT(a,b)     (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

#define PI     3.1415926535897931159979635

#define qEPS   1.e-14

#define EPS08  1.E-08
#define ANGCUT 2.9

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
#define DEG175 3.1241
#define QA0  0
#define QA1  100
#define QA2  10000
#define QA3  100000000
#define QACB 1000000000  // crosses domain bounds: SUPER INVALID



typedef struct{
  int verts[4], qadj[4], id;
} Quad;


typedef struct{
  int  *verts, *quads;
  int   nV, nQ; // nV = n + 1 =  origin(1) + peaks (n)
  int  *idxV, *idxQ, *areas;
  double *angles, *ratios;
} vStar;


typedef struct {
  int      fID, oriQ, oriV, plotcount, totQ, totV, pp,
           sizeV, sizeQ, *qIdx, *qAdj, **valence,
           *vInv, *vType, *remQ, *remV, invsteps, regBd, regBd0;
  ego      face;
  double   range[4],  *xyzs, *uvs, minArea, maxArea, avArea;
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

extern int  EG_createMeshMap(bodyQuad *bodydata);
extern int  EG_meshRegularization(meshMap *qm);
extern int  EG_makeQuadTess(bodyQuad bodydata, ego *quadTess);
extern void EG_destroyMeshMap(bodyQuad *bodydata);
