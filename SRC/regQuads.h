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

#define DEG160 2.79252680319093
#define DEG10  0.1745329251994329
#define DEG20  0.349065850398866
#define PI     3.1415926535897931159979635

#define EPS11  1.E-11

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


typedef struct{
  int verts[4], qadj[4], id;
} Quad;


typedef struct{
  int  *verts, *quads;
  int   nV, nQ; // nV = n + 1 =  origin(1) + peaks (n)
  int  *idxV, *idxQ;
} vStar;


typedef struct {
  int      fID, oriQ, oriV, plotcount, uvtype, estQ, totQ,
           totV, sizeV, sizeQ, *vFix, *qIdx, *qAdj,
	    **valence, *vType, *remQ, *remV, pp, invsteps;
  ego      face;
  double   minsize, avsize, range[4], minAngle, maxAngle;
  double  *xyzs, *uvs;
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


extern int  EG_outLevel(const egObject *object );

extern int  EG_createMeshMap(bodyQuad *bodydata, int uvtype);
extern int  EG_fullMeshRegularization(meshMap *qm);
extern int  EG_makeQuadTess(bodyQuad bodydata, ego *quadTess);
extern void EG_destroymeshMap(bodyQuad *bodydata);

