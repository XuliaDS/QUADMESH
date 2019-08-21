// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

//----------------------------------------------------------------------------//
// QuadratureArea<Triangle>
// area quadrature rules for reference triangle (x in [0, 1]; y in [0, 1-x])

#include <cmath> // sqrt
#include <string>
#include <iostream>

#include "tools/SANSnumerics.h"     // Real
#include "tools/SANSException.h"
#include "Topology/ElementTopology.h"
#include "QuadratureArea.h"

namespace SANS
{

// number of order indexes
template <>
const int QuadratureArea<Triangle>::nOrderIdx = 12;

template <>
int
QuadratureArea<Triangle>::getOrderIndex( const int orderReqd )
{
  if (orderReqd < 0)          // max quadrature
    return nOrderIdx-1;
  else if (orderReqd <= 1)
    return 0;
  else if (orderReqd <= 2)
    return 1;
  else if (orderReqd <= 3)
    return 2;
  else if (orderReqd <= 4)
    return 3;
  else if (orderReqd <= 5)
    return 4;
  else if (orderReqd <= 8)
    return 5;
  else if (orderReqd <= 10)
    return 6;
  else if (orderReqd <= 11)
    return 7;
  else if (orderReqd <= 12)
    return 8;
  else if (orderReqd <= 13)
    return 9;
  else if (orderReqd <= 14)
    return 10;
  else if (orderReqd <= 30)
    return 11;
  else
    SANS_DEVELOPER_EXCEPTION("Quadrature order=%d is not available", orderReqd);

  // suppress compiler warnings
  return -1;
}

template <>
int
QuadratureArea<Triangle>::getOrderFromIndex( const int orderidx )
{
  // returns the quadrature order for a given order index into the quadrature cache
  switch (orderidx)
  {
  case 0:
    return 1;
  case 1:
    return 2;
  case 2:
    return 3;
  case 3:
    return 4;
  case 4:
    return 5;
  case 5:
    return 8;
  case 6:
    return 10;
  case 7:
    return 11;
  case 8:
    return 12;
  case 9:
    return 13;
  case 10:
    return 14;
  case 11:
    return 30;
  default:
    SANS_DEVELOPER_EXCEPTION("Please add orderidx=%d mapping", orderidx);
  }

  // suppress compiler warnings
  return -1;
}

template <>
QuadratureArea<Triangle>::QuadratureArea( const int orderReqd ) :
    order_(0), orderidx_(getOrderIndex(orderReqd)), nQuad_(0), weight_(NULL), xcoord_(NULL), ycoord_(NULL)
{
  //orderReqd: order of the polynomial that needs to be integrated exactly

  int order;
  int orderMax;         // maximum allowed order
  int nQuad;

  // set quadrature rule order

  orderMax = 30;

  if (orderReqd < 0)          // max quadrature
    order = orderMax;
  else if (orderReqd <= 1)
    order = 1;
  else if (orderReqd <= 2)
    order = 2;
  else if (orderReqd <= 3)
    order = 3;
  else if (orderReqd <= 4)
    order = 4;
  else if (orderReqd <= 5)
    order = 5;
  else if (orderReqd <= 8)
    order = 8;
  else if (orderReqd <= 10)
    order = 10;
  else if (orderReqd <= 11)
    order = 11;
  else if (orderReqd <= 12)
    order = 12;
  else if (orderReqd <= 13)
    order = 13;
  else if (orderReqd <= 14)
    order = 14;
  else if (orderReqd <= 30)
    order = 30;
  else
    SANS_DEVELOPER_EXCEPTION("Quadrature order=%d is not available", orderReqd);

  // set weights, coordinates for reference triangle (x in [0, 1], y in [0, 1-x])
  // Note: weights sum to 1 instead of 0.5 (area of reference triangle).
  //The reference element area is accounted for by a member variable named areaRef in element topology class.
  // TODO: rationale behind this approach???

  // degree 1 polynomial; 1 point
  if (order == 1)
  {
    nQuad = 1;

    order_  = order;
    nQuad_  = nQuad;
    weight_ = new Real[nQuad];
    xcoord_ = new Real[nQuad];
    ycoord_ = new Real[nQuad];

    weight_[0] = 1;
    xcoord_[0] = 1./3.;
    ycoord_[0] = 1./3.;
  }

  // degree 2 polynomial; 3 points
  // ref: Strang-Fix formula #1
  else if (order == 2)
  {
    nQuad = 3;

    order_  = order;
    nQuad_  = nQuad;
    weight_ = new Real[nQuad];
    xcoord_ = new Real[nQuad];
    ycoord_ = new Real[nQuad];

    weight_[0] = 1./3.;
    xcoord_[0] = 1./6.;
    ycoord_[0] = 2./3.;

    weight_[1] = 1./3.;
    xcoord_[1] = 2./3.;
    ycoord_[1] = 1./6.;

    weight_[2] = 1./3.;
    xcoord_[2] = 1./6.;
    ycoord_[2] = 1./6.;
  }

#if 0
  // degree 3 polynomial; 4 points
  // ref: Strang-Fix formula #3
  else if (order == 3)
  {
    nQuad = 4;

    order_  = order;
    nQuad_  = nQuad;
    weight_ = new Real[nQuad];
    xcoord_ = new Real[nQuad];
    ycoord_ = new Real[nQuad];

    weight_[0] = -18./32.;     // NOTE: negative weight!!
    xcoord_[0] = 1./3.;
    ycoord_[0] = 1./3.;

    weight_[1] = 50./96.;
    xcoord_[1] = 0.2;
    ycoord_[1] = 0.2;

    weight_[2] = 50./96.;
    xcoord_[2] = 0.6;
    ycoord_[2] = 0.2;

    weight_[3] = 50./96.;
    xcoord_[3] = 0.2;
    ycoord_[3] = 0.6;
  }
#else
  // degree 3 polynomial; 6 points
  // ref: Strang-Fix formula #4; quadrature_26sep13.nb
  else if (order == 3)
  {
    nQuad = 6;

    order_  = order;
    nQuad_  = nQuad;
    weight_ = new Real[nQuad];
    xcoord_ = new Real[nQuad];
    ycoord_ = new Real[nQuad];

    Real w = 2./12.;
    Real a, b, c;
    a = 0.65902762237409221517838077125540;   // (1/3)(1 + cos((1/3) atan(3/4)))
    b = 0.23193336855303057249678456117469;
    c = 0.10903900907287721232483466756991;

    weight_[0] = w;
    xcoord_[0] = a;
    ycoord_[0] = b;

    weight_[1] = w;
    xcoord_[1] = b;
    ycoord_[1] = c;

    weight_[2] = w;
    xcoord_[2] = c;
    ycoord_[2] = a;

    weight_[3] = w;
    xcoord_[3] = a;
    ycoord_[3] = c;

    weight_[4] = w;
    xcoord_[4] = b;
    ycoord_[4] = a;

    weight_[5] = w;
    xcoord_[5] = c;
    ycoord_[5] = b;
  }
#endif

  // degree 4 polynomial; 6 points
  // ref: PXQuadTriangle (PXE_QuadRule3), quadrature_15oct14.nb
  else if (order == 4)
  {
    order_ = 4;
    nQuad_ = 6;
    weight_ = new Real[nQuad_];
    xcoord_ = new Real[nQuad_];
    ycoord_ = new Real[nQuad_];

    Real x, w;

    x = 4./9. - sqrt(10)/18. + sqrt(950 - 220*sqrt(10))/90.;
    w = 2.*(1./12. + 3*sqrt(95 - 22*sqrt(10))/496. - sqrt(950 - 220*sqrt(10))/7440.);
    xcoord_[0] = x;
    ycoord_[0] = x;
    weight_[0] = w;

    xcoord_[1] = x;
    ycoord_[1] = 1 - 2*x;
    weight_[1] = w;

    xcoord_[2] = 1 - 2*x;
    ycoord_[2] = x;
    weight_[2] = w;

    x = 4./9. - sqrt(10)/18. - sqrt(950 - 220*sqrt(10))/90.;
    w = 2.*(1./12. - 3*sqrt(95 - 22*sqrt(10))/496. + sqrt(950 - 220*sqrt(10))/7440.);
    xcoord_[3] = x;
    ycoord_[3] = x;
    weight_[3] = w;

    xcoord_[4] = x;
    ycoord_[4] = 1 - 2*x;
    weight_[4] = w;

    xcoord_[5] = 1 - 2*x;
    ycoord_[5] = x;
    weight_[5] = w;
  }

  // degree 5 polynomial; 7 points
  // ref: PXQuadTriangle (PXE_QuadRule4), quadrature_15oct14.nb
  else if (order == 5)
  {
    order_ = 5;
    nQuad_ = 7;
    weight_ = new Real[nQuad_];
    xcoord_ = new Real[nQuad_];
    ycoord_ = new Real[nQuad_];

    Real x, w;

    x = 1./3.;
    w = 2.*9./80.;
    xcoord_[0] = x;
    ycoord_[0] = x;
    weight_[0] = w;

    x = 2./7. + sqrt(15)/21.;
    w = 2.*(31./480. + sqrt(15)/2400.);
    xcoord_[1] = x;
    ycoord_[1] = x;
    weight_[1] = w;

    xcoord_[2] = x;
    ycoord_[2] = 1 - 2*x;
    weight_[2] = w;

    xcoord_[3] = 1 - 2*x;
    ycoord_[3] = x;
    weight_[3] = w;

    x = 2./7. - sqrt(15)/21.;
    w = 2.*(31./480. - sqrt(15)/2400.);
    xcoord_[4] = x;
    ycoord_[4] = x;
    weight_[4] = w;

    xcoord_[5] = x;
    ycoord_[5] = 1 - 2*x;
    weight_[5] = w;

    xcoord_[6] = 1 - 2*x;
    ycoord_[6] = x;
    weight_[6] = w;
  }

  // degree 8 polynomial; 16 points
  // ref: PXQuadTriangle (PXE_QuadRule7), quadrature_15oct14.nb
  else if (order == 8)
  {
    order_ = 8;
    nQuad_ = 16;
    weight_ = new Real[nQuad_];
    xcoord_ = new Real[nQuad_];
    ycoord_ = new Real[nQuad_];

    Real x, y, w;

    x = 1./3.;
    w = 2.*0.072157803838893584125545555244532;    // 2*729*(13764*sqrt(30) - 24517)/(513947840)
    xcoord_[0] = x;
    ycoord_[0] = x;
    weight_[0] = w;

    x = 0.17056930775176020662229350149146;
    w = 2.*0.051608685267359125140895775146065;
    xcoord_[1] = x;
    ycoord_[1] = x;
    weight_[1] = w;

    xcoord_[2] = x;
    ycoord_[2] = 1 - 2*x;
    weight_[2] = w;

    xcoord_[3] = 1 - 2*x;
    ycoord_[3] = x;
    weight_[3] = w;

    x = 0.45929258829272315602881551449417;
    w = 2.*0.047545817133642312396948052194292;
    xcoord_[4] = x;
    ycoord_[4] = x;
    weight_[4] = w;

    xcoord_[5] = x;
    ycoord_[5] = 1 - 2*x;
    weight_[5] = w;

    xcoord_[6] = 1 - 2*x;
    ycoord_[6] = x;
    weight_[6] = w;

    x = 0.050547228317030975458423550596599;
    w = 2.*0.016229248811599040155462964170890;
    xcoord_[7] = x;
    ycoord_[7] = x;
    weight_[7] = w;

    xcoord_[8] = x;
    ycoord_[8] = 1 - 2*x;
    weight_[8] = w;

    xcoord_[9] = 1 - 2*x;
    ycoord_[9] = x;
    weight_[9] = w;

    x = 0.72849239295540428124100037917606;
    y = 0.26311282963463811342178578628464;
    w = 2.*0.0136151570872174971324223450369545;
    xcoord_[10] = x;
    ycoord_[10] = y;
    weight_[10] = w;

    xcoord_[11] = y;
    ycoord_[11] = x;
    weight_[11] = w;

    xcoord_[12] = x;
    ycoord_[12] = 1 - x - y;
    weight_[12] = w;

    xcoord_[13] = 1 - x - y;
    ycoord_[13] = x;
    weight_[13] = w;

    xcoord_[14] = y;
    ycoord_[14] = 1 - x - y;
    weight_[14] = w;

    xcoord_[15] = 1 - x - y;
    ycoord_[15] = y;
    weight_[15] = w;
  }

  else if (order == 10)
  {
    // degree 10 polynomial; 25 points
    // ref: PXQuadTriangle (PXE_QuadRule9)
    order_ = 10;
    nQuad_ = 25;
    weight_ = new Real[nQuad_];
    xcoord_ = new Real[nQuad_];
    ycoord_ = new Real[nQuad_];

    xcoord_[  0] =   3.3333333333333331e-01;  ycoord_[  0] =   3.3333333333333331e-01;
    xcoord_[  1] =   4.2727317884677551e-01;  ycoord_[  1] =   1.4545364230644897e-01;
    xcoord_[  2] =   1.4545364230644897e-01;  ycoord_[  2] =   4.2727317884677551e-01;
    xcoord_[  3] =   4.2727317884677551e-01;  ycoord_[  3] =   4.2727317884677551e-01;
    xcoord_[  4] =   1.8309922244867502e-01;  ycoord_[  4] =   6.3380155510264991e-01;
    xcoord_[  5] =   6.3380155510264991e-01;  ycoord_[  5] =   1.8309922244867502e-01;
    xcoord_[  6] =   1.8309922244867502e-01;  ycoord_[  6] =   1.8309922244867502e-01;
    xcoord_[  7] =   4.9043401970113060e-01;  ycoord_[  7] =   1.9131960597738806e-02;
    xcoord_[  8] =   1.9131960597738806e-02;  ycoord_[  8] =   4.9043401970113060e-01;
    xcoord_[  9] =   4.9043401970113060e-01;  ycoord_[  9] =   4.9043401970113060e-01;
    xcoord_[ 10] =   1.2572445551580533e-02;  ycoord_[ 10] =   9.7485510889683891e-01;
    xcoord_[ 11] =   9.7485510889683891e-01;  ycoord_[ 11] =   1.2572445551580533e-02;
    xcoord_[ 12] =   1.2572445551580533e-02;  ycoord_[ 12] =   1.2572445551580533e-02;
    xcoord_[ 13] =   3.0804600168524771e-01;  ycoord_[ 13] =   3.7685330394686190e-02;
    xcoord_[ 14] =   3.7685330394686190e-02;  ycoord_[ 14] =   3.0804600168524771e-01;
    xcoord_[ 15] =   6.5426866792006610e-01;  ycoord_[ 15] =   3.7685330394686190e-02;
    xcoord_[ 16] =   3.7685330394686190e-02;  ycoord_[ 16] =   6.5426866792006610e-01;
    xcoord_[ 17] =   6.5426866792006610e-01;  ycoord_[ 17] =   3.0804600168524771e-01;
    xcoord_[ 18] =   3.0804600168524771e-01;  ycoord_[ 18] =   6.5426866792006610e-01;
    xcoord_[ 19] =   3.3371833739304788e-02;  ycoord_[ 19] =   8.4382358919213596e-01;
    xcoord_[ 20] =   8.4382358919213596e-01;  ycoord_[ 20] =   3.3371833739304788e-02;
    xcoord_[ 21] =   1.2280457706855927e-01;  ycoord_[ 21] =   8.4382358919213596e-01;
    xcoord_[ 22] =   8.4382358919213596e-01;  ycoord_[ 22] =   1.2280457706855927e-01;
    xcoord_[ 23] =   1.2280457706855927e-01;  ycoord_[ 23] =   3.3371833739304788e-02;
    xcoord_[ 24] =   3.3371833739304788e-02;  ycoord_[ 24] =   1.2280457706855927e-01;

    weight_[  0] =   4.0468714398811438e-02*2.0;
    weight_[  1] =   3.8649294001481559e-02*2.0;
    weight_[  2] =   3.8649294001481559e-02*2.0;
    weight_[  3] =   3.8649294001481559e-02*2.0;
    weight_[  4] =   3.9228819306185866e-02*2.0;
    weight_[  5] =   3.9228819306185866e-02*2.0;
    weight_[  6] =   3.9228819306185866e-02*2.0;
    weight_[  7] =   8.7345839979647427e-03*2.0;
    weight_[  8] =   8.7345839979647427e-03*2.0;
    weight_[  9] =   8.7345839979647427e-03*2.0;
    weight_[ 10] =   2.1461870924164139e-03*2.0;
    weight_[ 11] =   2.1461870924164139e-03*2.0;
    weight_[ 12] =   2.1461870924164139e-03*2.0;
    weight_[ 13] =   1.8734429105233822e-02*2.0;
    weight_[ 14] =   1.8734429105233822e-02*2.0;
    weight_[ 15] =   1.8734429105233822e-02*2.0;
    weight_[ 16] =   1.8734429105233822e-02*2.0;
    weight_[ 17] =   1.8734429105233822e-02*2.0;
    weight_[ 18] =   1.8734429105233822e-02*2.0;
    weight_[ 19] =   1.3474676295939981e-02*2.0;
    weight_[ 20] =   1.3474676295939981e-02*2.0;
    weight_[ 21] =   1.3474676295939981e-02*2.0;
    weight_[ 22] =   1.3474676295939981e-02*2.0;
    weight_[ 23] =   1.3474676295939981e-02*2.0;
    weight_[ 24] =   1.3474676295939981e-02*2.0;

  }
  else if (order == 11)
  {
    // degree 11 polynomial; 28 points
    // ref: PXQuadTriangle (PXE_QuadRule10)
    order_ = 11;
    nQuad_ = 28;
    weight_ = new Real[nQuad_];
    xcoord_ = new Real[nQuad_];
    ycoord_ = new Real[nQuad_];

    xcoord_[  0] =   3.3333333333333331e-01;  ycoord_[  0] =   3.3333333333333331e-01;
    xcoord_[  1] =   3.0938355245430784e-02;  ycoord_[  1] =   9.3812328950913848e-01;
    xcoord_[  2] =   9.3812328950913848e-01;  ycoord_[  2] =   3.0938355245430784e-02;
    xcoord_[  3] =   3.0938355245430784e-02;  ycoord_[  3] =   3.0938355245430784e-02;
    xcoord_[  4] =   4.3649818113412886e-01;  ycoord_[  4] =   1.2700363773174228e-01;
    xcoord_[  5] =   1.2700363773174228e-01;  ycoord_[  5] =   4.3649818113412886e-01;
    xcoord_[  6] =   4.3649818113412886e-01;  ycoord_[  6] =   4.3649818113412886e-01;
    xcoord_[  7] =   4.9898476370259326e-01;  ycoord_[  7] =   2.0304725948134816e-03;
    xcoord_[  8] =   2.0304725948134816e-03;  ycoord_[  8] =   4.9898476370259326e-01;
    xcoord_[  9] =   4.9898476370259326e-01;  ycoord_[  9] =   4.9898476370259326e-01;
    xcoord_[ 10] =   2.1468819795859434e-01;  ycoord_[ 10] =   5.7062360408281132e-01;
    xcoord_[ 11] =   5.7062360408281132e-01;  ycoord_[ 11] =   2.1468819795859434e-01;
    xcoord_[ 12] =   2.1468819795859434e-01;  ycoord_[ 12] =   2.1468819795859434e-01;
    xcoord_[ 13] =   1.1368310404211339e-01;  ycoord_[ 13] =   7.7263379191577319e-01;
    xcoord_[ 14] =   7.7263379191577319e-01;  ycoord_[ 14] =   1.1368310404211339e-01;
    xcoord_[ 15] =   1.1368310404211339e-01;  ycoord_[ 15] =   1.1368310404211339e-01;
    xcoord_[ 16] =   1.5974230459185018e-01;  ycoord_[ 16] =   1.4638929243286969e-02;
    xcoord_[ 17] =   1.4638929243286969e-02;  ycoord_[ 17] =   1.5974230459185018e-01;
    xcoord_[ 18] =   8.2561876616486285e-01;  ycoord_[ 18] =   1.4638929243286969e-02;
    xcoord_[ 19] =   1.4638929243286969e-02;  ycoord_[ 19] =   8.2561876616486285e-01;
    xcoord_[ 20] =   8.2561876616486285e-01;  ycoord_[ 20] =   1.5974230459185018e-01;
    xcoord_[ 21] =   1.5974230459185018e-01;  ycoord_[ 21] =   8.2561876616486285e-01;
    xcoord_[ 22] =   3.1178371570959901e-01;  ycoord_[ 22] =   4.7743974155535718e-02;
    xcoord_[ 23] =   4.7743974155535718e-02;  ycoord_[ 23] =   3.1178371570959901e-01;
    xcoord_[ 24] =   6.4047231013486527e-01;  ycoord_[ 24] =   4.7743974155535718e-02;
    xcoord_[ 25] =   4.7743974155535718e-02;  ycoord_[ 25] =   6.4047231013486527e-01;
    xcoord_[ 26] =   6.4047231013486527e-01;  ycoord_[ 26] =   3.1178371570959901e-01;
    xcoord_[ 27] =   3.1178371570959901e-01;  ycoord_[ 27] =   6.4047231013486527e-01;

    weight_[  0] =   4.0588980148433582e-02*2.0;
    weight_[  1] =   6.1620217534547467e-03*2.0;
    weight_[  2] =   6.1620217534547467e-03*2.0;
    weight_[  3] =   6.1620217534547467e-03*2.0;
    weight_[  4] =   3.1414004872205054e-02*2.0;
    weight_[  5] =   3.1414004872205054e-02*2.0;
    weight_[  6] =   3.1414004872205054e-02*2.0;
    weight_[  7] =   6.1101895246822649e-03*2.0;
    weight_[  8] =   6.1101895246822649e-03*2.0;
    weight_[  9] =   6.1101895246822649e-03*2.0;
    weight_[ 10] =   3.3850674476405752e-02*2.0;
    weight_[ 11] =   3.3850674476405752e-02*2.0;
    weight_[ 12] =   3.3850674476405752e-02*2.0;
    weight_[ 13] =   2.0109846814425847e-02*2.0;
    weight_[ 14] =   2.0109846814425847e-02*2.0;
    weight_[ 15] =   2.0109846814425847e-02*2.0;
    weight_[ 16] =   7.3811363588580509e-03*2.0;
    weight_[ 17] =   7.3811363588580509e-03*2.0;
    weight_[ 18] =   7.3811363588580509e-03*2.0;
    weight_[ 19] =   7.3811363588580509e-03*2.0;
    weight_[ 20] =   7.3811363588580509e-03*2.0;
    weight_[ 21] =   7.3811363588580509e-03*2.0;
    weight_[ 22] =   2.0363998229149520e-02*2.0;
    weight_[ 23] =   2.0363998229149520e-02*2.0;
    weight_[ 24] =   2.0363998229149520e-02*2.0;
    weight_[ 25] =   2.0363998229149520e-02*2.0;
    weight_[ 26] =   2.0363998229149520e-02*2.0;
    weight_[ 27] =   2.0363998229149520e-02*2.0;
  }
  else if (order == 12)
  {
    // degree 12 polynomial; 33 points
    // ref: PXQuadTriangle (PXE_QuadRule11)
    order_ = 12;
    nQuad_ = 33;
    weight_ = new Real[nQuad_];
    xcoord_ = new Real[nQuad_];
    ycoord_ = new Real[nQuad_];

    xcoord_[  0] =   2.1317350453210371e-02;  ycoord_[  0] =   9.5736529909357926e-01;
    xcoord_[  1] =   9.5736529909357926e-01;  ycoord_[  1] =   2.1317350453210371e-02;
    xcoord_[  2] =   2.1317350453210371e-02;  ycoord_[  2] =   2.1317350453210371e-02;
    xcoord_[  3] =   2.7121038501211592e-01;  ycoord_[  3] =   4.5757922997576816e-01;
    xcoord_[  4] =   4.5757922997576816e-01;  ycoord_[  4] =   2.7121038501211592e-01;
    xcoord_[  5] =   2.7121038501211592e-01;  ycoord_[  5] =   2.7121038501211592e-01;
    xcoord_[  6] =   1.2757614554158592e-01;  ycoord_[  6] =   7.4484770891682817e-01;
    xcoord_[  7] =   7.4484770891682817e-01;  ycoord_[  7] =   1.2757614554158592e-01;
    xcoord_[  8] =   1.2757614554158592e-01;  ycoord_[  8] =   1.2757614554158592e-01;
    xcoord_[  9] =   4.3972439229446025e-01;  ycoord_[  9] =   1.2055121541107949e-01;
    xcoord_[ 10] =   1.2055121541107949e-01;  ycoord_[ 10] =   4.3972439229446025e-01;
    xcoord_[ 11] =   4.3972439229446025e-01;  ycoord_[ 11] =   4.3972439229446025e-01;
    xcoord_[ 12] =   4.8821738977380486e-01;  ycoord_[ 12] =   2.3565220452390290e-02;
    xcoord_[ 13] =   2.3565220452390290e-02;  ycoord_[ 13] =   4.8821738977380486e-01;
    xcoord_[ 14] =   4.8821738977380486e-01;  ycoord_[ 14] =   4.8821738977380486e-01;
    xcoord_[ 15] =   2.8132558098993954e-01;  ycoord_[ 15] =   2.2838332222257063e-02;
    xcoord_[ 16] =   2.2838332222257063e-02;  ycoord_[ 16] =   2.8132558098993954e-01;
    xcoord_[ 17] =   6.9583608678780340e-01;  ycoord_[ 17] =   2.2838332222257063e-02;
    xcoord_[ 18] =   2.2838332222257063e-02;  ycoord_[ 18] =   6.9583608678780340e-01;
    xcoord_[ 19] =   6.9583608678780340e-01;  ycoord_[ 19] =   2.8132558098993954e-01;
    xcoord_[ 20] =   2.8132558098993954e-01;  ycoord_[ 20] =   6.9583608678780340e-01;
    xcoord_[ 21] =   1.1625191590759715e-01;  ycoord_[ 21] =   2.5734050548330167e-02;
    xcoord_[ 22] =   2.5734050548330167e-02;  ycoord_[ 22] =   1.1625191590759715e-01;
    xcoord_[ 23] =   8.5801403354407269e-01;  ycoord_[ 23] =   2.5734050548330167e-02;
    xcoord_[ 24] =   2.5734050548330167e-02;  ycoord_[ 24] =   8.5801403354407269e-01;
    xcoord_[ 25] =   8.5801403354407269e-01;  ycoord_[ 25] =   1.1625191590759715e-01;
    xcoord_[ 26] =   1.1625191590759715e-01;  ycoord_[ 26] =   8.5801403354407269e-01;
    xcoord_[ 27] =   2.7571326968551418e-01;  ycoord_[ 27] =   1.1534349453469800e-01;
    xcoord_[ 28] =   1.1534349453469800e-01;  ycoord_[ 28] =   2.7571326968551418e-01;
    xcoord_[ 29] =   6.0894323577978782e-01;  ycoord_[ 29] =   1.1534349453469800e-01;
    xcoord_[ 30] =   1.1534349453469800e-01;  ycoord_[ 30] =   6.0894323577978782e-01;
    xcoord_[ 31] =   6.0894323577978782e-01;  ycoord_[ 31] =   2.7571326968551418e-01;
    xcoord_[ 32] =   2.7571326968551418e-01;  ycoord_[ 32] =   6.0894323577978782e-01;

    weight_[  0] =   3.0831305257795088e-03*2.0;
    weight_[  1] =   3.0831305257795088e-03*2.0;
    weight_[  2] =   3.0831305257795088e-03*2.0;
    weight_[  3] =   3.1429112108942552e-02*2.0;
    weight_[  4] =   3.1429112108942552e-02*2.0;
    weight_[  5] =   3.1429112108942552e-02*2.0;
    weight_[  6] =   1.7398056465354472e-02*2.0;
    weight_[  7] =   1.7398056465354472e-02*2.0;
    weight_[  8] =   1.7398056465354472e-02*2.0;
    weight_[  9] =   2.1846272269019203e-02*2.0;
    weight_[ 10] =   2.1846272269019203e-02*2.0;
    weight_[ 11] =   2.1846272269019203e-02*2.0;
    weight_[ 12] =   1.2865533220227668e-02*2.0;
    weight_[ 13] =   1.2865533220227668e-02*2.0;
    weight_[ 14] =   1.2865533220227668e-02*2.0;
    weight_[ 15] =   1.1178386601151722e-02*2.0;
    weight_[ 16] =   1.1178386601151722e-02*2.0;
    weight_[ 17] =   1.1178386601151722e-02*2.0;
    weight_[ 18] =   1.1178386601151722e-02*2.0;
    weight_[ 19] =   1.1178386601151722e-02*2.0;
    weight_[ 20] =   1.1178386601151722e-02*2.0;
    weight_[ 21] =   8.6581155543294461e-03*2.0;
    weight_[ 22] =   8.6581155543294461e-03*2.0;
    weight_[ 23] =   8.6581155543294461e-03*2.0;
    weight_[ 24] =   8.6581155543294461e-03*2.0;
    weight_[ 25] =   8.6581155543294461e-03*2.0;
    weight_[ 26] =   8.6581155543294461e-03*2.0;
    weight_[ 27] =   2.0185778883190463e-02*2.0;
    weight_[ 28] =   2.0185778883190463e-02*2.0;
    weight_[ 29] =   2.0185778883190463e-02*2.0;
    weight_[ 30] =   2.0185778883190463e-02*2.0;
    weight_[ 31] =   2.0185778883190463e-02*2.0;
    weight_[ 32] =   2.0185778883190463e-02*2.0;

  }
  else if (order == 13)
  {
    // degree 13 polynomial; 37 points
    // ref: PXQuadTriangle (PXE_QuadRule12)
    order_ = 13;
    nQuad_ = 37;
    weight_ = new Real[nQuad_];
    xcoord_ = new Real[nQuad_];
    ycoord_ = new Real[nQuad_];

    xcoord_[  0] =   3.3333333333333331e-01;  ycoord_[  0] =   3.3333333333333331e-01;
    xcoord_[  1] =   4.2694141425980042e-01;  ycoord_[  1] =   1.4611717148039916e-01;
    xcoord_[  2] =   1.4611717148039916e-01;  ycoord_[  2] =   4.2694141425980042e-01;
    xcoord_[  3] =   4.2694141425980042e-01;  ycoord_[  3] =   4.2694141425980042e-01;
    xcoord_[  4] =   2.2137228629183289e-01;  ycoord_[  4] =   5.5725542741633416e-01;
    xcoord_[  5] =   5.5725542741633416e-01;  ycoord_[  5] =   2.2137228629183289e-01;
    xcoord_[  6] =   2.2137228629183289e-01;  ycoord_[  6] =   2.2137228629183289e-01;
    xcoord_[  7] =   2.1509681108843184e-02;  ycoord_[  7] =   9.5698063778231368e-01;
    xcoord_[  8] =   9.5698063778231368e-01;  ycoord_[  8] =   2.1509681108843184e-02;
    xcoord_[  9] =   2.1509681108843184e-02;  ycoord_[  9] =   2.1509681108843184e-02;
    xcoord_[ 10] =   4.8907694645253935e-01;  ycoord_[ 10] =   2.1846107094921297e-02;
    xcoord_[ 11] =   2.1846107094921297e-02;  ycoord_[ 11] =   4.8907694645253935e-01;
    xcoord_[ 12] =   4.8907694645253935e-01;  ycoord_[ 12] =   4.8907694645253935e-01;
    xcoord_[ 13] =   3.0844176089211778e-01;  ycoord_[ 13] =   6.8012243554206597e-02;
    xcoord_[ 14] =   6.8012243554206597e-02;  ycoord_[ 14] =   3.0844176089211778e-01;
    xcoord_[ 15] =   6.2354599555367562e-01;  ycoord_[ 15] =   6.8012243554206597e-02;
    xcoord_[ 16] =   6.8012243554206597e-02;  ycoord_[ 16] =   6.2354599555367562e-01;
    xcoord_[ 17] =   6.2354599555367562e-01;  ycoord_[ 17] =   3.0844176089211778e-01;
    xcoord_[ 18] =   3.0844176089211778e-01;  ycoord_[ 18] =   6.2354599555367562e-01;
    xcoord_[ 19] =   1.1092204280346339e-01;  ycoord_[ 19] =   2.4370186901093840e-02;
    xcoord_[ 20] =   2.4370186901093840e-02;  ycoord_[ 20] =   1.1092204280346339e-01;
    xcoord_[ 21] =   8.6470777029544277e-01;  ycoord_[ 21] =   2.4370186901093840e-02;
    xcoord_[ 22] =   2.4370186901093840e-02;  ycoord_[ 22] =   8.6470777029544277e-01;
    xcoord_[ 23] =   8.6470777029544277e-01;  ycoord_[ 23] =   1.1092204280346339e-01;
    xcoord_[ 24] =   1.1092204280346339e-01;  ycoord_[ 24] =   8.6470777029544277e-01;
    xcoord_[ 25] =   1.6359740106785048e-01;  ycoord_[ 25] =   8.7895483032197297e-02;
    xcoord_[ 26] =   8.7895483032197297e-02;  ycoord_[ 26] =   1.6359740106785048e-01;
    xcoord_[ 27] =   7.4850711589995222e-01;  ycoord_[ 27] =   8.7895483032197297e-02;
    xcoord_[ 28] =   8.7895483032197297e-02;  ycoord_[ 28] =   7.4850711589995222e-01;
    xcoord_[ 29] =   7.4850711589995222e-01;  ycoord_[ 29] =   1.6359740106785048e-01;
    xcoord_[ 30] =   1.6359740106785048e-01;  ycoord_[ 30] =   7.4850711589995222e-01;
    xcoord_[ 31] =   2.7251581777342965e-01;  ycoord_[ 31] =   5.1263891023823338e-03;
    xcoord_[ 32] =   5.1263891023823338e-03;  ycoord_[ 32] =   2.7251581777342965e-01;
    xcoord_[ 33] =   7.2235779312418802e-01;  ycoord_[ 33] =   5.1263891023823338e-03;
    xcoord_[ 34] =   5.1263891023823338e-03;  ycoord_[ 34] =   7.2235779312418802e-01;
    xcoord_[ 35] =   7.2235779312418802e-01;  ycoord_[ 35] =   2.7251581777342965e-01;
    xcoord_[ 36] =   2.7251581777342965e-01;  ycoord_[ 36] =   7.2235779312418802e-01;

    weight_[  0] =   3.3980018293415820e-02*2.0;
    weight_[  1] =   2.7800983765226665e-02*2.0;
    weight_[  2] =   2.7800983765226665e-02*2.0;
    weight_[  3] =   2.7800983765226665e-02*2.0;
    weight_[  4] =   2.9139242559599991e-02*2.0;
    weight_[  5] =   2.9139242559599991e-02*2.0;
    weight_[  6] =   2.9139242559599991e-02*2.0;
    weight_[  7] =   3.0261685517695858e-03*2.0;
    weight_[  8] =   3.0261685517695858e-03*2.0;
    weight_[  9] =   3.0261685517695858e-03*2.0;
    weight_[ 10] =   1.1997200964447365e-02*2.0;
    weight_[ 11] =   1.1997200964447365e-02*2.0;
    weight_[ 12] =   1.1997200964447365e-02*2.0;
    weight_[ 13] =   1.7320638070424187e-02*2.0;
    weight_[ 14] =   1.7320638070424187e-02*2.0;
    weight_[ 15] =   1.7320638070424187e-02*2.0;
    weight_[ 16] =   1.7320638070424187e-02*2.0;
    weight_[ 17] =   1.7320638070424187e-02*2.0;
    weight_[ 18] =   1.7320638070424187e-02*2.0;
    weight_[ 19] =   7.4827005525828338e-03*2.0;
    weight_[ 20] =   7.4827005525828338e-03*2.0;
    weight_[ 21] =   7.4827005525828338e-03*2.0;
    weight_[ 22] =   7.4827005525828338e-03*2.0;
    weight_[ 23] =   7.4827005525828338e-03*2.0;
    weight_[ 24] =   7.4827005525828338e-03*2.0;
    weight_[ 25] =   1.2089519905796910e-02*2.0;
    weight_[ 26] =   1.2089519905796910e-02*2.0;
    weight_[ 27] =   1.2089519905796910e-02*2.0;
    weight_[ 28] =   1.2089519905796910e-02*2.0;
    weight_[ 29] =   1.2089519905796910e-02*2.0;
    weight_[ 30] =   1.2089519905796910e-02*2.0;
    weight_[ 31] =   4.7953405017716316e-03*2.0;
    weight_[ 32] =   4.7953405017716316e-03*2.0;
    weight_[ 33] =   4.7953405017716316e-03*2.0;
    weight_[ 34] =   4.7953405017716316e-03*2.0;
    weight_[ 35] =   4.7953405017716316e-03*2.0;
    weight_[ 36] =   4.7953405017716316e-03*2.0;

  }
  else if (order == 14)
  {
    // degree 14 polynomial; 46 points
    // ref: PXQuadTriangle (PXE_QuadRule13)
    order_ = 14;
    nQuad_ = 46;
    weight_ = new Real[nQuad_];
    xcoord_ = new Real[nQuad_];
    ycoord_ = new Real[nQuad_];

    xcoord_[  0] =   3.3333333333333331e-01;  ycoord_[  0] =   3.3333333333333331e-01;
    xcoord_[  1] =   9.9797608064584320e-03;  ycoord_[  1] =   9.8004047838708308e-01;
    xcoord_[  2] =   9.8004047838708308e-01;  ycoord_[  2] =   9.9797608064584320e-03;
    xcoord_[  3] =   9.9797608064584320e-03;  ycoord_[  3] =   9.9797608064584320e-03;
    xcoord_[  4] =   4.7997789352118841e-01;  ycoord_[  4] =   4.0044212957623171e-02;
    xcoord_[  5] =   4.0044212957623171e-02;  ycoord_[  5] =   4.7997789352118841e-01;
    xcoord_[  6] =   4.7997789352118841e-01;  ycoord_[  6] =   4.7997789352118841e-01;
    xcoord_[  7] =   1.5381195917696691e-01;  ycoord_[  7] =   6.9237608164606623e-01;
    xcoord_[  8] =   6.9237608164606623e-01;  ycoord_[  8] =   1.5381195917696691e-01;
    xcoord_[  9] =   1.5381195917696691e-01;  ycoord_[  9] =   1.5381195917696691e-01;
    xcoord_[ 10] =   7.4023477116987813e-02;  ycoord_[ 10] =   8.5195304576602437e-01;
    xcoord_[ 11] =   8.5195304576602437e-01;  ycoord_[ 11] =   7.4023477116987813e-02;
    xcoord_[ 12] =   7.4023477116987813e-02;  ycoord_[ 12] =   7.4023477116987813e-02;
    xcoord_[ 13] =   1.3035468250332999e-01;  ycoord_[ 13] =   7.3929063499334002e-01;
    xcoord_[ 14] =   7.3929063499334002e-01;  ycoord_[ 14] =   1.3035468250332999e-01;
    xcoord_[ 15] =   1.3035468250332999e-01;  ycoord_[ 15] =   1.3035468250332999e-01;
    xcoord_[ 16] =   2.3061722602665313e-01;  ycoord_[ 16] =   5.3876554794669373e-01;
    xcoord_[ 17] =   5.3876554794669373e-01;  ycoord_[ 17] =   2.3061722602665313e-01;
    xcoord_[ 18] =   2.3061722602665313e-01;  ycoord_[ 18] =   2.3061722602665313e-01;
    xcoord_[ 19] =   4.2233208341914780e-01;  ycoord_[ 19] =   1.5533583316170441e-01;
    xcoord_[ 20] =   1.5533583316170441e-01;  ycoord_[ 20] =   4.2233208341914780e-01;
    xcoord_[ 21] =   4.2233208341914780e-01;  ycoord_[ 21] =   4.2233208341914780e-01;
    xcoord_[ 22] =   1.9061636003190091e-01;  ycoord_[ 22] =   2.3146254033438118e-02;
    xcoord_[ 23] =   2.3146254033438118e-02;  ycoord_[ 23] =   1.9061636003190091e-01;
    xcoord_[ 24] =   7.8623738593466097e-01;  ycoord_[ 24] =   2.3146254033438118e-02;
    xcoord_[ 25] =   2.3146254033438118e-02;  ycoord_[ 25] =   7.8623738593466097e-01;
    xcoord_[ 26] =   7.8623738593466097e-01;  ycoord_[ 26] =   1.9061636003190091e-01;
    xcoord_[ 27] =   1.9061636003190091e-01;  ycoord_[ 27] =   7.8623738593466097e-01;
    xcoord_[ 28] =   3.6232313774354713e-01;  ycoord_[ 28] =   7.1247185958454584e-03;
    xcoord_[ 29] =   7.1247185958454584e-03;  ycoord_[ 29] =   3.6232313774354713e-01;
    xcoord_[ 30] =   6.3055214366060741e-01;  ycoord_[ 30] =   7.1247185958454584e-03;
    xcoord_[ 31] =   7.1247185958454584e-03;  ycoord_[ 31] =   6.3055214366060741e-01;
    xcoord_[ 32] =   6.3055214366060741e-01;  ycoord_[ 32] =   3.6232313774354713e-01;
    xcoord_[ 33] =   3.6232313774354713e-01;  ycoord_[ 33] =   6.3055214366060741e-01;
    xcoord_[ 34] =   2.9077120588366739e-01;  ycoord_[ 34] =   8.2651464260026286e-02;
    xcoord_[ 35] =   8.2651464260026286e-02;  ycoord_[ 35] =   2.9077120588366739e-01;
    xcoord_[ 36] =   6.2657732985630632e-01;  ycoord_[ 36] =   8.2651464260026286e-02;
    xcoord_[ 37] =   8.2651464260026286e-02;  ycoord_[ 37] =   6.2657732985630632e-01;
    xcoord_[ 38] =   6.2657732985630632e-01;  ycoord_[ 38] =   2.9077120588366739e-01;
    xcoord_[ 39] =   2.9077120588366739e-01;  ycoord_[ 39] =   6.2657732985630632e-01;
    xcoord_[ 40] =   7.1165710877750768e-02;  ycoord_[ 40] =   1.4624304192623769e-02;
    xcoord_[ 41] =   1.4624304192623769e-02;  ycoord_[ 41] =   7.1165710877750768e-02;
    xcoord_[ 42] =   9.1420998492962546e-01;  ycoord_[ 42] =   1.4624304192623769e-02;
    xcoord_[ 43] =   1.4624304192623769e-02;  ycoord_[ 43] =   9.1420998492962546e-01;
    xcoord_[ 44] =   9.1420998492962546e-01;  ycoord_[ 44] =   7.1165710877750768e-02;
    xcoord_[ 45] =   7.1165710877750768e-02;  ycoord_[ 45] =   9.1420998492962546e-01;

    weight_[  0] =   2.9298142613014298e-02*2.0;
    weight_[  1] =   8.6757561486263378e-04*2.0;
    weight_[  2] =   8.6757561486263378e-04*2.0;
    weight_[  3] =   8.6757561486263378e-04*2.0;
    weight_[  4] =   1.3081891279307261e-02*2.0;
    weight_[  5] =   1.3081891279307261e-02*2.0;
    weight_[  6] =   1.3081891279307261e-02*2.0;
    weight_[  7] =   1.9598646212009145e-03*2.0;
    weight_[  8] =   1.9598646212009145e-03*2.0;
    weight_[  9] =   1.9598646212009145e-03*2.0;
    weight_[ 10] =   6.1236798784704335e-03*2.0;
    weight_[ 11] =   6.1236798784704335e-03*2.0;
    weight_[ 12] =   6.1236798784704335e-03*2.0;
    weight_[ 13] =   1.4099814251628980e-02*2.0;
    weight_[ 14] =   1.4099814251628980e-02*2.0;
    weight_[ 15] =   1.4099814251628980e-02*2.0;
    weight_[ 16] =   2.5443543592979744e-02*2.0;
    weight_[ 17] =   2.5443543592979744e-02*2.0;
    weight_[ 18] =   2.5443543592979744e-02*2.0;
    weight_[ 19] =   2.5226719950801800e-02*2.0;
    weight_[ 20] =   2.5226719950801800e-02*2.0;
    weight_[ 21] =   2.5226719950801800e-02*2.0;
    weight_[ 22] =   8.5318221061167262e-03*2.0;
    weight_[ 23] =   8.5318221061167262e-03*2.0;
    weight_[ 24] =   8.5318221061167262e-03*2.0;
    weight_[ 25] =   8.5318221061167262e-03*2.0;
    weight_[ 26] =   8.5318221061167262e-03*2.0;
    weight_[ 27] =   8.5318221061167262e-03*2.0;
    weight_[ 28] =   4.8417332127533002e-03*2.0;
    weight_[ 29] =   4.8417332127533002e-03*2.0;
    weight_[ 30] =   4.8417332127533002e-03*2.0;
    weight_[ 31] =   4.8417332127533002e-03*2.0;
    weight_[ 32] =   4.8417332127533002e-03*2.0;
    weight_[ 33] =   4.8417332127533002e-03*2.0;
    weight_[ 34] =   1.8192877964242501e-02*2.0;
    weight_[ 35] =   1.8192877964242501e-02*2.0;
    weight_[ 36] =   1.8192877964242501e-02*2.0;
    weight_[ 37] =   1.8192877964242501e-02*2.0;
    weight_[ 38] =   1.8192877964242501e-02*2.0;
    weight_[ 39] =   1.8192877964242501e-02*2.0;
    weight_[ 40] =   3.4823316867592063e-03*2.0;
    weight_[ 41] =   3.4823316867592063e-03*2.0;
    weight_[ 42] =   3.4823316867592063e-03*2.0;
    weight_[ 43] =   3.4823316867592063e-03*2.0;
    weight_[ 44] =   3.4823316867592063e-03*2.0;
    weight_[ 45] =   3.4823316867592063e-03*2.0;
  }
  // degree 30 polynomial; 175 points
  // ref: PXQuadTriangle (PXE_QuadRule14; Wandzura, Xiao)
  else if (order == 30)
  {
    order_ = 30;
    nQuad_ = 175;
    weight_ = new Real[nQuad_];
    xcoord_ = new Real[nQuad_];
    ycoord_ = new Real[nQuad_];

    xcoord_[  0] =  0.333333333333330;  ycoord_[  0] =  0.333333333333330;  weight_[  0] =  0.015579960202899;
    xcoord_[  1] =  0.007330116432770;  ycoord_[  1] =  0.496334941783620;  weight_[  1] =  0.003177233700534;
    xcoord_[  2] =  0.496334941783620;  ycoord_[  2] =  0.496334941783620;  weight_[  2] =  0.003177233700534;
    xcoord_[  3] =  0.496334941783620;  ycoord_[  3] =  0.007330116432770;  weight_[  3] =  0.003177233700534;
    xcoord_[  4] =  0.082995675802960;  ycoord_[  4] =  0.458502162098520;  weight_[  4] =  0.010483426635731;
    xcoord_[  5] =  0.458502162098520;  ycoord_[  5] =  0.458502162098520;  weight_[  5] =  0.010483426635731;
    xcoord_[  6] =  0.458502162098520;  ycoord_[  6] =  0.082995675802960;  weight_[  6] =  0.010483426635731;
    xcoord_[  7] =  0.150980956125410;  ycoord_[  7] =  0.424509521937290;  weight_[  7] =  0.013209459577744;
    xcoord_[  8] =  0.424509521937290;  ycoord_[  8] =  0.424509521937290;  weight_[  8] =  0.013209459577744;
    xcoord_[  9] =  0.424509521937290;  ycoord_[  9] =  0.150980956125410;  weight_[  9] =  0.013209459577744;
    xcoord_[ 10] =  0.235905859892170;  ycoord_[ 10] =  0.382047070053920;  weight_[ 10] =  0.014975006966271;
    xcoord_[ 11] =  0.382047070053920;  ycoord_[ 11] =  0.382047070053920;  weight_[ 11] =  0.014975006966271;
    xcoord_[ 12] =  0.382047070053920;  ycoord_[ 12] =  0.235905859892170;  weight_[ 12] =  0.014975006966271;
    xcoord_[ 13] =  0.438024308407850;  ycoord_[ 13] =  0.280987845796080;  weight_[ 13] =  0.014987904443384;
    xcoord_[ 14] =  0.280987845796080;  ycoord_[ 14] =  0.280987845796080;  weight_[ 14] =  0.014987904443384;
    xcoord_[ 15] =  0.280987845796080;  ycoord_[ 15] =  0.438024308407850;  weight_[ 15] =  0.014987904443384;
    xcoord_[ 16] =  0.545302048291930;  ycoord_[ 16] =  0.227348975854030;  weight_[ 16] =  0.013338864741022;
    xcoord_[ 17] =  0.227348975854030;  ycoord_[ 17] =  0.227348975854030;  weight_[ 17] =  0.013338864741022;
    xcoord_[ 18] =  0.227348975854030;  ycoord_[ 18] =  0.545302048291930;  weight_[ 18] =  0.013338864741022;
    xcoord_[ 19] =  0.650881776982540;  ycoord_[ 19] =  0.174559111508730;  weight_[ 19] =  0.010889171113902;
    xcoord_[ 20] =  0.174559111508730;  ycoord_[ 20] =  0.174559111508730;  weight_[ 20] =  0.010889171113902;
    xcoord_[ 21] =  0.174559111508730;  ycoord_[ 21] =  0.650881776982540;  weight_[ 21] =  0.010889171113902;
    xcoord_[ 22] =  0.753483145597130;  ycoord_[ 22] =  0.123258427201440;  weight_[ 22] =  0.008189440660893;
    xcoord_[ 23] =  0.123258427201440;  ycoord_[ 23] =  0.123258427201440;  weight_[ 23] =  0.008189440660893;
    xcoord_[ 24] =  0.123258427201440;  ycoord_[ 24] =  0.753483145597130;  weight_[ 24] =  0.008189440660893;
    xcoord_[ 25] =  0.839831542215610;  ycoord_[ 25] =  0.080084228892200;  weight_[ 25] =  0.005575387588608;
    xcoord_[ 26] =  0.080084228892200;  ycoord_[ 26] =  0.080084228892200;  weight_[ 26] =  0.005575387588608;
    xcoord_[ 27] =  0.080084228892200;  ycoord_[ 27] =  0.839831542215610;  weight_[ 27] =  0.005575387588608;
    xcoord_[ 28] =  0.904451065184200;  ycoord_[ 28] =  0.047774467407900;  weight_[ 28] =  0.003191216473412;
    xcoord_[ 29] =  0.047774467407900;  ycoord_[ 29] =  0.047774467407900;  weight_[ 29] =  0.003191216473412;
    xcoord_[ 30] =  0.047774467407900;  ycoord_[ 30] =  0.904451065184200;  weight_[ 30] =  0.003191216473412;
    xcoord_[ 31] =  0.956558970639720;  ycoord_[ 31] =  0.021720514680140;  weight_[ 31] =  0.001296715144327;
    xcoord_[ 32] =  0.021720514680140;  ycoord_[ 32] =  0.021720514680140;  weight_[ 32] =  0.001296715144327;
    xcoord_[ 33] =  0.021720514680140;  ycoord_[ 33] =  0.956558970639720;  weight_[ 33] =  0.001296715144327;
    xcoord_[ 34] =  0.990470644769130;  ycoord_[ 34] =  0.004764677615440;  weight_[ 34] =  0.000298262826135;
    xcoord_[ 35] =  0.004764677615440;  ycoord_[ 35] =  0.004764677615440;  weight_[ 35] =  0.000298262826135;
    xcoord_[ 36] =  0.004764677615440;  ycoord_[ 36] =  0.990470644769130;  weight_[ 36] =  0.000298262826135;
    xcoord_[ 37] =  0.000925371193350;  ycoord_[ 37] =  0.415295270913310;  weight_[ 37] =  0.000998905685079;
    xcoord_[ 38] =  0.415295270913310;  ycoord_[ 38] =  0.583779357893340;  weight_[ 38] =  0.000998905685079;
    xcoord_[ 39] =  0.583779357893340;  ycoord_[ 39] =  0.000925371193350;  weight_[ 39] =  0.000998905685079;
    xcoord_[ 40] =  0.415295270913310;  ycoord_[ 40] =  0.000925371193350;  weight_[ 40] =  0.000998905685079;
    xcoord_[ 41] =  0.583779357893340;  ycoord_[ 41] =  0.415295270913310;  weight_[ 41] =  0.000998905685079;
    xcoord_[ 42] =  0.000925371193350;  ycoord_[ 42] =  0.583779357893340;  weight_[ 42] =  0.000998905685079;
    xcoord_[ 43] =  0.001385925855560;  ycoord_[ 43] =  0.061189909785350;  weight_[ 43] =  0.000462850849173;
    xcoord_[ 44] =  0.061189909785350;  ycoord_[ 44] =  0.937424164359090;  weight_[ 44] =  0.000462850849173;
    xcoord_[ 45] =  0.937424164359090;  ycoord_[ 45] =  0.001385925855560;  weight_[ 45] =  0.000462850849173;
    xcoord_[ 46] =  0.061189909785350;  ycoord_[ 46] =  0.001385925855560;  weight_[ 46] =  0.000462850849173;
    xcoord_[ 47] =  0.937424164359090;  ycoord_[ 47] =  0.061189909785350;  weight_[ 47] =  0.000462850849173;
    xcoord_[ 48] =  0.001385925855560;  ycoord_[ 48] =  0.937424164359090;  weight_[ 48] =  0.000462850849173;
    xcoord_[ 49] =  0.003682415455910;  ycoord_[ 49] =  0.164908690136910;  weight_[ 49] =  0.001234451336382;
    xcoord_[ 50] =  0.164908690136910;  ycoord_[ 50] =  0.831408894407180;  weight_[ 50] =  0.001234451336382;
    xcoord_[ 51] =  0.831408894407180;  ycoord_[ 51] =  0.003682415455910;  weight_[ 51] =  0.001234451336382;
    xcoord_[ 52] =  0.164908690136910;  ycoord_[ 52] =  0.003682415455910;  weight_[ 52] =  0.001234451336382;
    xcoord_[ 53] =  0.831408894407180;  ycoord_[ 53] =  0.164908690136910;  weight_[ 53] =  0.001234451336382;
    xcoord_[ 54] =  0.003682415455910;  ycoord_[ 54] =  0.831408894407180;  weight_[ 54] =  0.001234451336382;
    xcoord_[ 55] =  0.003903223424160;  ycoord_[ 55] =  0.025035062232000;  weight_[ 55] =  0.000570719852243;
    xcoord_[ 56] =  0.025035062232000;  ycoord_[ 56] =  0.971061714343840;  weight_[ 56] =  0.000570719852243;
    xcoord_[ 57] =  0.971061714343840;  ycoord_[ 57] =  0.003903223424160;  weight_[ 57] =  0.000570719852243;
    xcoord_[ 58] =  0.025035062232000;  ycoord_[ 58] =  0.003903223424160;  weight_[ 58] =  0.000570719852243;
    xcoord_[ 59] =  0.971061714343840;  ycoord_[ 59] =  0.025035062232000;  weight_[ 59] =  0.000570719852243;
    xcoord_[ 60] =  0.003903223424160;  ycoord_[ 60] =  0.971061714343840;  weight_[ 60] =  0.000570719852243;
    xcoord_[ 61] =  0.003233248155010;  ycoord_[ 61] =  0.306064465151100;  weight_[ 61] =  0.001126946125878;
    xcoord_[ 62] =  0.306064465151100;  ycoord_[ 62] =  0.690702286693890;  weight_[ 62] =  0.001126946125878;
    xcoord_[ 63] =  0.690702286693890;  ycoord_[ 63] =  0.003233248155010;  weight_[ 63] =  0.001126946125878;
    xcoord_[ 64] =  0.306064465151100;  ycoord_[ 64] =  0.003233248155010;  weight_[ 64] =  0.001126946125878;
    xcoord_[ 65] =  0.690702286693890;  ycoord_[ 65] =  0.306064465151100;  weight_[ 65] =  0.001126946125878;
    xcoord_[ 66] =  0.003233248155010;  ycoord_[ 66] =  0.690702286693890;  weight_[ 66] =  0.001126946125878;
    xcoord_[ 67] =  0.006467432112240;  ycoord_[ 67] =  0.107073283730220;  weight_[ 67] =  0.001747866949407;
    xcoord_[ 68] =  0.107073283730220;  ycoord_[ 68] =  0.886459284157540;  weight_[ 68] =  0.001747866949407;
    xcoord_[ 69] =  0.886459284157540;  ycoord_[ 69] =  0.006467432112240;  weight_[ 69] =  0.001747866949407;
    xcoord_[ 70] =  0.107073283730220;  ycoord_[ 70] =  0.006467432112240;  weight_[ 70] =  0.001747866949407;
    xcoord_[ 71] =  0.886459284157540;  ycoord_[ 71] =  0.107073283730220;  weight_[ 71] =  0.001747866949407;
    xcoord_[ 72] =  0.006467432112240;  ycoord_[ 72] =  0.886459284157540;  weight_[ 72] =  0.001747866949407;
    xcoord_[ 73] =  0.003247475491330;  ycoord_[ 73] =  0.229957549345580;  weight_[ 73] =  0.001182818815032;
    xcoord_[ 74] =  0.229957549345580;  ycoord_[ 74] =  0.766794975163080;  weight_[ 74] =  0.001182818815032;
    xcoord_[ 75] =  0.766794975163080;  ycoord_[ 75] =  0.003247475491330;  weight_[ 75] =  0.001182818815032;
    xcoord_[ 76] =  0.229957549345580;  ycoord_[ 76] =  0.003247475491330;  weight_[ 76] =  0.001182818815032;
    xcoord_[ 77] =  0.766794975163080;  ycoord_[ 77] =  0.229957549345580;  weight_[ 77] =  0.001182818815032;
    xcoord_[ 78] =  0.003247475491330;  ycoord_[ 78] =  0.766794975163080;  weight_[ 78] =  0.001182818815032;
    xcoord_[ 79] =  0.008675090806750;  ycoord_[ 79] =  0.337036633305780;  weight_[ 79] =  0.001990839294675;
    xcoord_[ 80] =  0.337036633305780;  ycoord_[ 80] =  0.654288275887460;  weight_[ 80] =  0.001990839294675;
    xcoord_[ 81] =  0.654288275887460;  ycoord_[ 81] =  0.008675090806750;  weight_[ 81] =  0.001990839294675;
    xcoord_[ 82] =  0.337036633305780;  ycoord_[ 82] =  0.008675090806750;  weight_[ 82] =  0.001990839294675;
    xcoord_[ 83] =  0.654288275887460;  ycoord_[ 83] =  0.337036633305780;  weight_[ 83] =  0.001990839294675;
    xcoord_[ 84] =  0.008675090806750;  ycoord_[ 84] =  0.654288275887460;  weight_[ 84] =  0.001990839294675;
    xcoord_[ 85] =  0.015597026467310;  ycoord_[ 85] =  0.056256576182060;  weight_[ 85] =  0.001900412795036;
    xcoord_[ 86] =  0.056256576182060;  ycoord_[ 86] =  0.928146397350630;  weight_[ 86] =  0.001900412795036;
    xcoord_[ 87] =  0.928146397350630;  ycoord_[ 87] =  0.015597026467310;  weight_[ 87] =  0.001900412795036;
    xcoord_[ 88] =  0.056256576182060;  ycoord_[ 88] =  0.015597026467310;  weight_[ 88] =  0.001900412795036;
    xcoord_[ 89] =  0.928146397350630;  ycoord_[ 89] =  0.056256576182060;  weight_[ 89] =  0.001900412795036;
    xcoord_[ 90] =  0.015597026467310;  ycoord_[ 90] =  0.928146397350630;  weight_[ 90] =  0.001900412795036;
    xcoord_[ 91] =  0.017976721253690;  ycoord_[ 91] =  0.402451375212400;  weight_[ 91] =  0.004498365808817;
    xcoord_[ 92] =  0.402451375212400;  ycoord_[ 92] =  0.579571903533910;  weight_[ 92] =  0.004498365808817;
    xcoord_[ 93] =  0.579571903533910;  ycoord_[ 93] =  0.017976721253690;  weight_[ 93] =  0.004498365808817;
    xcoord_[ 94] =  0.402451375212400;  ycoord_[ 94] =  0.017976721253690;  weight_[ 94] =  0.004498365808817;
    xcoord_[ 95] =  0.579571903533910;  ycoord_[ 95] =  0.402451375212400;  weight_[ 95] =  0.004498365808817;
    xcoord_[ 96] =  0.017976721253690;  ycoord_[ 96] =  0.579571903533910;  weight_[ 96] =  0.004498365808817;
    xcoord_[ 97] =  0.017124245353890;  ycoord_[ 97] =  0.243654702010830;  weight_[ 97] =  0.003478719460275;
    xcoord_[ 98] =  0.243654702010830;  ycoord_[ 98] =  0.739221052635280;  weight_[ 98] =  0.003478719460275;
    xcoord_[ 99] =  0.739221052635280;  ycoord_[ 99] =  0.017124245353890;  weight_[ 99] =  0.003478719460275;
    xcoord_[100] =  0.243654702010830;  ycoord_[100] =  0.017124245353890;  weight_[100] =  0.003478719460275;
    xcoord_[101] =  0.739221052635280;  ycoord_[101] =  0.243654702010830;  weight_[101] =  0.003478719460275;
    xcoord_[102] =  0.017124245353890;  ycoord_[102] =  0.739221052635280;  weight_[102] =  0.003478719460275;
    xcoord_[103] =  0.022883405346580;  ycoord_[103] =  0.165389585614530;  weight_[103] =  0.004102399036724;
    xcoord_[104] =  0.165389585614530;  ycoord_[104] =  0.811727009038880;  weight_[104] =  0.004102399036724;
    xcoord_[105] =  0.811727009038880;  ycoord_[105] =  0.022883405346580;  weight_[105] =  0.004102399036724;
    xcoord_[106] =  0.165389585614530;  ycoord_[106] =  0.022883405346580;  weight_[106] =  0.004102399036724;
    xcoord_[107] =  0.811727009038880;  ycoord_[107] =  0.165389585614530;  weight_[107] =  0.004102399036724;
    xcoord_[108] =  0.022883405346580;  ycoord_[108] =  0.811727009038880;  weight_[108] =  0.004102399036724;
    xcoord_[109] =  0.032737597287770;  ycoord_[109] =  0.099301874495850;  weight_[109] =  0.004021761549744;
    xcoord_[110] =  0.099301874495850;  ycoord_[110] =  0.867960528216390;  weight_[110] =  0.004021761549744;
    xcoord_[111] =  0.867960528216390;  ycoord_[111] =  0.032737597287770;  weight_[111] =  0.004021761549744;
    xcoord_[112] =  0.099301874495850;  ycoord_[112] =  0.032737597287770;  weight_[112] =  0.004021761549744;
    xcoord_[113] =  0.867960528216390;  ycoord_[113] =  0.099301874495850;  weight_[113] =  0.004021761549744;
    xcoord_[114] =  0.032737597287770;  ycoord_[114] =  0.867960528216390;  weight_[114] =  0.004021761549744;
    xcoord_[115] =  0.033821012342340;  ycoord_[115] =  0.308478333069050;  weight_[115] =  0.006033164660795;
    xcoord_[116] =  0.308478333069050;  ycoord_[116] =  0.657700654588600;  weight_[116] =  0.006033164660795;
    xcoord_[117] =  0.657700654588600;  ycoord_[117] =  0.033821012342340;  weight_[117] =  0.006033164660795;
    xcoord_[118] =  0.308478333069050;  ycoord_[118] =  0.033821012342340;  weight_[118] =  0.006033164660795;
    xcoord_[119] =  0.657700654588600;  ycoord_[119] =  0.308478333069050;  weight_[119] =  0.006033164660795;
    xcoord_[120] =  0.033821012342340;  ycoord_[120] =  0.657700654588600;  weight_[120] =  0.006033164660795;
    xcoord_[121] =  0.035547614460020;  ycoord_[121] =  0.460668318592110;  weight_[121] =  0.003946290302130;
    xcoord_[122] =  0.460668318592110;  ycoord_[122] =  0.503784066947870;  weight_[122] =  0.003946290302130;
    xcoord_[123] =  0.503784066947870;  ycoord_[123] =  0.035547614460020;  weight_[123] =  0.003946290302130;
    xcoord_[124] =  0.460668318592110;  ycoord_[124] =  0.035547614460020;  weight_[124] =  0.003946290302130;
    xcoord_[125] =  0.503784066947870;  ycoord_[125] =  0.460668318592110;  weight_[125] =  0.003946290302130;
    xcoord_[126] =  0.035547614460020;  ycoord_[126] =  0.503784066947870;  weight_[126] =  0.003946290302130;
    xcoord_[127] =  0.050539790306870;  ycoord_[127] =  0.218815299453930;  weight_[127] =  0.006644044537680;
    xcoord_[128] =  0.218815299453930;  ycoord_[128] =  0.730644910239200;  weight_[128] =  0.006644044537680;
    xcoord_[129] =  0.730644910239200;  ycoord_[129] =  0.050539790306870;  weight_[129] =  0.006644044537680;
    xcoord_[130] =  0.218815299453930;  ycoord_[130] =  0.050539790306870;  weight_[130] =  0.006644044537680;
    xcoord_[131] =  0.730644910239200;  ycoord_[131] =  0.218815299453930;  weight_[131] =  0.006644044537680;
    xcoord_[132] =  0.050539790306870;  ycoord_[132] =  0.730644910239200;  weight_[132] =  0.006644044537680;
    xcoord_[133] =  0.057014714915730;  ycoord_[133] =  0.379209551560270;  weight_[133] =  0.008254305856078;
    xcoord_[134] =  0.379209551560270;  ycoord_[134] =  0.563775733523990;  weight_[134] =  0.008254305856078;
    xcoord_[135] =  0.563775733523990;  ycoord_[135] =  0.057014714915730;  weight_[135] =  0.008254305856078;
    xcoord_[136] =  0.379209551560270;  ycoord_[136] =  0.057014714915730;  weight_[136] =  0.008254305856078;
    xcoord_[137] =  0.563775733523990;  ycoord_[137] =  0.379209551560270;  weight_[137] =  0.008254305856078;
    xcoord_[138] =  0.057014714915730;  ycoord_[138] =  0.563775733523990;  weight_[138] =  0.008254305856078;
    xcoord_[139] =  0.064152806421200;  ycoord_[139] =  0.142960819418190;  weight_[139] =  0.006496056633406;
    xcoord_[140] =  0.142960819418190;  ycoord_[140] =  0.792886374160610;  weight_[140] =  0.006496056633406;
    xcoord_[141] =  0.792886374160610;  ycoord_[141] =  0.064152806421200;  weight_[141] =  0.006496056633406;
    xcoord_[142] =  0.142960819418190;  ycoord_[142] =  0.064152806421200;  weight_[142] =  0.006496056633406;
    xcoord_[143] =  0.792886374160610;  ycoord_[143] =  0.142960819418190;  weight_[143] =  0.006496056633406;
    xcoord_[144] =  0.064152806421200;  ycoord_[144] =  0.792886374160610;  weight_[144] =  0.006496056633406;
    xcoord_[145] =  0.080501148287630;  ycoord_[145] =  0.283731282105920;  weight_[145] =  0.009252778144147;
    xcoord_[146] =  0.283731282105920;  ycoord_[146] =  0.635767569606450;  weight_[146] =  0.009252778144147;
    xcoord_[147] =  0.635767569606450;  ycoord_[147] =  0.080501148287630;  weight_[147] =  0.009252778144147;
    xcoord_[148] =  0.283731282105920;  ycoord_[148] =  0.080501148287630;  weight_[148] =  0.009252778144147;
    xcoord_[149] =  0.635767569606450;  ycoord_[149] =  0.283731282105920;  weight_[149] =  0.009252778144147;
    xcoord_[150] =  0.080501148287630;  ycoord_[150] =  0.635767569606450;  weight_[150] =  0.009252778144147;
    xcoord_[151] =  0.104367068134530;  ycoord_[151] =  0.196737441004440;  weight_[151] =  0.009164920726294;
    xcoord_[152] =  0.196737441004440;  ycoord_[152] =  0.698895490861030;  weight_[152] =  0.009164920726294;
    xcoord_[153] =  0.698895490861030;  ycoord_[153] =  0.104367068134530;  weight_[153] =  0.009164920726294;
    xcoord_[154] =  0.196737441004440;  ycoord_[154] =  0.104367068134530;  weight_[154] =  0.009164920726294;
    xcoord_[155] =  0.698895490861030;  ycoord_[155] =  0.196737441004440;  weight_[155] =  0.009164920726294;
    xcoord_[156] =  0.104367068134530;  ycoord_[156] =  0.698895490861030;  weight_[156] =  0.009164920726294;
    xcoord_[157] =  0.113844894428750;  ycoord_[157] =  0.355889141211660;  weight_[157] =  0.011569524628098;
    xcoord_[158] =  0.355889141211660;  ycoord_[158] =  0.530265964359590;  weight_[158] =  0.011569524628098;
    xcoord_[159] =  0.530265964359590;  ycoord_[159] =  0.113844894428750;  weight_[159] =  0.011569524628098;
    xcoord_[160] =  0.355889141211660;  ycoord_[160] =  0.113844894428750;  weight_[160] =  0.011569524628098;
    xcoord_[161] =  0.530265964359590;  ycoord_[161] =  0.355889141211660;  weight_[161] =  0.011569524628098;
    xcoord_[162] =  0.113844894428750;  ycoord_[162] =  0.530265964359590;  weight_[162] =  0.011569524628098;
    xcoord_[163] =  0.145363487715520;  ycoord_[163] =  0.259818685351910;  weight_[163] =  0.011761116467609;
    xcoord_[164] =  0.259818685351910;  ycoord_[164] =  0.594817826932560;  weight_[164] =  0.011761116467609;
    xcoord_[165] =  0.594817826932560;  ycoord_[165] =  0.145363487715520;  weight_[165] =  0.011761116467609;
    xcoord_[166] =  0.259818685351910;  ycoord_[166] =  0.145363487715520;  weight_[166] =  0.011761116467609;
    xcoord_[167] =  0.594817826932560;  ycoord_[167] =  0.259818685351910;  weight_[167] =  0.011761116467609;
    xcoord_[168] =  0.145363487715520;  ycoord_[168] =  0.594817826932560;  weight_[168] =  0.011761116467609;
    xcoord_[169] =  0.189945652821980;  ycoord_[169] =  0.321923181231300;  weight_[169] =  0.013824702182165;
    xcoord_[170] =  0.321923181231300;  ycoord_[170] =  0.488131165946720;  weight_[170] =  0.013824702182165;
    xcoord_[171] =  0.488131165946720;  ycoord_[171] =  0.189945652821980;  weight_[171] =  0.013824702182165;
    xcoord_[172] =  0.321923181231300;  ycoord_[172] =  0.189945652821980;  weight_[172] =  0.013824702182165;
    xcoord_[173] =  0.488131165946720;  ycoord_[173] =  0.321923181231300;  weight_[173] =  0.013824702182165;
    xcoord_[174] =  0.189945652821980;  ycoord_[174] =  0.488131165946720;  weight_[174] =  0.013824702182165;
  }
}


template <>
QuadratureArea<Triangle>::~QuadratureArea()
{
  delete [] ycoord_;
  delete [] xcoord_;
  delete [] weight_;
  nQuad_ = 0;
}


template <>
void
QuadratureArea<Triangle>::weight( int n, Real& wght ) const
{
  SANS_ASSERT ((n >= 0) && (n < nQuad_));
  wght = weight_[n];
}


template <>
void
QuadratureArea<Triangle>::coordinates( int n, Real xy[] ) const
{
  SANS_ASSERT ((n >= 0) && (n < nQuad_));
  xy[0] = xcoord_[n];
  xy[1] = ycoord_[n];
}


template <>
void
QuadratureArea<Triangle>::coordinates( int n, Real& x, Real& y ) const
{
  SANS_ASSERT ((n >= 0) && (n < nQuad_));
  x = xcoord_[n];
  y = ycoord_[n];
}

template <>
SANS::DLA::VectorS<2,Real>
QuadratureArea<Triangle>::coordinates( int n ) const
{
  SANS_ASSERT ((n >= 0) && (n < nQuad_));
  return {xcoord_[n], ycoord_[n]};
}

template <>
SANS::QuadraturePoint<TopoD2>
QuadratureArea<Triangle>::coordinates_cache( int n ) const
{
  return SANS::QuadraturePoint<TopoD2>(SANS::QuadratureRule::eGauss, n, orderidx_, {xcoord_[n], ycoord_[n]});
}


template <>
void
QuadratureArea<Triangle>::dump( int indentSize, std::ostream& out ) const
{
  std::string indent(indentSize, ' ');
  out << indent << "QuadratureArea<Triangle>: total points = " << nQuad_ << std::endl;
  if (weight_ != 0)
  {
    out << indent << "  weights = ";
    for (int k = 0; k < nQuad_; k++)
      out << weight_[k] << " ";
    out << std::endl;
  }
  if (xcoord_ != 0)
  {
    out << indent << "  xcoords = ";
    for (int k = 0; k < nQuad_; k++)
      out << xcoord_[k] << " ";
    out << std::endl;
  }
  if (ycoord_ != 0)
  {
    out << indent << "  ycoords = ";
    for (int k = 0; k < nQuad_; k++)
      out << ycoord_[k] << " ";
    out << std::endl;
  }
}


// I/O
template <>
std::ostream&
operator<<( std::ostream& out, const QuadratureArea<Triangle>& quadrature )
{
  int nquad = quadrature.nQuadrature();
  Real w, xy[2];
  out << "QuadratureArea<Triangle> " << nquad << ": ";
  for (int k = 0; k < nquad; k++)
  {
    quadrature.weight(k, w);
    quadrature.coordinates(k, xy);
    out << "(" << w << ", " << xy[0] <<  ", " << xy[1] <<") ";
  }
  return out;
}

} // namespace SANS
