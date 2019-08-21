// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "BasisFunctionArea_Triangle_Lagrange.h"

#include "tools/SANSException.h"
#include "LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"

namespace SANS
{

#ifndef ONESIXTH
#define ONESIXTH 0.16666666666666666666666666666667
#endif

#ifndef ONETHIRD
#define ONETHIRD 0.3333333333333333333333333333333
#endif

// Node locations for (uniform) Lagrange basis

/*

Triangle Q1:            Triangle Q2:         Triangle Q3:          Triangle Q4:

t
^                                                                   2
|                                                                   | \
2                       2                    2                      6   5
|`\                     |`\                  | \                    |     \
|  `\                   |  `\                5   4                  7 (14)  4
|    `\                 4    `3              |     \                |         \
|      `\               |      `\            6  (9)  3              8 (12) (13) 3
|        `\             |        `\          |         \            |             \
0----------1 --> s      0-----5----1         0---7---8---1          0---9--10--11---1

*/

const std::vector<Real> BasisFunctionArea<Triangle,Lagrange,1>::coords_s_ = {0.0, 1.0, 0.0};
const std::vector<Real> BasisFunctionArea<Triangle,Lagrange,1>::coords_t_ = {0.0, 0.0, 1.0};

const std::vector<Real> BasisFunctionArea<Triangle,Lagrange,2>::coords_s_ = {0.0, 1.0, 0.0, 0.5, 0.0, 0.5};
const std::vector<Real> BasisFunctionArea<Triangle,Lagrange,2>::coords_t_ = {0.0, 0.0, 1.0, 0.5, 0.5, 0.0};

const std::vector<Real> BasisFunctionArea<Triangle,Lagrange,3>::coords_s_ = {0.0, 1.0, 0.0, 2./3., 1./3., 0.0, 0.0, 1./3., 2./3., 1./3.};
const std::vector<Real> BasisFunctionArea<Triangle,Lagrange,3>::coords_t_ = {0.0, 0.0, 1.0, 1./3., 2./3., 2./3., 1./3., 0.0, 0.0, 1./3.};

const std::vector<Real> BasisFunctionArea<Triangle,Lagrange,4>::coords_s_ = {0.0, 1.0, 0.0,
                                                                             0.75, 0.5, 0.25,
                                                                             0.0, 0.0, 0.0,
                                                                             0.25, 0.5, 0.75,
                                                                             0.25, 0.5, 0.25
                                                                            };
const std::vector<Real> BasisFunctionArea<Triangle,Lagrange,4>::coords_t_ = {0.0, 0.0, 1.0,
                                                                             0.25, 0.5, 0.75,
                                                                             0.75, 0.5, 0.25,
                                                                             0.0, 0.0, 0.0,
                                                                             0.25, 0.25, 0.5
                                                                            };

const std::vector<Real> BasisFunctionArea<Triangle,Lagrange,5>::coords_s_ = {0.0, 1.0, 0.0,
                                                                             0.8, 0.6, 0.4, 0.2,
                                                                             0.0, 0.0, 0.0, 0.0,
                                                                             0.2, 0.4, 0.6, 0.8,
                                                                             0.2, 0.4, 0.6, 0.2, 0.4, 0.2
                                                                            };
const std::vector<Real> BasisFunctionArea<Triangle,Lagrange,5>::coords_t_ = {0.0, 0.0, 1.0,
                                                                             0.2, 0.4, 0.6, 0.8,
                                                                             0.8, 0.6, 0.4, 0.2,
                                                                             0.0, 0.0, 0.0, 0.0,
                                                                             0.2, 0.2, 0.2, 0.4, 0.4, 0.6
                                                                            };
//const std::vector<Real> BasisFunctionArea<Triangle,Lagrange,6>::coords_s_ = {0.0, 1.0, 0.0, //vertices
//                                                                             5./6., 4./6., 3./6., 2./6., 1./6., // edge 0
//                                                                               0.0,   0.0,   0.0,   0.0,   0.0, // edge 1
//                                                                             1./6., 2./6., 3./6., 4./6., 5./6., // edge 2
//                                                                             1./6., 2./6., 3./6., 4./6., // row 1
//                                                                             1./6., 2./6., 3./6., // row 2
//                                                                             1./6., 2./6., // row 3
//                                                                             1./6. // row 4
//                                                                            };
//const std::vector<Real> BasisFunctionArea<Triangle,Lagrange,6>::coords_t_ = { 0.0, 0.0, 1.0, //vertices
//                                                                              1./6., 2./6., 3./6., 4./6., 5./6., // edge 0
//                                                                              5./6., 4./6., 3./6., 2./6., 1./6., // edge 1
//                                                                                0.0,   0.0,   0.0,   0.0,   0.0, // edge 2
//                                                                              1./6., 1./6., 1./6., 1./6., // row 1
//                                                                              2./6., 2./6., 2./6., // row 2
//                                                                              3./6., 3./6., // row 3
//                                                                              4./6. // row 4
//                                                                            };
//const std::vector<Real> BasisFunctionArea<Triangle,Lagrange,7>::coords_s_ = {0.0, 1.0, 0.0, //vertices
//                                                                             6./7., 5./7., 4./7., 3./7., 2./7., 1./7., // edge 0
//                                                                               0.0,   0.0,   0.0,   0.0,   0.0,   0.0, // edge 1
//                                                                             1./7., 2./7., 3./7., 4./7., 5./7., 6./7., // edge 2
//                                                                             1./7., 2./7., 3./7., 4./7., 5./7., // row 1
//                                                                             1./7., 2./7., 3./7., 4./7., // row 2
//                                                                             1./7., 2./7., 3./7., // row 3
//                                                                             1./7., 2./7., // row 4
//                                                                             1./7. // row 5
//                                                                            };
//const std::vector<Real> BasisFunctionArea<Triangle,Lagrange,7>::coords_t_ = { 0.0, 0.0, 1.0, //vertices
//                                                                              1./7., 2./7., 3./7., 4./7., 5./7., 6./7., // edge 0
//                                                                              6./7., 5./7., 4./7., 3./7., 2./7., 1./7., // edge 1
//                                                                                0.0,   0.0,   0.0,   0.0,   0.0, // edge 2
//                                                                              1./7., 1./7., 1./7., 1./7., 1./7., // row 1
//                                                                              2./7., 2./7., 2./7., 2./7., // row 2
//                                                                              3./7., 3./7., 3./7., // row 3
//                                                                              4./7., 4./7., // row 4
//                                                                              5./7. // row 5
//                                                                            };

void getLagrangeNodes_Triangle(const int order, std::vector<Real>& s, std::vector<Real>& t)
{
  switch (order)
  {
    case 1:
    {
      const BasisFunctionArea<Triangle,Lagrange,1>* basisP1 = BasisFunctionArea<Triangle,Lagrange,1>::self();
      basisP1->coordinates(s, t);
      break;
    }

    case 2:
    {
      const BasisFunctionArea<Triangle,Lagrange,2>* basisP2 = BasisFunctionArea<Triangle,Lagrange,2>::self();
      basisP2->coordinates(s, t);
      break;
    }

    case 3:
    {
      const BasisFunctionArea<Triangle,Lagrange,3>* basisP3 = BasisFunctionArea<Triangle,Lagrange,3>::self();
      basisP3->coordinates(s, t);
      break;
    }

    case 4:
    {
      const BasisFunctionArea<Triangle,Lagrange,4>* basisP4 = BasisFunctionArea<Triangle,Lagrange,4>::self();
      basisP4->coordinates(s, t);
      break;
    }

    case 5:
    {
      const BasisFunctionArea<Triangle,Lagrange,5>* basisP5 = BasisFunctionArea<Triangle,Lagrange,5>::self();
      basisP5->coordinates(s, t);
      break;
    }

//    case 6:
//    {
//      const BasisFunctionArea<Triangle,Lagrange,6>* basisP6 = BasisFunctionArea<Triangle,Lagrange,6>::self();
//      basisP6->coordinates(s, t);
//      break;
//    }
//
//    case 7:
//    {
//      const BasisFunctionArea<Triangle,Lagrange,7>* basisP7 = BasisFunctionArea<Triangle,Lagrange,7>::self();
//      basisP7->coordinates(s, t);
//      break;
//    }

    default:
      SANS_DEVELOPER_EXCEPTION( "getLagrangeNodes_Triangle - Unsupported order = %d.", order );
  }
}

//---------------------------------------------------------------------------//
const int LagrangeNodes<Triangle>::PMax;

//---------------------------------------------------------------------------//
void
LagrangeNodes<Triangle>::
get(const int order, std::vector<DLA::VectorS<TopoD2::D,Real>>& sRef)
{
  std::vector<Real> s, t;
  getLagrangeNodes_Triangle(order, s, t);
  sRef.resize(s.size());

  for (std::size_t i = 0; i < s.size(); i++)
  {
    sRef[i][0] = s[i];
    sRef[i][1] = t[i];
  }
}

//---------------------------------------------------------------------------//
void
BasisFunctionArea<Triangle,Lagrange,1>::evalBasis(
const Real& s, const Real& t, const Int3& sgn , Real phi[], int nphi) const
{
  SANS_ASSERT(nphi==3);

  // phi
  phi[0] =  -s-t+1.0;
  phi[1] =  s;
  phi[2] =  t;
}

void
BasisFunctionArea<Triangle,Lagrange,1>::evalBasisDerivative(
const Real& s, const Real& t, const Int3& sgn , Real phis[], Real phit[], int nphi) const
{
  SANS_ASSERT(nphi==3);

  // phis
  phis[0] =  -1.0;
  phis[1] =  1.0;
  phis[2] =  0.0;

  // phit
  phit[0] =  -1.0;
  phit[1] =  0.0;
  phit[2] =  1.0;
}

void
BasisFunctionArea<Triangle,Lagrange,1>::evalBasisHessianDerivative(
const Real& s, const Real& t, const Int3& sgn , Real phiss[], Real phist[], Real phitt[], int nphi) const
{
  SANS_ASSERT(nphi==3);

  // phiss
  phiss[0] =  0.0;
  phiss[1] =  0.0;
  phiss[2] =  0.0;

  // phist
  phist[0] =  0.0;
  phist[1] =  0.0;
  phist[2] =  0.0;

  // phitt
  phitt[0] =  0.0;
  phitt[1] =  0.0;
  phitt[2] =  0.0;
}

void
BasisFunctionArea<Triangle,Lagrange,1>::coordinates( std::vector<Real>& s, std::vector<Real>& t ) const
{
  s = coords_s_;
  t = coords_t_;
}

void
BasisFunctionArea<Triangle,Lagrange,2>::evalBasis(
const Real& s, const Real& t, const Int3& sgn , Real phi[], int nphi) const
{
  SANS_ASSERT(nphi==6);

  // phi
  phi[0] =  s*-3.0-t*3.0+s*t*4.0+(s*s)*2.0+(t*t)*2.0+1.0;
  phi[1] =  -s+(s*s)*2.0;
  phi[2] =  -t+(t*t)*2.0;
  phi[3] =  s*t*4.0;
  phi[4] =  t*(s+t-1.0)*-4.0;
  phi[5] =  -s*(s*4.0+t*4.0-4.0);
}

void
BasisFunctionArea<Triangle,Lagrange,2>::evalBasisDerivative(
const Real& s, const Real& t, const Int3& sgn , Real phis[], Real phit[], int nphi) const
{
  SANS_ASSERT(nphi==6);

  // phis
  phis[0] =  s*4.0+t*4.0-3.0;
  phis[1] =  s*4.0-1.0;
  phis[2] =  0.0;
  phis[3] =  t*4.0;
  phis[4] =  t*-4.0;
  phis[5] =  s*-8.0-t*4.0+4.0;

  // phit
  phit[0] =  s*4.0+t*4.0-3.0;
  phit[1] =  0.0;
  phit[2] =  t*4.0-1.0;
  phit[3] =  s*4.0;
  phit[4] =  s*-4.0-t*8.0+4.0;
  phit[5] =  s*-4.0;
}

void
BasisFunctionArea<Triangle,Lagrange,2>::evalBasisHessianDerivative(
const Real& s, const Real& t, const Int3& sgn , Real phiss[], Real phist[], Real phitt[], int nphi) const
{
  SANS_ASSERT(nphi==6);

  // phiss
  phiss[0] =  4.0;
  phiss[1] =  4.0;
  phiss[2] =  0.0;
  phiss[3] =  0.0;
  phiss[4] =  0.0;
  phiss[5] =  -8.0;

  // phist
  phist[0] =  4.0;
  phist[1] =  0.0;
  phist[2] =  0.0;
  phist[3] =  4.0;
  phist[4] =  -4.0;
  phist[5] =  -4.0;

  // phitt
  phitt[0] =  4.0;
  phitt[1] =  0.0;
  phitt[2] =  4.0;
  phitt[3] =  0.0;
  phitt[4] =  -8.0;
  phitt[5] =  0.0;
}

void
BasisFunctionArea<Triangle,Lagrange,2>::coordinates( std::vector<Real>& s, std::vector<Real>& t ) const
{
  s = coords_s_;
  t = coords_t_;
}

void
BasisFunctionArea<Triangle,Lagrange,3>::evalBasis(
const Real& s, const Real& t, const Int3& sgn , Real phi[], int nphi) const
{
  SANS_ASSERT(nphi==10);

  // phi
  phi[0] =  s*(-1.1E1/2.0)-t*(1.1E1/2.0)+s*t*1.8E1-s*(t*t)*(2.7E1/2.0)-(s*s)*t*(2.7E1/2.0)+(s*s)*9.0-
              (s*s*s)*(9.0/2.0)+(t*t)*9.0-(t*t*t)*(9.0/2.0)+1.0;

  phi[1] =  s-(s*s)*(9.0/2.0)+(s*s*s)*(9.0/2.0);

  phi[2] =  t-(t*t)*(9.0/2.0)+(t*t*t)*(9.0/2.0);

  phi[3] =  s*t*(-9.0/2.0)+(s*s)*t*(2.7E1/2.0);

  phi[4] =  s*t*(-9.0/2.0)+s*(t*t)*(2.7E1/2.0);

  phi[5] =  t*(-9.0/2.0)+s*t*(9.0/2.0)-s*(t*t)*(2.7E1/2.0)+(t*t)*1.8E1-(t*t*t)*(2.7E1/2.0);

  phi[6] =  t*9.0-s*t*(4.5E1/2.0)+s*(t*t)*2.7E1+(s*s)*t*(2.7E1/2.0)-(t*t)*(4.5E1/2.0)+(t*t*t)*(2.7E1/2.0);

  phi[7] =  s*9.0-s*t*(4.5E1/2.0)+s*(t*t)*(2.7E1/2.0)+(s*s)*t*2.7E1-(s*s)*(4.5E1/2.0)+(s*s*s)*(2.7E1/2.0);

  phi[8] =  s*(-9.0/2.0)+s*t*(9.0/2.0)-(s*s)*t*(2.7E1/2.0)+(s*s)*1.8E1-(s*s*s)*(2.7E1/2.0);

  phi[9] =  s*t*2.7E1-s*(t*t)*2.7E1-(s*s)*t*2.7E1;
}

void
BasisFunctionArea<Triangle,Lagrange,3>::evalBasisDerivative(
const Real& s, const Real& t, const Int3& sgn , Real phis[], Real phit[], int nphi) const
{
  SANS_ASSERT(nphi==10);

  // phis
  phis[0] =  s*1.8E1+t*1.8E1-s*t*2.7E1-(s*s)*(2.7E1/2.0)-(t*t)*(2.7E1/2.0)-1.1E1/2.0;

  phis[1] =  s*-9.0+(s*s)*(2.7E1/2.0)+1.0;

  phis[2] =  0.0;

  phis[3] =  t*(-9.0/2.0)+s*t*2.7E1;

  phis[4] =  t*(-9.0/2.0)+(t*t)*(2.7E1/2.0);

  phis[5] =  t*(9.0/2.0)-(t*t)*(2.7E1/2.0);

  phis[6] =  t*(-4.5E1/2.0)+s*t*2.7E1+(t*t)*2.7E1;

  phis[7] =  s*-4.5E1-t*(4.5E1/2.0)+s*t*5.4E1+(s*s)*(8.1E1/2.0)+(t*t)*(2.7E1/2.0)+9.0;

  phis[8] =  s*3.6E1+t*(9.0/2.0)-s*t*2.7E1-(s*s)*(8.1E1/2.0)-9.0/2.0;

  phis[9] =  t*2.7E1-s*t*5.4E1-(t*t)*2.7E1;


  // phit
  phit[0] =  s*1.8E1+t*1.8E1-s*t*2.7E1-(s*s)*(2.7E1/2.0)-(t*t)*(2.7E1/2.0)-1.1E1/2.0;

  phit[1] =  0.0;

  phit[2] =  t*-9.0+(t*t)*(2.7E1/2.0)+1.0;

  phit[3] =  s*(-9.0/2.0)+(s*s)*(2.7E1/2.0);

  phit[4] =  s*(-9.0/2.0)+s*t*2.7E1;

  phit[5] =  s*(9.0/2.0)+t*3.6E1-s*t*2.7E1-(t*t)*(8.1E1/2.0)-9.0/2.0;

  phit[6] =  s*(-4.5E1/2.0)-t*4.5E1+s*t*5.4E1+(s*s)*(2.7E1/2.0)+(t*t)*(8.1E1/2.0)+9.0;

  phit[7] =  s*(-4.5E1/2.0)+s*t*2.7E1+(s*s)*2.7E1;

  phit[8] =  s*(9.0/2.0)-(s*s)*(2.7E1/2.0);

  phit[9] =  s*2.7E1-s*t*5.4E1-(s*s)*2.7E1;
}

void
BasisFunctionArea<Triangle,Lagrange,3>::evalBasisHessianDerivative(
const Real& s, const Real& t, const Int3& sgn , Real phiss[], Real phist[], Real phitt[], int nphi) const
{
  SANS_ASSERT(nphi==10);

  // phiss
  phiss[0] =  s*-2.7E1-t*2.7E1+1.8E1;

  phiss[1] =  s*2.7E1-9.0;

  phiss[2] =  0.0;

  phiss[3] =  t*2.7E1;

  phiss[4] =  0.0;

  phiss[5] =  0.0;

  phiss[6] =  t*2.7E1;

  phiss[7] =  s*8.1E1+t*5.4E1-4.5E1;

  phiss[8] =  s*-8.1E1-t*2.7E1+3.6E1;

  phiss[9] =  t*-5.4E1;


  // phist
  phist[0] =  s*-2.7E1-t*2.7E1+1.8E1;

  phist[1] =  0.0;

  phist[2] =  0.0;

  phist[3] =  s*2.7E1-9.0/2.0;

  phist[4] =  t*2.7E1-9.0/2.0;

  phist[5] =  t*-2.7E1+9.0/2.0;

  phist[6] =  s*2.7E1+t*5.4E1-4.5E1/2.0;

  phist[7] =  s*5.4E1+t*2.7E1-4.5E1/2.0;

  phist[8] =  s*-2.7E1+9.0/2.0;

  phist[9] =  s*-5.4E1-t*5.4E1+2.7E1;


  // phitt
  phitt[0] =  s*-2.7E1-t*2.7E1+1.8E1;

  phitt[1] =  0.0;

  phitt[2] =  t*2.7E1-9.0;

  phitt[3] =  0.0;

  phitt[4] =  s*2.7E1;

  phitt[5] =  s*-2.7E1-t*8.1E1+3.6E1;

  phitt[6] =  s*5.4E1+t*8.1E1-4.5E1;

  phitt[7] =  s*2.7E1;

  phitt[8] =  0.0;

  phitt[9] =  s*-5.4E1;
}

void
BasisFunctionArea<Triangle,Lagrange,3>::coordinates( std::vector<Real>& s, std::vector<Real>& t ) const
{
  s = coords_s_;
  t = coords_t_;
}

void
BasisFunctionArea<Triangle,Lagrange,4>::evalBasis(
const Real& s, const Real& t, const Int3& sgn , Real phi[], int nphi) const
{
  SANS_ASSERT(nphi==15);

  // phi
  phi[0] =  (s*s)*(t*t)*6.4E1-ONETHIRD*s*2.5E1-ONETHIRD*t*2.5E1+ONETHIRD*(s*s)*7.0E1-ONETHIRD*(s*s*s)*8.0E1+
      ONETHIRD*(s*s*s*s)*3.2E1+ONETHIRD*(t*t)*7.0E1-ONETHIRD*(t*t*t)*8.0E1+ONETHIRD*(t*t*t*t)*3.2E1-s*(t*t)*
      8.0E1-(s*s)*t*8.0E1+ONETHIRD*s*t*1.4E2+ONETHIRD*s*(t*t*t)*1.28E2+ONETHIRD*(s*s*s)*t*1.28E2+1.0;

  phi[1] =  -s+ONETHIRD*(s*s)*2.2E1+ONETHIRD*(s*s*s*s)*3.2E1-(s*s*s)*1.6E1;

  phi[2] =  -t+ONETHIRD*(t*t)*2.2E1+ONETHIRD*(t*t*t*t)*3.2E1-(t*t*t)*1.6E1;

  phi[3] =  (s*s)*t*-3.2E1+ONETHIRD*s*t*1.6E1+ONETHIRD*(s*s*s)*t*1.28E2;

  phi[4] =  (s*s)*(t*t)*6.4E1+s*t*4.0-s*(t*t)*1.6E1-(s*s)*t*1.6E1;

  phi[5] =  s*(t*t)*-3.2E1+ONETHIRD*s*t*1.6E1+ONETHIRD*s*(t*t*t)*1.28E2;

  phi[6] =  ONETHIRD*t*1.6E1-ONETHIRD*(t*t)*1.12E2+ONETHIRD*(t*t*t)*2.24E2-ONETHIRD*(t*t*t*t)*1.28E2+s*
      (t*t)*3.2E1-ONETHIRD*s*t*1.6E1-ONETHIRD*s*(t*t*t)*1.28E2;

  phi[7] =  t*-1.2E1+(s*s)*(t*t)*6.4E1+s*t*2.8E1-s*(t*t)*1.44E2-(s*s)*t*1.6E1+s*(t*t*t)*1.28E2+(t*t)*7.6E1-
      (t*t*t)*1.28E2+(t*t*t*t)*6.4E1;

  phi[8] =  t*1.6E1-(s*s)*(t*t)*1.28E2-ONETHIRD*(t*t)*2.08E2-ONETHIRD*(t*t*t*t)*1.28E2+s*(t*t)*1.92E2+(s*
      s)*t*9.6E1-s*(t*t*t)*1.28E2+(t*t*t)*9.6E1-ONETHIRD*s*t*2.08E2-ONETHIRD*(s*s*s)*t*1.28E2;

  phi[9] =  s*1.6E1-(s*s)*(t*t)*1.28E2-ONETHIRD*(s*s)*2.08E2-ONETHIRD*(s*s*s*s)*1.28E2+s*(t*t)*9.6E1+(s*
      s)*t*1.92E2-(s*s*s)*t*1.28E2+(s*s*s)*9.6E1-ONETHIRD*s*t*2.08E2-ONETHIRD*s*(t*t*t)*1.28E2;

  phi[10] =  s*-1.2E1+(s*s)*(t*t)*6.4E1+s*t*2.8E1-s*(t*t)*1.6E1-(s*s)*t*1.44E2+(s*s*s)*t*1.28E2+(s*s)*7.6E1-
      (s*s*s)*1.28E2+(s*s*s*s)*6.4E1;

  phi[11] =  ONETHIRD*s*1.6E1-ONETHIRD*(s*s)*1.12E2+ONETHIRD*(s*s*s)*2.24E2-ONETHIRD*(s*s*s*s)*1.28E2+(s*
      s)*t*3.2E1-ONETHIRD*s*t*1.6E1-ONETHIRD*(s*s*s)*t*1.28E2;

  phi[12] =  (s*s)*(t*t)*2.56E2+s*t*9.6E1-s*(t*t)*2.24E2-(s*s)*t*2.24E2+s*(t*t*t)*1.28E2+(s*s*s)*t*1.28E2;

  phi[13] =  (s*s)*(t*t)*-1.28E2-s*t*3.2E1+s*(t*t)*3.2E1+(s*s)*t*1.6E2-(s*s*s)*t*1.28E2;

  phi[14] =  (s*s)*(t*t)*-1.28E2-s*t*3.2E1+s*(t*t)*1.6E2+(s*s)*t*3.2E1-s*(t*t*t)*1.28E2;

}

void
BasisFunctionArea<Triangle,Lagrange,4>::evalBasisDerivative(
const Real& s, const Real& t, const Int3& sgn , Real phis[], Real phit[], int nphi) const
{
  SANS_ASSERT(nphi==15);

  // phis
  phis[0] =  ONETHIRD*-2.5E1+ONETHIRD*s*1.4E2+ONETHIRD*t*1.4E2-s*t*1.6E2-ONETHIRD*(s*s)*2.4E2+ONETHIRD*
      (s*s*s)*1.28E2+ONETHIRD*(t*t*t)*1.28E2+s*(t*t)*1.28E2-(t*t)*8.0E1+ONETHIRD*(s*s)*t*3.84E2;

  phis[1] =  ONETHIRD*s*4.4E1+ONETHIRD*(s*s*s)*1.28E2-(s*s)*4.8E1-1.0;

  phis[2] =  0.0;

  phis[3] =  ONETHIRD*t*1.6E1-s*t*6.4E1+ONETHIRD*(s*s)*t*3.84E2;

  phis[4] =  t*4.0-s*t*3.2E1+s*(t*t)*1.28E2-(t*t)*1.6E1;

  phis[5] =  ONETHIRD*t*1.6E1+ONETHIRD*(t*t*t)*1.28E2-(t*t)*3.2E1;

  phis[6] =  ONETHIRD*t*-1.6E1-ONETHIRD*(t*t*t)*1.28E2+(t*t)*3.2E1;

  phis[7] =  t*2.8E1-s*t*3.2E1+s*(t*t)*1.28E2-(t*t)*1.44E2+(t*t*t)*1.28E2;

  phis[8] =  ONETHIRD*t*-2.08E2+s*t*1.92E2-s*(t*t)*2.56E2+(t*t)*1.92E2-(t*t*t)*1.28E2-ONETHIRD*(s*s)*t*
      3.84E2;

  phis[9] =  ONETHIRD*s*-4.16E2-ONETHIRD*t*2.08E2+s*t*3.84E2-ONETHIRD*(s*s*s)*5.12E2-ONETHIRD*(t*t*t)*1.28E2-
      s*(t*t)*2.56E2-(s*s)*t*3.84E2+(s*s)*2.88E2+(t*t)*9.6E1+1.6E1;

  phis[10] =  s*1.52E2+t*2.8E1-s*t*2.88E2+s*(t*t)*1.28E2+(s*s)*t*3.84E2-(s*s)*3.84E2+(s*s*s)*2.56E2-(t*
      t)*1.6E1-1.2E1;

  phis[11] =  ONETHIRD*1.6E1-ONETHIRD*s*2.24E2-ONETHIRD*t*1.6E1+s*t*6.4E1+ONETHIRD*(s*s)*6.72E2-ONETHIRD*
      (s*s*s)*5.12E2-ONETHIRD*(s*s)*t*3.84E2;

  phis[12] =  t*9.6E1-s*t*4.48E2+s*(t*t)*5.12E2+(s*s)*t*3.84E2-(t*t)*2.24E2+(t*t*t)*1.28E2;

  phis[13] =  t*-3.2E1+s*t*3.2E2-s*(t*t)*2.56E2-(s*s)*t*3.84E2+(t*t)*3.2E1;

  phis[14] =  t*-3.2E1+s*t*6.4E1-s*(t*t)*2.56E2+(t*t)*1.6E2-(t*t*t)*1.28E2;


  // phit
  phit[0] =  ONETHIRD*-2.5E1+ONETHIRD*s*1.4E2+ONETHIRD*t*1.4E2-s*t*1.6E2+ONETHIRD*(s*s*s)*1.28E2-ONETHIRD*
      (t*t)*2.4E2+ONETHIRD*(t*t*t)*1.28E2+(s*s)*t*1.28E2-(s*s)*8.0E1+ONETHIRD*s*(t*t)*3.84E2;

  phit[1] =  0.0;

  phit[2] =  ONETHIRD*t*4.4E1+ONETHIRD*(t*t*t)*1.28E2-(t*t)*4.8E1-1.0;

  phit[3] =  ONETHIRD*s*1.6E1+ONETHIRD*(s*s*s)*1.28E2-(s*s)*3.2E1;

  phit[4] =  s*4.0-s*t*3.2E1+(s*s)*t*1.28E2-(s*s)*1.6E1;

  phit[5] =  ONETHIRD*s*1.6E1-s*t*6.4E1+ONETHIRD*s*(t*t)*3.84E2;

  phit[6] =  ONETHIRD*1.6E1-ONETHIRD*s*1.6E1-ONETHIRD*t*2.24E2+s*t*6.4E1+ONETHIRD*(t*t)*6.72E2-ONETHIRD*
      (t*t*t)*5.12E2-ONETHIRD*s*(t*t)*3.84E2;

  phit[7] =  s*2.8E1+t*1.52E2-s*t*2.88E2+s*(t*t)*3.84E2+(s*s)*t*1.28E2-(s*s)*1.6E1-(t*t)*3.84E2+(t*t*t)*
      2.56E2-1.2E1;

  phit[8] =  ONETHIRD*s*-2.08E2-ONETHIRD*t*4.16E2+s*t*3.84E2-ONETHIRD*(s*s*s)*1.28E2-ONETHIRD*(t*t*t)*5.12E2-
      s*(t*t)*3.84E2-(s*s)*t*2.56E2+(s*s)*9.6E1+(t*t)*2.88E2+1.6E1;

  phit[9] =  ONETHIRD*s*-2.08E2+s*t*1.92E2-(s*s)*t*2.56E2+(s*s)*1.92E2-(s*s*s)*1.28E2-ONETHIRD*s*(t*t)*
      3.84E2;

  phit[10] =  s*2.8E1-s*t*3.2E1+(s*s)*t*1.28E2-(s*s)*1.44E2+(s*s*s)*1.28E2;

  phit[11] =  ONETHIRD*s*-1.6E1-ONETHIRD*(s*s*s)*1.28E2+(s*s)*3.2E1;

  phit[12] =  s*9.6E1-s*t*4.48E2+s*(t*t)*3.84E2+(s*s)*t*5.12E2-(s*s)*2.24E2+(s*s*s)*1.28E2;

  phit[13] =  s*-3.2E1+s*t*6.4E1-(s*s)*t*2.56E2+(s*s)*1.6E2-(s*s*s)*1.28E2;

  phit[14] =  s*-3.2E1+s*t*3.2E2-s*(t*t)*3.84E2-(s*s)*t*2.56E2+(s*s)*3.2E1;
}

void
BasisFunctionArea<Triangle,Lagrange,4>::evalBasisHessianDerivative(
const Real& s, const Real& t, const Int3& sgn , Real phiss[], Real phist[], Real phitt[], int nphi) const
{
  SANS_ASSERT(nphi==15);

  // phiss
  phiss[0] =  ONETHIRD*1.4E2-t*1.6E2-ONETHIRD*s*4.8E2+ONETHIRD*(s*s)*3.84E2+(t*t)*1.28E2+ONETHIRD*s*t*7.68E2;

  phiss[1] =  ONETHIRD*4.4E1-s*9.6E1+ONETHIRD*(s*s)*3.84E2;

  phiss[2] =  0.0;

  phiss[3] =  t*-6.4E1+ONETHIRD*s*t*7.68E2;

  phiss[4] =  t*-3.2E1+(t*t)*1.28E2;

  phiss[5] =  0.0;

  phiss[6] =  0.0;

  phiss[7] =  t*-3.2E1+(t*t)*1.28E2;

  phiss[8] =  t*1.92E2-(t*t)*2.56E2-ONETHIRD*s*t*7.68E2;

  phiss[9] =  ONETHIRD*-4.16E2+s*5.76E2+t*3.84E2-s*t*7.68E2-ONETHIRD*(s*s)*1.536E3-(t*t)*2.56E2;

  phiss[10] =  s*-7.68E2-t*2.88E2+s*t*7.68E2+(s*s)*7.68E2+(t*t)*1.28E2+1.52E2;

  phiss[11] =  ONETHIRD*-2.24E2+t*6.4E1+ONETHIRD*s*1.344E3-ONETHIRD*(s*s)*1.536E3-ONETHIRD*s*t*7.68E2;

  phiss[12] =  t*-4.48E2+s*t*7.68E2+(t*t)*5.12E2;

  phiss[13] =  t*3.2E2-s*t*7.68E2-(t*t)*2.56E2;

  phiss[14] =  t*6.4E1-(t*t)*2.56E2;


  // phist
  phist[0] =  ONETHIRD*1.4E2-s*1.6E2-t*1.6E2+s*t*2.56E2+ONETHIRD*(s*s)*3.84E2+ONETHIRD*(t*t)*3.84E2;

  phist[1] =  0.0;

  phist[2] =  0.0;

  phist[3] =  ONETHIRD*1.6E1-s*6.4E1+ONETHIRD*(s*s)*3.84E2;

  phist[4] =  s*-3.2E1-t*3.2E1+s*t*2.56E2+4.0;

  phist[5] =  ONETHIRD*1.6E1-t*6.4E1+ONETHIRD*(t*t)*3.84E2;

  phist[6] =  ONETHIRD*-1.6E1+t*6.4E1-ONETHIRD*(t*t)*3.84E2;

  phist[7] =  s*-3.2E1-t*2.88E2+s*t*2.56E2+(t*t)*3.84E2+2.8E1;

  phist[8] =  ONETHIRD*-2.08E2+s*1.92E2+t*3.84E2-s*t*5.12E2-ONETHIRD*(s*s)*3.84E2-(t*t)*3.84E2;

  phist[9] =  ONETHIRD*-2.08E2+s*3.84E2+t*1.92E2-s*t*5.12E2-ONETHIRD*(t*t)*3.84E2-(s*s)*3.84E2;

  phist[10] =  s*-2.88E2-t*3.2E1+s*t*2.56E2+(s*s)*3.84E2+2.8E1;

  phist[11] =  ONETHIRD*-1.6E1+s*6.4E1-ONETHIRD*(s*s)*3.84E2;

  phist[12] =  s*-4.48E2-t*4.48E2+s*t*1.024E3+(s*s)*3.84E2+(t*t)*3.84E2+9.6E1;

  phist[13] =  s*3.2E2+t*6.4E1-s*t*5.12E2-(s*s)*3.84E2-3.2E1;

  phist[14] =  s*6.4E1+t*3.2E2-s*t*5.12E2-(t*t)*3.84E2-3.2E1;


  // phitt
  phitt[0] =  ONETHIRD*1.4E2-s*1.6E2-ONETHIRD*t*4.8E2+ONETHIRD*(t*t)*3.84E2+(s*s)*1.28E2+ONETHIRD*s*t*7.68E2;

  phitt[1] =  0.0;

  phitt[2] =  ONETHIRD*4.4E1-t*9.6E1+ONETHIRD*(t*t)*3.84E2;

  phitt[3] =  0.0;

  phitt[4] =  s*-3.2E1+(s*s)*1.28E2;

  phitt[5] =  s*-6.4E1+ONETHIRD*s*t*7.68E2;

  phitt[6] =  ONETHIRD*-2.24E2+s*6.4E1+ONETHIRD*t*1.344E3-ONETHIRD*(t*t)*1.536E3-ONETHIRD*s*t*7.68E2;

  phitt[7] =  s*-2.88E2-t*7.68E2+s*t*7.68E2+(s*s)*1.28E2+(t*t)*7.68E2+1.52E2;

  phitt[8] =  ONETHIRD*-4.16E2+s*3.84E2+t*5.76E2-s*t*7.68E2-ONETHIRD*(t*t)*1.536E3-(s*s)*2.56E2;

  phitt[9] =  s*1.92E2-(s*s)*2.56E2-ONETHIRD*s*t*7.68E2;

  phitt[10] =  s*-3.2E1+(s*s)*1.28E2;

  phitt[11] =  0.0;

  phitt[12] =  s*-4.48E2+s*t*7.68E2+(s*s)*5.12E2;

  phitt[13] =  s*6.4E1-(s*s)*2.56E2;

  phitt[14] =  s*3.2E2-s*t*7.68E2-(s*s)*2.56E2;
}

void
BasisFunctionArea<Triangle,Lagrange,4>::coordinates( std::vector<Real>& s, std::vector<Real>& t ) const
{
  s = coords_s_;
  t = coords_t_;
}

void
BasisFunctionArea<Triangle,Lagrange,5>::evalBasis(
const Real& s, const Real& t, const Int3& sgn , Real phi[], int nphi) const
{
  SANS_ASSERT(nphi==21);

  // phi
  phi[0] =  s*(-1.37E2/1.2E1)-t*(1.37E2/1.2E1)+(s*s)*(t*t)*4.6875E2-(s*s)*(t*t*t)*2.604166666666667E2-(s*
      s*s)*(t*t)*2.604166666666667E2+s*t*(3.75E2/4.0)-s*(t*t)*2.65625E2-(s*s)*t*2.65625E2+s*(t*t*t)*(6.25E2/
          2.0)+(s*s*s)*t*(6.25E2/2.0)-s*(t*t*t*t)*1.302083333333333E2-(s*s*s*s)*t*1.302083333333333E2+(s*s)*(3.75E2/
              8.0)-(s*s*s)*8.854166666666667E1+(s*s*s*s)*(6.25E2/8.0)-(s*s*s*s*s)*(6.25E2/2.4E1)+(t*t)*(3.75E2/8.0)-
              (t*t*t)*8.854166666666667E1+(t*t*t*t)*(6.25E2/8.0)-(t*t*t*t*t)*(6.25E2/2.4E1)+1.0;

  phi[1] =  s-(s*s)*(1.25E2/1.2E1)+(s*s*s)*(8.75E2/2.4E1)-(s*s*s*s)*(6.25E2/1.2E1)+(s*s*s*s*s)*(6.25E2/
      2.4E1);

  phi[2] =  t-(t*t)*(1.25E2/1.2E1)+(t*t*t)*(8.75E2/2.4E1)-(t*t*t*t)*(6.25E2/1.2E1)+(t*t*t*t*t)*(6.25E2/
      2.4E1);

  phi[3] =  s*t*(-2.5E1/4.0)+(s*s)*t*5.729166666666667E1-(s*s*s)*t*(6.25E2/4.0)+(s*s*s*s)*t*1.302083333333333E2;

  phi[4] =  (s*s)*(t*t)*(-6.25E2/4.0)+(s*s*s)*(t*t)*2.604166666666667E2+(s*s)*t*(1.25E2/4.0)-(s*s*s)*t*
      (6.25E2/1.2E1)-ONESIXTH*s*t*2.5E1+ONESIXTH*s*(t*t)*1.25E2;

  phi[5] =  (s*s)*(t*t)*(-6.25E2/4.0)+(s*s)*(t*t*t)*2.604166666666667E2+s*(t*t)*(1.25E2/4.0)-s*(t*t*t)*
      (6.25E2/1.2E1)-ONESIXTH*s*t*2.5E1+ONESIXTH*(s*s)*t*1.25E2;

  phi[6] =  s*t*(-2.5E1/4.0)+s*(t*t)*5.729166666666667E1-s*(t*t*t)*(6.25E2/4.0)+s*(t*t*t*t)*1.302083333333333E2;

  phi[7] =  t*(-2.5E1/4.0)+s*t*(2.5E1/4.0)-s*(t*t)*5.729166666666667E1+s*(t*t*t)*(6.25E2/4.0)-s*(t*t*t*
      t)*1.302083333333333E2+(t*t)*6.354166666666667E1-(t*t*t)*2.135416666666667E2+(t*t*t*t)*2.864583333333333E2-
      (t*t*t*t*t)*1.302083333333333E2;

  phi[8] =  (s*s)*(t*t)*(-6.25E2/4.0)+(s*s)*(t*t*t)*2.604166666666667E2+ONETHIRD*t*5.0E1-s*t*(7.5E1/2.0)+
      s*(t*t)*3.229166666666667E2-s*(t*t*t)*7.8125E2-(t*t)*(3.25E2/2.0)+(t*t*t)*5.104166666666667E2-(t*t*t*
          t)*6.25E2+(t*t*t*t*t)*2.604166666666667E2+ONESIXTH*(s*s)*t*1.25E2+ONESIXTH*s*(t*t*t*t)*3.125E3;

  phi[9] =  t*-2.5E1+(s*s)*(t*t)*7.8125E2-(s*s)*(t*t*t)*7.8125E2-(s*s*s)*(t*t)*2.604166666666667E2+s*t*
      9.791666666666667E1-s*(t*t)*7.395833333333333E2-(s*s)*t*1.25E2+s*(t*t*t)*1.40625E3+(s*s*s)*t*(6.25E2/
          1.2E1)-s*(t*t*t*t)*7.8125E2+(t*t)*2.229166666666667E2-(t*t*t)*6.145833333333333E2+(t*t*t*t)*6.770833333333333E2-
          (t*t*t*t*t)*2.604166666666667E2;

  phi[10] =  t*2.5E1-(s*s)*(t*t)*1.09375E3+(s*s)*(t*t*t)*7.8125E2-s*t*1.604166666666667E2+s*(t*t)*7.395833333333333E2+
      (s*s)*t*3.697916666666667E2-s*(t*t*t)*1.09375E3-(s*s*s)*t*3.645833333333333E2+(s*s*s*s)*t*1.302083333333333E2-
      (t*t)*1.604166666666667E2+(t*t*t)*3.697916666666667E2-(t*t*t*t)*3.645833333333333E2+(t*t*t*t*t)*1.302083333333333E2+
      ONESIXTH*(s*s*s)*(t*t)*3.125E3+ONESIXTH*s*(t*t*t*t)*3.125E3;

  phi[11] =  s*2.5E1-(s*s)*(t*t)*1.09375E3+(s*s*s)*(t*t)*7.8125E2-s*t*1.604166666666667E2+s*(t*t)*3.697916666666667E2+
      (s*s)*t*7.395833333333333E2-s*(t*t*t)*3.645833333333333E2-(s*s*s)*t*1.09375E3+s*(t*t*t*t)*1.302083333333333E2-
      (s*s)*1.604166666666667E2+(s*s*s)*3.697916666666667E2-(s*s*s*s)*3.645833333333333E2+(s*s*s*s*s)*1.302083333333333E2+
      ONESIXTH*(s*s)*(t*t*t)*3.125E3+ONESIXTH*(s*s*s*s)*t*3.125E3;

  phi[12] =  s*-2.5E1+(s*s)*(t*t)*7.8125E2-(s*s)*(t*t*t)*2.604166666666667E2-(s*s*s)*(t*t)*7.8125E2+s*t*
      9.791666666666667E1-s*(t*t)*1.25E2-(s*s)*t*7.395833333333333E2+s*(t*t*t)*(6.25E2/1.2E1)+(s*s*s)*t*1.40625E3-
      (s*s*s*s)*t*7.8125E2+(s*s)*2.229166666666667E2-(s*s*s)*6.145833333333333E2+(s*s*s*s)*6.770833333333333E2-
      (s*s*s*s*s)*2.604166666666667E2;

  phi[13] =  (s*s)*(t*t)*(-6.25E2/4.0)+(s*s*s)*(t*t)*2.604166666666667E2+ONETHIRD*s*5.0E1-s*t*(7.5E1/2.0)+
      (s*s)*t*3.229166666666667E2-(s*s*s)*t*7.8125E2-(s*s)*(3.25E2/2.0)+(s*s*s)*5.104166666666667E2-(s*s*s*
          s)*6.25E2+(s*s*s*s*s)*2.604166666666667E2+ONESIXTH*s*(t*t)*1.25E2+ONESIXTH*(s*s*s*s)*t*3.125E3;

  phi[14] =  s*(-2.5E1/4.0)+s*t*(2.5E1/4.0)-(s*s)*t*5.729166666666667E1+(s*s*s)*t*(6.25E2/4.0)-(s*s*s*s)*
      t*1.302083333333333E2+(s*s)*6.354166666666667E1-(s*s*s)*2.135416666666667E2+(s*s*s*s)*2.864583333333333E2-
      (s*s*s*s*s)*1.302083333333333E2;

  phi[15] =  (s*s)*(t*t)*2.5E3-(s*s)*(t*t*t)*1.5625E3-(s*s*s)*(t*t)*1.5625E3+s*t*2.5E2+s*(t*t*t)*1.25E3+
      (s*s*s)*t*1.25E3-ONESIXTH*s*(t*t)*5.875E3-ONESIXTH*(s*s)*t*5.875E3-ONESIXTH*s*(t*t*t*t)*3.125E3-ONESIXTH*
      (s*s*s*s)*t*3.125E3;

  phi[16] =  (s*s)*(t*t)*(-1.71875E3)+(s*s)*(t*t*t)*7.8125E2+(s*s*s)*(t*t)*1.5625E3-s*t*1.25E2+s*(t*t)*
      2.8125E2+(s*s)*t*9.0625E2-s*(t*t*t)*(6.25E2/4.0)-(s*s*s)*t*1.5625E3+(s*s*s*s)*t*7.8125E2;

  phi[17] =  (s*s)*(t*t)*(6.25E2/2.0)-ONESIXTH*(s*s*s)*(t*t)*3.125E3+ONETHIRD*s*t*1.25E2-ONETHIRD*s*(t*
      t)*1.25E2+ONETHIRD*(s*s*s)*t*2.5E3-ONESIXTH*(s*s)*t*2.125E3-ONESIXTH*(s*s*s*s)*t*3.125E3;

  phi[18] =  (s*s)*(t*t)*(-1.71875E3)+(s*s)*(t*t*t)*1.5625E3+(s*s*s)*(t*t)*7.8125E2-s*t*1.25E2+s*(t*t)*
      9.0625E2+(s*s)*t*2.8125E2-s*(t*t*t)*1.5625E3-(s*s*s)*t*(6.25E2/4.0)+s*(t*t*t*t)*7.8125E2;

  phi[19] =  (s*s)*(t*t)*1.09375E3-(s*s)*(t*t*t)*7.8125E2-(s*s*s)*(t*t)*7.8125E2+s*t*(1.25E2/4.0)-s*(t*
      t)*(3.75E2/2.0)-(s*s)*t*(3.75E2/2.0)+s*(t*t*t)*(6.25E2/4.0)+(s*s*s)*t*(6.25E2/4.0);

  phi[20] =  (s*s)*(t*t)*(6.25E2/2.0)-ONESIXTH*(s*s)*(t*t*t)*3.125E3+ONETHIRD*s*t*1.25E2-ONETHIRD*(s*s)*
      t*1.25E2+ONETHIRD*s*(t*t*t)*2.5E3-ONESIXTH*s*(t*t)*2.125E3-ONESIXTH*s*(t*t*t*t)*3.125E3;
}

void
BasisFunctionArea<Triangle,Lagrange,5>::evalBasisDerivative(
const Real& s, const Real& t, const Int3& sgn , Real phis[], Real phit[], int nphi) const
{
  SANS_ASSERT(nphi==21);

  // phis
  phis[0] =  s*(3.75E2/4.0)+t*(3.75E2/4.0)-(s*s)*(t*t)*7.8125E2-s*t*5.3125E2+s*(t*t)*9.375E2+(s*s)*t*9.375E2-
      s*(t*t*t)*5.208333333333333E2-(s*s*s)*t*5.208333333333333E2-(s*s)*2.65625E2+(s*s*s)*(6.25E2/2.0)-(s*s*
          s*s)*1.302083333333333E2-(t*t)*2.65625E2+(t*t*t)*(6.25E2/2.0)-(t*t*t*t)*1.302083333333333E2-1.37E2/1.2E1;

  phis[1] =  s*(-1.25E2/6.0)+(s*s)*(8.75E2/8.0)-(s*s*s)*(6.25E2/3.0)+(s*s*s*s)*1.302083333333333E2+1.0;

  phis[2] =  0.0;

  phis[3] =  t*(-2.5E1/4.0)+s*t*1.145833333333333E2-(s*s)*t*4.6875E2+(s*s*s)*t*5.208333333333333E2;

  phis[4] =  (s*s)*(t*t)*7.8125E2-ONESIXTH*t*2.5E1+s*t*(1.25E2/2.0)+ONESIXTH*(t*t)*1.25E2-s*(t*t)*(6.25E2/
      2.0)-(s*s)*t*(6.25E2/4.0);

  phis[5] =  ONESIXTH*t*-2.5E1-s*(t*t)*(6.25E2/2.0)+s*(t*t*t)*5.208333333333333E2+(t*t)*(1.25E2/4.0)-(t*
      t*t)*(6.25E2/1.2E1)+ONESIXTH*s*t*2.5E2;

  phis[6] =  t*(-2.5E1/4.0)+(t*t)*5.729166666666667E1-(t*t*t)*(6.25E2/4.0)+(t*t*t*t)*1.302083333333333E2;

  phis[7] =  t*(2.5E1/4.0)-(t*t)*5.729166666666667E1+(t*t*t)*(6.25E2/4.0)-(t*t*t*t)*1.302083333333333E2;

  phis[8] =  t*(-7.5E1/2.0)+ONESIXTH*(t*t*t*t)*3.125E3-s*(t*t)*(6.25E2/2.0)+s*(t*t*t)*5.208333333333333E2+
      (t*t)*3.229166666666667E2-(t*t*t)*7.8125E2+ONESIXTH*s*t*2.5E2;

  phis[9] =  t*9.791666666666667E1-(s*s)*(t*t)*7.8125E2-s*t*2.5E2+s*(t*t)*1.5625E3+(s*s)*t*(6.25E2/4.0)-
      s*(t*t*t)*1.5625E3-(t*t)*7.395833333333333E2+(t*t*t)*1.40625E3-(t*t*t*t)*7.8125E2;

  phis[10] =  t*(-1.604166666666667E2)+s*t*7.395833333333333E2+ONESIXTH*(t*t*t*t)*3.125E3-s*(t*t)*2.1875E3-
      (s*s)*t*1.09375E3+s*(t*t*t)*1.5625E3+(s*s*s)*t*5.208333333333333E2+(t*t)*7.395833333333333E2-(t*t*t)*
      1.09375E3+ONESIXTH*(s*s)*(t*t)*9.375E3;

  phis[11] =  s*(-3.208333333333333E2)-t*1.604166666666667E2+(s*s)*(t*t)*2.34375E3+s*t*1.479166666666667E3-
      s*(t*t)*2.1875E3-(s*s)*t*3.28125E3+(s*s)*1.109375E3-(s*s*s)*1.458333333333333E3+(s*s*s*s)*6.510416666666667E2+
      (t*t)*3.697916666666667E2-(t*t*t)*3.645833333333333E2+(t*t*t*t)*1.302083333333333E2+ONESIXTH*s*(t*t*t)*
      6.25E3+ONESIXTH*(s*s*s)*t*1.25E4+2.5E1;

  phis[12] =  s*4.458333333333333E2+t*9.791666666666667E1-(s*s)*(t*t)*2.34375E3-s*t*1.479166666666667E3+
      s*(t*t)*1.5625E3+(s*s)*t*4.21875E3-s*(t*t*t)*5.208333333333333E2-(s*s*s)*t*3.125E3-(s*s)*1.84375E3+(s*
          s*s)*2.708333333333333E3-(s*s*s*s)*1.302083333333333E3-(t*t)*1.25E2+(t*t*t)*(6.25E2/1.2E1)-2.5E1;

  phis[13] =  ONETHIRD*5.0E1-s*3.25E2-t*(7.5E1/2.0)+(s*s)*(t*t)*7.8125E2+s*t*6.458333333333333E2+ONESIXTH*
      (t*t)*1.25E2-s*(t*t)*(6.25E2/2.0)-(s*s)*t*2.34375E3+(s*s)*1.53125E3-(s*s*s)*2.5E3+(s*s*s*s)*1.302083333333333E3+
      ONESIXTH*(s*s*s)*t*1.25E4;

  phis[14] =  s*1.270833333333333E2+t*(2.5E1/4.0)-s*t*1.145833333333333E2+(s*s)*t*4.6875E2-(s*s*s)*t*5.208333333333333E2-
      (s*s)*6.40625E2+(s*s*s)*1.145833333333333E3-(s*s*s*s)*6.510416666666667E2-2.5E1/4.0;

  phis[15] =  t*2.5E2-(s*s)*(t*t)*4.6875E3-ONESIXTH*(t*t)*5.875E3-ONESIXTH*(t*t*t*t)*3.125E3+s*(t*t)*5.0E3+
      (s*s)*t*3.75E3-s*(t*t*t)*3.125E3+(t*t*t)*1.25E3-ONESIXTH*s*t*1.175E4-ONESIXTH*(s*s*s)*t*1.25E4;

  phis[16] =  t*-1.25E2+(s*s)*(t*t)*4.6875E3+s*t*1.8125E3-s*(t*t)*3.4375E3-(s*s)*t*4.6875E3+s*(t*t*t)*1.5625E3+
      (s*s*s)*t*3.125E3+(t*t)*2.8125E2-(t*t*t)*(6.25E2/4.0);

  phis[17] =  ONETHIRD*t*1.25E2-ONETHIRD*(t*t)*1.25E2+s*(t*t)*6.25E2-ONESIXTH*(s*s)*(t*t)*9.375E3-ONESIXTH*
      s*t*4.25E3+ONETHIRD*(s*s)*t*7.5E3-ONESIXTH*(s*s*s)*t*1.25E4;

  phis[18] =  t*-1.25E2+(s*s)*(t*t)*2.34375E3+s*t*5.625E2-s*(t*t)*3.4375E3-(s*s)*t*4.6875E2+s*(t*t*t)*3.125E3+
      (t*t)*9.0625E2-(t*t*t)*1.5625E3+(t*t*t*t)*7.8125E2;

  phis[19] =  t*(1.25E2/4.0)-(s*s)*(t*t)*2.34375E3-s*t*3.75E2+s*(t*t)*2.1875E3+(s*s)*t*4.6875E2-s*(t*t*
      t)*1.5625E3-(t*t)*(3.75E2/2.0)+(t*t*t)*(6.25E2/4.0);

  phis[20] =  ONETHIRD*t*1.25E2+ONETHIRD*(t*t*t)*2.5E3-ONESIXTH*(t*t)*2.125E3-ONESIXTH*(t*t*t*t)*3.125E3+
      s*(t*t)*6.25E2-ONETHIRD*s*t*2.5E2-ONESIXTH*s*(t*t*t)*6.25E3;


  // phit
  phit[0] =  s*(3.75E2/4.0)+t*(3.75E2/4.0)-(s*s)*(t*t)*7.8125E2-s*t*5.3125E2+s*(t*t)*9.375E2+(s*s)*t*9.375E2-
      s*(t*t*t)*5.208333333333333E2-(s*s*s)*t*5.208333333333333E2-(s*s)*2.65625E2+(s*s*s)*(6.25E2/2.0)-(s*s*
          s*s)*1.302083333333333E2-(t*t)*2.65625E2+(t*t*t)*(6.25E2/2.0)-(t*t*t*t)*1.302083333333333E2-1.37E2/1.2E1;

  phit[1] =  0.0;

  phit[2] =  t*(-1.25E2/6.0)+(t*t)*(8.75E2/8.0)-(t*t*t)*(6.25E2/3.0)+(t*t*t*t)*1.302083333333333E2+1.0;

  phit[3] =  s*(-2.5E1/4.0)+(s*s)*5.729166666666667E1-(s*s*s)*(6.25E2/4.0)+(s*s*s*s)*1.302083333333333E2;

  phit[4] =  ONESIXTH*s*-2.5E1-(s*s)*t*(6.25E2/2.0)+(s*s*s)*t*5.208333333333333E2+(s*s)*(1.25E2/4.0)-(s*
      s*s)*(6.25E2/1.2E1)+ONESIXTH*s*t*2.5E2;

  phit[5] =  (s*s)*(t*t)*7.8125E2-ONESIXTH*s*2.5E1+s*t*(1.25E2/2.0)+ONESIXTH*(s*s)*1.25E2-s*(t*t)*(6.25E2/
      4.0)-(s*s)*t*(6.25E2/2.0);

  phit[6] =  s*(-2.5E1/4.0)+s*t*1.145833333333333E2-s*(t*t)*4.6875E2+s*(t*t*t)*5.208333333333333E2;

  phit[7] =  s*(2.5E1/4.0)+t*1.270833333333333E2-s*t*1.145833333333333E2+s*(t*t)*4.6875E2-s*(t*t*t)*5.208333333333333E2-
      (t*t)*6.40625E2+(t*t*t)*1.145833333333333E3-(t*t*t*t)*6.510416666666667E2-2.5E1/4.0;

  phit[8] =  ONETHIRD*5.0E1-s*(7.5E1/2.0)-t*3.25E2+(s*s)*(t*t)*7.8125E2+s*t*6.458333333333333E2+ONESIXTH*
      (s*s)*1.25E2-s*(t*t)*2.34375E3-(s*s)*t*(6.25E2/2.0)+(t*t)*1.53125E3-(t*t*t)*2.5E3+(t*t*t*t)*1.302083333333333E3+
      ONESIXTH*s*(t*t*t)*1.25E4;

  phit[9] =  s*9.791666666666667E1+t*4.458333333333333E2-(s*s)*(t*t)*2.34375E3-s*t*1.479166666666667E3+
      s*(t*t)*4.21875E3+(s*s)*t*1.5625E3-s*(t*t*t)*3.125E3-(s*s*s)*t*5.208333333333333E2-(s*s)*1.25E2+(s*s*
          s)*(6.25E2/1.2E1)-(t*t)*1.84375E3+(t*t*t)*2.708333333333333E3-(t*t*t*t)*1.302083333333333E3-2.5E1;

  phit[10] =  s*(-1.604166666666667E2)-t*3.208333333333333E2+(s*s)*(t*t)*2.34375E3+s*t*1.479166666666667E3-
      s*(t*t)*3.28125E3-(s*s)*t*2.1875E3+(s*s)*3.697916666666667E2-(s*s*s)*3.645833333333333E2+(s*s*s*s)*1.302083333333333E2+
      (t*t)*1.109375E3-(t*t*t)*1.458333333333333E3+(t*t*t*t)*6.510416666666667E2+ONESIXTH*s*(t*t*t)*1.25E4+
      ONESIXTH*(s*s*s)*t*6.25E3+2.5E1;

  phit[11] =  s*(-1.604166666666667E2)+s*t*7.395833333333333E2+ONESIXTH*(s*s*s*s)*3.125E3-s*(t*t)*1.09375E3-
      (s*s)*t*2.1875E3+s*(t*t*t)*5.208333333333333E2+(s*s*s)*t*1.5625E3+(s*s)*7.395833333333333E2-(s*s*s)*1.09375E3+
      ONESIXTH*(s*s)*(t*t)*9.375E3;

  phit[12] =  s*9.791666666666667E1-(s*s)*(t*t)*7.8125E2-s*t*2.5E2+s*(t*t)*(6.25E2/4.0)+(s*s)*t*1.5625E3-
      (s*s*s)*t*1.5625E3-(s*s)*7.395833333333333E2+(s*s*s)*1.40625E3-(s*s*s*s)*7.8125E2;

  phit[13] =  s*(-7.5E1/2.0)+ONESIXTH*(s*s*s*s)*3.125E3-(s*s)*t*(6.25E2/2.0)+(s*s*s)*t*5.208333333333333E2+
      (s*s)*3.229166666666667E2-(s*s*s)*7.8125E2+ONESIXTH*s*t*2.5E2;

  phit[14] =  s*(2.5E1/4.0)-(s*s)*5.729166666666667E1+(s*s*s)*(6.25E2/4.0)-(s*s*s*s)*1.302083333333333E2;

  phit[15] =  s*2.5E2-(s*s)*(t*t)*4.6875E3-ONESIXTH*(s*s)*5.875E3-ONESIXTH*(s*s*s*s)*3.125E3+s*(t*t)*3.75E3+
      (s*s)*t*5.0E3-(s*s*s)*t*3.125E3+(s*s*s)*1.25E3-ONESIXTH*s*t*1.175E4-ONESIXTH*s*(t*t*t)*1.25E4;

  phit[16] =  s*-1.25E2+(s*s)*(t*t)*2.34375E3+s*t*5.625E2-s*(t*t)*4.6875E2-(s*s)*t*3.4375E3+(s*s*s)*t*3.125E3+
      (s*s)*9.0625E2-(s*s*s)*1.5625E3+(s*s*s*s)*7.8125E2;

  phit[17] =  ONETHIRD*s*1.25E2+ONETHIRD*(s*s*s)*2.5E3-ONESIXTH*(s*s)*2.125E3-ONESIXTH*(s*s*s*s)*3.125E3+
      (s*s)*t*6.25E2-ONETHIRD*s*t*2.5E2-ONESIXTH*(s*s*s)*t*6.25E3;

  phit[18] =  s*-1.25E2+(s*s)*(t*t)*4.6875E3+s*t*1.8125E3-s*(t*t)*4.6875E3-(s*s)*t*3.4375E3+s*(t*t*t)*3.125E3+
      (s*s*s)*t*1.5625E3+(s*s)*2.8125E2-(s*s*s)*(6.25E2/4.0);

  phit[19] =  s*(1.25E2/4.0)-(s*s)*(t*t)*2.34375E3-s*t*3.75E2+s*(t*t)*4.6875E2+(s*s)*t*2.1875E3-(s*s*s)*
      t*1.5625E3-(s*s)*(3.75E2/2.0)+(s*s*s)*(6.25E2/4.0);

  phit[20] =  ONETHIRD*s*1.25E2-ONETHIRD*(s*s)*1.25E2+(s*s)*t*6.25E2-ONESIXTH*(s*s)*(t*t)*9.375E3-ONESIXTH*
      s*t*4.25E3+ONETHIRD*s*(t*t)*7.5E3-ONESIXTH*s*(t*t*t)*1.25E4;
}

void
BasisFunctionArea<Triangle,Lagrange,5>::evalBasisHessianDerivative(
const Real& s, const Real& t, const Int3& sgn , Real phiss[], Real phist[], Real phitt[], int nphi) const
{
  SANS_ASSERT(nphi==21);

  // phiss
  phiss[0] =  s*(-5.3125E2)-t*5.3125E2+s*t*1.875E3-s*(t*t)*1.5625E3-(s*s)*t*1.5625E3+(s*s)*9.375E2-(s*s*
      s)*5.208333333333333E2+(t*t)*9.375E2-(t*t*t)*5.208333333333333E2+3.75E2/4.0;

  phiss[1] =  s*(8.75E2/4.0)-(s*s)*6.25E2+(s*s*s)*5.208333333333333E2-1.25E2/6.0;

  phiss[2] =  0.0;

  phiss[3] =  t*1.145833333333333E2-s*t*9.375E2+(s*s)*t*1.5625E3;

  phiss[4] =  t*(1.25E2/2.0)-s*t*(6.25E2/2.0)+s*(t*t)*1.5625E3-(t*t)*(6.25E2/2.0);

  phiss[5] =  ONESIXTH*t*2.5E2-(t*t)*(6.25E2/2.0)+(t*t*t)*5.208333333333333E2;

  phiss[6] =  0.0;

  phiss[7] =  0.0;

  phiss[8] =  ONESIXTH*t*2.5E2-(t*t)*(6.25E2/2.0)+(t*t*t)*5.208333333333333E2;

  phiss[9] =  t*-2.5E2+s*t*(6.25E2/2.0)-s*(t*t)*1.5625E3+(t*t)*1.5625E3-(t*t*t)*1.5625E3;

  phiss[10] =  t*7.395833333333333E2-s*t*2.1875E3+(s*s)*t*1.5625E3-(t*t)*2.1875E3+(t*t*t)*1.5625E3+ONESIXTH*
      s*(t*t)*1.875E4;

  phiss[11] =  s*2.21875E3+t*1.479166666666667E3-s*t*6.5625E3+ONESIXTH*(t*t*t)*6.25E3+s*(t*t)*4.6875E3-
      (s*s)*4.375E3+(s*s*s)*2.604166666666667E3-(t*t)*2.1875E3+ONESIXTH*(s*s)*t*3.75E4-3.208333333333333E2;

  phiss[12] =  s*(-3.6875E3)-t*1.479166666666667E3+s*t*8.4375E3-s*(t*t)*4.6875E3-(s*s)*t*9.375E3+(s*s)*
      8.125E3-(s*s*s)*5.208333333333333E3+(t*t)*1.5625E3-(t*t*t)*5.208333333333333E2+4.458333333333333E2;

  phiss[13] =  s*3.0625E3+t*6.458333333333333E2-s*t*4.6875E3+s*(t*t)*1.5625E3-(s*s)*7.5E3+(s*s*s)*5.208333333333333E3-
      (t*t)*(6.25E2/2.0)+ONESIXTH*(s*s)*t*3.75E4-3.25E2;

  phiss[14] =  s*(-1.28125E3)-t*1.145833333333333E2+s*t*9.375E2-(s*s)*t*1.5625E3+(s*s)*3.4375E3-(s*s*s)*
      2.604166666666667E3+1.270833333333333E2;

  phiss[15] =  ONESIXTH*t*-1.175E4+s*t*7.5E3-s*(t*t)*9.375E3+(t*t)*5.0E3-(t*t*t)*3.125E3-ONESIXTH*(s*s)*
      t*3.75E4;

  phiss[16] =  t*1.8125E3-s*t*9.375E3+s*(t*t)*9.375E3+(s*s)*t*9.375E3-(t*t)*3.4375E3+(t*t*t)*1.5625E3;

  phiss[17] =  ONESIXTH*t*-4.25E3+(t*t)*6.25E2+ONETHIRD*s*t*1.5E4-ONESIXTH*s*(t*t)*1.875E4-ONESIXTH*(s*
      s)*t*3.75E4;

  phiss[18] =  t*5.625E2-s*t*9.375E2+s*(t*t)*4.6875E3-(t*t)*3.4375E3+(t*t*t)*3.125E3;

  phiss[19] =  t*-3.75E2+s*t*9.375E2-s*(t*t)*4.6875E3+(t*t)*2.1875E3-(t*t*t)*1.5625E3;

  phiss[20] =  ONETHIRD*t*-2.5E2-ONESIXTH*(t*t*t)*6.25E3+(t*t)*6.25E2;


  // phist
  phist[0] =  s*(-5.3125E2)-t*5.3125E2+s*t*1.875E3-s*(t*t)*1.5625E3-(s*s)*t*1.5625E3+(s*s)*9.375E2-(s*s*
      s)*5.208333333333333E2+(t*t)*9.375E2-(t*t*t)*5.208333333333333E2+3.75E2/4.0;

  phist[1] =  0.0;

  phist[2] =  0.0;

  phist[3] =  s*1.145833333333333E2-(s*s)*4.6875E2+(s*s*s)*5.208333333333333E2-2.5E1/4.0;

  phist[4] =  ONESIXTH*-2.5E1+s*(1.25E2/2.0)+ONESIXTH*t*2.5E2-s*t*6.25E2+(s*s)*t*1.5625E3-(s*s)*(6.25E2/
      4.0);

  phist[5] =  ONESIXTH*-2.5E1+t*(1.25E2/2.0)+ONESIXTH*s*2.5E2-s*t*6.25E2+s*(t*t)*1.5625E3-(t*t)*(6.25E2/
      4.0);

  phist[6] =  t*1.145833333333333E2-(t*t)*4.6875E2+(t*t*t)*5.208333333333333E2-2.5E1/4.0;

  phist[7] =  t*(-1.145833333333333E2)+(t*t)*4.6875E2-(t*t*t)*5.208333333333333E2+2.5E1/4.0;

  phist[8] =  t*6.458333333333333E2+ONESIXTH*s*2.5E2-s*t*6.25E2+ONESIXTH*(t*t*t)*1.25E4+s*(t*t)*1.5625E3-
      (t*t)*2.34375E3-7.5E1/2.0;

  phist[9] =  s*-2.5E2-t*1.479166666666667E3+s*t*3.125E3-s*(t*t)*4.6875E3-(s*s)*t*1.5625E3+(s*s)*(6.25E2/
      4.0)+(t*t)*4.21875E3-(t*t*t)*3.125E3+9.791666666666667E1;

  phist[10] =  s*7.395833333333333E2+t*1.479166666666667E3-s*t*4.375E3+ONESIXTH*(t*t*t)*1.25E4+s*(t*t)*
      4.6875E3-(s*s)*1.09375E3+(s*s*s)*5.208333333333333E2-(t*t)*3.28125E3+ONESIXTH*(s*s)*t*1.875E4-1.604166666666667E2;

  phist[11] =  s*1.479166666666667E3+t*7.395833333333333E2-s*t*4.375E3+ONESIXTH*(s*s*s)*1.25E4+(s*s)*t*
      4.6875E3-(s*s)*3.28125E3-(t*t)*1.09375E3+(t*t*t)*5.208333333333333E2+ONESIXTH*s*(t*t)*1.875E4-1.604166666666667E2;

  phist[12] =  s*(-1.479166666666667E3)-t*2.5E2+s*t*3.125E3-s*(t*t)*1.5625E3-(s*s)*t*4.6875E3+(s*s)*4.21875E3-
      (s*s*s)*3.125E3+(t*t)*(6.25E2/4.0)+9.791666666666667E1;

  phist[13] =  s*6.458333333333333E2+ONESIXTH*t*2.5E2-s*t*6.25E2+ONESIXTH*(s*s*s)*1.25E4+(s*s)*t*1.5625E3-
      (s*s)*2.34375E3-7.5E1/2.0;

  phist[14] =  s*(-1.145833333333333E2)+(s*s)*4.6875E2-(s*s*s)*5.208333333333333E2+2.5E1/4.0;

  phist[15] =  ONESIXTH*s*-1.175E4-ONESIXTH*t*1.175E4+s*t*1.0E4-ONESIXTH*(s*s*s)*1.25E4-ONESIXTH*(t*t*t)*
      1.25E4-s*(t*t)*9.375E3-(s*s)*t*9.375E3+(s*s)*3.75E3+(t*t)*3.75E3+2.5E2;

  phist[16] =  s*1.8125E3+t*5.625E2-s*t*6.875E3+s*(t*t)*4.6875E3+(s*s)*t*9.375E3-(s*s)*4.6875E3+(s*s*s)*
      3.125E3-(t*t)*4.6875E2-1.25E2;

  phist[17] =  ONETHIRD*1.25E2-ONESIXTH*s*4.25E3-ONETHIRD*t*2.5E2+s*t*1.25E3+ONETHIRD*(s*s)*7.5E3-ONESIXTH*
      (s*s*s)*1.25E4-ONESIXTH*(s*s)*t*1.875E4;

  phist[18] =  s*5.625E2+t*1.8125E3-s*t*6.875E3+s*(t*t)*9.375E3+(s*s)*t*4.6875E3-(s*s)*4.6875E2-(t*t)*4.6875E3+
      (t*t*t)*3.125E3-1.25E2;

  phist[19] =  s*-3.75E2-t*3.75E2+s*t*4.375E3-s*(t*t)*4.6875E3-(s*s)*t*4.6875E3+(s*s)*4.6875E2+(t*t)*4.6875E2+
      1.25E2/4.0;

  phist[20] =  ONETHIRD*1.25E2-ONETHIRD*s*2.5E2-ONESIXTH*t*4.25E3+s*t*1.25E3+ONETHIRD*(t*t)*7.5E3-ONESIXTH*
      (t*t*t)*1.25E4-ONESIXTH*s*(t*t)*1.875E4;


  // phitt
  phitt[0] =  s*(-5.3125E2)-t*5.3125E2+s*t*1.875E3-s*(t*t)*1.5625E3-(s*s)*t*1.5625E3+(s*s)*9.375E2-(s*s*
      s)*5.208333333333333E2+(t*t)*9.375E2-(t*t*t)*5.208333333333333E2+3.75E2/4.0;

  phitt[1] =  0.0;

  phitt[2] =  t*(8.75E2/4.0)-(t*t)*6.25E2+(t*t*t)*5.208333333333333E2-1.25E2/6.0;

  phitt[3] =  0.0;

  phitt[4] =  ONESIXTH*s*2.5E2-(s*s)*(6.25E2/2.0)+(s*s*s)*5.208333333333333E2;

  phitt[5] =  s*(1.25E2/2.0)-s*t*(6.25E2/2.0)+(s*s)*t*1.5625E3-(s*s)*(6.25E2/2.0);

  phitt[6] =  s*1.145833333333333E2-s*t*9.375E2+s*(t*t)*1.5625E3;

  phitt[7] =  s*(-1.145833333333333E2)-t*1.28125E3+s*t*9.375E2-s*(t*t)*1.5625E3+(t*t)*3.4375E3-(t*t*t)*
      2.604166666666667E3+1.270833333333333E2;

  phitt[8] =  s*6.458333333333333E2+t*3.0625E3-s*t*4.6875E3+(s*s)*t*1.5625E3-(s*s)*(6.25E2/2.0)-(t*t)*7.5E3+
      (t*t*t)*5.208333333333333E3+ONESIXTH*s*(t*t)*3.75E4-3.25E2;

  phitt[9] =  s*(-1.479166666666667E3)-t*3.6875E3+s*t*8.4375E3-s*(t*t)*9.375E3-(s*s)*t*4.6875E3+(s*s)*1.5625E3-
      (s*s*s)*5.208333333333333E2+(t*t)*8.125E3-(t*t*t)*5.208333333333333E3+4.458333333333333E2;

  phitt[10] =  s*1.479166666666667E3+t*2.21875E3-s*t*6.5625E3+ONESIXTH*(s*s*s)*6.25E3+(s*s)*t*4.6875E3-
      (s*s)*2.1875E3-(t*t)*4.375E3+(t*t*t)*2.604166666666667E3+ONESIXTH*s*(t*t)*3.75E4-3.208333333333333E2;

  phitt[11] =  s*7.395833333333333E2-s*t*2.1875E3+s*(t*t)*1.5625E3-(s*s)*2.1875E3+(s*s*s)*1.5625E3+ONESIXTH*
      (s*s)*t*1.875E4;

  phitt[12] =  s*-2.5E2+s*t*(6.25E2/2.0)-(s*s)*t*1.5625E3+(s*s)*1.5625E3-(s*s*s)*1.5625E3;

  phitt[13] =  ONESIXTH*s*2.5E2-(s*s)*(6.25E2/2.0)+(s*s*s)*5.208333333333333E2;

  phitt[14] =  0.0;

  phitt[15] =  ONESIXTH*s*-1.175E4+s*t*7.5E3-(s*s)*t*9.375E3+(s*s)*5.0E3-(s*s*s)*3.125E3-ONESIXTH*s*(t*
      t)*3.75E4;

  phitt[16] =  s*5.625E2-s*t*9.375E2+(s*s)*t*4.6875E3-(s*s)*3.4375E3+(s*s*s)*3.125E3;

  phitt[17] =  ONETHIRD*s*-2.5E2-ONESIXTH*(s*s*s)*6.25E3+(s*s)*6.25E2;

  phitt[18] =  s*1.8125E3-s*t*9.375E3+s*(t*t)*9.375E3+(s*s)*t*9.375E3-(s*s)*3.4375E3+(s*s*s)*1.5625E3;

  phitt[19] =  s*-3.75E2+s*t*9.375E2-(s*s)*t*4.6875E3+(s*s)*2.1875E3-(s*s*s)*1.5625E3;

  phitt[20] =  ONESIXTH*s*-4.25E3+(s*s)*6.25E2+ONETHIRD*s*t*1.5E4-ONESIXTH*s*(t*t)*3.75E4-ONESIXTH*(s*s)*
      t*1.875E4;
}

void
BasisFunctionArea<Triangle,Lagrange,5>::coordinates( std::vector<Real>& s, std::vector<Real>& t ) const
{
  s = coords_s_;
  t = coords_t_;
}

}
