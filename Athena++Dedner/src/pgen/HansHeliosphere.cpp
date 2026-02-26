// heliosphere.cpp:
// Spherically symmetric stellar wind with oncoming ISM flow
// 
// Requires Cartesian coordinates
// Fixes star at the origin
//
// Copyright(C) 2020-2021 Greg Szypko, Hans-Reinhard Mueller, and Noah Jerris

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
// NOTE: assumes this file has been moved to `src/pgen` in the Athena++ source tree
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

//========================================================================================
// Machine unit values
//========================================================================================

// FIXED VALUES:
// distance: 1 AU
// velocity: 52.483 km/s
// mass density: 2.8E13 kg/AU^3 (or 5 m_p/cm^3)
// mass: 2.8E13 kg
// DERIVED VALUES:
// time: 33.0 days
// pressure: 2.30E11 N/m^2
// energy density: 2.3E11 N/m^2
// magnetic field: 5.38E-9 T
// momentum: 262.4 (m_p cm^-3)(km s^-1)
// temperature: 111110 K

//========================================================================================
// Parameters
//========================================================================================

// Theoretical solar wind overwriting region:
Real rout;              // Outer radius of overwriting region
Real rramp;             // Radius to begin ramping the wind from theoretical to simulated
Real rin;               // Inner radius of overwriting region
Real rref;              // Reference radius (where most variables are given an initial value)
// NOTE: must satisfy `rin < rramp < rout`
Real d_wind;            // Density
Real e_wind;            // EXPLANATION FROM PROF. MUELLER HERE
Real b_rad_wind;        // EXPLANATION FROM PROF. MUELLER HERE
Real mom_rad_wind;      // EXPLANATION FROM PROF. MUELLER HERE

// ISM:
Real vx_ISM, vy_ISM, vz_ISM;  // Velocity components
Real d_ISM, p_ISM;            // Density and pressure
// NOTE: I think we may need to change the boundary conditions to allow inflow from up to three faces,
// depending on the ISM velocity

// Heliosphere Plasma Sheet:
bool plasma_sheet;      // If the plasma sheet is enabled
Real hps_height;        // EXPLANATION FROM PROF. MUELLER HERE
Real hps_factor;        // EXPLANATION FROM PROF. MUELLER HERE

// Magnetic fields:
bool split_monopole;    // If the split monopole condition is to be enforced
Real b0;                // EXPLANATION FROM PROF. MUELLER HERE
Real angle;             // EXPLANATION FROM PROF. MUELLER HERE
Real rin_mag;           // EXPLANATION FROM PROF. MUELLER HERE

// Miscellaneous:
Real amr_threshold;     // AMR threshold value (unused)

//========================================================================================
// Parker magnetic field
//========================================================================================

Real sign(Real x) {
  if (x > 0) {
    return 1.0;
  } else if (x < 0) {
    return -1.0;
  } else {
    return 0;
  }
}

std::vector<Real> SphToCart(Real x, Real y, Real z, Real *v) {
  std::vector<Real> vec;

  Real rxy = pow(x, 2) + pow(y, 2);
  if (rxy <= 0.05) {
    vec.push_back(sign(z) * v[1]);
    vec.push_back(v[2]);
    vec.push_back(sign(z) * v[0]);

    return vec;
  } else {
    Real rr = std::sqrt(rxy + pow(z, 2));
    rxy = std::sqrt(rxy);

    Real transform[3][3] = {
      {x, x * z / rxy, (-y) * rr / rxy},
      {y, y * z / rxy, x * rr / rxy},
      {z, -rxy, 0},
    };

    Real out[3];
    for (int i = 0; i < 3; i++) {
      out[i] = 0.0;
      for (int j = 0; j < 3; j++) {
        out[i] += transform[i][j] * v[j];
      }
    }

    for (int i = 0; i < 3; i++) {
      vec.push_back(out[i] / rr);
    }

    return vec;
  }
}

std::vector<Real> ParkerField(Real x, Real y, Real z) {
  Real omega = 2.7e-6;
  Real vm = 2.67e-6;
  Real b = 0.05 / 1.496;
  Real r = std::max(std::sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)), b);
  Real b0;
  if (z < 0) {
    b0 = -5;
  } else {
    b0 = 5;
  }

  Real sinth = std::sqrt(pow(x, 2) + pow(y, 2)) / r;

  Real magnitude = b0 * pow(b / r, 2);
  Real v[3];
  v[0] = magnitude;                                   // r
  v[1] = 0;                                           // theta
  v[2] = magnitude * -(omega / vm) * (r - b) * sinth; // rho

  return SphToCart(x, y, z, v);
}

//========================================================================================
// Athena forward declarations
//========================================================================================

int RefinementCondition(MeshBlock *pmb);

void WindSource(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
    const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar);

void ISMBoundary_ix1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh);

//========================================================================================
// Problem generator method implementations
//========================================================================================

// TODO: document
void Mesh::InitUserMeshData(ParameterInput *pin) {
  EnrollUserExplicitSourceFunction(WindSource);
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, ISMBoundary_ix1);
  // if (adaptive) {
  //   EnrollUserRefinementCondition(RefinementCondition);
  //   threshold = pin->GetReal("problem","thr");
  // }
}

// Validate input and initialize the simulation
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Check that coordinates are Cartesian
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in heliosphere.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM= " << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  // Check that `rin < rramp < rout` condition holds
  rout = pin->GetReal("problem", "rout");
  rramp = pin->GetReal("problem", "rramp");
  rin  = pin->GetReal("problem", "rin");
  if (!(rin < rramp && rramp < rout)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in heliosphere.cpp ProblemGenerator" << std::endl
        << "condition `rin < rramp < rout` not satisfied" << std::endl;
    ATHENA_ERROR(msg);
  }
  rref = pin->GetReal("problem", "rref");

  vx_ISM = pin->GetOrAddReal("problem", "vx_ISM", 0.0);
  vy_ISM = pin->GetOrAddReal("problem", "vy_ISM", 0.0);
  vz_ISM = pin->GetOrAddReal("problem", "vz_ISM", 0.0);
  d_ISM = pin->GetOrAddReal("problem", "d_ISM", 0.0);
  p_ISM = pin->GetOrAddReal("problem", "p_ISM", 0.0);
  rin_mag = pin->GetOrAddReal("problem", "rin_mag", 1.0);
  b_rad_wind = pin->GetOrAddReal("problem", "b_rad_wind", 0.0);
  e_wind = pin->GetOrAddReal("problem", "e_wind", 0.0);
  mom_rad_wind = pin->GetOrAddReal("problem", "mom_rad_wind", 0.0);
  d_wind = pin->GetOrAddReal("problem", "d_wind", 0.0);
  split_monopole = pin->GetOrAddBoolean("problem", "split_monopole", false);
  plasma_sheet = pin->GetOrAddBoolean("problem", "plasma_sheet", false);
  hps_height = pin->GetOrAddReal("problem", "hps_height", 0.01);
  hps_factor = pin->GetOrAddReal("problem", "hps_factor", 4.0);
  Real pa  = p_ISM;
  Real da  = d_ISM;

  if (MAGNETIC_FIELDS_ENABLED) {
    b0 = pin->GetReal("problem", "b0");
    angle = (PI/180.0)*pin->GetReal("problem", "angle");
  }
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;

  // setup uniform ambient medium
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        
        Real den = da;

        phydro->u(IDN,k,j,i) = d_ISM;

        // MOMENTUM-INITIALIZED ISM
        phydro->u(IM1,k,j,i) = d_ISM*vx_ISM;
        phydro->u(IM2,k,j,i) = d_ISM*vy_ISM;
        phydro->u(IM3,k,j,i) = d_ISM*vz_ISM;

        //BOUNDARY-ONLY ISM MOMENTUM
        // phydro->u(IM1,k,j,i) = 0.0;
        // phydro->u(IM2,k,j,i) = 0.0;
        // phydro->u(IM3,k,j,i) = 0.0;

        if (NON_BAROTROPIC_EOS) {
          Real pres = pa;
          phydro->u(IEN,k,j,i) = pres/gm1; //thermal energy
          phydro->u(IEN,k,j,i) += 0.5*d_ISM*(vx_ISM*vx_ISM+vy_ISM*vy_ISM+vz_ISM*vz_ISM); //bulk kinetic energy
          if (RELATIVISTIC_DYNAMICS)  // this should only ever be SR with this file
            phydro->u(IEN,k,j,i) += den;
        }
      }
    }
  }
              // Real falloff = rin*rin/(rad*rad);
              // Real curr_b_rad_wind = b_rad_wind * falloff;
              // if(y<0.0) curr_b_rad_wind *= -1.0; //enforce split monopole
              // cons(IB1,k,j,i) = curr_b_rad_wind * x / rad;
              // cons(IB2,k,j,i) = curr_b_rad_wind * y / rad;
              // cons(IB3,k,j,i) = curr_b_rad_wind * z / rad;

  // initialize solar wind magnetic field
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          Real x = pcoord->x1v(i);
          Real y = pcoord->x2v(j);
          Real z = pcoord->x3v(k);
          std::vector<Real> bfield = ParkerField(x, y, z);
          pfield->b.x1f(k,j,i) = bfield[0];
          pfield->b.x2f(k,j,i) = bfield[1];
          pfield->b.x3f(k,j,i) = bfield[2];
        }
      }
    }
  }
}

// Boundary condition definition
void ISMBoundary_ix1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt,
                int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set primitive variables in inlet ghost zones
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
          prim(IDN,k,j,il-i) = d_ISM;
          prim(IVX,k,j,il-i) = vx_ISM;
          prim(IVY,k,j,il-i) = vy_ISM;
          prim(IVZ,k,j,il-i) = vz_ISM;
          prim(IPR,k,j,il-i) = p_ISM;
      }
    }
  }

  // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
            // b.x1f(k,j,il-i) = b0;
            // b.x2f(k,j,il-i) = 0.0;
            // b.x3f(k,j,il-i) = 0.0;
            b.x1f(k,j,il-i) = b0*std::cos(angle);
            b.x2f(k,j,il-i) = b0*std::sin(angle);
            b.x3f(k,j,il-i) = 0.0;
        }
      }
    }
  }
}

// User-defined source term, here overwriting a spherical region around the origin to produce a
// solar wind
void WindSource(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
    const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar) {

  // Iterate over all cells in the mesh
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real x = pmb->pcoord->x1v(i);
        Real y = pmb->pcoord->x2v(j);
        Real z = pmb->pcoord->x3v(k);
        Real rad = std::sqrt(SQR(x) + SQR(y) + SQR(z));

        // Overwriting only occurs in cells such that the radius is less than the overwriting limit `rout`
        if (rad < rout) {
          if (MAGNETIC_FIELDS_ENABLED) {
            // Calculate the magnetic field vector
            std::vector<Real> bfield = ParkerField(x, y, z);
            pmb->pfield->b.x1f(k,j,i) = bfield[0];
            pmb->pfield->b.x2f(k,j,i) = bfield[1];
            pmb->pfield->b.x3f(k,j,i) = bfield[2];
          }

          if (rad < rin) {  // Within innermost overwriting region
            // Increase momentum linearly to its value at rref
            // Real md_factor = pow(rref/rin, 2);
            Real md_factor = SQR(rref/rin);
            Real curr_mom_rad_wind = mom_rad_wind * md_factor * rad / rin;

            // Keep density constant to product accelerating winds, but possibly account
            // for HPS
            Real idn = d_wind * md_factor;
            if(plasma_sheet && std::abs(z)<=hps_height){
              Real t = (hps_height - std::abs(z))/hps_height; //varies from 0 at outer edge to 1 at equator
              Real factor = (3.0*t*t - 2.0*t*t*t) * (hps_factor-1.0) + 1.0; //varies factor from 1 at outer edge to 4 at equator (cubic Hermite spline)
              curr_mom_rad_wind *= factor;
              idn *= factor;
            } 
            Real e_factor = pow(rref/rin, 3.333333333);
            cons(IDN,k,j,i) = idn;
            cons(IM1,k,j,i) = curr_mom_rad_wind * x / rad;
            cons(IM2,k,j,i) = curr_mom_rad_wind * y / rad;
            cons(IM3,k,j,i) = curr_mom_rad_wind * z / rad;
            cons(IEN,k,j,i) = e_wind * e_factor + 0.5*curr_mom_rad_wind*curr_mom_rad_wind/idn;
          } else {  // Within normal overwriting region
            // If `rad < rramp`, we want conserved variables to equal their theoretical values.
            // Otherwise, we want them to be the weighted average of their theoretical and simulated
            // values. Calculate the theoretical values first, then determine if those values or
            // the averaged values are needed.
            Real idn, im1, im2, im3, ien, ien2, ib1, ib2, ib3;

            // Compute theoretical values

            // Decrease momentum and density by factor rad^(-2)
            // Real md_factor = pow(rref/rad, 2);
            Real md_factor = SQR(rref/rad);
            idn = d_wind * md_factor;
            Real curr_mom_rad_wind = mom_rad_wind * md_factor;

            // Alter momentum and density for HPS if enabled
            if(plasma_sheet && std::abs(z)<=hps_height){
              Real t = (hps_height - std::abs(z))/hps_height; //varies from 0 at outer edge to 1 at equator
              Real factor = (3.0*t*t - 2.0*t*t*t) * (hps_factor-1.0) + 1.0; //varies factor from 1 at outer edge to 4 at equator (cubic Hermite spline)
              curr_mom_rad_wind *= factor;
              idn *= factor;
            }

            // Momentum components
            im1 = curr_mom_rad_wind * x / rad;
            im2 = curr_mom_rad_wind * y / rad;
            im3 = curr_mom_rad_wind * z / rad;

            // Decrease e by factor rad^(-10/3)
            // Real e_factor = pow(rref/rad, 10/3);
            // Real e_factor = pow( std::cbrt(rref/rad), 10);
            Real e_factor = pow(rref/rad, 3.333333333);
            ien  = e_wind * e_factor;
            ien2 = 0.5*curr_mom_rad_wind*curr_mom_rad_wind/idn;

            // Apply theoretical or hybrid values, depending on radius
            if ( (rad < rramp) || (time < 10.0)) {  // Within purely theoretical overwriting region
              cons(IDN,k,j,i) = idn;
              cons(IM1,k,j,i) = im1;
              cons(IM2,k,j,i) = im2;
              cons(IM3,k,j,i) = im3;
              cons(IEN,k,j,i) = ien + ien2;
            } else {            // Within ramped overwriting region
              Real ws = (rad - rramp)/(rout - rramp);   // Weight for simulated values
              Real wt = 1 - ws;                         // Weight for theoretical values
          //  Real eold = cons(IEN,k,j,i) - 0.5*
          //       ( pow(cons(IM1,k,j,i),2)+pow(cons(IM2,k,j,i),2)+pow(cons(IM3,k,j,i),2) )
          //         /cons(IDN,k,j,i);
          //  cons(IDN,k,j,i) = ws * cons(IDN,k,j,i) + wt * idn;
          //  cons(IM1,k,j,i) = ws * cons(IM1,k,j,i) + wt * im1;
          //  cons(IM2,k,j,i) = ws * cons(IM2,k,j,i) + wt * im2;
          //  cons(IM3,k,j,i) = ws * cons(IM3,k,j,i) + wt * im3;
              cons(IDN,k,j,i) = idn;
              cons(IM1,k,j,i) = im1;
              cons(IM2,k,j,i) = im2;
              cons(IM3,k,j,i) = im3;
              cons(IEN,k,j,i) = ien + ien2;
          //  cons(IEN,k,j,i) =(ws * eold            + wt * ien) + 0.5*
          //       ( pow(cons(IM1,k,j,i),2)+pow(cons(IM2,k,j,i),2)+pow(cons(IM3,k,j,i),2) )
          //         /cons(IDN,k,j,i);
            }
          }  // else case of rad<rin
        }  // if rad<rout
      }   // i
    }    // j
  }     // k
}      // end routine
