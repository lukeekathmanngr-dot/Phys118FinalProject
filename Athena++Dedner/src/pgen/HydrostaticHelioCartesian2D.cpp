//Problem initialization file for a hydrodynamic Heliosphere
//This serves to set the values of density, pressure, momentum, etc across the entire simulation grid at t=0
//Boundary conditions for the inner radius should be input later on in the IN file
//We need to define a custom boundary condition for inflowing ISM in -x diection, done here in this file
//ASSUME AN ISM OF ONLY PROTONS, NO MAGNETIC FIELDS, NO CHARGE EXCHANGE

//C++ Headers
#include <cmath>
#include <fstream>

//Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

//Physical Constants
//We will use normalized variables to simplify the computations. Here are the scale constants that will be used:
const double r0 = 1; //AU
const double v0 = 52.483; //km/s
const double t0 = 2.851e6; //s = 33 days
const double rho0 = 8.361e-27; //kg cm-3
const double P0 = 2.3e-11; //N/m^2
const double e0 = P0;

//Compute ISM, SW parameters
//Assume the LISM has the following plasma parameters: n=0.06 cm-3, v=26.24 km/s (-x direction), T=6530K
//Assume solar wind at 1AU has the following plasma parameters: n=5 cm-3, v=400 km s-1, T = 100000 K
//Assume Adibatic EOS with gamma=5/3
const double BOLTZMANN_CONSTANT = 1.38e-23; // J K-1
const double PROTON_MASS = 1.6726e-27; // kg
const double ISM_DENS = 0.1; //cm-3
const double ISM_TEMP = 6530.0; //K
const double SW_DENS1AU = 5.0; //cm-3
const double SW_TEMP1AU = 100000.0; //K

const double rho_ism = ISM_DENS * PROTON_MASS; //kg cm-3
const double thermPr_ism = 2 * ISM_DENS * BOLTZMANN_CONSTANT * ISM_TEMP * 1e6; // (J/m-3 = N/m^2), use ideal gas law P=2nkT (factor of 2 since both protons and electrons)
const double vel_ism = 26.53; //km s-1
const double vel_sw1AU = 400; //km s-1
const double rho_sw1AU = SW_DENS1AU * PROTON_MASS; //kg cm-3
const double thermPr_sw1AU = 2 * SW_DENS1AU * BOLTZMANN_CONSTANT * SW_TEMP1AU * 1e6; //(J/m^3 = N/m^2)

//Headers for user defined BC and solar wind source term
void Boundary_ox1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void WindInjection(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar);      

//Set our user defined BC at outer radius, set source term
void Mesh::InitUserMeshData(ParameterInput *pin){
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, Boundary_ox1);
    EnrollUserExplicitSourceFunction(WindInjection);
}

//Initialize mesh with SW values at 1AU, density decaying as 1/r^2 (spherical expansion), temperature as 1/r^4/3 (adiabatic, monoatomic gas decay), constant radial velocity
void MeshBlock::ProblemGenerator(ParameterInput *pin){

    //Note we need to normalize our quantities for dimensionless form of MHD eqns
    //Notes: phydro->u is the array of conserved quantities

    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
                //Get x,y positions
                Real x = pcoord->x1v(i);
                Real y = pcoord->x2v(j);
                Real z = pcoord->x2v(k);
                Real r = std::sqrt(x * x + y * y);
                
                //Initialize grid with solar wind values. Density falls as 1/r^2
                phydro->w(IDN,k,j,i) =  (rho_sw1AU / rho0) * std::pow(r,-2.0);
                //Initialize grid with velocity components
                phydro->w(IVX,k,j,i) =  (vel_sw1AU / v0) * (x / r);
                phydro->w(IVY,k,j,i) =  (vel_sw1AU / v0) * (y / r);
                phydro->w(IVZ,k,j,i) =  0.0;
                //Initialize energy density. Thermal pressure falls off as (1/r^2) * (1/r^4/3) since density times temperature.
                phydro->w(IPR,k,j,i) =  (thermPr_sw1AU / P0) * std::pow(r,-10.0/3.0);
            }
        }
    }
    peos->PrimitiveToConserved(phydro->w, phydro->w, phydro->u, pcoord, is, ie, js, je, ks, ke);
}

//Define our outer x BC function
//Note il= is, iu=ie, but need to use the ones defined in this function due to SMR
void Boundary_ox1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh){
                        for (int k=kl; k<=ku; ++k){
                            for (int j=jl; j<=ju; ++j){
                                //iu+1: first ghost cell beyond last active z1 cell. iu+ngh: last ghost cell
                                //Set primitive variables
                                for (int i=iu+1; i<=iu+ngh; ++i){
                                    prim(IDN,k,j,i) = rho_ism / rho0;
                                    prim(IVX,k,j,i) = -vel_ism / v0; //x-direction
                                    prim(IVY,k,j,i) = 0;
                                    prim(IVZ,k,j,i) = 0;
                                    prim(IPR,k,j,i) = thermPr_ism / P0;
                                }
                            }
                        }           
}

//Define our solar wind source term
void WindInjection(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar){
                Real r_inj = 2.0; // Injection radius in code units (AU)

                Coordinates *coords = pmb->pcoord;
                auto *peos = pmb->peos;

                for (int k = pmb->ks; k <= pmb->ke; ++k) {
                    for (int j = pmb->js; j <= pmb->je; ++j) {
                        for (int i = pmb->is; i <= pmb->ie; ++i) {
                            Real x = coords->x1v(i);
                            Real y = coords->x2v(j);
                            Real r = std::sqrt(x * x + y * y);
                            Real gamma = peos->GetGamma();
                            if (r <= r_inj) {
                                //Create injection primitive values 
                                Real dens = (rho_sw1AU / rho0) * std::pow(r,-2.0);
                                Real vx = (vel_sw1AU / v0) * (x / r);
                                Real vy = (vel_sw1AU / v0) * (y / r);
                                Real vz = 0;
                                Real pres = (thermPr_sw1AU / P0) * std::pow(r,-10.0/3.0);

                                //Convert primitive to conserved, update source terms in conserved variables by multiplying dt
                                cons(IDN,k,j,i) += dt * dens;
                                cons(IM1,k,j,i) += dt * vx * dens;
                                cons(IM2,k,j,i) += dt * vy * dens;
                                cons(IM3,k,j,i) += 0.0;
                                cons(IEN,k,j,i) += dt * ((pres / (gamma - 1)) + (0.5 * dens * (vx * vx + vy * vy)));
                            }
                        }
                    }
                }
}

/*
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin){
    std::ofstream output("SimValues.csv");
    output << "x,y,dens,pres" << "\n";
    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
                if (pcoord->x1v(i) < 1.0 && pcoord->x1v(i) > -1.0){
                    output << pcoord->x1v(i) << "," << pcoord->x2v(j) << "," << phydro->w(IDN,k,j,i) << "," << phydro->w(IPR,k,j,i) << '\n';
                }
            }
        }
    }
    output.close();
}
*/

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
    std::ofstream output("SimValues.csv");
    output << "x,y,dens,pres,velx,vely" << "\n";
    //Loop through all mesh blocks
    for (int b=0; b<nblocal; ++b) { 
        MeshBlock *pmb = my_blocks(b);
        for (int k=pmb->ks; k<=pmb->ke; k++) {
            for (int j=pmb->js; j<=pmb->je; j++) {
                for (int i=pmb->is; i<=pmb->ie; i++) {
                    Real x = pmb->pcoord->x1v(i);
                    Real y = pmb->pcoord->x2v(j);
                    Real vx = pmb->phydro->w(IVX,k,j,i);
                    Real vy = pmb->phydro->w(IVY,k,j,i);
                    if (x < 0.1 && x > -0.1){
                        output << x << "," << y << "," << pmb->phydro->w(IDN,k,j,i) << "," << pmb->phydro->w(IPR,k,j,i) << ',' << vx << ',' << vy << '\n';
                    }
                }
            }
        }
    }
    output.close();
}


/*
void MeshBlock::UserWorkInLoop() {
  Real gamma = peos->GetGamma();
  Real r_inj = 2.0; // Injection radius in code units (AU)
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real r_dist = std::sqrt(x * x + y * y);
        Real r = std::max(r_dist, 1e-6);
        if (r <= r_inj) {
            // Enforce solar wind primitives inside the injection radius
            phydro->w(IDN,k,j,i) =  (rho_sw1AU / rho0) * std::pow(r,-2.0);
            phydro->w(IVX,k,j,i) =  (vel_sw1AU / v0) * (x / r);
            phydro->w(IVY,k,j,i) =  (vel_sw1AU / v0) * (y / r);
            phydro->w(IVZ,k,j,i) =  0.0;
            phydro->w(IPR,k,j,i) =  (thermPr_sw1AU / P0) * std::pow(r,-10.0/3.0);
        }
      }
    }
  }
  peos->PrimitiveToConserved(phydro->w, phydro->w, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
*/