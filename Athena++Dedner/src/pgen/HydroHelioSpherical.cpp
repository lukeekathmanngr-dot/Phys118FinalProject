//Problem initialization file for a hydrodynamic Heliosphere
//This serves to set the values of density, pressure, momentum, etc across the entire simulation grid at t=0
//Boundary conditions for the inner radius should be input later on in the IN file
//We need to define a custom boundary condition for inflowing ISM in -x diection, done here in this file
//ASSUME AN ISM OF ONLY PROTONS, NO MAGNETIC FIELDS, NO CHARGE EXCHANGE

//C++ Headers
#include <cmath>
#include <fstream>
#include <string>

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
const double BOLTZMANN_CONSTANT = 1.38 * pow(10,-23); // J K-1
const double PROTON_MASS = 1.6726 * pow(10,-27); // kg
const double ISM_DENS = 0.06; //cm-3
const double ISM_TEMP = 6530.0; //K
const double SW_DENS1AU = 5.0; //cm-3
const double SW_TEMP1AU = 100000.0; //K

const double rho_ism = ISM_DENS * PROTON_MASS; //kg cm-3
const double thermPr_ism = 2 * ISM_DENS * BOLTZMANN_CONSTANT * ISM_TEMP * 1e6; // (J/m-3 = N/m^2), use ideal gas law P=2nkT (factor of 2 since both protons and electrons)
const double vel_ism = 26.24; //km s-1
const double vel_sw1AU = 400; //km s-1
const double rho_sw1AU = SW_DENS1AU * PROTON_MASS; //kg cm-3
const double thermPr_sw1AU = 2 * SW_DENS1AU * BOLTZMANN_CONSTANT * SW_TEMP1AU * 1e6; //(J/m^3 = N/m^2)

//Header for user defined BC at outer and inner radius
void Boundary_ox1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void Boundary_ix1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh);

//Set our user defined BC at outer radius
void Mesh::InitUserMeshData(ParameterInput *pin){
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, Boundary_ox1);
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, Boundary_ix1);
}

//Initialize mesh with SW values at 1AU, density decaying as 1/r^2 (spherical expansion), temperature as 1/r^4/3 (adiabatic, monoatomic gas decay), constant radial velocity
void MeshBlock::ProblemGenerator(ParameterInput *pin){

    Real gm1 = peos->GetGamma() - 1.0;

    //Note we need to normalize our quantities for dimensionless form of MHD eqns
    //Notes: phydro->u is the array of conserved quantities

    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            Real phi = pcoord->x2v(j); // polar angle
            for (int i=is; i<=ie; i++) {
                //Get radial coordinate
                Real r = pcoord->x1v(i);
                Real r_thresh = 10000; // set up ISM parameters past this value
                
                if (r<r_thresh){
                    //Initialize mesh
                    phydro->w(IDN,k,j,i) =  0; //rho_sw1AU / rho0) * pow(r,-2.0);
                    phydro->w(IVX,k,j,i) =  0; //(vel_sw1AU / v0);
                    phydro->w(IVY,k,j,i) =  0.0;
                    phydro->w(IVZ,k,j,i) =  0.0;
                    //Initialize energy density. Thermal pressure falls off as (1/r^2) * (1/r^4/3) since density times temperature.
                    phydro->w(IPR,k,j,i) =  0; //thermPr_sw1AU / P0 * pow(r,-10.0/3.0);
                }
                else{
                    phydro->w(IDN,k,j,i) =  (rho_ism / rho0);
                    phydro->w(IVX,k,j,i) =  -(vel_ism / v0) * std::cos(phi);
                    phydro->w(IVY,k,j,i) =  (vel_ism / v0) * std::sin(phi);
                    phydro->w(IVZ,k,j,i) =  0.0;
                    phydro->w(IPR,k,j,i) =  ((thermPr_ism)) / P0;
                }
            }
        }
    }
    peos->PrimitiveToConserved(phydro->w, phydro->w, phydro->u, pcoord, is, ie, js, je, ks, ke);
}


//Define our outer radius BC function
//Note il= is, iu=ie, but need to use the ones defined in this function due to SMR
//Note: only want pressure, mass density, etc from inflowing plasma if theta<pi/2, or else we artificially inject additional plasma pressure, mass on the left side of the mesh
void Boundary_ox1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh){
                        for (int k=kl; k<=ku; ++k){
                            for (int j=jl; j<=ju; ++j){
                                //Get phi (angular) coordinate (x2v is x2 'vector', cell centered)
                                Real phi = pco->x2v(j);
                                //Transform velocity in x direction to polar counterparts
                                Real v_r = -vel_ism * std::cos(phi);
                                Real v_phi = vel_ism * std::sin(phi);
                                //iu+1: first ghost cell beyond last active z1 cell. iu+ngh: last ghost cell
                                for (int i=iu+1; i<=iu+ngh; ++i){
                                    //Only assign values to points less than pi/2 in azimuth
                                    if (phi < M_PI/2){
                                        prim(IDN,k,j,i) = rho_ism / rho0;
                                        prim(IVX,k,j,i) = v_r / v0;
                                        prim(IVY,k,j,i) = v_phi / v0;
                                        prim(IVZ,k,j,i) = 0;
                                        prim(IPR,k,j,i) = thermPr_ism / P0;
                                    }
                                    else{
                                        prim(IDN,k,j,i) = prim(IDN,k,j,i-1);
                                        prim(IVX,k,j,i) = prim(IVX,k,j,i-1);
                                        prim(IVY,k,j,i) = prim(IVY,k,j,i-1);
                                        prim(IVZ,k,j,i) = prim(IVZ,k,j,i-1);
                                        prim(IPR,k,j,i) = prim(IPR,k,j,i-1);
                                    }
                                }
                            }
                        }           
}

//Define our inner radius BC function
void Boundary_ix1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh){
                        for (int k=kl; k<=ku; ++k){
                            for (int j=jl; j<=ju; ++j){
                                for (int i=il-ngh; i<=il-1; ++i){
                                    prim(IDN,k,j,i) = rho_sw1AU / rho0;
                                    prim(IVX,k,j,i) = vel_sw1AU / v0;
                                    prim(IVY,k,j,i) = 0;
                                    prim(IVZ,k,j,i) = 0;
                                    prim(IPR,k,j,i) = thermPr_sw1AU / P0;
                                }
                            }
                        }
                                   
}