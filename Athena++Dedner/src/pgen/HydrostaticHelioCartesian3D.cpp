//Problem initialization file for a hydrodynamic Heliosphere
//This serves to set the values of density, pressure, momentum, etc across the entire simulation grid at t=0
//Boundary conditions for the inner radius should be input later on in the IN file
//We need to define a custom boundary condition for inflowing ISM in -x diection, done here in this file
//ASSUME AN ISM OF ONLY PROTONS, NO MAGNETIC FIELDS, NO CHARGE EXCHANGE
//Written by Luke Kathmann, 2025

//C++ Headers
#include <cmath>
#include <fstream>
#include <vector>

//Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

//Physical Constants
//We will use normalized variables to simplify the computations. Here are the scale constants that will be used as outlined in Hans document:
const double r0 = 1; //AU
const double v0 = 52.483; //km/s
const double t0 = 2.851e6; //s = 33 days
const double rho0 = 8.361e-27; //kg cm-3
const double P0 = 2.3e-11; //N/m^2
const double e0 = P0;

//Compute ISM, SW parameters
//Assume the LISM has the following plasma parameters: n=0.06 cm-3, v=26.24 km/s (-x direction), T=6530K
//Assume solar wind at 1AU has the following plasma parameters at 1 AU: n=5 cm-3, v=400 km s-1, T = 100000 K
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

//INJECTION RADIUS (in AU)
const double r_inj = 50.0;

//Function Headers
std::vector<Real> SphericalTransform(Real x, Real y, Real z);
void Boundary_ox1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh);
int RefinementCondition(MeshBlock *pmb);

//Set our user defined BC at outer radius. Set up adapative mesh refinement condition if on
void Mesh::InitUserMeshData(ParameterInput *pin){
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, Boundary_ox1);
    if (adaptive==true){
        EnrollUserRefinementCondition(RefinementCondition);
    }
}

//Intakes Cartesian coordinates, converts to spherical coordinate parameters needed for primitive, conserved variable transforms from Cartesian to spherical
std::vector<Real> SphericalTransform(Real x, Real y, Real z){
    std::vector<Real> transformVec(6);
    transformVec[0] = std::sqrt(x * x + y * y); //Polar Radius
    transformVec[1] = std::sqrt(x * x + y * y + z * z); //Spherical Radius
    transformVec[2] = y / transformVec[0]; //Sin_phi
    transformVec[3] = x / transformVec[0]; //Cos_phi
    transformVec[4] = transformVec[0] / transformVec[1]; //Sin_theta
    transformVec[5] = z / transformVec[1]; //Cos_theta
    return transformVec;
}

//Initialize mesh with SW values at 1AU, density decaying as 1/r^2 (spherical expansion), temperature as 1/r^4/3 (adiabatic, monoatomic gas decay), constant radial velocity within injection radius. Set ISM values outside
void MeshBlock::ProblemGenerator(ParameterInput *pin){

    //Notes: phydro->u is the array of conserved quantities

    //dimensions of meshblock
    const int Nx = ie - is + 1;
    const int Ny = je - js + 1;
    const int Nz = ke - ks + 1; 

    AthenaArray<Real> b; //needed for PrimitiveToConserved()
    b.NewAthenaArray(Nz, Ny, Nx);

    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
                Real x = pcoord->x1v(i);
                Real y = pcoord->x2v(j);
                Real z = pcoord->x3v(k);
                std::vector<Real> spherical = SphericalTransform(x, y, z);

                //Initialize values within injection radius
                if (spherical[1] < r_inj){
                    phydro->w(IDN,k,j,i) =  (rho_sw1AU / rho0) * std::pow(spherical[1],-2.0);
                    phydro->w(IPR,k,j,i) =  (thermPr_sw1AU / P0) * std::pow(spherical[1],-3.3333333);
                    phydro->w(IVX,k,j,i) =  (vel_sw1AU / v0) * x / spherical[1]; 
                    phydro->w(IVY,k,j,i) =  (vel_sw1AU / v0) * y / spherical[1]; 
                    phydro->w(IVZ,k,j,i) =  (vel_sw1AU / v0) * z / spherical[1]; 
                }

                //Initialize remaining values to ISM values
                else {
                    phydro->w(IDN,k,j,i) =  rho_ism / rho0;
                    phydro->w(IPR,k,j,i) =  thermPr_ism / P0;
                    phydro->w(IVX,k,j,i) =  -vel_ism / v0; //-x direction
                    phydro->w(IVY,k,j,i) =  0;
                    phydro->w(IVZ,k,j,i) =  0;
                }
            }   
        }
    }
    peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord, is, ie, js, je, ks, ke);
}


//Define our outer x BC function for ISM inflow
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

//Implementation for overwriting values within injection radius to create wind. Userworkinloop is called after EVERY timestep.
void MeshBlock::UserWorkInLoop() {

    //dimensions of meshblock
    const int Nx = ie - is + 1;
    const int Ny = je - js + 1;
    const int Nz = ke - ks + 1; 

    AthenaArray<Real> b; //needed for PrimitiveToConserved()
    b.NewAthenaArray(Nz, Ny, Nx);

  for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
                Real x = pcoord->x1v(i);
                Real y = pcoord->x2v(j);
                Real z = pcoord->x3v(k);
                std::vector<Real> spherical = SphericalTransform(x, y, z);

                //Overwrite only values within injection radius to create a continuous wind
                if (spherical[1] < r_inj){
                    phydro->w(IDN,k,j,i) =  (rho_sw1AU / rho0) * std::pow(spherical[1],-2.0);
                    phydro->w(IPR,k,j,i) =  (thermPr_sw1AU / P0) * std::pow(spherical[1],-3.3333333);
                    phydro->w(IVX,k,j,i) =  (vel_sw1AU / v0) * x / spherical[1]; 
                    phydro->w(IVY,k,j,i) =  (vel_sw1AU / v0) * y / spherical[1]; 
                    phydro->w(IVZ,k,j,i) =  (vel_sw1AU / v0) * z / spherical[1]; 
                }
            }
        }
    }
    peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord, is, ie, js, je, ks, ke);
}


//Refine based on Temperature gradient - INCOMPLETE
//ks, ke, js, je, is, ie are the bounds of the active (non-ghost) cell-centered zones in a MeshBlock
int RefinementCondition(MeshBlock *pmb) {
    Real gamma = pmb->peos->GetGamma();
    AthenaArray<Real> &w = pmb->phydro->w;
    //Initialize density gradient to zero
    Real maxTempGrad = 0.0;
    //SET DENSITY GRADIENT THRESHOLD FOR REFINEMENT
    Real threshold = 100000;

    for (int k=pmb->ks; k<=pmb->ke; k++) {
        for (int j=pmb->js; j<=pmb->je; j++) {
            for (int i=pmb->is; i<=pmb->ie; i++) {
                //Compute numerical gradients. T = 111111 (e/rho) = 111111*((pres/(gamma-1))/rho). P=2nkT (electrons and protons -> factor of 2)
                Real tempGradX = 111111 * ((std::abs((w(IPR,k,j,i+1) / w(IDN,k,j,i+1))  - (w(IPR,k,j,i) / w(IDN,k,j,i)))) / (gamma - 1)); /// (pmb->pcoord->dx1v(i));
                Real tempGradY = 111111 * ((std::abs((w(IPR,k,j+1,i) / w(IDN,k,j+1,i))  - (w(IPR,k,j,i) / w(IDN,k,j,i)))) / (gamma - 1)); /// (pmb->pcoord->dx1v(j));
                Real tempGradZ = 111111 * ((std::abs((w(IPR,k+1,j,i) / w(IDN,k+1,j,i))  - (w(IPR,k,j,i) / w(IDN,k,j,i)))) / (gamma - 1)); /// (pmb->pcoord->dx3v(k));
                Real tempGrad = std::sqrt(tempGradX*tempGradX + tempGradY*tempGradY + tempGradZ*tempGradZ);
                //If the newly computed density gradient is larger than the previously stored max gradient, replace it
                maxTempGrad = std::max(maxTempGrad, tempGrad);
                } 
            }           
        }
    if (maxTempGrad > threshold){
        return 1;     
    }
    else if (maxTempGrad < 0.5 * threshold){
        return -1;
    } 
    else{
        return 0;
    }
}