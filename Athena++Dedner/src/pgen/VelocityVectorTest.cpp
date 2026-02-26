//Problem generator to test fluid solver/transport for velocity along an axis vs 45 degree velocity vector
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

//Function Headers
std::vector<Real> SphericalTransform(Real x, Real y, Real z);

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

                if(spherical[1]<5){
                    if(x>0){ 
                    //Set axial values here (along x-axis)
                    phydro->w(IDN,k,j,i) =  (rho_sw1AU / rho0) * std::pow(spherical[1],-2.0);
                    phydro->w(IPR,k,j,i) =  (thermPr_sw1AU / P0) * std::pow(spherical[1],-3.3333333);
                    phydro->w(IVX,k,j,i) =  (vel_sw1AU / v0);
                    phydro->w(IVY,k,j,i) =  0;
                    phydro->w(IVZ,k,j,i) =  0;
                    }

                    //Set 45deg values here (along z=0 plane)
                    else{ 
                        phydro->w(IDN,k,j,i) =  (rho_sw1AU / rho0) * std::pow(spherical[1],-2.0);
                        phydro->w(IPR,k,j,i) =  (thermPr_sw1AU / P0) * std::pow(spherical[1],-3.3333333);
                        phydro->w(IVX,k,j,i) =  (vel_sw1AU / v0) * std::cos(M_PI/4);
                        phydro->w(IVY,k,j,i) =  (vel_sw1AU / v0) * std::sin(M_PI/4);;
                        phydro->w(IVZ,k,j,i) =  0;
                    }   
                }
                else{
                    phydro->w(IDN,k,j,i) =  rho_ism / rho0;
                    phydro->w(IPR,k,j,i) =  thermPr_ism / P0;
                    phydro->w(IVX,k,j,i) =  0;
                    phydro->w(IVY,k,j,i) =  0;
                    phydro->w(IVZ,k,j,i) =  0;
                }
            }
        }
    }
    peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord, is, ie, js, je, ks, ke);
}