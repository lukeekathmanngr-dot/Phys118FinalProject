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
const double BOLTZMANN_CONSTANT = 1.38e-23; // J K-1
const double PROTON_MASS = 1.6726e-27; // kg
const double rho_inner = 1e-8; //MUST BE HIGH ENOUGH TO PREVENT NUMERICAL INSTABILITY, LIKELY TO DO WITH SHARP GRADIENTS
const double press_inner= 1e-12; //MUST BE HIGH ENOUGH TO PREVENT NUMERICAL INSTABILITY, LIKELY TO DO WITH SHARP GRADIENTS

//Create variables for params needed from input file
Real ism_dens;
Real ism_temp;
Real ism_vel;
Real radius_swBC;
Real sw_densBC; 
Real sw_tempBC;
Real sw_velBC;
Real r_Inner;

//Create variables for primitive values which will be derived from the input file values
Real rho_ism;
Real thermPr_ism;
Real rho_swBC;
Real thermPr_swBC;
Real pfloor;

//Function Headers
std::vector<Real> SphericalTransform(Real x, Real y, Real z);
Real FindRefinementRegionDxSpacing(MeshBlock *pmb, Coordinates *pco);
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

//Intakes Cartesian coordinates, converts to spherical coordinate parameters needed for primitive, conserved variable transforms
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

//Loop over all cells, find the smallest radial cell spacing value, which corresponds to the cell spacing in the region with the highest refinement level, which is where
//our inner radius boundary condition is located. Need to correctly set radial BC on our Cartesian grid
Real FindRefinementRegionDrSpacing(Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku){
    Real smallest_dr = 1000; //(AU), initialize with an arbitrarily large value
    for (int i=il; i<=iu; i++) {
        Real dx = pco->dx1f(i);
        /*for (int j=jl; j<=ju; j++) {
            Real dy = pco->dx2f(j);
            for (int k=il; k<=ku; k++) {
                Real dz = pco->dx3f(k);
                Real dr = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (dr < smallest_dr){
                    smallest_dr = dr;
                }
            }       
        }
            */
        if (dx < smallest_dr){
            smallest_dr = dx;
        }
    }
    return (smallest_dr);
}


void MeshBlock::ProblemGenerator(ParameterInput *pin){

    //Get physical parameters needed from input file
    ism_dens = pin->GetReal("problem", "ism_dens");
    ism_temp = pin->GetReal("problem", "ism_temp");
    ism_vel = pin->GetReal("problem", "ism_vel");

    radius_swBC = pin->GetReal("problem", "radius_swBC");
    sw_densBC = pin->GetReal("problem", "sw_densBC");
    sw_tempBC = pin->GetReal("problem", "sw_tempBC");
    sw_velBC = pin->GetReal("problem", "sw_velBC");
    r_Inner = pin->GetReal("problem", "r_inner");
    pfloor = pin->GetReal("hydro", "pfloor");

    //Compute prim values needed for initilization
    rho_ism = ism_dens * PROTON_MASS; //kg cm-3
    thermPr_ism = 2 * ism_dens * BOLTZMANN_CONSTANT * ism_temp * 1e6; // (J/m-3 = N/m^2), use ideal gas law P=2nkT (factor of 2 since both protons and electrons), 1e-6 factor for cm-3 to m-3
    rho_swBC = sw_densBC * PROTON_MASS * std::pow(radius_swBC / r_Inner, 2); //kg cm-3, assume 1/r^2 scaling of density for solar wind
    thermPr_swBC = 2 * sw_densBC * BOLTZMANN_CONSTANT * sw_tempBC * 1e6 * std::pow(radius_swBC / r_Inner, 10.0/3.0); //(J/m^3 = N/m^2), 1e-6: conversion factor from cm-3 to m-3. Assume 1/r^(10/3) scaling of therm press (adiabatic)

    //dimensions of meshblock
    const int Nx = ie - is + 1;
    const int Ny = je - js + 1;
    const int Nz = ke - ks + 1; 

    //Set a small dr such that cells with r>r_inner, r<r_inner+dr are assigned boundary values. Set to the value of spacing between cells in x,y,z direction including SMR/AMR refinement.
    Real dr = FindRefinementRegionDrSpacing(pcoord, is, ie, js, je, ks, ke);
    Real dr_withBlendZone = 3 * dr;

    AthenaArray<Real> b; //needed for PrimitiveToConserved()
    b.NewAthenaArray(Nz, Ny, Nx);

    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
                Real x = pcoord->x1v(i);
                Real y = pcoord->x2v(j);
                Real z = pcoord->x3v(k);
                std::vector<Real> spherical = SphericalTransform(x, y, z);

                //Set values less than inner radial boundary to zero; outside simulation domain
                if (spherical[1] < r_Inner){
                    phydro->w(IDN,k,j,i) =  rho_inner;
                    phydro->w(IPR,k,j,i) =  press_inner;
                    phydro->w(IVX,k,j,i) =  0; 
                    phydro->w(IVY,k,j,i) =  0; 
                    phydro->w(IVZ,k,j,i) =  0; 
                }
                //Set radial BC (should be ~1 cell thick radial boundary)
                else if(spherical[1] > r_Inner && spherical[1] < r_Inner + dr){
                    phydro->w(IDN,k,j,i) =  (rho_swBC / rho0);
                    phydro->w(IPR,k,j,i) =  (thermPr_swBC / P0);
                    phydro->w(IVX,k,j,i) =  (sw_velBC / v0) * spherical[4] * spherical[3]; //sin_theta*cos_phi
                    phydro->w(IVY,k,j,i) =  (sw_velBC / v0) * spherical[4] * spherical[2]; //sin_theta*sin_phi
                    phydro->w(IVZ,k,j,i) =  (sw_velBC / v0) * spherical[5]; //cos_theta
                }
                //Apply a smoothing layer, assuming normal expansions of density and pressure (should be ~2-3 cells thick)
                else if(spherical[1] > r_Inner + dr && spherical[1] < r_Inner + dr_withBlendZone){
                    phydro->w(IDN,k,j,i) =  (rho_swBC/ rho0) * std::pow(r_Inner / spherical[1], 2.0);
                    phydro->w(IPR,k,j,i) =  (thermPr_swBC / P0) * std::pow(r_Inner / spherical[1], 10.0/3.0);
                    phydro->w(IVX,k,j,i) =  (sw_velBC / v0) * spherical[4] * spherical[3]; 
                    phydro->w(IVY,k,j,i) =  (sw_velBC / v0) * spherical[4] * spherical[2]; 
                    phydro->w(IVZ,k,j,i) =  (sw_velBC / v0) * spherical[5];
                }
                //Initialize remaining values to ISM values
                else {
                    phydro->w(IDN,k,j,i) =  rho_ism / rho0;
                    phydro->w(IPR,k,j,i) =  thermPr_ism / P0;
                    phydro->w(IVX,k,j,i) =  -ism_vel / v0; //-x direction
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
                                    prim(IVX,k,j,i) = -ism_vel / v0; //x-direction
                                    prim(IVY,k,j,i) = 0;
                                    prim(IVZ,k,j,i) = 0;
                                    prim(IPR,k,j,i) = thermPr_ism / P0;
                                }
                            }
                        }           
}

void MeshBlock::UserWorkInLoop() {
    //dimensions of meshblock
    const int Nx = ie - is + 1;
    const int Ny = je - js + 1;
    const int Nz = ke - ks + 1; 

    Real dr = FindRefinementRegionDrSpacing(pcoord, is, ie, js, je, ks, ke);
    Real dr_withBlendZone = 3 * dr;

    AthenaArray<Real> b; //needed for PrimitiveToConserved()
    b.NewAthenaArray(Nz, Ny, Nx);

  for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
                Real x = pcoord->x1v(i);
                Real y = pcoord->x2v(j);
                Real z = pcoord->x3v(k);
                std::vector<Real> spherical = SphericalTransform(x, y, z);
                if(spherical[1] > r_Inner && spherical[1] < r_Inner + dr){
                    phydro->w(IDN,k,j,i) = (rho_swBC/ rho0);
                    phydro->w(IPR,k,j,i) = (thermPr_swBC / P0);
                    phydro->w(IVX,k,j,i) = (sw_velBC / v0) * spherical[4] * spherical[3]; //sin_theta*cos_phi
                    phydro->w(IVY,k,j,i) = (sw_velBC / v0) * spherical[4] * spherical[2]; //sin_theta*sin_phi
                    phydro->w(IVZ,k,j,i) = (sw_velBC / v0) * spherical[5] ; //cos_theta
                }
                ///Apply radial BC including a smoothing layer, assuming normal expansions of density and pressure (should be ~2-3 cells thick)
                else if(spherical[1] > r_Inner + dr && spherical[1] < r_Inner + dr_withBlendZone){
                    phydro->w(IDN,k,j,i) =  (rho_swBC/ rho0) * std::pow(r_Inner / spherical[1], 2.0);
                    phydro->w(IPR,k,j,i) =  (thermPr_swBC / P0) * std::pow(r_Inner / spherical[1], 10.0/3.0); //(2 * phydro->w(IDN,k,j,i) * BOLTZMANN_CONSTANT * sw_tempBC * 1e6 * std::pow(radius_swBC / spherical[1], 4.0/3.0)) / P0;
                    phydro->w(IVX,k,j,i) =  (sw_velBC / v0) * spherical[4] * spherical[3]; 
                    phydro->w(IVY,k,j,i) =  (sw_velBC / v0) * spherical[4] * spherical[2]; 
                    phydro->w(IVZ,k,j,i) =  (sw_velBC / v0) * spherical[5];
                }
                else if(spherical[1] < r_Inner){
                    phydro->w(IDN,k,j,i) =  rho_inner;
                    phydro->w(IPR,k,j,i) =  press_inner;
                    phydro->w(IVX,k,j,i) =  0; 
                    phydro->w(IVY,k,j,i) =  0; 
                    phydro->w(IVZ,k,j,i) =  0; 
                }
                //If pressure values become negative due to numerical errors, floor values will be enforced. The low temp at the BC will cause failure of pressure expansion without pickup ions.
                //For now, this only happens in the supersonic, kinetically dominated wind, so we will manually compute pressure values within termination shock and replace values according to adiabatic expansion.
                //This is the best way to ensure our simulation expands correctly with an adiabatic pressure/temp profile, the pressure expansion will naturally begin to fail if no values are updated.
                //If we try to replace pressure values only when they are below a certain value, the profile becomes strange and not spherically symmetric.
                Real vel_mag = std::sqrt(std::pow(phydro->w(IVX,k,j,i),2) + std::pow(phydro->w(IVY,k,j,i),2) + std::pow(phydro->w(IVZ,k,j,i),2));
                if ((spherical[1] > (r_Inner + dr_withBlendZone)) && (vel_mag > ((sw_velBC / v0) - 0.2))){
                    phydro->w(IPR,k,j,i) = (thermPr_swBC / P0) * std::pow(r_Inner / spherical[1], 10.0/3.0); //(2 * phydro->w(IDN,k,j,i) * BOLTZMANN_CONSTANT * sw_tempBC * 1e6 * std::pow(radius_swBC / spherical[1], 4.0/3.0)) / P0; //r^-4/3 factor: adiabatic temperature
                }
            }
        }
    }
    peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord, is, ie, js, je, ks, ke);
}

//Refinement based on Temperature gradient - INCOMPLETE
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