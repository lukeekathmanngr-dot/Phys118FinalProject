//========================================================================================
// HydrodynamicHelioSpherical3D.cpp
//========================================================================================
// Problem initialization file for a magnetized stellar wind interacting with a PLASMA ONLY, UNMAGNETIZED inflowing ISM
// This serves to set the primitive/conserved values across the entire simulation grid at t=0
// The simulation has an inner radial boundary, set in the in file (in AU) by the user
// The inflowing ISM is set as a boundary condition on the outer x1 (x coordinate in Cartesian), inflowing in the -x direction
// ASSUME AN ISM OF ONLY PROTONS, NO CHARGE EXCHANGE
// Written by Luke Kathmann, 2025

// C++ Headers
#include <cmath>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

// Ensure magnetic field are enabled
#if !MAGNETIC_FIELDS_ENABLED
#error "This code requires magnetic fields"
#endif

//================================================================================================================================================
// Normalization and Physical Constants as outlined in Hans document: https://drive.google.com/drive/u/2/folders/1RPFssIHfPF7QZJYQLEVFdbCg5tahYjXA
//================================================================================================================================================
const double r0 = 1; //AU
const double v0 = 52.483; //km/s
const double t0 = 2.851e6; //s = 33 days
const double rho0 = 8.361e-27; //kg cm-3
const double P0 = 2.3e-11; //N/m^2
const double e0 = P0;
const double BOLTZMANN_CONSTANT = 1.38e-23; // J K-1
const double PROTON_MASS = 1.6726e-27; // kg
const double b0 = 5.38e-9; //Tesla
const double rho_inner = 1e-6; //MUST BE HIGH ENOUGH TO PREVENT NUMERICAL INSTABILITY, LIKELY TO DO WITH SHARP GRADIENTS (floor value for inactive simulation region)
const double press_inner= 1e-8; //MUST BE HIGH ENOUGH TO PREVENT NUMERICAL INSTABILITY, LIKELY TO DO WITH SHARP GRADIENTS (floor value for inactive simulation region)
const double b_inner = 0; //pow factor for nT to T

//========================================================================================
// Varaiables needed for initialization and boundary conditions
//========================================================================================
Real ism_dens;
Real ism_temp;
Real ism_vel;
Real ism_B;
Real radius_swBC;
Real sw_densBC; 
Real sw_tempBC;
Real sw_velBC;
Real sw_bBC;
Real r_Inner;
Real angle_ism;
Real omega;
Real dFloor;
Real pFloor;
Real bCeiling;
Real bFloor;

//Create variables for primitive values (thermal pressure, mass density) which will be derived from the input file values
Real rho_ism;
Real thermPr_ism;
Real rho_swBC;
Real thermPr_swBC;

//========================================================================================
// Function Headers
//========================================================================================

std::vector<Real> ComputeParkerSpiralField(Real x, Real y, Real z);

void Boundary_ox1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void Boundary_ix1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh);
                    
int RefinementCondition(MeshBlock *pmb);

//Set our user defined BC for inflowing ISM. Set up adapative mesh refinement condition if enabled in input file
void Mesh::InitUserMeshData(ParameterInput *pin){
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, Boundary_ox1);
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, Boundary_ix1);
    //if (adaptive==true){
        //EnrollUserRefinementCondition(RefinementCondition);
    //}
}

//=====================================================================================================================================================
// Function: ComputeParkerSpiral
// Purpose: Take in Cartesian coordinates and find cooresponding Parker Spiral magnetic field as derived in Parker's 1958 paper
// Output: 3 component vector of the corresponding normalized Parker Spiral magnetic field for a point x,y,z transformed back to Cartesian coordinates
//=====================================================================================================================================================
std::vector<Real> ComputeParkerSpiralField(Real r, Real theta, Real phi){
    std::vector<Real> bField(3);

    Real bR = (sw_bBC / b0) * std::pow(radius_swBC / r, 2); //sphericalTrans[1]: spherical radius;
    Real bTheta = 0; //Always zero according to Parker derivation
    Real bPhi = -(sw_bBC / b0) * ((radius_swBC * radius_swBC) / (sw_velBC / v0)) * (1 / r) * (omega * t0) * std::sin(theta);
    
    bField[0] = bR; 
    bField[1] = 0; 
    bField[2] = bPhi;
    return bField;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin){
    //Get physical parameters needed from input file
    ism_dens = pin->GetReal("problem", "ism_dens");
    ism_temp = pin->GetReal("problem", "ism_temp");
    ism_vel = pin->GetReal("problem", "ism_vel");
    ism_B = pin->GetReal("problem", "ism_B");
    angle_ism = pin->GetReal("problem", "angle_ism");
    omega = pin->GetReal("problem", "omega");
    dFloor = pin->GetReal("problem", "dfloor");
    pFloor = pin->GetReal("problem", "pfloor");
    bCeiling = pin->GetReal("problem", "bCeiling");
    bFloor = pin->GetReal("problem", "bFloor");

    radius_swBC = pin->GetReal("problem", "radius_swBC");
    sw_densBC = pin->GetReal("problem", "sw_densBC");
    sw_tempBC = pin->GetReal("problem", "sw_tempBC");
    sw_velBC = pin->GetReal("problem", "sw_velBC");
    sw_bBC = pin->GetReal("problem", "sw_bBC");
    r_Inner = pin->GetReal("problem", "r_inner");
    //Compute prim values needed for initilization
    rho_ism = ism_dens * PROTON_MASS; //kg cm-3
    thermPr_ism = 2 * ism_dens * BOLTZMANN_CONSTANT * ism_temp * 1e6; // (J/m-3 = N/m^2), use ideal gas law P=2nkT (factor of 2 since both protons and electrons), 1e-6 factor for cm-3 to m-3
    rho_swBC = sw_densBC * PROTON_MASS * std::pow(radius_swBC / r_Inner, 2); //kg cm-3, assume 1/r^2 scaling of density for solar wind
    thermPr_swBC = 2 * sw_densBC * BOLTZMANN_CONSTANT * sw_tempBC * 1e6 * std::pow(radius_swBC / r_Inner, 10.0/3.0); //(J/m^3 = N/m^2), 1e-6: conversion factor from cm-3 to m-3. Assume 1/r^(10/3) scaling of therm press (adiabatic)
    sw_bBC = sw_bBC * pow(10,-9); //(nT to T)
    ism_B = ism_B * pow(10,-9); //(nT to T)
    bCeiling = bCeiling * pow(10,-9); //(nT to T)

    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
                Real r = pcoord->x1v(i);
                Real theta = pcoord->x2v(j);
                Real phi = pcoord->x3v(k);
                phydro->w(IDN,k,j,i) =  ism_dens / rho0;
                phydro->w(IPR,k,j,i) =  thermPr_ism / P0;
                phydro->w(IVX,k,j,i) =  (-ism_vel / v0) * std::sin(theta) * std::cos(phi); 
                phydro->w(IVY,k,j,i) =  (-ism_vel / v0) * std::cos(theta) * std::cos(phi); 
                phydro->w(IVZ,k,j,i) =  (-ism_vel / v0) * std::sin(phi); 
            }   
        }
    }
    
    //Set Interface Bx
    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie+1; i++) {
                Real rFace = pcoord->x1f(i);
                Real theta = pcoord->x2v(j);
                Real phi = pcoord->x3v(k);
                pfield->b.x1f(k,j,i) = 0; //Start with no magnetic field in ISM
            }
        }
    }
    //Set Interface By
    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
            for (int i=is; i<=ie; i++) {
                Real r = pcoord->x1v(i);
                Real thetaFace = pcoord->x2f(j);
                Real phi = pcoord->x3v(k);
                pfield->b.x2f(k,j,i) = 0; //Start with no magnetic field in ISM
            }
        }
    }
    //Set Interface Bz
    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
            for (int i=is; i<=ie; i++) {
                Real r = pcoord->x1v(i);
                Real theta = pcoord->x2v(j);
                Real phiFace = pcoord->x3f(k);
                pfield->b.x3f(k,j,i) = 0; //Start with no magnetic field in ISM
            }
        }
    }
    pfield->CalculateCellCenteredField(pfield->b, pfield->bcc, pcoord, is, ie, js, je, ks, ke);
    peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
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
                                    Real theta = pco->x2v(j);
                                    Real phi = pco->x3v(k);
                                    prim(IDN,k,j,i) = rho_ism / rho0;
                                    prim(IPR,k,j,i) = thermPr_ism / P0;
                                    prim(IVX,k,j,i) = -(ism_vel / v0) * std::sin(theta) * std::cos(phi); //x-direction
                                    prim(IVY,k,j,i) = -(ism_vel / v0) * std::cos(theta) * std::cos(phi);
                                    prim(IVZ,k,j,i) = -(ism_vel / v0) * (-std::sin(phi));
                                }
                            }
                        }
                        //Set x1 Bfield
                        for (int k=kl; k<=ku; ++k){
                            for (int j=jl; j<=ju; ++j){
                                for (int i=iu+2; i<=iu+ngh+1; ++i){
                                    Real theta = pco->x2v(j);
                                    Real phi = pco->x3v(k);
                                    b.x1f(k,j,i) = 0;
                                }
                            }
                        }
                        //Set x2 Bfield
                        for (int k=kl; k<=ku; ++k){
                            for (int j=jl; j<=ju+1; ++j){
                                for (int i=iu+1; i<=iu+ngh; ++i){
                                    Real theta = pco->x2v(j);
                                    Real phi = pco->x3v(k);
                                    b.x2f(k,j,i) = 0;
                                }
                            }
                        }
                        //Set x3 Bfield
                        for (int k=kl; k<=ku+1; ++k){
                            for (int j=jl; j<=ju; ++j){
                                for (int i=iu+1; i<=iu+ngh; ++i){
                                    Real theta = pco->x2v(j);
                                    Real phi = pco->x3v(k);
                                    b.x3f(k,j,i) = 0;
                                }
                            }
                        }       
}

//Define our inner radial solar wind BC
//Note il= is, iu=ie, but need to use the ones defined in this function due to SMR
void Boundary_ix1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh){
                        for (int k=kl; k<=ku; ++k){
                            for (int j=jl; j<=ju; ++j){
                                //iu+1: first ghost cell beyond last active z1 cell. iu+ngh: last ghost cell
                                //Set primitive variables
                                for (int i=il-ngh; i<=il-1; ++i){
                                    Real r = pco->x1v(i);
                                    prim(IDN,k,j,i) = (rho_swBC/ rho0) * std::pow(r_Inner / r, 2.0);
                                    prim(IPR,k,j,i) = (thermPr_swBC / P0) * std::pow(r_Inner / r, 10.0/3.0);
                                    prim(IVX,k,j,i) = (sw_velBC / v0);
                                    prim(IVY,k,j,i) = 0;
                                    prim(IVZ,k,j,i) = 0;
                                }
                            }
                        }
                        //Set x1 Bfield
                        for (int k=kl; k<=ku; ++k){
                            for (int j=jl; j<=ju; ++j){
                                for (int i=il-ngh; i<=il-1; ++i){
                                    Real rFace = pco->x1f(i);
                                    Real theta = pco->x2v(j);
                                    Real phi = pco->x3v(k);
                                    std::vector<Real> parker = ComputeParkerSpiralField(rFace, theta, phi);
                                    b.x1f(k,j,i) = parker[0];
                                }
                            }
                        }
                        //Set x2 Bfield
                        for (int k=kl; k<=ku; ++k){
                            for (int j=jl; j<=ju+1; ++j){
                                for (int i=il-ngh; i<=il-1; ++i){
                                    Real r = pco->x1v(i);
                                    Real thetaFace = pco->x2f(j);
                                    Real phi = pco->x3v(k);
                                    std::vector<Real> parker = ComputeParkerSpiralField(r, thetaFace, phi);
                                    b.x2f(k,j,i) = parker[1];
                                }
                            }
                        }
                        //Set x3 Bfield
                        for (int k=kl; k<=ku+1; ++k){
                            for (int j=jl; j<=ju; ++j){
                                for (int i=il-ngh; i<=il-1; ++i){
                                    Real r = pco->x1v(i);
                                    Real theta = pco->x2v(j);
                                    Real phiFace = pco->x3f(k);
                                    std::vector<Real> parker = ComputeParkerSpiralField(r, theta, phiFace);
                                    b.x2f(k,j,i) = parker[2];
                                }
                            }
                        }
}