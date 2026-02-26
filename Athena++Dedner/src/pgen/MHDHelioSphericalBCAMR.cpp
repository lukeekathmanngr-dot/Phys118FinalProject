//========================================================================================
// HydrodynamicHelioCarteisan3DMagnetizedWind.cpp
//========================================================================================
// Problem initialization file for a magnetized stellar wind interacting with a PLASMA ONLY, UNMAGNETIZED inflowing ISM
// This serves to set the primitive/conserved values across the entire simulation grid at t=0
// The simulation has an inner radial boundary, set in the in file (in AU) by the user
// The inflowing ISM is set as a boundary condition on the outer x1 (x coordinate in Cartesian), inflowing in the -x direction
// ASSUME AN ISM OF ONLY PROTONS, NO CHARGE EXCHANGE
// Written by Luke Kathmann, 2025

// NOTE: Source code was modified to implement Dedner cleaning. Description of Dedner cleaning and files that were changed can be found here:
// PASTE DEDNER DOCUMENT HERE!

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
const double b0 = 5.38e-9; //Tesla
const double BOLTZMANN_CONSTANT = 1.38e-23; // J K-1
const double PROTON_MASS = 1.6726e-27; // kg
const double rho_inner = 1e-6; //MUST BE HIGH ENOUGH TO PREVENT NUMERICAL INSTABILITY, LIKELY TO DO WITH SHARP GRADIENTS (floor value for inactive simulation region)
const double press_inner= 1e-8; //MUST BE HIGH ENOUGH TO PREVENT NUMERICAL INSTABILITY, LIKELY TO DO WITH SHARP GRADIENTS (floor value for inactive simulation region)
const double b_inner = 0;

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
Real ch_glm;
Real cp_glm;

//Create variables for primitive values (thermal pressure, mass density) which will be derived from the input file temperatures, densities, etc
Real rho_ism;
Real thermPr_ism;
Real rho_swBC;
Real thermPr_swBC;

//Create variables for AMR refinement
Real threshold_densGrad;
Real deref_threshold_densGrad;
Real threshold_tempGrad;
Real deref_threshold_tempGrad;
Real threshold_pressGrad;
Real deref_threshold_pressGrad;

//========================================================================================
// Function Headers
//========================================================================================
std::vector<Real> SphericalTransform(Real x, Real y, Real z);

std::vector<Real> ComputeParkerSpiralField(Real x, Real y, Real z);

Real FindRefinementRegionSpacing(Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku);

void Boundary_ox1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh);

int RefinementCondition(MeshBlock *pmb);

//Set our user defined BC for inflowing ISM. Set up adapative mesh refinement condition if enabled in input file
void Mesh::InitUserMeshData(ParameterInput *pin){
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, Boundary_ox1);
    if (adaptive==true){
        EnrollUserRefinementCondition(RefinementCondition);
    }
}

//==================================================================================================
// Function: SphericalTransform
// Purpose: Take in Cartesian coordinates and find cooresponding spherical transformation
// Output: 6 Component vector of variables needed for Cartesian -> spherical (rho, r, sintheta, etc)
//=================================================================================================
std::vector<Real> SphericalTransform(Real x, Real y, Real z){
    std::vector<Real> transformVec(6);
    transformVec[0] = std::sqrt(x * x + y * y); //Polar Radius
    transformVec[1] = std::sqrt(x * x + y * y + z * z); //Spherical Radius
    transformVec[2] = y / transformVec[0]; //Sin_phi, x/r
    transformVec[3] = x / transformVec[0]; //Cos_phi, y/r
    transformVec[4] = transformVec[0] / transformVec[1]; //Sin_theta, r/rho
    transformVec[5] = z / transformVec[1]; //Cos_theta, z/rho
    return transformVec;
}

//=====================================================================================================================================================
// Function: ComputeParkerSpiral
// Purpose: Take in Cartesian coordinates and find cooresponding Parker Spiral magnetic field as derived in Parker's 1958 paper
// Output: 3 component vector of the corresponding normalized Parker Spiral magnetic field for a point x,y,z transformed back to Cartesian coordinates
//=====================================================================================================================================================
std::vector<Real> ComputeParkerSpiralField(Real x, Real y, Real z){
    std::vector<Real> bField(3);
    std::vector<Real> sphericalTrans = SphericalTransform(x, y, z);
    Real bR;
    Real bTheta;
    Real bPhi;

    bR = (sw_bBC / b0) * std::pow(radius_swBC / sphericalTrans[1], 2); //sphericalTrans[1]: spherical radius
    bTheta = 0; //Always zero according to Parker derviation
    bPhi = -(sw_bBC / b0) * ((radius_swBC * radius_swBC) / (sw_velBC / v0)) * (1 / sphericalTrans[1]) * (omega * t0) * sphericalTrans[4];
    
    //Convert the solution in Spherical coordinates back to Cartesian
    bField[0] = bR * sphericalTrans[4] * sphericalTrans[3]  - bPhi * sphericalTrans[2]; //x-component, br * sintheta cosphi - bphi *sinphi
    bField[1] = bR * sphericalTrans[4] * sphericalTrans[2]  + bPhi * sphericalTrans[3]; //y-component, br * sintheta sinphi - bphi *cosphi
    bField[2] = bR * sphericalTrans[5]; //z-component, br * costheta

    return bField;
}

//========================================================================================================================================================
// Function: FindRefinementRegionSpacing
// Purpose: Find the smallest cell spacing distance on the entire mesh
// Output: Return the smallest cell spacing distance on our mesh, which corresponds to the inner refinement region. Needed in order to 
//         properly set the spherical boundary condition (overwrite SW values between inner radius set in input file and r_inner + integer * smallest_dr)
//========================================================================================================================================================
Real FindRefinementRegionSpacing(Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku){
    Real smallest_dr = 1000; //(AU), initialize with an arbitrarily large value
    for (int i=il; i<=iu; i++) { //Loop over all x1 in the meshblock (il: lower/starting index. ie: ending index)
        Real dx = pco->dx1v(i);
        if (dx < smallest_dr){
            smallest_dr = dx;
        }
    }
    for (int j=jl; j<=ju; j++) { //Loop over all x2
        Real dy = pco->dx2v(j);
        if (dy < smallest_dr){
            smallest_dr = dy;
        }
    }
    for (int k=kl; k<=ku; k++) { //Loop over all x3
        Real dz = pco->dx3v(k);
        if (dz< smallest_dr){
            smallest_dr = dz;
        }
    }
    return (smallest_dr);
}

//========================================================================================================================================================
// Function: ProblemGenerator
// Purpose: Initialize mesh with SW boundary conditions, rest of mesh inititalized to ambient ISM values. Called automatically at beginning of simulation
// Output: None
//========================================================================================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin){

    //Get physical parameters needed from input file
    ism_dens = pin->GetReal("problem", "ism_dens");
    ism_temp = pin->GetReal("problem", "ism_temp");
    ism_vel = pin->GetReal("problem", "ism_vel");
    ism_B = pin->GetReal("problem", "ism_B");
    angle_ism = pin->GetReal("problem", "angle_ism");
    omega = pin->GetReal("problem", "omega");
    radius_swBC = pin->GetReal("problem", "radius_swBC");
    sw_densBC = pin->GetReal("problem", "sw_densBC");
    sw_tempBC = pin->GetReal("problem", "sw_tempBC");
    sw_velBC = pin->GetReal("problem", "sw_velBC");
    sw_bBC = pin->GetReal("problem", "sw_bBC");
    r_Inner = pin->GetReal("problem", "r_inner");
    phydro->hsrc.flag_dedner_source = pin->GetOrAddBoolean("hydro", "flag_dedner_source", true); //TURNS ON DEDNER PSI FIELD SOURCE TERM/EVOLUTION
    cp_glm = pin->GetReal("hydro", "cp_glm"); // Parabolic Damping coefficient, for Dedner
    ch_glm = pin->GetReal("hydro", "ch_glm"); // Hyperbolic Wave coefficient, for Dedner

    //Get thresholds for refinement from input file
    Real threshold_densGrad = pin->GetReal("problem", "threshold_densGrad");
    Real deref_threshold_densGrad = pin->GetReal("problem", "deref_threshold_densGrad");
    Real threshold_TGrad = pin->GetReal("problem", "threshold_tempGrad");
    Real deref_threshold_TGrad = pin->GetReal("problem", "deref_threshold_tempGrad");
    Real threshold_pressGrad = pin->GetReal("problem", "threshold_pressGrad");
    Real deref_threshold_pressGrad = pin->GetReal("problem", "deref_threshold_pressGrad");

    //Compute prim values needed for initilization
    rho_ism = ism_dens * PROTON_MASS; //kg cm-3
    thermPr_ism = 2 * ism_dens * BOLTZMANN_CONSTANT * ism_temp * 1e6; // (J/m-3 = N/m^2), use ideal gas law P=2nkT (factor of 2 since both protons and electrons), 1e-6 factor for cm-3 to m-3
    rho_swBC = sw_densBC * PROTON_MASS * std::pow(radius_swBC / r_Inner, 2); //kg cm-3, assume 1/r^2 scaling of density for solar wind
    thermPr_swBC = 2 * sw_densBC * BOLTZMANN_CONSTANT * sw_tempBC * 1e6 * std::pow(radius_swBC / r_Inner, 10.0/3.0); //(J/m^3 = N/m^2), 1e-6: conversion factor from cm-3 to m-3. Assume 1/r^(10/3) scaling of therm press (adiabatic)
    sw_bBC = sw_bBC * pow(10,-9); //(nT to T)
    ism_B = ism_B * pow(10,-9); //(nT to T)

    Real dr = FindRefinementRegionSpacing(pcoord, is, ie, js, je, ks, ke); //smallest cell spacing on mesh
    Real r_outer = r_Inner + 1*dr; // Sets spherical overwrite region: overwrite values with SW values between r_Inner and r_Inner + r_outer
    Real r_blend = r_outer + 3*dr; // Sets blend region: overwrite values with blended expected SW values and values from the computational domain, acts to smooth solution initially
    Real dr_blend = 3*dr; // Thickness of blend region

    //Set the hydro class public ch and cp to initial input to be used by code (should already be set, but set again just to be safe)
    phydro->ch_glm = ch_glm;
    phydro->cp_glm = cp_glm;

    //Loop over current mesh block, initialize hydro variables
    for (int k=ks; k<=ke; k++) { //ks: lower x3 bound. ke: upper x3 bound.
        for (int j=js; j<=je; j++) { //js: lower x2 bound. je: upper x2 bound.
            for (int i=is; i<=ie; i++) { //is: lower x1 bound. ie: upper x1 bound.

                //Get coordinates at corresponding to indices location on meshblock
                Real x = pcoord->x1v(i);
                Real y = pcoord->x2v(j);
                Real z = pcoord->x3v(k);
                //Compute spherical parameters 
                std::vector<Real> spherical = SphericalTransform(x, y, z);

                //Initialize scalar field psi used for Dedner divergence cleaning to zero everywhere on mesh
                //Note: phydro->u corresponds to array of conserved quantities used in flux calculations (momentum, energy, etc), phydro->w corresponds to primitive values (temp, velocity, etc)
                phydro->u(IPSIC,k,j,i) = 0;
                phydro->w(IPSIW,k,j,i) = 0;

                //If radius (spherical[1]) of point on mesh is less than r_Inner, set hydro variables to zero/floor values (outside of computational domain)
                if (spherical[1] < r_Inner){
                    phydro->w(IDN,k,j,i) =  rho_inner;
                    phydro->w(IPR,k,j,i) =  press_inner;
                    phydro->w(IVX,k,j,i) =  0; 
                    phydro->w(IVY,k,j,i) =  0; 
                    phydro->w(IVZ,k,j,i) =  0; 
                }

                //If radius (spherical[1]) of point on mesh is greater than r_Inner and less than r_outer, this corresponds to our spherical shell BC overwrite region. SET SW VALUES HERE!!!
                else if(spherical[1] > r_Inner && spherical[1] < r_outer){
                    phydro->w(IDN,k,j,i) =  (rho_swBC/ rho0) * std::pow(r_Inner / spherical[1], 2.0); //Assume 1/r^2 expansion of density
                    phydro->w(IPR,k,j,i) =  (thermPr_swBC / P0) * std::pow(r_Inner / spherical[1], 10.0/3.0); //Assume 1/r^10/3 expansion of pressure corresponding to adiabtic expansion
                    //assume constant radial velocity, need to transform v_r rhat = v_r * (sintheta*cosphi xhat + sintheta * sinphi yhat + costheta zhat)
                    phydro->w(IVX,k,j,i) =  (sw_velBC / v0) * spherical[4] * spherical[3]; // sintheta * cosphi
                    phydro->w(IVY,k,j,i) =  (sw_velBC / v0) * spherical[4] * spherical[2]; //sintheta * sinphi
                    phydro->w(IVZ,k,j,i) =  (sw_velBC / v0) * spherical[5]; //costheta
                }

                //If radius (spherical[1]) of point on mesh is greater than r_outer and less than r_blend, this cirresponds to spherical blend region. SET VALUES WITH MIX OF SW AND ISM VALUES, 
                //WEIGHTED LINEARLY BASED ON PROXIMITY TO SPHERICAL BOUNDARY OVERWRITE REGION. This blend region will act to smooth sharp gradients initially.
                else if(spherical[1] > r_outer && spherical[1] < r_blend){
                    Real alpha = (spherical[1] - r_outer) / dr_blend; //Linear weight
                    //Overwrite with same assumption for SW values as before, multiply by weight. Add on ISM values, multipled by (1-weight) to smooth initial gradients from SW to ISM values
                    phydro->w(IDN,k,j,i) =  (1-alpha) * (rho_swBC/ rho0) * std::pow(r_Inner / spherical[1], 2.0) + alpha * (rho_ism / rho0);
                    phydro->w(IPR,k,j,i) =  (1-alpha) * (thermPr_swBC / P0) * std::pow(r_Inner / spherical[1], 10.0/3.0) + alpha * (thermPr_ism / P0);
                    phydro->w(IVX,k,j,i) =  (1-alpha) * (sw_velBC / v0) * spherical[4] * spherical[3] - (alpha) * ism_vel / v0; 
                    phydro->w(IVY,k,j,i) =  (1-alpha) * (sw_velBC / v0) * spherical[4] * spherical[2]; 
                    phydro->w(IVZ,k,j,i) =  (1-alpha) * (sw_velBC / v0) * spherical[5];

                }

                //Any points outside this radius, assign to ISM values
                else {
                    phydro->w(IDN,k,j,i) =  rho_ism / rho0;
                    phydro->w(IPR,k,j,i) =  thermPr_ism / P0;
                    phydro->w(IVX,k,j,i) =  -ism_vel / v0; //-x inflow
                    phydro->w(IVY,k,j,i) =  0; 
                    phydro->w(IVZ,k,j,i) =  0;
                }
            }   
        }
    }

    //Set Interface Bx
    //IMPORTANT NOTE: fields used by Athena++ are face centered. To properly assign all magnetic field values, need to loop over first over x1 running from is to ie+1 (face centered means the upper bound for B is extended by 1 vs volume centered hydro variables).
    //Need to use face centered x1 coordinates for first loop, volume centered for x2 and x3 to ensure x1 fields are set at the correct positions.
    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie+1; i++) {
                //Get coordinates corresponding to face centered positions in the x-direction
                Real xFace = pcoord->x1f(i);
                Real y = pcoord->x2v(j);
                Real z = pcoord->x3v(k);
                //Compute spherical parameters corresponding to these positions
                std::vector<Real> spherical = SphericalTransform(xFace, y, z);

                //Radius (spherical[1]) < r_Inner: outside computational domain, set to zero
                if (spherical[1] < r_Inner){
                    pfield->b.x1f(k,j,i) = 0;
                }

                //Radius (spherical[1]) > r_Inner and < r_outer: corresponds to spherical BC. SET TO PARKER SPIRAL SW VALUES. Parker[0] corresponds to x component of Parker spiral
                else if (spherical[1] < r_outer && spherical[1] > r_Inner){
                    std::vector<Real> parker = ComputeParkerSpiralField(xFace, y, z);
                    pfield->b.x1f(k,j,i) = parker[0];
                }

                //Radius (spherical[1]) > r_outer and < r_blend: corresponds to spherical blend region. SET TO MIX OF PARKER AND ISM VALUES, WITH CORRESPONDING LINEAR WEIGHTS TO SMOOTH INITIAL SHARP GRADIENTS
                else if (spherical[1] < r_blend && spherical[1] > r_outer){
                    Real alpha = (spherical[1] - r_outer) / dr_blend;
                    std::vector<Real> parker = ComputeParkerSpiralField(xFace, y, z);
                    pfield->b.x1f(k,j,i) = (1-alpha) * parker[0] + alpha * (ism_B / b0) * std::cos(angle_ism);
                }

                //Any points outside: set to ISM Bx
                else{
                    pfield->b.x1f(k,j,i) = (ism_B / b0) * std::cos(angle_ism);
                }
            }
        }
    }

    //Set Interface By
    //Exact same logic as Bx, except now we index the x2 variable from js to je+1 to properly set face centered x2 B values. Code in loop stays the same, except access face centered x2 positions, cell centered
    //x1 and x3 positions, set x2 B instead of x1 B, and use Parker[1] which corresponds to y component of Parker Spiral.
    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
            for (int i=is; i<=ie; i++) {
                Real x = pcoord->x1v(i);
                Real yFace = pcoord->x2f(j);
                Real z = pcoord->x3v(k);

                std::vector<Real> spherical = SphericalTransform(x, yFace, z);

                if (spherical[1] < r_Inner){
                    pfield->b.x2f(k,j,i) = 0;
                }

                else if (spherical[1] < r_outer && spherical[1] > r_Inner){
                    std::vector<Real> parker = ComputeParkerSpiralField(x, yFace, z);
                    pfield->b.x2f(k,j,i) = parker[1];
                }

                else if (spherical[1] < r_blend && spherical[1] > r_outer){
                    Real alpha = (spherical[1] - r_outer) / dr_blend;
                    std::vector<Real> parker = ComputeParkerSpiralField(x, yFace, z);
                    pfield->b.x2f(k,j,i) = (1-alpha) * parker[1] + alpha * (ism_B / b0) * std::sin(angle_ism);
                }

                else{
                    pfield->b.x2f(k,j,i) = (ism_B / b0) * std::sin(angle_ism);
                }
            }
        }
    }

    //Set Interface Bz
    //Exact same logic as Bx and By, except now we index the x3 variable from ks to ke+1 to properly set face centered x2 B values. Code in loop stays the same, except access face centered x3 positions, cell centered
    //x2 and x3 positions, set x3 B instead of x1/x2 B, and use Parker[2] which corresponds to z component of Parker Spiral.
    for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
                Real x = pcoord->x1v(i);
                Real y = pcoord->x2v(j);
                Real zFace = pcoord->x3f(k);

                std::vector<Real> spherical = SphericalTransform(x, y, zFace);

                if (spherical[1] < r_Inner){
                    pfield->b.x3f(k,j,i) = 0;
                }

                else if (spherical[1] < r_outer && spherical[1] > r_Inner){
                    std::vector<Real> parker = ComputeParkerSpiralField(x, y, zFace);
                    pfield->b.x3f(k,j,i) = parker[2];
                }

                else if (spherical[1] < r_blend && spherical[1] > r_outer){
                    Real alpha = (spherical[1] - r_outer) / dr_blend;
                    std::vector<Real> parker = ComputeParkerSpiralField(x, y, zFace);
                    pfield->b.x3f(k,j,i) = (1-alpha) * parker[2];
                }

                else{
                    pfield->b.x3f(k,j,i) = 0;
                }
            }
        }
    }
    //IMPORTANT: after modifying face centered fields, need to call CalculateCellCenteredField to properly set cell centered fields needed by Athena to calculate magnetic pressure, etc needed for primitive/conserved quantity updates.
    //ALSO: we set primitive variables, need to call PrimitiveToConserved to properly set the conserved quantities used in flux calculations/code evolution
    pfield->CalculateCellCenteredField(pfield->b, pfield->bcc, pcoord, is, ie, js, je, ks, ke); //Face centered B, Cell centered B, meshblock coords, lower/upper cell centered meshblock bounds
    peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke); //Primitive variables, cell centered fields, conserved variables, meshblock coords, lower/upper cell centered meshblock bounds
}

//========================================================================================================================================================
// Function: Boundary_ox1
// Purpose: Set inflow of ISM boundary condition, gets enrolled near top of code (void Mesh::InitUserMeshData(ParameterInput *pin))
// Output: None
//========================================================================================================================================================
void Boundary_ox1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh){
                        for (int k=kl; k<=ku; ++k){
                            for (int j=jl; j<=ju; ++j){
                                //iu+1: first ghost cell beyond last active z1 cell. iu+ngh: last ghost cell. Loop limits for the boundary cells that should be filled can be found here:
                                //https://github.com/PrincetonUniversity/athena/wiki/Boundary-Conditions
                                //Set primitive variables
                                for (int i=iu+1; i<=iu+ngh; ++i){
                                    prim(IDN,k,j,i) = rho_ism / rho0;
                                    prim(IPR,k,j,i) = thermPr_ism / P0;
                                    prim(IVX,k,j,i) = -ism_vel / v0; //x-direction
                                    prim(IVY,k,j,i) = 0;
                                    prim(IVZ,k,j,i) = 0;
                                    prim(IPSIW,k,j,i) = 0;
                                }
                            }
                        }
                        //Set x1 Bfield. Use bounds corresponding to the link provided above
                        for (int k=kl; k<=ku; ++k){
                            for (int j=jl; j<=ju; ++j){
                                for (int i=iu+2; i<=iu+ngh+1; ++i){
                                    b.x1f(k,j,i) = (ism_B / b0) * std::cos(angle_ism);
                                }
                            }
                        }
                        //Set x2 Bfield. Use bounds corresponding to the link provided above
                        for (int k=kl; k<=ku; ++k){
                            for (int j=jl; j<=ju+1; ++j){
                                for (int i=iu+1; i<=iu+ngh; ++i){
                                    b.x2f(k,j,i) = (ism_B / b0) * std::sin(angle_ism);
                                }
                            }
                        }
                        //Set x3 Bfield. Use bounds corresponding to the link provided above
                        for (int k=kl; k<=ku+1; ++k){
                            for (int j=jl; j<=ju; ++j){
                                for (int i=iu+1; i<=iu+ngh; ++i){
                                    b.x3f(k,j,i) = 0;
                                }
                            }
                        }    
}

//========================================================================================================================================================
// Function: UserWorkInLoop
// Purpose: This function is called at the end of Athena++ task list for a given timestep before moving on to next timestep and restarting task list.
//          Use to overwrite hydro/magnetic fields values correctly to set spherical boundary conditions on our Cartesian grid.
// Output: None
//========================================================================================================================================================
void MeshBlock::UserWorkInLoop() {

    Real dr = FindRefinementRegionSpacing(pcoord, is, ie, js, je, ks, ke); //smallest cell spacing on mesh
    Real r_outer = r_Inner + 1*dr; // Sets spherical overwrite region: overwrite values with SW values between r_Inner and r_Inner + r_outer
    Real r_blend = r_outer + 3*dr; // Sets blend region: overwrite values with blended expected SW values and values from the computational domain, acts to smooth solution initially
    Real dr_blend = 3*dr; // Thickness of blend region

    //Loop over current mesh block, Overwrite hydro variables
    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
                //Cell centered coords corresponding to meshblock indices
                Real x = pcoord->x1v(i);
                Real y = pcoord->x2v(j);
                Real z = pcoord->x3v(k);
                //Compute spherical parameters
                std::vector<Real> spherical = SphericalTransform(x, y, z);
                
                //if radius (spherical[1]) < r_Inner, set to zero/floor values (outside computational domain)
                if(spherical[1] < r_Inner){
                    phydro->w(IDN,k,j,i) = rho_inner;
                    phydro->w(IPR,k,j,i) = press_inner;
                    phydro->w(IVX,k,j,i) = 0; 
                    phydro->w(IVY,k,j,i) = 0; 
                    phydro->w(IVZ,k,j,i) = 0; 
                }

                //if radius (spherical[1]) > r_Inner and < r_outer, this corresponds to the spherical BC region. OVERWRITE WITH SW VALUES
                else if(spherical[1] > r_Inner && spherical[1] < r_outer){
                    phydro->w(IDN,k,j,i) = (rho_swBC/rho0) * std::pow(r_Inner/spherical[1], 2.0);
                    phydro->w(IPR,k,j,i) = (thermPr_swBC/P0) * std::pow(r_Inner/spherical[1], 10.0/3.0);
                    phydro->w(IVX,k,j,i) = (sw_velBC/v0) * spherical[4] * spherical[3]; 
                    phydro->w(IVY,k,j,i) = (sw_velBC/v0) * spherical[4] * spherical[2]; 
                    phydro->w(IVZ,k,j,i) = (sw_velBC/v0) * spherical[5];
                }

                //if radius (spherical[1]) > r_outer and < r_blend, this corresponds to the spherical blend region. OVERWRITE WITH LINEARLY WEIGHTED SW AND CURRENT SIMULATION VALUES TO SMOOTH INITIAL SHARP GRADIENTS
                else if(spherical[1] > r_outer && spherical[1] < r_blend){
                    Real alpha = (spherical[1] - r_outer) / dr_blend;
                    phydro->w(IDN,k,j,i) = (1-alpha) * (rho_swBC/rho0) * std::pow(r_Inner/spherical[1], 2.0) 
                                          + alpha * phydro->w(IDN,k,j,i);
                    phydro->w(IPR,k,j,i) = (1-alpha) * (thermPr_swBC/P0) * std::pow(r_Inner/spherical[1], 10.0/3.0) 
                                          + alpha * phydro->w(IPR,k,j,i);
                    phydro->w(IVX,k,j,i) = (1-alpha) * (sw_velBC/v0) * spherical[4] * spherical[3] 
                                          + alpha * phydro->w(IVX,k,j,i); 
                    phydro->w(IVY,k,j,i) = (1-alpha) * (sw_velBC/v0) * spherical[4] * spherical[2] 
                                          + alpha * phydro->w(IVY,k,j,i); 
                    phydro->w(IVZ,k,j,i) = (1-alpha) * (sw_velBC/v0) * spherical[5] 
                                          + alpha * phydro->w(IVZ,k,j,i);
                }
                
                //Due to low SW temperatures and densities near the termination shock, thermal pressure struggles to be reliably derived from conserved energy (P = e - B^2 / 2 - rho * v^2 / 2, e is energy density) and can lead to roundoff errors, negative values, etc.
                //For regions within the termination shock where thermal pressure expansion is purely ~1/r^10/3, we will manually reset the thermal pressure values every timestep to ensure a correct, physically expected profile.
                //This is fine since we are still following the expected physics and helping the solver maintain physically expected values within the termination shock.
                //CHECK IF CORRESPONDING POINT ON MESHBLOCK IS WITHIN EXPECTED MAX TERMINATION SHOCK BOUNDS IN SIMULATION AND THAT VELOCITY IS ROUGHLY THE SW BOUNDARY CONDITION VALUE (little/no forces act on the fluid within TS, so we expect the velocity to be constant within TS)
                Real vel_mag = std::sqrt(std::pow(phydro->w(IVX,k,j,i),2) + 
                                        std::pow(phydro->w(IVY,k,j,i),2) + 
                                        std::pow(phydro->w(IVZ,k,j,i),2));
                if ((spherical[1] > r_outer) && (spherical[1] < 250) && (vel_mag > ((sw_velBC/v0) - 0.3))){
                    phydro->w(IPR,k,j,i) = (thermPr_swBC/P0) * std::pow(r_Inner/spherical[1], 10.0/3.0); //reset cells within termination shock to expected adiabatic expansion values!
                }
            }
        }
    }
    
    // Overwrite Bx (x1-faces)
    //IMPORTANT NOTE: fields used by Athena++ are face centered. To properly assign all magnetic field values, need to loop over first over x1 running from is to ie+1 (face centered means the upper bound for B is extended by 1 vs volume centered hydro variables).
    //Need to use face centered x1 coordinates for first loop, volume centered for x2 and x3 to ensure x1 fields are set at the correct positions.
    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie+1; i++) {
                //Get positions centered on x-face where bx should be set
                Real xFace = pcoord->x1f(i);
                Real y = pcoord->x2v(j);
                Real z = pcoord->x3v(k);
                //Compute spherical radius
                Real r = std::sqrt(xFace*xFace + y*y + z*z);
                
                //If r <r_Inner, this is outside the computational domain, set to zero
                if (r < r_Inner) {
                    pfield->b.x1f(k,j,i) = 0.0;
                }

                //if r > r_inner and r < r_outer, this corresponds to the spherical BC shell. OVERWRITE WITH SW PARKER SPIRAL BX
                else if (r >= r_Inner && r <= r_outer) {
                    std::vector<Real> parker = ComputeParkerSpiralField(xFace, y, z);
                    pfield->b.x1f(k,j,i) = parker[0]; //Parker[0] -> x-component of Parker Spiral field
                }

                //if r (spherical[1]) > r_outer and < r_blend: corresponds to spherical blend region. OVERWRITE WITH MIX OF PARKER AND CURRENT SIM VALUES, WITH CORRESPONDING LINEAR WEIGHTS TO SMOOTH INITIAL SHARP GRADIENTS
                else if (r <= r_blend && r >= r_outer){
                    Real alpha = (r - r_outer) / dr_blend;
                    std::vector<Real> parker = ComputeParkerSpiralField(xFace, y, z);
                    pfield->b.x1f(k,j,i) = (1-alpha) * parker[0] + alpha * pfield->b.x1f(k,j,i);
                }
                //else, leave untouched and let the MHD solver evolve the field
            }
        }
    }
    
    // Overwrite By (x2-faces)
    // Same logic as Bx (x1-faces), but loop from js, je+1 (return x1 bounds to is, ie) and used face centered x2, cell centered x1/x3 to properly assign By on y faces.
    for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
            for (int i=is; i<=ie; i++) {
                //Get positions centered on y-face where by should be set
                Real x = pcoord->x1v(i);
                Real yFace = pcoord->x2f(j);
                Real z = pcoord->x3v(k);
                //Compute spherical radius
                Real r = std::sqrt(x*x + yFace*yFace + z*z);
                
                //Same overwrite logic used in Bx, except replace b2 (by) instead and use Parker[1] corresponding to y-component of Parker Spiral
                if (r < r_Inner) {
                    pfield->b.x2f(k,j,i) = 0.0;
                }
                
                else if (r >= r_Inner && r <= r_outer) {
                    std::vector<Real> parker = ComputeParkerSpiralField(x, yFace, z);
                    pfield->b.x2f(k,j,i) = parker[1];
                }

                else if (r <= r_blend && r >= r_outer){
                    Real alpha = (r - r_outer) / dr_blend;
                    std::vector<Real> parker = ComputeParkerSpiralField(x, yFace, z);
                    pfield->b.x2f(k,j,i) = (1-alpha) * parker[1] + alpha * pfield->b.x2f(k,j,i);
                }
            }
        }
    }
    
    // Overwrite Bz (x3-faces)
    // Same logic as Bx and By, but loop from ks, ke+1 (return x1,x2 bounds to is, ie, js, je) and used face centered x3, cell centered x1/x2 to properly assign Bz on z faces.
    for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
                //Get positions centered on z-face where bz should be set
                Real x = pcoord->x1v(i);
                Real y = pcoord->x2v(j);
                Real zFace = pcoord->x3f(k);
                //Compute spherical radius
                Real r = std::sqrt(x*x + y*y + zFace*zFace);
                
                //Same overwrite logic used in Bx/By, except replace b3 instead and use Parker[2] corresponding to z-component of Parker Spiral
                if (r < r_Inner) {
                    pfield->b.x3f(k,j,i) = 0.0;
                }

                else if (r >= r_Inner && r <= r_outer) {
                    std::vector<Real> parker = ComputeParkerSpiralField(x, y, zFace);
                    pfield->b.x3f(k,j,i) = parker[2];
                }

                else if (r <= r_blend && r >= r_outer){
                    Real alpha = (r - r_outer) / dr_blend;
                    std::vector<Real> parker = ComputeParkerSpiralField(x, y, zFace);
                    pfield->b.x3f(k,j,i) = (1-alpha) * parker[2] + alpha * pfield->b.x3f(k,j,i);
                }
            }
        }
    }
    
    //IMPORTANT: after modifying face centered fields, need to call CalculateCellCenteredField to properly set cell centered fields needed by Athena to calculate magnetic pressure, etc needed for primitive/conserved quantity updates.
    //ALSO: we set primitive variables, need to call PrimitiveToConserved to properly set the conserved quantities used in flux calculations/code evolution
    pfield->CalculateCellCenteredField(pfield->b, pfield->bcc, pcoord, is, ie, js, je, ks, ke);
    peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);     
}

//========================================================================================================================================================
// Function: RefinementCondition
// Purpose: Set up adaptive mesh refinement based on density, temperature gradients
// Output: None
//========================================================================================================================================================
int RefinementCondition(MeshBlock *pmb) {

  //Combine Temp, density, total pressure gradients to get a wholistic approach to refinement
  
  Real max_rhoGrad = 0.0;
  Real max_tempGrad = 0.0;
  Real max_pressGrad = 0.0;
  int counter = 0;
  
  //Loop over whole meshblock
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        
        //Safeguard and don't go over max/min bounds
        if (i == pmb->is || i == pmb->ie || 
            j == pmb->js || j == pmb->je || 
            k == pmb->ks || k == pmb->ke) {
          continue;
        }
        
        //Compute normalized dens gradient in each direction. Track the RELATIVE,FRACTIONAL change across cells. Fractional change = abs(current cell value - previous cell value) / current cell value
        Real dRhox = std::abs(pmb->phydro->w(IDN,k,j,i) - pmb->phydro->w(IDN,k,j,i-1)) / pmb->phydro->w(IDN,k,j,i);
        Real dRhoy = std::abs(pmb->phydro->w(IDN,k,j,i) - pmb->phydro->w(IDN,k,j-1,i)) / pmb->phydro->w(IDN,k,j,i);
        Real dRhoz = std::abs(pmb->phydro->w(IDN,k,j,i) - pmb->phydro->w(IDN,k-1,j,i)) / pmb->phydro->w(IDN,k,j,i);

        //Set max_rhoGrad to maximum gradient only if the gradient in a given direction is greater than current max_rhoGrad, check over all directions
        if (dRhox > max_rhoGrad) max_rhoGrad = dRhox;
        if (dRhoy > max_rhoGrad) max_rhoGrad = dRhoy;
        if (dRhoz > max_rhoGrad) max_rhoGrad = dRhoz;
    
        //Compute relative T gradient. T = (111111 / (gamma-1)) * (press/rho)
        Real tempCell = (pmb->phydro->w(IPR,k,j,i) / pmb->phydro->w(IDN,k,j,i)); //dont need the factor in front of (press/rho) since it will get divided out when we take relative/normalized change across cells
        Real dTx = (tempCell - (pmb->phydro->w(IPR,k,j,i-1) / pmb->phydro->w(IDN,k,j,i-1))) / tempCell;
        Real dTy = (tempCell - (pmb->phydro->w(IPR,k,j-1,i) / pmb->phydro->w(IDN,k,j-1,i))) / tempCell;
        Real dTz = (tempCell - (pmb->phydro->w(IPR,k-1,j,i) / pmb->phydro->w(IDN,k-1,j,i))) / tempCell;

        //Set max_tempGrad to maximum gradient only if the gradient in a given direction is greater than current max_tempGrad, check over all directions
        if (dTx > max_tempGrad) max_tempGrad = dTx;
        if (dTy > max_tempGrad) max_tempGrad = dTy;
        if (dTz > max_tempGrad) max_tempGrad = dTz;

        //Compute relative Total pressure gradient
        Real bMagCell = (std::pow(pmb->pfield->bcc(IB1,k,j,i),2) + std::pow(pmb->pfield->bcc(IB2,k,j,i),2) + std::pow(pmb->pfield->bcc(IB3,k,j,i),2));
        Real bMagCellxMinus1 = (std::pow(pmb->pfield->bcc(IB1,k,j,i-1),2) + std::pow(pmb->pfield->bcc(IB2,k,j,i-1),2) + std::pow(pmb->pfield->bcc(IB3,k,j,i-1),2));
        Real bMagCellyMinus1 = (std::pow(pmb->pfield->bcc(IB1,k,j-1,i),2) + std::pow(pmb->pfield->bcc(IB2,k,j-1,i),2) + std::pow(pmb->pfield->bcc(IB3,k,j-1,i),2));
        Real bMagCellzMinus1 = (std::pow(pmb->pfield->bcc(IB1,k-1,j,i),2) + std::pow(pmb->pfield->bcc(IB2,k-1,j,i),2) + std::pow(pmb->pfield->bcc(IB3,k-1,j,i),2));

        //Calculate total pressures of the cells we need
        Real pressCell =  bMagCell/2 + pmb->phydro->w(IPR,k,j,i);
        Real pressCellxMinus1 =  bMagCellxMinus1/2 + pmb->phydro->w(IPR,k,j,i-1);
        Real pressCellyMinus1 =  bMagCellyMinus1/2 + pmb->phydro->w(IPR,k,j-1,i);
        Real pressCellzMinus1 =  bMagCellzMinus1/2 + pmb->phydro->w(IPR,k-1,j,i);

        //Compute gradient in each direction
        Real dPressx = std::abs(pressCell - pressCellxMinus1) / pressCell;
        Real dPressy = std::abs(pressCell - pressCellyMinus1) / pressCell;
        Real dPressz = std::abs(pressCell - pressCellzMinus1) / pressCell;

        if (dPressx > max_pressGrad) max_pressGrad = dPressx;
        if (dPressy > max_pressGrad) max_pressGrad = dPressy;
        if (dPressz > max_pressGrad) max_pressGrad = dPressz;
      }
    }
  }

  //Tally up counter by seeing if the maximum gradients are larger than thresholds
  if (max_rhoGrad > threshold_densGrad) counter++;
  if (max_tempGrad > threshold_tempGrad) counter++;
  if (max_pressGrad > threshold_pressGrad) counter++;

  if (counter >= 2) { //Refine if two thresholds are satisfied
    return 1;
  } else if (counter == 0) { //Derefine if no thresholds satisfied
    return -1;
  } else {
    return 0; //Else do nothing
  }
}