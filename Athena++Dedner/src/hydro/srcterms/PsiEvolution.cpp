// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../field/field.hpp"
#include "../../mesh/mesh.hpp"
#include "../../orbital_advection/orbital_advection.hpp"
#include "../hydro.hpp"

// this class header
#include "hydro_srcterms.hpp"

void HydroSourceTerms::PsiEvolution(const Real dt, 
                                  const AthenaArray<Real> &prim,
                                  AthenaArray<Real> &cons) {
                       
  Hydro *ph = pmy_hydro_;
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  Field *pf = pmb->pfield;

  Coordinates *pcoord = pmb->pcoord;
  Real ch = ph->ch_glm;
  Real cp = ph->cp_glm;
  Real damping_coeff = ch*ch / (cp*cp);

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        
        // Compute div(B) using face-centered B field
        Real divB = 0;
        
        // x1-component
          Real dx_center = pcoord->x1f(i+1) - pcoord->x1f(i);
          divB += (pf->b.x1f(k,j,i+1) - pf->b.x1f(k,j,i)) / dx_center;
        
        // x2-component
          Real dy_center = pcoord->x2f(j+1) - pcoord->x2f(j);
          divB += (pf->b.x2f(k,j+1,i) - pf->b.x2f(k,j,i)) / dy_center;
        
        // x3-component
          Real dz_center = pcoord->x3f(k+1) - pcoord->x3f(k);
          divB += (pf->b.x3f(k+1,j,i) - pf->b.x3f(k,j,i)) / dz_center;
        
        Real psi = cons(IPSIC,k,j,i);
        
        cons(IPSIC,k,j,i) -= dt * ((damping_coeff * psi) + (ch * ch * divB));
      }
    }
  }
}