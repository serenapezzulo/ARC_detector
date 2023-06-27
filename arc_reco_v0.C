
// Framework include files
#include "DD4hep/Detector.h"
#include "DD4hep/DDTest.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/BitFieldCoder.h"
#include "DDRec/CellIDPositionConverter.h"


R__LOAD_LIBRARY(libDDCore) 

int ncell = 1;
int nphi = 1;

void arc_reco_v0()
{

auto lcdd = &(dd4hep::Detector::getInstance());
lcdd->fromCompact("./compact/arc_barrel_v0.xml");
auto barrelDE = lcdd->detector("ARCBARREL");
auto cellDE = barrelDE.child( Form("ARCBARREL_cell%d_ref0_phi%dDE", ncell, nphi) );
auto mirrorDE = cellDE.child(Form("ARCBARREL_mirror%d_ref0_phi%dDE", ncell, nphi) );

auto mirrorSolid = mirrorDE.solid();
auto m = (TGeoCompositeShape*)mirrorSolid.access();

// this is just the matrix for building the intersection, translation
auto matrix_boolean = m->GetBoolNode()->GetRightMatrix(); //->Print()
// matrix for placin the cell ( TranslationZ(ncell)*RotationZ(phin)*RotationY(-90) )
auto matrix_cellvol = cellDE.nominal().worldTransformation(); //.Print()

double local_coord [3] = {0.,0.,0.};
double global_coord[3] = {0.,0.,0.};

(matrix_cellvol**matrix_boolean).LocalToMaster(local_coord, global_coord);

double center_sphere_global_rho = sqrt(global_coord[0]*global_coord[0] + global_coord[1]*global_coord[1]);
double center_sphere_global_z   = global_coord[2];

std::cout << "center_sphere_global_rho = " << center_sphere_global_rho << " cm"  << std::endl;
std::cout << "center_sphere_global_z = "   << center_sphere_global_z   << " cm"  << std::endl;

TGeoSphere * spherical_mirror_shape = (TGeoSphere*)m->GetBoolNode()->GetRightShape()

double mirror_r_min = spherical_mirror_shape->GetRmin();
std::cout << "radius_sphere = "   << mirror_r_min << " cm" << std::endl;

return;
}
