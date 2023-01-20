//----------------------------------
//  pfRICH: Proximity Focusing RICH
//  Author: C. Dilks
//----------------------------------

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include <XML/Helper.h>

using namespace dd4hep;
using namespace dd4hep::rec;
using dd4hep::SubtractionSolid;

// create the detector
static Ref_t createDetector(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
  xml::Component dims = detElem.dimensions();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement det(detName, detID);
  sens.setType("tracker");


  ///----------->>> constant attributes
  // - vessel
  double cell_x = dims.attr<double>(_Unicode(cell_x));
  double cell_y = dims.attr<double>(_Unicode(cell_y));
  double cell_z = dims.attr<double>(_Unicode(cell_z));
  double cell_wall_thickness = dims.attr<double>(_Unicode(cell_wall_thickness));

  auto vesselMat = desc.material(detElem.attr<std::string>(_Unicode(material)));
  auto gasvolMat = desc.material(detElem.attr<std::string>(_Unicode(gas)));
  auto vesselVis = desc.visAttributes(detElem.attr<std::string>(_Unicode(vis_vessel)));
  auto gasvolVis = desc.visAttributes(detElem.attr<std::string>(_Unicode(vis_gas)));


  // - mirror
  auto   mirrorElem         = detElem.child(_Unicode(mirror)).child(_Unicode(module));
  auto   mirrorVis          = desc.visAttributes(mirrorElem.attr<std::string>(_Unicode(vis)));
  auto   mirrorMat          = desc.material(detElem.attr<std::string>(_Unicode(material)));
  auto   mirrorR            = mirrorElem.attr<double>(_Unicode(radius));
  auto   mirrorThickness    = mirrorElem.attr<double>(_Unicode(thickness));

  // - sensor module
  auto   sensorElem         = detElem.child(_Unicode(sensors)).child(_Unicode(module));
  auto   sensorVis          = desc.visAttributes(sensorElem.attr<std::string>(_Unicode(vis)));
  double sensorSide         = sensorElem.attr<double>(_Unicode(side));
  double sensorThickness    = sensorElem.attr<double>(_Unicode(thickness));
  auto   sensorSurf         = surfMgr.opticalSurface(sensorElem.attr<std::string>(_Unicode(surface)));

  ///----------->>> Define vessel and gas volumes
  /* - `vessel`: aluminum enclosure, mother volume
   * - `gasvol`: gas volume, which fills `vessel`; all other volumes defined below
   *   are children of `gasvol`
   */

  double gasUpperZlimit = cell_z / 2. - cell_wall_thickness;
  Box vesselSolid(cell_x / 2.,
                  cell_y / 2.,
                  cell_z / 2.);
  Box gasvolSolid(cell_x / 2. - cell_wall_thickness,
                  cell_y / 2. - cell_wall_thickness,
                  gasUpperZlimit);

  Volume vesselVol(detName, vesselSolid, vesselMat);
  Volume gasvolVol(detName + "_gas", gasvolSolid, gasvolMat);
  vesselVol.setVisAttributes(vesselVis);
  gasvolVol.setVisAttributes(gasvolVis);

  ///----------->>> Sensor
  {
    Box sensorShape( sensorSide/2.0,
                     sensorSide/2.0,
                    sensorThickness / 2.);
    Volume sensorVol(detName + "_sensor", sensorShape, desc.material("AirOptical"));

    sensorVol.setVisAttributes(sensorVis);
    sensorVol.setSensitiveDetector(sens);
    PlacedVolume sensorPV = gasvolVol.placeVolume(sensorVol, Position(0, 0, -cell_z / 2. + cell_wall_thickness + sensorThickness / 2.));
    sensorPV.addPhysVolID("module", 123);

    // Make sensor sensitive + define optical properties
    DetElement sensorDE(det, "ARC_sensor", 123);
    sensorDE.setPlacement(sensorPV);
    SkinSurface sensorSkin(desc, sensorDE, "sensor_optical_surface", sensorSurf, sensorVol); // FIXME: 3rd arg needs `imod`?
    sensorSkin.isValid();
  }
 ///----------->>> Mirror
  {
    // define "mirrorVolFull" as a hollow sphere of Aluminium 
    Sphere mirrorShapeFull( mirrorR - mirrorThickness,
                            mirrorR ,
                            0.,
                            3.14/2);

    // 3D transformation of mirrorVolFull in order to place it inside the gas volume
    Transform3D mirrorTr( RotationZYX(0.,0,0.), Translation3D( 0 , 0 , gasUpperZlimit - mirrorR ) );

    // Define the actual mirror as intersection of the mother volume and the hollow sphere just defined
    Solid mirrorSol = IntersectionSolid(gasvolSolid, mirrorShapeFull, mirrorTr);
    Volume mirrorVol(detName + "_mirror", mirrorSol, mirrorMat);
    mirrorVol.setVisAttributes(mirrorVis);
    /* PlacedVolume mirrorPV =*/ gasvolVol.placeVolume(mirrorVol);



  }



  // place gas volume
  /* PlacedVolume gasvolPV = */ vesselVol.placeVolume(gasvolVol, Position(0, 0, 0));

  // place mother volume (vessel)
  Volume motherVol = desc.pickMotherVolume(det);
  PlacedVolume vesselPV = motherVol.placeVolume(vesselVol, Position(0, 0, 0));
  vesselPV.addPhysVolID("system", detID);
  det.setPlacement(vesselPV);

  return det;
}

// clang-format off
DECLARE_DETELEMENT(ARCTYPE, createDetector)
