//----------------------------------
//         ARC detector v0
//----------------------------------

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include <XML/Helper.h>

using namespace dd4hep;
using namespace dd4hep::rec;
using dd4hep::SubtractionSolid;

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

/*! 
 *  \brief     Detector constructor of single cubic RICH cell
 *  \details   This code creates a minimal working example of a RICH cell. Evolved from the pfRICH example in DD4hep.
 *  \author    Alvaro Tolosa-Delgado alvaro.tolosa.delgado@cern.ch
 *  \author    Martin Tat            martin.tat@cern.ch
 *  \version   0
 *  \date      2023
 *  \pre       DD4hep compiled with Geant4+Qt
 *  \bug       Walls do not reflect optical photons
 */
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
  auto mirrorElem = detElem.child(_Unicode(mirror)).child(_Unicode(module));
  auto mirrorVis = desc.visAttributes(mirrorElem.attr<std::string>(_Unicode(vis)));
  auto mirrorMat = desc.material(mirrorElem.attr<std::string>(_Unicode(material)));
  auto mirrorR = mirrorElem.attr<double>(_Unicode(radius));
  auto mirrorThickness = mirrorElem.attr<double>(_Unicode(thickness));
  auto mirrorSurf = surfMgr.opticalSurface(mirrorElem.attr<std::string>(_Unicode(surface)));

  // - sensor module
  auto sensorElem = detElem.child(_Unicode(sensors)).child(_Unicode(module));
  auto sensorVis = desc.visAttributes(sensorElem.attr<std::string>(_Unicode(vis)));
  double sensorX = sensorElem.attr<double>(_Unicode(sensorX));
  double sensorY = sensorElem.attr<double>(_Unicode(sensorY));
  double sensorThickness = sensorElem.attr<double>(_Unicode(thickness));
  auto sensorSurf = surfMgr.opticalSurface(sensorElem.attr<std::string>(_Unicode(surface)));
  auto sensorMat = desc.material(sensorElem.attr<std::string>(_Unicode(material)));

  // - aerogel
  auto aerogelElem = detElem.child(_Unicode(aerogel)).child(_Unicode(module));
  auto aerogelVis = desc.visAttributes(aerogelElem.attr<std::string>(_Unicode(vis)));
  auto aerogelMat = desc.material(aerogelElem.attr<std::string>(_Unicode(material)));
  double aerogel_thickness = aerogelElem.attr<double>(_Unicode(thickness));

  // - cooling
  auto coolingElem = detElem.child(_Unicode(cooling)).child(_Unicode(module));
  auto coolingVis = desc.visAttributes(coolingElem.attr<std::string>(_Unicode(vis)));
  auto coolingMat = desc.material(coolingElem.attr<std::string>(_Unicode(material)));
  double cooling_thickness = coolingElem.attr<double>(_Unicode(thickness));

  ///----------->>> Define vessel and gas volumes
  /* - `vessel`: aluminum enclosure, mother volume
   * - `gasvol`: gas volume, which fills `vessel`; all other volumes defined below
   *   as children of `gasvol`. Sensor is placed inside Cooling.
   * vessel (cubic) -> Gasvol (cubic) -> Cooling (cubic) -> Sensor CCD (cubic)
   *                                 \-> Mirror (sphere intersection with gasvol)
   *                                 \-> Aerogel
   */

  // Vessel
  Box vesselSolid(cell_x / 2.,
                  cell_y / 2.,
                  cell_z / 2.);
  Volume vesselVol(detName, vesselSolid, vesselMat);
  vesselVol.setVisAttributes(vesselVis);

  // Gas
  // Thickness of gas volume (z-direction) if we ignore the mirror
  double gasThickness = cell_z - 2 * cell_wall_thickness;
  Box gasvolSolid(cell_x / 2. - cell_wall_thickness,
                  cell_y / 2. - cell_wall_thickness,
                  gasThickness / 2.);
  Volume gasvolVol(detName + "_gas", gasvolSolid, gasvolMat);
  gasvolVol.setVisAttributes(gasvolVis);
  /* PlacedVolume gasvolPV = */ vesselVol.placeVolume(gasvolVol, Position(0, 0, 0));

  ///----------->>> Aerogel (+sensor)
  {
    Box coolingSolid(cell_x / 2. - cell_wall_thickness,
                     cell_y / 2. - cell_wall_thickness,
                     cooling_thickness / 2.);
    Volume coolingVol(detName + "cooling", coolingSolid, coolingMat);
    coolingVol.setVisAttributes(coolingVis);

    // place cooling volume
    double coolingCentre = -gasThickness / 2. + cooling_thickness / 2.;
    PlacedVolume coolingPV = gasvolVol.placeVolume(coolingVol, Position(0, 0, coolingCentre));
    coolingPV.addPhysVolID("module", 63);

    ///----------->>> Sensor
    {
      Box sensorShape(sensorX / 2.,
                      sensorY / 2.,
                      sensorThickness / 2.);
      Volume sensorVol(detName + "_sensor", sensorShape, sensorMat);

      sensorVol.setVisAttributes(sensorVis);
      sensorVol.setSensitiveDetector(sens);
      double sensorCentre = cooling_thickness / 2. - sensorThickness / 2.;
      PlacedVolume sensorPV = coolingVol.placeVolume(sensorVol, Position(0, 0, sensorCentre));
      sensorPV.addPhysVolID("module", 127);

      // // Make sensor sensitive + define optical properties
      // DetElement sensorDE(aerogelDE, "ARC_sensor", 127);
      // sensorDE.setPlacement(sensorPV);
      // SkinSurface sensorSkin(desc, sensorDE, "sensor_optical_surface", sensorSurf, sensorVol); // FIXME: 3rd arg needs `imod`?
      // sensorSkin.isValid();
    }
  }
  ///----------->>> Aerogel (+sensor)
  {
    Box aerogelSolid(cell_x / 2. - cell_wall_thickness,
                     cell_y / 2. - cell_wall_thickness,
                     aerogel_thickness / 2.);
    Volume aerogelVol(detName + "_aerogel", aerogelSolid, aerogelMat);
    aerogelVol.setVisAttributes(aerogelVis);

    // place aerogel volume
    // z-position of gas volume
    double aerogelCentre = cooling_thickness - gasThickness / 2. + aerogel_thickness / 2.;
    /*PlacedVolume aerogelPV = */ gasvolVol.placeVolume(aerogelVol, Position(0, 0, aerogelCentre));
  }
  ///----------->>> Mirror
  {
    // define "mirrorVolFull" as a hollow sphere of Aluminium
    Sphere mirrorShapeFull(mirrorR - mirrorThickness,
                           mirrorR,
                           0.,
                           3.14 / 2);

    // 3D transformation of mirrorVolFull in order to place it inside the gas volume
    Transform3D mirrorTr(RotationZYX(0., 0, 0.), Translation3D(0, 0, gasThickness / 2. - mirrorR));

    // Define the actual mirror as intersection of the mother volume and the hollow sphere just defined
    Solid mirrorSol = IntersectionSolid(gasvolSolid, mirrorShapeFull, mirrorTr);
    Volume mirrorVol(detName + "_mirror", mirrorSol, mirrorMat);
    mirrorVol.setVisAttributes(mirrorVis);
    PlacedVolume mirrorPV = gasvolVol.placeVolume(mirrorVol);
    mirrorPV.addPhysVolID("module", 3);
    DetElement mirrorDE(det, "ARC_mirror", 3);
    mirrorDE.setPlacement(mirrorPV);
    SkinSurface mirrorSkin(desc, mirrorDE, "mirror_optical_surface", mirrorSurf, mirrorVol); // FIXME: 3rd arg needs `imod`?
    mirrorSkin.isValid();
  }

  // place mother volume (vessel)
  Volume motherVol = desc.pickMotherVolume(det);
  PlacedVolume vesselPV = motherVol.placeVolume(vesselVol, Position(0, 0, 0));
  vesselPV.addPhysVolID("system", detID);
  det.setPlacement(vesselPV);

  return det;
}

// clang-format off
DECLARE_DETELEMENT(ARCTYPE, createDetector)
