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

// create the detector
static Ref_t createDetector(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
  xml::Component dims = detElem.dimensions();
  // OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement det(detName, detID);

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
  auto sensorVis = desc.visAttributes(detElem.attr<std::string>(_Unicode(vis_sensor)));

  double sensor_thickness = 1.0 * mm;

  ///----------->>> Define vessel and gas volumes
  /* - `vessel`: aluminum enclosure, mother volume
   * - `gasvol`: gas volume, which fills `vessel`; all other volumes defined below
   *   are children of `gasvol`
   */

  Box vesselSolid(cell_x / 2.,
                  cell_y / 2.,
                  cell_z / 2.);
  Box gasvolSolid(cell_x / 2. - cell_wall_thickness,
                  cell_y / 2. - cell_wall_thickness,
                  cell_z / 2. - cell_wall_thickness);

  Volume vesselVol(detName, vesselSolid, vesselMat);
  Volume gasvolVol(detName + "_gas", gasvolSolid, gasvolMat);
  vesselVol.setVisAttributes(vesselVis);
  gasvolVol.setVisAttributes(gasvolVis);

  ///----------->>> Sensor
  {
    Box sensorShape(cell_x / 2. - cell_wall_thickness,
                    cell_y / 2. - cell_wall_thickness,
                    sensor_thickness / 2.);
    Volume sensorVol(detName + "_sensor", sensorShape, desc.material("Silicon"));
    sensorVol.setVisAttributes(sensorVis);
    gasvolVol.placeVolume(sensorVol, Position(0, 0, -cell_z / 2. + cell_wall_thickness + sensor_thickness / 2.));
  }

  // place gas volume
  PlacedVolume gasvolPV = vesselVol.placeVolume(gasvolVol, Position(0, 0, 0));

  // place mother volume (vessel)
  Volume motherVol = desc.pickMotherVolume(det);
  PlacedVolume vesselPV = motherVol.placeVolume(vesselVol, Position(0, 0, 0));
  vesselPV.addPhysVolID("system", detID);
  det.setPlacement(vesselPV);

  return det;
}

// clang-format off
DECLARE_DETELEMENT(ARCTYPE, createDetector)
