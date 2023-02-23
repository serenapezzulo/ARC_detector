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

// local stuff inside anonymous namespace to avoid collisions
namespace
{
  // store parameters by name
  std::map<std::string, double> cell_parameters_m;

  /// tokenize string by space as delimiter
  void mytokenizer(std::string &istring, std::vector<std::string> &tokens)
  {
    std::stringstream myline_ss(istring);
    std::string intermediate;
    while (getline(myline_ss, intermediate, ' '))
      tokens.push_back(intermediate);
  }

  // function to fill map with parameters cell_parameters_m
  void fill_cell_parameters_m()
  {
    // avoid calling this function twice
    if (cell_parameters_m.size())
      return;

    // hardcoded, to be developed later
    std::ifstream ifile("RadiatorCell_FinaOptimisation.txt");

    // prepare some counters for later sanity check
    int Curvature_counter(0);
    int XPosition_counter(0);
    int ZPosition_counter(0);
    int DetPosition_counter(0);
    int DetTilt_counter(0);
    int barrel_unique_cells(0);
    int endcap_unique_cells(0);

    while (ifile.good())
    {
      // read one line and tokenize by delimiter
      std::string myline("");
      std::getline(ifile, myline);
      // skip if empty line of line start with #
      if ((0 == myline.size()) || (myline.size() && '#' == myline[0]))
        continue;

      std::vector<std::string> tokens;
      mytokenizer(myline, tokens);
      // skip if not 2 elements are provided
      if (2 != tokens.size())
        continue;

      std::string &parname = tokens.at(0);
      cell_parameters_m[parname] = atof(tokens[1].c_str());

      // increase corresponding parameter counter
      // and calibrate parameter according to Martin units
      if (std::string::npos != parname.find("Curvature"))
      {
        ++Curvature_counter;
        cell_parameters_m[parname] = cell_parameters_m[parname] * m;
      }
      else if (std::string::npos != parname.find("XPosition"))
      {
        ++XPosition_counter;
        cell_parameters_m[parname] = cell_parameters_m[parname] * m;
      }
      else if (std::string::npos != parname.find("ZPosition"))
      {
        ++ZPosition_counter;
        cell_parameters_m[parname] = cell_parameters_m[parname] * m;
      }
      else if (std::string::npos != parname.find("DetPosition"))
      {
        ++DetPosition_counter;
        cell_parameters_m[parname] = cell_parameters_m[parname] * m;
      }
      else if (std::string::npos != parname.find("DetTilt"))
      {
        ++DetTilt_counter;
        cell_parameters_m[parname] = cell_parameters_m[parname] * rad;
      }
      if (std::string::npos != parname.find("EndCapRadiator"))
        ++endcap_unique_cells;
      else if (std::string::npos != parname.find("Radiator"))
        ++barrel_unique_cells;
    }
    ifile.close();

    // normalize to the number of parameters per cell
    endcap_unique_cells /= 5;
    barrel_unique_cells /= 5;

    // check if number of parameters is ok, if not, throw exception
    if (23 != endcap_unique_cells)
      throw std::runtime_error("Number of endcap cells different from expected (23)");
    if (18 != barrel_unique_cells)
      throw std::runtime_error("Number of barrel cells different from expected (18)");
    if (0 != Curvature_counter - endcap_unique_cells - barrel_unique_cells)
      throw std::runtime_error("Number of Curvature parameters different from expected (23+18)");
    if (0 != XPosition_counter - endcap_unique_cells - barrel_unique_cells)
      throw std::runtime_error("Number of XPosition parameters different from expected (23+18)");
    if (0 != ZPosition_counter - endcap_unique_cells - barrel_unique_cells)
      throw std::runtime_error("Number of ZPosition parameters different from expected (23+18)");
    if (0 != DetPosition_counter - endcap_unique_cells - barrel_unique_cells)
      throw std::runtime_error("Number of DetPosition parameters different from expected (23+18)");
    if (0 != DetTilt_counter - endcap_unique_cells - barrel_unique_cells)
      throw std::runtime_error("Number of DetTilt parameters different from expected (23+18)");

    return;
  } // end void fill_cell_parameters_m()

} // end anonymous namespace

// create the detector
static Ref_t create_barrel_cell(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
  xml::Component dims = detElem.dimensions();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement det(detName, detID);
  sens.setType("tracker");

  // read Martin file and store parameters by name in the map
  fill_cell_parameters_m();

  // mother volume corresponds to the world
  Volume motherVol = desc.pickMotherVolume(det);

  double vessel_outer_r = 210 * cm;
  double vessel_inner_r = 190 * cm;
  double vessel_length = 440 * cm;
  double vessel_wall_thickness = 0.1 * cm;
  double hexagon_side_length = 14.815 * cm;
  double thickness_sphere(10 * mm);
  double zstep = 2 * hexagon_side_length;
  double phistep = 13.333 * deg;

  if (vessel_outer_r <= vessel_inner_r)
    throw std::runtime_error("Ilegal parameters: vessel_outer_r <= vessel_inner_r");

  // Build cylinder for gas.
  Tube gasvolSolid(vessel_inner_r,
                   vessel_outer_r,
                   vessel_length);

  // Use regular polyhedra for endcaps cells
  // PolyhedraRegular cellS("aa",6,0,4*cm);

  // Use pyramid for barrel cells
  std::vector<double> zplanes = {0 * cm, vessel_outer_r};
  std::vector<double> rs = {0 * cm, hexagon_side_length};
  /// Hexagonal pyramid
  Polyhedra shape("mypyramid", 6, 30 * deg, 360 * deg, zplanes, rs);
  /// rotation of 90deg around Y axis, to align Z axis of pyramid with X axis of cylinder
  Transform3D pyramidTr(RotationZYX(0, 90. * deg, 0. * deg), Translation3D(0, 0, 0));



  // Build the mirror for ncell=1..18

  std::vector<int> ncell_vector = {-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13, -14,-15,-16,-17,-18,
                                  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18
                                  };
  // int ncell = 1;
  // for (int ncell = 1; ncell <= 18; ++ncell)
  for (auto ncell : ncell_vector)
  {
    /// cell shape. Coordinate system still the same as cylinder!
    Solid cellS = IntersectionSolid(gasvolSolid, shape, pyramidTr);
    Volume cellVol(detName + "_cell" + std::to_string(ncell), cellS, desc.material("C4F10_PFRICH"));
    cellVol.setVisAttributes(desc.visAttributes("gas_vis"));

    // there is no cell number 0, and cell number 1 do not need to be reflected
    if (0 == ncell || -1 == ncell)
      continue;
    bool reflect_parameters = false;
    if (0 > ncell)
    {
      ncell *= -1;
      reflect_parameters = true;
    }

    // initialize parameters for creating the mirror
    double center_of_sphere_x(-1);
    double center_of_sphere_z(-1);
    double radius_of_sphere(-1);

    int name_col = ncell / 2;
    int name_row = ncell % 2 ? 1 : 2;
    // retrieve stored parameters
    {
      std::string name_col_s = std::to_string(name_col);
      std::string name_row_s = std::to_string(name_row);
      radius_of_sphere = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_Curvature"];
      center_of_sphere_x = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_XPosition"];
      double zposition = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_ZPosition"];
      center_of_sphere_z = vessel_outer_r - vessel_wall_thickness - radius_of_sphere - zposition;

      if (radius_of_sphere <= thickness_sphere)
        throw std::runtime_error("Ilegal parameters: radius_of_sphere <= thickness_sphere");
    }

    if (reflect_parameters)
    {
      center_of_sphere_x *= -1.0;
    }

    // create the semi-sphere that will result in the mirror
    Sphere mirrorShapeFull(radius_of_sphere - thickness_sphere,
                           radius_of_sphere,
                           0.,
                           3.14 / 2);

    /// 3D transformation of mirrorVolFull in order to place it inside the gas volume
    Transform3D mirrorTr(RotationZYX(0, 0, 0), Translation3D(center_of_sphere_x, 0, center_of_sphere_z));

    /// Define the actual mirror as intersection of the mother volume and the hollow sphere just defined
    Solid mirrorSol = IntersectionSolid(shape, mirrorShapeFull, mirrorTr);
    Volume mirrorVol(detName + "_mirror" + std::to_string(ncell)+ "z" + std::to_string(reflect_parameters), mirrorSol, desc.material("Aluminum"));
    mirrorVol.setVisAttributes(desc.visAttributes(Form("mirror_vis%d", ncell)));
    cellVol.placeVolume(mirrorVol, pyramidTr);

    // if (bool build_full_barrel = false)
    // {

    //   for (int ringn = -8; ringn <= 8; ++ringn)
    //   {
    //     for (int phin = 0; phin < 27; ++phin)
    //     {
    //       PlacedVolume cellPV = motherVol.placeVolume(cellVol, RotationZ(phistep * phin) * Translation3D(0, 0, ringn * 29.63 * cm));
    //       cellPV.addPhysVolID("system", detID).addPhysVolID("module", 17 * phin + ringn);
    //       // create mirrors as separate detectors, so properties can be adjusted lated!
    //       det.setPlacement(cellPV);
    //     }
    //   }
    //   for (double ringn = -7.5; ringn <= 7.5; ++ringn)
    //   {
    //     for (int phin = 0; phin < 27; ++phin)
    //     {
    //       PlacedVolume cellPV = motherVol.placeVolume(cellVol, RotationZ(phistep * (0.5 + phin)) * Translation3D(0, 0, ringn * 29.63 * cm));
    //       cellPV.addPhysVolID("system", detID).addPhysVolID("module", 17 * phin + ringn);
    //       // create mirrors as separate detectors, so properties can be adjusted lated!
    //       det.setPlacement(cellPV);
    //     }
    //   }
    //   return det;
    // }
    // else
    {
      double mirror_abs_pos_z = name_col * zstep - 0.5 * zstep * (2 == name_row);
      if (reflect_parameters)
        mirror_abs_pos_z *= -1.0;
      double phi_offset = 0 + 0.5 * phistep * (2 == name_row);

      for (int phin = 0; phin < 27; ++phin)
      {
        PlacedVolume cellPV = motherVol.placeVolume(cellVol, RotationZ(phistep * phin + phi_offset) * Translation3D(0, 0, mirror_abs_pos_z));
        cellPV.addPhysVolID("system", detID).addPhysVolID("module", 17 * phin + name_col);
        // create mirrors as separate detectors, so properties can be adjusted lated!
        det.setPlacement(cellPV);
      }
    }
  }

  return det;
}
DECLARE_DETELEMENT(ARCBARREL_T, create_barrel_cell)

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
