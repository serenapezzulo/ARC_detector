//----------------------------------
//         ARC detector v0
//----------------------------------

/*!
 *  \brief     Detector constructor of barrel of ARC detector
 *  \details   This code creates full geometry of barrel
 *             Evolved from the pfRICH example in DD4hep.
 *  \author    Alvaro Tolosa-Delgado alvaro.tolosa.delgado@cern.ch
 *  \author    Martin Tat            martin.tat@cern.ch
 *  \version   0
 *  \date      2023
 *  \pre       DD4hep compiled with Geant4+Qt
 *  \bug       Walls do not reflect optical photons. Hard-coded values in many places.
 */

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include <XML/Helper.h>

using namespace dd4hep;
using namespace dd4hep::rec;
using dd4hep::SubtractionSolid;


#include "ARC_par_reader.hpp"

/**
 * create barrel as sum of single cells. DEPRECATED!!
 * next step is to place mirrors+sensors in a cylindral shape gas volume
 */
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

  // Vessel, cylindral
  double vessel_outer_r = 210 * cm;
  double vessel_inner_r = 190 * cm;
  double vessel_length = 440 * cm;
  double vessel_wall_thickness = 1.0 * cm;
  if (vessel_outer_r <= vessel_inner_r)
    throw std::runtime_error("Ilegal parameters: vessel_outer_r <= vessel_inner_r");

  // Cooling,
  double cooling_radial_thickness = 1 * cm;

  // Cell parameters
  /// Cell is intersection of hexagonal pyramid and the cylinder
  double hexagon_side_length = 14.815 * cm;
  /// Distance in x-direction
  double zstep = 2 * hexagon_side_length;
  /// Distance in phi angle between cells
  /// since cells are regular hexagons, this distance matches
  /// the angle that one cell covers
  double phistep = 13.333 * deg;
  /// number of repetition of unique cells around the barrel
  int phinmax = 1; // 27;

  // Mirror parameters
  double thickness_sphere(10 * mm);
  double mirror_diameter_safe_shrink = 1*mm;
  double mirror_z_safe_shrink = 1*mm;
  double mirror_z_origin_Martin = vessel_outer_r - vessel_wall_thickness - 37 * cm - mirror_z_safe_shrink;

  // Light sensor parameters
  double sensor_sidex = 8 * cm;
  double sensor_sidey = 8 * cm;
  double sensor_thickness = 0.2 * cm;
  double sensor_z_origin_Martin = vessel_inner_r + vessel_wall_thickness + cooling_radial_thickness;

  // Build cylinder for gas.
  Tube gasvolSolid(vessel_inner_r + vessel_wall_thickness,
                   vessel_outer_r - vessel_wall_thickness,
                   vessel_length);

  // Use regular polyhedra for endcaps cells
  // PolyhedraRegular cellS("aa",6,0,4*cm);

  // Use pyramid for barrel cells
  std::vector<double> zplanes = {0 * cm, vessel_outer_r - vessel_wall_thickness - mirror_z_safe_shrink};
  std::vector<double> rs = {0 * cm, hexagon_side_length - mirror_diameter_safe_shrink};
  /// Hexagonal pyramid
  Polyhedra shape("mypyramid", 6, 30 * deg, 360 * deg, zplanes, rs);
  /// rotation of 90deg around Y axis, to align Z axis of pyramid with X axis of cylinder
  Transform3D pyramidTr(RotationZYX(0, -90. * deg, 0. * deg), Translation3D(0, 0, 0));

  // Build sensor shape
  Box sensorSol(sensor_sidex / 2, sensor_sidey / 2, sensor_thickness / 2);
  Volume sensorVol(detName + "_sensor", sensorSol, desc.material("Aluminum"));

  // Build the mirror for ncell=1..18
  std::vector<int> ncell_vector = {-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18,
                                   1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
  // std::vector<int> ncell_vector = { /*-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13, -14,-15,-16,-17,-18,*/
  //                                 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18
  //                                 };
  ncell_vector = {7};
  for (auto ncell : ncell_vector)
  {
    // The following line skips even number cells
    // if (!(ncell % 2))
    //   continue;

    // The following line skips odd number cells
    // if ((ncell % 2))
    // continue;

    /// cell shape. Coordinate system still the same as cylinder!
    Solid cellS = IntersectionSolid(gasvolSolid, shape, pyramidTr);
    Volume cellVol(detName + "_cell" + std::to_string(ncell), cellS, desc.material("C4F10_PFRICH"));
    cellVol.setVisAttributes(desc.visAttributes("gas_vis"));

    // there is no cell number 0, and cell number 1 do not need to be reflected
    if (0 == ncell || -1 == ncell)
      continue;

    // cells at z
    bool reflect_parameters = false;
    if (0 > ncell)
    {
      ncell *= -1;
      reflect_parameters = true;
    }

    // initialize parameters for creating the mirror
    double center_of_sphere_x(-999.);
    double center_of_sphere_z(-999.);
    double radius_of_sphere(-999.);

    double center_of_sensor_x(-999.);
    double angle_of_sensor(-999.);

    // convert Roger nomenclature (one cell number) to Martin nomenclature (row and col numbers)
    int name_col = ncell / 2;
    int name_row = ncell % 2 ? 1 : 2;
    // retrieve stored parameters
    {
      std::string name_col_s = std::to_string(name_col);
      std::string name_row_s = std::to_string(name_row);
      radius_of_sphere = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_Curvature"];
      center_of_sphere_x = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_XPosition"];
      double zposition = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_ZPosition"];

      center_of_sensor_x = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_DetPosition"];
      angle_of_sensor = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_DetTilt"];
      center_of_sphere_z = mirror_z_origin_Martin + zposition;

      // check if parameters are ok
      if (-999. == center_of_sphere_x)
        throw std::runtime_error("Ilegal parameters: center_of_sphere_x not provided");
      if (-999. == center_of_sphere_z)
        throw std::runtime_error("Ilegal parameters: center_of_sphere_z not provided");
      if (-999. == radius_of_sphere)
        throw std::runtime_error("Ilegal parameters: radius_of_sphere not provided");
      if (radius_of_sphere <= thickness_sphere)
        throw std::runtime_error("Ilegal parameters: radius_of_sphere <= thickness_sphere");

      if (-999. == center_of_sensor_x)
        throw std::runtime_error("Ilegal parameters: center_of_sensor_x not provided");
      if (-999. == angle_of_sensor)
        throw std::runtime_error("Ilegal parameters: angle_of_sensor not provided");
    }

    if (reflect_parameters)
    {
      center_of_sphere_x *= -1.0;
      center_of_sensor_x *= -1.0;
      angle_of_sensor *= -1.0;
    }

    // create the semi-sphere that will result in the mirror
    Sphere mirrorShapeFull(radius_of_sphere - thickness_sphere,
                           radius_of_sphere,
                           0.,
                           3.14 / 2);
    /// 3D transformation of mirrorVolFull in order to place it inside the gas volume
    Transform3D mirrorTr(RotationZYX(0, 0, 0), Translation3D(center_of_sphere_x, 0, center_of_sphere_z - mirror_z_safe_shrink));

    // TODO: cell 18 corresponds to half a pyramid, currently is full pyramid
    /// Define the actual mirror as intersection of the mother volume and the hollow sphere just defined
    Solid mirrorSol = IntersectionSolid(shape, mirrorShapeFull, mirrorTr);
    Volume mirrorVol(detName + "_mirror" + std::to_string(ncell) + "z" + std::to_string(reflect_parameters), mirrorSol, desc.material("Aluminum"));
    mirrorVol.setVisAttributes(desc.visAttributes(Form("mirror_vis%d", ncell)));
    cellVol.placeVolume(mirrorVol, pyramidTr);

    // Place detector in cell
    Transform3D sensorTr(RotationZYX(0, 90 * deg - angle_of_sensor, 0), Translation3D(-sensor_z_origin_Martin, 0, center_of_sensor_x));
    cellVol.placeVolume(sensorVol, sensorTr);

    // position of mirror in cylinder coordinate system
    double mirror_abs_pos_z = name_col * zstep - 0.5 * zstep * (2 == name_row);
    if (reflect_parameters)
      mirror_abs_pos_z *= -1.0;

    // row 2 is shifted half step size
    double phi_offset = 0 + 0.5 * phistep * (2 == name_row);

    for (int phin = 0; phin < phinmax; ++phin)
    {
      PlacedVolume cellPV = motherVol.placeVolume(cellVol, RotationZ(phistep * phin + phi_offset) * Translation3D(0, 0, mirror_abs_pos_z));
      cellPV.addPhysVolID("system", detID).addPhysVolID("module", 17 * phin + name_col);
      // create mirrors as separate detectors, so properties can be adjusted lated!
      det.setPlacement(cellPV);
    }
  }

  return det;
}

static Ref_t create_barrel(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
  xml::Component dims = detElem.dimensions();
  DetElement det(detName, detID);
  sens.setType("tracker");
  // Move to begining of detector constructor?
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  auto sensorSurf = surfMgr.opticalSurface("SensorSurface_PFRICH");
  auto mirrorSurf = surfMgr.opticalSurface("MirrorSurface");

  // Vessel, cylindral
  double vessel_outer_r = 210 * cm;
  double vessel_inner_r = 190 * cm;
  double vessel_length = 440. * cm;
  double vessel_wall_thickness = 1.0 * cm;
  if (vessel_outer_r <= vessel_inner_r)
    throw std::runtime_error("Ilegal parameters: vessel_outer_r <= vessel_inner_r");

  // Cooling,
  double cooling_radial_thickness = 1 * cm;

  // Cell parameters
  /// Cell is intersection of hexagonal pyramid and the cylinder
  double hexagon_side_length = 14.815 * cm;
  double hexagon_apothem = hexagon_side_length * cos(30 * deg);
  /// Distance in x-direction
  double zstep = 2 * hexagon_apothem;
  /// Distance in phi angle between cells
  /// since cells are regular hexagons, this distance matches
  /// the angle that one cell covers
  double phistep = 13.333 * deg;
  /// number of repetition of unique cells around the barrel
  int phinmax = 1; // 27;

  // Mirror parameters
  double thickness_sphere(1 * mm);
  double mirror_diameter_safe_shrink = 1*mm;
  double mirror_z_safe_shrink = 2.6*mm;
  double mirror_z_origin_Martin = vessel_outer_r - vessel_wall_thickness - 37 * cm;

  // Light sensor parameters
  double sensor_sidex = 8 * cm;
  double sensor_sidey = 8 * cm;
  double sensor_thickness = 0.2 * cm;
  double sensor_z_safe_shrink = 1.5*mm;
  double sensor_z_origin_Martin = vessel_inner_r + vessel_wall_thickness + cooling_radial_thickness + sensor_z_safe_shrink;

  // read Martin file and store parameters by name in the map
  fill_cell_parameters_m();

  ///----------->>> Define vessel and gas volumes
  /* - `vessel`: aluminum enclosure, mother volume of the barrel
   * - `gasvol`: gas volume, which fills `vessel`; all other volumes defined below
   *   as children of `gasvol`. Sensor is placed inside Cooling.
   * vessel (cylind) -> Gasvol (cylind) -> Cooling (cylind) -> Sensor CCD (cylind)
   *                                   \-> Mirror (sphere intersection with cell volume)
   *                                   \-> Aerogel (cylind)
   */
  // Build cylinder vol for vessel.
  Tube vesselSolid(vessel_inner_r,
                   vessel_outer_r,
                   vessel_length / 2.);
  Volume vesselVol(detName + "_vessel", vesselSolid, desc.material("CarbonFibStr"));
  vesselVol.setVisAttributes(desc.visAttributes("vessel_vis"));

  // Build cylinder vol for gas. It is placed inside vessel
  Tube gasvolSolid(vessel_inner_r + vessel_wall_thickness,
                   vessel_outer_r - vessel_wall_thickness,
                   vessel_length / 2.);
  Volume gasVol(detName + "_gas", gasvolSolid, desc.material("C4F10_PFRICH"));
  gasVol.setVisAttributes(desc.visAttributes("gas_vis"));
#ifdef __CREATE_COOLING__
  // Build cylinder vol for cooling.It is placed inside gas
  Tube coolinSolid(vessel_inner_r + vessel_wall_thickness,
                   vessel_inner_r + vessel_wall_thickness + cooling_radial_thickness,
                   vessel_length / 2.);
  Volume coolingVol(detName + "_cooling", coolinSolid, desc.material("Copper"));
  coolingVol.setVisAttributes(desc.visAttributes("cooling_vis"));
#endif

  //----->> Place mirrors and sensors
  {

    // Use pyramid for barrel cells
    std::vector<double> zplanes = {0 * cm, vessel_outer_r - vessel_wall_thickness - mirror_z_safe_shrink};
    std::vector<double> rs = {0 * cm, hexagon_apothem};
    /// Hexagonal pyramid
    Polyhedra shape("mypyramid", 6, 30 * deg, 360 * deg, zplanes, rs);
    /// rotation of 90deg around Y axis, to align Z axis of pyramid with X axis of cylinder
    Transform3D pyramidTr(RotationZYX(0, -90. * deg, 0. * deg), Translation3D(0, 0, 0));

    // Build sensor shape
    Box sensorSol(sensor_sidex / 2, sensor_sidey / 2, sensor_thickness / 2);
    Volume sensorVol(detName + "_sensor", sensorSol, desc.material("AirOptical"));
    sensorVol.setSensitiveDetector(sens);
    SkinSurface sensorSkin(desc, det, "sensor_optical_surface", sensorSurf, sensorVol); // FIXME: 3rd arg needs `imod`?
    sensorSkin.isValid();

    // Build the mirror for ncell=1..18
    std::vector<int> ncell_vector = {-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18,
                                     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
    // std::vector<int> ncell_vector = { /*-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13, -14,-15,-16,-17,-18,*/
    //                                 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18
    //                                 };
    ncell_vector = {17};

    /// Dummy counter to place elements inside the barrel
    int sensorcounter(0);
    for (auto ncell : ncell_vector)
    {
      // The following line skips even number cells
      // if (!(ncell % 2))
      //   continue;

      // The following line skips odd number cells
      // if ((ncell % 2))
      // continue;

      // there is no cell number 0, and cell number 1 do not need to be reflected
      if (0 == ncell || -1 == ncell)
        continue;

      // cells at z
      bool reflect_parameters = false;
      if (0 > ncell)
      {
        ncell *= -1;
        reflect_parameters = true;
      }
      // initialize parameters for creating the mirror
      double center_of_sphere_x(-999.);
      double center_of_sphere_z(-999.);
      double radius_of_sphere(-999.);

      double center_of_sensor_x(-999.);
      double angle_of_sensor(-999.);

      // convert Roger nomenclature (one cell number) to Martin nomenclature (row and col numbers)
      int name_col = ncell / 2;
      int name_row = ncell % 2 ? 1 : 2;
      // retrieve stored parameters
      {
        std::string name_col_s = std::to_string(name_col);
        std::string name_row_s = std::to_string(name_row);
        radius_of_sphere = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_Curvature"];
        center_of_sphere_x = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_XPosition"];
        double zposition = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_ZPosition"];

        center_of_sensor_x = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_DetPosition"];
        angle_of_sensor = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_DetTilt"];
        center_of_sphere_z = mirror_z_origin_Martin + zposition;

        // check if parameters are ok
        if (-999. == center_of_sphere_x)
          throw std::runtime_error("Ilegal parameters: center_of_sphere_x not provided");
        if (-999. == center_of_sphere_z)
          throw std::runtime_error("Ilegal parameters: center_of_sphere_z not provided");
        if (-999. == radius_of_sphere)
          throw std::runtime_error("Ilegal parameters: radius_of_sphere not provided");
        if (radius_of_sphere <= thickness_sphere)
          throw std::runtime_error("Ilegal parameters: radius_of_sphere <= thickness_sphere");

        if (-999. == center_of_sensor_x)
          throw std::runtime_error("Ilegal parameters: center_of_sensor_x not provided");
        if (-999. == angle_of_sensor)
          throw std::runtime_error("Ilegal parameters: angle_of_sensor not provided");
      }

      if (reflect_parameters)
      {
        center_of_sphere_x *= -1.0;
        center_of_sensor_x *= -1.0;
        angle_of_sensor *= -1.0;
      }

      // create the semi-sphere that will result in the mirror
      Sphere mirrorShapeFull(radius_of_sphere - thickness_sphere,
                             radius_of_sphere,
                             0.,
                             3.14 / 2.6);
      /// 3D transformation of mirrorVolFull in order to place it inside the gas volume
      Transform3D mirrorTr(RotationZYX(0, 0, 0), Translation3D(center_of_sphere_x, 0, center_of_sphere_z- mirror_z_safe_shrink));

      // TODO: cell 18 corresponds to half a pyramid, currently is full pyramid
      /// Define the actual mirror as intersection of the mother volume and the hollow sphere just defined
      Solid mirrorSol = IntersectionSolid(shape, mirrorShapeFull, mirrorTr);
      Volume mirrorVol(detName + "_mirror" + std::to_string(ncell) + "z" + std::to_string(reflect_parameters), mirrorSol, desc.material("Aluminum"));
      mirrorVol.setVisAttributes(desc.visAttributes(Form("mirror_vis%d", ncell)));

      // cellVol.placeVolume(mirrorVol, pyramidTr);

      // // Place detector in cell
      Transform3D sensorTr(RotationZYX(0, 90 * deg - angle_of_sensor, 0), Translation3D(-sensor_z_origin_Martin, 0, center_of_sensor_x));
      // cellVol.placeVolume(sensorVol, sensorTr);

      // position of mirror in cylinder coordinate system
      double mirror_abs_pos_z = name_col * zstep - 0.5 * zstep * (2 == name_row);
      if (reflect_parameters)
        mirror_abs_pos_z *= -1.0;

      // row 2 is shifted half step size
      double phi_offset = 0 + 0.5 * phistep * (2 == name_row);

      for (int phin = 0; phin < phinmax; ++phin)
      {
        auto cellTr = RotationZ(phistep * phin + phi_offset) * Translation3D(0, 0, mirror_abs_pos_z);
        PlacedVolume mirrorPV = gasVol.placeVolume(mirrorVol, cellTr * pyramidTr);
#ifdef __CREATE_COOLING__
        PlacedVolume sensorPV = coolingVol.placeVolume(sensorVol, cellTr * sensorTr);
#else
        PlacedVolume sensorPV = gasVol.placeVolume(sensorVol, cellTr * sensorTr);
#endif
        // sensorPV.addPhysVolID("phin", phin );
        // sensorPV.addPhysVolID("cellrow", name_row );
        // sensorPV.addPhysVolID("cellcolumn", name_col );

        // sensorPV.addPhysVolID("uniquecellid", ncell );
        // sensorPV.addPhysVolID("side", reflect_parameters?  0 : 1);
        sensorPV.addPhysVolID("cellnumber", 2*sensorcounter);

        std::cout << sensorPV.volIDs().str() << std::endl;

        // create mirrors as separate detectors, so properties can be adjusted later!
        DetElement mirrorDE(det, Form("ARC_DEmirror%d", sensorcounter), 2 * sensorcounter +1);
        mirrorDE.setPlacement(mirrorPV);
        SkinSurface mirrorSkin(desc, mirrorDE, Form("mirror_optical_surface%d", sensorcounter), mirrorSurf, mirrorVol); // FIXME: 3rd arg needs `imod`?
        mirrorSkin.isValid();

        // increase counter
        sensorcounter++;
      }

    } //---> End for loop over cell number vector

  } //---- End placing mirrors and sensors
#ifdef __CREATE_COOLING__
  gasVol.placeVolume(coolingVol);
#endif
  vesselVol.placeVolume(gasVol);
  // place vessel in mother volume
  // mother volume corresponds to the world
  Volume motherVol = desc.pickMotherVolume(det);
  PlacedVolume vesselPV = motherVol.placeVolume(vesselVol);
  vesselPV.addPhysVolID("system", detID);
  vesselPV.addPhysVolID("barrel", 0);
  // create mirrors as separate detectors, so properties can be adjusted lated!
  det.setPlacement(vesselPV);
  return det;
}

DECLARE_DETELEMENT(ARCBARREL_T, create_barrel)
