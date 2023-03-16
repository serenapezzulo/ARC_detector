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
 *  \details   This code creates a minimal working example of a RICH cell.
 *             Evolved from the pfRICH example in DD4hep.
 *  \author    Alvaro Tolosa-Delgado alvaro.tolosa.delgado@cern.ch
 *  \author    Martin Tat            martin.tat@cern.ch
 *  \version   0
 *  \date      2023
 *  \pre       DD4hep compiled with Geant4+Qt
 *  \bug       Walls do not reflect optical photons. Positions of unique cells in endcap are hardcoded.
 */

// Interface to Martins file inside anonymous namespace to avoid collisions
namespace
{
  // TODO: move parameters per cell to struct
  // store parameters by name
  std::map<std::string, double> cell_parameters_m;

  /// tokenize string by space as delimiter
  void mytokenizer(std::string &istring, std::vector<std::string> &tokens, char delimiter = ' ')
  {
    std::stringstream myline_ss(istring);
    std::string intermediate;
    while (getline(myline_ss, intermediate, delimiter))
      tokens.push_back(intermediate);
  }

  // function to fill map with parameters cell_parameters_m
  void fill_cell_parameters_m()
  {
    // avoid calling this function twice
    if ( not cell_parameters_m.empty() )
      return;

    // Read Martins File
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
      double parvalue = atof(tokens[1].c_str());
      cell_parameters_m.emplace( parname, parvalue );

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

/**
 * create barrel as sum of single cells
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
  int phinmax = 27; // 27;

  // Mirror parameters
  double thickness_sphere(10 * mm);
  double mirror_z_origin_Martin = vessel_outer_r - vessel_wall_thickness - 37 * cm;

  // Light sensor parameters
  double sensor_sidex = 8 * cm;
  double sensor_sidey = 8 * cm;
  double sensor_thickness = 0.2 * cm;
  double sensor_z_origin_Martin = vessel_inner_r + vessel_wall_thickness + 0.5 * cooling_radial_thickness;

  // Build cylinder for gas.
  Tube gasvolSolid(vessel_inner_r + vessel_wall_thickness,
                   vessel_outer_r - vessel_wall_thickness,
                   vessel_length);

  // Use regular polyhedra for endcaps cells
  // PolyhedraRegular cellS("aa",6,0,4*cm);

  // Use pyramid for barrel cells
  std::vector<double> zplanes = {0 * cm, vessel_outer_r - vessel_wall_thickness};
  std::vector<double> rs = {0 * cm, hexagon_side_length};
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
  // ncell_vector = {1};
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
    Transform3D mirrorTr(RotationZYX(0, 0, 0), Translation3D(center_of_sphere_x, 0, center_of_sphere_z));

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
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement det(detName, detID);
  sens.setType("tracker"); 

  auto sensorSurf = surfMgr.opticalSurface( "MirrorSurface");

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
  int phinmax = 27; // 27;

  // Mirror parameters
  double thickness_sphere(10 * mm);
  double mirror_z_origin_Martin = vessel_outer_r - vessel_wall_thickness - 37 * cm;

  // Light sensor parameters
  double sensor_sidex = 8 * cm;
  double sensor_sidey = 8 * cm;
  double sensor_thickness = 0.2 * cm;
  double sensor_z_origin_Martin = vessel_inner_r + vessel_wall_thickness + 0.5 * cooling_radial_thickness;
  
  // read Martin file and store parameters by name in the map
  fill_cell_parameters_m();

  ///----------->>> Define vessel and gas volumes
  /* - `vessel`: aluminum enclosure, mother volume of the barrel
   * - `gasvol`: gas volume, which fills `vessel`; all other volumes defined below
   *   as children of `gasvol`. Sensor is placed inside Cooling.
   * vessel (cylind) -> Gasvol (cylind) -> Cooling (cylind) -> Sensor CCD (cylind)
   *                                 \-> Mirror (sphere intersection with cell volume)
   *                                 \-> Aerogel (cylind)
   */
  // Build cylinder vol for vessel.
  Tube vesselSolid(vessel_inner_r,
                   vessel_outer_r,
                   vessel_length/2.);
  Volume vesselVol(detName + "_vessel", vesselSolid, desc.material("Aluminum"));
  vesselVol.setVisAttributes(desc.visAttributes( "gas_vis")); 

  // Build cylinder vol for gas. It is placed inside vessel
  Tube gasvolSolid(vessel_inner_r + vessel_wall_thickness,
                   vessel_outer_r - vessel_wall_thickness,
                   vessel_length/2.);
  Volume gasVol(detName + "_gas", gasvolSolid, desc.material("C4F10_PFRICH"));
  gasVol.setVisAttributes(desc.visAttributes("gas_vis"));

  // Build cylinder vol for cooling.It is placed inside gas
  Tube coolinSolid(vessel_inner_r + vessel_wall_thickness,
                   vessel_inner_r + vessel_wall_thickness + cooling_radial_thickness,
                   vessel_length/2.);
  Volume coolingVol(detName + "_cooling", coolinSolid, desc.material("Copper"));
  coolingVol.setVisAttributes(desc.visAttributes("gas_vis"));

  //----->> Place mirrors and sensors
  {
    
    // Use pyramid for barrel cells
    std::vector<double> zplanes = {0 * cm, vessel_outer_r - vessel_wall_thickness};
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
    // ncell_vector = {1};
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
                            3.14 / 2);
      /// 3D transformation of mirrorVolFull in order to place it inside the gas volume
      Transform3D mirrorTr(RotationZYX(0, 0, 0), Translation3D(center_of_sphere_x, 0, center_of_sphere_z));

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
        gasVol.placeVolume(mirrorVol, cellTr*pyramidTr);
        PlacedVolume sensorPV = coolingVol.placeVolume(sensorVol, cellTr*sensorTr);
        sensorPV.addPhysVolID("barrel", sensorcounter++ );
        // sensorPV.addPhysVolID("barrel", 1);
        // sensorPV.addPhysVolID("side", reflect_parameters?  0 : 1);
        // // sensorPV.addPhysVolID("sector", phin); //sector:5,
        // // sensorPV.addPhysVolID("uniquecell", ncell); // uniquecell:5,
        // sensorPV.addPhysVolID("isreflected", reflect_parameters);



        
        // create mirrors as separate detectors, so properties can be adjusted later!
      }


    } //---> End for loop over cell number vector


  } //---- End placing mirrors and sensors


  gasVol.placeVolume( coolingVol );
  vesselVol.placeVolume( gasVol );
  // place vessel in mother volume
  // mother volume corresponds to the world
  Volume motherVol = desc.pickMotherVolume(det);
  PlacedVolume vesselPV = motherVol.placeVolume( vesselVol );
  vesselPV.addPhysVolID("system", detID);
  // create mirrors as separate detectors, so properties can be adjusted lated!
  det.setPlacement(vesselPV);
  return det;





  // for (auto ncell : ncell_vector)
  // {
  //   // // The following line skips even number cells
  //   // // if (!(ncell % 2))
  //   // //   continue;

  //   // // The following line skips odd number cells
  //   // // if ((ncell % 2))
  //   // // continue;

  //   // /// cell shape. Coordinate system still the same as cylinder!
  //   // Solid cellS = IntersectionSolid(gasvolSolid, shape, pyramidTr);
  //   // Volume cellVol(detName + "_cell" + std::to_string(ncell), cellS, desc.material("C4F10_PFRICH"));
  //   // cellVol.setVisAttributes(desc.visAttributes("gas_vis"));

  //   // // there is no cell number 0, and cell number 1 do not need to be reflected
  //   // if (0 == ncell || -1 == ncell)
  //   //   continue;

  //   // // cells at z
  //   // bool reflect_parameters = false;
  //   // if (0 > ncell)
  //   // {
  //   //   ncell *= -1;
  //   //   reflect_parameters = true;
  //   // }

  //   // // initialize parameters for creating the mirror
  //   // double center_of_sphere_x(-999.);
  //   // double center_of_sphere_z(-999.);
  //   // double radius_of_sphere(-999.);

  //   // double center_of_sensor_x(-999.);
  //   // double angle_of_sensor(-999.);

  //   // // convert Roger nomenclature (one cell number) to Martin nomenclature (row and col numbers)
  //   // int name_col = ncell / 2;
  //   // int name_row = ncell % 2 ? 1 : 2;
  //   // // retrieve stored parameters
  //   // {
  //   //   std::string name_col_s = std::to_string(name_col);
  //   //   std::string name_row_s = std::to_string(name_row);
  //   //   radius_of_sphere = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_Curvature"];
  //   //   center_of_sphere_x = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_XPosition"];
  //   //   double zposition = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_ZPosition"];

  //   //   center_of_sensor_x = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_DetPosition"];
  //   //   angle_of_sensor = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_DetTilt"];
  //   //   center_of_sphere_z = mirror_z_origin_Martin + zposition;

  //   //   // check if parameters are ok
  //   //   if (-999. == center_of_sphere_x)
  //   //     throw std::runtime_error("Ilegal parameters: center_of_sphere_x not provided");
  //   //   if (-999. == center_of_sphere_z)
  //   //     throw std::runtime_error("Ilegal parameters: center_of_sphere_z not provided");
  //   //   if (-999. == radius_of_sphere)
  //   //     throw std::runtime_error("Ilegal parameters: radius_of_sphere not provided");
  //   //   if (radius_of_sphere <= thickness_sphere)
  //   //     throw std::runtime_error("Ilegal parameters: radius_of_sphere <= thickness_sphere");

  //   //   if (-999. == center_of_sensor_x)
  //   //     throw std::runtime_error("Ilegal parameters: center_of_sensor_x not provided");
  //   //   if (-999. == angle_of_sensor)
  //   //     throw std::runtime_error("Ilegal parameters: angle_of_sensor not provided");
  //   // }

  //   // if (reflect_parameters)
  //   // {
  //   //   center_of_sphere_x *= -1.0;
  //   //   center_of_sensor_x *= -1.0;
  //   //   angle_of_sensor *= -1.0;
  //   // }

  //   // create the semi-sphere that will result in the mirror
  //   Sphere mirrorShapeFull(radius_of_sphere - thickness_sphere,
  //                          radius_of_sphere,
  //                          0.,
  //                          3.14 / 2);
  //   /// 3D transformation of mirrorVolFull in order to place it inside the gas volume
  //   Transform3D mirrorTr(RotationZYX(0, 0, 0), Translation3D(center_of_sphere_x, 0, center_of_sphere_z));

  //   // TODO: cell 18 corresponds to half a pyramid, currently is full pyramid
  //   /// Define the actual mirror as intersection of the mother volume and the hollow sphere just defined
  //   Solid mirrorSol = IntersectionSolid(shape, mirrorShapeFull, mirrorTr);
  //   Volume mirrorVol(detName + "_mirror" + std::to_string(ncell) + "z" + std::to_string(reflect_parameters), mirrorSol, desc.material("Aluminum"));
  //   mirrorVol.setVisAttributes(desc.visAttributes(Form("mirror_vis%d", ncell)));
  //   cellVol.placeVolume(mirrorVol, pyramidTr);

  //   // Place detector in cell
  //   Transform3D sensorTr(RotationZYX(0, 90 * deg - angle_of_sensor, 0), Translation3D(-sensor_z_origin_Martin, 0, center_of_sensor_x));
  //   cellVol.placeVolume(sensorVol, sensorTr);

  //   // position of mirror in cylinder coordinate system
  //   double mirror_abs_pos_z = name_col * zstep - 0.5 * zstep * (2 == name_row);
  //   if (reflect_parameters)
  //     mirror_abs_pos_z *= -1.0;

  //   // row 2 is shifted half step size
  //   double phi_offset = 0 + 0.5 * phistep * (2 == name_row);

  //   for (int phin = 0; phin < phinmax; ++phin)
  //   {
  //     PlacedVolume cellPV = motherVol.placeVolume(cellVol, RotationZ(phistep * phin + phi_offset) * Translation3D(0, 0, mirror_abs_pos_z));
  //     cellPV.addPhysVolID("system", detID).addPhysVolID("module", 17 * phin + name_col);
  //     // create mirrors as separate detectors, so properties can be adjusted lated!
  //     det.setPlacement(cellPV);
  //   }
  // }

  return det;
}

DECLARE_DETELEMENT(ARCBARREL_T, create_barrel)

/**
 * create endcap as sum of single cells
 * next step is to place mirrors+sensors in a cylindral shape gas volume
 */
static Ref_t create_endcap_cell(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
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

  // Vessel, endcap
  double vessel_outer_r = 190 * cm;
  double vessel_inner_r = 30.2 * cm;
  double vessel_length = 20 * cm;
  double vessel_wall_thickness = 1.0 * cm;
  if (vessel_outer_r <= vessel_inner_r)
    throw std::runtime_error("Ilegal parameters: vessel_outer_r <= vessel_inner_r");

  // Cooling,
  double cooling_thickness = 1 * cm;

  // Cell parameters
  /// Cell is an hexagonal prysm
  double hexagon_side_length = 14.815 * cm;
  double hexagon_apothem = hexagon_side_length * cos(30 * deg);
  /// each cell corresponds to one object of the following class
  /// which gathers all the important parameters
  /// this info can be moved later to compact file
  struct mycell_t
  {
    /// Martin number for row
    int row = {0};
    /// Martin number for column
    int col = {0};
    /// Roger ID number
    int RID = {-1};
    /// x position of the cell
    double x = {0.};
    /// y position of the cell
    double y = {0.};
    /// if reflected
    bool isReflected = {false};
  };

  /// vector with cell geometric parameters
  std::vector<mycell_t> mycell_v(21);
  {
    double hx_u = hexagon_apothem;
    double hx_x = hexagon_side_length;
    mycell_v[0] = {1, 2, 2, 0, 4 * hx_u};
    mycell_v[1] = {1, 3, 4, 0, 6 * hx_u};
    mycell_v[2] = {1, 4, 7, 0, 8 * hx_u};
    mycell_v[3] = {1, 5, 10, 0, 10 * hx_u};
    mycell_v[4] = {1, 6, 14, 0, 12 * hx_u};
    mycell_v[5] = {1, 7, 18, 0, 14 * hx_u};
    mycell_v[6] = {2, 2, 1, -1.5 * hx_x, 3 * hx_u};
    mycell_v[7] = {2, 3, 3, -1.5 * hx_x, 5 * hx_u, true};
    mycell_v[8] = {2, 4, 6, -1.5 * hx_x, 7 * hx_u, true};
    mycell_v[9] = {2, 5, 9, -1.5 * hx_x, 9 * hx_u, true};
    mycell_v[10] = {2, 6, 13, -1.5 * hx_x, 11 * hx_u, true};
    mycell_v[11] = {2, 7, 17, -1.5 * hx_x, 13 * hx_u, true};
    mycell_v[12] = {3, 3, 5, -3.0 * hx_x, 6 * hx_u};
    mycell_v[13] = {3, 4, 8, -3.0 * hx_x, 8 * hx_u, true};
    mycell_v[14] = {3, 5, 12, -3.0 * hx_x, 10 * hx_u, true};
    mycell_v[15] = {3, 6, 16, -3.0 * hx_x, 12 * hx_u, true};
    mycell_v[16] = {3, 7, 21, -3.0 * hx_x, 14 * hx_u, true};
    mycell_v[17] = {4, 5, 11, -4.5 * hx_x, 9 * hx_u};
    mycell_v[18] = {4, 6, 15, -4.5 * hx_x, 11 * hx_u, true};
    mycell_v[19] = {4, 7, 20, -4.5 * hx_x, 13 * hx_u, true};
    mycell_v[20] = {5, 6, 19, -6.0 * hx_x, 12 * hx_u};
  }

  /// Distance in phi angle between cells
  /// since cells are regular hexagons, this distance matches
  /// the angle that one cell covers
  double phistep = 60 * deg;
  /// number of repetition of unique cells around the endcap
  int phinmax = 6; // 6;

  // Mirror parameters
  double thickness_sphere(10 * mm);
  double mirror_z_origin_Martin = vessel_length / 2. - vessel_wall_thickness - 37 * cm;

  // Light sensor parameters
  double sensor_sidex = 8 * cm;
  double sensor_sidey = 8 * cm;
  double sensor_thickness = 0.2 * cm;
  double sensor_z_origin_Martin = -vessel_length / 2. + vessel_wall_thickness + 0.5 * cooling_thickness;

  // Build cylinder for gas.
  Tube gasvolSolid(vessel_inner_r + vessel_wall_thickness,
                   vessel_outer_r - vessel_wall_thickness,
                   vessel_length);

  // Use regular polyhedra for endcaps cells
  PolyhedraRegular cellS(6, 0, 0., hexagon_apothem, vessel_length);

  // Build sensor shape
  Box sensorSol(sensor_sidex / 2, sensor_sidey / 2, sensor_thickness / 2);
  Volume sensorVol(detName + "_sensor", sensorSol, desc.material("Aluminum"));

  // Build the mirror for ncell=1..21
  // auto ncell = mycell_v[0];
  for (auto &ncell : mycell_v)
  {
    if (-1 == ncell.RID)
      continue;
    // if( 5 != ncell.row )
    // continue;
    for (int phin = 0; phin < phinmax; phin++)
    {

      std::string volname = detName + "_cell" + std::to_string(ncell.RID);
      volname += "_phin" + std::to_string(phin);
      Volume cellV(volname, cellS, desc.material("C4F10_PFRICH"));

      // The following line skips even number cells
      // if ( 1 != ncell.row )
      //   continue;

      // // initialize parameters for creating the mirror
      double center_of_sphere_x(-999.);
      double center_of_sphere_z(-999.);
      double radius_of_sphere(-999.);

      double center_of_sensor_x(-999.);
      double angle_of_sensor(-999.);

      // convert Roger nomenclature (one cell number) to Martin nomenclature (row and col numbers)
      int name_col = ncell.col;
      int name_row = ncell.row;
      // retrieve stored parameters
      {
        std::string name_col_s = std::to_string(name_col);
        std::string name_row_s = std::to_string(name_row);
        radius_of_sphere = cell_parameters_m.at("EndCapRadiator_c" + name_col_s + "_r" + name_row_s + "_Curvature");
        center_of_sphere_x = cell_parameters_m.at("EndCapRadiator_c" + name_col_s + "_r" + name_row_s + "_XPosition");
        double zposition = cell_parameters_m.at("EndCapRadiator_c" + name_col_s + "_r" + name_row_s + "_ZPosition");
        center_of_sphere_z = mirror_z_origin_Martin + zposition;

        center_of_sensor_x = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_DetPosition"];
        angle_of_sensor = cell_parameters_m["Radiator_c" + name_col_s + "_r" + name_row_s + "_DetTilt"];

        // check if parameters are ok
        if (-999. == center_of_sphere_x)
          throw std::runtime_error("Ilegal parameters: center_of_sphere_x not provided");
        if (-999. == center_of_sphere_z)
          throw std::runtime_error("Ilegal parameters: center_of_sphere_z not provided");
        if (-999. == radius_of_sphere)
          throw std::runtime_error("Ilegal parameters: radius_of_sphere not provided");
        if (radius_of_sphere <= thickness_sphere)
          throw std::runtime_error(Form("Ilegal parameters cell %d: %g <= %g", ncell.RID, radius_of_sphere, thickness_sphere));

        if (-999. == center_of_sensor_x)
          throw std::runtime_error("Ilegal parameters: center_of_sensor_x not provided");
        if (-999. == angle_of_sensor)
          throw std::runtime_error("Ilegal parameters: angle_of_sensor not provided");
      }

      // create the semi-sphere that will result in the mirror
      Sphere mirrorShapeFull(radius_of_sphere - thickness_sphere,
                             radius_of_sphere,
                             0.,
                             3.14 / 2);
      /// 3D transformation of mirrorVolFull in order to place it inside the gas volume
      double alpha = atan(ncell.y / ncell.x) * rad;
      if (0 > alpha)
        alpha += 180 * deg;
      double dx = center_of_sphere_x * cos(alpha);
      double dy = center_of_sphere_x * sin(alpha);
      std::cout << ncell.RID << '\t' << center_of_sphere_x << '\t' << alpha / rad * deg << std::endl;

      Transform3D mirrorTr(RotationZYX(0, 0, 0), Translation3D(dx, dy, center_of_sphere_z));

      /// Define the actual mirror as intersection of the hex cell volume and the hollow sphere just defined
      Solid mirrorSol = IntersectionSolid(cellS, mirrorShapeFull, mirrorTr);
      std::string mirrorVolName = detName + "_mirror" + std::to_string(ncell.RID) + "z" + std::to_string(ncell.isReflected);
      Volume mirrorVol(mirrorVolName, mirrorSol, desc.material("Aluminum"));
      mirrorVol.setVisAttributes(desc.visAttributes(Form("mirror_vis%d", ncell.RID)));
      cellV.placeVolume(mirrorVol);

      // // Place detector in cell
      Transform3D sensorTr(RotationZYX(alpha - 90 * deg, 0 /*90*deg-angle_of_sensor*/, angle_of_sensor /*ncell.row*20*deg*/),
                           Translation3D(0, center_of_sensor_x, sensor_z_origin_Martin));
      cellV.placeVolume(sensorVol, sensorTr);

      // // position of mirror in cylinder coordinate system
      // double mirror_abs_pos_z = name_col * zstep - 0.5 * zstep * (2 == name_row);
      // if (reflect_parameters)
      //   mirror_abs_pos_z *= -1.0;

      // // row 2 is shifted half step size
      // double phi_offset = 0 + 0.5 * phistep * (2 == name_row);

      // for (int phin = 0; phin < phinmax; ++phin)
      // {
      //   PlacedVolume cellPV = motherVol.placeVolume(cellVol, RotationZ(phistep * phin + phi_offset) * Translation3D(0, 0, mirror_abs_pos_z));
      //   cellPV.addPhysVolID("system", detID).addPhysVolID("module", 17 * phin + name_col);
      //   // create mirrors as separate detectors, so properties can be adjusted lated!
      //   det.setPlacement(cellPV);
      // }

      cellV.setVisAttributes(desc.visAttributes("gas_vis"));
      PlacedVolume cellPV = motherVol.placeVolume(cellV, RotationZ(phistep * phin) * Translation3D(ncell.x, ncell.y, 0));
      cellPV.addPhysVolID("system", detID).addPhysVolID("module", ncell.RID);
      // create mirrors as separate detectors, so properties can be adjusted lated!
      det.setPlacement(cellPV);
      if (ncell.isReflected)
      {
        Volume cellV_reflected(volname + "_z1", cellS, desc.material("C4F10_PFRICH"));
        cellV_reflected.setVisAttributes(desc.visAttributes("gas_vis"));
        Transform3D mirrorTr_reflected(RotationZYX(0, 0, 0), Translation3D(-dx, dy, center_of_sphere_z));

        // TODO: cell 18 corresponds to half a pyramid, currently is full pyramid
        /// Define the actual mirror as intersection of the mother volume and the hollow sphere just defined
        Solid mirrorSol_reflected = IntersectionSolid(cellS, mirrorShapeFull, mirrorTr_reflected);
        std::string mirror_name = detName + "_mirror" + std::to_string(ncell.RID) + "z" + std::to_string(ncell.isReflected);
        Volume mirrorVol_reflected(mirror_name, mirrorSol_reflected, desc.material("Aluminum"));
        mirrorVol_reflected.setVisAttributes(desc.visAttributes(Form("mirror_vis%d", ncell.RID)));
        cellV_reflected.placeVolume(mirrorVol_reflected);
        Transform3D sensorTr_reflected(RotationZYX(-alpha + 90 * deg, 0 /*90*deg-angle_of_sensor*/, angle_of_sensor),
                                       Translation3D(0, center_of_sensor_x, sensor_z_origin_Martin));
        cellV_reflected.placeVolume(sensorVol, sensorTr_reflected);

        motherVol.placeVolume(cellV_reflected, RotationZ(phistep * phin) * Translation3D(-ncell.x, ncell.y, 0));
      }
    } //-- end loop for sector
  }   //-- end loop for endcap

  return det;
}
DECLARE_DETELEMENT(ARCENDCAP_T, create_endcap_cell)

/// This fcn just build the individual cell volume, without elements
static Ref_t create_endcap_cell_volumes(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
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

  // Vessel, endcap
  double vessel_outer_r = 190 * cm;
  double vessel_inner_r = 30.2 * cm;
  double vessel_length = 20 * cm;
  double vessel_wall_thickness = 1.0 * cm;
  if (vessel_outer_r <= vessel_inner_r)
    throw std::runtime_error("Ilegal parameters: vessel_outer_r <= vessel_inner_r");

  // Cooling,
  double cooling_thickness = 1 * cm;

  // Cell parameters
  /// Cell is an hexagonal prysm
  double hexagon_side_length = 14.815 * cm;
  double hexagon_apothem = hexagon_side_length * cos(30 * deg);
  /// each cell corresponds to one object of the following class
  /// which gathers all the important parameters
  /// this info can be moved later to compact file
  struct mycell_t
  {
    /// Martin number for row
    int row = {0};
    /// Martin number for column
    int col = {0};
    /// Roger ID number
    int RID = {0};
    /// x position of the cell
    double x = {0.};
    /// y position of the cell
    double y = {0.};
    /// if reflected
    bool isReflected = {false};
  };

  /// vector with cell geometric parameters
  std::vector<mycell_t> mycell_v(21);
  {
    double hx_u = hexagon_apothem;
    double hx_x = hexagon_side_length;
    mycell_v[0] = {1, 1, 2, 0, 4 * hx_u};
    mycell_v[1] = {1, 2, 4, 0, 6 * hx_u};
    mycell_v[2] = {1, 3, 7, 0, 8 * hx_u};
    mycell_v[3] = {1, 4, 10, 0, 10 * hx_u};
    mycell_v[4] = {1, 5, 14, 0, 12 * hx_u};
    mycell_v[5] = {1, 6, 18, 0, 14 * hx_u};
    mycell_v[6] = {2, 1, 1, -1.5 * hx_x, 3 * hx_u};
    mycell_v[7] = {2, 2, 3, -1.5 * hx_x, 5 * hx_u, true};
    mycell_v[8] = {2, 3, 6, -1.5 * hx_x, 7 * hx_u, true};
    mycell_v[9] = {2, 4, 9, -1.5 * hx_x, 9 * hx_u, true};
    mycell_v[10] = {2, 5, 13, -1.5 * hx_x, 11 * hx_u, true};
    mycell_v[11] = {2, 6, 17, -1.5 * hx_x, 13 * hx_u, true};
    mycell_v[12] = {3, 1, 5, -3.0 * hx_x, 6 * hx_u};
    mycell_v[13] = {3, 2, 8, -3.0 * hx_x, 8 * hx_u, true};
    mycell_v[14] = {3, 3, 12, -3.0 * hx_x, 10 * hx_u, true};
    mycell_v[15] = {3, 4, 16, -3.0 * hx_x, 12 * hx_u, true};
    mycell_v[16] = {3, 5, 21, -3.0 * hx_x, 14 * hx_u, true};
    mycell_v[17] = {4, 1, 11, -4.5 * hx_x, 9 * hx_u};
    mycell_v[18] = {4, 2, 15, -4.5 * hx_x, 11 * hx_u, true};
    mycell_v[19] = {4, 3, 20, -4.5 * hx_x, 13 * hx_u, true};
    mycell_v[20] = {5, 1, 19, -6.0 * hx_x, 12 * hx_u};
  }

  /// Distance in phi angle between cells
  /// since cells are regular hexagons, this distance matches
  /// the angle that one cell covers
  double phistep = 60 * deg;
  /// number of repetition of unique cells around the endcap
  int phinmax = 6; // 27;

  // Mirror parameters
  double thickness_sphere(10 * mm);
  double mirror_z_origin_Martin = vessel_length / 2. - vessel_wall_thickness - 37 * cm;

  // Light sensor parameters
  double sensor_sidex = 8 * cm;
  double sensor_sidey = 8 * cm;
  double sensor_thickness = 0.2 * cm;
  double sensor_z_origin_Martin = -vessel_length / 2. + vessel_wall_thickness + 0.5 * cooling_thickness;

  // Build cylinder for gas.
  Tube gasvolSolid(vessel_inner_r + vessel_wall_thickness,
                   vessel_outer_r - vessel_wall_thickness,
                   vessel_length);

  // Use regular polyhedra for endcaps cells
  PolyhedraRegular cellS(6, 0, 0., hexagon_apothem, vessel_length);

  // Build the mirror for ncell=1..21

  for (auto ncell : mycell_v)
  {
    for (int phin = 0; phin < 6; phin++)
    {

      std::string volname = detName + "_cell" + std::to_string(ncell.RID);
      volname += "_phin" + std::to_string(phin);
      Volume cellV(volname, cellS, desc.material("Aluminum"));
      cellV.setVisAttributes(desc.visAttributes(Form("mirror_vis%d", ncell.RID)));
      Transform3D cellTr(RotationZ(60 * phin * deg), Translation3D(ncell.x, ncell.y, 0));
      PlacedVolume cellPV = motherVol.placeVolume(cellV, RotationZ(phistep * phin) * Translation3D(ncell.x, ncell.y, 0));
      cellPV.addPhysVolID("system", detID).addPhysVolID("module", ncell.RID);
      // create mirrors as separate detectors, so properties can be adjusted lated!
      det.setPlacement(cellPV);
      if (ncell.isReflected)
      {
        Transform3D cellTrReflected(RotationZ(60 * phin * deg), Translation3D(-ncell.x, ncell.y, 0));
        motherVol.placeVolume(cellV, RotationZ(phistep * phin) * Translation3D(-ncell.x, ncell.y, 0));
      }

    } //-- end loop for sector
  }   //-- end loop for endcap

  return det;
}
