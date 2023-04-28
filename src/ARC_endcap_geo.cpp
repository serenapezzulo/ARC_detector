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

/// Function to build one ARC endcap
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

  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //          VESSEL PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  double vessel_outer_r = 190 * cm;
  double vessel_inner_r = 30.2 * cm;
  double vessel_length = 20 * cm;
  double vessel_wall_thickness = 1.0 * cm;
  if (vessel_outer_r <= vessel_inner_r)
    throw std::runtime_error("Ilegal parameters: vessel_outer_r <= vessel_inner_r");
  // // //-------------------------------------------------------------// // //


  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //         AEROGEL PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  double aerogel_thickness = 1.0 * cm;
  auto aerogelMat = desc.material("Aerogel_PFRICH");
  // // //-------------------------------------------------------------// // //

  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //         COOLING PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  double cooling_thickness = 1 * cm;
  // // //-------------------------------------------------------------// // //


  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //           CELL PARAMETERS           // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
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

// // // // // // // // // // // // // // // // // // // // // // //
//   SCHEME OF UNIQUE CELLS INSIDE A SECTOR
//   CELL NUMBERING CORRESPONDS TO ROGERS
//   MARTIN NUMBERING SPECIFIED BY ROW/COLUMN
//   THIS SECTOR MUST BE MIRRORED, AND THEN
//   BOTH (ORIGINAL AND MIRRORED) REPEATED 6 TIMES
//
//                         _____         _____
//                        /     \       /     \
//    7             _____/  21   \_____/  18   \
//                 /     \       /     \       /
//    7           /  20   \_____/  17   \_____/
//                \       /     \       /     \
//    6            \_____/  16   \_____/  14   \
//                 /     \       /     \       /
//    6           /  15   \_____/  13   \_____/
//                \       /     \       /     \
//    5            \_____/  12   \_____/  10   \
//                 /     \       /     \       /
//    5           /  11   \_____/   9   \_____/
//                \       /     \       /     \
//    4            \_____/   8   \_____/   7   \
//                       \       /     \       /
//    4                   \_____/   6   \_____/
//                        /     \       /     \
//    3                  /   5   \_____/   4   \
//                       \       /     \       /
//    3                   \_____/   3   \_____/
//                              \       /     \
//    2                          \_____/   2   \
//                               /     \       /
//    2                         /   1   \_____/
//                              \       /
//    COLUMN ^                   \_____/
//    ROW->       4        3        2     1
//
//   Y axis = column
//   X axis = row
// // // // // // // // // // // // // // // // // // // // // // //

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


  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //          MIRROR PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  double thickness_sphere(10 * mm);
  double mirror_z_origin_Martin = vessel_length / 2. - vessel_wall_thickness - 37 * cm;
  //   auto mirrorSurf = surfMgr.opticalSurface("MirrorSurface");
  auto mirrorElem = detElem.child(_Unicode(mirror)).child(_Unicode(module));
  double mirrorThickness = mirrorElem.attr<double>(_Unicode(thickness));
  auto mirrorSurf = surfMgr.opticalSurface(mirrorElem.attr<std::string>(_Unicode(surface)));
  auto mirrorMat = desc.material(mirrorElem.attr<std::string>(_Unicode(material)));
  // // //-------------------------------------------------------------// // //

  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // //          LIGHT SENSOR PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  double sensor_sidex = 8 * cm;
  double sensor_sidey = 8 * cm;
  double sensor_thickness = 0.2 * cm;
  double sensor_z_origin_Martin = -vessel_length / 2. + vessel_wall_thickness + 0.5 * cooling_thickness;
  // - sensor module
  auto sensorElem = detElem.child(_Unicode(sensors)).child(_Unicode(module));
  auto sensorVis = desc.visAttributes(sensorElem.attr<std::string>(_Unicode(vis)));
  //   double sensorX = sensorElem.attr<double>(_Unicode(sensorX));
  //   double sensorY = sensorElem.attr<double>(_Unicode(sensorY));
  //   double sensorThickness = sensorElem.attr<double>(_Unicode(thickness));
  auto sensorSurf = surfMgr.opticalSurface(sensorElem.attr<std::string>(_Unicode(surface)));
  auto sensorMat = desc.material(sensorElem.attr<std::string>(_Unicode(material)));

  // // //-------------------------------------------------------------// // //



  // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //
  // // //++++++++++++  BUILD VESSEL, CELL AND SENOR VOLUMES ++++++++++// // //
  // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //

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

/// Deprecated
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
