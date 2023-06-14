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

using namespace dd4hep;

// #define DUMP_SENSOR_POSITIONS

/// Function to build one ARC endcap
static Ref_t create_endcap_cell(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement det(detName, detID);
  sens.setType("tracker");

  double zpos_endcap = detElem.attr<double>(_Unicode(zpos));

  auto gasElem    = detElem.child(_Unicode(radiatorgas));
  auto gasvolMat  = desc.material(gasElem.attr<std::string>(_Unicode(material)));
  auto gasvolVis  = desc.visAttributes(gasElem.attr<std::string>(_Unicode(vis)));

  auto vesselElem = detElem.child(_Unicode(vessel));
  auto vesselSkinMat  = desc.material(vesselElem.attr<std::string>(_Unicode(skinMaterial)));
  auto vesselSkinVis  = desc.visAttributes(vesselElem.attr<std::string>(_Unicode(skin_vis)));

  auto vesselBulkMat  = desc.material(vesselElem.attr<std::string>(_Unicode(bulk_material)));
  auto vesselBulkVis  = desc.visAttributes(vesselElem.attr<std::string>(_Unicode(bulk_vis)));

    double bulk_skin_ratio = vesselElem.attr<double>(_Unicode(bulk_skin_ratio));

    if( 0 > bulk_skin_ratio || 1 < bulk_skin_ratio )
        throw std::runtime_error("ARC: bulk_skin_ratio must be a number between 0 and 1");


  // read Martin file and store parameters by name in the map
//   fill_cell_parameters_m();

  // mother volume corresponds to the world
  Volume motherVol = desc.pickMotherVolume(det);

  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //          VESSEL PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  double vessel_outer_r = 210 * cm; // 190 * cm;
  double vessel_inner_r = 25 * cm;  // 30.2 * cm;
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
  double cooling_thickness = 2 * mm;
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
//                        /     \       /     \    .
//    7             _____/  21   \_____/  18   \   .
//                 /     \       /     \       /   .
//    7           /  20   \_____/  17   \_____/    .
//                \       /     \       /     \    .
//    6            \_____/  16   \_____/  14   \   .
//                 /     \       /     \       /   .
//    6           /  15   \_____/  13   \_____/    .
//                \       /     \       /     \    .
//    5            \_____/  12   \_____/  10   \   .
//                 /     \       /     \       /   .
//    5           /  11   \_____/   9   \_____/    .
//                \       /     \       /     \    .
//    4            \_____/   8   \_____/   7   \   .
//                       \       /     \       /   .
//    4                   \_____/   6   \_____/    .
//                        /     \       /     \    .
//    3                  /   5   \_____/   4   \   .
//                       \       /     \       /   .
//    3                   \_____/   3   \_____/    .
//                              \       /     \    .
//    2                          \_____/   2   \   .
//                               /     \       /   .
//    2                         /   1   \_____/    .
//                              \       /          .
//    COLUMN ^                   \_____/           .
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

  /// Distance in phi angle between complete sectors
  double phistep = 60 * deg;
  /// number of repetition of sectors
  int phinmax = 6; // 6;


  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //          MIRROR PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  double mirror_z_origin_Martin = vessel_length / 2. - vessel_wall_thickness - 37 * cm;
  auto mirrorElem = detElem.child(_Unicode(mirror));
  double mirrorThickness = mirrorElem.attr<double>(_Unicode(thickness));
  auto mirrorSurf = surfMgr.opticalSurface(mirrorElem.attr<std::string>(_Unicode(surface)));
  auto mirrorMat = desc.material(mirrorElem.attr<std::string>(_Unicode(material)));
  // // //-------------------------------------------------------------// // //


    // // // // // // // // // // // // // // // // // // // // // // // // // //
    // // // // // //          LIGHT SENSOR PARAMETERS          // // // // // //
    // // // // // // // // // // // // // // // // // // // // // // // // // //
    //default values
    double sensor_sidex     = 8 * cm;
    double sensor_sidey     = 8 * cm;
    double sensor_thickness = 0.2 * cm;
    // empirical distance to keep the sensor inside the cell volume
    double sensor_z_origin_Martin = -vessel_length / 2. + vessel_wall_thickness + 0.5 * cooling_thickness;
    auto sensorMat = desc.material("SiliconOptical");
    auto sensorVis = desc.visAttributes("no_vis");

    // Read from xml the parameters for the sensor module
    {
        auto sensorElem  = detElem.child(_Unicode(sensors));
        sensor_sidex     = sensorElem.attr<double>(_Unicode(sensor_side_X));
        sensor_sidey     = sensorElem.attr<double>(_Unicode(sensor_side_Y));
        sensor_thickness = sensorElem.attr<double>(_Unicode(thickness));
        sensorMat        = desc.material(sensorElem.attr<std::string>(_Unicode(material)));
        sensorVis        = desc.visAttributes(sensorElem.attr<std::string>(_Unicode(vis)));
    }
    // // //-------------------------------------------------------------// // //


    // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //
    // // //+++++++++++  BUILD VESSEL, CELL AND SENSOR VOLUMES ++++++++++// // //
    // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //

    // Build cylinder for gas, and the vessel for the gas
    Tube gasenvelopeS(  vessel_inner_r + vessel_wall_thickness,
                        vessel_outer_r - vessel_wall_thickness,
                        vessel_length/2.);
    Volume endcap_cells_gas_envelope (detName+"_gasEnvelope", gasenvelopeS, gasvolMat );
    endcap_cells_gas_envelope.setVisAttributes( desc.visAttributes("arc_envelope_vis") );

    Tube vesselEnvelopeSolid(  vessel_inner_r,
                               vessel_outer_r,
                               vessel_length/2. + vessel_wall_thickness);
    Volume endcap_cells_vessel_envelope (detName+"_vesselEnvelope", vesselEnvelopeSolid, vesselSkinMat );
    endcap_cells_vessel_envelope.setVisAttributes( vesselSkinVis );

    // if 0==bulk_skin_ratio do not create bulk at all
    if(0<bulk_skin_ratio)
    {
      // build bulk for inner wall
      double vessel_bulk_inner_r_ini = vessel_inner_r + (1 - bulk_skin_ratio)*0.5*vessel_wall_thickness;
      double vessel_bulk_inner_r_fin = vessel_inner_r + (1 + bulk_skin_ratio)*0.5*vessel_wall_thickness;

      Tube vesselInnerBulkSolid( vessel_bulk_inner_r_ini,
                            vessel_bulk_inner_r_fin,
                            vessel_length/2. + vessel_wall_thickness - (1-bulk_skin_ratio)*0.5*vessel_wall_thickness);
      Volume vessel_innerbulk_vol (detName+"_vesselInnerBulk", vesselInnerBulkSolid, vesselBulkMat );
      vessel_innerbulk_vol.setVisAttributes( vesselBulkVis );
      endcap_cells_vessel_envelope.placeVolume(vessel_innerbulk_vol);

      // build bulk for outer wall
      double vessel_bulk_outer_r_ini = vessel_outer_r - (1 + bulk_skin_ratio)*0.5*vessel_wall_thickness;
      double vessel_bulk_outer_r_fin = vessel_outer_r - (1 - bulk_skin_ratio)*0.5*vessel_wall_thickness;

      Tube vesselOuterBulkSolid( vessel_bulk_outer_r_ini,
                                vessel_bulk_outer_r_fin,
                                vessel_length/2. + vessel_wall_thickness -  (1-bulk_skin_ratio)*0.5*vessel_wall_thickness);
      Volume vessel_outerbulk_vol (detName+"_vesselOuterBulk", vesselOuterBulkSolid, vesselBulkMat );
      vessel_outerbulk_vol.setVisAttributes( vesselBulkVis );
      endcap_cells_vessel_envelope.placeVolume(vessel_outerbulk_vol);

      Tube vesselBaseBulkSolid(  vessel_bulk_inner_r_fin,
                                vessel_bulk_outer_r_ini,
                                bulk_skin_ratio*0.5*vessel_wall_thickness);
      Volume vessel_base_bulk_vol (detName+"_vesselBaseBulk", vesselBaseBulkSolid, vesselBulkMat );
      vessel_base_bulk_vol.setVisAttributes( vesselBulkVis );
      auto posZPositive = Position(0, 0, vessel_length/2. + 0.5*vessel_wall_thickness);
      endcap_cells_vessel_envelope.placeVolume(vessel_base_bulk_vol,posZPositive);

      auto posZNegative = Position(0, 0, -vessel_length/2. - 0.5*vessel_wall_thickness);
      endcap_cells_vessel_envelope.placeVolume(vessel_base_bulk_vol,posZNegative);
    }


    // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //

  // Use regular polyhedra for endcaps cells
  PolyhedraRegular cellS(6, 0, 0., hexagon_apothem, vessel_length);

  // Build sensor shape
  Box sensorSol(sensor_sidex / 2, sensor_sidey / 2, sensor_thickness / 2);
  Volume sensorVol(detName + "_sensor", sensorSol, sensorMat);
  sensorVol.setSensitiveDetector(sens);
  sensorVol.setVisAttributes( sensorVis );

  // Build cooling plate
  double cooling_z_offset =   sensor_thickness  + cooling_thickness + 0.5*mm;
  Tube coolingSol_tube(0, 1.5*hexagon_side_length, cooling_thickness);

  // Build aerogel plate
  double aerogel_z_offset =   sensor_thickness  + aerogel_thickness + 0.5*mm;
  Tube aerogelSol_tube(0, 1.5*hexagon_side_length, aerogel_thickness);

  // Build cells of a sector
//   mycell_v = {mycell_v[16], mycell_v[19]};
  phinmax = 1;
  int cellCounter = 0;
  int physicalVolumeCounter = 0;
  auto createPhysVolID = [&](){return physicalVolumeCounter++;};

#ifdef DUMP_SENSOR_POSITIONS
  std::ofstream ofile_sensor_pos("ofile_sensor_pos_endcap.txt");
#endif
  for (auto &ncell : mycell_v)
  {
    // sanity check, skip non initialized cells
    if (-1 == ncell.RID)
      continue;

    // The following line skips even number cells
    // if ( 1 != ncell.row )
    //   continue;

    /// repeat the sector 6 times
    for (int phin = 0; phin < phinmax; phin++, cellCounter++)
    {

      /// function to create names in a systematic way
      /// final name = detName + part + cell parameters
      auto create_part_name_ff = [ncell,detName,phin](std::string  partName){
          std::string fullName = detName + "_" + partName;
          fullName += std::to_string(ncell.RID);
          fullName += "_phi" +  std::to_string(phin);
          dd4hep::printout(dd4hep::DEBUG,"ARCENDCAP_T", "+++ New name:%s",fullName.c_str());
          return fullName;
        };

      /// cell volume, hex prism
      /// the elements must be placed inside
      std::string cellName = create_part_name_ff("cell");
      Volume cellV(cellName, cellS, gasvolMat);
      cellV.setVisAttributes( gasvolVis );
      /// Detector element that will contain cellVol later
      /// there are 3 elements with ID:
      /// the cell, ID= 6 * cellCounter
      /// its mirror. ID = 6 * cellCounter +1
      /// and its sensor, ID = 6 * cellCounter +2
      DetElement cellDE(det, cellName+"DE", 6 * cellCounter + 0);



      // // initialize parameters for creating the mirror
      double center_of_sphere_x(-999.);
      double center_of_sphere_z(-999.);
      double radius_of_sphere(-999.);

      double center_of_sensor_x(-999.);
      double angle_of_sensor(-999.);
      double zoffset_of_sensor(0);


      // retrieve cell parameters
      // if parameter not present, exception is thrown and not catched
      {
        // convert Roger nomenclature (one cell number) to Martin nomenclature (row and col numbers)
        std::string name_col_s = std::to_string(ncell.col);
        std::string name_row_s = std::to_string(ncell.row);
        std::string MartinCellName = "EndCapRadiator_c" + name_col_s + "_r" + name_row_s;

        radius_of_sphere = desc.constantAsDouble(MartinCellName + "_Curvature");

        center_of_sphere_x = desc.constantAsDouble(MartinCellName + "_XPosition");

        double zposition = desc.constantAsDouble(MartinCellName + "_ZPosition");

        center_of_sphere_z = mirror_z_origin_Martin + zposition;

        center_of_sensor_x = desc.constantAsDouble(MartinCellName + "_DetPosition");

        angle_of_sensor = desc.constantAsDouble(MartinCellName + "_DetTilt");

        std::string ZOffsetSensorParName = MartinCellName + "_DetZOffset";


        if( desc.constants().count(MartinCellName + "_DetPositionZ") )
            zoffset_of_sensor = desc.constantAsDouble(MartinCellName + "_DetPositionZ");
        else
          dd4hep::printout(dd4hep::WARNING,"ARCENDCAP_T", "+++ Constant %s is missing in xml file, default is 0",ZOffsetSensorParName.c_str());

        if (radius_of_sphere <= mirrorThickness)
          throw std::runtime_error(Form("Ilegal parameters cell %d: %g <= %g", ncell.RID, radius_of_sphere, mirrorThickness));

      }
      double sensor_z_pos = zoffset_of_sensor + sensor_z_origin_Martin;

      // create the semi-sphere that will result in the mirror
      Sphere mirrorShapeFull(radius_of_sphere - mirrorThickness,
                             radius_of_sphere,
                             0.,
                             3.14 / 3);
      /// alpha: angle of position vector of first sector n-cell with respect to x-axis
      double alpha = atan(ncell.y / ncell.x) * rad;
      if (0 > alpha)
        alpha += 180 * deg;
      double dx = center_of_sphere_x * cos(alpha);
      double dy = center_of_sphere_x * sin(alpha);

      /// 3D transformation of mirrorVolFull in order to place it inside the gas volume
      Transform3D mirrorTr(RotationZYX(0, 0, 0), Translation3D(dx, dy, center_of_sphere_z));

      /// Define the actual mirror as intersection of the hex cell volume and the hollow sphere just defined
      Solid mirrorSol = IntersectionSolid(cellS, mirrorShapeFull, mirrorTr);
      std::string mirrorVolName = create_part_name_ff("mirror");
      Volume mirrorVol(mirrorVolName, mirrorSol, mirrorMat);
      mirrorVol.setVisAttributes(desc.visAttributes(Form("arc_mirror_vis%d", ncell.RID)));
      PlacedVolume mirrorPV = cellV.placeVolume(mirrorVol);

      DetElement mirrorDE(cellDE, mirrorVolName + "DE", 6 * cellCounter+1 );
      mirrorDE.setPlacement(mirrorPV);
      SkinSurface mirrorSkin(desc, mirrorDE, Form("mirror_optical_surface%d", cellCounter), mirrorSurf, mirrorVol); // FIXME: 3rd arg needs `imod`?
      mirrorSkin.isValid();


      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  COOLING PLATE  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
      auto coolingTrCell = RotationZYX(0, 0, angle_of_sensor ) *
                           Translation3D(0, center_of_sensor_x, sensor_z_pos-cooling_z_offset);

      Solid coolingSol = IntersectionSolid(cellS, coolingSol_tube, coolingTrCell);
      std::string coolingName = create_part_name_ff("cooling");
      /// TODO: change material
      Volume coolingVol( coolingName , coolingSol, mirrorMat );
      coolingVol.setVisAttributes( desc.visAttributes("arc_cooling_vis") );
      cellV.placeVolume(coolingVol);

      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  AEROGEL PLATE  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
      // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
      auto aerogelTrCell = RotationZYX(0, 0, angle_of_sensor ) *
                           Translation3D(0, center_of_sensor_x, sensor_z_pos+aerogel_z_offset);

      Solid aerogelSol = IntersectionSolid(cellS, aerogelSol_tube, aerogelTrCell);
      std::string aerogelName = create_part_name_ff("aerogel");
      Volume aerogelVol( aerogelName , aerogelSol, aerogelMat );
      aerogelVol.setVisAttributes( desc.visAttributes("arc_aerogel_vis") );
      cellV.placeVolume(aerogelVol);

      auto sensorTr = RotationZYX(alpha - 90 * deg, 0 , angle_of_sensor )*
                           Translation3D(0, center_of_sensor_x, sensor_z_pos );


      PlacedVolume sensorPV = cellV.placeVolume(sensorVol, sensorTr);
//       sensorPV.addPhysVolID("cellnumber", 6 * cellCounter+2);
#ifdef DUMP_SENSOR_POSITIONS
      ofile_sensor_pos  << 6 * cellCounter+2 << '\t'
                         << ncell.RID << '\t'
                         << ncell.isReflected << '\t'
                         << ncell.x << '\t'
                         << ncell.y << '\t'
                         << phin << '\n';
#endif
      DetElement sensorDE(cellDE, create_part_name_ff("sensor") + "DE", 6 * cellCounter+2 );
      sensorDE.setType("tracker");
      sensorDE.setPlacement(sensorPV);

      PlacedVolume cellPV = endcap_cells_gas_envelope.placeVolume(cellV, RotationZ(phistep * phin) * Translation3D(ncell.x, ncell.y, 0));
      cellPV.addPhysVolID("cellnumber", createPhysVolID() );//6*cellCounter + 0);
      cellDE.setPlacement( cellPV );

      if ( ncell.isReflected)
      {
        std::string cellRefName = create_part_name_ff("cell_ref");
        Volume cellV_reflected(cellRefName, cellS, gasvolMat);
        cellV_reflected.setVisAttributes( gasvolVis );
        DetElement cell_reflected_DE(det, cellRefName+"DE", 6 * cellCounter + 3);
        Transform3D mirrorTr_reflected(RotationZYX(0, 0, 0), Translation3D(-dx, dy, center_of_sphere_z));

        /// Define the actual mirror as intersection of the mother volume and the hollow sphere just defined
        Solid mirrorSol_reflected = IntersectionSolid(cellS, mirrorShapeFull, mirrorTr_reflected);
        Volume mirrorVol_reflected(mirrorVolName + "_ref1", mirrorSol_reflected, mirrorMat);
        mirrorVol_reflected.setVisAttributes(desc.visAttributes(Form("arc_mirror_vis%d", ncell.RID)));
        PlacedVolume mirror_ref_PV = cellV_reflected.placeVolume(mirrorVol_reflected);

        DetElement mirror_ref_DE(cell_reflected_DE, mirrorVolName + "_ref1" + "DE", 6 * cellCounter+4 );
        mirror_ref_DE.setPlacement(mirror_ref_PV);
        SkinSurface mirror_ref_Skin(desc, mirror_ref_DE, Form("mirror_ref_optical_surface%d", cellCounter), mirrorSurf, mirrorVol_reflected); // FIXME: 3rd arg needs `imod`?
        mirror_ref_Skin.isValid();

//         Transform3D sensorTr_reflected(RotationZYX(-alpha + 90 * deg, 0 /*90*deg-angle_of_sensor*/, angle_of_sensor),
//                                        Translation3D(0, center_of_sensor_x, sensor_z_origin_Martin));
        auto sensorTr_reflected = RotationZYX(-alpha + 90 * deg, 0 /*90*deg-angle_of_sensor*/, angle_of_sensor)*
                                       Translation3D(0, center_of_sensor_x, sensor_z_origin_Martin);
        PlacedVolume sensor_ref_PV = cellV_reflected.placeVolume(sensorVol, sensorTr_reflected);
//         sensor_ref_PV.addPhysVolID("cellnumber", 6 * cellCounter+5);
        DetElement sensor_ref_DE(cell_reflected_DE, create_part_name_ff("sensor") + "_ref_DE", 6 * cellCounter+5 );
        sensor_ref_DE.setPlacement(sensor_ref_PV);

#ifdef DUMP_SENSOR_POSITIONS
        ofile_sensor_pos  << 6 * cellCounter+5 << '\t'
                    << ncell.RID << '\t'
                    << ncell.isReflected << '\t'
                    << ncell.x << '\t'
                    << ncell.y << '\t'
                    << phin << '\n';
#endif

        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  COOLING PLATE  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  AEROGEL PLATE  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
          cellV_reflected.placeVolume(coolingVol);
          cellV_reflected.placeVolume(aerogelVol);


        PlacedVolume cell_ref_PV = endcap_cells_gas_envelope.placeVolume(cellV_reflected, RotationZ(phistep * phin) * Translation3D(-ncell.x, ncell.y, 0));
        cell_ref_PV.addPhysVolID("cellnumber", createPhysVolID() );//6*cellCounter + 3);
        cell_reflected_DE.setPlacement( cell_ref_PV );
      }
    } //-- end loop for sector
  }   //-- end loop for endcap

  endcap_cells_vessel_envelope.placeVolume(endcap_cells_gas_envelope);



  Assembly endcaps_assemblyV("endcaps_assemblyV");

  Transform3D endcapZPos_Tr(RotationZYX(0,0,0), Translation3D(0, 0, zpos_endcap));
  PlacedVolume endcapZPos_PV = endcaps_assemblyV.placeVolume(endcap_cells_vessel_envelope, endcapZPos_Tr);
  endcapZPos_PV.addPhysVolID("barrel", 1);

  DetElement endcapZPos_DE(det, "endcapZPos_DE", 0 );
  endcapZPos_DE.setPlacement(endcapZPos_PV);

/*
  Transform3D envelope_zreflected_Tr(RotationZYX( 0 ,0,180*deg), Translation3D(0, 0, -zpos_endcap));
  PlacedVolume endcapZNeg_PV = endcaps_assemblyV.placeVolume(endcap_cells_vessel_envelope, envelope_zreflected_Tr);
  endcapZNeg_PV.addPhysVolID("barrel", 2);

  DetElement endcapZNeg_DE(det, "endcapZNeg_DE", 2 );
  endcapZNeg_DE.setPlacement(endcapZNeg_PV);*/


  PlacedVolume endcaps_PV = motherVol.placeVolume(endcaps_assemblyV);
  endcaps_PV.addPhysVolID("system", detID);
  det.setPlacement(endcaps_PV);

  return det;
}
DECLARE_DETELEMENT(ARCENDCAP_T, create_endcap_cell)

/// Deprecated
static Ref_t create_endcap_cell_volumes(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
//   xml::Component dims = detElem.dimensions();
//   OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement det(detName, detID);
  sens.setType("tracker");

  // read Martin file and store parameters by name in the map
//   fill_cell_parameters_m();

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
//   double cooling_thickness = 1 * cm;

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
//   double thickness_sphere(10 * mm);
//   double mirror_z_origin_Martin = vessel_length / 2. - vessel_wall_thickness - 37 * cm;

  // Light sensor parameters
//   double sensor_sidex = 8 * cm;
//   double sensor_sidey = 8 * cm;
//   double sensor_thickness = 0.2 * cm;
//   double sensor_z_origin_Martin = -vessel_length / 2. + vessel_wall_thickness + 0.5 * cooling_thickness;

  // Build cylinder for gas.
  Tube gasvolSolid(vessel_inner_r + vessel_wall_thickness,
                   vessel_outer_r - vessel_wall_thickness,
                   vessel_length);

  // Use regular polyhedra for endcaps cells
  PolyhedraRegular cellS(6, 0, 0., hexagon_apothem, vessel_length);

  // Build the mirror for ncell=1..21

  for (auto ncell : mycell_v)
  {
    for (int phin = 0; phin < phinmax; phin++)
    {

      std::string volname = detName + "_cell" + std::to_string(ncell.RID);
      volname += "_phin" + std::to_string(phin);
      Volume cellV(volname, cellS, desc.material("Aluminum"));
      cellV.setVisAttributes(desc.visAttributes(Form("arc_mirror_vis%d", ncell.RID)));
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
