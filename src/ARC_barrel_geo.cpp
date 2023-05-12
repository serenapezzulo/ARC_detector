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


#include "ARC_par_reader.hpp"

/**
 * create barrel as sum of single cells.
 * next step is to place mirrors+sensors in a cylindral shape gas volume
 */
static Ref_t create_barrel_cell(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
//   xml::Component dims = detElem.dimensions();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement det(detName, detID);
  sens.setType("tracker");

  auto vesselMat = desc.material(detElem.attr<std::string>(_Unicode(vessel_material)));
  auto gasvolMat = desc.material(detElem.attr<std::string>(_Unicode(gas_material)));
  auto vesselVis = desc.visAttributes(detElem.attr<std::string>(_Unicode(vessel_vis)));
  auto gasvolVis = desc.visAttributes(detElem.attr<std::string>(_Unicode(gas_vis)));

  // read Martin file and store parameters by name in the map
  fill_cell_parameters_m();

  // mother volume corresponds to the world
  Volume motherVol = desc.pickMotherVolume(det);


  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //          VESSEL PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  double vessel_outer_r = 210 * cm;
  double vessel_inner_r = 190 * cm;
  double vessel_length = 440 * cm;
  double vessel_wall_thickness = 1.0 * cm;
  if (vessel_outer_r <= vessel_inner_r)
    throw std::runtime_error("Ilegal parameters: vessel_outer_r <= vessel_inner_r");
  // // //-------------------------------------------------------------// // //


  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //         AEROGEL PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  double aerogel_radial_thickness = 1.0 * cm;
  auto aerogelMat = desc.material("Aerogel_PFRICH");

  // // //-------------------------------------------------------------// // //

  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //         COOLING PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  double cooling_radial_thickness = 1.0 * cm;
  // // //-------------------------------------------------------------// // //


  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //           CELL PARAMETERS           // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  /// Distance in x-direction
  double hexagon_side_length = 14.815 * cm;
  double hex_apothem_length = hexagon_side_length*cos( M_PI / 6. ); //
  double zstep = 2 * hex_apothem_length;
  /// Distance in phi angle between cells
  /// since cells are regular hexagons, this distance matches
  /// the angle that one cell covers
  double phistep = 13.333 * deg;
  /// number of repetition of unique cells around the barrel
  int phinmax = 27; // 27;
  // // //-------------------------------------------------------------// // //


  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // // // //          MIRROR PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  auto mirrorElem = detElem.child(_Unicode(mirror)).child(_Unicode(module));
  auto mirrorSurf = surfMgr.opticalSurface(mirrorElem.attr<std::string>(_Unicode(surface)));
  auto mirrorMat = desc.material(mirrorElem.attr<std::string>(_Unicode(material)));
  double mirrorThickness = mirrorElem.attr<double>(_Unicode(thickness));
  // if this z shrink is not applied, the upper tip of the mirrors are cut
  // TODO: crosscheck with Martin distances between mirrors and sensors
  double mirror_z_safe_shrink = 6*mm;
  double mirror_z_origin_Martin = vessel_outer_r - vessel_wall_thickness - 37 * cm  - mirrorThickness - mirror_z_safe_shrink;
  // // //-------------------------------------------------------------// // //



  // // // // // // // // // // // // // // // // // // // // // // // // // //
  // // // // // //          LIGHT SENSOR PARAMETERS          // // // // // //
  // // // // // // // // // // // // // // // // // // // // // // // // // //
  //default values
  double sensor_sidex     = 8 * cm;
  double sensor_sidey     = 8 * cm;
  double sensor_thickness = 0.2 * cm;
  // empirical distance to keep the sensor inside the cell volume
  double sensor_z_safe_distance = 1.25*mm;
  double sensor_z_origin_Martin = vessel_inner_r + vessel_wall_thickness + cooling_radial_thickness + sensor_z_safe_distance;
  auto sensorMat = desc.material("SiliconOptical");
  auto sensorVis = desc.visAttributes("no_vis");
  // auto sensorSurf = surfMgr.opticalSurface(sensorElem.attr<std::string>(_Unicode(surface)));

  // Read from xml the parameters for the sensor module
  {
    auto sensorElem  = detElem.child(_Unicode(sensors)).child(_Unicode(module));
    sensor_sidex     = sensorElem.attr<double>(_Unicode(sensor_side_Phi));
    sensor_sidey     = sensorElem.attr<double>(_Unicode(sensor_side_Z));
    sensor_thickness = sensorElem.attr<double>(_Unicode(thickness));
    sensorMat        = desc.material(sensorElem.attr<std::string>(_Unicode(material)));
    sensorVis        = desc.visAttributes(sensorElem.attr<std::string>(_Unicode(vis)));
  }
  // // //-------------------------------------------------------------// // //



  // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //
  // // //+++++++++++  BUILD VESSEL, CELL AND SENSOR VOLUMES ++++++++++// // //
  // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //

  // Build cylinder for gas.
  Tube gasenvelopeS(  vessel_inner_r + vessel_wall_thickness,
                   vessel_outer_r - vessel_wall_thickness,
                   vessel_length/2.);
  Volume barrel_cells_gas_envelope (detName+"_gasEnvelope", gasenvelopeS, gasvolMat );
  barrel_cells_gas_envelope.setVisAttributes( desc.visAttributes("envelope_vis") );
  Tube vesselEnvelopeSolid(  vessel_inner_r ,
                             vessel_outer_r ,
                             vessel_length/2. + vessel_wall_thickness);
  Volume barrel_cells_vessel_envelope (detName+"_vesselEnvelope", vesselEnvelopeSolid, vesselMat );
  barrel_cells_vessel_envelope.setVisAttributes( vesselVis );
  // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //

  // Use truncated pyramid for barrel cells
  double angle_hex = 2*asin( hex_apothem_length / vessel_outer_r );
  std::vector<double> zplanes = { vessel_inner_r + vessel_wall_thickness,
                                 (vessel_outer_r - vessel_wall_thickness)*cos(angle_hex) };
  std::vector<double> rs = { hex_apothem_length -1*mm, hex_apothem_length-1*mm };
  /// Hexagonal truncated pyramid
  Polyhedra cell_shape("cellShape", 6, 30 * deg, 360 * deg, zplanes, rs);
  /// rotation of 90deg around Y axis, to align Z axis of pyramid with X axis of cylinder
  Transform3D pyramidTr(RotationZYX(0, -90. * deg, 0. * deg), Translation3D(0, 0, 0));

  // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //


  // // //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++// // //

  // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
  // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
  // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
  // Build the mirror for ncell=1..18
  std::vector<int> ncell_vector = {-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15, -16, /*-17, -18,*/
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 /*, 17, 18 */};

    int cellCounter(0);
    // WARNING for developping purposes
//         ncell_vector = {1};

//         phinmax = 1;

    // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
    // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
    // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
    for (auto ncell_original : ncell_vector)
    {
      // Activate flag to reflect parameters later
      int ncell = ncell_original;
      bool reflect_parameters = false;
      if (0 > ncell)
      {
        ncell *= -1;
        reflect_parameters = true;
      }
      for (int phin = 0; phin < phinmax; ++phin)
      {

        // WARNING for developping purposes
        // The following line skips even number cells
        // if (!(ncell % 2))
        //   continue;

        // The following line skips odd number cells
        // if ((ncell % 2))
        // continue;
        // there is no cell number 0, and cell number 1 do not need to be reflected
        if (0 == ncell || -1 == ncell)
          continue;


        auto create_part_name_ff = [ncell,detName,phin,reflect_parameters](std::string  partName){
          std::string fullName = detName + "_" + partName;
          fullName += std::to_string(ncell);
          fullName += "_ref" + std::to_string(reflect_parameters);
          fullName += "_phi" +  std::to_string(phin);
          dd4hep::printout(dd4hep::DEBUG,"ARCBARREL_T", "+++ New name:%s",fullName.c_str());
          return fullName;
        };






        std::string cellName = create_part_name_ff( "cell");

        /// Volume that contains gas and other stuff
        Volume cellVol(cellName, cell_shape, gasvolMat);
        cellVol.setVisAttributes( gasvolVis );
        /// Detector element that will contain cellVol later
        /// there are 3 elements with ID:
        /// the cell, ID= 3 * cellCounter
        /// its mirror. ID = 3 * cellCounter +1
        /// and its sensor, ID = 3 * cellCounter +2
        DetElement cellDE(det, cellName+"DE", 3 * cellCounter);



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
          if (radius_of_sphere <= mirrorThickness)
            throw std::runtime_error("Ilegal parameters: radius_of_sphere <= mirrorThickness");

          if (-999. == center_of_sensor_x)
            throw std::runtime_error("Ilegal parameters: center_of_sensor_x not provided");
          if (-999. == angle_of_sensor)
            throw std::runtime_error("Ilegal parameters: angle_of_sensor not provided");
        }

        // reflect parameters for cells with z<0,
        // ie, ncell < 0
        if (reflect_parameters)
        {
          center_of_sphere_x *= -1.0;
          center_of_sensor_x *= -1.0;
          angle_of_sensor *= -1.0;
        }

        // create the semi-sphere that will result in the mirror
        Sphere mirrorShapeFull(radius_of_sphere - mirrorThickness,
                               radius_of_sphere,
                               0.,
                               3.14 / 2);
        /// 3D transformation of mirrorVolFull in order to place it inside the gas volume
        Transform3D mirrorTr(RotationZYX(0, 0, 0), Translation3D(center_of_sphere_x, 0, center_of_sphere_z - mirror_z_safe_shrink));

        // TODO: cell 18 corresponds to half a pyramid, currently is full pyramid
        Solid mirrorSol = IntersectionSolid(cell_shape, mirrorShapeFull, mirrorTr);

        std::string mirrorName = create_part_name_ff("mirror"); // detName + "_mirror" + std::to_string(ncell) + "z" + std::to_string(reflect_parameters)

        Volume mirrorVol( mirrorName , mirrorSol, mirrorMat );
        mirrorVol.setVisAttributes(desc.visAttributes(Form("mirror_vis%d", ncell)));
        PlacedVolume mirrorPV = cellVol.placeVolume(mirrorVol);

        DetElement mirrorDE(cellDE, mirrorName + "DE", 3 * cellCounter+1 );
        mirrorDE.setPlacement(mirrorPV);
        SkinSurface mirrorSkin(desc, mirrorDE, Form("mirror_optical_surface%d", cellCounter), mirrorSurf, mirrorVol); // FIXME: 3rd arg needs `imod`?
        mirrorSkin.isValid();

        // Place detector in cell
        // Build sensor shape
        // TODO dummy while developping
//         angle_of_sensor = 45*deg;
        Box sensorSol(sensor_sidex / 2, sensor_sidey / 2, sensor_thickness / 2);
        std::string sensorName = create_part_name_ff("sensor");
        Volume sensorVol( sensorName , sensorSol, sensorMat );
        sensorVol.setVisAttributes( sensorVis );
        sensorVol.setSensitiveDetector(sens);

        Transform3D sensorTr(RotationZYX(0, 90 * deg - angle_of_sensor, 0), Translation3D(-sensor_z_origin_Martin, 0, center_of_sensor_x));
        PlacedVolume sensorPV = cellVol.placeVolume(sensorVol, RotationZYX(0, 90. * deg, 0. * deg)*sensorTr);
        sensorPV.addPhysVolID("cellnumber", 3 * cellCounter+2);
        DetElement sensorDE(cellDE, sensorName + "DE", 3 * cellCounter+2 );
        sensorDE.setType("tracker");
        sensorDE.setPlacement(sensorPV);

        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  COOLING PLATE  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
        {
          double cooling_z_offset =   sensor_thickness  + cooling_radial_thickness;
          Tube coolingSol_tube(0, 1.5*hexagon_side_length, cooling_radial_thickness/2.);
          Transform3D coolingTr(RotationZYX(0, 90 * deg - angle_of_sensor, 0), Translation3D(-sensor_z_origin_Martin+cooling_z_offset, 0, center_of_sensor_x));
          auto coolingTrCell = RotationZYX(0, 90. * deg, 0. * deg)*coolingTr;
          Solid coolingSol = IntersectionSolid(cell_shape, coolingSol_tube, coolingTrCell);
          std::string coolingName = create_part_name_ff("cooling");
          /// TODO: change material
          Volume coolingVol( coolingName , coolingSol, mirrorMat );
          coolingVol.setVisAttributes( desc.visAttributes("cooling_vis") );
          cellVol.placeVolume(coolingVol );
        }
        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //

        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  AEROGEL PLATE  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
        double aerogel_z_offset =   sensor_thickness + aerogel_radial_thickness;
        Tube aerogelSol_tube(0, 1.5*hexagon_side_length, cooling_radial_thickness);
        Transform3D aerogelTr(RotationZYX(0, 90 * deg - angle_of_sensor, 0), Translation3D(-sensor_z_origin_Martin - aerogel_z_offset, 0, center_of_sensor_x));
        auto aerogelTrCell = RotationZYX(0, 90. * deg, 0. * deg)*aerogelTr;
        Solid aerogelSol = IntersectionSolid(cell_shape, aerogelSol_tube, aerogelTrCell);
        std::string aerogelName = create_part_name_ff("aerogel");
        /// TODO: change material
        Volume aerogelVol( aerogelName , aerogelSol, aerogelMat );
        aerogelVol.setVisAttributes( desc.visAttributes("aerogel_vis") );
        cellVol.placeVolume(aerogelVol );
        // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //


        // position of mirror in cylinder coordinate system
        double mirror_abs_pos_z = name_col * zstep - 0.5 * zstep * (2 == name_row);
        if (reflect_parameters)
          mirror_abs_pos_z *= -1.0;

        // row 2 is shifted half step size
        double phi_offset = 0 + 0.5 * phistep * (2 == name_row);


        auto cellTr = RotationZ(phistep * phin + phi_offset) * Translation3D(0, 0, mirror_abs_pos_z)*pyramidTr ;
        PlacedVolume cellPV = barrel_cells_gas_envelope.placeVolume(cellVol, cellTr);
        cellDE.setPlacement(cellPV);


        // increase counter
        cellCounter++;

      }
    }
    // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
    // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
    // // // ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> ~> // // //
    barrel_cells_vessel_envelope.placeVolume(barrel_cells_gas_envelope);
    PlacedVolume assemblyPV = motherVol.placeVolume(barrel_cells_vessel_envelope);
    assemblyPV.addPhysVolID("system", detID);
    det.setPlacement(assemblyPV);
    return det;
}

DECLARE_DETELEMENT(ARCBARREL_T, create_barrel_cell)

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
  //   double mirror_diameter_safe_shrink = 0*mm;
  //   double mirror_z_safe_shrink = 0*mm;
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
    double pyramid_height = vessel_outer_r - vessel_wall_thickness /*- mirror_z_safe_shrink*/;
    double pyramid_base_apothem = hexagon_apothem;
    double myfactor = cos( phistep /2 );
    std::vector<double> zplanes = {0 * cm, pyramid_height*myfactor };
    std::vector<double> rs = {0 * cm, pyramid_base_apothem*myfactor };
    /// Hexagonal pyramid
    Polyhedra shape("mypyramid", 6, 30 * deg, 360 * deg, zplanes, rs);
    /// rotation of 90deg around Y axis, to align Z axis of pyramid with X axis of cylinder
    Transform3D pyramidTr(RotationZYX(0, -90. * deg, 0. * deg), Translation3D(0, 0, 0));

    // Build sensor shape
    Box sensorSol(sensor_sidex / 2, sensor_sidey / 2, sensor_thickness / 2);
    Volume sensorVol(detName + "_sensor", sensorSol, desc.material("SiliconOptical"));
    sensorVol.setSensitiveDetector(sens);
//     SkinSurface sensorSkin(desc, det, "sensor_optical_surface", sensorSurf, sensorVol); // FIXME: 3rd arg needs `imod`?
//     sensorSkin.isValid();

    // Build the mirror for ncell=1..18
    std::vector<int> ncell_vector = {-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18,
      1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
      // std::vector<int> ncell_vector = { /*-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13, -14,-15,-16,-17,-18,*/
      //                                 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18
      //                                 };
//       ncell_vector = {17};

      /// Dummy counter to place elements inside the barrel
      int cellCounter(0);
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
        Transform3D mirrorTr(RotationZYX(0, 0, 0), Translation3D(center_of_sphere_x , 0, center_of_sphere_z /*- mirror_z_safe_shrink*/));

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
          sensorPV.addPhysVolID("cellnumber", 2*cellCounter);

          std::cout << sensorPV.volIDs().str() << std::endl;

          // // Make sensor sensitive + define optical properties
          DetElement sensorDE(det, Form("ARC_DEsensor%d", cellCounter), 2 * cellCounter );
          sensorDE.setPlacement(sensorPV);
          SkinSurface sensorSkin(desc, sensorDE, Form("sensor_optical_surface%d", cellCounter), sensorSurf, sensorVol); // FIXME: 3rd arg needs `imod`?
          sensorSkin.isValid();

          // create mirrors as separate detectors, so properties can be adjusted later!
          DetElement mirrorDE(det, Form("ARC_DEmirror%d", cellCounter), 2 * cellCounter +1);
          mirrorDE.setPlacement(mirrorPV);
          SkinSurface mirrorSkin(desc, mirrorDE, Form("mirror_optical_surface%d", cellCounter), mirrorSurf, mirrorVol); // FIXME: 3rd arg needs `imod`?
          mirrorSkin.isValid();

          // increase counter
          cellCounter++;
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

DECLARE_DETELEMENT(ARCBARREL_LEGACY_T, create_barrel)
