#include "ARCalg.h"
#include "ArcQuarticEquation.h"
#include "Reconstruction.h"
#include "test.h"
#include <TVector3.h>
#include <map>
#include <random>

// ROOT
#include "TGeoBoolNode.h"
#include <TGeoManager.h>
#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"



// STL
#include <iostream>
#include <sstream>

#include <TH1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TApplication.h>


std::map<double, double> getPDETable(const std::string& tableName, dd4hep::Detector* detector) {
    std::map<double, double> pdeTable;



    const TObjArray* c = detector->manager().GetListOfGDMLMatrices();
    TObjArrayIter arr(c);
    //std::cout << "GDML tables " << std::endl;
    for (TObject* i = arr.Next(); i; i = arr.Next()) {
      TGDMLMatrix* gdmlMat = (TGDMLMatrix*)i;
      if (gdmlMat->GetName() == tableName) {
       // std::cout << "Trovata la matrice: " << gdmlMat->GetName() << std::endl;
	int rows = gdmlMat->GetRows();
        int cols = gdmlMat->GetCols();
	
	

	for (size_t r = 0; r < rows; ++r) {
	  double energy = gdmlMat->Get(r, 0);      
	  double efficiency = gdmlMat->Get(r, 1);
          
	  pdeTable[energy] = efficiency;
	  //std::cout << "energy: " << energy << " eV, efficiency: " << efficiency << std::endl;
        }
	break;
	
      }
    }

        
    return pdeTable;
}

double interpolatePDE(double energy, const std::map<double, double>& pdeTable) {
       if (pdeTable.empty()) {
        std::cerr << "PDE table is empty." << std::endl;
        return 0.0; // or another appropriate value
    }

    // Print the PDE table for debugging
    //std::cout << "PDE Table:" << std::endl;
    for (const auto& pair : pdeTable) {
        //std::cout << "Energy: " << pair.first << ", PDE: " << pair.second << std::endl;
    }

    // Find two closest points
    auto upper = pdeTable.upper_bound(energy);
    if (upper == pdeTable.begin()) {
        // Energy is less than the smallest value in the table
        return upper->second;
    }
    if (upper == pdeTable.end()) {
        // Energy is greater than the largest value in the table
        return std::prev(upper)->second;
    }
    
    auto lower = std::prev(upper);
    double x0 = lower->first;
    double x1 = upper->first;
    double y0 = lower->second;
    double y1 = upper->second;

    // Linear interpolation
    return y0 + (energy - x0) * (y1 - y0) / (x1 - x0);
}


ROOT::Math::XYZPoint IntersectSphereWithTrack(
					      const ROOT::Math::XYZPoint& shootingPoint,//mm
					      const ROOT::Math::XYZVector& trace,
					      const ROOT::Math::XYZPoint& mirrorCenter,
					      double mirrorRadius)
{
  // Componenti della direzione del vettore traccia
  double vx = trace.x();
  double vy = trace.y();
  double vz = trace.z();
  
  // Punto iniziale della traccia (shooting point)
  double x0 = shootingPoint.x();
  double y0 = shootingPoint.y();
  double z0 = shootingPoint.z();

  std::cout<<"shooting x "<<shootingPoint.x()<<std::endl;
  std::cout<<"shooting y "<<shootingPoint.y()<<std::endl;
  std::cout<<"shooting z "<<shootingPoint.z()<<std::endl;
  
  
  // Coordinate del centro dello specchio
  double x_mirror = mirrorCenter.x()*10;
  double y_mirror = mirrorCenter.y()*10;
  double z_mirror = mirrorCenter.z()*10;

  std::cout<<"x_mirror "<<x_mirror<<std::endl;
  std::cout<<"y_mirror "<<y_mirror<<std::endl;
  std::cout<<"z_mirror "<<z_mirror<<std::endl;
    
  
  std::cout<<"mirrorRadius "<<mirrorRadius<<std::endl;
  mirrorRadius *=10;
  
  double A = vx * vx + vy * vy + vz * vz;
  double B = 2 * ((x0 - x_mirror) * vx + (y0 - y_mirror) * vy + (z0 - z_mirror) * vz);
  double C = (x0 - x_mirror) * (x0 - x_mirror) +
    (y0 - y_mirror) * (y0 - y_mirror) +
    (z0 - z_mirror) * (z0 - z_mirror) - mirrorRadius * mirrorRadius;
  
  
  double discriminant = B * B - 4 * A * C;
  
  
  if (discriminant < 0) {
    
    std::cerr << "La traccia non interseca lo specchio sferico." << std::endl;
    return ROOT::Math::XYZPoint(); 
  }

  
  double t1 = (-B + sqrt(discriminant)) / (2 * A);
  double t2 = (-B - sqrt(discriminant)) / (2 * A);
  
  
  ROOT::Math::XYZPoint intersection1(x0 + t1 * vx, y0 + t1 * vy, z0 + t1 * vz);
  ROOT::Math::XYZPoint intersection2(x0 + t2 * vx, y0 + t2 * vy, z0 + t2 * vz);
  
  
  ROOT::Math::XYZPoint intersection = (intersection1.z() > intersection2.z()) ? intersection1 : intersection2;
  std::cout<<"intersezione1 specchio "<<intersection1.x()<<" "<<intersection1.y()<<" "<<intersection1.z()<<std::endl;
  std::cout<<"intersezione2 specchio "<<intersection2.x()<<" "<<intersection2.y()<<" "<<intersection2.z()<<std::endl;

  //double piano=-87.0;
  double piano = -77.0;
  
  double t_plane = (piano - z0) / vz;
  ROOT::Math::XYZPoint intersection_plane(x0 + t_plane * vx, y0 + t_plane * vy, piano);
  std::cout<<"intersezione piano a -90 "<<intersection_plane<<std::endl;
  
  
  ROOT::Math::XYZPoint midpoint(
				(intersection_plane.x() + intersection.x()) / 2.0,
				(intersection_plane.y() + intersection.y()) / 2.0,
				(intersection_plane.z() + intersection.z()) / 2.0
				);
  
  
  return midpoint;
}


ROOT::Math::XYZPoint IntersectPlanesWithTrack(
                const ROOT::Math::XYZPoint& shootingPoint,
                const ROOT::Math::XYZVector& trace)
{

  // Componenti della direzione del vettore traccia
  double vx = trace.x();
  double vy = trace.y();
  double vz = trace.z();
  std::cout<<"trace "<<trace<<std::endl;
  
  // Punto iniziale della traccia (shooting point)
  double x0 = shootingPoint.x();
  double y0 = shootingPoint.y();
  double z0 = shootingPoint.z();

  double plane1 = -87.0;
  double plane2 = -77.0;

  ROOT::Math::XYZPoint intersection1(x0 + (plane1 - z0) * vx / vz, y0 + (plane1 - z0) * vy / vz, plane1);
  ROOT::Math::XYZPoint intersection2(x0 + (plane2 - z0) * vx / vz, y0 + (plane2 - z0) * vy / vz, plane2);
  std::cout<<"intersezione1 piano "<<intersection1.x()<<" "<<intersection1.y()<<" "<<intersection1.z()<<std::endl;
  std::cout<<"intersezione2 piano "<<intersection2.x()<<" "<<intersection2.y()<<" "<<intersection2.z()<<std::endl;


  ROOT::Math::XYZPoint midpoint(
      (intersection1.x() + intersection2.x()) / 2.0,
      (intersection1.y() + intersection2.y()) / 2.0,
      (intersection1.z() + intersection2.z()) / 2.0
  );

  std::cout<<"midpoint "<<midpoint<<std::endl;
  return midpoint;


}

    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       ARCalg constructor       ////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    // -- KeyValues("name of the variable that holds the name of the collection exposed in the python steering file", {"default name for the collection"}),
    ARCalg::ARCalg(const std::string& name, ISvcLocator* svcLoc)
        : MultiTransformer(name, svcLoc,
                            {
                                KeyValues("ARC_simhits", {""}),
                                KeyValues("MCParticles", {""}),
                                KeyValues("HeaderName", {"EventHeader"}),
                            },
                            {
                                KeyValues("ARC_PID", {"ARC_PID"})
                            }
                        )
    {
    m_geoSvc = serviceLocator()->service(m_geoSvcName);
    }


    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       initialize       ////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    StatusCode ARCalg::initialize() {


        //-----------------
        // Retrieve the subdetector
      std::string ARC_name(m_ARC_name.value());
        if ( 0 == m_geoSvc->getDetector()->detectors().count(ARC_name) )
        {
            ThrowException( "Detector <<" + ARC_name + ">> does not exist." );
        }
	
        dd4hep::DetElement ARC_det = m_geoSvc->getDetector()->detectors().at(ARC_name);

        ///////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////  retrieve mirror properties  //////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////
        auto cellDE = ARC_det; // for case of cubic single cell of ARC
        dd4hep::DetElement mirrorDE = ARC_det.child("ARC_mirror");
        auto mirrorSolid = mirrorDE.solid();
        auto m = (TGeoCompositeShape*)mirrorSolid.access();

        /// Retrieve the center of the sphere
        //-- this is just the matrix for building the intersection, translation
        auto matrix_boolean = m->GetBoolNode()->GetRightMatrix();
        //-- matrix for placing the cell
        auto matrix_cellvol = cellDE.nominal().worldTransformation();
        double local_coord [3] = {0.,0.,0.};
        double global_coord[3] = {0.,0.,0.};
        (matrix_cellvol**matrix_boolean).LocalToMaster(local_coord, global_coord);
        /// global_coord contains now the XYZ global coordinates of the center of the mirror
        mirrorCenter = dd4hep::Position( global_coord[0], global_coord[1], global_coord[2]);
	      mirrorCenterXYZ = ROOT::Math::XYZPoint(  global_coord[0], global_coord[1], global_coord[2]);
        /// Retrieve the inner radius of the spherical mirror
        TGeoSphere * spherical_mirror_shape = (TGeoSphere*)m->GetBoolNode()->GetRightShape();
        mirrorRadius = spherical_mirror_shape->GetRmin();

        ///////////////////////////////////////////////////////////////////////////////////

        //-----------------
        // Retrieve the readout associated with the detector element (subdetector)
        dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector(ARC_name);
        if(not dch_sd.isValid() )
            ThrowException("No valid Sensitive Detector was found for detector <<" + ARC_name + ">>.");
	
        dd4hep::Readout dch_readout = dch_sd.readout();
        // set the cellID decoder
        m_decoder = dch_readout.idSpec().decoder();
        //-----------------

	/*
	auto description = m_geoSvc->getDetector();
	const TObjArray* c = description->manager().GetListOfGDMLMatrices();
	TObjArrayIter arr(c);
	std::cout<<"GDML tables "<<std::endl;
	//printout(INFO,"Dump_GDMLTables","+++ Dumping known GDML tables from TGeoManager.");
	for(TObject* i = arr.Next(); i; i=arr.Next())   {
	  TGDMLMatrix* gdmlMat = (TGDMLMatrix*)i;
	  //num_elements += (gdmlMat->GetRows()*gdmlMat->GetCols());
	  //++num_tables;
	  gdmlMat->Print();
	}
	*/
        std::stringstream ss;
        PrintConfiguration(ss);
        info() << ss.str().c_str() <<endmsg;

        if(m_create_debug_histos.value()){
            hThetaRecoTrue = new TH1D("hThetaRecoTrue", "Reconstructed Theta Truth", 500, 0, 0.3);
            hThetaRecoTrue->SetDirectory(0);
            hThetaReco = new TH1D("hThetaReco", "Reconstructed Theta", 500, 0, 0.3);
            hThetaReco->SetDirectory(0);
            hThetaRecoEm = new TH1D("hThetaRecoEm", "Reconstructed Theta", 500, 0, 0.3);
            hThetaRecoEm->SetDirectory(0);
            hThetaErrorEm = new TH1D("hThetaErrorEm", "Reconstructed Theta Error", 500, -0.015, 0.015);
            hThetaErrorEm->SetDirectory(0);
            hPixelID = new TH2F("hPixelID", "Pixel ID", 100,0,100,100,0,100 );
            hPixelID->SetDirectory(0);
            hTrueHit = new TH2F("hTrueHit", "Pixel Position", 160, -80, 80, 160, -80, 80);
            hTrueHit->SetDirectory(0);
            hPixelHit = new TH2F("hPixelHit", "Pixel Hit (center)", 160, -80, 80, 160, -80, 80);
            hPixelHit->SetDirectory(0);
            hPixelError = new TH1D("hPixelError", "Pixel Error", 500, -0.05, 0.05);
            hPixelError->SetDirectory(0);
            hThetaRecoPixel = new TH1D("hThetaRecoPixel", "Reconstructed Theta", 100, 0, 0.07);
            hThetaRecoPixel->SetDirectory(0);
            hPhotonYield = new TH1D("hPhotonYield", "Photon Yield", 100, 0, 150);
            hPhotonYield->SetDirectory(0);
	          hEmissionZTrue = new TH1D("hEmissionZTrue", "hEmissionZTrue", 100, -100, 100);
            hEmissionZTrue->SetDirectory(0);
	          hEmissionXTrue = new TH1D("hEmissionXTrue", "hEmissionXTrue", 1000, -200, 200);
            hEmissionXTrue->SetDirectory(0);
	          hEmissionYTrue = new TH1D("hEmissionYTrue", "hEmissionYTrue", 1000, -200, 200);
            hEmissionYTrue->SetDirectory(0);
            hPhotonYieldGas = new TH1D("hPhotonYieldGas", "Photon Yield Gas", 60, 0, 60);
            hPhotonYieldGas->SetDirectory(0);
            hPhotonYieldAerogel = new TH1D("hPhotonYieldAerogel", "Photon Yield Aerogel", 60, 0, 60);
            hPhotonYieldAerogel->SetDirectory(0);
            hThetaGas = new TH1D("hThetaGas", "Reconstructed em Theta Gas", 500, 0, 0.3);
            hThetaGas->SetDirectory(0);
            hThetaAerogel = new TH1D("hThetaAerogel", "Reconstructed Theta Aerogel", 500, 0, 0.3);
            hThetaAerogel->SetDirectory(0);
            hThetaGen = new TH1D("hThetaGen", "Generated Theta", 500, 0, 0.3);
            hThetaGen->SetDirectory(0);
            hRecoGen = new TH1D("hRecoGen", "Reco-Gen", 500, -0.2, 0.2);
            hRecoGen->SetDirectory(0);	   
            hEMerrorGas = new TH1D("hEMerrorGas", "EM error Gas", 500, -0.2, 0.2);
            hEMerrorGas->SetDirectory(0); 
            hEMerrorAerogel = new TH1D("hEMerrorAerogel", "EM error Aerogel", 500, -0.2, 0.2);
            hEMerrorAerogel->SetDirectory(0);
        }

      
        return StatusCode::SUCCESS;
    }

    

    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       operator()       ////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    std::tuple<colltype_out>
    ARCalg::operator()(const colltype_in& input_sim_hits, const colltype_mc& input_mc_particles,
                            const edm4hep::EventHeaderCollection&  /*headers*/) const {


    debug() << "Input Sim Hit collection size: " << input_sim_hits.size() << endmsg;

    // Create the collections we are going to return
    colltype_out output_digi_hits;

    int photonCounter = 0;
    int gamma_gas = 0;
    int gamma_aerogel = 0;

    //loop over hit collection
    for (const auto& input_sim_hit : input_sim_hits) {


    const edm4hep::MCParticle& photon = input_sim_hit.getMCParticle();
    //std::cout<<"PDG Photon "<<photon.getPDG()<<std::endl;
    double photonEnergy = photon.getEnergy();
    //double photonMomx = photon.getMomentum().x;
    TVector3 photonMom(photon.getMomentum().x, photon.getMomentum().y, photon.getMomentum().z);
    std::cout << "Momentum Photon" << photonMom.x() << " " << photonMom.y() << " " << photonMom.z() << std::endl;
    std::cout << "Energy Photon " << photonEnergy << std::endl;


    
    std::string qeTableName = "ARC_SiPM_QuantumEfficiency";
    std::string ARC_nome(m_ARC_name.value());

    
    
    dd4hep::Detector* detector = m_geoSvc->getDetector();
      
    
          dd4hep::DetElement element = m_geoSvc->getDetector()->detectors().at(ARC_nome);
    
    
    std::map<double, double> qeTable = getPDETable(qeTableName, detector);

    double pde = interpolatePDE(photonEnergy, qeTable);
    std::cout << "PDE for photon energy " << photonEnergy << " eV: " << pde << std::endl;

    //Random number
      std::random_device rd;  
      std::mt19937 gen(rd()); // Seed 
      std::uniform_real_distribution<> dis(0.0, 1.0); //Range

      double randomNumber = dis(gen);
      std::cout << "Generated random number: " << randomNumber << std::endl;

      // if random number < PDE, keep the photon
      if (randomNumber < pde) {
        
          photonCounter++;
          std::cout << "Photon kept. Current photon count: " << photonCounter << std::endl;
      
          dd4hep::DDSegmentation::CellID id = input_sim_hit.getCellID();
          //pixels ID -> integer numbers
          double pixel_x = m_decoder->get(id, "x");
          double pixel_y = m_decoder->get(id, "y");
          std::cout << "pixel_x: " << pixel_x << " pixel_y: " << pixel_y << std::endl;
          
          //auto pixelcenter_x = -3.96 + pixel_x*0.08;
          //auto pixelcenter_y = -3.96 + pixel_y*0.08;
          //auto pixelcenter_x = -3.975 + pixel_x*0.05;
          //auto pixelcenter_y = -3.975 + pixel_y*0.05;

          auto truehit_x = input_sim_hit.getPosition().x;
          auto truehit_y = input_sim_hit.getPosition().y;

          //  pixel indices (pixel 0.5mm)
          auto pixelnewx = floor(truehit_x / 0.5);  
          auto pixelnewy = floor(truehit_y / 0.5);  

        //  pixel center
          auto pixelcenter_x = (pixelnewx + 0.5) * 0.5; 
          auto pixelcenter_y = (pixelnewy + 0.5) * 0.5; 
          std::cout << "pixelcenter_x: " << pixelcenter_x << " pixelcenter_y: " << pixelcenter_y << std::endl;

         
          
          
          



          //auto mother = photon.getParents();
          int pdgid;
          for (auto parents : photon.getParents())
            {
              pdgid = parents.getPDG();
              //std::cout<<"PDG Mother "<<pdgid<<std::endl;
            
            }
          //if (TMath::Abs(pdgid)!=321) continue;

          

          
                // Get the emission point (vertex) of the photon
          auto emissionPoint = photon.getVertex();


          
          //std::cout<<"PDG Photon Mother "<<input_mc_particles.PDG<<std::endl;
          
          ROOT::Math::XYZPoint EmissionPointTrue(emissionPoint.x, emissionPoint.y, emissionPoint.z);
    //std::cout<<"EmissionPointTrue.z() "<<EmissionPointTrue.z()<<std::endl;
    //if(EmissionPointTrue.z()>(-0.0) || EmissionPointTrue.z()<(-20.0)) continue;
    
    
          // Convert edm4hep::Vector3d to ROOT::Math::XYZVector
          ROOT::Math::XYZPoint DetectionPointTrue(input_sim_hit.getPosition().x,
                                    input_sim_hit.getPosition().y,
                                    input_sim_hit.getPosition().z);

          ROOT::Math::XYZVector Trace(input_mc_particles.momentum().at(0).x,
                              input_mc_particles.momentum().at(0).y,
                              input_mc_particles.momentum().at(0).z);

         


          ROOT::Math::XYZPoint DetectionPointPixel = ROOT::Math::XYZPoint(pixelcenter_x, pixelcenter_y, input_sim_hit.getPosition().z);
          

          double TrueAngle = reconstruct(mirrorRadius, EmissionPointTrue, DetectionPointTrue, Trace);
          //std::cout<<"TrueAngle "<<TrueAngle<<std::endl;
          
        //Calculate the approx emission point and the emission point error

          // Define a dummy TGeoBBox (surface) that cuts the cell in half in the z direction
          double dx = 135.0;  
          double dy = 135.0;  
          double dz = 1e-7;  
          TGeoBBox* dummyBox = new TGeoBBox("dummyBox", dx, dy, dz);
          
          //get the shooting point (vertex) of the parent particle
          const edm4hep::MCParticle& mcParticle = input_sim_hit.getMCParticle();
          const edm4hep::MCParticle& parentParticle = mcParticle.getParents(0);
          auto shootingPoint = parentParticle.getVertex();
          //std::cout<<"shooting point "<<shootingPoint.x<<" "<<shootingPoint.y<<" "<<shootingPoint.z<<std::endl;
          ROOT::Math::XYZPoint ShootingPointTrue(shootingPoint.x, shootingPoint.y, shootingPoint.z);

    
    
          
                // Calculate the intersection between the track and the dummy box
          ROOT::Math::XYZPoint EmissionPointGas = IntersectSphereWithTrack(ShootingPointTrue, Trace, mirrorCenterXYZ, mirrorRadius);
          ROOT::Math::XYZPoint EmissionPointAerogel = IntersectPlanesWithTrack(ShootingPointTrue, Trace);
          
          //ROOT::Math::XYZPoint EmissionPointGas = CalculateIntersectionPoint(dummyBox, ShootingPointTrue, Trace);
          
          
          //std::cout << "Emission Point Approximation: " << EmissionPointGas.x() << " " << EmissionPointGas.y() << " " << EmissionPointGas.z() << std::endl;

      


          double Angle_em_gas = reconstruct(mirrorRadius, EmissionPointGas, DetectionPointTrue, Trace);
          double Angle_em_aerogel = reconstruct(mirrorRadius, EmissionPointAerogel, DetectionPointTrue, Trace);
          //std::cout<<"Angle_em_gas "<<Angle_em_gas<<std::endl;
          //std::cout<<"Angle_em_aerogel "<<Angle_em_aerogel<<std::endl;

          //for now the angles are hard coded, I will copy paste it from the reconstruction.cpp
          std::array<double, 2> delta_em = {
          TMath::Abs(Angle_em_gas - 0.051),
          TMath::Abs(Angle_em_aerogel - 0.221)};


          // Find the minimum difference and its index
          auto min_it = std::min_element(delta_em.begin(), delta_em.end());
          int min_index = std::distance(delta_em.begin(), min_it);
          std::cout<<"diff[i]"<<delta_em[0]<< ", "<< delta_em[1] <<std::endl;
          std::cout<<"min index "<<min_index<<std::endl;
          
          ROOT::Math::XYZPoint EmissionPoint;
          double Angle_em;
          double ThetaGas_em;
          double ThetaAerogel_em;
          if (min_index == 0) {
            Angle_em = Angle_em_gas;
            EmissionPoint = EmissionPointGas;
            gamma_gas ++;
            ThetaGas_em = reconstruct(mirrorRadius, EmissionPointGas, DetectionPointTrue, Trace);
            hThetaGas->Fill(ThetaGas_em);
            hEMerrorGas->Fill(ThetaGas_em - TrueAngle);

          } else {
            Angle_em = Angle_em_aerogel;
            EmissionPoint = EmissionPointAerogel;
            gamma_aerogel ++;
            ThetaAerogel_em = reconstruct(mirrorRadius, EmissionPointAerogel, DetectionPointTrue, Trace);
            hThetaAerogel->Fill(ThetaAerogel_em);
            hEMerrorAerogel->Fill(ThetaAerogel_em - TrueAngle);
          }


          double GenAngle;
          TVector3 v1(Trace.x(), Trace.y(), Trace.z());
          GenAngle = photonMom.Angle(v1);
          //std::cout<<"GenAngle "<<GenAngle<<std::endl;
          hThetaGen->Fill(GenAngle);
         
          
          double Angle_pixel = reconstruct(mirrorRadius, EmissionPointTrue, DetectionPointPixel, Trace);
          double AngleReco = reconstruct(mirrorRadius, EmissionPoint, DetectionPointPixel, Trace);

          if(m_create_debug_histos.value()){
            hThetaRecoTrue->Fill(TrueAngle);
            hRecoGen->Fill(TrueAngle - GenAngle);
            hThetaReco->Fill(AngleReco);
            hThetaRecoEm->Fill(Angle_em);
            hThetaErrorEm->Fill(Angle_em - TrueAngle);
            hPixelID->Fill(pixel_x, pixel_y);
            hTrueHit->Fill(truehit_x, truehit_y);
            hPixelHit->Fill(pixelcenter_x, pixelcenter_y);
            hThetaRecoPixel->Fill(Angle_pixel);
            hPixelError->Fill(Angle_pixel - TrueAngle);
            hEmissionZTrue->Fill(EmissionPointTrue.z());
            hEmissionXTrue->Fill(EmissionPointTrue.x());
            hEmissionYTrue->Fill(EmissionPointTrue.y());
            
            }

      
      } else {
              // Discard the photon
              std::cout << "Photon discarded." << std::endl;
              // Continue with the next photon
      }


        
      
        
          
          // emplace output objects
          // output_digi_hits.create(....);

      }// end loop over hit collection

      if (m_create_debug_histos.value()) {
          hPhotonYield->Fill(photonCounter);
          hPhotonYieldGas->Fill(gamma_gas);
          hPhotonYieldAerogel->Fill(gamma_aerogel);
      }

          

          /////////////////////////////////////////////////////////////////
          return std::make_tuple<colltype_out>(std::move(output_digi_hits));
      }



    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       finalize       //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    StatusCode ARCalg::finalize() { 
        
        

        if(m_create_debug_histos.value()){
            std::unique_ptr<TFile> f(TFile::Open("Reconstruction.root", "RECREATE"));

          

            hThetaRecoTrue->Write();
            hThetaReco->Write();
            hThetaRecoEm->Write();
            hThetaErrorEm->Write();
            hPixelID->Write();
            hTrueHit->Write();
            hPixelHit->Write();
            hPixelError->Write();
            hThetaRecoPixel->Write();
            hPhotonYield->Write();
            hEmissionZTrue->Write();
            hEmissionXTrue->Write();
            hEmissionYTrue->Write();
            hPhotonYieldGas->Write();
            hPhotonYieldAerogel->Write();
            hThetaGas->Write();
            hThetaAerogel->Write();
            hThetaGen->Write();
            hRecoGen->Write();
            hEMerrorGas->Write();
            hEMerrorAerogel->Write();
            
            }
        
        
        return StatusCode::SUCCESS; }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       ThrowException       ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    void ARCalg::ThrowException(std::string s) const {
        error() << s.c_str()  << endmsg;
        throw std::runtime_error(s);
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       PrintConfiguration       ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    void ARCalg::PrintConfiguration(std::ostream& io)
    {
    io << "ARCalg will use the following components:\n";
    io << "\tGeometry Service: "                  << m_geoSvcName.value().c_str()           << "\n";
    io << "\tDetector name: "                     << m_ARC_name.value().c_str()             << "\n";
    io << "\t\t|--Volume bitfield: "              << m_decoder->fieldDescription().c_str()  << "\n";
    io << "\t\t|--Mirror center (cm): "           << Form("(%g, %g, %g)",
                                                            mirrorCenter.x(),
                                                                mirrorCenter.y(),
                                                                    mirrorCenter.z())       << "\n";
    io << "\t\t|--Mirror radius (cm): "           << mirrorRadius                           << "\n";
    return;
    }
