    #include "ARCalg.h"
    #include "ArcQuarticEquation.h"
    #include "Reconstruction.h"
    //#include "Intersection.h"
    #include "test.h"
    #include <TVector3.h>

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


        std::stringstream ss;
        PrintConfiguration(ss);
        info() << ss.str().c_str() <<endmsg;

        if(m_create_debug_histos.value()){
            hThetaRecoTrue = new TH1D("hThetaRecoTrue", "Reconstructed Theta", 100, 0, 0.07);
            hThetaRecoTrue->SetDirectory(0);
            hThetaRecoEm = new TH1D("hThetaRecoEm", "Reconstructed Theta", 100, 0, 0.07);
            hThetaRecoEm->SetDirectory(0);
            hThetaErrorEm = new TH1D("hThetaErrorEm", "Reconstructed Theta Error", 100, -0.05, 0.05);
            hThetaErrorEm->SetDirectory(0);
            hPixelID = new TH2F("hPixelID", "Pixel ID", 100,0,100,100,0,100 );
            hPixelID->SetDirectory(0);
            hTrueHit = new TH2F("hTrueHit", "Pixel Position", 100, -5, 5, 100, -5, 5);
            hTrueHit->SetDirectory(0);
            hPixelHit = new TH2F("hPixelHit", "Pixel Hit (center)", 100, -5, 5, 100, -5, 5);
            hPixelHit->SetDirectory(0);
            hPixelError = new TH1D("hPixelError", "Pixel Error", 100, -0.05, 0.05);
            hPixelError->SetDirectory(0);
            hThetaRecoPixel = new TH1D("hThetaRecoPixel", "Reconstructed Theta", 100, 0, 0.07);
            hThetaRecoPixel->SetDirectory(0);
            hPhotonYield = new TH1D("hPhotonYield", "Photon Yield", 100, 0, 150);
            hPhotonYield->SetDirectory(0);

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

        //loop over hit collection
        for (const auto& input_sim_hit : input_sim_hits) {
        dd4hep::DDSegmentation::CellID id = input_sim_hit.getCellID();
         //pixels ID -> integer numbers
        double pixel_x = m_decoder->get(id, "x");
        double pixel_y = m_decoder->get(id, "y");
        auto pixelcenter_x = -3.96 + pixel_x*0.08;
        auto pixelcenter_y = -3.96 + pixel_y*0.08;
       
        
        auto truehit_x = input_sim_hit.getPosition().x;
        auto truehit_y = input_sim_hit.getPosition().y;


       const edm4hep::MCParticle& photon = input_sim_hit.getMCParticle();

        // Get the emission point (vertex) of the photon
        auto emissionPoint = photon.getVertex();

        ROOT::Math::XYZPoint EmissionPointTrue(emissionPoint.x, emissionPoint.y, emissionPoint.z);

     
        // Convert edm4hep::Vector3d to ROOT::Math::XYZVector
        ROOT::Math::XYZPoint DetectionPointTrue(input_sim_hit.getPosition().x,
                                  input_sim_hit.getPosition().y,
                                  input_sim_hit.getPosition().z);

        ROOT::Math::XYZVector Trace(input_mc_particles.momentum().at(0).x,
                            input_mc_particles.momentum().at(0).y,
                            input_mc_particles.momentum().at(0).z);


        ROOT::Math::XYZPoint DetectionPointPixel = ROOT::Math::XYZPoint(pixelcenter_x*10, pixelcenter_y*10, input_sim_hit.getPosition().z);
        

        double TrueAngle = reconstruct(mirrorRadius, EmissionPointTrue, DetectionPointTrue, Trace);
        
        
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
        ROOT::Math::XYZPoint ShootingPointTrue(shootingPoint.x, shootingPoint.y, shootingPoint.z);


        // Calculate the intersection between the track and the dummy box
        ROOT::Math::XYZPoint EmissionPointApprox = CalculateIntersectionPoint(dummyBox, ShootingPointTrue, Trace);

        
        double Angle_em = reconstruct(mirrorRadius, EmissionPointApprox, DetectionPointTrue, Trace);
        
    
        std::cout << "Trace: " << Trace << std::endl;
        //std::cout << "entering point" << enteringpoint << std::endl;
        std::cout << "EmissionPointApprox: " << EmissionPointApprox << std::endl;
        std::cout << "EmissionPointTrue: " << EmissionPointTrue << std::endl;
        
        double Angle_pixel = reconstruct(mirrorRadius, EmissionPointTrue, DetectionPointPixel, Trace);
        




        if(m_create_debug_histos.value()){
            hThetaRecoTrue->Fill(TrueAngle);
            hThetaRecoEm->Fill(Angle_em);
            hThetaErrorEm->Fill(Angle_em - TrueAngle);
            hPixelID->Fill(pixel_x, pixel_y);
            hTrueHit->Fill(truehit_x*0.1, truehit_y*0.1);
            hPixelHit->Fill(pixelcenter_x, pixelcenter_y);
            hThetaRecoPixel->Fill(Angle_pixel);
            hPixelError->Fill(Angle_pixel - TrueAngle);
            hPhotonYield->Fill(input_sim_hits.size());

        

            
        }
    
       
        
        // emplace output objects
        // output_digi_hits.create(....);

        }// end loop over hit collection

        

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
            hThetaRecoEm->Write();
            hThetaErrorEm->Write();
            hPixelID->Write();
            hTrueHit->Write();
            hPixelHit->Write();
            hPixelError->Write();
            hThetaRecoPixel->Write();
            hPhotonYield->Write();
        
            
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
