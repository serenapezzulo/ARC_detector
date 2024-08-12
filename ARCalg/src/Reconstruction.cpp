#include "Reconstruction.h"
#include <iostream>
#include <cmath> 
#include <gsl/gsl_complex.h>
#include "ArcQuarticEquation.h"
   // ROOT
    #include "TGeoBoolNode.h"
    #include <TGeoManager.h>


    // Math
    #include "Math/Math.h"
    #include "Math/Rotation3D.h"
    #include "Math/RotationX.h"
    #include "Math/RotationY.h"
    #include "Math/RotationZ.h"
    #include "Math/Transform3D.h"
    #include "Math/Vector3D.h"
    #include "Math/AxisAngle.h"
 

    // STL
    #include <iostream>
    #include <sstream>



    ArcQuarticEquation* m_ArcQuarticEquation;

    double mirrorRadius(367); 
    ROOT::Math::XYZPoint MirrorCoC(0, 0, -279);

    double reconstruct(double radius, ROOT::Math::XYZPoint EmissionPoint, ROOT::Math::XYZPoint DetectionPoint, ROOT::Math::XYZVector Trace){


    //auto mirror_ptr = get_mirror_properties(tr, geoSvc);
    ROOT::Math::XYZPoint CoC=MirrorCoC;
    //std::cout<<"CoC "<<CoC<<endl;
    ROOT::Math::XYZVector ev = (EmissionPoint - CoC);
    ROOT::Math::XYZVector dv = (DetectionPoint - CoC);
    double d2 = dv.Mag2();
    double e2 = ev.Mag2();
    double d = std::sqrt(d2);
    double e = std::sqrt(e2);
    

    double r=mirrorRadius;
    double r2=r*r;
    //double gamma = (ev).Angle(dv); 
    double gamma = acos(ev.Dot(dv) / (e * d));
    double dy = d*sin(gamma);
    double dx = d*cos(gamma);
    
    //std::cout << "gamma:" << gamma << std::endl; 
    //std::cout << "dy: " << dy << std::endl;
    //std::cout << "dx: " << dy << std::endl;
    //std::cout << "e: " << e << std::endl;
    double den=4*e2*d2;
    double coeff[4]={(-4 * e2 * dy * r), (dy * dy * r2 + ( e + dx ) * ( e + dx ) * r2 - den),
        (2 * e * dy * ( e - dx ) * r), (( e2 - r2 ) * dy * dy)};
    
    if( den !=0){
        for(int i=0;i<4;i++){
        coeff[i]/=den;
        }
    }
    gsl_complex solq[4];
    m_ArcQuarticEquation->gsl_poly_complex_solve_quartic(coeff[0],coeff[1],coeff[2],coeff[3], &solq[0], &solq[1], &solq[2], &solq[3]);

    //store real solutions in an array
    double real_sol[4]; //me ne aspetto solo 2 ma creo array per tenerle tutte nel caso non funzionasse
    int n_real_sol=0;
    
    
    for(int i=0;i<4;i++){
        if(GSL_IMAG(solq[i])==0 && GSL_REAL(solq[i]) <=1.0){
        real_sol[n_real_sol]=GSL_REAL(solq[i]);
        n_real_sol++;
        }
    }
    if(n_real_sol==0){
        //std::cout << "No real solutions" << std::endl;
        return -1;
    }

        //std::cout << "N. real solution" << n_real_sol <<std::endl;

    for (int i = 0; i < n_real_sol; ++i) {
            //std::cout << "Real sol:"<< real_sol[i] << " ";
        }
        std::cout << std::endl;
        ROOT::Math::XYZVector delta[2]={ROOT::Math::XYZVector(0.0,0.0,0.0), ROOT::Math::XYZVector(0.0,0.0,0.0)};
        ROOT::Math::XYZVector delta_pi[2]={ROOT::Math::XYZVector(0.0,0.0,0.0), ROOT::Math::XYZVector(0.0,0.0,0.0)};

    //tra le 2 soluzioni reali ne prendo una in base alla geometria 
    //proviamo con quella più vicina a (0,0,0)
    std::vector<ROOT::Math::XYZPoint> Mv; //vettore di punti M

    for(int i=0;i<n_real_sol;i++){
        
        double beta = asin(real_sol[i]);
        ROOT::Math::XYZVector u=ev.Unit();
        u *= (r);

        ROOT::Math::XYZVector nvec = ev.Cross(dv); //vettore normale al piano individuato dai vettori ec/dc
        ROOT::Math::Rotation3D rotn = ROOT::Math::Rotation3D(ROOT::Math::AxisAngle(nvec, beta));
        ROOT::Math::Rotation3D rotn_pi = ROOT::Math::Rotation3D(ROOT::Math::AxisAngle(nvec, ROOT::Math::Pi()-beta));
     
        
        ROOT::Math::XYZVector bb=rotn(u);
        ROOT::Math::XYZVector bb_pi = rotn_pi(u);
        
        delta[i]=bb;
        ROOT::Math::XYZPoint M = CoC + delta[i];
        Mv.push_back(M);
        
        delta_pi[i]=bb_pi;
        ROOT::Math::XYZPoint M_pi = CoC + delta_pi[i];
       

    
        ROOT::Math::XYZVector EM = (M - EmissionPoint);
        ROOT::Math::XYZVector p = EM.Unit();
        double    pmag     = pow( EM.Mag2(), 0.5 );
        double    tmag     = pow( Trace.Mag2(), 0.5 );
        double    aAngle   = ( ( pmag * tmag ) != 0.0 ) ? acos( ( EM.Dot( Trace ) ) / ( pmag * tmag ) ) : 0.0;
    
        
        ROOT::Math::XYZVector EM_pi = M_pi - EmissionPoint;
        ROOT::Math::XYZVector p_pi = EM_pi.Unit();

        double    pmag_pi     = pow( EM_pi.Mag2(), 0.5 );
        double    tmag_pi     = pow( Trace.Mag2(), 0.5 );
        double    aAngle_pi   = ( ( pmag_pi * tmag_pi ) != 0.0 ) ? acos( ( EM_pi.Dot( Trace ) ) / ( pmag_pi * tmag_pi ) ) : 0.0;
       
        
    }


    ROOT::Math::XYZPoint punto(0.0, 0.0, 0.0);
    if(n_real_sol >=0){
        ROOT::Math::XYZVector delta_fin;
        delta_fin = ( ( n_real_sol == 0 ) || ( delta[0].z() > delta[1].z() ) ) ? delta[0] : delta[1];
        punto = CoC + delta_fin;
    }


    ROOT::Math::XYZVector EM_fin = (punto - EmissionPoint);
    ROOT::Math::XYZVector p_fin = EM_fin.Unit();
    double    pmag_fin     = pow( EM_fin.Mag2(), 0.5 );
    double    tmag_fin     = pow( Trace.Mag2(), 0.5 );
    double    aAngle_fin   = ( ( pmag_fin * tmag_fin ) != 0.0 ) ? acos( ( EM_fin.Dot( Trace ) ) / ( pmag_fin * tmag_fin ) ) : 0.0;
    
    
    

    return aAngle_fin;//theta;
    }


