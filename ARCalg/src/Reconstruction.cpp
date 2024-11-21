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
	double angle[2];
	double angle_pi[2];
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
	//std::cout<<"aAngle "<<aAngle<<std::endl;
	angle[i]=aAngle;
	
        ROOT::Math::XYZVector EM_pi = M_pi - EmissionPoint;
        ROOT::Math::XYZVector p_pi = EM_pi.Unit();

        double    pmag_pi     = pow( EM_pi.Mag2(), 0.5 );
        double    tmag_pi     = pow( Trace.Mag2(), 0.5 );
        double    aAngle_pi   = ( ( pmag_pi * tmag_pi ) != 0.0 ) ? acos( ( EM_pi.Dot( Trace ) ) / ( pmag_pi * tmag_pi ) ) : 0.0;
	//std::cout<<"aAngle_pi "<<aAngle_pi<<std::endl;
        angle_pi[i]=aAngle_pi;
    }


//QUI c'è la condizione hard coded con 0.05
    double n_gas = 1.00135; //C4F10 refractive index
    double n_aerogel = 1.025; //aerogel refractive index
    double mass = 0.4937; //kaon mass in GeV
    double momentum = 50; //da passare in maniera più efficiente
    double beta = momentum/sqrt(pow(mass,2)+pow(momentum,2));
    std::cout<<"beta "<<beta<<std::endl;
    double angle_gas = acos(1/(beta*n_gas)); //già in radianti
    double angle_aerogel = acos(1/(beta*n_aerogel));
    std::cout<<"angle_gas "<<angle_gas<<std::endl;
    std::cout<<"angle_aerogel "<<angle_aerogel<<std::endl;
    
    std::cout<<"angle[0] "<<angle[0]<<std::endl;
    std::cout<<"angle[1] "<<angle[1]<<std::endl;
    
    ROOT::Math::XYZPoint punto(0.0, 0.0, 0.0);
    if (n_real_sol >=0){
    std::array<double, 4> differences = {
        TMath::Abs(angle[0] - angle_gas),
        TMath::Abs(angle[1] - angle_gas),
        TMath::Abs(angle[0] - angle_aerogel),
        TMath::Abs(angle[1] - angle_aerogel)
    };

    // Find the minimum difference and its index
    auto min_it = std::min_element(differences.begin(), differences.end());
    int min_index = std::distance(differences.begin(), min_it);
    std::cout<<"diff[i]"<<differences[1]<< ", "<< differences[2] <<", "<<  differences[3] <<", "<<  differences[4]<<std::endl;
    std::cout<<"min index "<<min_index<<std::endl;

    //mettere if con theta true se la differenza è maggiore di tot flag (1.5 mrad)

    //counter per vedere da dove vengono i fotoni
    //però non può stare qui perchè ne facciamo uno alla volta
/*
    int gamma_gas = 0; 
    int gamma_aerogel = 0;
    for (int i = 0; i < 10000; i++){
        if (min_index == 0 || min_index == 1) {
            gamma_gas += 1;
        } else {
            gamma_aerogel += 1;
        }
    }
*/

    // Determine corresponding delta vector
    ROOT::Math::XYZVector delta_fin;
    if (min_index < 2) {
        delta_fin = delta[min_index];
    } else {
        delta_fin = delta[min_index - 2];
    }
        punto = CoC + delta_fin;
    }
/*
    double angle_fin;
    ROOT::Math::XYZPoint punto(0.0, 0.0, 0.0);
    if(n_real_sol >=0){
        ROOT::Math::XYZVector delta_fin;
	double tmp1,tmp2;
	tmp1 = TMath::Abs(angle[0]-angle_gas);
	tmp2 = TMath::Abs(angle[1]-angle_gas);
    pippo1 = TMath::Abs(angle[0]-angle_aerogel);
    pippo2 = TMath::Abs(angle[1]-angle_aerogel);
    
	if(tmp1<tmp2){
	  //delta_fin = ( ( n_real_sol == 0 ) || ( delta[0].z() > delta[1].z() ) ) ? delta[0] : delta[1];
	  delta_fin=delta[0];
	}else{
	  delta_fin=delta[1];
	}
        punto = CoC + delta_fin;
    }
*/

    ROOT::Math::XYZVector EM_fin = (punto - EmissionPoint);
    ROOT::Math::XYZVector p_fin = EM_fin.Unit();
    double    pmag_fin     = pow( EM_fin.Mag2(), 0.5 );
    double    tmag_fin     = pow( Trace.Mag2(), 0.5 );
    double    aAngle_fin   = ( ( pmag_fin * tmag_fin ) != 0.0 ) ? acos( ( EM_fin.Dot( Trace ) ) / ( pmag_fin * tmag_fin ) ) : 0.0;
    
    
    
    //std::cout<<"angle fin "<<aAngle_fin<<std::endl;
    return aAngle_fin;//theta;
    }


