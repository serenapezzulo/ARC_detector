#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include "Math/Point3D.h"
#include "Math/Vector3D.h"

double reconstruct(double radius, ROOT::Math::XYZPoint EmissionPoint, 
                   ROOT::Math::XYZPoint DetectionPoint, ROOT::Math::XYZVector Trace);

#endif