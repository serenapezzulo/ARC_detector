#ifndef TEST_H
#define TEST_H

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "TGeoBBox.h"

// Function to calculate the intersection point between a track and the surface of a TGeoBBox
ROOT::Math::XYZPoint CalculateIntersectionPoint(TGeoBBox* box, const ROOT::Math::XYZPoint& localPoint, const ROOT::Math::XYZVector& localDir);

#endif // TEST_H
