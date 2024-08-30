#include <iostream>
#include <limits>  // For std::numeric_limits
#include "TGeoBBox.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"

// Function to calculate the intersection point between a track and the surface of a TGeoBBox
ROOT::Math::XYZPoint CalculateIntersectionPoint(TGeoBBox* box, const ROOT::Math::XYZPoint& localPoint, const ROOT::Math::XYZVector& localDir) {
    ROOT::Math::XYZPoint intersectionPoint;
    double minDistance = std::numeric_limits<double>::max();

    // Get half-dimensions of the TGeoBBox
    double dx = box->GetDX();
    double dy = box->GetDY();
    double dz = box->GetDZ();

    // Check for intersection with each face of the bounding box
    for (int i = 0; i < 3; i++) {
        double d1, d2;
        double distance1 = std::numeric_limits<double>::max();
        double distance2 = std::numeric_limits<double>::max();

        if (i == 0) {  // X-axis planes
            d1 = dx;  // x = dx
            d2 = -dx; // x = -dx
            if (localDir.X() != 0) {
                distance1 = (d1 - localPoint.X()) / localDir.X();
                distance2 = (d2 - localPoint.X()) / localDir.X();
            }
        } else if (i == 1) {  // Y-axis planes
            d1 = dy;  // y = dy
            d2 = -dy; // y = -dy
            if (localDir.Y() != 0) {
                distance1 = (d1 - localPoint.Y()) / localDir.Y();
                distance2 = (d2 - localPoint.Y()) / localDir.Y();
            }
        } else if (i == 2) {  // Z-axis planes
            d1 = dz;  // z = dz
            d2 = -dz; // z = -dz
            if (localDir.Z() != 0) {
                distance1 = (d1 - localPoint.Z()) / localDir.Z();
                distance2 = (d2 - localPoint.Z()) / localDir.Z();
            }
        }

        // For each intersection, check if it is within the bounds of the other two dimensions
        std::vector<std::pair<double, ROOT::Math::XYZPoint>> candidates = {
            {distance1, localPoint + distance1 * localDir},
            {distance2, localPoint + distance2 * localDir}
        };

        for (const auto& [distance, point] : candidates) {
            if (distance > 0 && distance < minDistance) {
                // Check if the intersection point is within the box bounds
                if (std::abs(point.X()) <= dx && std::abs(point.Y()) <= dy && std::abs(point.Z()) <= dz) {
                    minDistance = distance;
                    intersectionPoint = point;
                }
            }
        }
    }

    return intersectionPoint;
}
