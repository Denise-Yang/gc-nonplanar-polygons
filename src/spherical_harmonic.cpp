#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "polyscope/point_cloud.h"


#include "geometrycentral/surface/direction_fields.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"
#include "polyscope/curve_network.h"
#include "harnack.hpp"


#include "args/args.hxx"
#include "imgui.h"
#include <cmath>

using namespace geometrycentral;
using namespace geometrycentral::surface;


double sphereicalHarmonic(Eigen::Vector3d pt){
    double x = pt[0];
    double y = pt[1];
    double z = pt[2];

    double res = pow(x,2.f) + pow(y,2.f) - z;
    return res;

}

double traceSphericalHarmonic(Eigen::Vector3d o, Eigen::Vector3d dir, double r){


}