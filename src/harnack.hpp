#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "polyscope/point_cloud.h"


#include "geometrycentral/surface/direction_fields.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"
#include "polyscope/curve_network.h"

#include "args/args.hxx"
#include "imgui.h"
#include <cmath>
#include <float.h>

class HSphere{
    public:
        double R;
        double r;
        Eigen::Vector3d pt;
    HSphere(double r, double R,Eigen::Vector3d pt ){
        this->R = R;
        this->r = r;
        this->pt = pt;
    }
};
double getMaxStep(double fx, double R, double lo_bound, double up_bound);
double getMinDis(Eigen::Vector3d pt, Eigen::MatrixXd curvePoints);
double calculateSA(Eigen::VectorXd x, Eigen::MatrixXd curvePoints);
std::vector<HSphere> harnackTrace(Eigen::Vector3d origin, Eigen::Vector3d dir, Eigen::MatrixXd curvePoints);

