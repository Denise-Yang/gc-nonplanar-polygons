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
#include <float.h>


//Do a solid angle function as harmonic function
//calculate harnack's bound in 3d
// raytrace using harnack's bound --> let's get that oworking for next week
// add the 100th point to the thing and plot

double getMaxStep(double fx, double R, double lo_bound, double up_bound){
    double v = up_bound/fx;
    double w = lo_bound/fx;
    double lo_r = (-(2w+1) + sqrt(8*w + 1)) * R/(2*w);
    double up_r = (2v+1 - sqrt(8*v + 1)) * R/(2*v);
    return std::min(up_r, lo_r)
}

double getMinDis(Eigen::Vector3d pt, Eigen::MatrixXd curvePoints){
    Eigen::Vector3d minP;
    double minDis = DBL_MAX;
    for (Eigen::Index i = 0; i < curvePoints.rows();i++){
        
        Eigen::Vector3d p0 = curvePoints.row(i);
        double distance = sqrt(pow(pt[0]-p0[0],2)  + pow(pt[1]-p0[1],2) + pow(pt[2]-p0[2],2))
        if (dis < minDis){
            minP = p0;
            minDis = dis;
        }
    }
    return minP;

}

double calculateSA(Eigen::Vector3d x, Eigen::MatrixXd curvePoints){
    double exteriorAngleSum = 0.;
    float dir = 0.f;
    for (Eigen::Index i = 0; i < curvePoints.rows();i++){

        Eigen::Vector3d p0, p1, p2, vX0, vX1, vX2, v, u, c;
        float d, dir;

        if (i - 1 < 0){
            p0 = curvePoints.row(curvePoints.rows()-1);
        }else{
            p0 = curvePoints.row(i - 1);
        }
        p1 = curvePoints.row(i);
        p2 = curvePoints.row((i+1)%curvePoints.rows());

        vX0 = p0-x;
        vX1 = p1-x;
        vX2 = p2-x; 
        v = vX0.cross(vX1).normalized();
        u = vX1.cross(vX2).normalized();

        d = u.dot(v);
        c = u.cross(v);
        dir = c.dot(vX1)/vX1.norm();
        double theta =  atan2(dir, d);
        // if (theta < 0) theta += 2*M_PI;
        exteriorAngleSum += M_PI - theta;
    }
    // double solidAng = (2*M_PI  + exteriorAngleSum);
    // printf("solid anlge %f\n",solidAng);
    double solidAng = ((2.f-numPoints)*M_PI+exteriorAngleSum);
    //   if (solidAng < 0.){
    //     return 4*M_PI + solidAng;
    //   }
    return 4*M_PI  + solidAng;
}

std::vector<HSphere> harnackTrace(Eigen::Vector3d origin, Eigen::Vector3d dir, Eigen::MatrixXd curvePoints){
    std::vector<HSphere> spheres;
    int curr_sample = 0;
    int max_sample = 50;
    double theta = atan2(sqrt(math.pow(dir[0],2) + math.pow(dir[1],2))/dir[2]);
    double ro = atan2(dir[1]/dir[0]);

    Eigen::Vector3d v = origin; 
    double r = 0.f;
    while (cur_sample < max_samples){
        double fx = calculateSA(v, curvePoints);
        double R = getMinDis(v, curvePoints)
        r = get_r(fx, R, 4*M_PI, 6*M_PI);
        double xStep = r * sin(theta) * cos(ro);
        double yStep = r * sin(theta) * sin(ro);
        double zStep = r * cos(theta);
        v[0] += xStep;
        v[1] += yStep;
        v[2] += zStep;
        spheres.push_back(HSphere(r,R,v));
        curr_sample++;
    }
    return spheres;

}