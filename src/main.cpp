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

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
polyscope::VolumeMeshVertexScalarQuantity* scalarQ;
polyscope::CurveNetworkNodeScalarQuantity* rayScalarQ;
polyscope::PointCloud* ptCloud;
polyscope::PointCloud* harnackPtCloud;
polyscope::PointCloud* spherePtCloud;
polyscope::PointCloud* boundingSphere;


// Polyscope visualization handle, to quickly add data to the surface

polyscope::VolumeMesh *psMesh;
// Some algorithm parameters
float param1 = 42.0;

// Example computation function -- this one computes and registers a scalar
// quantity
void doWork() {
  polyscope::warning("Computing Gaussian curvature.\nalso, parameter value = " +
                     std::to_string(param1));

  geometry->requireVertexGaussianCurvatures();
 
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
}



// numbered like in this diagram
// https://polyscope.run/media/tet_element_orderings.jpg
// clang-format off
const std::array<std::array<size_t, 4>, 6> stencilHex =
  {
   std::array<size_t, 4>{2,1,0,3},
   std::array<size_t, 4>{4,0,1,5},
   std::array<size_t, 4>{5,1,2,6},
   std::array<size_t, 4>{7,3,0,4},
   std::array<size_t, 4>{6,2,3,7},
   std::array<size_t, 4>{7,4,5,6}
  };
// clang-format on

// == Geometry data
double wx = 2, wy = 2, wz = 2;
size_t N = 20;
size_t numPoints = 4;
Eigen::MatrixXd gridPoints;
Eigen::MatrixXi gridHexes;
Eigen::MatrixXd planePoints;
Eigen::MatrixXi planeHexes;

Eigen::MatrixXd curvePoints;

class SolidAngleResults{
  public:
    double shifted;   //range [-2π, 2π]
    double unshifted; //range [0, 4π]
};

double getMaxStep(double fx, double R, double lo_bound, double up_bound){
    double v = up_bound/fx;
    double w = lo_bound/fx;
    double lo_r = (-(2*w+1) + sqrt(8*w + 1))*R/(2*w);
    double up_r = (2*v+1 - sqrt(8*v + 1))*R/(2*v);
    // printf("\tup bound: %f lo bound: %f\n", v, w);
    return std::min(up_r, lo_r);
}

double getDis(Eigen::Vector3d pt1, Eigen::Vector3d pt2){
    return sqrt(pow(pt1[0]-pt2[0],2)  + pow(pt2[1]-pt2[1],2) + pow(pt1[2]-pt2[2],2));
}

double getMinDis(Eigen::Vector3d pt, Eigen::MatrixXd curvePoints){
    double minDis = DBL_MAX;
    for (Eigen::Index i = 0; i < curvePoints.rows();i++){
        
        Eigen::Vector3d p0 = curvePoints.row(i);
        double dis = getDis(pt, p0);
        if (dis < minDis){
            minDis = dis;
        }
    }
    return minDis;

}


double sphereicalHarmonic(Eigen::Vector3d pt){
    double x = pt[0];
    double y = pt[1];
    double z = pt[2];

    double res = pow(x,2.f) - pow(y,2.f) + z + 1;
    return res;

}

/*Manually adjusted the range to be between [-2pi,2pi]*/
SolidAngleResults calculateSolidAngle(Eigen::Vector3d x, Eigen::MatrixXd curvePoints){
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
    exteriorAngleSum += theta;
  }
  SolidAngleResults angles = SolidAngleResults();

  double solidAng = (2*M_PI  + exteriorAngleSum) - 2*M_PI;
  angles.unshifted = solidAng + 2*M_PI; 
  // printf("solid anlge %f\n",solidAng);
  // double solidAng = ((2.f-numPoints)*M_PI+exteriorAngleSum);
  if (solidAng < 0.){
    angles.shifted =  4*M_PI + solidAng;
  }else{
    angles.shifted = solidAng;
  }
 
  return angles;
}

void initCurve(){
    float scale = .5;
    curvePoints = Eigen::MatrixXd(numPoints, 3);
    double x,y,z;
    for (size_t i = 0; i < numPoints; i++){
      if (i%2) {
        z = 1.f;
        x = 1.f;
        y = 1.f;
      }
      else {
        z = 0;
        x = -1.f;
        y = 1.f;}
      if (i >= numPoints/2){
        x *= -1.f;
        y *= -1.f;
      }
      curvePoints(i,0) = x*scale;
      curvePoints(i,1) = y*scale;
      curvePoints(i,2) = z*scale;
      // printf("created point %f %f %f\n",x,y,z);
    }
}




    
// double traceSphericalHarmonic(Eigen::Vector3d o, Eigen::Vector3d dir, double r){

float intersect(Eigen::Vector3d origin, Eigen::Vector3d dir) 
{
        float t0, t1,t; // solutions for t if the ray intersects
        // analytic solution
        Eigen::Vector3d L = origin - Eigen::Vector3d(0.f,0.f,1.f);
        float a = dir.dot(dir);
        float b = 2 * dir.dot(L);
        float c = L.dot(L) - 1.f;
        float discr = b * b - 4 * a * c;
        printf("discr %f\n", discr);
        if (discr < 0.f) return -1.f;
        else if (discr == 0.f) t0 = t1 = - 0.5 * b / a;
        else {
          printf("hit this case\n");
          float q = (b > 0) ?
            -0.5 * (b + sqrt(discr)) :
            -0.5 * (b - sqrt(discr));
            t0 = q / a;
            t1 =  c / q;
        }
        if (t0 > t1) std::swap(t0, t1);
    
        if (t0 < 0) {
            t0 = t1; // if t0 is negative, let's use t1 instead
        }

        t = t0;
        printf("final t0 %f  (%f, %f)\n", t, t0, t1);

        return t+.1;
}

// }
std::vector<HSphere> harnackTrace(Eigen::Vector3d origin, Eigen::Vector3d dir, Eigen::MatrixXd curvePoints){
    std::vector<HSphere> spheres;
    int curr_sample = 0;
    int max_samples = 100;
    double loBound = .0001;
    double hiBound = 2.f;
    const double mag = dir.norm();
    dir /= mag;
    Eigen::Vector3d v = origin; 
    double r = 0.f;
    float t = intersect(v, dir);
    if (t<0.f){
      return spheres;
    }
    v += dir*t;
    printf("sH: %f\n", sphereicalHarmonic(v));
    if (sphereicalHarmonic(v) > 2.f){
      loBound = 2.f;
      hiBound = 4.f;
      
    }
    while (curr_sample < max_samples){
       
        Eigen::Vector3d shift = Eigen::Vector3d(0,0,1);
        Eigen::Vector3d center = Eigen::Vector3d(0,0,0);
        double fx = sphereicalHarmonic(v);
        double R = 1.f - getDis(v, center + shift);
        // printf("point (%f %f,%f): val %f", v[0],v[1],v[2], mag);
        r = getMaxStep(fx, R, loBound, hiBound);      
        v += r*dir;
        spheres.push_back(HSphere(r,R,v));
        curr_sample++;
    }
    return spheres;

}

void initGeom() {
    gridPoints = Eigen::MatrixXd(N * N * N, 3);
    gridHexes  = Eigen::MatrixXi((N - 1) * (N - 1) * (N - 1), 8);

    double dx = wx / (N - 1);
    double dy = wy / (N - 1);
    double dz = wz / (N - 1);
    for (size_t iX = 0; iX < N; iX++) {
        for (size_t iY = 0; iY < N; iY++) {
            for (size_t iZ = 0; iZ < N; iZ++) {
                double x = dx * iX - wx / 2.;
                double y = dy * iY - wx / 2.;
                double z = dz * iZ - wx / 2. + 1;

                size_t iV         = N * N * iX + N * iY + iZ;
                gridPoints(iV, 0) = x;
                gridPoints(iV, 1) = y;
                gridPoints(iV, 2) = z;

                if (iX + 1 >= N || iY + 1 >= N || iZ + 1 >= N) continue;

                size_t iC        = (N - 1) * (N - 1) * iX + (N - 1) * iY + iZ;
                gridHexes(iC, 0) = N * N * (iX + 0) + N * (iY + 0) + (iZ + 0);
                gridHexes(iC, 1) = N * N * (iX + 1) + N * (iY + 0) + (iZ + 0);
                gridHexes(iC, 2) = N * N * (iX + 1) + N * (iY + 1) + (iZ + 0);
                gridHexes(iC, 3) = N * N * (iX + 0) + N * (iY + 1) + (iZ + 0);
                gridHexes(iC, 4) = N * N * (iX + 0) + N * (iY + 0) + (iZ + 1);
                gridHexes(iC, 5) = N * N * (iX + 1) + N * (iY + 0) + (iZ + 1);
                gridHexes(iC, 6) = N * N * (iX + 1) + N * (iY + 1) + (iZ + 1);
                gridHexes(iC, 7) = N * N * (iX + 0) + N * (iY + 1) + (iZ + 1);
            }
        }
    }
}

int main(int argc, char **argv) {

  // Initialize polyscope
  initGeom();
  initCurve();

   // Initialize polyscope
  polyscope::init();

  Eigen::Vector3d origin = Eigen::Vector3d(0,0,3);
  Eigen::Vector3d dir = Eigen::Vector3d(0,0,-1.f);
 
  // Eigen::Vector3d origin = Eigen::Vector3d(.45,0,.75);
  // Eigen::Vector3d dir = Eigen::Vector3d(-3,0,2);
 

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // std::vector<HSphere> harnackPoints = harnackTrace(origin,dir,curvePoints);
  std::vector<Eigen::VectorXd> samplePoints;
  std::vector<double> harnackRadius;
  std::vector<double> sphereRadius;
  std::vector<Eigen::VectorXd> bounds;


  // for (auto &h : harnackPoints){
  //   samplePoints.push_back(h.pt);
  //   harnackRadius.push_back(h.r);
  //   sphereRadius.push_back(h.R);

  // }

  
  // ptCloud = polyscope::registerPointCloud("rayPoints", samplePoints);
  // harnackPtCloud = polyscope::registerPointCloud("harnackSpheres", samplePoints);

  //enable visualization for max harnack step. 
  // spherePtCloud = polyscope::registerPointCloud("sphereRadii", samplePoints);
  // auto sr = spherePtCloud->addScalarQuantity("sphereRadius", sphereRadius); // add the quantity
  // spherePtCloud->setPointRadiusQuantity(sr, false);

  // harnackPtCloud->setTransparency(0.25);
  // auto hr = harnackPtCloud->addScalarQuantity("harnackRadius", harnackRadius); // add the quantity
  // harnackPtCloud->setPointRadiusQuantity(hr, false);

  //create visual spherical bound for harmmonic function
  // bounds.push_back(Eigen::Vector3d(0,0,1));
  // boundingSphere = polyscope::registerPointCloud("boundingSphere", bounds);
  
  // boundingSphere->setPointRadius(1.f,false);
  // boundingSphere->setTransparency(0.5);
  


  Eigen::VectorXd f(gridPoints.rows());
  for (Eigen::Index iV = 0; iV < gridPoints.rows(); iV++) {
      //f(iV) = std::fmod(calculateSolidAngle( gridPoints.row(iV), tri0) + calculateSolidAngle( gridPoints.row(iV), tri1),(4*M_PI));
      // SolidAngleResults angles = calculateSolidAngle( gridPoints.row(iV), curvePoints);
      double res = sphereicalHarmonic( gridPoints.row(iV));
      f(iV) = res;
  }

  // Register the curve with polyscope
  // polyscope::registerCurveNetworkLoop("my network", curvePoints);
  psMesh = polyscope::registerHexMesh(
      "my mesh",
      gridPoints,gridHexes);
 
  scalarQ = psMesh->addVertexScalarQuantity("scalar Q", f);
  scalarQ->setEnabledLevelSet(true);
  scalarQ->setLevelSetValue(1.f);

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}


/*
//HARNACK TRACE FOR SOLID ANGLE
std::vector<HSphere> harnackTrace(Eigen::Vector3d origin, Eigen::Vector3d dir, Eigen::MatrixXd curvePoints){
    std::vector<HSphere> spheres;
    int curr_sample = 0;
    int max_samples = 50;
    const double mag = dir.norm();
    dir /= mag;
    // double theta = atan2(sqrt(pow(dir[0],2) + pow(dir[1],2)),dir[2])+M_PI;
    // double ro = atan2(dir[1], dir[0]) + M_PI;
    // printf("atan2(0): %f, atan2(1): %f\n",atan2(-0,1),atan2(-1,-1) );

    Eigen::Vector3d v = origin; 
    double r = 0.f;
    while (curr_sample < max_samples){
        //Add 2π to make solid angle positive
        //Add 4π to shift the function after the discontinuity above 0
        SolidAngleResults angles = calculateSolidAngle(v, curvePoints);
        double fx = angles.shifted + 4*M_PI;
        double loBound = 4*M_PI;
        double hiBound = 8*M_PI;
        \\R =  min dis between point and curvePoints
        r = getMaxStep(fx, R, loBound, hiBound);      
        v += r*dir;
        spheres.push_back(HSphere(r,R,v));
        curr_sample++;
    }
    return spheres;

}
*/