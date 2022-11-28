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
polyscope::PointCloud* psCloud;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::VolumeMesh *psMesh;

// Some algorithm parameters
float param1 = 42.0;

/*TODO
    1. Find point of intersection along ray and 2pi level set
    2. Find the gradient via Bios Savart
    3. MARIOOOOOOO
*/




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
double wx = 2, wy = 2, wz = 6;
size_t N = 20;
size_t numPoints = 4;
Eigen::MatrixXd gridPoints;
Eigen::MatrixXi gridHexes;
Eigen::MatrixXd planePoints;
Eigen::MatrixXi planeHexes;

Eigen::MatrixXd curvePoints;


double calculateSolidAngle(Eigen::Vector3d x, Eigen::MatrixXd curvePoints){
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
  if (solidAng < 0.){
    return 4*M_PI + solidAng;
  }
  return solidAng;
}

std::vector<Eigen::Vector3d> traceRay(Eigen::Vector3d v, Eigen::Vector3d offset, float xrange, float numSamples,Eigen::MatrixXd curvePoints){
  std::vector<Eigen::Vector3d > points;
  float step = xrange/numSamples; 
  for (float t  = 0; t < xrange; t += step){
    Eigen::Vector3d x  = v*t + offset;
    // glm::vec2 point(x);
    points.push_back(x);
  }
  return points;
}


void initCurve(){
    curvePoints = Eigen::MatrixXd(numPoints, 3);
    double x,y,z;
    for (size_t i = 0; i < numPoints; i++){
      if (i%2) {
        z = .5;
        x = .5;
        y = .5;
      }
      else {
        z = 0;
        x = -.5;
        y = .5;}
      if (i >= numPoints/2){
        x *= -1.;
        y *= -1.;
      }
      curvePoints(i,0) = x;
      curvePoints(i,1) = y;
      curvePoints(i,2) = z;
      printf("created point %f %f %f\n",x,y,z);
    }
}

Eigen::MatrixXd initTriangulatedCurve(bool doLeft){
    size_t pointsNum = 3;
    int flipTri = 1;
    size_t start = 0;
    if (doLeft) 
    {
      flipTri = -1;
      start = 1;
    }
    Eigen::MatrixXd tri = Eigen::MatrixXd(pointsNum, 3);
    double x,y,z;
    for (size_t i = start; i < pointsNum +start; i++){
      if (i == start){
        z = .0;
        x = .5*flipTri;
        y = -.5*flipTri;
      }else if(i%2 == 0) {
        z = .5;
        x = .5;
        y = .5;
      }
      else if (i%2==1){
        z = .5;
        x = -.5;
        y = -.5;
      }
      tri(i-start,0) = x;
      tri(i-start,1) = y;
      tri(i-start,2) = z;
      printf("created point %f %f %f\n",x,y,z);
    }
    return tri;
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
                double z = dz * iZ - wx / 2.;

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

void initPlane() {
    planePoints = Eigen::MatrixXd(N * N, 3);

    double dx = wx / (N - 1);
    double dy = wy / (N - 1);
    for (size_t iY = 0; iY < N; iY++) {
        for (size_t iX = 0; iX < N; iX++) {
              double x = dx * iX - wx / 2.;
              double y = dy * iY - wx / 2.;
              double z = 0;

              size_t iV         = N * iY + iX;
              planePoints(iV, 0) = x;
              planePoints(iV, 1) = y;
              planePoints(iV, 2) = z;
              if (iX + 1 >= N || iY + 1 >= N) continue;
        }
    }
}


// Example computation function -- this one computes and registers a scalar
// quantity
void doWork() {
  polyscope::warning("Computing Gaussian curvature.\nalso, parameter value = " +
                     std::to_string(param1));

  geometry->requireVertexGaussianCurvatures();
  psMesh->addVertexScalarQuantity("curvature",
                                  geometry->vertexGaussianCurvatures,
                                  polyscope::DataType::SYMMETRIC);
}


// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    // static bool showLevelSet = true;
    // if (!showLevelSet && ImGui::Button("Show level set")) {
    //     showLevelSet = true;
    //     scalarQ->setEnabledLevelSet(showLevelSet);
    // } else if (showLevelSet && ImGui::Button("Hide level set")) {
    //     showLevelSet = false;
    //     scalarQ->setEnabledLevelSet(showLevelSet);
    // }

    // static float val = 0.3f;
    // if (ImGui::SliderFloat("Level set val", &val, 0.0f, 0.5f)) {
    //     scalarQ->setLevelSetValue(val);
    // }
}

Eigen::Vector3d traceLevelSet(std::vector<Eigen::Vector3d> points, Eigen::MatrixXd curvePoints, float levelSet){
  double SA;
  double prevSA = calculateSolidAngle(points.at(0), curvePoints);
  for (size_t i = 1; i < points.size()-1; i++) {
      SA = calculateSolidAngle(points.at(i), curvePoints);
      //f(iV) = std::fmod(calculateSolidAngle( gridPoints.row(iV), tri0) + calculateSolidAngle( gridPoints.row(iV), tri1),(4*M_PI));
      if ((prevSA >= levelSet && SA <= levelSet) || (SA >= levelSet && prevSA <= levelSet)){
        return (points.at(i-1) + points.at(i))/2;
      }
    }
  printf("didnt intersect\n");
  return Eigen::Vector3d (0,0,-1);

}



Eigen::Vector3d binarySearch(float lo, float hi, float To_Find, Eigen::Vector3d v, Eigen::Vector3d offset, Eigen::MatrixXd curvePoints)
{
    float mid;
    // This below check covers all cases , so need to check
    // for mid=lo-(hi-lo)/2
    float precision = 1.f;
    while (hi - lo > 0) {
        mid = (hi + lo) / 2.f;
        Eigen::Vector3d x  = v*mid + offset;
        float sa = calculateSolidAngle(x, curvePoints); 
        if (sa < To_Find) {
            lo = mid + precision;
        }
        else {
            hi = mid;
        }
    }
    Eigen::Vector3d loX  = v*lo + offset;
    Eigen::Vector3d hiX  = v*hi + offset;
    float saLo = calculateSolidAngle(loX, curvePoints); 
    float saHi = calculateSolidAngle(hiX, curvePoints); 

    if (saLo <= To_Find + precision && saLo >= To_Find - precision) {
        std::cout << "Found"
             << " At Index " << lo << std::endl;
        return loX;

    }
    else if (saHi <= To_Find + precision && saHi >= To_Find - precision) {
        std::cout << "Found"
             << " At Index " << hi << std::endl;
        return hiX;
    }
    else {
        std::cout << "Not Found" << std::endl;
        return Eigen::Vector3d (0,0,-1);
    }
}

int main(int argc, char **argv) {

  // Configure the argument parser
  // args::ArgumentParser parser("geometry-central & Polyscope example project");
  // args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

  // // Parse args
  // try {
  //   parser.ParseCLI(argc, argv);
  // } catch (args::Help &h) {
  //   std::cout << parser;
  //   return 0;
  // } catch (args::ParseError &e) {
  //   std::cerr << e.what() << std::endl;
  //   std::cerr << parser;
  //   return 1;
  // }

  // // Make sure a mesh name was given
  // if (!inputFilename) {
  //   std::cerr << "Please specify a mesh file as argument" << std::endl;
  //   return EXIT_FAILURE;
  // }
  initGeom();
  initCurve();
  printf("registering plane\n");
  initPlane();
  printf("registered plane\n");
  Eigen::MatrixXd tri0 = initTriangulatedCurve(true);
  Eigen::MatrixXd tri1 = initTriangulatedCurve(false);

  


  // Initialize polyscope
  polyscope::init();

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  std::string ray1Name = "ray1";
  Eigen::Vector3d v = Eigen::Vector3d(0,0,.5);
  // Eigen::Vector3d offset = Eigen::Vector3d(0,0,-1);
  // std::vector<Eigen::Vector3d> points = traceRay(v,offset,100,5,curvePoints);
  // polyscope::registerCurveNetworkLine(ray1Name, points);
  std::vector<Eigen::Vector3d> plane2Pi; 

  polyscope::registerCurveNetworkLoop("my network", curvePoints);
  // polyscope::registerCurveNetworkLoop("tri0", tri0);
  // polyscope::registerCurveNetworkLoop("tri1", tri1);


  // std::vector<double> saQ(points.size());
  // for (size_t i = 0; i < points.size(); i++) {
  //   saQ[i] = calculateSolidAngle(points[i], curvePoints); 
  // }

  // visualize
  // rayScalarQ = polyscope::getCurveNetwork(ray1Name)->addNodeScalarQuantity("SolidAngle", saQ);
  // rayScalarQ->setMapRange(std::pair<double,double> (0.0, 12.56637));

  

  // Register the mesh with polyscope
  // psMesh = polyscope::registerSurfaceMesh(
  //     polyscope::guessNiceNameFromPath(args::get(inputFilename)),
  //     geometry->inputVertexPositions, mesh->getFaceVertexList(),
  //     polyscopePermutations(*mesh));

// std::vector<Eigen::Vector3d> plane2Pi; 
//     for (Eigen::Index iV = 0; iV < planePoints.rows(); iV++) {
//       //f(iV) = std::fmod(calculateSolidAngle( gridPoints.row(iV), tri0) + calculateSolidAngle( gridPoints.row(iV), tri1),(4*M_PI));
//       Eigen::Vector3d point2Pi_push.back(binarySearch(0.f, 10.f, 2*M_PI, v, offset, curvePoints));
//       // plane2Pi(iV, 0) = point2Pi[0];
//       // plane2Pi(iV, 1) = point2Pi[1];
//       // plane2Pi(iV, 2) = point2Pi[2];
//     }
    unsigned long counter =0;
    std::vector<Eigen::Vector3d> points;
    std::vector<std::array<unsigned long, 2>> edges;
    for (Eigen::Index iV = 0; iV < planePoints.rows(); iV++) {
      //f(iV) = std::fmod(calculateSolidAngle( gridPoints.row(iV), tri0) + calculateSolidAngle( gridPoints.row(iV), tri1),(4*M_PI));
      std::vector<Eigen::Vector3d> rayPoints = traceRay(v,planePoints.row(iV),2,4,curvePoints);
      // points.insert(std::end(points), std::begin(rayPoints), std::end(rayPoints));
      // std::array<unsigned long, 2>  ptIdx = {counter, counter+1};
      // counter += 2;
      // edges.push_back(ptIdx);
      plane2Pi.push_back(binarySearch(0, 10.f, 2*M_PI, v, planePoints.row(iV), curvePoints));
      // plane2Pi.push_back(traceLevelSet(points, curvePoints, 2*M_PI));

    }
    // polyscope::CurveNetwork* psCurve = polyscope::registerCurveNetwork("ray intersects", points, edges);
    // psCurve->setRadius(.0005);
    psCloud = polyscope::registerPointCloud("2pi points", plane2Pi);


  // Visualize the Solid Angle
  psMesh = polyscope::registerHexMesh(
      "my mesh",
      gridPoints,gridHexes);

  Eigen::VectorXd f(gridPoints.rows());
    for (Eigen::Index iV = 0; iV < gridPoints.rows(); iV++) {
      //f(iV) = std::fmod(calculateSolidAngle( gridPoints.row(iV), tri0) + calculateSolidAngle( gridPoints.row(iV), tri1),(4*M_PI));
      f(iV) = std::fmod(calculateSolidAngle( gridPoints.row(iV), curvePoints),(4*M_PI));
    }

  scalarQ = psMesh->addVertexScalarQuantity("scalar Q", f);
  scalarQ->setEnabledLevelSet(true);
  scalarQ->setLevelSetValue(.3f);

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
