#include "../EZGeometry.hpp"
#include <iostream>

using namespace ezgeometry;
using namespace ezgeometry::geometry2d;

int main() {
  auto *myPose = new Pose2d(Translation2d(), Rotation2d());
  myPose->Print();
  std::cout << "Hello, World!" << std::endl;
  return 0;
}
