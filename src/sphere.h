#ifndef _SPHERE_H_
#define _SPHERE_H_

#include <vector>

#include "argparser.h"
#include "boundingbox.h"
#include "vectors.h"
#include "spring.h"
#include "balloon.h"
#include "face.h"
#include "balloon_particle.h"

class Sphere {

public:
  Sphere(ArgParser *args);
    ~Sphere() { delete [] particles; }

  // ACCESSORS
  const BoundingBox& getBoundingBox() const { return box; }

  // PAINTING & ANIMATING
  void PackMesh();
  void PackSphereSurface(float* &current);
  void Animate();

private:

  // PRIVATE ACCESSORS
  //Vec3f computeGouraudNormal(int i, int j) const;
    Vec3f computeGouraudNormal(int i) const;
  // HELPER FUNCTION
  void computeBoundingBox();
    void computeFaceNormals();

  // HELPER FUNCTIONS FOR ANIMATION
  void AddWireFrameTriangle(float* &current,
                            const Vec3f &apos, const Vec3f &bpos, const Vec3f &cpos,
                            const Vec3f &anormal, const Vec3f &bnormal, const Vec3f &cnormal,
                            const Vec3f &abcolor, const Vec3f &bccolor, const Vec3f &cacolor);

  // REPRESENTATION
  ArgParser *args;
  // grid data structure
  //int nx, ny;
  BoundingBox box;
    std::vector<Face> mesh_faces;
    std::vector<Vec3f> mesh_vertices;
    Vec3f position;
    float radius;
    BalloonParticle* particles;
};

// ========================================================================

#endif
