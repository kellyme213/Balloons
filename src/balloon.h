#ifndef _BALLOON_H_
#define _BALLOON_H_

#include <vector>

#include "argparser.h"
#include "boundingbox.h"
#include "vectors.h"
#include "spring.h"
#include "face.h"
#include "sphere.h"
#include "balloon_particle.h"

class Sphere;

// =====================================================================================
// Cloth Particles
// =====================================================================================


// =====================================================================================
// Cloth System
// =====================================================================================

class Balloon {

public:
  Balloon(ArgParser *args);
  ~Balloon() { delete [] particles; }

  // ACCESSORS
  const BoundingBox& getBoundingBox() const { return box; }

  // PAINTING & ANIMATING
  void PackMesh();
  void PackBalloonSurface(float* &current);
  void PackBalloonVelocities(float* &current);
  void PackBalloonForces(float* &current);
  void Animate();

//private:

  // PRIVATE ACCESSORS
    /*
  const BalloonParticle& getParticle(int i, int j) const {
    assert (i >= 0 && i < nx && j >= 0 && j < ny);
    return particles[i + j*nx]; }
  BalloonParticle& getParticle(int i, int j) {
    assert (i >= 0 && i < nx && j >= 0 && j < ny);
    return particles[i + j*nx]; }
*/
  //Vec3f computeGouraudNormal(int i, int j) const;
    Vec3f computeGouraudNormal(int i) const;
  // HELPER FUNCTION
  void computeBoundingBox();
  void computeFaceNormals();
  void collisionDetection();
  Vec3f isStretched(BalloonParticle& p1, BalloonParticle& p2, double k_constant);
  void ProvotCorrection();
  void Correct(BalloonParticle& p1, BalloonParticle& p2, double constraint);
  Vec3f angularSpring(AngularSpring& spring);

  void angleCorrect(double constraint);

  // HELPER FUNCTIONS FOR ANIMATION
  void AddWireFrameTriangle(float* &current,
                            const Vec3f &apos, const Vec3f &bpos, const Vec3f &cpos,
                            const Vec3f &anormal, const Vec3f &bnormal, const Vec3f &cnormal,
                            const Vec3f &abcolor, const Vec3f &bccolor, const Vec3f &cacolor);

  // REPRESENTATION
  ArgParser *args;
  // grid data structure
  //int nx, ny;
    BalloonParticle* particles;
  BoundingBox box;
    std::vector<Face> mesh_faces;
    std::vector<Vec3f> mesh_vertices;
  // simulation parameters
  double damping;
  // spring constants
  double k_structural;
  double k_shear;
  double k_bend;
  // correction thresholds
  double provot_structural_correction;
  double provot_shear_correction;
    double provot_flexion_correction;
    double provot_angular_correction;
    double k_string;
    double string_stretch;
    double k_normal = 1.0;
    std::vector<Sphere> spheres;
    Vec3f string_pos;
    int string_id = -1;
    bool use_string = false;
    bool use_provot = false;
    BalloonParticle string_particle;
    //int num_spheres = 1;
};

// ========================================================================

#endif
