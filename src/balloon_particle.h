#ifndef _BALLOON_PARTICLE_H_
#define _BALLOON_PARTICLE_H_

#include <vector>

#include "argparser.h"
#include "boundingbox.h"
#include "vectors.h"
#include "spring.h"
#include "face.h"
#include "sphere.h"

// =====================================================================================
// Cloth Particles
// =====================================================================================

struct ShearSpring;
struct StructuralSpring;
struct FlexionSpring;
struct AngularSpring;


class BalloonParticle {
public:
  // ACCESSORS
  const Vec3f& getOriginalPosition() const{ return original_position; }
  const Vec3f& getPosition() const{ return position; }
  const Vec3f& getVelocity() const{ return velocity; }
  const Vec3f& getAcceleration() const { return acceleration; }
  Vec3f getForce() const { return float(mass)*acceleration; }
  double getMass() const { return mass; }
  bool isFixed() const { return fixed; }
  // MODIFIERS
  void setOriginalPosition(const Vec3f &p) { original_position = p; }
  void setPosition(const Vec3f &p) { position = p; }
  void setVelocity(const Vec3f &v) { velocity = v; }
  void setAcceleration(const Vec3f &a) { acceleration = a; }
  void setMass(double m) { mass = m; }
  void setFixed(bool b) { fixed = b; }
  // REPRESENTATION
  Vec3f original_position;
  Vec3f position;
  Vec3f velocity;
  Vec3f acceleration;
    
  //Vec3f new_position;
  //Vec3f new_velocity;
  Vec3f new_acceleration;
  double mass;
  bool fixed = false;
  
  std::vector<ShearSpring> shear_springs;
  std::vector<StructuralSpring> structural_springs;
  std::vector<FlexionSpring> flexion_springs;
  std::vector<AngularSpring> angular_springs;
  Balloon* balloon;
  int particle_id;
  std::vector<int> nearest_faces;
  std::vector<int> nearest_particles;
  bool valid_cache = false;
  Vec3f cached_normal;

};

// ========================================================================

#endif
