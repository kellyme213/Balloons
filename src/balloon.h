#ifndef _BALLOON_H_
#define _BALLOON_H_

#include <vector>

#include "argparser.h"
#include "boundingbox.h"
#include "vectors.h"
#include "spring.h"

// =====================================================================================
// Cloth Particles
// =====================================================================================

struct ShearSpring;
struct StructuralSpring;
struct FlexionSpring;
struct AngularSpring;



struct Face
{
    int v[4];
    
    bool contains(int n)
    {
        return (n == v[0]) || (n == v[1]) || (n == v[2]) || (n == v[3]);
    }
    
    int connectivity(int n1, int n2)
    {
        if (!contains(n1) || !contains(n2))
        {
            return -1;
        }
        
        if (n1 == n2)
        {
            return 2;
        }
        
        int a, b;
        
        for (int x = 0; x < 4; x++)
        {
            if (n1 == v[x])
            {
                a = x;
            }
            if (n2 == v[x])
            {
                b = x;
            }
        }
        
        //1 = structural spring
        //0 = shear spring
        return (abs(a - b) % 2);
    }
    
    bool equals(Face& f)
    {
        return (v[0] == f.v[0]) && (v[1] == f.v[1]) &&
        (v[2] == f.v[2]) && (v[3] == f.v[3]);
    }
};

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
  bool fixed;
    
    std::vector<ShearSpring> shear_springs;
    std::vector<StructuralSpring> structural_springs;
    std::vector<FlexionSpring> flexion_springs;
    std::vector<AngularSpring> angular_springs;

};

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

private:

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

  // HELPER FUNCTIONS FOR ANIMATION
  void AddWireFrameTriangle(float* &current,
                            const Vec3f &apos, const Vec3f &bpos, const Vec3f &cpos,
                            const Vec3f &anormal, const Vec3f &bnormal, const Vec3f &cnormal,
                            const Vec3f &abcolor, const Vec3f &bccolor, const Vec3f &cacolor);

  // REPRESENTATION
  ArgParser *args;
  // grid data structure
  //int nx, ny;
  BalloonParticle *particles;
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
    
    double isStretched(BalloonParticle& p1, BalloonParticle& p2, double k_constant);
};

// ========================================================================

#endif
