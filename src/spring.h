// ================================================================
// Parse the command line arguments and the input file
// ================================================================

#ifndef __SPRING_H__
#define __SPRING_H__

#include <string>
#include <random>
#include "cloth.h"
#include "vectors.h"

//these might be needed in the future, but who knows?
//class MeshData;
//class Mesh;
//class Cloth;
//class BoundingBox;

// ====================================================================
// ====================================================================

struct Spring {
    
    ClothParticle* leftParticle;
    ClothParticle* rightParticle;
    float k_constant;
    
    virtual Vec3f calculateForce() = 0;
};


struct StructuralSpring: Spring {
    
    
    Vec3f calculateForce();
};

struct ShearSpring: Spring {
    
    Vec3f calculateForce();
};

struct FlexionSpring: Spring {
    
    Vec3f calculateForce();
};

struct AngularSpring: Spring {
    
    ClothParticle* middleParticle;
    Vec3f calculateForce();
};

#endif
