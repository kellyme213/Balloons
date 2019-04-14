// ================================================================
// Parse the command line arguments and the input file
// ================================================================

#ifndef __SPRING_H__
#define __SPRING_H__

#include <string>
#include <random>
#include "balloon.h"
#include "vectors.h"


//these might be needed in the future, but who knows?
class MeshData;
class Mesh;
class Balloon;
class BalloonParticle;
//class Cloth;
class BoundingBox;

// ====================================================================
// ====================================================================

struct Spring {
    
    BalloonParticle* leftParticle;
    BalloonParticle* rightParticle;
    float k_constant;
    
    virtual Vec3f calculateForce() = 0;
    
    virtual ~Spring()
    {
        //HAHAHA WHO CARES ABOUT MEMEORY DEALLOCATION
    }
    
};


struct StructuralSpring: Spring {
    
    bool equals(StructuralSpring s)
    {
        return leftParticle == s.leftParticle && rightParticle == s.rightParticle;
    }
    
    Vec3f calculateForce();
};

struct ShearSpring: Spring {
    
    bool equals(ShearSpring& s)
    {
        return leftParticle == s.leftParticle && rightParticle == s.rightParticle;
    }
    
    Vec3f calculateForce();
};

struct FlexionSpring: Spring {
    
    bool equals(FlexionSpring& s)
    {
        return leftParticle == s.leftParticle && rightParticle == s.rightParticle;
    }
    
    Vec3f calculateForce();
};

struct AngularSpring: Spring {
    
    bool equals(AngularSpring& s)
    {
        return leftParticle == s.leftParticle && rightParticle == s.rightParticle && middleParticle == s.middleParticle;
    }
    
    BalloonParticle* middleParticle;
    Vec3f calculateForce();
};

#endif
