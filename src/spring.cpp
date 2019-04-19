// ================================================================
// Parse the command line arguments and the input file
// ================================================================


#include "spring.h"
#include "vectors.h"
#include "cloth.h"
#define PI 3.14159265

// ====================================================================
// ====================================================================

Vec3f StructuralSpring::calculateForce()
{
    return Vec3f();
}

Vec3f ShearSpring::calculateForce()
{
    return Vec3f();
}

Vec3f FlexionSpring::calculateForce()
{
    return Vec3f();
}

Vec3f AngularSpring::calculateForce()
{
    return Vec3f();
}

void AngularSpring::setCross(){
    Vec3f::Cross3(og_cross, leftParticle->getOriginalPosition(), rightParticle->getOriginalPosition());
}

void AngularSpring::setAngle(){
    double og_mag = og_cross.Length();
    double left_mag = leftParticle->getOriginalPosition().Length();
    double right_mag = rightParticle->getOriginalPosition().Length();

    og_angle = asin(og_mag/(left_mag*right_mag)) * (180.0/PI);
}

double AngularSpring::calculateAngle(){
    double og_mag = og_cross.Length();
    double left_mag = leftParticle->position.Length();
    double right_mag = rightParticle->position.Length();

    return asin(og_mag/(left_mag*right_mag)) * (180.0/PI);
}