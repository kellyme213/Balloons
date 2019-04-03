// ================================================================
// Parse the command line arguments and the input file
// ================================================================


#include "spring.h"
#include "vectors.h"
#include "cloth.h"

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
