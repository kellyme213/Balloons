#include <fstream>
#include "balloon.h"
#include "argparser.h"
#include "utils.h"
#include "meshdata.h"

extern MeshData *mesh_data;

// ================================================================================
// ================================================================================

//TODO
Balloon::Balloon(ArgParser *_args) {
  args =_args;

  // open the file
  std::ifstream istr(std::string(args->path+'/'+args->cloth_file).c_str());
  assert (istr.good());
  std::string token;

  // read in the simulation parameters
  istr >> token >> k_structural; assert (token == "k_structural");  // (units == N/m)  (N = kg*m/s^2)
  istr >> token >> k_shear; assert (token == "k_shear");
  istr >> token >> k_bend; assert (token == "k_bend");
  istr >> token >> damping; assert (token == "damping");
  // NOTE: correction factor == .1, means springs shouldn't stretch more than 10%
  //       correction factor == 100, means don't do any correction
  istr >> token >> provot_structural_correction; assert (token == "provot_structural_correction");
  istr >> token >> provot_shear_correction; assert (token == "provot_shear_correction");

  // the cloth dimensions
  istr >> token >> nx >> ny; 
  assert (token == "m");
  assert (nx >= 2 && ny >= 2);

  // the corners of the cloth
  // (units == meters)
  Vec3f a,b,c,d;
  double x,y,z;
  istr >> token >> x >> y >> z; assert (token == "p");
  a.set(x,y,z);
  istr >> token >> x >> y >> z; assert (token == "p");
  b.set(x,y,z);
  istr >> token >> x >> y >> z; assert (token == "p");
  c.set(x,y,z);
  istr >> token >> x >> y >> z; assert (token == "p");
  d.set(x,y,z);
  
  // fabric weight  (units == kg/m^2)
  // denim ~300 g/m^2
  // silk ~70 g/m^2
  double fabric_weight;
  istr >> token >> fabric_weight; assert (token == "fabric_weight");
  double area = AreaOfTriangle(a,b,c) + AreaOfTriangle(a,c,d);

  // create the particles
  particles = new BalloonParticle[nx*ny];
  double mass = area*fabric_weight / double(nx*ny);
  for (int i = 0; i < nx; i++) {
    double x = i/double(nx-1);
    Vec3f ab = float(1-x)*a + float(x)*b;
    Vec3f dc = float(1-x)*d + float(x)*c;
    for (int j = 0; j < ny; j++) {
      double y = j/double(ny-1);
      BalloonParticle &p = getParticle(i,j);
      Vec3f abdc = float(1-y)*ab + float(y)*dc;
      p.setOriginalPosition(abdc);
      p.setPosition(abdc);
      p.setVelocity(Vec3f(0,0,0));
      p.setMass(mass);
      p.setFixed(false);
    }
  }

  // the fixed particles
  while (istr >> token) {
    assert (token == "f");
    int i,j;
    double x,y,z;
    istr >> i >> j >> x >> y >> z;
    BalloonParticle &p = getParticle(i,j);
    p.setPosition(Vec3f(x,y,z));
    p.setFixed(true);
  }

  computeBoundingBox();
}

// ================================================================================

void Balloon::computeBoundingBox() {
  box = BoundingBox(getParticle(0,0).getPosition());
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      box.Extend(getParticle(i,j).getPosition());
      box.Extend(getParticle(i,j).getOriginalPosition());
    }
  }
}

// ================================================================================

Vec3f Balloon::calculateForce(int i1, int j1, int i2, int j2, double k_constant)
{
    if (i1 < 0 || i1 >= nx || i2 < 0 || i2 >= nx || j1 < 0 || j1 >= ny || j2 < 0 || j2 >= ny)
    {
        return Vec3f(0, 0, 0);
    }
    BalloonParticle p1 = getParticle(i1, j1);
    BalloonParticle p2 = getParticle(i2, j2);
    Vec3f originalLength = p2.original_position - p1.original_position;
    Vec3f currentLength = p2.position - p1.position;
    double dist = currentLength.Length() - originalLength.Length();
    currentLength.Normalize();
    return k_constant * dist * currentLength;
}

bool Balloon::particleExists(int i, int j)
{
    return i >= 0 && i < nx && j >= 0 && j < ny;
}

double Balloon::isStretched(BalloonParticle& p1, BalloonParticle& p2, double k_constant)
{
    Vec3f originalLength = p2.original_position - p1.original_position;
    Vec3f currentLength = p2.position - p1.position;
    
    return (currentLength.Length() / originalLength.Length());
}
static int num = 0;
void Balloon::provotCorrection(int i1, int j1, int i2, int j2, double k_constant)
{
    
    if (particleExists(i2, j2))
    {
        double stretchAmount = isStretched(getParticle(i1, j1), getParticle(i2, j2), k_constant);
        if (stretchAmount > 1.0 + k_constant)
        {
            num++;
            if (num > 10000000)
            {
            //std::cerr << stretchAmount << std::endl;
            }
            BalloonParticle p1 = getParticle(i1, j1);
            BalloonParticle p2 = getParticle(i2, j2);
            Vec3f originalLength = p2.original_position - p1.original_position;
            //std::cout << originalLength.Length();
            
            Vec3f direction = p2.position - p1.position;
            //direction.Normalize();
            direction *= originalLength.Length();
            direction *= stretchAmount - 1.0 - k_constant;
            
            //std::cout << direction.Length() << " " << stretchAmount * originalLength.Length() << std::endl;
            if (p1.fixed)
            {
                getParticle(i2, j2).position -= direction;
            }
            else if (p2.fixed)
            {
                getParticle(i1, j1).position += direction;
            }
            else
            {
                getParticle(i1, j1).position += 0.5 * direction;
                getParticle(i2, j2).position -= 0.5 * direction;
            }
            
            Vec3f p = getParticle(i2, j2).position - getParticle(i1, j1).position;
            //std::cout << p.Length() / originalLength.Length() << " " << stretchAmount << std::endl;
            
        }
        else
        {
            //std::cout << "boo" << rand() << std::endl;
        }
    }
    
    // a = total length
    // b = extra stretch
    // n = rest length
    // b = a - n
    // b + n -> b
    // * b / b + n
    // * a - n / n
    
}

static float blah = 0;
void Balloon::Animate() {


  // *********************************************************************  
  // ASSIGNMENT:
  //
  // Compute the forces on each particle, and update the state
  // (position & velocity) of each particle.
  //
  // Also, this is where you'll put the Provot correction for super-elasticity
  //
  // *********************************************************************    

    /*
    Vec3f gravity(rand() % 500 - 250,
                  rand() % 500 - 250,
                  rand() % 500 - 250);
*/
    Vec3f gravity(args->mesh_data->gravity.data[0],
                  args->mesh_data->gravity.data[1],// * 50 * cos(blah),
                  args->mesh_data->gravity.data[2]);
    
    double smallest = 10000;
    double largest = -1;
    
    blah += 0.001;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            Vec3f structuralTotal;
            Vec3f shearTotal;
            Vec3f flexTotal;
            
            structuralTotal += calculateForce(i, j, i, j + 1, k_structural);
            structuralTotal += calculateForce(i, j, i, j - 1, k_structural);
            structuralTotal += calculateForce(i, j, i + 1, j, k_structural);
            structuralTotal += calculateForce(i, j, i - 1, j, k_structural);
            
            shearTotal += calculateForce(i, j, i + 1, j + 1, k_shear);
            shearTotal += calculateForce(i, j, i - 1, j - 1, k_shear);
            shearTotal += calculateForce(i, j, i + 1, j - 1, k_shear);
            shearTotal += calculateForce(i, j, i - 1, j + 1, k_shear);

            flexTotal += calculateForce(i, j, i, j + 2, k_bend);
            flexTotal += calculateForce(i, j, i, j - 2, k_bend);
            flexTotal += calculateForce(i, j, i + 2, j, k_bend);
            flexTotal += calculateForce(i, j, i - 2, j, k_bend);

            Vec3f totalForce = structuralTotal + shearTotal + flexTotal;
            totalForce += -damping * getParticle(i, j).velocity;
            totalForce += getParticle(i, j).mass * gravity;
            
            Vec3f acceleration = (1.0 / getParticle(i, j).mass) * totalForce;
            getParticle(i, j).new_acceleration = acceleration;
            
            if (acceleration.Length() > largest)
            {
                largest = acceleration.Length();
            }
            
            if (acceleration.Length() < smallest)
            {
                smallest = acceleration.Length();
            }
            
            if (acceleration.Length() * args->mesh_data->timestep > 5.0)
            {
                args->mesh_data->timestep /= 2.0;
                std::cout << "Halving timestep to " << args->mesh_data->timestep << std::endl;
            }
        }
    }
    
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            Vec3f acceleration = getParticle(i, j).new_acceleration;
            
            Vec3f velocity = args->mesh_data->timestep * (0.5 * (acceleration + getParticle(i, j).acceleration)) + getParticle(i, j).velocity;
            Vec3f position = args->mesh_data->timestep * (0.5 * (velocity + getParticle(i, j).velocity)) + getParticle(i, j).position;
            
            getParticle(i, j).acceleration = acceleration;
            getParticle(i, j).velocity = velocity;
            if (!getParticle(i, j).fixed)
            {
                getParticle(i, j).position = position;
            }
        }
    }
    
    
    for (int n = 0; n < 4; n++)
    {
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            provotCorrection(i, j, i, j + 1, provot_structural_correction);
            //provotCorrection(i, j, i, j - 1, provot_structural_correction);
            provotCorrection(i, j, i + 1, j, provot_structural_correction);
            //provotCorrection(i, j, i - 1, j, provot_structural_correction);
            
            provotCorrection(i, j, i + 1, j + 1, provot_shear_correction);
            //provotCorrection(i, j, i - 1, j - 1, provot_shear_correction);
            provotCorrection(i, j, i + 1, j - 1, provot_shear_correction);
            //provotCorrection(i, j, i - 1, j + 1, provot_shear_correction);
        }
    }
    }

}
