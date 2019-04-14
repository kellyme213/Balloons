#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <set>
#include <cstdlib>
#include <algorithm>
#include "balloon.h"
#include "argparser.h"
#include "utils.h"
#include "meshdata.h"

extern MeshData *mesh_data;

// ================================================================================
// ================================================================================

struct Face;

std::vector<Face> getFacesWithVertex(int n, std::vector<Face>& faces)
{
    std::vector<Face> ret;
    for (int x = 0; x < faces.size(); x++)
    {
        if (faces[x].contains(n))
        {
            ret.push_back(faces[x]);
        }
    }
    
    return ret;
}

#define MAX_CHAR_PER_LINE 200

Balloon::Balloon(ArgParser *_args) {
  args =_args;

  // open the file
  std::ifstream istr(std::string(args->path+'/'+args->balloon_file).c_str());
    
  assert (istr.good());
  std::string token;
    char line[MAX_CHAR_PER_LINE];

    std::vector<Face> faces;
    std::vector<Vec3f> vertices;
    while (istr.getline(line,MAX_CHAR_PER_LINE))
    {
        std::stringstream ss;
        ss << line;
        
        // check for blank line
        token = "";
        ss >> token;
        if (token == "") continue;
        if (token == "#")
        {
            continue;
        }
        
        if (token == "v")
        {
            float x;
            float y;
            float z;
            ss >> x >> y >> z;
            
            vertices.push_back(Vec3f(x, y, z));
        }
        else if (token == "f")
        {
            Face f;
            
            ss >> f.v[0] >> f.v[1] >> f.v[2] >> f.v[3];
            f.v[0] -= 1;
            f.v[1] -= 1;
            f.v[2] -= 1;
            f.v[3] -= 1;

            faces.push_back(f);
        }
    }
    
    
    particles = new BalloonParticle[vertices.size()];
    
    for (int x = 0; x < vertices.size(); x++)
    {
        BalloonParticle& p = particles[x];
        
        p.setOriginalPosition(vertices[x]);
        p.setPosition(vertices[x]);
        p.setVelocity(Vec3f(0,0,0));
        p.setMass(1.0);
    }
    

    for (int x = 0; x < vertices.size(); x++)
    {
        std::vector<Face> collectedFaces = getFacesWithVertex(x, faces);
        
        std::set<int> uniqueVerts;
        for (int y = 0; y < collectedFaces.size(); y++)
        {
            uniqueVerts.insert(collectedFaces[y].v[0]);
            uniqueVerts.insert(collectedFaces[y].v[1]);
            uniqueVerts.insert(collectedFaces[y].v[2]);
            uniqueVerts.insert(collectedFaces[y].v[3]);
        }
        uniqueVerts.erase(x);
        
        bool containsStructural = false;
        for (std::set<int>::iterator itr = uniqueVerts.begin();
             itr != uniqueVerts.end(); itr++)
        {
            int testVert = *itr;
            for (int y = 0; y < collectedFaces.size(); y++)
            {
                int connectivity = collectedFaces[y].connectivity(x, testVert);
                
                if (connectivity == 0) //shear
                {
                    ShearSpring spring;
                    spring.leftParticle = &particles[x];
                    spring.rightParticle = &particles[testVert];
                    particles[x].shear_springs.push_back(spring);
                    std::cout << x << " " << testVert << std::endl;
                }
                else if (connectivity == 1) //structural
                {
                    StructuralSpring spring;
                    spring.leftParticle = &particles[x];
                    spring.rightParticle = &particles[testVert];
                    particles[x].structural_springs.push_back(spring);
                    containsStructural = true;
                }
            }
            
            if (!containsStructural)
            {
                continue;
            }
            
            std::vector<Face> collectedFaces2 = getFacesWithVertex(testVert, faces);
            
            std::set<int> uniqueVerts2;
            for (int y = 0; y < collectedFaces2.size(); y++)
            {
                uniqueVerts2.insert(collectedFaces2[y].v[0]);
                uniqueVerts2.insert(collectedFaces2[y].v[1]);
                uniqueVerts2.insert(collectedFaces2[y].v[2]);
                uniqueVerts2.insert(collectedFaces2[y].v[3]);
            }
            uniqueVerts2.erase(testVert);
            
            //ugh
            for (std::set<int>::iterator itr2 = uniqueVerts2.begin();
                 itr2 != uniqueVerts2.end(); itr2++)
            {
                int testVert2 = *itr2;
                if (testVert2 == x)
                {
                    continue;
                }
                
                for (int y = 0; y < collectedFaces.size(); y++)
                {
                    int connectivity = collectedFaces[y].connectivity(x, testVert);
                    if (connectivity == 1)
                    {

                        for (int z = 0; z < collectedFaces2.size(); z++)
                        {
                            int connectivity2 = collectedFaces2[z].connectivity(testVert, testVert2);
                            if (connectivity2 == 1)
                            {
                                bool badSpring = false;
                                for (int n = 0; n < collectedFaces.size() && !badSpring; n++)
                                {
                                    if (collectedFaces[n].connectivity(x, testVert2) >= 0)
                                    {
                                        badSpring = true;
                                    }
                                }
                                
                                for (int n = 0; n < collectedFaces2.size() && !badSpring; n++)
                                {
                                    if (collectedFaces2[n].connectivity(x, testVert2) >= 0)
                                    {
                                        badSpring = true;
                                    }
                                }
                                if (badSpring)
                                {
                                    continue;
                                }
                                
                                AngularSpring spring1;
                                spring1.leftParticle = &particles[x];
                                spring1.rightParticle = &particles[testVert2];
                                spring1.middleParticle = &particles[testVert];
                                particles[x].angular_springs.push_back(spring1);
                                std::cout << x << " " << testVert << " " << testVert2 << std::endl;
                                
                                FlexionSpring spring2;
                                spring1.leftParticle = &particles[x];
                                spring1.rightParticle = &particles[testVert2];
                                particles[x].flexion_springs.push_back(spring2);
                            }
                        }
                    }
                }
            }
        }
    }
    
    this->mesh_faces = faces;

    std::cout << std::endl;
    for (int x = 0; x < vertices.size(); x++)
    {
        std::cout << x << std::endl;
        std::cout << particles[x].structural_springs.size() << " ";
        std::cout << particles[x].shear_springs.size() << " ";
        std::cout << particles[x].angular_springs.size() << " ";
        std::cout << particles[x].flexion_springs.size() << std::endl;

    }
/*
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
 */
}

// ================================================================================

void Balloon::computeBoundingBox() {
    /*
  box = BoundingBox(getParticle(0,0).getPosition());
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      box.Extend(getParticle(i,j).getPosition());
      box.Extend(getParticle(i,j).getOriginalPosition());
    }
  }*/
}

// ================================================================================


double Balloon::isStretched(BalloonParticle& p1, BalloonParticle& p2, double k_constant)
{
    Vec3f originalLength = p2.original_position - p1.original_position;
    Vec3f currentLength = p2.position - p1.position;
    
    return (currentLength.Length() / originalLength.Length());
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

   
}
