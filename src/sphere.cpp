#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <set>
#include <cstdlib>
#include <algorithm>
#include "argparser.h"
#include "utils.h"
#include "meshdata.h"
#include "sphere.h"
#if _WIN32
#include <windows.h>
#endif

extern MeshData *mesh_data;

// ================================================================================
// ================================================================================


void Sphere::computeFaceNormals()
{
    for (Face& f: this->mesh_faces)
    {
        f.normal = 0.5 * (computeNormal(mesh_vertices[f.v[0]],
                                        mesh_vertices[f.v[1]],
                                        mesh_vertices[f.v[2]]) +
                          computeNormal(mesh_vertices[f.v[1]],
                                        mesh_vertices[f.v[2]],
                                        mesh_vertices[f.v[3]]));
    }
}

#define MAX_CHAR_PER_LINE 200

Sphere::Sphere(ArgParser *_args) {
  args =_args;

  // open the file
  std::ifstream istr(std::string(args->path+'/'+"collision_sphere.obj").c_str());
    
  assert (istr.good());
  std::string token;
    char line[MAX_CHAR_PER_LINE];

    std::vector<Face> faces;
    std::vector<Vec3f> vertices;
    int trianglesInMesh = 0;
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
            f.v[0] = f.v[1] = f.v[2] = f.v[3] = 0;
            ss >> f.v[0] >> f.v[1] >> f.v[2] >> f.v[3];
            f.v[0] -= 1;
            f.v[1] -= 1;
            f.v[2] -= 1;
            f.v[3] -= 1;
            if (f.v[3] < 0)
            {
                trianglesInMesh++;
                f.v[3] = f.v[2];
            }

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
        //p.balloon = this;
        p.particle_id = x;
    }
    
    this->mesh_faces = faces;
    this->mesh_vertices = vertices;
    
    for (int x = 0; x < vertices.size(); x++)
    {
        BalloonParticle& p = particles[x];
        p.nearest_faces = getFacesIDWithVertex(x, faces);
        p.nearest_particles = getClosestParticles(x, p.nearest_faces, faces);
    }
    
    if (trianglesInMesh > 0)
    {
        #if _WIN32
        
            HANDLE  hConsole;
            hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
            SetConsoleTextAttribute(hConsole, 244);
            std::cout << "WARNING: " << trianglesInMesh << " TRIANGLES IN MESH." << std::endl;
            std::cout << "BAD BAD BAD!!!!!!" << std::endl;
            std::cout << "ONLY USE QUADS PLEASE." << std::endl;
            SetConsoleTextAttribute(hConsole, 15);
        
        #endif
        
        #if __APPLE__
        
            std::cout << "\033[1;4;31mWARNING: " << trianglesInMesh << " TRIANGLES IN MESH." << std::endl;
            std::cout << "BAD BAD BAD!!!!!!" << std::endl;
            std::cout << "ONLY USE QUADS PLEASE.\033[0m" << std::endl;
        
        #endif
        
            std::cout << "thank you." << std::endl << std::endl;
    }
    
    computeBoundingBox();
}

// ================================================================================

void Sphere::computeBoundingBox() {
    
    box = BoundingBox(particles[0].getPosition());
    for (int i = 0; i < mesh_vertices.size(); i++) {
        //for (int i = 0; i < nx; i++) {
        //for (int j = 0; j < ny; j++) {
        box.Extend(particles[i].getPosition());
        box.Extend(particles[i].getOriginalPosition());
    }
  //}
}

// ================================================================================


void Sphere::Animate() {


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
