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
    //particles = new BalloonParticle[vertices.size()];
    
    for (int x = 0; x < vertices.size(); x++)
    {
        BalloonParticle p;
        
        p.setOriginalPosition(vertices[x]);
        p.setPosition(vertices[x]);
        p.setVelocity(Vec3f(0,0,0));
        p.setMass(1.0);
        //p.balloon = this;
        p.particle_id = x;
        
        particles.push_back(p);
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

    
    position -= Vec3f(0, 0.0001, 0);
}

void Sphere::PackSphereSurface(float* &current) {
    // like the last assignment...  to make wireframe edges...
    //
    //   a-----------------------b
    //   |\                     /|
    //   |  \                 /  |
    //   |    \             /    |
    //   |      \         /      |
    //   |        \     /        |
    //   |          \ /          |
    //   |           x           |
    //   |          / \          |
    //   |        /     \        |
    //   |      /         \      |
    //   |    /             \    |
    //   |  /                 \  |
    //   |/                     \|
    //   d-----------------------c
    //
    // mesh surface positions & normals
    computeFaceNormals();
    
    
    //reset the cached normals
    for (int i = 0; i < mesh_vertices.size(); i++)
    {
        particles[i].valid_cache = false;
    }
    
    for (int i = 0; i < mesh_faces.size(); i++) {
        //for (int i = 0; i < nx-1; i++) {
        //for (int j = 0; j < ny-1; j++) {
        
        Face f = mesh_faces[i];
        const BalloonParticle &a = particles[f.v[0]];//getParticle(i,j);
        const BalloonParticle &b = particles[f.v[1]];//getParticle(i,j+1);
        const BalloonParticle &c = particles[f.v[2]];//getParticle(i+1,j+1);
        const BalloonParticle &d = particles[f.v[3]];//getParticle(i+1,j);
        
        const Vec3f &a_pos = radius * a.getPosition() + position;
        const Vec3f &b_pos = radius * b.getPosition() + position;
        const Vec3f &c_pos = radius * c.getPosition() + position;
        const Vec3f &d_pos = radius * d.getPosition() + position;
        
        Vec3f x_pos = (a_pos+b_pos+c_pos+d_pos) * 0.25f;
        
        Vec3f a_normal = computeGouraudNormal(f.v[0]);
        Vec3f b_normal = computeGouraudNormal(f.v[1]);
        Vec3f c_normal = computeGouraudNormal(f.v[2]);
        Vec3f d_normal = computeGouraudNormal(f.v[3]);
        if (!mesh_data->gouraud) {
            // compute normals at each corner and average them
            Vec3f top = b_pos-a_pos;
            Vec3f bottom = c_pos-d_pos;
            Vec3f horiz = (top+bottom); horiz.Normalize();
            Vec3f left = d_pos-a_pos;
            Vec3f right = c_pos-b_pos;
            Vec3f vert = (left+right); vert.Normalize();
            Vec3f normal;
            Vec3f::Cross3(normal,horiz,vert);
            normal.Normalize();
            a_normal = b_normal = c_normal = d_normal = normal;
        }
        Vec3f x_normal = (a_normal+b_normal+c_normal+d_normal);
        x_normal.Normalize();
        //a_normal = b_normal = c_normal = d_normal = x_normal;
        
        Vec3f ab_color = Vec3f(1,0,0);
        Vec3f bc_color = Vec3f(1,0,0);
        Vec3f cd_color = Vec3f(1,0,0);
        Vec3f da_color = Vec3f(1,0,0);
        
        Vec3f ac_color = Vec3f(1,0,0);
        Vec3f bd_color = Vec3f(1,0,0);
        
        AddWireFrameTriangle(current,
                             a_pos,b_pos,x_pos,
                             a_normal,b_normal,x_normal,
                             ab_color,bd_color,ac_color);
        AddWireFrameTriangle(current,
                             b_pos,c_pos,x_pos,
                             b_normal,c_normal,x_normal,
                             bc_color,ac_color,bd_color);
        AddWireFrameTriangle(current,
                             c_pos,d_pos,x_pos,
                             c_normal,d_normal,x_normal,
                             cd_color,bd_color,ac_color);
        AddWireFrameTriangle(current,
                             d_pos,a_pos,x_pos,
                             d_normal,a_normal,x_normal,
                             da_color,ac_color,bd_color);
    }
    //}
}

void Sphere::AddWireFrameTriangle(float* &current,
                                   const Vec3f &apos, const Vec3f &bpos, const Vec3f &cpos,
                                   const Vec3f &anormal, const Vec3f &bnormal, const Vec3f &cnormal,
                                   const Vec3f &_abcolor, const Vec3f &_bccolor, const Vec3f &_cacolor) {
    
    Vec3f white = Vec3f(1,1,1);
    Vec3f xpos = (apos+bpos+cpos) * (1/3.0f);
    Vec3f xnormal = (anormal+bnormal+cnormal);
    xnormal.Normalize();
    
    Vec3f abcolor=_abcolor;
    Vec3f bccolor=_bccolor;
    Vec3f cacolor=_cacolor;
    if (!mesh_data->wireframe) {
        abcolor=white;
        bccolor=white;
        cacolor=white;
    }
    //std::cout << apos << bpos << cpos << std::endl;
    //std::cout << anormal << bnormal << cnormal << std::endl;
    //std::cout << apos << bpos << cpos << std::endl;
    
    // Draw 3 triangles if wireframe is active
    float12 ta = { float(apos.x()),float(apos.y()),float(apos.z()),1, float(anormal.x()),float(anormal.y()),float(anormal.z()),0, float(abcolor.r()),float(abcolor.g()),float(abcolor.b()),1 };
    float12 tb = { float(bpos.x()),float(bpos.y()),float(bpos.z()),1, float(bnormal.x()),float(bnormal.y()),float(bnormal.z()),0, float(abcolor.r()),float(abcolor.g()),float(abcolor.b()),1 };
    float12 tc = { float(xpos.x()),float(xpos.y()),float(xpos.z()),1, float(xnormal.x()),float(xnormal.y()),float(xnormal.z()),0, float(white.r()),float(white.g()),float(white.b()),1 };
    memcpy(current, &ta, sizeof(float)*12); current += 12;
    memcpy(current, &tb, sizeof(float)*12); current += 12;
    memcpy(current, &tc, sizeof(float)*12); current += 12;
    
    ta = { float(bpos.x()),float(bpos.y()),float(bpos.z()),1, float(bnormal.x()),float(bnormal.y()),float(bnormal.z()),0, float(bccolor.r()),float(bccolor.g()),float(bccolor.b()),1 };
    tb = { float(cpos.x()),float(cpos.y()),float(cpos.z()),1, float(cnormal.x()),float(cnormal.y()),float(cnormal.z()),0, float(bccolor.r()),float(bccolor.g()),float(bccolor.b()),1 };
    tc = { float(xpos.x()),float(xpos.y()),float(xpos.z()),1, float(xnormal.x()),float(xnormal.y()),float(xnormal.z()),0, float(white.r()),float(white.g()),float(white.b()),1 };
    memcpy(current, &ta, sizeof(float)*12); current += 12;
    memcpy(current, &tb, sizeof(float)*12); current += 12;
    memcpy(current, &tc, sizeof(float)*12); current += 12;
    
    ta = { float(cpos.x()),float(cpos.y()),float(cpos.z()),1, float(cnormal.x()),float(cnormal.y()),float(cnormal.z()),0, float(cacolor.r()),float(cacolor.g()),float(cacolor.b()),1 };
    tb = { float(apos.x()),float(apos.y()),float(apos.z()),1, float(anormal.x()),float(anormal.y()),float(anormal.z()),0, float(abcolor.r()),float(cacolor.g()),float(cacolor.b()),1 };
    tc = { float(xpos.x()),float(xpos.y()),float(xpos.z()),1, float(xnormal.x()),float(xnormal.y()),float(xnormal.z()),0, float(white.r()),float(white.g()),float(white.b()),1 };
    memcpy(current, &ta, sizeof(float)*12); current += 12;
    memcpy(current, &tb, sizeof(float)*12); current += 12;
    memcpy(current, &tc, sizeof(float)*12); current += 12;
}

Vec3f Sphere::computeGouraudNormal(int i) {
    
    //return Vec3f(0,0,0);
    //std::cout << i << std::endl;
    //std::cout << mesh_faces.size() << std::endl;
    //if normal has already been calculated this iteration
    if (particles[i].valid_cache)
    {
        return particles[i].cached_normal;
    }
    
    Vec3f normal;
    
    for (int n: particles[i].nearest_faces)
    {
        //std::cout << n << std::endl;
        normal += mesh_faces[n].normal;
    }
    
    normal.Normalize();
    
    //cache the normal
    particles[i].valid_cache = true;
    particles[i].cached_normal = normal;
    
    return normal;
}

void Sphere::fastCollide(std::vector<std::pair<int, float>>& collision_ids)
{
    
}

void Sphere::slowCollide(std::vector<std::pair<int, float>>& collision_ids)
{
    for (int x = 0; x < balloon->mesh_vertices.size(); x++)
    {
        float dist = (balloon->particles[x].position - position).Length();
        if (dist < radius)
        {
            collision_ids.push_back(std::make_pair(x, dist));
        }
    }
}

void Sphere::collide(std::vector<std::pair<int, float>>& collision_ids)
{
    slowCollide(collision_ids);
}

