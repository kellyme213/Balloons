#include <fstream>
#include <cstring>
#include "balloon.h"
#include "argparser.h"
#include "utils.h"
#include "meshdata.h"

extern MeshData *mesh_data;

// ================================================================================


void Balloon::AddWireFrameTriangle(float* &current,
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


Vec3f super_elastic_color(const BalloonParticle &a, const BalloonParticle &b, double correction) {

    return Vec3f(0,0,0);
    
  Vec3f a_o, b_o, a_c, b_c;
  a_c = a.getPosition();
  b_c = b.getPosition();
  a_o = a.getOriginalPosition();
  b_o = b.getOriginalPosition();
  double length_o,length;
  length = (a_c-b_c).Length();
  length_o = (a_o-b_o).Length();

  if (length >= (1+0.99*correction) * length_o){
    // spring is too long, make it cyan
    return Vec3f(0,1,1);
  } else if (length <= (1-0.99*correction) * length_o) {
    // spring is too short, make it yellow
    return Vec3f(1,1,0);
  } else {
    return Vec3f(0,0,0);
  }
}


void Balloon::PackMesh() {

  int new_balloon_tri_count = 0;

  if (mesh_data->surface) {
      new_balloon_tri_count += 12 * 4 * mesh_faces.size();//(nx-1) * (ny-1);
      if (spheres.size() > 0)
      {
          new_balloon_tri_count += 12 * 4 * spheres[0].mesh_faces.size() * spheres.size();
      }

  }
  if (mesh_data->velocity) {
    //new_balloon_tri_count += 12 * (nx) * (ny);
  }
  if (mesh_data->force) {
    // *********************************************************************  
    // ASSIGNMENT:
    //
    // Insert the # of forces you will visualize...
    //
    // new_balloon_tri_count += ????
    // *********************************************************************
  }
  if (mesh_data->bounding_box) {
    new_balloon_tri_count += 12 * 12;
  }
  if (mesh_data->balloonTriCount != new_balloon_tri_count) {
    delete [] mesh_data->balloonTriData;
    mesh_data->balloonTriCount = new_balloon_tri_count;
    // allocate space for the new data
    if (mesh_data->balloonTriCount == 0) {
      mesh_data->balloonTriData = 0;
    } else {
      mesh_data->balloonTriData = new float[12*3* mesh_data->balloonTriCount];
    }
  }
  
  // Loop over all of the triangles
  float* current = mesh_data->balloonTriData;

  if (mesh_data->surface) {
    PackBalloonSurface(current);
      
      for (int x = 0; x < spheres.size(); x++)
      {
          spheres[x].PackSphereSurface(current);
      }
  }
  if (mesh_data->velocity) {
    PackBalloonVelocities(current);
  }
  if (mesh_data->force) {
    PackBalloonForces(current);
  }
  if (mesh_data->bounding_box) {
    PackBoundingBox(current,getBoundingBox());
  }
}


void Balloon::PackBalloonSurface(float* &current) {
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

      const Vec3f &a_pos = a.getPosition();
      const Vec3f &b_pos = b.getPosition();
      const Vec3f &c_pos = c.getPosition();
      const Vec3f &d_pos = d.getPosition();

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

      Vec3f ab_color = super_elastic_color(a,b,provot_structural_correction);
      Vec3f bc_color = super_elastic_color(b,c,provot_structural_correction);
      Vec3f cd_color = super_elastic_color(c,d,provot_structural_correction);
      Vec3f da_color = super_elastic_color(d,a,provot_structural_correction);

      Vec3f ac_color = super_elastic_color(a,c,provot_shear_correction);
      Vec3f bd_color = super_elastic_color(b,d,provot_shear_correction);

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



void Balloon::PackBalloonVelocities(float* &current) {
/*
  float thickness = 0.005 * mesh_data->bb_max_dim;
  float dt = mesh_data->timestep;
      
  // velocity & force visualization
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      const BalloonParticle &p = getParticle(i,j);
      const Vec3f &pos = p.getPosition();
      const Vec3f &vel = p.getVelocity();

      addEdgeGeometry(current,
                      pos,pos+dt*100*vel,
                      Vec3f(1,0,0),Vec3f(1,1,0),thickness,thickness);
      
    }
  }*/
}



void Balloon::PackBalloonForces(float* &current) {


  // *********************************************************************  
  // ASSIGNMENT:
  //
  // Visualize the forces
  //
  // *********************************************************************    
  
}


Vec3f Balloon::computeGouraudNormal(int i) const {
    
    //return Vec3f(0,0,0);
    
    //if normal has already been calculated this iteration
    if (particles[i].valid_cache)
    {
        return particles[i].cached_normal;
    }
    
    Vec3f normal;
    
    for (int n: particles[i].nearest_faces)
    {
        normal += mesh_faces[n].normal;
    }
    
    normal.Normalize();
    
    //cache the normal
    particles[i].valid_cache = true;
    particles[i].cached_normal = normal;
    
    return normal;
    
    /*
    //Face f = mesh_faces[i];
  //assert (i >= 0 && i < nx && j >= 0 && j < ny);

  Vec3f pos = particles[f.v[0]].getPosition();
  Vec3f north = pos;
  Vec3f south = pos;
  Vec3f east = pos;
  Vec3f west = pos;
  
  north = particles[f.v[0]].getPosition();
  south = particles[f.v[1]].getPosition();
  east = particles[f.v[2]].getPosition();
  west = particles[f.v[3]].getPosition();

  Vec3f vns = north - south;
  Vec3f vwe = west - east;
  vns.Normalize();
  vwe.Normalize();

  // compute normals at each corner and average them
  Vec3f normal;
  Vec3f::Cross3(normal,vns,vwe);
  normal.Normalize();
  return normal;
    */
}

// ================================================================================

