#include <fstream>
#include <iomanip>
#include <algorithm>
#include <iostream>

#include "fluid.h"
#include "argparser.h"
#include "boundingbox.h"
#include "marching_cubes.h"
#include "utils.h"
#include "meshdata.h"

extern MeshData *mesh_data;

#define BETA_0 1.7
#define EPSILON 0.0001

// ==============================================================
// ==============================================================
// CONSTRUCTOR
// ==============================================================
// ==============================================================

Fluid::Fluid(ArgParser *_args) {
  args = _args;
  Load();
  marchingCubes = new MarchingCubes(nx+1,ny+1,nz+1,dx,dy,dz);
  SetEmptySurfaceFull();
}

Fluid::~Fluid() { 
  delete [] cells; 
  delete marchingCubes; 
}

// ==============================================================

void Fluid::Load() {    

  // open the file
  std::ifstream istr(std::string(args->path+'/'+args->fluid_file).c_str());
  assert (istr.good());
  std::string token, token2, token3;

  // load in the grid size & dimensions
  istr >> token >> nx >> ny >> nz;  assert (token=="grid");
  assert (nx > 0 && ny > 0 && nz > 0);
  istr >> token >> dx >> dy >> dz; assert (token=="cell_dimensions");
  cells = new Cell[(nx+2)*(ny+2)*(nz+2)];

  // simulation parameters
  istr >> token >> token2;  assert (token=="flow");
  if (token2 == "compressible") compressible = true;
  else { assert (token2 == "incompressible"); compressible = false; }
  istr >> token >> token2;  assert (token=="xy_boundary");
  if (token2 == "free_slip") xy_free_slip = true;
  else { assert  (token2 == "no_slip"); xy_free_slip = false; }
  istr >> token >> token2;  assert (token=="yz_boundary");
  if (token2 == "free_slip") yz_free_slip = true;
  else { assert  (token2 == "no_slip"); yz_free_slip = false; }
  istr >> token >> token2;  assert (token=="zx_boundary");
  if (token2 == "free_slip") zx_free_slip = true;
  else { assert  (token2 == "no_slip"); zx_free_slip = false; }
  istr >> token >> viscosity;  assert (token=="viscosity");
  double gravity;
  istr >> token >> gravity;  assert (token=="gravity");
  mesh_data->gravity.data[0] = 0;
  mesh_data->gravity.data[1] = float(-9.8 * gravity);
  mesh_data->gravity.data[2] = 0;
  
  // initialize marker particles 
  istr >> token >> token2 >> token3;  assert (token=="initial_particles");
  istr >> token >> density;  assert (token=="density");
  GenerateParticles(token2,token3);

  // initialize velocities
  istr >> token >> token2;  assert (token=="initial_velocity");
  if (token2 == "zero") {
    // default is zero
  } else {
    assert (token2 == "random");
    int i,j,k;
    double max_dim = std::max(dx,std::max(dy,dz));
    for (i = -1; i <= nx; i++) {
      for (j = -1; j <= ny; j++) {
        for (k = -1; k <= nz; k++) {
          getCell(i,j,k)->set_u_plus((2*args->rand()-1)*max_dim);
	  getCell(i,j,k)->set_v_plus((2*args->rand()-1)*max_dim);
	  getCell(i,j,k)->set_w_plus((2*args->rand()-1)*max_dim);
        }
      }
    }
  }
  // read in custom velocities
  while(istr >> token) {
    int i,j,k;
    double velocity;
    assert (token == "u" || token == "v" || token == "w");
    istr >> i >> j >> k >> velocity;
    assert(i >= 0 && i < nx);
    assert(j >= 0 && j < ny);
    assert(k >= 0 && k < nz);
    if      (token == "u") getCell(i,j,k)->set_u_plus(velocity);
    else if (token == "v") getCell(i,j,k)->set_v_plus(velocity);
    else if (token == "w") getCell(i,j,k)->set_w_plus(velocity);
    else assert(0);
  }
  SetBoundaryVelocities();
}

// ==============================================================

bool Fluid::inShape(Vec3f &pos, const std::string &shape) {
  // return true if this point is inside the "shape"
  // defined procedurally (using an implicit surface)
  if (shape == "everywhere") {
    return true;
  } else if (shape == "left") {
    // a blob of particles on the lower left (for the dam)
    return (pos.x() < 0.2*nx*dx && pos.y() < 0.5*ny*dy);
  } else if (shape == "drop") {
    // a shallow pool of particles on the bottom
    double h = ny*dy/6.0;
    if (pos.y() < 2*h) return true;
    // and a sphere of particles above
    Vec3f center = Vec3f(nx*dx*0.5, 5*h,nz*dz*0.5);
    double length = (center-pos).Length();
    if (length < 0.8*h) return true;
    return false;
  } else {
    std::cout << "unknown shape: " << shape << std::endl;
    exit(0);
  }
}

// ==============================================================

void Fluid::GenerateParticles(const std::string &shape, const std::string &placement) {
  // create a set of points according to the "placement" token,
  // then check whether they are inside of the "shape"
  if (placement == "uniform") {
    int dens = (int)pow(density,0.334);
    assert (dens*dens*dens == density);
    // the uniform grid spacing
    double spacing = 1/double(dens);
    for (double x = 0.5*spacing*dx; x < nx*dx; x += spacing*dx) {
      for (double y = 0.5*spacing*dy; y < ny*dy; y += spacing*dy) {
        for (double z = 0.5*spacing*dz; z < nz*dz; z += spacing*dz) {
          Vec3f pos = Vec3f(x,y,z);
          if (inShape(pos,shape)) {
            Cell *cell = getCell(int(x/dx),int(y/dy),int(z/dz));
            FluidParticle *p = new FluidParticle();
            p->setPosition(pos);
            cell->addParticle(p);
          }
        }
      }
    }
  } else {
    assert (placement == "random");
    // note: we don't necessarily have the same number of particles in each cell
    for (int n = 0; n < nx*ny*nz*density; n++) {
      Vec3f pos = Vec3f(args->rand()*nx*dx,
                                args->rand()*ny*dy,
                                args->rand()*nz*dz);
      if (inShape(pos,shape)) {      
        Cell *cell = getCell(int(pos.x()/dx),int(pos.y()/dy),int(pos.z()/dz));
        FluidParticle *p = new FluidParticle();
        p->setPosition(pos);
        cell->addParticle(p);
      }
    }
  }
}

// ==============================================================
// ==============================================================
// ANIMATION
// ==============================================================
// ==============================================================

void Fluid::Animate() {

  // the animation manager:  this is what gets done each timestep!

  ComputeNewVelocities();
  SetBoundaryVelocities();
  
  // compressible / incompressible flow
  if (compressible == false) {
    for (int iters = 0; iters < 20; iters++) {
      double max_divergence = AdjustForIncompressibility();
      SetBoundaryVelocities();
      if (max_divergence < EPSILON) break;
    }
  }

  UpdatePressures();
  CopyVelocities();

  // advanced the particles through the fluid
  MoveParticles();
  ReassignParticles();
  SetEmptySurfaceFull();
}

// ==============================================================

void Fluid::ComputeNewVelocities() {
  double dt = mesh_data->timestep;
  int i,j,k;

  // using the formulas from Foster & Metaxas

  for (i = 0; i < nx-1; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
        double new_u_plus =
          get_u_plus(i,j,k) +            
          dt * ((1/dx) * (square(get_u_avg(i,j,k)) - square(get_u_avg(i+1,j,k))) +
                (1/dy) * (get_uv_plus(i,j-1,k) - get_uv_plus(i,j,k)) + 
                (1/dz) * (get_uw_plus(i,j,k-1) - get_uw_plus(i,j,k)) +
                mesh_data->gravity.data[0] +
                (1/dx) * (getPressure(i,j,k)-getPressure(i+1,j,k)) +
                (viscosity/square(dx)) * (get_u_plus(i+1,j  ,k  ) - 2*get_u_plus(i,j,k) + get_u_plus(i-1,j  ,k  )) +
                (viscosity/square(dy)) * (get_u_plus(i  ,j+1,k  ) - 2*get_u_plus(i,j,k) + get_u_plus(i  ,j-1,k  )) +
                (viscosity/square(dz)) * (get_u_plus(i  ,j  ,k+1) - 2*get_u_plus(i,j,k) + get_u_plus(i  ,j  ,k-1)) );
        cell->set_new_u_plus(new_u_plus);
      }
    }
  }

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny-1; j++) {
      for (k = 0; k < nz; k++) {	
        Cell *cell = getCell(i,j,k);
        double new_v_plus =
          get_v_plus(i,j,k) +
          dt * ((1/dx) * (get_uv_plus(i-1,j,k) - get_uv_plus(i,j,k)) +
                (1/dy) * (square(get_v_avg(i,j,k)) - square(get_v_avg(i,j+1,k))) +
                (1/dz) * (get_vw_plus(i,j,k-1) - get_vw_plus(i,j,k)) +
                mesh_data->gravity.data[1] +
                (1/dy) * (getPressure(i,j,k)-getPressure(i,j+1,k)) +
                (viscosity/square(dx)) * (get_v_plus(i+1,j  ,k  ) - 2*get_v_plus(i,j,k) + get_v_plus(i-1,j  ,k  )) +
                (viscosity/square(dy)) * (get_v_plus(i  ,j+1,k  ) - 2*get_v_plus(i,j,k) + get_v_plus(i  ,j-1,k  )) +
                (viscosity/square(dz)) * (get_v_plus(i  ,j  ,k+1) - 2*get_v_plus(i,j,k) + get_v_plus(i  ,j  ,k-1)) );
        cell->set_new_v_plus(new_v_plus);
      }
    }
  }

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz-1; k++) {
        Cell *cell = getCell(i,j,k);
        double new_w_plus =
          get_w_plus(i,j,k) +
          dt * ((1/dx) * (get_uw_plus(i-1,j,k) - get_uw_plus(i,j,k)) +
                (1/dy) * (get_vw_plus(i,j-1,k) - get_vw_plus(i,j,k)) +
                (1/dz) * (square(get_w_avg(i,j,k)) - square(get_w_avg(i,j,k+1))) +
                mesh_data->gravity.data[2] +
                (1/dz) * (getPressure(i,j,k)-getPressure(i,j,k+1)) +
                (viscosity/square(dx)) * (get_w_plus(i+1,j  ,k  ) - 2*get_w_plus(i,j,k) + get_w_plus(i-1,j  ,k  )) +
                (viscosity/square(dy)) * (get_w_plus(i  ,j+1,k  ) - 2*get_w_plus(i,j,k) + get_w_plus(i  ,j-1,k  )) +
                (viscosity/square(dz)) * (get_w_plus(i  ,j  ,k+1) - 2*get_w_plus(i,j,k) + get_w_plus(i  ,j  ,k-1)) );
        cell->set_new_w_plus(new_w_plus);
      }
    }
  }
}


// ==============================================================

void Fluid::SetBoundaryVelocities() {

  // zero out flow perpendicular to the boundaries (no sources or sinks)
  for (int j = -1; j <= ny; j++) {
    for (int k = -1; k <= nz; k++) {
      getCell(-1  ,j,k)->set_u_plus(0);
      getCell(nx-1,j,k)->set_u_plus(0);
      getCell(nx  ,j,k)->set_u_plus(0);
    }
  }
  for (int i = -1; i <= nx; i++) {
    for (int k = -1; k <= nz; k++) {
      getCell(i,-1  ,k)->set_v_plus(0);
      getCell(i,ny-1,k)->set_v_plus(0);
      getCell(i,ny  ,k)->set_v_plus(0);
    }
  }
  for (int i = -1; i <= nx; i++) {
    for (int j = -1; j <= ny; j++) {
      getCell(i,j,-1  )->set_w_plus(0);
      getCell(i,j,nz-1)->set_w_plus(0);
      getCell(i,j,nz  )->set_w_plus(0);
    }
  }

  // free slip or no slip boundaries (friction with boundary)
  double xy_sign = (xy_free_slip) ? 1 : -1;
  double yz_sign = (yz_free_slip) ? 1 : -1;
  double zx_sign = (zx_free_slip) ? 1 : -1;
  for (int i = 0; i < nx; i++) {
    for (int j = -1; j <= ny; j++) {
      getCell(i,j,-1)->set_u_plus(xy_sign*getCell(i,j,0)->get_u_plus());
      getCell(i,j,nz)->set_u_plus(xy_sign*getCell(i,j,nz-1)->get_u_plus());
    }
    for (int k = -1; k <= nz; k++) {
      getCell(i,-1,k)->set_u_plus(zx_sign*getCell(i,0,k)->get_u_plus());
      getCell(i,ny,k)->set_u_plus(zx_sign*getCell(i,ny-1,k)->get_u_plus());
    }
  }
  for (int j = 0; j < ny; j++) {
    for (int i = -1; i <= nx; i++) {
      getCell(i,j,-1)->set_v_plus(xy_sign*getCell(i,j,0)->get_v_plus());
      getCell(i,j,nz)->set_v_plus(xy_sign*getCell(i,j,nz-1)->get_v_plus());
    }
    for (int k = -1; k <= nz; k++) {
      getCell(-1,j,k)->set_v_plus(yz_sign*getCell(0,j,k)->get_v_plus());
      getCell(nx,j,k)->set_v_plus(yz_sign*getCell(nx-1,j,k)->get_v_plus());
    }
  }
  for (int k = 0; k < nz; k++) {
    for (int i = -1; i <= nx; i++) {
      getCell(i,-1,k)->set_w_plus(zx_sign*getCell(i,0,k)->get_w_plus());
      getCell(i,ny,k)->set_w_plus(zx_sign*getCell(i,ny-1,k)->get_w_plus());
    }
    for (int j = -1; j <= ny; j++) {
      getCell(-1,j,k)->set_w_plus(yz_sign*getCell(0,j,k)->get_w_plus());
      getCell(nx,j,k)->set_w_plus(yz_sign*getCell(nx-1,j,k)->get_w_plus());
    }
  }
}

// ==============================================================

void Fluid::EmptyVelocities(int i, int j, int k) {
  Cell *c = getCell(i,j,k);
  if (c->getStatus() != CELL_EMPTY) return;
  Cell *ciplus = getCell(i+1,j,k);
  Cell *cjplus = getCell(i,j+1,k);
  Cell *ckplus = getCell(i,j,k+1);
  if (ciplus->getStatus() == CELL_EMPTY)
    c->set_new_u_plus(0);
  if (cjplus->getStatus() == CELL_EMPTY)
    c->set_new_v_plus(0);
  if (ckplus->getStatus() == CELL_EMPTY)
    c->set_new_w_plus(0);
}


// move to new timestep
void Fluid::CopyVelocities() {
  double dt = mesh_data->timestep;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
	Cell *c = getCell(i,j,k);
	
	EmptyVelocities(i,j,k);

	c->copyVelocity();
	if (fabs(c->get_u_plus()) > 0.5*dx/dt ||
	    fabs(c->get_v_plus()) > 0.5*dy/dt ||
	    fabs(c->get_w_plus()) > 0.5*dz/dt) {
	  // velocity has exceeded reasonable threshhold
	  std::cout << "velocity has exceeded reasonable threshhold, stopping animation" << std::endl;
	  mesh_data->animate=false;
	}
      }
    }
  }
}

// ==============================================================

double Fluid::AdjustForIncompressibility() {


  // *********************************************************************  
  // ASSIGNMENT:
  //
  // This is not a complete implementation of the Marker and Cell (MAC) method.
  // Additional boundary velocities should be equalized as described in the references
  // depending on whether the boundaries are free-slip or no-slip.
  //
  // Also play around with compressible flow!
  //
  // *********************************************************************    


  // return the maximum divergence
  // (will iterate for specified # of iterations or until divergence is near zero)

  // placeholder...
  return 0;
}

// ==============================================================

void Fluid::UpdatePressures() {
  for (int i = -1; i <= nx; i++) {
    for (int j = -1; j <= ny; j++) {
      for (int k = -1; k <= nz; k++) {
	Cell *c = getCell(i,j,k);
	if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz) {
	  // compute divergence and increment/decrement pressure
	  double pressure = c->getPressure();
	  double divergence = 
	    - ( (1/dx) * (get_new_u_plus(i,j,k) - get_new_u_plus(i-1,j,k)) +
		(1/dy) * (get_new_v_plus(i,j,k) - get_new_v_plus(i,j-1,k)) +
		(1/dz) * (get_new_w_plus(i,j,k) - get_new_w_plus(i,j,k-1)) );
	  double dt = mesh_data->timestep;
	  double beta = BETA_0/((2*dt) * (1/square(dx) + 1/square(dy) + 1/square(dz)));
	  double dp = beta*divergence;
	  c->setPressure(pressure + dp);
	} else {
	  // zero out boundary cells (just in case)
	  c->setPressure(0);
	}

	
	// =======================================
	// HACK? From Foster 2001 paper?
	// zero out empty cells
	if (c->getStatus() == CELL_EMPTY) {
	  c->setPressure(0);
	}
	// ========================================

      }
    }
  }
}

// ==============================================================

void Fluid::MoveParticles() {
  double dt = mesh_data->timestep;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
	std::vector<FluidParticle*> &particles = cell->getParticles();
        for (unsigned int iter = 0; iter < particles.size(); iter++) {
          FluidParticle *p = particles[iter];
          Vec3f pos = p->getPosition();
          Vec3f vel = getInterpolatedVelocity(pos);
          Vec3f pos2 = pos + float(dt)*vel;
          // euler integration
          p->setPosition(pos2);
        }
      }
    }
  }
}

// ==============================================================

void Fluid::ReassignParticles() {
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
	std::vector<FluidParticle*> &particles = cell->getParticles();
        for (unsigned int iter = 0; iter < particles.size(); iter++) {
          FluidParticle *p = particles[iter];
          Vec3f pos = p->getPosition();
          int i2 = (int)std::min(double(nx-1),std::max(0.0,floor(pos.x()/dx)));
          int j2 = (int)std::min(double(ny-1),std::max(0.0,floor(pos.y()/dy)));
          int k2 = (int)std::min(double(nz-1),std::max(0.0,floor(pos.z()/dz)));
          // if the particle has crossed one of the cell faces 
          // assign it to the new cell
          if (i != i2 || j != j2 || k != k2) {
            cell->removeParticle(p);
            getCell(i2,j2,k2)->addParticle(p);
          } 
        }
      }
    }
  }
}

// ==============================================================

void Fluid::SetEmptySurfaceFull() {
  int i,j,k;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
        if (cell->numParticles() == 0)
          cell->setStatus(CELL_EMPTY);
        else 
          cell->setStatus(CELL_FULL);
      }
    }
  }

  // pick out the boundary cells
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
        if (cell->getStatus() == CELL_FULL &&
            (getCell(i-1,j,k)->getStatus() == CELL_EMPTY ||
             getCell(i+1,j,k)->getStatus() == CELL_EMPTY ||
             getCell(i,j-1,k)->getStatus() == CELL_EMPTY ||
             getCell(i,j+1,k)->getStatus() == CELL_EMPTY ||
             getCell(i,j,k-1)->getStatus() == CELL_EMPTY ||
             getCell(i,j,k+1)->getStatus() == CELL_EMPTY)) {
          cell->setStatus(CELL_SURFACE);
        }
      }
    }
  }
}




Vec3f Fluid::calculateWeights(const Vec3f& pos, const Vec3f& offset, int i, int j, int k) const
{
    
    //a, b, c
    //c = u * a + (1 - u) * b
    //c = u * (a - b) + b
    //u = (c - b) / (a - b)
    
    
    double a;
    double b;
    double c;
    
    c = pos.x()/dx;
    //b = floor(c) + offset[0];
    //a = ceil(c) + offset[0];
    b = (i - 1) + offset[0];
    a = i + offset[0];
    //std::cout << c << " " << i << std::endl;
    
    double u = (c - b) / (a - b);
    //std::cout << u << std::endl;
    //std::cout << u << " " << c - i + 1 - offset[0] << std::endl;
    if (u < 0 || u > 1)
    {
        //std::cout << pos;
        //std::cout << offset[0] << std::endl;
        //std::cout << a << " " << b << " " << c << std::endl;
        //std::cout << u << std::endl << std::endl;
        
    }
    //assert(u >= 0);
    //u = c - i + 1 - offset
    //c >= i
    //c - i >=0
    //1 - offset >= 0
    //u >= 0
    
    b = (j - 1) + offset[1];
    a = j + offset[1];
    c = pos.y()/dy;
    
    double v = (c - b) / (a - b);
    
    b = (k - 1) + offset[2];
    a = k + offset[2];
    c = pos.z()/dz;
    
    double w = (c - b) / (a - b);
    
    //u -= offset[0];
    //v -= offset[1];
    //w -= offset[2];
    
    
    //if (u < 0) u = 0;
    //if (u > 1) u = 1;
    //if (v < 0) v = 0;
    //if (v > 1) v = 1;
    //if (w < 0) w = 0;
    //if (w > 1) w = 1;
    
    Vec3f ret(u, v, w);

    return ret;
}

// ==============================================================

Vec3f Fluid::getInterpolatedVelocity(const Vec3f &pos) const {


  // *********************************************************************  
  // ASSIGNMENT:
  //
  // Here is the naive velocity interpolation.
  // (use a simple average of the face velocity throughout the cell)
  // Do it right, as described in the papers.
  //
  int i = int(floor(pos.x()/dx)); if (i < 0) i = 0; if (i >= nx) i = nx-1;
  int j = int(floor(pos.y()/dy)); if (j < 0) j = 0; if (j >= ny) j = ny-1;
  int k = int(floor(pos.z()/dz)); if (k < 0) k = 0; if (k >= nz) k = nz-1;
  

    
    double u, v, w;

    double xyz[2][2][2];
    double yz[2][2];
    double a[2];
    Vec3f blendVals;
    Vec3f offset;
    
    
    int roundI;
    int roundJ;
    int roundK;
    Cell* cells[2][2][2];
    //Vec3f velocities[2][2][2];
    roundI = int(floor(pos.x()/dx));
    if (roundI < 0) {roundI = 0;}
    if (roundI >= nx) {roundI = nx-1;}
    
    roundJ = int(round(pos.y()/dy));
    if (roundJ < 0) {roundJ = 0;}
    if (roundJ >= ny) {roundJ = ny-1;}
    
    roundK = int(round(pos.z()/dz));
    if (roundK < 0) {roundK = 0;}
    if (roundK >= nz) {roundK = nz-1;}
    
    for (int x = 0; x <= 1; x++)
    {
        for (int y = 0; y <= 1; y++)
        {
            for (int z = 0; z <= 1; z++)
            {
                cells[1 - x][1 - y][1 - z] = getCell(roundI - x, roundJ - y, roundK - z);
            }
        }
    }
    
    
    
    
    
    for (int x = 0; x <= 1; x++)
    {
        for (int y = 0; y <= 1; y++)
        {
            for (int z = 0; z <= 1; z++)
            {
                Cell* c = cells[x][y][z];
                xyz[x][y][z] = c->get_u_plus();
            }
        }
    }
    
    offset = Vec3f(1, 0.5, 0.5);
    blendVals = calculateWeights(pos, offset, roundI, roundJ, roundK);

    
    for (int y = 0; y <= 1; y++)
    {
        for (int z = 0; z <= 1; z++)
        {
            yz[y][z] = (1 - blendVals.x()) * xyz[0][y][z] + (blendVals.x()) * xyz[1][y][z];
        }
    }
    
    for (int z = 0; z <= 1; z++)
    {
        a[z] = (1 - blendVals.y()) * yz[0][z] + (blendVals.y()) * yz[1][z];
    }
    
    u = (1 - blendVals.z()) * a[0] + (blendVals.z()) * a[1];
    
    
    
    
    
    
    
    
    roundI = int(round(pos.x()/dx));
    if (roundI < 0) {roundI = 0;}
    if (roundI >= nx) {roundI = nx-1;}
    
    roundJ = int(floor(pos.y()/dy));
    if (roundJ < 0) {roundJ = 0;}
    if (roundJ >= ny) {roundJ = ny-1;}
    
    roundK = int(round(pos.z()/dz));
    if (roundK < 0) {roundK = 0;}
    if (roundK >= nz) {roundK = nz-1;}
    
    for (int x = 0; x <= 1; x++)
    {
        for (int y = 0; y <= 1; y++)
        {
            for (int z = 0; z <= 1; z++)
            {
                cells[1 - x][1 - y][1 - z] = getCell(roundI - x, roundJ - y, roundK - z);
            }
        }
    }
    
    
    for (int x = 0; x <= 1; x++)
    {
        for (int y = 0; y <= 1; y++)
        {
            for (int z = 0; z <= 1; z++)
            {
                Cell* c = cells[x][y][z];
                xyz[x][y][z] = c->get_v_plus();
            }
        }
    }
    
    offset = Vec3f(0.5, 1, 0.5);
    blendVals = calculateWeights(pos, offset, roundI, roundJ, roundK);

    
    for (int y = 0; y <= 1; y++)
    {
        for (int z = 0; z <= 1; z++)
        {
            yz[y][z] = (1 - blendVals.x()) * xyz[0][y][z] + (blendVals.x()) * xyz[1][y][z];
        }
    }
    
    for (int z = 0; z <= 1; z++)
    {
        a[z] = (1 - blendVals.y()) * yz[0][z] + (blendVals.y()) * yz[1][z];
    }
    
    v = (1 - blendVals.z()) * a[0] + (blendVals.z()) * a[1];
    
    
    
    
    
    
    
    
    
    
    roundI = int(round(pos.x()/dx));
    if (roundI < 0) {roundI = 0;}
    if (roundI >= nx) {roundI = nx-1;}
    
    roundJ = int(round(pos.y()/dy));
    if (roundJ < 0) {roundJ = 0;}
    if (roundJ >= ny) {roundJ = ny-1;}
    
    roundK = int(floor(pos.z()/dz));
    if (roundK < 0) {roundK = 0;}
    if (roundK >= nz) {roundK = nz-1;}
    
    for (int x = 0; x <= 1; x++)
    {
        for (int y = 0; y <= 1; y++)
        {
            for (int z = 0; z <= 1; z++)
            {
                cells[1 - x][1 - y][1 - z] = getCell(roundI - x, roundJ - y, roundK - z);
            }
        }
    }
    
    
    for (int x = 0; x <= 1; x++)
    {
        for (int y = 0; y <= 1; y++)
        {
            for (int z = 0; z <= 1; z++)
            {
                Cell* c = cells[x][y][z];
                xyz[x][y][z] = c->get_w_plus();
            }
        }
    }
    
    offset = Vec3f(0.5, 0.5, 1);
    blendVals = calculateWeights(pos, offset, roundI, roundJ, roundK);

    
    for (int y = 0; y <= 1; y++)
    {
        for (int z = 0; z <= 1; z++)
        {
            yz[y][z] = (1 - blendVals.x()) * xyz[0][y][z] + (blendVals.x()) * xyz[1][y][z];
        }
    }
    
    for (int z = 0; z <= 1; z++)
    {
        a[z] = (1 - blendVals.y()) * yz[0][z] + (blendVals.y()) * yz[1][z];
    }
    
    w = (1 - blendVals.z()) * a[0] + (blendVals.z()) * a[1];
    
    
    
    
    
    /*
    if (i == 0 || i == nx - 1)
    {
        u = get_u_avg(i,j,k);
    }
    else
    {
        //u = getU(i, j, k, pos);
    }
    
    if (j == 0 || j == ny - 1)
    {
        v = get_v_avg(i,j,k);
    }
    else
    {
        //v = getU(i, j, k, pos);
    }
    */
    //u = getU(i, j, k, pos);
    //v = getV(i, j, k, pos);

    //u = get_u_avg(i,j,k);
    //v = get_v_avg(i,j,k);
    //w = get_w_avg(i,j,k);
    
    return Vec3f(u, v, w);
  //return Vec3f(u,v,get_w_avg(i,j,k));
  //
  // *********************************************************************  

}

