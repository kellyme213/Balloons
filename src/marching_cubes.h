#ifndef _MARCHING_CUBES_
#define _MARCHING_CUBES_

#include <vector>
#include "vectors.h"

// ==================================================================================
// A helper class for the Marching Cubes algorithm

class GridValue {
public:
  GridValue() {}
  GridValue(const Vec3f &p, const Vec3f &n, double v) : position(p), normal(n), value(v) { } 
  Vec3f position;
  Vec3f normal;
  double value;
};



class MCTri {
public:
  Vec3f a;
  Vec3f b;
  Vec3f c;
  Vec3f na;
  Vec3f nb;
  Vec3f nc;
};

// ==================================================================================
// The marching cubes algorithm is used to render an isosurface
// of a signed distance field.

class MarchingCubes {

public:
  MarchingCubes(int _nx, int _ny, int _nz, double _dx, double _dy, double _dz) : 
    nx(_nx), ny(_ny), nz(_nz), dx(_dx), dy(_dy), dz(_dz) { values = new double[nx*ny*nz]; }
  ~MarchingCubes() { delete [] values; }

  // position
  double get(int x, int y, int z) const {
    assert (x >= 0 && x < nx);
    assert (y >= 0 && y < ny);
    assert (z >= 0 && z < nz);
    return values[z*nx*ny + y*nx + x];
  }
  // normal
  Vec3f getNormal(int x, int y, int z) const {
    assert (x >= 0 && x < nx);
    assert (y >= 0 && y < ny);
    assert (z >= 0 && z < nz);
    double dx;
    double dy;
    double dz;
    if (x == 0)         dx = get(x+1,y,z)-get(x  ,y,z);
    else if (x == nx-1) dx = get(x  ,y,z)-get(x-1,y,z);
    else                dx = get(x+1,y,z)-get(x-1,y,z);
    if (y == 0)         dy = get(x,y+1,z)-get(x,y  ,z);
    else if (y == ny-1) dy = get(x,y  ,z)-get(x,y-1,z);
    else                dy = get(x,y+1,z)-get(x,y-1,z);
    if (z == 0)         dz = get(x,y,z+1)-get(x,y,z  );
    else if (z == nz-1) dz = get(x,y,z  )-get(x,y,z-1);
    else                dz = get(x,y,z+1)-get(x,y,z-1);
    Vec3f norm(-dx,-dy,-dz);
    norm.Normalize();
    return norm;
  }
  void set(int x, int y, int z, double v) const {
    assert (x >= 0 && x < nx);
    assert (y >= 0 && y < ny);
    assert (z >= 0 && z < nz);
    values[z*nx*ny + y*nx + x] = v;
  }

  // =============
  // THE DRAW CODE
  void GenerateTris();
  int numMCTris() const { return mc_tris.size(); }
  const std::vector<MCTri>& getMCTris() const { return mc_tris; }
  
private:

  // private helper functions
  void PaintTetra(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, double isosurface); 
  void PaintTetraHelper4(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, double isosurface); 
  void PaintTetraHelper3(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, double isosurface);
  void PaintTetraHelper2(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, double isosurface); 
  void PaintTetraHelper1(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, double isosurface); 
  void drawIfBoundary(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3);
  void drawTriangleWithNormals(const Vec3f &n1, const Vec3f &p1, 
			       const Vec3f &n2, const Vec3f &p2, 
			       const Vec3f &n3, const Vec3f &p3);

  // ==============
  // REPRESENTATION
  int nx, ny, nz;
  double dx, dy, dz;
  double* values;

  std::vector<MCTri> mc_tris;
};

// ==================================================================================

#endif

