// ================================================================
// Parse the command line arguments and the input file
// ================================================================

#include <iostream>
#include <fstream>

#include "balloon.h"
#include "cloth.h"
#include "fluid.h"
#include "argparser.h"
#include "meshdata.h"

#if __APPLE__
#include "matrix.h"
#else
#include <glm/gtc/type_ptr.hpp>
#endif


ArgParser *GLOBAL_args;

// The command line arguments
ArgParser::ArgParser(int argc, const char *argv[], MeshData *_mesh_data) {

  // set some default values
  mesh_data = _mesh_data;
  path = ".";
  cloth = NULL;
  fluid = NULL;
    balloon = NULL;
  
  // parse the command line arguments
  for (int i = 1; i < argc; i++) {
    if (argv[i] == std::string("-cloth")) {
      i++; assert (i < argc); 
      separatePathAndFile(argv[i],path,cloth_file);
    } else if (argv[i] == std::string("-fluid")) {
      i++; assert (i < argc); 	
      separatePathAndFile(argv[i],path,fluid_file);
    } else if (argv[i] == std::string("-balloon")) {
        i++; assert (i < argc);
        separatePathAndFile(argv[i],path,balloon_file);
    } else if (argv[i] == std::string("-size")) {
      i++; assert (i < argc); 
      mesh_data->width = mesh_data->height = atoi(argv[i]);
    } else if (argv[i] == std::string("-timestep")) {
      i++; assert (i < argc); 
      mesh_data->timestep = atof(argv[i]);
      assert (mesh_data->timestep > 0);
    } else {
      std::cout << "ERROR: unknown command line argument " 
                << i << ": '" << argv[i] << "'" << std::endl;
      exit(1);
    }
  }

  Load();
  GLOBAL_args = this;
  packMesh(mesh_data,cloth,fluid,balloon);
}

void ArgParser::Load() {
  delete cloth;
  delete fluid;
  if (cloth_file != "") {
    cloth = new Cloth(this);
  } else {
    cloth = NULL;
  }
  if (fluid_file != "") {
    fluid = new Fluid(this);
  } else {
    fluid = NULL;
  }
    if (balloon_file != "")
    {
        balloon = new Balloon(this);
    }
    else
    {
        balloon = NULL;
    }
}

void ArgParser::separatePathAndFile(const std::string &input, std::string &path, std::string &file) {
  // we need to separate the filename from the path
  // (we assume the vertex & fragment shaders are in the same directory)
  // first, locate the last '/' in the filename
  size_t last = std::string::npos;  
  while (1) {
    int next = input.find('/',last+1);
    if (next != (int)std::string::npos) { 
      last = next;
      continue;
    }
    next = input.find('\\',last+1);
    if (next != (int)std::string::npos) { 
      last = next;
      continue;
    }
    break;
  }
  if (last == std::string::npos) {
    // if there is no directory in the filename
    file = input;
    path = ".";
  } else {
    // separate filename & path
    file = input.substr(last+1,input.size()-last-1);
    path = input.substr(0,last);
  }
}



void packMesh(MeshData *mesh_data, Cloth *cloth, Fluid *fluid, Balloon *balloon) {
  if (cloth != NULL)
    cloth->PackMesh();

  if (fluid != NULL)
    fluid->PackMesh();

    if (balloon != NULL)
    {
        balloon->PackMesh();
    }
    
  BoundingBox bbox;
  if (cloth != NULL) {
    bbox = cloth->getBoundingBox();
  } else if (fluid != NULL) {
    bbox = fluid->getBoundingBox();
  }

  // the boundingbox center and size will be used to adjust the camera
  Vec3f center;
  bbox.getCenter(center);
  mesh_data->bb_center.data[0] = center.x();
  mesh_data->bb_center.data[1] = center.y();
  mesh_data->bb_center.data[2] = center.z();
  mesh_data->bb_max_dim = bbox.maxDim();
  mesh_data->bb_scale = 1.8 / float(bbox.maxDim());
}
