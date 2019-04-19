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
#if _WIN32
#include <windows.h>
#endif

using namespace std;
extern MeshData *mesh_data;

// ================================================================================
// ================================================================================


void Balloon::computeFaceNormals()
{
    for (Face& f: this->mesh_faces)
    {
        f.normal = 0.5 * (computeNormal(particles[f.v[0]].position,
                                        particles[f.v[1]].position,
                                        particles[f.v[2]].position) +
                          computeNormal(particles[f.v[1]].position,
                                        particles[f.v[2]].position,
                                        particles[f.v[3]].position));
    }
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
        p.setMass(0.001);
        p.balloon = this;
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
    
    for (int x = 0; x < vertices.size(); x++)
    {
        std::vector<Face> collectedFaces = idsToFaces(particles[x].nearest_faces, this->mesh_faces);//getFacesWithVertex(x, faces);
        /*
        std::set<int> uniqueVerts;
        for (int y = 0; y < collectedFaces.size(); y++)
        {
            uniqueVerts.insert(collectedFaces[y].v[0]);
            uniqueVerts.insert(collectedFaces[y].v[1]);
            uniqueVerts.insert(collectedFaces[y].v[2]);
            uniqueVerts.insert(collectedFaces[y].v[3]);
        }
        uniqueVerts.erase(x);
        */
        bool containsStructural = false;
        //for (std::set<int>::iterator itr = uniqueVerts.begin();
        //     itr != uniqueVerts.end(); itr++)
        for (int testVert: particles[x].nearest_particles)
        {
            //int testVert = *itr;
            for (int y = 0; y < collectedFaces.size(); y++)
            {
                int connectivity = collectedFaces[y].connectivity(x, testVert);
                
                if (connectivity == 0) //shear
                {
                    ShearSpring spring;
                    spring.leftParticle = &particles[x];
                    spring.rightParticle = &particles[testVert];
                    bool found = false;
                    
                    for (ShearSpring s: particles[x].shear_springs)
                    {
                        if (s.equals(spring))
                        {
                            found = true;
                        }
                    }
                    
                    if (!found)
                    {
                        particles[x].shear_springs.push_back(spring);
                    }
                }
                else if (connectivity == 1) //structural
                {
                    StructuralSpring spring;
                    spring.leftParticle = &particles[x];
                    spring.rightParticle = &particles[testVert];
                    
                    bool found = false;
                    
                    for (StructuralSpring s: particles[x].structural_springs)
                    {
                        if (s.equals(spring))
                        {
                            found = true;
                        }
                    }
                    
                    if (!found)
                    {
                        particles[x].structural_springs.push_back(spring);
                        //std::cout << x << " " << testVert << std::endl;
                    }
                    
                    containsStructural = true;
                }
            }
            
            if (!containsStructural)
            {
                continue;
            }
            
            std::vector<Face> collectedFaces2 = idsToFaces(particles[testVert].nearest_faces, this->mesh_faces);//getFacesWithVertex(testVert, faces);
            
            /*std::set<int> uniqueVerts2;
            for (int y = 0; y < collectedFaces2.size(); y++)
            {
                uniqueVerts2.insert(collectedFaces2[y].v[0]);
                uniqueVerts2.insert(collectedFaces2[y].v[1]);
                uniqueVerts2.insert(collectedFaces2[y].v[2]);
                uniqueVerts2.insert(collectedFaces2[y].v[3]);
            }
            uniqueVerts2.erase(testVert);
            */
            //ugh
            //for (std::set<int>::iterator itr2 = uniqueVerts2.begin();
            //     itr2 != uniqueVerts2.end(); itr2++)
            for (int testVert2: particles[testVert].nearest_particles)
            {
                //int testVert2 = *itr2;
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
                                spring1.setCross();
                                spring1.setAngle();

                                
                                bool found = false;
                                
                                for (AngularSpring s: particles[x].angular_springs)
                                {
                                    if (s.equals(spring1))
                                    {
                                        found = true;
                                    }
                                }
                                
                                if (!found)
                                {
                                    particles[x].angular_springs.push_back(spring1);
                                }
                                

                                FlexionSpring spring2;
                                spring2.leftParticle = &particles[x];
                                spring2.rightParticle = &particles[testVert2];
                                //particles[x].flexion_springs.push_back(spring2);
                                found = false;
                                
                                for (FlexionSpring s: particles[x].flexion_springs)
                                {
                                    if (s.equals(spring2))
                                    {
                                        found = true;
                                    }
                                }
                                
                                if (!found)
                                {
                                    particles[x].flexion_springs.push_back(spring2);

                                }
                            }
                        }
                    }
                }
            }
        }
    }
    

    int shearTotal = 0;
    int structuralTotal = 0;
    int angularTotal = 0;
    int flexionTotal = 0;
    
    
    std::cout << std::endl;
    for (int x = 0; x < vertices.size(); x++)
    {
        /*
        std::cout << "Particle " << x << std::endl;
        std::cout << particles[x].structural_springs.size() << " structural, ";
        std::cout << particles[x].shear_springs.size() << " shear, ";
        std::cout << particles[x].angular_springs.size() << " angular, ";
        std::cout << particles[x].flexion_springs.size() << " flexion, " << std::endl;
        */
        shearTotal += particles[x].shear_springs.size();
        structuralTotal += particles[x].structural_springs.size();
        angularTotal += particles[x].angular_springs.size();
        flexionTotal += particles[x].flexion_springs.size();
    }
    
    
    std::cout << std::endl;
    std::cout << "Loaded " << vertices.size() << " particles, " << this->mesh_faces.size() << " faces, " << std::endl;
    std::cout << structuralTotal << " structural, ";
    std::cout << shearTotal << " shear, ";
    std::cout << angularTotal << " angular, ";
    std::cout << flexionTotal << " flexion." << std::endl << std::endl;
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
    

    
    std::ifstream istr2(std::string(args->path+'/'+args->input_file).c_str());
    //Sphere s(_args);
    
    assert (istr2.good());
    
    while (istr2.getline(line,MAX_CHAR_PER_LINE))
    {
        std::stringstream ss;
        ss << line;
        
        // check for blank line
        token = "";
        ss >> token;
        if (token == "") continue;
        
        if (token == "k_shear")
        {
            float k;
            ss >> k;
            
            for (int x = 0; x < mesh_vertices.size(); x++)
            {
                for (ShearSpring& s: particles[x].shear_springs)
                {
                    s.k_constant = k;
                }
            }
        }
        else if (token == "k_structural")
        {
            float k;
            ss >> k;
            
            for (int x = 0; x < mesh_vertices.size(); x++)
            {
                for (StructuralSpring& s: particles[x].structural_springs)
                {
                    s.k_constant = k;
                }
            }
        }
        else if (token == "k_angular")
        {
            float k;
            ss >> k;
            
            for (int x = 0; x < mesh_vertices.size(); x++)
            {
                for (AngularSpring& s: particles[x].angular_springs)
                {
                    s.k_constant = k;
                }
            }
        }
        else if (token == "k_flexion")
        {
            float k;
            ss >> k;
            
            for (int x = 0; x < mesh_vertices.size(); x++)
            {
                for (FlexionSpring& s: particles[x].flexion_springs)
                {
                    s.k_constant = k;
                }
            }
        }
        else if (token == "string_pos")
        {
            float x, y, z;
            ss >> x >> y >> z;
            string_pos = Vec3f(x, y, z);
        }
        else if (token == "string_id")
        {
            ss >> string_id;
        }
        else if (token == "s")
        {
            float x, y, z, r, m, vx, vy, vz;
            ss >> x >> y >> z >> r >> m;
            ss >> vx >> vy >> vz;
            
            Sphere s(_args);
            s.position = Vec3f(x, y, z);
            s.velocity = Vec3f(vx, vy, vz);
            s.radius = r;
            s.mass = m;
            s.balloon = this;
            spheres.push_back(s);
        }
        else if (token == "use_string")
        {
            use_string = true;
        }
        else if (token == "damping")
        {
            ss >> this->damping;
        }
        else if (token == "f")
        {
            int n;
            ss >> n;
            particles[n].fixed = true;
        }
        else if (token == "m")
        {
            float n;
            ss >> n;
            
            for (int x = 0; x < this->mesh_vertices.size(); x++)
            {
                particles[x].mass = n;
            }
        }
        else if (token == "struct_cor")
        {
            ss >> this->provot_structural_correction;
        }
        else if (token == "shear_cor")
        {
            ss >> this->provot_shear_correction;
        }
        else if (token == "flex_cor")
        {
            ss >> this->provot_flexion_correction;
        }
        else if (token == "ang_cor")
        {
            ss >> this->provot_angular_correction;
        }
        else if (token == "use_provot")
        {
            _args->mesh_data->use_provot = true;
        }
        else if (token == "k_string")
        {
            ss >> this->k_string;
        }
        else if (token == "string_stretch")
        {
            ss >> this->string_stretch;
        }
        else if (token == "k_normal")
        {
            ss >> _args->mesh_data->k_normal;
        }
        
    }
    
    if (string_id == -1 && use_string)
    {
        float distance = (particles[0].getOriginalPosition() - string_pos).Length();
        
        for (int x = 0; x < mesh_vertices.size(); x++)
        {
            float new_distance = (particles[x].getOriginalPosition() - string_pos).Length();
            if (distance > new_distance)
            {
                distance = new_distance;
                string_id = x;
            }
        }
    }
    
    string_particle.setOriginalPosition(string_pos);
    string_particle.setPosition(string_pos);

}

// ================================================================================

void Balloon::computeBoundingBox() {
    
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

void Balloon::collisionDetection()
{
    for (int x = 0; x < spheres.size(); x++)
    {
        spheres[x].Animate();
    }
}

void Balloon::Correct(BalloonParticle& p1, BalloonParticle& p2, double constraint){
      //p0 is the particle we are "looking at"
    //p1 is the particle connected to it
    Vec3f p0pos = p1.getPosition();
    Vec3f p0orig = p1.getOriginalPosition();
    Vec3f p1pos = p2.getPosition();
    Vec3f p1orig = p2.getOriginalPosition();

    double stretch = (p0pos - p1pos).Length();
    double rest = (p0orig - p1orig).Length();
    double ratio = (stretch/rest) - 1.0;

    if(ratio > constraint){
        Vec3f diff = p0pos - p1pos;
        rest *= 1.0 + constraint;
        double half = ((diff.Length() - rest)/diff.Length());
        
        if(p1.isFixed() == false){
        p1.setPosition(p0pos - (diff * (half * 0.5)));
        }else{
        p2.setPosition(p1pos + (diff * half));
        }
        if(p2.isFixed() == false){
        p2.setPosition(p1pos + (diff * (half * 0.5)));
        }else{
        p1.setPosition(p0pos - (diff * half));
        }
    }
}

void Balloon::ProvotCorrection(){
    int iterations = 6;
    for(int its = 0; its < iterations; its++){  
        for(int i = 0; i < mesh_vertices.size(); i++){
            Vec3f springforces(0.0, 0.0, 0.0);
            for(int j = 0; j < particles[i].shear_springs.size(); j++){
                Correct(*particles[j].shear_springs[j].leftParticle, *particles[j].shear_springs[j].rightParticle, provot_shear_correction);
            }
            for(int k = 0; k < particles[i].structural_springs.size(); k++){
                Correct(*particles[k].structural_springs[k].leftParticle, *particles[k].structural_springs[k].rightParticle, provot_structural_correction);
            }
            
            for(int l = 0; l < particles[i].flexion_springs.size(); l++){
                Correct(*particles[l].flexion_springs[l].leftParticle, *particles[l].flexion_springs[l].rightParticle, provot_flexion_correction);
            }
        }
    }
}

Vec3f Balloon::isStretched(BalloonParticle& p1, BalloonParticle& p2, double k_constant)
{
    //p1 is the particle we are "looking at"
    //p2 is the particle connected to it
    
    if (k_constant <= 0.001)
    {
        return Vec3f(0,0,0);
    }
    
    Vec3f originalLength = p1.original_position - p2.original_position;
    Vec3f currentLength = p1.position - p2.position;

    Vec3f p0pos = p1.position;
    Vec3f p0orig = p1.original_position;
    Vec3f p1pos = p2.position;
    Vec3f p1orig = p2.original_position;

    double stretch = (p0pos - p1pos).Length();
    double rest = (p0orig - p1orig).Length();

    Vec3f fvec = ((p0pos - p1pos) *(1/stretch))*(stretch-rest);
    return fvec * k_constant;
}

void Balloon::angleCorrect(double constraint){
    for(int i = 0; i < mesh_vertices.size(); i++){
        for(int j = 0; j < particles[i].angular_springs.size(); j++){
            double cur_angle = particles[i].angular_springs[j].calculateAngle();
            double ratio = (cur_angle/particles[i].angular_springs[i].og_angle) - 1.0;
            
            if(ratio > constraint){
                


                /*
                if(p1.isFixed() == false){
                    p1.setPosition(p0pos - (diff * (half * 0.5)));
                }else{
                    p2.setPosition(p1pos + (diff * half));
                }
                   
                if(p2.isFixed() == false){
                    p2.setPosition(p1pos + (diff * (half * 0.5)));
                }else{
                    p1.setPosition(p0pos - (diff * half));
                }*/
            }
        }
    }
}

Vec3f Balloon::angularSpring(AngularSpring& spring)
{
    
    if (spring.k_constant <= 0.001)
    {
        return Vec3f(0,0,0);
    }
    
    Vec3f lm = spring.leftParticle->position - spring.middleParticle->position;
    Vec3f rm = spring.rightParticle->position - spring.middleParticle->position;

    Vec3f lmo = spring.leftParticle->original_position - spring.middleParticle->original_position;
    Vec3f rmo = spring.rightParticle->original_position - spring.middleParticle->original_position;
    
    lm.Normalize();
    rm.Normalize();
    
    lmo.Normalize();
    rmo.Normalize();
    
    
    float angle = acos(lm.Dot3(rm));
    float o_angle = acos(lmo.Dot3(rmo));
    
    float d_angle = angle - o_angle;
    
    Vec3f dir = spring.leftParticle->position - spring.rightParticle->position;
    dir.Normalize();
    
    return spring.k_constant * d_angle * dir;
}

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
   
    //collisionDetection();
    
    this->use_provot = mesh_data->use_provot;
    this->k_normal = mesh_data->k_normal;
    float timestep = mesh_data->timestep;
    float g[3];
    g[0] = mesh_data->gravity.data[0];
    g[1] = mesh_data->gravity.data[1];
    g[2] = mesh_data->gravity.data[2];
    Vec3f gravity(g[0],g[1],g[2]);
    Vec3f helium(0.0, 200 + 5 * k_normal, 0.0);
    
    for(int i = 0; i < mesh_vertices.size(); i++){
        if(particles[i].fixed == false){
            Vec3f inflate = particles[i].cached_normal;
            inflate *= k_normal;
            BalloonParticle p = particles[i];
            Vec3f springforces(0.0, 0.0, 0.0);
            for(int j = 0; j < p.shear_springs.size(); j++){
                springforces += isStretched(*p.shear_springs[j].leftParticle, *p.shear_springs[j].rightParticle, p.shear_springs[j].k_constant);
            }
            for(int k = 0; k < p.structural_springs.size(); k++){
                springforces += isStretched(*p.structural_springs[k].leftParticle, *p.structural_springs[k].rightParticle, p.structural_springs[k].k_constant);
            }
            for(int l = 0; l < p.flexion_springs.size(); l++){
                springforces += isStretched(*p.flexion_springs[l].leftParticle, *p.flexion_springs[l].rightParticle, p.flexion_springs[l].k_constant);
            }
            /*
            for (AngularSpring& s: p.angular_springs)
            {
                springforces += angularSpring(s);
            }
            */
            if (i == string_id && use_string)
            {
                double dist = (string_pos - particles[i].position).Length();
                double orig_dist = (string_pos - particles[i].original_position).Length();
                if (string_stretch * orig_dist < dist)
                {
                    springforces += isStretched(particles[i], string_particle, k_string);
                }
            }
            Vec3f gravforces = (gravity + helium) * particles[i].getMass();
            Vec3f dampforces = damping * particles[i].getVelocity();
            Vec3f totforces = gravforces - (springforces + dampforces);
            totforces += inflate;
            Vec3f acc = totforces*(1/particles[i].getMass());
            Vec3f nvel = particles[i].getVelocity() + timestep*(acc);
            Vec3f npos  = particles[i].getPosition() + (timestep*nvel);
            particles[i].setVelocity(nvel);
            particles[i].setAcceleration(acc);
            particles[i].setPosition(npos);
        }
    }
    
    if (use_provot)
    {
        ProvotCorrection();
    }
    
    for (Sphere& s: spheres)
    {
        std::vector<std::pair<int, float>> collision_ids;
        s.collide(collision_ids);
        for (int x = 0; x < collision_ids.size(); x++)
        {
            std::pair<int, float> pair = collision_ids[x];
            float dist = pair.second - s.radius;
            int c_id = pair.first;
            
            Vec3f dir = particles[c_id].position - s.position;
            dir.Normalize();
            
            particles[c_id].position = s.radius * dir + s.position;
            
            s.force -= particles[c_id].structural_springs[0].k_constant * dir;
        }
        
        s.force += s.mass * gravity- (0.1 * damping * s.velocity);
        s.acceleration = (1.0 / s.mass) * s.force;
        s.velocity += s.acceleration * timestep;
        s.position += s.velocity * timestep;
        s.force = Vec3f(0.0, 0.0, 0.0);
    }
}
