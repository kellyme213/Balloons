
#include <vector>
#include <set>
#include "face.h"



std::vector<Face> idsToFaces(std::vector<int>& ids, std::vector<Face>& faces)
{
    std::vector<Face> ret;
    for (int id: ids)
    {
        ret.push_back(faces[id]);
    }
    
    return ret;
};

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
};

std::vector<int> getClosestParticles(int n, std::vector<int>& face_ids, std::vector<Face>& faces)
{
    std::set<int> verts;
    
    for (int x = 0; x < face_ids.size(); x++)
    {
        Face f = faces[face_ids[x]];
        
        verts.insert(f.v[0]);
        verts.insert(f.v[1]);
        verts.insert(f.v[2]);
        verts.insert(f.v[3]);
    }
    
    if (verts.find(n) != verts.end())
    {
        verts.erase(verts.find(n));
    }
    
    std::vector<int> ret;
    for (std::set<int>::iterator itr = verts.begin(); itr != verts.end(); itr++)
    {
        ret.push_back(*itr);
    }
    
    return ret;
};

std::vector<int> getFacesIDWithVertex(int n, std::vector<Face>& faces)
{
    std::vector<int> ret;
    for (int x = 0; x < faces.size(); x++)
    {
        if (faces[x].contains(n))
        {
            ret.push_back(x);
        }
    }
    
    return ret;
};
