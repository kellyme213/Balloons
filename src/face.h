#ifndef _FACE_H_
#define _FACE_H_

#include <vector>
#include <set>
#include "vectors.h"

struct Face
{
    int v[4];
    Vec3f normal;
    float area;
    
    bool contains(int n)
    {
        return (n == v[0]) || (n == v[1]) || (n == v[2]) || (n == v[3]);
    }
    
    int connectivity(int n1, int n2)
    {
        if (!contains(n1) || !contains(n2))
        {
            return -1;
        }
        
        if (n1 == n2)
        {
            return 2;
        }
        
        int a, b;
        
        for (int x = 0; x < 4; x++)
        {
            if (n1 == v[x])
            {
                a = x;
            }
            if (n2 == v[x])
            {
                b = x;
            }
        }
        
        //1 = structural spring
        //0 = shear spring
        return (abs(a - b) % 2);
    }
    
    bool equals(Face& f)
    {
        return (v[0] == f.v[0]) && (v[1] == f.v[1]) &&
        (v[2] == f.v[2]) && (v[3] == f.v[3]);
    }
};
std::vector<Face> idsToFaces(std::vector<int>& ids, std::vector<Face>& faces);
std::vector<Face> getFacesWithVertex(int n, std::vector<Face>& faces);
std::vector<int> getClosestParticles(int n, std::vector<int>& face_ids, std::vector<Face>& faces);
std::vector<int> getFacesIDWithVertex(int n, std::vector<Face>& faces);


#endif
