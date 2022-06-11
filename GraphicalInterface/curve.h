#ifndef CURVE_H
#define CURVE_H
#include "vector.h"
#include <vector>

class Curve
{
public:
    std::vector<std::vector<Vector>> polygons;

    Curve(size_t numberOfPolygons = 6)
        : polygons{std::vector<std::vector<Vector>>(numberOfPolygons)} {
    }

    void setNumberOfPolygons(size_t n){
        polygons.resize(n);
    }

    void clean(){
        polygons.clear();
    }

    void addPoint(size_t index, double x, double y, double z) {
        polygons[index].push_back(Vector{x, y, z});
    }
};

#endif // CURVE_H
