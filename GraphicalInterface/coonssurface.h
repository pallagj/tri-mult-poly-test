#ifndef COONS_H
#define COONS_H

#include "bsplinecurve.h"
#include "memory"

class CoonsSurface
{
public:
    std::vector<BSplineCurve> curves;

    std::vector<std::shared_ptr<Point>> controlPoints;
    std::vector<size_t> controlPointsSizes;

    CoonsSurface(
            std::vector<std::shared_ptr<Point>> controlPoints = std::vector<std::shared_ptr<Point>>{},
            std::vector<size_t> controlPointsSizes = std::vector<size_t>{},
            std::vector<BSplineCurve> curves = std::vector<BSplineCurve>{})
        : controlPoints{controlPoints},
          controlPointsSizes{controlPointsSizes},
          curves{curves}
    {

    }

    Point evaluate(double u, double v) const;
    double calculateMean(double u, double v);

    Vector derivateU(double u, double v);
    Vector derivateV(double u, double v);

    Vector derivateUU(double u, double v);
    Vector derivateVV(double u, double v);
    Vector derivateUV(double u, double v);

    Vector getNormal(double u, double v);

    static CoonsSurface createCoonsSurface3(std::vector<Point> controlPoints, std::vector<size_t> controlPointsSizes);
};

#endif // COONS_H
