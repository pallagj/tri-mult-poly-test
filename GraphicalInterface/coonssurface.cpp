#include "coonssurface.h"

Point CoonsSurface::evaluate(double u, double v) const {
    Point out;

    std::vector<double> uv = {u, v};

    for(int i=0; i<2; i++) {
        out += this->curves[i*2].evaluate(uv[i])*(1-uv[1-i]) + this->curves[i*2+1].evaluate(uv[i])*(uv[1-i]);
    }

    Point S00 = this->curves[0].evaluate(0);
    Point S01 = this->curves[1].evaluate(0);
    Point S10 = this->curves[0].evaluate(1);
    Point S11 = this->curves[1].evaluate(1);

    Point corr = S00*(1-u)*(1-v) + S01*(1-u)* v +
           S10*  u  *(1-v) + S11*  u  * v;

    return out - corr;
}

Vector CoonsSurface::derivateU(double u, double v){
    Point out;

    out += this->curves[0].derivativesByControlPoints(u, 1)*(1-v) + this->curves[1].derivativesByControlPoints(u, 1)*(v);
    out += this->curves[2].evaluate(v)*(-1) + this->curves[3].evaluate(v);

    Point S00 = this->curves[0].evaluate(0);
    Point S01 = this->curves[1].evaluate(0);
    Point S10 = this->curves[0].evaluate(1);
    Point S11 = this->curves[1].evaluate(1);

    Point corr = S00*(-1)*(1-v) + S01*(-1)* v +
                 S10*  1  *(1-v) + S11*  1  * v;

    return out - corr;
}

Vector CoonsSurface::derivateV(double u, double v){
    Point out;

    out += this->curves[0].evaluate(u)*(-1) + this->curves[1].evaluate(u)*(1);
    out += this->curves[2].derivativesByControlPoints(v)*(1-u) + this->curves[3].derivativesByControlPoints(v)*(u);

    Point S00 = this->curves[0].evaluate(0);
    Point S01 = this->curves[1].evaluate(0);
    Point S10 = this->curves[0].evaluate(1);
    Point S11 = this->curves[1].evaluate(1);

    Point corr = S00*(1-u)*(-1) + S01*(1-u)* 1 +
                 S10*  u  *(-1) + S11*  u  * 1;

    return out - corr;
}

Vector CoonsSurface::derivateUU(double u, double v){
    Point out;

    out += this->curves[0].derivativesByControlPoints(u, 2)*(1-v) + this->curves[1].derivativesByControlPoints(u, 2)*(v);
    out += this->curves[2].evaluate(v)*(0) + this->curves[3].evaluate(v)*(0);

    Point S00 = this->curves[0].evaluate(0);
    Point S01 = this->curves[1].evaluate(0);
    Point S10 = this->curves[0].evaluate(1);
    Point S11 = this->curves[1].evaluate(1);

    Point corr = S00*0;

    return out - corr;
}

Vector CoonsSurface::derivateVV(double u, double v){
    Point out;

    out += this->curves[0].evaluate(u)*(0) + this->curves[1].evaluate(u)*(0);
    out += this->curves[2].derivativesByControlPoints(v, 2)*(1-u) + this->curves[3].derivativesByControlPoints(v, 2)*(u);

    Point S00 = this->curves[0].evaluate(0);
    Point S01 = this->curves[1].evaluate(0);
    Point S10 = this->curves[0].evaluate(1);
    Point S11 = this->curves[1].evaluate(1);

    Point corr = S00*0;

    return out - corr;
}

Vector CoonsSurface::derivateUV(double u, double v){
    Point out;

    out += this->curves[0].derivativesByControlPoints(u)*(-1) + this->curves[1].derivativesByControlPoints(u)*(1);
    out += this->curves[2].derivativesByControlPoints(v)*(-1) + this->curves[3].derivativesByControlPoints(v)*(1);

    Point S00 = this->curves[0].evaluate(0);
    Point S01 = this->curves[1].evaluate(0);
    Point S10 = this->curves[0].evaluate(1);
    Point S11 = this->curves[1].evaluate(1);

    Point corr = S00*(-1)*(-1) + S01*(-1)* 1 +
                 S10*  1  *(-1) + S11*  1  * 1;

    return out - corr;
}

Vector CoonsSurface::getNormal(double u, double v) {
    Vector normal = derivateU(u, v) ^ derivateV(u, v);
    normal.normalize();
    return normal;
}


double CoonsSurface::calculateMean(double u, double v) {
    double E = derivateU(u, v) * derivateU(u, v);
    double F = derivateU(u, v) * derivateV(u, v);
    double G = derivateV(u, v) * derivateV(u, v);

    double L = getNormal(u, v) * derivateUU(u, v);
    double M = getNormal(u, v) * derivateUV(u, v);
    double N = getNormal(u, v) * derivateVV(u, v);

    return (N*E - 2*M*F + L*G) / (2*(E*G - F*F));
}

CoonsSurface CoonsSurface::createCoonsSurface3(std::vector<Point> controlPoints, std::vector<size_t> controlPointsSizes) {
    std::vector<std::shared_ptr<Point>> controlPointsPointers;

    for(Point cp : controlPoints) {
        controlPointsPointers.push_back(std::make_shared<Point>(cp));
    }

    std::vector<BSplineCurve> curves;

    int index = 0;
    for(int i=0; i < 4; i++) {
        int n = controlPointsSizes[i];
        int endIndex = index + n;

        std::vector<std::shared_ptr<Point>> cp;

        for(; index < endIndex; index++) {
            cp.push_back(controlPointsPointers[index]);
        }

        cp.push_back(controlPointsPointers[index % controlPointsPointers.size()]);

        if(i >= 2){
            std::reverse(cp.begin(), cp.end());
        }

        curves.push_back(BSplineCurve::create3BSplineCurve(cp));

    }

    BSplineCurve temp = curves[1];
    curves[1] = curves[2];
    curves[2] = curves[3];
    curves[3] = temp;


    return CoonsSurface{controlPointsPointers, controlPointsSizes, curves};
}

