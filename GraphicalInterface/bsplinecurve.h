#ifndef BSPLINECURVE_H
#define BSPLINECURVE_H

#include "vector.h"
#include <algorithm>
#include <memory>

class BSplineCurve
{
public:

  size_t p;                     // fokszám
  size_t n;                     // n + 1 = cp.size()
  DoubleVector knots;           // első és utolsó p+1 érték megegyezik (*)
  std::vector<std::shared_ptr<Point>> cp;               // knots.size() = cp.size() + p + 1

  BSplineCurve(size_t p = 0, DoubleVector knots=DoubleVector{}, std::vector<std::shared_ptr<Point>> cp = std::vector<std::shared_ptr<Point>>{}) : p{p}, n{cp.size()-1}, knots{knots}, cp{cp} {}

  size_t findSpan(double u) const;
  void basisFunctions(size_t i, double u, DoubleVector &coeff) const;
  Point evaluate(double u) const;
  size_t findSpanWithMultiplicity(double u, size_t &multi) const;
  void derivativeControlPoints(size_t d, size_t r1, size_t r2, PointMatrix &dcp) const;
  void basisFunctionsAll(size_t i, double u, DoubleMatrix &coeff) const;
  Point derivativesByControlPoints(double u, size_t d, VectorVector &der) const;
  Point derivativesByControlPoints(double u, size_t d = 1) const;
  Point evaluate2DRational(double u) const;
  Point derivatives2DRational(double u, size_t d, VectorVector &der);
  size_t binomial(size_t n, size_t k);


  static BSplineCurve create3BSplineCurve(std::vector<std::shared_ptr<Point>> cp){
      DoubleVector knots{};

      for(int i = 0; i<4; i++) {
          knots.push_back(0);
      }

      int s = cp.size() - 3;
      for(int i=1; i<=s-1; i++) {
          knots.push_back(i/(double)s);
      }


      for(int i = 0; i<4; i++) {
          knots.push_back(1);
      }

      return BSplineCurve{3, knots, cp};
  }
};






#endif // BSPLINECURVE_H
