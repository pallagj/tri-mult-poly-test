#include "bsplinecurve.h"

// (*) Biztosítja a végpont-interpolációt


// [10] Az első nem 0 bázisfüggvény megkeresése

size_t BSplineCurve::findSpan(double u) const
{
  // Külön kell kezelni, ha az intervallum végén vagyunk
  if (u == knots[n+1])
    return n;
  return (std::upper_bound(knots.begin() + p + 1, knots.end(), u)
          - knots.begin()) - 1;
}


// [11] A bázisfüggvények kiszámítása

// knots = {0,..,0,1,...,1} => bernstein
void BSplineCurve::
basisFunctions(size_t i, double u, DoubleVector &coeff) const
{
  coeff.clear(); coeff.reserve(p + 1);
  coeff.push_back(1.0);
  DoubleVector left(p + 1), right(p + 1); // 0. elem nincs kihasználva
  for (size_t j = 1; j <= p; ++j) {
    left[j]  = u - knots[i+1-j];
    right[j] = knots[i+j] - u;
    double saved = 0.0;
    for (size_t r = 0; r < j; ++r) {
      double tmp = coeff[r] / (right[r+1] + left[j-r]);
      coeff[r] = saved + tmp * right[r+1];
      saved = tmp * left[j-r];
    }
    coeff.push_back(saved);
  }
}

// vö. [3]

// [12] Görbepont kiértékelés bázisfüggvényekkel

Point BSplineCurve::evaluate(double u) const
{
  double span = findSpan(u);
  DoubleVector coeff; basisFunctions(span, u, coeff);
  Point point(0.0, 0.0, 0.0);
  for (size_t i = 0; i <= p; ++i)
    point += (*cp[span - p + i]) * coeff[i];
  return point;
}

// vö. [4]


// [13] Az első nem 0 bázisfüggvény megkeresése + multiplicitás

size_t BSplineCurve::findSpanWithMultiplicity(double u, size_t &multi) const
{
  auto range = std::equal_range(knots.begin(), knots.end(), u);
  multi = range.second - range.first;

  if (u == knots[n+1])
    return n;
  return (range.second - knots.begin()) - 1;
}


// [14] deBoor algoritmus



// vö. [5]


// [15] BSpline deriváltak kontrollpontjainak kiszámítása

// Feltételezi, hogy d <= p
// Csak az [r1, r2] kontrollpont-intervallumra számolja ki
void BSplineCurve::derivativeControlPoints(size_t d, size_t r1, size_t r2,
                                           PointMatrix &dcp) const
{
  dcp.clear(); dcp.resize(d + 1);
  size_t r = r2 - r1;
  dcp[0].reserve(r + 1);
  for (size_t i = 0; i <= r; ++i)
    dcp[0].push_back(*cp[r1+i]);
  for (size_t k = 1; k <= d; ++k) {
    dcp[k].reserve(r + 1 - k);
    size_t tmp = p - k + 1;
    for (size_t i = 0; i <= r - k; ++i)
      dcp[k].push_back((dcp[k-1][i+1] - dcp[k-1][i]) * tmp /
                       (knots[r1+i+p+1] - knots[r1+i+k]));
  }
}

// vö. [6]


// [16] Az összes bázisfügvény (minden kisebb fokszámra is)

void BSplineCurve::basisFunctionsAll(size_t i, double u,
                                     DoubleMatrix &coeff) const
{
  coeff.clear(); coeff.resize(p + 1);
  coeff[0].push_back(1.0);
  DoubleVector left(p + 1), right(p + 1);
  for (size_t j = 1; j <= p; ++j) {
    coeff[j].reserve(j + 1);
    left[j]  = u - knots[i+1-j];
    right[j] = knots[i+j] - u;
    double saved = 0.0;
    for (size_t r = 0; r < j; ++r) {
      double tmp = coeff[j-1][r] / (right[r+1] + left[j-r]);
      coeff[j].push_back(saved + tmp * right[r+1]);
      saved = tmp * left[j-r];
    }
    coeff[j].push_back(saved);
  }
}


// [17] Deriváltszámítás

Point BSplineCurve::derivativesByControlPoints(double u, size_t d,
                                               VectorVector &der) const
{
  size_t du = std::min(d, p);
  der.clear();
  size_t span = findSpan(u);
  DoubleMatrix coeff; basisFunctionsAll(span, u, coeff);
  PointMatrix dcp; derivativeControlPoints(du, span - p, span, dcp);
  for (size_t k = 0; k <= du; ++k) {
    der.emplace_back(0.0, 0.0, 0.0);
    for (size_t j = 0; j <= p - k; ++j)
      der[k] += dcp[k][j] * coeff[p-k][j];
  }
  for (size_t k = p + 1; k <= d; ++k)
    der.emplace_back(0.0, 0.0, 0.0);
  return der[0];
}

Point BSplineCurve::derivativesByControlPoints(double u, size_t d) const
{
  VectorVector der;
  size_t du = std::min(d, p);
  der.clear();
  size_t span = findSpan(u);
  DoubleMatrix coeff; basisFunctionsAll(span, u, coeff);
  PointMatrix dcp; derivativeControlPoints(du, span - p, span, dcp);
  for (size_t k = 0; k <= du; ++k) {
    der.emplace_back(0.0, 0.0, 0.0);
    for (size_t j = 0; j <= p - k; ++j)
      der[k] += dcp[k][j] * coeff[p-k][j];
  }
  for (size_t k = p + 1; k <= d; ++k)
    der.emplace_back(0.0, 0.0, 0.0);
  return der[d];
}

// vö. [8]


// [18] 2D NURBS görbe kiértékelés

Point BSplineCurve::evaluate2DRational(double u) const
{
  Point p = evaluate(u);
  return Point(p.x / p.z, p.y / p.z, 1.0);
}


// [19] 2D NURBS görbe deriváltjainak kiértékelése

Point BSplineCurve::derivatives2DRational(double u, size_t d,
                                          VectorVector &der)
{
  der.clear(); der.reserve(d + 1);
  VectorVector der3d; derivativesByControlPoints(u, d, der3d);
  for (size_t k = 0; k <= d; ++k) {
    Vector v = der3d[k];
    for (size_t i = 1; i <= k; ++i)
      v = v - der[k-i] * der3d[i].z * binomial(k, i);
    der.push_back(v / der3d[0].z);
  }
  return der[0];
}


// [20] Binomiális együtthatók hatékony kiszámítása

size_t BSplineCurve::binomial(size_t n, size_t k)
{
  if (k > n)
    return 0;
  size_t result = 1;
  for (size_t d = 1; d <= k; ++d, --n)
    result = result * n / d;
  return result;
}

