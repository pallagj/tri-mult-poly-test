#ifndef GORBEALGORITMUSOK_H
#define GORBEALGORITMUSOK_H
/******************************************************************************







                    Más hasznos görbe-algoritmusok

                        a "NURBS Book" alapján


                    Dr. Várady Tamás, Salvi Péter
               BME Villamosmérnöki és Informatikai Kar
               Irányítástechnika és Informatika Tanszék






******************************************************************************/


// [1] Csomóbeszúrás

// k: span, s: a csomó eredeti multiplicitása, r: beszúrások száma
// feltételezi, hogy s + r <= p
BSplineCurve BSplineCurve::
insertKnot(double u, size_t k, size_t s, size_t r) const
{
  BSplineCurve result;
  result.p = p; result.n = n + r;

  result.knots.reserve(knots.size() + r);
  std::copy_n(knots.begin(), k + 1, std::back_inserter(result.knots));
  std::fill_n(std::back_inserter(result.knots), r, u);
  std::copy(knots.begin() + k + 1, knots.end(),
            std::back_inserter(result.knots));

  result.cp.resize(cp.size() + r);
  std::copy_n(cp.begin(), k - p + 1, result.cp.begin());
  std::copy(cp.begin() + k - s, cp.end(), result.cp.begin() + r + k - s);

  // ↓


  // ↑

  PointVector tmp; tmp.reserve(p - s + 1);

  std::copy_n(cp.begin() + k - p, p - s + 1, std::back_inserter(tmp));

  size_t L = k - p + 1;         // result.cp következő indexe
  for (size_t j = 1; j <= r; ++j, ++L) {
    for (size_t i = 0; i <= p - j - s; ++i) {
      double alpha = (u - knots[L+i]) / (knots[i+k+1] - knots[L+i]);
      tmp[i] = tmp[i+1] * alpha + tmp[i] * (1.0 - alpha);
    }
    result.cp[L] = tmp[0];
    result.cp[k+r-j-s] = tmp[p-j-s];
  }
  std::copy_n(tmp.begin() + 1, p - s - 1 - r, result.cp.begin() + L);

  return result;
}


// [2] Többszörös csomóbeszúrás

// pl. new_knots = {2,2,4,4,4,5}
BSplineCurve BSplineCurve::refineKnots(DoubleVector new_knots) const
{
  size_t r = new_knots.size();
  size_t a = findSpan(new_knots[0]);
  size_t b = findSpan(new_knots[r-1]) + 1;

  BSplineCurve result;
  result.p = p; result.n = n + r;
  result.knots.resize(knots.size() + r);
  result.cp.resize(cp.size() + r);

  std::copy_n(knots.begin(), a + 1, result.knots.begin());
  std::copy(knots.begin() + b + p, knots.end(),
            result.knots.begin() + b + p + r);
  std::copy_n(cp.begin(), a - p + 1, result.cp.begin());
  std::copy(cp.begin() + b - 1, cp.end(), result.cp.begin() + b - 1 + r);

  // ↓

  // ↑
  size_t i = b+p-1, k = b+p+r, j = r; // i: knots, k: result.knots, j: new_knots
  do { --j; --k;
    for (; new_knots[j] <= knots[i] && i > a; --i, --k) {
      result.knots[k] = knots[i];
      result.cp[k-p-1] = cp[i-p-1];
    }
    result.cp[k-p-1] = result.cp[k-p];
    for (size_t l = 1; l <= p; ++l) {
      size_t index = k - p + l;
      double alpha = result.knots[k+l] - new_knots[j];
      if (fabs(alpha) == 0.0)
        result.cp[index-1] = result.cp[index];
      else {
        alpha /= result.knots[k+l] - knots[i-p+l];
        result.cp[index-1] = result.cp[index-1] * alpha +
                             result.cp[index] * (1.0 - alpha);
      }
    }
    result.knots[k] = new_knots[j];
  } while (j > 0);
  return result;
}

// [3] Projekció

Point BSplineCurve::
projectPoint(const Point &point, double &u, double &distance,
             size_t resolution, double distance_tol, double cosine_tol) const
{
  // Ha a pont a görbén van, elég volna a pont körüli konvex burkot nézni...
  double span_min = knots[p], span_max = knots[n+1];
  distance = std::numeric_limits<double>::max();
  for (size_t i = 0; i < resolution; ++i) {
    double v = span_min +
      (span_max - span_min) * (double)i / (double)(resolution - 1);
    double d = (evaluate(v) - point).norm();
    if (d < distance) {
      distance = d;
      u = v;
    }
  }

  VectorVector der;
  derivativesByControlPoints(u, 2, der);
  Vector deviation = der[0] - point;
  // ↓

  // ↑

  while (distance > distance_tol) {
    double scaled_error = der[1] * deviation;
    double cosine_err = fabs(scaled_error) / (der[1].norm() * distance);
    if (cosine_err < cosine_tol)
      break;

    double old = u;
    u -= scaled_error / (der[2] * deviation + der[1] * der[1]);
    u = std::min(std::max(u, span_min), span_max);

    if ((der[1] * (u - old)).norm() < distance_tol)
      break;

    derivativesByControlPoints(u, 2, der);
    deviation = der[0] - point;
    distance = deviation.norm();
  }

  return der[0];
}


// [4] Interpoláció

BSplineCurve BSplineCurve::interpolate(PointVector points, size_t p)
{
  size_t n = points.size() - 1;
  size_t m = n + p + 2;
  DoubleVector u; u.reserve(n + 1);
  u.push_back(0.0);
  double total = 0.0;
  for (size_t i = 1; i <= n; ++i) {
    u.push_back(std::sqrt((points[i] - points[i-1]).norm()));
    total += u.back();
  }
  for (size_t i = 1; i < n; ++i)
    u[i] = u[i-1] + u[i] / total;
  u.back() = 1.0;

  BSplineCurve result;
  result.p = p; result.n = n;
  result.knots.reserve(m + 1);
  result.cp.reserve(n + 1);

  // ↓

  // ↑

  std::fill_n(std::back_inserter(result.knots), p + 1, 0.0);
  for (size_t j = p + 1; j <= n; ++j) {
    double t = 0.0;
    for (size_t i = j - p; i <= j - 1; ++i)
      t += u[i];
    result.knots.push_back(t / p);
  }
  std::fill_n(std::back_inserter(result.knots), p + 1, 1.0);

  // ↓

  // ↑

  Eigen::MatrixXd A(n + 1, n + 1); A.setZero();
  for (size_t i = 0; i <= n; ++i) {
    size_t span = result.findSpan(u[i]);
    DoubleVector coeff; result.basisFunctions(span, u[i], coeff);
    for (size_t j = 0; j <= p; ++j)
      A(i, span - p + j) = coeff[j];
  }
  Eigen::MatrixXd b(n + 1, 3);
  for (size_t i = 0; i <= n; ++i) {
    b(i, 0) = points[i].x;
    b(i, 1) = points[i].y;
    b(i, 2) = points[i].z;
  }
  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  for (size_t i = 0; i <= n; ++i)
    result.cp.emplace_back(x(i, 0), x(i, 1), x(i, 2));

  return result;
}


// [5] Approximáció

BSplineCurve BSplineCurve::approximate(PointVector points, size_t p, size_t n)
{
  size_t m = points.size() - 1;
  assert(n <= m);

  BSplineCurve result;
  result.p = p; result.n = n;
  result.knots.reserve(n + p + 2);
  result.cp.reserve(n + 1);

  DoubleVector u; u.reserve(m + 1);
  u.push_back(0.0);
  double total = 0.0;
  for (size_t i = 1; i <= m; ++i) {
    u.push_back(std::sqrt((points[i] - points[i-1]).norm()));
    total += u.back();
  }
  for (size_t i = 1; i < m; ++i)
    u[i] = u[i-1] + u[i] / total;
  u.back() = 1.0;
  // ↓

  // ↑
  std::fill_n(std::back_inserter(result.knots), p + 1, 0.0);
  for (size_t j = 1; j <= n - p; ++j) {
    double d = (double)(m + 1) / (n - p + 1);
    size_t i = d * j;
    double alpha = d * j - i;
    double knot = (1.0 - alpha) * u[i-1] + alpha * u[i];
    result.knots.push_back(knot);
  }
  std::fill_n(std::back_inserter(result.knots), p + 1, 1.0);

  Eigen::MatrixXd N(m + 1, n + 1); N.setZero();
  Eigen::MatrixXd S(m + 1, 3);
  for (size_t i = 0; i <= m; ++i) {
    size_t span = result.findSpan(u[i]);
    DoubleVector coeff; result.basisFunctions(span, u[i], coeff);
    for (size_t j = 0; j <= p; ++j)
      N(i, span - p + j) = coeff[j];
    S(i, 0) = points[i].x;
    S(i, 1) = points[i].y;
    S(i, 2) = points[i].z;
  }
  // ↓

  // ↑

  Eigen::MatrixXd x = N.fullPivLu().solve(S);
  for (size_t i = 0; i <= n; ++i)
    result.cp.emplace_back(x(i, 0), x(i, 1), x(i, 2));

  return result;
}


// Vége

#endif // GORBEALGORITMUSOK_H
