// -*- mode: c++ -*-
#pragma once

#include <string>

#include <QGLViewer/qglviewer.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include "coonssurface.h"
#include "curve.h"

using qglviewer::Vec;
enum class ModelType { NONE, MESH, BEZIER_SURFACE, COONS_SURFACE, CURVE };

class MyViewer : public QGLViewer {
  Q_OBJECT

public:
  explicit MyViewer(QWidget *parent);
  virtual ~MyViewer();

  inline double getCutoffRatio() const;
  inline void setCutoffRatio(double ratio);
  inline double getMeanMin() const;
  inline void setMeanMin(double min);
  inline double getMeanMax() const;
  inline void setMeanMax(double max);
  inline const double *getSlicingDir() const;
  inline void setSlicingDir(double x, double y, double z);
  inline double getSlicingScaling() const;
  inline void setSlicingScaling(double scaling);
  bool openMesh(const std::string &filename, bool update_view = true);
  bool openBezier(const std::string &filename, bool update_view = true);
  bool saveBezier(const std::string &filename);

  //Coons
  bool openCoons(const std::string &filename, bool update_view = true);
  bool saveCoons(const std::string &filename);

  //Curve
  bool openCurve(const std::string &filename, bool update_view = true);

  ModelType getModelType() const;

signals:
  void startComputation(QString message);
  void midComputation(int percent);
  void endComputation();

protected:
  virtual void init() override;
  virtual void draw() override;
  virtual void drawWithNames() override;
  virtual void postSelection(const QPoint &p) override;
  virtual void keyPressEvent(QKeyEvent *e) override;
  virtual void mouseMoveEvent(QMouseEvent *e) override;
  virtual QString helpString() const override;

private:
  struct MyTraits : public OpenMesh::DefaultTraits {
      using Point  = OpenMesh::Vec3d; // the default would be Vec3f
      using Normal = OpenMesh::Vec3d;
      VertexTraits {
        double mean;              // approximated mean curvature
      };
    };
    using MyMesh = OpenMesh::TriMesh_ArrayKernelT<MyTraits>;
    using Vector = OpenMesh::VectorT<double,3>;

    // Mesh
    void updateMesh(bool update_mean_range = true);
    void updateVertexNormals();
  #ifdef USE_JET_FITTING
    void updateWithJetFit(size_t neighbors);
  #endif
    void localSystem(const Vector &normal, Vector &u, Vector &v);
    double voronoiWeight(MyMesh::HalfedgeHandle in_he);
    void updateMeanMinMax();
    void updateMeanCurvature();
    void updateAccurateMeanCurvature();
    void updateAlmostBoundary();

    // Bezier
    static void bernsteinAll(size_t n, double u, std::vector<double> &coeff);
    void generateMesh(size_t resolution);

    //Coons
    void generateCoonsMesh(size_t resolution);

    // Visualization
    void setupCamera();
    Vec meanMapColor(double d) const;
    Vec accurateMeanMapColor(double u, double v) const;
    void drawControlNet() const;
    void drawCurves() const;
    void drawCoonsControlNet() const;
    void drawAxes() const;
    void drawAxesWithNames() const;
    static Vec intersectLines(const Vec &ap, const Vec &ad, const Vec &bp, const Vec &bd);

    // Other
    void fairMesh();

    //////////////////////
    // Member variables //
    //////////////////////

    ModelType model_type;

    // Mesh
    MyMesh mesh;
    std::map<int, std::pair<double, double>> uvs;

    // Bezier
    size_t degree[2];
    std::vector<Vec> control_points;

    //Coons
    CoonsSurface coons;

    //Curve
    Curve curve;



  // Visualization
  double mean_min, mean_max, cutoff_ratio;
  bool show_control_points, show_solid, show_wireframe;
  bool vertexBoundary = true; /* Faces with edge at boundary means boundary face */
  std::map<OpenMesh::SmartFaceHandle, bool> almostBoundaryFaces;
  enum class Visualization { PLAIN, MEAN, ACCURATE_MEAN, SLICING, ISOPHOTES, ALMOSTBOUNDARY } visualization;
  GLuint isophote_texture, environment_texture, current_isophote_texture, slicing_texture;
  Vector slicing_dir;
  double slicing_scaling;
  int selected_vertex;
  struct ModificationAxes {
    bool shown;
    float size;
    int selected_axis;
    Vec position, grabbed_pos, original_pos;
  } axes;
  std::string last_filename;
};

#include "MyViewer.hpp"
