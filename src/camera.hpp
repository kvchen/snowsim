#ifndef CAMERA_H
#define CAMERA_H

#include <Eigen/Dense>
#include <iostream>

#include "math.h"

using namespace Eigen;

namespace SnowSimulator {

#define PI (3.14159265358979323)
#define EPS_D (0.00000000001)
#define EPS_F (0.00001f)
#define INF_D (std::numeric_limits<double>::infinity())
#define INF_F (std::numeric_limits<float>::infinity())

/*
  Takes any kind of number and converts from degrees to radians.
*/
template <typename T> inline T radians(T deg) { return deg * (PI / 180); }

/*
  Takes any kind of number and converts from radians to degrees.
*/
template <typename T> inline T degrees(T rad) { return rad * (180 / PI); }

/*
  Takes any kind of number, as well as a lower and upper bound, and clamps the
  number to be within the bound.
  NOTE: x, lo, and hi must all be the same type or compilation will fail. A
        common mistake is to pass an int for x and size_ts for lo and hi.
*/
template <typename T> inline T clamp(T x, T lo, T hi) {
  return std::min(std::max(x, lo), hi);
}

class Camera {
public:
  /*
    Sets the field of view to match screen screenW/H.
    NOTE: data and screenW/H will almost certainly disagree about the aspect
          ratio. screenW/H are treated as the source of truth, and the field
          of view is expanded along whichever dimension is too narrow.
    NOTE2: info.hFov and info.vFov are expected to be in DEGREES.
  */
  void configure(double nearClip, double farClip, double hFov, double vFov,
                 size_t screenW, size_t screenH);

  /*
    Phi and theta are in RADIANS.
  */
  void place(const Vector3d &targetPos, const double phi, const double theta,
             const double r, const double minR, const double maxR);

  std::string param_string() { return ""; }

  /*
    Copies just placement data from the other camera.
  */
  void copy_placement(const Camera &other);

  /*
    Updates the screen size to be the specified size, keeping screenDist
    constant.
  */
  void set_screen_size(const size_t screenW, const size_t screenH);

  /*
    Translates the camera such that a value at distance d directly in front of
    the camera moves by (dx, dy). Note that dx and dy are in screen coordinates,
    while d is in world-space coordinates (like pos/dir/up).
  */
  void move_by(const double dx, const double dy, const double d);

  /*
    Move the specified amount along the view axis.
  */
  void move_forward(const double dist);

  /*
    Rotate by the specified amount around the target.
  */
  void rotate_by(const double dPhi, const double dTheta);

  Vector3d position() const { return pos; }
  Vector3d view_point() const { return targetPos; }
  Vector3d up_dir() const { return c2w.col(1); }

  double v_fov() const { return vFov; }
  double aspect_ratio() const { return ar; }

  double near_clip() const { return nearClip; }
  double far_clip() const { return farClip; }

  // virtual void dump_settings(std::string filename);
  // virtual void load_settings(std::string filename);

private:
  // Computes pos, screenXDir, screenYDir from target, r, phi, theta.
  void compute_position();

  // Field of view aspect ratio, clipping planes.
  double hFov, vFov, ar, nearClip, farClip;

  // Current position and target point (the point the camera is looking at).
  Vector3d pos, targetPos;

  // Orientation relative to target, and min & max distance from the target.
  double phi, theta, r, minR, maxR;

  // camera-to-world rotation matrix (note: also need to translate a
  // camera-space point by 'pos' to perform a full camera-to-world
  // transform)
  Matrix3d c2w;

  // Info about screen to render to; it corresponds to the camera's full field
  // of view at some distance.
  size_t screenW, screenH;
  double screenDist;
};

} // namespace SnowSimulator

#endif // CAMERA_H
