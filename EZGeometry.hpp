#ifndef EZGEOMETRY_EZGEOMETRY_HPP
#define EZGEOMETRY_EZGEOMETRY_HPP

#include "cfloat"
#include "cmath"
#include "cstdlib"
#include "iostream"

namespace ezgeometry {
inline double Min(double a, double b) { return (a > b ? b : a); }
inline double Max(double a, double b) { return (a > b ? a : b); }

inline bool EpsilonEquals(double a, double b, double epsilon) {
  return fabs(a - b) < epsilon;
}

inline double ToDegrees(double radians) { return 180.0 * radians / M_PI; }
inline double ToRadians(double degrees) { return degrees * M_PI / 180.0; }

inline double Signum(double value) {
  if (EpsilonEquals(value, 0.0, DBL_EPSILON)) {
    return 0.0;
  } else if (value > 0) {
    return 1.0;
  } else {
    return -1.0;
  }
}

namespace geometry2d {
class Twist2d {
public:
  explicit Twist2d() : dx_{0.0}, dy_{0.0}, dtheta_{0.0} {}

  explicit Twist2d(double dx, double dy, double dtheta)
      : dx_{dx}, dy_{dy}, dtheta_{dtheta} {}

  Twist2d(const Twist2d &other) = default;

  static Twist2d Identity() { return Twist2d(0.0, 0.0, 0.0); }

  Twist2d Scaled(double scale) const {
    return Twist2d(dx_ * scale, dy_ * scale, dtheta_ * scale);
  }

  double Norm() const {
    // Common case of dy == 0
    if (dy_ == 0.0)
      return fabs(dx_);
    return hypot(dx_, dy_);
  }

  double Curvature() const {
    if (fabs(dtheta_) < EPSILON_VALUE && Norm() < EPSILON_VALUE)
      return 0.0;
    return dtheta_ / Norm();
  }

  double Dx() const { return dx_; }

  double Dy() const { return dy_; }

  double Dtheta() const { return dtheta_; }

  void Print() const {
    std::cout << "Twist2d:" << std::endl
              << "  "
              << "dx->" << dx_ << "  dy->" << dy_ << std::endl;
  }

private:
  static constexpr double EPSILON_VALUE = DBL_EPSILON;

  const double dx_;
  const double dy_;
  const double dtheta_; // Radians
};

class Rotation2d {
public:
  explicit Rotation2d() {
    cos_angle_ = 1.0;
    sin_angle_ = 0.0;
    theta_degrees_ = 0.0;
  }

  explicit Rotation2d(double x, double y, bool normalize) {
    if (normalize) {
      // From trig, we know that sin^2 + cos^2 == 1, but as we do math on
      // this object we might accumulate rounding errors. Normalizing
      // forces us to re-scale the sin and cos to reset rounding errors.
      double magnitude = hypot(x, y);
      if (magnitude > EPSILON_VALUE) {
        sin_angle_ = y / magnitude;
        cos_angle_ = x / magnitude;
      } else {
        sin_angle_ = 0;
        cos_angle_ = 1;
      }
    } else {
      cos_angle_ = x;
      sin_angle_ = y;
    }
    theta_degrees_ = ToDegrees(atan2(sin_angle_, cos_angle_));
  }

  explicit Rotation2d(const double theta_degrees) {
    cos_angle_ = cos(ToRadians(theta_degrees));
    sin_angle_ = sin(ToRadians(theta_degrees));
    theta_degrees_ = theta_degrees;
  }

  Rotation2d(const Rotation2d &other) {
    cos_angle_ = other.cos_angle_;
    sin_angle_ = other.sin_angle_;
    theta_degrees_ = ToDegrees(atan2(sin_angle_, cos_angle_));
  }

  static Rotation2d Identity() { return Rotation2d(1.0, 0.0, false); }

  static Rotation2d FromRadians(double angle_radians) {
    return Rotation2d(cos(angle_radians), sin(angle_radians), false);
  }

  static Rotation2d FromDegrees(double angle_degrees) {
    return Rotation2d(angle_degrees);
  }

  double Cos() const { return cos_angle_; }

  double Sin() const { return sin_angle_; }

  double Tan() const {
    if (fabs(cos_angle_) < EPSILON_VALUE) {
      if (sin_angle_ >= 0.0) {
        return 1e9;
      } else {
        return -1e9;
      }
    }
    return sin_angle_ / cos_angle_;
  }

  double GetRadians() const { return atan2(sin_angle_, cos_angle_); }

  double GetDegrees() const { return ToDegrees(GetRadians()); }

  double GetUnboundedDegrees() const { return theta_degrees_; }

  /**
   * We can rotate this Rotation2d by adding together the effects of it and
   * another rotation.
   *
   * @param other The other rotation. See:
   *              https://en.wikipedia.org/wiki/Rotation_matrix
   * @return This rotation rotated by other.
   */
  Rotation2d RotateBy(const Rotation2d &other) const {
    return Rotation2d(
        cos_angle_ * other.cos_angle_ - sin_angle_ * other.sin_angle_,
        cos_angle_ * other.sin_angle_ + sin_angle_ * other.cos_angle_, true);
  }

  Rotation2d Normal() const {
    return Rotation2d(-sin_angle_, cos_angle_, false);
  }

  /**
   * The inverse of a Rotation2d "undoes" the effect of this rotation.
   *
   * @return The opposite of this rotation.
   */
  Rotation2d Inverse() const {
    return Rotation2d(cos_angle_, -sin_angle_, false);
  }

  /**
   * @return The pole nearest to this rotation.
   */
  Rotation2d NearestPole() const {
    double pole_sin = 0.0;
    double pole_cos = 0.0;
    if (fabs(cos_angle_) > fabs(sin_angle_)) {
      pole_cos = Signum(cos_angle_);
      pole_sin = 0.0;
    } else {
      pole_cos = 0.0;
      pole_sin = Signum(sin_angle_);
    }
    return Rotation2d(pole_cos, pole_sin, false);
  }

  Rotation2d Interpolate(const Rotation2d &other, double x) const {
    if (x <= 0) {
      return *this;
    } else if (x >= 1) {
      return other;
    }
    double angle_diff = Inverse().RotateBy(other).GetRadians();
    return RotateBy(Rotation2d::FromRadians(angle_diff * x));
  }

  double Distance(const Rotation2d &other) const {
    return Inverse().RotateBy(other).GetRadians();
  }

  bool EpsilonEquals(const Rotation2d &other, const double epsilon) const {
    return Distance(other) < epsilon;
  }

  bool Equals(const Rotation2d &other) const {
    return Distance(other) < EPSILON_VALUE;
  }

  void Print() const {
    std::cout << "Rotation2d:" << std::endl
              << "  "
              << "cos_angle->" << cos_angle_ << "  sin_angle->" << sin_angle_
              << "  theta_degrees->" << theta_degrees_ << std::endl;
  }

private:
  static constexpr double EPSILON_VALUE = DBL_EPSILON;

  double cos_angle_;
  double sin_angle_;
  double theta_degrees_;
};

class Translation2d {
public:
  explicit Translation2d() {
    x_ = 0;
    y_ = 0;
  }

  explicit Translation2d(double x, double y) {
    x_ = x;
    y_ = y;
  }

  explicit Translation2d(const Translation2d &start, const Translation2d &end) {
    x_ = end.x_ - start.x_;
    y_ = end.y_ - start.y_;
  }

  Translation2d(const Translation2d &other) = default;

  static Translation2d Identity() {
    static Translation2d identity = Translation2d(0.0, 0.0);
    return identity;
  }

  static Translation2d FromPolar(const Rotation2d &direction,
                                 double magnitude) {
    return Translation2d(direction.Cos() * magnitude,
                         direction.Sin() * magnitude);
  }

  /**
   * The "norm" of a transform is the Euclidean distance in x and y.
   *
   * @return sqrt(x ^ 2 + y ^ 2)
   */
  double Norm() const { return hypot(x_, y_); }

  /**
   * Normalizing a vector scales it so that its norm is 1 while maintaining
   * its direction. If input is a zero vector, return a zero vector.
   *
   * @return r / norm(r) or (0,0)
   */
  Translation2d Normalize() const {
    if (EpsilonEquals(Identity(), EPSILON_VALUE))
      return *this;
    return Scale(1.0 / Norm());
  }

  double X() const { return x_; }

  double Y() const { return y_; }

  void SetX(double x) { x_ = x; }

  void SetY(double y) { y_ = y; }

  /**
   * We can compose Translation2d's by adding together the x and y shifts.
   *
   * @param other The other translation to add.
   * @return The combined effect of translating by this object and the other.
   */
  Translation2d TranslateBy(const Translation2d &other) const {
    return Translation2d(x_ + other.x_, y_ + other.y_);
  }

  /**
   * We can also rotate Translation2d's. See:
   * https://en.wikipedia.org/wiki/Rotation_matrix
   *
   * @param rotation The rotation to apply.
   * @return This translation rotated by rotation.
   */
  Translation2d RotateBy(const Rotation2d &rotation) const {
    return Translation2d(x_ * rotation.Cos() - y_ * rotation.Sin(),
                         x_ * rotation.Sin() + y_ * rotation.Cos());
  }

  Rotation2d Direction() const { return Rotation2d(x_, y_, true); }

  /**
   * The inverse simply means a Translation2d that "undoes" this object.
   *
   * @return Translation by -x and -y.
   */
  Translation2d Inverse() const { return Translation2d(-x_, -y_); }

  Translation2d Interpolate(const Translation2d &other, double x) const {
    if (x <= 0) {
      return *this;
    } else if (x >= 1) {
      return other;
    }
    return Extrapolate(other, x);
  }

  Translation2d Extrapolate(const Translation2d &other, double x) const {
    return Translation2d(x * (other.x_ - x_) + x_, x * (other.y_ - y_) + y_);
  }

  Translation2d Scale(double s) const { return Translation2d(x_ * s, y_ * s); }

  bool EpsilonEquals(const Translation2d &other, const double epsilon) const {
    return ezgeometry::EpsilonEquals(X(), other.X(), epsilon) &&
           ezgeometry::EpsilonEquals(Y(), other.Y(), epsilon);
  }

  static double Dot(const Translation2d &a, const Translation2d &b) {
    return a.x_ * b.x_ + a.y_ * b.y_;
  }

  /**
   * The scalar projection of a vector u onto a vector v is the length of the
   * "shadow" cast by u onto v under a "light" that is placed on a line normal
   * to v and containing the endpoint of u, given that u and v share a
   * starting point. tl;dr: _* u /| / | / | v *---+---------------->* \___/ |
   * scal_v(u) u.scal(v)
   *
   * @return (u . v) / norm(v)
   */
  double Scal(const Translation2d &v) const { return Dot(*this, v) / v.Norm(); }

  /**
   * The projection of a vector u onto a vector v is the vector in the
   * direction of v with the magnitude u.scal(v).
   *
   * @return u.scal(v) * v / norm(v)
   */
  Translation2d Proj(const Translation2d &v) const {
    return v.Normalize().Scale(Scal(v));
  }

  bool IsWithinAngle(Translation2d A, const Translation2d &B,
                     const Translation2d &C, bool vertical) const {
    Translation2d M = A.Interpolate(C, 0.5);               // midpoint
    Translation2d m = Translation2d(B, M).Normalize();     // mid-vector
    Translation2d a = Translation2d(B, A).Normalize();     // side vector
    Translation2d d = Translation2d(B, *this).Normalize(); // vector to here
    if (vertical) {
      m = m.Inverse();
      a = a.Inverse();
    }
    return Translation2d::Dot(d, m) > Translation2d::Dot(a, m);
  }

  bool IsWithinAngle(const Translation2d &A, const Translation2d &B,
                     const Translation2d &C) const {
    return IsWithinAngle(A, B, C, false);
  }

  /** Assumes an angle centered at the origin. */
  bool IsWithinAngle(const Translation2d &A, const Translation2d &C,
                     bool vertical) const {
    return IsWithinAngle(A, Identity(), C, vertical);
  }

  bool IsWithinAngle(const Translation2d &A, const Translation2d &C) const {
    return IsWithinAngle(A, C, false);
  }

  static Rotation2d GetAngle(const Translation2d &a, const Translation2d &b) {
    double cos_angle = Dot(a, b) / (a.Norm() * b.Norm());
    return Rotation2d::FromRadians(acos(Min(1.0, Max(cos_angle, -1.0))));
  }

  static double Cross(const Translation2d &a, const Translation2d &b) {
    return a.x_ * b.y_ - a.y_ * b.x_;
  }

  /**
   * The distance between a point and a line can be computed as a scalar
   * projection.
   *
   * @param Translation2d a One point on the line.
   * @param Translation2d b Another point on the line.
   */
  double DistanceToLine(const Translation2d &a, const Translation2d &b) const {
    Translation2d point = Translation2d(a, *this);
    Translation2d line = Translation2d(a, b);
    Translation2d perpLine = line.RotateBy(Rotation2d(90.0));
    return fabs(point.Scal(perpLine));
  }

  double Distance(const Translation2d &other) const {
    return Inverse().TranslateBy(other).Norm();
  }

  bool Equals(const Translation2d &other) const {
    return Distance(other) < EPSILON_VALUE;
  }

  void Print() const {
    std::cout << "Translation2d:" << std::endl
              << "  "
              << "x->" << x_ << "  y->" << y_ << std::endl;
  }

private:
  static constexpr double EPSILON_VALUE = DBL_EPSILON;

  double x_ = 0.0;
  double y_ = 0.0;
};

class Pose2d {
public:
  explicit Pose2d() : translation_{}, rotation_{} {}

  explicit Pose2d(double x, double y, const Rotation2d &rotation)
      : translation_{x, y}, rotation_{rotation} {}

  explicit Pose2d(const Translation2d &translation, const Rotation2d &rotation)
      : translation_{translation}, rotation_{rotation} {}

  Pose2d(const Pose2d &other) = default;

  static Pose2d Identity() { return Pose2d(); }

  static Pose2d FromTranslation(const Translation2d &translation) {
    return Pose2d(translation, Rotation2d());
  }

  static Pose2d FromRotation(const Rotation2d &rotation) {
    return Pose2d(Translation2d(), rotation);
  }

  /**
   * Obtain a new Pose2d from a (constant curvature) velocity. See:
   * https://github.com/strasdat/Sophus/blob/master/sophus/se2.hpp
   */
  static Pose2d Exp(const Twist2d &delta) {
    double sin_theta = sin(delta.Dtheta());
    double cos_theta = cos(delta.Dtheta());
    double s, c;
    if (fabs(delta.Dtheta()) < EPSILON_VALUE) {
      s = 1.0 - 1.0 / 6.0 * delta.Dtheta() * delta.Dtheta();
      c = .5 * delta.Dtheta();
    } else {
      s = sin_theta / delta.Dtheta();
      c = (1.0 - cos_theta) / delta.Dtheta();
    }
    return Pose2d(Translation2d(delta.Dx() * s - delta.Dy() * c,
                                delta.Dx() * c + delta.Dy() * s),
                  Rotation2d(cos_theta, sin_theta, false));
  }

  /**
   * Logical inverse of the above.
   */
  static Twist2d Log(const Pose2d &transform) {
    const double dtheta = transform.GetRotation().GetRadians();
    const double half_dtheta = 0.5 * dtheta;
    const double cos_minus_one = transform.GetRotation().Cos() - 1.0;
    double half_theta_by_tan_of_half_dtheta;
    if (fabs(cos_minus_one) < EPSILON_VALUE) {
      half_theta_by_tan_of_half_dtheta = 1.0 - 1.0 / 12.0 * dtheta * dtheta;
    } else {
      half_theta_by_tan_of_half_dtheta =
          -(half_dtheta * transform.GetRotation().Sin()) / cos_minus_one;
    }
    const Translation2d translation_part = transform.GetTranslation().RotateBy(
        Rotation2d(half_theta_by_tan_of_half_dtheta, -half_dtheta, false));
    return Twist2d(translation_part.X(), translation_part.Y(), dtheta);
  }

  Translation2d GetTranslation() const { return translation_; }

  Rotation2d GetRotation() const { return rotation_; }

  /**
   * Transforming this RigidTransform2d means first translating by
   * other.translation and then rotating by other.rotation
   *
   * @param other The other transform.
   * @return This transform * other
   */
  Pose2d TransformBy(const Pose2d &other) const {
    return Pose2d(
        translation_.TranslateBy(other.translation_.RotateBy(rotation_)),
        rotation_.RotateBy(other.rotation_));
  }

  /**
   * The inverse of this transform "undoes" the effect of translating by this
   * transform.
   *
   * @return The opposite of this transform.
   */
  Pose2d Inverse() const {
    Rotation2d rotation_inverted = rotation_.Inverse();
    return Pose2d(translation_.Inverse().RotateBy(rotation_inverted),
                  rotation_inverted);
  }

  Pose2d Normal() const { return Pose2d(translation_, rotation_.Normal()); }

  /**
   * Finds the point where the heading of this pose intersects the heading of
   * another. Returns (+INF, +INF) if parallel.
   */
  Translation2d Intersection(const Pose2d &other) const {
    const Rotation2d other_rotation = other.GetRotation();
    if (rotation_.EpsilonEquals(other_rotation, EPSILON_VALUE)) {
      // Lines are parallel.
      return Translation2d(1e9, 1e9);
    }
    if (fabs(rotation_.Cos()) < fabs(other_rotation.Cos())) {
      return IntersectionInternal(*this, other);
    } else {
      return IntersectionInternal(other, *this);
    }
  }

  bool EpsilonEquals(const Pose2d &other, const double epsilon) const {
    return GetTranslation().EpsilonEquals(other.GetTranslation(), epsilon) &&
           GetRotation().EpsilonEquals(other.GetRotation(), epsilon);
  }

  /**
   * Do twist interpolation of this pose assuming constant curvature.
   */
  Pose2d Interpolate(const Pose2d &other, double x) const {
    if (x <= 0) {
      return *this;
    } else if (x >= 1) {
      return other;
    }
    const Twist2d twist = Pose2d::Log(Inverse().TransformBy(other));
    return TransformBy(Pose2d::Exp(twist.Scaled(x)));
  }

  double Distance(const Pose2d &other) const {
    return Pose2d::Log(Inverse().TransformBy(other)).Norm();
  }

  bool Equals(const Pose2d &other) const {
    return EpsilonEquals(other, EPSILON_VALUE);
  }

  Pose2d Mirror() const {
    return Pose2d(Translation2d(GetTranslation().X(), -GetTranslation().Y()),
                  GetRotation().Inverse());
  }

  void Print() const {
    std::cout << "Pose2d:" << std::endl;
    std::cout << "  ";
    translation_.Print();
    std::cout << "  ";
    rotation_.Print();
  }

private:
  static constexpr double EPSILON_VALUE = DBL_EPSILON;

  Translation2d translation_;
  Rotation2d rotation_;

  static Translation2d IntersectionInternal(const Pose2d &a, const Pose2d &b) {
    const Rotation2d a_r = a.GetRotation();
    const Rotation2d b_r = b.GetRotation();
    const Translation2d a_t = a.GetTranslation();
    const Translation2d b_t = b.GetTranslation();

    const double tan_b = b_r.Tan();
    const double t = ((a_t.X() - b_t.X()) * tan_b + b_t.Y() - a_t.Y()) /
                     (a_r.Sin() - a_r.Cos() * tan_b);
    if (ezgeometry::EpsilonEquals(t, 0.0, EPSILON_VALUE)) {
      return Translation2d(1e9, -1e9);
    }
    return a_t.TranslateBy(Translation2d::FromPolar(a_r, 1.0));
  }
};
} // namespace geometry2d
} // namespace ezgeometry

#endif // EZGEOMETRY_EZGEOMETRY_HPP
