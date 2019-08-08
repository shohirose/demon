#ifndef DEMON_VECTOR2D_HPP
#define DEMON_VECTOR2D_HPP

#include <cmath>

/// @brief 2D vector
struct Vector2d {
  /// @brief X component
  double x;
  /// @brief Y component
  double y;

  /// @name Constructors
  /// @{
  Vector2d() = default;
  Vector2d(double tx, double ty) : x{tx}, y{ty} {}
  Vector2d(const Vector2d&) = default;
  Vector2d(Vector2d&&) = default;
  /// @}

  Vector2d& operator=(const Vector2d&) = default;
  Vector2d& operator=(Vector2d&&) = default;

  /// @brief Returns squared norm
  auto squared_norm() const noexcept { return x * x + y * y; }

  /// @brief Returns norm
  auto norm() const noexcept { return std::sqrt(this->squared_norm()); }

  /// @brief Returns dot product
  auto dot(const Vector2d& other) const noexcept {
    return x * other.x + y * other.y;
  }

  /// @brief Returns zero vector
  static Vector2d zero() { return {0.0, 0.0}; }

  /// @brief Returns one vector
  static Vector2d one() { return {1.0, 1.0}; }

  /// @brief Returns constant vector
  static Vector2d constant(double c) { return {c, c}; }
};

inline Vector2d operator+(const Vector2d& v1, const Vector2d& v2) {
  return {v1.x + v2.x, v1.y + v2.y};
}

inline Vector2d operator-(const Vector2d& v1, const Vector2d& v2) {
  return {v1.x - v2.x, v1.y - v2.y};
}

inline Vector2d operator*(const Vector2d& v, double a) {
  return {v.x * a, v.y * a};
}

inline Vector2d operator*(double a, const Vector2d& v) {
  return {v.x * a, v.y * a};
}

inline Vector2d operator/(const Vector2d& v, double a) {
  return {v.x / a, v.y / a};
}

#endif  // DEMON_VECTOR2D_HPP