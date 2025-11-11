#include <iostream>

#include "rdr/accel.h"
#include "rdr/canary.h"
#include "rdr/interaction.h"
#include "rdr/math_aliases.h"
#include "rdr/platform.h"
#include "rdr/shape.h"

RDR_NAMESPACE_BEGIN

/* ===================================================================== *
 *
 * AABB Implementations
 *
 * ===================================================================== */

bool 
AABB::isOverlap (const AABB &other) const 
{
  return ((other.low_bnd[0] >= this->low_bnd[0] &&
           other.low_bnd[0] <= this->upper_bnd[0]) ||
           (this->low_bnd[0] >= other.low_bnd[0] &&
            this->low_bnd[0] <= other.upper_bnd[0])) &&
         ((other.low_bnd[1] >= this->low_bnd[1] &&
           other.low_bnd[1] <= this->upper_bnd[1]) ||
           (this->low_bnd[1] >= other.low_bnd[1] &&
            this->low_bnd[1] <= other.upper_bnd[1])) &&
         ((other.low_bnd[2] >= this->low_bnd[2] &&
           other.low_bnd[2] <= this->upper_bnd[2]) ||
           (this->low_bnd[2] >= other.low_bnd[2] &&
            this->low_bnd[2] <= other.upper_bnd[2]));
}

bool
AABB::intersect (const Ray &ray, Float *t_in, Float *t_out) const 
{
  // TODO(HW3): implement ray intersection with AABB.
  // ray distance for two intersection points are returned by pointers.
  //
  // This method should modify t_in and t_out as the "time"
  // when the ray enters and exits the AABB respectively.
  //
  // And return true if there is an intersection, false otherwise.
  //
  // Useful Functions:
  // @see Ray::safe_inverse_direction
  //    for getting the inverse direction of the ray.
  // @see Min/Max/ReduceMin/ReduceMax
  //    for vector min/max operations.
  
  Vec3f t_min, t_max;
  t_min = (low_bnd - ray.origin) / ray.direction;
  t_max = (upper_bnd - ray.origin) / ray.direction;

  if (t_min.x > t_max.x) std::swap(t_min.x, t_max.x);
  if (t_min.y > t_max.y) std::swap(t_min.y, t_max.y);
  if (t_min.z > t_max.z) std::swap(t_min.z, t_max.z);
  
  *t_in  = Max(t_min.x, t_min.y, t_min.z);
  *t_out = Min(t_max.x, t_max.y, t_max.z);

  if (*t_in > *t_out) return false;

  return true;
}

/* ===================================================================== *
 *
 * Accelerator Implementations
 *
 * ===================================================================== */

bool 
TriangleIntersect (Ray &ray, const uint32_t &triangle_index, 
                   const ref<TriangleMeshResource> &mesh, 
                   SurfaceInteraction &interaction) 
{
  using InternalScalarType = Double;
  using InternalVecType    = Vec<InternalScalarType, 3>;

  AssertAllValid(ray.direction, ray.origin);
  AssertAllNormalized(ray.direction);

  const auto &vertices = mesh->vertices;
  const Vec3u v_idx(&mesh->v_indices[3 * triangle_index]);
  assert(v_idx.x < mesh->vertices.size());
  assert(v_idx.y < mesh->vertices.size());
  assert(v_idx.z < mesh->vertices.size());

  InternalVecType d   = Cast<InternalScalarType>(ray.direction);
  InternalVecType o   = Cast<InternalScalarType>(ray.origin);
  InternalVecType v0  = Cast<InternalScalarType>(vertices[v_idx[0]]);
  InternalVecType v1  = Cast<InternalScalarType>(vertices[v_idx[1]]);
  InternalVecType v2  = Cast<InternalScalarType>(vertices[v_idx[2]]);

  // TODO(HW3): implement ray-triangle intersection test.
  // You should compute the u, v, t as InternalScalarType
  //
  //   InternalScalarType u = ...;
  //   InternalScalarType v = ...;
  //   InternalScalarType t = ...;
  //
  // And exit early with `return false` if there is no intersection.
  //
  // The intersection points is denoted as:
  // (1 - u - v) * v0 + u * v1 + v * v2 == ray.origin + t * ray.direction
  // where the left side is the barycentric interpolation of the triangle
  // vertices, and the right side is the parametric equation of the ray.
  //
  // You should also make sure that:
  // u >= 0, v >= 0, u + v <= 1, and, ray.t_min <= t <= ray.t_max
  //
  // Useful Functions:
  // You can use @see Cross and @see Dot for determinant calculations.

  InternalVecType e1 = v1 - v0;
  InternalVecType e2 = v2 - v0;
  InternalVecType s = o - v0;
  InternalVecType s1 = Cross(d, e2);
  InternalVecType s2 = Cross(s, e1);
  InternalScalarType coef = 1 / Dot(s1, e1);
  InternalScalarType t = coef * Dot(s2, e2);
  InternalScalarType u = coef * Dot(s1, s);
  InternalScalarType v = coef * Dot(s2, d);

  if (!(u >= 0 && v >= 0 && u + v <= 1 && t >= ray.t_min && t <= ray.t_max))
    return false;

  // We will reach here if there is an intersection

  CalculateTriangleDifferentials(interaction,
                                 {static_cast<Float>(1 - u - v), 
                                  static_cast<Float>(u),
                                  static_cast<Float>(v)},
                                 mesh, triangle_index);
  AssertNear(interaction.p, ray(t));
  assert(ray.withinTimeRange(t));
  ray.setTimeMax(t);
  return true;
}

void Accel::setTriangleMesh(const ref<TriangleMeshResource> &mesh) {
  // Build the bounding box
  AABB bound(Vec3f(Float_INF, Float_INF, Float_INF),
      Vec3f(Float_MINUS_INF, Float_MINUS_INF, Float_MINUS_INF));
  for (auto &vertex : mesh->vertices) {
    bound.low_bnd   = Min(bound.low_bnd, vertex);
    bound.upper_bnd = Max(bound.upper_bnd, vertex);
  }

  this->mesh  = mesh;   // set the pointer
  this->bound = bound;  // set the bounding box
}

void Accel::build() {}

AABB Accel::getBound() const {
  return bound;
}

bool Accel::intersect(Ray &ray, SurfaceInteraction &interaction) const {
  bool success = false;
  for (int i = 0; i < mesh->v_indices.size() / 3; i++)
    success |= TriangleIntersect(ray, i, mesh, interaction);
  return success;
}

RDR_NAMESPACE_END
