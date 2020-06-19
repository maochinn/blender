#pragma once

#include "DRW_render.h"

#include "BLI_math.h"
#include "DNA_listBase.h"

#include "custom_ray_tracer.h"

/* local function */

/* https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates */
void barycentric(float *u,
                 float *v,
                 float *w,
                 const float p[3],
                 const float a[3],
                 const float b[3],
                 const float c[3])
{
  // Compute barycentric coordinates (u, v, w) for
  // point p with respect to triangle (a, b, c)

  float v0[3], v1[3], v2[3];
  // Vector v0 = b - a, v1 = c - a, v2 = p - a;
  sub_v3_v3v3(v0, b, a);
  sub_v3_v3v3(v1, c, a);
  sub_v3_v3v3(v2, p, a);

  float d00 = dot_v3v3(v0, v0);
  float d01 = dot_v3v3(v0, v1);
  float d11 = dot_v3v3(v1, v1);
  float d20 = dot_v3v3(v2, v0);
  float d21 = dot_v3v3(v2, v1);
  float denom = d00 * d11 - d01 * d01;
  *v = (d11 * d20 - d01 * d21) / denom;
  *w = (d00 * d21 - d01 * d20) / denom;
  *u = 1.0f - *v - *w;
}

inline float schlick(float cosine, float ref_idx)
{
  float r0 = (1 - ref_idx) / (1 + ref_idx);
  r0 = r0 * r0;
  return r0 + (1 - r0) * pow((1 - cosine), 5);
}

bool refract(float refracted[3], const float v[3], const float n[3], float ni_over_nt)
{
  float uv[3];
  normalize_v3_v3(uv, v);
  float dt = dot_v3v3(uv, n);
  float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
  if (discriminant > 0) {
    // refracted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
    madd_v3_v3v3fl(refracted, uv, n, -dt);
    mul_v3_fl(refracted, ni_over_nt);
    madd_v3_v3fl(refracted, n, -sqrt(discriminant));
    return true;
  }
  else
    return false;
}

void getSphereUV(float *u, float *v, const CUSTOM_Vector *p)
{
  float phi = atan2(p->z, p->x);
  float theta = asin(p->y);
  *u = 1.0f - (phi + M_PI) / (2.0f * M_PI);
  *v = (theta + M_PI / 2) / M_PI;
}
/******************************************************************
 * get a point on a ray by parameter
 * Output:
 *  p:      CUSTOM_Vector vec3
 * Input:
 *  ray:    pointer of const CUSTOM_Ray*
 *  sphere: pointer of const CUSTOM_Sphere*
 *  t_min:  float, min of t parameter of ray
 *  t_max:  float, max of t parameter of ray
 *  rec:    pointer of CUSTOM_HitRecord*
 */
CUSTOM_Vector pointAtParameter(const CUSTOM_Ray *ray, float t)
{
  CUSTOM_Vector p;
  // A + t * B;
  madd_v3_v3v3fl(p.vec3, ray->origin.vec3, ray->direction.vec3, t);
  return p;
}
/******************************************************************
 * check hit by a ray through whole world
 * Output:
 *  hit:  bool
 * Input:
 *  ray:    pointer of const CUSTOM_Ray*
 *  sphere:  pointer of const CUSTOM_Sphere*
 *  t_min:  float, min of t parameter of ray
 *  t_max:  float, max of t parameter of ray
 *  rec:    pointer of CUSTOM_HitRecord*
 */

bool CUSTOM_HittableHit(CUSTOM_HitRecord *rec,
                        const CUSTOM_Ray *ray,
                        const CUSTOM_Hittable *hittable,
                        float t0,
                        float t1)
{
  if (hittable->type == CUSTOM_SPHERE) {
    return CUSTOM_SphereHit(rec, ray, (CUSTOM_Sphere *)hittable, t0, t1);
  }
  else if (hittable->type == CUSTOM_RECT_XY) {
    return CUSTOM_RectXYHit(rec, ray, (CUSTOM_RectXY *)hittable, t0, t1);
  }
  else if (hittable->type == CUSTOM_RECT_XZ) {
    return CUSTOM_RectXZHit(rec, ray, (CUSTOM_RectXZ *)hittable, t0, t1);
  }
  else if (hittable->type == CUSTOM_RECT_YZ) {
    return CUSTOM_RectYZHit(rec, ray, (CUSTOM_RectYZ *)hittable, t0, t1);
  }
  else if (hittable->type == CUSTOM_BOX) {
    return CUSTOM_BoxHit(rec, ray, (CUSTOM_Box *)hittable, t0, t1);
  }
  else if (hittable->type == CUSTOM_TRIANGLE) {
    return CUSTOM_TriangleHit(rec, ray, (CUSTOM_Triangle *)hittable, t0, t1);
  }
  else if (hittable->type == CUSTOM_BVH_NODE) {
    return CUSTOM_BvhNodeHit(rec, ray, (CUSTOM_BvhNode *)hittable, t0, t1);
  }

  return false;
}

bool CUSTOM_HittableBoundingBox(CUSTOM_AABB *output_box, const CUSTOM_Hittable *hittable)
{
  if (hittable->type == CUSTOM_SPHERE) {
    return CUSTOM_SphereBoundingBox(output_box, (CUSTOM_Sphere *)hittable);
  }
  else if (hittable->type == CUSTOM_BOX) {
    return CUSTOM_BoxBoundingBox(output_box, (CUSTOM_Box *)hittable);
  }
  else if (hittable->type == CUSTOM_TRIANGLE) {
    return CUSTOM_TriangleBoundingBox(output_box, (CUSTOM_Triangle *)hittable);
  }
  else if (hittable->type == CUSTOM_BVH_NODE) {
    return CUSTOM_BvhNodeBoundingBox(output_box, (CUSTOM_BvhNode *)hittable);
  }
  return false;
}

/******************************************************************
 * Generate ramdom vector, whose length < 1.0f
 * Output:
 *  p:  CUSTOM_Vector vec3
 */
CUSTOM_Vector randomInUnitSphere(void)
{
  CUSTOM_Vector p;
  do {
    // p = 2.0 * vec3(random_float(), random_float(), random_float()) - vec3(1, 1, 1);
    mul_v3_v3fl(p.vec3, (float[3]){randomFloat(), randomFloat(), randomFloat()}, 2.0f);
    sub_v3_v3(p.vec3, (float[3]){1.0f, 1.0f, 1.0f});
  } while (len_squared_v3(p.vec3) >= 1.0);
  return p;
}

/**/
CUSTOM_Vector randomInUnitHemisphere(const CUSTOM_Vector *dir)
{
  CUSTOM_Vector p = randomInUnitSphere();
  if (dot_v3v3(p.vec3, dir->vec3) < 0.0f) {
    negate_v3(p.vec3);
  }
  return p;
}

CUSTOM_Vector randomInUnitDisk(void)
{
  CUSTOM_Vector p;
  do {
    // p = 2.0 * vec3(random_double(), random_double(), 0) - vec3(1, 1, 0);
    mul_v3_v3fl(p.vec3, (float[3]){randomFloat(), randomFloat(), 0.0f}, 2.0f);
    sub_v3_v3(p.vec3, (float[3]){1.0f, 1.0f, 0.0f});
  } while (dot_v3v3(p.vec3, p.vec3) >= 1.0);
  return p;
}

/******************************************************************
 * check hit by a ray through whole world
 * Output:
 *  hit_anything:  bool
 * Input:
 *  ray:    pointer of const CUSTOM_Ray*
 *  world:  pointer of const ListBase*
 *  t_min:  float, min of t parameter of ray
 *  t_max:  float, max of t parameter of ray
 *  rec:    pointer of CUSTOM_HitRecord*
 */
bool list_hit(
    CUSTOM_HitRecord *rec, const CUSTOM_Ray *ray, const ListBase *list, float t_min, float t_max)
{
  CUSTOM_HitRecord temp_rec;
  bool hit_anything = false;
  double closest_so_far = t_max;
  for (CUSTOM_Hittable *hittable = (CUSTOM_Hittable *)list->first; hittable;
       hittable = hittable->next) {
    if (CUSTOM_HittableHit(&temp_rec, ray, hittable, t_min, closest_so_far)) {
      hit_anything = true;
      closest_so_far = temp_rec.t;
      *rec = temp_rec;
    }
  }
  return hit_anything;
}

CUSTOM_AABB surroundingBox(const CUSTOM_AABB *box0, const CUSTOM_AABB *box1)
{
  // point3 small(fmin(box0.min().x(), box1.min().x()),
  //             fmin(box0.min().y(), box1.min().y()),
  //             fmin(box0.min().z(), box1.min().z()));

  // point3 big(fmax(box0.max().x(), box1.max().x()),
  //           fmax(box0.max().y(), box1.max().y()),
  //           fmax(box0.max().z(), box1.max().z()));

  CUSTOM_Vector small = {fmin(box0->min.x, box1->min.x),
                         fmin(box0->min.y, box1->min.y),
                         fmin(box0->min.z, box1->min.z)};
  CUSTOM_Vector big = {fmax(box0->max.x, box1->max.x),
                       fmax(box0->max.y, box1->max.y),
                       fmax(box0->max.z, box1->max.z)};

  return (CUSTOM_AABB){small, big};
}
bool list_boundingBox(CUSTOM_AABB *output_box, const ListBase *list, float t_min, float t_max)
{
  if (BLI_listbase_is_empty(list))
    return false;

  CUSTOM_AABB temp_box;
  bool first_box = true;

  for (CUSTOM_Hittable *hittable = (CUSTOM_Hittable *)list->first; hittable;
       hittable = hittable->next) {
    if (!CUSTOM_HittableBoundingBox(&temp_box, hittable, t_min, t_max)) {
      return false;
    }
    *output_box = first_box ? temp_box : (surroundingBox(output_box, &temp_box));
    first_box = false;
  }
  return true;
}

bool AABB_hit(const CUSTOM_AABB *aabb, const CUSTOM_Ray *ray, float t_min, float t_max)
{
  for (int a = 0; a < 3; a++) {
    // float t0 = fmin((aabb->min.vec3[a] - ray->origin.vec3[a]) / ray->direction.vec3[a],
    //          (aabb->max.vec3[a] - ray->origin.vec3[a]) / ray->direction.vec3[a]);
    // float t1 = fmax((aabb->min.vec3[a] - ray->origin.vec3[a]) / ray->direction.vec3[a],
    //          (aabb->max.vec3[a] - ray->origin.vec3[a]) / ray->direction.vec3[a]);

    // t_min = fmax(t0, t_min);
    // t_max = fmin(t1, t_max);

    float invD = 1.0f / ray->direction.vec3[a];
    float t0 = (aabb->min.vec3[a] - ray->origin.vec3[a]) * invD;
    float t1 = (aabb->max.vec3[a] - ray->origin.vec3[a]) * invD;
    if (invD < 0.0f) {
      SWAP(float, t0, t1);
      // float temp = t0;
      // t0 = t1;
      // t1 = temp;
    };
    t_min = t0 > t_min ? t0 : t_min;
    t_max = t1 < t_max ? t1 : t_max;

    // if it is plane, t_max == t_min will wrong
    // if (t_max <= t_min)
    if (t_max < t_min)
      return false;
  }
  return true;
}

int boxCompareX(const void *a, const void *b)
{
  return boxCompare(*(CUSTOM_Hittable **)a, *(CUSTOM_Hittable **)b, 0);
}

int boxCompareY(const void *a, const void *b)
{
  return boxCompare(*(CUSTOM_Hittable **)a, *(CUSTOM_Hittable **)b, 1);
}

int boxCompareZ(const void *a, const void *b)
{
  return boxCompare(*(CUSTOM_Hittable **)a, *(CUSTOM_Hittable **)b, 2);
}

CUSTOM_Vector textureValue(const CUSTOM_Texture *texture,
                           float u,
                           float v,
                           const CUSTOM_Vector *pos)
{
  if (texture->type == CUSTOM_CONSTANT_TEXTURE)
    return CUSTOM_ConstantTextureValue((CUSTOM_ConstantTexture *)texture);
  else if (texture->type == CUSTOM_CHECKER_TEXTURE)
    return CUSTOM_CheckerTextureValue((CUSTOM_CheckerTexture *)texture, u, v, pos);
  else if (texture->type == CUSTOM_IMAGE_TEXTURE)
    return CUSTOM_ImageTextureValue((CUSTOM_ImageTexture *)texture, u, v, pos);
  else
    return (CUSTOM_Vector){0.0f, 0.0f, 0.0f};
}

bool materialScatter(CUSTOM_Ray *scattered,
                     CUSTOM_Vector *attenuation,
                     const CUSTOM_Ray *ray_in,
                     const CUSTOM_HitRecord *rec,
                     const CUSTOM_Material *mat)
{
  if (rec->mat_ptr->type == CUSTOM_LAMBERTIAN)
    return CUSTOM_LambertianScatter(
        scattered, attenuation, rec, (CUSTOM_Lambertian *)rec->mat_ptr);
  else if (rec->mat_ptr->type == CUSTOM_METAL)
    return CUSTOM_MetalScatter(scattered, attenuation, ray_in, rec, (CUSTOM_Metal *)rec->mat_ptr);
  else if (rec->mat_ptr->type == CUSTOM_DIELECTRIC)
    return CUSTOM_DielectricScatter(
        scattered, attenuation, ray_in, rec, (CUSTOM_Dielectric *)rec->mat_ptr);
  else if (rec->mat_ptr->type == CUSTOM_DIFFUSE_LIGHT)
    return CUSTOM_DiffuseLightScatter((CUSTOM_DiffuseLight *)rec->mat_ptr);
  else
    return false;  // error
}
CUSTOM_Vector materialEmitted(float u, float v, const CUSTOM_Vector *p, const CUSTOM_Material *mat)
{
  if (mat->type < CUSTOM_MATERIAL_MAX)
    return (CUSTOM_Vector){0.0f, 0.0f, 0.0f};
  else if (mat->type == CUSTOM_DIFFUSE_LIGHT)
    return CUSTOM_DiffuseLightEmitted(u, v, p, (CUSTOM_DiffuseLight *)mat);
}

/**/

bool CUSTOM_SphereHit(CUSTOM_HitRecord *rec,
                      const CUSTOM_Ray *ray,
                      const CUSTOM_Sphere *sphere,
                      float t_min,
                      float t_max)
{
  float co[3];  // sphere's center to ray's origin
  sub_v3_v3v3(co, ray->origin.vec3, sphere->center.vec3);

  float a = dot_v3v3(ray->direction.vec3, ray->direction.vec3);
  float b = dot_v3v3(co, ray->direction.vec3);
  float c = dot_v3v3(co, co) - sphere->radius * sphere->radius;
  float discriminant = b * b - a * c;
  if (discriminant > 0) {
    float temp = (-b - sqrt(discriminant)) / a;
    if (temp < t_max && temp > t_min) {
      rec->t = temp;
      rec->position = pointAtParameter(ray, rec->t);
      sub_v3_v3v3(rec->normal.vec3, rec->position.vec3, sphere->center.vec3);
      // mul_v3_fl(rec->normal.vec3, 1.0f / sphere->radius);
      normalize_v3(rec->normal.vec3);
      rec->mat_ptr = sphere->mat_ptr;
      getSphereUV(&rec->u, &rec->v, &rec->normal);

      return true;
    }
    temp = (-b + sqrt(discriminant)) / a;
    if (temp < t_max && temp > t_min) {
      rec->t = temp;
      rec->position = pointAtParameter(ray, rec->t);
      sub_v3_v3v3(rec->normal.vec3, rec->position.vec3, sphere->center.vec3);
      // mul_v3_fl(rec->normal.vec3, 1.0f / sphere->radius);
      normalize_v3(rec->normal.vec3);
      rec->mat_ptr = sphere->mat_ptr;
      getSphereUV(&rec->u, &rec->v, &rec->normal);

      return true;
    }
  }
  return false;
}
bool CUSTOM_RectXYHit(
    CUSTOM_HitRecord *rec, const CUSTOM_Ray *ray, const CUSTOM_RectXY *rect, float t0, float t1)
{
  float t = (rect->k - ray->origin.z) / ray->direction.z;
  if (t < t0 || t > t1)
    return false;
  float x = ray->origin.x + t * ray->direction.x;
  float y = ray->origin.y + t * ray->direction.y;
  if (x < rect->x0 || x > rect->x1 || y < rect->y0 || y > rect->y1)
    return false;
  rec->u = (x - rect->x0) / (rect->x1 - rect->x0);
  rec->v = (y - rect->y0) / (rect->y1 - rect->y0);
  rec->t = t;
  rec->mat_ptr = rect->mat_ptr;
  rec->position = pointAtParameter(ray, rec->t);
  rec->normal = rect->normal;
  return true;
}
bool CUSTOM_RectXZHit(
    CUSTOM_HitRecord *rec, const CUSTOM_Ray *ray, const CUSTOM_RectXZ *rect, float t0, float t1)
{
  float t = (rect->k - ray->origin.y) / ray->direction.y;
  if (t < t0 || t > t1)
    return false;
  float x = ray->origin.x + t * ray->direction.x;
  float z = ray->origin.z + t * ray->direction.z;
  if (x < rect->x0 || x > rect->x1 || z < rect->z0 || z > rect->z1)
    return false;
  rec->u = (x - rect->x0) / (rect->x1 - rect->x0);
  rec->v = (z - rect->z0) / (rect->z1 - rect->z0);
  rec->t = t;
  rec->mat_ptr = rect->mat_ptr;
  rec->position = pointAtParameter(ray, rec->t);
  rec->normal = rect->normal;
  return true;
}
bool CUSTOM_RectYZHit(
    CUSTOM_HitRecord *rec, const CUSTOM_Ray *ray, const CUSTOM_RectYZ *rect, float t0, float t1)
{
  float t = (rect->k - ray->origin.x) / ray->direction.x;
  if (t < t0 || t > t1)
    return false;
  float y = ray->origin.y + t * ray->direction.y;
  float z = ray->origin.z + t * ray->direction.z;
  if (y < rect->y0 || y > rect->y1 || z < rect->z0 || z > rect->z1)
    return false;
  rec->u = (y - rect->y0) / (rect->y1 - rect->y0);
  rec->v = (z - rect->z0) / (rect->z1 - rect->z0);
  rec->t = t;
  rec->mat_ptr = rect->mat_ptr;
  rec->position = pointAtParameter(ray, rec->t);
  rec->normal = rect->normal;
  return true;
}
bool CUSTOM_BoxHit(
    CUSTOM_HitRecord *rec, const CUSTOM_Ray *ray, const CUSTOM_Box *box, float t0, float t1)
{
  if (list_hit(rec, ray, &box->faces, t0, t1))
    return true;
  else
    return false;
  // return list_hit(rec, ray, &box->faces, t0, t1);
}
bool CUSTOM_TriangleHit(
    CUSTOM_HitRecord *rec, const CUSTOM_Ray *ray, const CUSTOM_Triangle *tri, float t0, float t1)
{
  /*Möller–Trumbore intersection algorithm*/
  /* https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm */

  const float EPSILON = 0.0000001;
  // Vector3D vertex0 = inTriangle->vertex0;
  // Vector3D vertex1 = inTriangle->vertex1;
  // Vector3D vertex2 = inTriangle->vertex2;
  // Vector3D edge1, edge2, h, s, q;
  float edge1[3], edge2[3], h[3], s[3], q[3];
  // edge1 = vertex1 - vertex0;
  // edge2 = vertex2 - vertex0;
  sub_v3_v3v3(edge1, tri->v1.vec3, tri->v0.vec3);
  sub_v3_v3v3(edge2, tri->v2.vec3, tri->v0.vec3);
  // h = rayVector.crossProduct(edge2);
  // a = edge1.dotProduct(h);
  cross_v3_v3v3(h, ray->direction.vec3, edge2);
  float a = dot_v3v3(edge1, h);
  if (a > -EPSILON && a < EPSILON)
    return false;  // This ray is parallel to this triangle.
  float f = 1.0 / a;
  // s = rayOrigin - vertex0;
  // u = f * s.dotProduct(h);
  sub_v3_v3v3(s, ray->origin.vec3, tri->v0.vec3);
  float u = f * dot_v3v3(s, h);
  if (u < 0.0 || u > 1.0)
    return false;
  // q = s.crossProduct(edge1);
  // v = f * rayVector.dotProduct(q);
  cross_v3_v3v3(q, s, edge1);
  float v = f * dot_v3v3(ray->direction.vec3, q);
  if (v < 0.0 || (double)u + (double)v > 1.0)
    return false;
  // At this stage we can compute t to find out where the intersection point is on the line.
  // float t = f * edge2.dotProduct(q);
  float t = f * dot_v3v3(edge2, q);

  if (t < t0 || t > t1) {
    return false;
  }
  // outIntersectionPoint = rayOrigin + rayVector * t;
  CUSTOM_Vector P = pointAtParameter(ray, t);

  /* https://codeplea.com/triangular-interpolation */
  float uv[2];
  //float w[3] = {
  //    1.0f / len_v3v3(P.vec3, tri->v0.vec3),
  //    1.0f / len_v3v3(P.vec3, tri->v1.vec3),
  //    1.0f / len_v3v3(P.vec3, tri->v2.vec3),
  //};
  //float w012 = w[0] + w[1] + w[2];
  //w[0] /= w012;
  //w[1] /= w012;
  //w[2] /= w012;
  float w[3];

  barycentric(&w[0], &w[1], &w[2], P.vec3, tri->v0.vec3, tri->v1.vec3, tri->v2.vec3);

  float normal[3];
  interp_v2_v2v2v2(uv, tri->uv0.vec2, tri->uv1.vec2, tri->uv2.vec2, w);
  interp_v3_v3v3v3(normal, tri->n0.vec3, tri->n1.vec3, tri->n2.vec3, w);

  normalize_v3(normal);

  rec->u = uv[0];
  rec->v = uv[1];
  rec->t = t;
  rec->mat_ptr = tri->mat_ptr;
  rec->position = P;
  rec->normal = (CUSTOM_Vector){normal[0], normal[1], normal[2]};

  return true;  // this ray hits the triangle
}
bool CUSTOM_BvhNodeHit(CUSTOM_HitRecord *rec,
                       const CUSTOM_Ray *ray,
                       const CUSTOM_BvhNode *node,
                       double t_min,
                       double t_max)
{
  if (!AABB_hit(&node->box, ray, t_min, t_max))
    return false;

  bool hit_left = CUSTOM_HittableHit(rec, ray, node->left, t_min, t_max);
  bool hit_right = CUSTOM_HittableHit(rec, ray, node->right, t_min, hit_left ? rec->t : t_max);

  return hit_left || hit_right;
}

/******************************************************************
 * Sampling Color by a ray
 * Return:
 *  color:  pointer of CUSTOM_Vector vec3 rgb
 * Input:
 *  ray:    pointer of const CUSTOM_Ray
 *  world:  pointer of const ListBase
 *  depth:  int, now recursive depth
 */
CUSTOM_Vector CUSTOM_sampleColor(const CUSTOM_Ray *ray, const ListBase *world, int depth)
{
  CUSTOM_HitRecord rec;
  if (list_hit(&rec, ray, world, 0.01f, FLT_MAX)) {
    // return (CUSTOM_Vector){1.0f, 0.0f, 0.0f};
    CUSTOM_Ray scattered;
    CUSTOM_Vector attenuation;

    CUSTOM_Vector emitted = materialEmitted(rec.u, rec.v, &rec.position, rec.mat_ptr);

    if (depth < CUSTOM_RECURSIVE_MAX &&
        materialScatter(&scattered, &attenuation, ray, &rec, rec.mat_ptr)) {
      CUSTOM_Vector color = CUSTOM_sampleColor(&scattered, world, depth + 1);
      mul_v3_v3(color.vec3, attenuation.vec3);
      add_v3_v3(color.vec3, emitted.vec3);
      return color;
    }
    else {
      return emitted;
    }
  }
  else {
    return (CUSTOM_Vector){0.0f, 0.0f, 0.0f};
  }
}
/******************************************************************
 * Initialize Camera
 * Return:
 *  camera:  CUSTOM_Camera
 * Input:
 *  lookfrom[3]:  float, camera's origin
 *  lookat[3]:    float, camera's focus point
 *  vup[3]:       float, view up
 *  hfov:         float, horizontal field of view(radians)
 *  aspect:       float, width / height
 *  aperture:     float,
 *  focus_dist:   float, focus distance
 */
CUSTOM_Camera *CUSTOM_CameraCreate(const float lookfrom[3],
                                   const float lookat[3],
                                   const float vup[3],
                                   const float hfov,
                                   const float aspect,
                                   const float aperture,
                                   const float focus_dist)
{
  float half_width = tan(hfov / 2);
  float half_height = half_width / aspect;

  CUSTOM_Vector u, v, w;
  {
    sub_v3_v3v3(w.vec3, lookfrom, lookat);
    normalize_v3(w.vec3);

    cross_v3_v3v3(u.vec3, vup, w.vec3);
    normalize_v3(u.vec3);

    cross_v3_v3v3(v.vec3, w.vec3, u.vec3);
  }

  CUSTOM_Vector lower_left_corner;
  {
    madd_v3_v3v3fl(lower_left_corner.vec3, lookfrom, u.vec3, -half_width * focus_dist);
    madd_v3_v3fl(lower_left_corner.vec3, v.vec3, -half_height * focus_dist);
    madd_v3_v3fl(lower_left_corner.vec3, w.vec3, -focus_dist);
  }
  CUSTOM_Vector horizontal, vertical;
  mul_v3_v3fl(horizontal.vec3, u.vec3, 2 * half_width * focus_dist);
  mul_v3_v3fl(vertical.vec3, v.vec3, 2 * half_height * focus_dist);

  CUSTOM_Camera *camera = (CUSTOM_Camera *)MEM_callocN(sizeof(CUSTOM_Camera), "CUSTOM_Camera");
  camera->origin = (CUSTOM_Vector){lookfrom[0], lookfrom[1], lookfrom[2]};
  camera->lower_left_corner = lower_left_corner;
  camera->horizontal = horizontal;
  camera->vertical = vertical;
  camera->u = u;
  camera->v = v;
  camera->w = w;
  camera->aspect = aspect;
  camera->lens_radius = aperture / 2.0f;

  return camera;
}

/******************************************************************
 * Get a ray from camera with (u, v)
 * Return   :pointer of CUSTOM_Camera, Camera
 * Input:
 *  u       :float
 *  v       :float
 */

CUSTOM_Ray CUSTOM_CameraGetRay(const CUSTOM_Camera *cam, float s, float t)
{
  // vec3 rd = lens_radius * random_in_unit_disk();
  float rd[3], offset[3];
  mul_v3_v3fl(rd, randomInUnitDisk().vec3, cam->lens_radius);
  // vec3 offset = u * rd.x() + v * rd.y();
  mul_v3_v3fl(offset, cam->u.vec3, rd[0]);
  madd_v3_v3fl(offset, cam->v.vec3, rd[1]);
  // return ray(origin + offset, lower_left_corner + s * horizontal + t * vertical - origin -
  // offset);
  CUSTOM_Vector ray_origin;
  add_v3_v3v3(ray_origin.vec3, cam->origin.vec3, offset);
  CUSTOM_Vector ray_dir;
  madd_v3_v3v3fl(ray_dir.vec3, cam->lower_left_corner.vec3, cam->horizontal.vec3, s);
  madd_v3_v3fl(ray_dir.vec3, cam->vertical.vec3, t);
  sub_v3_v3(ray_dir.vec3, ray_origin.vec3);
  return (CUSTOM_Ray){ray_origin, ray_dir};
}
/******************************************************************
 * Get a ray from camera with (u, v)
 * Return     :pointer of CUSTOM_Sphere, sphere
 * Input:
 *  center[3]     :float, center of sphere
 *  r             :float, radius of sphere
 *  material      :pointer of CUSTOM_Material, material of sphere
 */
CUSTOM_Sphere *CUSTOM_SphereCreate(const float center[3], float r, const CUSTOM_Material *material)
{
  CUSTOM_Sphere *sphere = (CUSTOM_Sphere *)MEM_callocN(sizeof(CUSTOM_Sphere), "CUSTOM_Sphere");

  sphere->type = CUSTOM_SPHERE;
  sphere->prev = NULL;
  sphere->next = NULL;

  sphere->center = (CUSTOM_Vector){center[0], center[1], center[2]};
  sphere->radius = r;
  sphere->mat_ptr = material;

  return sphere;
}
/**/
CUSTOM_RectXY *CUSTOM_RectXYCreate(
    float x0, float x1, float y0, float y1, float k, const CUSTOM_Material *mat, bool flip_normal)
{
  CUSTOM_RectXY *rect = (CUSTOM_RectXY *)MEM_callocN(sizeof(CUSTOM_RectXY), "CUSTOM_RectXY");
  rect->type = CUSTOM_RECT_XY;
  rect->prev = NULL;
  rect->next = NULL;

  rect->x0 = x0;
  rect->x1 = x1;
  rect->y0 = y0;
  rect->y1 = y1;
  rect->k = k;
  rect->mat_ptr = mat;
  rect->normal = flip_normal ? (CUSTOM_Vector){0.0f, 0.0f, -1.0f} :
                               (CUSTOM_Vector){0.0f, 0.0f, 1.0f};
  return rect;
}
/**/
CUSTOM_RectXZ *CUSTOM_RectXZCreate(
    float x0, float x1, float z0, float z1, float k, const CUSTOM_Material *mat, bool flip_normal)
{
  CUSTOM_RectXZ *rect = (CUSTOM_RectXZ *)MEM_callocN(sizeof(CUSTOM_RectXZ), "CUSTOM_RectXZ");
  rect->type = CUSTOM_RECT_XZ;
  rect->prev = NULL;
  rect->next = NULL;

  rect->x0 = x0;
  rect->x1 = x1;
  rect->z0 = z0;
  rect->z1 = z1;
  rect->k = k;
  rect->mat_ptr = mat;
  rect->normal = flip_normal ? (CUSTOM_Vector){0.0f, -1.0f, 0.0f} :
                               (CUSTOM_Vector){0.0f, 1.0f, 0.0f};
  return rect;
}
/**/
CUSTOM_RectYZ *CUSTOM_RectYZCreate(
    float y0, float y1, float z0, float z1, float k, const CUSTOM_Material *mat, bool flip_normal)
{
  CUSTOM_RectYZ *rect = (CUSTOM_RectYZ *)MEM_callocN(sizeof(CUSTOM_RectYZ), "CUSTOM_RectYZ");
  rect->type = CUSTOM_RECT_YZ;
  rect->prev = NULL;
  rect->next = NULL;

  rect->y0 = y0;
  rect->y1 = y1;
  rect->z0 = z0;
  rect->z1 = z1;
  rect->k = k;
  rect->mat_ptr = mat;
  rect->normal = flip_normal ? (CUSTOM_Vector){-1.0f, 0.0f, 0.0f} :
                               (CUSTOM_Vector){1.0f, 0.0f, 0.0f};
  return rect;
}

CUSTOM_Box *CUSTOM_BoxCreate(const float p_min[3],
                             const float p_max[3],
                             const CUSTOM_Material *mat)
{
  CUSTOM_Box *box = (CUSTOM_Box *)MEM_callocN(sizeof(CUSTOM_Box), "CUSTOM_Box");
  box->type = CUSTOM_BOX;
  box->prev = NULL;
  box->next = NULL;

  box->pos_min = (CUSTOM_Vector){p_min[0], p_min[1], p_min[2]};
  box->pos_max = (CUSTOM_Vector){p_max[0], p_max[1], p_max[2]};

  CUSTOM_RectYZ *front = CUSTOM_RectYZCreate(
      box->pos_min.y, box->pos_max.y, box->pos_min.z, box->pos_max.z, box->pos_max.x, mat, false);
  CUSTOM_RectYZ *back = CUSTOM_RectYZCreate(
      box->pos_min.y, box->pos_max.y, box->pos_min.z, box->pos_max.z, box->pos_min.x, mat, true);
  CUSTOM_RectXY *up = CUSTOM_RectXYCreate(
      box->pos_min.x, box->pos_max.x, box->pos_min.y, box->pos_max.y, box->pos_max.z, mat, false);
  CUSTOM_RectXY *down = CUSTOM_RectXYCreate(
      box->pos_min.x, box->pos_max.x, box->pos_min.y, box->pos_max.y, box->pos_min.z, mat, true);
  CUSTOM_RectXZ *right = CUSTOM_RectXZCreate(
      box->pos_min.x, box->pos_max.x, box->pos_min.z, box->pos_max.z, box->pos_max.y, mat, false);
  CUSTOM_RectXZ *left = CUSTOM_RectXZCreate(
      box->pos_min.x, box->pos_max.x, box->pos_min.z, box->pos_max.z, box->pos_min.y, mat, true);

  box->faces.first = box->faces.last = NULL;
  BLI_addtail(&box->faces, front);
  BLI_addtail(&box->faces, back);
  BLI_addtail(&box->faces, up);
  BLI_addtail(&box->faces, down);
  BLI_addtail(&box->faces, right);
  BLI_addtail(&box->faces, left);

  return box;
}
CUSTOM_Triangle *CUSTOM_TriangleCreate(const float v0[3],
                                       const float v1[3],
                                       const float v2[3],
                                       const float n0[3],
                                       const float n1[3],
                                       const float n2[3],
                                       const float uv0[3],
                                       const float uv1[3],
                                       const float uv2[3],
                                       const CUSTOM_Material *mat)
{
  CUSTOM_Triangle *tri = (CUSTOM_Triangle *)MEM_callocN(sizeof(CUSTOM_Triangle),
                                                        "CUSTOM_Triangle");
  tri->type = CUSTOM_TRIANGLE;
  tri->prev = NULL;
  tri->next = NULL;

  tri->v0 = (CUSTOM_Vector){v0[0], v0[1], v0[2]};
  tri->v1 = (CUSTOM_Vector){v1[0], v1[1], v1[2]};
  tri->v2 = (CUSTOM_Vector){v2[0], v2[1], v2[2]};
  tri->n0 = (CUSTOM_Vector){n0[0], n0[1], n0[2]};
  tri->n1 = (CUSTOM_Vector){n1[0], n1[1], n1[2]};
  tri->n2 = (CUSTOM_Vector){n2[0], n2[1], n2[2]};
  tri->uv0 = (CUSTOM_Vector){uv0[0], uv0[1], uv0[2]};
  tri->uv1 = (CUSTOM_Vector){uv1[0], uv1[1], uv1[2]};
  tri->uv2 = (CUSTOM_Vector){uv2[0], uv2[1], uv2[2]};
  tri->mat_ptr = mat;

  return tri;
}

CUSTOM_BvhNode *CUSTOM_BvhNodeCreate(const CUSTOM_Hittable **objects, size_t start, size_t end)
{
  CUSTOM_BvhNode *node = (CUSTOM_BvhNode *)MEM_callocN(sizeof(CUSTOM_BvhNode), "CUSTOM_BvhNode");

  int axis = randomInt(0, 2);
  int (*comparator)(const void *, const void *) = (axis == 0) ?
                                                      boxCompareX :
                                                      (axis == 1) ? boxCompareY : boxCompareZ;

  // int amount = BLI_listbase_count(objects);
  size_t object_span = end - start;

  if (object_span == 1) {
    node->left = node->right = objects[start];
  }
  else if (object_span == 2) {
    if (comparator(&objects[start], &objects[start + 1]) < 0) {
      node->left = objects[start];
      node->right = objects[start + 1];
    }
    else {
      node->left = objects[start + 1];
      node->right = objects[start];
    }
  }
  else {
    // BLI_listbase_sort(objects, comparator);
    qsort(objects, end - start, sizeof(CUSTOM_Hittable *), comparator);

    int mid = start + object_span / 2;

    node->left = (CUSTOM_Hittable *)CUSTOM_BvhNodeCreate(objects, start, mid);
    node->right = (CUSTOM_Hittable *)CUSTOM_BvhNodeCreate(objects, mid, end);
    // ListBase left_list = {NULL, NULL};
    // ListBase right_list = {NULL, NULL};
    // CUSTOM_Hittable *hittable = (CUSTOM_Hittable *)objects->first;
    // for (int i = 0; hittable; hittable = hittable->next, i++) {
    //  if (i < mid) {
    //    BLI_addtail(&left_list, hittable);
    //  }
    //  else {
    //    BLI_addtail(&right_list, hittable);
    //  }
    //}
    // left = make_shared<bvh_node>(objects, start, mid, time0, time1);
    // right = make_shared<bvh_node>(objects, mid, end, time0, time1);
    // node->left = (CUSTOM_Hittable*)CUSTOM_BvhNodeCreate(&left_list);
    // node->right = (CUSTOM_Hittable *)CUSTOM_BvhNodeCreate(&right_list);
  }

  CUSTOM_AABB box_left, box_right;

  // if (!left->bounding_box(time0, time1, box_left) || !right->bounding_box(time0, time1,
  // box_right))
  //  std::cerr << "No bounding box in bvh_node constructor.\n";
  if (!CUSTOM_HittableBoundingBox(&box_left, node->left) ||
      !CUSTOM_HittableBoundingBox(&box_right, node->right)) {
    puts("No bounding box in bvh_node constructor.");
  }

  node->box = surroundingBox(&box_left, &box_right);
  node->type = CUSTOM_BVH_NODE;
  node->prev = NULL;
  node->next = NULL;
  return node;
}

/******************************************************************
 * Create a image with 4 channel float rgba
 * Output:
 *  image[wdt*hgt]  :CUSTOM_Vector array
 * Input:
 *  wdt             :float, center of sphere
 *  hgt             :float, radius of sphere
 */
void CUSTOM_ImageCreate(CUSTOM_Vector **image, const int wdt, const int hgt)
{
  *image = (CUSTOM_Vector *)MEM_calloc_arrayN(wdt * hgt, sizeof(CUSTOM_Vector), "CUSTOM_Image");
}

/******************************************************************
 * Return:
 *
 * Output:
 *
 * Input:
 *
 */
CUSTOM_Material *CUSTOM_LambertianCreate(const float albedo[3], const CUSTOM_Texture *texture)
{
  CUSTOM_Lambertian *lamb = (CUSTOM_Lambertian *)MEM_callocN(sizeof(CUSTOM_Lambertian),
                                                             "CUSTOM_Lambertian");
  lamb->type = CUSTOM_LAMBERTIAN;
  lamb->albedo = (CUSTOM_Vector){albedo[0], albedo[1], albedo[2]};
  lamb->texture = texture;
  return (CUSTOM_Material *)lamb;
}

/******************************************************************
 * Return:
 *
 * Output:
 *
 * Input:
 *
 */
CUSTOM_Material *CUSTOM_MetalCreate(const float albedo[3], const float f)
{
  CUSTOM_Metal *metal = (CUSTOM_Metal *)MEM_callocN(sizeof(CUSTOM_Metal), "CUSTOM_Metal");
  metal->type = CUSTOM_METAL;
  metal->albedo = (CUSTOM_Vector){albedo[0], albedo[1], albedo[2]};
  metal->fuzzy = f < 1.0f ? f : 1.0f;
  metal->texture = NULL;
  return (CUSTOM_Material *)metal;
}

/******************************************************************
 * Return:
 *
 * Output:
 *
 * Input:
 *
 */
CUSTOM_Material *CUSTOM_DielectricCreate(const float ri)
{
  CUSTOM_Dielectric *diele = (CUSTOM_Dielectric *)MEM_callocN(sizeof(CUSTOM_Dielectric),
                                                              "CUSTOM_Dielectric");
  diele->type = CUSTOM_DIELECTRIC;
  diele->ref_idx = ri;
  diele->texture = NULL;
  return (CUSTOM_Material *)diele;
}

CUSTOM_Material *CUSTOM_DiffuseLightCreate(const CUSTOM_Texture *texture)
{
  CUSTOM_DiffuseLight *light = (CUSTOM_DiffuseLight *)MEM_callocN(sizeof(CUSTOM_DiffuseLight),
                                                                  "CUSTOM_DiffuseLight");
  light->type = CUSTOM_DIFFUSE_LIGHT;
  light->emit = texture;
  return (CUSTOM_Material *)light;
}

/******************************************************************
 * Return:
 *  true
 * Output:
 *  scattered     : CUSTOM_Ray*, scattered ray
 *  attenuation   : CUSTOM_Vector* vec3, attenuation ratio
 * Input:
 *  rec           : const CUSTOM_HitRecord*, record
 *  lamb          : const CUSTOM_Lambertian*, be scattered lambertian
 */
bool CUSTOM_LambertianScatter(CUSTOM_Ray *scattered,
                              CUSTOM_Vector *attenuation,
                              const CUSTOM_HitRecord *rec,
                              const CUSTOM_Lambertian *lamb)
{
  float target[3];
  CUSTOM_Vector random = randomInUnitSphere();
  add_v3_v3v3(target, rec->position.vec3, rec->normal.vec3);
  add_v3_v3(target, random.vec3);
  // CUSTOM_Vector random = randomInUnitHemisphere(&rec->normal);
  // add_v3_v3v3(target, rec->position.vec3, random.vec3);

  copy_v3_v3(scattered->origin.vec3, rec->position.vec3);
  sub_v3_v3v3(scattered->direction.vec3, target, rec->position.vec3);

  if (lamb->texture)
    *attenuation = textureValue(lamb->texture, rec->u, rec->v, &rec->position);
  else
    *attenuation = lamb->albedo;
  return true;
}
/******************************************************************
 * Return:
 *  scatter result: bool, between scattered ray and normal < 90 angle = true, else false
 * Output:
 *  scattered     : CUSTOM_Ray*, scattered ray
 *  attenuation   : CUSTOM_Vector* vec3, attenuation ratio
 * Input:
 *  ray_in        " const CUSTOM_Ray*, inserted ray
 *  rec           : const CUSTOM_HitRecord*, record
 *  lamb          : const CUSTOM_Lambertian*, be scattered lambertian
 */
bool CUSTOM_MetalScatter(CUSTOM_Ray *scattered,
                         CUSTOM_Vector *attenuation,
                         const CUSTOM_Ray *ray_in,
                         const CUSTOM_HitRecord *rec,
                         const CUSTOM_Metal *metal)
{
  float reflected[3], in_unit[3];
  normalize_v3_v3(in_unit, ray_in->direction.vec3);
  reflect_v3_v3v3(reflected, in_unit, rec->normal.vec3);
  CUSTOM_Vector random = randomInUnitSphere();

  copy_v3_v3(scattered->origin.vec3, rec->position.vec3);
  madd_v3_v3v3fl(scattered->direction.vec3, reflected, random.vec3, metal->fuzzy);

  *attenuation = metal->albedo;

  return (dot_v3v3(scattered->direction.vec3, rec->normal.vec3) > 0.0f);
}
/******************************************************************
 * Return:
 *  true
 * Output:
 *  scattered     : CUSTOM_Ray*, scattered ray
 *  attenuation   : CUSTOM_Vector* vec3, attenuation ratio
 * Input:
 *  rec           : const CUSTOM_HitRecord*, record
 *  lamb          : const CUSTOM_Lambertian*, be scattered lambertian
 */
bool CUSTOM_DielectricScatter(CUSTOM_Ray *scattered,
                              CUSTOM_Vector *attenuation,
                              const CUSTOM_Ray *ray_in,
                              const CUSTOM_HitRecord *rec,
                              const CUSTOM_Dielectric *diele)
{
  float outward_normal[3], reflected[3], refracted[3];
  // vec3 reflected = reflect(r_in.direction(), rec.normal);
  reflect_v3_v3v3(reflected, ray_in->direction.vec3, rec->normal.vec3);
  *attenuation = (CUSTOM_Vector){1.0f, 1.0f, 1.0f};

  float reflect_prob, cosine, ni_over_nt;
  float ray_dot_normal = dot_v3v3(ray_in->direction.vec3, rec->normal.vec3);
  if (ray_dot_normal > 0) {  // insert
    // outward_normal = -rec.normal;
    mul_v3_v3fl(outward_normal, rec->normal.vec3, -1.0f);
    ni_over_nt = diele->ref_idx;
    cosine = diele->ref_idx * ray_dot_normal / len_v3(ray_in->direction.vec3);
  }
  else {  // out
    // outward_normal = rec.normal;
    copy_v3_v3(outward_normal, rec->normal.vec3);
    ni_over_nt = 1.0f / diele->ref_idx;
    cosine = -ray_dot_normal / len_v3(ray_in->direction.vec3);
  }

  if (refract(&refracted, ray_in->direction.vec3, outward_normal, ni_over_nt)) {
    reflect_prob = schlick(cosine, diele->ref_idx);
  }
  else {
    reflect_prob = 1.0;
  }

  if (randomFloat() < reflect_prob) {
    // scattered = ray(rec->position, reflected);
    copy_v3_v3(scattered->origin.vec3, rec->position.vec3);
    copy_v3_v3(scattered->direction.vec3, reflected);
  }
  else {
    // scattered = ray(rec->position, refracted);
    copy_v3_v3(scattered->origin.vec3, rec->position.vec3);
    copy_v3_v3(scattered->direction.vec3, refracted);
  }
  return true;
}
/**/
bool CUSTOM_DiffuseLightScatter(const CUSTOM_DiffuseLight *light)
{
  return false;
}

CUSTOM_Vector CUSTOM_DiffuseLightEmitted(float u,
                                         float v,
                                         const CUSTOM_Vector *p,
                                         const CUSTOM_DiffuseLight *light)
{
  return textureValue(light->emit, u, v, p);
}

/**/
CUSTOM_Texture *CUSTOM_ConstantTextureCreate(const float color[3])
{
  CUSTOM_ConstantTexture *texture = (CUSTOM_ConstantTexture *)MEM_callocN(
      sizeof(CUSTOM_ConstantTexture), "CUSTOM_ConstantTexture");

  texture->type = CUSTOM_CONSTANT_TEXTURE;
  texture->color = (CUSTOM_Vector){color[0], color[1], color[2]};

  return (CUSTOM_Texture *)texture;
}
/**/
CUSTOM_Texture *CUSTOM_CheckerTextureCreate(const CUSTOM_Texture *odd, const CUSTOM_Texture *even)
{
  CUSTOM_CheckerTexture *texture = (CUSTOM_CheckerTexture *)MEM_callocN(
      sizeof(CUSTOM_CheckerTexture), "CUSTOM_CheckerTexture");

  texture->type = CUSTOM_CHECKER_TEXTURE;
  texture->odd = odd;
  texture->even = even;

  return (CUSTOM_Texture *)texture;
}

CUSTOM_Texture *CUSTOM_ImageTextureCreate(const float *rgba, int w, int h)
{
  CUSTOM_ImageTexture *texture = (CUSTOM_ImageTexture *)MEM_callocN(sizeof(CUSTOM_ImageTexture),
                                                                    "CUSTOM_ImageTexture");

  texture->type = CUSTOM_IMAGE_TEXTURE;
  texture->img = (CUSTOM_Vector *)rgba;
  texture->wdt = w;
  texture->hgt = h;

  return (CUSTOM_Texture *)texture;
}

/**/
CUSTOM_Vector CUSTOM_ConstantTextureValue(const CUSTOM_ConstantTexture *texture)
{
  return texture->color;
}

/**/
CUSTOM_Vector CUSTOM_CheckerTextureValue(const CUSTOM_CheckerTexture *texture,
                                         float u,
                                         float v,
                                         const CUSTOM_Vector *pos)
{
  float sines = sin(10.0f * pos->x) * sin(10.0f * pos->y) * sin(10.0f * pos->z);
  if (sines < 0)
    return textureValue(texture->odd, u, v, pos);
  else
    return textureValue(texture->even, u, v, pos);
}

/**/
CUSTOM_Vector CUSTOM_ImageTextureValue(const CUSTOM_ImageTexture *texture,
                                       float u,
                                       float v,
                                       const CUSTOM_Vector *pos)
{
  int w = texture->wdt;
  int h = texture->hgt;

  int i = u * w;
  // int j = (1 - v) * texture->hgt;
  int j = v * h;

  i = clamp_i(i, 0, w - 1);
  j = clamp_i(j, 0, h - 1);
  // i = max(0, min(i, w - 1));
  // j = max(0, min(j, h - 1));

  float r = texture->img[i + w * j].r;
  float g = texture->img[i + w * j].g;
  float b = texture->img[i + w * j].b;
  return (CUSTOM_Vector){r, g, b};
}

bool CUSTOM_BoxBoundingBox(CUSTOM_AABB *output_box, const CUSTOM_Box *box)
{
  output_box->min = box->pos_min;
  output_box->max = box->pos_max;
  return true;
}
bool CUSTOM_SphereBoundingBox(CUSTOM_AABB *output_box, const CUSTOM_Sphere *sphere)
{
  sub_v3_v3v3(output_box->min.vec3,
              sphere->center.vec3,
              (float[3]){sphere->radius, sphere->radius, sphere->radius});
  add_v3_v3v3(output_box->max.vec3,
              sphere->center.vec3,
              (float[3]){sphere->radius, sphere->radius, sphere->radius});

  return true;
}
bool CUSTOM_BvhNodeBoundingBox(CUSTOM_AABB *output_box, const CUSTOM_BvhNode *node)
{
  *output_box = node->box;
  return true;
}
bool CUSTOM_TriangleBoundingBox(CUSTOM_AABB *output_box, const CUSTOM_Triangle *tri)
{
  output_box->min.x = min(tri->v0.x, min(tri->v1.x, tri->v2.x));
  output_box->min.y = min(tri->v0.y, min(tri->v1.y, tri->v2.y));
  output_box->min.z = min(tri->v0.z, min(tri->v1.z, tri->v2.z));

  output_box->max.x = max(tri->v0.x, max(tri->v1.x, tri->v2.x));
  output_box->max.y = max(tri->v0.y, max(tri->v1.y, tri->v2.y));
  output_box->max.z = max(tri->v0.z, max(tri->v1.z, tri->v2.z));

  return true;
}
