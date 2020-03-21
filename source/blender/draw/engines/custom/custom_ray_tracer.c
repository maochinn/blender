#pragma once

#include "DRW_render.h"

#include "BLI_math.h"
#include "DNA_listBase.h"

#include "custom_ray_tracer.h"

/* local function */

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
    //refracted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
    madd_v3_v3v3fl(refracted, uv, n, -dt);
    mul_v3_fl(refracted, ni_over_nt);
    madd_v3_v3fl(refracted, n, -sqrt(discriminant));
    return true;
  }
  else
    return false;
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
bool Sphere_hit(const CUSTOM_Ray *ray,
                const CUSTOM_Sphere *sphere,
                float t_min,
                float t_max,
                CUSTOM_HitRecord *rec)
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
      //mul_v3_fl(rec->normal.vec3, 1.0f / sphere->radius);
      normalize_v3(rec->normal.vec3);
      rec->mat_ptr = sphere->mat_ptr;
      return true;
    }
    temp = (-b + sqrt(discriminant)) / a;
    if (temp < t_max && temp > t_min) {
      rec->t = temp;
      rec->position = pointAtParameter(ray, rec->t);
      sub_v3_v3v3(rec->normal.vec3, rec->position.vec3, sphere->center.vec3);
      //mul_v3_fl(rec->normal.vec3, 1.0f / sphere->radius);
      normalize_v3(rec->normal.vec3);
      rec->mat_ptr = sphere->mat_ptr;
      return true;
    }
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

CUSTOM_Vector randomInUnitDisk(void)
{
  CUSTOM_Vector p;
  do {
    //p = 2.0 * vec3(random_double(), random_double(), 0) - vec3(1, 1, 0);
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
bool world_hit(
    const CUSTOM_Ray *ray, const ListBase *world, float t_min, float t_max, CUSTOM_HitRecord *rec)
{
  CUSTOM_HitRecord temp_rec;
  bool hit_anything = false;
  double closest_so_far = t_max;
  for (CUSTOM_Hittable *hittable = (CUSTOM_Hittable *)world->first; hittable;
       hittable = hittable->next) {
    if (hittable->type == CUSTOM_SPHERE) {
      CUSTOM_Sphere *sphere = (CUSTOM_Sphere *)hittable;
      if (Sphere_hit(ray, sphere, t_min, closest_so_far, &temp_rec)) {
        hit_anything = true;
        closest_so_far = temp_rec.t;
        *rec = temp_rec;
      }
    }
    else
      return false;
  }
  return hit_anything;
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
  if (world_hit(ray, world, 0.0, FLT_MAX, &rec)) {
    CUSTOM_Ray scattered;
    CUSTOM_Vector attenuation;

    if (depth < CUSTOM_RECURSIVE_MAX) {
      bool scatter_result = false;
      if (rec.mat_ptr->type == CUSTOM_LAMBERTIAN)
        scatter_result = CUSTOM_LambertianScatter(
            &scattered, &attenuation, &rec, (CUSTOM_Lambertian *)rec.mat_ptr);
      else if (rec.mat_ptr->type == CUSTOM_METAL)
        scatter_result = CUSTOM_MetalScatter(
            &scattered, &attenuation, ray, &rec, (CUSTOM_Metal *)rec.mat_ptr);
      else if (rec.mat_ptr->type == CUSTOM_DIELECTRIC)
        scatter_result = CUSTOM_DielectricScatter(
            &scattered, &attenuation, ray, &rec, (CUSTOM_Dielectric *)rec.mat_ptr);
      else
        return (CUSTOM_Vector){0.0f, 0.0f, 0.0f};  // error

      if (scatter_result) {
        CUSTOM_Vector color = CUSTOM_sampleColor(&scattered, world, depth + 1);
        mul_v3_v3(color.vec3, attenuation.vec3);
        return color;
      }
      else
        return (CUSTOM_Vector){0.0f, 0.0f, 0.0f};

    }
    else {
      return (CUSTOM_Vector){0.0f, 0.0f, 0.0f};
    }

    // float target[3];
    //// vec3 target = rec.p + rec.normal + random_in_unit_sphere();
    // CUSTOM_Vector random = randomInUnitSphere();
    // add_v3_v3v3(target, rec.position.vec3, rec.normal.vec3);
    // add_v3_v3(target, random.vec3);

    // CUSTOM_Ray new_ray;
    // copy_v3_v3(new_ray.origin.vec3, rec.position.vec3);
    // sub_v3_v3v3(new_ray.direction.vec3, target, rec.position.vec3);

    //// return 0.5 * CUSTOM_sampleColor(color, &new_ray, world);
    // CUSTOM_Vector color = CUSTOM_sampleColor(&new_ray, world);
    // mul_v3_fl(color.vec3, 0.5f);
    // return color;

    // return (CUSTOM_Vector){1.0f, 0.0f, 0.0f};

    // 0.5 * vec3(N.x() + 1, N.y() + 1, N.z() + 1)
    // return (CUSTOM_Vector){
    //    0.5f * (rec.normal.x + 1.0f),
    //    0.5f * (rec.normal.y + 1.0f),
    //    0.5f * (rec.normal.z + 1.0f),
    //    1.0f};
  }
  else {
    float unit_direction[3];
    normalize_v3_v3(unit_direction, ray->direction.vec3);
    float t = 0.5f * (unit_direction[2] + 1.0f);

    CUSTOM_Vector background;
    mul_v3_v3fl(background.vec3, (float[3]){1.0f, 1.0f, 1.0f}, 1.0f - t);
    madd_v3_v3fl(background.vec3, (float[3]){0.5f, 0.7f, 1.0f}, t);

    // mul_v3_v3(color->vec3, background.vec3);
    return background;
    // return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
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
  //vec3 rd = lens_radius * random_in_unit_disk();
  float rd[3], offset[3];
  mul_v3_v3fl(rd, randomInUnitDisk().vec3, cam->lens_radius);
  //vec3 offset = u * rd.x() + v * rd.y();
  mul_v3_v3fl(offset, cam->u.vec3, rd[0]);
  madd_v3_v3fl(offset, cam->v.vec3, rd[1]);
  //return ray(origin + offset, lower_left_corner + s * horizontal + t * vertical - origin - offset);
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

CUSTOM_Sphere *CUSTOM_SphereCreate(const float center[3], float r, CUSTOM_Material *material)
{
  CUSTOM_Sphere *sphere = (CUSTOM_Sphere *)MEM_callocN(sizeof(CUSTOM_Sphere), "CUSTOM_Sphere");
  sphere->type = CUSTOM_SPHERE;
  sphere->center = (CUSTOM_Vector){center[0], center[1], center[2]};
  sphere->radius = r;
  sphere->mat_ptr = material;

  return sphere;
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
CUSTOM_Lambertian *CUSTOM_LambertianCreate(const float albedo[3])
{
  CUSTOM_Lambertian *lamb = (CUSTOM_Lambertian *)MEM_callocN(sizeof(CUSTOM_Lambertian),
                                                             "CUSTOM_Lambertian");
  lamb->type = CUSTOM_LAMBERTIAN;
  lamb->albedo = (CUSTOM_Vector){albedo[0], albedo[1], albedo[2]};
  return lamb;
}
/******************************************************************
 * Return:
 *
 * Output:
 *
 * Input:
 *
 */
CUSTOM_Metal *CUSTOM_MetalCreate(const float albedo[3], const float f)
{
  CUSTOM_Metal *metal = (CUSTOM_Metal *)MEM_callocN(sizeof(CUSTOM_Metal), "CUSTOM_Metal");
  metal->type = CUSTOM_METAL;
  metal->albedo = (CUSTOM_Vector){albedo[0], albedo[1], albedo[2]};
  metal->fuzzy = f < 1.0f ? f : 1.0f;
  return metal;
}
/******************************************************************
 * Return:
 *
 * Output:
 *
 * Input:
 *
 */
CUSTOM_Dielectric *CUSTOM_DielectricCreate(const float ri)
{
  CUSTOM_Dielectric *diele = (CUSTOM_Dielectric *)MEM_callocN(sizeof(CUSTOM_Dielectric),
                                                              "CUSTOM_Dielectric");
  diele->type = CUSTOM_DIELECTRIC;
  diele->ref_idx = ri;
  return diele;
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

  copy_v3_v3(scattered->origin.vec3, rec->position.vec3);
  sub_v3_v3v3(scattered->direction.vec3, target, rec->position.vec3);

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
  if (ray_dot_normal > 0) { //insert
    // outward_normal = -rec.normal;
    mul_v3_v3fl(outward_normal, rec->normal.vec3, -1.0f);
    ni_over_nt = diele->ref_idx;
    cosine = diele->ref_idx * ray_dot_normal /
             len_v3(ray_in->direction.vec3);
  }
  else {  //out
    //outward_normal = rec.normal;
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
    //scattered = ray(rec->position, reflected);
    copy_v3_v3(scattered->origin.vec3, rec->position.vec3);
    copy_v3_v3(scattered->direction.vec3, reflected);
  }
  else {
    //scattered = ray(rec->position, refracted);
    copy_v3_v3(scattered->origin.vec3, rec->position.vec3);
    copy_v3_v3(scattered->direction.vec3, refracted);
  }
  return true;
}
