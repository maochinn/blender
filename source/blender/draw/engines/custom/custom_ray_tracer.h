#pragma once

#define CUSTOM_RECURSIVE_MAX 10
#define AO_RAY_AMOUNT 10
#define LIGHT_RAY_AMOUNT 10

// vector 4
typedef union CUSTOM_Vector {
  struct {
    float x, y, z, w;
  };
  struct {
    float s, t, p, q;
  };
  struct {
    float r, g, b, a;
  };
  float vec2[2];
  float vec3[3];
  float vec4[4];

  // CUSTOM_Vector(float x = 0, float y = 0, float z = 0, float w = 0) : x(x), y(y), z(z), w(w){}
} CUSTOM_Vector;

/*Texture type*/
enum {
  CUSTOM_CONSTANT_TEXTURE = 0,
  CUSTOM_CHECKER_TEXTURE,
  CUSTOM_IMAGE_TEXTURE,
  CUSTOM_TEXTURE_MAX
};

typedef struct CUSTOM_Texture {
  short type;
} CUSTOM_Texture;
typedef struct CUSTOM_ConstantTexture {
  short type;
  CUSTOM_Vector color;
} CUSTOM_ConstantTexture;
typedef struct CUSTOM_CheckerTexture {
  short type;
  CUSTOM_Texture *odd;
  CUSTOM_Texture *even;
} CUSTOM_CheckerTexture;
typedef struct CUSTOM_ImageTexture {
  short type;

  CUSTOM_Vector *img;
  int wdt;
  int hgt;
} CUSTOM_ImageTexture;

/*Material Type*/
enum {
  CUSTOM_LAMBERTIAN = 0,
  CUSTOM_METAL,
  CUSTOM_DIELECTRIC,
  CUSTOM_MATERIAL_MAX = 10,
  CUSTOM_DIFFUSE_LIGHT,
  CUSTOM_LIGHT_MAX
};
typedef struct CUSTOM_Material {
  short type;
  CUSTOM_Texture *texture;
} CUSTOM_Material;
typedef struct CUSTOM_Lambertian {
  short type;
  CUSTOM_Texture *texture;

  CUSTOM_Vector albedo;  // vec3
} CUSTOM_Lambertian;
typedef struct CUSTOM_Metal {
  short type;
  CUSTOM_Texture *texture;

  CUSTOM_Vector albedo;  // vec3
  float fuzzy;
} CUSTOM_Metal;
typedef struct CUSTOM_Dielectric {
  short type;
  CUSTOM_Texture *texture;

  float ref_idx;  // Refractive index
} CUSTOM_Dielectric;
typedef struct CUSTOM_DiffuseLight {
  short type;
  CUSTOM_Texture *emit;

} CUSTOM_DiffuseLight;

/*Hittable type*/
enum {
  CUSTOM_SPHERE = 0,
  CUSTOM_RECT_XY,
  CUSTOM_RECT_XZ,
  CUSTOM_RECT_YZ,
  CUSTOM_BOX,
  CUSTOM_TRIANGLE,
  CUSTOM_BVH_NODE,
  CUSTOM_HITTABLE_MAX
};

typedef struct CUSTOM_AABB {
  CUSTOM_Vector min, max;
} CUSTOM_AABB;

typedef struct CUSTOM_Hittable {
  struct CUSTOM_Hittable *next, *prev;
  short type;
} CUSTOM_Hittable;

typedef struct CUSTOM_BvhNode {
  struct CUSTOM_Hittable *next, *prev;
  short type;

  struct CUSTOM_Hittable *left, *right;
  CUSTOM_AABB box;
} CUSTOM_BvhNode;

typedef struct CUSTOM_Sphere {
  struct CUSTOM_Hittable *next, *prev;
  short type;

  CUSTOM_Material *mat_ptr;
  CUSTOM_Vector center;  // vec3
  float radius;
} CUSTOM_Sphere;
typedef struct CUSTOM_RectXY {
  struct CUSTOM_Hittable *next, *prev;
  short type;

  CUSTOM_Material *mat_ptr;
  CUSTOM_Vector normal;
  float x0, x1, y0, y1, k;
} CUSTOM_RectXY;
typedef struct CUSTOM_RectXZ {
  struct CUSTOM_Hittable *next, *prev;
  short type;

  CUSTOM_Material *mat_ptr;
  CUSTOM_Vector normal;
  float x0, x1, z0, z1, k;
} CUSTOM_RectXZ;
typedef struct CUSTOM_RectYZ {
  struct CUSTOM_Hittable *next, *prev;
  short type;

  CUSTOM_Material *mat_ptr;
  CUSTOM_Vector normal;
  float y0, y1, z0, z1, k;
} CUSTOM_RectYZ;
typedef struct CUSTOM_Box {
  struct CUSTOM_Hittable *next, *prev;
  short type;

  CUSTOM_Vector pos_min, pos_max;
  // CUSTOM_RectXY *up, *down;
  // CUSTOM_RectXZ *right, *left;
  // CUSTOM_RectYZ *front, *back;
  ListBase faces;
} CUSTOM_Box;
typedef struct CUSTOM_Triangle {
  struct CUSTOM_Hittable *next, *prev;
  short type;

  CUSTOM_Material *mat_ptr;
  CUSTOM_Vector v0, v1, v2;
  CUSTOM_Vector n0, n1, n2;
  CUSTOM_Vector uv0, uv1, uv2;
} CUSTOM_Triangle;

typedef struct CUSTOM_HitRecord {
  float t;
  CUSTOM_Vector position;  // vec3
  CUSTOM_Vector normal;    // vec3
  float u, v;

  CUSTOM_Material *mat_ptr;
} CUSTOM_HitRecord;

typedef struct CUSTOM_Ray {
  CUSTOM_Vector origin;     // vec3
  CUSTOM_Vector direction;  // vec3
} CUSTOM_Ray;

typedef struct CUSTOM_Camera {
  CUSTOM_Vector origin;             // vec3
  CUSTOM_Vector lower_left_corner;  // vec3
  CUSTOM_Vector horizontal;         // vec3
  CUSTOM_Vector vertical;           // vec3
  CUSTOM_Vector u, v, w;            // vec3
  float aspect;
  float lens_radius;
} CUSTOM_Camera;

CUSTOM_Vector CUSTOM_sampleColor(const CUSTOM_Ray *ray, const ListBase *world, int depth);

CUSTOM_Camera *CUSTOM_CameraCreate(const float lookfrom[3],
                                   const float lookat[3],
                                   const float vup[3],
                                   const float hfov,
                                   const float aspect,
                                   const float aperture,
                                   const float focus_dist);
CUSTOM_Ray CUSTOM_CameraGetRay(const CUSTOM_Camera *cam, float s, float t);

CUSTOM_Sphere *CUSTOM_SphereCreate(const float center[3],
                                   float r,
                                   const CUSTOM_Material *material);
CUSTOM_RectXY *CUSTOM_RectXYCreate(
    float x0, float x1, float y0, float y1, float k, const CUSTOM_Material *mat, bool flip_normal);
CUSTOM_RectXZ *CUSTOM_RectXZCreate(
    float x0, float x1, float z0, float z1, float k, const CUSTOM_Material *mat, bool flip_normal);
CUSTOM_RectYZ *CUSTOM_RectYZCreate(
    float y0, float y1, float z0, float z1, float k, const CUSTOM_Material *mat, bool flip_normal);
CUSTOM_Box *CUSTOM_BoxCreate(const float p_min[3],
                             const float p_max[3],
                             const CUSTOM_Material *mat);
CUSTOM_Triangle *CUSTOM_TriangleCreate(const float v0[3],
                                       const float v1[3],
                                       const float v2[3],
                                       const float n0[3],
                                       const float n1[3],
                                       const float n2[3],
                                       const float uv0[3],
                                       const float uv1[3],
                                       const float uv2[3],
                                       const CUSTOM_Material *mat);
CUSTOM_BvhNode *CUSTOM_BvhNodeCreate(const CUSTOM_Hittable **objects, size_t start, size_t end);

CUSTOM_Material *CUSTOM_LambertianCreate(const float albedo[3], const CUSTOM_Texture *texture);
CUSTOM_Material *CUSTOM_MetalCreate(const float albedo[3], const float f);
CUSTOM_Material *CUSTOM_DielectricCreate(const float ri);
CUSTOM_Material *CUSTOM_DiffuseLightCreate(const CUSTOM_Texture *texture);

bool CUSTOM_HittableHit(CUSTOM_HitRecord *rec,
                        const CUSTOM_Ray *ray,
                        const CUSTOM_Hittable *hittable,
                        float t0,
                        float t1);
bool CUSTOM_SphereHit(CUSTOM_HitRecord *rec,
                      const CUSTOM_Ray *ray,
                      const CUSTOM_Sphere *sphere,
                      float t_min,
                      float t_max);
bool CUSTOM_RectXYHit(
    CUSTOM_HitRecord *rec, const CUSTOM_Ray *ray, const CUSTOM_RectXY *rect, float t0, float t1);
bool CUSTOM_RectXZHit(
    CUSTOM_HitRecord *rec, const CUSTOM_Ray *ray, const CUSTOM_RectXZ *rect, float t0, float t1);
bool CUSTOM_RectYZHit(
    CUSTOM_HitRecord *rec, const CUSTOM_Ray *ray, const CUSTOM_RectYZ *rect, float t0, float t1);
bool CUSTOM_BoxHit(
    CUSTOM_HitRecord *rec, const CUSTOM_Ray *ray, const CUSTOM_Box *box, float t0, float t1);
bool CUSTOM_TriangleHit(
    CUSTOM_HitRecord *rec, const CUSTOM_Ray *ray, const CUSTOM_Triangle *tri, float t0, float t1);
bool CUSTOM_BvhNodeHit(CUSTOM_HitRecord *rec,
                       const CUSTOM_Ray *ray,
                       const CUSTOM_BvhNode *node,
                       double t_min,
                       double t_max);

bool CUSTOM_HittableBoundingBox(CUSTOM_AABB *output_box, const CUSTOM_Hittable *hittable);
bool CUSTOM_SphereBoundingBox(CUSTOM_AABB *output_box, const CUSTOM_Sphere *sphere);
bool CUSTOM_BoxBoundingBox(CUSTOM_AABB *output_box, const CUSTOM_Box *box);
bool CUSTOM_TriangleBoundingBox(CUSTOM_AABB *output_box, const CUSTOM_Triangle *tri);
bool CUSTOM_BvhNodeBoundingBox(CUSTOM_AABB *output_box, const CUSTOM_BvhNode *node);

void CUSTOM_ImageCreate(CUSTOM_Vector **image, const int wdt, const int hgt);

bool CUSTOM_LambertianScatter(CUSTOM_Ray *scattered,
                              CUSTOM_Vector *attenuation,
                              const CUSTOM_HitRecord *rec,
                              const CUSTOM_Lambertian *lamb);
bool CUSTOM_MetalScatter(CUSTOM_Ray *scattered,
                         CUSTOM_Vector *attenuation,
                         const CUSTOM_Ray *ray_in,
                         const CUSTOM_HitRecord *rec,
                         const CUSTOM_Metal *metal);
bool CUSTOM_DielectricScatter(CUSTOM_Ray *scattered,
                              CUSTOM_Vector *attenuation,
                              const CUSTOM_Ray *ray_in,
                              const CUSTOM_HitRecord *rec,
                              const CUSTOM_Dielectric *diele);
bool CUSTOM_DiffuseLightScatter(const CUSTOM_DiffuseLight *light);

CUSTOM_Vector CUSTOM_DiffuseLightEmitted(float u,
                                         float v,
                                         const CUSTOM_Vector *p,
                                         const CUSTOM_DiffuseLight *light);

CUSTOM_Texture *CUSTOM_ConstantTextureCreate(const float color[3]);
CUSTOM_Texture *CUSTOM_CheckerTextureCreate(const CUSTOM_Texture *odd, const CUSTOM_Texture *even);
CUSTOM_Texture *CUSTOM_ImageTextureCreate(const float *rgba, int w, int h);

CUSTOM_Vector CUSTOM_ConstantTextureValue(const CUSTOM_ConstantTexture *texture);
CUSTOM_Vector CUSTOM_CheckerTextureValue(const CUSTOM_CheckerTexture *texture,
                                         float u,
                                         float v,
                                         const CUSTOM_Vector *pos);
CUSTOM_Vector CUSTOM_ImageTextureValue(const CUSTOM_ImageTexture *texture,
                                       float u,
                                       float v,
                                       const CUSTOM_Vector *pos);

inline float randomFloat()
{
  return rand() / (RAND_MAX + 1.0);
}
inline int randomInt(int min, int max)
{
  return min + rand() % (max - min + 1);
}
inline int boxCompare(const CUSTOM_Hittable *a, const CUSTOM_Hittable *b, int axis)
{
  CUSTOM_AABB box_a;
  CUSTOM_AABB box_b;

  if (!CUSTOM_HittableBoundingBox(&box_a, a, 0.0, 0.0) ||
      !CUSTOM_HittableBoundingBox(&box_b, b, 0.0, 0.0)) {
    puts("No bounding box in bvh_node constructor.");
  }

  return box_a.min.vec3[axis] < box_b.min.vec3[axis] ? -1 : 1;
}
