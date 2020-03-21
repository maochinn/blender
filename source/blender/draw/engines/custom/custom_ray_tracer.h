#pragma once

#define CUSTOM_RECURSIVE_MAX 10

//vector 4
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

/*Material Type*/
enum { CUSTOM_LAMBERTIAN = 0, CUSTOM_METAL, CUSTOM_DIELECTRIC, CUSTOM_MATERIAL_MAX };
typedef struct CUSTOM_Material {
  short type;
} CUSTOM_Material;
typedef struct CUSTOM_Lambertian {
  short type;
  CUSTOM_Vector albedo; //vec3
} CUSTOM_Lambertian;
typedef struct CUSTOM_Metal {
  short type;
  CUSTOM_Vector albedo;  // vec3
  float fuzzy;
} CUSTOM_Metal;
typedef struct CUSTOM_Dielectric {
  short type;
  float ref_idx;  //Refractive index
} CUSTOM_Dielectric;

/*Hittable type*/
enum { CUSTOM_SPHERE = 0, CUSTOM_HITTABLE_MAX };

typedef struct CUSTOM_Hittable {
  struct CUSTOM_Hittable *next, *prev;
  short type;
} CUSTOM_Hittable;
typedef struct CUSTOM_Sphere {
  struct CUSTOM_Hittable *next, *prev;
  short type;
  CUSTOM_Vector center; //vec3
  float radius;

  CUSTOM_Material *mat_ptr;
  // CUSTOM_Sphere() : type(CUSTOM_SPHERE){ }
} CUSTOM_Sphere;

typedef struct CUSTOM_HitRecord {
  float t;
  CUSTOM_Vector position; //vec3
  CUSTOM_Vector normal;   //vec3

  CUSTOM_Material *mat_ptr;
} CUSTOM_HitRecord;

typedef struct CUSTOM_Ray {
  CUSTOM_Vector origin;   //vec3
  CUSTOM_Vector direction;//vec3
} CUSTOM_Ray;

typedef struct CUSTOM_Camera {
  CUSTOM_Vector origin;             //vec3
  CUSTOM_Vector lower_left_corner;  //vec3
  CUSTOM_Vector horizontal;         //vec3
  CUSTOM_Vector vertical;           //vec3
  CUSTOM_Vector u, v, w;            //vec3
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

CUSTOM_Sphere *CUSTOM_SphereCreate(const float center[3], float r, CUSTOM_Material *material);

CUSTOM_Lambertian *CUSTOM_LambertianCreate(const float albedo[3]);
CUSTOM_Metal *CUSTOM_MetalCreate(const float albedo[3], const float f);
CUSTOM_Dielectric *CUSTOM_DielectricCreate(const float ri);

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

inline float randomFloat() {  return rand() / (RAND_MAX + 1.0);}

