#pragma once

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

typedef struct CUSTOM_Material {
  int temp;
} CUSTOM_Material;

/*Hittable type*/
enum { CUSTOM_SPHERE = 0, CUSTOM_TYPE_MAX };

typedef struct CUSTOM_Hittable {
  struct CUSTOM_Hittable *next, *prev;
  const short type;
} CUSTOM_Hittable;
typedef struct CUSTOM_Sphere {
  struct CUSTOM_Hittable *next, *prev;
  const short type;
  CUSTOM_Vector center;
  float radius;

  CUSTOM_Material *mat_ptr;
  // CUSTOM_Sphere() : type(CUSTOM_SPHERE){ }
} CUSTOM_Sphere;
typedef struct CUSTOM_HitRecord {
  float t;
  CUSTOM_Vector position;
  CUSTOM_Vector normal;

  CUSTOM_Material *mat_ptr;
} CUSTOM_HitRecord;

typedef struct CUSTOM_Ray {
  CUSTOM_Vector origin;
  CUSTOM_Vector direction;
} CUSTOM_Ray;

typedef struct CUSTOM_Camera {
  CUSTOM_Vector origin;
  CUSTOM_Vector lower_left_corner;
  CUSTOM_Vector horizontal;
  CUSTOM_Vector vertical;
  CUSTOM_Vector u, v, w;
  float aspect;
  float lens_radius;
} CUSTOM_Camera;

CUSTOM_Vector CUSTOM_sampleColor(const CUSTOM_Ray *ray, const ListBase *world);
CUSTOM_Camera *CUSTOM_CameraCreate(float lookfrom[3],
                                   float lookat[3],
                                   float vup[3],
                                   float hfov,
                                   float aspect,
                                   float aperture,
                                   float focus_dist);
CUSTOM_Ray CUSTOM_CameraGetRay(const CUSTOM_Camera *cam, float u, float v);

CUSTOM_Sphere *CUSTOM_SphereCreate(float center[3], float r);
void CUSTOM_ImageCreate(CUSTOM_Vector **image, const int wdt, const int hgt);

inline float random_float() {  return rand() / (RAND_MAX + 1.0);}
