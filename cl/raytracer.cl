typedef float4 Sphere;
typedef struct Ray { float3 origin; float3 direction; } Ray;
typedef struct Material { int type; float4 attributes; } Material;
typedef struct Hitrec { float t; float3 p; float3 normal; Material mat;} Hitrec;
typedef struct Camera {
  float3 origin;
  float3 low_left;
  float3 horizontal;
  float3 vertical;
  float3 u, v, w;
  float lens_radius;
} Camera;

float3
Unit(float3 v)
{
  return v / fast_length(v);
}

float3
RandomInUnitDisk(constant float3 *uvs)
{
  return (2.0f * (float3)(uvs[get_global_id(2)].x, uvs[get_global_id(2)].z, 0.0)) - (float3)(1.0, 1.0, 0.0);
}

Camera
SetupCamera(const float3 lookfrom, const float3 lookat, const float3 vup, const float vfov,
            const float aspect, const float aperture, const float focus_dist)
{
  Camera camera;
  camera.lens_radius = aperture * 0.5;
  camera.origin = lookfrom;
  float theta = vfov * M_PI_F / 180.0;
  float half_height = tan(theta/2.0f);
  float half_width = half_height * aspect;
  camera.w = Unit(lookfrom - lookat);
  camera.u = Unit(cross(vup, camera.w));
  camera.v = cross(camera.w, camera.u);
  camera.low_left = camera.origin - half_width * focus_dist * camera.u - 
                    half_height * focus_dist * camera.v - focus_dist * camera.w;
  camera.horizontal = 2 * half_width * focus_dist * camera.u;
  camera.vertical = 2 * half_height * focus_dist * camera.v;

  return camera;
}

Ray
GetRay(const Camera camera, const float s, const float t, constant float3 *uvs)
{
  float3 rd = camera.lens_radius * RandomInUnitDisk(uvs);
  float3 offset = camera.u * rd.x + camera.v * rd.y;
  Ray ray;
  ray.origin = camera.origin + offset;
  ray.direction = camera.low_left + s * camera.horizontal + t * camera.vertical - camera.origin - offset;

  return ray;
}

float3
Scale(float3 v, float s)
{
  float3 result = {v.x * s, v.y * s, v.z * s};
  return result;
}

float
Rand(constant float2 *uvs)
{
  // return (abs(rand()) % 445) / 444.0;
  // return float(abs(rand()) % 100001) / 100000.0;
  switch (rand() % 2) {
    case 0:
      return uvs[get_global_id(2)].x;
    case 1:
      return uvs[get_global_size(2) - get_global_id(2)].y;
  }
}

float3
Sky(const Ray ray)
{
  float3 uv = ray.direction / fast_length(ray.direction);
  float t = 0.5 * (uv.y + 1.0);
  return Scale((float3)(1.0, 1.0, 1.0), (1.0 - t)) + Scale((float3)(0.5, 0.7, 1.0), t);
}

float3
PointAt(Ray ray, float t)
{
  return ray.origin + Scale(ray.direction, t);
}

float
Schlick(float cosine, float ref_idx)
{
  float r0 = pow(float(1.0 - ref_idx) / float(1.0 + ref_idx), 2);
  return r0 + (1 - r0) * pow((1 - cosine), 5);
}

float3
Reflect(const float3 v, const float3 n)
{
  return v - 2.0f * dot(v, n) * n;
}

bool
Refract(const float3 v, const float3 n, float ni_over_nt, float3 *refracted)
{
  float3 uv = Unit(v);
  float dt = dot(uv, n);
  float disc = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
  if (disc > 0.0)
  {
    float3 p = ni_over_nt * (uv - n * dt) - n * sqrt(disc);
    refracted->x = p.x;
    refracted->y = p.y;
    refracted->z = p.z;
    return true;
  }
  return false;
}

bool
InSphere(Sphere s, Material mat, Ray ray, float t_min, float t_max, Hitrec *rec)
{
  float3 oc = ray.origin - s.xyz;
  float a = dot(ray.direction, ray.direction);
  float b = dot(oc, ray.direction);
  float c = dot(oc, oc) - s.w * s.w;
  float disc = b * b - a * c;
  if (disc > 0.0)
  {
    float tmp = (-b - sqrt(disc)) / a;
    if (t_min < tmp && tmp < t_max)
    {
      rec->t = tmp;
      rec->p = PointAt(ray, tmp);
      rec->normal = (rec->p - s.xyz) / s.w;
      rec->mat = mat;
      return true;
    }
    tmp = (-b + sqrt(disc)) / a;
    if (t_min < tmp && tmp < t_max)
    {
      rec->t = tmp;
      rec->p = PointAt(ray, tmp);
      rec->normal = (rec->p - s.xyz) / s.w;
      rec->mat = mat;
      return true;
    }
  }
  return false;
}

float
SphereDistance(Sphere s, Ray ray)
{
  float3 oc = ray.origin - s.xyz;
  float a = dot(ray.direction, ray.direction);
  float b = 2.0 * dot(oc, ray.direction);
  float c = dot(oc, oc) - s.w * s.w;
  float disc = b * b - 4 * a * c;
  if (disc < 0.0)
  {
    return -1.0;
  }
  return (-b - sqrt(disc)) / (2.0 * a);
}

float3
RandomInUnitSphere(constant float3 *uvs)
{
  // float3 p = (float3)(1.0, 1.0, 1.0);
  // while (p.x * p.x + p.y * p.y + p.z * p.z >= 1.0) {
  //   p = (2.0f * (float3)(Rand(uvs), Rand(uvs), Rand(uvs))) - (float3)(1.0, 1.0, 1.0);
  // }
  // return p;
  int id0 = get_global_id(0)+ get_global_id(1) + get_global_id(2);
  int id2 = get_global_id(2);
  int id = id0 % get_global_size(2);
  return (float)(uvs[id].x, uvs[id].y, uvs[id].z);
}

bool
Scatter(Ray ray, Hitrec rec, float3 *attenuation, Ray *scattered, constant float3 *uvs)
{
  switch (rec.mat.type) {
    case 0: {
      float3 target = rec.p + rec.normal + RandomInUnitSphere(uvs);
      scattered->origin = rec.p;
      scattered->direction = target - rec.p;
      attenuation->x = rec.mat.attributes.x;
      attenuation->y = rec.mat.attributes.y;
      attenuation->z = rec.mat.attributes.z;
      return true;
      }
    case 1:{
      float3 reflected = Reflect(Unit(ray.direction), rec.normal);
      scattered->origin = rec.p;
      scattered->direction = reflected + rec.mat.attributes.w * RandomInUnitSphere(uvs);
      attenuation->x = rec.mat.attributes.x;
      attenuation->y = rec.mat.attributes.y;
      attenuation->z = rec.mat.attributes.z;

      return dot(scattered->direction, rec.normal) > 0.0;
      }
    case 2: {
      float3 outward_normal;
      float3 reflected = Reflect(ray.direction, rec.normal);
      float ni_over_nt = 1.0;
      attenuation->x = 1.0; attenuation->y = 1.0; attenuation->z = 1.0;
      float3 refracted;
      float reflect_prob = 1.0;
      float cosine = 1.0;
      if (dot(ray.direction, rec.normal) > 0.0)
      {
        outward_normal = -rec.normal;
        ni_over_nt = rec.mat.attributes.x;
        cosine = rec.mat.attributes.x * dot(ray.direction, rec.normal) / fast_length(ray.direction);
      }
      else
      {
        outward_normal = rec.normal;
        ni_over_nt = 1.0 / rec.mat.attributes.x;
        cosine = -dot(ray.direction, rec.normal) / fast_length(ray.direction);
      }
      if (Refract(ray.direction, outward_normal, ni_over_nt, &refracted))
      {
        reflect_prob = Schlick(cosine, rec.mat.attributes.x);
      }
      if (uvs[get_global_id(2)].x < reflect_prob)
      {
        scattered->origin = rec.p;
        scattered->direction = reflected;
      }
      else
      {
        scattered->origin = rec.p;
        scattered->direction = refracted;
      }
      return true;
    }
  }
}


float3
Color(Ray ray, constant Sphere *s, int n_sphere, constant int *types, constant float4 *mats, constant float3 *uvs, int depth)
{
  bool is_hit = true;
  float3 color = (float3)(1.0, 1.0, 1.0);
  while (is_hit)
  {
    Hitrec rec;
    rec.t = 999999.9;
    is_hit = false;
    for (int i = 0; i < n_sphere; ++i)
    {
      Material mat;
      mat.type = types[i];
      mat.attributes = mats[i];
      is_hit = InSphere(s[i], mat, ray, 0.001, rec.t, &rec) || is_hit;
    }
    if (is_hit)
    {
      Ray scattered;
      float3 attenuation;
      if (depth < 50 && Scatter(ray, rec, &attenuation, &scattered, uvs))
      {
        color *= attenuation;
        ray.origin = scattered.origin;
        ray.direction = scattered.direction;
        ++depth;
      }
      else
      {
        color *= attenuation;
        return color;
      }
    }
  }
  return color * Sky(ray);
}

kernel
void draw_sky(global const float *d_pixels, global float *d_results, int2 size, float3 low, float3 hor, float3 ver, float3 ori)
{
  const int x = get_global_id(0);
  const int y = get_global_id(1);
  const float u = float(x) / size.x;
  const float v = float(y) / size.y;
  
  const float3 dir = low + u * hor + v * ver;
  struct Ray ray;
  ray.origin = ori;
  ray.direction = dir;
  float3 color = Sky(ray);

  const int idx = (y * size.x + x) * 3;
  d_results[idx] = color.x;
  d_results[idx + 1] = color.y;
  d_results[idx + 2] = color.z;
}

// kernel
// void draw_sphere(global const float *d_pixels, global float *d_results, int2 size, float3 low, float3 hor, float3 ver, float3 ori, constant Sphere *s)
// {
//   const int x = get_global_id(0);
//   const int y = get_global_id(1);
//   const float u = float(x) / size.x;
//   const float v = float(y) / size.y;
  
//   const float3 dir = low + u * hor + v * ver;
//   struct Ray ray;
//   ray.origin = ori;
//   ray.direction = dir;
//   float3 color = Color(ray, s, 1);

//   const int idx = (y * size.x + x) * 3;
//   d_results[idx] = color.x;
//   d_results[idx + 1] = color.y;
//   d_results[idx + 2] = color.z;
// }

// kernel
// void draw_spheres(global const float *d_pixels, global float *d_results, int2 size, float3 low, float3 hor, float3 ver, float3 ori, constant Sphere *s, int n_sphere)
// {
//   const int x = get_global_id(0);
//   const int y = get_global_id(1);
//   const float u = float(x) / size.x;
//   const float v = float(y) / size.y;
  
//   const float3 dir = low + u * hor + v * ver;
//   struct Ray ray;
//   ray.origin = ori;
//   ray.direction = dir;
//   float3 color = Color(ray, s, n_sphere);

//   const int idx = (y * size.x + x) * 3;
//   d_results[idx] = color.x;
//   d_results[idx + 1] = color.y;
//   d_results[idx + 2] = color.z;
// }

// kernel
// void diffuse(global float *d_results, int2 size, float3 low, float3 hor, float3 ver, float3 ori, constant Sphere *s, int n_sphere)
// {
//   const int x = get_global_id(0);
//   const int y = get_global_id(1);
//   const int sn = get_global_id(2);

//   const float u = float(x + Rand()) / size.x;
//   const float v = float(y + Rand()) / size.y;

//   const float3 dir = low + u * hor + v * ver - ori;
//   Ray ray;
//   ray.origin = ori;
//   ray.direction = dir;
//   float3 color = Color(ray, s, n_sphere);
  
//   const int idx = ((y * size.x + x) * get_global_size(2) + sn) * 3;
//     // const int idx = (y * size.x + x) * 3;
//   d_results[idx] = color.x;
//   d_results[idx + 1] = color.y;
//   d_results[idx + 2] = color.z;

//   // colors[t_num] = Color(ray, s, n_sphere, t_num);
//   // // reduction
//   // for (int offset = 1; offset < num_items; offset *= 2)
//   // {
//   //  int mask = 2 * offset - 1;
//   //  barrier(CLK_LOCAL_MEM_FENCE);
//   //  if ((t_num & mask) == 0)
//   //  {
//   //    colors[t_num] += colors[t_num + offset];
//   //  }
//   // }

//   // barrier(CLK_LOCAL_MEM_FENCE);
//   // if (t_num == 0)
//   // {
//   //   const int idx = (y * size.x + x) * 3;
//   //   d_results[idx] = colors[0].x / num_items;
//   //   d_results[idx + 1] = colors[0].y / num_items;
//   //   d_results[idx + 2] = colors[0].z / num_items;
//   // }
// }

// kernel
// void Materials(global float *d_results, float3 low, float3 hor, float3 ver, float3 ori,
//    constant Sphere *s, int n_sphere, constant int *types, constant float4 *materials, constant float3 *uvs)//,
//  //  local Material *mats)
// {
//   const int x = get_global_id(0);
//   const int y = get_global_id(1);
//   const int sn = get_global_id(2);

//   const int idx = ((y * get_global_size(0) + x) * get_global_size(2) + sn) * 3;

//   const float u = float(x + uvs[sn].x) / get_global_size(0);
//   const float v = float(y + uvs[sn].y) / get_global_size(1);

//   const float3 dir = low + u * hor + v * ver - ori;
//   Ray ray;
//   ray.origin = ori;
//   ray.direction = dir;

//   // for (int i = 0; i < n_sphere; ++i)
//   // {
//   //   Material mat;
//   //   mat.type = types[i];
//   //   mat.attributes = (float4)(materials[i].x, materials[i].y, materials[i].z, materials[i].w);
//   //   mats[i] = mat;
//   // }

//   float3 color = Color(ray, s, n_sphere, types, materials, uvs, 0);
  
//     // const int idx = (y * size.x + x) * 3;
//   d_results[idx] = color.x;
//   d_results[idx + 1] = color.y;
//   d_results[idx + 2] = color.z;
// }

kernel
void Draw(global float *d_results, constant Sphere *s, int n_sphere, constant int *types,
                constant float4 *materials, constant float3 *uvs)
{
  const int x = get_global_id(0);
  const int y = get_global_id(1);
  const int sn = get_global_id(2);

  const int idx = ((y * get_global_size(0) + x) * get_global_size(2) + sn) * 3;

  float3 lookfrom = (float3)(11.0, 2.0, 3.0);
  float3 lookat = (float3)(0.0, 0.0, 0.0);
  float dist_to_focus = fast_length(lookfrom - lookat);
  float aperture = 0.1;

  const Camera camera = SetupCamera(lookfrom, lookat, (float3)(0.0, 1.0, 0.0), 20.0,
                                    float(get_global_size(0))/ float(get_global_size(1)), aperture, dist_to_focus);

  const float u = float(x + uvs[sn].x) / get_global_size(0);
  const float v = float(y + uvs[sn].y) / get_global_size(1);

  Ray ray = GetRay(camera, u, v, uvs);

  float3 color = Color(ray, s, n_sphere, types, materials, uvs, 0);
  
    // const int idx = (y * size.x + x) * 3;
  d_results[idx] = color.x;
  d_results[idx + 1] = color.y;
  d_results[idx + 2] = color.z;
}