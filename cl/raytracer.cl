typedef float4 Sphere;
typedef struct Hitrec { float t; float3 p; float3 normal; } Hitrec;
typedef struct Ray { float3 origin; float3 direction; } Ray;

kernel
void gradation(global const float *d_pixels, global float *d_results, int2 size)
{
  const int x = get_global_id(0);
  const int y = get_global_id(1);
  const int idx = (y * size.x + x) * 3;
  d_results[idx] = float(x) / 200.0;
  d_results[idx + 1] = float(y) / 100.0;
  d_results[idx + 2] = 0.5;
}

float3
Scale(float3 v, float s)
{
  float3 result = {v.x * s, v.y * s, v.z * s};
  return result;
}

float3
Unit(float3 v)
{
  return v / fast_length(v);
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

bool
InSphere(Sphere s, Ray ray, float t_min, float t_max, Hitrec *rec)
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
      return true;
    }
    tmp = (-b + sqrt(disc)) / a;
        if (t_min < tmp && tmp < t_max)
    {
      rec->t = tmp;
      rec->p = PointAt(ray, tmp);
      rec->normal = (rec->p - s.xyz) / s.w;
      return true;
    }
    return false;
  }
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
Color(Ray ray, constant Sphere *s, int n_sphere)
{
  struct Hitrec rec;
  bool is_hit = false;
  rec.t = 999999.9;
  for (int i = 0; i < n_sphere; ++i)
  {
    is_hit = InSphere(s[i], ray, 0.0, rec.t, &rec) || is_hit;
  }
  if (is_hit)
  {
    return Scale((float3)(rec.normal.x + 1.0, rec.normal.y + 1.0, rec.normal.z + 1.0), 0.5);
  }
  return Sky(ray);
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

kernel
void draw_sphere(global const float *d_pixels, global float *d_results, int2 size, float3 low, float3 hor, float3 ver, float3 ori, constant Sphere *s)
{
  const int x = get_global_id(0);
  const int y = get_global_id(1);
  const float u = float(x) / size.x;
  const float v = float(y) / size.y;
  
  const float3 dir = low + u * hor + v * ver;
  struct Ray ray;
  ray.origin = ori;
  ray.direction = dir;
  float3 color = Color(ray, s, 1);

  const int idx = (y * size.x + x) * 3;
  d_results[idx] = color.x;
  d_results[idx + 1] = color.y;
  d_results[idx + 2] = color.z;
}

kernel
void draw_spheres(global const float *d_pixels, global float *d_results, int2 size, float3 low, float3 hor, float3 ver, float3 ori, constant Sphere *s, int n_sphere)
{
  const int x = get_global_id(0);
  const int y = get_global_id(1);
  const float u = float(x) / size.x;
  const float v = float(y) / size.y;
  
  const float3 dir = low + u * hor + v * ver;
  struct Ray ray;
  ray.origin = ori;
  ray.direction = dir;
  float3 color = Color(ray, s, n_sphere);

  const int idx = (y * size.x + x) * 3;
  d_results[idx] = color.x;
  d_results[idx + 1] = color.y;
  d_results[idx + 2] = color.z;
}

kernel
void antialias(global const float *d_pixels, global float *d_results, int2 size, float3 low, float3 hor, float3 ver, float3 ori, constant Sphere *s, int n_sphere)
{
  const int x = get_global_id(0);
  const int y = get_global_id(1);
  const int z = get_global_id(2);
  const float u = float(x) / size.x;
  const float v = float(y) / size.y;

  const float u = float(x + drand48()) / size.x;
  const float v = float(y + drand48()) / size.y;
  
  const float3 dir = low + u * hor + v * ver - ori;
  struct Ray ray;
  ray.origin = ori;
  ray.direction = dir;
  colors[t_num] = Color(ray, s, n_sphere);
  // reduction
  //for (int offset = 1; offset < num_items; offset *= 2)
  //{
  //  int mask = 2 * offset - 1;
  //  barrier(CLK_LOCAL_MEM_FENCE);
  //  if ((t_num & mask) == 0)
  //  {
  //    colors[t_num] += colors[t_num + offset];
  //  }
  //}

  //barrier(CLK_LOCAL_MEM_FENCE);
  //if (t_num == 0)
  //{
    const int idx = (y * size.x + x) * 3;
    d_results[idx] = colors[0].x / num_items;
    d_results[idx + 1] = colors[0].y / num_items;
    d_results[idx + 2] = colors[0].z / num_items;
  //}
}