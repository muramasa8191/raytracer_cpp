typedef float4 sphere;

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
Sky(float3 origin, float3 dir)
{
  float3 uv = dir / fast_length(dir);
  float t = 0.5 * (uv.y + 1.0);
  return Scale((float3)(1.0, 1.0, 1.0), (1.0 - t)) + Scale((float3)(0.5, 0.7, 1.0), t);
}

bool
InSphere(sphere s, float3 ori, float3 dir)
{
  float3 oc = ori - s.xyz;
  float a = dot(dir, dir);
  float b = 2.0 * dot(oc, dir);
  float c = dot(oc, oc) - s.w * s.w;
  float disc = b * b - 4 * a * c;
  return disc > 0.0;
}

float3
Color(sphere s, float3 ori, float3 dir)
{
  if (InSphere(s, ori, dir))
  {
    return (float3)(1.0, 0.0, 0.0);
  }
  return Sky(ori, dir);
}

kernel
void draw_sky(global const float *d_pixels, global float *d_results, int2 size, float3 low, float3 hor, float3 ver, float3 ori)
{
  const int x = get_global_id(0);
  const int y = get_global_id(1);
  const float u = float(x) / size.x;
  const float v = float(y) / size.y;
  
  const float3 dir = low + u * hor + v * ver;
  float3 color = Sky(ori, dir);

  const int idx = (y * size.x + x) * 3;
  d_results[idx] = color.x;
  d_results[idx + 1] = color.y;
  d_results[idx + 2] = color.z;
}

kernel
void draw_sphere(global const float *d_pixels, global float *d_results, int2 size, float3 low, float3 hor, float3 ver, float3 ori, sphere s)
{
  const int x = get_global_id(0);
  const int y = get_global_id(1);
  const float u = float(x) / size.x;
  const float v = float(y) / size.y;
  
  const float3 dir = low + u * hor + v * ver;
  float3 color = Color(s, ori, dir);

  const int idx = (y * size.x + x) * 3;
  d_results[idx] = color.x;
  d_results[idx + 1] = color.y;
  d_results[idx + 2] = color.z;
}