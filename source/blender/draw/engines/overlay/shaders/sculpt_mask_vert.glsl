
uniform float maskOpacity;

in vec3 pos;
in float msk;

out vec4 finalColor;

void main()
{
  vec3 world_pos = point_object_to_world(pos);
  gl_Position = point_world_to_ndc(world_pos);

  finalColor = vec4(0.0, 0.0, 0.0, msk * maskOpacity);

#ifdef USE_WORLD_CLIP_PLANES
  world_clip_planes_calc_clip_distance(world_pos);
#endif
}
