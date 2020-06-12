out vec4 FragColor;

#if defined(DEBUG_MAO)
//uniform sampler2D occlusionBuffer;  //define in ambient_occlusion.glsl
uniform sampler2D curvatureBuffer;
uniform sampler2D normalBuffer;

uniform sampler1D colormapBuffer;
void main()
{
  vec2 uvs = saturate(gl_FragCoord.xy / vec2(textureSize(depthBuffer, 0).xy));
  float depth = textureLod(depthBuffer, uvs, 0.0).r;

  if (depth == 1.0) {
    /* Do not trace for background */
    FragColor = vec4(0.0);
    return;
  }

  /* Avoid self shadowing. */
  depth = saturate(depth - 3e-6); /* Tweaked for 24bit depth buffer. */

  vec3 viewPosition = get_view_space_from_depth(uvs, depth);
  vec3 V = viewCameraVec;
  vec3 normal = normal_decode(texture(normalBuffer, uvs).rg, V);
  vec3 N = normalize(normal);

  float visiblity = textureLod(occlusionBuffer, uvs, 0.0).r;
  float curvature = textureLod(curvatureBuffer, uvs, 0.0).r;
  float D_curvature = abs(dFdx(curvature)) + abs(dFdy(curvature));

  float contour = dot(N, V) < 0.2 ? 1.0 - dot(N, V) : 0.0;
  float suggestive_contour = abs(curvature) <= 0.01 && D_curvature > 0.0 ? 1.0 : 0.0;

  vec3 color = vec3(visiblity);
  color = contour > 0 ? contour * vec3(1.0, 0.0, 0.0) : color;
  color = suggestive_contour > 0 ? vec3(0.0, 0.0, 1.0) : color;
  //color += contour * vec3(0.5, 0.0, 0.0);
  //color += suggestive_contour * vec3(0.0, 0.5, 0.0);

  FragColor = vec4(color, 1.0);
  FragColor = vec4(texture(colormapBuffer, clamp(curvature, 0.0, 1.0)).rgb, 1.0);
}

#elif defined(SSAO)
layout(location = 0) out float occlusion_factor;

uniform sampler2D normalBuffer;
void main()
{
  vec2 uvs = saturate(gl_FragCoord.xy / vec2(textureSize(depthBuffer, 0).xy));
  float depth = textureLod(depthBuffer, uvs, 0.0).r;

  if (depth == 1.0) {
    /* Do not trace for background */
    //FragColor = vec4(0.0);
    occlusion_factor = 0.0;
    return;
  }

  /* Avoid self shadowing. */
  depth = saturate(depth - 3e-6); /* Tweaked for 24bit depth buffer. */

  vec3 viewPosition = get_view_space_from_depth(uvs, depth);
  vec3 V = viewCameraVec;
  vec3 normal = normal_decode(texture(normalBuffer, uvs).rg, V);

  vec3 N, T, B;
  N = normalize(normal);
  make_orthonormal_basis(N, T, B); /* Generate tangent space */

  // Iterate over the sample kernel and calculate occlusion factor
  float occlusion = 0.0;
  for (float i = 0; i < sampleCount; i++) {
    // sample
    vec3 L = sample_hemisphere(i, N, T, B); /* Microfacet normal */
    vec3 viewSample = viewPosition + L * aoDistance;

    vec4 offset = vec4(viewSample, 1.0);
    offset = ProjectionMatrix * offset;   // from view to clip-space
    offset.xyz /= offset.w;               // perspective divide
    offset.xyz = offset.xyz * 0.5 + 0.5;  // transform to range 0.0 - 1.0

    // get sample depth
    float sampleDepth = textureLod(depthBuffer, offset.xy, 0.0).r;
    float sampleDepthZ = get_view_z_from_depth(sampleDepth);

    // range check & accumulate
    float rangeCheck = smoothstep(0.0, 1.0, aoDistance / abs(viewPosition.z - sampleDepthZ));
    occlusion += (sampleDepthZ >= viewSample.z ? 1.0 : 0.0) * rangeCheck;
  }
  occlusion = occlusion / sampleCount;

  occlusion_factor = (1.0 - occlusion);
}
#elif defined(SSC)
layout(location = 0) out float curvature_factor;

uniform sampler2D normalBuffer;

uniform vec2 viewportSize;
void main()
{
  vec2 uvs = saturate(gl_FragCoord.xy / vec2(textureSize(depthBuffer, 0).xy));
  float depth = textureLod(depthBuffer, uvs, 0.0).r;

  if (depth == 1.0) {
    /* Do not trace for background */
    FragColor = vec4(0.0);
    return;
  }

  /* Avoid self shadowing. */
  depth = saturate(depth - 3e-6); /* Tweaked for 24bit depth buffer. */

  /*code from http://madebyevan.com/shaders/curvature/ */

  vec3 viewPosition = get_view_space_from_depth(uvs, depth);
  vec3 V = viewCameraVec;
  vec3 normal = normal_decode(textureLod(normalBuffer, uvs, 0.0).rg, V);
  vec3 N = normalize(normal);

  vec3 dNdx = dFdx(N);
  vec3 dNdy = dFdy(N);
  vec3 xneg = N - dNdx;
  vec3 xpos = N + dNdx;
  vec3 yneg = N - dNdy;
  vec3 ypos = N + dNdy;
  float curvature = (cross(xneg, xpos).y - cross(yneg, ypos).x) * 0.5 / length(viewPosition);
  //float x = dot(dNdx, dNdx);
  //float y = dot(dNdy, dNdy);
  //float  curvature = max(x, y);
  curvature_factor = curvature + 0.5;
  
  //float Cx = 2.0 / (viewportSize.x * tan(0.7853981625));
  //float Cy = 2.0 / (viewportSize.y * tan(0.7853981625));

  //float depth_c = texelFetch(depthBuffer, ivec2(gl_FragCoord.xy), 0).r;
  //float depth_d = texelFetch(depthBuffer, ivec2(gl_FragCoord.xy) + ivec2(0, -1), 0).r;
  //float depth_l = texelFetch(depthBuffer, ivec2(gl_FragCoord.xy) + ivec2(-1, 0), 0).r;
  //float depth_r = texelFetch(depthBuffer, ivec2(gl_FragCoord.xy) + ivec2(1, 0), 0).r;
  //float depth_u = texelFetch(depthBuffer, ivec2(gl_FragCoord.xy) + ivec2(0, 1), 0).r;

  ////depth_c = get_view_z_from_depth(depth_c);
  ////depth_d = get_view_z_from_depth(depth_d);
  ////depth_l = get_view_z_from_depth(depth_l);
  ////depth_r = get_view_z_from_depth(depth_r);
  ////depth_u = get_view_z_from_depth(depth_u);

  //float dx = (0.5f * (depth_r - depth_l));
  //float dy = (0.5f * (depth_u - depth_d));

  //float dxx = (depth_l - 2.0f * depth + depth_r);
  //float dyy = (depth_d - 2.0f * depth + depth_u);

  //float dx2 = dx * dx;
  //float dy2 = dy * dy;

  //float Cx2 = Cx * Cx;
  //float Cy2 = Cy * Cy;

  //float D = Cy2 * dx2 + Cx2 * dy2 + Cx2 * Cy2 * depth_c * depth_c;
  //float H = Cy * dxx * D - Cy * dx * (Cy2 * dx * dxx + Cx2 * Cy2 * depth_c * dx) + Cx * dyy * D -
  //          Cx * dy * (Cx2 * dy * dyy + Cx2 * Cy2 * depth_c * dy);
  //H /= pow(D, 3.0f / 2.0f);

  //curvature_factor = H * 0.5 + 0.5;
  //curvature_factor = viewportSize.x;
}
#endif
