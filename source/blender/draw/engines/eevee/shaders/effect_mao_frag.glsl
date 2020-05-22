out vec4 FragColor;

#ifdef DEBUG_MAO
uniform sampler2D occlusionBuffer;

void main()
{
  vec2 texel_size = 1.0 / vec2(textureSize(occlusionBuffer, 0)).xy;
  vec2 uvs = saturate(gl_FragCoord.xy * texel_size);

  float occlusion = textureLod(occlusionBuffer, uvs, 0.0).r;
  FragColor = vec4(vec3(occlusion), 1.0);
}

#else
layout(location = 0) out float occlusion_factor;

uniform sampler2D normalBuffer;
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
  // occlusion = 1.0 - (occlusion / sampleCount);

  // FragColor = vec4(occlusion, occlusion, occlusion, 1.0);
  occlusion_factor = (1.0 - occlusion);
}
#endif
