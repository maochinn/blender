/*
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Copyright 2016, Blender Foundation.
 */

/** \file
 * \ingroup draw_engine
 *
 * Implementation of the screen space Ground Truth Ambient Occlusion.
 */

#include "DRW_render.h"

#include "BLI_string_utils.h"
#include "BLI_rand.h"

#include "DEG_depsgraph_query.h"

#include "BKE_global.h" /* for G.debug_value */

#include "eevee_private.h"

#include "GPU_extensions.h"
#include "GPU_platform.h"
#include "GPU_state.h"

static struct {
  /* Screen Space Ambient Occlusion */
  struct GPUShader *ssao_sh;
  /* Screen Space Curvature*/
  struct GPUShader *ssc_sh;

  struct GPUShader *mao_debug_sh;

  struct GPUTexture *hammersley;
  struct GPUTexture *colormap;
  struct GPUTexture *dummy;  // must declare a texture, avoid texture == NULL
} e_data = {NULL};           /* Engine data */

extern char datatoc_ambient_occlusion_lib_glsl[];
extern char datatoc_common_view_lib_glsl[];
extern char datatoc_common_uniforms_lib_glsl[];
extern char datatoc_bsdf_common_lib_glsl[];
extern char datatoc_bsdf_sampling_lib_glsl[];
extern char datatoc_effect_mao_frag_glsl[];

static void eevee_mao_create_shader(void)
{
  char *mao_frag_str = BLI_string_joinN(datatoc_common_view_lib_glsl,
                                        datatoc_common_uniforms_lib_glsl,
                                        datatoc_bsdf_common_lib_glsl,
                                        datatoc_ambient_occlusion_lib_glsl,
                                        datatoc_bsdf_sampling_lib_glsl,
                                        datatoc_effect_mao_frag_glsl);

  e_data.ssao_sh = DRW_shader_create_fullscreen(
      mao_frag_str,
      "#define SSAO\n"
      "#define HAMMERSLEY_SIZE " STRINGIFY(HAMMERSLEY_SIZE) "\n");
  e_data.ssc_sh = DRW_shader_create_fullscreen(mao_frag_str, "#define SSC\n");
  e_data.mao_debug_sh = DRW_shader_create_fullscreen(mao_frag_str, "#define DEBUG_MAO\n");

  MEM_freeN(mao_frag_str);
}

static struct GPUTexture *create_hammersley_sample_texture(int samples)
{
  struct GPUTexture *tex;
  float(*texels)[2] = MEM_mallocN(sizeof(float[2]) * samples, "hammersley_tex");
  int i;

  for (i = 0; i < samples; i++) {
    double dphi;
    BLI_hammersley_1d(i, &dphi);
    float phi = (float)dphi * 2.0f * M_PI;
    texels[i][0] = cosf(phi);
    texels[i][1] = sinf(phi);
  }

  tex = DRW_texture_create_1d(samples, GPU_RG16F, DRW_TEX_WRAP, (float *)texels);
  MEM_freeN(texels);
  return tex;
}

static struct GPUTexture *create_colormap_texture()
{
  float color[] = {
    0.000000, 0.000000, 0.517825,
    0.000000, 0.000000, 0.553476,
    0.000000, 0.000000, 0.589127,
    0.000000, 0.000000, 0.624777,
    0.000000, 0.000000, 0.660428,
    0.000000, 0.000000, 0.696078,
    0.000000, 0.000000, 0.731729,
    0.000000, 0.000000, 0.767380,
    0.000000, 0.000000, 0.803030,
    0.000000, 0.000000, 0.820856,
    0.000000, 0.000000, 0.874332,
    0.000000, 0.000000, 0.892157,
    0.000000, 0.000000, 0.945633,
    0.000000, 0.000000, 0.963458,
    0.000000, 0.000000, 1.000000,
    0.000000, 0.000000, 1.000000,
    0.000000, 0.017647, 1.000000,
    0.000000, 0.049020, 1.000000,
    0.000000, 0.064706, 1.000000,
    0.000000, 0.111765, 1.000000,
    0.000000, 0.143137, 1.000000,
    0.000000, 0.174510, 1.000000,
    0.000000, 0.190196, 1.000000,
    0.000000, 0.237255, 1.000000,
    0.000000, 0.268627, 1.000000,
    0.000000, 0.300000, 1.000000,
    0.000000, 0.315686, 1.000000,
    0.000000, 0.362745, 1.000000,
    0.000000, 0.394118, 1.000000,
    0.000000, 0.425490, 1.000000,
    0.000000, 0.441176, 1.000000,
    0.000000, 0.488235, 1.000000,
    0.000000, 0.519608, 1.000000,
    0.000000, 0.550980, 1.000000,
    0.000000, 0.582353, 1.000000,
    0.000000, 0.613725, 1.000000,
    0.000000, 0.645098, 1.000000,
    0.000000, 0.660784, 1.000000,
    0.000000, 0.707843, 1.000000,
    0.000000, 0.739216, 1.000000,
    0.000000, 0.770588, 1.000000,
    0.000000, 0.801961, 1.000000,
    0.000000, 0.833333, 1.000000,
    0.000000, 0.864706, 0.996205,
    0.000000, 0.896078, 0.970904,
    0.009488, 0.911765, 0.958254,
    0.047438, 0.958824, 0.920304,
    0.072739, 0.990196, 0.895003,
    0.098039, 1.000000, 0.869703,
    0.123340, 1.000000, 0.844402,
    0.148640, 1.000000, 0.819102,
    0.173941, 1.000000, 0.793801,
    0.199241, 1.000000, 0.768501,
    0.211891, 1.000000, 0.755851,
    0.249842, 1.000000, 0.717900,
    0.275142, 1.000000, 0.692600,
    0.300443, 1.000000, 0.667299,
    0.325743, 1.000000, 0.641999,
    0.351044, 1.000000, 0.616698,
    0.376344, 1.000000, 0.591398,
    0.401645, 1.000000, 0.566097,
    0.414295, 1.000000, 0.553447,
    0.452245, 1.000000, 0.515497,
    0.477546, 1.000000, 0.490196,
    0.502846, 1.000000, 0.464896,
    0.528147, 1.000000, 0.439595,
    0.553447, 1.000000, 0.414295,
    0.578748, 1.000000, 0.388994,
    0.604048, 1.000000, 0.363694,
    0.629349, 1.000000, 0.338393,
    0.654649, 1.000000, 0.313093,
    0.679949, 1.000000, 0.287793,
    0.705250, 1.000000, 0.262492,
    0.730550, 1.000000, 0.237192,
    0.755851, 1.000000, 0.211891,
    0.768501, 1.000000, 0.199241,
    0.806452, 1.000000, 0.161290,
    0.831752, 1.000000, 0.135990,
    0.857052, 1.000000, 0.110689,
    0.882353, 1.000000, 0.085389,
    0.907653, 1.000000, 0.060089,
    0.932954, 1.000000, 0.034788,
    0.958254, 0.973856, 0.009488,
    0.983555, 0.944808, 0.000000,
    1.000000, 0.915759, 0.000000,
    1.000000, 0.886710, 0.000000,
    1.000000, 0.857662, 0.000000,
    1.000000, 0.828613, 0.000000,
    1.000000, 0.799564, 0.000000,
    1.000000, 0.770516, 0.000000,
    1.000000, 0.741467, 0.000000,
    1.000000, 0.726943, 0.000000,
    1.000000, 0.683370, 0.000000,
    1.000000, 0.654321, 0.000000,
    1.000000, 0.625272, 0.000000,
    1.000000, 0.596224, 0.000000,
    1.000000, 0.567175, 0.000000,
    1.000000, 0.538126, 0.000000,
    1.000000, 0.509078, 0.000000,
    1.000000, 0.480029, 0.000000,
    1.000000, 0.450980, 0.000000,
    1.000000, 0.421932, 0.000000,
    1.000000, 0.392883, 0.000000,
    1.000000, 0.363834, 0.000000,
    1.000000, 0.334786, 0.000000,
    1.000000, 0.305737, 0.000000,
    1.000000, 0.276688, 0.000000,
    1.000000, 0.262164, 0.000000,
    1.000000, 0.218591, 0.000000,
    1.000000, 0.189542, 0.000000,
    1.000000, 0.160494, 0.000000,
    1.000000, 0.131445, 0.000000,
    1.000000, 0.102397, 0.000000,
    0.999109, 0.073348, 0.000000,
    0.963458, 0.044299, 0.000000,
    0.927807, 0.015251, 0.000000,
    0.892157, 0.000000, 0.000000,
    0.856506, 0.000000, 0.000000,
    0.820856, 0.000000, 0.000000,
    0.785205, 0.000000, 0.000000,
    0.749554, 0.000000, 0.000000,
    0.713904, 0.000000, 0.000000,
    0.678253, 0.000000, 0.000000,
    0.660428, 0.000000, 0.000000,
    0.606952, 0.000000, 0.000000,
    0.571301, 0.000000, 0.000000,
    0.535651, 0.000000, 0.000000,
    0.500000, 0.000000, 0.000000,
  };
  struct GPUTexture *tex;
  tex = DRW_texture_create_1d(128, GPU_RGB16F, DRW_TEX_FILTER, (float *)color);
  return tex;
}

int EEVEE_mao_init(EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata)
{
  EEVEE_CommonUniformBuffer *common_data = &sldata->common_data;
  EEVEE_FramebufferList *fbl = vedata->fbl;
  EEVEE_StorageList *stl = vedata->stl;
  EEVEE_EffectsInfo *effects = stl->effects;

  const DRWContextState *draw_ctx = DRW_context_state_get();
  const Scene *scene_eval = DEG_get_evaluated_scene(draw_ctx->depsgraph);

  if (!e_data.dummy) {
    float pixel[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    e_data.dummy = DRW_texture_create_2d(1, 1, GPU_RGBA8, DRW_TEX_WRAP, pixel);
  }

  if (!e_data.hammersley) {
    e_data.hammersley = create_hammersley_sample_texture(HAMMERSLEY_SIZE);
    e_data.colormap = create_colormap_texture();
    eevee_mao_create_shader();
  }

  if (scene_eval->eevee.flag & SCE_EEVEE_GTAO_ENABLED) {
    const float *viewport_size = DRW_viewport_size_get();
    const int fs_size[2] = {(int)viewport_size[0], (int)viewport_size[1]};

    common_data->ao_dist = scene_eval->eevee.gtao_distance;

    effects->ssao = DRW_texture_pool_query_2d(
        fs_size[0], fs_size[1], GPU_R8, &draw_engine_eevee_type);
    GPU_framebuffer_ensure_config(&fbl->ssao_fb,
                                  {GPU_ATTACHMENT_NONE, GPU_ATTACHMENT_TEXTURE(effects->ssao)});
    effects->ssc = DRW_texture_pool_query_2d(
        fs_size[0], fs_size[1], GPU_R8, &draw_engine_eevee_type);
    GPU_framebuffer_ensure_config(&fbl->ssc_fb,
                                  {GPU_ATTACHMENT_NONE, GPU_ATTACHMENT_TEXTURE(effects->ssc)});

    if (G.debug_value == 33) {
      effects->mao_debug = DRW_texture_pool_query_2d(
          fs_size[0], fs_size[1], GPU_RGBA8, &draw_engine_eevee_type);
      GPU_framebuffer_ensure_config(
          &fbl->mao_debug_fb, {GPU_ATTACHMENT_NONE, GPU_ATTACHMENT_TEXTURE(effects->mao_debug)});
    }
    else {
      effects->mao_debug = NULL;
    }
    return EFFECT_GTAO | EFFECT_NORMAL_BUFFER;
  }

  /* Cleanup */
  effects->ssao = e_data.dummy;
  effects->ssc = e_data.dummy;
  GPU_FRAMEBUFFER_FREE_SAFE(fbl->ssao_fb);
  GPU_FRAMEBUFFER_FREE_SAFE(fbl->ssc_fb);
  common_data->ao_dist = 0.0f;

  return 0;
}

void EEVEE_mao_cache_init(EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata)
{
  EEVEE_PassList *psl = vedata->psl;
  EEVEE_StorageList *stl = vedata->stl;
  // EEVEE_TextureList *txl = vedata->txl;
  EEVEE_EffectsInfo *effects = stl->effects;
  // DefaultTextureList *dtxl = DRW_viewport_texture_list_get();

  struct GPUBatch *quad = DRW_cache_fullscreen_quad_get();

  if ((effects->enabled_effects & EFFECT_GTAO) != 0) {
    DRW_PASS_CREATE(psl->ssao, DRW_STATE_WRITE_COLOR);
    DRWShadingGroup *grp = DRW_shgroup_create(e_data.ssao_sh, psl->ssao);
    DRW_shgroup_uniform_texture_ref(grp, "depthBuffer", &effects->ao_src_depth);
    DRW_shgroup_uniform_texture_ref(grp, "normalBuffer", &effects->ssr_normal_input);
    DRW_shgroup_uniform_block(grp, "common_block", sldata->common_ubo);
    DRW_shgroup_uniform_float_copy(grp, "sampleCount", 64.0f);
    DRW_shgroup_uniform_float_copy(grp, "invSampleCount", 1 / 64.0f);
    DRW_shgroup_uniform_texture(grp, "texHammersley", e_data.hammersley);
    DRW_shgroup_call(grp, quad, NULL);

    DRW_PASS_CREATE(psl->ssc, DRW_STATE_WRITE_COLOR);
    grp = DRW_shgroup_create(e_data.ssc_sh, psl->ssc);
    DRW_shgroup_uniform_texture_ref(grp, "depthBuffer", &effects->ao_src_depth);
    DRW_shgroup_uniform_texture_ref(grp, "normalBuffer", &effects->ssr_normal_input);
    DRW_shgroup_uniform_vec2_copy(grp, "viewportSize", DRW_viewport_size_get(), 2);
    DRW_shgroup_uniform_block(grp, "common_block", sldata->common_ubo);
    DRW_shgroup_call(grp, quad, NULL);

    if (G.debug_value == 33) {
      DRW_PASS_CREATE(psl->mao_debug, DRW_STATE_WRITE_COLOR);
      grp = DRW_shgroup_create(e_data.mao_debug_sh, psl->mao_debug);
      DRW_shgroup_uniform_texture_ref(grp, "depthBuffer", &effects->ao_src_depth);
      DRW_shgroup_uniform_texture_ref(grp, "normalBuffer", &effects->ssr_normal_input);
      DRW_shgroup_uniform_block(grp, "common_block", sldata->common_ubo);
      DRW_shgroup_uniform_texture_ref(grp, "occlusionBuffer", &effects->ssao);
      DRW_shgroup_uniform_texture_ref(grp, "curvatureBuffer", &effects->ssc);
      DRW_shgroup_uniform_texture(grp, "colormapBuffer", e_data.colormap);
      DRW_shgroup_call(grp, quad, NULL);
    }
  }
}

void EEVEE_mao_compute(EEVEE_ViewLayerData *UNUSED(sldata),
                       EEVEE_Data *vedata,
                       struct GPUTexture *depth_src)
{
  EEVEE_PassList *psl = vedata->psl;
  EEVEE_FramebufferList *fbl = vedata->fbl;
  EEVEE_StorageList *stl = vedata->stl;
  EEVEE_EffectsInfo *effects = stl->effects;

  if ((effects->enabled_effects & EFFECT_GTAO) != 0) {
    DRW_stats_group_start("SSAO");
    effects->ao_src_depth = depth_src;

    GPU_framebuffer_bind(fbl->ssao_fb);

    DRW_draw_pass(psl->ssao);

    if (GPU_mip_render_workaround() ||
        GPU_type_matches(GPU_DEVICE_INTEL_UHD, GPU_OS_WIN, GPU_DRIVER_ANY)) {
      /* Fix dot corruption on intel HD5XX/HD6XX series. */
      GPU_flush();
    }

    /* Restore */
    GPU_framebuffer_bind(fbl->main_fb);

    DRW_stats_group_end();

    //
    DRW_stats_group_start("SSC");
    effects->ao_src_depth = depth_src;

    GPU_framebuffer_bind(fbl->ssc_fb);

    DRW_draw_pass(psl->ssc);

    if (GPU_mip_render_workaround() ||
        GPU_type_matches(GPU_DEVICE_INTEL_UHD, GPU_OS_WIN, GPU_DRIVER_ANY)) {
      /* Fix dot corruption on intel HD5XX/HD6XX series. */
      GPU_flush();
    }

    /* Restore */
    GPU_framebuffer_bind(fbl->main_fb);

    DRW_stats_group_end();
  }
}

void EEVEE_mao_draw_debug(EEVEE_ViewLayerData *UNUSED(sldata), EEVEE_Data *vedata)
{
  EEVEE_PassList *psl = vedata->psl;
  EEVEE_FramebufferList *fbl = vedata->fbl;
  EEVEE_StorageList *stl = vedata->stl;
  EEVEE_EffectsInfo *effects = stl->effects;

  if (((effects->enabled_effects & EFFECT_GTAO) != 0) && (G.debug_value == 33)) {
    DRW_stats_group_start("MAO Debug");

    GPU_framebuffer_bind(fbl->mao_debug_fb);
    DRW_draw_pass(psl->mao_debug);

    /* Restore */
    GPU_framebuffer_bind(fbl->main_fb);

    DRW_stats_group_end();
  }
}

void EEVEE_mao_free(void)
{
  DRW_SHADER_FREE_SAFE(e_data.ssao_sh);
  DRW_SHADER_FREE_SAFE(e_data.ssc_sh);
  DRW_SHADER_FREE_SAFE(e_data.mao_debug_sh);
  DRW_TEXTURE_FREE_SAFE(e_data.hammersley);
}
