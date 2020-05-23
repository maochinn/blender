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
} e_data = {NULL}; /* Engine data */

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
      mao_frag_str, "#define SSAO\n" "#define HAMMERSLEY_SIZE " STRINGIFY(HAMMERSLEY_SIZE) "\n");
  e_data.ssc_sh = DRW_shader_create_fullscreen(
      mao_frag_str, "#define SSC\n");
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

int EEVEE_mao_init(EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata)
{
  EEVEE_CommonUniformBuffer *common_data = &sldata->common_data;
  EEVEE_FramebufferList *fbl = vedata->fbl;
  EEVEE_StorageList *stl = vedata->stl;
  EEVEE_EffectsInfo *effects = stl->effects;

  const DRWContextState *draw_ctx = DRW_context_state_get();
  const Scene *scene_eval = DEG_get_evaluated_scene(draw_ctx->depsgraph);

  if (!e_data.hammersley) {
    e_data.hammersley = create_hammersley_sample_texture(HAMMERSLEY_SIZE);
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
          &fbl->mao_debug_fb,
          {GPU_ATTACHMENT_NONE, GPU_ATTACHMENT_TEXTURE(effects->mao_debug)});
    }
    else {
      effects->mao_debug = NULL;
    }
    return EFFECT_GTAO | EFFECT_NORMAL_BUFFER;
  }

  /* Cleanup */
  GPU_FRAMEBUFFER_FREE_SAFE(fbl->ssao_fb);
  //common_data->ao_settings = 0.0f;

  return 0;
}

void EEVEE_mao_cache_init(EEVEE_ViewLayerData *sldata, EEVEE_Data *vedata)
{
  EEVEE_PassList *psl = vedata->psl;
  EEVEE_StorageList *stl = vedata->stl;
  //EEVEE_TextureList *txl = vedata->txl;
  EEVEE_EffectsInfo *effects = stl->effects;
  //DefaultTextureList *dtxl = DRW_viewport_texture_list_get();

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
