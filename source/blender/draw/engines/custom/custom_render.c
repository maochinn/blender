
#include "DRW_engine.h"
#include "DRW_render.h"

#include "DNA_node_types.h"
#include "DNA_object_types.h"

#include "BKE_global.h"
#include "BKE_object.h"

#include "BLI_rand.h"
#include "BLI_rect.h"

#include "DEG_depsgraph_query.h"

#include "GPU_framebuffer.h"
#include "GPU_extensions.h"
#include "GPU_state.h"
#include "GPU_shader.h"

#include "RE_pipeline.h"

#include "custom_private.h"

/*shader*/
// lib
extern char datatoc_common_view_lib_glsl[];
extern char datatoc_common_globals_lib_glsl[];
extern char datatoc_gpu_shader_common_obinfos_lib_glsl[];

extern char datatoc_custom_uniform_lib_glsl[];
    // vertex shader
extern char datatoc_custom_vert_glsl[];
// fragment shadet
extern char datatoc_custom_frag_glsl[];

/* *********** STATIC *********** */

// engine data
static struct {
  struct GPUShader *sh_data;
} e_data = {NULL};

/* *********** FUNCTION *********** */

void custom_render_update_passes(RenderEngine *engine, Scene *scene, ViewLayer *view_layer)
{
  RE_engine_register_pass(engine, scene, view_layer, RE_PASSNAME_COMBINED, 4, "RGBA", SOCK_RGBA);
}

// shader

GPUShader *CUSTOM_shader(void)
{
  const DRWContextState *draw_ctx = DRW_context_state_get();
  const GPUShaderConfigData *sh_cfg = &GPU_shader_cfg_data[draw_ctx->sh_cfg];
  struct GPUShader *sh_data = e_data.sh_data;
  if (!sh_data) {
    sh_data = GPU_shader_create_from_arrays({
        .vert = (const char *[]){sh_cfg->lib,
                                 datatoc_common_view_lib_glsl,
                                 datatoc_common_globals_lib_glsl,
                                 datatoc_gpu_shader_common_obinfos_lib_glsl,
                                 datatoc_custom_vert_glsl,
                                 NULL},
        .frag = (const char *[]){datatoc_common_view_lib_glsl,
                                 datatoc_custom_uniform_lib_glsl,
                                 datatoc_custom_frag_glsl,
                                 NULL},
        .defs = (const char *[]){sh_cfg->def, NULL},
    });
  }
  return sh_data;
}

void CUSTOM_shader_free(void)
{
  DRW_SHADER_FREE_SAFE(e_data.sh_data);
}
