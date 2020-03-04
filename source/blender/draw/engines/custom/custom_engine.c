#include "DRW_render.h"

#include "BLI_rand.h"

#include "BKE_object.h"
#include "BKE_global.h" /* for G.debug_value */

#include "DEG_depsgraph_query.h"

#include "DNA_world_types.h"

#include "custom_private.h"
#include "custom_engine.h" /* own include */

#define CUSTOM_ENGINE "MAOCHINN_CUSTOM"

static void custom_engine_init(void *vedata)
{
  CUSTOM_Data *data = vedata;
  CUSTOM_StorageList *stl = data->stl;
  const DRWContextState *draw_ctx = DRW_context_state_get();
  const View3D *v3d = draw_ctx->v3d;

  if (!stl->pd)
    stl->pd = MEM_callocN(sizeof(*stl->pd), __func__);

  CUSTOM_PrivateData *pd = stl->pd;

  pd->hide_overlays = (v3d->flag2 & V3D_HIDE_OVERLAYS) != 0;
  if (!stl->pd->hide_overlays)
    stl->pd->overlay = v3d->overlay;
  else
    memset(&stl->pd->overlay, 0, sizeof(stl->pd->overlay));

  if (v3d->shading.type == OB_WIRE)
    stl->pd->overlay.flag |= V3D_OVERLAY_WIREFRAMES;
}

static void custom_engine_free(void)
{
  CUSTOM_shader_free();
}

static void custom_cache_init(void *vedata)
{
  CUSTOM_Data *data = vedata;
  CUSTOM_PassList *psl = data->psl;
  CUSTOM_PrivateData *pd = data->stl->pd;
  const DRWContextState *draw_ctx = DRW_context_state_get();

  View3DShading *shading = &draw_ctx->v3d->shading;

  pd->wire_step_param = pd->overlay.wireframe_threshold - 254.0f / 255.0f;

  bool is_wire_shmode = (shading->type == OB_WIRE);
  bool is_material_shmode = (shading->type > OB_SOLID);
  bool is_object_color = is_wire_shmode && (shading->wire_color_type == V3D_SHADING_OBJECT_COLOR);
  bool is_random_color = is_wire_shmode && (shading->wire_color_type == V3D_SHADING_RANDOM_COLOR);

  GPUShader *shader = CUSTOM_shader();

  DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL |
                   DRW_STATE_STENCIL_EQUAL | DRW_STATE_FIRST_VERTEX_CONVENTION;
  // DRW_STATE_STENCIL_EQUAL 會執行stencil test

  DRW_PASS_CREATE(psl->custom_default_ps, state);
  DRWPass *pass = psl->custom_default_ps;

  DRWShadingGroup *grp = pd->shgrp = DRW_shgroup_create(shader, pass);
  DRW_shgroup_uniform_block_persistent(grp, "globalsBlock", G_draw.block_ubo);
  DRW_shgroup_uniform_float_copy(grp, "wireStepParam", pd->wire_step_param);
  DRW_shgroup_uniform_bool_copy(grp, "useColoring", true);
  DRW_shgroup_uniform_bool_copy(grp, "isTransform", (G.moving & G_TRANSFORM_OBJ) != 0);
  DRW_shgroup_uniform_bool_copy(grp, "isObjectColor", false);
  DRW_shgroup_uniform_bool_copy(grp, "isRandomColor", true);
  DRW_shgroup_stencil_mask(grp, 0xFF);  // stencil buffer 可寫
}

static void custom_cache_populate(void *vedata, Object *ob)
{
  CUSTOM_Data *data = vedata;
  CUSTOM_PrivateData *pd = data->stl->pd;
  struct GPUBatch *geom = NULL;

  if (ob->type == OB_MESH)
    geom = DRW_cache_mesh_face_wireframe_get(ob);
  else if (ob->type == OB_CURVE)
    geom = DRW_cache_curve_edge_wire_get(ob);

  if (geom) {
    DRW_shgroup_call(pd->shgrp, geom, ob);
  }
}

static void custom_cache_finish(void *vedata)
{

}

static void custom_draw_scene(void *vedata)
{
  CUSTOM_Data* data = vedata;
  CUSTOM_PassList *psl = data->psl;
  CUSTOM_PrivateData *pd = data->stl->pd;

  DRW_draw_pass(psl->custom_default_ps);
}

static const DrawEngineDataSize custom_data_size = DRW_VIEWPORT_DATA_SIZE(CUSTOM_Data);

DrawEngineType draw_engine_custom_type = {
    NULL,
    NULL,
    N_("Custom"),
    &custom_data_size,
    &custom_engine_init,
    &custom_engine_free,
    &custom_cache_init,
    &custom_cache_populate,
    &custom_cache_finish,
    NULL,
    &custom_draw_scene,
    NULL,
    NULL,
    NULL,
};

RenderEngineType DRW_engine_custom_type = {
    NULL,
    NULL,
    CUSTOM_ENGINE,
    N_("Custom"),
    RE_INTERNAL | RE_USE_PREVIEW,
    NULL,
    &DRW_render_to_image,
    NULL,
    NULL,
    NULL,
    NULL,
    &custom_render_update_passes,
    &draw_engine_custom_type,
    {NULL, NULL, NULL},
};

#undef CUSTOM_ENGINE
