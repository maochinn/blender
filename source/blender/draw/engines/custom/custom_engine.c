#include "DRW_render.h"

#include "BLI_rand.h"

#include "BKE_object.h"
#include "BKE_global.h" /* for G.debug_value */
#include "BKE_camera.h"

#include "DEG_depsgraph_query.h"

#include "DNA_world_types.h"

#include "custom_private.h"
#include "custom_engine.h" /* own include */

///
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_camera_types.h"

#include "custom_ray_tracer.h"
///

#define CUSTOM_ENGINE "PBRT_CUSTOM"

static void custom_engine_init(void *vedata)
{
  CUSTOM_Data *data = vedata;
  CUSTOM_TextureList *txl = data->txl;
  CUSTOM_StorageList *stl = data->stl;
  CUSTOM_ViewLayerData *vldata = CUSTOM_view_layer_data_ensure();
  const DRWContextState *draw_ctx = DRW_context_state_get();
  const View3D *v3d = draw_ctx->v3d;

  if (!stl->pd)
    stl->pd = MEM_callocN(sizeof(*stl->pd), __func__);

  if (vldata->custom_ubo == NULL) {
    vldata->custom_ubo = DRW_uniformbuffer_create(sizeof(vldata->custom_data),
                                                  &vldata->custom_data);
  }
}

static void custom_engine_free(void)
{
  //CUSTOM_shader_free();
}

static void custom_cache_init(void *vedata)
{
  //CUSTOM_Data *data = vedata;
  //CUSTOM_PassList *psl = data->psl;
  //CUSTOM_PrivateData *pd = data->stl->pd;
  //CUSTOM_ViewLayerData *vldata = CUSTOM_view_layer_data_ensure();
  //const DRWContextState *draw_ctx = DRW_context_state_get();

  //GPUShader *shader = CUSTOM_shader();

  //DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL |
  //                 DRW_STATE_STENCIL_EQUAL | DRW_STATE_FIRST_VERTEX_CONVENTION;
  // DRW_STATE_STENCIL_EQUAL 會執行stencil test

  //DRW_PASS_CREATE(psl->custom_ps, state);
  //DRWPass *pass = psl->custom_ps;

  //DRWShadingGroup *grp = pd->shgrp = DRW_shgroup_create(shader, pass);

  //DRW_shgroup_uniform_block(grp, "custom_block", vldata->custom_ubo);
  //DRW_shgroup_uniform_block_persistent(grp, "globalsBlock", G_draw.block_ubo);
  //DRW_shgroup_uniform_float_copy(grp, "wireStepParam", 1.0f);
  //DRW_shgroup_uniform_bool_copy(grp, "useColoring", true);
  //DRW_shgroup_uniform_bool_copy(grp, "isTransform", (G.moving & G_TRANSFORM_OBJ) != 0);
  //DRW_shgroup_uniform_bool_copy(grp, "isObjectColor", false);
  //DRW_shgroup_uniform_bool_copy(grp, "isRandomColor", true);
  //DRW_shgroup_stencil_mask(grp, 0xFF);  // stencil buffer 可寫
}

static void custom_cache_populate(void *vedata, Object *ob)
{
  CUSTOM_Data *data = vedata;
  CUSTOM_PrivateData *pd = data->stl->pd;

  if (ob->type == OB_MESH) {
    BLI_addtail(&pd->world,
                CUSTOM_SphereCreate(
                    ob->loc,
                    0.5f,
                    (CUSTOM_Material *)CUSTOM_LambertianCreate((float[3]){0.8f, 0.3f, 0.3f})));
    BLI_addtail(&pd->world,
                CUSTOM_SphereCreate(
                    (float[3]){ob->loc[0], ob->loc[1], ob->loc[2] - 100.5f},
                    100.0f,
                    (CUSTOM_Material *)CUSTOM_LambertianCreate((float[3]){0.8f, 0.8f, 0.0f}))); 
    BLI_addtail(&pd->world,
                CUSTOM_SphereCreate(
                    (float[3]){ob->loc[0] + 1.0f, ob->loc[1], ob->loc[2]},
                    0.5f,
                    (CUSTOM_Material *)CUSTOM_MetalCreate((float[3]){0.8f, 0.6f, 0.2f}, 0.3f)));
    BLI_addtail(&pd->world,
                CUSTOM_SphereCreate((float[3]){ob->loc[0] - 1.0f, ob->loc[1], ob->loc[2]},
                                    0.5f,
                                    (CUSTOM_Material *)CUSTOM_DielectricCreate(1.5f)));
    BLI_addtail(&pd->world,
                CUSTOM_SphereCreate((float[3]){ob->loc[0] - 1.0f, ob->loc[1], ob->loc[2]},
                                    -0.45f,
                                    (CUSTOM_Material *)CUSTOM_DielectricCreate(1.5f)));
    BLI_addtail(&pd->world,
                CUSTOM_SphereCreate(
                    (float[3]){ob->loc[0], ob->loc[1] + 1.0f, ob->loc[2]},
                    0.5f,
                    (CUSTOM_Material *)CUSTOM_MetalCreate((float[3]){0.8f, 0.8f, 0.8f}, 1.0f)));
    BLI_addtail(&pd->world,
                CUSTOM_SphereCreate((float[3]){ob->loc[0], ob->loc[1]-1.0f, ob->loc[2]},
                                    0.5f,
                                    (CUSTOM_Material *)CUSTOM_DielectricCreate(1.5f)));
  }
  else if (ob->type == OB_CAMERA) {
    Camera *cam = ob->data;
    const float *size = DRW_viewport_size_get();
    CameraParams params;
    float sensor = BKE_camera_sensor_size(cam->sensor_fit, cam->sensor_x, cam->sensor_y);

    float look_at[3];
    sub_v3_v3v3(look_at, ob->loc, ob->obmat[2]);

    const DRWContextState *draw_ctx = DRW_context_state_get();
    const Scene *scene = draw_ctx->scene;
    int resolution_x = scene->r.xsch;
    int resolution_y = scene->r.ysch;

    float aspect = (float)resolution_x / (float)resolution_y;
    //float aspect = cam->sensor_x / cam->sensor_y;
    //float aspect = size[0] / size[1];

    pd->camera = CUSTOM_CameraCreate(ob->loc,
                                     look_at,
                                     (float[3]){0.0f, 0.0f, 1.0f},
                                     focallength_to_fov(cam->lens, sensor),
                                     aspect,
                                     cam->dof.aperture_fstop,
                                     cam->dof.focus_distance);
  }
}

static void custom_cache_finish(void *vedata)
{
}

static void custom_draw_scene(void *vedata)
{
  CUSTOM_Data *data = vedata;
  CUSTOM_TextureList *txl = data->txl;
  //CUSTOM_PassList *psl = data->psl;
  CUSTOM_PrivateData *pd = data->stl->pd;

  CUSTOM_ViewLayerData *vldata = CUSTOM_view_layer_data_ensure();



  const float *size = DRW_viewport_size_get();
  CUSTOM_Vector *frame = NULL;
  int wdt = (int)size[0]/4;
  int hgt = (int)size[1]/4;
  CUSTOM_ImageCreate(&frame, wdt, hgt);

  const DRWContextState *draw_ctx = DRW_context_state_get();
  const Scene *scene = draw_ctx->scene;

  if (pd->world.first != NULL) {
    int nx = (float)wdt / (float)hgt < pd->camera->aspect ? wdt : (int)((float)hgt * pd->camera->aspect);
    int ny = (int)((float)nx / pd->camera->aspect);
    //int ns = 10;

    int ns = scene->custom.viewport_samples;

    int offset_x = (wdt - nx)/2;
    int offset_y = (hgt - ny)/2;

    for (int j = ny - 1; j >= 0; j--)
      for (int i = 0; i < nx; i++) {
        float color[3];
        zero_v3(color);
        for (int s = 0; s < ns; s++) {
          float u = (float)(i + randomFloat()) / (float)nx;
          float v = (float)(j + randomFloat()) / (float)ny;
          CUSTOM_Ray ray = CUSTOM_CameraGetRay(pd->camera, u, v);
          add_v3_v3(color, CUSTOM_sampleColor(&ray, &pd->world, 0).vec3);
        }
        mul_v3_fl(color, 1.0f / (float)ns);
        copy_v3_fl3(color, sqrt(color[0]), sqrt(color[1]), sqrt(color[2]));

        int idx_x = i + offset_x;
        int idx_y = j + offset_y;
        frame[idx_x + idx_y * wdt] = (CUSTOM_Vector){color[0], color[1], color[2], 1.0f};
      }
  }

  txl->custom_tx = DRW_texture_create_2d(wdt, hgt, GPU_RGBA16F, 0, (float *)frame);

  DRW_transform_to_display(txl->custom_tx, false, false);

  // free
  BLI_freelistN(&pd->world);
  pd->world.first = NULL;
  pd->world.last = NULL;
}
static void custom_view_update(void *vedata)
{
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
    &custom_view_update,
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
