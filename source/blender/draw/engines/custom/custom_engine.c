#include "DRW_render.h"

#include "BLI_rand.h"

#include "BKE_object.h"
#include "BKE_global.h" /* for G.debug_value */
#include "BKE_camera.h"
#include "BKE_mesh.h"
#include "BKE_node.h"
#include "BKE_image.h"

#include "DEG_depsgraph_query.h"
#include "IMB_imbuf_types.h"

#include "DNA_world_types.h"

#include "custom_private.h"
#include "custom_engine.h" /* own include */

///
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_camera_types.h"
#include "DNA_node_types.h"

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

  if (!stl->pd) {
    stl->pd = MEM_callocN(sizeof(*stl->pd), __func__);
  }


  if (vldata->custom_ubo == NULL) {
    vldata->custom_ubo = DRW_uniformbuffer_create(sizeof(vldata->custom_data),
                                                  &vldata->custom_data);
  }
}

static void custom_engine_free(void)
{
  // CUSTOM_shader_free();
}

static void custom_cache_init(void *vedata)
{
  // CUSTOM_Data *data = vedata;
  // CUSTOM_PassList *psl = data->psl;
  // CUSTOM_PrivateData *pd = data->stl->pd;
  // CUSTOM_ViewLayerData *vldata = CUSTOM_view_layer_data_ensure();
  // const DRWContextState *draw_ctx = DRW_context_state_get();

  // GPUShader *shader = CUSTOM_shader();

  // DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL |
  //                 DRW_STATE_STENCIL_EQUAL | DRW_STATE_FIRST_VERTEX_CONVENTION;
  // DRW_STATE_STENCIL_EQUAL 會執行stencil test

  // DRW_PASS_CREATE(psl->custom_ps, state);
  // DRWPass *pass = psl->custom_ps;

  // DRWShadingGroup *grp = pd->shgrp = DRW_shgroup_create(shader, pass);

  // DRW_shgroup_uniform_block(grp, "custom_block", vldata->custom_ubo);
  // DRW_shgroup_uniform_block_persistent(grp, "globalsBlock", G_draw.block_ubo);
  // DRW_shgroup_uniform_float_copy(grp, "wireStepParam", 1.0f);
  // DRW_shgroup_uniform_bool_copy(grp, "useColoring", true);
  // DRW_shgroup_uniform_bool_copy(grp, "isTransform", (G.moving & G_TRANSFORM_OBJ) != 0);
  // DRW_shgroup_uniform_bool_copy(grp, "isObjectColor", false);
  // DRW_shgroup_uniform_bool_copy(grp, "isRandomColor", true);
  // DRW_shgroup_stencil_mask(grp, 0xFF);  // stencil buffer 可寫
}

static void custom_cache_populate(void *vedata, Object *ob)
{
  CUSTOM_Data *data = vedata;
  CUSTOM_PrivateData *pd = data->stl->pd;
  const DRWContextState *draw_ctx = DRW_context_state_get();
  const Scene *scene = draw_ctx->scene;
  float unit_scale = scene->unit.scale_length;

  float loc[3] = {ob->loc[0] * unit_scale, ob->loc[1] * unit_scale, ob->loc[2] * unit_scale};

  if (ob->type == OB_MESH) {
    const Mesh *mesh = ob->data;

    const CUSTOM_Material *material = CUSTOM_LambertianCreate(
        (float[3]){0.5f, 0.5f, 0.5f}, CUSTOM_ConstantTextureCreate((float[3]){0.5f, 0.5f, 0.5f}));

    if (mesh->mat[0] && mesh->mat[0]->use_nodes) {
      const ListBase nodes = mesh->mat[0]->nodetree->nodes;
      for (const bNode *node = nodes.first; node; node = node->next) {
        if (node->type == SH_NODE_TEX_IMAGE) {
          NodeTexImage *tex = node->storage;
          Image *ima = (Image *)node->id;
          ImBuf *ibuf;
          void *lock;
          int i, size;

          ibuf = BKE_image_acquire_ibuf(ima, NULL, &lock);

          if (ibuf) {
            size = ibuf->x * ibuf->y * ibuf->channels;
            float *pixels = (float *)MEM_callocN(sizeof(float) * size, "float array");

            if (ibuf->rect_float) {
              memcpy(pixels, ibuf->rect_float, sizeof(float) * size);
            }
            else {
              for (i = 0; i < size; i++) {
                pixels[i] = ((unsigned char *)ibuf->rect)[i] * (1.0f / 255.0f);
              }
            }
            material = CUSTOM_LambertianCreate(
                (float[3]){0.5f, 0.5f, 0.5f}, CUSTOM_ImageTextureCreate(pixels, ibuf->x, ibuf->y));
          }

          BKE_image_release_ibuf(ima, ibuf, lock);
        }
        else if (node->type == SH_NODE_BSDF_DIFFUSE) {
          float color[3] = {1.0f, 0.0f, 0.0f};
          for (const bNodeSocket *iosock = node->inputs.first; iosock; iosock = iosock->next) {
            if (strcmp("Color", iosock->name) == 0 && iosock->type == SOCK_RGBA)
              copy_v3_v3(color, ((bNodeSocketValueRGBA *)iosock->default_value)->value);
          }
          material = CUSTOM_LambertianCreate(
              (float[3]){0.5f, 0.5f, 0.5f},
              CUSTOM_ConstantTextureCreate((float[3]){color[0], color[1], color[2]}));
        }
        else if (node->type == SH_NODE_BSDF_GLOSSY) {
          float color[3] = {1.0f, 0.0f, 0.0f};
          float roughness = 0.0f;
          for (const bNodeSocket *iosock = node->inputs.first; iosock; iosock = iosock->next) {
            if (strcmp("Color", iosock->name) == 0 && iosock->type == SOCK_RGBA)
              copy_v3_v3(color, ((bNodeSocketValueRGBA *)iosock->default_value)->value);
            else if (strcmp("Roughness", iosock->name) == 0 && iosock->type == SOCK_FLOAT)
              roughness = ((bNodeSocketValueFloat *)iosock->default_value)->value;
          }
          material = CUSTOM_MetalCreate((float[3]){0.5f, 0.5f, 0.5f}, roughness);
        }
      }
    }
    if (strncmp("OBSphere", ob->id.name, 8) == 0) {
      BLI_addtail(
          &pd->bvh_nodes,
          CUSTOM_SphereCreate((float[3]){loc[0], loc[1], loc[2]}, 0.5f * ob->scale[0], material));
    }
    else if (strncmp("OBCube", ob->id.name, 6) == 0) {
      float min[3], max[3];
      madd_v3_v3v3v3(min, loc, (float[3]){-0.5f, -0.5f, -0.5f}, ob->scale);
      madd_v3_v3v3v3(max, loc, (float[3]){0.5f, 0.5f, 0.5f}, ob->scale);
      BLI_addtail(&pd->bvh_nodes, CUSTOM_BoxCreate(min, max, material));
    }
    else if (strncmp("OBPlane", ob->id.name, 7) == 0) {
      float dimension[3];
      BKE_object_dimensions_get(ob, dimension);
      const MPoly *mpoly = mesh->mpoly;
      float normal[3];

      BKE_mesh_calc_poly_normal(mpoly, mesh->mloop + mpoly->loopstart, mesh->mvert, normal);
      if (dimension[0] <= 0.001f)
        BLI_addtail(&pd->world,
                    CUSTOM_RectYZCreate(loc[1] - 0.5f * ob->scale[1],
                                        loc[1] + 0.5f * ob->scale[1],
                                        loc[2] - 0.5f * ob->scale[2],
                                        loc[2] + 0.5f * ob->scale[2],
                                        loc[0],
                                        material,
                                        normal[0] < 0.0f ? true : false));
      else if (dimension[1] <= 0.001f)
        BLI_addtail(&pd->world,
                    CUSTOM_RectXZCreate(loc[0] - 0.5f * ob->scale[0],
                                        loc[0] + 0.5f * ob->scale[0],
                                        loc[2] - 0.5f * ob->scale[2],
                                        loc[2] + 0.5f * ob->scale[2],
                                        loc[1],
                                        material,
                                        normal[1] < 0.0f ? true : false));
      else if (dimension[2] <= 0.001f)
        BLI_addtail(&pd->world,
                    CUSTOM_RectXYCreate(loc[0] - 0.5f * ob->scale[0],
                                        loc[0] + 0.5f * ob->scale[0],
                                        loc[1] - 0.5f * ob->scale[1],
                                        loc[1] + 0.5f * ob->scale[1],
                                        loc[2],
                                        material,
                                        normal[2] < 0.0f ? true : false));
    }
  }
  else if (ob->type == OB_LAMP) {
    const Light *light = ob->data;
    if (light->type == LA_AREA) {
      const CUSTOM_Material *material = CUSTOM_DiffuseLightCreate(
          CUSTOM_ConstantTextureCreate((float[3]){4.0f, 4.0f, 4.0f}));
      BLI_addtail(&pd->world,
                  CUSTOM_RectXYCreate(loc[0] - 0.5f * ob->scale[0],
                                      loc[0] + 0.5f * ob->scale[0],
                                      loc[1] - 0.5f * ob->scale[1],
                                      loc[1] + 0.5f * ob->scale[1],
                                      loc[2],
                                      material,
                                      true));
    }
  }
  else if (ob->type == OB_CAMERA) {
    const Camera *cam = ob->data;
    const float *size = DRW_viewport_size_get();
    // CameraParams params;
    float sensor = BKE_camera_sensor_size(cam->sensor_fit, cam->sensor_x, cam->sensor_y);

    float look_at[3];
    sub_v3_v3v3(look_at, loc, ob->obmat[2]);

    const DRWContextState *draw_ctx = DRW_context_state_get();
    const Scene *scene = draw_ctx->scene;
    int resolution_x = scene->r.xsch;
    int resolution_y = scene->r.ysch;

    float aspect = (float)resolution_x / (float)resolution_y;
    // float aspect = cam->sensor_x / cam->sensor_y;
    // float aspect = size[0] / size[1];

    if (cam->dof.flag) {
      // use depth of field
      float focus_distance;
      if (cam->dof.focus_object)
        focus_distance = len_v3v3(loc, cam->dof.focus_object->loc);
      else
        focus_distance = cam->dof.focus_distance;

      pd->camera = CUSTOM_CameraCreate(loc,
                                       look_at,
                                       (float[3]){0.0f, 0.0f, 1.0f},
                                       focallength_to_fov(cam->lens, sensor),
                                       aspect,
                                       cam->dof.aperture_fstop,
                                       focus_distance);
    }
    else {
      pd->camera = CUSTOM_CameraCreate(loc,
                                       look_at,
                                       (float[3]){0.0f, 0.0f, 1.0f},
                                       focallength_to_fov(cam->lens, sensor),
                                       aspect,
                                       0.0f,
                                       1.0f);
    }
  }
}

static void custom_cache_finish(void *vedata)
{
}

static void custom_draw_scene(void *vedata)
{
  CUSTOM_Data *data = vedata;
  CUSTOM_TextureList *txl = data->txl;
  // CUSTOM_PassList *psl = data->psl;
  CUSTOM_PrivateData *pd = data->stl->pd;

  CUSTOM_ViewLayerData *vldata = CUSTOM_view_layer_data_ensure();

  const DRWContextState *draw_ctx = DRW_context_state_get();
  const Scene *scene = draw_ctx->scene;

  const float *size = DRW_viewport_size_get();
  const float precentage = (float)scene->r.size / 100.0f;
  CUSTOM_Vector *frame = NULL;
  int wdt = (int)size[0] * precentage;
  int hgt = (int)size[1] * precentage;
  CUSTOM_ImageCreate(&frame, wdt, hgt);

  if (!BLI_listbase_is_empty(&pd->bvh_nodes)) {
    int count = BLI_listbase_count(&pd->bvh_nodes);
    struct CUSTOM_Hittable **nodes = MEM_calloc_arrayN(count, sizeof(CUSTOM_Hittable*), "CUSTOM_Hittable");
    int i = 0;
    for (CUSTOM_Hittable *hittable = (CUSTOM_Hittable *)pd->bvh_nodes.first; hittable;
         hittable = hittable->next) {
      nodes[i++] = hittable;
    }
    BLI_addtail(&pd->world, CUSTOM_BvhNodeCreate(nodes, 0, count));
  }

  if (!BLI_listbase_is_empty(&pd->world)) {
    int nx = (float)wdt / (float)hgt < pd->camera->aspect ? wdt :
                                                            (int)((float)hgt * pd->camera->aspect);
    int ny = (int)((float)nx / pd->camera->aspect);
    // int ns = 10;

    int ns = scene->custom.viewport_samples;

    int offset_x = (wdt - nx) / 2;
    int offset_y = (hgt - ny) / 2;

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
  BLI_freelistN(&pd->bvh_nodes);
  pd->world.first = pd->world.last = NULL;
  pd->bvh_nodes.first = pd->bvh_nodes.last = NULL;
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
