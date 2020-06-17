#pragma once

#ifndef __CUSTOM_PRIVATE_H__
#define __CUSTOM_PRIVATE_H__

#include "custom_ray_tracer.h"

extern struct DrawEngineType draw_engine_custom_type;

typedef struct CUSTOM_UniformBuffer {
  float matrix[4][4];   /* mat4 */
  float vector[4];    /* vec4 */
} CUSTOM_UniformBuffer;

typedef struct CUSTOM_ViewLayerData {
  struct CUSTOM_UniformBuffer custom_data;
  struct GPUUniformBuffer *custom_ubo;
} CUSTOM_ViewLayerData;

typedef struct CUSTOM_FramebufferList {
	struct GPUFrameBuffer* custom_fb;
} CUSTOM_FramebufferList;

typedef struct CUSTOM_TextureList {
	struct GPUTexture* custom_tx;
} CUSTOM_TextureList;

typedef struct CUSTOM_PassList {
	DRWPass* custom_ps;
} CUSTOM_PassList;

typedef struct CUSTOM_PrivateData {
	DRWShadingGroup* shgrp;

  
  ListBase world;
  ListBase bvh_nodes;
  CUSTOM_Camera *camera;
} CUSTOM_PrivateData;

typedef struct CUSTOM_StorageList {
	struct CUSTOM_PrivateData* pd;
} CUSTOM_StorageList;

typedef struct CUSTOM_Data {
	void* engine_type;
  CUSTOM_FramebufferList *fbl;
	CUSTOM_TextureList* txl;
	CUSTOM_PassList* psl;
	CUSTOM_StorageList* stl;
} CUSTOM_Data;


void custom_render_update_passes(struct RenderEngine *engine,
                                 struct Scene *scene,
                                 struct ViewLayer *view_layer);

GPUShader *CUSTOM_shader(void);
void CUSTOM_shader_free(void);

void CUSTOM_view_layer_data_free(void *sldata);
CUSTOM_ViewLayerData *CUSTOM_view_layer_data_ensure(void);


#endif /* __CUSTOM_PRIVATE_H__ */
