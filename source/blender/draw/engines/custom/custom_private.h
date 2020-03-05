#pragma once

#ifndef __CUSTOM_PRIVATE_H__
#define __CUSTOM_PRIVATE_H__

typedef struct CUSTOM_FramebufferList {
	struct GPUFrameBuffer* custom_default_fb;
} CUSTOM_FramebufferList;

typedef struct CUSTOM_TextureList {
	struct GPUTexture* custom_default_tx;
} CUSTOM_TextureList;

typedef struct CUSTOM_PassList {
	DRWPass* custom_default_ps;
} CUSTOM_PassList;

typedef struct CUSTOM_PrivateData {
	DRWShadingGroup* shgrp;
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


void custom_render_update_passes(struct RenderEngine* engine,
	struct Scene* scene,
	struct ViewLayer* view_layer);

GPUShader *CUSTOM_shader(void);
void CUSTOM_shader_free(void);

#endif /* __CUSTOM_PRIVATE_H__ */
