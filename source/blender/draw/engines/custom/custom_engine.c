#include "DRW_render.h"

#include "BLI_rand.h"

#include "BKE_object.h"
#include "BKE_global.h" /* for G.debug_value */

#include "DEG_depsgraph_query.h"

#include "DNA_world_types.h"

#include "custom_private.h"
#include "custom_engine.h" /* own include */

#define CUSTOM_ENGINE "MAOCHINN_CUSTOM"

static void custom_engine_init(void* vedata)
{

}

static void custom_engine_free(void)
{

}

static void custom_cache_init(void* vedata)
{

}

static void custom_cache_populate(void* vedata, Object* ob)
{

}

static void custom_cache_finish(void* vedata)
{

}

static void custom_draw_scene(void* vedata)
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
