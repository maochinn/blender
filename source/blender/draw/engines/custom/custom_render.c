
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

#include "RE_pipeline.h"

#include "custom_private.h"

void custom_render_update_passes(RenderEngine* engine, Scene* scene, ViewLayer* view_layer)
{
	RE_engine_register_pass(engine, scene, view_layer, RE_PASSNAME_COMBINED, 4, "RGBA", SOCK_RGBA);
}