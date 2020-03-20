#include "DRW_render.h"

#include "custom_private.h"

CUSTOM_ViewLayerData* CUSTOM_view_layer_data_ensure(void)
{
	CUSTOM_ViewLayerData** vldata = (CUSTOM_ViewLayerData**)DRW_view_layer_engine_data_ensure(
		&draw_engine_custom_type, &CUSTOM_view_layer_data_free);

	if (*vldata == NULL) {
    *vldata = (CUSTOM_ViewLayerData *)MEM_callocN(sizeof(**vldata), "CUSTOM_ViewLayerData");
	}

	return *vldata;
}

void CUSTOM_view_layer_data_free(void* storage)
{
	CUSTOM_ViewLayerData* vldata = (CUSTOM_ViewLayerData*)storage;

	DRW_UBO_FREE_SAFE(vldata->custom_ubo);
}
