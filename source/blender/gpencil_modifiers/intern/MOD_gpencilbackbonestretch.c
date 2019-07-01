/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
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
 * along with this program; if not, write to the Free Software  Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2017, Blender Foundation
 * This is a new part of Blender
 *
 * Contributor(s): Antonio Vazquez, Joshua Leung, Yiming Wu
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

/** \file blender/gpencil_modifiers/intern/MOD_gpencilstrokes.c
 *  \ingroup modifiers
 */

#include <stdio.h>

#include "MEM_guardedalloc.h"

#include "DNA_scene_types.h"
#include "DNA_object_types.h"
#include "DNA_gpencil_types.h"
#include "DNA_gpencil_modifier_types.h"

#include "BLI_blenlib.h"
#include "BLI_rand.h"
#include "BLI_math.h"
#include "BLI_utildefines.h"
#include "BLI_linklist.h"
#include "BLI_alloca.h"

#include "BKE_gpencil.h"
#include "BKE_gpencil_modifier.h"
#include "BKE_modifier.h"
#include "BKE_context.h"
#include "BKE_global.h"
#include "BKE_object.h"
#include "BKE_main.h"
#include "BKE_scene.h"
#include "BKE_layer.h"
#include "BKE_library_query.h"
#include "BKE_collection.h"
#include "BKE_mesh.h"
#include "BKE_mesh_mapping.h"

#include "bmesh.h"
#include "bmesh_tools.h"

#include "DEG_depsgraph.h"
#include "DEG_depsgraph_build.h"
#include "DEG_depsgraph_query.h"

#include "MOD_gpencil_util.h"
#include "MOD_gpencil_modifiertypes.h"

static void initData(GpencilModifierData *md)
{
  BackboneGpencilModifierData *gpmd = (BackboneGpencilModifierData *)md;
}

static void copyData(const GpencilModifierData *md, GpencilModifierData *target)
{
  BKE_gpencil_modifier_copyData_generic(md, target);
}

static void bakeModifier(Main *UNUSED(bmain),
                         Depsgraph *depsgraph,
                         GpencilModifierData *md,
                         Object *ob)
{

  bGPdata *gpd = ob->data;

  for (bGPDlayer *gpl = gpd->layers.first; gpl; gpl = gpl->next) {
    for (bGPDframe *gpf = gpl->frames.first; gpf; gpf = gpf->next) {
      BackboneGpencilModifierData *lmd = (BackboneGpencilModifierData *)md;
      bGPDstroke *gps;
      for (gps = gpf->strokes.first; gps; gps = gps->next) {
        BKE_gpencil_stretch_stroke(gps, lmd->length);
      }
      return;
    }
  }
}

/* -------------------------------- */

/* Generic "generateStrokes" callback */
static void generateStrokes(
    GpencilModifierData *md, Depsgraph *depsgraph, Object *ob, bGPDlayer *gpl, bGPDframe *gpf)
{
  BackboneGpencilModifierData *lmd = (BackboneGpencilModifierData *)md;
  bGPDstroke *gps;
  for (gps = gpf->strokes.first; gps; gps = gps->next) {
    BKE_gpencil_stretch_stroke(gps, lmd->length);
  }
}

static void updateDepsgraph(GpencilModifierData *md, const ModifierUpdateDepsgraphContext *ctx)
{
  BackboneGpencilModifierData *lmd = (BackboneGpencilModifierData *)md;
}

static void foreachObjectLink(GpencilModifierData *md,
                              Object *ob,
                              ObjectWalkFunc walk,
                              void *userData)
{
  BackboneGpencilModifierData *mmd = (BackboneGpencilModifierData *)md;
}

GpencilModifierTypeInfo modifierType_Gpencil_Backbone = {
    /* name */ "Backbone Stretch",
    /* structName */ "BackboneGpencilModifierData",
    /* structSize */ sizeof(BackboneGpencilModifierData),
    /* type */ eGpencilModifierTypeType_Gpencil,
    /* flags */ 0,

    /* copyData */ copyData,

    /* deformStroke */ NULL,
    /* generateStrokes */ generateStrokes,
    /* bakeModifier */ bakeModifier,
    /* remapTime */ NULL,

    /* initData */ initData,
    /* freeData */ NULL,
    /* isDisabled */ NULL,
    /* updateDepsgraph */ updateDepsgraph,
    /* dependsOnTime */ NULL,
    /* foreachObjectLink */ foreachObjectLink,
    /* foreachIDLink */ NULL,
    /* foreachTexLink */ NULL,
    /* getDuplicationFactor */ NULL,
};