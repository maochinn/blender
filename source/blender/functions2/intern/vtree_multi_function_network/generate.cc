#include "FN_vtree_multi_function_network_builder.h"
#include "FN_vtree_multi_function_network_generation.h"
#include "FN_multi_functions.h"

#include "BLI_math_cxx.h"
#include "BLI_string_map.h"

namespace FN {

using BLI::float3;
using BLI::StringMap;

static MFDataType get_type_by_socket(const VSocket &vsocket)
{
  StringRef idname = vsocket.idname();

  if (idname == "fn_FloatSocket") {
    return MFDataType::ForSingle<float>();
  }
  else if (idname == "fn_VectorSocket") {
    return MFDataType::ForSingle<float3>();
  }
  else if (idname == "fn_IntegerSocket") {
    return MFDataType::ForSingle<int32_t>();
  }
  else if (idname == "fn_BooleanSocket") {
    return MFDataType::ForSingle<bool>();
  }
  else if (idname == "fn_ObjectSocket") {
    return MFDataType::ForSingle<Object *>();
  }
  else if (idname == "fn_TextSocket") {
    return MFDataType::ForSingle<std::string>();
  }
  else if (idname == "fn_FloatListSocket") {
    return MFDataType::ForVector<float>();
  }
  else if (idname == "fn_VectorListSocket") {
    return MFDataType::ForVector<float3>();
  }
  else if (idname == "fn_IntegerListSocket") {
    return MFDataType::ForVector<int32_t>();
  }
  else if (idname == "fn_BooleanListSocket") {
    return MFDataType::ForVector<bool>();
  }
  else if (idname == "fn_ObjectListSocket") {
    return MFDataType::ForVector<Object *>();
  }
  else if (idname == "fn_TextListSocket") {
    return MFDataType::ForVector<std::string>();
  }

  return MFDataType();
}

static const CPPType &get_cpp_type_by_name(StringRef name)
{
  if (name == "Float") {
    return FN::GET_TYPE<float>();
  }
  else if (name == "Vector") {
    return FN::GET_TYPE<float3>();
  }
  else if (name == "Integer") {
    return FN::GET_TYPE<int32_t>();
  }
  else if (name == "Boolean") {
    return FN::GET_TYPE<bool>();
  }
  else if (name == "Object") {
    return FN::GET_TYPE<Object *>();
  }
  else if (name == "Text") {
    return FN::GET_TYPE<std::string>();
  }

  BLI_assert(false);
  return FN::GET_TYPE<float>();
}

using InsertVNodeFunction = std::function<void(
    VTreeMFNetworkBuilder &builder, OwnedResources &resources, const VNode &vnode)>;
using InsertUnlinkedInputFunction = std::function<MFBuilderOutputSocket &(
    VTreeMFNetworkBuilder &builder, OwnedResources &resources, const VSocket &vsocket)>;
using InsertImplicitConversionFunction =
    std::function<std::pair<MFBuilderInputSocket *, MFBuilderOutputSocket *>(
        VTreeMFNetworkBuilder &builder, OwnedResources &resources)>;

template<typename T, typename... Args>
T &allocate_resource(const char *name, OwnedResources &resources, Args &&... args)
{
  std::unique_ptr<T> value = BLI::make_unique<T>(std::forward<Args>(args)...);
  T &value_ref = *value;
  resources.add(std::move(value), name);
  return value_ref;
}

static void INSERT_vector_math(VTreeMFNetworkBuilder &builder,
                               OwnedResources &resources,
                               const VNode &vnode)
{
  const MultiFunction &fn = allocate_resource<FN::MF_AddFloat3s>("vector math function",
                                                                 resources);
  builder.add_function(fn, {0, 1}, {2}, vnode);
}

static const MultiFunction &get_vectorized_function(
    const MultiFunction &base_function,
    OwnedResources &resources,
    PointerRNA *rna,
    ArrayRef<const char *> is_vectorized_prop_names)
{
  Vector<bool> input_is_vectorized;
  for (const char *prop_name : is_vectorized_prop_names) {
    char state[5];
    RNA_string_get(rna, prop_name, state);
    BLI_assert(STREQ(state, "BASE") || STREQ(state, "LIST"));

    bool is_vectorized = STREQ(state, "LIST");
    input_is_vectorized.append(is_vectorized);
  }

  if (input_is_vectorized.contains(true)) {
    return allocate_resource<FN::MF_SimpleVectorize>(
        "vectorized function", resources, base_function, input_is_vectorized);
  }
  else {
    return base_function;
  }
}

static void INSERT_float_math(VTreeMFNetworkBuilder &builder,
                              OwnedResources &resources,
                              const VNode &vnode)
{
  const MultiFunction &base_fn = allocate_resource<FN::MF_AddFloats>("float math function",
                                                                     resources);
  const MultiFunction &fn = get_vectorized_function(
      base_fn, resources, vnode.rna(), {"use_list__a", "use_list__b"});

  builder.add_function(fn, {0, 1}, {2}, vnode);
}

static void INSERT_combine_vector(VTreeMFNetworkBuilder &builder,
                                  OwnedResources &resources,
                                  const VNode &vnode)
{
  const MultiFunction &base_fn = allocate_resource<FN::MF_CombineVector>("combine vector function",
                                                                         resources);
  const MultiFunction &fn = get_vectorized_function(
      base_fn, resources, vnode.rna(), {"use_list__x", "use_list__y", "use_list__z"});
  builder.add_function(fn, {0, 1, 2}, {3}, vnode);
}

static void INSERT_separate_vector(VTreeMFNetworkBuilder &builder,
                                   OwnedResources &resources,
                                   const VNode &vnode)
{
  const MultiFunction &base_fn = allocate_resource<FN::MF_SeparateVector>(
      "separate vector function", resources);
  const MultiFunction &fn = get_vectorized_function(
      base_fn, resources, vnode.rna(), {"use_list__vector"});
  builder.add_function(fn, {0}, {1, 2, 3}, vnode);
}

static void INSERT_list_length(VTreeMFNetworkBuilder &builder,
                               OwnedResources &resources,
                               const VNode &vnode)
{
  char *type_name = RNA_string_get_alloc(vnode.rna(), "active_type", nullptr, 0);
  const CPPType &type = get_cpp_type_by_name(type_name);
  MEM_freeN(type_name);

  const MultiFunction &fn = allocate_resource<FN::MF_ListLength>(
      "list length function", resources, type);
  builder.add_function(fn, {0}, {1}, vnode);
}

static void INSERT_get_list_element(VTreeMFNetworkBuilder &builder,
                                    OwnedResources &resources,
                                    const VNode &vnode)
{
  char *type_name = RNA_string_get_alloc(vnode.rna(), "active_type", nullptr, 0);
  const CPPType &type = get_cpp_type_by_name(type_name);
  MEM_freeN(type_name);

  const MultiFunction &fn = allocate_resource<FN::MF_GetListElement>(
      "get list element function", resources, type);
  builder.add_function(fn, {0, 1, 2}, {3}, vnode);
}

static MFBuilderOutputSocket &build_pack_list_node(VTreeMFNetworkBuilder &builder,
                                                   OwnedResources &resources,
                                                   const VNode &vnode,
                                                   const CPPType &base_type,
                                                   const char *prop_name,
                                                   uint start_index)
{
  Vector<bool> input_is_list;
  RNA_BEGIN (vnode.rna(), itemptr, prop_name) {
    int state = RNA_enum_get(&itemptr, "state");
    if (state == 0) {
      /* single value case */
      input_is_list.append(false);
    }
    else if (state == 1) {
      /* list case */
      input_is_list.append(true);
    }
    else {
      BLI_assert(false);
    }
  }
  RNA_END;

  uint input_amount = input_is_list.size();
  uint output_param_index = (input_amount > 0 && input_is_list[0]) ? 0 : input_amount;

  const MultiFunction &fn = allocate_resource<FN::MF_PackList>(
      "pack list function", resources, base_type, input_is_list);
  MFBuilderFunctionNode &node = builder.add_function(
      fn, IndexRange(input_amount).as_array_ref(), {output_param_index});

  for (uint i = 0; i < input_amount; i++) {
    builder.map_sockets(vnode.input(start_index + i), *node.inputs()[i]);
  }

  return *node.outputs()[0];
}

static void INSERT_pack_list(VTreeMFNetworkBuilder &builder,
                             OwnedResources &resources,
                             const VNode &vnode)
{
  char *type_name = RNA_string_get_alloc(vnode.rna(), "active_type", nullptr, 0);
  const CPPType &type = get_cpp_type_by_name(type_name);
  MEM_freeN(type_name);

  MFBuilderOutputSocket &packed_list_socket = build_pack_list_node(
      builder, resources, vnode, type, "variadic", 0);
  builder.map_sockets(vnode.output(0), packed_list_socket);
}

static void INSERT_object_location(VTreeMFNetworkBuilder &builder,
                                   OwnedResources &resources,
                                   const VNode &vnode)
{
  const MultiFunction &fn = allocate_resource<FN::MF_ObjectWorldLocation>(
      "object location function", resources);
  builder.add_function(fn, {0}, {1}, vnode);
}

static void INSERT_text_length(VTreeMFNetworkBuilder &builder,
                               OwnedResources &resources,
                               const VNode &vnode)
{
  const MultiFunction &fn = allocate_resource<FN::MF_TextLength>("text length function",
                                                                 resources);
  builder.add_function(fn, {0}, {1}, vnode);
}

static StringMap<InsertVNodeFunction> get_node_inserters()
{
  StringMap<InsertVNodeFunction> inserters;
  inserters.add_new("fn_FloatMathNode", INSERT_float_math);
  inserters.add_new("fn_VectorMathNode", INSERT_vector_math);
  inserters.add_new("fn_CombineVectorNode", INSERT_combine_vector);
  inserters.add_new("fn_SeparateVectorNode", INSERT_separate_vector);
  inserters.add_new("fn_ListLengthNode", INSERT_list_length);
  inserters.add_new("fn_PackListNode", INSERT_pack_list);
  inserters.add_new("fn_GetListElementNode", INSERT_get_list_element);
  inserters.add_new("fn_ObjectTransformsNode", INSERT_object_location);
  inserters.add_new("fn_TextLengthNode", INSERT_text_length);
  return inserters;
}

static MFBuilderOutputSocket &INSERT_vector_socket(VTreeMFNetworkBuilder &builder,
                                                   OwnedResources &resources,
                                                   const VSocket &vsocket)
{
  float3 value;
  RNA_float_get_array(vsocket.rna(), "value", value);

  const MultiFunction &fn = allocate_resource<FN::MF_ConstantValue<float3>>(
      "vector socket", resources, value);
  MFBuilderFunctionNode &node = builder.add_function(fn, {}, {0});
  return *node.outputs()[0];
}

static MFBuilderOutputSocket &INSERT_float_socket(VTreeMFNetworkBuilder &builder,
                                                  OwnedResources &resources,
                                                  const VSocket &vsocket)
{
  float value = RNA_float_get(vsocket.rna(), "value");

  const MultiFunction &fn = allocate_resource<FN::MF_ConstantValue<float>>(
      "float socket", resources, value);
  MFBuilderFunctionNode &node = builder.add_function(fn, {}, {0});
  return *node.outputs()[0];
}

static MFBuilderOutputSocket &INSERT_int_socket(VTreeMFNetworkBuilder &builder,
                                                OwnedResources &resources,
                                                const VSocket &vsocket)
{
  int value = RNA_int_get(vsocket.rna(), "value");

  const MultiFunction &fn = allocate_resource<FN::MF_ConstantValue<int>>(
      "int socket", resources, value);
  MFBuilderFunctionNode &node = builder.add_function(fn, {}, {0});
  return *node.outputs()[0];
}

static MFBuilderOutputSocket &INSERT_object_socket(VTreeMFNetworkBuilder &builder,
                                                   OwnedResources &resources,
                                                   const VSocket &vsocket)
{
  Object *value = (Object *)RNA_pointer_get(vsocket.rna(), "value").data;

  const MultiFunction &fn = allocate_resource<FN::MF_ConstantValue<Object *>>(
      "object socket", resources, value);
  MFBuilderFunctionNode &node = builder.add_function(fn, {}, {0});
  return *node.outputs()[0];
}

static MFBuilderOutputSocket &INSERT_text_socket(VTreeMFNetworkBuilder &builder,
                                                 OwnedResources &resources,
                                                 const VSocket &vsocket)
{
  char *value = RNA_string_get_alloc(vsocket.rna(), "value", nullptr, 0);
  std::string text = value;
  MEM_freeN(value);

  const MultiFunction &fn = allocate_resource<FN::MF_ConstantValue<std::string>>(
      "text socket", resources, text);
  MFBuilderFunctionNode &node = builder.add_function(fn, {}, {0});
  return *node.outputs()[0];
}

template<typename T>
static MFBuilderOutputSocket &INSERT_empty_list_socket(VTreeMFNetworkBuilder &builder,
                                                       OwnedResources &resources,
                                                       const VSocket &UNUSED(vsocket))
{
  const MultiFunction &fn = allocate_resource<FN::MF_EmptyList<T>>("empty list socket", resources);
  MFBuilderFunctionNode &node = builder.add_function(fn, {}, {0});
  return *node.outputs()[0];
}

static StringMap<InsertUnlinkedInputFunction> get_unlinked_input_inserter()
{
  StringMap<InsertUnlinkedInputFunction> inserters;
  inserters.add_new("fn_VectorSocket", INSERT_vector_socket);
  inserters.add_new("fn_FloatSocket", INSERT_float_socket);
  inserters.add_new("fn_IntegerSocket", INSERT_int_socket);
  inserters.add_new("fn_ObjectSocket", INSERT_object_socket);
  inserters.add_new("fn_TextSocket", INSERT_text_socket);
  inserters.add_new("fn_VectorListSocket", INSERT_empty_list_socket<float3>);
  inserters.add_new("fn_FloatListSocket", INSERT_empty_list_socket<float>);
  inserters.add_new("fn_IntegerListSocket", INSERT_empty_list_socket<int32_t>);
  inserters.add_new("fn_ObjectListSocket", INSERT_empty_list_socket<Object *>);
  inserters.add_new("fn_TextListSocket", INSERT_empty_list_socket<std::string>);
  return inserters;
}

template<typename FromT, typename ToT>
static std::pair<MFBuilderInputSocket *, MFBuilderOutputSocket *> INSERT_convert(
    VTreeMFNetworkBuilder &builder, OwnedResources &resources)
{
  const MultiFunction &fn = allocate_resource<FN::MF_Convert<FromT, ToT>>("converter function",
                                                                          resources);
  MFBuilderFunctionNode &node = builder.add_function(fn, {0}, {1});
  return {node.inputs()[0], node.outputs()[0]};
}

template<typename FromT, typename ToT>
static std::pair<MFBuilderInputSocket *, MFBuilderOutputSocket *> INSERT_convert_list(
    VTreeMFNetworkBuilder &builder, OwnedResources &resources)
{
  const MultiFunction &fn = allocate_resource<FN::MF_ConvertList<FromT, ToT>>(
      "convert list function", resources);
  MFBuilderFunctionNode &node = builder.add_function(fn, {0}, {1});
  return {node.inputs()[0], node.outputs()[0]};
}

template<typename T>
static std::pair<MFBuilderInputSocket *, MFBuilderOutputSocket *> INSERT_element_to_list(
    VTreeMFNetworkBuilder &builder, OwnedResources &resources)
{
  const MultiFunction &fn = allocate_resource<FN::MF_SingleElementList<T>>(
      "single element list function", resources);
  MFBuilderFunctionNode &node = builder.add_function(fn, {0}, {1});
  return {node.inputs()[0], node.outputs()[0]};
}

static Map<std::pair<std::string, std::string>, InsertImplicitConversionFunction>
get_conversion_inserters()
{
  Map<std::pair<std::string, std::string>, InsertImplicitConversionFunction> inserters;

  inserters.add_new({"fn_IntegerSocket", "fn_FloatSocket"}, INSERT_convert<int, float>);
  inserters.add_new({"fn_FloatSocket", "fn_IntegerSocket"}, INSERT_convert<float, int>);

  inserters.add_new({"fn_FloatSocket", "fn_BooleanSocket"}, INSERT_convert<float, bool>);
  inserters.add_new({"fn_BooleanSocket", "fn_FloatSocket"}, INSERT_convert<bool, float>);

  inserters.add_new({"fn_IntegerSocket", "fn_BooleanSocket"}, INSERT_convert<int, bool>);
  inserters.add_new({"fn_BooleanSocket", "fn_IntegerSocket"}, INSERT_convert<bool, int>);

  inserters.add_new({"fn_IntegerListSocket", "fn_FloatListSocket"},
                    INSERT_convert_list<int, float>);
  inserters.add_new({"fn_FloatListSocket", "fn_IntegerListSocket"},
                    INSERT_convert_list<float, int>);

  inserters.add_new({"fn_FloatListSocket", "fn_BooleanListSocket"},
                    INSERT_convert_list<float, bool>);
  inserters.add_new({"fn_BooleanListSocket", "fn_FloatListSocket"},
                    INSERT_convert_list<bool, float>);

  inserters.add_new({"fn_IntegerListSocket", "fn_BooleanListSocket"},
                    INSERT_convert_list<int, bool>);
  inserters.add_new({"fn_BooleanListSocket", "fn_IntegerListSocket"},
                    INSERT_convert_list<bool, int>);

  inserters.add_new({"fn_IntegerSocket", "fn_IntegerListSocket"}, INSERT_element_to_list<int>);
  inserters.add_new({"fn_FloatSocket", "fn_FloatListSocket"}, INSERT_element_to_list<float>);
  inserters.add_new({"fn_BooleanSocket", "fn_BooleanListSocket"}, INSERT_element_to_list<bool>);

  return inserters;
}

static bool insert_nodes(VTreeMFNetworkBuilder &builder, OwnedResources &resources)
{
  const VirtualNodeTree &vtree = builder.vtree();
  auto inserters = get_node_inserters();

  for (const VNode *vnode : vtree.nodes()) {
    StringRef idname = vnode->idname();
    InsertVNodeFunction *inserter = inserters.lookup_ptr(idname);

    if (inserter != nullptr) {
      (*inserter)(builder, resources, *vnode);
      BLI_assert(builder.data_sockets_of_vnode_are_mapped(*vnode));
    }
    else if (builder.has_data_sockets(*vnode)) {
      builder.add_dummy(*vnode);
    }
  }

  return true;
}

static bool insert_links(VTreeMFNetworkBuilder &builder, OwnedResources &resources)
{
  auto conversion_inserters = get_conversion_inserters();

  for (const VInputSocket *to_vsocket : builder.vtree().all_input_sockets()) {
    ArrayRef<const VOutputSocket *> origins = to_vsocket->linked_sockets();
    if (origins.size() != 1) {
      continue;
    }

    if (!builder.is_data_socket(*to_vsocket)) {
      continue;
    }

    const VOutputSocket *from_vsocket = origins[0];
    if (!builder.is_data_socket(*from_vsocket)) {
      return false;
    }

    MFBuilderOutputSocket &from_socket = builder.lookup_socket(*from_vsocket);
    MFBuilderInputSocket &to_socket = builder.lookup_socket(*to_vsocket);

    if (from_socket.type() == to_socket.type()) {
      builder.add_link(from_socket, to_socket);
    }
    else {
      InsertImplicitConversionFunction *inserter = conversion_inserters.lookup_ptr(
          {from_vsocket->idname(), to_vsocket->idname()});
      if (inserter == nullptr) {
        return false;
      }
      auto new_sockets = (*inserter)(builder, resources);
      builder.add_link(from_socket, *new_sockets.first);
      builder.add_link(*new_sockets.second, to_socket);
    }
  }

  return true;
}

static bool insert_unlinked_inputs(VTreeMFNetworkBuilder &builder, OwnedResources &resources)
{
  Vector<const VInputSocket *> unlinked_data_inputs;
  for (const VInputSocket *vsocket : builder.vtree().all_input_sockets()) {
    if (builder.is_data_socket(*vsocket)) {
      if (!builder.is_input_linked(*vsocket)) {
        unlinked_data_inputs.append(vsocket);
      }
    }
  }

  auto inserters = get_unlinked_input_inserter();

  for (const VInputSocket *vsocket : unlinked_data_inputs) {
    InsertUnlinkedInputFunction *inserter = inserters.lookup_ptr(vsocket->idname());

    if (inserter == nullptr) {
      return false;
    }
    MFBuilderOutputSocket &from_socket = (*inserter)(builder, resources, *vsocket);
    MFBuilderInputSocket &to_socket = builder.lookup_socket(*vsocket);
    builder.add_link(from_socket, to_socket);
  }

  return true;
}

std::unique_ptr<VTreeMFNetwork> generate_vtree_multi_function_network(const VirtualNodeTree &vtree,
                                                                      OwnedResources &resources)
{
  Vector<MFDataType> type_by_vsocket{vtree.socket_count()};
  for (const VSocket *vsocket : vtree.all_sockets()) {
    MFDataType data_type = get_type_by_socket(*vsocket);
    type_by_vsocket[vsocket->id()] = data_type;
  }

  VTreeMFNetworkBuilder builder(vtree, std::move(type_by_vsocket));
  if (!insert_nodes(builder, resources)) {
    BLI_assert(false);
  }
  if (!insert_links(builder, resources)) {
    BLI_assert(false);
  }
  if (!insert_unlinked_inputs(builder, resources)) {
    BLI_assert(false);
  }

  auto vtree_network = builder.build();
  return vtree_network;
}

}  // namespace FN