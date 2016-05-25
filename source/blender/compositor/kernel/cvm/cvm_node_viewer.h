CVM_float4_node_start(viewer)
  float4 color_result = make_float4(0.0,0.0,0.0,0.0);
  float add_sample = 1.0f/global.subpixel_samples_xy;
  float total_weight = 0.;
  float inset = 1.0f/(2*global.subpixel_samples_xy);
  for (float addx = inset; addx < 1.0 ; addx += add_sample ) {
    for (float addy = inset; addy < 1.0 ; addy += add_sample ) {
      float weight = 1.;
      float2 sample_coord = make_float2(xy.x+addx, xy.y+addy);
      float4 new_result = CVM_float4_node_call(0, sample_coord, __cvm_in_ref weight);
      color_result += new_result*weight;
      total_weight += weight;
    }
  }

  return color_result / total_weight;

CVM_float4_node_end
