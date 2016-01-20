#ifdef __nope__
#define __kernel
#define __global
#define __local
#define CLK_GLOBAL_MEM_FENCE
#define CLK_LOCAL_MEM_FENCE
#define float4
#define INFINITY 0
#endif

__kernel
void reduce_sum(__global float* buffer,
            __local float* scratch,
            __global unsigned int* pLength,
            __global float* result) {

  __const int length = *pLength;
  int global_index = get_global_id(0);
  float accumulator = 0;
  // Loop sequentially over chunks of input vector
  while (global_index < length) {
    float element = buffer[global_index];
    accumulator += element;
    global_index += get_global_size(0);
  }

  // Perform parallel reduction
  int local_index = get_local_id(0);
  scratch[local_index] = accumulator;
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int offset = get_local_size(0) / 2;
      offset > 0;
      offset = offset / 2) {
    if (local_index < offset) {
      float other = scratch[local_index + offset];
      float mine = scratch[local_index];
      scratch[local_index] += other;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  if (local_index == 0) {
    result[get_group_id(0)] = scratch[0];
  }
}

__kernel
void reduce_max(__global float* buffer,
            __local float* scratch,
            __global unsigned int* pLength,
            __global float* result) {

  __const int length = *pLength;
  int global_index = get_global_id(0);
  float accumulator = 0;
  // Loop sequentially over chunks of input vector
  while (global_index < length) {
    float element = buffer[global_index];
    accumulator = (accumulator > element) ? accumulator : element;
    global_index += get_global_size(0);
  }

  // Perform parallel reduction
  int local_index = get_local_id(0);
  scratch[local_index] = accumulator;
  barrier(CLK_LOCAL_MEM_FENCE);
  for(int offset = get_local_size(0) / 2;
      offset > 0;
      offset = offset / 2) {
    if (local_index < offset) {
      float other = scratch[local_index + offset];
      float mine = scratch[local_index];
      scratch[local_index] = (mine > other) ? mine : other;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  if (local_index == 0) {
    result[get_group_id(0)] = scratch[0];
  }
}
