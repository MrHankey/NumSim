__kernel void poisson_jacobi(	__global const float *oldGrid,
								__global const float *rhs,
								__global float *newGrid,
								__global float *dx
						)
 {
 
 
    // Get the index of the current element to be processed
    int i = get_global_id(0) + 1;
    int j = get_global_id(1) + 1;
    
    int index_p = j*get_global_size(0) + i;
    int index_p_d = (j-1)*get_global_size(0) + i;
    int index_p_u = (j+1)*get_global_size(0) + i;
    int index_p_r = j*get_global_size(0) + i + 1;
    int index_p_l = j*get_global_size(0) + i - 1;
    
    float h = *dx;
    float newPressure = 0.25*(oldGrid[index_p_r] + oldGrid[index_p_l] + oldGrid[index_p_u] + oldGrid[index_p_d] - (h*h)*rhs[index_p]);
    
    //if ( i == 2 && j == 2)
    	//printf("dx: %lf, rhs: %lf, newP: %lf \n", h, rhs[index_p], newPressure);
    	  
    newGrid[index_p] = newPressure;
}

__kernel void reduction_vector(__global float4* data,__local float4* partial_sums, __global float* output) 
{
    int lid = get_local_id(0);
    int group_size = get_local_size(0);
    partial_sums[lid] = data[get_global_id(0)];
    barrier(CLK_LOCAL_MEM_FENCE);

    for(int i = group_size/2; i>0; i >>= 1) {
        if(lid < i) {
            partial_sums[lid] += partial_sums[lid + i];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if(lid == 0) {
        output[get_group_id(0)] = dot(partial_sums[0], (float4)(1.0f));
    }
}