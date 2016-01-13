#ifdef __nope__
#define __kernel
#define __global
#define __local
#define CLK_GLOBAL_MEM_FENCE
#define CLK_LOCAL_MEM_FENCE
#define float4
#endif

__kernel void poisson_jacobi(	__global const float *oldGrid,
								__global const float *rhs,
								__global float *newGrid,
								__global float *dx,
								__global float *resGrid
						)
 {
 
 
 	int strideSize = 4;
 
    // Get the index of the current element to be processed
    int i = get_global_id(0)*strideSize + 1;
    int j = get_global_id(1)*strideSize + 1;
    
    
    
    float h = *dx;
    float h_square = h*h;
    float h_square_inv = 1/h_square;
    int strideFac = (get_global_size(0)*strideSize+2);
    int base_index = 0;
    int base_index_d = 0;
    int base_index_u = 0;
    int index_p = 0;
    int index_p_u = 0;
    int index_p_d = 0;
    int index_p_r = 0;
    int index_p_l = 0;
    float newPressure = 0.0f;
    float dxx = 0.0f;
    float dyy = 0.0f;
    
#pragma unroll
    for ( int j_off = 0; j_off < strideSize; j_off++ )
    {

    	base_index = (j + j_off)*strideFac + i;
    	base_index_d = base_index - strideFac;
    	base_index_u = base_index + strideFac;
#pragma unroll
    	for ( int i_off = 0; i_off < strideSize; i_off++ )
    	{
    
		    index_p = base_index + i_off;
		    index_p_d = base_index_d + i_off;
		    index_p_u = base_index_u + i_off;
		    index_p_r = index_p + 1;
		    index_p_l = index_p - 1;
		    
		    newPressure = 0.25f*(oldGrid[index_p_r] + oldGrid[index_p_l] + oldGrid[index_p_u] + oldGrid[index_p_d] - h_square*rhs[index_p]);
		    
		    //if ( i == 5 && j == 5)
		    	//printf("dx: %lf, rhs: %lf, newP: %lf \n", h, rhs[index_p], newPressure);
		    	  
		    newGrid[index_p] = newPressure;
		    
		    /*float dxx = (newGrid[index_p_r]-2*newGrid[index_p]+newGrid[index_p_l])*h_square_inv;
			float dyy = (newGrid[index_p_u]-2*newGrid[index_p]+newGrid[index_p_d])*h_square_inv;
    
    		resGrid[index_p] = fabs( dxx + dyy - rhs[index_p]);*/
		}
	}
    
    barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
    
#pragma unroll
    for ( int j_off = 0; j_off < strideSize; j_off++ )
    {
    	base_index = (j + j_off)*strideFac + i;
    	base_index_d = base_index - strideFac;
    	base_index_u = base_index + strideFac;

#pragma unroll
    	for ( int i_off = 0; i_off < strideSize; i_off++ )
    	{
    		index_p = base_index + i_off;
		    index_p_d = base_index_d + i_off;
		    index_p_u = base_index_u + i_off;
		    index_p_r = index_p + 1;
		    index_p_l = index_p - 1;
    
		    dxx = (newGrid[index_p_r]-2*newGrid[index_p]+newGrid[index_p_l])*h_square_inv;
			dyy = (newGrid[index_p_u]-2*newGrid[index_p]+newGrid[index_p_d])*h_square_inv;
    
    		resGrid[index_p] = fabs( dxx + dyy - rhs[index_p]);
		}
	}
    		 
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
