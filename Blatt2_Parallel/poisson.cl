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
								__global float *resGrid,
								__global float *h_square,
								__global float *h_square_inv,
								__local float* scratch
						)
 {
 
 
 	int strideSize = 1;
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
 
    // Get the index of the current element to be processed
    int i = get_global_id(0)*strideSize + 1;
    int j = get_global_id(1)*strideSize + 1;
    
    for ( int i_loc = 0; i_loc < strideSize+2; i_loc++)
    {
    	for (int j_loc = 0; j_loc < strideSize+2; j_loc++)
    	{
    		//scratch[j_loc*(strideSize+2) + i_loc] = oldGrid[(j-1+j_loc)*strideFac + i - 1 + i_loc];
    		scratch[j_loc*(strideSize+2) + i_loc] = oldGrid[(j-1+j_loc)*strideFac + (i-1) + i_loc];
    	}
    }

    
    
    //float h_square = h*h;
    //float h_square_inv = 1/h_square;
    

#pragma unroll
    for ( int j_off = 0; j_off < strideSize; j_off++ )
    {

    	base_index = (j + j_off)*strideFac + i;
    	/*base_index_d = base_index - strideFac;
    	base_index_u = base_index + strideFac;*/
#pragma unroll
    	for ( int i_off = 0; i_off < strideSize; i_off++ )
    	{
    
		    index_p = (j_off+1)*(strideSize + 2) + i_off + 1;
		    index_p_d = index_p - strideSize - 2;
		    index_p_u = index_p + strideSize + 2;
		    index_p_r = index_p + 1;
		    index_p_l = index_p - 1;
		    
		    newPressure = 0.25f*(scratch[index_p_r] + scratch[index_p_l] + scratch[index_p_u] + scratch[index_p_d] - (*h_square)*rhs[base_index+i_off]);
		    
		    //if ( i == 5 && j == 5)
		    	//printf("dx: %lf, rhs: %lf, newP: %lf \n", h, rhs[index_p], newPressure);
		    	  
		    newGrid[base_index+i_off] = newPressure;
		    
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
    
		    dxx = (newGrid[index_p_r]-2*newGrid[index_p]+newGrid[index_p_l])*(*h_square_inv);
			dyy = (newGrid[index_p_u]-2*newGrid[index_p]+newGrid[index_p_d])*(*h_square_inv);
    
    		resGrid[index_p] = fabs( dxx + dyy - rhs[index_p]);
		}
	}
    		 
}
