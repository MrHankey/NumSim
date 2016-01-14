#ifdef __nope__
#define __kernel
#define __global
#define __local
#define CLK_GLOBAL_MEM_FENCE
#define CLK_LOCAL_MEM_FENCE
#define float4
#endif

__kernel void sor(	__global const float *grid,
					__global const float *rhs,
					__global float *resGrid,
					__global float *h_square,
					__global float *h_square_inv,
					__local float* localGrid,
					__local float* localRHS
				)
 {
	float omega = 1.7f;

    float fac = 0.25f; // 1/4;

    int g_x = get_global_id(0) + 1;
    int g_y = get_global_id(1) + 1;
    int gridSize = get_global_size(0) + 2;
    int locGridSize = get_local_size(0) + 2;
    int l_x = get_local_id(0) + 1;
    int l_y = get_local_id(1) + 1;
    int idx = g_y*gridSize + g_x;
    int localIdx = l_y*locGridSize + l_x;

    int wg_x = get_group_id(0)*get_local_size(0)+1;
    int wg_y = get_group_id(1)*get_local_size(1)+1;

    for (int i = 0; i < locGridSize; i++)
    {
    	for (int j = 0; j < locGridSize; j++)
    	{
    		localGrid[j*locGridSize + i] = grid	[(wg_y-1 + j)*gridSize + wg_x-1 + i];
    		localRHS[j*locGridSize + i] = rhs	[(wg_y-1 + j)*gridSize + wg_x-1 + i];
    	}
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    //red cycle
    if ((g_x + g_y)&1 == 0)
    {
    	float p_star = fac*(localGrid[localIdx-locGridSize] + localGrid[localIdx-1] + localGrid[localIdx+1] + localGrid[localIdx+locGridSize] - h_square*localRHS[localIdx]);
    	localGrid[localIdx] = (1-omega)*localGrid[localIdx] + omega*p_star;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    //black cycle
    if ((g_x + g_y)&1 != 0)
    {
    	float p_star = fac*(localGrid[localIdx-locGridSize] + localGrid[localIdx-1] + localGrid[localIdx+1] + localGrid[localIdx+locGridSize] - h_square*localRHS[localIdx]);
    	localGrid[localIdx] = (1-omega)*localGrid[localIdx] + omega*p_star;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

	resGrid[idx] = localGrid[localIdx-locGridSize] + localGrid[localIdx-1] + localGrid[localIdx+1] + localGrid[localIdx+locGridSize] - 4*localGrid[localIdx] - h_square*localRHS[localIdx];

    for (int i = 0; i < locGridSize; i++)
    {
    	for (int j = 0; j < locGridSize; j++)
    	{
    		grid[(g_y-1 + j)*gridSize + (g_x-1) + i] = localGrid[j*locGridSize + i];
    	}
    }
}
