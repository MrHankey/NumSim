#ifdef __nope__
#define __kernel
#define __global
#define __local
#define CLK_GLOBAL_MEM_FENCE
#define CLK_LOCAL_MEM_FENCE
#define float4
#endif

__kernel void sor(	__global float *grid,
					__global const float *rhs,
					__global float *resGrid,
					__global float *h_square,
					__global float *h_square_inv
					//__local float* localGrid,
					//__local float* localRHS
				)
 {
	//printf("hi");
	float omega = 1.0f;
	float hs = (*h_square);
	float hsi = (*h_square_inv);

	float fac = 0.25f;
    //float fac2 = 0.5f*((hs*hs)/(hs + hs));

    int g_x = get_global_id(0) + 1;
    int g_y = get_global_id(1) + 1;

    int gridSize = get_global_size(0) + 2;
    int gridSizeInterior = get_global_size(0);
    int idx = g_y*gridSize + g_x;
    float res = 1.0f;
    float eps = 0.001f;
    int it = 1;
    int itMax = 100;

    /*if (g_x == 1 && g_y == 1)
    	printf("while \n");*/

	while(res>eps && it <= itMax){

		//if (g_x == 1 && g_y == 1)
			//printf("it: %d \n", it);

		barrier(0);


	    float p = grid[idx];
	    float p_l = grid[idx-1];
	    float p_r = grid[idx+1];
	    float p_u = grid[idx+gridSize];
	    float p_d = grid[idx-gridSize];
	    float p_rhs = rhs[idx];


	    //red cycle
	    if ((((g_x + g_y)) & 1) == 0)
	    {

	    	float p_star = fac*(p_d + p_l + p_r + p_u - hs*p_rhs);
	    	//float p_star = fac2*((p_d + p_u)*hsi + (p_l + p_r)*hsi - p_rhs);
	    	grid[idx] = (1.0f-omega)*p + omega*p_star;
	    	//printf("red: %d %d pstar: %f result: %f \n", g_x, g_y, p_star, grid[idx]);
	    }

	    barrier(CLK_LOCAL_MEM_FENCE);

	    //black cycle
	    if (((g_x + g_y) & 1) != 0)
	    {
	    	float p_star = fac*(p_d + p_l + p_r + p_u - hs*p_rhs);
	    	//float p_star = fac2*((p_d + p_u)*hsi + (p_l + p_r)*hsi - p_rhs);
	    	grid[idx] = (1.0f-omega)*p + omega*p_star;
	    	//printf("black: %d %d pstar: %f result: %f \n", g_x, g_y, p_star, grid[idx]);
	    }

	    p = grid[idx];
	    p_l = grid[idx-1];
	    p_r = grid[idx+1];
	    p_u = grid[idx+gridSize];
	    p_d = grid[idx-gridSize];

	    //float dxx = (p_l - 2*p + p_r) * hsi;
		//float dyy = (p_u - 2*p + p_d) * hsi;

		//resGrid[idx] = fabs( dxx + dyy - p_rhs);
	    resGrid[idx] = fabs(p_l + p_r + p_u + p_d - 6.0f*p - hs*p_rhs);

	    //if (resGrid[idx] != 0.0)
			//		printf("idx: %d res_start: %f \n", g_x, resGrid[idx]);

		//barrier(0);


		/*int teiler = 2;
		for (int i = 0; i < log2((float)gridSize); i++){
			if( ((g_x % teiler) == 1) && (g_x + teiler) < (gridSize)){
				resGrid[idx]=resGrid[idx]+resGrid[idx+(teiler/2)];

				//if (g_x == 1 && g_y == 100)
					//printf("res_x: %f \n", resGrid[idx]);

				teiler = teiler*2;
			}
			else
			{
				break;
			}
			barrier(0);
		}

		teiler = 2;


		for (int i = 0; i < log2((float)gridSize); i++){
			if( (g_x == 1 && (g_y % teiler) == 1) && (g_y+teiler) < (gridSize)){
				resGrid[idx] += resGrid[idx+((teiler/2)*gridSize)];
				if (g_x == 1 && g_y == 1)
					//printf("teiler: %d res: %f \n", teiler, resGrid[idx]);

				teiler=teiler*2;
			}
			else
			{
				break;
			}
			barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
		}*/



		// Update Boundaries
		if (g_x == 1){
			grid[idx-1] = grid[idx];
		}
		if (g_x == gridSizeInterior){
			grid[idx+1] = grid[idx];
		}
		if (g_y == 1){
			grid[idx-gridSize] = grid[idx];
		}
		if (g_y == gridSizeInterior){
			grid[idx+gridSize] = grid[idx];
		}

		//barrier(0);
		//res = resGrid[gridSize+1]/(gridSize*gridSize);

		//if (g_x == 1 && g_y == 1)
			//printf("res: %f \n", resGrid[0]);

		it++;
	}

}
