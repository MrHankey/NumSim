#ifdef __nope__
#define __kernel
#define __global
#define __local
#define CLK_GLOBAL_MEM_FENCE
#define CLK_LOCAL_MEM_FENCE
#define float4
#endif

float d_l(float left, float mid, float hi) {
	//return (_data[it]-_data[it.Left()])/_geom->Mesh()[0];
	return (mid-left)*hi;
}

float d_r(float right, float mid, float hi) {
	//return (_data[it.Right()]-_data[it])/_geom->Mesh()[0];
	return (right-mid)*hi;
}

float ddx(float right, float left, float mid, float hsi) {
	//return (_data[it.Top()]-2*_data[it]+_data[it.Down()])/_geom->Mesh()[1]/_geom->Mesh()[1];
	return (right-2*mid+left)*hsi;
}

float udu_x(float u_r, float u_l, float u, float hi) {
	//real_t A = (_data[it]+_data[it.Right()])/2.0;
	//real_t B = (_data[it.Left()]+_data[it])/2.0;
	//real_t C = (_data[it]-_data[it.Right()])/2.0;
	//real_t D = (_data[it.Left()]-_data[it])/2.0;
	//return (( A*A - B*B) + alpha*(fabs(A)*C-fabs(B)*D) )/_geom->Mesh()[0];
	float A = (u+u_r)*0.5f;
	float B = (u_l+u)*0.5f;
	float C = (u-u_r)*0.5f;
	float D = (u_l-u)*0.5f;

	//TODO ALPHA
	return (( A*A - B*B) + 0.9f*(fabs(A)*C-fabs(B)*D) )*hi;
}

float vdv_y(float v_t, float v_d, float v, float hi) {
	//real_t A = (_data[it]+_data[it.Right()])/2.0;
	//real_t B = (_data[it.Left()]+_data[it])/2.0;
	//real_t C = (_data[it]-_data[it.Right()])/2.0;
	//real_t D = (_data[it.Left()]-_data[it])/2.0;
	//return (( A*A - B*B) + alpha*(fabs(A)*C-fabs(B)*D) )/_geom->Mesh()[0];
	float A = (v+v_t)*0.5f;
	float B = (v_d+v)*0.5f;
	float C = (v-v_t)*0.5f;
	float D = (v_d-v)*0.5f;

	//TODO ALPHA
	return (( A*A - B*B) + 0.9f*(fabs(A)*C-fabs(B)*D) )*hi;
}

float vdu_y(float u_t, float u_d, float u, float v_r, float v_d, float v_dr, float v, float hi) {
	float A = (v + v_r)*0.5f;
	float B = (u + u_t)*0.5f;
	float C = (v_d + v_dr)*0.5f;
	float D = (u + u_d)*0.5f;
	float E = (u - u_t)*0.5f;
	float F = (u_d - u)*0.5f;

	//TODO ALPHA
	return (( A * B - C * D) + 0.9f * ( fabs(A)*E-fabs(C)*F )) * hi;

	//TODO CHECK SYMMETRY
}

float udv_x(float v_r, float v_l, float v, float u_t, float u_l, float u_lt, float u, float hi) {
	float A = (u + u_t)*0.5f;
	float B = (v + v_r)*0.5f;
	float C = (u_l + u_lt)*0.5f;
	float D = (v + v_l)*0.5f;
	float E = (v - v_r)*0.5f;
	float F = (v_l - v)*0.5f;

	//TODO ALPHA
	return (( A * B - C * D) + 0.9f * ( fabs(A)*E-fabs(C)*F )) * hi;

	//TODO CHECK SYMMETRY
}

__kernel void sor(	__global float *grid,
					__global const float *rhs,
					__global float *resGrid,
					__global float *h_square,
					__global float *h_square_inv
					//__local float* localGrid,
					//__local float* localRHS
				)
 {
	float omega = 1.7f;
	float hs = (*h_square);
	float hsi = (*h_square_inv);

	float fac = 0.25f;
    //float fac2 = 0.5f*((hs*hs)/(hs + hs));

    int g_x = get_global_id(0) + 1;
    int g_y = get_global_id(1) + 1;

    int gridSizeInterior = get_global_size(0);
    int gridSize = gridSizeInterior + 2;

    int idx = g_y*gridSize + g_x;


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

    p = grid[idx];
    p_l = grid[idx-1];
    p_r = grid[idx+1];
    p_u = grid[idx+gridSize];
    p_d = grid[idx-gridSize];

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

    float dxx = (p_l - 2.0f*p + p_r) * hsi;
	float dyy = (p_u - 2.0f*p + p_d) * hsi;

	resGrid[idx] = fabs( dxx + dyy - p_rhs);
    //resGrid[idx] = fabs(p_l + p_r + p_u + p_d - 6.0f*p - hs*p_rhs);

    if (g_x == 1){
		grid[idx-1] = p;
	}
	if (g_x == gridSizeInterior){
		grid[idx+1] = p;
	}
	if (g_y == 1){
		grid[idx-gridSize] = p;
	}
	if (g_y == gridSizeInterior){
		grid[idx+gridSize] = p;
	}
}

__kernel void newvel(	__global const float *FGrid,
					__global const float *GGrid,
					__global const float *pGrid,
					__global float *uGrid,
					__global float *vGrid,
					__global float *h_inv,
					__global float *deltaT
				)
{

	float hi = (*h_inv);
	float dt = (*deltaT);

	//TODO DONT HARDCODE
	float velocity = 0.0f;

	int g_x = get_global_id(0) + 1;
	int g_y = get_global_id(1) + 1;

	int gridSize = get_global_size(0) + 2;
	int idx = g_y*gridSize + g_x;

	float F = FGrid[idx];
	float G = GGrid[idx];

	float p = pGrid[idx];
	float p_r = pGrid[idx+1];
	float p_u = pGrid[idx+gridSize];

	uGrid[idx] = F - dt * d_r(p_r, p, hi);//(p_r - p)*hi;
	vGrid[idx] = G - dt * d_r(p_u, p, hi);//(p_u - p)*hi;

	//UPDATE u v
	//left
	if (g_x == 1){
		uGrid[idx-1] = 0;
		vGrid[idx-1] = -1.0f*vGrid[idx];
	}
	//right
	if (g_x == gridSize-2){
		uGrid[idx+1] = 0;
		uGrid[idx] = 0;
		vGrid[idx+1] = -1.0f*vGrid[idx];
	}
	//bottom
	if (g_y == 1){
		uGrid[idx-gridSize] = -1.0f*uGrid[idx];
		vGrid[idx-gridSize] = 0;
	}
	//top
	if (g_y == gridSize-2){
		uGrid[idx+gridSize] = 2.0f*velocity - uGrid[idx];
		vGrid[idx+gridSize] = 0;
		vGrid[idx] = 0;
	}


}

__kernel void rhs(	__global const float *FGrid,
					__global const float *GGrid,
					__global float *rhs,
					__global const float *h_inv,
					__global const float *deltaTInv
				)
{

	float hi = (*h_inv);
	float dti = (*deltaTInv);

	int g_x = get_global_id(0) + 1;
	int g_y = get_global_id(1) + 1;

	int gridSize = get_global_size(0) + 2;
	int idx = g_y*gridSize + g_x;

	float F = FGrid[idx];
	float G = GGrid[idx];

	float F_l = FGrid[idx-1];
	float G_d = GGrid[idx-gridSize];

	rhs[idx] = (d_l(F_l, F, hi) + d_l(G_d, G, hi)) * dti;

}

__kernel void momentumeq(	__global float *FGrid,
					__global float *GGrid,
					__global float *uGrid,
					__global float *vGrid,
					__global const float *h_inv,
					__global const float *h_inv_squared,
					__global const float *deltaT,
					__global const float *RE_inv,
					__global const float *TGrid
				)
{

	float hi = (*h_inv);
	float hsi = (*h_inv_squared);
	float dt = (*deltaT);
	float rei = (*RE_inv);

	//TODO DONT HARDCODE
	float velocity = 0.0f;
	float gravity = -9.81f;
	float beta = 10.0f;//0.5f;

	int g_x = get_global_id(0) + 1;
	int g_y = get_global_id(1) + 1;

	int gridSize = get_global_size(0) + 2;
	int idx = g_y*gridSize + g_x;

	float u = uGrid[idx];
	float v = vGrid[idx];

	float u_r = uGrid[idx+1];
	float v_r = vGrid[idx+1];

	float u_l = uGrid[idx-1];
	float v_l = vGrid[idx-1];

	float u_t = uGrid[idx+gridSize];
	float v_t = vGrid[idx+gridSize];

	float u_d = uGrid[idx-gridSize];
	float v_d = vGrid[idx-gridSize];

	float u_lt = uGrid[idx-1+gridSize];
	float v_dr = vGrid[idx+1-gridSize];


	float A = rei * ( ddx(u_r, u_l, u, hsi) + ddx(u_t, u_d, u, hsi) ) - udu_x(u_r, u_l, u, hi) - vdu_y(u_t, u_d, u, v_r, v_d, v_dr, v, hi);
	float B = rei * ( ddx(v_r, v_l, v, hsi) + ddx(v_t, v_d, v, hsi) ) - vdv_y(v_t, v_d, v, hi) - udv_x(v_r, v_l, v, u_t, u_l, u_lt, u, hi) + gravity;

	FGrid[idx] = u + dt*A;
	GGrid[idx] = v + dt*B - dt * beta * gravity * ( TGrid[idx] + TGrid[idx + gridSize] )*0.5f;


	//UPDATE F G
	//left
	if (g_x == 1){
		FGrid[idx-1] = 0;
		GGrid[idx-1] = -1.0f*GGrid[idx];
	}
	//right
	if (g_x == gridSize-2){
		FGrid[idx+1] = 0;
		FGrid[idx] = 0;
		GGrid[idx+1] = -1.0f*GGrid[idx];
	}
	//bottom
	if (g_y == 1){
		FGrid[idx-gridSize] = -1.0f*FGrid[idx];
		GGrid[idx-gridSize] = 0;
	}
	//top
	if (g_y == gridSize-2){
		FGrid[idx+gridSize] = 2.0f*velocity - FGrid[idx];
		GGrid[idx+gridSize] = 0;
		GGrid[idx] = 0;
	}
}

__kernel void compT(	__global float *TGrid,
					__global float *uGrid,
					__global float *vGrid,
					__global const float *h_inv,
					__global const float *h_inv_squared,
					__global const float *deltaT,
					__global const float *RE_inv,
					__global const float *PR_inv
				)
{



	float hi = (*h_inv);
	float hsi = (*h_inv_squared);
	float dt = (*deltaT);
	float rei = (*RE_inv);
	float pri = (*PR_inv);

	//TODO DONT HARDCODE
	float heat = 1.0f;
	float Q = 0.0f;
	float gamma = 0.9f;

	int g_x = get_global_id(0) + 1;
	int g_y = get_global_id(1) + 1;

	int gridSize = get_global_size(0) + 2;
	int idx = g_y*gridSize + g_x;

	float u = uGrid[idx];
	float v = vGrid[idx];
	float u_l = uGrid[idx-1];
	float v_d = vGrid[idx-gridSize];

	float T = TGrid[idx];
    float T_l = TGrid[idx-1];
    float T_r = TGrid[idx+1];
    float T_u = TGrid[idx+gridSize];
    float T_d = TGrid[idx-gridSize];

    float duTdx = hi * ( (u * (T_r + T)*0.5f - u_l * (T + T_l)*0.5f) + gamma * (fabs(u) * (T - T_r)*0.5 - fabs(u_l) * (T_l - T)*0.5f ) );

 	float dvTdy = hi * ( (v * (T_u + T)*0.5f - v_d * (T + T_d)*0.5f) + gamma * (fabs(v) * (T - T_u)*0.5 - fabs(v_d) * (T_d - T)*0.5f ) );

 	float delT = rei*pri * (ddx(T_r, T_l, T, hsi) + ddx(T_u, T_d, T, hsi)) - duTdx - dvTdy + Q;

 	TGrid[idx] += dt * delT;

 	//UPDATE T
	//left
	if (g_x == 1){
		TGrid[idx-1] = TGrid[idx];
	}
	//right
	if (g_x == gridSize-2){
		TGrid[idx+1] = TGrid[idx];
	}
	//bottom
	if (g_y == 1){
		TGrid[idx-gridSize] = heat;
	}
	//top
	if (g_y == gridSize-2){
		TGrid[idx+gridSize] = TGrid[idx];
	}
}
