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
                  __global float *h_square_inv,
                  //__local float* localGrid,
                  //__local float* localRHS
				  __global float *h,
				  __global float *alpha,
				  __global float *uGrid,
				  __global float *vGrid,
                  )
{
    int g_x = get_global_id(0) + 1;
    int g_y = get_global_id(1) + 1;

    int gridSize = get_global_size(0) + 2;
    int idx 	 = g_y*gridSize + g_x;
	
    
    float p = grid;
    float v = vGrid;
    float u = uGrid;
    
    
    
    int down (int it) {return it-gridSize;}
    int top  (int it) {return it+gridSize;}
    int left (int it) {return it-1;}
    int right(int it) {return it+1;}
    
    
    //BSP  uGrid[idx] = F - dt * (p_r - p)*hi;
	
	
	float dx_l() const {
		//return (_data[it]-_data[it.Left()])/_geom->Mesh()[0];
		return (p[idx]-p[left(idx)])/h[0];
	}
	
	float dx_r() const {
		//return (_data[it.Right()]-_data[it])/_geom->Mesh()[0];
		return (p[right(idx)]-p[idx])/h[0];
	}
	
	float dy_l() const {
		//return (_data[it]-_data[it.Down()])/_geom->Mesh()[1];
		return (p[idx]-p[down(idx)])/h[1];
	}
	
	float dy_r() const {
		//return (_data[it.Top()]-_data[it])/_geom->Mesh()[1];
		return (p[top(idx)]-p[idx])/h[1];
	}
	
	
	
	
	float Grid::dxx() const {
		//return (_data[it.Right()]-2*_data[it]+_data[it.Left()])/_geom->Mesh()[0]/_geom->Mesh()[0];
		return (p[right(idx)]-2*p[idx]+p[left(idx)])/h[0]/h[0];
	}
	
	float dyy() const {
		//return (_data[it.Top()]-2*_data[it]+_data[it.Down()])/_geom->Mesh()[1]/_geom->Mesh()[1];
		return (p[top(idx)]-2*p[idx]+p[down(idx)])/h[1]/h[1];
	}
	
	
	
	
	float DC_udu_x() const {
		//real_t A = (_data[it]+_data[it.Right()])/2.0;
		//real_t B = (_data[it.Left()]+_data[it])/2.0;
		//real_t C = (_data[it]-_data[it.Right()])/2.0;
		//real_t D = (_data[it.Left()]-_data[it])/2.0;
		//return (( A*A - B*B) + alpha*(fabs(A)*C-fabs(B)*D) )/_geom->Mesh()[0];
		float A = (p[idx]+p[right(idx)])/2.0;
		float B = (p[left(idx)]+p[idx])/2.0;
		float C = (p[idx]-p[right(idx)])/2.0;
		float D = (p[left(idx)]-p[idx])/2.0;
		return (( A*A - B*B) + alpha*(fabs(A)*C-fabs(B)*D) )/h[0];
	}
	
	float DC_vdv_y() const {
		//real_t A = (_data[it]+_data[it.Top()])/2.0;
		//real_t B = (_data[it.Down()]+_data[it])/2.0;
		//real_t C = (_data[it]-_data[it.Top()])/2.0;
		//real_t D = (_data[it.Down()]-_data[it])/2.0;
		//return ( (A*A - B*B) + alpha*( fabs(A)*C-fabs(B)*D) )/_geom->Mesh()[1];
		float A = (p[idx]+p[top(idx)])/2.0;
		float B = (p[down(idx)]+p[idx])/2.0;
		float C = (p[idx]-p[top(idx)])/2.0;
		float D = (p[down(idx)]-p[idx])/2.0;
		return (( A*A - B*B) + alpha*(fabs(A)*C-fabs(B)*D) )/h[1];
	}
	
	
	
	float DC_vdu_y() const {
		//real_t A = (v->Cell(it)+v->Cell(it.Right()))/2.0;
		//real_t B = (_data[it]+_data[it.Top()])/2.0;
		//real_t C = (v->Cell(it.Down())+ v->Cell(it.Down().Right()))/2.0;
		//real_t D = (_data[it]+_data[it.Down()] )/2.0;
		//real_t E = (_data[it]-_data[it.Top()])/2.0;
		//real_t F = (_data[it.Down()]-_data[it])/2.0;
		//return (( A*B - C * D) + alpha*( fabs(A)*E-fabs(C)*F ))  /_geom->Mesh()[1];
		float A = (v[idx]+v[right(idx])/2.0;
		float B = (p[idx]+p[top(idx)])/2.0;
		float C = (v[down(idx)]+v[right(down(idx))])/2.0;
		float D = (p[idx]+p[down(idx)])/2.0;
		float E = (p[idx]-p[top(idx)])/2.0;
		float F = (p[down(idx)]-p[idx])/2.0;
		return (( A*B - C * D) + alpha*( fabs(A)*E-fabs(C)*F )) /h[1];
	}

	float DC_udv_x() const {
		//real_t A = (u->Cell(it)+u->Cell(it.Top()))/2.0;
		//real_t B = (_data[it]+_data[it.Right()])/2.0;
		//real_t C = (u->Cell(it.Left())+ u->Cell(it.Left().Top()))/2.0;
		//real_t D = (_data[it]+_data[it.Left()] )/2.0;
		//real_t E = (_data[it]-_data[it.Right()])/2.0;
		//real_t F = (_data[it.Left()]-_data[it] )/2.0;
		//return (( A*B - C * D) + alpha*( fabs(A)*E-fabs(C)*F ))  /_geom->Mesh()[0];
		float A = (u[idx]+u[top(idx])/2.0;
		float B = (p[idx]+p[right(idx)])/2.0;
		float C = (u[left(idx)]+u[left(top(idx))])/2.0;
		float D = (p[idx]+p[left(idx)])/2.0;
		float E = (p[idx]-p[right(idx)])/2.0;
		float F = (p[left(idx)]-p[idx])/2.0;
		return (( A*B - C * D) + alpha*( fabs(A)*E-fabs(C)*F )) /h[0];
	}
	
	
}