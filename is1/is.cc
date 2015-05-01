#include "is.h"
#include "vector.h"
#include<iostream>
 #include <omp.h>
using namespace std;
/*double4_t calObj(double4_t* inter,int nx,int x0,int y0,int x1,int y1,double4_t VPC,int X, int Y)
{
	double4_t VXC = double4_0;
	//double4_t ret = double4_0;
	//cout<<"Inside"<<endl;
	for(int y=y0; y<=y1; y++)
		for(int x=x0; x<=x1; x++)
		{
			VXC += inter[y*nx+x];
			//cout<<inter[y*nx+x][0]<<" "<<inter[y*nx+x][1]<<" "<<inter[y*nx+x][2]<<endl;
		}
	//cout<<"End"<<endl;
	ret = ((VXC*VXC)/X) ;
	if(Y>0)
		ret += (((VPC - VXC) * (VPC - VXC))/ Y);
	//cout<<ret[0]<<" "<<ret[1]<<" "<<ret[2]<<endl;
	//cout<<"End"<<endl;  
	return VXC;
}  */

Result segment(int ny, int nx, const float* data) {
    // FIXME
	double4_t *inter= double4_alloc(ny*nx);
	double4_t VPC = double4_0;
	double res[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	double4_t *s= double4_alloc((ny+1)*(nx+1));
	#pragma omp parallel for 
	for(int y=0; y<ny; y++)
		for(int x=0; x<nx; x++)
		{
			inter[y*nx+x] =  double4_0;
			for(int c = 0;c<3;c++)
			{
				inter[y*nx+x][c]= data[c + 3 * x + 3 * nx * y];
				//cout<<data[c + 3 * x + 3 * nx * y]<<" ";
			}
			
		}
	//cout<<endl<<"VPC"<<VPC[0]<<" "<<VPC[1]<<" "<<VPC[2]<<endl;
	for(int y =0;y<=ny;y++)
		s[y*(nx+1)] = double4_0;
	for(int x=0; x<=nx; x++)
		 s[x]= double4_0;
	for(int y=1; y<=ny; y++)
		for(int x=1; x<=nx; x++)
		{
				s[y*(nx+1)+(x)] = s[(y-1)*(nx+1)+x] +  s[y*(nx+1)+x-1] - s[(y-1)*(nx+1)+x-1] + inter[(y-1)*nx+(x-1)];
				VPC = VPC + inter[(y-1)*nx+(x-1)];
		}
	int opt_x0_t[10],opt_y0_t[10],opt_x1_t[10],opt_y1_t[10];
	#pragma omp parallel for schedule(static,10)
	for(int hy = 1;hy <= ny;hy++)
		for(int wx = 1;wx <=nx;wx++)
		{
			int t1 = wx*hy;
			int t2 = (nx*ny);
			double X = 1.0/t1;
			double Y = 1.0;
			if(t1< t2)
			{
				Y = 1.0/(t2-t1);
				asm ("#dummy");
			}
			for(int y0 = 0;y0 < ny-hy+1;y0++)
				for(int x0 = 0;x0 < nx-wx+1;x0++)
				{
					int x1 = x0 + wx;
					int y1 = y0 + hy;
					double4_t VXC = double4_0;
					VXC = s[(y1*(nx+1))+(x1)]+ s[(y0*(nx+1))+(x0)] - s[(y0*(nx+1))+(x1)] - s[(y1*(nx+1))+(x0)];
					//cout<<"X=="<<X<<" "<<Y<<endl;	
					//cout<<"VXC=="<<VXC[0]<<" "<<VXC[1]<<" "<<VXC[2]<<endl;
					double4_t temp = ((VXC*VXC)*X) +(((VPC - VXC) * (VPC - VXC))* (Y)); 

					//cout<<"VXC===="<<(VXC[0]*VXC[0])*X<<" "<<((VPC[0] - VXC[0]) * (VPC[0] - VXC[0]) )*Y<<endl;
					double t =temp[0]+temp[1]+temp[2];
					//cout<<endl<<"Test=="<<x0<<" "<<y0<<" "<<x1<<" "<<y1<<" "<<t<<endl;
					int i = omp_get_thread_num(); 
					if( t > res[i] )
					{
						opt_x0_t[i] = x0;
						opt_y0_t[i] = y0;
						opt_x1_t[i] = x1;
						opt_y1_t[i] = y1;	
						res[i]= t;
						asm ("#dummy");
						////cout<<"Yes"<<endl;
					}
				}
		}
	int opt_x0=0,opt_y0=0,opt_x1=0,opt_y1=0;
	double max = 0.0;
	for(int i = 0;i<10;i++)
		if(res[i] > max)
		{	
			opt_x0=opt_x0_t[i] ;
			opt_y0=opt_y0_t[i] ;
			opt_x1=opt_x1_t[i];
			opt_y1=opt_y1_t[i] ;	
			max = res[i];
		}
	double4_t a_star = double4_0;
	double4_t b_star = double4_0;
	for(int y=opt_y0; y<opt_y1; y++)
		for(int x=opt_x0; x<opt_x1; x++)
		{
			a_star += inter[y*nx+x];
			////cout<<"Inter "<<inter[y*nx+x][0]<<" "<<inter[y*nx+x][1]<<" "<<inter[y*nx+x][2]<<endl;
			////cout<<"No"<<endl;
		}
	int X = (opt_x1-opt_x0)*(opt_y1-opt_y0);
	if(((nx*ny) - X) > 0 )	
		b_star = (VPC-a_star)/((nx*ny) - X) ;
	a_star = a_star/X;
	float a[3] = { 0.0,0.0,0.0};
	float b[3] = {0.0, 0.0,0.0};
	for(int i = 0 ;i <3; i++)
		a[i] = a_star[i],b[i]=b_star[i];
	////cout<<X << " " <<((nx*ny) - X)<<endl;
	////cout<<"A "<<a_star[0]<<" "<<a_star[1]<<" "<<a_star[2]<<endl;
	////cout<<"B "<<b_star[0]<<" "<<b_star[1]<<" "<<b_star[2]<<endl;
	Result result { opt_y0,opt_x0,opt_y1,opt_x1, {b[0], b[1], b[2]}, {a[0], a[1],a[2]} };
  //Result result { ny/3, nx/3, 2*ny/3, 2*nx/3, {0.0f, 0.0f, 1.0f}, {1.0f, 0.0f, 0.0f} };
	return result;
}
