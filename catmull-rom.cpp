#include <stdlib.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

// ja desenha o teapot


float camX = 0, camY, camZ = 5;
int startX, startY, tracking = 0;
//float y[3]={0,-1,0};
int alpha = 0, beta = 0, r = 5;
float y[3]={0,0,1};

float res[3];
float res1[3];
float res2[3];
float res3[3];

#define POINT_COUNT 5
// Points that make up the loop for catmull-rom interpolation
float p[POINT_COUNT][3] = {{-1,0,0},{-1,1,0},{1,1,0},{1,0,0},{1,-1,0}};
float p2[POINT_COUNT][3] = {{1,0,0},{1,-1,0},{-1,-1,0},{-1,0,0},{1,-1,0}};
float p3[POINT_COUNT][3] = {{-0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5},{0.5,-0.5,0.5},{1,-1,0}};
float p4[POINT_COUNT][3] = {{-0.25,-0.25,0.75},{-0.25,0.25,0.75},{0.25,0.25,0.75},{0.25,-0.25,0.75},{1,-1,0}};

// vai gerar 
void getCatmullRomPoint(float t, int *indices, float *res,float *deriv, int pn) {

	// catmull-rom matrix
	float m[4][4] = {	{-1.0f,  3.0f, -3.0f,  1.0f},
						{ 3.0f, -6.0f,  3.0f,  0.0f},
						{-3.0f,  3.0f,  0.0f,  0.0f},
						{ 1.0f,  0.0f,  0.0f,  0.0f}};
						
	res[0] = 0.0; res[1] = 0.0; res[2] = 0.0;
	float ax[4];
	float ay[4];
	float az[4];	

	if(pn==0){
	ax[0]=m[0][0]*p[indices[0]][0]+m[0][1]*p[indices[1]][0]+m[0][2]*p[indices[2]][0]+m[0][3]*p[indices[3]][0];
	ax[1]=m[1][0]*p[indices[0]][0]+m[1][1]*p[indices[1]][0]+m[1][2]*p[indices[2]][0]+m[1][3]*p[indices[3]][0];
	ax[2]=m[2][0]*p[indices[0]][0]+m[2][1]*p[indices[1]][0]+m[2][2]*p[indices[2]][0]+m[2][3]*p[indices[3]][0];
	ax[3]=m[3][0]*p[indices[0]][0]+m[3][1]*p[indices[1]][0]+m[3][2]*p[indices[2]][0]+m[3][3]*p[indices[3]][0];

	ay[0]=m[0][0]*p[indices[0]][1]+m[0][1]*p[indices[1]][1]+m[0][2]*p[indices[2]][1]+m[0][3]*p[indices[3]][1];
	ay[1]=m[1][0]*p[indices[0]][1]+m[1][1]*p[indices[1]][1]+m[1][2]*p[indices[2]][1]+m[1][3]*p[indices[3]][1];
	ay[2]=m[2][0]*p[indices[0]][1]+m[2][1]*p[indices[1]][1]+m[2][2]*p[indices[2]][1]+m[2][3]*p[indices[3]][1];
	ay[3]=m[3][0]*p[indices[0]][1]+m[3][1]*p[indices[1]][1]+m[3][2]*p[indices[2]][1]+m[3][3]*p[indices[3]][1];

	az[0]=m[0][0]*p[indices[0]][2]+m[0][1]*p[indices[1]][2]+m[0][2]*p[indices[2]][2]+m[0][3]*p[indices[3]][2];
	az[1]=m[1][0]*p[indices[0]][2]+m[1][1]*p[indices[1]][2]+m[1][2]*p[indices[2]][2]+m[1][3]*p[indices[3]][2];
	az[2]=m[2][0]*p[indices[0]][2]+m[2][1]*p[indices[1]][2]+m[2][2]*p[indices[2]][2]+m[2][3]*p[indices[3]][2];
	az[3]=m[3][0]*p[indices[0]][2]+m[3][1]*p[indices[1]][2]+m[3][2]*p[indices[2]][2]+m[3][3]*p[indices[3]][2];
}

	else if(pn==1){
	ax[0]=m[0][0]*p2[indices[0]][0]+m[0][1]*p2[indices[1]][0]+m[0][2]*p2[indices[2]][0]+m[0][3]*p2[indices[3]][0];
	ax[1]=m[1][0]*p2[indices[0]][0]+m[1][1]*p2[indices[1]][0]+m[1][2]*p2[indices[2]][0]+m[1][3]*p2[indices[3]][0];
	ax[2]=m[2][0]*p2[indices[0]][0]+m[2][1]*p2[indices[1]][0]+m[2][2]*p2[indices[2]][0]+m[2][3]*p2[indices[3]][0];
	ax[3]=m[3][0]*p2[indices[0]][0]+m[3][1]*p2[indices[1]][0]+m[3][2]*p2[indices[2]][0]+m[3][3]*p2[indices[3]][0];

	ay[0]=m[0][0]*p2[indices[0]][1]+m[0][1]*p2[indices[1]][1]+m[0][2]*p2[indices[2]][1]+m[0][3]*p2[indices[3]][1];
	ay[1]=m[1][0]*p2[indices[0]][1]+m[1][1]*p2[indices[1]][1]+m[1][2]*p2[indices[2]][1]+m[1][3]*p2[indices[3]][1];
	ay[2]=m[2][0]*p2[indices[0]][1]+m[2][1]*p2[indices[1]][1]+m[2][2]*p2[indices[2]][1]+m[2][3]*p2[indices[3]][1];
	ay[3]=m[3][0]*p2[indices[0]][1]+m[3][1]*p2[indices[1]][1]+m[3][2]*p2[indices[2]][1]+m[3][3]*p2[indices[3]][1];

	az[0]=m[0][0]*p2[indices[0]][2]+m[0][1]*p2[indices[1]][2]+m[0][2]*p2[indices[2]][2]+m[0][3]*p2[indices[3]][2];
	az[1]=m[1][0]*p2[indices[0]][2]+m[1][1]*p2[indices[1]][2]+m[1][2]*p2[indices[2]][2]+m[1][3]*p2[indices[3]][2];
	az[2]=m[2][0]*p2[indices[0]][2]+m[2][1]*p2[indices[1]][2]+m[2][2]*p2[indices[2]][2]+m[2][3]*p2[indices[3]][2];
	az[3]=m[3][0]*p2[indices[0]][2]+m[3][1]*p2[indices[1]][2]+m[3][2]*p2[indices[2]][2]+m[3][3]*p2[indices[3]][2];
}


	else if(pn==2){
	ax[0]=m[0][0]*p3[indices[0]][0]+m[0][1]*p3[indices[1]][0]+m[0][2]*p3[indices[2]][0]+m[0][3]*p3[indices[3]][0];
	ax[1]=m[1][0]*p3[indices[0]][0]+m[1][1]*p3[indices[1]][0]+m[1][2]*p3[indices[2]][0]+m[1][3]*p3[indices[3]][0];
	ax[2]=m[2][0]*p3[indices[0]][0]+m[2][1]*p3[indices[1]][0]+m[2][2]*p3[indices[2]][0]+m[2][3]*p3[indices[3]][0];
	ax[3]=m[3][0]*p3[indices[0]][0]+m[3][1]*p3[indices[1]][0]+m[3][2]*p3[indices[2]][0]+m[3][3]*p3[indices[3]][0];

	ay[0]=m[0][0]*p3[indices[0]][1]+m[0][1]*p3[indices[1]][1]+m[0][2]*p3[indices[2]][1]+m[0][3]*p3[indices[3]][1];
	ay[1]=m[1][0]*p3[indices[0]][1]+m[1][1]*p3[indices[1]][1]+m[1][2]*p3[indices[2]][1]+m[1][3]*p3[indices[3]][1];
	ay[2]=m[2][0]*p3[indices[0]][1]+m[2][1]*p3[indices[1]][1]+m[2][2]*p3[indices[2]][1]+m[2][3]*p3[indices[3]][1];
	ay[3]=m[3][0]*p3[indices[0]][1]+m[3][1]*p3[indices[1]][1]+m[3][2]*p3[indices[2]][1]+m[3][3]*p3[indices[3]][1];

	az[0]=m[0][0]*p3[indices[0]][2]+m[0][1]*p3[indices[1]][2]+m[0][2]*p3[indices[2]][2]+m[0][3]*p3[indices[3]][2];
	az[1]=m[1][0]*p3[indices[0]][2]+m[1][1]*p3[indices[1]][2]+m[1][2]*p3[indices[2]][2]+m[1][3]*p3[indices[3]][2];
	az[2]=m[2][0]*p3[indices[0]][2]+m[2][1]*p3[indices[1]][2]+m[2][2]*p3[indices[2]][2]+m[2][3]*p3[indices[3]][2];
	az[3]=m[3][0]*p3[indices[0]][2]+m[3][1]*p3[indices[1]][2]+m[3][2]*p3[indices[2]][2]+m[3][3]*p3[indices[3]][2];
}


	else if(pn==3){
	ax[0]=m[0][0]*p4[indices[0]][0]+m[0][1]*p4[indices[1]][0]+m[0][2]*p4[indices[2]][0]+m[0][3]*p4[indices[3]][0];
	ax[1]=m[1][0]*p4[indices[0]][0]+m[1][1]*p4[indices[1]][0]+m[1][2]*p4[indices[2]][0]+m[1][3]*p4[indices[3]][0];
	ax[2]=m[2][0]*p4[indices[0]][0]+m[2][1]*p4[indices[1]][0]+m[2][2]*p4[indices[2]][0]+m[2][3]*p4[indices[3]][0];
	ax[3]=m[3][0]*p4[indices[0]][0]+m[3][1]*p4[indices[1]][0]+m[3][2]*p4[indices[2]][0]+m[3][3]*p4[indices[3]][0];

	ay[0]=m[0][0]*p4[indices[0]][1]+m[0][1]*p4[indices[1]][1]+m[0][2]*p4[indices[2]][1]+m[0][3]*p4[indices[3]][1];
	ay[1]=m[1][0]*p4[indices[0]][1]+m[1][1]*p4[indices[1]][1]+m[1][2]*p4[indices[2]][1]+m[1][3]*p4[indices[3]][1];
	ay[2]=m[2][0]*p4[indices[0]][1]+m[2][1]*p4[indices[1]][1]+m[2][2]*p4[indices[2]][1]+m[2][3]*p4[indices[3]][1];
	ay[3]=m[3][0]*p4[indices[0]][1]+m[3][1]*p4[indices[1]][1]+m[3][2]*p4[indices[2]][1]+m[3][3]*p4[indices[3]][1];

	az[0]=m[0][0]*p4[indices[0]][2]+m[0][1]*p4[indices[1]][2]+m[0][2]*p4[indices[2]][2]+m[0][3]*p4[indices[3]][2];
	az[1]=m[1][0]*p4[indices[0]][2]+m[1][1]*p4[indices[1]][2]+m[1][2]*p4[indices[2]][2]+m[1][3]*p4[indices[3]][2];
	az[2]=m[2][0]*p4[indices[0]][2]+m[2][1]*p4[indices[1]][2]+m[2][2]*p4[indices[2]][2]+m[2][3]*p4[indices[3]][2];
	az[3]=m[3][0]*p4[indices[0]][2]+m[3][1]*p4[indices[1]][2]+m[3][2]*p4[indices[2]][2]+m[3][3]*p4[indices[3]][2];
	}




	// resultado final do produto de matrizes do pdf
	res[0]=pow(t,3)*ax[0]+pow(t,2)*ax[1]+t*ax[2]+ax[3];
	res[1]=pow(t,3)*ay[0]+pow(t,2)*ay[1]+t*ay[2]+ay[3];
	res[2]=pow(t,3)*az[0]+pow(t,2)*az[1]+t*az[2]+az[3];

	// Compute point res = T * M * P
	// where Pi = p[indices[i]]
	// ...
}


void cross(float *a, float *b, float *res) { 
	res[0] = a[1]*b[2] - a[2]*b[1]; 
	res[1] = a[2]*b[0] - a[0]*b[2]; 
	res[2] = a[0]*b[1] - a[1]*b[0];
}



// given  global t, returns the point in the curve
// vai dar um ponto da curva a partir do getCatmullPoint
// usa os 5 pontos acima para gerar os outros pontos
void getGlobalCatmullRomPoint(float gt, float *res, float *deriv, int pn) {

	float t = gt * POINT_COUNT; // this is the real global t
	int index = floor(t);  // which segment
	t = t - index; // where within  the segment

	// indices store the points
	int indices[4]; 
	indices[0] = 0;	
	indices[1] = 1;
	indices[2] = 2; 
	indices[3] = 3;

	getCatmullRomPoint(t, indices, res, deriv, pn);
}


void changeSize(int w, int h) {

	// Prevent a divide by zero, when window is too short
	// (you cant make a window with zero width).
	if(h == 0)
		h = 1;

	// compute window's aspect ratio 
	float ratio = w * 1.0 / h;

	// Reset the coordinate system before modifying
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
	// Set the viewport to be the entire window
    glViewport(0, 0, w, h);

	// Set the correct perspective
	gluPerspective(45,ratio,1,1000);

	// return to the model view matrix mode
	glMatrixMode(GL_MODELVIEW);
}


void renderCatmullRomCurve() {

// desenhar a curva usando segmentos de reta - GL_LINE_LOOP
	float res[3];
	float deriv[3];
	
	//gera a linha a partir dos pontos gerados acima
	for(int j=0;j<2;j++){
	glBegin(GL_LINE_LOOP);
	for(int i=0; i<100; i++){
		getGlobalCatmullRomPoint(i/100.0,res,deriv,j);
		glVertex3f(res[0],res[1],res[2]);
	}

	glEnd();
	}
	

}

void buildRotMatrix(float *x, float *y, float *z, float *m) {
	m[0] = x[0];
	m[1] = x[1]; 
	m[2] = x[2]; 
	m[3] = 0; 
	m[4] = y[0]; 
	m[5] = y[1];
	m[6] = y[2];
	m[7] = 0;
	m[8] = z[0]; 
	m[9] = z[1];
	m[10] = z[2];
	m[11] = 0; 
	m[12] = 0;
	m[13] = 0; 
	m[14] = 0; 
	m[15] = 1;
}

void normalize(float *a) { 
	float l = sqrt(a[0]*a[0] + a[1] * a[1] + a[2] * a[2]); 
	a[0] = a[0]/l; 
	a[1] = a[1]/l; 
	a[2] = a[2]/l;
}

float* getPoints(int *patchs, int n_patch, float *vertices, int n_vertices, int detail, float t){

	float Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,Ex,Ey,Ez,Fx,Fy,Fz,Gx,Gy,Gz,Hx,Hy,Hz,Ix,Iy,Iz,Jx,Jy,Jz,Kx,Ky,Kz,Lx,Ly,Lz,Mx,My,Mz,Nx,Ny,Nz,Ox,Oy,Oz,Px,Py,Pz;
	float change = 1.0 / detail, *points=(float*)malloc(n_patch*(3*(detail+1)*(detail+1))*sizeof(float));
	float a,b,c,d;
	int v=0,i,j,k;

	float m[4][4] = {	{-1.0f,  3.0f, -3.0f,  1.0f},
						{ 3.0f, -6.0f,  3.0f,  0.0f},
						{-3.0f,  3.0f,  0.0f,  0.0f},
						{ 1.0f,  0.0f,  0.0f,  0.0f}};

	float ax[4];
	float ay[4];
	float az[4];

    for(k=n_patch;k<n_patch+1;k++){
	//for(k=n_patch;k<n_patch+1;k++){

        Ax = vertices[3*patchs[(16*k)]];         Ay = vertices[3*patchs[(16*k)]+1];     Az = vertices[3*patchs[(16*k)]+2];
        Bx = vertices[3*patchs[(16*k+1)]];       By = vertices[3*patchs[(16*k+1)]+1];   Bz = vertices[3*patchs[(16*k)+1]+2];
        Cx = vertices[3*patchs[(16*k+2)]];       Cy = vertices[3*patchs[(16*k+2)]+1];   Cz = vertices[3*patchs[(16*k)+2]+2];
        Dx = vertices[3*patchs[(16*k+3)]];       Dy = vertices[3*patchs[(16*k+3)]+1];   Dz = vertices[3*patchs[(16*k)+3]+2];
        Ex = vertices[3*patchs[(16*k+4)]];       Ey = vertices[3*patchs[(16*k+4)]+1];   Ez = vertices[3*patchs[(16*k)+4]+2];
        Fx = vertices[3*patchs[(16*k+5)]];       Fy = vertices[3*patchs[(16*k+5)]+1];   Fz = vertices[3*patchs[(16*k)+5]+2];
        Gx = vertices[3*patchs[(16*k+6)]];       Gy = vertices[3*patchs[(16*k+6)]+1];   Gz = vertices[3*patchs[(16*k)+6]+2];
        Hx = vertices[3*patchs[(16*k+7)]];       Hy = vertices[3*patchs[(16*k+7)]+1];   Hz = vertices[3*patchs[(16*k)+7]+2];
        Ix = vertices[3*patchs[(16*k+8)]];       Iy = vertices[3*patchs[(16*k+8)]+1];   Iz = vertices[3*patchs[(16*k)+8]+2];
        Jx = vertices[3*patchs[(16*k+9)]];       Jy = vertices[3*patchs[(16*k+9)]+1];   Jz = vertices[3*patchs[(16*k)+9]+2];
        Kx = vertices[3*patchs[(16*k+10)]];      Ky = vertices[3*patchs[(16*k+10)]+1];  Kz = vertices[3*patchs[(16*k)+10]+2];
        Lx = vertices[3*patchs[(16*k+11)]];      Ly = vertices[3*patchs[(16*k+11)]+1];  Lz = vertices[3*patchs[(16*k)+11]+2];
        Mx = vertices[3*patchs[(16*k+12)]];      My = vertices[3*patchs[(16*k+12)]+1];  Mz = vertices[3*patchs[(16*k)+12]+2];
        Nx = vertices[3*patchs[(16*k+13)]];      Ny = vertices[3*patchs[(16*k+13)]+1];  Nz = vertices[3*patchs[(16*k)+13]+2];
        Ox = vertices[3*patchs[(16*k+14)]];      Oy = vertices[3*patchs[(16*k+14)]+1];  Oz = vertices[3*patchs[(16*k)+14]+2];
        Px = vertices[3*patchs[(16*k+15)]];      Py = vertices[3*patchs[(16*k+15)]+1];  Pz = vertices[3*patchs[(16*k)+15]+2];

        /*for(i=0; i<=detail; i++){
            for(j=0; j<=detail; j++){
                
                // First get the vertices
                points[v++] = Ax*a*a*a*c*c*c   + Bx*3*a*a*a*c*c*d
                            + Cx*3*a*a*a*c*d*d + Dx*a*a*a*d*d*d
                            + Ex*3*a*a*b*c*c*c + Fx*9*a*a*b*c*c*d
                            + Gx*9*a*a*b*c*d*d + Hx*3*a*a*b*d*d*d
                            + Ix*3*a*b*b*c*c*c + Jx*9*a*b*b*c*c*d
                            + Kx*9*a*b*b*c*d*d + Lx*3*a*b*b*d*d*d
                            + Mx*b*b*b*c*c*c   + Nx*3*b*b*b*c*c*d
                            + Ox*3*b*b*b*c*d*d + Px*b*b*b*d*d*d;
                                      
                points[v++] = Ay*a*a*a*c*c*c   + By*3*a*a*a*c*c*d
                            + Cy*3*a*a*a*c*d*d + Dy*a*a*a*d*d*d
                            + Ey*3*a*a*b*c*c*c + Fy*9*a*a*b*c*c*d
                            + Gy*9*a*a*b*c*d*d + Hy*3*a*a*b*d*d*d
                            + Iy*3*a*b*b*c*c*c + Jy*9*a*b*b*c*c*d
                            + Ky*9*a*b*b*c*d*d + Ly*3*a*b*b*d*d*d
                            + My*b*b*b*c*c*c   + Ny*3*b*b*b*c*c*d
                            + Oy*3*b*b*b*c*d*d + Py*b*b*b*d*d*d;
            
                points[v++] = Az*a*a*a*c*c*c   + Bz*3*a*a*a*c*c*d
                            + Cz*3*a*a*a*c*d*d + Dz*a*a*a*d*d*d
                            + Ez*3*a*a*b*c*c*c + Fz*9*a*a*b*c*c*d
                            + Gz*9*a*a*b*c*d*d + Hz*3*a*a*b*d*d*d
                            + Iz*3*a*b*b*c*c*c + Jz*9*a*b*b*c*c*d
                            + Kz*9*a*b*b*c*d*d + Lz*3*a*b*b*d*d*d
                            + Mz*b*b*b*c*c*c   + Nz*3*b*b*b*c*c*d
                            + Oz*3*b*b*b*c*d*d + Pz*b*b*b*d*d*d;
                
                
                //change the c-variable within the inner loop
                c += change;
                d  = 1.0 - c;
            }
            //change the a-variable outside the inner loop
            a += change;
            b  = 1.0 - a;
                                  
            // Reset the c-variable to make it ready for the inner loop again
            c = 0.0;
            d = 1.0 - c;
        }
    }

	return points;*/


        ax[0]=m[0][0]*Ax+m[0][1]*Bx+m[0][2]*Cx+m[0][3]*Dx;
		ax[1]=m[1][0]*Ax+m[1][1]*Bx+m[1][2]*Cx+m[1][3]*Dx;
		ax[2]=m[2][0]*Ax+m[2][1]*Bx+m[2][2]*Cx+m[2][3]*Dx;
		ax[3]=m[3][0]*Ax+m[3][1]*Bx+m[3][2]*Cx+m[3][3]*Dx;
		
		ay[0]=m[0][0]*Ay+m[0][1]*By+m[0][2]*Cy+m[0][3]*Dy;
		ay[1]=m[1][0]*Ay+m[1][1]*By+m[1][2]*Cy+m[1][3]*Dy;
		ay[2]=m[2][0]*Ay+m[2][1]*By+m[2][2]*Cy+m[2][3]*Dy;
		ay[3]=m[3][0]*Ay+m[3][1]*By+m[3][2]*Cy+m[3][3]*Dy;

		az[0]=m[0][0]*Az+m[0][1]*Bz+m[0][2]*Cz+m[0][3]*Dz;
		az[1]=m[1][0]*Az+m[1][1]*Bz+m[1][2]*Cz+m[1][3]*Dz;
		az[2]=m[2][0]*Az+m[2][1]*Bz+m[2][2]*Cz+m[2][3]*Dz;
		az[3]=m[3][0]*Az+m[3][1]*Bz+m[3][2]*Cz+m[3][3]*Dz;

		res[0]=pow(t,3)*ax[0]+pow(t,2)*ax[1]+t*ax[2]+ax[3];
		res[1]=pow(t,3)*ay[0]+pow(t,2)*ay[1]+t*ay[2]+ay[3];
		res[2]=pow(t,3)*az[0]+pow(t,2)*az[1]+t*az[2]+az[3];

		ax[0]=m[0][0]*Ex+m[0][1]*Fx+m[0][2]*Gx+m[0][3]*Hx;
		ax[1]=m[1][0]*Ex+m[1][1]*Fx+m[1][2]*Gx+m[1][3]*Hx;
		ax[2]=m[2][0]*Ex+m[2][1]*Fx+m[2][2]*Gx+m[2][3]*Hx;
		ax[3]=m[3][0]*Ex+m[3][1]*Fx+m[3][2]*Gx+m[3][3]*Hx;
		
		ay[0]=m[0][0]*Ey+m[0][1]*Fy+m[0][2]*Gy+m[0][3]*Hy;
		ay[1]=m[1][0]*Ey+m[1][1]*Fy+m[1][2]*Gy+m[1][3]*Hy;
		ay[2]=m[2][0]*Ey+m[2][1]*Fy+m[2][2]*Gy+m[2][3]*Hy;
		ay[3]=m[3][0]*Ey+m[3][1]*Fy+m[3][2]*Gy+m[3][3]*Hy;

		az[0]=m[0][0]*Ez+m[0][1]*Fz+m[0][2]*Gz+m[0][3]*Hz;
		az[1]=m[1][0]*Ez+m[1][1]*Fz+m[1][2]*Gz+m[1][3]*Hz;
		az[2]=m[2][0]*Ez+m[2][1]*Fz+m[2][2]*Gz+m[2][3]*Hz;
		az[3]=m[3][0]*Ez+m[3][1]*Fz+m[3][2]*Gz+m[3][3]*Hz;

		res1[0]=pow(t,3)*ax[0]+pow(t,2)*ax[1]+t*ax[2]+ax[3];
		res1[1]=pow(t,3)*ay[0]+pow(t,2)*ay[1]+t*ay[2]+ay[3];
		res1[2]=pow(t,3)*az[0]+pow(t,2)*az[1]+t*az[2]+az[3];

		ax[0]=m[0][0]*Ix+m[0][1]*Jx+m[0][2]*Kx+m[0][3]*Lx;
		ax[1]=m[1][0]*Ix+m[1][1]*Jx+m[1][2]*Kx+m[1][3]*Lx;
		ax[2]=m[2][0]*Ix+m[2][1]*Jx+m[2][2]*Kx+m[2][3]*Lx;
		ax[3]=m[3][0]*Ix+m[3][1]*Jx+m[3][2]*Kx+m[3][3]*Lx;
		
		ay[0]=m[0][0]*Iy+m[0][1]*Jy+m[0][2]*Ky+m[0][3]*Ly;
		ay[1]=m[1][0]*Iy+m[1][1]*Jy+m[1][2]*Ky+m[1][3]*Ly;
		ay[2]=m[2][0]*Iy+m[2][1]*Jy+m[2][2]*Ky+m[2][3]*Ly;
		ay[3]=m[3][0]*Iy+m[3][1]*Jy+m[3][2]*Ky+m[3][3]*Ly;

		az[0]=m[0][0]*Iz+m[0][1]*Jz+m[0][2]*Kz+m[0][3]*Lz;
		az[1]=m[1][0]*Iz+m[1][1]*Jz+m[1][2]*Kz+m[1][3]*Lz;
		az[2]=m[2][0]*Iz+m[2][1]*Jz+m[2][2]*Kz+m[2][3]*Lz;
		az[3]=m[3][0]*Iz+m[3][1]*Jz+m[3][2]*Kz+m[3][3]*Lz;

		res2[0]=pow(t,3)*ax[0]+pow(t,2)*ax[1]+t*ax[2]+ax[3];
		res2[1]=pow(t,3)*ay[0]+pow(t,2)*ay[1]+t*ay[2]+ay[3];
		res2[2]=pow(t,3)*az[0]+pow(t,2)*az[1]+t*az[2]+az[3];

		ax[0]=m[0][0]*Mx+m[0][1]*Nx+m[0][2]*Ox+m[0][3]*Px;
		ax[1]=m[1][0]*Mx+m[1][1]*Nx+m[1][2]*Ox+m[1][3]*Px;
		ax[2]=m[2][0]*Mx+m[2][1]*Nx+m[2][2]*Ox+m[2][3]*Px;
		ax[3]=m[3][0]*Mx+m[3][1]*Nx+m[3][2]*Ox+m[3][3]*Px;
		
		ay[0]=m[0][0]*My+m[0][1]*Ny+m[0][2]*Oy+m[0][3]*Py;
		ay[1]=m[1][0]*My+m[1][1]*Ny+m[1][2]*Oy+m[1][3]*Py;
		ay[2]=m[2][0]*My+m[2][1]*Ny+m[2][2]*Oy+m[2][3]*Py;
		ay[3]=m[3][0]*My+m[3][1]*Ny+m[3][2]*Oy+m[3][3]*Py;

		az[0]=m[0][0]*Mz+m[0][1]*Nz+m[0][2]*Oz+m[0][3]*Pz;
		az[1]=m[1][0]*Mz+m[1][1]*Nz+m[1][2]*Oz+m[1][3]*Pz;
		az[2]=m[2][0]*Mz+m[2][1]*Nz+m[2][2]*Oz+m[2][3]*Pz;
		az[3]=m[3][0]*Mz+m[3][1]*Nz+m[3][2]*Oz+m[3][3]*Pz;

		res3[0]=pow(t,3)*ax[0]+pow(t,2)*ax[1]+t*ax[2]+ax[3];
		res3[1]=pow(t,3)*ay[0]+pow(t,2)*ay[1]+t*ay[2]+ay[3];
		res3[2]=pow(t,3)*az[0]+pow(t,2)*az[1]+t*az[2]+az[3];
	}

}


void read_Patch(FILE *f_patch, FILE *f, int detail){
    
	int i=0,j=0,v=0,avanco,k,n_patch,n_vertices, *patchs=NULL;
	float *vertices=NULL,*points=NULL,x,y,z;
    
    
	fscanf(f_patch,"%d\n",&n_patch);
	patchs=(int*)malloc(16*n_patch*sizeof(int));
    
	while(i<n_patch){
		for(j=0;j<15;j++){
			fscanf(f_patch,"%d, ",&k);
			patchs[v++]=k;
		}
		fscanf(f_patch,"%d\n",&k);
		patchs[v++]=k;
		i++;
	}
    
	fscanf(f_patch,"%d\n",&n_vertices);
	vertices=(float*)malloc(3*n_vertices*sizeof(float));
	v=0;
	while(fscanf(f_patch," %f, %f, %f\n",&x,&y,&z)!=EOF){
		vertices[v++]=x;
		vertices[v++]=y;
		vertices[v++]=z;
	}

	/*

	points = getPoints(patchs,n_patch,vertices,n_vertices,detail,i/100.0);

	glBegin(GL_LINE_LOOP);
	for(int y=0;y<n_patch*(3*(detail+1)*(detail+1));y+=3){
		glVertex3f(points[y],points[y+1],points[y+2]);
	}
	glEnd();*/

	float a;

	for(int co = n_patch-1;co<n_patch;co++){
	for(int i=0; i<20; i++){
	glBegin(GL_LINE_STRIP);
		getPoints(patchs,co,vertices,n_vertices,detail,i/20.0);
		glVertex3f(res[0],res[1],res[2]);
		glVertex3f(res1[0],res1[1],res1[2]);
		glVertex3f(res2[0],res2[1],res2[2]);
		glVertex3f(res3[0],res3[1],res3[2]);
	glEnd();
	}

	for(int j=0;j<4;j++){
	glBegin(GL_LINE_STRIP);
		for(int k=0;k<5;k++){
			if(j==0){
				glVertex3f(res[0]-(res[0]*k/5), )
			}
		}
	}


	//}

/*	glBegin(GL_LINE_STRIP);
		for(int j=0;j<20;j++){
			glVertex3f(res[0]*)
		glVertex3f(res[0],res[1],res[2]);
		glVertex3f(res1[0],res1[1],res1[2]);
		glVertex3f(res2[0],res2[1],res2[2]);
		glVertex3f(res3[0],res3[1],res3[2]);*/



//}

	/*for(int i=0; i<20; i++){
	glBegin(GL_LINE_STRIP);
		glVertex3f(res[0],res[1],res[2]);
	glEnd();
}*/
/*
	glBegin(GL_LINE_STRIP);
		glVertex3f(res1[0],res1[1],res1[2]);
	glEnd();

	glBegin(GL_LINE_STRIP);
		glVertex3f(res2[0],res2[1],res2[2]);
	glEnd();

	glBegin(GL_LINE_STRIP);
		glVertex3f(res3[0],res3[1],res3[2]);
	glEnd();*/
//}
}


//glTranslatef(6.0,0.0,0.0);
//glutSolidTeapot(2);



	/*glTranslatef(6.0,0.0,0.0);
	glutSolidTeapot(2.0);*/

	//points=getPoints(patchs,n_patch,vertices,n_vertices,detail);
    
	//Imprimir Vertices
	/*n_vertices=n_patch*(3*(detail+1)*(detail+1));
	fprintf(f, "%d\n",n_vertices);
    for(i=0;i<n_vertices;i+=3)
        fprintf(f, "%f %f %f\n",points[i],points[i+1],points[i+2]);
    
    
    //Imprimir Indices
    n_vertices=n_patch * detail * detail *3*2;
    fprintf(f, "%d\n",n_vertices);
    
    avanco=(detail+1)*(detail+1);
    
    for(i=0;i<n_patch;i++){
        for (j=0; j<detail; j++) {
            for (v=0; v<detail; v++) {
                
                
                fprintf(f, "%d %d %d\n",i*avanco + j*(detail+1) + v, i*avanco + j*(detail+1) + v+1 ,i*avanco + (j+1)*(detail+1) + v );
                fprintf(f, "%d %d %d\n", i*avanco + j*(detail+1) + v+1 , i*avanco + (j+1)*(detail+1) + v+1 ,i*avanco + (j+1)*(detail+1) + v  );
              
            }
        }
        
    }*/
}


void renderScene(void) {

	FILE *f=NULL, *patch=NULL;

	patch=fopen("teapot.patch","r");
	f=fopen("teap.txt","w");

	static float t = 0;
	float res[3];
	float deriv[3];
	
	//yo =(1,0,0)
	glClearColor(0.0f,0.0f,0.0f,0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();
	gluLookAt(camX, camY, camZ, 
		      0.0,0.0,0.0,
			  0.0f,1.0f,0.0f);

	// apply transformations here
	// ...

	float m[16]; 
	
	float z[3];

	read_Patch(patch,f,1);

	fclose(patch);
	fclose(f);

	/*getGlobalCatmullRomPoint(t,res,deriv,0);
	getGlobalCatmullRomPoint(t,res,deriv,1);
	getGlobalCatmullRomPoint(t,res,deriv,2);
	getGlobalCatmullRomPoint(t,res,deriv,3);*/
	//normalize(deriv);
	//cross(deriv,y,z);
	//normalize(z);
	//cross(z,deriv,y);
	//glTranslatef(res[0],res[1],res[2]);


	//buildRotMatrix(deriv,y,z,m);
	
	//glMultMatrixf(m); 
	
	//glutWireTeapot(0.1);


	glutSwapBuffers();
	t+=0.001;
}


void processMouseButtons(int button, int state, int xx, int yy) 
{
	if (state == GLUT_DOWN)  {
		startX = xx;
		startY = yy;
		if (button == GLUT_LEFT_BUTTON)
			tracking = 1;
		else if (button == GLUT_RIGHT_BUTTON)
			tracking = 2;
		else
			tracking = 0;
	}
	else if (state == GLUT_UP) {
		if (tracking == 1) {
			alpha += (xx - startX);
			beta += (yy - startY);
		}
		else if (tracking == 2) {
			
			r -= yy - startY;
			if (r < 3)
				r = 3.0;
		}
		tracking = 0;
	}
}

void processKeys(unsigned char c, int xx, int yy) {
  switch (c){
  case 'q': r=r+1.0; break;
  case 'e': r=r-1.0; break;
  case 'w': if((beta+0.1<=M_PI/2) && (beta+0.1>=-M_PI/2)) beta+=0.1; break;
  case 's': if((beta-0.1<=M_PI/2) && (beta-0.1>=-M_PI/2)) beta-=0.1; break;
  case 'm': alpha+=0.1; break;
  case 'n': alpha-=0.1; break;
  default: break;
  }
  glutPostRedisplay();
}


void processMouseMotion(int xx, int yy)
{
	int deltaX, deltaY;
	int alphaAux, betaAux;
	int rAux;

	if (!tracking)
		return;

	deltaX = xx - startX;
	deltaY = yy - startY;

	if (tracking == 1) {

		alphaAux = alpha + deltaX;
		betaAux = beta + deltaY;

		if (betaAux > 85.0)
			betaAux = 85.0;
		else if (betaAux < -85.0)
			betaAux = -85.0;

		rAux = r;
	}
	else if (tracking == 2) {

		alphaAux = alpha;
		betaAux = beta;
		rAux = r - deltaY;
		if (rAux < 3)
			rAux = 3;
	}
	camX = rAux * sin(alphaAux * 3.14 / 180.0) * cos(betaAux * 3.14 / 180.0);
	camZ = rAux * cos(alphaAux * 3.14 / 180.0) * cos(betaAux * 3.14 / 180.0);
	camY = rAux *							     sin(betaAux * 3.14 / 180.0);
}


int main(int argc, char **argv) {

// inicialization
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGBA);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(320,320);
	glutCreateWindow("CG@DI-UM");
		
// callback registration 
	glutDisplayFunc(renderScene);
	glutIdleFunc(renderScene);
	glutReshapeFunc(changeSize);

// mouse callbacks
	glutMouseFunc(processMouseButtons);
	glutMotionFunc(processMouseMotion);
    glutKeyboardFunc(processKeys);

// OpenGL settings
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

// enter GLUT's main cycle 
	glutMainLoop();
}

