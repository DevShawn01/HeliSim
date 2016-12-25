
// graphics functions examples

#include <cmath>   // math functions
#include <cstdio>  // standard I/O functions
#include <cstring> // string manipulation functions

#include <iostream>  // console stream I/O
#include <fstream>   // file stream I/O
#include <strstream> // string stream I/0
#include <iomanip>   // I/O manipulators

#include "3D_graphics.h"
#include "rotation.h"
#include "graphics.h"

// array size for vertex buffer
#define BUFF_MAX 40000

// global variable to store keyboard input
extern unsigned char Key_input[256];

// arrays for storing D3D vertices
static double X[BUFF_MAX],Y[BUFF_MAX],Z[BUFF_MAX];
static double R[BUFF_MAX],G[BUFF_MAX],B[BUFF_MAX];
static double U[BUFF_MAX],V[BUFF_MAX];

const double PI = 4*atan(1.0);

void draw_plane2(double R1[3+1][3+1], double rc[3+1])
// R1 - rotation matrix
// rc - translation vector
{
	int i,n;
	// construct a list of vertices ///////////////////
	i = 0;
	// point #1
	X[i] = -1.0; // x
	Y[i] = -1.0; // y
	Z[i] = 0.01; // z
	R[i] = 1.0; // red component varies from 0 to 1
	G[i] = 0.0; // green component varies from 0 to 1
	B[i] = 1.0; // blue component varies from 0 to 1

	// point #2
	i++;
	X[i] = 2.0; // x
	Y[i] = 0.0; // y
	Z[i] = 0.01; // z
	R[i] = 1.0; // red component varies from 0 to 1
	G[i] = 1.0; // green component varies from 0 to 1
	B[i] = 0.0; // blue component varies from 0 to 1

	// point #3
	i++;
	X[i] = -1.0; // x
	Y[i] = 1.0; // y
	Z[i] = 0.0; // z
	R[i] = 1.0; // red component varies from 0 to 1
	G[i] = 0.0; // green component varies from 0 to 1
	B[i] = 1.0; // blue component varies from 0 to 1

	// point #4
	i++;
	X[i] = -1.0; // x
	Y[i] = -1.0; // y
	Z[i] = 0.0; // z
	R[i] = 0.0; // red component varies from 0 to 1
	G[i] = 0.0; // green component varies from 0 to 1
	B[i] = 1.0; // blue component varies from 0 to 1

	// point #5
	i++;
	X[i] = -1.0; // x
	Y[i] = 1.0; // y
	Z[i] = 0.01; // z
	R[i] = 0.0; // red component varies from 0 to 1
	G[i] = 0.0; // green component varies from 0 to 1
	B[i] = 1.0; // blue component varies from 0 to 1

	// point #6
	i++;
	X[i] = 2.0; // x
	Y[i] = 0.0; // y
	Z[i] = 0.0; // z
	R[i] = 0.0; // red component varies from 0 to 1
	G[i] = 1.0; // green component varies from 0 to 1
	B[i] = 1.0; // blue component varies from 0 to 1

	n = 6; // set number of verticies

	// transform the points to global coord
	// -> rotatate an translate the triangles
	rt_transform(R1,rc,X,Y,Z,X,Y,Z,n);

	draw_triangle_list(X,Y,Z,R,G,B,n);
}


void draw_cube(double R1[3+1][3+1], double rc[3+1], double size,
				  double r, double g, double b)
// R1 - rotation matrix
// rc - translation vector
// size - size
// r,g,b - colour
{
	int n = 36,i;

	// input the 8 vertices for the cube
	const double x1 = 0.0, y1 = 0.0, z1 = 1.0;
	const double x2 = 1.0, y2 = 0.0, z2 = 1.0;
	const double x3 = 0.0, y3 = 1.0, z3 = 1.0;
	const double x4 = 1.0, y4 = 1.0, z4 = 1.0;
	const double x5 = 0.0, y5 = 0.0, z5 = 0.0;
	const double x6 = 1.0, y6 = 0.0, z6 = 0.0;
	const double x7 = 0.0, y7 = 1.0, z7 = 0.0;
	const double x8 = 1.0, y8 = 1.0, z8 = 0.0;

	// construct a vertex list for a cube with size 1
	static double xc[36] = {
	 x1,x2,x4,x1,x4,x3,x2,x6,x8,x2,x8,x4,x3,x4,x8,x3,x8,x7,
	 x1,x6,x2,x1,x5,x6,x5,x8,x6,x5,x7,x8,x1,x7,x5,x1,x3,x7};

	static double yc[36] = { 
	 y1,y2,y4,y1,y4,y3,y2,y6,y8,y2,y8,y4,y3,y4,y8,y3,y8,y7,
	 y1,y6,y2,y1,y5,y6,y5,y8,y6,y5,y7,y8,y1,y7,y5,y1,y3,y7};

	static double zc[36] = {
	 z1,z2,z4,z1,z4,z3,z2,z6,z8,z2,z8,z4,z3,z4,z8,z3,z8,z7,
	 z1,z6,z2,z1,z5,z6,z5,z8,z6,z5,z7,z8,z1,z7,z5,z1,z3,z7};
		
	// process the vertex list
	for(i=0;i<36;i++) {
		// make the colour green
		R[i] = r;
		G[i] = g;
		B[i] = b;
		// scale and translate the vertices
		X[i] = size*(xc[i]-0.5);
		Y[i] = size*(yc[i]-0.5);
		Z[i] = size*(zc[i]-0.5);
	}

	// transform the points to global coord
	// -> rotatate an translate the triangles
	rt_transform(R1,rc,X,Y,Z,X,Y,Z,n);
	draw_triangle_list(X,Y,Z,R,G,B,n);
}

void draw_lines(double R1[3 + 1][3 + 1], double rc[3 + 1]){
	int i, n;

	// construct a list of vertices ///////////////////
	i = 0;
	// point #1
	X[i] = -1.0; // x
	Y[i] = -1.0; // y
	Z[i] = 0.0; // z
	R[i] = 1.0; // red component varies from 0 to 1
	G[i] = 0.0; // green component varies from 0 to 1
	B[i] = 1.0; // blue component varies from 0 to 1

	// point #2
	i++;
	X[i] = 2.0; // x
	Y[i] = 0.0; // y
	Z[i] = 0.0; // z
	R[i] = 1.0; // red component varies from 0 to 1
	G[i] = 1.0; // green component varies from 0 to 1
	B[i] = 0.0; // blue component varies from 0 to 1
	/*
	// point #3
	i++;
	X[i] = -1.0; // x
	Y[i] = 1.0; // y
	Z[i] = 0.0; // z
	R[i] = 1.0; // red component varies from 0 to 1
	G[i] = 0.0; // green component varies from 0 to 1
	B[i] = 1.0; // blue component varies from 0 to 1

	// point #4
	i++;
	X[i] = -1.0; // x
	Y[i] = -1.0; // y
	Z[i] = 0.0; // z
	R[i] = 0.0; // red component varies from 0 to 1
	G[i] = 0.0; // green component varies from 0 to 1
	B[i] = 1.0; // blue component varies from 0 to 1

	// point #5
	i++;
	X[i] = -1.0; // x
	Y[i] = 1.0; // y
	Z[i] = 0.01; // z
	R[i] = 0.0; // red component varies from 0 to 1
	G[i] = 0.0; // green component varies from 0 to 1
	B[i] = 1.0; // blue component varies from 0 to 1

	// point #6
	i++;
	X[i] = 2.0; // x
	Y[i] = 0.0; // y
	Z[i] = 0.0; // z
	R[i] = 0.0; // red component varies from 0 to 1
	G[i] = 1.0; // green component varies from 0 to 1
	B[i] = 1.0; // blue component varies from 0 to 1
	*/
	n = 2; // set number of verticies

	// transform the points to global coord
	// -> rotatate an translate the triangles
	rt_transform(R1, rc, X, Y, Z, X, Y, Z, n);
	//draw_line_list(double *x, double *y, double *z, double *R, double *G, double *B, int n);
	//draw_triangle_list(X, Y, Z, R, G, B, n);
	draw_line_list(X, Y, Z, R, G, B, n);
}

void draw_TriangleStrip(double d, int width, int height, float x0, float z0) {
	static double h[1000][1000];
	static double X[2], Y[2], Z[2];
	static double R[BUFF_MAX], G[BUFF_MAX], B[BUFF_MAX];
	static double R1[] = { 1.0,1.0 }, B1[] = { 1.0,1.0 }, G1[] = { 1.0,1.0 };

	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			// j = row number
			// i = column number

			// horizontal (X,Y,Z) -> (X+1,Y,Z)
			X[0] = d*i;
			X[1] = d*(i + 1);

			Y[0] = h[j][i] * d;
			Y[1] = h[j][i+1] * d;

			Z[0] = d*j;
			Z[1] = d*j;

			// vertical
			X[0] = d*i;
			X[1] = d*i;

			Y[0] = h[j][i] * d;
			Y[1] = h[j+1][i] * d;

			Z[0] = d*j;
			Z[1] = d*(j + 1);

			//diagonal
			X[0] = d*i;
			X[1] = d*(i + 1);

			Y[0] = h[i][j] * d;
			Y[1] = h[i + 1][j + 1] * d;

			Z[0] = d*j;
			Z[1] = d*(j + 1);

		}
	}
}
double linearVal(double y0, double y1, double x, double xdelta) {
	//assumes x0 = 0
	return ((x)*(y1 - y0) / (xdelta)) + y0;
}

double RGBheight(int option, float y, float yo) {
	//option 0 - R
	//option 1 - G
	//option 2 - B

	float pos = (y - yo) / 40.0;
	//gradient( #0008fc 0 % , #207cca 9 % , #f9e52f 18 % , #72d829 26 % , #187f18 44 % , #bc8223 65 % , #7c420c 79 % , #777465 90 % , #ffffff 100 % ); 
	if (pos <= 0.09) {
		if (option == 0) {
			return linearVal(0, 32.0, pos, 0.09) / 255.0;
		}
		else if (option == 1) {
			return linearVal(8.0, 124.0, pos, 0.09) / 255.0;
		}
		else if (option == 2) {
			return linearVal(252.0, 202.0, pos, 0.09) / 255.0;
		}
	}
	else if (pos > 0.09 && pos <= 0.18) {
		pos = pos - 0.09;
		if (option == 0) {
			return linearVal(32.0, 249.0, pos, 0.09) / 255.0;
		}
		else if (option == 1) {
			return linearVal(124.0, 229.0, pos, 0.09) / 255.0;
		}
		else if (option == 2) {
			return linearVal(202.0, 47.0, pos, 0.09) / 255.0;
		}
	}
	else if (pos > 0.18 && pos <= 0.26) {
		pos = pos - 0.18;
		if (option == 0) {
			return linearVal(249.0, 114.0, pos, 0.08) / 255.0;
		}
		else if (option == 1) {
			return linearVal(229.0, 216.0, pos, 0.08) / 255.0;
		}
		else if (option == 2) {
			return linearVal(47.0, 41.0, pos, 0.08) / 255.0;
		}

	}
	else if (pos > 0.26 && pos <= 0.44) {
		pos = pos - 0.26;
		if (option == 0) {
			return linearVal(114.0, 24.0, pos, 0.18) / 255.0;
		}
		else if (option == 1) {
			return linearVal(216.0, 127.0, pos, 0.18) / 255.0;
		}
		else if (option == 2) {
			return linearVal(41.0, 24.0, pos, 0.18) / 255.0;
		}
	}
	else if (pos > 0.44 && pos <= 0.65) {
		pos = pos - 0.44;
		if (option == 0) {
			return linearVal(24.0, 188.0, pos, 0.21) / 255.0;
		}
		else if (option == 1) {
			return linearVal(127.0, 130.0, pos, 0.21) / 255.0;
		}
		else if (option == 2) {
			return linearVal(24.0, 35.0, pos, 0.21) / 255.0;
		}
	}
	else if (pos > 0.65 && pos <= 0.79) {
		pos = pos - 0.65;
		if (option == 0) {
			return linearVal(188.0, 124.0, pos, 0.14) / 255.0;
		}
		else if (option == 1) {
			return linearVal(130.0, 66.0, pos, 0.14) / 255.0;
		}
		else if (option == 2) {
			return linearVal(35.0, 12.0, pos, 0.14) / 255.0;
		}
	}
	else if (pos > 0.79 && pos <= 0.90) {
		pos = pos - 0.79;
		if (option == 0) {
			return linearVal(124.0, 119.0, pos, 0.11) / 255.0;
		}
		else if (option == 1) {
			return linearVal(66.0, 116.0, pos, 0.11) / 255.0;
		}
		else if (option == 2) {
			return linearVal(12.0, 101.0, pos, 0.14) / 255.0;
		}
	}
	else if (pos>0.9 && pos<1.0) {
		pos = pos - 0.9;
		if (option == 0) {
			return linearVal(119.0, 255.0, pos, 0.10) / 255.0;
		}
		else if (option == 1) {
			return linearVal(116.0, 255.0, pos, 0.10) / 255.0;
		}
		else if (option == 2) {
			return linearVal(101.0, 255.0, pos, 0.10) / 255.0;
		}
	}
	else return 1.0;

}
/*

void drawStars(int n, double *aX, double *aY, double *aR, double *aG, double *aB) {
static double X[4], Y[4], Z[4];
static double R[4], G[4], B[4];

for (int i = 0; i < n; i++) {
X[0] = aX[i];
Y[0] = aY[i];
Z[0] = 10.0;
R[0] = aR[i];
G[0] = aG[i];
B[0] = aB[i];

X[1] = X[0] - 2.0;
Y[1] = Y[0];
Z[1] = 10.0;
R[1] = R[0];
G[1] = G[0];
B[1] = B[0];

X[2] = X[0] - 2.0;
Y[2] = Y[0] - 2.0;
Z[2] = 10.0;
R[2] = R[0];
G[2] = G[0];
B[2] = B[0];

draw_triangle_list(X, Y, Z, R, G, B, 3);
}

}

void drawTriangles(double d, int nWidth, int nHeight) {
static double X[3], Y[3], Z[3];
static double R[3], G[3], B[3];
static double h[1000][1000];
float exp = 2.36;

// populate height values
for (int i = 0; i < nHeight + 1; i++) {
for (int j = 0; j < nWidth + 1; j++) {
double nx = (float)j / nWidth - flyingy,
ny = (float)i / nHeight - flying;
double e = noiseOrig(nx, ny) + noiseOrig(2 * nx, 2 * ny) + noiseOrig(4 * nx, 4 * ny);
h[i][j] = pow(e, exp);
}

}

// draw triangles
for (int z = 0; z < nHeight; z++) {
for (int x = 0; x < nWidth; x++) {

X[0] = x*d + xoff;
Z[0] = h[z][x] * d / 2.0 + yoff;
Y[0] = z*d + zoff;
R[0] = RGBheight(0, Y[0], yoff);
G[0] = RGBheight(1, Y[0], yoff);
B[0] = RGBheight(2, Y[0], yoff);

//vertical 0->1
X[1] = x*d + xoff;
Z[1] = h[z + 1][x] * d / 2.0 + yoff;
Y[1] = (z + 1)*d + zoff;
R[1] = RGBheight(0, Y[1],yoff);
G[1] = RGBheight(1, Y[1],yoff);
B[1] = RGBheight(2, Y[1],yoff);

//horizontal 0->2
X[2] = (x + 1)*d + xoff;
Z[2] = h[z][x + 1] * d / 2.0 + yoff;
Y[2] = z*d + zoff;
R[2] = RGBheight(0, Y[2],yoff);
G[2] = RGBheight(1, Y[2],yoff);
B[2] = RGBheight(2, Y[2],yoff);
draw_triangle_list(X, Y, Z, R, G, B, 3);

X[0] = X[1];
Z[0] = Z[1];
Y[0] = Y[1];
R[0] = R[1];
G[0] = G[1];
B[0] = B[1];

X[1] = (x + 1)*d + xoff;
Z[1] = h[z + 1][x + 1] * d / 2.0 + yoff;
Y[1] = (z + 1)*d + zoff;
R[1] = RGBheight(0, Y[0],yoff);
G[1] = RGBheight(1, Y[0],yoff);
B[1] = RGBheight(2, Y[0],yoff);
draw_triangle_list(X, Y, Z, R, G, B, 3);

}
}
}
void draw_TriangleStrip(double d, int width, int height, float x0, float z0) {
	static double verts[BUFF_MAX][3 + 1];
	static double X[2], Y[2], Z[2];
	static double R[BUFF_MAX], G[BUFF_MAX], B[BUFF_MAX];
	static double R1[] = { 1.0,1.0 }, B1[] = { 1.0,1.0 }, G1[] = { 1.0,1.0 };
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			//horizontal 
			X[0] = verts[i + j*width][1] = x0 + d*i;
			Y[0] = verts[i + j*width][2] = terr[i][j];
			Z[0] = verts[i + j*width][3] = z0 + d*j;

			X[1] = verts[i + j*width + 1][1] = x0 + d*(i + 1);
			Y[1] = verts[i + j*width + 1][2] = terr[i + 1][j];
			Z[1] = verts[i + j*width + 1][3] = z0 + d*j;
			draw_line_list(X, Y, Z, R1, G1, B1, 2);
			//vertical
			X[0] = verts[i + j*width][1] = x0 + d*i;
			Y[0] = verts[i + j*width][2] = terr[i][j];;
			Z[0] = verts[i + j*width][3] = z0 + d*j;

			X[1] = verts[i + j*width + 2][1] = x0 + d*i;
			Y[1] = verts[i + j*width + 2][2] = terr[i][j + 1];;
			Z[1] = verts[i + j*width + 2][3] = z0 + d*(j + 1);
			draw_line_list(X, Y, Z, R1, G1, B1, 2);
			//diagonal
			X[0] = verts[i + j*width + 1][1] = x0 + d*i;
			Y[0] = verts[i + j*width + 1][2] = terr[i][j + 1];;
			Z[0] = verts[i + j*width + 1][3] = z0 + d*j;

			X[1] = verts[i + j*width + 1][1] = x0 + d*i;
			Y[1] = verts[i + j*width + 1][2] = terr[i + 1][j];;
			Z[1] = verts[i + j*width + 1][3] = z0 + d*(j + 1);
			draw_line_list(X, Y, Z, R1, G1, B1, 2);
		}
	}
}

*/


