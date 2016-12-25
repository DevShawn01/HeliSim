
#include <cmath>   // math functions
#include <cstdio>  // standard I/O functions
#include <cstring> // string manipulation functions
#include <conio.h>

#include <iostream>  // console stream I/O
#include <fstream>   // file stream I/O
#include <strstream> // string stream I/0
#include <iomanip>   // I/O manipulators
#include <windows.h> // for keyboard input

// user defined functions
#include "timer.h" // for measuring time
#include "rotation.h" // for computing rotation functions
#include "3D_graphics.h" // for DirectX 3D graphics

#include "ran.h"
#include "graphics.h"
#include "Heli.h"
#include "noise\noise.h"
noise::module::Perlin pModule;

#define BUFF_MAX 40000
#define KEY(c) ( GetAsyncKeyState((int)(c)) & (SHORT)0x8000 )
RECT desktop;
BOOL ifDp = GetWindowRect(GetDesktopWindow(), &desktop);

// 3D graphics window size
int WIDTH_MIN = 0.0;
int HEIGHT_MIN = 0.0;
int WIDTH_MAX = 1000; // increase this to increase the window width
int HEIGHT_MAX = 800; // increase this to increase the window height

// background colour for the scene
float BACK_R = 0.0f; // red colour component (0 to 1)
float BACK_G = 0.0f; // green colour component (0 to 1)
float BACK_B = 0.0f; // blue colour component (0 to 1)
double flying = 0.0, flyingy = 0.0;
// default min and max viewing distances.
// objects closer than VMIN and farther than VMAX are not drawn (ie cannot be seen).
// note the ratio of VMAX/VMIN should be less than 10,000 for most graphics cards.
double VMIN = 1.0; // units of m (or whatever units you draw your object in)
double VMAX = 1000.0; // units of m                                        190

const double PI = 4*atan(1.0);

// global variable for keyboard input
extern unsigned char Key_input[256];

using namespace std;

int mheight = 0, mwidth = 0;
double xoff = 0.0, yoff=0.0, zoff = 0.0;
double scl;


void helpMenu() {
	text_xy("Increase Thrust: Space", 10, 10, 20);
	text_xy("Decrease Thrust: L Control", 10, 50, 20);
	text_xy("Roll CCW: A", 10, 100, 20);
	text_xy("Roll CW: D", 10, 150, 20);
	text_xy("Pitch Down: W", 10, 200, 20);
	text_xy("Pitch Up: S", 10, 250, 20);
	text_xy("Yaw CCW: L Arrow", 10, 300, 20);
	text_xy("Yaw CW: R Arrow", 10, 350, 20);
}

double noiseOrig(double nx, double ny) {
	return pModule.GetValue(nx, ny, 0) / 2.0 + 0.5;
}

double drawTerrain(double d, int nWidth, int nHeight){
	static int n=2;
	static double X[2], Y[2], Z[2];
	static double R[2], G[2], B[2];
	static double h[101][101];
	float exp = 5.36; // increase this for taller/sharper mountaints

		// populate height values
		for (int i = 0; i < nHeight+1; i++) {
			for (int j = 0; j < nWidth+1; j++) {
				// Height using Perlin noise
				double nx = (float)j / nWidth - flyingy,
					ny = (float)i / nHeight - flying;
				double e = noiseOrig(nx, ny) + 0.4*noiseOrig(2 * nx, 2 * ny) + 0.6*noiseOrig(4 * nx, 4 * ny); // change 0.4 amd 0.6 for different terrain variation
				h[i][j] = pow(e, exp);
			}
		}

		// draw triangles
		for (int z = 0; z < nHeight; z++) {
			for (int x = 0; x < nWidth; x++) {

				X[0] = x*d + xoff;
				Z[0] = h[z][x] * d / 2.0 + yoff;
				Y[0] = z*d + zoff;
				R[0] = RGBheight(0, Z[0],yoff);
				G[0] = RGBheight(1, Z[0],yoff);
				B[0] = RGBheight(2, Z[0],yoff);

				//vertical 0->1
				X[1] = x*d + xoff;
				Z[1] = h[z+1][x]*d/2.0 + yoff;
				Y[1] = (z + 1)*d + zoff;
				R[1] = RGBheight(0, Z[1],yoff);
				G[1] = RGBheight(1, Z[1],yoff);
				B[1] = RGBheight(2, Z[1],yoff);
				draw_line_list(X, Y, Z, R, G, B, n);

				//horizontal 0->2
				X[1] = (x + 1)*d + xoff;
				Z[1] = h[z][x+1]*d/2.0 + yoff;
				Y[1] = z*d + zoff;
				R[1] = RGBheight(0, Z[1], yoff);
				G[1] = RGBheight(1, Z[1], yoff);
				B[1] = RGBheight(2, Z[1], yoff);
				draw_line_list(X, Y, Z, R, G, B, n);

				//diagonal
				X[1] = (x + 1)*d + xoff;
				Z[1] = h[z+1][x+1]*d/2.0 + yoff;
				Y[1] = (z + 1)*d + zoff;
				R[1] = RGBheight(0, Z[1], yoff);
				G[1] = RGBheight(1, Z[1], yoff);
				B[1] = RGBheight(2, Z[1], yoff);
				draw_line_list(X, Y, Z, R, G, B, n);

				// Make top and right side of nWidth x nHeight square straight - without this it would appear jagged
				if (x == nWidth - 1) { //end of row - draw a vertical line
					X[0] = (x+1)*d + xoff;
					Z[0] = h[z][x+1] * d / 2.0 + yoff;
					Y[0] = z*d + zoff;
					R[0] = RGBheight(0, Z[0], yoff);
					G[0] = RGBheight(1, Z[0], yoff);
					B[0] = RGBheight(2, Z[0], yoff);

					X[1] = (x+1)*d + xoff;
					Z[1] = h[z + 1][x+1] * d / 2.0 + yoff;
					Y[1] = (z + 1)*d + zoff;
					R[1] = RGBheight(0, Z[1], yoff);
					G[1] = RGBheight(1, Z[1], yoff);
					B[1] = RGBheight(2, Z[1], yoff);
					draw_line_list(X, Y, Z, R, G, B, n);
				}
				if (z == nHeight - 1) { //top of column - draw horizontal line
					X[0] = x*d + xoff;
					Z[0] = h[z+1][x] * d / 2.0 + yoff;
					Y[0] = (z+1)*d + zoff;
					R[0] = RGBheight(0, Z[0], yoff);
					G[0] = RGBheight(1, Z[0], yoff);
					B[0] = RGBheight(2, Z[0], yoff);

					X[1] = (x + 1)*d + xoff;
					Z[1] = h[z + 1][x + 1] * d / 2.0 + yoff; 
					Y[1] = (z + 1)*d + zoff;
					R[1] = RGBheight(0, Z[1], yoff);
					G[1] = RGBheight(1, Z[1], yoff);
					B[1] = RGBheight(2, Z[1], yoff);
					draw_line_list(X, Y, Z, R, G, B, n);
				}
			}
	}
	
		return h[nWidth / 2][nHeight / 2];
}				

void draw_3D_graphics()
{
	static bool isOnStart = true, isOnHelp = false, collision = false;

	static double stX[100], stY[100], stR[100], stG[100], stB[100];

	static double X0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0  };
	static double U0[] = {0.0, 0.0 };
	static double la_pt[4] = { 0, 5, 1, 5 };
	static double eye_pt[4] = { 0, 0, 0, 0.0 };
	static double up_pt[4] = {0, 0.0, 0.0, 1.0 }; 
	static int init = 0; // initialization flag
	static double t, dt;
	static double yaw = 0.0, pitch = 0.0, roll = 0.0, roll0 = 0.0, yaw0 = 0.0;
	static double r;
	static double q[10]; 
	static Heli *heli;
	static double Ttr = 0, Tmr = 5.0, d_lon = 0.0, d_lat = 0.0, y0;
	static double dot_p, Cx, Cy, Cz, vx, vy, vz;

	//****** initalization section ******
	if (!init) {
		heli = new Heli(X0, 9, U0, 2, 0.0);
		mwidth = 100;
		mheight = 100;
		scl = 10.0;
		xoff = -0.5*mwidth*scl;
		yoff = -5.0*scl;
		zoff = -0.5*mheight*scl;
		flyingy = 0.0;
		flying = 0.0;

		for (int i = 0; i < 100; i++) {
			stX[i] = 2.0*xoff + abs(xoff)*4.0*ran();
			stY[i] = 2.0*zoff + abs(zoff)*4.0*ran();
			stR[i] = 1.0 - 0.1*ran();
			stG[i] = 1.0 - 0.1*ran();
			stB[i] = 1.0 - 0.1*ran();
		}
		
		init = 1;
	} // end of initialization section
	
	if (KEY('H')) {
		isOnHelp = true;
	}
	else isOnHelp = false;

	if (!isOnStart && !isOnHelp) {
		text_xy("Help - H", 10, 10, 10);
		//drawStars(100, stX, stY, stR, stG, stB);
		//****** Key presses handling ******

		// Thrust
		if (KEY(VK_SPACE)) {
			Tmr = (Tmr + 1.0 > 10.0) ? Tmr : Tmr + 1.0;
		}
		else if (KEY(VK_LCONTROL)) {

			Tmr = (Tmr - 0.5 < 0.0) ? Tmr : Tmr - 0.5;
		}
		else Tmr = (abs(Tmr) <= 5.00) ? 5.00 : Tmr / 1.01;

		// Pitch
		if (KEY('W')) {
			d_lon = (d_lon - 0.1 < -PI / 2) ? -PI / 2 : d_lon - 0.1;
		}
		else if (KEY('S')) {
			d_lon = (d_lon + 0.1 > PI / 2) ? PI / 2 : d_lon + 0.1;
		}
		else d_lon = (abs(d_lon) <= 0.01) ? 0.000 : d_lon / 1.2;

		// Spin Left - Right
		if (KEY(VK_LEFT)) {
			Ttr = (Ttr - 0.1 < -1.0) ? Ttr : Ttr - 0.1;
		}
		else if (KEY(VK_RIGHT)) {
			Ttr = (Ttr + 0.1 > 1.0) ? Ttr : Ttr + 0.1;
		}
		else Ttr = (abs(Ttr) <= 0.01) ? 0.000 : Ttr / 1.2; //decel of Ttr when no key pressed

		// Roll Left - Right
		if (KEY('A')) {
			d_lat = (d_lat + 0.1 > PI) ? PI : d_lat + 0.1;
		}
		else if (KEY('D')) {
			d_lat = (d_lat - 0.1 < -PI) ? -PI : d_lat - 0.1;
		}
		else d_lat = (abs(d_lat) <= 0.001) ? 0.000 : d_lat / 1.2;

		//Update time
		dt = high_resolution_time() - t;
		t = high_resolution_time();

		//Simulate
		heli->sim_step(dt, Ttr, Tmr, d_lat, d_lon, q, collision);

		//****** Update view params ******

		yaw = q[5];
		pitch = PI / 2 - q[4];
		roll = q[3];
		yoff += scl*q[2] * dt;

		// Orient look-at point based on yaw and pitch
		la_pt[1] = 1.0*cos(yaw)*sin(pitch);
		la_pt[2] = 1.0*sin(yaw)*sin(pitch);
		la_pt[3] = 1.0*cos(pitch);

		// Rotate up-direction using Rodrigue's formula (this will apply any roll to the view)
		dot_p = la_pt[1] * up_pt[1] + la_pt[2] * up_pt[2] + la_pt[3] * up_pt[3];
		Cx = up_pt[2] * la_pt[3] - up_pt[3] * la_pt[2];
		Cy = up_pt[3] * la_pt[1] - up_pt[1] * la_pt[3];
		Cz = up_pt[1] * la_pt[2] - up_pt[2] * la_pt[1];

		up_pt[1] = up_pt[1] * cos(roll - roll0) + Cx*sin(roll - roll0) + dot_p*la_pt[1] * (1 - cos(roll - roll0));
		up_pt[2] = up_pt[2] * cos(roll - roll0) + Cy*sin(roll - roll0) + dot_p*la_pt[2] * (1 - cos(roll - roll0));
		up_pt[3] = up_pt[3] * cos(roll - roll0) + Cz*sin(roll - roll0) + dot_p*la_pt[3] * (1 - cos(roll - roll0));
		roll0 = roll;

		// Rotate up-direction about Z-axis by delta yaw (yaw - yaw0) - necessary step or up-direction will not rotate for yaw changes
		vx = up_pt[1];
		vy = up_pt[2];
		vz = up_pt[3];

		up_pt[1] = vx * cos(yaw - yaw0) - vy*sin(yaw - yaw0);
		up_pt[2] = vx*sin(yaw - yaw0) + vy*cos(yaw - yaw0);
		up_pt[3] = vz;
		yaw0 = yaw;

		// Move terrain to look like we are moving (in reality we are always stuck at (0,0,0))
		flyingy += (q[0] * dt*cos(yaw) - q[1] * dt*sin(yaw))*0.05;
		flying += (q[0] * dt*sin(yaw) + q[1] * dt*cos(yaw))*0.05;

		//draw
		set_view(eye_pt, la_pt, up_pt);
		y0 = drawTerrain(scl, mwidth, mheight);
		//check collision
		if (y0*scl/2.0 + yoff + 3.0 > 0.0) collision = true;
		else collision = false;
	}
	else if(isOnStart && !isOnHelp){
		// start menu

		text_xy("Press Enter to Begin", WIDTH_MAX/2-120, HEIGHT_MAX/2-50, 20);
		text_xy("Help - H", 10, 10, 10);
		up_pt[1] = 0.0;
		up_pt[2] = 0.0;
		up_pt[3] = 1.0;

		dt = high_resolution_time() - t;
		t = high_resolution_time();

		flyingy -= 0.008;

		yaw = 0.0;
		pitch = PI / 2;
		roll = 0.0;

		la_pt[1] = 1.0;
		la_pt[2] = 0.0;
		la_pt[3] = 0.0;

		set_view(eye_pt, la_pt, up_pt);
		y0 = drawTerrain(scl, mwidth, mheight);
		if (y0*scl / 2.0 + yoff > -1.0) {
			yoff -= 2.0;
		}
		if (y0*scl / 2.0 + yoff > -5.0) {
			yoff -= 8.0;
		}


		if (y0)
		if (KEY(VK_RETURN)) {
			isOnStart = false;
		}
	}
	else {
		dt = high_resolution_time() - t;
		t = high_resolution_time();
		// help menu
		helpMenu();
	}

}

