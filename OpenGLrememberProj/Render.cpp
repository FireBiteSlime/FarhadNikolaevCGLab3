

#include "Render.h"

#include <Windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>
#include <corecrt_math.h>
#include <iostream>




#define PI 3.14159265

double Calculating_Vector_Angle(double a1, double b1, double c1, double a2, double b2, double c2) {
	double vector_sum = (a1 * a2) + (b1 * b2) + (c1 * c2);
	double mod1 = sqrt( a1 * a1 + b1 * b1 + c1 * c1);
	double mod2 = sqrt( a2 * a2 + b2 * b2 + c2 * c2);
	double cos_angle = vector_sum / (mod1 * mod2);
	double angle = acos(cos_angle) * 180.0 / PI;

	return angle;
}


double Calculating_Length_Points(double a1, double b1, double c1, double a2, double b2, double c2)
{
	double res;
	res = sqrt(((a2 - a1) * (a2 - a1)) + ((b2 - b1) * (b2 - b1)) + ((c2 - c1) * (c2 - c1)));
	return res;
}

double* Calculating_Normal(double a1, double b1, double c1)
{
	double *normal = new double[3];
	double locLength = Calculating_Length_Points(0, 0, 0, a1, b1, c1);
	double inv_length = (1 / locLength);
	normal[0] = a1 * inv_length;
	normal[1] = b1 * inv_length;
	normal[2] = c1 * inv_length;
	return normal;
}



double Calculating_First_Formula(double p1, double p2, double p3, double p4, double t)
{
	return (1 - t) * (1 - t) * (1 - t) * p1 + 3 * t * (1 - t) * (1 - t) * p2 + 3 * t * t * (1 - t) * p3 + t * t * t * p4;
}


double* Calculating_Second_Formula(double* p1, double* p2, double* p3, double* p4, double t){
	double* Trace = new double[3];
	Trace[0] = -3 * p1[0] + 6 * p1[0] * t - 3 * p1[0] * t * t + 3 * p2[0] - 12 * p2[0] * t + 9 * p2[0] * t * t + 6 * p3[0] * t - 9 * p3[0] * t * t - 3 * p4[0] * t * t;
	Trace[1] = -3 * p1[1] + 6 * p1[1] * t - 3 * p1[1] * t * t + 3 * p2[1] - 12 * p2[1] * t + 9 * p2[1] * t * t + 6 * p3[1] * t - 9 * p3[1] * t * t - 3 * p4[1] * t * t;
	Trace[2] = -3 * p1[2] + 6 * p1[2] * t - 3 * p1[2] * t * t + 3 * p2[2] - 12 * p2[2] * t + 9 * p2[2] * t * t + 6 * p3[2] * t - 9 * p3[2] * t * t - 3 * p4[2] * t * t;
	return Trace;
}

double Calculating_Hermite_Formula(double p1, double p4, double r1, double r4, double t)
{
	return p1 * (2 * t * t * t * -3 * t * t + 1) + p4 * (-2 * t * t * t + 3 * t * t) + r1 * (t * t * t - 2 * t * t + t) + r4*(t * t * t - t * t);
}
double Factorial(double num) {
	double res=1;
	for (int i = 1; i <= num; i++) {
		res *= i;
	}
	return res;
}

double Calculating_Bernoulli_Formula(double n, double i , double u ) {
	return ((Factorial(n)/(Factorial(i)* Factorial(n-i)))*pow(u,i)*pow((1-u),(n-i)));
}

double* Calculating_Third_Formula( double P[4][4][3], double u, double v) {

	double* Trace = new double[3];
	Trace[0] = 0;
	Trace[1] = 0;
	Trace[2] = 0;
	double P2[3] = { P[1][3][0],P[1][3][1],P[1][3][2] };
	for (int i = 0; i <= 3; i++) {
		for (int j = 0; j <= 3; j++) {
			Trace[0] += Calculating_Bernoulli_Formula(3, i, u) * Calculating_Bernoulli_Formula(3, j, v) * P[i][j][0];
			Trace[1] += Calculating_Bernoulli_Formula(3, i, u) * Calculating_Bernoulli_Formula(3, j, v) * P[i][j][1];
			Trace[2] += Calculating_Bernoulli_Formula(3, i, u) * Calculating_Bernoulli_Formula(3, j, v) * P[i][j][2];
		}
	}
	return Trace;
}



double t_max =0;
int flag = 0;

void Shuriken() {
	double A1[] = { 0,1,0 };
	double A2[] = { -0.2,0.8,0 };
	double A3[] = { -0.2,-0.8,0 };
	double A4[] = { 0,-1,0 };
	double A5[] = { 0.2,-0.8,0 };
	double A6[] = { 0.2,0.8,0 };

	double B1[] = { 0.2,0.2,0 };
	double B2[] = { 0.2,-0.2,0 };
	double B3[] = { 1,-0.2,0 };
	double B4[] = { 1,0,0 };

	double B5[] = { -0.2,0.2,0 };
	double B6[] = {-0.2,-0.2,0 };
	double B7[] = { -1,-0.2,0 };
	double B8[] = { -1,0,0 };

	double C1[] = { 0,-1,0.4 };
	double C2[] = { 0,-0.8,0.4 };
	double C3[] = { 0,-0.6,0 };
	double C4[] = { 0,-1,0 };

	
	glColor3d(0.3, 0.3, 0.3);
	glBegin(GL_TRIANGLES);

	glVertex3d(0, 0, 0);
	glVertex3d(-0.6, 0, 0);
	glVertex3d(-0.6, 0.5, 0);

	glVertex3d(0, 0, 0);
	glVertex3d(0, -0.6, 0);
	glVertex3d(-0.5, -0.6, 0);

	glVertex3d(0, 0, 0);
	glVertex3d(0.6, 0, 0);
	glVertex3d(0.6, -0.5, 0);

	glVertex3d(0, 0, 0);
	glVertex3d(0, 0.6, 0);
	glVertex3d(0.5, 0.6, 0);

	glEnd();
	
}

void Animation(double Point1[], double Point2[], double Point3[], double Point4[], double t_max) {

	double* Trace = new double[3];
	double* Vector = new double[3];
	Trace = Calculating_Second_Formula(Point1, Point2, Point3, Point4, t_max);

	Vector = Calculating_Normal(Trace[0], Trace[1], Trace[2]);


	Trace[0] = Calculating_First_Formula(Point1[0], Point2[0], Point3[0], Point4[0], t_max);
	Trace[1] = Calculating_First_Formula(Point1[1], Point2[1], Point3[1], Point4[1], t_max);
	Trace[2] = Calculating_First_Formula(Point1[2], Point2[2], Point3[2], Point4[2], t_max);
	int angle_key = 0;


	double angle = Calculating_Vector_Angle(0, 1, 0, Vector[0], Vector[1], 0);
	if (Vector[0] < 0) {
		if (angle_key == 0)
			angle_key = 1;
		else if (angle_key == 1)
			angle_key = 0;
	};
	if (angle_key == 1) angle = angle * 1;
	if (angle_key == 0) angle = angle * -1;
	glPushMatrix();
	glTranslated(Trace[0], Trace[1], Trace[2]);
	glRotated(angle, 0, 0, 1);

	angle_key = 1;
	angle = Calculating_Vector_Angle(0, 1, 0, 0, Vector[1], Vector[2]);
	if (Vector[2] < 0) {
		if (angle_key == 0)
			angle_key = 1;
		else if (angle_key == 1)
			angle_key = 0;
	};
	if (angle_key == 1) angle = angle * 1;
	if (angle_key == 0) angle = angle * -1;
	glRotated(angle, 1, 0, 0);
	glRotated(angle * 100, 0, 0, 1);
	Shuriken();
	
	glPopMatrix();
}

void Bezier_Curve(double delta_time, double P1[], double P2[], double P3[], double P4[], bool animation_option) {
	if (flag == 0) t_max += delta_time / 10;
	if (t_max > 1) flag = 1;
	if (flag == 1) t_max -= delta_time / 10; 
	if (t_max < 0) flag = 0;



	double* Point1 = new double[3];
	double* Point2 = new double[3];
	double* Point3 = new double[3];
	double* Point4 = new double[3];

	Point1 = P1;
	Point2 = P2;
	Point3 = P3;
	Point4 = P4;

	double t_res = 0.001;

	double* Trace = new double[3];
	double Vector_Points[3];
	

	Vector_Points[0] = Calculating_First_Formula(Point1[0], Point2[0], Point3[0], Point4[0], t_max + t_res);
	Vector_Points[1] = Calculating_First_Formula(Point1[1], Point2[1], Point3[1], Point4[1], t_max + t_res);
	Vector_Points[2] = Calculating_First_Formula(Point1[2], Point2[2], Point3[2], Point4[2], t_max + t_res);

	Trace[0] = Calculating_First_Formula(Point1[0], Point2[0], Point3[0], Point4[0], t_max);
	Trace[1] = Calculating_First_Formula(Point1[1], Point2[1], Point3[1], Point4[1], t_max);
	Trace[2] = Calculating_First_Formula(Point1[2], Point2[2], Point3[2], Point4[2], t_max);

	

	Vector_Points[0] = Calculating_First_Formula(Point1[0], Point2[0], Point3[0], Point4[0], t_max + t_res);
	Vector_Points[1] = Calculating_First_Formula(Point1[1], Point2[1], Point3[1], Point4[1], t_max + t_res);
	Vector_Points[2] = Calculating_First_Formula(Point1[2], Point2[2], Point3[2], Point4[2], t_max + t_res);

	while (Calculating_Length_Points(Trace[0], Trace[1], Trace[2], Vector_Points[0], Vector_Points[1], Vector_Points[2]) <= 1) {
		Vector_Points[0] = Calculating_First_Formula(Point1[0], Point2[0], Point3[0], Point4[0], t_max + t_res);
		Vector_Points[1] = Calculating_First_Formula(Point1[1], Point2[1], Point3[1], Point4[1], t_max + t_res);
		Vector_Points[2] = Calculating_First_Formula(Point1[2], Point2[2], Point3[2], Point4[2], t_max + t_res);
		t_res += 0.001;
	}


	glColor3d(0, 0, 1);
	glBegin(GL_LINE_STRIP);
	glVertex3dv(Point1);
	glVertex3dv(Point2);
	glVertex3dv(Point3);
	glVertex3dv(Point4);
	glEnd();
	glLineWidth(3);
	glColor3d(0, 1, 0);
	glBegin(GL_LINE_STRIP);

	double animation_time = 1;

	if (animation_option) {
		animation_time = t_max;
	}
	for (double t = 0; t <= 1; t += 0.001)
	{
		Trace[0] = Calculating_First_Formula(Point1[0], Point2[0], Point3[0], Point4[0], t * animation_time );
		Trace[1] = Calculating_First_Formula(Point1[1], Point2[1], Point3[1], Point4[1], t * animation_time );
		Trace[2] = Calculating_First_Formula(Point1[2], Point2[2], Point3[2], Point4[2], t * animation_time );
		glVertex3dv(Trace);
	}
	glEnd();


	glColor3d(1, 0, 1);
	glLineWidth(1);
	glPointSize(10);
	glColor3d(1, 0, 0);
	glBegin(GL_POINTS);
	glVertex3dv(Point1);
	glVertex3dv(Point2);
	glVertex3dv(Point3);
	glVertex3dv(Point4);
	glEnd();

	if (animation_option) {
		Animation(Point1, Point2, Point3, Point4, t_max);
	}
}



void Hermite_Curve(double P1[], double P2[], double P3[], double P4[]) {
	
	glPopMatrix();
	glPushMatrix();
	

	double* Ermit_Point1 = new double[3];
	double* Ermit_Point2 = new double[3];
	double* Ermit_Point3 = new double[3];
	double* Ermit_Point4 = new double[3];
	Ermit_Point1 = P1;
	Ermit_Point2 = P2;
	Ermit_Point3 = P3;
	Ermit_Point4 = P4;
	
	double Trace[3];
	glColor3d(1, 0, 0);

	glBegin(GL_LINE_STRIP);

	for (double t = 0; t <= 1; t += 0.001)
	{
		Trace[0] = Calculating_Hermite_Formula(Ermit_Point1[0], Ermit_Point2[0], Ermit_Point3[0], Ermit_Point4[0], t);
		Trace[1] = Calculating_Hermite_Formula(Ermit_Point1[1], Ermit_Point2[1], Ermit_Point3[1], Ermit_Point4[1], t);
		Trace[2] = Calculating_Hermite_Formula(Ermit_Point1[2], Ermit_Point2[2], Ermit_Point3[2], Ermit_Point4[2], t);
		glVertex3dv(Trace);
	}
	glEnd();


	glColor3d(0, 0, 0);
	glBegin(GL_LINE_STRIP);
	glVertex3dv(Ermit_Point1);
	glVertex3dv(Ermit_Point3);
	glEnd();

	glBegin(GL_LINE_STRIP);
	glVertex3dv(Ermit_Point2);
	glVertex3dv(Ermit_Point4);
	glEnd();
	glColor3d(0, 1, 0);
	glBegin(GL_POINTS);
	glVertex3dv(Ermit_Point1);
	glVertex3dv(Ermit_Point2);
	glEnd();
	glPopMatrix();
}

void Bezier_Surface() {

	double* Trace = new double[3];
	glPushMatrix();
	glTranslated(-18, 0, 3);
	double Points[4][4][3] = {
		{{0,9,1}, {3,11,0}, {6,9,0}, {10,9,1}},
		{{0,7,0}, {3,6,0}, {6,6,4}, {9,6,0}},
		{{4,3,0}, {3,4,3}, {7,3,1}, {9,4,0}},
		{{0,0,1}, {3,0,5}, {9,0,0}, {12,0,1}},
	};
	int i = 0;
	int j = 0;
	glPointSize(6);
	glColor3d(1, 0, 0);
	glBegin(GL_POINTS);
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			glVertex3dv(Points[i][j]);
		}
		glVertex3dv(Points[i][j]);
	}
	glEnd();

	i = 0;
	j = 0;
	glColor3d(0, 0, 0);
	
	while (i<4) {
		glBegin(GL_LINE_STRIP);
		while (j<4){
			glVertex3dv(Points[i][j]);
			j += 1;
		}
		j = 0;
		i += 1;
		glEnd();
	}
	i = 0;
	j = 0;
	while (i < 4) {
		glBegin(GL_LINE_STRIP);
		while (j < 4) {
			glVertex3dv(Points[j][i]);
			j += 1;
		}
		j = 0;
		i += 1;
		glEnd();
	}



	glPointSize(4);
	glColor3d(1, 0, 0);
	glBegin(GL_POINTS);
	
	 i = 0;
	 j = 0;
	for (double u = 0; u <= 1; u += 0.1) {
		j = 0;
		for (double v = 0; v <= 1; v += 0.1) {
			Trace = Calculating_Third_Formula(Points,u,v);
			
			glVertex3dv(Trace);

			j += 1;
		}
		i += 1;
	}
	glEnd();
	glPopMatrix();

	glTranslated(-18, 0, 0.5);
	glColor3d(1, 1, 0);


	double u = 0;
	double v = 0;
	glColor3d(0,1, 1);

	while (u <= 1) {
		glBegin(GL_LINE_STRIP);
		while (v <= 1) {
			Trace = Calculating_Third_Formula(Points, u, v);
			
			glVertex3dv(Trace); 
			v += 0.1;
		}
		v = 0;
		u += 0.1;
		glEnd();
	}
	u = 0;
	v = 0;
	while (u <= 1) {
		glBegin(GL_LINE_STRIP);
		while (v <= 1) {
			Trace = Calculating_Third_Formula(Points, v,u);
			
			glVertex3dv(Trace); 
			v += 0.1;
		}
		v = 0;
		u += 0.1;
		glEnd();
	}

}


void Render(double delta_time)
{
	double P1[] = { -18,21,0 };
	double P2[] = { -2,7,10 };
	double P3[] = { -4,9,5 };
	double P4[] = { -8,20,6 };


	double P11[] = { 0,0,0 };
	double P22[] = { 3,15,9 };
	double P33[] = { 14,12,4 };
	double P44[] = { -3,10,4 };


	double Ermit_Point1[] = { 0,0,0 };
	double Ermit_Point2[] = { 3,7,4 };
	double Ermit_Point3[] = { 2,1,1 };
	double Ermit_Point4[] = { 3,-12,4 };

	double Ermit_Point11[] = { 0,0,0 };
	double Ermit_Point22[] = { 7,4,1 };
	double Ermit_Point33[] = { 4,3,1.8 };
	double Ermit_Point44[] = { 5,-7,3.12 };

	

	Bezier_Curve(delta_time, P1, P2, P3, P4, true);
	
	Bezier_Curve(delta_time, P11, P22, P33, P44, true);
	


	Bezier_Surface();
	glTranslated(34, 0, 0);
	Hermite_Curve(Ermit_Point1, Ermit_Point2, Ermit_Point3, Ermit_Point4);
	glPopMatrix();
	
	

	glTranslated(14, 0, 0);
	Hermite_Curve(Ermit_Point11, Ermit_Point22, Ermit_Point33, Ermit_Point44);
	glPopMatrix();
	
	
}
	



