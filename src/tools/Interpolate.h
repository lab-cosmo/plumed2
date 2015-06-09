/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#ifndef __PLUMED_tools_Interpolate_h
#define __PLUMED_tools_Interpolate_h
#include <algorithm>
//Following codes have been used from Numerical Recipes in C
//~ 
//~ struct Base_interp{
	//~ //Abstract base class used by all interpolation routines in this chapter. Only the routine interp
	//~ //is called directly by the user.
 //~ int n, mm, jsav, cor, dj;
 //~ const double *xx, *yy;
 //~ Base_interp(vector<double> &x, const double *y, int m) : n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y) {
 //~ dj = MIN(1,(int)pow((double)n,0.25));
 //~ }
//~ 
//~ double interp(double x) {
//~ //Given a value x, return an interpolated value, using data pointed to by xx and yy.
 //~ int jlo = cor ? hunt(x) : locate(x);
 //~ return rawinterp(jlo,x);
//~ }
//~ 
 //~ int locate(const double x);
 //~ int hunt(const double x);
 //~ double virtual rawinterp(int jlo, double x) = 0;
 //~ //Derived classes provide this as the actual interpolation method.
//~ 
//~ };
//~ 
//~ int Base_interp::locate(const double x)
//~ //Given a value x, return a value j such that x is (insofar as possible) centered in the subrange
//~ //xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either
//~ //increasing or decreasing. The returned value is not less than 0, nor greater than n-1.
//~ {
 //~ int ju,jm,jl;
 //~ bool ascnd = (xx[n-1] >= xx[0]);
 //~ jl=0;
 //~ ju=n-1;
 //~ 
//~ while (ju-jl > 1){
 //~ jm = (ju+jl) >> 1;
 //~ if (x >= xx[jm] == ascnd)
 //~ jl=jm;
 //~ else
 //~ ju=jm;
//~ }
 //~ cor = abs(jl-jsav) > dj ? 0 : 1;
 //~ jsav = jl;
//~ return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
//~ }
//~ 
//~ int Base_interp::hunt(const double x)
//~ Given a value x, return a value j such that x is (insofar as possible) centered in the subrange
//~ xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either
//~ increasing or decreasing. The returned value is not less than 0, nor greater than n-1.
//~ {
	//~ int jl=jsav, jm, ju, inc=1;
	//~ //if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
	//~ bool ascnd=(xx[n-1] >= xx[0]);
	//~ if (jl < 0 || jl > n-1) {
	//~ jl=0;
	//~ ju=n-1;
	//~ } 
	//~ else {
	//~ if (x >= xx[jl] == ascnd) {
		//~ for (;;) {
			//~ ju = jl + inc;
			//~ if (ju >= n-1) { ju = n-1; break;}
			//~ else if (x < xx[ju] == ascnd) break;
			//~ else {
				//~ jl = ju;
				//~ inc += inc;
			//~ }
		//~ }
	//~ } 
	//~ else {
	//~ ju = jl;
			//~ for (;;) {
				//~ jl = jl - inc;
				//~ if (jl <= 0) { jl = 0; break;}
				//~ else if (x >= xx[jl] == ascnd) break;
				//~ else {
				//~ ju = jl;
				//~ inc += inc;
				//~ }
			//~ }
		//~ }
	//~ }
	//~ 
	//~ while (ju-jl > 1) {
//~ //Hunt is done, so begin the final bisection phase:
		//~ jm = (ju+jl) >> 1;
		//~ if (x >= xx[jm] == ascnd)
			//~ jl=jm;
		//~ else
			//~ ju=jm;
	//~ }
	//~ cor = abs(jl-jsav) > dj ? 0 : 1;
	//~ Decide whether to use hunt or locate next
	//~ jsav = jl;
	//~ time.
	//~ return max(0,min(n-mm,jl-((mm-2)>>1)));
//~ }
//~ 
//~ 
//~ struct Poly_interp : Base_interp
 //~ //Polynomial interpolation object. Construct with x and y vectors, and the number M of points
 //~ //to be used locally (polynomial order plus one), then call interp for interpolated values.
//~ {
	//~ double dy;
	//~ Poly_interp(vector<double> &xv, VecDoub_I &yv, int m) : Base_interp(xv,&yv[0],m), dy(0.) {}
	//~ double rawinterp(int jl, double x);
//~ };
//~ 
//~ 
//~ double Poly_interp::rawinterp(int jl, double x)
//~ //Given a value x, and using pointers to data xx and yy, this routine returns an interpolated
//~ //value y, and stores an error estimate dy. The returned value is obtained by mm-point polynomial
//~ //interpolation on the subrange xx[jl..jl+mm-1].
//~ {
	//~ int i,m,ns=0;
	//~ double  y,den,dif,dift,ho,hp,w;
	//~ const double *xa = &xx[jl], *ya = &yy[jl];
	//~ vecot<double> c(mm),d(mm);
	//~ dif=abs(x-xa[0]);
	//~ for(i=0;i<mm;i++) {
			//~ Here we find the index ns of the closest table entry,
			//~ if ((dift=abs(x-xa[i])) < dif) {
			//~ ns=	i;
			//~ dif=dift;
			//~ }
		//~ c[i]=ya[i];
		//~ and initialize the tableau of c’s and d’s.
		//~ d[i]=ya[i];
	//~ }	
	//~ y=ya[ns--];
	//~ for (m=1;m<mm;m++) {
		//~ for (i=0;i<mm-m;i++) {
			//~ ho=xa[i]-x;
			//~ hp=xa[i+m]-x;
			//~ w=c[i+1]-d[i];
			//~ den=w/den;
			//~ d[i]=hp*den;
			//~ c[i]=ho*den;
		//~ }
	 //~ y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
	//~ }
//~ 
 //~ return y;
//~ }




void bcucof(std::vector<double> &y, std::vector<double> &y1, std::vector<double> &y2, std::vector<double> &y12,const double d1, const double d2, std::vector<std::vector <double> > &c) {
//~ Given arrays y[0..3], y1[0..3], y2[0..3], and y12[0..3], containing the function, gradients,
//~ and cross-derivative at the four grid points of a rectangular grid cell (numbered counterclockwise
//~ from the lower left), and given d1 and d2, the length of the grid cell in the 1 and 2 directions, this
//~ routine returns the table c[0..3][0..3] that is used by routine bcuint for bicubic interpolation.
	//static Int wt_d[16*16]=
	int l,k,j,i;
	double xx,d1d2=d1*d2;
	std::vector<double> cl(16),x(16);
	//static MatInt wt(16,16,wt_d);
	
	static int wt[16][16] = {
	{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
	{-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0},
	{2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0},
	{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
	{0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1},
	{0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1},
	{-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0},
	{9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2},
	{-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2},
	{2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0},
	{-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1},
	{4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1}};
	for (i=0;i<4;i++) {
		x[i]=y[i];
		x[i+4]=y1[i]*d1;
		x[i+8]=y2[i]*d2;
		x[i+12]=y12[i]*d1d2;
	}
	for (i=0;i<16;i++){
	//Matrix-multiply by the stored table.
	xx=0.0;
	for (k=0;k<16;k++) xx += wt[i][k]*x[k];
	cl[i]=xx;
	}
	l=0;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++) c[i][j]=cl[l++];
}

void bcuint(std::vector<double> &y, std::vector<double> &y1, std::vector<double> &y2, std::vector<double> &y12,const double x1l, const double x1u, const double x2l, const double x2u,const double x1, const double x2, double &ansy, double &ansy1, double &ansy2){
	int i;
	double t,u,d1=x1u-x1l,d2=x2u-x2l;
    std::vector<std::vector<double> > c(4,std::vector<double>(4));
    bcucof(y,y1,y2,y12,d1,d2,c);
    if (x1u == x1l || x2u == x2l)
       printf("Bad input in routine bcuint");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	ansy=ansy2=ansy1=0.0;
	for (i=3;i>=0;i--) {
		//~ Equation (3.6.6).
		ansy=t*ansy+((c[i][3]*u+c[i][2])*u+c[i][1])*u+c[i][0];
		ansy2=t*ansy2+(3.0*c[i][3]*u+2.0*c[i][2])*u+c[i][1];
		ansy1=u*ansy1+(3.0*c[3][i]*t+2.0*c[2][i])*t+c[1][i];
	}	
	ansy1 /= d1;
	ansy2 /= d2;
}

#endif


