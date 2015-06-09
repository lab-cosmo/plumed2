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
#ifndef __PLUMED_tools_GridSearch_h
#define __PLUMED_tools_GridSearch_h

#include "MinimiseBase.h"
#include "ConjugateGradient.h"
#include <iostream>
#include "Interpolate.h"
#include <math.h>
namespace PLMD{

template <class FCLASS>
class GridSearch : public MinimiseBase<FCLASS>{
private:
/// This is the pointer to the member funciton in the energy
/// calculating class that calculates the energy
	typedef double(FCLASS::*engf_pointer)( const std::vector<double>& p, std::vector<double>& der );
	const unsigned ITMAX;
	const double EPS;
	FCLASS* myclass_func;
public:
	GridSearch( FCLASS* funcc ) : MinimiseBase<FCLASS>(funcc), ITMAX(200), EPS(1E-10) {}
	void minimise( std::vector<double>& p, engf_pointer myfunc,std::vector< std::vector<double> > &stressgrid,const int &nfgrid,double &min_eng,std::vector<double> &ptsinx,std::vector<double> &ptsiny  ); 
};

template <class FCLASS>

void GridSearch<FCLASS>::minimise( std::vector<double> &p, engf_pointer myfunc,std::vector< std::vector<double> > &stressgrid,const int &nfgrid,double &min_eng,std::vector<double> &ptsinx,std::vector<double> &ptsiny ){

	std::vector<double> der( p.size() );
	std::vector<double> temp( p.size() );
	std::vector<double> grid_point( p.size() );
	double eng_pt;
	int xlimit = stressgrid.size();
	int ylimit = stressgrid[0].size();
	std::vector<double> y(4);
	std::vector<double> y1(4);
	std::vector<double> y2(4);
	std::vector<double> y12(4);
	for(int i=1;i<xlimit-2;i++){
	  for(int j=1;j<ylimit -2;j++){
		  y[0] = stressgrid[i][j];
		  y[1] = stressgrid[i+1][j];
		  y[2] = stressgrid[i+1][j+1];
		  y[3] = stressgrid[i][j+1];
		  
		  y1[0] = (stressgrid[i+1][j] - stressgrid[i-1][j])/(ptsinx[i+1] - ptsinx[i-1]);
		  y1[1] = (stressgrid[i+2][j] - stressgrid[i][j])/(ptsinx[i+2] - ptsinx[i]);
		  y1[2] = (stressgrid[i+2][j+1] - stressgrid[i][j])/(ptsinx[i+2] - ptsinx[i]);
		  y1[3] = (stressgrid[i+1][j+1] - stressgrid[i-1][j+1])/(ptsinx[i+1] - ptsinx[i-1]);

		  y2[0] = (stressgrid[i][j+1] - stressgrid[i][j-1])/(ptsiny[j+1] - ptsiny[j-1]);
		  y2[1] = (stressgrid[i+1][j+1] - stressgrid[i+1][j-1])/(ptsiny[j+1] - ptsiny[j-1]);
		  y2[2] = (stressgrid[i+1][j+2] - stressgrid[i+1][j])/(ptsiny[j+2] - ptsiny[j]);
		  y2[3] = (stressgrid[i][j+2] - stressgrid[i][j])/(ptsiny[j+2] - ptsiny[j]);

		  y12[0] = (stressgrid[i+1][j+1] - stressgrid[i+1][j-1] - stressgrid[i-1][j+1] + stressgrid[i-1][j-1] )/(ptsinx[i+1] - ptsinx[i-1])*(ptsiny[j+1] - ptsiny[j-1]);
		  y12[1] = (stressgrid[i+2][j+1] - stressgrid[i+2][j-1] - stressgrid[i][j+1] + stressgrid[i][j-1] )/(ptsinx[i+2] - ptsinx[i])*(ptsiny[j+1] - ptsiny[j-1]);
		  y12[2] = (stressgrid[i+2][j+2] - stressgrid[i+2][j] - stressgrid[i][j+2] + stressgrid[i][j] )/(ptsinx[i+2] - ptsinx[i])*(ptsiny[j+2] - ptsiny[j]);
		  y12[3] = (stressgrid[i+1][j+2] - stressgrid[i+1][j] - stressgrid[i-1][j+2] + stressgrid[i-1][j] )/(ptsinx[i+1] - ptsinx[i-1])*(ptsiny[j+2] - ptsiny[j]);
		  

		  double x1l = ptsinx[i];
		  double x1u = ptsinx[i+1];
		  double x2l = ptsiny[j];
		  double x2u = ptsiny[j+1];
		  int ncgrid=xlimit;
		  int nmult=nfgrid/ncgrid;
		  double dx=(x1l-x1u)/(double)nmult;
		  double dy=(x2l-x2u)/(double)nmult;
		  for (unsigned ii=0;ii<nmult;ii++){
		    double x1 = x1l+ii*dx;
		    for (unsigned jj=0;jj<nmult;jj++){
		      double x2 = x2l+jj*dy;
		      double ansy,ansy1,ansy2;
			  bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,ansy2); //Can be turned off by just commenting out this line.
		      if(ansy < min_eng && ansy >0){
			     min_eng = ansy;
			     p[0] = x1;p[1] = x2;
		         }
		      }
		    
	      }
	  }

	}


		  
//Flower Search (Bruteforce search in fixed directions):

/*   
	for(unsigned i=0;i<p.size();i++) {
		grid_point[i] = p[i];
		temp[i] = p[i];
	}
    
    min_eng = this->calcDerivatives(grid_point,der,myfunc);
	for(int i=-5;i<5;i++){	  
		for(unsigned k=0;k<p.size();k++) grid_point[k] = temp[k]+0.25*i;
		double engy = this->calcDerivatives(grid_point,der,myfunc); 
		if(engy<min_eng){
			min_eng = engy;
			//for(unsigned j=0;j<p.size();j++) p[j] = grid_point[j];
		}	  
		for(unsigned k=0;k<p.size();k+=2) grid_point[k] = temp[k]+0.25*i;
		engy = this->calcDerivatives( grid_point, der, myfunc );
		if(engy<min_eng){
			min_eng = engy;
			//for(unsigned j=0;j<p.size();j++) p[j] = grid_point[j];
		}	  
		for(unsigned k=1;k<p.size();k+=2) grid_point[k] = temp[k]+0.25*i;
		engy = this->calcDerivatives( grid_point, der, myfunc );
		if(engy<min_eng){
			min_eng = engy;
			//for(unsigned j=0;j<p.size();j++) p[j] = grid_point[j];
		}	  
	}
	
*/	

}
  
}
#endif

