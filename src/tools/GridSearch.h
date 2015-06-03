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
	void minimise( std::vector<double>& p, engf_pointer myfunc,const double minx,const double miny,const double maxx,const double maxy ); 
};

template <class FCLASS>

void GridSearch<FCLASS>::minimise( std::vector<double> &p, engf_pointer myfunc,const double minx,const double miny,const double maxx,const double maxy ){
	
	std::vector<double> der( p.size() );
	std::vector<double> temp( p.size() );
	std::vector<double> grid_point( p.size() );
	double min_eng,eng_pt;

  //~ //minx,maxx and miny,maxy define the boundary of the grid we want and the step size.
  //~ //Bicubic Interpolation <Grid Search>
  //std::cout<< "minx "<<minx<<" miny "<<miny<<" maxx "<<maxx<<" maxy "<<maxy<<"\n";
  double stepx = fabs(maxx-minx)/10.0;
  double stepy = fabs(maxy-miny)/10.0;
  //std::cout<< "here stepy "<<stepy<<"stepx "<<stepx<<"\n";
  //~ double curr_x = minx;
  //~ double curr_y = miny;
  //~ std::vector<double> ptsinx;  
  //~ std::vector<double> ptsiny;
  //~ std::vector< std::vector<double> > engmatrix(11,std::vector<double>(11));
  //~ grid_point[0] = curr_x;grid_point[1] = curr_y;
  //~ min_eng = this->calcDerivatives(grid_point,der,myfunc) ;
  
  //~ for(unsigned i=0;curr_x <= maxx;i++){
	  //~ curr_x = curr_x + stepx*i;
	  //~ ptsinx.push_back(curr_x);
	  //~ curr_y = miny;
	  //~ for(unsigned j=0;curr_y<=maxy;j++){
		  //~ curr_y = curr_y + stepy*j;
		  //~ if(i==0) ptsiny.push_back(curr_y);
		  //~ grid_point[0] = curr_x;grid_point[1] = curr_y;
		  //~ eng_pt = this->calcDerivatives(grid_point,der,myfunc);
		  //~ engmatrix[i][j] = eng_pt;
		  //~ if(eng_pt < min_eng){
			  //~ min_eng = eng_pt;
			  //~ p[0] = curr_x;p[1] = curr_y;
		  //~ }		  		  
	  //~ }
	//~ }
  //~ std::cout<<p[0]<<" "<<p[1]<<" "<<min_eng<<"\n";
  
  //std::cout<<"Minimum energy on the coarse grid is "<<min_eng<<"\n";
  //~ // We have stress values on the coarse grid now. We need to implement the fine interpolation to find the true minimum.
  //~ //Definitions of required vairables for bicubic routine.
  //~ std::vector<double> y(4),y1(4),y2(4),y12(4);
  //~ 
  //~ 
  //~ for(int i=1;i<ptsinx.size()-2;i++){
	  //~ for(int j=1;j<ptsiny.size()-2;j++){
		  //~ 
		  //~ y[0] = engmatrix[i][j];
		  //~ y[1] = engmatrix[i+1][j];
		  //~ y[2] = engmatrix[i+1][j+1];
		  //~ y[3] = engmatrix[i][j+1];
		  //~ 
		  //~ y1[0] = (engmatrix[i+1][j] - engmatrix[i-1][j])/(ptsinx[i+1] - ptsinx[i-1]);
		  //~ y1[1] = (engmatrix[i+2][j] - engmatrix[i][j])/(ptsinx[i+2] - ptsinx[i]);
		  //~ y1[2] = (engmatrix[i+2][j+1] - engmatrix[i][j])/(ptsinx[i+2] - ptsinx[i]);
		  //~ y1[3] = (engmatrix[i+1][j+1] - engmatrix[i-1][j+1])/(ptsinx[i+1] - ptsinx[i-1]);
		  //~ 
		  //~ y2[0] = (engmatrix[i][j+1] - engmatrix[i][j-1])/(ptsiny[j+1] - ptsiny[j-1]);
		  //~ y2[1] = (engmatrix[i+1][j+1] - engmatrix[i+1][j-1])/(ptsiny[j+1] - ptsiny[j-1]);
		  //~ y2[2] = (engmatrix[i+1][j+2] - engmatrix[i+1][j])/(ptsiny[j+2] - ptsiny[j]);
		  //~ y2[3] = (engmatrix[i][j+2] - engmatrix[i][j])/(ptsiny[j+2] - ptsiny[j]);
		  //~ 
		  //~ y12[0] = (engmatrix[i+1][j+1] - engmatrix[i+1][j-1] - engmatrix[i-1][j+1] + engmatrix[i-1][j-1] )/(ptsinx[i+1] - ptsinx[i-1])*(ptsiny[j+1] - ptsiny[j-1]);
		  //~ y12[1] = (engmatrix[i+2][j+1] - engmatrix[i+2][j-1] - engmatrix[i][j+1] + engmatrix[i][j-1] )/(ptsinx[i+2] - ptsinx[i])*(ptsiny[j+1] - ptsiny[j-1]);
		  //~ y12[2] = (engmatrix[i+2][j+2] - engmatrix[i+2][j] - engmatrix[i][j+2] + engmatrix[i][j] )/(ptsinx[i+2] - ptsinx[i])*(ptsiny[j+2] - ptsiny[j]);
		  //~ y12[3] = (engmatrix[i+1][j+2] - engmatrix[i+1][j] - engmatrix[i-1][j+2] + engmatrix[i-1][j] )/(ptsinx[i+1] - ptsinx[i-1])*(ptsiny[j+2] - ptsiny[j]);
		  //~ 
		  //~ double x1l = ptsinx[i];
		  //~ double x1u = ptsinx[i+1];
		  //~ double x2l = ptsiny[j];
		  //~ double x2u = ptsiny[j+1];
	      //~ double x1 = (x1l+x1u)/2;
	      //~ double x2 = (x2l+x2u)/2;
	      //~ double ansy,ansy1,ansy2;
	      //~ bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,ansy2);
	      //~ //std::cout<<"hereisit\n";
		  //~ if(ansy < min_eng){
			//~ min_eng = ansy;
			//~ p[0] = x1;p[1] = x2;
		  //~ }
	      //~ 
	  //~ }
  //~ }
		  
//Flower Search :   
	//~ for(unsigned i=0;i<p.size();i++) {
		//~ grid_point[i] = p[i];
		//~ temp[i] = p[i];
	//~ }
    //~ 
    //~ min_eng = this->calcDerivatives(grid_point,der,myfunc);
    //~ 
   //~ // std::cout<<"Minimum energy for this frame is "<<min_eng<<"\n";
	//~ for(int i=-5;i<5;i++){	  
		//~ for(unsigned k=0;k<p.size();k++) grid_point[k] = temp[k]+0.25*i;
		//~ double engy = this->calcDerivatives(grid_point,der,myfunc); 
		//~ if(engy<min_eng){
			//~ min_eng = engy;
			//~ //for(unsigned j=0;j<p.size();j++) p[j] = grid_point[j];
		//~ }	  
		//~ for(unsigned k=0;k<p.size();k+=2) grid_point[k] = temp[k]+0.25*i;
		//~ engy = this->calcDerivatives( grid_point, der, myfunc );
		//~ if(engy<min_eng){
			//~ min_eng = engy;
			//~ //for(unsigned j=0;j<p.size();j++) p[j] = grid_point[j];
		//~ }	  
		//~ for(unsigned k=1;k<p.size();k+=2) grid_point[k] = temp[k]+0.25*i;
		//~ engy = this->calcDerivatives( grid_point, der, myfunc );
		//~ if(engy<min_eng){
			//~ min_eng = engy;
			//~ //for(unsigned j=0;j<p.size();j++) p[j] = grid_point[j];
		//~ }	  
	//~ }
	//std::cout<<"Minimum energy for moved frame is "<<min_eng<<"\n";
	

}
  
}
#endif

