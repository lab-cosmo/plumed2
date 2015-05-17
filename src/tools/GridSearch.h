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

namespace PLMD{

template <class FCLASS>
class GridSearch : public MinimiseBase<FCLASS>  {
private:
/// This is the pointer to the member funciton in the energy
/// calculating class that calculates the energy
  typedef double(FCLASS::*engf_pointer)( const std::vector<double>& p, std::vector<double>& der );
  const unsigned ITMAX;
  const double EPS;
  FCLASS* myclass_func;
public:
  GridSearch( FCLASS* funcc ) : MinimiseBase<FCLASS>(funcc), ITMAX(200), EPS(1E-10) {}
  void minimise( std::vector<double>& p, engf_pointer myfunc ); 
};

template <class FCLASS>

void GridSearch<FCLASS>::minimise( std::vector<double> &p, engf_pointer myfunc){
  std::vector<double> der( p.size() );
  std::vector<double> temp( p.size() );
  std::vector<double> grid_point( p.size() );
  
  for(unsigned i=0;i<p.size();i++) {
	  grid_point[i] = p[i];
	  temp[i] = p[i];
  }
      double min_eng = this->calcDerivatives( grid_point, der, myfunc );
  
  for(int i=-5;i<5;i++){
	  
	  
	  for(unsigned k=0;k<p.size();k++) grid_point[k] = temp[k]+0.25*i;
	  double engy = this->calcDerivatives( grid_point, der, myfunc ); 
	  if(engy<min_eng){
		  min_eng = engy;
		  for(unsigned j=0;j<p.size();j++) p[j] = grid_point[j];
	  }	  
	  for(unsigned k=0;k<p.size();k+=2) grid_point[k] = temp[k]+0.25*i;
	  engy = this->calcDerivatives( grid_point, der, myfunc );
	  if(engy<min_eng){
		  min_eng = engy;
		  for(unsigned j=0;j<p.size();j++) p[j] = grid_point[j];
	  }	  
	  for(unsigned k=1;k<p.size();k+=2) grid_point[k] = temp[k]+0.25*i;
	  engy = this->calcDerivatives( grid_point, der, myfunc );
	  if(engy<min_eng){
		  min_eng = engy;
		  for(unsigned j=0;j<p.size();j++) p[j] = grid_point[j];
	  }
	  
  }
	//double cgtol = 0.0001;
	//ConjugateGradient<FCLASS>::minimise(cgtol,p,myfunc )
  	//ConjugateGradient<DimensionalityReductionBase> myminimiser2( this );
	//myminimiser2.minimise( cgtol, p, &DimensionalityReductionBase::calculateStress );
  //std::cout<<" Min engy "<<min_eng<<"\n";
}

}

#endif
