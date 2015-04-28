/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
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
#include "LandmarkSelectionBase.h"
#include "LandmarkRegister.h"
#include "tools/Random.h"
#include <iostream>
#include <math.h>

#define ESP (1e-9)

//#define gamma 0.5
namespace PLMD {
namespace analysis {

class StagedSampling : public LandmarkSelectionBase {
private:
  unsigned seed;
  double gamma;
public:
  StagedSampling( const LandmarkSelectionOptions& lo );
  void select( MultiReferenceBase* );
};

PLUMED_REGISTER_LANDMARKS(StagedSampling,"STAGED")

StagedSampling::StagedSampling( const LandmarkSelectionOptions& lo ):
LandmarkSelectionBase(lo)
{
  parse("SEED",seed);
  parse("GAMMA",gamma);
}

void StagedSampling::select( MultiReferenceBase* myframes ){
  unsigned int n = getNumberOfLandmarks(); //Only one call to such functions.
  unsigned int N = getNumberOfFrames();
  unsigned int m = (int)sqrt(n*N);   // this should be the default but perhaps we should have and option
  std::vector<unsigned> fpslandmarks(m);
  for (unsigned i=0;i<m;++i) fpslandmarks[i]=0;
  // Select first point at random
  Random random; random.setSeed(-seed); double rand=random.RandU01();
  fpslandmarks[0] = std::floor( N*rand );
  //using FPS we want to find m landmarks where m = sqrt(nN)
 
  double sum_weight=0.0;
  bool flag=false;
  std::vector<double> mdlist(N); 
  std::vector<double> wt_frames(N); 
  for(unsigned i=0;i<N;i++){ 
	  mdlist[i] = getDistanceBetweenFrames(fpslandmarks[0],i);
	  wt_frames[i] = getWeightOfFrame(i);
	  sum_weight=sum_weight+wt_frames[i];
  }
  
  if(sum_weight - 0.0 > ESP) flag = true;
  
  unsigned maxj;
  double maxd;
  for(unsigned i=1;i<m;i++){
	  maxd = 0.0;
	  for(unsigned j=0;j<N;j++) if(mdlist[j] > maxd) {maxd = mdlist[j]; maxj = j;}
	  fpslandmarks[i] = maxj;
	  for(unsigned j=0;j<N;j++){
		  double dij = getDistanceBetweenFrames(maxj,j);
		  if(mdlist[j] > dij) mdlist[j] = dij;
	  }
  }
 
 
  for(unsigned int i=0;i<N;i++) wt_frames[i] = pow(wt_frames[i],gamma);
   
 // std::cout << "after FPS"<<std::endl;
  std::vector<std::vector<int> > lneighbours(m);
  std::vector<double> weights(m);  
  for(unsigned i=0;i<m;i++) {
	  weights[i] = (flag == false) ? 1 : (wt_frames[fpslandmarks[i]]);   //!todo: probably frames can have weights, so these should be included
//	  weights[i] = pow(weights[i],gamma);
	  lneighbours[i].push_back(fpslandmarks[i]); // Each element has itself atleast in the neighbourhood.
  }
  int mind_index=0;
  //Now after selecting m landmarks we calculate vornoiweights.
  for(unsigned i=0;i<N;i++){
      double mind_vor=getDistanceBetweenFrames(fpslandmarks[0],i);
	  for(unsigned j=1;j<m;j++){
		if (i != fpslandmarks[j]){
		  double tempd = getDistanceBetweenFrames(fpslandmarks[j],i);
		  if(tempd < mind_vor){
			  mind_vor = tempd;
			  mind_index = j;
		        }
                 }
	  }
	  weights[mind_index] += (flag==false) ? 1 : (wt_frames[fpslandmarks[mind_index]]);
	  lneighbours[mind_index].push_back(i);
  }
  //for(int i=0;i<m;i++) std::cout<<fpslandmarks[i]<< " ";
  //Calulate cumulative weights.
  std::vector<double> cum_weights(m);
  cum_weights[0] = weights[0];
  for(unsigned i=1;i<m;i++){
	  cum_weights[i] = cum_weights[i-1] + weights[i];
  }
  
  //Calculate unique n random sampling from this .
  
  std::vector<bool> selected(m);
  for (unsigned i=0;i<m;i++) selected[i]=false;
  unsigned ncount=0;
  while ( ncount<n){
//  std::cout << "here"<<std::endl;
// generate random weight and check which point it belongs to. select only it was not selected before
      double rand=random.RandU01();
	  int rand_ind = std::floor( cum_weights[m-1]*rand );
	  for(int j=m-2;j>=0;j--){
		  if(rand_ind - cum_weights[j] > 0 && !selected[j+1] ) {
			  unsigned k = lneighbours[j+1].size();
			  
			  double rand=random.RandU01();
			  unsigned in = std::floor(k*rand);
			  unsigned isel = lneighbours[j+1][in];
			  //std::cout<<" isel is"<<isel<<" ";
			  selectFrame(isel,myframes);
              selected[j+1]=true;
              ncount++;
		      break;	
		  }
	  }
  } 
}
}
}
