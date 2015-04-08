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
namespace PLMD {
namespace analysis {

class StagedSampling : public LandmarkSelectionBase {
private:
  unsigned seed;
public:
  StagedSampling( const LandmarkSelectionOptions& lo );
  void select( MultiReferenceBase* );
};

PLUMED_REGISTER_LANDMARKS(StagedSampling,"STAGED")

StagedSampling::StagedSampling( const LandmarkSelectionOptions& lo ):
LandmarkSelectionBase(lo)
{
  parse("SEED",seed);
}

void StagedSampling::select( MultiReferenceBase* myframes ){
  unsigned int n = getNumberOfLandmarks(); //Only one call to such functions.
  unsigned int N = getNumberOfFrames();
  unsigned int m = (int)sqrt(n*N);
  std::vector<unsigned> fpslandmarks(m);
  // Select first point at random
  Random random; random.setSeed(-seed); double rand=random.RandU01();
  fpslandmarks[0] = std::floor( N*rand );

  //using FPS we want to find m landmarks where m = sqrt(nN)
  // Now find distance to all other points
  Matrix<double> distances( m, N );
  for(unsigned int i=0;i<N;++i) {
         distances(0,i) = getDistanceBetweenFrames( fpslandmarks[0], i );
}
  // Now find all other landmarks
  for(unsigned i=1;i<m;++i){
      // Find point that has the largest minimum distance from the landmarks selected thus far
      double maxd=0;
      for(unsigned j=0;j<N;++j){
          double mind=distances(0,j);
          for(unsigned k=1;k<i;++k){
              if( distances(k,j)<mind ){ mind=distances(k,j); }
          }
          if( mind>maxd ){ maxd=mind; fpslandmarks[i]=j; }
      }
      for(unsigned k=0;k<getNumberOfFrames();++k) distances(i,k) = getDistanceBetweenFrames( fpslandmarks[i], k );
         
  }
   
  int weights[m];
  for(int i=0;i<m;i++) weights[i] = 1;
  double mind_vor=999999999999;
  int mind_index=0;
  //Now after selecting m landmarks we calculate vornoiweights.
  for(int i=0;i<N;i++){
	  for(int j=0;j<m;j++){
                  if (i != j) {
		  double tempd = getDistanceBetweenFrames(i,j);
		  if(tempd < mind_vor){
			  mind_vor = tempd;
			  mind_index = j;
		   }
                   }
	  }
	  weights[mind_index]++;
	  mind_vor = 999999999999;
  }
  //Calulate cumulative weights.
  int cum_weights[m];
  cum_weights[0] = weights[0];
  for(int i=1;i<m;i++){
	  cum_weights[i] = cum_weights[i-1] + weights[i];
  }
  //Calculate unique n random sampling from this .
  bool selected[m];
  for (int i=0;i<m;i++) selected[i]=false;
  int ncount=0;
  while ( ncount<n){
  int flag =0;
// generate random weight and check which point it belongs to. select only it was not selected before
          double rand=random.RandU01();
	  int rand_ind = std::floor( cum_weights[m-1]*rand );
	  for(int j=m-2;j>=0;j--){
		  if(rand_ind - cum_weights[j] > 0 && !selected[j+1] ) {
			  selectFrame(j+1,myframes);
 //                         std::cout<<"selecting"<<j+1<<" frame"<<std::endl;
                          selected[j+1]=true;
                          ncount++;
                          flag=1;
			  break;
		  }
//  in case our random weight falls in the border of two blocks

		  if(rand_ind - cum_weights[j] == 0 && !selected[j]) {
			  selectFrame(j,myframes);
                          selected[j]=true;
                          ncount++;
                          flag=1;
			  break;
		  }
	  }
// if our random weight belongs to first sample
         if(flag==0){
			  selectFrame(0,myframes);
                          selected[0]=true;
                          ncount++;
                     }
                       
         
  } 
  
	  
		   
}

}
}
