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
  unsigned int m = (int)sqrt(n*N);   // this should be the default but perhaps we should have and option
  std::vector<unsigned> fpslandmarks(m);
  for (unsigned i=0;i<m;++i) fpslandmarks[i]=0;
  // Select first point at random
  Random random; random.setSeed(-seed); double rand=random.RandU01();
  fpslandmarks[0] = std::floor( N*rand );
  //using FPS we want to find m landmarks where m = sqrt(nN)
 
 //Michele's method for finding landmarks.
 
     //~ else if (smode=="minmax") 
    //~ {
        //~ // farthest point sampling selection of the points. no risk on getting duplicates here
        //~ 
        //~ for (unsigned long i=1; i<n; ++i)
        //~ {
//~ 
            //~ maxd=0.;  for (unsigned long j=0; j<N; ++j) if (mdlist[j]>maxd) {maxd=mdlist[j]; maxj=j;}
            //~ std::cerr<<"selecting point "<<i<<" : "<<maxj<<"("<<maxd<<")\n";
//~ 
            //~ isel[i]=maxj;
            //~ LP.row(i)=HP.row(maxj);
            //~ for (unsigned long j=0; j<N; ++j) 
            //~ { dij=metric->dist(&LP(i,0),&HP(j,0),D); if (mdlist[j]>dij) mdlist[j]=dij; }
        //~ }
    //~ }
 
  std::vector<double> mdlist(N);
//  m=30;  
  for(unsigned i=0;i<N;i++){ mdlist[i] = getDistanceBetweenFrames(fpslandmarks[0],i);}
  unsigned maxj;
  double maxd;
  for(unsigned i=1;i<m;i++){
	  maxd = 0.0;
	  for(unsigned j=0;j<N;j++) if(mdlist[j] > maxd) {maxd = mdlist[j]; maxj = j;}
//	  std::cout<<"maxj ="<<maxj<< " maxd="<<maxd<<std::endl;
	  fpslandmarks[i] = maxj;
	  for(unsigned j=0;j<N;j++){
		  double dij = getDistanceBetweenFrames(maxj,j);
		  if(mdlist[j] > dij) mdlist[j] = dij;
	  }
  }
 
std::cout<<"out"<<std::endl; 
 
 
 
 
 
  //~ // Now find all other landmarks
    //~ Matrix<double> distances( m, N );
  //~ for(unsigned i=0;i<N;++i) distances(0,i) = getDistanceBetweenFrames( fpslandmarks[0], i );
//~ 
  //~ // Now find all other landmarks
  //~ for(unsigned i=1;i<m;++i){
      //~ // Find point that has the largest minimum distance from the landmarks selected thus far
      //~ double maxd=0;
      //~ for(unsigned j=0;j<N;++j){
          //~ double mind=distances(0,j);
          //~ for(unsigned k=1;k<i;++k){
              //~ if( distances(k,j)<mind ){ mind=distances(k,j); }
          //~ }
          //~ if( mind>maxd ){ maxd=mind; fpslandmarks[i]=j; }
      //~ }
      //~ //selectFrame( landmarks[i], myframes );
      //~ for(unsigned k=0;k<getNumberOfFrames();++k) distances(i,k) = getDistanceBetweenFrames( fpslandmarks[i], k );
  //~ }
  
 
   
 // std::cout << "after FPS"<<std::endl;
  std::vector<std::vector<int> > lneighbours(m);
  std::vector<int> weights(m);  
  for(unsigned i=0;i<m;i++) {
	  weights[i] = 1;   //!todo: probably frames can have weights, so these should be included
	  lneighbours[i].push_back(i); // Each element has itself atleast in the neighbourhood.
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
	  weights[mind_index]++;
	  lneighbours[mind_index].push_back(i);
  }
  
  //Calulate cumulative weights.
  std::vector<int> cum_weights(m);
  cum_weights[0] = weights[0];
  for(unsigned i=1;i<m;i++){
	  cum_weights[i] = cum_weights[i-1] + weights[i];
  }
  
std::cout<<"out 2"<<std::endl; 
  //Calculate unique n random sampling from this .
  std::vector<bool> selected(N);
  for (unsigned i=0;i<m;i++) selected[i]=false;
  unsigned ncount=0;
  while ( ncount<n){
//  std::cout << "here"<<std::endl;
// generate random weight and check which point it belongs to. select only it was not selected before
      double rand=random.RandU01();
	  int rand_ind = std::floor( cum_weights[m-1]*rand );
	  for(int j=m-2;j>=0;j--){
		  if(rand_ind - cum_weights[j] > 0 && !selected[j+1] ) {
			  int k = lneighbours[fpslandmarks[j+1]].size();
			  double rand=random.RandU01();
			  int index_of_selection = std::floor(k*rand);
			  selectFrame(index_of_selection,myframes);
              selected[index_of_selection]=true;
              ncount++;
		      break;	
           //!todo: reaally select one random point from the selected voronoi polyhedron
		  }
	  }
  } 
}
}
}
