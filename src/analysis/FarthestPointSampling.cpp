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

namespace PLMD {
namespace analysis {

class FarthestPointSampling : public LandmarkSelectionBase {
private:
  unsigned seed;
public:
  FarthestPointSampling( const LandmarkSelectionOptions& lo );
  void select( MultiReferenceBase* );
};

PLUMED_REGISTER_LANDMARKS(FarthestPointSampling,"FPS")

FarthestPointSampling::FarthestPointSampling( const LandmarkSelectionOptions& lo ):
LandmarkSelectionBase(lo)
{
  parse("SEED",seed);
}

void FarthestPointSampling::select( MultiReferenceBase* myframes ){
  

  //~ fpslandmarks[0] = std::floor( N*rand );
  //~ //using FPS we want to find m landmarks where m = sqrt(nN)
 //~ 
//~ 
  //~ std::vector<double> mdlist(N); 
  //~ for(unsigned i=0;i<N;i++){ mdlist[i] = getDistanceBetweenFrames(fpslandmarks[0],i);}
  //~ unsigned maxj;
  //~ double maxd;
  //~ for(unsigned i=1;i<m;i++){
	  //~ maxd = 0.0;
	  //~ for(unsigned j=0;j<N;j++) if(mdlist[j] > maxd) {maxd = mdlist[j]; maxj = j;}
	  //~ fpslandmarks[i] = maxj;
	  //~ for(unsigned j=0;j<N;j++){
		  //~ double dij = getDistanceBetweenFrames(maxj,j);
		  //~ if(mdlist[j] > dij) mdlist[j] = dij;
	  //~ }
  //~ }
  
  unsigned N = getNumberOfFrames();
  unsigned n = getNumberOfLandmarks();
  
  std::vector<unsigned> landmarks(n);
  Random random; random.setSeed(-seed); double rand=random.RandU01();
  landmarks[0] = std::floor( N*rand );
  selectFrame( landmarks[0], myframes );
  

  std::vector<double>mdlist(N);
  for(unsigned i=0;i<N;i++){ mdlist[i] = getDistanceBetweenFrames(landmarks[0],i);}
  unsigned maxj;
  double maxd;
  for(unsigned i=1;i<n;i++){
	  maxd = 0.0;
	  for(unsigned j=0;j<N;j++) if(mdlist[j] > maxd) {maxd = mdlist[j]; maxj = j;}
	  landmarks[i] = maxj;
	  selectFrame(landmarks[i],myframes);
	  for(unsigned j=0;j<N;j++){
		  double dij = getDistanceBetweenFrames(maxj,j);
		  if(mdlist[j] > dij) mdlist[j] = dij;
	  }
  }



  //~ // Select first point at random
  //~ Random random; random.setSeed(-seed); double rand=random.RandU01();
  //~ landmarks[0] = std::floor( getNumberOfFrames()*rand );
  //~ selectFrame( landmarks[0], myframes );
//~ 
  //~ // Now find distance to all other points
  //~ Matrix<double> distances( getNumberOfLandmarks(), getNumberOfFrames() );
  //~ for(unsigned i=0;i<getNumberOfFrames();++i) distances(0,i) = getDistanceBetweenFrames( landmarks[0], i );
//~ 
  //~ // Now find all other landmarks
  //~ for(unsigned i=1;i<getNumberOfLandmarks();++i){
      //~ // Find point that has the largest minimum distance from the landmarks selected thus far
      //~ double maxd=0;
      //~ for(unsigned j=0;j<getNumberOfFrames();++j){
          //~ double mind=distances(0,j);
          //~ for(unsigned k=1;k<i;++k){
              //~ if( distances(k,j)<mind ){ mind=distances(k,j); }
          //~ }
          //~ if( mind>maxd ){ maxd=mind; landmarks[i]=j; }
      //~ }
      //~ selectFrame( landmarks[i], myframes );
      //~ for(unsigned k=0;k<getNumberOfFrames();++k) distances(i,k) = getDistanceBetweenFrames( landmarks[i], k );
  //~ } 
}

}
}
