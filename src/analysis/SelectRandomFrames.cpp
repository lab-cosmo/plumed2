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
#include <fstream>
namespace PLMD {
namespace analysis {

class SelectRandomFrames : public LandmarkSelectionBase {
private:
  unsigned seed;
public:
  SelectRandomFrames( const LandmarkSelectionOptions& lo );
  void select( MultiReferenceBase* );
};

PLUMED_REGISTER_LANDMARKS(SelectRandomFrames,"RANDOM")

SelectRandomFrames::SelectRandomFrames( const LandmarkSelectionOptions& lo ):
LandmarkSelectionBase(lo)
{
parse("SEED",seed);
}

void SelectRandomFrames::select( MultiReferenceBase* myframes ){
//  std::ofstream outfile;
//  outfile.open ("random1.dat");
  Random random; random.setSeed(-seed);
  unsigned nframe= getNumberOfFrames();
  unsigned nlandmarks= getNumberOfLandmarks();
  std::vector<bool>selected(nframe);
  for (unsigned i=0;i<nframe;++i) selected[i]=false;
 
  
//  outfile << "selecting  " << nlandmarks << "From " << getNumberOfFrames() << std::endl ;
  unsigned fcount=0;
  while (fcount<nlandmarks) 
   {
   double rand=random.RandU01();
   unsigned iframe= std::floor( nframe*rand );
   if (!selected[iframe])
     {
 //    outfile << iframe <<"  "<< fcount << std::endl; 
     selected[iframe]=true;
     selectFrame( iframe, myframes );
     ++fcount;
     } 
  }
// outfile.close();
}

}
}
