/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
#include "AnalysisWithLandmarks.h"
#include "ClassicalScaling.h"
#include "SMACOF.h"
#include "reference/PointWiseMapping.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace analysis {

//+PLUMEDOC ANALYSIS LANDMARK_SELECT
/*
Just select landmarks and outputs them to a file.


\par Examples

The following command instructs plumed to construct a classical multidimensional scaling projection of a trajectory.
The RMSD distance between atoms 1-256 have moved is used to measure the distances in the high-dimensional space.

\verbatim
LANDMARK_SELECT ...
  ATOMS=1-256
  METRIC=OPTIMAL-FAST 
  USE_ALL_DATA
  OUTPUT_FILE=rmsd-embed
... LANDMARK_SELECT
\endverbatim
*/
//+ENDPLUMEDOC

class LandmarkSelector : public AnalysisWithLandmarks {
private:  
  MultiReferenceBase *mydata;
  std::string ofilename;
public:
  static void registerKeywords( Keywords& keys );
  LandmarkSelector( const ActionOptions& ao );
  ~LandmarkSelector();
  void analyzeLandmarks();
};

PLUMED_REGISTER_ACTION(LandmarkSelector,"LANDMARK_SELECT")

void LandmarkSelector::registerKeywords( Keywords& keys ){
  AnalysisWithLandmarks::registerKeywords( keys );
  keys.add("compulsory","OUTPUT_FILE","file on which to output the final embedding coordinates");  
}

LandmarkSelector::LandmarkSelector( const ActionOptions& ao ):
Action(ao),
AnalysisWithLandmarks(ao)
{  
  mydata = new MultiReferenceBase( getMetricName(), false );
  setDataToAnalyze( mydata );
  parseOutputFile("OUTPUT_FILE",ofilename);
}

LandmarkSelector::~LandmarkSelector(){
}

void LandmarkSelector::analyzeLandmarks(){
  OFile gfile; gfile.link(*this); 
  gfile.setBackupString("analysis");
  gfile.fmtField(getOutputFormat()+" ");
  gfile.open( ofilename.c_str() );
  gfile <<this->getNumberOfLandmarks()<< " hello\n";
  unsigned M=mydata->getNumberOfReferenceFrames(); 
  for(unsigned i=0; i<M; ++i){
     mydata->getFrame(i)->print(gfile, "random string");
     }  
  gfile.close();
}

}
}
