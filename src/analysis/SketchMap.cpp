/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "DimensionalityReductionBase.h"
#include "ClassicalScaling.h"
#include "SMACOF.h"
#include "tools/ConjugateGradient.h"
#include "tools/GridSearch.h"
#include "reference/PointWiseMapping.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"
#include <iostream>
#include <algorithm> 
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
//+PLUMEDOC ANALYSIS SKETCHMAP
/*
Perform a dimensionality reduction using the sketch-map algorithm

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class SketchMap : public DimensionalityReductionBase {
private:
  bool nosmacof,doglobal;
  double smactol, smaptol, regulariser,ncgrid,nfgrid;
  SwitchingFunction lowdf, highdf;
  unsigned niter;
  double recalculateWeights( const Matrix<double>& Distances, const Matrix<double>& F, PointWiseMapping* mymap, Matrix<double>& Weights );
public:
  static void registerKeywords( Keywords& keys );
  SketchMap( const ActionOptions& ao );
  std::string getAlgorithmName() const { return "sketch-map"; }
  void calculateAllDistances( PointWiseMapping* mymap, Matrix<double>& targets );
  void generateProjections( PointWiseMapping* mymap );
  double transformHD( const double& val, double& df ) const ;
  double transformLD( const double& val, double& df ) const ;
};

PLUMED_REGISTER_ACTION(SketchMap,"SKETCHMAP")

void SketchMap::registerKeywords( Keywords& keys ){
  DimensionalityReductionBase::registerKeywords( keys );
  keys.add("compulsory","HIGH_DIM_FUNCTION","the parameters of the switching function in the high dimensional space");
  keys.add("compulsory","LOW_DIM_FUNCTION","the parameters of the switching function in the low dimensional space");
  keys.add("compulsory","SMACOF_TOL","1E-4","the tolerance for each SMACOF cycle");
  keys.add("compulsory","SMAP_TOL","1E-4","the tolerance for sketch-map");
  keys.add("compulsory","REGULARISE_PARAM","0.001","this is used to ensure that we don't divide by zero when updating weights");
 // keys.add("compulsory","DOGLOBAL","false","logical. to do it or not");
  keys.add("compulsory","CGrid","10","Coarse Grid dimensions ");
  keys.add("compulsory","FGrid","100","Fine Grid dimensions");
  keys.add("compulsory","NITER","1","Number of times to repeat Grid Search");
}

SketchMap::SketchMap( const ActionOptions& ao ):
Action(ao),
DimensionalityReductionBase(ao)
{
  // Read in the switching functions
  std::string linput,hinput, errors;
  parse("HIGH_DIM_FUNCTION",hinput);
  highdf.set(hinput,errors);
  if(errors.length()>0) error(errors);
  parse("LOW_DIM_FUNCTION",linput);
  lowdf.set(hinput,errors);
  if(errors.length()>0) error(errors);
//  parse("DOGLOBAL",doglobal);
  // Read tolerances
  parse("SMACOF_TOL",smactol);
  parse("SMAP_TOL",smaptol);
  parse("REGULARISE_PARAM",regulariser);
  parse("CGrid",ncgrid);
  parse("FGrid",nfgrid);
  parse("NITER",niter);
}

void SketchMap::calculateAllDistances( PointWiseMapping* mymap, Matrix<double>& targets ){
  // Calculate matrix of dissimilarities (High dimensional space) 
  mymap->calculateAllDistances( getPbc(), getArguments(), comm, mymap->modifyDmat(), false );
  double dr; unsigned M = mymap->getNumberOfReferenceFrames(); 
  for(unsigned i=1; i<M; ++i){
      for(unsigned j=0; j<i; ++j) {
          targets(i,j) = targets(j,i) = 1.0 - highdf.calculate( mymap->modifyDmat()(i,j), dr ); // high dim space
      }
  }
}


void SketchMap::generateProjections( PointWiseMapping* mymap ){
	Matrix<double> Distances( mymap->modifyDmat() ); //making a copy
	// Calculates the first guess of projections in LD space 
	ClassicalScaling::run( mymap );

	// Calculate the value of sigma and the weights
	unsigned M = mymap->getNumberOfReferenceFrames();
	Matrix<double> Weights(M,M); double filt = recalculateWeights( Distances, getTargets(), mymap, Weights );

	unsigned MAXSTEPS=100; double newsig;
	for(unsigned i=0;i<MAXSTEPS;++i){
	  // Run the smacof algorithm
	  SMACOF::run( Weights, mymap, smactol );
	  // Recalculate weights matrix and sigma
	  newsig = recalculateWeights( Distances, getTargets(), mymap, Weights );
	  printf("HELLO GARETH AND RACHEL %d %f %f %f \n",i, newsig, filt, fabs( newsig - filt ) );
	  // Test whether or not the algorithm has converged
	  if( fabs( newsig - filt )<smaptol ) break;
	  // Make initial sigma into new sigma so that the value of new sigma is used every time so that the error can be reduced
	  filt=newsig;
	}

	//targets matrix contains the distances of each frame with other frames
	Matrix<double> targets(mymap->modifyDmat());
	targets = getTargets();

	std::vector<double> smacof_error(M);
	std::vector<double> cg_error(M);
	std::vector<double> grid_error(M);

	double smacof_totalerror=0;

	//Calculates error just after Smacof
	for(unsigned i=0;i<M;i++){
	   std::vector<double> p(mymap->getNumberOfProperties());
	   std::vector<double> der(mymap->getNumberOfProperties());
	   //Modify fframes for calculateStress routine	   
	   setTargetVectorForPointwiseGlobalMinimisation(i,targets);
	   for(unsigned k=0;k<mymap->getNumberOfProperties();k++) p[k] = mymap->getProjectionCoordinate(i,k);
	   double error = calculateStress(p,der);
	   smacof_error[i] = error;
	   smacof_totalerror+=error;
	}

	std::ofstream myfile;
	/* Routine to write Projected coordinates after Smacof */
	myfile.open("Smacof.dat");
	for(unsigned i=0;i<M;i++){
		for(unsigned j=0;j<mymap->getNumberOfProperties();j++){
			myfile<<mymap->getProjectionCoordinate(i,j)<<" ";
		}
		myfile<<smacof_error[i]/(M*(M-1));
	myfile<<"\n";
	}
	myfile.close();
	std::cout<<"Total error after Smacof "<< smacof_totalerror/(M*(M-1)) <<"\n";	
/////////////////////////////////////////////////////////////////////////////   

	double cgtol = 1E-6;
	double totalerror=0;
	   for(unsigned i=0;i<M;i++){
		   std::vector<double> p(mymap->getNumberOfProperties());
		   std::vector<double> der(mymap->getNumberOfProperties());
	           setTargetVectorForPointwiseGlobalMinimisation(i,targets);
		   //Copy frame coordinates to p.
		   for(unsigned j=0;j<mymap->getNumberOfProperties();j++) p[j] = mymap->getProjectionCoordinate(i,j); 
		   ConjugateGradient<DimensionalityReductionBase> myminimiser2( this );
		   myminimiser2.minimise( cgtol, p, &DimensionalityReductionBase::calculateStress );
		   for(unsigned j=0;j<p.size();++j) mymap->setProjectionCoordinate(i,j,p[j]); //Set the projected point in the mymap object.
	           double error = calculateStress(p,der);
                   totalerror+=error;
            }
	std::cout<<"Total error after Smacof+CG "<<totalerror/(M*(M-1)) <<"\n";	
	/* Routine to write Projected coordinates after Smacof */
	myfile.open("Smacof+CG.dat");
	for(unsigned i=0;i<M;i++){
		for(unsigned j=0;j<mymap->getNumberOfProperties();j++){
			myfile<<mymap->getProjectionCoordinate(i,j)<<" ";
		}
	myfile<<"\n";
	}
	myfile.close();
//Routine for grid search.
   
	std::cout<<"Starting Pointwise Global Optimization-->";	
        std::cout<<" Coarse Grid "<< ncgrid<<" X " << ncgrid ;
        std::cout<<" --> Fine Grid "<< nfgrid<<" X " << nfgrid <<"\n";
	for(int cnt=0;cnt<niter;cnt++){
		double minx,miny,maxx,maxy;
		maxy = miny = mymap->getProjectionCoordinate(0,1);
		maxx = minx = mymap->getProjectionCoordinate(0,0);
	  
		for(unsigned i=1;i<mymap->getNumberOfReferenceFrames();i++){
			if(mymap->getProjectionCoordinate(i,0) < minx) minx = mymap->getProjectionCoordinate(i,0);
			if(mymap->getProjectionCoordinate(i,1) < miny) miny = mymap->getProjectionCoordinate(i,1);
			if(mymap->getProjectionCoordinate(i,0) > maxx) maxx = mymap->getProjectionCoordinate(i,0);
			if(mymap->getProjectionCoordinate(i,1) > maxy) maxy = mymap->getProjectionCoordinate(i,1);
		} 
// padding for the grid 
		minx=minx-0.1*minx;
		miny=miny-0.1*miny;
		maxx=maxx+0.1*minx;
		maxy=maxy+0.1*miny;
		double stepx = fabs(maxx-minx)/(double)ncgrid;
		double stepy = fabs(maxy-miny)/(double)ncgrid;
       //Run Grid Search for each frame and get the correct projection.
	   for(unsigned i=0;i<M;i++){
		   std::vector<double> p(mymap->getNumberOfProperties());
		   std::vector<double> der(mymap->getNumberOfProperties());
	           setTargetVectorForPointwiseGlobalMinimisation(i,targets);
		   //Copy frame coordinates to p.
		   for(unsigned j=0;j<mymap->getNumberOfProperties();j++) p[j] = mymap->getProjectionCoordinate(i,j); 
		   std::vector<double> temp(p.size());
		   std::vector<double> grid_point(p.size());
		   double min_eng,eng_pt,curr_x,curr_y;
		   //To hold stress at each grid point.
		   std::vector< std::vector<double> > stressgrid(ncgrid,std::vector<double>(ncgrid));
		   //To hold x and y coordinates of grid points.
		   std::vector<double> ptsinx;
		   std::vector<double> ptsiny;	
		   grid_point[0] = curr_x = minx;
		   grid_point[1] = curr_y = miny;		   
		   //Calculate stress at the frame point.
		   min_eng = calculateStress(p,der) ;		
       //            std::cout<<"min_eng= "<<min_eng<<" "<<i<<"\n"; 
       //            std::cout<<"min_eng= "<<min_eng<<" "<<i<<"\n"; 
	//           std::cout<<"\n Optimizing frame: " << i <<" Coarse Grid " ;	   
		   for(int l=0;l<ncgrid;l++){
				curr_x = minx + stepx*l;
				ptsinx.push_back(curr_x);
				grid_point[0] = curr_x;
				for(int j=0;j<ncgrid;j++){					
					curr_y = miny + stepy*j;
					if(l==0) ptsiny.push_back(curr_y);
					grid_point[1] = curr_y;
					eng_pt = calculateStress(grid_point,der);
					stressgrid[l][j] = eng_pt; //Store stress in the stress grid.
					if(eng_pt < min_eng){
	  //                                      std::cout<<"FOUND " ;	   
						min_eng = eng_pt;			
						p[0] = curr_x;p[1] = curr_y;	//Update the point if stress is less here.	
					}		  		  
				}
		   }
		   GridSearch<DimensionalityReductionBase> myminimiser( this );
            //       double temp_en=min_eng;    
		   myminimiser.minimise(p,&DimensionalityReductionBase::calculateStress,stressgrid,nfgrid,min_eng,ptsinx,ptsiny);
              //     if(min_eng < temp_en) std::cout<<"Fine Grid: Found" ;	
		   grid_error[i] = calculateStress(p,der);
		   ConjugateGradient<DimensionalityReductionBase> myminimiser2( this );
		   myminimiser2.minimise( cgtol, p, &DimensionalityReductionBase::calculateStress );
		   for(unsigned j=0;j<p.size();++j) mymap->setProjectionCoordinate(i,j,p[j]); //Set the projected point in the mymap object.

	   }
	   
/*	    myfile.open("gridsearch.dat");
		for(unsigned i=0;i<M;i++){
			for(unsigned j=0;j<mymap->getNumberOfProperties();j++){
			myfile<<mymap->getProjectionCoordinate(i,j)<<" ";
			}
			myfile<<grid_error[i]<<"\n";
		}
		myfile.close();		  
*/	   
		double grid_totalerror = 0;
		for(unsigned q=0;q<M;q++){
			std::vector<double> p(mymap->getNumberOfProperties());
			std::vector<double> deri(mymap->getNumberOfProperties());	   
			setTargetVectorForPointwiseGlobalMinimisation(q,targets);
			for(unsigned k=0;k<mymap->getNumberOfProperties();k++) p[k] = mymap->getProjectionCoordinate(q,k);
			double error = calculateStress(p,deri);
			grid_error[q] = error;
			grid_totalerror+=grid_error[q];
		}
	 
	 std::cout<<"Total error after Pointwise Global Optimization after iteration number "<<cnt+1<<" : " <<grid_totalerror/(M*(M-1)) <<"\n";
	 std::ostringstream fn;
	 fn << "global_" << cnt << ".dat";
	 std::ofstream out(fn.str().c_str(),std::ios_base::binary);
	 for(unsigned i=0;i<M;i++){
		for(unsigned j=0;j<mymap->getNumberOfProperties();j++){
			out<<mymap->getProjectionCoordinate(i,j)<<" ";
		}
		out<<grid_error[i]/(M*(M-1));
		out<<"\n";
	 }
	 out.close();

   }
		  //~ 


}

double SketchMap::recalculateWeights( const Matrix<double>& Distances, const Matrix<double>& F, PointWiseMapping* mymap, Matrix<double>& Weights ){
  double filt=0, totalWeight=0.;; double dr;
  for(unsigned i=1; i<Weights.nrows(); ++i){
      for(unsigned j=0; j<i; ++j){
          double tempd=0;
          for(unsigned k=0;k<mymap->getNumberOfProperties();++k){
             double tmp = mymap->getProjectionCoordinate( i, k ) - mymap->getProjectionCoordinate( j, k );
             tempd += tmp*tmp;
          }
          double ninj=mymap->getWeight(i)*mymap->getWeight(j);
          totalWeight += ninj;

          double dij=sqrt(tempd);
          double fij = 1.0 - lowdf.calculate( dij, dr );
          double filter=F(i,j)-fij;
          double diff=Distances(i,j) - dij;
          if( fabs(diff)<regulariser ) Weights(i,j)=Weights(j,i)=0.0;
          else Weights(i,j)=Weights(j,i) = ( -ninj*filter*dij*dr ) / diff;
          filt += ninj*filter*filter;
      }
  }
  return filt / totalWeight;
}

double SketchMap::transformHD( const double& val, double& df ) const { 
  double vv=1.0 - highdf.calculate( val, df );
  df=-df; return vv;
}
double SketchMap::transformLD( const double& val, double& df ) const {
  double vv=1.0 - lowdf.calculate( val, df );
  df=-df; return vv;
}

}
}
