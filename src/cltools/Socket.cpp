/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
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
#include "CLTool.h"
#include "CLToolRegister.h"
#include "tools/Tools.h"
#include "wrapper/Plumed.h"
#include "tools/Communicator.h"
#include "tools/Random.h"
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
#include "tools/Units.h"



#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <netdb.h>


using namespace std;

namespace PLMD {
namespace cltools{

//+PLUMEDOC TOOLS socket
/*
Runs PLUMED in driver mode fetching atomic configurations from a socket.
Data exchange is implemented based on the i-PI protocol. 

\par Examples

\verbatim
plumed socket --plumed plumed.dat  --host host [ --port port  | --unix ]
\endverbatim


*/
//+ENDPLUMEDOC
//

namespace sockets {
void open(int *psockfd, int* inet, int* port, const char* host)
/* Opens a socket.

Note that fortran passes an extra argument for the string length, but this is
ignored here for C compatibility.

Args:
   psockfd: The id of the socket that will be created.
   inet: An integer that determines whether the socket will be an inet or unix
      domain socket. Gives unix if 0, inet otherwise.
   port: The port number for the socket to be created. Low numbers are often
      reserved for important channels, so use of numbers of 4 or more digits is
      recommended.
   host: The name of the host server.
*/

{
   int sockfd, ai_err;

   if (*inet>0)
   {  // creates an internet socket
      
      // fetches information on the host      
      struct addrinfo hints, *res;  
      char service[256];
   
      memset(&hints, 0, sizeof(hints));
      hints.ai_socktype = SOCK_STREAM;
      hints.ai_family = AF_UNSPEC;
      hints.ai_flags = AI_PASSIVE;

      sprintf(service,"%d",*port); // convert the port number to a string
      ai_err = getaddrinfo(host, service, &hints, &res); 
      if (ai_err!=0) { perror("Error fetching host data. Wrong host name?"); exit(-1); }

      // creates socket
      sockfd = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
      if (sockfd < 0) { perror("Error opening socket"); exit(-1); }
    
      // makes connection
      if (connect(sockfd, res->ai_addr, res->ai_addrlen) < 0) 
      { perror("Error opening INET socket: wrong port or server unreachable"); exit(-1); }
      freeaddrinfo(res);
   }
   else
   {  
      struct sockaddr_un serv_addr;

      // fills up details of the socket addres
      memset(&serv_addr, 0, sizeof(serv_addr));
      serv_addr.sun_family = AF_UNIX;
      strcpy(serv_addr.sun_path, "/tmp/ipi_");
      strcpy(serv_addr.sun_path+9, host);
      // creates a unix socket
  
      // creates the socket
      sockfd = socket(AF_UNIX, SOCK_STREAM, 0);

      // connects
      if (connect(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) 
      { perror("Error opening UNIX socket: path unavailable, or already existing"); exit(-1); }
   }


   *psockfd=sockfd;
  }
  
  void writebuffer(int *psockfd, const char *data, int* plen)
/* Writes to a socket.

Args:
   psockfd: The id of the socket that will be written to.
   data: The data to be written to the socket.
   plen: The length of the data in bytes.
*/

{
   int n;
   int sockfd=*psockfd;
   int len=*plen;

   n = write(sockfd,data,len);
   if (n < 0) { perror("Error writing to socket: server has quit or connection broke"); exit(-1); }
}


void readbuffer(int *psockfd, char *data, int* plen)
/* Reads from a socket.

Args:
   psockfd: The id of the socket that will be read from.
   data: The storage array for data read from the socket.
   plen: The length of the data in bytes.
*/

{
   int n, nr;
   int sockfd=*psockfd;
   int len=*plen;

   n = nr = read(sockfd,data,len);

   while (nr>0 && n<len )
   {  nr=read(sockfd,&data[n],len-n); n+=nr; }

   if (n == 0) { perror("Error reading from socket: server has quit or connection broke"); exit(-1); }
}
}


template<typename real>
class Socket : public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  explicit Socket(const CLToolOptions& co );
  int main(FILE* in, FILE*out, Communicator& pc);
  string description() const;
};

template<typename real>
void Socket<real>::registerKeywords( Keywords& keys ){
  CLTool::registerKeywords( keys ); keys.isDriver();
  keys.add("compulsory","--plumed","plumed.dat","specify the name of the plumed input file");
  keys.add("compulsory","--host","localhost","the name/IP address/UNIX socket for the host");
  keys.add("optional","--port","the port number (for Internet sockets)");
  keys.addFlag("--unix",false,"uses a UNIX domain socket");
}
template<typename real>
Socket<real>::Socket(const CLToolOptions& co ):
CLTool(co)
{
 inputdata=commandline;
}
template<typename real>
string Socket<real>::description()const{ return "runs PLUMED over a socket"; }

template<typename real>
int Socket<real>::main(FILE* in,FILE*out,Communicator& pc){

  // variables needed for storage and communication stuff
  int port; int inet, master, hasdata, sockerr;
  int ipisock, me; double *buffer; long bsize;
// Parse everything  
// Read the plumed input file name  
  string plumedFile; parse("--plumed",plumedFile);
// the stride
  port=0; parse("--port",port);  
  bool unix; parseFlag("--unix",unix); if (unix) inet = 0; else inet = 1;
  if (!unix && port==0) error("You should either use a UNIX domain socket or the port to connect to.");
// the hostname
  string hostname; parse("--host",hostname);
  
  std::cerr << "Setting up internals\n";
  
  Plumed p;
  int rr=sizeof(real);
  p.cmd("setRealPrecision",&rr);
  real bohr2nm = 0.052917721; p.cmd("setMDLengthUnits", &bohr2nm);
  real ha2kjmol = 2625.4996; p.cmd("setMDEnergyUnits", &ha2kjmol);
  string timeunits="ps"; p.cmd("setMDEnergyUnits", &timeunits);
  p.cmd("setMDEngine","socket");
  real timestep = 1.0; p.cmd("setTimestep", &timestep);
  p.cmd("setPlumedDat",plumedFile.c_str());
  p.cmd("setLog",out);
  std::cerr << "setting up the communicator\n";  
  p.cmd("setMPIComm", pc.Get_comm());
  
  int natoms;
  

  std::cerr << "GOT TO OPENIG THE SOCKET\n";
  if (pc.Get_rank()==0) sockets::open(&ipisock, &inet, &port, hostname.c_str());
  else ipisock=0;
  
  std::string line;
  std::vector<real> coordinates;
  std::vector<real> forces;
  std::vector<real> masses;
  std::vector<real> charges;
  std::vector<real> cell;
  std::vector<real> virial;

/*
  while(true){
    if(!noatoms){

    bool first_step=false;
    if(!noatoms){
      if(use_molfile==false && (trajectory_fmt=="xyz" || trajectory_fmt=="gro")){
        if(trajectory_fmt=="gro") if(!Tools::getline(fp,line)) error("premature end of trajectory file");
        sscanf(line.c_str(),"%100d",&natoms);
      }
    }
    if(checknatoms<0 && !noatoms){
      pd_nlocal=natoms;
      pd_start=0;
      first_step=true;
      masses.assign(natoms,real(1.0));
      charges.assign(natoms,real(0.0));
//case pdb: structure
      if(pdbfile.length()>0){
        for(unsigned i=0;i<pdb.size();++i){
          AtomNumber an=pdb.getAtomNumbers()[i];
          unsigned index=an.index();
          if( index>=unsigned(natoms) ) error("atom index in pdb exceeds the number of atoms in trajectory");
          masses[index]=pdb.getOccupancy()[i];
          charges[index]=pdb.getBeta()[i];
        }
      }
      if(mcfile.length()>0){
        IFile ifile;
        ifile.open(mcfile);
        int index; double mass; double charge;
        while(ifile.scanField("index",index).scanField("mass",mass).scanField("charge",charge).scanField()){
          masses[index]=mass;
          charges[index]=charge;
        }
      }
    } else if( checknatoms<0 && noatoms ){ 
      natoms=0; 
    }
    if( checknatoms<0 ){
      checknatoms=natoms;
      p.cmd("setNatoms",&natoms);
      p.cmd("init");
    }
    if(checknatoms!=natoms){
       std::string stepstr; Tools::convert(step,stepstr);
       error("number of atoms in frame " + stepstr + " does not match number of atoms in first frame");
    }

    coordinates.assign(3*natoms,real(0.0));
    forces.assign(3*natoms,real(0.0));
    cell.assign(9,real(0.0));
    virial.assign(9,real(0.0));

    if( first_step || rnd.U01()>0.5){
      if(debug_pd){
        int npe=intracomm.Get_size();
        vector<int> loc(npe,0);
        vector<int> start(npe,0);
        for(int i=0;i<npe-1;i++){
          int cc=(natoms*2*rnd.U01())/npe;
          if(start[i]+cc>natoms) cc=natoms-start[i];
          loc[i]=cc;
          start[i+1]=start[i]+loc[i];
        }
        loc[npe-1]=natoms-start[npe-1];
        intracomm.Bcast(loc,0);
        intracomm.Bcast(start,0);
        pd_nlocal=loc[intracomm.Get_rank()];
        pd_start=start[intracomm.Get_rank()];
        if(intracomm.Get_rank()==0){
          fprintf(out,"\nDRIVER: Reassigning particle decomposition\n");
          fprintf(out,"DRIVER: "); for(int i=0;i<npe;i++) fprintf(out,"%d ",loc[i]); printf("\n");
          fprintf(out,"DRIVER: "); for(int i=0;i<npe;i++) fprintf(out,"%d ",start[i]); printf("\n");
        }
        p.cmd("setAtomsNlocal",&pd_nlocal);
        p.cmd("setAtomsContiguous",&pd_start);
      } else if(debug_dd){
        int npe=intracomm.Get_size();
        int rank=intracomm.Get_rank();
        dd_charges.assign(natoms,0.0);
        dd_masses.assign(natoms,0.0);
        dd_gatindex.assign(natoms,-1);
        dd_g2l.assign(natoms,-1);
        dd_coordinates.assign(3*natoms,0.0);
        dd_forces.assign(3*natoms,0.0);
        dd_nlocal=0;
        for(int i=0;i<natoms;++i){
          double r=rnd.U01()*npe;
          int n; for(n=0;n<npe;n++) if(n+1>r)break;
          plumed_assert(n<npe);
          if(n==rank){
            dd_gatindex[dd_nlocal]=i;
            dd_g2l[i]=dd_nlocal;
            dd_charges[dd_nlocal]=charges[i];
            dd_masses[dd_nlocal]=masses[i];
            dd_nlocal++;
          }
        }
        if(intracomm.Get_rank()==0){
          fprintf(out,"\nDRIVER: Reassigning particle decomposition\n");
        }
        p.cmd("setAtomsNlocal",&dd_nlocal);
        p.cmd("setAtomsGatindex",&dd_gatindex[0]);
      }
    }

    int plumedStopCondition=0;
    p.cmd("setStep",&step);
    p.cmd("setStopFlag",&plumedStopCondition);
    if(!noatoms){
       if(use_molfile){
#ifdef __PLUMED_HAS_MOLFILE
    	   if(pbc_cli_given==false) {
                 if(ts_in.A>0.0){ // this is negative if molfile does not provide box
    		   // info on the cell: convert using pbcset.tcl from pbctools in vmd distribution
    		   real cosBC=cos(ts_in.alpha*pi/180.);
    		   //double sinBC=sin(ts_in.alpha*pi/180.);
    		   real cosAC=cos(ts_in.beta*pi/180.);
    		   real cosAB=cos(ts_in.gamma*pi/180.);
    		   real sinAB=sin(ts_in.gamma*pi/180.);
    		   real Ax=ts_in.A;
    		   real Bx=ts_in.B*cosAB;
    		   real By=ts_in.B*sinAB;
                   real Cx=ts_in.C*cosAC;
                   real Cy=(ts_in.C*ts_in.B*cosBC-Cx*Bx)/By;
                   real Cz=sqrt(ts_in.C*ts_in.C-Cx*Cx-Cy*Cy);
    		   cell[0]=Ax/10.;cell[1]=0.;cell[2]=0.;
    		   cell[3]=Bx/10.;cell[4]=By/10.;cell[5]=0.;
    		   cell[6]=Cx/10.;cell[7]=Cy/10.;cell[8]=Cz/10.;
                 } else {
                   cell[0]=0.0; cell[1]=0.0; cell[2]=0.0;
                   cell[3]=0.0; cell[4]=0.0; cell[5]=0.0;
                   cell[6]=0.0; cell[7]=0.0; cell[8]=0.0;
                 }
    	   }else{
    		   for(unsigned i=0;i<9;i++)cell[i]=pbc_cli_box[i];
    	   }
    	   // info on coords
    	   // the order is xyzxyz...
    	   for(unsigned i=0;i<3*natoms;i++){
    		   coordinates[i]=real(ts_in.coords[i]/10.); //convert to nm
    		   //cerr<<"COOR "<<coordinates[i]<<endl;
    	   }
#endif
       }else if(trajectory_fmt=="xdr-xtc" || trajectory_fmt=="xdr-trr"){
#ifdef __PLUMED_HAS_XDRFILE
         int step;
         float time;
         matrix box;
         rvec* pos=new rvec[natoms];
         float prec,lambda;
         int ret;
         if(trajectory_fmt=="xdr-xtc") ret=read_xtc(xd,natoms,&step,&time,box,pos,&prec);
         if(trajectory_fmt=="xdr-trr") ret=read_trr(xd,natoms,&step,&time,&lambda,box,pos,NULL,NULL);
         if(ret==exdrENDOFFILE) break;
         if(ret!=exdrOK) break;
         for(unsigned i=0;i<3;i++) for(unsigned j=0;j<3;j++) cell[3*i+j]=box[i][j];
         for(unsigned i=0;i<natoms;i++) for(unsigned j=0;j<3;j++)
                 coordinates[3*i+j]=real(pos[i][j]);
         delete [] pos;
#endif
       }else{
       if(trajectory_fmt=="xyz"){
         if(!Tools::getline(fp,line)) error("premature end of trajectory file");

         std::vector<double> celld(9,0.0);
         if(pbc_cli_given==false) {
           std::vector<std::string> words;
           words=Tools::getWords(line);
           if(words.size()==3){
             sscanf(line.c_str(),"%100lf %100lf %100lf",&celld[0],&celld[4],&celld[8]);
           } else if(words.size()==9){
             sscanf(line.c_str(),"%100lf %100lf %100lf %100lf %100lf %100lf %100lf %100lf %100lf",
                    &celld[0], &celld[1], &celld[2],
                    &celld[3], &celld[4], &celld[5],
                    &celld[6], &celld[7], &celld[8]);
           } else error("needed box in second line of xyz file");
         } else {			// from command line
           celld=pbc_cli_box;
         }
         for(unsigned i=0;i<9;i++)cell[i]=real(celld[i]);
       }
  	   int ddist=0;
       // Read coordinates
       for(int i=0;i<natoms;i++){
         bool ok=Tools::getline(fp,line);
         if(!ok) error("premature end of trajectory file");
         double cc[3];
         if(trajectory_fmt=="xyz"){
           char dummy[1000];
           int ret=std::sscanf(line.c_str(),"%999s %100lf %100lf %100lf",dummy,&cc[0],&cc[1],&cc[2]);
           if(ret!=4) error("cannot read line"+line);
         } else if(trajectory_fmt=="gro"){
           // do the gromacs way
           if(!i){
        	   //
        	   // calculate the distance between dots (as in gromacs gmxlib/confio.c, routine get_w_conf )
        	   //
        	   const char      *p1, *p2, *p3;
        	   p1 = strchr(line.c_str(), '.');
        	   if (p1 == NULL) error("seems there are no coordinates in the gro file");
        	   p2 = strchr(&p1[1], '.');
        	   if (p2 == NULL) error("seems there is only one coordinates in the gro file");
        	   ddist = p2 - p1;
        	   p3 = strchr(&p2[1], '.');
        	   if (p3 == NULL)error("seems there are only two coordinates in the gro file");
        	   if (p3 - p2 != ddist)error("not uniform spacing in fields in the gro file");
           }
           Tools::convert(line.substr(20,ddist),cc[0]);
           Tools::convert(line.substr(20+ddist,ddist),cc[1]);
           Tools::convert(line.substr(20+ddist+ddist,ddist),cc[2]);
         } else plumed_error();
         if(!debug_pd || ( i>=pd_start && i<pd_start+pd_nlocal) ){
           coordinates[3*i]=real(cc[0]);
           coordinates[3*i+1]=real(cc[1]);
           coordinates[3*i+2]=real(cc[2]);
         }
       }
       if(trajectory_fmt=="gro"){
         if(!Tools::getline(fp,line)) error("premature end of trajectory file");
         std::vector<string> words=Tools::getWords(line);
         if(words.size()<3) error("cannot understand box format");
         Tools::convert(words[0],cell[0]);
         Tools::convert(words[1],cell[4]);
         Tools::convert(words[2],cell[8]);
         if(words.size()>3) Tools::convert(words[3],cell[1]);
         if(words.size()>4) Tools::convert(words[4],cell[2]);
         if(words.size()>5) Tools::convert(words[5],cell[3]);
         if(words.size()>6) Tools::convert(words[6],cell[5]);
         if(words.size()>7) Tools::convert(words[7],cell[6]);
         if(words.size()>8) Tools::convert(words[8],cell[7]);
       }

     }

       if(debug_dd){
         for(int i=0;i<dd_nlocal;++i){
           int kk=dd_gatindex[i];
           dd_coordinates[3*i+0]=coordinates[3*kk+0];
           dd_coordinates[3*i+1]=coordinates[3*kk+1];
           dd_coordinates[3*i+2]=coordinates[3*kk+2];
         }
         p.cmd("setForces",&dd_forces[0]);
         p.cmd("setPositions",&dd_coordinates[0]);
         p.cmd("setMasses",&dd_masses[0]);
         p.cmd("setCharges",&dd_charges[0]);
       } else {
         p.cmd("setForces",&forces[3*pd_start]);
         p.cmd("setPositions",&coordinates[3*pd_start]);
         p.cmd("setMasses",&masses[pd_start]);
         p.cmd("setCharges",&charges[pd_start]);
       }
       p.cmd("setBox",&cell[0]);
       p.cmd("setVirial",&virial[0]);
   }
   p.cmd("calc");

// this is necessary as only processor zero is adding to the virial:
   intracomm.Bcast(virial,0);
   if(debug_pd) intracomm.Sum(forces);
   if(debug_dd){
     for(int i=0;i<dd_nlocal;i++){
       forces[3*dd_gatindex[i]+0]=dd_forces[3*i+0];
       forces[3*dd_gatindex[i]+1]=dd_forces[3*i+1];
       forces[3*dd_gatindex[i]+2]=dd_forces[3*i+2];
     }
     dd_forces.assign(3*natoms,0.0);
     intracomm.Sum(forces);
   }
   if(debug_grex &&step%grex_stride==0){
     p.cmd("GREX savePositions");
     if(intracomm.Get_rank()>0){
       p.cmd("GREX prepare");
     } else {
       int r=intercomm.Get_rank();
       int n=intercomm.Get_size();
       int partner=r+(2*((r+step/grex_stride)%2))-1;
       if(partner<0)partner=0;
       if(partner>=n) partner=n-1;
       p.cmd("GREX setPartner",&partner);
       p.cmd("GREX calculate");
       p.cmd("GREX shareAllDeltaBias");
       for(int i=0;i<n;i++){
         string s; Tools::convert(i,s);
         real a; s="GREX getDeltaBias "+s; p.cmd(s.c_str(),&a);
         if(grex_log) fprintf(grex_log," %f",a);
       }
       if(grex_log) fprintf(grex_log,"\n");
     }
   }


   if(fp_forces){
     fprintf(fp_forces,"%d\n",natoms);
     string fmtv=dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+"\n";
     string fmt=dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+"\n";
     if(dumpfullvirial){
       fprintf(fp_forces,fmtv.c_str(),virial[0],virial[1],virial[2],virial[3],virial[4],virial[5],virial[6],virial[7],virial[8]);
     } else {
       fprintf(fp_forces,fmt.c_str(),virial[0],virial[4],virial[8]);
     }
     fmt="X "+fmt;
     for(int i=0;i<natoms;i++)
       fprintf(fp_forces,fmt.c_str(),forces[3*i],forces[3*i+1],forces[3*i+2]);
   }

    if(noatoms && plumedStopCondition) break;

    step+=stride;
  }
  p.cmd("runFinalJobs");

  if(fp_forces) fclose(fp_forces);
  if(fp && fp!=in)fclose(fp);
#ifdef __PLUMED_HAS_XDRFILE
  if(xd) xdrfile_close(xd);
#endif
#ifdef __PLUMED_HAS_MOLFILE
  if(h_in) api->close_file_read(h_in);
  if(ts_in.coords) delete [] ts_in.coords;
#endif
  if(grex_log) fclose(grex_log);
*/
  return 0;
}

typedef Socket<double> SocketDouble;

/// Specialized version
template<>
string Socket<double>::description()const{ return "run plumed over a socket"; }

PLUMED_REGISTER_CLTOOL(SocketDouble,"socket")

}
}
