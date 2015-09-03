/*

flex_roc.cc

Code for flexible ROC area calculation.
USAGE: ./fex_roc <datafile> <width> [opt: roc output file]
A positive prediction +/- width registers as success.
Note: setting width=0 is standard ROC calculation

dependency: gsl (for numerical integration)
compile command: g++ -O3 flex_roc.cc -o flex_roc -lgsl

Copyright 2011 Ishanu Chattopadhyay <ishanu@uchicago.edu>

*/
//------------------------------------------------------------
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.
//      
//-------------------------------------------------------------


#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <stdlib.h>
#include <set>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#define DEBUG_ 0

using namespace std;

const double TOL=1e-3;

//---------------------------------------------
//---------------------------------------------
unsigned int INF_RNG=0;
//---------------------------------------------
//---------------------------------------------

ostream& operator << (ostream &out, 
		      map<double,double> H)
{
  for(map<double,double>::iterator itr=H.begin();
      itr!=H.end();
      ++itr)
    out << itr->first 
	 << " " << itr->second 
	 <<endl;
  return out;
}
//---------------------------------------------
//---------------------------------------------

class prf__
{
private:
  map<string,double> perf;

public:
  prf__(){};
  void set(string ST, double val)
  {
    perf[ST] = val;
  };

  double get(string ST)
  {
    if(perf.find(ST)!=perf.end())
      return perf[ST];
    else
      return 0.0;
  };

  pair<double,double> get()
  {
    return make_pair(perf["FP"]/(perf["TN"]+perf["FP"]+0.0),
		     perf["TP"]/(perf["TP"]+perf["FN"]+0.0) );
  };

  friend ostream& operator << (ostream &out, prf__& P);
};
//---------------------------------------------
//---------------------------------------------
class numInt_; //forward declaration
double getx(double xi, void* N);
//---------------------------------------------
//---------------------------------------------
//---------------------------------------------
//---------------------------------------------
class numInt_
{
  vector <double> x,y;
  gsl_interp_accel *acc;
  gsl_spline *spline;

public:
  numInt_(map<double,double> H)
  {
    for(map<double,double>::iterator itr=H.begin();
	itr!=H.end();
	++itr)
      {
	x.push_back(itr->first);
	y.push_back(itr->second);
      }
    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc (gsl_interp_cspline, H.size());
    gsl_spline_init (spline, &x[0], &y[0], H.size());
  };

  double get(double xi) const
  {
    return gsl_spline_eval (spline, xi, acc);
  };

  double area(double tol=TOL)
  {
    gsl_integration_workspace * w 
      = gsl_integration_workspace_alloc (1000);
    double result, error;
    gsl_function F;
    F.function = &getx;
    F.params = (void*) this;
    gsl_integration_qags (&F, 0, 1, 0, tol, 1000,
			  w, &result, &error);
    return result;
  }
};
//---------------------------------------------
//---------------------------------------------

double getx(double xi, void* N)
{
  return ((numInt_*)N)->get(xi);
}
//---------------------------------------------
//---------------------------------------------
ostream& operator << (ostream &out, prf__& P)
{
  out << " TN " <<  P.perf["TN"]  
      << " FN " <<  P.perf["FN"]  
      << " TP " << P.perf["TP"] 
      << " FP " << P.perf["FP"] 
      <<  " FPR " << P.perf["FP"]/(P.perf["TN"]+P.perf["FP"]+0.0) 
      << " TPR " <<  P.perf["TP"]/(P.perf["TP"]+P.perf["FN"]+0.0)   
      << endl;
  return out;
}

//---------------------------------------------
//---------------------------------------------
prf__&  getdata(map<unsigned int, 
		pair<unsigned int, double> > H,
		double thr,
		prf__& PERF)
{
  vector <pair <unsigned int,double> > rocdata;

  set<set<unsigned int> >  rSets;
  set<unsigned int> S;

  unsigned int TN=0, FN=0, FP=0, TP=0;

  for(map<unsigned int, 
	pair<unsigned int,
	double> >::iterator itr=H.begin();
      itr!=H.end();
      ++itr)
    {
      if(itr->second.second > thr)
	{
	  bool ev_=false;
	  for(unsigned int j=itr->first-INF_RNG;
	      j<=itr->first+INF_RNG;
	      ++j)
	    if(H.find(j)!=H.end())
	      if(H[j].first == 1)
		ev_=true;
	  if(ev_)
	    TP+=2*INF_RNG+1;
	  else
	    FP++;
	  if(ev_)
	    for(unsigned int j=itr->first-INF_RNG;
		j<=itr->first+INF_RNG;
		++j)
	      S.insert(j);
	  else
	    S.insert(itr->first);
	}
    }
  for(map<unsigned int, 
	pair<unsigned int,
	double> >::iterator itr=H.begin();
      itr!=H.end();
      ++itr)
    if(S.find(itr->first)==S.end())
      {
	if(itr->second.first==1)
	  FN++;
	else
	  TN++;
      }
  PERF.set("TN", TN);
  PERF.set("FN", FN);
  PERF.set("TP", TP);
  PERF.set("FP", FP);

  return PERF;
};
//---------------------------------------------
//---------------------------------------------

int main(int argc, char* argv[])
{
  string rocfile="roc.dat";
  if(argc <=1)
    {
      cout << "missing datafile" << endl;
      exit(1);
    }
  if(argc>2)
    INF_RNG=atoi(argv[2]);
  if(argc>4)
    rocfile = argv[3];

  string datafile=argv[1];
  ifstream IN(datafile.c_str());
  string line;
  vector<double> vals;
  map<unsigned int, 
      pair<unsigned int, double> > H;

  while(getline(IN,line))
    {
      stringstream ss(line);
      unsigned int id,event;
      double val;
      while(ss>>id>>event>>val)
	H[id]=make_pair(event,val);
    }

  prf__ PERF;
  map<double,double> HH;
  ofstream RC(rocfile.c_str());

  for(double th=0.0001;th<0.99;th+=0.0001)
    {
      pair <double,double> valpair = getdata(H,th,PERF).get();
      if(HH.find(valpair.first)==HH.end())
	HH[valpair.first] = valpair.second;
      else
	{
	  if(HH[valpair.first] < valpair.second)
	    HH[valpair.first] = valpair.second;
	}
    }
  RC << HH;
  RC.close();

  numInt_ N(HH);
  cout << "ROC AREA: " << N.area() << endl;

  return 0;
}
//---------------------------------------------
//---------------------------------------------
