/***************************************************************************
 *   Copyright (C) 2007 by Bjorn Harpe,,,   *
 *   bjorn@Ouelong.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/* Modified for GP Shape estimation, 2020
 * Author: Sudharshan Suresh <suddhu@cmu.edu> */

#ifndef CONTOURS_H
#define CONTOURS_H
#include <vector>
#include <algorithm>

#define DIFFERENCE 0.000001
#define EQ(_x_,_y_) (((_x_-_y_<DIFFERENCE)&&(_y_-_x_<DIFFERENCE))?1:0)
#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])
#define MIN(x,y) (x<y?x:y)
#define MAX(x,y) (x>y?x:y)

using namespace std;

using cnt = std::vector<std::vector<double>>;  // added
using cnts = std::vector<cnt>;  // added

inline std::ostream& operator<< (std::ostream& stream, const cnt& c) // added
{   
    for(auto& pt:c){
        std::cout<<pt[0]<<" "<<pt[1]<<std::endl;
    }
    return stream;
}

inline std::ostream& operator<< (std::ostream& stream, const cnts& cs) // added
{   
    int i = 1;
    for(auto& c:cs){
        std::cout<<"Contour #"<<i<<" : \n"; 
        std::cout<<"-----------------------\n";
        std::cout<<c;
        i++;
    }
    return stream;
}

struct SVector
{
   double dx,dy;
};

struct SPoint
{
   SPoint(double x,double y){this->x=x;this->y=y;}
   SPoint(){}
   double x,y;
};

struct SPair
{
   SPair(SPoint p1, SPoint p2){this->p1=p1;this->p2=p2;}
   SPair reverse(){return(SPair(p2,p1));}
   SPoint p1,p2;
};

bool operator <(SPoint p1, SPoint p2);
bool operator <(SPair p1,SPair p2);
bool operator ==(SPoint p1,SPoint p2);
bool operator !=(SPoint p1,SPoint p2);
SPoint operator +=(SPoint p, SVector v);

class CContour
{
   public:
      CContour(){contour=NULL;}
      ~CContour();
      int merge(CContour *c);
      int reverse();
      int add_vector(SPoint start,SPoint end);
      int condense(double difference = 0.000000001);
      int dump();
      void getShapeContour(cnt& cont); // added
      bool closed(){return(_start==_end);}
      SPoint start(){return(_start);}
      SPoint end(){return(_end);}
      vector<SVector> *contour;
   private:
      SPoint _start,_end;
};

class CRaster
{
   public:
      CRaster(){};
      virtual double value(double, double){return(0);};
      virtual SPoint upper_bound(){return(SPoint(0,0));};
      virtual SPoint lower_bound(){return(SPoint(0,0));};
      virtual ~CRaster(){};
};

class CContourLevel
{
   public:
     CContourLevel(){contour_lines=NULL;raw=NULL;};
     int dump();
     int getShapeContour(cnts& conts); // added
     int merge();
     int consolidate();
     vector<CContour*> *contour_lines;
     vector<SPair> *raw;
     ~CContourLevel();
}; 

class CContourMap
{
   public:
      CContourMap();
      int generate_levels(double min,double max, int num);
      int add_segment(SPair t,int level);
      int dump();
      int getShapeContour(cnts& conts); // added
      void getBestShapeContour(cnt& cont,  const double& threshold); //added
      int contour(Eigen::MatrixXd& d,
	double *x,
	double *y);
      int consolidate();
      CContourLevel* level(int i){return((*contour_level)[i]);}
      double alt(int i){return(levels[i]);}
      ~CContourMap();
   private:
      vector<CContourLevel*> *contour_level;
      int n_levels;
      double *levels;      
};
#endif