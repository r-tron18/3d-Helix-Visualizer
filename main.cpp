#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
#include<fstream>

//  for random values
#include<ctime>
#include<cstdlib>

//  define preprocessor values

#define PI  (3.14159265358979323846)
#define M_ERROR (0.00001)
#define VERBOSE 0
#define INFO    1

int height_mul = 1;

using namespace std;
ofstream fout("point.txt",ios::out);

// ----------------  UTILITY FUNCTIONS ------------------

void fix(int &a,int &b)  
	{
		int mn = min(a,b);
		int mx = max(a,b);
		a=mn; b=mx;
	}

double rad_to_deg(double angle)
	{
		return (180 * angle)/PI;
	}

double deg_to_rad(double angle)
	{
		return (angle * PI) / 180 ;
	}

bool is_solvable(double a, double b)  // check whether solve function is callable based on angles a and b
	{
		if ((a>0 && b<0) || (a<0 && b>0)) return false;
		return true;
	}

// -------------- Point ----------

class Point
 {
	double x;double y;double z;
	
		public:

	Point() {}
	Point(double a,double b) : x(a),y(b),z(0) {}
	Point(double a,double b,double c) : x(a),y(b),z(c) {}

	double _x() {return x;}
	double _y() {return y;}
	double _z() {return z;}

	void write()
		{
			fout.setf(ios::fixed,ios::floatfield);  
			fout<<x<<","<<y<<","<<z<<"\n"; 		
		}

	void display()
		{
			cout.setf(ios::fixed,ios::floatfield);  
			cout<<"("<<x<<", "<<y<<", "<<z<<")\n";  
		}

	double magnitude() // magnitude is the square of distanceFrom value
		{
			return sqrt(x*x + y*y + z*z);
		}

	double mod()
		{
			return (x*x + y*y + z*z);
		}

	double distanceFrom(Point &other)   //  distance of current point to given point
		{
			return sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y) + (z - other.z) * (z - other.z)) ;
		}

	Point operator -()                  // unary minus
		{
			return Point(-x,-y,-z);
		}

	Point operator -(Point &other)      // binary minus
		{
			return Point(this->x - other.x,this->y - other.y, this->z - other.z);
		}
	
	Point operator +(Point &other)      // binary plus 
		{
			return Point(this->x + other.x,this->y + other.y, this->z + other.z);
		}

	bool operator == (Point & other)
		{
			if ( abs(x-other.x) < M_ERROR && abs(y-other.y) < M_ERROR && abs(z-other.z) < M_ERROR ) 
				{return true; }
			else
	   		{return false;}
		}

	double operator *(Point &other)   //  scalar product
		{
			return this->x * other.x + this->y * other.y + this->z * other.z;
		}
	
	friend Point operator*(Point &other,double l);   // multiplication of vector with scalar
	friend Point operator*(double l,Point &other);   // multiplication of scalar with vector
	friend Point vec_pro(Point a,Point b);           // vector product
};

Point vec_pro(Point a,Point b)
	{
		return Point(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
	}

Point operator*(Point &other,double l)
	{
		return Point(other.x * l, other.y * l, other.z * l);
	}

Point operator*(double l,Point &other)
	{
		return Point(l*other.x, l * other.y, l * other.z);
	}

typedef Point Vector;
//--------------------------------


// -------------- Axis -------------

class Axis
 {
	Point p;
	Vector v;

		public:
	Axis() : p(Point(0,0,0)), v(Vector(0,0,0)) {}
	Axis(Point P,Point V) : p(P), v(V) {}

	Point per_point_on_axis(Point b)
		{
			double lambda = ((b-p)*v)/(v.mod());
			Point temp =   lambda * v;
			temp = temp + p;
			return temp;
		}
	Axis operator -()
		{
			Point _v = -v;
			return Axis(p,_v);
		}
	void display()
		{        
			cout<<"\t The AXIS is \n";
			cout<<"Axis Point ";p.display();
			cout<<"Axis Vector ("<<v._x()<<", "<<v._y()<<", "<<v._z()<<")\n";
		}
	
	Point leg_of(Point B,Point A) //  here B is a point in space while A is the point in plane
		{
			double d = A * v;
			double lambda = (d - (B*v)) / (v*v);

			return lambda * v + B;
		}

	Point _getP() {return p;}
  Vector _getV() {return v;}
};

// ----------------------------------------------------

Point getPoint(Point a,double len_b,Axis axis,double angle,double height);

// ----------------------------------------------------


double simple_angle(Vector a,Vector b)  // Takes vector a and b and returns an angle in radian in range 0 to PI
{
		double mag = (a*b)/(a.magnitude() * b.magnitude());
		return acos(mag);
}

//  Takes 2 points a and b in space an axis and returns angle between points with respect to axis

double angle_between(Point b, Point a,Axis _ax) // angle returned is in radians in range 0 to 2*PI
{
 	Point b_leg       = _ax.leg_of(b,a); 
	Point start_perp  = _ax.per_point_on_axis(a);
  Point end_perp    = _ax.per_point_on_axis(b);
	Vector start_vec  = a  - start_perp;
	Vector end_vec    = b  - end_perp;
	
 	double mag = (start_vec*end_vec)/(start_vec.magnitude() * end_vec.magnitude());
		if (mag<-1) mag =-1;
		if (mag>1)  mag = 1;

	//double height = b.distanceFrom(b_leg);
	Point temp = getPoint(a,end_vec.magnitude(),_ax,acos(mag),0); //either send height here and change b_leg below by b

	if (temp == b_leg) 
		{
		 return acos(mag) ; 			
		}
	else{
		double angle = rad_to_deg( acos(mag ));
		angle = 360 - angle;
		return deg_to_rad(angle);
			}

// always angle is returned such that Point a is moved counter clockwise by angle radians to reach Point b
}


Point getPoint(Point a,double len_b,Axis axis,double angle,double height)
{
	Vector on = axis._getV();
	Point o = axis.per_point_on_axis(a);
	Vector oa = a - o;

	Vector j = vec_pro(on,oa);
	double c = (double) 1/(j.magnitude());
	double lambda = height/on.magnitude();	
  j = j * c;
	
	Point ans = o ;
	Point _tmp1 =   len_b * cos(angle) / oa.magnitude()  * oa ;
	Point _tmp2 =   len_b * sin(angle) * j  ;

	ans = ans + _tmp1;
	ans = ans + _tmp2;
	ans =  lambda * on + ans ;
	//Point ans_perp = axis.per_point_on_axis(ans);
	return ans;
}

// -------------------------------------------------------

void solve(Point ax_point,Vector ax_vec,double move_angle,Point start, Point end, double rotation_step)
 {

		//  define axis and other necessary points
		Axis axis(ax_point,ax_vec);
		Point start_perp = axis.per_point_on_axis(start);
		Point end_perp   = axis.per_point_on_axis(end);
		Point end_leg    = axis.leg_of(end,start);
		Vector start_to_end_vector = end - start;
		Vector start_vec = start - start_perp;
		Vector end_vec   = end   - end_perp;
		// ------------------------------------------
		

		//  The function requires axis_vec such that it goes in direction of end point
		
		double _angle = simple_angle(start_to_end_vector,axis._getV()); // _angle is in radian
			if ( rad_to_deg(_angle) > 90 ) {height_mul = -1;}


		//   abs_angle is the abs angle between end and start in the ccw or cw direction range (0, 2*PI)
		
 		double abs_angle = angle_between(end,start,axis);       // abs_angle is the angle in counter clockwose direction
			if (abs_angle < M_ERROR) abs_angle = 0;
			if ( move_angle < 0 ) {  abs_angle = 2*PI - abs_angle;} //	abs_angle in clock wise direction is 2*PI - angle

		
		// ------- approximate total angle
		double new_rotation_step;
		int k;		
			if (move_angle > 0)  k = (rotation_step  - abs_angle) /(2*PI);  // positive for ccw 
			else                 k = (rotation_step + abs_angle) /(2*PI);  // negative for cw 
 
		new_rotation_step = 2*PI*k;
			if (move_angle<0) new_rotation_step = new_rotation_step - abs_angle;  //cw  angle is negative
			else              new_rotation_step = new_rotation_step + abs_angle;  //ccw angle is positive

		//  Find radius and pitch
		double radius= start.distanceFrom(start_perp);
		double pitch = end.distanceFrom(end_leg);
	    if (pitch<M_ERROR) pitch = 0;
   	
		
		//  Find unit_height and unit_radius that should be changed per radian angle
		double unit_height =   pitch / new_rotation_step; // always positive
		double unit_radius =  ( end_vec.magnitude() - start_vec.magnitude() ) / new_rotation_step ;
		

		// ----------------- Print Informative Values --------------
			#ifndef VERBOSE
			#define VERBOSE 0
			#endif

			#ifndef INFO
			#define INFO 0
			#endif

			#if INFO
				{
					cout<<"\n\tGiven Rotational step : "<<rad_to_deg(rotation_step)<<endl;
    			cout<<"\tUsing Rotational step : "<<rad_to_deg(new_rotation_step)<<endl;
					cout<<"\tradius = "<<radius<<", pitch = "<<pitch<<endl;
    			cout<<"\n";
				}
			#endif
		// --------------

		//  while loop variables		
		int i = 0;
		double current_angle  = 0 ;
		double current_radius = radius;
		double current_height = 0 ;

		// the above is the properties of start

		/*
			Using getPoint While loop generates points based on 
					current_radius
					current_angle
					current_height
					axis
		*/
		
		while ( abs(current_angle) <= abs(new_rotation_step) )
			{
	
				Point temp = getPoint(start,current_radius,axis,current_angle,height_mul * current_height);
					
					#if VERBOSE
						{
							if (i) 
								cout<<"Intermediate -> " ; 
							else   
								cout<<"Start        -> " ;
							temp.display();
						}
					#endif

					temp.write();   //writes in a file
					i+=1;
				
				current_angle  = move_angle*i ;
				current_radius = radius + unit_radius * move_angle * i;
				current_height = unit_height * move_angle * i ;
			}
					#if VERBOSE
						{
							cout<<"End          -> " ;
							end.display();
							cout<<"------------------ done ------------------ \n";
						}
					#else
						{
							if ( !( VERBOSE<0 ))
								cout<<"Verbose mode off.\n";
							#undef VERBOSE
							#define VERBOSE -1			
						}
					#endif
		fout.close();
 }


//------*****--------------****------------****-----------*****-------


int main()
	{
		/* Axis paramters
    Point p(0,0,0);
    Vector vector(0,0,1);
    
    // start and end point in space
    Point start(1, 0, 0);
    Point end(-2, 0, 1);
    
    // angle     (the angle to be moved to find next point)
    // tot_angle  (The total angle between start and end which needs to be approximated to be correct)
    double angle = -3.0000, tot_angle = -180.0000 - 360.0000*4;
    
    // solve it
    if (is_solvable(angle,tot_angle)) 
    	solve(p,vector,deg_to_rad(angle),start,end,deg_to_rad(tot_angle)); 
    else 
      cout<<"Wrong Input\n";
    
    cout<<"\n------------END---------------\n";
*/

		srand(time(0));
    cout<<"  Arbitrary Case (Random Values) \n\n";
    Point p4(rand()%40,rand()%40,rand()%40);
    Point vector4(rand()%40,rand()%40,rand()%40);
    Axis tmp(p4,vector4); // not needed for solution ; only for display function purpose
    tmp.display(); cout<<"\n";
    
    Point start4(rand()%40,rand()%40,rand()%40);
    Point end4(rand()%40,rand()%40,rand()%40);
    
    double angle4=1, tot_angle4=360*3;
    
    cout<<"Start Point ";start4.display();
    cout<<"End Point   ";end4.display();
    cout<<"angle moved per step: "<< angle4 <<", tot_angle given: "<< tot_angle4 <<endl;
    
    solve(p4,vector4,deg_to_rad(angle4),start4,end4,deg_to_rad(tot_angle4)); 
    cout<<"\tTesting Done\n";
 
    return 0;

	}


