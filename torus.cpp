#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
#include<fstream>
#include<utility>
#include<map>

//  define preprocessor values

#define PI  (3.14159265358979323846)
#define M_ERROR (0.00001)


int height_mul = 1;

using namespace std;
ofstream fout("point.txt",ios::out); // for plotting purposes only

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

Point operator *(Point &other,double l)
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
	Point temp = getPoint(a,end_vec.magnitude(),_ax,acos(mag),0)	; //either send height here and change b_leg below by b

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

vector<Point> solve(Point ax_point,Vector ax_vec,double move_angle,Point start, Point end, double rotation_step)
 {

		vector<Point> ArcPoints;
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

		//  while loop variables		
		int i = 0;
		double current_angle  = 0 ;
		double current_radius = radius;
		double current_height = 0 ;
		
		while ( abs(current_angle) <= abs(new_rotation_step) )
			{
	
				Point temp = getPoint(start,current_radius,axis,current_angle,height_mul * current_height);
					
				ArcPoints.push_back(temp);
				temp.write();   //writes in a file
				i+=1;
				
				current_angle  = move_angle*i ;
				current_radius = radius + unit_radius * move_angle * i;
				current_height = unit_height * move_angle * i ;
			}
		//fout.close();
	return ArcPoints;
 }


//------*****--------------****------------****-----------*****-------

void disp_ans(map<string, pair<double, Point> > Map)
{
	cout<< "Min X " ; Map["min_x"].second.display() ;
	cout<< "Max X " ; Map["max_x"].second.display() ;
	cout<< "Min Y " ; Map["min_y"].second.display() ;
	cout<< "Max Y " ; Map["max_y"].second.display() ;
	cout<< "Min Z " ; Map["min_z"].second.display() ;
	cout<< "Max Z " ; Map["max_z"].second.display() ;
}

void update_ans(map<string, pair<double, Point> >& Map, Point p, int ch=0)
{
	double	x = p._x(), y = p._y(), z = p._z();

	if(x < Map["min_x"].first)
		Map["min_x"] = make_pair(x, p);
	if(x > Map["max_x"].first)
		Map["max_x"] = make_pair(x, p);

	if(y < Map["min_y"].first)
		Map["min_y"] = make_pair(y, p);
	if(y > Map["max_y"].first)
		Map["max_y"] = make_pair(y, p);

	if(z < Map["min_z"].first)
		Map["min_z"] = make_pair(z, p);
	else if(z > Map["max_z"].first)
		Map["max_z"] = make_pair(z, p);
}


void print_6_points(Point ax_point, Vector ax_vector, Point cir_point, Vector cir_vector, Point start, Point end)
{
		vector<Point> ArcPoints;

		Axis arcAxis(cir_point, cir_vector);
    Axis torusAxis(ax_point, ax_vector);

    // start and end point in space
    
    double Alpha = angle_between(end, start, arcAxis);
		cout<< rad_to_deg(Alpha )<< "\n";
    
		ArcPoints = solve(cir_point, cir_vector,deg_to_rad(5),start, end, Alpha);
		vector<Point> TorusPoints;
		
		for(auto apoint: ArcPoints)
			{
				Point ax_p =  torusAxis.per_point_on_axis(apoint) ;
				double len = apoint.distanceFrom(ax_p);
				Point en = getPoint(apoint, len, torusAxis, deg_to_rad(360), 0);
				vector<Point> TMP = solve(ax_point, ax_vector, deg_to_rad(5), apoint, en, deg_to_rad(360)); 
				TorusPoints.insert(TorusPoints.end(), TMP.begin(), TMP.end());
			}
		
		map<string, pair<double, Point> > Map;
		Point fp = TorusPoints[0];
		double x = fp._x(), y = fp._y(), z = fp._z();
		Map["min_x"] = make_pair(x, fp);
		Map["max_x"] = make_pair(x, fp);
		
		Map["min_y"] = make_pair(y, fp);
		Map["max_y"] = make_pair(y, fp);

		Map["min_z"] = make_pair(z, fp);
		Map["max_z"] = make_pair(z, fp);
		
		for(int i=1;i<TorusPoints.size();i++)
			{
				Point tmp = TorusPoints[i];
				update_ans(Map, tmp);
			}
		
		update_ans(Map, start, 1);
		update_ans(Map, end);
		
		disp_ans(Map);

	fout.close(); // Only for getting points in file so it can be plotted
}

// driver function
int main()
{
 Point axisPoint(0,0,0);
 Vector axisVector(0,0,1);

 Point CirPoint(1,0,0);
 Vector CirVector(0,1,0);

 Point start(1,0,1);
 Point end(1,0,-1);

	print_6_points(axisPoint, axisVector, CirPoint, CirVector, start, end);
return 0;
}


