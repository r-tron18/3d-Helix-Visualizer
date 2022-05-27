#include<iostream>
#include<cmath>

#define M_ERROR (0.00001)
using namespace std;

class Point
{
double x;double y;double z;

public:
	Point(): x(0),y(0),z(0) {}
	Point(double a,double b) : x(a),y(b),z(0) {}
	Point(double a,double b,double c) : x(a),y(b),z(c) {}

	double _x() {return x;}
	double _y() {return y;}
	double _z() {return z;}

	void display()
		{
		//cout.precision(5);
		cout.setf(ios::fixed,ios::floatfield);
		cout<<"The point is ("<<x<<", "<<y<<", "<<z<<")\n";
		}
	
	double distanceFrom(Point &other)
		{
			return sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y) + (z - other.z) * (z - other.z)) ;
		}
	double magnitude()
		{
			return sqrt(x*x + y*y + z*z);
		}
	
	double mod()
		{
			return (x*x + y*y + z*z);
		}

	Point operator -()  // unary minus
		{
			return Point(-x,-y,-z);
		}

	Point operator -(Point &other)
		{
			return Point(this->x - other.x,this->y - other.y, this->z - other.z);
		}
	
	Point operator +(Point &other)
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

	double operator*(Point &other)   //  scalar product
		{
			return this->x * other.x + this->y * other.y + this->z * other.z;
		}
	
	friend Point operator*(Point &other,double l);
	friend Point operator*(double l,Point &other);
};

Point operator *(Point &other,double l)
{
	return Point(other.x * l, other.y * l, other.z * l);
}

Point operator*(double l,Point &other)
{
	return Point(l*other.x, l * other.y, l * other.z);
}

// -------------- Axis -------------
class Axis
{
Point p;
Point v;

public:
	Axis()
		{
		 p = Point(0,0,0);
		 v = Point(0,0,0);
		}

	Axis(Point &P,Point &V)
		{
			p=P; v=V;
		}
	
	Point _getP() {return p;}
  Point _getV() {return v;}

	Point per_point_on_axis(Point b)
		{
			double lambda = ((b-p)*v)/(v.mod());
			Point temp =   lambda * v + p;
			return temp;
		}
	
	Point leg_of(Point B,Point A) //  here B is a point in space while A is the point in plane
		{
			double d = A * v;
			double lambda = (d - (B*v)) / (v*v);

			return lambda * v + B;
		}

	void display()
		{        
			cout<<"\t The AXIS is \n";
			p.display();
			cout<<"Vector is ("<<v._x()<<", "<<v._y()<<", "<<v._z()<<")\n";
		}
};

int main()
  {  
    cout<<"\t\t  Testing Axis \n\n";
    
    cout<<"Creating Axis without values \n";
    Axis ax;
    ax.display();
    cout<<endl;
    
    cout<<"Creating Axis with Point(1,0,0) and vec(0,0,1) \n";
    Point o1(1,0,0),vec1(0,0,1);
    Axis ax1(o1,vec1);
    ax1.display();
    cout<<endl;
    
    cout<<"Finding Point X on previous axis using Point(1,1,4) such that ...\n";
    cout<<"\t a line passing through Point X and Point(1,1,4) is perpendicular to axis vector\n";
    Point p(1,1,4);
    Point perp_p = ax1.per_point_on_axis(p);
    perp_p.display();
    cout<<endl;
    
    cout<<"Creating Axis with Point(1,1,2) and vec(1,3,1) \n";
    Point o2(1,1,2),vec2(1,3,1);
    Axis ax2(o2, vec2);
    cout<<"Creating Point (2,2,2) \n";
    Point st(2,2,2);
    cout<<" The Plane in consideration is the one with Point (2,2,2) lying on plane and plane_vector (1,3,1)\n";
    cout<<" We create anothe point (10,20,30) \n";
    Point st2(10,20,30);
    cout<<" The perpendicular point x of (10,20,30) is such that distance between x and plane is minimum \n";
    Point perp_st2  = ax2.leg_of(st2,st);
    perp_st2.display();
    cout<<endl;
    
    return 0;    
  }
