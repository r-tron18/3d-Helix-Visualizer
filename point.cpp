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


int main()
  {
    cout<<"\t\t  Testing Point \n\n";
    
    cout<<"Creating Point without values\n";
    Point a;
    a.display();
    cout<<endl;
    
    cout<<"Creating Point with 2 values 1, 2.34 \n";
    Point a1(1,2.34);
    a1.display();
    cout<<endl;
    
    cout<<"Creating Point with 3 values 11.23, 6, 7.64 \n";
    Point a2(11,6,7.64);
    a2.display();
    cout<<endl;
    
    cout<<"Creating two points (1,2,3)  (2,4,5) \n";
    cout<<"Distance between points is ";
    Point a3(1,2,3);
    Point a4(2,4,5);
    cout<<a3.distanceFrom(a4)<<"\n\n";
    
    cout<<"The magnitude of (1,2,3) is "<<a3.magnitude()<<endl;
    cout<<"The mod of (1,2,3) is "<<a3.mod()<<endl<<"\n";
    
    cout<<"The point is (1.234,0,0.234) \n"; Point a5(1.234,0,0.234);
    cout<<"Testing unary operator - \n"; a5 = -a5 ; a5.display(); 
    cout<<endl;
    
    cout<<"Testing binary - for two points (1,2,3)  (2,4,5) \n";
    Point a6 = a3-a4;
    a6.display();
    cout<<endl;
    
    cout<<"Testing binary + for two points (1,2,3)  (2,4,5) \n";
    Point a7 = a3+a4;
    a7.display();
    cout<<endl;
    
    cout<<"Testing binary == for two points (1,2,3)  (2,4,5) \n";
    if (a3==a4)cout<<"Equal Points\n";
    else       cout<<"Unequal Points\n";
    cout<<endl;
    
    
    cout<<"Testing binary == for two points (1,2,3)  (2,4,5) \n";
    Point a8(1,2,3);
    if (a3==a8)cout<<"Equal Points\n";
    else       cout<<"Unequal Points\n";
    cout<<endl;
    
    cout<<"Testing multiplication of scalar 5 with vector (11,.23,34.5) \n";
    Point vec(11,.23,34.5);
    Point new_vec = 5 * vec;
    new_vec.display();
    cout<<endl;
    
    cout<<"Testing multiplication of vector (7.9,1.3,4.15) with scalar 6\n";
    Point vec1(7.9,1.3,4.15);
    Point new_vec1 =vec1 * 6;
    new_vec1.display();
    cout<<endl;
    
    cout<<"Testing scalar Product of vector (1,2.3,4) and vector (12.9,6,1) \n";
    Point vec2(1,2.3,4);
    Point vec3(12.9,6,1);
    
    cout<<"The scalar Product is "<<vec2 * vec1<<endl;
    cout<<endl;
    
    return 0;

  }
