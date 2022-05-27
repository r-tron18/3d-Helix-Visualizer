#include<iostream>
#include<algorithm>
using namespace std;

void fix(int &a,int &b)
{
	int mn = min(a,b);
	int mx = max(a,b);
	a=mn; b=mx;
}

class solution
{
double x,y,z;
	public:
	solution(): x(0),y(0),z(0) {}
	solution(double a,double b,double c): x(a),y(b),z(c) {}
	
	void display()
		{
			cout<<"The Solution is\n";
			cout<<" X = "<<x<<endl;
			cout<<" Y = "<<y<<endl;
			cout<<" Z = "<<z<<endl;
			cout<<endl;
		}
};

class matrix
{
int row,col;
double Matrix[3][3];
public:
	
	matrix()
		{
			row=3;col=3;
			for(int i=0;i<row;i++)
				for(int j=0;j<col;j++)
						Matrix[i][j] = 0;
		}

	
	matrix(double arr[][3])
		{
			row=3;col=3;
			for(int i=0;i<row;i++)
				for(int j=0;j<col;j++)
						Matrix[i][j] = arr[i][j];
		}

	bool isinverse()
		{
			if (det()) return true;
			else return false;
		}

	matrix inverse()
		{
			double t = 1/this->det();
			matrix temp = this->adjoint() * t  ;	// adjoint of matrix A  /  |A|
			return temp;
		}
	
	matrix adjoint()
		{
			double arr[3][3];
			for(int i=0;i<3;i++)
				for(int j=0;j<3;j++)
						arr[i][j] = cofactor(i,j);
			matrix temp(arr);
			return temp.transpose();
		}
	friend matrix operator* (matrix a, double C);
	friend matrix operator* (double C, matrix a);
	
	matrix transpose()
		{
			matrix temp;
			for(int i=0;i<3;i++)
				for(int j=0;j<3;j++)
						temp.Matrix[i][j] = this->Matrix[j][i];

			return temp;
		}
	double det()
		{
			double DET=0;
			for(int j=0;j<3;j++) 
					DET+= Matrix[0][j]*cofactor(0,j);		
			return DET; 
		}

	double _minor(int i,int j)
		{
			int r1=(i+1)%3,r2=(i+2)%3;
			int c1=(j+1)%3,c2=(j+2)%3;

			fix(r1,r2);
			fix(c1,c2);

			//cout<<"r1= "<<r1<<" r2= "<<r2<<" c1= "<<c1<<" c2= "<<c2<<endl;
			return Matrix[r1][c1] * Matrix[r2][c2] - Matrix[r2][c1] * Matrix[r1][c2];
		}

	double cofactor(int i,int j)
		{
			if( (i+j)%2==0 ) {return   _minor(i,j); } 
			else          {return  -_minor(i,j); }		
		}

	void display()
		{
			for(int i=0;i<row;i++)
			{
				for(int j=0;j<col;j++)
					{
						cout<<Matrix[i][j]<<" ";
					}
			cout<<endl;			
			}
		}

	friend solution linear_system_sol(matrix a,double arr[3]);

};

solution linear_system_sol(matrix m,double arr[3])
{
double X=0,Y=0,Z=0;
matrix mat = m.inverse();
for(int j=0;j<mat.col;j++)
		X += mat.Matrix[0][j] * arr[j];				

for(int j=0;j<mat.col;j++)
		Y += mat.Matrix[1][j] * arr[j];				

for(int j=0;j<mat.col;j++)
		Z += mat.Matrix[2][j] * arr[j];				

	return solution(X,Y,Z);
}

matrix operator* (matrix a, double C)
{
	matrix temp;
	for(int i=0;i<a.row;i++)
		for(int j=0;j<a.col;j++)
			temp.Matrix[i][j] = a.Matrix[i][j] * C;

	return temp;
}

matrix operator* (double C,matrix a)
{
	matrix temp;
	for(int i=0;i<a.row;i++)
		for(int j=0;j<a.col;j++)
			temp.Matrix[i][j] = a.Matrix[i][j] * C;

	return temp;
}

int main()
{
// test code
cout<<"\t Testing matrix class\n";
char var[3]= {'x','y','z'};
double arr[3][3] ={
										{1,1,5},
										{0,2,5},
										{2,5,-1}
		};
double d[3] = {6,-4,27};

matrix temp;
cout<<"Default matrix made \n";
cout<<"\t Matrix is\n";
temp.display();
cout<<endl;

cout<<"Matrix made with arr value \n";
matrix mat1(arr);
cout<<"\t Matrix is\n";
mat1.display();
cout<<endl;

cout<<"\t Determinant is\n";
cout<<mat1.det();
cout<<endl;

matrix trans = mat1.transpose();
cout<<"\t Transpose Matrix is\n";
trans.display();
cout<<endl;

cout<<"\tminors are \n";
for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
		cout<<"Minor ["<<i<<"]["<<j<<"] : "<<mat1._minor(i,j)<<endl;
cout<<endl;

cout<<"\t Cofactors are \n";
for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
		cout<<"Cofactor ["<<i<<"]["<<j<<"] : "<<mat1.cofactor(i,j)<<endl;
cout<<endl;

matrix adj = mat1.adjoint();
cout<<"\t Adjoint Matrix is\n";
adj.display();
cout<<endl;

matrix inv = mat1.inverse();
cout<<"\t Inverse Matrix is\n";
inv.display();
cout<<endl;


cout<<"\t The solution for system of linear equations \n";
for(int i=0;i<3;i++)
	{
		for(int j=0;j<2;j++)
			cout<<arr[i][j]<<var[j]<<" + ";
		cout<<arr[i][2]<<var[2]<<" = "; 
		cout<<d[i]<<endl;
	}
cout<<endl;
solution s = linear_system_sol(mat1,d);
s.display();

return 0;
}
