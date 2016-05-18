#include "../algebra/Exception.hpp"
#include "../algebra/Vector.hpp"
#include "../algebra/Matrix.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>



using namespace std;



unsigned int assertTrue(bool b, string description)
{
	if(b){
		cout << "\033[1;32mTest OK\033[0m\t\t" << description << endl;
		return 0;
	} else {
		cout << "\033[1;31mTest FAILED\033[0m\t" << description << endl;
		return 1;
	}
}



template<typename T, typename U>
unsigned int assertEqual(T a, U b, string description)
{
	if(a==b){
		cout << "\033[1;32mTest OK\033[0m\t\t" << description << " -> " << a << " == " << b << endl;
		return 0;
	} else {
		cout << "\033[1;31mTest FAILED\033[0m\t" << description << " -> " << a << " != " << b<< endl;
		return 1;
	}
}



template<typename T, typename U>
unsigned int assertDifferent(T a, U b, string description)
{
	if(!(a==b)){
		cout << "\033[1;32mTest OK\033[0m\t\t" << description << " -> " << a << " != " << b << endl;
		return 0;
	} else {
		cout << "\033[1;31mTest FAILED\033[0m\t" << description << " -> " << a << " == " << b<< endl;
		return 1;
	}
}



template<typename T, typename U, typename V>
unsigned int assertEqual(T a, U b, V diff, string description)
{
	if(abs(a-b)<diff){
		cout << "\033[1;32mTest OK\033[0m\t\t" << description << " -> " << a << " = " << b << endl;
		return 0;
	} else {
		cout << "\033[1;31mTest FAILED\033[0m\t" << description << " -> " << a << " != " << b<< endl;
		return 1;
	}
}



unsigned int test_vector_constructors(bool verbose)
{
	cout << "\033[1;36mRun Test: test_vector_constructors\033[0m" << endl;

	default_random_engine gen(chrono::high_resolution_clock::now().time_since_epoch().count());
    uniform_real_distribution<double> rdist(-100,100);
	unsigned int failcount = 0;

	// standard constructor
	Vector v;
	failcount += assertEqual(v[0],0.,"standard constructor at 0");
	failcount += assertEqual(v[1],0.,"standard constructor at 1");
	failcount += assertEqual(v[2],0.,"standard constructor at 2");		
	failcount += assertEqual(v.dim(),3,"standard constructor dimension");	
	failcount += assertTrue(v.isNullvector(),"standard constructor is null");
	failcount += assertEqual(v.norm(),0,"standard constructor norm");	
	try { v[4]; } catch( Exception& e ){ failcount += assertTrue(true,e.desc()); }
	try { v[-1]; } catch( Exception& e ){ failcount += assertTrue(true,e.desc()); }		
	
	// 2d constructor
	double x0 = rdist(gen);
	double x1 = rdist(gen);
	Vector v2d(x0,x1);
	failcount += assertEqual(v2d[0],x0,"2d constructor at 0");
	failcount += assertEqual(v2d[1],x1,"2d constructor at 1");
	failcount += assertEqual(v2d.dim(),2,"2d constructor dimension");	
	failcount += assertTrue(!v2d.isNullvector(),"2d constructor is not null");
	failcount += assertEqual(v2d.norm(),sqrt(x0*x0+x1*x1),"2d constructor norm");		

	// 3d constructor
	double x2 = rdist(gen);
	Vector v3d(x0,x1,x2);
	failcount += assertEqual(v3d[0],x0,"3d constructor at 0");
	failcount += assertEqual(v3d[1],x1,"3d constructor at 1");
	failcount += assertEqual(v3d[2],x2,"3d constructor at 2");	
	failcount += assertEqual(v3d.dim(),3,"3d constructor dimension");	
	failcount += assertTrue(!v3d.isNullvector(),"3d constructor is not null");
	failcount += assertEqual(v3d.norm(),sqrt(x0*x0+x1*x1+x2*x2),"3d constructor norm");	
	
	// 6d constructor
	double x3 = rdist(gen);
	double x4 = rdist(gen);
	double x5 = rdist(gen);
	Vector v6d(x0,x1,x2,x3,x4,x5);
	failcount += assertEqual(v6d[0],x0,"6d constructor at 0");
	failcount += assertEqual(v6d[1],x1,"6d constructor at 1");
	failcount += assertEqual(v6d[2],x2,"6d constructor at 2");	
	failcount += assertEqual(v6d[3],x3,"6d constructor at 3");
	failcount += assertEqual(v6d[4],x4,"6d constructor at 4");	
	failcount += assertEqual(v6d[5],x5,"6d constructor at 5");		
	failcount += assertEqual(v6d.dim(),6,"6d constructor dimension");	
	failcount += assertTrue(!v6d.isNullvector(),"6d constructor is not null");
	failcount += assertEqual(v6d.norm(),sqrt(x0*x0+x1*x1+x2*x2+x3*x3+x4*x4+x5*x5),"6d constructor norm");	
	
	// copy constructor
	Vector v6d_copy(v6d);
	failcount += assertEqual(v6d[0],v6d_copy[0],"copy constructor at 0");
	failcount += assertEqual(v6d[1],v6d_copy[1],"copy constructor at 1");
	failcount += assertEqual(v6d[2],v6d_copy[2],"copy constructor at 2");	
	failcount += assertEqual(v6d[3],v6d_copy[3],"copy constructor at 3");
	failcount += assertEqual(v6d[4],v6d_copy[4],"copy constructor at 4");	
	failcount += assertEqual(v6d[5],v6d_copy[5],"copy constructor at 5");		
	failcount += assertEqual(v6d_copy.dim(),6,"copy constructor dimension");	
	failcount += assertTrue(!v6d_copy.isNullvector(),"copy constructor is not null");
	failcount += assertEqual(v6d_copy.norm(),sqrt(x0*x0+x1*x1+x2*x2+x3*x3+x4*x4+x5*x5),"6d constructor norm");	
		
	if( verbose ){
		cout << v << endl;
		cout << v2d << endl;	
		cout << v3d << endl;			
		cout << v6d << endl;			
		cout << v6d_copy << endl;
	}		
	
	return failcount;
}



unsigned int test_vector_operators(bool verbose)
{
	cout << "\033[1;36mRun Test: test_vector_operators\033[0m" << endl;

	default_random_engine gen(chrono::high_resolution_clock::now().time_since_epoch().count());
    uniform_real_distribution<double> rdist(-100,100);
	unsigned int failcount = 0;
	
	double x0 = rdist(gen);
	double x1 = rdist(gen);
	double x2 = rdist(gen);		
	double y0 = rdist(gen);	
	double y1 = rdist(gen);
	double y2 = rdist(gen);		
		
	Vector x(x0,x1,x2);
	Vector u(-x0,-x1,-x2);	
	Vector y(y0,y1,y2);
	Vector z;

	// get and set
	Vector t(x);
	failcount += assertEqual(t[0],x0,"operator [] get");
	t[0] = x1;
	failcount += assertEqual(t[0],x1,"operator [] set");
		
	// negation operator
	Vector minus_x = -x;
	failcount += assertEqual(minus_x[0],-x0,"operator -x[0]");	
	failcount += assertEqual(minus_x[1],-x1,"operator -x[1]");	
	failcount += assertEqual(minus_x[2],-x2,"operator -x[2]");
	failcount += assertEqual(minus_x.dim(),x.dim(),"operator -x dimension");			
	
	// comparison operators
	failcount += assertTrue(x==x,"operator x==x");
	failcount += assertTrue(u==u,"operator u==u");	
	failcount += assertTrue(x!=u,"operator x!=u");
	failcount += assertEqual(y,y,"operator y==y");
	failcount += assertDifferent(y,z,"operator y!=z");	
	failcount += assertEqual(-y,-y,"operator -y==-y");
	failcount += assertEqual(x,-u,"operator x==-u");
	failcount += assertEqual(u,-x,"operator u==-x");	
	failcount += assertEqual(z,-z,"operator z==-z");
	failcount += assertDifferent(x,y,"operator x!=y");	
	
	// copy operator
	Vector q=x;
	failcount += assertEqual(q,x,"copy operator q=x");
		
	// add and subtract
	Vector sum = x+y;	
	Vector dif = x-y;		
	failcount += assertEqual(sum[0],x[0]+y[0],"add vectors 0");
	failcount += assertEqual(sum[1],x[1]+y[1],"add vectors 1");
	failcount += assertEqual(sum[2],x[2]+y[2],"add vectors 2");	
	failcount += assertEqual(z,x+u,"add vectors to 0");	
	failcount += assertEqual(dif[0],x[0]-y[0],"subtract vectors 0");
	failcount += assertEqual(dif[1],x[1]-y[1],"subtract vectors 1");
	failcount += assertEqual(dif[2],x[2]-y[2],"subtract vectors 2");			
	failcount += assertEqual(z,x-x,"subtract vectors to 0");
	failcount += assertEqual(x+y,y+x,"commutation");
	failcount += assertEqual(x-y,-(y-x),"commutation 2");
	try { Vector(2,1)+x; } catch( Exception& e ){ failcount += assertTrue(true,e.desc()); }	
	
	// factor and division
	double f = rdist(gen);
	double d = rdist(gen); if (d==0){d=1;}
	Vector fx = f*x;
	Vector xd = x/d;
	failcount += assertEqual(fx[0],f*x[0],"mult factor 0");
	failcount += assertEqual(fx[1],f*x[1],"mult factor 1");
	failcount += assertEqual(fx[2],f*x[2],"mult factor 2");		
	failcount += assertEqual(-f*x,-(f*x),"factor in braces");	
	failcount += assertEqual(0*x,z,"factor 0");	
	failcount += assertEqual(xd[0],x[0]/d,0.00001,"divide 0");
	failcount += assertEqual(xd[1],x[1]/d,0.00001,"divide 1");
	failcount += assertEqual(xd[2],x[2]/d,0.00001,"divide 2");
	failcount += assertEqual(2*x/2,x,"factor and division");	
	failcount += assertEqual(f*-x,-f*(2*x/2),"some braces");	
	
	// algebra
	Vector n(2,-1,5);
	Vector m(-1,0,4);
	double f2 = 2;
	failcount += assertEqual(f2*(n+m),f2*n+f2*m,"algebra 1");
	failcount += assertEqual((3*(2*x-4*y-x))[0],(3*x-12*y)[0],0.00001,"algebra 2_0");
	failcount += assertEqual((3*(2*x-4*y-x))[1],(3*x-12*y)[1],0.00001,"algebra 2_1");
	failcount += assertEqual((3*(2*x-4*y-x))[2],(3*x-12*y)[2],0.00001,"algebra 2_2");	
	
	// scalar product and norm
	failcount += assertEqual(z*z,0,"scalar product null vector");
	failcount += assertEqual(x*x,u*u,"scalar product x*x and -x*-x");	
	failcount += assertEqual(x*x,x0*x0+x1*x1+x2*x2,"scalar product x*x");
	failcount += assertEqual(x*y,x0*y0+x1*y1+x2*y2,"scalar product x*y");	
	failcount += assertEqual(x*y,y*x,"scalar product x*y == y*x");
	failcount += assertEqual(abs(x),sqrt(x0*x0+x1*x1+x2*x2),"norm x");
	failcount += assertEqual(abs(y),sqrt(y0*y0+y1*y1+y2*y2),"norm y");
	failcount += assertEqual(x.norm(),u.norm(),"norm x = norm -x");		
	failcount += assertEqual(x.norm(),sqrt(x*x),"norm x 2");		
	failcount += assertEqual((2*x).norm(),2*abs(x),"norm(2x) = 2abs(x)");
	failcount += assertEqual(abs(x.normalize()),1.,0.00001,"normalization");	
	try { z.normalize(); } catch( Exception& e ){ failcount += assertTrue(true,e.desc()); }

	if( verbose ){
		cout << x << endl;
		cout << u << endl;	
		cout << y << endl;
		cout << z << endl;
	}	
	
	return failcount;
}



unsigned int test_vector_angle_and_cross(bool verbose)
{
	cout << "\033[1;36mRun Test: test_vector_angle_and_cross\033[0m" << endl;

	default_random_engine gen(chrono::high_resolution_clock::now().time_since_epoch().count());
    uniform_real_distribution<double> rdist(-100,100);
	unsigned int failcount = 0;
	
	double x0 = rdist(gen);
	double x1 = rdist(gen);
	double x2 = rdist(gen);		
	double y0 = rdist(gen);	
	double y1 = rdist(gen);
	double y2 = rdist(gen);	
		
	Vector x(x0,x1);	
	Vector y(y0,y1);
	Vector z(0,0);
	Vector e1(1,0);
	Vector e2(0,1);	
	Vector v45(1,1);
	Vector xx(x0,x1,x2);
	Vector yy(y0,y1,y2);	
	
	// angles
	failcount += assertEqual(angle(e1,e2),0.5*M_PI,0.0000001,"angle 90°");	
	failcount += assertEqual(angle(x,x),0,0.0000001,"angle self");
	failcount += assertEqual(angle(e1,-e1),M_PI,0.0000001,"angle 180°");	
	failcount += assertEqual(angle(e1,-e2),0.5*M_PI,0.0000001,"angle 270°");	
	failcount += assertEqual(angle(e1,v45),0.25*M_PI,0.0000001,"angle 45°");		

	// cross product
	Vector c = xx.cross(yy);
	failcount += assertEqual(c[0],xx[1]*yy[2]-xx[2]*yy[1],0.0000001,"cross 0");
	failcount += assertEqual(c[1],xx[2]*yy[0]-xx[0]*yy[2],0.0000001,"cross 1");
	failcount += assertEqual(c[2],xx[0]*yy[1]-xx[1]*yy[0],0.0000001,"cross 2");
	double o = abs(c)/(abs(xx)*abs(yy));
	if( o>1 ){ o=1; }	if( o<-1 ){ o=-1; }	
	double p = (xx*yy)/(abs(xx)*abs(yy));
	if( p>1 ){ p=1; }	if( p<-1 ){ p=-1; }		
	failcount += assertEqual(abs(cos(asin(o))),abs(p),0.0000001,"criss cross ");	

	if( verbose ){
		cout << x << endl;
		cout << y << endl;
		cout << z << endl;
		cout << e1 << endl;
		cout << e2 << endl;
		cout << v45 << endl;		
	}	
	
	return failcount;
}



unsigned int test_matrix_constructors(bool verbose)
{
	cout << "\033[1;36mRun Test: test_matrix_constructors\033[0m" << endl;

	default_random_engine gen(chrono::high_resolution_clock::now().time_since_epoch().count());
    uniform_real_distribution<double> rdist(-100,100);
	unsigned int failcount = 0;

	// standard constructor
	Matrix m;
	failcount += assertEqual(m.rows(),3,"standard constructor rows");
	failcount += assertEqual(m.cols(),3,"standard constructor columns");	
	failcount += assertEqual(m[0][0],0,"standard constructor element 0,0");
	failcount += assertEqual(m[0][1],0,"standard constructor element 0,1");	
	failcount += assertEqual(m[0][2],0,"standard constructor element 0,2");
	failcount += assertEqual(m[1][0],0,"standard constructor element 1,0");
	failcount += assertEqual(m[1][1],0,"standard constructor element 1,1");
	failcount += assertEqual(m[1][2],0,"standard constructor element 1,2");
	failcount += assertEqual(m[2][0],0,"standard constructor element 2,0");
	failcount += assertEqual(m[2][1],0,"standard constructor element 2,1");
	failcount += assertEqual(m[2][2],0,"standard constructor element 2,2");
	try { m[2][3]; } catch( Exception& e ){ failcount += assertTrue(true,e.desc()); }
	try { m[3][0]; } catch( Exception& e ){ failcount += assertTrue(true,e.desc()); }
	
	// rows-columns constructor
	Matrix m2(8,2);
	failcount += assertEqual(m2.rows(),8,"rows-columns constructor rows");
	failcount += assertEqual(m2.cols(),2,"rows-columns constructor columns");		
	failcount += assertEqual(m2[3][0],0,"rows-columns constructor element 3,0");
	failcount += assertEqual(m2[0][1],0,"rows-columns constructor element 0,1");	
	failcount += assertEqual(m2[4][0],0,"rows-columns constructor element 4,0");
	failcount += assertEqual(m2[1][0],0,"rows-columns constructor element 1,0");
	failcount += assertEqual(m2[5][1],0,"rows-columns constructor element 5,1");
	failcount += assertEqual(m2[1][0],0,"rows-columns constructor element 1,0");
	failcount += assertEqual(m2[6][0],0,"rows-columns constructor element 6,0");
	failcount += assertEqual(m2[2][1],0,"rows-columns constructor element 2,1");
	failcount += assertEqual(m2[7][1],0,"rows-columns constructor element 7,1");	
	try { m2[2][3]; } catch( Exception& e ){ failcount += assertTrue(true,e.desc()); }
	try { m2[-1][0]; } catch( Exception& e ){ failcount += assertTrue(true,e.desc()); }
	
	// copy constructor
	Matrix m3(2,3);
	m3[0][0] = rdist(gen);
	m3[0][1] = rdist(gen);
	m3[0][2] = rdist(gen);
	m3[1][0] = rdist(gen);
	m3[1][1] = rdist(gen);
	m3[1][2] = rdist(gen);
	Matrix m3_copy(m3);
	failcount += assertEqual(m3.rows(),m3_copy.rows(),"copy constructor rows");
	failcount += assertEqual(m3.cols(),m3_copy.cols(),"copy constructor columns");		
	failcount += assertEqual(m3[0][0],m3_copy[0][0],"copy constructor element [0][0]");
	failcount += assertEqual(m3[0][1],m3_copy[0][1],"copy constructor element [0][1]");
	failcount += assertEqual(m3[0][2],m3_copy[0][2],"copy constructor element [0][2]");
	failcount += assertEqual(m3[1][0],m3_copy[1][0],"copy constructor element [1][0]");
	failcount += assertEqual(m3[1][1],m3_copy[1][1],"copy constructor element [1][1]");
	failcount += assertEqual(m3[1][2],m3_copy[1][2],"copy constructor element [1][2]");

	if( verbose ){
		cout << m << endl;
		cout << m2 << endl;		
		cout << m3 << endl;		
	}		
	
	return failcount;
}



unsigned int test_matrix_operators(bool verbose)
{
	cout << "\033[1;36mRun Test: test_matrix_operators\033[0m" << endl;

	default_random_engine gen(chrono::high_resolution_clock::now().time_since_epoch().count());
    uniform_real_distribution<double> rdist(-100,100);
	unsigned int failcount = 0;

	Matrix m;
	Matrix m1(2,3);
	m1[0][0] = rdist(gen);
	m1[0][1] = rdist(gen);
	m1[0][2] = rdist(gen);
	m1[1][0] = rdist(gen);
	m1[1][1] = rdist(gen);
	m1[1][2] = rdist(gen);
	Matrix m2(2,3);
	m2[0][0] = rdist(gen);
	m2[0][1] = rdist(gen);
	m2[0][2] = rdist(gen);	
	m2[1][0] = rdist(gen);
	m2[1][1] = rdist(gen);
	m2[1][2] = rdist(gen);		

	// operator == and !=
	failcount += assertEqual(m1,m1,"operator == 1");
	failcount += assertEqual(m2,m2,"operator == 2");
	failcount += assertDifferent(m,m1,"operator != 1");	
	failcount += assertDifferent(m,m2,"operator != 2");	
	failcount += assertDifferent(m1,m2,"operator != 2");
	
	// operator =
	m = m1;	
	failcount += assertEqual(m,m1,"operator = 1");		
	m = m2;
	failcount += assertEqual(m,m2,"operator = 2");
	
	// operator + and -
	Matrix sum = m1+m2;
	failcount += assertEqual(sum[0][0],m1[0][0]+m2[0][0],"standard sum element [0][0]");
	failcount += assertEqual(sum[0][1],m1[0][1]+m2[0][1],"standard sum element [0][1]");
	failcount += assertEqual(sum[0][2],m1[0][2]+m2[0][2],"standard sum element [0][2]");	
	failcount += assertEqual(sum[1][0],m1[1][0]+m2[1][0],"standard sum element [1][0]");
	failcount += assertEqual(sum[1][1],m1[1][1]+m2[1][1],"standard sum element [1][1]");
	failcount += assertEqual(sum[1][2],m1[1][2]+m2[1][2],"standard sum element [1][2]");	
	Matrix sum2(2,3);
	sum2+=m1;
	sum2+=m2;
	failcount += assertEqual(sum2,sum,"standard sum 2");	
	Matrix diff = m1-m2;
	failcount += assertEqual(diff[0][0],m1[0][0]-m2[0][0],"standard diff element [0][0]");
	failcount += assertEqual(diff[0][1],m1[0][1]-m2[0][1],"standard diff element [0][1]");
	failcount += assertEqual(diff[0][2],m1[0][2]-m2[0][2],"standard diff element [0][2]");	
	failcount += assertEqual(diff[1][0],m1[1][0]-m2[1][0],"standard diff element [1][0]");
	failcount += assertEqual(diff[1][1],m1[1][1]-m2[1][1],"standard diff element [1][1]");
	failcount += assertEqual(diff[1][2],m1[1][2]-m2[1][2],"standard diff element [1][2]");	
	Matrix diff2(2,3);
	diff2+=m1;
	diff2-=m2;
	failcount += assertEqual(diff2,diff,"standard diff 2");
	
	// factor
	double f=rdist(gen);
	failcount += assertEqual(2*m1,m1*2,"factor 2");	
	failcount += assertEqual(-m2+m1,diff,"factor -1");	
	failcount += assertEqual(m2,-(-(m2)),"factor -1*-1");
	failcount += assertEqual(f*(m1-m2)[0][0],(f*m1-f*m2)[0][0],0.000001,"algebra [0][0]");
	failcount += assertEqual(f*(m1-m2)[0][1],(f*m1-f*m2)[0][1],0.000001,"algebra [0][1]");
	failcount += assertEqual(f*(m1-m2)[0][2],(f*m1-f*m2)[0][2],0.000001,"algebra [0][2]");
	failcount += assertEqual(f*(m1-m2)[1][0],(f*m1-f*m2)[1][0],0.000001,"algebra [1][0]");
	failcount += assertEqual(f*(m1-m2)[1][1],(f*m1-f*m2)[1][1],0.000001,"algebra [1][1]");
	failcount += assertEqual(f*(m1-m2)[1][2],(f*m1-f*m2)[1][2],0.000001,"algebra [1][2]");					
	failcount += assertEqual(f*2*m1[0][0],f*(m1+m1+m2-m2)[0][0],0.000001,"algebra2 [0][0]");
	failcount += assertEqual(f*2*m1[0][1],f*(m1+m1+m2-m2)[0][1],0.000001,"algebra2 [0][1]");
	failcount += assertEqual(f*2*m1[0][2],f*(m1+m1+m2-m2)[0][2],0.000001,"algebra2 [0][2]");
	failcount += assertEqual(f*2*m1[1][0],f*(m1+m1+m2-m2)[1][0],0.000001,"algebra2 [1][0]");
	failcount += assertEqual(f*2*m1[1][1],f*(m1+m1+m2-m2)[1][1],0.000001,"algebra2 [1][1]");
	failcount += assertEqual(f*2*m1[1][2],f*(m1+m1+m2-m2)[1][2],0.000001,"algebra2 [1][2]");
	
	if( verbose ){
		cout << m1 << endl;	
		cout << m2 << endl;			
		cout << sum << endl;	
		cout << diff << endl;			
		cout << diff2 << endl;					
		cout << f << endl;				
	}		
	
	return failcount;
}



unsigned int test_matrix_multiplication(bool verbose)
{
	cout << "\033[1;36mRun Test: test_matrix_multiplication\033[0m" << endl;

	default_random_engine gen(chrono::high_resolution_clock::now().time_since_epoch().count());
    uniform_int_distribution<double> rdist(-9,9);
    uniform_int_distribution<double> idist(1,20);    
	unsigned int failcount = 0;

	Matrix m1(idist(gen),idist(gen));
	Matrix m2(m1.cols(),m1.rows());
	for(int i=0; i<m1.rows(); i++){
		for(int j=0; j<m1.cols(); j++){
			m1[i][j] = rdist(gen);
			m2[j][i] = rdist(gen);
		}
	}	
	Matrix m = m1*m2;
	Matrix n(m1.rows(),m2.cols());
	
	// matrix multiplication
	failcount += assertEqual(n.cols(),m.cols(),"matrix multiplication cols");	
	failcount += assertEqual(n.rows(),m.rows(),"matrix multiplication rows");		

	for(int i=0; i<m1.rows(); i++){
		for(int j=0; j<m1.cols(); j++){
			for(int r=0; r<m2.cols(); r++){	
				n[i][r] += m1[i][j] * m2[j][r];
			}
		}
	}		
	
	for(int i=0; i<m.rows(); i++){
		for(int j=0; j<m.cols(); j++){
			failcount += assertEqual(m[i][j],n[i][j],0.000001,"matrix multiplication element");
		}
	}		

	// another matrix multiplication	
	Matrix a(2,3);
	Matrix b(3,2);
	Matrix c(2,2);
	a[0][0] = 3;
	a[0][1] = 2;
	a[0][2] = 1;
	a[1][0] = 1;
	a[1][1] = 0;
	a[1][2] = 2;
	
	b[0][0] = 1;
	b[0][1] = 2;
	b[1][0] = 0;
	b[1][1] = 1;
	b[2][0] = 4;
	b[2][1] = 0;
	
	c[0][0] = 7;
	c[0][1] = 8;
	c[1][0] = 9;
	c[1][1] = 2;

	failcount += assertEqual(a*b,c,"matrix multiplication wikipedia");		
	failcount += assertEqual(2*(a*b),(a)*(2*b),"matrix multiplication and braces");	
	failcount += assertEqual((a*b)+(a*b),(2*a*b),"matrix multiplication and braces and +");	

	if( verbose ){
		cout << m1 << endl;	
		cout << m2 << endl;				
		cout << m << endl;			
		cout << n << endl;					
	}		
	
	return failcount;
}



int main(int argc, char* argv[])
{
	double failcount = 0;
	double verbose = false;
	
	if( argc > 1 ){
		if ( string(argv[1]) == string("-v") ){
			verbose = true;
		}
	}

	try {
		failcount+=test_vector_constructors(verbose);	
		failcount+=test_vector_operators(verbose);	
		failcount+=test_vector_angle_and_cross(verbose);
		failcount+=test_matrix_constructors(verbose);
		failcount+=test_matrix_operators(verbose);	
		failcount+=test_matrix_multiplication(verbose);
	} catch(Exception& e){ e.print(); }
	
	
	if(failcount==0){
		cout << "\033[1;32mAll tests passed!\033[0m\t\t" << endl;
	} else {
		cout << "\033[1;31m" << failcount << " tests failed\033[0m\t\t" << endl;
	}

	return 0;
}
