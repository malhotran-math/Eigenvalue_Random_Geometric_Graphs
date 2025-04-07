#include <iostream>
#include <conio.h>
#include <random>
#include<math.h>
#include<chrono>
#include <Eigen/Dense>
#include<time.h>
#define pi atan(1) * 4       //defining pi as arctan(1)*1 , as tan(pi/4)=1
using namespace std;
using namespace Eigen;
void clear_up(MatrixXd A, int dim)
{
    for(int i=0;i<dim;i++)
    {   for(int j=0;j<dim;j++)
        {   A(i,j)=0;   }   }
}
int sumrow(MatrixXd A,int row, int dim)   //adding up elements of a particular row
{
    int sum=0;
    for(int i=0;i<dim;i++)
    {
        sum=sum+A(row,i);
    }
    return sum;
}
MatrixXd generate_D(MatrixXd A, int dim)   //directly computing D^(-1/2)
{
    MatrixXd D(dim,dim);
    clear_up(D,dim);
    for(int x=0;x<dim;x++)
    {
        for(int y=0;y<dim;y++)
        {
            if(x==y)
            {
                D(x,y)=1/sqrt(sumrow(A,x,dim));  /*diagonal elements are one over square root of sum of
                                                    corresponding row elements of A */
            }
            else
                { D(x,y) =0;  }
        }
    }
    return D;
}
MatrixXd sqmatrix_multiply(MatrixXd A,MatrixXd B,int dim)
{
    MatrixXd C(dim,dim);
    clear_up(C,dim);
    for(int i=0;i<dim;i++)
    {
        for(int j=0;j<dim;j++)
        {
            for(int k=0;k<dim;k++)
            {
                C(i,j) = C(i,j) + (A(i,k)*B(k,j));
            }
        }
    }
    return C;
}
MatrixXd generate_B(MatrixXd A, int dim)
{
    MatrixXd B(dim,dim), D(dim,dim);
    MatrixXd inter1(dim,dim), inter2(dim,dim);
    D=generate_D(A,dim);
    inter1=sqmatrix_multiply(A,D,dim);  //multiplying A and D^(-1/2)
    inter2=sqmatrix_multiply(D,inter1,dim); //multiplying D^(-1/2) and the  previous product 
    for(int i=0;i<dim;i++)
    {
        for(int j=0;j<dim;j++)
        {
            if(i==j)
            {
                    B(i,j) = 1 - inter2(i,j);  //diagonal element of identity matrix is 1
            }
            else
            {
                B(i,j) = 0 - inter2(i,j); //non diagonal elements of identity matrix are 0
            }

        }
    }
    return B;
}
int main()
{
    srand(time(NULL));  //randomizing the seed
    double X,Y,x,y,a,b;
    MatrixXd L(5000,5000),M(5000,5000);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();    //setting a seed for Poisson distribution
    default_random_engine g (seed);
    int N, I;
    cout<<"Enter your number and iterations";
    cin>>N>>I;
    int P[5000];
    poisson_distribution<int> distribution(N);  //using the Poisson operator to make a Poisson distribution
    for(int i=0; i<I; i++)
    {
        P[i]=distribution(g);       //Putting Poisson distribution values in an Array
    }
    for(int k=0; k<I; k++)
    {
        for (int j=0;j<P[k];j++)
        {
          do

            {   X=(float)rand()/RAND_MAX;    //X takes float type random values in [0,1]
                Y=(float)rand()/RAND_MAX;
                x=sqrt(X);    /*we take the root of the random number generated. This follows from the
                        distribution of distance given as F(r) = r^2/R^2 (R is the radius, which in our
                        problem is 1 => F(r) = r^2 . Say u = F(r) => r = sqrt(u). We generated u (which is X)
                        and thus found r (which is x) as root of u. */
                y = 2*(pi)*Y;        /* the distribution for the angle is given as F(a)=a/2*pi . If we put b = F(a), we get
                               a = 2*pi*b. We generate b randomly (which is Y) and obtain a from it (which is y). */
                a=x*cos(y);        // X polar coordinate is rcos(theta) where r=distance and theta = angle from X axis
                b=x*sin(y);
            }while((x*x)+(y*y)<1);
           L(k,j)=a; M(k,j)=b;  //storing X and Y values in L and M for each k
        }
    }
    for(int a=0;a<I;a++)
    {
        MatrixXd A(P[a],P[a]), B(P[a],P[a]);
        clear_up(A,P[a]);   clear_up(B,P[a]);
        for(int a1=0;a1<P[a];a1++)
        {
            for(int a2=0;a2<P[a];a2++)
            {
                if(a1==a2)
                {
                    A(a1,a2)=0;     //no edge between a point and itself
                }
                else
                {
                    if(((L(a,a1)-L(a,a2))*(L(a,a1)-L(a,a2)))+((M(a,a1)-M(a,a2))*(M(a,a1)-M(a,a2)))< 1)  //edge exists
                    {
                        A(a1,a2)=1;
                    }
                    else        //edge does not exist
                    {
                        A(a1,a2)=0;
                    }
                }
            }
        }
        B=generate_B(A,P[a]);   //generating B for a matrix A with P[a] number of points
        EigenSolver<MatrixXd> b(B,false);  //"false" means we obtain only eigen values and not eigen vectors
        cout<<a+1<<" Eigen values are: \n "<<b.eigenvalues()<<"\n"; /*calculating eigen values for B as (a,b) where
			a is the real part and b is coeffictient of iota */
    }
    return 0;

}

code.txt
