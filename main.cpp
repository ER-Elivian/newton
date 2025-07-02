#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>


using namespace std;

//Вероятно, стоит разделить это для упрощения чтения на класс-файлы (но одним куском проще навигация)
class Matrix
{
  friend class SLAU_solve_Gauss;
public:
  Matrix(vector<vector<double> > i_matrix)
  {
    matrix=i_matrix;
  }

  Matrix(vector<double> i_vector)
  {//Сразу сделаем из него столбец
    matrix=vector<vector<double> >(1,vector<double>(1));
    for(int i=0;i<i_vector.size()-1;++i)
      matrix.push_back(vector<double>(1));

    for(int i=0;i<i_vector.size();++i)
      matrix[i][0]=i_vector[i];
  }

  matadd(Matrix mat_A)
  {
    if(matrix.size()!=mat_A.matrix.size() || matrix[0].size()!=mat_A.matrix[0].size())
      throw "Matrices have different shape";
    for(int i=0;i<matrix.size();++i)
      for(int j=0;j<matrix[0].size();++j)
        matrix[i][j]+=mat_A.matrix[i][j];
  }

  matmul(Matrix mat_A)
  {
    if(matrix[0].size()!=mat_A.matrix.size())
      throw "1st width != 2nd height";

    vector < vector<double> > res_mat(matrix.size(), vector<double>(mat_A.matrix[0].size()));
    for(int i=0;i<matrix.size();++i)
      for(int j=0;j<mat_A.matrix[0].size();++j)
        for(int k=0;k<matrix[0].size();++k)
        {
          res_mat[i][j]+=matrix[i][k]*mat_A.matrix[k][j];
        }
    matrix=res_mat;
  }

private:
  vector<vector<double> > matrix;
};

class Eq
  {
public:
  Eq(vector<vector<double> > i_polynom) // состоит из векторов мономов: coef,pow_x1,pow_x2...

  {
    polynom=i_polynom;
  }

  double calc(vector<double> arg)
  {
    double res=0;

    for(int i=0; i<polynom.size(); ++i)
    {
      double mon_mul=polynom[i][0];
      for(int j=1; j<polynom[i].size(); ++j)
      {
        mon_mul*=pow(arg[j-1],polynom[i][j]);
      }
      res+=mon_mul;
    }
    return res;
  }

  double deriv(vector <double> arg, int num, double eps)//можно было сделать и производную по направлению, но тут нужен будет всего лишь якобиан
  {
    vector <double> arg2=arg;
    arg2[num]+=eps;
    return (calc(arg2)-calc(arg))/eps;//Возможно, стоит использовать  двустороннюю производную
  }

//private:
  vector<vector<double> > polynom;
  };

class SLAU_solve_Gauss //обращать матрицу - ещё хуже
{
public:
  SLAU_solve_Gauss(Matrix i_m, vector<double> i_v)
  {
  m=i_m.matrix;
  for(int i=0;i<m.size();++i)
    m[i].push_back(i_v[i]);
  }

  vector<double> solve(void)
  {
  //forvard
    double thresh_val=0.0001;
    for(int i = 0;i<m.size()-1;++i)
    {
      if(m[i][i]<thresh_val && m[i][i]> -thresh_val)//на диагонали ноль
      {
        int j = i+1;
        for(; m[i][j]<thresh_val && m[i][j]> -thresh_val;++j)//ищем не ноль ниже
        {
          if(j>=m.size())
            throw "SLAU unsolveable";//не нашли
        }
        addc(i,j,1);//нашли
      }
      for(int j = i+1; j<m.size();++j)
      {
        double rel=m[j][i]/m[i][i];
        addc(j,i,-rel);
      }
    }
  //reverse
    for(int i = m.size()-1; i>=0;--i)
    {
    //на диагонали нуля быть уже не может
      for(int j = i-1; j>=0;--j)
      {
        double rel=m[j][i]/m[i][i];
        addc(j,i,-rel);
        }
      mul_line(i,1/m[i][i]);
    }
  //result
    vector<double> res(0);
    for(int i = 0; i < m.size(); ++i)
      res.push_back(m[i][m[i].size()-1]);
    return res; //возможно, стоит передавать указатель на вектор, но эта функция всё равно одноразовая, так что потом надо будет копировать
  }

  void printmatrix(void)
  {
  for(int i=0;i<m.size();++i)
    {
    for(int j=0;j<m[0].size();++j)
      printf("%lf \t",m[i][j]);
    printf("\n");
    }
  printf("\n");
  }

private:
  void mul_line (int num, double coef)
  {
    if(num>=m.size())
      throw "num > matrix height";

    for(int i=0;i<m[0].size();++i)
    {
      m[num][i]*=coef;
    }
  }

  void addc(int num_a, int num_b, double coef) //a=a+coef*b, a,b - номера строк. Эта операция только в Гауссе нужна, а не в матрицах
  {
    if(num_a>=m.size())
      throw "num_a > matrix height";
    if(num_b>=m.size())
      throw "num_b > matrix height";
    if(num_a==num_b)
      throw "num_a = num_b";

    for(int i=0;i<m[0].size();++i)
    {
      m[num_a][i]+=coef*m[num_b][i];
    }
  }

  vector<vector<double> > m;
};

class Newton
{
public:
  Newton(vector<Eq> i_to_find, vector<double> i_start_point, double i_eps, int i_maxstep)
  {
    to_find=i_to_find;
    last_point=i_start_point;
    eps=i_eps;//погрешность
    maxstep=i_maxstep;
  }

  vector<double> solve(void)
  {
  calc_F();
  F_norm_val=calc_norm_val(F);
  printf ("dF %lf\n", F_norm_val);
  while (F_norm_val>eps)
    {
    make_step();
    F_norm_val=calc_norm_val(F);
    d_x_norm_val=calc_norm_val(delta_x);
    printf ("dx %lf\n", d_x_norm_val);
    printf ("dF %lf\n", F_norm_val);

    if(step>maxstep)
      throw "solution point unrechable";
    }

    return last_point;
  }

  void make_step(void)
  {
    calc_F();
    calc_J();
    Matrix m(J);
    SLAU_solve_Gauss d_x(m,F);
    delta_x=d_x.solve();
    xadd(delta_x);
    step+=1;
  }

  void printJ(void)
  {
  printf("J:\n");
  for(int i=0;i<J.size();++i)
    {
    for(int j=0;j<J[0].size();++j)
      printf("%lf \t",J[i][j]);
    printf("\n");
    }
  printf("\n");
  }

private:
  double calc_norm_val(vector<double> v)
  {
    double res;
    for(int i=0;i<v.size();++i)
      res+=v[i]*v[i];

    return sqrt(res);
  }

  void xadd(vector<double> delta_x)
  {
    if(delta_x.size()!=last_point.size())
      throw "vectors have different dim";
    for(int i=0;i<last_point.size();++i)
      last_point[i]-=delta_x[i];
  }

  void calc_F(void)
  {
    F=vector<double>(last_point.size());
    for(int i=0;i<last_point.size();++i)
      F[i]=to_find[i].calc(last_point);
  }

  void calc_J(void)
  {
    J=vector<vector<double> >(to_find.size(), vector<double>(last_point.size()));

    for(int i=0; i<to_find.size();++i)
    {
      for(int j=0; j<last_point.size();++j)
      {
        J[i][j]=to_find[i].deriv(last_point,j,0.0001);//! Это стоит делать как-то более оптимально, используя значение функции в текущей точке
      }
    }
  }

  vector<double> delta_x;
  vector<Eq> to_find;
  vector<double> last_point;
  vector<vector<double> > J; //якобиан
  vector<double> F;// F = to_find(last_point)
  double d_x_norm_val;
  double F_norm_val;
  double eps;
  int step=0,maxstep;
};

int main()
{
  double coefs[6];
  ifstream c_file("coefs.txt"); // окрываем файл для чтения
    if (c_file.is_open())
    {
    for(int i=0;i<6;++i)//6 чисел в файле
      c_file >> coefs[i];
    }


  double eq_arr1[3][3]={{1,3,0},{1,0,2},{1,0,0}};
  eq_arr1[0][0]=coefs[0];
  eq_arr1[1][0]=coefs[1];
  eq_arr1[2][0]=-coefs[2];
  vector<vector<double> > Eq1_v(3,(vector<double>(3)));

  for(int i=0;i<3;++i)
    for(int j=0;j<3;++j)
      Eq1_v[i][j]=eq_arr1[i][j];

  double eq_arr2[3][3]={{1,1,3},{1,0,1},{1,0,0}};
  eq_arr2[0][0]=coefs[3];
  eq_arr2[1][0]=coefs[4];
  eq_arr2[2][0]=-coefs[5];
  vector<vector<double> > Eq2_v(3,(vector<double>(3)));

  for(int i=0;i<3;++i)
    for(int j=0;j<3;++j)
      Eq2_v[i][j]=eq_arr2[i][j];

  Eq Eq1(Eq1_v);
  Eq Eq2(Eq2_v);

  vector<Eq> to_find;
  to_find.push_back(Eq1);
  to_find.push_back(Eq2);

  vector<double> start_point;
  start_point.push_back(-3);
  start_point.push_back(-7);

  double eps=0.001;

  Newton sys(to_find,start_point,eps,1000);
  vector<double> solv=sys.solve();

  printf("x_vector:\n");
  for(int i=0;i<solv.size();++i)
  {
    printf("%lf ",solv[i]);
  }
  return 0;
}


/*
вычисление "полиноминальной" функции и производных

#define NUM_VARS 2
#define NUM_MONOMS 3

  vector<vector<double> > polynom(NUM_MONOMS,vector<double>(NUM_VARS+1));
  double polynom_arr[NUM_MONOMS][NUM_VARS+1]={{2,2,0},{2,0,2},{3,2,2}};

  for(int i=0;i<NUM_MONOMS;++i)
    for(int j=0;j<NUM_VARS+1;++j)
      polynom[i][j]=polynom_arr[i][j];

  Eq eq1(polynom);
  vector<double> p1 (2);
  p1[0]=0,p1[1]=0;
  vector<double> p2 (2);
  p2[0]=1,p2[1]=0;
  vector<double> p3 (2);
  p3[0]=0,p3[1]=1;
  vector<double> p4 (2);
  p4[0]=1,p4[1]=1;


  printf("%lf\n",eq1.calc(p1));
  printf("%lf\n",eq1.calc(p2));
  printf("%lf\n",eq1.calc(p3));
  printf("%lf\n",eq1.calc(p4));

  printf("%lf\n",eq1.deriv(p1,0,0.001));
  printf("%lf\n",eq1.deriv(p2,0,0.001));
  printf("%lf\n",eq1.deriv(p3,1,0.001));
  printf("%lf\n",eq1.deriv(p4,0,0.001));

*** для перемножения матриц

#define NUM_ROWS 3
#define NUM_COLUMNS 3

  vector<vector<double> > matrix1(NUM_ROWS,vector<double>(NUM_COLUMNS));
  double matrix_arr1[NUM_ROWS][NUM_COLUMNS]={{1,1,1},{0,1,1},{0,0,1}};

  for(int i=0;i<NUM_ROWS;++i)
    for(int j=0;j<NUM_COLUMNS;++j)
      matrix1[i][j]=matrix_arr1[i][j];

#define NUM 3

  vector<double> matrix2(vector<double>(NUM));
  double matrix_arr2[NUM]={1,2,3};

  for(int i=0;i<NUM;++i)
    matrix2[i]=matrix_arr2[i];

  Matrix matr(matrix1);
  Matrix vec(matrix2);
  matr.matmul(vec);
  printf ("%d ", matr.matrix.size());
  printf ("%d\n", matr.matrix[0].size());
  for(int i=0;i<matr.matrix.size();++i)
    {
    for(int j=0;j<matr.matrix[0].size();++j)
      printf("%lf \t",matr.matrix[i][j]);
    printf("\n");
    }

*** для метода Гаусса

#define NUM_ROWS 3
#define NUM_COLUMNS 3

  vector<vector<double> > matrix1(NUM_ROWS,vector<double>(NUM_COLUMNS));
  double matrix_arr1[NUM_ROWS][NUM_COLUMNS]={{0,-4,-2},{2,-1,3},{1,2,-1}};

  for(int i=0;i<NUM_ROWS;++i)
    for(int j=0;j<NUM_COLUMNS;++j)
      matrix1[i][j]=matrix_arr1[i][j];

#define NUM 3

  vector<double> vec(vector<double>(NUM));
  double matrix_arr2[NUM]={-28,13,9};

  for(int i=0;i<NUM;++i)
    vec[i]=matrix_arr2[i];

  Matrix matr(matrix1);
  SLAU_solve_Gauss SLAU(matr,vec);
  vector<double>vector_solve=SLAU.solve();

  for(int i=0;i<vector_solve.size();++i)
    printf("%lf\n",vector_solve[i]);
*/
