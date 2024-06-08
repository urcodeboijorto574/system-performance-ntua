#include <iostream>
#include <algorithm>
// #include <vector>

using namespace std;

enum TypeOfStation
{
  DELAY,
  LD,
  LI
};

const int n[2 + 1] = {0, 304, 240};
const int N = n[1] + n[2]; // αριθμός εργασιών (k)

const int M = 15;       // αριθμός σταθμών (i)
const int C = 2;        // αριθμός κατηγοριών (j)
double Q[M + 1][C + 1]; // see Bard-Schweitzer approximation (algorithm 5.2)
double X[C + 1];
const double D[M + 1][C + 1] = {
    {0, 0, 0},
    {0, 25.0, 29.0},
    {0, 0.032, 0.0},
    {0, 0.0, 0.090},
    {0, 0.048, 0.134},
    {0, 0.025, 0.0},
    {0, 0.0, 0.016},
    {0, 0.059, 0.0},
    {0, 0.070, 0.0},
    {0, 0.0, 0.28},
    {0, 0.0, 0.025},
    {0, 0.048, 0.058},
    {0, 0.054, 0.066},
    {0, 0.069, 0.106}, /* TODO: different */
    {0, 0.067, 0.072},
    {0, 0.088, 0.096}};
const TypeOfStation type_of_station[M + 1] = {/*DummyValue*/ DELAY,
                                              /* 1 */ DELAY,
                                              /* 2 */ LI,
                                              /* 3 */ LI,
                                              /* 4 */ LI,
                                              /* 5 */ DELAY,
                                              /* 6 */ DELAY,
                                              /* 7 */ LI,
                                              /* 8 */ LI,
                                              /* 9 */ LI,
                                              /* 10 */ LI,
                                              /* 11 */ LI,
                                              /* 12 */ LI,
                                              /* 13 */ LD,
                                              /* 14 */ LI,
                                              /* 15 */ LI};
double a[N + 1];
double p[N + 1];
double R[M + 1][C + 1];
double U[M + 1];

double maxD(int j)
{
  return D[13][j];
}

double sumD(int j)
{
  double sum = 0.0;
  for (int i = 1; i <= M; ++i)
    sum += D[i][j];
  return sum;
}

void calc_p(int i = 13)
{
  double product[N + 1];

  // Calculating p13(0|N)

  double sum = 0;
  for (int k = 1; k <= N; ++k)
  {
    product[k] = 1;
    for (int l = 1; l <= k; ++l)
    {
      double sum2 = 0;
      for (int j = 1; j <= C; ++j)
      {
        sum2 += (X[j] * D[i][j]);
      }
      product[k] *= sum2 / a[l];
    }
    sum += product[k];
  }

  p[0] = 1.0 / (1.0 + sum);

  // Calculating p13(k|N)

  for (int k = 1; k <= N; ++k)
  {
    p[k] = p[0] * product[k];
  }
}

int main()
{

  /* Initialize Qij */
  for (int j = 1; j <= C; ++j)
  {
    for (int i = 1; i <= M; ++i)
    {
      if (type_of_station[i] == LI)
        Q[i][j] = n[j] / M;
    }
  }

  /* Initialize Xij (for i=13 only) */
  for (int j = 1; j <= C; ++j)
  {
    X[j] = min(1 / maxD(j), n[j] / sumD(j));
  }

  /* Initialization of a */
  for (int k = 1; k <= N; ++k)
  {
    a[k] = (k <= 64) ? 0.40 + 0.60 * k : 38.80;
  }

  /* Main algorithm */

  for (int iteration = 0; iteration < 100; ++iteration)
  {
    calc_p();
    for (int k = 1; k <= N; ++k)
    {
      for (int i = 1; i <= M; ++i)
      {
        for (int j = 1; j <= C; ++j)
        {
          double sum = 0;
          switch (type_of_station[i])
          {
          case DELAY:
            R[i][j] = D[j][i];
            break;
          case LI:
            for (int l = 1; l <= C; ++l)
              sum += Q[i][l];
            R[i][j] = D[j][i] * sum;
            break;
          case LD:
            for (int j = 1; j <= C; ++j)
            {
              for (int k = 1; k <= N; ++k)
              {
                sum += k * p[k - 1] / a[k];
              }
              R[i][j] = D[j][i] * sum;
            }
            break;
          }
        }
      }

      // printf("X[j]:\n");
      for (int j = 1; j <= C; ++j)
      {
        double sum = 0.0;
        for (int i = 1; i <= M; ++i)
        {
          sum += R[i][j];
        }
        X[j] = k / sum;
        // printf("\tX[%d]: %f\n", j, X[j]);
      }

      for (int i = 1; i <= M; ++i)
      {
        for (int j = 1; j <= C; ++j)
        {
          Q[i][j] = X[j] * R[i][j];
        }
      }
    }
  }

  for (int j = 1; j <= C; ++j)
  {
    for (int i = 1; i <= M; ++i)
    {
      U[i] = X[j] * D[i][j];
    }
  }

  /* TODO: Print results */
  // Ζητούμενα: X,R,Q,U ανά κατηγορία και συνολικά

  printf("Ρυθμός απόδοσης X ανά κατηγορία:\n\tX_1: %f\n\tX_2: %f\n", X[1], X[2]);
  printf("Συνολικός ρυθμός απόδοσης X: %f\n", X[1] + X[2]);

  double R_total, R_j1, R_j2;
  for (int j = 1; j <= C; ++j)
  {
    double sum[j] = {0.0, 0.0};
    for (int i = 1; i <= M; ++i)
    {
      sum[j] += R[i][j];
    }
    if (j == 1)
      R_j1 = sum[j];
    else
      R_j2 = sum[j];
  }
  R_total = R_j1 + R_j2;

  printf("Χρόνος απόκρισης ανά κατηγορία:\n\tR_1: %f\n\tR_2: %f\n", R_j1, R_j2);
  printf("Συνολικός χρόνος απόκρισης συστήματος R: %f\n", R_total);

  printf("Μέσος αριθμός εργασιών ανά σταθμό Q_i:\n");
  double Q_mean[M + 1];
  for (int i = 1; i <= M; ++i)
  {
    Q_mean[i] = (Q[i][1] + Q[i][2]) / 2;
    printf("\tQ_%d: %f\n", i, Q_mean[i]);
  }

  printf("Βαθμός χρησιμοποίησης ανά σταθμό U_j:\n");
  for (int i = 1; i <= M; ++i)
  {
    printf("\tU_%d: %f\n", i, U[i]);
  }
}
