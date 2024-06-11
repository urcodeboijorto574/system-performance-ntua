#include <iostream>
#include <algorithm> // min()
#include <cmath>     // abs()

using namespace std;

#define DumVal 0 // The 1st element of all the arrays is a dummy value
                 // so the iteration starts at 1 and ends at its size.

enum TypeOfStation
{
  DELAY,
  LD,
  LI
};

const int M = 15; // αριθμός σταθμών (i)
const int C = 2;  // αριθμός κατηγοριών (j)
const int N1 = 304, N2 = 240;
const int N = N1 + N2; // αριθμός εργασιών (k)
const int n[C + 1] = {DumVal, N1, N2};

double Q[M + 1][C + 1];     // see Bard-Schweitzer approximation (algorithm 5.2)
double Q_old[M + 1][C + 1]; // used for calculating the diversion
double X[C + 1];
const double D[M + 1][C + 1] = {
    {DumVal, DumVal, DumVal},
    /* 1 */ {DumVal, 25.0, 29.0},
    /* 2 */ {DumVal, 0.032, 0.0},
    /* 3 */ {DumVal, 0.0, 0.090},
    /* 4 */ {DumVal, 0.048, 0.134},
    /* 5 */ {DumVal, 0.025, 0.0},
    /* 6 */ {DumVal, 0.0, 0.016},
    /* 7 */ {DumVal, 0.059, 0.0},
    /* 8 */ {DumVal, 0.070, 0.0},
    /* 9 */ {DumVal, 0.0, 0.028},
    /* 10 */ {DumVal, 0.0, 0.025},
    /* 11 */ {DumVal, 0.048, 0.058},
    /* 12 */ {DumVal, 0.054, 0.066},
    /* 13 */ {DumVal, 0.069, 0.106}, // Dij(1) (never explicitly used)
    /* 14 */ {DumVal, 0.067, 0.072},
    /* 15 */ {DumVal, 0.088, 0.096}};
double D_13[C + 1][N + 1];
const TypeOfStation type_of_station[M + 1] = {/* DumVal */ DELAY,
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
double m[C + 1][N + 1];
double p[N + 1];
double R[M + 1][C + 1];
double U[M + 1][C + 1];

double maxD(int j)
{
  // TODO: check correctness of code below
  double max = D[1][j];
  for (int i = 2; i <= n[j]; ++i)
    if (max < D[i][j])
      max = D[i][j];
  return max;
}

double sumD(int j)
{ // TODO: check correctness
  double sum = 0.0;
  for (int i = 1; i <= M; ++i)
    sum += D[i][j];
  return sum;
}

void calc_p13() // depends on Xj, Dij, ak
{
  /* Calculating p13(0|N) */
  double sum = 0.0, product[N + 1];
  for (int k = 1; k <= N; ++k)
  {
    product[k] = 1;
    for (int l = 1; l <= k; ++l)
    {
      double sum_nominator = 0;
      for (int j = 1; j <= C; ++j)
        sum_nominator += (X[j] * D[13][j]); // TODO: check correctness of D
      product[k] *= sum_nominator / a[l];
    }
    sum += product[k];
  }
  p[0] = 1.0 / (1.0 + sum);

  /* Calculating p13(k|N) */
  for (int k = 1; k <= N; ++k)
    p[k] = p[0] * product[k];
}

int main()
{
  /* Initialization of D_13jk, ak, mjk */
  D_13[1][1] = D[13][1];
  D_13[2][1] = D[13][2];
  for (int k = 1; k <= N; ++k)
  {
    a[k] = (k <= 64) ? 0.40 + 0.60 * k : 38.80;
    for (int j = 1; j <= C; ++j)
    {
      D_13[j][k] = D_13[j][1] / a[k];
      m[j][k] = 1 / D_13[j][k]; // TODO: check later for validity
    }
  }

  /* Initialization of Qij (and Q_old) */
  for (int i = 1; i <= M; ++i)
    for (int j = 1; j <= C; ++j)
      Q_old[i][1] = Q[i][j] =
          type_of_station[i] != DELAY ? n[j] / M : DumVal;

  /* Initialization of Xj */
  for (int j = 1; j <= C; ++j)
    X[j] = min(1 / maxD(j), n[j] / sumD(j));

  bool is_precision_achieved;
  double precision = 0.001;
  do
  {
    is_precision_achieved = true;

    calc_p13();

    // R
    for (int i = 1; i <= M; ++i)
      for (int j = 1; j <= C; ++j)
      {
        double sum = 0.0;
        switch (type_of_station[i])
        {
        case DELAY:
          R[i][j] = D[i][j];
          break;
        case LI:
          R[i][j] = D[i][j] * (1.0 + (n[j] - 1) / n[j] * Q[i][j] + (j == 1 ? Q[i][2] : Q[i][1]));
          break;
        case LD:
          for (int k = 1; k <= N; ++k)
            sum += k / a[k] * p[k - 1];
          R[i][j] = D[i][j] * sum; // TODO: check correctness
          break;
        }
      }

    // X
    for (int j = 1; j <= C; ++j)
    {
      double sum = 0.0;
      for (int i = 1; i <= M; ++i)
        sum += R[i][j];
      X[j] = n[j] / sum;
    }

    // Q
    for (int i = 1; i <= M; ++i)
      for (int j = 1; j <= C; ++j)
        Q[i][j] = X[j] * R[i][j];

    /* Precision check */
    for (int i = 1; i <= M; ++i)
      for (int j = 1; j <= C; ++j)
        if (abs(Q[i][j] - Q_old[i][j]) > precision)
        {
          is_precision_achieved = false;
          break;
        }

  } while (is_precision_achieved);

  /* Calculation of Uij */
  for (int j = 1; j <= C; ++j)
    for (int i = 1; i <= M; ++i)
      U[i][j] = X[j] * D[i][j];

  /* Display Results */ // X,R,Q,U ανά κατηγορία και συνολικά
  printf("Ρυθμός απόδοσης Xj:\n\tX_1: %f\n\tX_2: %f\n", X[1], X[2]);
  printf("Συνολικός ρυθμός απόδοσης X: %f\n", X[1] + X[2]);

  printf("Χρόνος απόκρισης Rj:\n");
  double R_[C + 1] = {DumVal, 0.0, 0.0};
  for (int j = 1; j <= C; ++j)
  {
    for (int i = 1; i <= M; ++i)
      R_[j] += R[i][j];
    printf("\tR_%d: %f\n", j, R_[j]);
  }
  printf("Συνολικός χρόνος απόκρισης R: %f\n", R_[1] + R_[2]);

  printf("Αριθμός εργασιών Qij:\n");
  for (int j = 1; j <= C; ++j)
    for (int i = 1; i <= M; ++i)
      printf("\tQ_%d_%d: %f\n", i, j, Q[i][j]);

  printf("Αριθμός εργασιών Qi (sum):\n");
  double Q_sum[M + 1];
  for (int i = 1; i <= M; ++i)
  {
    Q_sum[i] = Q[i][1] + Q[i][2];
    printf("\tQ_%d: %f\n", i, Q_sum[i]);
  }

  printf("Βαθμός χρησιμοποίησης Uij:\n");
  for (int j = 1; j <= C; ++j)
    for (int i = 1; i <= M; ++i)
      printf("\tU_%d_%d: %f\n", i, j, U[i][j]);

  printf("Βαθμός χρησιμοποίησης Ui (sum):\n");
  double U_sum[M + 1];
  for (int i = 1; i <= M; ++i)
  {
    U_sum[i] = 0.0;
    for (int j = 1; j <= C; ++j)
      U_sum[i] += U[i][j];
    printf("\tU_%d: %f\n", i, U_sum[i]);
  }
}
