#include <iostream>
#include <algorithm> // min()

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
const int n[C + 1] = {DumVal, N1, N2};
const int N = n[1] + n[2]; // αριθμός εργασιών (k)

double Q[M + 1][C + 1][N1 + 1][N2 + 1]; // see Bard-Schweitzer approximation (algorithm 5.2)
double X[C + 1];
const double D[M + 1][C + 1] = {
    {DumVal, DumVal, DumVal},
    {DumVal, 25.0, 29.0},
    {DumVal, 0.032, 0.0},
    {DumVal, 0.0, 0.090},
    {DumVal, 0.048, 0.134},
    {DumVal, 0.025, 0.0},
    {DumVal, 0.0, 0.016},
    {DumVal, 0.059, 0.0},
    {DumVal, 0.070, 0.0},
    {DumVal, 0.0, 0.28},
    {DumVal, 0.0, 0.025},
    {DumVal, 0.048, 0.058},
    {DumVal, 0.054, 0.066},
    {DumVal, 0.069, 0.106}, // Dij(1) (never explicitly used)
    {DumVal, 0.067, 0.072},
    {DumVal, 0.088, 0.096}};
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
double p[N + 1][N1 + 1][N2 + 1];
double R[M + 1][C + 1];
double U[M + 1][C + 1];
double u[M + 1][C + 1];

int main()
{
  /* Initialization of D_13jk */
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

  /* Initialization of p(k|N) */
  for (int k = 1; k <= N; ++k)
  {
    int n1 = (k >= N2 ? k - N2 : 0),
        n2 = (k >= N2 ? N2 : k);
    for (; n1 <= N1 && n2 >= 0; ++n1, --n2)
    {
      double sum = 0.0;
      for (int j = 1; j <= C; ++j)
      {
        int pos1, pos2;
        if (n1 == 0 && n2 != 0)
        {
          pos1 = 0;
          pos2 = n2 - 1;
        }
        else if (n1 != 0 && n2 == 0)
        {
          pos1 = n1 - 1;
          pos2 = 0;
        }
        else if (n1 == 0 && n2 == 0)
        {
          pos1 = 0;
          pos2 = 0;
        }
        else
        {
          if (j == 1)
          {
            pos1 = n1 - 1;
            pos2 = n2;
          }
          else
          {
            pos1 = n1;
            pos2 = n2 - 1;
          }
        }
        sum += X[j] / m[j][k] * p[k - 1][pos1][pos2];
      }
      p[k][n1][n2] = sum;
    }
  }
  double sum = 0.0;
  for (int l = 1; l <= N; ++l)
    sum += p[l][N1][N2];
  p[0][N1][N2] = 1 - sum;

  /* MVA: Main algorithm */

  /* Initialize Qij */
  for (int i = 1; i <= M; ++i)
    if (type_of_station[i] == LI)
      for (int j = 1; j <= C; ++j)
        Q[i][j][0][0] = 0.0;

  for (int iteration = 0; iteration < 1; ++iteration)
  {
    int k[3];
    for (k[1] = 0; k[1] <= n[1]; ++k[1])
      for (k[2] = (k[1] ? 0 : 1); k[2] <= n[2]; ++k[2])
      {
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
              for (int l = 1, pos1, pos2; l <= C; ++l)
              {
                pos1 = (j == 1 && k[1] > 0 ? k[1] - 1 : k[1]),
                pos2 = (j == 2 && k[2] > 0 ? k[2] - 1 : k[2]);
                sum += Q[i][l][pos1][pos2];
              }
              R[i][j] = D[i][j] * (1.0 + sum);
              break;
            case LD: // TODO: i think that this is false (see 4.56)
              // This can be calculated only if uij is known, which
              // in our case it isn't. So, we change to the approximating
              // techniques for calculating Rij. Goodbyee! ^_^
              // (Change branch :) )
              break;
            }
          }

        // X
        for (int j = 1; j <= C; ++j)
        {
          double sum = 0.0;
          for (int i = 1; i <= M; ++i)
            sum += R[i][j];
          X[j] = k[j] / sum;
        }

        // Q
        for (int i = 1; i <= M; ++i)
          for (int j = 1; j <= C; ++j)
            Q[i][j][k[1]][k[2]] = X[j] * R[i][j];
      }
  }

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
      printf("\tQ_%d_%d: %f\n", i, j, Q[i][j][n[1]][n[2]]);

  printf("Αριθμός εργασιών Qi (sum):\n");
  double Q_sum[M + 1];
  for (int i = 1; i <= M; ++i)
  {
    Q_sum[i] = 0.0;
    for (int j = 1; j <= C; ++j)
      Q_sum[i] += Q[i][j][n[1]][n[2]];
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
