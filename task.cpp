#include<iostream>
#include<fstream>
#include<cmath>
#include <ctime>
using namespace std;
int main() {

	const int n = 101;
	int iter = 0;
	double ro = 0.9, K = 1.0, alfa = 0.05, betta = 1.0;
	double dx = 1.0 / (n - 1), dy = 1.0 / (n - 1), dt = dx*dx*0.25, eps = pow(10, -8), dif, dif2;

	double gamma = 20, s0 = 22, a0 = 440;

	double A[n][n], A0[n][n];

	double S[n][n], S0[n][n];

	double g[n][n], f[n][n], F[n][n];
	srand(time(0));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {

			A0[i][j] = rand() % 440;
			S0[i][j] = rand() % 22;

		}
	}
	double a = 50.0, s = 80.0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			F[i][j] = S0[i][j] * A0[i][j] / (1.0 + S0[i][j] + K*S0[i][j] * S0[i][j]);
			g[i][j] = gamma*(s0 - S0[i][j] - ro*F[i][j]);
			f[i][j] = gamma*(alfa*(a0 - A0[i][j]) - ro*F[i][j]);

		}
	}

	do {
		for (int i = 0; i <n; i++) {
			for (int j = 0; j <n; j++) {
				A0[0][j] = A0[1][j];

				A0[n - 1][j] = A0[n-2][j];

				A0[i][n - 1] = A0[i][n-2];

				A0[i][0] = A0[i][1];

				///////////////////////////////

				S0[0][j] = S0[1][j];

				S0[n - 1][j] = S0[n - 2][j];

				S0[i][n - 1] = S0[i][n - 2];

				S0[i][0] = S0[i][1];

			}
		}

		for (int i = 0; i <= n - 1; i++) {
			for (int j = 0; j <= n - 1; j++) {

				S[i][j] = (g[i][j] * (dx*dx) + ((S0[i + 1][j] - 2.0*S0[i][j] + S0[i - 1][j])) + ((S0[i][j + 1] - 2.0*S0[i][j] + S0[i][j - 1])))*dt / (dx*dx) + S0[i][j];
				A[i][j] = (f[i][j] * (dx*dx) + betta*(((A0[i + 1][j] - 2.0*A0[i][j] + A0[i - 1][j])) + ((A0[i][j + 1] - 2.0*A0[i][j] + A0[i][j - 1]))))*dt / (dx*dx) + A0[i][j];

			}
		}


		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				A0[i][j] = A[i][j];
				S0[i][j] = S[i][j];


			}
		}
		iter++;

	} while (iter < 100);

	fstream fout("task4.dat", ios::out);
	fout << "VARIABLES = \"X\",\"Y\",\"A\",\"S\"" << endl;
	fout << "ZONE I=" << n - 2 << ",J=" << n - 2 << ",F=POINT" << endl;
	for (int i = 1; i < n - 1; i++) {
		for (int j = 1; j < n - 1; j++) {
			fout << i * dx << '\t' << j * dy << '\t' << A0[i][j] << '\t' << S0[i][j] << endl;
		}
	}


	cout << "iterations = " << iter << endl;
	system("pause");
	return 0;
}
