#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

// Function to calculate values
double function(double x) {
    return 2 * pow(x, 7) + 3 * pow(x, 4) + 2 * pow(x, 2) + 2;
}

// Print the system of equations (matrix A and vector b)
void printSystem(const vector<vector<double>>& A, const vector<double>& b) {
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].size(); ++j) {
            cout << setw(10) << A[i][j] << " ";
        }
        cout << "| " << setw(10) << b[i];
        cout << endl;
    }
    cout << endl;
}

void forwardElimination(vector<vector<double>>& A, vector<double>& b) {
    int n = A.size();
    for (int i = 0; i < n; ++i) {
        cout << "Normalizing row " << i + 1 << ":\n";

        double diagElement = A[i][i];
        for (int j = 0; j < n; ++j) {
            A[i][j] /= diagElement;
        }
        b[i] /= diagElement;

        printSystem(A, b);

        for (int k = i + 1; k < n; ++k) {
            double factor = A[k][i];
            cout << "Eliminating element in row " << k + 1 << ", column " << i + 1 << ":\n";
            for (int j = 0; j < n; ++j) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
            printSystem(A, b);
        }
    }
}

vector<double> backSubstitution(const vector<vector<double>>& A, const vector<double>& b) {
    int n = A.size();
    vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A[i][j] * x[j];
        }
    }
    return x;
}

int main() {
    vector<double> nodes = { 0, 1, 2, 3, 4 };
    int n = nodes.size();

    cout << "1. Define the polynomial:\n";
    cout << "L(x) = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4\n\n";

    cout << "2. Calculate L(x) values at nodes:\n";
    vector<double> values;
    for (double x : nodes) {
        double value = function(x);
        values.push_back(value);
        cout << "L(" << x << ") = " << value << endl;
    }

    vector<vector<double>> A(n, vector<double>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i][j] = pow(nodes[i], j);
        }
    }

    cout << "\n3. Form the matrix A and vector b:\n";
    printSystem(A, values);

    cout << "4. Perform forward elimination (Gauss method):\n";
    forwardElimination(A, values);

    cout << "5. Find coefficients using back substitution:\n";
    vector<double> coefficients = backSubstitution(A, values);

    cout << "\nPolynomial coefficients:\n";
    for (int i = 0; i < coefficients.size(); ++i) {
        cout << "c" << i << " = " << coefficients[i] << endl;
    }

    cout << "\nFinal polynomial:\n";
    cout << "L(x) = ";
    for (int i = 0; i < coefficients.size(); ++i) {
        if (i > 0) cout << " + ";
        cout << coefficients[i] << "*x^" << i;
    }
    cout << endl;

    return 0;
}