#include <iostream>
#include <random>
#include <string>
using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Please input a positive interger" << std::endl;
        exit(1);
    }
    int k = stoi(argv[1]);
    double* a = new double[2 * k + 1];
    mt19937 e(time(nullptr));
    uniform_real_distribution<double> urd(0, 1);
    for (int i = 0; i <= 2 * k; i++) a[i] = urd(e);

    double* b = a + k;

    for (int i = -k; i <= k; i++) {
        cout << "a[" << i + k << "] = " << a[i + k] << ", b[" << i << "] = " << b[i]
             << endl;
    }
    delete[] a;
}