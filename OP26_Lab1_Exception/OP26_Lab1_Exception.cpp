#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

class ErrorNoFile
{
    string str;
public:
    ErrorNoFile(string s) : str(s) { }
    void Message()
    {
        cout << "ErrorNoFile: " << str << " - switching to Algorithm 4." << endl;
    }
};

class SignalAlg2 {};
class SignalAlg3 {};

void CreateDataFiles() {
    ofstream f1("dat_X_1_1.dat");
    f1 << -1.000 << " " << -4.935 << " " << 1.935 << endl;
    f1 << -0.500 << " " << -1.097 << " " << 2.951 << endl;
    f1 << 0.000 << " " << 0.000 << " " << 4.571 << endl;
    f1 << 0.500 << " " << 0.548 << " " << 4.596 << endl;
    f1 << 1.000 << " " << 0.000 << " " << 3.000 << endl;
    f1.close();

    ofstream f2("dat_X1_00.dat");
    f2 << 0.000 << " " << -4.935 << " " << 1.935 << endl;
    f2 << 0.500 << " " << 3.356 << " " << 1.388 << endl;
    f2 << 1.000 << " " << 5.890 << " " << 0.377 << endl;
    f2.close();

    ofstream f3("dat_X00_1.dat");
    f3 << 0.000 << " " << -4.935 << " " << 1.935 << endl;
    f3 << -0.500 << " " << -0.141 << " " << 0.976 << endl;
    f3 << -1.000 << " " << 3.480 << " " << 0.252 << endl;
    f3.close();
}

double GetTableVal(double x, int mode)
{
    string filename;
    double search_x = x;

    if (abs(x) <= 1) {
        filename = "dat_X_1_1.dat";
    }
    else if (x < -1) {
        search_x = 1.0 / x;
        filename = "dat_X00_1.dat";
    }
    else {
        search_x = 1.0 / x;
        filename = "dat_X1_00.dat";
    }

    ifstream is(filename.c_str());
    if (!is) throw ErrorNoFile(filename);

    double xi, ti, ui, xi1, ti1, ui1;
    double res = 0;
    bool found = false;

    is >> xi1 >> ti1 >> ui1;

    if (abs(xi1 - search_x) < 1e-9) {
        return (mode == 0) ? ti1 : ui1;
    }

    while (is >> xi >> ti >> ui)
    {
        double min_x = (xi1 < xi) ? xi1 : xi;
        double max_x = (xi1 < xi) ? xi : xi1;

        if (search_x >= min_x && search_x <= max_x) {
            double val_prev = (mode == 0) ? ti1 : ui1;
            double val_curr = (mode == 0) ? ti : ui;

            res = val_prev + (val_curr - val_prev) * (search_x - xi1) / (xi - xi1);
            found = true;
            break;
        }

        xi1 = xi; ti1 = ti; ui1 = ui;
    }

    is.close();

    if (!found) res = (mode == 0) ? ti1 : ui1;

    return res;
}

double Srz(double x, double y, double z)
{
    double Tx = GetTableVal(x, 0);
    double Ty = GetTableVal(y, 0);
    double Uy = GetTableVal(y, 1);
    double Uz = GetTableVal(z, 1);

    if (x > y)
        return Tx + Uz - Ty;
    else
        return Ty + Uy - Uz;
}

double Algorithm4_fun(double x, double y, double z)
{
    return 1.3498 * x + 2.2362 * y * z - 2.348 * x * y;
}

double Srs1_Alg2(double x, double y, double z)
{
    if (z > y) return Srz(x, y, z) + 1.44 * y * z;
    else return y + 1.44 * Srz(z, x, y);
}

double Qrz1_Alg2(double x, double y)
{
    if (abs(y) < 1) return x * Srs1_Alg2(x, y, x);
    else return y * Srs1_Alg2(y, x, y);
}

double Rrz_Alg2(double x, double y, double z)
{
    if (x > y) return x * y * Qrz1_Alg2(y, z);
    else return x * z * Qrz1_Alg2(x, y);
}

double Srs2_Alg3(double x, double y, double z)
{
    if (z > y) return Srz(x, y, z) + y * x;
    else return y * z + Srz(z, x, y);
}

double Qrz2_Alg3(double x, double y)
{
    if (abs(x) < 1) return x * Srs2_Alg3(x, y, x);
    else return y * Srs2_Alg3(y, x, y);
}

double Rrz_Alg3(double x, double y, double z)
{
    if (x > y) return x * y * Qrz2_Alg3(y, z);
    else return y * z * Qrz2_Alg3(x, y);
}

double Srs_Alg1(double x, double y, double z)
{
    if ((z * z + x * y) <= 0) throw SignalAlg2();
    if ((x * x + z * y) <= 0) throw SignalAlg3();

    if (z > y)
        return Srz(x, y, z) + y * sqrt(z * z + x * y);
    else
        return y + Srz(z, x, y) * sqrt(x * x + z * y);
}

double Qrz_Alg1(double x, double y)
{
    if (abs(x) < 1)
        return x * Srs_Alg1(x, y, x);
    else
        return y * Srs1_Alg2(y, x, y);
}

double Rrz_Alg1_Internal(double x, double y, double z)
{
    if (x > y)
        return x * z * Qrz_Alg1(y, z);
    else
        return y * x * x * Qrz_Alg1(x, y);
}

double Grs(double x, double y, double z)
{
    double term1;
    try {
        term1 = Rrz_Alg1_Internal(x, y, y);
    }
    catch (SignalAlg2) {
        term1 = Rrz_Alg2(x, y, y);
    }
    catch (SignalAlg3) {
        term1 = Rrz_Alg3(x, y, y);
    }

    double term2;
    try {
        term2 = Rrz_Alg1_Internal(x - y, z, y);
    }
    catch (SignalAlg2) {
        term2 = Rrz_Alg2(x - y, z, y);
    }
    catch (SignalAlg3) {
        term2 = Rrz_Alg3(x - y, z, y);
    }

    return 0.1389 * term1 + 1.8389 * term2;
}

double fun(double x, double y, double z)
{
    return x * Grs(x, y, z) + y * Grs(x, z, y);
}

int main()
{
    CreateDataFiles();

    double x, y, z, f;

    cout << "Input x y z: ";
    cin >> x >> y >> z;

    try {
        f = fun(x, y, z);
    }
    catch (ErrorNoFile& e)
    {
        e.Message();
        f = Algorithm4_fun(x, y, z);
    }
    catch (...)
    {
        cout << "\n Unknown error catch...";
        f = Algorithm4_fun(x, y, z);
    }

    cout << "\n fun = " << f << endl;

    return 0;
}