//
//  main.cpp
//  sekibun
//
//  Created by イグノタスペベレル on 2023/05/01.
//

#include <bits/stdc++.h>
//#include "MT.h"

using namespace std;
float inf = numeric_limits<float>::infinity();
using ll = long long;
const double radius = 100;
const double GAMMA = pow(10, 6.0 / 10.0);
const double end_time = 1000;
constexpr double PI = 3.14159265358979323846264338;
const double theta = pow(10, 0.4);
//const double delta = 2.0 / pass_loss_exponent;
const double lambda_BS = 1.;
const double POWER = 1.0;
const double Euler_constant = 0.5772156649;


double urand(){
    double m, a;
    m = RAND_MAX + 1.0;
    a = (rand() + 0.5)/m;
    a = (rand() + a)/m;
    return (rand() + a)/m;
}

double exp_dist(double lambda) {
    //double g = genrand_real3();
    double g = urand();
    double tau = - log(1 - g) / lambda;
    return tau;
}

double my_func(double x, double w, double alpha, double lambda) {
    double delta = 2.0 / alpha;
    double A = pow(theta, delta) * PI * delta / sin(PI * delta);
    double C = 1.0 / x - 1.0;
    double a = C / (x * x);
    double B = -w * theta * pow(C, alpha) - (1 + lambda * A) * PI * C * C;
    return  a * exp(B);
}

double func2(double y, double w, double alpha, double lambda) {
    double delta = 2.0 / alpha;
    double A = pow(theta, delta) * PI * delta / sin(PI * delta);
    double a = - log(y);
    double b = pow(sqrt(a / (1. + lambda * A) / PI), alpha);
    return exp(-w * theta * b);
}

double NOMA_L1_pr1(double phi, double k, double t, double lambda, double alpha) {
    double ans;
    ans = pow(tan(phi), alpha + 2.) * pow(tan(k), alpha + 1.) / (1. + pow(1./lambda/PI, alpha / 2.) * pow(log(1./t), alpha / 2.) * cos(k) * cos(k));
    return ans;
}

double NOMA_L1_pr2(double phi, double k, double tmp) {
    return tan(phi) / cos(k) / cos(k) - tmp;
}

double NOMA_L1_pr3(double phi, double tmp, double lambda) {
    return 2. * PI * lambda * sin(phi) / cos(phi) / cos(phi) * exp(-PI * lambda * tan(phi) * tan(phi)) * exp(-2. * PI * lambda *  tmp);
    //return exp(-PI * lambda * tan(phi) * tan(phi)) * exp(-2. * PI * lambda *  tmp);
}


int main() {
    double dx, dphi, dt, dk, n = 1000.0;
    double s = 0;
    double x1, x2, f1, f2;
    double w = 10, alpha = 4.5, lambda = 3.5;
    ofstream outputfile;
    string filename = "pr_noma_L1_theo.txt";
    outputfile.open(filename);
    
    
    
    for (lambda = 0; lambda <= 5.0; lambda += 0.1) {
        outputfile << lambda << " ";
        dphi = PI/ 2. / n; dk = PI / 4. / n; dt = 1. / n;
        s = 0;
        
        for (double i = 0; i < n-1; i++) {
            double p1 = dphi * i;
            double p2 = dphi * (i+1);
            double tmp2 = 0;
            for (double j = n-1; j >= 0; j--) {
                double k1 = dk * j;
                double k2 = dk * (j+1);
                double tmp = 0;
                for (double l = 0; l < n-1; l++) {
                    double t1 = dt * l;
                    double t2 = dt * (l+1);
                    f1 = NOMA_L1_pr1(p1, k1, t1, lambda, alpha);
                    f2 = NOMA_L1_pr1(p2, k2, t2, lambda, alpha);
                    tmp += (f1 + f2) * dt / 2.0;
                }
                //cout << "tmp = " << tmp << endl;
                double g1 = NOMA_L1_pr2(p1, k1, tmp);
                double g2 = NOMA_L1_pr2(p2, k2, tmp);
                tmp2 += (g1 + g2) * dk / 2.;
            }
            //cout << "tmp2 = " << tmp2 << endl;
            double h1 = NOMA_L1_pr3(p1, tmp2, lambda);
            double h2 = NOMA_L1_pr3(p2, tmp2, lambda);
            s += (h1 + h2) * dphi / 2.;
        }
        cout << lambda << " " << s << endl;
        outputfile << s << endl;
    }
    
    
    /*
    //outputfile << "#d= none 1000 500 200" << endl;
    for (lambda = 0; lambda <= 5.0; lambda += 0.1) {
        //w = 1.0;
        outputfile << lambda << " ";
        for (alpha = 3.5; alpha <= 4.5; alpha += 1.0) {
//        for (int i = 0; i < 4; i++) {
//            double x = 0;
//            if (i == 1) x = 1000;
//            else if (i == 2) x = 500;
//            else if (i == 3) x = 200;
//            if (x == 0) w = 0;
//            else w = pow(x / 1000., -alpha);
//            cout << alpha << endl;
            n = 1000000.0;
            dx = 1.0 / n;
            
            s = 0;
            for (double i = dx; i < n - 1; i++) {
                x1 = dx * i;
                x2 = dx * (i+1);
                f1 = my_func(x1, w, alpha, lambda);
                f2 = my_func(x2, w, alpha, lambda);
                
                
                s += (f1 + f2) * dx / 2.0;
                //cout << s << endl;
            }
            
            cout << 2 * PI * s * lambda << endl;
            outputfile << 2 * PI * s * lambda << " ";
        }
        
        cout << endl;
        outputfile << endl;
    }
//    }
     */
    
}
