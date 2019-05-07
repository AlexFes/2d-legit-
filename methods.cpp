#include "window.h"
#include <stdio.h>
#include <QDebug>
#include <QMessageBox>
#include <math.h>

void Window::use_all_methods () {
    int N = n;

    if (n < 4)
        N = 4;
    else if (n > 3000)
        N = 3000;

    update_arrays (N+1);
    double delta_x = (b-a)/(N-1);

    for (int i = 0; i < N; i++) {
        x[i] = a + i*delta_x;
        func[i] = f(x[i]);
        f_array[i] = func[i];
    }

    if (idx_delta_function < N) {
        func[idx_delta_function] += val_delta_function;
        f_array[idx_delta_function] += val_delta_function;
    }

    AkimaNumbers(N);

    for (int i = 0; i < N; i++) {
        x[i] = a + i * delta_x;
        func[i] = f(x[i]);
        f_array[i] = func[i];
    }

    func[idx_delta_function] += val_delta_function;
    f_array[idx_delta_function] += val_delta_function;

    SplineNumbers(N);
}

void Window::AkimaNumbers(int n) {
    double* d = new double[n];

    d[0] = df(x[0]);
    d[1] = df(x[1]);
    d[n-2] = df(x[n-2]);
    d[n-1] = df(x[n-1]);

    for (int i = 2; i < n - 2; ++i) {
        double w_1 = fabs((func[i] - func[i-1])/(x[i] - x[i-1]) - (func[i-1] - func[i-2])/(x[i-1] - x[i-2]));
        double w_2 = fabs((func[i+2] - func[i+1])/(x[i+2] - x[i+1]) - (func[i+1] - func[i])/(x[i+1] - x[i]));

        if (w_1*w_1 + w_2*w_2 < 1e-12) {
            d[i] = ((x[i+1] - x[i])*(func[i] - func[i-1])/(x[i] - x[i-1]) + (x[i] - x[i-1])*(func[i+1] - func[i])/(x[i+1] - x[i]))/(x[i+1] - x[i-1]);
        } else {
            d[i] = (w_2*(func[i] - func[i-1])/(x[i] - x[i-1]) + w_1*(func[i+1] - func[i])/(x[i+1] - x[i]))/(w_2 - w_1);
        }
    }

    for (int i = 0; i < n - 1; ++i) {
        ak_1[i] = func[i];
        ak_2[i] = d[i];
        ak_3[i] = (3*(func[i+1] - func[i])/(x[i+1] - x[i]) - 2*d[i] - d[i+1])/(x[i+1] - x[i]);
        ak_4[i] = (d[i] + d[i+1] - 2*(func[i+1] - func[i])/(x[i+1] - x[i]))/((x[i+1] - x[i])*(x[i+1] - x[i]));
    }

    delete[]d;
}

double Window::AkimaSolve(double x0, int n) {
    double dx = (b - a)/(n - 1);
    int i = (int) trunc((x0 - a)/dx);
    if (i == n-1) {
        i = n - 2;
    }

    // printf("\n i=%d  appr=%f\n", i, ak_1[i] + ak_2[i]*(x0 - x[i]) + ak_3[i]*(x0 - x[i])*(x0 - x[i]) + ak_4[i]*(x0 - x[i])*(x0 - x[i])*(x0 - x[i]));

    return ak_1[i] + ak_2[i]*(x0 - x[i]) + ak_3[i]*(x0 - x[i])*(x0 - x[i]) + ak_4[i]*(x0 - x[i])*(x0 - x[i])*(x0 - x[i]);
}

void Window::SplineNumbers (int n) {
    double tmp1;
    double *arr = new double[(n+1)*3];

    for (int i = 1; i < n; ++i)
        ksi[i] = 0.5 * (x[i - 1] + x[i]);

    ksi[0] = a - (ksi[2] - ksi[1])/2;
    ksi[n] = b + (ksi[n - 1] - ksi[n - 2])/2;

    for (int i = 1; i < n; ++i) {
        arr[i*3] = 1/(x[i-1]-ksi[i-1]) - 1/(ksi[i]-ksi[i-1]);
        arr[i*3 + 2] = 1/(ksi[i+1]-x[i]) - 1/(ksi[i+1]-ksi[i]);
        arr[i*3 + 1] = 1/(ksi[i]-x[i-1]) + 1/(ksi[i]-ksi[i-1]) +
                       1/(x[i]-ksi[i]) + 1/(ksi[i+1]-ksi[i]);
        v[i] = (1/(x[i-1]-ksi[i-1]) + 1/(ksi[i]-x[i-1]))*func[i-1] +
               (1/(x[i]-ksi[i]) + 1/(ksi[i+1]-x[i]))*func[i];
    }

    arr[0] = (1/(ksi[1]-ksi[0]))*(1/(x[0]-ksi[0]));         //two extra equations
    arr[2] = -1*(1/(ksi[2]-ksi[1]))*(1/(ksi[2]-x[1]));      //c(3,0)==c(3,1)
    arr[1] = (1/(ksi[1]-ksi[0]))*(1/(ksi[1]-x[0])) -
             (1/(ksi[2]-ksi[1]))*(1/(x[1]-ksi[1]));
    v[0] = (1/(ksi[1]-ksi[0]))*(1/(x[0]-ksi[0])+1/(ksi[1]-x[0]))*func[0] -
           (1/(ksi[2]-ksi[1]))*(1/(x[1]-ksi[1])+1/(ksi[2]-x[1]))*func[1];

    arr[n*3] = (1/(ksi[n-1]-ksi[n-2]))*(1/(x[n-2]-ksi[n-2])); //c(3,n-2)==c(3,n-1)
    arr[n*3 + 2] = -1*(1/(ksi[n]-ksi[n-1]))*(1/(ksi[n]-x[n-1]));
    arr[n*3 + 1] = (1/(ksi[n-1]-ksi[n-2]))*(1/(ksi[n-1]-x[n-2])) -
                   (1/(ksi[n]-ksi[n-1]))*(1/(x[n-1]-ksi[n-1]));
    v[n] = (1/(ksi[n-1]-ksi[n-2]))*(1/(x[n-2]-ksi[n-2])+1/(ksi[n-1]-x[n-2]))*func[n-2] -
           (1/(ksi[n]-ksi[n-1]))*(1/(x[n-1]-ksi[n-1])+1/(ksi[n]-x[n-1]))*func[n-1];

    arr[2]/=arr[0];
    arr[1]/=arr[0];
    v[0]/=arr[0];
//    if (fabs(arr[0])<1e-10) qDebug() << "Triggered i = " << 0;

    arr[4]-=arr[3]*arr[1];
    arr[5]-=arr[3]*arr[2];
    v[1]-=arr[3]*v[0];

    for (int i = 1; i < n-1; ++i) {
        arr[i*3 + 2]/=arr[i*3 + 1];
        v[i]/=arr[i*3 + 1];
//        if (fabs(arr[i*3 + 1])<1e-10) qDebug() << "Triggered i = " << i*3 + 1;
        arr[(i+1)*3 + 1]-=arr[(i+1)*3]*arr[i*3 + 2];
        v[i+1]-=arr[(i+1)*3]*v[i];
    }

    arr[n*3 + 1]-=arr[n*3]*arr[(n-2)*3 + 2];
    v[n]-=arr[n*3]*v[n-2];

    arr[(n-1)*3 + 2]/=arr[(n-1)*3 + 1];
    v[n-1]/=arr[(n-1)*3 + 1];
//    if (fabs(arr[(n-1)*3 + 1])<1e-10) qDebug() << "Triggered i = " << (n-1)*3 + 1;

    arr[n*3 + 2]-=arr[n*3 + 1]*arr[(n-1)*3 + 2];
    v[n]-=arr[n*3 + 1]*v[n-1];
    v[n]/=arr[n*3 + 2];

//    if (fabs(arr[n*3 + 2])<1e-10) qDebug() << "Triggered i = " << n*3 + 2;

    for (int i = n; i > 1; --i)
        v[i-1]-=arr[(i-1)*3 + 2]*v[i];

    v[0]-=arr[2]*v[2];
    v[0]-=arr[1]*v[1];

    for (int i = 0; i < n; i ++) {
        c_1[i] = v[i];
        tmp1 = ((v[i+1] - func[i])/(ksi[i+1] - x[i]) - (func[i] - v[i])/
               (x[i] - ksi[i]))/(ksi[i+1] - ksi[i]);
        c_2[i] = (func[i] - v[i])/(x[i] - ksi[i]) - (x[i] - ksi[i])*tmp1;
        c_3[i] = tmp1;
    }

    delete[]arr;
}

double Window::SplineSolve (double x_0, int n) {
    double dx = (b - a)/(n - 1);
    int i = (int) trunc((x_0 - ksi[0])/dx);

    return c_1[i] + c_2[i]*(x_0 - ksi[i]) + c_3[i]*(x_0 - ksi[i])*(x_0 - ksi[i]);
}
