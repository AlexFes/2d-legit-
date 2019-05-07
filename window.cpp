#include <QObject>
#include <QPainter>
#include <stdio.h>
#include <math.h>

#include "window.h"
#include <QLineEdit>
#include <QLabel>
#include <QHBoxLayout>
#include <QSpacerItem>
#include <QPushButton>
#include <QDebug>
#include <QMessageBox>
#include <QRadioButton>

#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10

#define DELTA 0.1

static inline double f_0 (double x) {
  // return x;
  //return exp(x);
  return 3*x*x+2*x+1;
  //return x*x*x+3*x*x+2*x+1;
}

static inline double df_0 (double x) {
  // return 1 + x - x;
  //return exp(x);
  return 6*x+2;
  //return 3*x*x+6*x+2;
}
/*
static inline double d_1 (double x)
{
  (void) x;
  return 1 + cos (x);
}

static inline double dd_1 (double x)
{
  (void) x;
  return -sin (x);
}

static inline double d_2 (double x)
{
  (void) x;
  return -sin (x);
}

static inline double dd_2 (double x)
{
  (void) x;
  return -cos (x);
}
*/
Window::Window (QWidget *parent)
  : QWidget (parent) {
    a = DEFAULT_A;
    b = DEFAULT_B;
    n = DEFAULT_N;
    parent_save = parent;

    f = f_0;
    df = df_0;
//  d1 = d_1;
//  d2 = d_2;

    x          = 0;
    f_array    = 0;
//  method_1_y = 0;
//  method_2_y = 0;
//  workspace_array = 0;

//  m_coef_a = 0;
//  m_coef_b = 0;
//  m_coef_c = 0;
//  m_coef_d = 0;

    ksi  = 0;
    v    = 0;
    func = 0;

//  coeff_v_left   = 0;
//  coeff_v_center = 0;
//  coeff_v_right  = 0;
//  rhs            = 0;

    c_1 = 0;
    c_2 = 0;
    c_3 = 0;

    ak_1 = 0;
    ak_2 = 0;
    ak_3 = 0;
    ak_4 = 0;

    // res = 0;

    change_func ();
}

Window::~Window() {
    delete [] x;
    delete [] f_array;
//  delete [] method_1_y;
//  delete [] method_2_y;
//  delete [] workspace_array;

//  delete [] m_coef_a;
//  delete [] m_coef_b;
//  delete [] m_coef_c;
//  delete [] m_coef_d;

    delete [] ksi;
    delete [] v;
    delete [] func;

//  delete [] coeff_v_left;
//  delete [] coeff_v_center;
//  delete [] coeff_v_right;
//  delete [] rhs;

    delete [] c_1;
    delete [] c_2;
    delete [] c_3;

    delete[]ak_1;
    delete[]ak_2;
    delete[]ak_3;
    delete[]ak_4;

    // delete[]res;
}

void Window::exit_all () {
    parent_save->close();
    delete this;
}

QSize Window::minimumSizeHint () const {
    return QSize (100, 100);
}

QSize Window::sizeHint () const {
    return QSize (1000, 1000);
}


int Window::parse_command_line (int argc, char *argv[]) {
    if (argc == 1)
        return -1;
    if (argc == 2)
        return -2;

    if (sscanf (argv[1], "%lf", &a) != 1
        || sscanf (argv[2], "%lf", &b) != 1
        || b - a < MIN_FOR_COMPARE
        || (argc > 3 && sscanf (argv[3], "%d", &n) != 1)
        || n <= 0)
        return -3;

    scale_coeff = 1.0;
    func_id = 0;
    idx_delta_function = 0;
    val_delta_function = 0;

    return 1;
}

void Window::init_drawing_plane (QPainter &painter, int N) {
    double w = width(), h = height();
    f_name.clear();

    switch (func_id) {
        case 0:
            f_name.append("Akima; n=");
            break;
        case 1:
            f_name.append("Spline; n=");
            break;
        case 2:
            f_name.append("Residual; n=");
    }

    f_name.append (QString::number(N));

    painter.setPen ("black");
    painter.drawText (0, 20, f_name);
    painter.translate (0.5 * w, 0.5 * h);

    painter.scale (scale_coeff, -scale_coeff);

    painter.setPen ("red");
    painter.drawLine (-w/2, 0., w/2, 0.);
    painter.drawLine (0., -h/2, 0., h/2);
}

void Window::paintEvent (QPaintEvent*) {
    QPainter painter (this);
    painter.save ();

    int i, N = n;
    double x1, x2, y1, y2;

    if (n < 4)
        N = 4;
    else if (n > 3000)
        N = 3000;

    double step = (b-a)/width();
    init_drawing_plane (painter, N);

    /// draw real graph
    if (func_id < 3) {
        painter.setPen ("blue");
        x1 = a;
        y1 = f (a);
        func[idx_delta_function] += val_delta_function;

        for (i = 1, x2 = x1 + step; i < width(); i++, x2 += step) {
            y2 = f (x2);

            if (fabs (x[idx_delta_function] - x2) <= step)
                y2 += val_delta_function;

            painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
            x1 = x2, y1 = y2;
        }

        x2 = b;
        y2 = f(x2);
    }

    /// draw Akima
    if (func_id == 0) {
        painter.setPen ("green");
        x1 = a;
        y1 = f(a);

        for (i = 0; i < width(); i++) {
            x2 = x1 + step;
            y2 = AkimaSolve(x2, N);

            painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
            x1 = x2, y1 = y2;
        }

        x2 = b;
        y2 = f(x2);

        painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
    }

    /// draw Spline
    if (func_id == 1) {
        painter.setPen ("black");
        x1 = a;
        y1 = f (a);

        for (i = 0; i < width(); i++) {
            x2 = x1 + step;
            y2 = SplineSolve(x2, N);

            painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
            x1 = x2, y1 = y2;
        }

        x2 = b;
        y2 = f(x2);

        painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
    }

    /// draw Residual
    if (func_id == 2) {
        painter.setPen ("green");

        x1 = a;
        y1 = 0.;//f (x1);

        for (i = 1; i < width (); i++) {
            x2 = x1 + step;
            y2 = f(x2) - AkimaSolve(x2, N);

            painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
            x1 = x2, y1 = y2;
        }

        x2 = b;
        y2 = 0.;
        painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
        painter.setPen ("black");
        x1 = a;
        y1 = 0.;//f (x1);

        for (i = 1; i < width (); i++) {
            x2 = x1 + step;
            y2 = f(x2) - SplineSolve(x2, N);

            painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
            x1 = x2, y1 = y2;
        }

        x2 = b;
        y2 = 0.;
        painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
    }

    painter.restore ();
}


static inline void free_array (double *arr) {
    if (arr)
        delete [] arr;

    arr = 0;
}

void Window::update_arrays (int n) {
    free_array (x);
    free_array (f_array);

    free_array (ksi);
    free_array (v);
    free_array (func);

    free_array (c_1);
    free_array (c_2);
    free_array (c_3);

    free_array (ak_1);
    free_array (ak_2);
    free_array (ak_3);
    free_array (ak_4);
    // free_array (res);

    x          = new double [n];
    f_array    = new double [n];

    ksi  = new double [n];
    v    = new double [n];
    func = new double [n];

    c_1 = new double [n];
    c_2 = new double [n];
    c_3 = new double [n];

    ak_1 = new double [n];
    ak_2 = new double [n];
    ak_3 = new double [n];
    ak_4 = new double [n];
    // res = new double [2*n];
}
