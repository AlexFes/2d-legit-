#ifndef WINDOW_H
#define WINDOW_H

#include <math.h>
#include <QWidget>
#include <QPushButton>
#include <QPainter>
#include <QLabel>
#include <QRadioButton>

#define MIN_N 50
#define MIN_FOR_DRAW 1e-6
#define MIN_FOR_COMPARE 1e-12

class Window : public QWidget {
    Q_OBJECT

    int func_id;
    QString f_name;
    double a;
    double b;
    double res1;
    double res2;
    int n;
    double (*f) (double);
    double (*df) (double);
    int idx_delta_function;
    double val_delta_function;

    QWidget *parent_save;

    double *x;
    double *f_array;

    double *v;
    double *ksi;
    double *func;

    double *c_1,
           *c_2,
           *c_3;

   double *ak_1,
          *ak_2,
          *ak_3,
          *ak_4;

public:
    Window (QWidget *parent);
    ~Window ();

    QSize minimumSizeHint () const;
    QSize sizeHint () const;

    double scale_coeff;

    void use_all_methods ();
    void AkimaNumbers(int n);
    double AkimaSolve(double x0, int n);
    void SplineNumbers (int n);
    double SplineSolve (double x_0, int n);

    void init_drawing_plane (QPainter &painter);

public:
    int parse_command_line (int argc, char *argv[]);

public slots:
    void change_func ();
    void get_closer();
    void get_further();
    void increase_n();
    void decrease_n();
    void delta_function_plus();
    void delta_function_minus();
    void delta_function_up();
    void delta_function_down();
    void exit_all();
    void expand_region();
    void shrink_region();

protected:
    void paintEvent (QPaintEvent *event);

private:
    void update_arrays (int n);
};

#endif
