#include "window.h"
#include <QObject>
#include <QPainter>
#include <stdio.h>
#include <QMessageBox>
#include <QDebug>

#define DELTA_FUNC_STEP 1

void Window::get_closer () {

    scale_coeff *= 1.2;
    update();
}

void Window::get_further () {

    scale_coeff /= 1.2;
    update();
}

void Window::change_func () {

    scale_coeff = 1.0;

    func_id = (func_id + 1) % 3;
    f_name.clear();

    switch (func_id) {
        case 0:
            f_name.append("Newton; n=");
            break;
        case 1:
            f_name.append("Spline; n=");
            break;
//    case 2:
//      f_name.append("both methods; n=");
//      break;
        case 2:
            f_name.append("Residual; n=");
    }

    f_name.append (QString::number(n));

    update ();
}

void Window::increase_n () {

    n *= 2;
    idx_delta_function *= 2;
  //func_id = 0;
  //memset (x, 0, 40000 * sizeof (double));

    use_all_methods ();

    f_name.clear();

    switch (func_id) {
        case 0:
            f_name.append("Newton; n=");
            break;
        case 1:
            f_name.append("Spline; n=");
            break;
//    case 2:
//      f_name.append("both methods; n=");
//      break;
        case 2:
            f_name.append("Residual; n=");
    }

    f_name.append (QString::number(n));

    update();
}

void Window::decrease_n () {

    if (n < 3) {
        n = 3;
        //qDebug () << n;
    }

    else {
        idx_delta_function /= 2;
        n /= 2;
    }

    use_all_methods ();

    f_name.clear();

    switch (func_id) {
        case 0:
            f_name.append("Newton; n=");
            break;
        case 1:
            f_name.append("Spline; n=");
            break;
        case 2:
            f_name.append("Residual; n=");
    }

    f_name.append (QString::number(n));

    update();
}

void Window::delta_function_plus() {

    if (idx_delta_function < n)
        idx_delta_function++;

    update();
}

void Window::delta_function_minus() {

    if (idx_delta_function > 0)
        idx_delta_function--;

    update();
}

void Window::delta_function_up() {

    val_delta_function += DELTA_FUNC_STEP;

    update();
}

void Window::delta_function_down() {

    val_delta_function -= DELTA_FUNC_STEP;

    update();
}

