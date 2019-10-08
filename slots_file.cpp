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
//    scale_coeff = 1.0;
    func_id = (func_id + 1) % 3;
    update ();
}

void Window::increase_n () {
    n *= 2;
    idx_delta_function *= 2;
  //func_id = 0;
  //memset (x, 0, 40000 * sizeof (double));
    use_all_methods();
    update();
}

void Window::decrease_n() {
    idx_delta_function /= 2;
    n /= 2;

    if (n < 4) {
        n = 4;
        //qDebug () << n;
    }


    use_all_methods();
    update();
}

void Window::expand_region() {
    a -= 10;
    a = (a < -300 ? -300 : a);
    b += 10;
    b = (b > 300 ? 300 : b);

    use_all_methods();
    update();
}

void Window::shrink_region() {
    a += 10;
    a = (a > -10 ? -10 : a);
    b -= 10;
    b = (b < 10 ? 10 : b);

    use_all_methods();
    update();
}

void Window::delta_function_plus() {
    if (idx_delta_function < n - 1)
        idx_delta_function++;

    use_all_methods();
    update();
}

void Window::delta_function_minus() {
    if (idx_delta_function > 0)
        idx_delta_function--;

    use_all_methods();
    update();
}

void Window::delta_function_up() {
    val_delta_function += DELTA_FUNC_STEP;

    use_all_methods();
    update();
}

void Window::delta_function_down() {
    val_delta_function -= DELTA_FUNC_STEP;

    use_all_methods();


    update();
}
