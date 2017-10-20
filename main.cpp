#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include <QAction>
#include <QDebug>
#include <QMenuBar>
#include <QMessageBox>
#include <stdio.h>
#include <stdlib.h>
#include "window.h"

int main (int argc, char *argv[]) {

    QApplication app (argc, argv);

    QMainWindow *window = new QMainWindow;
    QMenu *Menu = new QMenu ("Menu");
    QMenuBar *tool_bar = new QMenuBar (window);
    Window *graph_area = new Window (window);
    QAction *action;

    if (graph_area->parse_command_line (argc, argv) != 1) {
        QMessageBox::warning (0, "Wrong input arguments!",
                                "Wrong input arguments!");

        return -1;
    }
      
    action = Menu->addAction ("&Change function", graph_area, SLOT (change_func ()));
    action->setShortcut (Qt::Key_1);

    action = Menu->addAction ("Double n", graph_area, SLOT (increase_n()));
    action->setShortcut (Qt::Key_2);

    action = Menu->addAction ("Reduce n", graph_area, SLOT (decrease_n ()));
    action->setShortcut (Qt::Key_3);

    action = Menu->addAction ("Zoom In", graph_area, SLOT (get_closer()));
    action->setShortcut (Qt::Key_4);

    action = Menu->addAction ("Zoom Out", graph_area, SLOT (get_further ()));
    action->setShortcut (Qt::Key_5);

    action = Menu->addAction ("E&xit", graph_area, SLOT (exit_all ()));
    action->setShortcut (QString ("Ctrl+X"));

    action = Menu->addAction ("Delta func &u&p", graph_area,
                                  SLOT (delta_function_up ()));
    action->setShortcut (Qt::Key_Up);

    action = Menu->addAction ("Delta func &d&o&w&n", graph_area,
                                  SLOT (delta_function_down ()));
    action->setShortcut (Qt::Key_Down);

    action = Menu->addAction ("Delta func &l&e&f&t", graph_area,
                                  SLOT (delta_function_minus()));
    action->setShortcut (Qt::Key_Left);

    action = Menu->addAction ("Delta func &r&i&g&h&t", graph_area,
                                  SLOT (delta_function_plus ()));
    action->setShortcut (Qt::Key_Right);


    graph_area->use_all_methods ();

    tool_bar->setMaximumHeight (30);
    tool_bar->addMenu(Menu);

    window->setMenuBar (tool_bar);
    window->setCentralWidget (graph_area);
    window->setWindowTitle ("Graph");

    window->show ();

    app.exec ();

    delete tool_bar;
    delete window;

    return 0;
}
