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
    // QMenu *Menu = new QMenu ("Menu");
    QMenuBar *tool_bar = new QMenuBar (window);
    Window *graph_area = new Window (window);
    QAction *action;

    if (graph_area->parse_command_line (argc, argv) != 1) {
        QMessageBox::warning (0, "Wrong input arguments!",
                                "Wrong input arguments!");
        graph_area->exit_all();
        delete tool_bar;
        delete window;
        return -1;
    }

    action = tool_bar->addAction ("&Change method: 1", graph_area, SLOT (change_func ()));
    action->setShortcut (Qt::Key_1);
    action = tool_bar->addAction ("Double n: 2", graph_area, SLOT (increase_n()));
    action->setShortcut (Qt::Key_2);
    action = tool_bar->addAction ("Reduce n: 3", graph_area, SLOT (decrease_n ()));
    action->setShortcut (Qt::Key_3);
    action = tool_bar->addAction ("Zoom In: 4", graph_area, SLOT (get_closer()));
    action->setShortcut (Qt::Key_4);
    action = tool_bar->addAction ("Zoom Out: 5", graph_area, SLOT (get_further ()));
    action->setShortcut (Qt::Key_5);
    action = tool_bar->addAction ("E&xit: Ctrl+X", graph_area, SLOT (exit_all ()));
    action->setShortcut (QString ("Ctrl+X"));
    action = tool_bar->addAction ("up", graph_area,
                                  SLOT (delta_function_up ()));
    action->setShortcut (Qt::Key_Up);
    action = tool_bar->addAction ("down", graph_area,
                                  SLOT (delta_function_down ()));
    action->setShortcut (Qt::Key_Down);
    action = tool_bar->addAction ("left", graph_area,
                                  SLOT (delta_function_minus()));
    action->setShortcut (Qt::Key_Left);
    action = tool_bar->addAction ("right", graph_area,
                                  SLOT (delta_function_plus ()));
    action->setShortcut (Qt::Key_Right);
    action = tool_bar->addAction("shrink: 9", graph_area, SLOT(shrink_region()));
    action->setShortcut(Qt::Key_9);
    action = tool_bar->addAction("expand: 0", graph_area, SLOT(expand_region()));
    action->setShortcut(Qt::Key_0);

    graph_area->use_all_methods ();
    tool_bar->setMaximumHeight (30);

    window->setMenuBar (tool_bar);
    window->setCentralWidget (graph_area);
    window->setWindowTitle ("Graph");
    window->show ();
    app.exec ();

    delete tool_bar;
    delete window;
    return 0;
}
