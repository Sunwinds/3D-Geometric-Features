#ifndef DDTT_WIDGET_H
#define DDTT_WIDGET_H

#include <QWidget>

namespace Ui {
class ddtt_widget;
}

class ddtt_widget : public QWidget
{
    Q_OBJECT

public:
    explicit ddtt_widget(QWidget *parent = 0);
    ~ddtt_widget();

public:
    Ui::ddtt_widget *ui;
};

#endif // DDTT_WIDGET_H
