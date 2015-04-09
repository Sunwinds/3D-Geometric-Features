#pragma warning(disable:4267)

#include "ddtt_widget.h"
#include "ui_ddtt_widget.h"

ddtt_widget::ddtt_widget(QWidget *parent) : QWidget(parent), ui(new Ui::ddtt_widget)
{
    ui->setupUi(this);
}

ddtt_widget::~ddtt_widget()
{
    delete ui;
}
