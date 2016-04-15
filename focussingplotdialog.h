#ifndef FOCUSSINGPLOTDIALOG_H
#define FOCUSSINGPLOTDIALOG_H

#include <QDialog>


#include "ui_focussingCurvePlot.h"

class FocussingPlotDialog : public QDialog
{
public:
    FocussingPlotDialog(QWidget *parent = nullptr);

private:
    Ui::FocussingPlotForm ui;
};

#endif // FOCUSSINGPLOTDIALOG_H
