#include "expparamsdialog.h"
#include <QDebug>

ExpParamsDialog::ExpParamsDialog(QWidget *parent): QDialog(parent)
//  ,areaValidator(EXPPARAMSDIALOG_AREA_MIN,EXPPARAMSDIALOG_AREA_MAX,this)
{
    ui.setupUi(this);

//    ui.xstartLineEdit->setValidator(&areaValidator);
//    ui.xendLineEdit->setValidator(&areaValidator);
//    ui.ystartLineEdit->setValidator(&areaValidator);
//    ui.yendLineEdit->setValidator(&areaValidator);

    ui.xminSpinBox->setRange(EXPPARAMSDIALOG_AREA_MIN,EXPPARAMSDIALOG_AREA_MAX);
    ui.xmaxSpinBox->setRange(EXPPARAMSDIALOG_AREA_MIN,EXPPARAMSDIALOG_AREA_MAX);

    ui.yminSpinBox->setRange(EXPPARAMSDIALOG_AREA_MIN,EXPPARAMSDIALOG_AREA_MAX);
    ui.ymaxSpinBox->setRange(EXPPARAMSDIALOG_AREA_MIN,EXPPARAMSDIALOG_AREA_MAX);

    ui.xbinSpinBox->setRange(EXPPARAMSDIALOG_BIN_MIN,EXPPARAMSDIALOG_BIN_MAX);
    ui.ybinSpinBox->setRange(EXPPARAMSDIALOG_BIN_MIN,EXPPARAMSDIALOG_BIN_MAX);

    ui.expTimeSpinBox->setRange(EXPPARAMSDIALOG_EXPTIME_MIN,EXPPARAMSDIALOG_EXPTIME_MAX);

    QStringList rate_list;
    rate_list.append("Normal");
    QString filename;

    init(filename,rate_list,0,EXPPARAMSDIALOG_BIN_MIN,EXPPARAMSDIALOG_BIN_MIN);

    setArea(EXPPARAMSDIALOG_AREA_MIN,EXPPARAMSDIALOG_AREA_MIN,EXPPARAMSDIALOG_AREA_MAX,EXPPARAMSDIALOG_AREA_MAX);

    connect(ui.xminSpinBox, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged),
          [=](int i){ if (ui.xmaxSpinBox->value() < i) ui.xmaxSpinBox->setValue(i);});
    connect(ui.xmaxSpinBox, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged),
          [=](int i){ if (ui.xminSpinBox->value() > i) ui.xminSpinBox->setValue(i);});
    connect(ui.yminSpinBox, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged),
          [=](int i){ if (ui.ymaxSpinBox->value() < i) ui.ymaxSpinBox->setValue(i);});
    connect(ui.ymaxSpinBox, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged),
          [=](int i){ if (ui.yminSpinBox->value() > i) ui.yminSpinBox->setValue(i);});
}


void ExpParamsDialog::init(QString &filename, QStringList &items, int index, int xbin, int ybin)
{
    ui.filenameLineEdit->setText(filename);
    qDebug() << "EXPPARAM INIT: " << filename;
    ui.rateComboBox->clear();
    ui.rateComboBox->addItems(items);
    ui.rateComboBox->setCurrentIndex(index);

    ui.xbinSpinBox->setValue(xbin);
    ui.ybinSpinBox->setValue(ybin);
}


void ExpParamsDialog::setArea(int xmin, int ymin, int xmax, int ymax)
{
    ui.xminSpinBox->setValue(xmin);

    ui.yminSpinBox->setValue(ymin);

    ui.xmaxSpinBox->setValue(xmax);

    ui.ymaxSpinBox->setValue(ymax);
}


void ExpParamsDialog::setExptimeRange(double min, double max)
{
    ui.expTimeSpinBox->setRange(min,max);
}


void ExpParamsDialog::setBinRange(int min, int max)
{
    ui.xbinSpinBox->setRange(min,max);
    ui.ybinSpinBox->setRange(min,max);
}


void ExpParamsDialog::setAreaRange(int min, int max)
{
    ui.xminSpinBox->setRange(min,max);
    ui.xmaxSpinBox->setRange(min,max);
    ui.yminSpinBox->setRange(min,max);
    ui.ymaxSpinBox->setRange(min,max);
}


void ExpParamsDialog::getArea(int *xmin, int *ymin, int *xmax, int *ymax) const
{
    if ( ui.frameGeometryGroupBox->isChecked() ) {
        *xmin = ui.xminSpinBox->value();
        *xmax = ui.xmaxSpinBox->value();
        *ymin = ui.yminSpinBox->value();
        *ymax = ui.ymaxSpinBox->value();
    } else { // full frame
        *xmin = ui.xminSpinBox->minimum();
        *xmax = ui.xmaxSpinBox->maximum();
        *ymin = ui.yminSpinBox->minimum();
        *ymax = ui.ymaxSpinBox->maximum();
    }
}


void ExpParamsDialog::getBin(int *xbin, int *ybin) const
{
    *xbin = ui.xbinSpinBox->value();
    *ybin = ui.ybinSpinBox->value();
}


double ExpParamsDialog::getExptime() const
{
    return ui.expTimeSpinBox->value();
}


int ExpParamsDialog::getExpNum() const
{
    return ui.expNumSpinBox->value();
}


QString ExpParamsDialog::getFilename() const
{
    return ui.filenameLineEdit->text();
}


void ExpParamsDialog::setFilename(QString filename)
{
    ui.filenameLineEdit->setText(filename);
}
