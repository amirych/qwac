#include "expparamsdialog.h"

ExpParamsDialog::ExpParamsDialog(QWidget *parent): QDialog(parent),
    areaValidator(EXPPARAMSDIALOG_AREA_MIN,EXPPARAMSDIALOG_AREA_MAX,this)
{
    ui.setupUi(this);

    ui.xstartLineEdit->setValidator(&areaValidator);
    ui.xendLineEdit->setValidator(&areaValidator);
    ui.ystartLineEdit->setValidator(&areaValidator);
    ui.yendLineEdit->setValidator(&areaValidator);

    ui.xbinSpinBox->setRange(EXPPARAMSDIALOG_BIN_MIN,EXPPARAMSDIALOG_BIN_MAX);
    ui.ybinSpinBox->setRange(EXPPARAMSDIALOG_BIN_MIN,EXPPARAMSDIALOG_BIN_MAX);

    ui.expTimeSpinBox->setRange(EXPPARAMSDIALOG_EXPTIME_MIN,EXPPARAMSDIALOG_EXPTIME_MAX);

    QStringList rate_list;
    rate_list.append("Normal");
    QString filename;

    init(filename,rate_list,0,EXPPARAMSDIALOG_BIN_MIN,EXPPARAMSDIALOG_BIN_MIN);

    setArea(EXPPARAMSDIALOG_AREA_MIN,EXPPARAMSDIALOG_AREA_MIN,EXPPARAMSDIALOG_AREA_MAX,EXPPARAMSDIALOG_AREA_MAX);
}


void ExpParamsDialog::init(QString &filename, QStringList &items, int index, int xbin, int ybin)
{
    ui.filenameLineEdit->setText(filename);

    ui.rateComboBox->clear();
    ui.rateComboBox->addItems(items);
    ui.rateComboBox->setCurrentIndex(index);

    ui.xbinSpinBox->setValue(xbin);
    ui.ybinSpinBox->setValue(ybin);
}


void ExpParamsDialog::setArea(int xmin, int ymin, int xmax, int ymax)
{
    QString val;

    val = QString::number(xmin);
    ui.xstartLineEdit->setText(val);

    val = QString::number(ymin);
    ui.ystartLineEdit->setText(val);

    val = QString::number(xmax);
    ui.xendLineEdit->setText(val);

    val = QString::number(ymax);
    ui.yendLineEdit->setText(val);
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
    areaValidator.setRange(min,max);
}


void ExpParamsDialog::getArea(int *xmin, int *ymin, int *xmax, int *ymax) const
{
    *xmin = ui.xstartLineEdit->text().toInt();
    *xmax = ui.xendLineEdit->text().toInt();
    *ymin = ui.ystartLineEdit->text().toInt();
    *ymax = ui.yendLineEdit->text().toInt();
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


QString ExpParamsDialog::getFilename() const
{
    return ui.filenameLineEdit->text();
}
