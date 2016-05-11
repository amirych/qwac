#ifndef EXPPARAMSDIALOG_H
#define EXPPARAMSDIALOG_H

#include <QDialog>
#include <QString>
#include <QList>
#include <QStringList>
#include <QIntValidator>

#include "ui_expParamsForm.h"

#define EXPPARAMSDIALOG_EXPTIME_MIN 0.0
#define EXPPARAMSDIALOG_EXPTIME_MAX 1.0E6
#define EXPPARAMSDIALOG_BIN_MIN 1
#define EXPPARAMSDIALOG_BIN_MAX 100
#define EXPPARAMSDIALOG_AREA_MIN 1      // by defaults follow the FITS coordinate system standard (origin is in [1,1])
#define EXPPARAMSDIALOG_AREA_MAX 100000



class ExpParamsDialog : public QDialog
{
    Q_OBJECT
public:
    ExpParamsDialog(QWidget *parent = nullptr);

    void init(QString &filename, QStringList &items, int index, int xbin, int ybin);

    void setExptimeRange(double min, double max);
    void setBinRange(int min, int max);
    void setAreaRange(int min, int max);

    void setArea(int xmin, int ymin, int xmax, int ymax);

    void getArea(int *xmin, int *ymin, int *xmax, int *ymax) const;
    void getBin(int *xbin, int *ybin) const;
    double getExptime() const;
    int getExpNum() const;
    QString getFilename() const;

public slots:
    void setFilename(QString filename);

private:
    Ui::ExpParamsForm ui;

//    QIntValidator areaValidator;
};

#endif // EXPPARAMSDIALOG_H
