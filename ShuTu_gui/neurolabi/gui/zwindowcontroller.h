#ifndef ZWINDOWCONTROLLER_H
#define ZWINDOWCONTROLLER_H

#include <QObject>

class ZWindowController : public QObject
{
  Q_OBJECT
public:
  explicit ZWindowController(QObject *parent = nullptr);

signals:

public slots:
};

#endif // ZWINDOWCONTROLLER_H