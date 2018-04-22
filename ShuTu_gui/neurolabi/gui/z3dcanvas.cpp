#include "z3dcanvas.h"

#include <QPainter>

#include "z3dnetworkevaluator.h"
#include <algorithm>
#include "z3dcanvaseventlistener.h"
#include "QsLog/QsLog.h"
#ifdef _QT5_
#include <QWindow>
#endif
#include "zpainter.h"
#include "zstackdrawable.h"

#ifdef _QT5_
Z3DCanvas::Z3DCanvas(const QString &title, int width, int height, QWidget* parent, Qt::WindowFlags f)
#else
Z3DCanvas::Z3DCanvas(const QString &title, int width, int height, const QGLFormat &format,
                     QWidget* parent, const QGLWidget *shareWidget, Qt::WindowFlags f)
#endif
  : QGraphicsView(parent)
  , m_fullscreen(false)
  , m_glWidget(NULL)
  , m_3dScene(NULL)
  , m_networkEvaluator(NULL)
  , m_fakeStereoOnce(false)
{
  setAlignment(Qt::AlignLeft | Qt::AlignTop);
  resize(width, height);

#ifdef _QT5_
  m_glWidget = new QOpenGLWidget(nullptr, f);
#else
  m_glWidget = new QGLWidget(format, NULL, shareWidget, f);
  m_glWidget->makeCurrent();
#endif
  m_isStereoScene = m_glWidget->format().stereo();
  m_3dScene = new QGraphicsScene(0, 0, width, height, this);

  setViewport(m_glWidget);
  setViewportUpdateMode(FullViewportUpdate);
  setScene(m_3dScene);

  setWindowTitle(title);
  setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform);
  setOptimizationFlags(QGraphicsView::DontSavePainterState |
                       QGraphicsView::DontAdjustForAntialiasing);

  setAcceptDrops(true);
  setFocusPolicy(Qt::StrongFocus);

#ifndef __APPLE__
  setStyleSheet("border-style: none;");
#endif
  setMouseTracking(true);

//#if defined(_FLYEM_)
  connect(&m_interaction, SIGNAL(decorationUpdated()),
          this->viewport(), SLOT(update()));
  connect(&m_interaction, SIGNAL(strokePainted(ZStroke2d*)),
          this, SIGNAL(strokePainted(ZStroke2d*)));
//#endif
}

Z3DCanvas::~Z3DCanvas() {}

void Z3DCanvas::toggleFullScreen()
{
  if (m_fullscreen) {
    m_fullscreen = false;
    showNormal();
  } else {
    showFullScreen();
    m_fullscreen = !m_fullscreen;
  }
}

void Z3DCanvas::enterEvent(QEvent* e)
{
  broadcastEvent(e, width(), height());
}

void Z3DCanvas::leaveEvent(QEvent* e)
{
  broadcastEvent(e, width(), height());
}

void Z3DCanvas::mousePressEvent(QMouseEvent* e)
{
  broadcastEvent(e, width(), height());

  m_interaction.processMousePressEvent(e);
}

bool Z3DCanvas::suppressingContextMenu() const
{
//#if defined(_FLYEM_)
  if (m_interaction.isStateOn(ZInteractionEngine::STATE_DRAW_STROKE) ||
      m_interaction.isStateOn(ZInteractionEngine::STATE_DRAW_RECT)) {
    return true;
  }
//#endif

  return false;
}

void Z3DCanvas::mouseReleaseEvent (QMouseEvent* e)
{
  broadcastEvent(e, width(), height());
  m_interaction.processMouseReleaseEvent(e);
  setCursor(m_interaction.getCursorShape());
}

void Z3DCanvas::mouseMoveEvent(QMouseEvent*  e)
{
  m_interaction.processMouseMoveEvent(e);

  if (!m_interaction.lockingMouseMoveEvent()) {
    broadcastEvent(e, width(), height());
  }
}

void Z3DCanvas::mouseDoubleClickEvent(QMouseEvent* e)
{
  broadcastEvent(e, width(), height());
}

void Z3DCanvas::wheelEvent(QWheelEvent* e)
{
  broadcastEvent(e, width(), height());
}

void Z3DCanvas::keyPressEvent(QKeyEvent* event)
{
  if (!m_interaction.processKeyPressEvent(event)) {
    broadcastEvent(event, width(), height());
  }

  setCursor(m_interaction.getCursorShape());
}

void Z3DCanvas::keyReleaseEvent(QKeyEvent* event)
{
  broadcastEvent(event, width(), height());
}

void Z3DCanvas::resizeEvent(QResizeEvent *event)
{
#ifdef _QT5_
  glGetError(); // opengl error from qt5?
#endif
  getGLFocus();
  QGraphicsView::resizeEvent(event);
  if (m_3dScene)
    m_3dScene->setSceneRect(QRect(QPoint(0, 0), event->size()));

  emit canvasSizeChanged(event->size().width() * getDevicePixelRatio(),
                         event->size().height() * getDevicePixelRatio());
}

void Z3DCanvas::paintEvent(QPaintEvent *event)
{
  getGLFocus();
  QGraphicsView::paintEvent(event);
}

void Z3DCanvas::dragEnterEvent(QDragEnterEvent *event)
{
  event->ignore();
}

void Z3DCanvas::dropEvent(QDropEvent *event)
{
  event->ignore();
}

void Z3DCanvas::drawBackground(QPainter *painter, const QRectF &)
{
#ifdef _QT5_
  glGetError(); // opengl error from qt5?
#endif

  if (!m_networkEvaluator) {
    return;
  }

  // QPainter set glclearcolor to white, we set it back
  glClearColor(0.f, 0.f, 0.f, 0.f);

  m_networkEvaluator->process(m_isStereoScene || m_fakeStereoOnce);
  m_fakeStereoOnce = false;


#if 1
  QList<ZStackObject*> drawableList = m_interaction.getDecorationList();

  foreach (ZStackObject *drawable, drawableList) {
    //drawable->setVisible(true);
    drawable->display(painter, 0, ZStackObject::NORMAL,
                      ZStackObject::DISPLAY_SLICE_SINGLE, NeuTube::Z_AXIS);
  }
#else
  UNUSED_PARAMETER(painter);
#endif

#ifdef _DEBUG_2
  painter->setPen(QColor(255, 0, 0));
  painter->drawRect(QRect(10, 10, 40, 60));
#endif

#if 1
  QStringList text;
  if (m_interaction.getKeyMode() == ZInteractionEngine::KM_SWC_SELECTION) {
      text.append("Selection mode on:");
      text.append("  1: downstream;");
      text.append("  2: upstream;");
      text.append("  3: neighboring nodes;");
      text.append("  4: host branch;");
      text.append("  5: connected nodes;");
      text.append("  6: inverse selection;");
      text.append("  7: select small trees;");
  }

  drawText(*painter, text);
#endif

  //ZPainter painter()
  //painter->drawRect(QRect(10, 10, 40, 60));
  CHECK_GL_ERROR;
}

void Z3DCanvas::drawText(QPainter &painter, const QStringList &text)
{
  if (!text.empty()) {
    int maxLength = 0;
    QString compText;
    foreach (const QString &str, text) {
      if (str.length() > maxLength) {
        maxLength = str.length();
      }
      compText += str + "\n";
    }

    maxLength = std::min(maxLength, 80);

    int width = maxLength * 8;
    int height = text.length() * 18;

    if (width > 0 && height > 0) {
      QPixmap pixmap(width, height);
      pixmap.fill(QColor(0, 0, 0, 128));

      QPainter bufferPainter(&pixmap);
      bufferPainter.setPen(QColor(255, 255, 255));

      //    bufferPainter.fillRect(pixmap.rect(), QColor(0, 0, 0, 0));
      bufferPainter.drawText(QRectF(10, 1, width, height), compText);
      painter.save();
      painter.setCompositionMode(QPainter::CompositionMode_SourceOver);
      painter.drawPixmap(0, 0, pixmap);
      painter.restore();
    }
  }
}

void Z3DCanvas::drawText(QPainter &painter, const QString &text)
{
  if (!text.isEmpty()) {
    QStringList textList;
    textList.append(text);
    drawText(painter, text);
  }
}

void Z3DCanvas::timerEvent(QTimerEvent* e)
{
  broadcastEvent(e, width(), height());
}

void Z3DCanvas::setNetworkEvaluator(Z3DNetworkEvaluator *n)
{
  //m_3dScene->setNetworkEvaluator(n);
  m_networkEvaluator = n;
  if (n)
    n->setOpenGLContext(this);
}

void Z3DCanvas::setFakeStereoOnce()
{
  m_fakeStereoOnce = true;
  //m_3dScene->setFakeStereoOnce();
}

void Z3DCanvas::addEventListenerToBack(Z3DCanvasEventListener *e)
{
  if (e)
    m_listeners.push_back(e);
}

void Z3DCanvas::addEventListenerToFront(Z3DCanvasEventListener *e)
{
  if (e)
    m_listeners.push_front(e);
}

void Z3DCanvas::removeEventListener(Z3DCanvasEventListener *e)
{
  std::deque<Z3DCanvasEventListener*>::iterator pos;
  pos = std::find(m_listeners.begin(), m_listeners.end(), e);

  if (pos != m_listeners.end())
    m_listeners.erase(pos);
}

void Z3DCanvas::clearEventListeners()
{
  m_listeners.clear();
}

void Z3DCanvas::broadcastEvent(QEvent *e, int w, int h)
{
  getGLFocus();
  for (size_t i = 0 ; i < m_listeners.size() ; ++i) {
    m_listeners[i]->onEvent(e, w, h);
    if (e->isAccepted()) {
      break;
    }
  }
}

void Z3DCanvas::setKeyMode(ZInteractionEngine::EKeyMode mode)
{
  m_interaction.setKeyMode(mode);
  viewport()->update();
//  update(QRect(QPoint(0, 0), size()));
}

double Z3DCanvas::getDevicePixelRatio()
{
#ifdef _QT5_
  return (window() && window()->windowHandle()) ?
        window()->windowHandle()->devicePixelRatio() : 1.0;
#else
  return 1.0;
#endif
}

void Z3DCanvas::disableKeyEvent()
{
  m_interaction.setKeyEventEnabled(false);
}

void Z3DCanvas::set3DInteractionHandler(Z3DTrackballInteractionHandler *handler)
{
  m_interaction.set3DInteractionHandler(handler);
}

void Z3DCanvas::updateCursor()
{
  setCursor(m_interaction.getCursorShape());
}

