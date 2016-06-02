#include "mainwindow.hh"
#include "neteditwidget.hh"
#include "network.hh"
#include "nodes.hh"
#include "assembler.hh"

#include <QTabWidget>
#include <QMenuBar>
#include <QMenu>
#include <QFileDialog>
#include <QToolBar>
#include <QToolButton>
#include <QMessageBox>


MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent), _netedit(0)
{
  setMinimumSize(640, 480);
  setWindowTitle(tr("StochBB"));

  _netedit = new NetEditWidget();
  setCentralWidget(_netedit);

  QMenuBar *menu_bar = new QMenuBar();
  QMenu *file_menu = menu_bar->addMenu(tr("File"));
  QAction *newnet = file_menu->addAction(tr("New"), this, SLOT(onNewNetwork()));
  newnet->setShortcut(Qt::CTRL + Qt::Key_N);
  QAction *open = file_menu->addAction(tr("Open ..."), this, SLOT(onLoad()));
  open->setShortcut(Qt::CTRL + Qt::Key_O);
  QAction *save = file_menu->addAction(tr("Save"), this, SLOT(onSave()));
  save->setShortcut(Qt::CTRL + Qt::Key_S);
  QAction *save_as = file_menu->addAction(tr("Save as ..."), this, SLOT(onSave()));
  save_as->setShortcut(Qt::SHIFT + Qt::CTRL + Qt::Key_S);
  file_menu->addSeparator();
  QAction *check_action = file_menu->addAction(
        QIcon("://icons/check_64.png"), tr("Verify"), this, SLOT(onCheck()));
  QAction *run_action = file_menu->addAction(
        QIcon("://icons/play_64.png"), tr("Run"), this, SLOT(onRun()));
  file_menu->addSeparator();
  QAction *quit = file_menu->addAction(tr("Quit"), this, SLOT(close()));
  quit->setShortcut(Qt::CTRL + Qt::Key_Q);

  QMenu *edit_menu = menu_bar->addMenu(tr("Edit"));
  QMenu *proc_menu = edit_menu->addMenu(QIcon("://icons/proc_64.png"), tr("Add process"));
  proc_menu->addAction(tr("Fixed delay"), _netedit, SLOT(addDelay()));
  proc_menu->addAction(tr("Random delay"), _netedit, SLOT(addRandomDelay()));
  proc_menu->addAction(tr("Gamma process"), _netedit, SLOT(addGammaProc()));
  proc_menu->addAction(tr("compound Gamma process"), _netedit, SLOT(addCompoundGammaProc()));
  proc_menu->addAction(tr("inverse Gamma process"), _netedit, SLOT(addInvGammaProc()));
  proc_menu->addAction(tr("compound inverse Gamma process"), _netedit, SLOT(addCompoundInvGammaProc()));
  QMenu *comb_menu = edit_menu->addMenu(QIcon("://icons/join_64.png"), tr("Add combine"));
  comb_menu->addAction(tr("Minimum"), _netedit, SLOT(addMinimum()));
  comb_menu->addAction(tr("Maximum"), _netedit, SLOT(addMaximum()));
  comb_menu->addAction(tr("Inhibition"), _netedit, SLOT(addInhibition()));
  QMenu *traf_menu = edit_menu->addMenu(QIcon("://icons/trafo_64.png"), tr("Add transform"));
  traf_menu->addAction(tr("Affine"), _netedit, SLOT(addAffine()));
  QMenu *rand_menu = edit_menu->addMenu(QIcon("://icons/var_64.png"), tr("Add variable"));
  rand_menu->addAction(tr("Stimulus"), _netedit, SLOT(addStimulus()));
  rand_menu->addSeparator();
  rand_menu->addAction(tr("Constant"), _netedit, SLOT(addConstant()));
  rand_menu->addAction(tr("Gamma variable"), _netedit, SLOT(addGammaVar()));
  rand_menu->addAction(tr("compound Gamma variable"), _netedit, SLOT(addGammaCompoundVar()));
  QMenu *out_menu = edit_menu->addMenu(QIcon("://icons/output_64.png"), tr("Add output"));
  out_menu->addAction(tr("Marginal plot"), _netedit, SLOT(addMarginalPlot()));

  edit_menu->addSeparator();
  QAction *rm_action = edit_menu->addAction(
        QIcon("://icons/trash_64.png"), tr("Delete selected"), _netedit, SLOT(removeSelected()));
  rm_action->setShortcut(Qt::CTRL + Qt::Key_Backspace);

  QMenu *view_menu = menu_bar->addMenu(tr("View"));
  QAction *zoom_in_action = view_menu->addAction(
        QIcon("://icons/zoom-in_64.png"), tr("Zoom in"), _netedit, SLOT(zoomIn()));
  QAction *zoom_out_action = view_menu->addAction(
        QIcon("://icons/zoom-out_64.png"), tr("Zoom out"), _netedit, SLOT(zoomOut()));

  QMenu *help_menu = menu_bar->addMenu(tr("Help"));
  help_menu->addAction(tr("User Manual"));
  help_menu->addAction(tr("About StochBB"));

  QToolBar *toolbar = new QToolBar();
  toolbar->addAction(check_action);
  toolbar->addAction(run_action);

  toolbar->addSeparator();
  QToolButton* add_proc = new QToolButton();
  add_proc->setText(tr("Add Process"));
  add_proc->setIcon(QIcon("://icons/proc_64.png"));
  add_proc->setToolButtonStyle(Qt::ToolButtonFollowStyle);
  add_proc->setPopupMode(QToolButton::InstantPopup);
  add_proc->setMenu(proc_menu);
  toolbar->addWidget(add_proc);

  QToolButton* add_comb = new QToolButton();
  add_comb->setText(tr("Add Combine"));
  add_comb->setIcon(QIcon("://icons/join_64.png"));
  add_comb->setToolButtonStyle(Qt::ToolButtonFollowStyle);
  add_comb->setPopupMode(QToolButton::InstantPopup);
  add_comb->setMenu(comb_menu);
  toolbar->addWidget(add_comb);

  QToolButton* add_traf = new QToolButton();
  add_traf->setText(tr("Add Transform"));
  add_traf->setIcon(QIcon("://icons/trafo_64.png"));
  add_traf->setToolButtonStyle(Qt::ToolButtonFollowStyle);
  add_traf->setPopupMode(QToolButton::InstantPopup);
  add_traf->setMenu(traf_menu);
  toolbar->addWidget(add_traf);

  QToolButton* add_rand = new QToolButton();
  add_rand->setText(tr("Add Variable"));
  add_rand->setIcon(QIcon("://icons/var_64.png"));
  add_rand->setToolButtonStyle(Qt::ToolButtonFollowStyle);
  add_rand->setPopupMode(QToolButton::InstantPopup);
  add_rand->setMenu(rand_menu);
  toolbar->addWidget(add_rand);

  QToolButton* add_out = new QToolButton();
  add_out->setText(tr("Add Output"));
  add_out->setIcon(QIcon("://icons/output_64.png"));
  add_out->setToolButtonStyle(Qt::ToolButtonFollowStyle);
  add_out->setPopupMode(QToolButton::InstantPopup);
  add_out->setMenu(out_menu);
  toolbar->addWidget(add_out);

  toolbar->addSeparator();
  toolbar->addAction(zoom_in_action);
  toolbar->addAction(zoom_out_action);

  toolbar->addSeparator();
  toolbar->addAction(rm_action);

  addToolBar(toolbar);
}

void
MainWindow::onNewNetwork() {
  _netedit->network()->clear();
  setWindowTitle(tr("StochBB - new"));
}

void
MainWindow::onLoad() {
  QString filename = QFileDialog::getOpenFileName(0, "Load network", "", "*.xml");
  if (filename.isEmpty())
    return;
  if (_netedit->network()->load(filename))
    setWindowTitle(tr("StochBB - %0").arg(_netedit->network()->filename()));
}

void
MainWindow::onSave() {
  if (_netedit->network()->hasFilename())
    _netedit->network()->save();
  else
    onSaveAs();
}

void
MainWindow::onSaveAs() {
  QString filename = QFileDialog::getSaveFileName(0, "Save network", "", "*.xml");
  if (filename.isEmpty())
    return;
  if (_netedit->network()->save(filename))
    setWindowTitle(tr("StochBB - %0").arg(_netedit->network()->filename()));
}

void
MainWindow::onCheck() {
  QHash<Socket *, stochbb::Var> varTable;
  Messages messages;
  if (! Assembler::assemble(_netedit->network(), varTable, messages)) {
    QStringList tmp;
    foreach (Message msg, messages) {
      tmp.append(msg.text());
    }
    QMessageBox::critical(0, tr("Results"), tmp.join("\n"));
  } else {
    QMessageBox::information(0, tr("Success"), "The network is consistent.");
  }
}

void
MainWindow::onRun() {
  QHash<Socket *, stochbb::Var> varTable;
  if (! Assembler::assemble(_netedit->network(), varTable)) {
    QMessageBox::critical(
          0, tr("Can not run network."),
          tr("There was an error during the analysis step: Run 'Check Network'."));
    return;
  }

  // Collect output nodes from network
  Network::nodeIterator node = _netedit->network()->nodesBegin();
  for (; node != _netedit->network()->nodesEnd(); node++) {
    if (OutputNode *out = dynamic_cast<OutputNode *>(*node)) {
      out->execute(varTable);
    }
  }
}
