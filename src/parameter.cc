#include "parameter.hh"
#include "nodes.hh"
#include <QTextStream>
#include <QCheckBox>
#include <QLineEdit>
#include <QTextEdit>
#include <QFormLayout>
#include <QVBoxLayout>
#include <QIntValidator>
#include <QDoubleValidator>
#include <QDialogButtonBox>
#include <QLabel>



/* ********************************************************************************************* *
 * Implementation of ParserInfo
 * ********************************************************************************************* */
ParserInfo::ParserInfo()
  : _state(OK), _messages()
{
  // pass...
}

ParserInfo::State
ParserInfo::state() const {
  return _state;
}

QStringList
ParserInfo::messages() const {
  QStringList msgs;
  QList< QPair<State, QString> >::const_iterator item = _messages.begin();
  for(; item != _messages.end(); item++) {
    switch (item->first) {
      case OK:
        msgs.push_back("OK:   " + item->second);
        break;
      case WARNING:
        msgs.push_back("WARN: " + item->second);
        break;
      case ERROR:
        msgs.push_back("ERR:  " + item->second);
        break;
    }
  }
  return msgs;
}

void
ParserInfo::addMessage(State state, const QString &msg) {
  _state = std::max(_state, state);
  _messages.append(QPair<State, QString>(state, msg));
}

void
ParserInfo::addInfo(const QString &msg) {
  addMessage(OK, msg);
}

void
ParserInfo::addWarning(const QString &msg) {
  addMessage(WARNING, msg);
}

void
ParserInfo::addError(const QString &msg) {
  addMessage(ERROR, msg);
}


/* ********************************************************************************************* *
 * Implementation of Parameter
 * ********************************************************************************************* */
Parameter::Parameter()
  : _type(BOOL), _value(false)
{
  // pass...
}

Parameter::Parameter(bool val)
  : _type(BOOL), _value(val)
{
  // pass...
}

Parameter::Parameter(int val)
  : _type(INTEGER), _value(val)
{
  // pass...
}

Parameter::Parameter(double val)
  : _type(FLOAT), _value(val)
{
  // pass...
}

Parameter::Parameter(const QString &val)
  : _type(STRING), _value(val)
{
  // pass...
}

Parameter::Parameter(const Parameter &other)
  : _type(other._type), _value(other._value)
{
  // pass...
}

Parameter &
Parameter::operator =(const Parameter &other) {
  _type = other._type;
  _value = other._value;
  return *this;
}

Parameter::Type
Parameter::type() const {
  return _type;
}

bool
Parameter::isBool() const {
  return BOOL == _type;
}

bool
Parameter::asBool() const {
  return _value.toBool();
}

bool
Parameter::isInt() const {
  return INTEGER == _type;
}

int
Parameter::asInt() const {
  return _value.toInt();
}

bool
Parameter::isFloat() const {
  return FLOAT == _type;
}

double
Parameter::asFloat() const {
  return _value.toDouble();
}

bool
Parameter::isString() const {
  return STRING == _type;
}

QString
Parameter::asString() const {
  return _value.toString();
}

QDomElement
Parameter::serialize(QDomDocument &doc) const {
  QDomElement node = doc.createElement("parameter");
  switch (_type) {
    case BOOL:
      node.setAttribute("type", "bool");
      break;
    case INTEGER:
      node.setAttribute("type", "int");
      break;
    case FLOAT:
      node.setAttribute("type", "float");
      break;
    case STRING:
      node.setAttribute("type", "string");
      break;
  }
  node.appendChild(doc.createTextNode(_value.toString()));

  return node;
}

Parameter
Parameter::fromXml(const QDomElement &node, ParserInfo &info) {
  if (! node.hasAttribute("type")) {
    QString text; QTextStream msg(&text);
    msg << "@ line " << node.lineNumber() << ": Parameter node has no 'type' attribute.";
    info.addWarning(text);
    return Parameter();
  }

  if ("bool" == node.attribute("type")) {
    if ("true" == node.text())
      return Parameter(true);
    return Parameter(false);
  } else if ("int" == node.attribute("type")) {
    return Parameter(node.text().toInt());
  } else if ("float" == node.attribute("type")) {
    return Parameter(node.text().toDouble());
  } else if ("string" == node.attribute("type")) {
    return Parameter(node.text());
  }

  QString text; QTextStream msg(&text);
  msg << "@ line " << node.lineNumber() << ": Parameter node has unknown type '"
      << node.attribute("type") << ".";
  info.addWarning(text);
  return Parameter();
}


/* ********************************************************************************************* *
 * Implementation of ParameterEdit
 * ********************************************************************************************* */
ParameterEdit::ParameterEdit(const Parameter &param, QWidget *parent)
  : QWidget(parent), _parameter(param), _edit(0)
{
  if (Parameter::BOOL == _parameter.type()) {
    QCheckBox *obj = new QCheckBox();
    obj->setChecked(_parameter.asBool());
    _edit = obj;
  } else if (Parameter::INTEGER == _parameter.type()) {
    QLineEdit *obj = new QLineEdit(QString::number(_parameter.asInt()));
    obj->setValidator(new QIntValidator());
    _edit = obj;
  } else if (Parameter::FLOAT == _parameter.type()) {
    QLineEdit *obj = new QLineEdit(QString::number(_parameter.asFloat()));
    obj->setValidator(new QDoubleValidator());
    _edit = obj;
  } else if (Parameter::STRING == _parameter.type()) {
    _edit = new QLineEdit(_parameter.asString());
  }

  QVBoxLayout *layout = new QVBoxLayout();
  layout->addWidget(_edit);
  layout->setSpacing(0);
  layout->setMargin(0);
  setLayout(layout);
}

Parameter
ParameterEdit::parameter() const {
  if (Parameter::BOOL == _parameter.type()) {
    QCheckBox *obj = dynamic_cast<QCheckBox *>(_edit);
    return Parameter(obj->isChecked());
  } else if (Parameter::INTEGER == _parameter.type()) {
    QLineEdit *obj = dynamic_cast<QLineEdit *>(_edit);
    return Parameter(obj->text().toInt());
  } else if (Parameter::FLOAT == _parameter.type()) {
    QLineEdit *obj = dynamic_cast<QLineEdit *>(_edit);
    return Parameter(obj->text().toDouble());
  } else if (Parameter::STRING == _parameter.type()) {
    QLineEdit *obj = dynamic_cast<QLineEdit *>(_edit);
    return Parameter(obj->text());
  }

  return Parameter();
}


/* ********************************************************************************************* *
 * Implementation of NodeConfigDialog
 * ********************************************************************************************* */
NodeConfigDialog::NodeConfigDialog(NodeBase *node)
  : QDialog(), _node(node)
{
  QFormLayout *nodeprops = new QFormLayout();
  nodeprops->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);
  _label = new QLineEdit(_node->label());
  _label->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Preferred);
  nodeprops->addRow(tr("Label"), _label);

  _description = new QTextEdit(_node->description());
  nodeprops->addRow(tr("Description"), _description);

  QFormLayout *param_layout = 0;
  if (node->hasParameters()) {
    param_layout = new QFormLayout();
    param_layout->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);
  }
  const QHash<QString, Parameter> &params = node->parameters();
  QHash<QString, Parameter>::const_iterator param = params.begin();
  for (; param != params.end(); param++) {
    ParameterEdit *p = new ParameterEdit(param.value());
    _params.insert(param.key(), p);
    param_layout->addRow(param.key(), p);
  }

  QVBoxLayout *layout = new QVBoxLayout();
  layout->addWidget(new QLabel(tr("<h2>Configure <i>%1</i> node.</h2>").arg(node->type())));
  layout->addSpacing(15);
  layout->addLayout(nodeprops);
  if (param_layout) {
    layout->addSpacing(15);
    layout->addWidget(new QLabel(tr("<h3>Parameters:</h3>")));
    layout->addLayout(param_layout);
  }

  QDialogButtonBox *bb = new QDialogButtonBox(QDialogButtonBox::Cancel | QDialogButtonBox::Ok);
  connect(bb, SIGNAL(accepted()), this, SLOT(apply()));
  connect(bb, SIGNAL(rejected()), this, SLOT(reject()));
  layout->addWidget(bb);

  setLayout(layout);
}

void
NodeConfigDialog::apply() {
  // Update label
  _node->setLabel(_label->text());
  // update description
  _node->setDescription(_description->toPlainText());
  // update parameters
  QHash<QString, ParameterEdit *>::iterator item = _params.begin();
  for (; item != _params.end(); item++) {
    _node->setParameter(item.key(), item.value()->parameter());
  }

  accept();
}

