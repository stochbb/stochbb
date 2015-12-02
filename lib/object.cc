#include "object.hh"
#include <list>
#include <iostream>
#include <typeinfo>
#include "logger.hh"

using namespace sbb;


/* ********************************************************************************************* *
 * Implementation of GC
 * ********************************************************************************************* */
GC *GC::_instance = 0;

GC::GC()
  : _objects()
{
  // pass...
}

GC::~GC() {
  // pass...
}

void
GC::run()
{
  // Mark all boxed objects:
  std::unordered_map<Object *, size_t>::iterator item = _objects.begin();
  for (; item!=_objects.end(); item++) {
    if (item->second)
      item->first->mark();
  }
  // Delete all unmarked objects:
  std::list<Object *> del;
  for (std::unordered_map<Object *, size_t>::iterator item=_objects.begin(); item!=_objects.end(); ) {
    if (item->first->isMarked()) {
      item->first->unmark(); item++;
    } else {
      Object *obj = item->first;
      item = _objects.erase(item);
      delete obj;
    }
  }
  // Free objects
  for (std::list<Object *>::iterator item=del.begin(); item!=del.end(); item++) {
    delete *item;
  }
}

GC&
GC::get() {
  if (0 == GC::_instance) {
    GC::_instance = new GC();
  }
  return *(GC::_instance);
}


/* ********************************************************************************************* *
 * Implementation of Object
 * ********************************************************************************************* */
Object::Object()
  : _marked(false)
{
  GC::get().box(this);
}

Object::~Object() {
  // pass...
}

void
Object::mark() {
  if (_marked) { return; }
  _marked = true;
}

void
Object::unmark() {
  _marked = false;
}
