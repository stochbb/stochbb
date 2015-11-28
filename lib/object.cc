#include "object.hh"
#include <list>
#include <iostream>
#include <typeinfo>

using namespace sbb;


/* ********************************************************************************************* *
 * Implementation of GC
 * ********************************************************************************************* */
GC *GC::_instance = 0;

GC::GC()
  : _objects(), _boxed(), _lockCount(0)
{
  // pass...
}

GC::~GC() {
  // pass...
}

void
GC::run()
{
  if (0 != _lockCount) { return; }
  _lockCount++;  // <- avoid recursive calls of run()

  // Mark all boxed objects:
  std::unordered_map<Object *, size_t>::iterator item = _boxed.begin();
  for (; item!=_boxed.end(); item++) {
    item->first->mark();
  }
  // Delete all unmarked objects:
  std::list<Object *> del;
  for (std::unordered_set<Object *>::iterator item=_objects.begin(); item!=_objects.end(); ) {
    if ((*item)->isMarked()) {
      (*item)->unmark(); item++;
    } else {
      Object *obj = *item;
      item = _objects.erase(item);
      delete obj;
    }
  }
  // Free objects
  for (std::list<Object *>::iterator item=del.begin(); item!=del.end(); item++) {
    _objects.erase(*item);
    delete *item;
  }
  _lockCount--;
}

void
GC::box(Object *obj) {
  std::unordered_map<Object *, size_t>::iterator item = _boxed.find(obj);
  if (_boxed.end() != item) {
    // If the object is already boxed -> increment counter
    //std::cerr << "GC: Increment counter for " << obj << " " << item->second;
    item->second += 1;
    //std::cerr << " -> " << item->second << std::endl;
  } else {
    //std::cerr << "GC: Box " << obj << " (" << typeid(obj).name() << ")" << std::endl;
    // Other wise add to table of boxed objects
    _boxed[obj] = 1;
  }
}

void
GC::unbox(Object *obj) {
  std::unordered_map<Object *, size_t>::iterator item = _boxed.find(obj);
  // If object is not boxed -> skip
  if (_boxed.end() == item) {
    std::cerr << "GC: Warning: Unbox unboxed object: " <<  obj << std::endl;
    return;
  }
  if (0 == item->second) {
    std::cerr << "GC: Warning: Object already unboxed!" <<  obj << std::endl;
    return;
  }
  // Decrease counter
  //std::cerr << "GC: Decrement counter for " << item->first << " " << item->second;
  item->second--;
  //std::cerr << " -> " << item->second << std::endl;
  // If no container refers to this object anymore -> add to set of unboxed objects
  if (0 == item->second) {
    //std::cerr << "GC: Unbox " << item->first << std::endl;
    _boxed.erase(item);
    // Collect unreachables
    run();
  }
}

void
GC::add(Object *obj) {
  if (0 == obj) { return; }
  _objects.insert(obj);
}

void
GC::lock() {
  _lockCount++;
}

void
GC::unlock() {
  _lockCount--;
  if (0 == _lockCount) { run(); }
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
  GC::get().add(this);
}

Object::~Object() {
  GC::get().remove(this);
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


/* ********************************************************************************************* *
 * Implementation of Container
 * ********************************************************************************************* */
Container::Container()
  : _object(0)
{
  // pass...
}

Container::Container(Object *obj)
  : _object(obj)
{
  // Register object with GC
  if (0 != _object) {
    GC::get().box(_object);
  }
}

Container::Container(const Container &other)
  : _object(other._object)
{
  // Register object with GC
  if (0 != _object) {
    GC::get().box(_object);
  }
}

Container::~Container() {
  // Unregister object with GC
  if (0 != _object) {
    GC::get().unbox(_object);
  }
  _object = 0;
}

const Container &
Container::operator =(const Container &other)
{
  // Unregister current object:
  if (0 != _object) { GC::get().unbox(_object); }
  _object = other._object;
  // Register new object
  if (0 != _object) { GC::get().box(_object); }
  // done.
  return *this;
}

bool
Container::isNull() const {
  return 0 == _object;
}

void
Container::box(Object *obj) {
  if (0 != _object) { GC::get().unbox(_object); }
  _object = obj;
  if (0 != _object) { GC::get().box(_object); }
}
