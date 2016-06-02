/** \page mem Memory Management
 * The classes in StochBB are separated into two groups API are all classes intended to be
 * used by other software using libstochbb, while all remainnig classes should
 * only be used within libstochbb. Dividing them into two groups provides two advantages: (a) The
 * API is well defined and clean. It does not contain any interface that is intended to be
 * used by thrid-party-software. And (b), it allows to implement everything in terms of actual
 * objects and containers. The "objects" do the voodoo while containers provide the intended API and
 * manage "objects". This allows to implement some automatc memmory management and the user does not
 * need to deal with it.

 * Almost all classes of the API are derived from the @c Container class. This class implements a
 * simple reference counting which ensures that any @c Object being held in a container by the
 * user will never be freed. As objects may reference other object, a simple mark and sweep garbage
 * collector is implemented too, which searches for unreachable objects that are not held in a
 * container and frees them.
 *
 * Therefore, when extending libstochbb, two basic rules must be considered:
 *   -# Never pass around an @c Object directly. Only pass around @c Container. This ensures that
 *      the object referenced by the Container is freed although being used.
 *   -# Never store containers in an @c Object. This rule ensures that there are no circular
 *      references which will result in a memory leak.
 *
 * Following these rules will ensure that now @c Object is freed while beeing used and all
 * unreachable Objects are freed.
 */

/** @defgroup internal Internal used objects and functions.
 * This group collects all classes and functions that are not directly accessible through the API.
 * These objects are derived from the basic @c stochbb::Object class wich implements the reference
 * counting and interfaces for the mark & sweep grabage collector. All classes documented here are only of
 * interest if you want to understand the implementation of StochBB or if you want to extend it. */

#ifndef __SBB_OBJECT_HH__
#define __SBB_OBJECT_HH__

#include <unordered_set>
#include <unordered_map>
#include <cstddef>
#include <typeinfo>



namespace stochbb {

// Forward declaration
class Object;

/** Garbage collector.
 * @ingroup internal */
class GC
{
protected:
  /** Hidden constructor. Use the factory method @c get to obtain the instance. */
  GC();

public:
  /** Destructor. */
  virtual ~GC();
  /** Runs the garbage collector. */
  void run();
  /** Factory method. */
  static GC &get();

protected:
  /** Registers an object with the GC. */
  inline void add(Object *obj) { _objects.insert(obj); }

protected:
  /** The set of all objects. */
  std::unordered_set<Object *> _objects;
  /** The singleton instance. */
  static GC *_instance;

  friend class Object;
};


/** Base class of all managed objects.
 * @ingroup internal */
class Object
{
protected:
  /** Hidden constructor. */
  Object();

public:
  /** Destructor. */
  virtual ~Object();
  /** Returns true if the object is marked by the GC. */
  inline bool isMarked() const { return _marked; }
  /** Marks the object (and all objects referenced). */
  virtual void mark();
  /** Unmark the object. */
  void unmark();
  /** Increments the reference counter. */
  inline void ref() { _refcount++; }
  /** Decrement the reference counter. */
  inline void unref() { if (_refcount) _refcount--; if (!_refcount) GC::get().run(); }
  /** Retruns the reference count. */
  inline size_t refcount() const { return _refcount; }
  /** Prints a textual representation of the object. */
  virtual void print(std::ostream &stream) const;

protected:
  /** If @c true, the object is marked. */
  bool _marked;
  /** The reference counter. */
  size_t _refcount;
};

}

#endif
