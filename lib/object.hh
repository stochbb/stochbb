/** @defgroup internal Internal used objects
 */
#ifndef __SBB_OBJECT_HH__
#define __SBB_OBJECT_HH__

#include <unordered_set>
#include <unordered_map>
#include <cstddef>
#include <typeinfo>



namespace sbb {

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
  /** Registers an object with the GC. */
  inline void add(Object *obj) { _objects.insert(obj); }

protected:
  /** The set of all objects. */
  std::unordered_map<Object *, size_t> _objects;
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

protected:
  /** If @c true, the object is marked. */
  bool _marked;
};

}

#endif
