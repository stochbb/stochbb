#ifndef __SBB_OBJECT_HH__
#define __SBB_OBJECT_HH__

#include <unordered_set>
#include <unordered_map>
#include <cstddef>
#include <typeinfo>



namespace sbb {

// Forward declaration
class Object;

/** Garbage collector. */
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
  /** Box an object. This method implements the reference counting. */
  void box(Object *obj);
  /** Unbox an object. This method implements the reference counting. */
  void unbox(Object *obj);
  /** Adds an object to the GC. */
  void add(Object *obj);
  /** Lock the GC, this prevents the collection of unreachable objects until @c unlock()
   * is called. */
  void lock();
  /** Unlocks the GC. */
  void unlock();

protected:
  /** Removes an object from the GC. */
  inline void remove(Object *obj) { _objects.erase(obj); }

protected:
  /** The set of all objects. */
  std::unordered_set<Object *> _objects;
  /** The set of objects which are referenced by containers. */
  std::unordered_map<Object *, size_t> _boxed;
  /** The counter for blocking the GC. */
  size_t _lockCount;
  /** The singleton instance. */
  static GC *_instance;

  friend class Object;
};


/** Base class of all managed objects. */
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


/** Base class of all container classes.
 * @ingroup api */
class Container
{
public:
  /** The object type of the container. */
  typedef Object ObjectType;

protected:
  /** Hidden constructor. */
  Container();
  /** Packs the given objects. */
  explicit Container(Object *obj);
  /** Copy constructor. */
  Container(const Container &other);

public:
  /** Destructor. */
  virtual ~Container();

  /** Assignment operator. */
  const Container &operator=(const Container &other);

  /** Returns @c true if the container is empty. */
  bool isNull() const;

  /** Returns @c true if the object can be casted to the given container. */
  template <class T>
  bool is() const {
    return 0 != dynamic_cast<typename T::ObjectType *>(this->_object);
  }

  /** Casts the object to the given container. */
  template <class T>
  T as() const {
    return T(dynamic_cast<typename T::ObjectType *>(this->_object));
  }

protected:
  /** Boxes the given object. */
  void box(Object *obj);

protected:
  /** Holds the object. */
  Object *_object;
};

}

#endif
