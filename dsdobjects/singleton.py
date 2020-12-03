import logging
log = logging.getLogger(__name__)

import sys
from weakref import WeakValueDictionary

class SingletonError(Exception):
    def __init__(self, message, existing = None):
        super(SingletonError, self).__init__(message)
        self.existing = existing

def clear_singletons(cls):
    """ Hard reset of singleton dictionary.

    This metod is for debugging only, using this in production 
    might cause unexpected behavior or segfaults... well, well.
    """
    cls._instanceNames.clear()
    cls._instanceCanon.clear()
    cls._instanceNames = WeakValueDictionary() # [name] = canon
    cls._instanceCanon = WeakValueDictionary() # [canon] = obj

def show_singletons(cls):
    """ Yields data of all registered singleton objects of the given class.

    Yields:
        string: name, repr(obj), id(obj), refcount(obj)
    """
    for name, obj in cls._instanceNames.items():
        yield (f'name = {name}, obj = {repr(obj)}, id = {id(obj)}, sys.getrefcount(obj) = {sys.getrefcount(obj)}')

class Singleton(type):
    """ A singleton metaclass. 
    
    All classes in this metaclass have the property that only one instance is
    created per canonical form. The classes using the Singleton type must
    provide a @staticmethod to derive the canonical form from the arguments.
    Because some classes may be using a "name" that is independent from the
    canonical form, this Singleton Metaclass also checks that the names are
    unique for each canonical form. The objects using the singleton type 
    should ensure that both name and canonical form are IMMUTABLE. 

    Examples for Domains:
        a = Domain('a', 10)
        assert a is Domain('a')
        del a
        assert Domain('a', 15) is ~Domain('a*', 15)

    Returns:
        obj: Either the existing or a new object.
    """
    def __init__(cls, name, bases, dict):
        super(Singleton, cls).__init__(cls, bases, dict)
        cls._instanceNames = WeakValueDictionary() # [name] = canon
        cls._instanceCanon = WeakValueDictionary() # [canon] = obj

    def __call__(cls, *args, **kwargs):
        """ Returns the object as defined by the agruments. 

        Arguments may be such that a canonical form can be derived, or they
        provide a name to find the object in the instance cache. In order to
        use names for lookup, the cache has to be initialized by providing both
        the the canonical form and the corresponding name. 

        Raises:
            SingletonError: If aruments are insufficient. 
            SingletonError: If name and canonical form point to different objects.
        """
        # Note that cls.identifiers() may return extra arguments that should be
        # passed on the object initialization. That may look odd, but sometimes
        # finding the canonical form is expensive and one would like to provide
        # a new object with the aquired data. 
        canon, name, kwadd = cls.identifiers(*args, **kwargs)
        assert not any((arg in kwargs and kwargs[arg] is not None) for arg in kwadd.keys())
        kwargs.update(kwadd)

        Sobj = None
        if name and canon:
            if name in cls._instanceNames or canon in cls._instanceCanon:
                objN = cls._instanceNames.get(name, None)
                objC = cls._instanceCanon.get(canon, None)
                if objN is None:
                    raise SingletonError(f'Duplicate Singleton {cls.__name__}({name} vs. {objC.name}).', existing = objC)
                elif objC is None:
                    # Let's not return exisiting here, although we could.  It
                    # might not be clear that this can happen when there are
                    # automated naming problems. E.g. one expects the Error
                    # because the canonical form object exists already, but instead
                    # you get a different complex back just because of the name ...
                    raise SingletonError(f'Duplicate Singleton {cls.__name__}({name}).')
                if not (objN is objC):
                    raise SingletonError(f'Duplicate Singleton {cls.__name__}: name ({name}) and canonical form match different objects!')
                Sobj = objN
        elif name:
            if name not in cls._instanceNames:
                raise SingletonError(f'Cannot instantiate Singleton {cls.__name__} from: name {name} only.')
            Sobj = cls._instanceNames[name]
        else:
            if canon not in cls._instanceCanon:
                raise SingletonError(f'Cannot instantiate Singleton {cls.__name__} from: canonical form only.')
            Sobj = cls._instanceCanon[canon]

        if Sobj is None:
            Sobj = super(Singleton, cls).__call__(*args, **kwargs)
            cls._instanceNames[name] = Sobj
            cls._instanceCanon[canon] = Sobj
        else:
            log.debug(f'Returning exisiting Singleton {cls.__name__} with name {name}: {canon}!')
        return Sobj

