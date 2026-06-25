import weakref


class RegistryMeta(type):
    """Metaclass that adds instance registry to any class."""
    _registry: weakref.WeakValueDictionary
    
    def __new__(mcs, name, bases, namespace):
        cls = super().__new__(mcs, name, bases, namespace)
        # Use weak references to allow garbage collection
        cls._registry = weakref.WeakValueDictionary()
        return cls
