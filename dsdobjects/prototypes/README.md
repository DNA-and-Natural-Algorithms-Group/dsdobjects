# dsdobject.prototypes

Sometimes it is difficult to recognize the benefits of base classes, before you
run into problems of using prebuild classes. Sometimes, you don't need more
than the bare backbone of a Domain, Complex, StrandOrder, ..., and you just
want to prototype an algorithm quickly. Well, this is folder is what you are
looking for:

```
    from dsdobjects.prototypes import init_from_pil
    from dsdobjects.prototypes import LogicDomain, Domain, Complex, ...
```

If you realize that you are missing functionality, feel free to copy the folder
and adapt the files to your needs.
