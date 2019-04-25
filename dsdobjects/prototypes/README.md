# dsdobject.prototypes

Sometimes it is difficult to recognize the benefits of base classes, before you
run into problems of using prebuild classes. Sometimes, you don't need more
than the bare backbone of a Domain, Complex, StrandOrder, ..., and you just
want to prototype an algorithm quickly. Well, this is folder is what you are
looking for:

```
    from dsdobjects.prototypes import read_pil
    from dsdobjects.prototypes import LogicDomain, Domain, Complex, ...
```

If you realize that you are missing functionality, feel free to copy the folder
and adapt the files to your needs or, better, conribute some code and make a
pull request!

### Note: 
As of dsdobjects v0.7 this is a preliminary piece of code, expect updates and
issues with backward compatibility.
