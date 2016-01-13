import inspect
import testmethod
import ldinvert

allmethods = inspect.getmembers(testmethod, inspect.isclass) + \
        inspect.getmembers(ldinvert, inspect.isclass)

def find_method(name):
    for n, m in allmethods:
        if n == name:
            return m
    return None
