import inspect
import testmethod
import ldinvert
import ldinvert2
import ldscore
import truth

allmethods = inspect.getmembers(testmethod, inspect.isclass) + \
    inspect.getmembers(ldinvert, inspect.isclass) + \
    inspect.getmembers(truth, inspect.isclass) + \
    inspect.getmembers(ldscore, inspect.isclass) + \
    inspect.getmembers(ldinvert2, inspect.isclass)

def find_method(name):
    for n, m in allmethods:
        if n == name:
            return m
    return None
