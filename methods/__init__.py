import inspect
import testmethod

allmethods = inspect.getmembers(testmethod, inspect.isclass)

def find_method(name):
    for n, m in allmethods:
        if n == name:
            return m
    return None
