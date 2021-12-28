import os
import sys

from 

SHIFTX2 = ['shiftx2.py']
def find_executable(names):
    for possible in names:
        result = _find_executable(possible)
        if result is not None:
            return result
    return None

print(SHIFTX2_HOME)