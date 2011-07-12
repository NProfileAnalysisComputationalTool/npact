import sys
import os.path
import os

#http://docs.python.org/install/index.html#inst-alt-install

dirpath = os.path.dirname(os.path.realpath(__file__))

for d in os.listdir(dirpath) :
    p = os.path.join(dirpath,d)
    if os.path.isdir(p) :
        sys.path.append(p)
