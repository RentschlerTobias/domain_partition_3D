#!/usr/bin/env python3
import sys
sys.path.insert(0, '/dtOO-install/tools')
import dtOOPythonSWIG as dtOO

print("dtOO attributes:")
for attr in dir(dtOO):
    if not attr.startswith('_'):
        print(f"  {attr}")
