#!/usr/bin/python
import pandas as pd
import sys

print pd.read_pickle(sys.argv[1])
