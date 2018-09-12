#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

with open(sys.argv[1]) as fin:
    for i, line in enumerate(fin):
        if i % 4 == 1:
            line = line.upper()
        sys.stdout.write(line.upper())
