#!/usr/bin/env python

# http://stackoverflow.com/questions/1729824/transpose-a-file-in-bash/1730698#1730698

import sys

if sys.stdin.isatty():
    msg = "error: no stdin"
    msg += "\nusage: cat myfile.txt | transpose.py > myfiletransposed.txt"
    msg += "\nexample: echo -e \"1 a\\n2 b\\n3 c\"; echo; echo -e \"1 a\\n2 b\\n3 c\" | transpose.py"
    sys.stderr.write("%s\n" % msg)
    sys.exit(1)

for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip())):
    print(' '.join(c))
