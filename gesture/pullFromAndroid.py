#!/usr/bin/python

import os
import subprocess

proc = subprocess.Popen(['adb', 'shell', 'ls', 'sdcard/wear*.csv'], stdout=subprocess.PIPE)
res = proc.communicate()[0].decode('ascii')
files = [n.strip() for n in res.split('\n') if len(n) > 5]

for fname in files:
    dst = os.path.join('data/waiting_room', os.path.basename(fname))
    subprocess.call(['adb', 'pull', fname, dst])
    subprocess.call(['adb', 'shell', 'rm', fname])
       
history = set()
for fname in files:
    prefix = os.path.basename(fname)[:20]
    if name not in history:
        print(prefix)
    history.add(prefix)
