#!/usr/bin/python

import os
import subprocess

files = []
proc = subprocess.Popen(['adb', 'shell', 'ls', 'sdcard/wear*.csv'], stdout=subprocess.PIPE)
res = proc.communicate()[0].decode('ascii')
if not 'No' in res:
    files.extend([n.strip() for n in res.split('\n') if len(n) > 14])

proc = subprocess.Popen(['adb', 'shell', 'ls', 'sdcard/pointlog*.csv'], stdout=subprocess.PIPE)
res = proc.communicate()[0].decode('ascii')
if not 'No' in res:
    files.extend([n.strip() for n in res.split('\n') if len(n) > 14])


for fname in files:
    dst = os.path.join('../external_data/waiting_room', os.path.basename(fname))
    if subprocess.call(['adb', 'pull', fname, dst]) != 0:
        exit()
    subprocess.call(['adb', 'shell', 'rm', fname])
       
history = set()
for fname in files:
    prefix = os.path.basename(fname)[:20]
    if prefix not in history:
        print(prefix)
    history.add(prefix)
