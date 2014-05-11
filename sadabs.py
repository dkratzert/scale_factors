#-*- encoding: utf-8 -*-
#möp
#
#
from __future__ import print_function

import time
import subprocess

sad = subprocess.Popen('sadabs.exe',stdin=subprocess.PIPE, stdout=subprocess.PIPE)
time.sleep(0.5)
time.sleep(0.5)
print(Popen.poll())
#sad.communicate('\n')
#time.sleep(0.5)
#sad.communicate('\n')
#time.sleep(0.5)
#sad.communicate('\n')
#time.sleep(0.5)
sad.terminate()
sad.kill()
