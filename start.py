import subprocess
import time

for i in range(5):
    print "Running Process " + str(i)
    p = subprocess.Popen("python Main.py -l NAD -p 1ENY", shell=True)
    p.communicate()
    time.sleep(10)

print "Finish" 