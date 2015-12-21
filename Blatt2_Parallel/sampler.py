import subprocess
import sys
from time import sleep

numProcs = int(sys.argv[1])
numSamples = int(sys.argv[2])
deterministic = bool(sys.argv[3])

print "Starting " + str(numSamples) + " samples on " + str(numProcs) + " cores."

procs = []

args= ["./numsim"]

start = 1000.0
end = 2000.0
stepsize = (end-start)/(numSamples-1)

for p in range (0, numProcs):
    if deterministic:
        process = subprocess.Popen(args + [str(start)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        start += stepsize
    else:
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    procs.append(process)
    print "process sarted"

    
finished = 0

while (finished < numSamples):
    for index, proc in enumerate(procs):
        if proc.poll() is not None:
            finished += 1
            print "Processes finished: " + str(finished)
            if deterministic:
                procs[index] = subprocess.Popen(args + [str(start)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                start += stepsize
            else:
                procs[index] = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
    sleep(0.050)
            
print "All processes finished "

for proc in procs:
    proc.kill()
            

    
