import numpy as np
import matplotlib.pylab as plt
import sys
if len(sys.argv)==1:
    print "Enter these arguments infront of wfplot.py"
    print "1=number of files,2=file names"
    exit()
if sys.argv[1] in ("--help"):
    print "Enter these arguments infront of wfplot.py"
    print "1=number of files,2=file names"
    exit()
for i in range(0,int(sys.argv[1])):
   a=np.loadtxt(sys.argv[i+2])
   plt.plot(a[:,0],a[:,1],label=sys.argv[i+2])
   plt.legend()
plt.xlabel('r',fontsize=20)
plt.show()
plt.savefig('wfplot')
