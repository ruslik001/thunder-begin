#!/usr/bin/env python

# a bar plot with errorbars
import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy import loadtxt
import matplotlib.ticker as mtick
from matplotlib.ticker import AutoMinorLocator
import numpy as np

ele = sys.argv[1]
states = sys.argv[2]

s_x = None
p_x = None
d_x = None

mini=100
maxi=-100
maxx = -100

if 's' in states:
	s = loadtxt(ele+".wf-s0.dat")
	sext = loadtxt(ele+".wf-s1.dat")

	s_x = s[:,0]
	s_y = s[:,1]

	sext_x = sext[:,0]
	sext_y = sext[:,1]

	if maxi<np.max(np.concatenate((s_y, sext_y))):
                maxi=np.max(np.concatenate((s_y, sext_y)))
        if mini>np.min(np.concatenate((s_y, sext_y))):
		mini=np.min(np.concatenate((s_y, sext_y)))
        if maxx<np.max(np.concatenate((s_x, sext_x))):
		maxx=np.max(np.concatenate((s_x, sext_x)))

if 'p' in states:
	p = loadtxt(ele+".wf-p0.dat")
	pext = loadtxt(ele+".wf-p1.dat")

	p_x = p[:,0]
	p_y = p[:,1]

	pext_x = pext[:,0]
	pext_y = pext[:,1]


	if maxi<np.max(np.concatenate((p_y, pext_y))):
                maxi=np.max(np.concatenate((p_y, pext_y)))
        if mini>np.min(np.concatenate((p_y, pext_y))):
		mini=np.min(np.concatenate((p_y, pext_y)))
        if maxx<np.max(np.concatenate((p_x, pext_x))):
		maxx=np.max(np.concatenate((p_x, pext_x)))

if 'd' in states:
	d = loadtxt(ele+".wf-d0.dat")
	dext = loadtxt(ele+".wf-d1.dat")

	d_x = d[:,0]
	d_y = d[:,1]

	dext_x = dext[:,0]
	dext_y = dext[:,1]


	if maxi<np.max(np.concatenate((d_y, dext_y))):
                maxi=np.max(np.concatenate((d_y, dext_y)))
        if mini>np.min(np.concatenate((d_y, dext_y))):
		mini=np.min(np.concatenate((d_y, dext_y)))
        if maxx<np.max(np.concatenate((d_x, dext_x))):
		maxx=np.max(np.concatenate((d_x, dext_x)))

colors = ['#000000', '#ff0000', '#00ff00', '#0000ff', '#8c0733', '#460319', '#afc94c', '#51caf2', '#376180', '#660066' ]

fig = plt.figure()

ax = fig.add_subplot(111)


if s_x is not None:
	rects0 = plt.plot(s_x, s_y, color=colors[0], lw =2, label='s-state')
	rects2 = plt.plot(sext_x, sext_y, color=colors[2], lw =2, label='s-state (excited)')
 

if p_x is not None:
	rects1 = plt.plot(p_x, p_y, color=colors[1], lw =2, label='p-state')
	rects3 = plt.plot(pext_x, pext_y, color=colors[3], lw =2, label='p-state (excited)')

if d_x is not None:

	rects4 = plt.plot(d_x, d_y, color=colors[4], lw =2, label='d-state')
	rects5 = plt.plot(dext_x, dext_y, color=colors[5], lw =2, label='d-state (excited)')

plt.legend(loc=1)

plt.xlabel(r'$\mathbf{Radial}$'+" "+r'$\mathbf{Distance}$'+" "+r'$\mathbf{(a_{B})}$',fontsize=20)
plt.ylabel('Amplitude (arbitrary units)',fontsize=20)

minor_locator_x = AutoMinorLocator(4)
minor_locator_y = AutoMinorLocator(2)
ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
ax.xaxis.set_minor_locator(minor_locator_x)
ax.yaxis.set_minor_locator(minor_locator_y)

plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4)
plt.tight_layout()
rng=maxi-mini
rngx=maxx - 0 

plt.ylim([mini-.1*rng,maxi+.1*rng])
plt.xlim([0,maxx+.05*rngx])
plt.savefig(ele+'-states.png')  

print maxi, mini

maxi=-100
plt.show()

