import matplotlib
import matplotlib.lines as lines
import pylab as p
import matplotlib.animation as animation
from numpy import *
import sys
import config as cfg
dreload(cfg)
from config import *
import sys
import imp
#p.rc('font', family='serif', serif='cm10', weight='extra bold', size='8')
#p.rc('text', usetex=True)
#p.rc('axes', linewidth=1.0)
#p.rc('axes.grid', linewidth=2.0)
#p.rcParams['text.latex.preamble'] = [r'\usepackage{siunitx}',r'\boldmath',r'\sisetup{detect-weight=true, detect-family=true}']
width_A4 = 21/2.54
width_thesis = 15.55528/2.54
height_A4 = 29.7/2.54
height_thesis = 22.37076/2.54
hx=lx/(nx-1)
hy=ly/(ny-1)
x = linspace(0, lx, nx)
y = linspace(0, ly, ny)
xx,yy = meshgrid(x, y)

numBeg = int(sys.argv[1])
prefixBeg = resumePrefix
switchInd = resumeIndex
numEnd = int(sys.argv[2])
prefixEnd = savePrefix
fps = int(sys.argv[3])
if len(sys.argv) > 4:
	widthpart = double(sys.argv[4])
else:
	widthpart = 1.0

if numBeg < resumeIndex:
	filename2D = 'data/plotData/%sy_%s.bin' % (prefixBeg,numBeg)
else:
	filename2D = 'data/plotData/%sy_%s.bin' % (prefixEnd,numBeg)

#2D-Plot
fignum2D = 1

fig = p.figure(fignum2D, figsize=(8.35, 8), facecolor='w', edgecolor='k')
#fig = p.figure(fignum2D)

p.clf()
ax1 = p.subplot(1,1,1)
#ax1.set_xlabel(r'\bf $x$ in $\si{\milli \metre}$')
#ax1.set_ylabel(r'\bf $y$ in $\si{\milli \metre}$')
#title = r'\bf $t=\SI{%s}{\milli \second}$' % (tp0+numBeg*dtp)
ax1.tick_params(axis='both', which='both', bottom='off', top='off', right='off', left='off', labelbottom='off', labelleft='off')
#ax1.set_aspect('equal')
title = 't=%sms' % (tp0+numBeg*dtp)
ax1title = ax1.set_title(title)
#ax1title.set_text('time_template%(i*dt)')
Y = fromfile(filename2D)
YY = Y.reshape(len(Y)/(nx*ny),ny,nx)
zz = YY[0,:];
quad1 = ax1.pcolormesh(xx, yy, zz)
quad1.set_clim((-120.0,60.0))
fig.tight_layout()
#cbar = p.colorbar(quad1)
#cbar.set_label('U in mV')

def animate(i):
	if i < resumeIndex:
		filename2D = 'data/plotData/%sy_%s.bin' % (prefixBeg,i)
	else:
		filename2D = 'data/plotData/%sy_%s.bin' % (prefixEnd,i)
	Y = fromfile(filename2D)
	YY = Y.reshape(len(Y)/(nx*ny),ny,nx)
	zz = YY[0,:]
	#quad1 = ax1.pcolormesh(xx, yy, zz)
	#quad1.set_clim((-120.0,60.0))
	quad1.set_array(zz[0:999,0:999].ravel())
	#title = r'\bf $t=\SI{%s}{\milli \second}$' % (tp0+i*dtp)
	title = 't=%sms' % (tp0+i*dtp)
	ax1title.set_text(title)
	print i
	return quad1, ax1title

ani = animation.FuncAnimation(fig, animate, arange(numBeg, numEnd+1),
	interval=500, blit=False)
ani.save('test.mp4', fps=fps, bitrate=4096, writer="avconv")
p.ion()
p.show()
fig.tight_layout()
p.figure(fignum2D)
#print zz
#plotFrame(xx,yy,YY[0,:],1)
#1D-Plot
#fignum1D = 2
#fig1D = p.figure(fignum1D)
#p.clf()
#data = fromfile(filename1D)
#N = data[0]
#t = data[1:N+1];
#E = data[N+1:2*N+1]
#axE = fig1D.add_subplot(2,1,1)
#axE.plot(t, E)
#p.xlabel('t in ms')
#p.ylabel('E in mV/cm')
#if obs_Feedback:
#	U = data[2*N+1:3*N+1]
#	axU = fig1D.add_subplot(2,1,2)
#	axU.plot(t, U)
#	p.xlabel('t in ms')
#	p.ylabel('U in mV/cm')
#p.ion()
#p.show()
#p.figure(fignum1D)







