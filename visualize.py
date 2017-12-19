import matplotlib
#matplotlib.use("GTKAGG")
import pylab as p
from numpy import *
import sys
import config as cfg
dreload(cfg)
from config import *
#dreload(config_PB_160)
#from config_PB_160 import *
#import config
#dreload(config)
#from config import *
p.rc('font', family='serif', serif='cm10', weight='extra bold', size='16')
p.rc('text', usetex=True)
#p.rc('axes', linewidth=1.0)
#p.rc('axes.grid', linewidth=2.0)
p.rcParams['text.latex.preamble'] = [r'\usepackage{siunitx}',r'\boldmath',r'\sisetup{detect-weight=true, detect-family=true}']
#p.rcParams['image.cmap'] = 'jet'

def plotFrame(xx,yy,zz,fig=1):
	p.figure(fig)
	p.clf()
	p.pcolormesh(xx, yy, zz)
	p.xlabel('x in mm')
	p.ylabel('y in mm')
	p.title('test')
	p.clim(-0.5*85.7-84,1.4*85.7-84)
	cbar=p.colorbar()
	cbar.set_label('U in mV')
	p.ion()
	p.show()
	p.figure(fig)
	return 0

hx=lx/(nx-1.0)
hy=ly/(ny-1.0)
x = linspace(0, lx, nx)
y = linspace(0, ly, ny)
xx,yy = meshgrid(x, y)

if len(sys.argv) > 1:
	num = int(sys.argv[1])
	if len(sys.argv) > 2:
		prefix = sys.argv[2]
	elif num < resumeIndex:
		prefix = resumePrefix
	else:
		prefix = savePrefix	
	filename2D = 'data/plotData/%sy_%s.bin' % (prefix,num)
else:
	import os
	import glob
	import re
	filename2D = max(glob.iglob('data/plotData/*y_*.bin'), key=os.path.getctime)
	nums=re.findall('([0-9]+)', filename2D)
	num = int(nums[len(nums)-1]);
	prefix = savePrefix
filename1D = 'data/plotData/%s1DPlot.bin' % (prefix)
filenameHet = 'data/plotData/%shet.bin' % (prefix)

#2D-Plot
fignum2D = 1
fig2D = p.figure(fignum2D)
p.clf()
Y = fromfile(filename2D)
print Y.shape
YY = Y.reshape(len(Y)/(nx*ny),ny,nx)
zz = YY[0,:];
title = 'Plot%s t=%sms' % (num,tp0+num*dtp)
p.pcolormesh(xx, yy, zz)
p.xlabel('x in mm')
p.ylabel('y in mm')
p.title(title)
p.clim(-120.0,60.0)
cbar=p.colorbar()
cbar.set_label('U in mV')
p.ion()
p.show()
p.figure(fignum2D)
#print zz
#plotFrame(xx,yy,YY[0,:],1)

#1D-Plot
#fignum1D = 2
#fig1D = p.figure(fignum1D)
#p.clf()
#data = fromfile(filename1D)
#N = data.shape[0]/2;
#t = data[0:N];
#E = data[N:2*N]
#axE = fig1D.add_subplot(2,1,1)
#axE.plot(t, E)
#p.xlabel('t in ms')
#p.ylabel('E in mV/cm')
#p.ion()
#p.show()
#p.figure(fignum1D)

#Het-Plot
#data = fromfile(filenameHet)
#if data[0] > 0.0:
#	fignumHet = 3
#	fig1D = p.figure(fignumHet)
#	p.clf()
#	Nhet=data[0]
#	Rhet=data[3:4*Nhet:4]
#	R=[Rhet[0]]
#	nR=[1]
#	j=0
#	for i in range(1,len(Rhet)):
#		if Rhet[i] == R[j]:
#			nR[j] = nR[j]+1
#		else:			
#			R.append(Rhet[i])
#			nR.append(1)
#			j = j+1
#
#	p.loglog(R, nR/(data[0]*(het_rmin/2)),'rd')
#	R=linspace(het_rmin,het_rmax,1000)
#	p.loglog(R, pow(R,het_alpha)/((het_alpha+1)*(pow(het_rmax,het_alpha+1)-pow(het_rmin,het_alpha+1))))
#	title = 'A=%smm^2 rho=%s rmin=%smm rmax=%smm' % (lx*ly,het_rho,het_rmin,het_rmax)
#	p.title(title)
#	p.xlabel('R in mm')
#	p.ylabel('p in mm^-1')
#	p.grid(True)
#	p.ion()
#	p.show()
#	p.figure(fignum1D)









