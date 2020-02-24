import math
from fractions import Fraction
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import numpy as np
import matplotlib.pyplot as plt
#
nx    = int(3)                       # number of mesh cells in x
ny    = int(3)                       # number of mesh cells in y
lx    = 5.*nx                        # domain length in x
ly    = 5.*ny                        # domain length in y
dx    = lx/(1.0*nx)                  # mesh spacing in x
dy    = ly/(1.0*ny)                  # mesh spacing in y
xp     = np.arange(dx/2.,lx+dx/2.,dx)
yp     = np.arange(dy/2.,ly+dy/2.,dy)
xu     = np.arange(dx   ,lx+dx   ,dx)
yu     = np.arange(dy/2.,ly+dy/2.,dy)
xv     = np.arange(dx/2.,lx+dx/2.,dx)
yv     = np.arange(dy   ,ly+dy   ,dy)
XP,YP = np.meshgrid(xp,yp)
XP = np.transpose(XP)
YP = np.transpose(YP)
XU,YU = np.meshgrid(xu,yu)
XU = np.transpose(XU)
YU = np.transpose(YU)
XV,YV = np.meshgrid(xv,yv)
XV = np.transpose(XV)
YV = np.transpose(YV)
#
plt.plot([0.,lx], [0.,0.], '-k', lw=4)
plt.plot([0.,0.], [0.,ly], '-k', lw=4)
plt.plot([0.,lx], [0.,0.], '-k', lw=4)
plt.plot([0.,lx], [ly,ly], '-k', lw=4)
plt.plot([lx,lx], [0.,ly], '-k', lw=4)
for i in range(0,nx+1):
    plt.plot([(i)*dx, (i)*dx],[0., ly],'-k',color='0.2',lw=2.0)
for j in range(0,ny+1):
    plt.plot([0., lx],[(j)*dy, (j)*dy],'-k',color='0.2',lw=2.0)
for i in range(0,nx+2):
    plt.plot([(i-0.5)*dx, (i-0.5)*dx],[0., ly],'--r',lw=1.0)
for j in range(0,ny+1):
    plt.plot([-dx/2., lx+dx/2.],[(j-0.0)*dy, (j-0.0)*dy],'--r',lw=1.0)
for i in range(0,nx+1):
    plt.plot([(i-0.0)*dx, (i-0.0)*dx],[-dy/2., ly+dy/2.],'--g',lw=1.0)
for j in range(0,ny+2):
    plt.plot([0., lx],[(j-0.5)*dy, (j-0.5)*dy],'--g',lw=1.0)
eps = 0.15
for i in range(0,nx+2):
    xp = (i-0.5)*dx
    xu = xp + dx/2.
    xv = xp
    for j in range(0,ny+2):
        yp = (j-0.5)*dy
        yv = yp + dy/2.
        yu = yp
        plt.plot([xp],[yp], '.k')
        if( i < nx+1 ): plt.plot([xu],[yu], '>r')
        if( j < ny+1 ): plt.plot([xv],[yv], '^g')
        x0p = (xp-lx/2)/dx
        y0p = (yp-ly/2)/dy
        x0u = (xu-lx/2)/dx
        y0u = (yu-ly/2)/dy
        x0v = (xv-lx/2)/dx
        y0v = (yv-ly/2)/dy
        sxp = int(np.sign(x0p))
        syp = int(np.sign(y0p))
        sxu = int(np.sign(x0u))
        syu = int(np.sign(y0u))
        sxv = int(np.sign(x0v))
        syv = int(np.sign(y0v))
        if(    sxp == 1.  ): fxp = '+'
        elif(  sxp == -1. ): fxp = '-'
        else:                fxp = ''
        if(    syp == 1.  ): fyp = '+'
        elif(  syp == -1. ): fyp = '-'
        else:                fyp = ''
        if(    sxu == 1.  ): fxu = '+'
        elif(  sxu == -1. ): fxu = '-'
        else:                fxu = ''
        if(    syu == 1.  ): fyu = '+'
        elif(  syu == -1. ): fyu = '-'
        else:                fyu = ''
        if(    sxv == 1.  ): fxv = '+'
        elif(  sxv == -1. ): fxv = '-'
        else:                fxv = ''
        if(    syv == 1.  ): fyv = '+'
        elif(  syv == -1. ): fyv = '-'
        else:                fyv = ''
        if(sxp != 0.): fxp += str(Fraction(abs(x0p)))
        if(syp != 0.): fyp += str(Fraction(abs(y0p)))
        if(sxu != 0.): fxu += str(Fraction(abs(x0u)))
        if(syu != 0.): fyu += str(Fraction(abs(y0u)))
        if(sxv != 0.): fxv += str(Fraction(abs(x0v)))
        if(syv != 0.): fyv += str(Fraction(abs(y0v)))
        plt.text(xp+eps,yp+eps,'i'+fxp+','+'j'+fyp,fontsize=6,color='k')
        if( i < nx+1 ): plt.text(xu+eps,yu+eps,'i'+fxu+','+'j'+fyu,fontsize=6,color='r')
        if( j < ny+1 ): plt.text(xv+eps,yv+eps,'i'+fxv+','+'j'+fyv,fontsize=6,color='g')
epsx = dx/4.
plt.gca().set_xlim([-dx/2-epsx, lx+dx/2.+epsx])
epsy = dy/4.
plt.gca().set_ylim([-dy/2.-epsy, ly+dy/2.+epsy])
plt.gca().axis('off')
w, h = plt.figaspect(ly/lx)
plt.gcf().set_size_inches(w*1.3,h*1.3)
plt.savefig('staggered.pdf')
plt.show()
plt.clf()
plt.cla()
exit()
