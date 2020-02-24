import math
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import numpy as np
import matplotlib.pyplot as plt
#
nx    = int(5)                       # number of mesh cells in x
ny    = int(5)                       # number of mesh cells in y
lx    = 1.*nx                        # domain length in x
ly    = 1.*ny                        # domain length in y
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
    plt.plot([0., lx],[(j)*dy, (j)*dy],'-k',color='0.2',lw=1.0)
for i in range(0,nx+2):
    plt.plot([(i-0.5)*dx, (i-0.5)*dx],[0., ly],'--r',lw=1.0)
for j in range(0,ny+1):
    plt.plot([-dx/2., lx+dx/2.],[(j-0.0)*dy, (j-0.0)*dy],'--r',lw=1.0)
for i in range(0,nx+1):
    plt.plot([(i-0.0)*dx, (i-0.0)*dx],[-dy/2., ly+dy/2.],'--g',lw=1.0)
for j in range(0,ny+2):
    plt.plot([0., lx],[(j-0.5)*dy, (j-0.5)*dy],'--g',lw=1.0)
for i in range(0,nx+2):
    xp = (i-0.5)*dx
    xu = xp + dx/2.
    xv = xp
    for j in range(0,ny+2):
        yp = (j-0.5)*dy
        yv = yp + dy/2.
        yu = yp
        plt.plot([xp],[yp], '.k')
        plt.plot([xu],[yu], '>r')
        plt.plot([xv],[yv], '^g')
#plt.gca().set_xticks(np.arange(-dx/2.,lx+dx/2.,dx))
#plt.gca().set_yticks(np.arange(-dy/2.,ly+dx/2.,dy))
epsx = dx/4.
plt.gca().set_xlim([-dx/2.-epsx, lx+dx/2.+epsx])
epsy = dy/4.
plt.gca().set_ylim([-dy/2.-epsy, ly+dy/2.+epsy])
plt.gca().set_yticklabels([])#np.arange(1,nx+1))
plt.gca().set_xticklabels([])#np.arange(1,ny+1))
plt.gca().xaxis.set_ticks_position('none')
plt.gca().yaxis.set_ticks_position('none')
#plt.gca().set_xlabel(r'$\mathbf{i}$')
#plt.gca().set_ylabel(r'$\mathbf{j}$')
#plt.legend(loc=1)
w, h = plt.figaspect(ly/lx)
plt.gcf().set_size_inches(w,h)
plt.show()
plt.clf()
plt.cla()
exit()
