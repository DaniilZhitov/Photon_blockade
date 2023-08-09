import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from matplotlib import animation
def rho_from_vector(r=1, theta=0, phi=0):
    x, y, z = r * np.sin(theta) * np.cos(phi), r * np.sin(theta) * np.sin(phi), r * np.cos(theta)
    rho = 0.5 * (qeye(2) + x * sigmax() + y * sigmay() + z * sigmaz())
    return rho

x_min, x_max, y_min, y_max = -2, 2, -2, 2
xvec = np.linspace(x_min, x_max, 60)
yvec = np.linspace(y_min, y_max, 60)
if False:
    fig, axes = plt.subplots(1,2)
    r, theta, phi = 0.6, 0, 0.0
    theta = np.arccos(1/r-np.sqrt((1/r)**2-1))+0.1
    cx, cy, cz = r * np.sin(theta) * np.cos(phi), r * np.sin(theta) * np.sin(phi), r * np.cos(theta)
    rho = rho_from_vector(r, theta, phi)
    d_crit = -np.sqrt(2)*r*np.sin(theta)/(2*(1-cz))
    print("d_crit="+str(d_crit))

    contours0 = axes[0].contourf(xvec, yvec, np.real(wigner(rho_from_vector(0.5, 2.3, 0.0), xvec, yvec, method='laguerre')))
    x_mesh, y_mesh = np.meshgrid(xvec, yvec)
    d2 = x_mesh**2+y_mesh**2
    f00 = 1/np.pi * np.exp(-d2)
    f11 = -1/np.pi * (1-2*d2)*np.exp(-d2)
    f10 = 1/np.pi * np.exp(1.0j*np.arctan2(y_mesh, x_mesh)) * np.sqrt(2*d2)*np.exp(-d2)
    f01 = np.conj(f10)

    f_rho = rho.data[0,0]*f00+rho.data[0,1]*f10+rho.data[1,0]*f01+rho.data[1,1]*f11
    f_rho = f_rho*np.exp(d2)
    W_min = ((1 - cz) * d_crit ** 2 + np.sqrt(2) * r * np.sin(theta) * d_crit + cz)/np.pi
    print('W_min_theor=' + str(W_min))
    print('W_min_plot='+str(np.min(f_rho)))
    contours1 = axes[1].contourf(xvec, yvec, np.real(f_rho), levels=np.linspace(-0.2,0.2,20))

    cbar = fig.colorbar(contours0)
    cbar = fig.colorbar(contours1)
    plt.show()
    quit()

if True:
    fig = plt.figure(figsize=(4,6))
    plt.xlabel('$c_x$')
    plt.ylabel('$c_z$')
  #  plt.x_label('$c_x$')
 #   ax.set_ylabel('$c_z$')
    for r in np.linspace(0, 1, 30):
        print("current r=", str(r))
        for theta in np.linspace(0, np.pi, 41):
            phi = 0
            wig = wigner(rho_from_vector(r=r, theta=theta, phi=phi), xvec, yvec, method='laguerre')
            x, y, z = r * np.sin(theta) * np.cos(phi), r * np.sin(theta) * np.sin(phi), r * np.cos(theta)
            flag = np.any(np.array(wig) < -.0001)
            if flag:
                plt.scatter(x, z, color='blue', s=1)
            else:
                plt.scatter(x, z, color='orange', s=1)
    plt.tight_layout()
    plt.savefig('Bloch_sphere_cross_section_negativity.png')
    quit()


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
for r in np.linspace(0, 1, 4):
    print("current r=", str(r))
    for theta in np.linspace(0, np.pi, 21):
        for phi in np.linspace(0, 2*np.pi, 5):
            wig = wigner(rho_from_vector(r=r, theta=theta, phi=phi), xvec, yvec, method='laguerre')
            x, y, z = r * np.sin(theta) * np.cos(phi), r * np.sin(theta) * np.sin(phi), r * np.cos(theta)
            flag = np.any(np.array(wig) < -.0001)
            if flag:
                ax.scatter(x, y, z, color='blue', s=1)
            else:
                ax.scatter(x, y, z, color='orange', s=1)
plt.show()
quit()

fig, ax = plt.subplots()

#contours = ax.contourf(xvec, yvec, , xvec, yvec, method='laguerre'))
#cbar = fig.colorbar(contours)
ax.set_title('r=1')
plt.xlabel('Q')
plt.ylabel('P')
def update(frame):
    print(frame)
    r_curr = 1 - frame/100
    return ax.contourf(xvec, yvec, wigner(rho_from_vector(r=r_curr), xvec, yvec))

#ani = animation.FuncAnimation(fig=fig, func=update, frames=100, interval=60)
#ani.save(filename='wigner_animation_r.gif', writer='pillow', fps=15)
#plt.show()