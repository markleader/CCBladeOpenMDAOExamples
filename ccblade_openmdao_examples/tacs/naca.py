import matplotlib.pyplot as plt
import numpy as np

# TODO:
# - moments of inertia
# - region of negative space is not differentiable -> thickness is not differentiable analytically?
#   - instead of trimming the points that violate the centerline, use a differentiable max function and collapse them to the centerline so they have 0 area contribution?

def naca(x, m, p, tt):
    """
    Return coordinates for a NACA 4-digit airfoil (mptt) at normalized coordinate x
    """

    m /= 100.0
    p /= 10.0
    tt /= 100.0

    yt = 5.0*tt*(0.2969*x**0.5 - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)

    if (m > 0.0) and (p > 0.0):
        if x <= p:
            yc = (m/p**2)*(2*p*x - x**2)
            dycdx = 2*m*(p - x)/(p**2)
        else:
            yc = m*((1 - 2.0*p) + 2*p*x - x**2)/((1 - p)**2)
            dycdx = 2*m*(p - x)/((1 - p)**2)

        theta = np.arctan(dycdx)
        xl = x + yt*np.sin(theta)
        xu = x - yt*np.sin(theta)
        yl = yc - yt*np.cos(theta)
        yu = yc + yt*np.cos(theta)

    else:
        yc = 0.0
        xl = x
        xu = x
        yl = -yt
        yu = yt

    return xl, xu, yl, yu, yc

def compute_inner_points(m, p, tt, tw, npts=100):
    """
    Compute the inner points for a hollow airfoil with wall-thickness tw
    """

    x = np.linspace(0.0, 1.0, npts)
    xl = np.zeros(len(x))
    xu = np.zeros(len(x))
    yl = np.zeros(len(x))
    yu = np.zeros(len(x))
    yc = np.zeros(len(x))
    for i in range(len(x)):
        xl[i], xu[i], yl[i], yu[i], yc[i] = naca(x[i], m, p, tt)

    xli = np.zeros(len(x))
    xui = np.zeros(len(x))
    yli = np.zeros(len(x))
    yui = np.zeros(len(x))

    # Compute the slope of the normal for the leading edge
    theta = np.arctan2(xu[1]-xl[1], yu[1]-yl[1])
    xli[0] = xl[0] + tw*np.cos(theta)
    xui[0] = xu[0] + tw*np.cos(theta)
    yli[0] = yl[0] - tw*np.sin(theta)
    yui[0] = yu[0] + tw*np.sin(theta)

    # Compute the inner points by projecting from the normal on the upper and lower surfaces with a centered finite difference
    for i in range(1, len(x)-1):
        # Upper surface
        theta = np.arctan2(yu[i+1]-yu[i-1], xu[i+1]-xu[i-1])
        xui[i] = xu[i] + tw*np.sin(theta)
        yui[i] = yu[i] - tw*np.cos(theta)

        # Lower surface
        theta = np.arctan2(yl[i+1]-yl[i-1], xl[i+1]-xl[i-1])
        xli[i] = xl[i] - tw*np.sin(theta)
        yli[i] = yl[i] + tw*np.cos(theta)

    # Compute the trailing edge slope with a single backwards finite difference
    theta = np.arctan2(yu[-1]-yu[-2], xu[-1]-xu[-2])
    xui[-1] = xu[-1] + tw*np.sin(theta)
    yui[-1] = yu[-1] - tw*np.cos(theta)

    theta = np.arctan2(yl[-1]-yl[-2], xl[-1]-xl[-2])
    xli[-1] = xl[-1] - tw*np.sin(theta)
    yli[-1] = yl[-1] + tw*np.cos(theta)

    # Trim the interior points to be on the correct side of the centerline

    # First, interpolate yc onto the xui and xli grids
    yc_xui_interp = np.interp(xui, x, yc)
    yc_xli_interp = np.interp(xli, x, yc)

    xui = xui[(yui >= yc_xui_interp+1e-12)]
    yui = yui[(yui >= yc_xui_interp+1e-12)]
    xli = xli[(yli <= yc_xli_interp-1e-12)]
    yli = yli[(yli <= yc_xli_interp-1e-12)]

    return xli, xui, yli, yui

def quad_area(pl_i, pl_ip1, pu_i, pu_ip1):
    """
    Return the area of a quadrilateral defined by points using the shoelace formula
    """

    A = 0.0
    A += np.abs(0.5*((pl_i[0] - pl_ip1[0])*(pu_ip1[1] - pl_i[1]) - (pl_i[0] - pu_ip1[0])*(pl_ip1[1] - pl_i[1])))
    A += np.abs(0.5*((pl_i[0] - pu_i[0])*(pu_ip1[1] - pl_i[1]) - (pl_i[0] - pu_ip1[0])*(pu_i[1] - pl_i[1])))

    return A

def quad_centroid(pl_i, pl_ip1, pu_i, pu_ip1):
    """
    Return the centroid of a quadrilateral defined by points
    """

    A = quad_area(pl_i, pl_ip1, pu_i, pu_ip1)
    x_centroid = 0.0
    y_centroid = 0.0

    x_centroid += (pl_i[0] + pu_ip1[0] + pl_ip1[0])*np.abs(0.5*((pl_i[0] - pl_ip1[0])*(pu_ip1[1] - pl_i[1]) - (pl_i[0] - pu_ip1[0])*(pl_ip1[1] - pl_i[1])))/3.0
    x_centroid += (pl_i[0] + pu_ip1[0] + pu_i[0])*np.abs(0.5*((pl_i[0] - pl_ip1[0])*(pu_ip1[1] - pl_i[1]) - (pl_i[0] - pu_ip1[0])*(pl_ip1[1] - pl_i[1])))/3.0
    y_centroid += (pl_i[1] + pu_ip1[1] + pl_ip1[1])*np.abs(0.5*((pl_i[0] - pl_ip1[0])*(pu_ip1[1] - pl_i[1]) - (pl_i[0] - pu_ip1[0])*(pl_ip1[1] - pl_i[1])))/3.0
    y_centroid += (pl_i[1] + pu_ip1[1] + pu_i[1])*np.abs(0.5*((pl_i[0] - pl_ip1[0])*(pu_ip1[1] - pl_i[1]) - (pl_i[0] - pu_ip1[0])*(pl_ip1[1] - pl_i[1])))/3.0

    x_centroid /= A
    y_centroid /= A

    return x_centroid, y_centroid

def area(m, p, tt, tw, npts=100):
    """
    Compute the area of a hollow NACA airfoil with wall-thickness tw
    """

    x = np.linspace(0.0, 1.0, npts)

    # Compute the area of the solid cross section
    A = 0.0
    for i in range(npts-1):
        # Get the points
        xl_i,   xu_i,   yl_i,   yu_i,   yc_i   = naca(x[i],   m, p, tt)
        xl_ip1, xu_ip1, yl_ip1, yu_ip1, yc_ip1 = naca(x[i+1], m, p, tt)

        # Add the area of this quadrilateral
        dA = quad_area([xl_i, yl_i], [xl_ip1, yl_ip1], [xu_i, yu_i], [xu_ip1, yu_ip1])
        A += dA

    # Subtract the area of the internal negative space
    xli, xui, yli, yui = compute_inner_points(m, p, tt, tw, npts=npts)
    for i in range(len(xui)-1):
        dx = xui[i+1] - xui[i]
        _, _, _, _, yc_i = naca(xui[i], m, p, tt)
        _, _, _, _, yc_ip1 = naca(xui[i+1], m, p, tt)
        h1 = yui[i] - yc_i
        h2 = yui[i+1] - yc_ip1
        A -= 0.5*dx*(h1 + h2)

    for i in range(len(xli)-1):
        dx = xli[i+1] - xli[i]
        _, _, _, _, yc_i = naca(xli[i], m, p, tt)
        _, _, _, _, yc_ip1 = naca(xli[i+1], m, p, tt)
        h1 = yc_i - yli[i]
        h2 = yc_ip1 - yli[i+1]
        A -= 0.5*dx*(h1 + h2)

    return A

def centroid(m, p, tt, tw, npts=100):

    """
    Compute the centroid of a hollow NACA airfoil with wall-thickness tw
    """

    x = np.linspace(0.0, 1.0, npts)

    # Compute the centroid of the solid cross section
    A_solid = 0.0
    x_centroid_solid = 0.0
    y_centroid_solid = 0.0

    for i in range(npts-1):
        # Get the points
        xl_i, xu_i, yl_i, yu_i, yc_i = naca(x[i], m, p, tt)
        xl_ip1, xu_ip1, yl_ip1, yu_ip1, yc_ip1 = naca(x[i+1], m, p, tt)

        # Compute the area of this quadrilateral
        dA = quad_area([xl_i, yl_i], [xl_ip1, yl_ip1], [xu_i, yu_i], [xu_ip1, yu_ip1])
        dxc, dyc = quad_centroid([xl_i, yl_i], [xl_ip1, yl_ip1], [xu_i, yu_i], [xu_ip1, yu_ip1])
        A_solid += dA
        x_centroid_solid += dxc*dA
        y_centroid_solid += dyc*dA

    x_centroid_solid /= A_solid
    y_centroid_solid /= A_solid

    # Compute the centroid of the negative cross section
    A_negative = 0.0
    x_centroid_negative = 0.0
    y_centroid_negative = 0.0

    xli, xui, yli, yui = compute_inner_points(m, p, tt, tw, npts=npts)
    for i in range(len(xui)-1):
        _, _, _, _, yc_i = naca(xui[i], m, p, tt)
        _, _, _, _, yc_ip1 = naca(xui[i+1], m, p, tt)
        dA = quad_area([xui[i], yc_i], [xui[i+1], yc_ip1], [xui[i], yui[i]], [xui[i+1], yui[i+1]])
        dxc, dyc = quad_centroid([xui[i], yc_i], [xui[i+1], yc_ip1], [xui[i], yui[i]], [xui[i+1], yui[i+1]])
        A_negative += dA
        x_centroid_negative += dxc*dA
        y_centroid_negative += dyc*dA

    for i in range(len(xli)-1):
        _, _, _, _, yc_i = naca(xli[i], m, p, tt)
        _, _, _, _, yc_ip1 = naca(xli[i+1], m, p, tt)
        dA = quad_area([xli[i], yc_i], [xli[i+1], yc_ip1], [xli[i], yli[i]], [xli[i+1], yli[i+1]])
        dxc, dyc = quad_centroid([xli[i], yc_i], [xli[i+1], yc_ip1], [xli[i], yli[i]], [xli[i+1], yli[i+1]])
        A_negative += dA
        x_centroid_negative += dxc*dA
        y_centroid_negative += dyc*dA

    if A_negative > 0.0:
        x_centroid_negative /= A_negative
        y_centroid_negative /= A_negative

    x_centroid = (A_solid*x_centroid_solid - A_negative*x_centroid_negative)/(A_solid - A_negative)
    y_centroid = (A_solid*y_centroid_solid - A_negative*y_centroid_negative)/(A_solid - A_negative)

    return x_centroid, y_centroid

def moments_of_intertia(m, p, tt, tw, npts=100):
    """
    Compute the moments of inertia of a hollow airfoil with wall-thickness tw
    """

    # First, compute the moments of intertia of the solid section about the origin
    x = np.linspace(0.0, 1.0, npts)
    Iyy_solid = 0.0
    Izz_solid = 0.0
    Iyz_solid = 0.0

    for i in range(npts-1):
        # Get the points
        xl_i, xu_i, yl_i, yu_i, yc_i = naca(x[i], m, p, tt)
        xl_ip1, xu_ip1, yl_ip1, yu_ip1, yc_ip1 = naca(x[i+1], m, p, tt)

        # Compute the area of this quadrilateral
        dA = quad_area([xl_i, yl_i], [xl_ip1, yl_ip1], [xu_i, yu_i], [xu_ip1, yu_ip1])
        dxc, dyc = quad_centroid([xl_i, yl_i], [xl_ip1, yl_ip1], [xu_i, yu_i], [xu_ip1, yu_ip1])

        Iyy_solid += (1./12.)*

    # Next, compute the moments of inertis of the negative section about the origin


    # Add them together and shift the moments of inertia to the centroid


    return Ixx, Iyy, Iyz

m = 2
p = 4
tt = 12
tw = 0.05

x = np.linspace(0.0, 1.0, 1000)
xl = np.zeros(len(x))
xu = np.zeros(len(x))
yl = np.zeros(len(x))
yu = np.zeros(len(x))
for i in range(len(x)):
    xl[i], xu[i], yl[i], yu[i], yc = naca(x[i], m, p, tt)

xli, xui, yli, yui = compute_inner_points(m, p, tt, tw, npts=1000)

fig = plt.figure()
ax = fig.add_subplot(111)

plt.scatter(xl, yl, color="tab:blue")
plt.scatter(xu, yu, color="tab:blue")
plt.scatter(xli, yli, color="tab:orange")
plt.scatter(xui, yui, color="tab:green")
plt.xlim(-0.05,1.05)
plt.ylim(-0.15,0.15)

ax.set_aspect('equal', adjustable='box')

# plt.xlabel("x")
# plt.ylabel("sinx")

plt.show()

A = area(m, p, tt, tw, npts=200)
print(A)

# npts = [50, 100, 150, 200, 250, 300, 350, 500, 750, 1000]
# A = np.zeros(len(npts))
# for i in range(len(npts)):
#     A[i] = area(m, p, tt, tw, npts=npts[i])

# fig = plt.figure()
# ax = fig.add_subplot(111)

# plt.plot(npts, 100.0*(A-A[-1])/A[-1])
# plt.show()
print(centroid(m, p, tt, tw))
