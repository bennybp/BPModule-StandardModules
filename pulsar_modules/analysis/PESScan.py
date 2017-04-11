import pulsar as psr
import math

def rot_mat_(norm,theta):
    cost,sint,cos1=math.cos(theta),math.sin(theta),1-math.cos(theta)
    nx,ny,nz=norm[0],norm[1],norm[2]
    nx2,ny2,nz2=nx**2,ny**2,nz**2
    nxny,nxnz,nynz=nx*ny,nx*nz,ny*nz
    return [cost+nx2*cos1,     nxny*cos1-nz*sint, nxnz*cos1+ny*sint,
            nxny*cos1+nz*sint, cost+ny2*cos1,     nynz*cos1-nx*sint,
            nxnz*cos1-ny*sint, nynz*cos1+nx*sint, cost+nz2*cos1]

def vec_diff_(v1,v2):
    return [v2[i]-v1[i] for i in range(3)]

def vec_sum_(v1,v2):
    return [v1[i]+v2[i] for i in range(3)]

def vec_dot_(v1,v2):
    return sum([v1[i]*v2[i] for i in range(3)])

def vec_mag_(v1):
    return math.sqrt(vec_dot_(v1,v1))

def vec_scale_(v1,x):
    return [v1[i]*x for i in range(3)]

def vec_unit_(v1):
    return vec_scale_(v1,1.0/vec_mag_(v1))

def vec_cross_(v1,v2):
    return [v1[1]*v2[2]-v1[2]*v2[1],
            v1[2]*v2[0]-v1[0]*v2[2],
            v1[0]*v2[1]-v1[1]*v2[0]]

def mid_point(pt1,pt2):
    """Wrapper function for clarity, returns the point between the two points"""
    return vec_scale_(vec_sum_(pt1,pt2),0.5)


def angle_scan(un,t,sys_in,orig):
    to_origin=[orig[i]*-1 for i in range(3)]
    from_origin=[orig[i] for i in range(3)]
    R=rot_mat(un,t)
    trans_sys=psr.translate(sys_in[0],to_origin)
    rot_sys=psr.rotate(trans_sys,R)
    origin_sys=psr.translate(rot_sys,from_origin)
    return psr.System(origin_sys.as_universe()+sys_in[1].as_universe(),True)

def make_pes_range(start,stop,increment,to_au=1.0):
    """Computes a scan range for you.

    The Python range function won't work with integers, which is not what we
    usually want to scan over for a PES scan.  This little function will map
    those integers to a range for you.

    Note: range is inclusive, i.e. [start,stop], not [start,stop)

    Args:
    start : float the first displacement to consider
    stop  : float the last displacement to consider
    increment : float the step size between displacements
    to_au : the conversion from your input units to atomic units (or radians)

    Return: a list of displacements
    """
    nsteps=int((stop-start)/increment)
    return [(start+i*increment)*to_au for i in range(0,nsteps+1)]


def pes_scan(init_sys,points,scan_range):
    """Performs a potential energy scan (PES) over a range of values.

    Given a system split in two, a set of points defining the scan coordinate,
    and a range of coordinate values, this generator function will yield the
    system at the requested points.  This function always moves system 2.


    For a bond scan the first point is assumed to be part of system 1 and the
    second part of system 2.  For an angle scan the first point is assumed to be
    part of system 1, the second is the vertex (can be in either system), and the
    third is assumed associated with system2.

    Args:
    init_sys : pair of pulsar::System instance at the initial geometry the
               second will be moved
    points : tuple of points defining the coordinates, it is assumed that the
             first point is associated
    range : list of separations (in a.u.) or angles (in radians) to scan over

    Yields: The system you should compute the energy of
    """

    #Helps fit yield calls on one line and binds arguments
    def update_subsys(z):
        return psr.update_subsystem(init_sys[0]+init_sys[1],init_sys[1],z)

    if len(points)==2:
        uv=vec_unit_(vec_diff_(points[0],points[1]))
        for x in scan_range:
            new_sys=psr.translate(init_sys[1],vec_scale_(uv,x))
            yield update_subsys(new_sys)
    elif len(points)==3:
        r21=vec_unit_(vec_diff_(points[1],points[0]))
        r23=vec_unit_(vec_diff_(points[1],points[2]))
        un=vec_unit_(vec_cross_(r23,r21))
        for t in scan_range:
            yield angle_scan(un,t,init_sys,points[1])
    elif len(points)==4: #and len(init_sys)==2:
        r32=vec_unit_(vec_diff_(points[2],points[1]))
        for t in scan_range:
            yield angle_scan(r32,t,init_sys,points[2])
