import pulsar as psr
import math

def rot_mat(norm,theta):
    cost,sint,cos1=math.cos(theta),math.sin(theta),1-math.cos(theta)
    nx,ny,nz=norm[0],norm[1],norm[2]
    nx2,ny2,nz2=nx**2,ny**2,nz**2
    nxny,nxnz,nynz=nx*ny,nx*nz,ny*nz
    return [cost+nx2*cos1,     nxny*cos1-nz*sint, nxnz*cos1+ny*sint,
            nxny*cos1+nz*sint, cost+ny2*cos1,     nynz*cos1-nx*sint,
            nxnz*cos1-ny*sint, nynz*cos1+nx*sint, cost+nz2*cos1]

def vec_diff(v1,v2):
    return [v2[i]-v1[i] for i in range(3)]

def vec_dot(v1,v2):
    return sum([v1[i]*v2[i] for i in range(3)])

def vec_mag(v1):
    return math.sqrt(vec_dot(v1,v1))

def vec_scale(v1,x):
    return [v1[i]*x for i in range(3)]

def vec_unit(v1):
    return vec_scale(v1,1.0/vec_mag(v1))

def vec_cross(v1,v2):
    return [v1[1]*v2[2]-v1[2]*v2[1],
            v1[2]*v2[0]-v1[0]*v2[2],
            v1[0]*v2[1]-v1[1]*v2[0]]

def angle_scan(un,t,sys_in,orig):
    to_origin=[orig[i]*-1 for i in range(3)]
    from_origin=[orig[i] for i in range(3)]
    R=rot_mat(un,t)
    trans_sys=psr.translate(sys_in[0],to_origin)
    rot_sys=psr.rotate(trans_sys,R)
    origin_sys=psr.translate(rot_sys,from_origin)
    return psr.System(origin_sys.as_universe()+sys_in[1].as_universe(),True)

def pes_scan(init_sys,points,scan_range):
    """Performs a potential energy scan (PES) over a range of values.

    Given a system split in two, a set of points defining the scan coordinate,
    and a range of coordinate values, this generator function will yield the
    system at the requested points.  This function recognizes a few types for
    values of points: the string "COM" to represent the center of mass, a
    pulsar::Atom instance, or a list of three floats

    Args:
    init_sys : pair of pulsar::System instance at the initial geometry
    points : tuple of points defining the coordinates
    range : list of separations (in a.u.) or angles (in radians) to scan over

    Yields: The system you should compute the energy of
    """

    actual_points=[]
    for pt in points:
        if pt=="COM":
            com=psr.center_of_mass(init_sys[len(actual_points)])
            actual_points.append(com)
        else:
            actual_points.append(pt)

    if len(actual_points)==2:
        uv=vec_unit(vec_diff(actual_points[0],actual_points[1]))
        for x in scan_range:
            new_sys=psr.translate(init_sys[1],vec_scale(uv,x))
            yield psr.System(init_sys[0].as_universe()+new_sys.as_universe(),True)
    if len(actual_points)==3:
        r21=vec_unit(vec_diff(actual_points[1],actual_points[0]))
        r23=vec_unit(vec_diff(actual_points[1],actual_points[2]))
        un=vec_unit(vec_cross(r23,r21))
        for t in scan_range:
            yield angle_scan(un,t,init_sys,actual_points[1])
    if len(actual_points)==4: #and len(init_sys)==2:
        r32=vec_unit(vec_diff(actual_points[2],actual_points[1]))
        for t in scan_range:
            yield angle_scan(r32,t,init_sys,actual_points[2])
    elif len(actual_points)==4 and len(init_sys)==4:
        for t in scan_range:
            #Central system doesn't move
            new_sys=psr.System(init_sys[3].as_universe(),True)
            for shift in range(3):
                #Need to put points in backwards to maintain order
                #Cyclic permute points that move
                # actual_point idx -> temp_point_idx
                # 1st cycle) 2->0, 1->1, 0->2
                # 2nd cycle) 1->0, 0->1, 2->2
                # 3rd cycle) 0->0, 2->1, 1->2
                temp_points=[actual_points[2-(i+shift)%3] for i in range(3)]
                #Angle is temp_points[0]-actual_points[3]temp_points[1]-temp_points[2]
                r32=vec_unit(vec_diff(temp_points[1],actual_points[3]))
                new_sys=angle_scan(r32,t/2.0,[init_sys[2-shift],new_sys],actual_points[3])
            yield new_sys
