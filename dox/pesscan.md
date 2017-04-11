PES Analysis Tools                                                    {#pesscan}
==================

One of the most common tasks in quantum chemistry is scanning along some
coordinate of a potential energy surface (PES).  In the olden days, this was
usually done via elaborate Z-matrices, but now we want to automate much of this.
To that end, we introduce the PES Analysis tools.  The various functions are
described below.

Scans come in a variety of flavors.  Probably the most common type of scan is to
partition the system into two parts somehow (*e.g.* the two monomers of a dimer,
a surface and an adsorbed molecule, *etc.*) and to then stretch the system along
some vector between the two parts.  The vector can be defined by specifying two
points, one from each part of the system.  Arguably the two most common choices
for the location of these points are an atom or the center of mass.

Another common scan is to scan along an angle, torsion angle, or improper
torsion angle.  These can be specified respectively by three and four points.
Despite involving more than two points these scans also involve moving two parts
of the system away from each other (just no longer along a vector).

##Geometric Background

###Distance Scan

Given two rigid systems, one that contains the point \f$\vec{r}_1\f$ and
another that contains the point \f$\vec{r}_2\f$ our goal is to sweep along the
vector \f$\vec{r}_{12}\f$ that points from the first point to the second point.
We arbitraily choose to hold the system containing point 1 still and move the
other system.  Say we want to move system 2 a displacement \f$\Delta r\f$
that in general can be either positive (moving away from system 1) or negative
(moving towards system 1).  Either way, the new position of \f$\vec{r}_2\f$,
which we denote \f${\vec{r}_2}^\prime\f$, is now given by:
\f[
{\vec{r}_2}^\prime=\vec{r}_2+\Delta r\frac{\vec{r}_{12}}{||\vec{r}_{12}||}
\f]
All atoms in system 2 are shifted by the same vector.

###Angle Scan

Given three points located at \f$\vec{r}_1\f$, \f$\vec{r}_2\f$, and
\f$\vec{r}_3\f$,  such that \f$\vec{r}_1\f$ is associated with system 1,
\f$\vec{r}_{3}\f$ is associated with system 2, and \f$\vec{r}_2\f$ is the vertex
of the angle (association is arbitrary as it won't move) our goal is to sweep
over the angle \f$\theta\f$ formed between \f$\vec{r}_{21}\f$ and
\f$\vec{r}_{23}\f$.  Consequentially:
\f[
\vec{n}=\frac{\vec{r}_{23}\times\vec{r}_{21}}
{||\vec{r}_{23}||\ ||\vec{r}_{21}||}
\f]
is a unit norm perpendicular to the angle consistent with a right-handed
coordinate system in which \f$\vec{r}_{23}\f$ is the x-axis, \f$\vec{r}_{21}\f$
is the y-axis, and \f$\vec{n}\f$ is the z-axis.  If we hold \f$\vec{r}_{23}\f$
still a rotation matrix \f$\mathbf{R}\f$ about \f$\vec{n}\f$ by
\f$\Delta \theta\f$  (positive \f$\Delta\theta\f$ being counterclockwise
rotation of \f$\vec{r}_{21}\f$ away from \f$\vec{r}_{23}\f$ and negative
\f$\Delta\theta\f$ being clockwise rotation towards \f$\vec{r}_{23}\f$) is given
by:
\f[
\mathbf{R}=\begin{bmatrix}
\cos(\Delta\theta)+n_x^2[1-\cos(\Delta\theta)]&
n_xn_y[1-\cos(\Delta\theta)]-n_z\sin(\Delta\theta)&
n_xn_z[1-\cos(\Delta\theta)]+n_y\sin(\Delta\theta)\\
n_xn_y[1-\cos(\Delta\theta)]+n_z\sin(\Delta\theta)&
\cos(\Delta\theta)+n_y^2[1-\cos(\Delta\theta)]&
n_yn_z[1-\cos(\Delta\theta)]-n_x\sin(\Delta\theta) \\
n_xn_z[1-\cos(\Delta\theta)]-n_y\sin(\Delta\theta)&
n_yn_z[1-\cos(\Delta\theta)]+n_x\sin(\Delta\theta)&
\cos(\Delta\theta)+n_z^2[1-\cos(\Delta\theta)]
\end{bmatrix}
\f]
Multiplication of each point in the system associated with \f$\vec{r}_{1}\f$ by
\f$\mathbf{R}\f$ generates the appropriately rotated system.

###Torsion and Improper Torsion Scan

The proper torsion case is similar to the angle scenario.  Assume our torsion
is made of four points located at \f$\vec{r}_1\f$, \f$\vec{r}_2\f$,
\f$\vec{r}_3\f$, and \f$\vec{r}_4\f$, such that system 2 is associated with
\f$\vec{r}_{34}\f$ and system 1 is associated with \f$\vec{r}_{12}\f$ (assuming
disjoint systems this means \f$\vec{r}_{23}\f$ is split between the two
systems).  The goal of a torsion scan is to rotate about the \f$\vec{r}_{23}\f$
axis by an angle \f$\Delta\theta\f$.  A positive \f$\Delta\theta\f$ moves system
1 away from system 2 and a negative moves it towards it.

An improper torsion scan is a bit more complicated owing to the poor definition
of how to measure it.  Assume we have three orbital points located at
\f$\vec{r}_1\f$, \f$\vec{r}_3\f$, and \f$\vec{r}_4\f$ and the central point is
\f$\vec{r}_2\f$.  The value of the improper torsion is then taken as the torsion
angle between the plane formed from the first three points and the plane formed
from the last three points.  There exists 4! ways to define these planes given
by the permutations of the 4 points.  Half of these definitions correspond
to a definition of planarity for which the improper torsion angle is 180 degrees
and the other half are for a planarity definition of 0 degrees.  The first
definition occurs when the central point is either the second or third point and
the second definition arises when it is either the first or last point.  Within
either definition there are only six unique angles (the same set is obtained for
both unique placements of the central atom); furthermore, among these six angles
there are only three unique magnitudes (cyclic permutations have the same sign).

Pulsar adopts the philosophy that you can scan an improper torsion angle in the
same way you scan a torsion angle, giving four points.  The angle that will be
scanned is the torsion between \f$\vec{r}_{21}\f$ and \f$\vec{r}_{34}\f$ along
the \f$\vec{r}_{32}\f$ axis.  Note that this will only rotate the
\f$\vec{r}_{21}\f$ and not the other two legs of the pyramid.  In order to
scan such that all three angles change by the same amount (the scenario I think
most people would be interested in) requires me to work out some geometry that I
don't feel like doing at the moment...
