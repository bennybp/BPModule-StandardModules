Many-Body Expansion                                                       {#mbe}
===================

## Theory

Quantum mechanics is an inherently short-ranged theory.  Consequentially, a 
variety of methods have been suggested that are capable of exploiting spatial 
sparsity.  Within the many-body expansion (MBE) spatial sparsity is exploited by
 literally fragmenting one system of interest (called the supersystem) into a 
series of smaller systems.  Each of these smaller systems is then treated as if 
it was its own supersystem (i.e. it is run as an individual electronic structure 
computation).  Letting \f$E_I\f$ be the energy of subsystem \f$I\f$, then we can
approximate the energy of a supersystem that has been decomposed into \f$N\f$
subsystems as:
\f[
E\approx\sum_{I=1}^NE_I.
\f]
Assuming the smaller systems are much smaller than the supersystem this will be
a poor approximation to the supersystem energy.  To improve upon this we can
consider the so-called two-body interactions, the two-body interaction between
subsystems \f$I\f$ and \f$J\f$ being denoted, \f$\Delta E_{IJ}\f$, and given by:
\f[
\Delta E_{IJ}=E_{IJ}-E_{I}-E_{J}.
\f]
Summing over all possible unique two-body interactions we get:
\f[
E\approx\sum_{I=1}^NE_I+\sum_{I<J}^{_NC_2}\Delta E_{IJ}.
\f]
This is expected to be a much better approximation; however, for high accuracy
applications this is still not enough.  The next most important set of 
interactions are the three-body interactions, where the three-body interaction
among subsystems \f$I\f$, \f$J\f$, and \f$K\f$, denoted \f$\Delta E_{IJK}\f$ is
given by:
\f[
\Delta E_{IJK}=E_{IJK}-\Delta E_{IJ}-\Delta E_{IK}-\Delta E_{JK}-E_I-E_J-E_K.
\f]
This in turn leads to an approximate energy:
\f[
E\approx\sum_{I=1}^NE_I+\sum_{I<J}^{_NC_2}\Delta E_{IJ}+
         \sum_{I<J<K}^{_NC3}\Delta E_{IJK}
\f]
If we further include the four-body, the five-body, and up to the \f$N\f$-body
interactions we recover an exact expansion:
\f[
E=\sum_{I=1}^NE_I+\sum_{I<J}^{_NC_2}\Delta E_{IJ}+
         \sum_{I<J<K}^{_NC3}\Delta E_{IJK}\cdots\Delta E_{IJK\cdots N}.
\f]
In practice it is not the exactness of this expansion that is useful, but rather
the fact that it can be truncated at some \f$n\f$-body interact such that we 
ignore all higher contributions to the energy.  To the extent that 
\f$n\approx3\f$ this has the potential to afford huge savings.  Furthermore, all
of these individual computations are independent and can be run in an 
embarassingly paralell fashion.

With a truncation point of \f$n\f$ it can be shown that the approximate energy
can be written:
\f[
E\approx\sum_{m=1}^n(-1)^{n-m}{N-m-1\choose n-m}
        \sum_{I<J<\cdots<m}^{_NC_m}E_{IJ\cdots m}
\f]
The point is it ultimately looks like nothing more than a linear combination
of subsystem energies (if we pretend that each pair, triple, etc. of subsystems
is a single subsystem).

## Implementations

To some extent there is only one MBE; however, the literature is filled with a
large number of variants.  For the most part these variants only differ in how
the subsystems are defined.  Once we know the subsytems and their weights they
are all the same.  As coded here, the MBE deffers to an EnergyMethod and a
SystemFragmenter module to handle these variations.  Consequentially, it should
be possible to code any variant simply making a new SystemFragmenter.

Here are is a table of the SystemFragmenter, MBE variant pairs:


| MBE Variant                 | SystemFragmenter  |
|-----------------------------|-------------------|
| Plain MBE                   | NMerizer          |
| MBE in supersystem basis    | CPGhoster         |
| Valiron-Mayer Fucntional CP | VMFCGhoster       |

## Options

The Pulsar implementation of the MBE recognizes two options

- "METHOD_KEY" : The module key associated with the EnergyMethod to be used.
  The energy method will be called on each fragment
- "SYSTEM_FRAGMENTER_KEY" : The module key associated with the SystemFragmenter
  to be used to create the fragments

## Performance

The Pulsar implementation of the MBE is designed to be used with MPI
parallelizing the dispatching of the energy calls and threading parallelizing
the calls themselves.  
