NMerizer                                                             {#nmerizer}
========

In methods such as the many-body expansion we often start with a set of
fragments, called the monomers, and then want to compute the interaction between
each pair of fragments, each triple, etc.  In principle, we could get at these
interactions directly by using methods like SAPT, but in practice for three-body
and higher interactions the only sustainable way to compute these are via so
called supermolecular approaches.  In supermolecular approaches we compute the
interaction between two monomers by subtracting from the dimer (the system
formed from the union of the two monomers) the energy of each dimer.  For the
three-body energy, it is a linear combination of the energy of the trimer, the
dimers, and the monomers.  The point is we need to be able to form unions of
fragments and that's what this class does.

It is straightforward to show that the number of dimers you can form from 
\f$N\f$ fragments is given by the binomial coefficient \f$N\f$ choose 2, which
we denote \f$_NC_2\f$.  In general, you can form \f$_NC_m\f$, \f$m\f$-mers from
\f$N\f$ fragments.  This in turn means that the number of \f$m\f$-mers grows as
\f$\mathcal{O}(N^m)\f$, which gets big quickly.  Hence it is customary to try to
thin this number out somewhat.  The usual way is by distance thresholds.  It
follows from the theory of intermolecular interactions that the strength of an
interaction falls off as an inverse power of the distance between fragments with
higher powers for higher-body interactions.  Consequentially, we expect that
there is some critical distance at which an interaction is essentially zero and
can be ignored.  This class allows for truncating the number of \f$N\f$-mers in
this manner.

## Options
The recognized options for NMerizer are:
- SYSTEM_FRAGMENTER_KEY: The key for another system fragmenter that will
  generate the first set of fragments.
- TRUNCATION_KEY: The maximum number of fragments involved in a union
- DISTANCE_THRESHOLDS: A map from an integer, \f$m\f$, to a double, \f$r_0\f$, 
  such that if \f$r_0>(\prod_{i=1}^mr_i)^m}\f$ is true the \f$m\f$-mer is
  considered significant, where \f$r_i\f$ is the distance from the center of
  mass of the \f$m\f$-mer to the \f$i\f$-th monomer
