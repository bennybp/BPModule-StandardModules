Bondizer                                                             {#bondizer}
========

Especially for large covalently bonded molecules, fragments are usually made by
considering the bonding connectivity of the molecule.  Typically this involves
creating fragments such that the \f$i\f$-th fragment is atom \f$i\f$ and some
set of atoms that are up to \f$n\f$ bonds away.  There is quite a bit 
of ambiguity in this description, namely:
- Are only non-branching paths considered? (a branching path is a set of nodes
  and edges such that one of the nodes has 3 or more of its edges included, e.g.
  to grab all of the ammonia molecule in a single path requires allowing
  branched paths)
- Is each path eminating from atom \f$i\f$ an individual fragment? Or is the 
  union of all paths originating from the \f$i\f$-th atom the fragment?
- Are there emperical rules (e.g. \f$n\f$ bonds, but continue if in the middle
  of an aromatic system)

The \f$i\f$-th fragment resulting from the Bondizer class will be atom \f$i\f$
and all atoms that are within \f$n\f$ bonds of it.  Thus branching paths are
considered and the fragment is the union of all paths eminating from that atom.


## Options

The recognized options for Bondizer are:
- MAX_NBONDS: The maximum path length in C counting (MAX_NBONDS=1 will give you
  each pair of bonded atoms as a fragment, MAX_NBONDS=2 will give you all
  triples of bonded atoms, etc.)
  - By default is 2,147,483,647 (the maximum value of a 32-bit integer), which
    should be sufficiently large to ensure all bonded atoms are in a fragment.

