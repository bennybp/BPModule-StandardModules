CPGhoster                                                           {#cpghoster}
=========

Arguably the most common means of correcting for BSSE is to simply ensure each
fragment has the same basis set as the supersystem.  This is the esscence of the
Boys and Bernardi counterpoise scheme, or CP for short. This class takes a
set of fragments and returns the union of that fragment and its complement,
where the latter are converted to ghost atoms.

General notes:
- Setting the internal fragmenter to something like bondzier will result in a
  normal CP correction
- Setting the interanl fragmenter to an NMerizer will result in the many-body
  expansion being run in the super system basis set


##Options

- SYSTEM_FRAGMENTER_KEY : A key to a fragmenter that will provide the set of
  fragments to ghost
