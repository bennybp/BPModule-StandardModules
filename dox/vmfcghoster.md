VMFCGhoster                                                       {#vmfcghoster}
===========

Arguably the right way to perform a basis set superposition correction is with
the method known as Valiron-Mayer functional counterpoise correction.  If we
take \f$E_X(Y)\f$ to be the energy of some system \f$X\f$, that has been 
computed in the \f$Y\f$ basis, VMFC computes the energy of an \f$N\f$ fragment
system as:
\f[
E=\sum_{I=1}^N E_I(I)+\sum_{I<J=1}^{_NC_2}\left[E_{IJ}(IJ)-E_{I}(I)-E_J(IJ)+\cdots
\f]
that is all the terms needed to compute the one-body energy are in the one-body
basis, all the terms needed to compute the two-body energy are in the dimer
basis, etc.

This class is designed to be used after an NMerizer.  Specifically it calls an
NMerizer and then adds to the resulting systems the necessary ghost plus real
system computations necessary to compute VMFC.  For example if the input
fragmenter provides a set of monomers, dimers, and trimers, you will get back
your original set, plus the monomers in the dimer basis, the monomers in the
trimer basis, and the dimers in the trimer basis.  Of course the coefficients
will be setup correctly for you as well.  If you don't use this class after an
NMerizer it will not do anything.

##Options

- SYSTEM_FRAGMENTER_KEY : A key to a fragmenter that will provide the set of
  fragments to ghost
