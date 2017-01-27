MBE Helper Analysis Tool                                            {#mbehelper}
========================

There are two main uses of a many-body expansion (MBE): obtaining a
supersystem's energy (and energy related properties) at a reduced cost and to
study the physics of the underlying interactions.  The MBE Helper module is
designed to aid in the latter endeavor (the main MBE EnergyMethod can be used,
as is, for the former).  The remainder of this page describes the functions
contained in the MBE Helper suite.

MBE_wrapper
-----------

Assume you want the total one-body, the total two-body, the total three-body,
etc. energies of a system, this is the wrapper for you.  The MBE EnergyMethod
will only return one of these quantities at a time.  This wrapper will loop over
all orders of the MBE, up to the one requested by the user, store the results in
a dictionary and then return them for your analysis.  This will lead to many
duplicate calculations and thus you should only use this wrapper if your
underlying EnergyMethod runs quickly or uses the cache.  To use the MBE_wrapper:

~~~.py
mm=#The ModuleAdministrator instance
wfn=#The Wavefunction of the system you are studying
n=#The truncation order of the MBE to go up to
mbe_key=#The module key for the MBE EnergyMethod
nmerizer_key=#The key for the SystemFragmenter that makes the n-mers
egy_keys=#A list of EnergyMethods keys to call the MBE on
egys=MBE_wrapper(mm,wfn,n,mbe_key,nmerizer_key,egy_keys)
#egys will now be something like (assuming here n=2):
#{"egy_key_1":{1:-217.1234,2:-218.1234},
# "egy_key_2":{1:-219.1234,2:-220.1234}
#}
~~~
