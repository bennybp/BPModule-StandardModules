import copy
import pulsar as psr
def mbe_helper(order,wfn,mbe_key,frag_key,n,mm,parent=0,HasGhosts=False,**kwargs):
    """
    A function designed to quickly compute useful MBE quantities

        Applying BSSE corrections to the MBE is essential.  This function allows
        for that using the flag HasGhosts.  If this is true then we will be
        creating a set of fragments that conatins ghost atoms.  In all cases,
        frag_key is the key of the fragmenter that generates the non-ghosted
        fragments.  If HasGhost is true then you must specify the keyword
        argument MAIN_FRAGGER, which is the key to the fragmenter that will make
        the final ghosted fragments.  Several of these flavors are convenience
        wrappers around ghosters and so we optionally allow you to specify the
        ghoster via the GHOST_FRAGGER keyword. 

        order (int): The derivative order of the MBE to compute
        wfn (psr.datastore.wavefunction): The wavefunction to give to the MBE 
        mbe_key (string): The key for the MBE to run
        frag_key (string): The key for the fragmenter
        n (int): The highest truncation order to generate statistics for
        mm (psr.modulemanager.ModuleManager): The module manager to use
        parent (int): The ID of the parent module defaults to 0
        HasGhosts (bool) : Indicates that this is actually for some form of BSSE
                           corrected MBE.  See above for more details
 
        returns a dictionary of names of energetic quantities indexed by their
                names
    """
    
    #Run the MBE
    my_mod=mm.get_module(mbe_key,parent)
    NewWfn,newegy=my_mod.deriv(order,wfn)#Egy reomputed as byproduct below
    print(newegy)
    
    #Get the name of the method the MBE called and make instance
    egy_meth_name=my_mod.options().get("METHOD")
    print(egy_meth_name)
    egy_meth=mm.get_module(egy_meth_name,parent)
    
    Returns={}#Will be the results
    
    prefix="" if "PREFIX" not in kwargs else kwargs["PREFIX"] #For printing
    
    def egy_title(n_in):
        return prefix+" "+str(n_in)+"-body energy"
    
    ckco=psr.copy_key_change_options
    
    for i in range(n,0,-1):
        #This is the fragmenter that generates the non-ghosted fragments
        fragger=ckco(frag_key,parent,mm,{"TRUNCATION_ORDER":i})
        
        #Pass it to the fragger that makes ghosts
        if HasGhosts:
            fragger=ckco(kwargs["MAIN_FRAGGER"],parent,mm,{"SYSTEM_FRAGMENTER_KEY":fragger.key()})
        mmers=fragger.fragmentize(wfn.system)
        egyname=egy_title(i)
        Returns[egyname]=0.0
        print("------------------------")
        print(wfn.system)
        for key,val in mmers.items():
            oldsys=wfn.system
            wfn.system=val.NMer
            tempwfn,tempEgy=egy_meth.deriv(order,wfn)
            Returns[str(key)]=tempEgy
            Returns[egyname]+=(tempEgy[0]*val.Weight)
            wfn.system=oldsys
    print(Returns)
    exit()

    return NewWfn,Returns

def vmfc_helper(order,wfn,mbe_key,frag_key,sub_frag_key,n,mm,
                parent=0,ghoster_key="PSR_GHOST_FRAG"):
    """A function designed to quickly compute useful VMFC quantities
    
       A normal MBE is a subset of VMFC and is included in the results

        order (int): The derivative order of the MBE to compute
        wfn (psr.datastore.wavefunction): The wavefunction to give to the MBE 
        mbe_key (string): The key for the MBE to run
        frag_key (string): The key for the VMFC fragger
        sub_frag_key (string) : The key for the fragger that makes the initial
                                set of fragments for VMFC
        n (int): The highest truncation order to generate statistics for
        mm (psr.modulemanager.ModuleManager): The module manager to use
        parent (int): The ID of the parent module defaults to 0
        ghoster_key (string): The key used to make ghosted atoms
        
        returns a dictionary of names of energetic quantities indexed by their
                names           
        
    """
    #First run the VMFC computations
    psr.print_global_debug("VMFC(m) energies and interactions")
    BestWfn,Returns=mbe_helper(order,wfn,mbe_key,sub_frag_key,n,mm,parent,True,
        PREFIX="BSSE-Corrected ",MAIN_FRAGGER=frag_key,GHOST_FRAGGER=ghoster_key)
    
    #Now the normal MBE computations, we have to be careful or we end up running
    #the computations again and not using the cached results.  This is because
    #the universes in the supersystem won't be the same.  To get around this
    #we make a null ghoster
    ghoster_name=mm.generate_unique_key()
    mm.duplicate_key(ghoster_key,ghoster_name)
    mm.change_option(ghoster_name,"GHOST_TRUNCATION_ORDERS",{i:0 for i in range(0,n)})
    psr.print_global_debug("In VMFC(m), but MBE energies and interactions")
    okwfn,okReturns=mbe_helper(order,wfn,mbe_key,sub_frag_key,n,mm,parent,True,
        MAIN_FRAGGER=ghoster_name)
    Returns.update(okReturns)
    return BestWfn,Returns
