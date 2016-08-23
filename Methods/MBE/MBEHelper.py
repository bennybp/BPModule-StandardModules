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
    Returns={}
    prefix="" if "PREFIX" not in kwargs else kwargs["PREFIX"]
    def TotalEName(i):
        return "Total "+prefix+str(i)+"-body derivative"
    def IntEName(i):
        return "Total of "+str(i)+"-body "+prefix+"interactions"
    
    for i in range(n,0,-1):
        temp_mbe_name=mm.generate_unique_key()
        temp_frag_name=mm.generate_unique_key()
        mm.duplicate_key(mbe_key,temp_mbe_name)
        mm.duplicate_key(frag_key,temp_frag_name)
        mm.change_option(temp_frag_name,"TRUNCATION_ORDER",i)
        frag_name=temp_frag_name
        if HasGhosts:
            bsse_name=mm.generate_unique_key()
            mm.duplicate_key(kwargs["MAIN_FRAGGER"],bsse_name)
            frag_name=bsse_name #Will actually call this guy
            if "GHOST_FRAGGER" in kwargs:
                ghoster_name=mm.generate_unique_key()
                mm.duplicate_key(kwargs["GHOST_FRAGGER"],ghoster_name)
                mm.change_option(bsse_name,"GHOSTER_KEY",ghoster_name)
                bsse_name=ghoster_name #Ghoster actually gets the fragmenter
            mm.change_option(bsse_name,"SYSTEM_FRAGMENTER_KEY",temp_frag_name)
        mm.change_option(temp_mbe_name,"FRAGMENTIZER",frag_name)
        mymod=mm.get_module(temp_mbe_name,parent)
        oldwfn,Result=mymod.deriv(order,wfn)
        if i==n:BestWfn=oldwfn
        Returns[TotalEName(i)]=Result
        if i<n:
            Returns[IntEName(i+1)]=[Returns[TotalEName(i+1)][j]-Result[j] 
                for j in range(0,len(Result))]
    return BestWfn,Returns

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
    BestWfn,Returns=mbe_helper(order,wfn,mbe_key,sub_frag_key,n,mm,parent,True,
        PREFIX="BSSE-Corrected ",MAIN_FRAGGER=frag_key,GHOST_FRAGGER=ghoster_key)
    
    #Now the normal MBE computations, we have to be careful or we end up running
    #the computations again and not using the cached results.  This is because
    #the universes in the supersystem won't be the same.  To get around this
    #we make a null ghoster
    ghoster_name=mm.generate_unique_key()
    mm.duplicate_key(ghoster_key,ghoster_name)
    mm.change_option(ghoster_name,"GHOST_TRUNCATION_ORDERS",{i:0 for i in range(0,n)})
    okwfn,okReturns=mbe_helper(order,wfn,mbe_key,sub_frag_key,n,mm,parent,True,
        MAIN_FRAGGER=ghoster_name)
    Returns.update(okReturns)
    return BestWfn,Returns