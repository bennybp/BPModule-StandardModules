def mbe_helper(order,wfn,mbe_key,frag_key,n,mm,parent=0):
    """
    A function designed to quickly compute useful MBE quantities

        order (int): The derivative order of the MBE to compute
        wfn (psr.datastore.wavefunction): The wavefunction to give to the MBE 
        mbe_key (string): The key for the MBE to run
        frag_key (string): The key for the fragmenter
        n (int): The highest truncation order to generate statistics for
        mm (psr.modulemanager.ModuleManager): The module manager to use
        parent (int): The ID of the parent module defaults to 0
 
        returns a dictionary of names of energetic quantities indexed by their
                names
    """
    Returns={}
    def TotalEName(i):
        return "Total "+str(i)+"-body derivative"
    def IntEName(i):
        return "Total of "+str(i)+"-body interactions"
    
    for i in range(n,0,-1):
        temp_mbe_name=mm.generate_unique_key()
        temp_frag_name=mm.generate_unique_key()
        mm.duplicate_key(mbe_key,temp_mbe_name)
        mm.duplicate_key(frag_key,temp_frag_name)
        mm.change_option(temp_frag_name,"TRUNCATION_ORDER",i)
        mm.change_option(temp_mbe_name,"FRAGMENTIZER",temp_frag_name)
        mymod=mm.get_module(temp_mbe_name,parent)
        oldwfn,Result=mymod.deriv(order,wfn)
        if i==n:BestWfn=oldwfn
        Returns[TotalEName(i)]=Result
        if i<n:
            Returns[IntEName(i+1)]=[Returns[TotalEName(i+1)][j]-Result[j] 
                for j in range(0,len(Result))]
    return BestWfn,Returns

MP2D={'Total 3-body derivative': [-225.02122921895318], 'Total 1-body derivative': [-224.9963830035146], 'Total 2-body derivative': [-225.01705899123075], 'Total of 2-body interactions': [-0.020675987716145983], 'Total of 3-body interactions': [-0.004170227722426034]}

HFD={'Total 2-body derivative': [-224.90851105000854], 'Total 3-body derivative': [-224.91252968717606], 'Total of 2-body interactions': [-0.018281847793275574], 'Total of 3-body interactions': [-0.004018637167519046], 'Total 1-body derivative': [-224.89022920221527]}
