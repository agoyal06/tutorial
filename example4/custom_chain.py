
######################################
class CustomChain(object):

    # Defining the folder tree structure
    def __init__(self,vaspobj=None):
        self.names=['Relax','wfn']
        self.vasp=vaspobj

    # Defining the Extraction object
    def Extract(self, jobdir):
        from pylada.vasp import MassExtract
        extract = MassExtract(jobdir)
        success={}
        for name in self.names:
            success[name]=MassExtract(jobdir+'/'+name).success
        success=all(success.values())
        extract.success=success
        return extract

    # Creating the workflow
    def __call__(self, structure, outdir=None, **kwargs ):

        from copy import deepcopy
        from os import getcwd
        from os.path import join
        from pylada.misc import RelativePath
        from pylada.error import ExternalRunFailed
        from pylada.vasp.extract import Extract
        from pylada.vasp import Vasp
        from pylada.vasp.relax import Relax

        # make this function stateless.
        structure_ = structure.copy()
        outdir = getcwd() if outdir is None else RelativePath(outdir).path

        ############ Calc 1 ###############
        name  = self.names[0]

        ## functional for Calc 1
        relax = Relax(copy=deepcopy(self.vasp))
        relax.relaxation = "volume ionic cellshape"
        relax.maxiter = 10
        relax.keep_steps = True
        relax.first_trial = { "kpoints": "\n0\nAuto\n10", "encut": 0.9 }
        ## end of the functional
        
        params = deepcopy(kwargs)
        fulldir = join(outdir, name)
        
        ## if this calculation has not been done run it
        output = relax(structure_, outdir=fulldir, **params)
        if not output.success: 
            raise ExternalRunFailed("VASP calculation did not succeed.") 

        ############ Calc 2 ###############  
        name  = self.names[1]

        ## functional for Calc 2
        wfn = Vasp(copy=deepcopy(self.vasp))
        wfn.isym    = 1
        wfn.ismear  = -5
        wfn.nbands=24*len(structure_)
        wfn.kpoints="\n0\nGamma\n4 4 4\n0. 0. 0.\n"
        ## end of the functional
        
        params = deepcopy(kwargs)
        fulldir = join(outdir, name)
        
        ## if this calculation has not been done, run it
        output = wfn(structure_, outdir=fulldir, restart=output, **params)
        if not output.success: 
            raise ExternalRunFailed("VASP calculation did not succeed.") 

        return self.Extract(fulldir)
#########################################
