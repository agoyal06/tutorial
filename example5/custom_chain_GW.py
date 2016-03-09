
######################################

class CustomChain(object):

    # Defining the folder tree structure
    def __init__(self):
        self.names=['Relax','wfn','gw_0','gw_1']

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

    def __call__(self, structure, outdir=None, vasp=None, **kwargs ):

        from copy import deepcopy
        from os import getcwd
        from os.path import join
        from pylada.misc import RelativePath
        from pylada.error import ExternalRunFailed
        from pylada.vasp.extract import Extract, MassExtract
        from pylada.vasp import Vasp
        from pylada.vasp.relax import Relax

        # make this function stateless.
        structure_ = structure.copy()
        outdir = getcwd() if outdir is None else RelativePath(outdir).path

        ############ Calc 1 ###############
        name  = self.names[0]

        ## functional for Calc 1
        relax = Relax(copy=vasp)
        relax.relaxation = "volume ionic cellshape"
        relax.maxiter = 10
        relax.keep_steps = True
        relax.first_trial = { "kpoints": "\n0\nAuto\n10", "encut": 0.9 }
        ## end of the functional

        params = deepcopy(kwargs)
        fulldir = join(outdir, name)

        ## if this calculation has not been done run it
        output = relax(structure_, outdir=fulldir, restart=None, **params)
        if not output.success:
            raise ExternalRunFailed("VASP calculation did not succeed.") 

        ############ Calc 2 ###############  
        name  = self.names[1]

        ## functional for Calc 2
        final = Vasp(copy=vasp)
        final.nbands=24*len(structure_)
        final.kpoints="\n0\nGamma\n2 2 2\n0. 0. 0.\n"
        final.loptics=True
        final.relaxation="static"
        ## end of the functional

        params = deepcopy(kwargs)
        fulldir = join(outdir, name)

        ## if this calculation has not been done, run it
        output = final(structure_, outdir=fulldir, restart=output, **params)
        if not output.success:
            raise ExternalRunFailed("VASP calculation did not succeed.") 

        ############## GW Loop ########################
        for name in self.names[2:]:

            gw = Vasp(copy=vasp)

            gw.kpoints    ="\n0\nGamma\n2 2 2\n0. 0. 0.\n"
            gw.nbands     =24*len(structure_)
            gw.lcharg     = True

            gw.add_keyword('nelm',1)
            gw.add_keyword('algo','gw')
            gw.add_keyword('LMAXFOCKAE',4)
            gw.add_keyword('nomega',64)
            gw.add_keyword('precfock','fast')
            gw.add_keyword('encutgw',50)
            gw.add_keyword('encutlf',50)
            gw.add_keyword('lrpa',False)
            gw.add_keyword('nkred',2)

            params = deepcopy(kwargs)
            fulldir = join(outdir, name)

            ## if this calculation has not been done, run it
            output = gw(structure_, outdir=fulldir, restart=output, **params)
            if not output.success:
                raise ExternalRunFailed("VASP calculation did not succeed.") 
        #########################

        return self.Extract(fulldir)
