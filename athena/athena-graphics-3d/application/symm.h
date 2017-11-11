! header file declaring arrays used to test symmetry properties.
	integer nvarsym
	parameter(nvarsym=1)
	common/symm/symmvar(nvarsym,1:iu_bnd,1:ju_bnd,1:ku_bnd,2)

        integer lengthsym
        parameter(lengthsym=nvarsym*iu_bnd*ju_bnd*ku_bnd)
