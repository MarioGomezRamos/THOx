# THOx
 
** THOx **  is a CDCC code for two-body projectiles, with the possibility of including core excitations. The code has been developed by the Nuclear Theory Group based at the University of Seville. It works with standard fortran compilers (ifort, gfortran,etc)

The code assumes a reaction of the form

  a(=c+v) + A -> c + v + A
  
where "c" and "v" are the fragment constituents (denoted, for convenience, core and valence particles, respectively). 

## Short input description 
###  SYSTEM  namelist: Zv, Av,  Zc, Ac, sn
- Zv, Zc: charge of valence and core particles
- Av, Ac: mass number of valence and core particles
- sn    : intrinsic spin of valence particle

### CORESTATES namelist: spin, parity, ex
 - spin, parity, ex: intrinsic spin, parity and excitation energy of this core state
 
   This namelist can be repeated if more than one core state is present. After the last core state is included, an empty CORESTATE namelist must be included
 
### OUTPUT namelist: wfout, cdcc, verb, solapout. 
    Controls the information printed out in stdout and auxiliary files

  - wfout(:): array containing the indeces of the projectile eigenstates whose wavefunctions will be printed out.
  - solapout(:): array containing the indexes of the projectile eigenstates whose overlap with the scattering states will be printed. 
  - cdcc: 
  - verb: controls the amount of output information
  
   
### GRID namlist: ng, rmin, rint, rmax, dr, rlast
 - rmin, rmax, dr: radial grid for core+valence relative motion
 - rint: when calculating the c+v scattering states, maximum radius up to which the equation will be integrated numerically. For r>rint, the asymptotic form will be assumed

### &POTENTIAL namelist: ptype, ap, at, Vl0(:), r0, a0, Vso, rso, aso, pcmodel, lambda, kband, lpot, cptype, Vcp0, rcp0, acp, delta 
- ptype: potential type
- ap, at: core and valence masses used in the radii conversion from reduced to physical radii
- Vl0(:): array for potential depths for each partial wave of central potential
- r0, a0: reduced radius and diffuseness of central potential
- a0, diffuseness paramter Vso, rso, aso, pcmodel, lambda, kband, lpot, cptype, Vcp0, rcp0, acp, delta 

### PAULI namelist: 
          
### JPSET namelist: bastype, mlst, gamma, bosc, nho, nsp, exmin, exmax, bas2, JTOT, PARITY, lmax, 
- bastype= index to specify the basis type: 
  0 = HO
  1 = THO
  2 = bins
- bosc: oscillator parameters used in the HO and THO bases.
- mlst, gamma: parameters for the local scale transformation (LST) in the THO basis  
- eta: 
- nsp: number of single-particle eigenvalues to keep in the diagonalization of the full core+valence Hamiltonian
- nbins: number of bins (bastype=2)
- nk : number or scattering states within each bin interval used to construct the bin wavefunction
- tres: included T-matrix in the weight function of the bin wfs
- inc : incoming channel for bins
- ehat (T/F): to use mean bin energies (otherwise midpoint), 
- filewf: external file for wfs
- wcut: mininum weight per channel to be retained (default 1) 
- vscale: scaling factor for central potential for this particular J/Pi set 

       
### SCATWF namelist: ifcont, emin, emax, nk, inc, jset
 - ifcont (T/F): true, compute scattering states for the core+valence system
 - emin, emax, nk: minimum energy, maximum energy and number of scattering functions to be computed
 - inc: for multichannel states, index of incoming channel
 - jset: index of JPSET to specify the angular momentum and spin of the scattering states to be computed


### BELAMBDA namelist: ifbel, uwfgsfile, lambda, jset, emin, emax, nk

### REACTION namelist: elab, namep, mp, mt, namet, zt, jt 

### TRANS namelist: skip, rcc, writeff 
- skip  (T/F): if true, skip calculation of transition potentials
- rcc: radius for Coulom potential
- writeff: if true, print coupling potential in ff.fr file (Fresco format) and states information in states.fr (to be read with FRESCO using readstates variable)

### NUMEROV namelist: hcm, rmaxcc, hort, method, jtmin, jtmax, skip
  - method: method of solution of the CC equations. Available options are:
  
         0=PC-numerov, 
         1=ENA with 5 terms in Cosh[Sqrt[T]]
         2=ENA with 5 terms in Cosh[Sqrt[T]], only diagonal 
         3=Raynal
         4=Modified Numerov used in Fresco
         
   - hcm: radial step for projectile-target coordinate for solving the CC equations
   - jtmin, jtmax: min, max total angular momentum for solving the CC equations
   - hort: if nonzero, uses a stabilization procedure of the CC equations (see long text description)
   - skip (T/F): if true, skips this section
   
### XSECTIONS namelist:  fileamp, thmin, thmax, dth, thcut, doublexs, triplexs, phixs, icore, ner, ermin, ermax, jsets(:) 
- fileamp: if defined, filename of file containing scattering amplitudes for the computation of double and triple differential cross sections. If no defined, these cross sections are calculated with the amplitudes previously calculated. 
- thmin, thmax, dth: angular grid for printing scattering amplitudes and cross sections
- thcut: angular cutoff for computation of relative energy distribution
- doublexs (T/F): if true, evaluates and prints double differential cross sections as a function of core-valence relative energy
- triplexs (T/F): if true, evaluates and prints triple differential cross sections. The energy and angular grids for the fragments are later read using the framework, gridnener, gridthetac, gridthetav, gridphi namelists. 
- phixs (T/F): if true, prints five-fold differential cross sections without phi angle integration
- icore: index of core state for the calculation of double and triple differential cross sections
- ner, ermin, ermax
- jsets(:): array of jpsets indexes to be included in the computation of double and triple differential cross sections 
   
   
### FRAMEWORK namelist: sys, idet 
- sys=lab/com
- idet: specifies if the energy of triple differential cross section correspods to the core (idet=1) or valence (idet=2) particles

### GRIDENER namelist: Enlow, Enup, dEn     
- Energy grid of detected particle

### GRIDTHETAC namelist tcl, tcu, dtc      
 - Theta angle grid of core particle
 
### GRIDTHETAV namelist: tvl, tvu, dtv
 - Theta angle grid of valence particle
 
### GRIDPHI namelist: phil, phiu, dphi  
 - Phi angle grid. Computed cross sections assume that phi=0 for the core. 
 
