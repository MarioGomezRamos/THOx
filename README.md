# THOx
 v1.0
** THOx **  is a CDCC code for two-body projectiles, with the possibility of including core excitations. The code has been developed by the Nuclear Theory Group based at the University of Seville.

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

  - wfout(:): array containing the index of the projectile eigenstates whose wavefunctions will be printed out 
  - cdcc: 
  - verb: 
  - solapout: 
   
### GRID namlist: ng, rmin, rint, rmax, dr, rlast

### &POTENTIAL namelist: ptype, ap, at, Vl0(:), r0, a0, Vso, rso, aso, pcmodel, lambda, kband, lpot, cptype, Vcp0, rcp0, acp, delta 

### PAULI namelist: 
          
###  JPSET namelist: bastype, mlst, gamma, bosc, nho, nsp, exmin, exmax, bas2, JTOT, PARITY, lmax, 
       
### SCATWF namelist: ifcont, emin, emax, nk, inc, jset       

### BELAMBDA namelist: ifbel, uwfgsfile, lambda, jset, emin, emax, nk

### REACTION namelist: elab, namep, mp, mt, namet, zt, jt 

### TRANS namelist: skip, rcc, writeff 

### NUMEROV namelist: hcm, rmaxcc, hort, method, jtmin, jtmax, skip
  - method: method of solution of the CC equations. Available options are:
  
         0 PC-numerov, 
         1=ENA with 5 terms in Cosh[Sqrt[T]]
         2=ENA with 5 terms in Cosh[Sqrt[T]], only diagonal 
         3=Raynal
         4=Modified Numerov used in Fresco
         
   - hcm: radial step for projectile-target coordinate for solving the CC equations
   - jtmin, jtmax: min, max total angular momentum for solving the CC equations
   - hort: if nonzero, uses a stabilization procedure of the CC equations (see long text description)
   - skip (T/F): if true, skips this section
   
### XSECTIONS namelist:  fileamp, thmin, thmax, dth, thcut, doublexs, triplexs, phixs, icore, ner, ermin, ermax, jsets(:)  
   
   
### FRAMEWORK namelist: sys, idet 

### GRIDENER namelist: Enlow, Enup, dEn      

### GRIDTHETAC namelist tcl, tcu, dtc      

### GRIDTHETAV namelist: tvl=0.0, tvu, dtv

### GRIDPHI namelist: phil, phiu, dphi  
   
[TO BE COMPLETED] 
