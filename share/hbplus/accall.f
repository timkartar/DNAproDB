      program accall
c
c     --  AIM:-   
c     --  Input a Brookhaven entry file and output an 
c     --  PDB format file after filtering/cleaning, including Van der
c     --  Waal radii, contained in an external file "vdw.radii". 
c     
c     --  INPUT:- 
c     --  PDB format file, van der Waal radii file
c     
c     --  OPTIONS:- 
c     --  Inclusion of non-standard amino acids, het-atoms, waters
c     --  etc. Flagging of missing residues, chain-breaks,
c     --  non-standard atoms names, missing atoms.
c     
c     --  AUTHOR: S. Hubbard 3/92. EMBL.
c     
c     --  PARAMETERS: maxa = max atoms per residue
c     --              maxr = max number of "standard" residue types 
c     --              maxs = max number of atoms in PDB file
c
      parameter    (
     -               maxa=100,   
     -               maxr=100,
     -               maxs=20000,
     -               maxx=2000
     -             )
c
c     -- functions
c
      integer        readstring, parse, fopen, what_atom
      real           readfloat
c
c     -- variables
c
      integer        i, j, k, 
     -               l(256), len, 
     -               flen, vlen,
     -               numats(maxr), 
     -               resindex(maxs),
     -               atomtype(maxs),
     -               rty(maxx),
     -               num_chains,
     -               num_res, 
     -               nats,
     -               atype,
     -               nacids

      real           vradii(maxr,maxa),
     -               xyz(maxs,3),
     -               rads(maxs),
     -               accs(maxs),
     -               vdw,
     -               hyrad,
     -               probe,
     -               zslice,
     -               ressums(maxx,5,2), 
     -               tsums(5)

      character*1    b, alt, firsta, ch_nam(100)
      character*3    res,  aacids(maxr), rlab(3)
      character*4    atom, anames(maxr,maxa)
      character*6    reslabel(50)
      character*10   last, resnam(maxx)
      character*30   label(maxs)
      character*256  fname, vname     
      character*256  card, c(256)

      logical        hetas, falt, hydro, hgen,
     -               ok, start, wwaters, resok, aok

      data hyrad/1.00/
      data rlab/'RES','HEM','HOH'/
c
c     -- defaults
c      
      hetas=.false.
      hydro=.false.
      wwaters=.false.
c
c     -- Get USER directives
c
      do while ( readstring(5,card,len).ge.0 )
         n = parse(card,len,' ',c,l)
         call tolow(c(1),l(1))
         if(c(1)(1:4).eq.'pdbf')then
            fname=c(2)
            flen=l(2)
         elseif(c(1)(1:4).eq.'vdwf')then
            vname=c(2)
            vlen=l(2)
         elseif(c(1)(1:4).eq.'prob')then
            probe=readfloat(c(2),l(2))
         elseif(c(1)(1:4).eq.'zsli')then
            zslice=readfloat(c(2),l(2))
         elseif(c(1)(1:4).eq.'heta')then
            hetas=.true.
         elseif(c(1)(1:4).eq.'hydr')then
            hydro=.true.
         elseif(c(1)(1:4).eq.'wate')then
            wwaters=.true.
         endif
      enddo
c
c     --   open files
c
c
c     -- Read in VDW radii for all residues/atoms
c
      call vanin ( 
     -     vname,
     -     vlen,
     -     nacids, 
     -     aacids, 
     -     anames, 
     -     numats, 
     -     vradii,
     -     maxr,
     -     maxa
     -     )     

      i = index(fname,'.')-1
      k = 1
      ok=.false.
      do j = i, 1, -1
         if( fname(j:j).eq.'/' .and. .not.ok )then
            k=j+1
            ok=.true.
         endif
      enddo

      if( fopen(1,fname,flen,'old').eq.0 ) then
         STOP'ERROR: unable to open PDB file'
      endif

      open ( 
     -      unit = 2,                                                 
     -      file = fname(k:i)//'.asa',
     -      status = 'unknown'
     -     )
      open ( 
     -      unit = 3,       
     -      file = fname(k:i)//'.rsa',
     -      status = 'unknown'
     -     )
      open ( 
     -      unit = 4,                                                 
     -      file = fname(k:i)//'.log',
     -      status = 'unknown'
     -     )

      write(4,'(a)')' ACCALL - Accessibility calculations'
      write(4,'(2a)')' PDB FILE INPUT ',fname(1:flen)
      write(4,'(a,f6.2)')' PROBE SIZE     ',probe
      write(4,'(a,f6.3)')' Z-SLICE WIDTH  ',zslice
      if(hetas)then
         write(4,'(a)')' INCL HETATOMS'
      else
         write(4,'(a)')' EXCL HETATOMS'
      endif
      if(hydro)then
         write(4,'(a)')' INCL HYDROGENS'
      else
         write(4,'(a)')' EXCL HYDROGENS'
      endif
      if(wwaters)then
         write(4,'(a)')' INCL WATERS'
      else
         write(4,'(a)')' EXCL WATERS'         
      endif
      write(4,'(a,i3,a)')' READVDW ',nacids,' residues input'
c
c --  Initialise variables/logicals
c
      falt=.false.
      start=.true.                       
      last='          '                                               
c
c -- Read data & decode
c
      nats = 0

      do while ( readstring(1,card,len).ge.0 )
         atype=0
         if(card(1:4).eq.'ATOM')atype=1
         if(card(1:6).eq.'HETATM')atype=2
         if(card(18:20).eq.'HOH')atype=3
         if(atype.eq.1.or.
     -        (atype.eq.2.and.hetas).or.
     -        (atype.eq.3.and.wwaters) )then
c
c     -- Ignore Alternate positions, other than blanks or 1st
c     -- encountered
c
            alt=card(17:17)
            if(alt.ne.' ')then
               if(.not.falt)then
                  firsta=alt
                  falt=.true.
               endif
               if(alt.ne.firsta)goto5
            endif            
c     
c     -- Ignore hydrogens & deuteriums (unless flagged)
c     
            if(card(14:14).eq.'H'.or.card(14:14).eq.'D')then
               if(.not.hydro)goto5
               vdw=hyrad
               goto6
            endif
c
c     -- Next atom
c
            nats = nats + 1
c
c     -- First residue ?
c
            if(start)then
               start = .false.
               ch_nam(1) = card(22:22)
               num_chains = 1
            endif
c
c     -- New residue ?
c
            if ( last.ne.card(18:27) ) then
               last = card(18:27)
               num_res = num_res + 1
               res = card(18:20)
               resnam(num_res)=last
               call which3(res,aacids,i,nacids,maxr,resok)                              
               if(.not.resok)write(4,102)card(18:27)               
               rty(num_res)=atype
            endif
c
c     -- Get atom type 
c
            atom=card(13:16)
            atomtype(nats) = what_atom(atom)
c
c     -- Assign radius to atom
c     -- Special case(s)
c
            if(atom.eq.'OXT')then
               vdw=1.40
               goto6
            endif
c
c     -- known residue type
c            
            if(resok)then
               vdw = 0.0
               call ratom(atom,anames,i,numats,maxr,maxa,aok,j)
            endif
c
c     -- Not OK, then try atoms in all residue types
c
            if( .not.resok .or. .not.aok )then
               call gatom(atom,anames,nacids,numats,maxr,maxa,aok,i,j)               
               write(4,104)atom,card(18:27)
               if(aok)write(4,106)atom,card(18:27),
     -                vradii(i,j),aacids(i)
            endif
c
c     -- Still not OK, make a guess
c
            if(.not.aok)then
               call vguess(atom,vdw)
               write(4,108)atom,card(18:27),vdw
            else
               vdw=vradii(i,j)
            endif
c
c     -- Store data
c
 6          rads(nats) = vdw
            label(nats) = card(1:30)
            resindex(nats)= num_res
            read(card,'(30x,3f8.3)')(xyz(nats,k),k=1,3)
         endif
 5       continue
      enddo
c
c     -- output
c
      close(1)
      write(4,'(a)')' ADDED VDW RADII'
      write(4,110)num_chains,num_res,nats
c
c     -- calculate atomic accessibilities
c
      call solva(nats,xyz,rads,accs,probe,zslice,maxs)

      do i = 1, nats
C         write(2,'(a30,3f8.3,f8.3,1x,f5.2)')
C     -        label(i),
C     -        (xyz(i,j),j=1,3),
C     -        accs(i),
C     -        rads(i)
         write(2,'(a30,3f8.3,10x,f8.3)')
     -        label(i),
     -        (xyz(i,j),j=1,3),
     -        accs(i)
      enddo
      
      write(4,'(a)')' CALCULATED ATOMIC ACCESSIBILITES'

      call summer(
     -     nats,
     -     accs,
     -     maxs,
     -     maxx,
     -     atomtype,
     -     resindex,
     -     resnam,
     -     ressums,
     -     tsums
     -     )
      
      write(4,'(a)')' SUMMED ACCESSIBILITIES OVER RESIDUES'

      write(3,120)
      write(3,125)
      write(3,130)
      do i = 1, resindex(nats)
         write(3,150)
     -        rlab(rty(i)),
     -        resnam(i),
     -        (ressums(i,j,1),ressums(i,j,2),j=1,5)
      enddo
      write(3,160)(tsums(i),i=1,5)
      
c----------------------------------------------------------       
 102  format(' UNKNOWN residue type.............> ',a10)
 104  Format(' NON-STANDARD atom.',a4,' in residue> ',a10)
 106  format(' ASSUMED vdw of ',a4,' in ',a10,' = ',f5.2,
     -       ' (same as ',a3,')')
 108  format(' GUESSED vdw of ',a4,' in ',a10,' = ',f5.2)
 110  format(' CHAINS   ',i5,/,
     -       ' RESIDUES ',i5,/,
     -       ' ATOMS    ',i5)
 120  format('REM  File of summed (Sum) and % (per.)',
     -     ' accessibilities for ',a)
 125  format('REM RES _ NUM      All atoms   Non P side  ',
     -     ' Polar Side   Total Side   Main Chain')
 130  format('REM                ABS   REL    ABS   REL',
     -     '    ABS   REL    ABS   REL    ABS   REL')
 150  format(a3,1x,a10,1x,5(f7.2,f6.1))
 160  format('END  Absolute sums over accessible surface ',/,
     -     'TOTAL',8x,5(f8.1,5x))	
c----------------------------------------------------------
      end

      subroutine vguess ( atom, vdw )
      real         vdw
      character*4  atom

      vdw=1.80
c
c -- Make a guess then !
c
      if(atom(2:2).eq.'C')vdw=1.80
      if(atom(2:2).eq.'N')vdw=1.60
      if(atom(2:2).eq.'S')vdw=1.85
      if(atom(2:2).eq.'O')vdw=1.40
      if(atom(2:2).eq.'P')vdw=1.90
      if(atom(1:2).eq.'CA')vdw=2.07
      if(atom(1:2).eq.'FE')vdw=1.47
      if(atom(1:2).eq.'CU')vdw=1.78
      if(atom(1:2).eq.'ZN')vdw=1.39
      if(atom(1:2).eq.'MG')vdw=1.73

      return
      end

      subroutine gatom(atom,anames,nres,nats,maxr,maxa,ok,ir,ia)
      integer      maxr, maxa
      integer      nres, nats(maxr), ir, ia
      character*4  atom, anames(maxr,maxa)
      logical      ok
      
      ok=.false.
      ir=0
      do while ( ir.lt.nres .and. .not.ok )
         ir=ir+1
         call ratom(atom,anames,ir,nats,maxr,maxa,ok,ia)
      enddo

      return
      end

      subroutine ratom(atom,anames,ires,nats,maxr,maxa,ok,find)
c
c --  Checks to see if a "standard" atom name has been read in.
c --  Standard atom names are read from the file "vdw.radii", for
c --  the defined residue types therein. OK=.true. if found, and
c --  residue type = ires, and find = integer identifier of atom.
c
      integer      ires, nats(maxr), find
      character*4  atom, anames(maxr,maxa)
      logical      ok

      ok=.false.
      find=0
      if(ires.eq.0.or.ires.gt.maxr)return

      do i = 1, nats(ires)
         if(atom.eq.anames(ires,i))then
            find=i
            ok=.true.
            return
         endif
      enddo

      return
      end

      subroutine which3(res,acids,ires,nacids,maxr,ok)
c
c -- Search array "acids" for existence of residue "res".
c -- OK = .true. if found, and index of res returned in "ires".
c
      integer      ires, nacids
      character*3  acids(maxr), res
      logical      ok
      
      ires=nacids+1
      do i = 1, nacids
         if(res.eq.acids(i))then
            ires=i
            ok=.true.
            return
         endif
      enddo

      ok=.false.
      return
      end

      integer function fopen( iochan, filename, len, fstat )
      integer       iochan, len
      character*(*) filename, fstat

      open (
     -      unit   = iochan,
     -      file   = filename(1:len),
     -      status = fstat,
     -      err    = 100
     -     )

      fopen = 1
      return
 100  fopen = 0
      return
      end

      integer function readstring(file,card,len)
      integer       file, len
      character*(*) card

      if(file.gt.200)STOP'ERROR: file number too large'
      read(file,'(q,a)',err=100,end=100)len,card(1:len)
      card(len+1:len+1)=char(0)
      readstring=len
      return
 100  readstring=-1
      return
      end

      integer function parse(card,length,separator,chars,clen)
      integer       charpos(256,2), clen(256), length
      character*1   separator
      character*(*) card, chars(256)
      logical       search

      parse=0
      i=0
      search=.false.
      do while(i.lt.length)
         i = i + 1
         if(.not.search)then
            if(card(i:i).ne.separator)then
               parse = parse + 1
               charpos(parse,1) = i
               search=.true.
            endif
         else
            if(card(i:i).eq.separator)then
               charpos(parse,2) = i-1
               search=.false.
            endif
         endif
      enddo
      if(search)charpos(parse,2)=length
      do i = 1, parse
         chars(i)=card(charpos(i,1):charpos(i,2))
         clen(i)=charpos(i,2)-charpos(i,1)+1
      enddo
      return
      end

      real function readfloat(card,len)
      integer       len
      real          value
      character*(*) card

      read(card(1:len),*,err=100)value
      readfloat=value
      return
 100  readfloat=-999.9
      return
      end

      subroutine tolow(text,len)
      integer       len
      character*(*) text
      character*1   alph(26,2)
      logical       ok

      data alph/'A','B','C','D','E','F','G','H','I','J',
     -          'K','L','M','N','O','P','Q','R','S','T',
     -          'U','V','W','X','Y','Z',
     -          'a','b','c','d','e','f','g','h','i','j',
     -          'k','l','m','n','o','p','q','r','s','t',
     -          'u','v','w','x','y','z'/

      do i = 1, len
         j=0
         ok=.true.
         do while( j.lt.26.and.ok ) 
            j = j + 1
            if ( alph(j,1).eq.text(i:i) ) then
               ok=.false.
               text(i:i)=alph(j,2)
            endif
         enddo
      enddo

      return
      end

      subroutine solva(nats,xyz,rads,accs,probe,zslice,maxs)
c**
c*************************************************************
c**  SOLVA - LEE & RICHARDS TYPE ACCESSIBLITY CALCULATIONS **
c*************************************************************
c**
c** Calculate accessible surface area for a group of atoms.
c** The accessible area for a given atom is calculated by the 
c** formula:
c**     (arcsum) x (atom radius+probe radius) x (deltaz)
c** Numerical integration is carried out over z. in each z-
c** section, the arcsum for a given atom is the arclength of 
c** the circle (intersection of the atom sphere with the z-
c** section) that is not interior to any other atom circles
c** in the same z-section.
c**
c*************************************************************
c**
c**  error parameter  - this gives accuracy of calculation
c**                   - suitable values are 0.01 (high 
c**                     accuracy) to 0.1 (low accuracy)
c**                   - in detail the z sections are spaced
c**                     at about error*diameter of atom
c**
c**  probe size       - radius of probe in angstroms
c**                   - suitable value for water = 1.4
c**
c*************************************************************
c=============================================================
c
c PARAMETERS:  
c ncube = maximum number of cubes allowed for placing of atoms
c nac   = maximum number of atoms per cube
c nint  = maximum number of sphere intersections
c smax  = maximum number of atoms
c                   
      integer     ncube,       nac,     nint,      smax,  maxs
      parameter ( ncube=10000, nac=150, nint=2000, smax=20000 )  
c
c     the following are dimensioned to the max no of atoms (maxs)
c
      integer cube(smax)
      real    xyz(maxs,3),
     -        rads(maxs),
     -        rad(smax),
     -        radsq(smax) 
c     
c     the following are dimensioned to the max no of intersections
c     of neighbouring spheres (nint)
c     
      dimension
     -        inov(nint),
     -        tag(nint),
     -        arci(nint),
     -        arcf(nint),
     -        dx(nint),
     -        dy(nint),
     -        d(nint),
     -        dsq(nint)

      integer    natm(nac,ncube), 
     -           itab(ncube)
      integer    nats, tag    
      real       probe, zslice, accs(maxs), b, c
      real       trig_test

      data xmin,ymin,zmin,xmax,ymax,zmax/3*9999.,3*-9999./
c
c     initialise variables, constants
c
      ict=nint
      pi=acos(-1.0)
      pix2=2.0*pi
c
c     -- Radius of an atom sphere = atom radius + probe radius
c     -- Find maxima and minima
c
      rmax=0.0
      karc=ict
      do i = 1, nats
         rad(i)   = rads(i) + probe
         radsq(i) = rad(i)**2
         if (rad(i).gt.rmax)rmax = rad(i)
         if (xmin.gt.xyz(i,1))  xmin = xyz(i,1)
         if (ymin.gt.xyz(i,2))  ymin = xyz(i,2)
         if (zmin.gt.xyz(i,3))  zmin = xyz(i,3)
         if (xmax.lt.xyz(i,1))  xmax = xyz(i,1)
         if (ymax.lt.xyz(i,2))  ymax = xyz(i,2)
         if (zmax.lt.xyz(i,3))  zmax = xyz(i,3)
      enddo
c
c     rmax = max diameter
c
      rmax = rmax*2.
c
c     -- Cubicals containing the atoms are setup. 
c     -- The dimension of an edge equals the largest atom sphere radius
c     -- The cubes have a single index
c     -- Minimum of 3 by 3 cubic grid
c     -- EXIT if max cubes exceeded
c
      idim=(xmax-xmin)/rmax+1.
      if(idim.lt.3)idim=3
      jidim=(ymax-ymin)/rmax+1.
      if(jidim.lt.3)jidim=3
      jidim=idim*jidim
      kjidim=(zmax-zmin)/rmax+1.
      if(kjidim.lt.3)kjidim=3
      kjidim=jidim*kjidim
      if(kjidim.gt.ncube)STOP'SOLVA_ERROR: max cubes exceeded'
c
c     -- Prepare upto ncube cubes each containing upto nac atoms. The cube index
c     -- is kji. The atom index for each cube is in itab
c
      do l = 1, ncube
         itab(l)=0
      enddo
      do l = 1, nats
         i = (xyz(l,1)-xmin)/rmax+1.
         j = (xyz(l,2)-ymin)/rmax
         k = (xyz(l,3)-zmin)/rmax
         kji = k*jidim + j*idim+i
         n = itab(kji)+1
         if(n.gt.nac)STOP'SOLVA_ERROR: max atoms per cube exceeded'
         itab(kji) = n
         natm(n,kji) = l
         cube(l) = kji
      enddo
c
c     -- Process each atom in turn
c
      nzp=1./zslice+0.5
      do  ir = 1, nats
         kji=cube(ir)
         io=0
         area=0.
         xr=xyz(ir,1)
         yr=xyz(ir,2)
         zr=xyz(ir,3)
         rr=rad(ir)
         rrx2=rr*2.
         rrsq=radsq(ir)
c
c     -- Find the 'mkji' cubes neighboring the kji cube
c
         do k = -1, 1, 1
            do j = -1, 1, 1
               do i = -1, 1, 1
                  mkji=kji+k*jidim+j*idim+i
                  if(mkji.ge.1)then
                     if(mkji.gt.kjidim)goto14
                     nm=itab(mkji)
                     if(nm.ge.1)then
c     
c     -- record the atoms in inov that neighbor atom ir
c
                        do m = 1, nm
                           in=natm(m,mkji)
                           if (in.ne.ir)then
                              io=io+1
                              if (io.gt.ict)then
                                 STOP'SOLVA_ERROR: intrsctns > max'
                              endif
                              dx(io)=xr-xyz(in,1)
                              dy(io)=yr-xyz(in,2)
                              dsq(io)=dx(io)**2+dy(io)**2
                              d(io)=sqrt(dsq(io))
                              inov(io)=in
                           endif
                        enddo
                     endif
                  endif
               enddo
            enddo
         enddo

 14      if(io.ge.1)then
c     
c     z resolution determined
c     
            zres=rrx2/nzp
            zgrid=xyz(ir,3)-rr-zres/2.            
         else
            area=pix2*rrx2
            goto 18
         endif
c     
c     section atom spheres perpendicular to the z axis
c
         do i = 1, nzp
            zgrid=zgrid+zres
c
c     find the radius of the circle of intersection of 
c     the ir sphere on the current z-plane
c
            rsec2r=rrsq-(zgrid-zr)**2
            rsecr=sqrt(rsec2r)
            do k = 1, karc
               arci(k)=0.0
            enddo
            karc=0
            do j = 1, io
               in=inov(j)
c     
c     find radius of circle locus
c     
               rsec2n=radsq(in)-(zgrid-xyz(in,3))**2
               if (rsec2n.le.0.0) goto10
               rsecn=sqrt(rsec2n)
c
c     find intersections of n.circles with ir circles in section
c
               if (d(j).ge.rsecr+rsecn) goto10
c
c     do the circles intersect, or is one circle completely inside the other?
c
               b=rsecr-rsecn
               if (d(j).gt.abs(b)) goto20
               if (b.le.0.0) goto9
               goto10
c
c     if the circles intersect, find the points of intersection
c
 20            karc=karc+1
               if(karc.ge.ict)then
                  STOP'SOLVA_ERROR: max intersections exceeded2'
               endif
c
c     Initial and final arc endpoints are found for the ir circle intersected
c     by a neighboring circle contained in the same plane. The initial endpoint
c     of the enclosed arc is stored in arci, and the final arc in arcf
c     law of cosines
c     
               trig_test=(dsq(j)+rsec2r-rsec2n)/(2.*d(j)*rsecr)
               if(trig_test.ge.1.0)trig_test=0.99999
               if(trig_test.le.-1.0)trig_test=-0.99999
               alpha=acos(trig_test)
c     
c     alpha is the angle between a line containing a point of intersection and
c     the reference circle center and the line containing both circle centers
c     
               beta=atan2(dy(j),dx(j))+pi
c     
c     beta is the angle between the line containing both circle centers and the x-axis
c     
               ti=beta-alpha
               tf=beta+alpha
               if(ti.lt.0.0)ti=ti+pix2
               if(tf.gt.pix2)tf=tf-pix2
               arci(karc)=ti
               if(tf.ge.ti)go to 3
c     
c     if the arc crosses zero, then it is broken into two segments.
c     the first ends at pix2 and the second begins at zero
c     
               arcf(karc)=pix2
               karc=karc+1
 3             arcf(karc)=tf
 10         enddo
c
c     find the accessible surface area for the sphere ir on this section
c     
            if(karc.ne.0)goto19
            arcsum=pix2
            go to 25
c     
c     The arc endpoints are sorted on the value of the initial arc endpoint
c
 19         call sortag(arci(1),karc,tag)
c
c***************************************
c     calculate the accessible area
c***************************************
c     
            arcsum=arci(1)
            t=arcf(tag(1))
            if(karc.eq.1) go to 11
            do k = 2, karc
               if(t.lt.arci(k))arcsum=arcsum+arci(k)-t
               tt=arcf(tag(k))
               if(tt.gt.t)t=tt
            enddo
 11         arcsum=arcsum+pix2-t
c     
c     The area/radius is equal to the accessible arc length x the section thickness.
c     
 25         parea=arcsum*zres
c     
c     Add the accessible area for this atom in this section to the area for this
c     atom for all the section encountered thus far
c     
            area=area+parea
 9       enddo
c     
c     scale area to vdw shell
c     
 18      b=area*rr
         accs(ir)=b
c------------------------------------------------------------------
c The following line converts from accessible to contact surface
c         c=(b*(rad(ir)-probe)**2)/(rad(ir)**2)
c------------------------------------------------------------------     
      enddo

      write(4,'(a)')' SOLVA: PROGRAM ENDS CORRECTLY'
      return
      end

      subroutine sortag(a,n,tag)
      integer tag,tg
      dimension a(n),iu(16),il(16),tag(n)
      do i = 1, n
         tag(i)=i
      enddo
      m=1
      i=1
      j=n
 5    if(i.ge.j) go to 70
 10   k=i
      ij=(j+i)/2
      t=a(ij)
      if(a(i).le.t) go to 20
      a(ij)= a(i)
      a(i)=t
      t=a(ij)
      tg=tag(ij)
      tag(ij)=tag(i)
      tag(i)=tg
 20   l=j
      if(a(j).ge.t) go to 40
      a(ij)=a(j)
      a(j)=t
      t=a(ij)
      tg=tag(ij)
      tag(ij)=tag(j)
      tag(j)=tg
      if(a(i).le.t) go to 40
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      tg=tag(ij)
      tag(ij)=tag(i)
      tag(i)=tg
      go to 40
 30   a(l)=a(k)
      a(k)=tt
      tg=tag(l)
      tag(l)=tag(k)
      tag(k)=tg
 40   l=l-1
      if(a(l).gt.t) go to 40
      tt=a(l)
 50   k=k+1
      if(a(k).lt.t) go to 50
      if(k.le.l) go to 30
      if(l-i.le.j-k) go to 60
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 80
 60   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 80
 70   m=m-1
      if(m.eq.0) return
      i=il(m)
      j=iu(m)
 80   if(j-i.ge.1) go to 10
      if(i.eq.1) go to 5
      i=i-1
 90   i=i+1
      if(i.eq.j) go to 70
      t=a(i+1)
      if(a(i).le.t) go to 90
      tg=tag(i+1)
      k=i
 100  a(k+1)=a(k)
      tag(k+1)=tag(k)
      k=k-1
      if(t.lt.a(k)) go to 100
      a(k+1)=t
      tag(k+1)=tg
      go to 90
      end

      subroutine summer(
     -     nats,
     -     accs,
     -     m1,m2,
     -     atomtype,
     -     resindex,
     -     resnam,
     -     ressums,
     -     tsums
     -     )

c     -- 	program to sum atomic accessibilities by residue.
c     --	copes with atom and hetatom records
c     --	produces relative accessibilities for the 20 common aminos
c     --	ouput written to .rsa file (channel 4)      

C RUPP rmaxr is implicitly typed as real - needs explicit declaration as integer in 
C RUPP statement below to avoid conflict with calls to subroutine which3( ).  
      integer       fopen, readstring,rmaxr
      parameter    (rmaxr=100)
      integer       nats, 
     -              atomtype(m1), 
     -              resindex(m1),
     -              rindex(2000),
     -              nacids
      real	    ressums(m2,5,2),
     -              tsums(5), 
     -              accs(m1),
     -		    standarea(rmaxr,5)
      character*256 line
      character*10  resnam(m2)
      character*3   res, acids(rmaxr)
      logical       stand, ok
c
c     if "standard.data" exists in CWD, read them in.
c
      stand=.false.
      if( fopen(1,'standard.data',13,'old').ne.0 )then
         write(3,'(a,a)')'REM  Relative accessibilites read from',
     -                  ' external file "standard.data"'
         stand=.true.
         i=0
         do while( readstring(1,line,len).ge.0.and.i.lt.rmaxr )
            if(line(1:4).eq.'ATOM')then
               i=i+1
               acids(i)=line(13:15)
               read(line(17:23),'(f7.2)')standarea(i,1)
               read(line(30:36),'(f7.2)')standarea(i,2)
               read(line(43:49),'(f7.2)')standarea(i,3)
               read(line(56:62),'(f7.2)')standarea(i,4)
               read(line(69:75),'(f7.2)')standarea(i,5)
            endif
         enddo
         write(4,'(a,i3,a)')
     -        ' RELATIVE (STANDARD) ACCESSIBILITIES READFOR ',
     -        i,' AMINO ACIDS'
      else
         write(4,'(a)')' NO STANDARD VALUES INPUT'
      endif

      nacids=i
      do i = 1, resindex(nats)
         rindex(i)=0
         if(stand)then
            res=resnam(i)(1:3)
            call which3(res,acids,ires,nacids,rmaxr,ok)
            if(ok)rindex(i)=ires
         endif
      enddo
c
c     -- sum the values
c     
      do i = 1, nats
         ir=resindex(i)
         tsums(1)=tsums(1)+accs(i)
         ressums(ir,1,1)=ressums(ir,1,1)+accs(i)
         if(rindex(ir).ne.0)then
            if(atomtype(i).eq.3)then
               ressums(ir,5,1)=ressums(ir,5,1)+accs(i)
               tsums(5)=tsums(5)+accs(i)
            elseif(atomtype(i).ne.0)then
               ressums(ir,4,1)=ressums(ir,4,1)+accs(i)
               tsums(4)=tsums(4)+accs(i)
               if(atomtype(i).eq.1)then
                  ressums(ir,2,1)=ressums(ir,2,1)+accs(i)
                  tsums(2)=tsums(2)+accs(i)
               elseif(atomtype(i).eq.2)then
                  ressums(ir,3,1)=ressums(ir,3,1)+accs(i)
                  tsums(3)=tsums(3)+accs(i)
               endif
            endif
         endif
      enddo
c
c     -- calculate realtive accessibilities
c
      do i = 1, resindex(nats)
         ires=rindex(i)
         if(stand.and.ires.ne.0)then
            do j = 1, 5
               if(standarea(ires,j).gt.0.0)then
                  ressums(i,j,2)=100.0*ressums(i,j,1)/standarea(ires,j)
               else
                  ressums(i,j,2)=0.0
               endif
            enddo
         else
            do j = 1, 5
               ressums(i,j,2)=-99.9
            enddo               
         endif
      enddo
c     
c     -- COMPLETED
c
      
      return
      end
c     
c     --what atom is it ? ie non-polar/polar, main/sidechain.
c	
      integer function what_atom(atom)
      character*4 atom, phobs(20), phils(20), mc(4)
c
c	atoms classed as non-polar in sidechains
c
      DATA PHOBS/' CA ',' CB ',' CD ',' CD1',' CD2',' CG ', 
     -     ' CG1',' CG2',' CE ',' CE1',' CE2',' CE3',
     -     ' CH2',' CZ ',' CZ2',' CZ3',' SD ',' SG ',
     -     '****','****'/
C     
C     ATOMS CLASSED AS POLAR IN SIDECHAINS
C
      DATA PHILS/' AD1',' AD2',' AE1',' AE2',' ND1',' ND2',
     -     ' NE ',' NE1',' NE2',' NH1',' NH2',' NZ ',
     -     ' OD1',' OD2',' OE3',' OE2',' OE1',' OG ',
     -     ' OG1',' OH '/
C     
C       MAIN CHAIN ATOMS
C     
      DATA MC/' N  ',' C  ',' O  ',' OXT'/  
      
      what_atom=3
      do i = 1, 4
         if(atom.eq.mc(i))return
      enddo
      
      what_atom=1
      do i = 1, 20
         if(atom.eq.phobs(i))return
      enddo
      
      what_atom=2
      do i = 1, 20
         if(atom.eq.phils(i))return
      enddo
      
      what_atom=0
      if(atom(2:2).eq.'O'.or.atom(2:2).eq.'N')then
         what_atom=2
      else
         what_atom=3
      endif
      
      return
      end
      
      subroutine vanin ( 
     -     vname,
     -     vlen,
     -     nacids, 
     -     aacids, 
     -     anames, 
     -     numats, 
     -     vradii,
     -     maxr,
     -     maxa
     -     )
c     
c -- Read in van der Waal radii from external file "vfile"
c -- nacids = number of residues read in
c -- aacids = *3 character array containing amino acid names
c -- anames = *4 character array containing atom names for each residue
c -- numats = array containing number of atoms for each residue
c -- radii  = vdw radii for each atom, indexed identically to anames
c
      integer       readstring, parse, fopen, addlab
      real          readfloat

      integer       nacids, 
     -              numats(maxr),
     -              l(256),
     -              nlabs,
     -              vlen
      real          vradii(maxr,maxa)
      character*3   aacids(maxr)
      character*4   anames(maxr,maxa)
      character*256 card, c(256), vname
      logical       ok

      if( fopen(1,vname,vlen,'old').eq.0 )then
         STOP'ERROR: unable to open "vdw.radii"'
      endif

      nacids = 0
      nlabs = 0

      do while ( readstring(1,card,len).ge.0 )
         n = parse(card,len,' ',c,l)
         if( c(1).eq.'RESIDUE' )then
            nacids = nacids + 1
            if( nacids.gt.maxr )STOP'ERROR: increase max_r'
            aacids(nacids) = c(3)(1:3)
            numats(nacids) = 0
         endif
         if(c(1).eq.'ATOM')then
            numats(nacids) = numats(nacids) + 1
            if( numats(nacids).gt.maxa )STOP'ERROR: increase max_a'
            anames(nacids,numats(nacids)) = card(6:9)
            vradii(nacids,numats(nacids)) = readfloat(c(n),l(n))
         endif
      enddo
      close(1)
      return
      end


