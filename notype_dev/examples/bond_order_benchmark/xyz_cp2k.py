# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under
# the GNU General Public License.

# xyz tool

oneline = "Convert LAMMPS snapshots to XYZ format"

docstr = """
x = xyz(d)		d = object containing atom coords (dump, data)

x.one()                 write all snapshots to tmp.xyz
x.one("new")            write all snapshots to new.xyz
x.many()                write snapshots to tmp0000.xyz, tmp0001.xyz, etc
x.many("new")           write snapshots to new0000.xyz, new0001.xyz, etc
x.single(N)             write snapshot for timestep N to tmp.xyz
x.single(N,"file")      write snapshot for timestep N to file.xyz
"""

# History
#   8/05, Steve Plimpton (SNL): original version

# ToDo list

# Variables
#   data = data file to read from

# Imports and external programs

import sys

# Class definition

class xyz:

  # --------------------------------------------------------------------

  def __init__(self,data):
    self.data = data

  # --------------------------------------------------------------------

  def one(self,*args):
    if len(args) == 0: file = "tmp.xyz"
    elif args[0][-4:] == ".xyz": file = args[0]
    else: file = args[0] + ".xyz"

    f = open(file,"w")
    n = flag = 0
    while 1:
      which,time,flag = self.data.iterator(flag)
      if flag == -1: break
      time,box,atoms,bonds,tris,lines = self.data.viz(which)

      f.write( "%d\n"%len(atoms) )
      f.write( f"i=\t{time}\n" )
      for atom in atoms:
        itype = int(atom[1])
        f.write(f"{itype} {atom[2]} {atom[3]} {atom[4]}\n")

      print(time),
      sys.stdout.flush()
      n += 1

    f.close()
    print("\nwrote %d snapshots to %s in XYZ format" % (n,file) )
    print("\nbox [xlo,ylo,zlo,xhi,yhi,zhi,xy,xz,yz] =",box)
    xlob, ylob, zlob, xhib, yhib, zhib, xy, xz, yz = box
    xlo = xlob - min(min(0.0,xy), min(xz,xy+xz) );
    xhi = xhib - max( max(0.0,xy), max(xz,xy+xz) );
    ylo = ylob - min(0.0,yz);
    yhi = yhib - max(0.0,yz);
    zlo = zlob
    zhi = zhib

    BOX = [xhi - xlo, xy, xz, yhi - ylo, yz, zhi - zlo]
    print("\nBOX [ax,bx,cx,by,cy,cz] =",BOX)
    return BOX

  # --------------------------------------------------------------------

  def many(self,*args):
    if len(args) == 0: root = "tmp"
    else: root = args[0]

    n = flag = 0
    while 1:
      which,time,flag = self.data.iterator(flag)
      if flag == -1: break
      time,box,atoms,bonds,tris,lines = self.data.viz(which)

      if n < 10:
        file = root + "000" + str(n)
      elif n < 100:
        file = root + "00" + str(n)
      elif n < 1000:
        file = root + "0" + str(n)
      else:
        file = root + str(n)
      file += ".xyz"
      f = open(file,"w")
      f.write( "%d\n"%len(atoms) )
      f.write( "Atoms\n" )
      for atom in atoms:
        itype = int(atom[1])
        f.write(f"{itype} {atom[2]} {atom[3]} {atom[4]}\n")
      print(time),
      sys.stdout.flush()
      f.close()
      n += 1

    print("\nwrote %s snapshots in XYZ format" % n)

  # --------------------------------------------------------------------

  def single(self,time,*args):
    if len(args) == 0: file = "tmp.xyz"
    elif args[0][-4:] == ".xyz": file = args[0]
    else: file = args[0] + ".xyz"

    which = self.data.findtime(time)
    time,box,atoms,bonds,tris,lines = self.data.viz(which)
    f = open(file,"w")
    f.write( "%d\n"%len(atoms) )
    f.write( "Atoms\n" )
    for atom in atoms:
      itype = int(atom[1])
      f.write(f"{itype} {atom[2]} {atom[3]} {atom[4]}\n")
    f.close()
