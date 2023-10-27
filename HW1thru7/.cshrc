###############################################################
#
#         .cshrc file
#
#################################################################

# echo "entering .cshrc"

limit coredumpsize 0
umask 022             # prevent write access for group and world.

# put this here for parallel (still needed?)
setenv OvertureGridDirectories "/home/henshw/Overture.g/sampleGrids:/home/henshw/grids"

# setenv MPI_ROOT /usr/local/mpich2
# setenv MPI_ROOT /home/henshw/software/mpich-3.2-install
# Ubuntu -- mpiexec --version --> 3.3a2 
# setenv MPI_ROOT /usr
setenv MPI_ROOT /home/henshw/software/mpich-3.3.1-install

set path=(. $HOME/bin  $MPI_ROOT/bin  $path)

#         skip remaining setup if not an interactive shell
if ($?USER == 0 || $?prompt == 0) exit

#  ------  settings  for interactive shells -------

set history=50
set ignoreeof   # must type "exit" or "logout" to get out of a shell
set savehist=20
# set prompt="`hostname`{$LOGNAME}\!: "
set notify

# emacs needs libotf.so.0 which is found here:
## setenv LD_LIBRARY_PATH /usr/lib64/compat-openmpi-psm/lib

alias cd            'cd \!*;echo $cwd'
alias pwd           'echo $cwd'

# locale settings
setenv LANGUAGE "en_US.UTF-8"
setenv LC_ALL "en_US.UTF-8"
setenv LANG "en_US.UTF-8"

alias ls 'ls -F'
# alias lt 'ls -lth'
alias lm 'ls -lt | more'

# alias emacs /usr/bin/emacs -l /home/henshw/.emacs -geometry 120x70
# alias emacs /usr/bin/emacs -l /home/henshw/.emacs -geometry 120x80
# pipe errors into /dev/null from running remotely from a Mac: 
# alias emacs '/usr/bin/emacs -l /home/henshw/.emacs -geometry 120x80 >& /dev/null'
# alias emacs '/usr/bin/emacs -l /home/henshw/.emacs -font 9x15 -geometry 120x80 >& /dev/null'
## alias emacs '/usr/bin/emacs -l /home/henshw/.emacs -geometry 120x80+1200+50 >& /dev/null'

# alias emacs '/usr/bin/emacs -l /home/henshw/.emacs -geometry 120x70+1200+50 >& /dev/null'
# for class: 
alias emacs '/usr/bin/emacs -l /home/henshw/.emacs -geometry 100x70+100+50 >& /dev/null'

# higher emacs window to right: 
alias emacsr '/usr/bin/emacs -l /home/henshw/.emacs -geometry 120x100+1200+50 >& /dev/null'
# higher emacs window to left
# alias emacsl '/usr/bin/emacs -l /home/henshw/.emacs  --font="Monospace-8" -geometry 120x80+100+50 >& /dev/null'
alias emacsl '/usr/bin/emacs -l /home/henshw/.emacs  --font="Monospace-10" -geometry 120x80+100+50 >& /dev/null'
# use small font for remote bug 
## alias emacss /usr/bin/emacs -l /home/henshw/.emacs --font="DejaVu Sans Mono-8" -geometry 120x80
## alias emacss /usr/bin/emacs -l /home/henshw/.emacs --font="Monospace-11" -geometry 130x80+1200+50
alias emacss /usr/bin/emacs -l /home/henshw/.emacs -geometry 130x80+1200+50

setenv VENDOR intel

# for doxygen: 
setenv DOXYGEN_VERSION "Version 25"

setenv CG_DOXYGEN "/home/henshaw.0/cg.v25" 
setenv CGDOC_DOXYGEN "/home/henshaw.0/cgDoc/doxygen" 

setenv OV_DOXYGEN "/home/henshaw.0/Overture.v25.d" 
setenv OVDOC_DOXYGEN "/home/henshaw.0/Overture/doxygen" 

# setenv TEXINPUTS ".:~henshaw/latex2e::/home/henshaw/latex2e/seminar/inputs:/usr/local/latex2html/texinputs"
setenv TEXINPUTS ".::$HOME/latex2e:$HOME/latex2e/pgfplots//"
# setenv TEXINPUTS ".::$HOME/latex2e"

alias module 'set args = (\!*); source $HOME/bin/module'

alias df 'df -k'
alias ps2ppm $HOME/Overture.g/bin/ps2ppm

set ovg = $HOME/Overture.g
set ovp = $HOME/Overture.p

set ofs = /data2/henshw/overtureFramework
set ovs = /data2/henshw/Overture.s
set cgs = /data2/henshw/cg.s

setenv OVHOME $HOME/overtureFramework
set of = $HOME/overtureFramework

set ov = $OVHOME/Overture
set cg = $OVHOME/cg
set mapping = $OVHOME/Overture/mapping 
set gf     = $OVHOME/Overture/gf
set grid   = $OVHOME/Overture/grid
set primer = $OVHOME/Overture/primer
set ogshow = $OVHOME/Overture/ogshow
set ogmg   = $OVHOME/Overture/ogmg
set ogmgp  = $OVHOME/Overture/ogmg/ogmgt

set ogmgs  = $OVHOME/Overture/ogmgs
set ogen   = $OVHOME/Overture/ogen
set ogenp  = $OVHOME/Overture/ogenp
set ogens  = $OVHOME/Overture/ogens
set tosx   = $OVHOME/Overture/oges/results/tos
set tcm3x   = $OVHOME/Overture/op/tests/tcm3
set ugen   = $OVHOME/Overture/ugen
set oges   = $OVHOME/Overture/oges
set hype   = $OVHOME/Overture/hype
set rap    = $OVHOME/Overture/rap
set op     = $OVHOME/Overture/op
set os     = $OVHOME/Overture/otherStuff
set gui    = $OVHOME/Overture/gui
set parallel = $OVHOME/Overture/parallel


alias hype $OVHOME/Overture/hype/hype
set hypex = $OVHOME/Overture/hype/hype

set articles = $HOME/articles
set papers = $HOME/papers

set memo   = $HOME/memo
set admin  = $HOME/admin
set talks  = $HOME/talks
# set book   = $res/book
set mgp = $papers/mgp
set runs   = $HOME/runs

set db   = $HOME/runs/mx/dielectricBodies

set cgDoc = $OVHOME/cgDoc
set insDoc = $OVHOME/cgDoc/ins
set cnsDoc = $OVHOME/cgDoc/cns
set mxDoc = $OVHOME/cgDoc/mx
set smDoc = $OVHOME/cgDoc/sm
set mpDoc = $OVHOME/cgDoc/mp

set inses = $HOME/codes/inses

setenv CGWAVE $HOME/Dropbox/research/cgwave
set wave = $CGWAVE
set cgwhd = $CGWAVE
set cgwh = $CGWAVE/bin/cgwh
alias cgWave $CGWAVE/bin/cgWave
set cgWave = $CGWAVE/bin/cgWave
set cgWavex = $CGWAVE/bin/cgWave

alias cgwh $CGWAVE/bin/cgwh
set cgwh = $CGWAVE/bin/cgwh

set eig = /home/henshw/Dropbox/research/eig
alias genEigs $eig/bin/genEigs
set genEigs = $eig/bin/genEigs

alias eveSolver $eig/bin/eveSolver
set eveSolver = $eig/eveSolver

alias genEigsILE $eig/bin/genEigsILE
set genEigsILE = $eig/bin/genEigsILE

# new marching grid generator:
setenv STRIDER $HOME/Dropbox/research/strider
set stride = $STRIDER
alias strider $STRIDER/bin/strider
set strider = $STRIDER/bin/strider

# Test Tom's radiation BC's 
alias trad /home/henshw/cg.g/mx/bin/trad
set tradx = /home/henshw/cg.g/mx/bin/trad

alias tradp /home/henshw/cg.p/mx/bin/trad
set tradpx = /home/henshw/cg.p/mx/bin/trad

set mxsosup = $HOME/papers/mxsosup
# set fibr = $HOME/papers/fibr
# set fibr = $HOME/Dropbox/CG/fibr

set emh = "$HOME/Dropbox/EM_Homogenization"
set bah = "$HOME/Dropbox/EM_Homogenization/notes"

# set se = $HOME/runs/mp/shockEllipse

# set ib = $HOME/runs/mx/interfaceBump

set dmx = ~/Dropbox/CG/DMX/DMX_ADE_notes
set ssmx = ~/Dropbox/AMP/ssmx
# set ssmx = ~/Dropbox/CG/ssmx
set fins = ~/Dropbox/CG/fins

set adei = ~/Dropbox/DARPA/RPI/adePapers/adegdmi

# Maxwell-Bloch paper
set mbe = ~/Dropbox/DARPA/RPI/Maxwell-Bloch/papers/mbe

set ism = ~/Dropbox/AMP/ism
set waho = ~/Dropbox/AMP/waho
set champ4 = ~/Dropbox/AMP/champ4


set r2p = ~/Dropbox/r2P

# papers: (working copy on /data1,  bare repo on /data2) symbolic links to :
set papers = ~/papers
set talks = ~/talks

set beamins = ~/papers/beamins

set flunsi = $HOME/tux/papers/flunsi

set ad     = $cg/ad
set ins    = $cg/ins
set cns    = $cg/cns
set asf    = $cg/asf
set sm     = $cg/sm
set mx     = $cg/mx
set mp     = $cg/mp

set tcmc = $op/tests/tcmConstraint

alias trdmf /home/henshw/cg.g/mx/bin/trdmf
alias texact /home/henshw/cg.g/mx/bin/texact

set cgb = ~/Dropbox/CG6backup

setenv BAMX $HOME/Dropbox/research/bamxSolver
set bamxd = $HOME/Dropbox/research/bamxSolver
alias bamx "$HOME/Dropbox/research/bamxSolver/bin/bamx"
set bamxx = $HOME/Dropbox/research/bamxSolver/bin/bamx
# paper:
set bamxp = $HOME/Dropbox/DARPA/BAMX/papers/bamx
set bagdmp = $HOME/Dropbox/DARPA/BAMX/papers/bagdm

# LLNL:
alias s2 "ssh -X henshaw1@tux291"
alias sc "ssh -X henshaw@cab"

# Don's machine: Michael S. Schwendeman
# alias sm "ssh -Y henshaw@mss.math.rpi.edu"

# Don's machine Claudia L. S.
alias sm "ssh -Y henshaw@cls.math.rpi.edu"

alias srtg1 "ssh -Y henshw@rtg1.math.rpi.edu"

alias s3 "ssh -Y henshw@cg3.math.rpi.edu"

# web pages are here:
alias sr "ssh -X henshw@rcs.rpi.edu"

# CCI Blue Gene Ft...
alias sc "ssh -Y OVTRhnsh@lp01.ccni.rpi.edu"
alias sfc "sftp OVTRhnsh@lp01.ccni.rpi.edu"

alias sbg "ssh -Y OVTRhnsh@lp01.ccni.rpi.edu"
alias sfbg "sftp OVTRhnsh@lp01.ccni.rpi.edu"

alias svnc 'svn commit -m"wdh"'
alias svnr 'svn resolved'

alias gita 'git add -u -v'
alias gitc 'git commit -m"wdh"'
alias gits 'git status --untracked-files=no'
alias gitpush 'git push origin master'
# This can cause trouble with loosing commits: 
# alias gitpull 'git pull --rebase'
# Do this instead from now on - June 26, 2021
alias gitpull 'git pull'

# Use 'git fetch' to update local copies of the remote repo on sourceforge -- does not change branches
# See all differences with master: but first do 'git fetch'
alias gitdm 'git diff --name-only origin/master'

alias kpdf 'okular'

# --track-origins=yes : to see where uninitialized values come from
alias valgrindsuppress "/usr/bin/valgrind --gen-suppressions=all --tool=memcheck --suppressions=/home/henshw/valgrind.supp --num-callers=20 --error-limit=no --db-attach=yes --max-stackframe=3146824 --track-origins=yes"

# try this 
alias valgrindebug "/usr/bin/valgrind --tool=memcheck --suppressions=/home/henshw/valgrind.supp --num-callers=20 --error-limit=no --max-stackframe=3146824"

alias valgrind "/usr/bin/valgrind --tool=memcheck --suppressions=/home/henshw/valgrind.supp --num-callers=20 --error-limit=no --max-stackframe=3146824"
set valgrind = "/usr/bin/valgrind --tool=memcheck --suppressions=/home/henshw/valgrind.supp --num-callers=20 --error-limit=no --max-stackframe=3146824"


# old: 
# valgrindebug "/usr/bin/valgrind --tool=memcheck --suppressions=/home/henshw/valgrind.supp --num-callers=20 --error-limit=no --db-attach=yes  --max-stackframe=3146824"

# alias valgrindebug "/usr/bin/valgrind --tool=memcheck --suppressions=/home/henshw/valgrind.supp --num-callers=20 --error-limit=no --vgdb-error=1  --max-stackframe=3146824"


# alias valgrindsuppress "/usr/bin/valgrind --gen-suppressions=yes --tool=memcheck --suppressions=/home/henshw/valgrind.supp --num-callers=20 --error-limit=no --vgdb-error=1  --max-stackframe=3146824"

module g

# dropbox start

