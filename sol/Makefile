PROGRAM = ../orcwa.exe

OMPOPT  = /openmp

OBJS    =           \
Main.obj            \
input_data.obj      \
input2.obj          \
monitor.obj         \
memallocfree1.obj   \
memallocfree2.obj   \
memallocfree3.obj   \
setup.obj           \
setup_vector.obj    \
setupId.obj         \
setupSize.obj       \
setupDispersion.obj \
setupDft.obj        \
setupMurHx.obj      \
setupMurHy.obj      \
setupMurHz.obj      \
setupPmlEx.obj      \
setupPmlEy.obj      \
setupPmlEz.obj      \
setupPmlHx.obj      \
setupPmlHy.obj      \
setupPmlHz.obj      \
setupPml.obj        \
solve.obj           \
initfield.obj       \
updateEx.obj        \
updateEy.obj        \
updateEz.obj        \
updateHx.obj        \
updateHy.obj        \
updateHz.obj        \
dispersionEx.obj    \
dispersionEy.obj    \
dispersionEz.obj    \
vfeed.obj           \
efeed.obj           \
eload.obj           \
vpoint.obj          \
murH.obj            \
pmlEx.obj           \
pmlEy.obj           \
pmlEz.obj           \
pmlHx.obj           \
pmlHy.obj           \
pmlHz.obj           \
pbcx.obj            \
pbcy.obj            \
pbcz.obj            \
dftNear3d.obj       \
nearfield_c.obj     \
average.obj         \
outputChars.obj     \
outputZin.obj       \
outputSpara.obj     \
outputCoupling.obj  \
outputCross.obj     \
ingeometry.obj      \
geomlines.obj       \
farfield.obj        \
writeout.obj        \
utils.obj           \
cputime.obj         \
plot3dGeom.obj      \
ev3d.obj

CC      = cl.exe
CFLAGS  = /O2 /nologo /I../include /W4 /wd4996
CFLAGS  = $(CFLAGS) $(OMPOPT)
#CFLAGS  = $(CFLAGS) /D_DOUBLE
#CFLAGS  = $(CFLAGS) /D_ID32
LDFLAGS = /O2 /nologo
LIBS    =

$(PROGRAM): $(OBJS)
	@echo "Loading $(PROGRAM) ... "
	@$(CC) $(LDFLAGS) /Fe$(PROGRAM) $(OBJS) $(LIBS)

.c.obj:
	$(CC) /c $(CFLAGS) $<

clean:
	del *.obj
