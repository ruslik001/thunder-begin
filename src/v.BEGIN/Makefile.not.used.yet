BEGIN = M_psi.o M_rcatms_Coulomb.o M_rcatms_Coulomb_DMOL_excited.o M_rcatms_DMOL_excited.o M_rcatms.o orthogonalize_psi.o M_vnn.o
BEGIN_LOOPRC = M_psi.o M_rcatms_Coulomb.o M_rcatms.o M_looprc.o


all: libthd_BEGIN.a libthd_BEGIN_LOOPRC.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.f90 *~

libthd_BEGIN.a: $(BEGIN)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

libthd_BEGIN_LOOPRC.a: $(BEGIN_LOOPRC)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
