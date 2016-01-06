FC = gfortran
EXE = main
default:
	$(FC) -c m_definitions.f90
	$(FC) -c m_tools.f90
	$(FC) -c lebedev_quadrature.f90
	$(FC) -c m_atoms.f90	
	$(FC) -c m_shells.f90
	$(FC) -c m_basis.f90
	$(FC) -c m_rysquad.f90
	$(FC) -c m_rysroot.f90
	$(FC) -c m_eri.f90
	$(FC) -c m_hamilton.f90
	$(FC) -c m_dft_grid.f90
	$(FC) -c m_dft_xc.f90
	$(FC) m_definitions.f90 m_tools.f90 lebedev_quadrature.f90 m_atoms.f90 m_shells.f90 m_basis.f90 m_rysquad.f90 m_rysroot.f90 m_eri.f90 m_hamilton.f90 m_dft_grid.f90 m_dft_xc.f90 main.f90 -o $(EXE)
clean:
	rm -rf *.o *.mod $(EXE)




