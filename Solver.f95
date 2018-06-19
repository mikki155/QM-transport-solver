!Welcome to the electron transport solver!

!In order to run this program, simply use the shell script named 'execute.sh'.

!This program currently has two options:
!1) Schrodinger-Poisson, and 2) NEGF. 
!The first solves eigenvectors and eigenenergies
!using Schrodinger, and then uses these in Poisson to solve for a new potential.
!This potential is then used to solve for new eigenvectors and eigenenergies,
!and the overall process is repeated until convergence.

!In 2), the program first solves the left contact self-energy, then the block diagonals,
!and then the right contact self-energy. Scroll further down in order to learn more,
!about how the NEGF process is calculated. The electron density that is calculated
!is inserted into the Poisson routine in order to gain a new potential.

!The program only calculates potentials and states for the AlGaAs/GaAs layered
!structure, but the user is welcomed to extend this to include other structures
!as well that can be used as input.

!The code is by no means optimized nor finished, and the user is welcomed to improve
!and complete it. One issue is with regards to the variables, a lot of them take
!values that might be completely wrong. More information on this will be given
!further down where it is relevant.

!After the user has read all of the comments in this code, the challenges which
!could or must be solved a summed up as follows (ranked from easiest to hardest):

!1. The dimensions of variables and states must be corrected.
!2. Correct expressions for the various self-energies should be implemented.
!3. Different phonon modes could/should be introduced.
!4. As of now, the self-energies do not converge, solving the problems above should correct this.
!5. It's thought that the algorithm calculating GFs blocks could be improved, although it may be very challenging.

!NEGF can be frustratingly confusing, so feel free to contact me regarding the code
!at mi_hau94@hotmail.com


	program Solver
		implicit none
		!The following parameters are respectivaly the dimension of the number of layers (Nldim), 
		!the thickness of the structure (Nthick), the convergence number for the Schrodinger-Poisson loop (conv1), 
		!the convergence number for the NEGF loop (conv2), energy levels (n), 
		!the resolution of the structure for NEGF calculations (Ngf),
		!and the number of iterations for NEGF-Schrodinger-Poisson (endloop).
		!Feel free to change these parameters, but keep in mind that Ngf should not be adjusted much higher
		!because this will drastically increase computational times.
		!Keep in mind that if Nldim/Nthick are adjusted, then "layers" below must be adjusted accordingly.
		integer, parameter		:: Nldim = 9, Nthick = 600, conv1 = 50, conv2 = 10, n = 5, &
						 Ngf = 100, endloop = 5
		!Most of the variables below are explained in the subroutines they appear.
		!Note however that percent(conv1), percentGF(endloop), and perit are used to display how
		!much of the process is completed percentage wise.
		integer			:: percent(conv1), percentGF(conv2), LDZ, INFO, INFOTWO, &
					iter, true, perit, input, calcnd, horn, loop, fileunit, respsi
		integer, dimension(Nldim)		:: layers
		double precision, dimension(max(1, 2*Nthick-2))	:: WORK
		double precision			:: Vgf(Ngf), nd(Ngf, n), Dpel, hw, delta, & 
							te, gma, md, svel, Da, Dpinel, Vtgf(Ngf), V0gf(Ngf), &
							psi(Nthick, Nthick), psigf(Ngf, Ngf), hbar, eps0, me, &
							 kB, E0, el0, meff, Vb, x, T, pi, rhod, mu, interm, en, pm
		double precision, dimension(Nthick)		:: V0, Vt, V, D, E
		double complex			:: gzz(Ngf, Ngf), SigC(Ngf, Ngf), grl(Ngf, Ngf, Ngf), Gretl(Ngf, Ngf, Ngf), &
						 SigCr(Ngf, Ngf), gal(Ngf, Ngf, Ngf), Sigma(Ngf, Ngf, Ngf), &
						gnl(Ngf, Ngf, Ngf), Gecordiag(Ngf, Ngf, Ngf), gpl(Ngf, Ngf, Ngf), &
						Gecorminhw(Ngf, Ngf, Ngf), Gecormaxhw(Ngf, Ngf, Ngf), Gpdiag(Ngf, Ngf, Ngf), &
						 Sigin(Ngf, Ngf, Ngf), Sigout(Ngf, Ngf, Ngf), Sigelret(Ngf, Ngf, Ngf), &
						 Gpmaxhw(Ngf, Ngf, Ngf), Gpminhw(Ngf, Ngf, Ngf), &
						 SiginCr(Ngf, Ngf), SigoutC(Ngf, Ngf), SigoutCr(Ngf, Ngf), SiginC(Ngf, Ngf)
		psi = 0
		psigf = 0
		true = 0
		interm = 0
		loop = 0
		percent = 0
		mu = -8.6E-21 !Assumed constant throughout the structure. This value may be fixed 
		el0 = 1.60217646E-19 !Elementary charge
		E0 = 1E6 !Applied bias
		md = 5.32 !Mass density, used to find phonon self-energy
		svel = 5.29 !Sound speed
		hw = 36E-20 !LO phonon energy
		x = 15E-2 !Mole fraction
		Vb = el0*0.62*(1.594*x + x*(1-x)*(0.127-1.31*x)) !Conduction band offset factor
		Vt = 0 !The new potential which is found in Poisson.
		T = 200 !Temperature in Kelvin
		rhod = 1.9E22 !Doping density
		pi = 3.1415926535897932
		hbar = 1.05457168E-34 !Planck constant
		me = 9.10938188E-31 !Electron mass
		kB = 1.3806503E-23 !Boltzmann constant
		Da = 10 !Acoustic deformation pot (eV)
		gma = 6.545E9 !Wave number in contact
		meff = 0.067*me
		eps0 = 8.8541878176E-12
		layers = (/ 57, 80, 25, 65, 41, 155, 30, 90, 57 /) !Each layer of GaAs/AlGaAs in Angstrom.
		delta = 1E-10
		te = (hbar*hbar)/(2*0.01*me*delta*delta)
		Dpel = (Da*Da*kB*T)/(2*md*svel*svel) !Elastic deformation factor
		Dpinel = (Da*Da*hbar*hbar)/(2*md*hw) !Feel free to adjust these
		Dpinel = Dpinel*6E10 !Correct for order
		Sigma = 0 !Total self-energy
		Sigin = 0
		Sigout = 0
		Sigelret = 0 !Elastic self-energy
		SiginC = 0
		SigoutC = 0
		SiginCr = 0
		SigoutCr = 0
		calcnd = 0
		fileunit = 1
		horn = 0 !If hole or electron inelastic self-en. is to be calculated
		pm = 0 !Used in gzero to indicate sign of energy
		print*,"**********************************************"
		print*,"*                                            *"
		print*,"* Welcome to the electron transport solver   *"
		print*,"* for layered semiconductor devices!         *"
		print*,"*                                            *"
		print*,"**********************************************"
		print*,"                                	      "
		print*,"                                	      "
		print*,"This program solves energies, states and densities for the GaAs/AlGaAs layered structure."
		print*,"											 "
		print*,"											 "
		print*,"For preliminary results, choose Schrödinger-Poisson. To include scattering events and photons, choose NEGF:"
		print*," "
		print*,"1) Schrödinger-Poisson"
		print*,"2) NEGF"
		read(*,*) input
		call initprofile(Nthick, Nldim, layers, V0, Vb, E0, el0) !Calculates the potential profile of the system given the applied bias.
		!The code below plots first the initial profile, V0. It then updates the fileunit which is used
		!to create a new file (see 'plotprofile' comment furhter down). Then, the total potential V is
		!updated where Vt is the potential involved in the Poisson calculations.
		!It then iterates between Schrodinger and Poisson, and gives out how far the program has come in percentage.
		!Finally, it plots each of the states for the n levels chosen.
		if (input.EQ.1) then
			print*,"Initiating Schrödinger-Poisson calculations..."
			call plotprofile(V0, Nthick, fileunit)
			fileunit = fileunit + 1
			V = V0 + Vt
			do iter = 1, conv1
				call schrodinger(Nthick, D, E, V, LDZ, WORK, INFO, hbar, meff, psi)
				call poisson(mu, Nthick, n, INFOTWO, input, meff, pi, hbar, el0, kb, rhod, eps0, T, psi, D, Vt, nd)
				if (mod(iter, 10).EQ.0) then
					interm = real(iter)/real(conv1)*100
					interm = int(interm)
					percent(iter) = interm
					if (percent(iter).NE.0) then
						do perit = 1, conv1
							if ((percent(perit).EQ.percent(iter)).AND.(perit.NE.iter)) then
								true = 1
							end if
						end do
						if (true.NE.1) then
							print*, percent(iter),"% complete."
						end if
					end if
				true = 0
				end if
				V = V0 + Vt
			end do
			psi = psi/1E5
			psi = abs(psi)**2
			print*, "Self-consistent solution converged!"
			do iter = 1, n
				call plotprofile(psi(:, iter), Nthick, fileunit)
				fileunit = fileunit + 1
			end do
			call plotprofile(V, Nthick, fileunit)
			fileunit = fileunit + 1
			call plotprofile(D, Nthick, fileunit)

		!This is the NEGF part, which may be a bit confusing. As in 1), the total potential V is first calculated.
		!Then, it first calculates g_{0, 0}, resizes the potential, and starts to calculate the GFs for the blocks.
		!Now, the inelastic self-energy described further down in 'inelselfen' depends on the electron correlation
		!function at two different energies; E - hw and E + hw. So first, all the GFs at energy E - hw must be
		!calculated, then the total self-energy Sigma must be reset so that all the GFs may be calculated again
		!for the energy E + hw. This is the reason why pm is included; it ensures that the correct energy in
		!'gzero' is calculated. 

		!After the new Sigin/out have been calculated using G^n/G^p, the program goes
		!into a new loop. In this loop, new GFs are calculated using the new self-energies, and at the end
		!of the loop new self-energies Sigin/out are calculated from these GFs. This process repeats until
		!some convergence is reached (conv2). After the loop is done, the energy is reset to the original
		!energy, E. The corresponding GFs are calculated, along with the electron density. 

		!This density
		!is then inserted into the Poisson loop further down. The new potential is calculated, and the entire
		!process of calculating the electron density repeats using now the new potential. Note that it is not
		!necessary to calculate schrodinger in this case, but it can be nice to see the corresponding
		!states and compare them to 1). Once the integer loop has reached endloop,
		!the NEGF calculations are done. The new profiles are finally plotted,
		!and option 2) ends.

		!A challenge for the user is that 'elselfen' is only included once, while this must be included
		!for each of the energy cases. Another challenge is that the program only calculates the energy
		!for a single phonon mode, and probably wrongfully so. A more correct implementation is then
		!welcomed.

		else if (input.EQ.2) then
			print*,"Initiating NEGF calculations..."
			V = V0 + Vt
			do loop = 1, endloop
				call gzero(gzz, SigC, Ngf, hbar, me, pm)
				if (loop.EQ.1) then
					call negfPot(V, Vgf, Nthick, Ngf)
					call negfPot(V0, V0gf, Nthick, Ngf)
				end if
				en = te + 0.1*Vgf(1) - 2*te*cos(gma*delta) !Assuming positive potential, -hw for inelastic self-en.
				pm = 0
				call gblocks(gzz, grl, gal, gnl, gpl, SigC, SigCr, Sigma, Sigin, &
						 Sigout, Ngf, meff, me, hbar, en, pm, Vgf)
				call continselfen(SigC, SigCr, Vgf, Ngf, mu, kB, T, SiginC, SiginCr, SigoutC, SigoutCr)
				do iter = 1, Ngf !Initial guess for Sigin/out
					Sigin(iter, :, :) = SiginC(:, :) + SiginCr(:, :)
					Sigout(iter, :, :) = SigoutCr(:, :) + SigoutC(:, :)
				end do
				call GFs(grl, gal, gnl, gpl, Gretl, Gecordiag, Gpdiag, Ngf, hbar, meff, pi, n, nd, calcnd)
				call elselfen(Gretl, Sigelret, Dpel, Ngf) !Calculates elastic self-energy
				en = en - hw !Starts process of calculating inelastic self-en.
				pm = -hw
				Sigma = 0
				Sigma = Sigelret
				call gzero(gzz, SigC, Ngf, hbar, me, pm)
				call gblocks(gzz, grl, gal, gnl, gpl, SigC, SigCr, Sigma, Sigin, &
						 Sigout, Ngf, meff, me, hbar, en, pm, Vgf)
				call GFs(grl, gal, gnl, gpl, Gretl, Gecorminhw, Gpminhw, Ngf, hbar, meff, pi, n, nd, calcnd)
				en = en + 2*hw
				pm = hw
				Sigma = 0
				Sigma = Sigelret
				call gzero(gzz, SigC, Ngf, hbar, me, pm)
				call gblocks(gzz, grl, gal, gnl, gpl, SigC, SigCr, Sigma, Sigin, &
						 Sigout, Ngf, meff, me, hbar, en, pm, Vgf)
				call GFs(grl, gal, gnl, gpl, Gretl, Gecormaxhw, Gpmaxhw, Ngf, hbar, meff, pi, n, nd, calcnd)
				call inelselfen(Gecorminhw, Gecormaxhw, Sigin, Dpinel, hw, kB, T, Ngf, horn) !Inelastic in-scattering is calculated
				horn = 1
				call inelselfen(Gpminhw, Gpmaxhw, Sigout, Dpinel, hw, kB, T, Ngf, horn) !Inelastic out-scattering is calculated
				horn = 0
				Sigma = 0
				Sigma = Sigin + Sigout + Sigelret
				do iter = 1, conv2
					en = en - 2*hw
					pm = - hw
					call gzero(gzz, SigC, Ngf, hbar, me, pm)
					call gblocks(gzz, grl, gal, gnl, gpl, SigC, SigCr, Sigma, Sigin, &
							 Sigout, Ngf, meff, me, hbar, en, pm, Vgf)
					call GFs(grl, gal, gnl, gpl, Gretl, Gecorminhw, Gpminhw, Ngf, hbar, meff, pi, n, nd, calcnd)
					en = en + 2*hw
					pm = hw
					Sigma(1, :, :) = Sigma(1, :, :) - SigC !To avoid counting contact self-energies twice
					Sigma(Ngf, :, :) = Sigma(Ngf, :, :) - SigCr
					call gzero(gzz, SigC, Ngf, hbar, me, pm)
					call gblocks(gzz, grl, gal, gnl, gpl, SigC, SigCr, Sigma, Sigin, &
							 Sigout, Ngf, meff, me, hbar, en, pm, Vgf)
					call GFs(grl, gal, gnl, gpl, Gretl, Gecormaxhw, Gpmaxhw, Ngf, hbar, meff, pi, n, nd, calcnd)
					call inelselfen(Gecorminhw, Gecormaxhw, Sigin, Dpinel, hw, kB, T, Ngf, horn)
					horn = 1
					call inelselfen(Gpminhw, Gpmaxhw, Sigout, Dpinel, hw, kB, T, Ngf, horn)
					horn = 0
					Sigma = 0
					Sigma = Sigin + Sigout + Sigelret
					interm = (real(iter)/real(conv2))*100
					percentGF(iter) = interm !Need this vector to avoid fractional part of percentage
				print*, percentGF(iter), "% complete at loop", loop
				end do
				en = en - hw
				pm = 0
				Sigma(1, :, :) = Sigma(1, :, :) - SigC
				Sigma(Ngf, :, :) = Sigma(Ngf, :, :) - SigCr
				call gzero(gzz, SigC, Ngf, hbar, me, pm)
				call gblocks(gzz, grl, gal, gnl, gpl, SigC, SigCr, Sigma, Sigin, &
						Sigout, Ngf, meff, me, hbar, en, pm, Vgf)
				calcnd = 1
				call GFs(grl, gal, gnl, gpl, Gretl, Gecordiag, Gpdiag, Ngf, hbar, meff, pi, n, nd, calcnd)
				do iter = 1, conv1
					call schrodinger(Nthick, D, E, V, LDZ, WORK, INFO, hbar, meff, psi)
!					call schrodinger(Ngf, D, E, Vgf, LDZ, WORK, INFO, hbar, meff, psigf)
					do respsi = 1, n
						call negfPot(psi(:, respsi), psigf(:, respsi), Nthick, Ngf)
					end do
					call poisson(mu, Ngf, n, INFOTWO, input, meff, pi, hbar, el0, kb, &
						 rhod, eps0, T, psi, D, Vtgf, nd)
					Vgf = V0gf + Vtgf
					call revPot(V, Vgf, Nthick, Ngf)
				end do
			end do
			psigf = psigf/1E5
			psigf = abs(psigf)**2
			do iter = 1, n
				call plotprofile(psigf(:, iter), Ngf, fileunit)
				fileunit = fileunit + 1
			end do
			call plotprofile(Vgf, Ngf, fileunit)
			fileunit = fileunit + 1
			do iter = 1, n
				call plotprofile(nd(:, iter), Ngf, fileunit)
				fileunit = fileunit + 1
			end do			
			print*, "Calculations complete."
		else
			print*,"Error: Not an option."
		end if
	end program Solver

!The subroutine 'initprofile' takes in the layered structure, and calculates the conduction band profile
!given a bias E0 and the offset factor Vb. The resulting profile V0 should be "bent" in a degree
!proportional to the applied bias. The higher E0 is, the more "bending" should be obtained.

	subroutine initprofile(N, Nldim, layers, V0, Vb, E0, el0)
		implicit none
		integer				:: i, k1, k2
		double precision, intent(in)		:: Vb, E0, el0
		integer, intent(in)			:: N, Nldim, layers(Nldim)
		double precision, dimension(N), intent(out)		:: V0
		double precision, dimension(N)				:: z
		V0 = 0
		z = 0
		k1 = 1
		k2 = layers(1)
		do i=1,(N-1)
			if (mod(k1, 2).NE.0) then
				V0(i) = Vb
			end if
			z(i) = el0*i*E0*1E-10
			V0(i) = V0(i) - z(i)
			if (i.EQ.k2) then
				k1 = k1 + 1
				k2 = k2 + layers(k1)
			end if
		end do
		if (mod(k1, 2).NE.0) then
			V0(N) = Vb
		end if
		z(N) = el0*N*E0*1E-10
		V0(N) = V0(N) - z(N)
	end subroutine initprofile

!The subroutine 'plotprofile' is a simple code which plots x-coordinates up to Nvdim
!together with the y-coordinates V. Note that 'fileunit' is an integer that allows for
!'plotprofile' to be used up to 100 times in 1 simulation. I do realize that Fortran
!allows you to reopen a file and continue to write to it, but this gives you one file
!for each profile you find in the simulation. I found this to be a bit more convenient,
!as I stored each of the profiles in seperate vectors in MATLAB in order to plot them.

	subroutine plotprofile(V, Nvdim, fileunit)
		implicit none
		integer, intent(in)	:: Nvdim
		double precision, dimension(Nvdim), intent(in)	:: V
		integer, dimension(Nvdim)		:: z
		integer		:: i, fileunit
		character(13)		:: filename
		character(7)		:: intermid1
		character(4)		:: intermid2
		character(2)		:: num
		intermid1 = 'profile'
		intermid2 = '.dat'
		z = 0
		do i=1,Nvdim
			z(i) = z(i) + i
		end do
		if ((fileunit.LT.0).OR.(fileunit.GT.99)) then
			print*, "Error: could not write to file."
		else
			if (fileunit.LT.10) then
				write(num, '(I1)') fileunit
				filename = trim(intermid1)//trim(num)
				filename = trim(filename)//trim(intermid2)			
				open(unit = fileunit, file=filename)
					do i=1,Nvdim
						write(fileunit,*) z(i), V(i)
					end do
				close(unit=fileunit)
			else
				write(num, '(I2)') fileunit
				filename = trim(intermid1)//trim(num)
				filename = trim(filename)//trim(intermid2)
				open(unit = fileunit, file=filename)
					do i=1,Nvdim
						write(fileunit,*) z(i), V(i)
					end do
				close(unit=fileunit)
			end if
		end if
	end subroutine plotprofile

!The subroutine 'eigs' takes in matrices D, E, and Z, and gives out eigenvalues in E
!and eigenvectors in Z. It uses the subroutine 'dstev' from the lapack library, so
!I would recommend looking this up if it is unknown to the user.

	subroutine eigs(N, D, E, Z, LDZ, WORK, INFO)
		implicit none
		integer, intent(in)		:: N
		double precision	:: D(N), E(N-1), Z(N, N), WORK(max(1, 2*N-2))
		integer				:: LDZ, INFO
		external dstev
		LDZ = N
		call DSTEV('V', N, D, E, Z, LDZ, WORK, INFO)
		if (INFO.NE.0) then
			print*,"Error: Could not calculate eigenvalues/vectors."
		end if
	end subroutine eigs

!This uses the FDM Schrodinger equation, and accordingly calculates the eigenvalues and eigenenergies.
!The spacing is delta, and may be adjusted as the user sees fit.
!Note that psi is multiplied by a large factor. This is because when psi is sent into 'eigs',
!the 'dstev' subroutine normalizes the eigenvector that it calculates. I have not found any
!appropriate solution to this, other than multiplying psi with 1E5 which seems to give back
!the correct dimensions. If this is not done, then the Schrodinger-Poisson loop gives
!completely wrong results.

	subroutine schrodinger(Nthick, D, E, V, LDZ, WORK, INFO, hbar, meff, psi)
		implicit none	
		integer, intent(in)		:: Nthick
		double precision, dimension(Nthick, Nthick)	:: psi
		double precision, intent(in)		:: hbar, meff
		double precision		:: delta, s, diag
		double precision, dimension(Nthick)		:: D(Nthick), V(Nthick), E(Nthick-1)
		double precision, dimension(max(1, 2*Nthick-2))	:: WORK
		integer		:: LDZ, INFO
		external dstev
		delta = 1E-10
		s = hbar*hbar/(2*delta*delta*meff)
		diag = hbar*hbar/(meff*delta*delta)
		D = diag + V
		E = s
		call eigs(Nthick, D, E, psi, LDZ, WORK, INFO)
		if (INFO.NE.0) then
			print*,"Warning: Wrong eigenvalues were calculated."
		end if
		psi(1, :) = 0
		psi(Nthick, :) = 0
		psi = psi*1E5
	end subroutine schrodinger

!The star of the show, namely the 'Poisson' subroutine. Currently it has two options:
!If input = 1, then that means the user chose to run only the Schrodinger-Poisson loop.
!If input = 2, then the user chose to calculate NEGF-Poisson.
!The difference between them is that in 1), only an approximate value for the electron
!density is calculated. However in 2), 'Poisson' takes in the electron density that was
!calculated in the NEGF loop.

!In any case, when using 'layers' I defined in 'Solver', the variables tempnum, tempsize etc.
!make sure that the largest quantum well is doped with 'rhod*e0'. The reason I did it like this
!is because in 1), Nthick = 600. But in 2), Nthick = 100 because a lower resolution is used.

!Nevertheless, in 1) the Fermi-Dirac factor is used, and I have accounted for huge variables.
!If the exponent is extremely large, then the term '1' in the logarithm is negligible.
!So the expression is simplified to avoid infinity terms in the electron density.
!If the exponent is not very large, then the expression is calculated normally.

!The last part involves using the expression for rho, and then using rho to calculate
!the new potential with FDM. This requires matrix inversion, which is why 'dgesv' is
!used from lapack. This takes in rho, and gives out the potential in rho. So the
!new potential Vt is then equal to rho.

	subroutine poisson(mu, Nthick, n, INFOTWO, input, meff, pi, hbar, el0, kb, rhod, eps, T, psi, D, Vt, nd)
		implicit none
		real*16			:: infinity
		double precision		:: var, tempnum, tempiter
		integer, intent(in)		:: Nthick, n
		integer			:: it, LDA, LDB, NRHS, IPIV(Nthick), INFOTWO, jt, &
					iter, input, tempsize, tempend, zt, a
		double precision, intent(in)	:: meff, eps, hbar, kb, T, mu, el0, rhod, pi
		double precision		:: D(Nthick), rhotemp(Nthick, 1), ni(n), rho(Nthick), & 
						nd(Nthick, n), integral, psi(Nthick, Nthick), &
						Mat(Nthick, Nthick), ndtemp(Nthick), nitemp(Nthick)
		double precision, intent(out)		:: Vt(Nthick)
		external dgesv, trapez
		integral = 0
		a = 1
		ndtemp = 0
		rho = 0
		rhotemp = 0
		Mat = 0
		tempnum = dble(Nthick)*dble(0.45)
		tempsize = floor(tempnum)
		tempiter = dble(Nthick)*dble(0.7)
		tempend = floor(tempiter)
		LDA = Nthick
		LDB = Nthick
		NRHS = 1
		do zt = tempsize, tempend
			rhotemp(tempsize, 1) = el0*rhod
		end do
		if (input.EQ.1) then
			do jt = 1, n
				var = log(1+exp((mu-D(jt))/kb*T))
				infinity = huge(var)
				if (var.GT.infinity) then
					ni(jt) = (meff/(pi*hbar*hbar))*(mu-D(jt))
				else
					ni(jt) = (meff/(pi*hbar*hbar))*kb*T*log(1+exp((mu-D(jt))/(kb*T)))
				end if
			end do
			do it = 1, n
				nitemp(:) = ni(it)*abs(psi(:, it))**2
				call trapez(nitemp, a, Nthick, integral, Nthick)
				rhotemp(:, 1) = rhotemp(:, 1) + el0*integral
				nitemp = 0
				integral = 0
			end do
		else if (input.EQ.2) then
			do it = 1, n
				ndtemp(:) = nd(:, it)*0.1
				call trapez(ndtemp, a, Nthick, integral, Nthick)
				rhotemp(:, 1) = rhotemp(:, 1) + el0*integral
				ndtemp = 0
				integral = 0
			end do

		end if
		do iter = 1, Nthick
			Mat(iter, iter) = 2
			if (iter.NE.Nthick) then
				Mat(iter, iter+1) = -1
				Mat(iter+1, iter) = -1
			end if
		end do
		rhotemp = el0*rhotemp/dble(12.9)*eps
		rho(:) = rhotemp(:, 1)
		call dgesv(Nthick, NRHS, Mat, LDA, IPIV, rho, LDB, INFOTWO)
		if (INFOTWO.NE.0) then
			print*,"Warning: Linear system in Poisson was calculated incorrectly."
		end if
		Vt = rho
	end subroutine poisson

!The subroutine 'gzero' calculates the g_{0,0} matrix described in my report. The mass of the left contact
!is completely arbitrary, so I would suggest correctly calculating this more precisely. This goes for
!the wave number in the contact as well. Furthermore, eta is used to ensure that the energies will converge.
!I recommend calculating the elements of g_{0, 0}, as I'm not too sure if they are completely correct.
!Two subroutines from lapack are used, 'zgetrf' and 'zgetri'. The first will factorize g_{0, 0} so that
!it can be used in 'zgetri' in order to find the inverse matrix. Note that these routines will only
!take in double complex matrices.

	subroutine gzero(gzz, SigC, N, hbar, me, pm) !Left connected gf near emitter
		implicit none
		integer, intent(in)		:: N
		integer				:: it, info, ipiv(N)
		double complex			:: gzz(N, N), work(N) !gzz = g_0,0, SigC = self-en. for contacts
		double complex, intent(out)		:: SigC(N, N)
		double precision		:: t, m, delta, hbar, gma, me, pm, eta
		external zgetrf
		external zgetri
		delta = 1E-10
		m = 0.01*me !Mass of contact
		eta = 1E-10 !Small imaginary constant to shift poles
		eta = tiny(eta)
		gma = 6.545E9 !Wave number in contact
		t = (hbar*hbar)/(2*m*delta*delta)
		gzz = 0
		SigC = 0
		gzz(1, 1) = cmplx(-t*cos(gma*delta) + pm, -t*sin(-gma*delta) + eta)
		do it = 1, N
			if (it.EQ.1) then
				gzz(1, 2) = t
				gzz(2, 1) = t
			else if (it.EQ.N) then
				gzz(N, N) = cmplx(-2*t*cos(gma*delta) + pm, eta)
			else
				gzz(it, it) = cmplx(-2*t*cos(gma*delta) + pm, eta)
				gzz(it+1, it) = t
				gzz(it, it+1) = t
			end if
		end do
		call zgetrf(N, N, gzz, N, ipiv, info)
		if (info.EQ.0) then
			call zgetri(N, gzz, N, ipiv, work, N, info)
			if (info.NE.0) then
				print*,"Matrix inversion could not be computed. (gzero)"
			end if
		else
			print*,"Matrix factorization was not computed correctly. (gzero)"
		end if
		SigC = t*gzz*t

	end subroutine gzero

!Subroutine 'gblocks', the beast or heart of this program. It will do the following:
!Take in g_{0,0}, left contact self-energy, calculate blocks of the small retarded
!green function, calculate the right contact self-energy using the next last block of
!grl by the same method as for g_{0, 0} ('zgetrf' and 'zgetri'),
!take in a total self-energy Sigma and add to it the contact self-energies,
!and then finally calculate the advanced, electron and hole (small)
!correlation Green functions. The expressions for the latter two should
!be correct, but I advise checking it a few times in order to be sure.

!It is left as a challenge for the user to find a more efficient algorithm
!that calculates the blocks throughout the system. This would allow for much
!higher resolution, and more precise solutions to the NEGF calculations.

	subroutine gblocks(gzero, grl, gal, gnl, gpl, SigC, SigCr, Sigma, Sigin, &
				Sigout, N, meff, me, hbar, en, pm, V) !Blocks of left connected gfs
		implicit none
		integer, intent(in)		:: N
		double complex			:: gzero(N, N), SigCr(N, N), grl(N, N, N), tempgl(N, N), work(N), grr(N, N), &
						 SigC(N, N), gal(N, N, N), Sigma(N, N, N), gnl(N, N, N), gpl(N, N, N), &
						Sigin(N, N, N), Sigout(N, N, N)
		double precision		:: t, meff, hbar, delta, V(N), en, gma, me, te, pm
		integer				:: it, info, ipiv(N)
		external zgetrf
		external zgetri
		delta = 1E-10
		t = (hbar*hbar)/(2*meff*delta*delta)
		te = (hbar*hbar)/(2*0.01*me*delta*delta)
		grl = 0 !Ret gf
		gnl = 0 !Electron correlation gf
		gal = 0 !Adv. gf
		gpl = 0 !Hole correlation gf
		tempgl = 0
		SigCr = 0
		Sigma(1, :, :) = Sigma(1, :, :) + SigC(:, :)
		tempgl(:, :) = en - Sigma(1, :, :) - V(1) - t - t*gzero(:, :)*t
		call zgetrf(N, N, tempgl, N, ipiv, info)
		if (info.EQ.0) then
			call zgetri(N, tempgl, N, ipiv, work, N, info)
			if (info.NE.0) then
				print*, "Matrix inversion could not be computed. (gblocks, 1st)"
			end if
		else
			print*, "Matrix factorization was not computed correctly. (gblocks, 1st)"
		end if
		grl(1, :, :) = tempgl(:, :)
		do it = 2, N-1
			tempgl = 0
			tempgl(:, :) = en - V(it) - t - t*grl(it-1, :, :)*t - Sigma(it, :, :)
			call zgetrf(N, N, tempgl, N, ipiv, info)
			if (info.EQ.0) then
				call zgetri(N, tempgl, N, ipiv, work, N, info)
				if (info.NE.0) then
					print*, "Matrix inversion could not be computed. (gblocks, 2nd)"
				end if
			else
				print*, "Matrix factorization was not computed correctly. (gblocks, 2nd)"
			end if
			grl(it, :, :) = tempgl(:, :)
		end do
		call grcont(grr, SigCr, V, N, hbar, me, pm) !To get self-en. for collector
		tempgl = 0
		Sigma(N, :, :) = Sigma(N, :, :) + SigCr(:, :)
		tempgl(:, :) = en - Sigma(N, :, :) - V(N) - t - t*grl(N-1, :, :)*t
		call zgetrf(N, N, tempgl, N, ipiv, info)
		if (info.EQ.0) then
			call zgetri(N, tempgl, N, ipiv, work, N, info)
			if (info.NE.0) then
				print*, "Matrix inversion could not be computed. (gblocks, 3rd)"
			end if
		else
			print*, "Matrix factorization was not computed correctly. (gblocks, 3rd)"
		end if
		grl(N, :, :) = tempgl(:, :)
		gal(:, :, :) = conjg(grl)
		gnl(1, :, :) = matmul(grl(1, :, :), matmul(Sigin(1, :, :), gal(1, :, :)))
		do it = 2, N
			gnl(it, :, :) = matmul(grl(it, :, :), matmul(Sigin(it, :, :), gal(it, :, :))) + &
				       & t*matmul(grl(it, :, :), matmul(gnl(it, :, :), gal(it, :, :)))*t
		end do
		gpl(1, :, :) = matmul(grl(1, :, :), matmul(Sigout(1, :, :), gal(1, :, :)))
		do it = 2, N
			gpl(it, :, :) = matmul(grl(it, :, :), matmul(Sigout(it, :, :), gal(it, :, :))) + &
					& t*matmul(grl(it, :, :), matmul(gpl(it, :, :), gal(it, :, :)))*t
		end do
	end subroutine gblocks

!Subroutine 'grcont' calculates the right contact self-energy. The difference between this and 'gzero' is that it relies
!on the conduction band difference (dV = V(1) - V(N)). This is probably a bit wrong, so I would suggest to look at the
!expression for the right contact self-energy. 

	subroutine grcont(grr, SigCr, V, N, hbar, me, pm)
	implicit none
	integer, intent(in)		:: N
	double complex			:: grr(N, N), SigCr(N, N), work(N)
	double precision, intent(in)		:: V(N), hbar, me
	double precision			:: m, dV, delta, gma, t, pm, eta
	integer				:: it, info, ipiv(N)
	external zgetrf
	external zgetri
	eta = 1E-10
	eta = tiny(eta)
	m = 0.01*me
	gma = 6.545E9
	delta = 1E-20
	t = (hbar*hbar)/(2*m*delta*delta)
	grr = 0
	SigCr = 0
	dV = V(1) - V(N)
	do it = 1, N-1
		grr(it, it) = cmplx(dV - 2*t*cos(gma*delta) + pm, eta)
		grr(it+1,it) = t
		grr(it, it+1) = t
	end do	
	grr(N, N) = cmplx(dV - t*cos(gma*delta) + pm, -t*sin(-gma*delta) + eta)
	call zgetrf(N, N, grr, N, ipiv, info)
	if (info.EQ.0) then
		call zgetri(N, grr, N, ipiv, work, N, info)
		if (info.NE.0) then
			print*,"Matrix inversion could not be computed. (grcont)"
		end if
	else
		print*,"Matrix factorization was not computed correctly. (grcont)"
	end if
	SigCr = t*grr*t

	end subroutine grcont

!The subroutine 'GFs' calculates the fully-connected Green's functions. Additionally, if calcnd = 1 then
!it will calculate the electron density.

	subroutine GFs(grl, gal, gnl, gpl, Gretl, Gecorl, Gpdiag, N, hbar, meff, pi, nlvl, nd, calcnd) !Calculates right connected big GFs
	implicit none
	integer, intent(in)		:: N, nlvl
	integer				:: it, jt, calcnd
	double complex			:: grl(N, N, N), gal(N, N, N), Gretl(N, N, N)
	double complex, intent(in)			:: gnl(N, N, N), gpl(N, N, N)
	double complex			:: Gretu(N-1, N, N), Gretd(N-1, N, N), &
					 Gadvl(N, N, N), Gadvd(N-1, N, N), Gadvu(N-1, N, N)
	double complex, intent(out)		:: Gecorl(N, N, N), Gpdiag(N, N, N)
	double precision		:: t, delta, hbar, meff, nd(N, nlvl), pi, a
	Gretl = 0
	Gadvl = 0
	Gadvd = 0
	Gadvu = 0
	Gretd = 0
	Gretu = 0
	Gecorl = 0
	Gpdiag = 0
	nd = 0
	delta = 1E-10
	a = 1E-2 !Because larger steps are used for nd
	Gretl(N, :, :) = grl(N, :, :)
	t = (hbar*hbar)/(2*meff*delta*delta)
	do it = N-1, 1, -1
	Gretl(it, :, :) = grl(it, :, :) + t*matmul(grl(it, :, :),(matmul(Gretl(it+1, :, :),grl(it, :, :))))*t !Diagonal elements of ret GF
	Gretd(it, :, :) = -t*matmul(Gretl(it+1, :, :), grl(it, :, :)) !Sub-diagonal of ret GF
	Gretu(it, :, :) = -t*matmul(grl(it, :, :), Gretl(it+1, :, :)) !Super-diagonal of ret GF
	end do
	Gadvl = conjg(Gretl) !Diagonal of adv GF
	Gadvd = conjg(Gretu) !Sub-diagonal of adv GF
	Gadvu = conjg(Gretd) !Super-diagonal of adv GF
	Gecorl(N, :, :) = gnl(N, :, :)
	do it = N-1, 1, -1
		Gecorl(it, :, :) = gnl(it, :, :) + t*matmul(grl(it, :, :), matmul(Gecorl(it+1, :, :), gal(it, :, :)))*t &
			& - t*matmul(gnl(it, :, :), Gadvl(it, :, :)) - t*matmul(Gretu(it, :, :), gnl(it, :, :)) !Diagonal of Gn/G<
		Gpdiag(it, :, :) = gpl(it, :, :) + t*matmul(grl(it, :, :), matmul(Gpdiag(it+1, :, :), gal(it, :, :)))*t &
			& - t*matmul(gpl(it, :, :), Gadvl(it, :, :)) - t*matmul(Gretu(it, :, :), gpl(it, :, :)) !Diagonal of Gp/G>
	end do
	if (calcnd.EQ.1) then
		do it = 1, N
			do jt = 1, nlvl
				nd(it, jt) = aimag(Gecorl(it, jt, jt))/(pi*a)
			end do
		end do
	end if
	end subroutine GFs

!The subroutine 'negfPot' simply resizes the potential initially calculated from
!layers, so that the resolution of the GFs is respected. It will not calculate
!a new potential if the GF dimensions are larger.

	subroutine negfPot(V, Vgf, N, Ngf)
	implicit none
	integer, intent(in)				:: N, Ngf
	double precision, intent(in)		:: V(N)
	double precision, intent(out)		:: Vgf(Ngf)
	double precision			:: tempnum
	integer					:: it, jt, multnum, store
	if (Ngf.GT.N) then
		print*, "Error: Wrong dimensions in negfPot."
		return
	else
		store = 0
		jt = 1
		Vgf = 0
		tempnum = dble(N)/dble(Ngf)
		multnum = floor(tempnum)
		Vgf(1) = V(1)
		do it = 1, Ngf
			Vgf(it) = V(multnum*it)
		end do
		do it = 1, Ngf
			if (Vgf(it).EQ.0) then
				store = it
				exit
			end if
		end do
		if (store.NE.0) then
			do it = store, Ngf
				Vgf(it) = V(multnum*Ngf + jt)
				jt = jt + 1
			end do
		end if
	end if
	end subroutine negfPot

!The following subroutine does the opposite of negfPot.
!It takes in the GF potential, and resizes it into high res.

	subroutine revPot(V, Vgf, N, Ngf)
	implicit none
	integer, intent(in)				:: N, Ngf
	double precision, intent(out)		:: V(N)
	double precision, intent(in)		:: Vgf(Ngf)
	double precision			:: tempnum
	integer					:: it, jt, addnum, store
	store = 0
	V = 0
	tempnum = dble(N)/dble(Ngf)
	addnum = floor(tempnum)
	do it = 1, Ngf
		do jt = 1, addnum
			V(jt + store) = Vgf(it)
		end do
		store = store + addnum
	end do
	return
	end subroutine revPot

!This subroutine 'inelselfen' calculates the inelastic phonon self-energy as given in the report.
!It takes in electron/hole correlation GFs for the two energy levels E-hw and E+hw.
!The expression is most likely wrong, and I encourage the user to correct this and improve it.

	subroutine inelselfen(Gecorminhw, Gecormaxhw, Sigma, Dp, hw, kb, T, N, horn) !Needs to be completely revised
	implicit none
	integer, intent(in)		:: N
	double complex			:: Gecorminhw(N, N, N), Gecormaxhw(N, N, N), Sigma(N, N, N)
	double precision		:: hw, Dp, nb, kb, T
	integer				:: it, horn
	nb = 1/(exp(hw/kb*T) - 1) !Bose distribution
	Sigma = 0
	if (horn.EQ.0) then
		do it = 1, N
			Sigma(it, :, :) = Sigma(it, :, :) + 2*Dp*(nb*Gecorminhw(it, :, :) + (nb+1)*Gecormaxhw(it, :, :))
		end do
	else
		do it = 1, N
			Sigma(it, :, :) = Sigma(it, :, :) + 2*Dp*((nb+1)*Gecorminhw(it, :, :) + nb*Gecormaxhw(it, :, :))
		end do	
	end if	

	end subroutine inelselfen

!The subroutine 'elselfen' calculates the elastic self-energy, and this expression
!is also probably not entirely correct.

	subroutine elselfen(G, Sigma, Dp, N)
		implicit none
		integer, intent(in)		:: N
		double complex			:: Sigma(N, N, N), G(N, N, N)
		double precision		:: Dp
		integer			:: it
		Sigma = 0
		do it = 1, N
			Sigma(it, :, :) = Dp*G(it, :, :)
		end do
		Sigma = Sigma*1E-16
		
	end subroutine elselfen

!This subroutine 'continselfen' calculates the scattering in and out
!self-energies from/to the contacts. It is dependent on the fermi
!levels in the contacts, and I have assumed these values.

	subroutine continselfen(SigC, SigCr, V, N, mu, kb, T, SiginC, SiginCr, SigoutC, SigoutCr) !Scattering in/out from contacts
	!Assumes constant Fermi factor in contacts
	implicit none
	integer, intent(in)			:: N
	double complex, dimension(N, N), intent(in)		:: SigC, SigCr
	double complex, dimension(N, N)		:: SiginC, SiginCr, SigoutC, SigoutCr
	double precision			:: fr, fl, kb, T, mu, er, el, V(N), var1, var2
	el = 0.1*V(1)
	er = 0.1*V(N)
	SiginC = 0
	SiginCr = 0
	SigoutC = 0
	SigoutCr = 0
	var1 = exp(el/kb*T)
	var2 = exp(er/kb*T)
	fl = 1/(var1 + 1)
	fr = 1/(var2 + 1)
!	var1 = exp((el)/kb*T)
!	infty = huge(var1)
!	if (var1.GT.infty) then
!		fl = 0
!	else
!		fl = 1/(var1 + 1)
!	end if
!	var2 = exp((er)/kb*T)
!	infty = huge(var2)
!	if (var2.GT.infty) then
!		fr = 0
!	else
!		fr = 1/(var2 + 1)
!	end if
	SiginC(:, :) = -2*aimag(SigC(:, :))*fl
	SiginCr(:, :) = -2*aimag(SigCr(:, :))*fr
	SigoutC(:, :) = -2*aimag(SigC(:, :))*(1-fl)
	SigoutCr(:, :) = -2*aimag(SigCr(:, :))*(1-fr)
	end subroutine continselfen

	!I wrote this custom trapezoidal integration method which assumes the
	!limits to be integrals, so that for an array of x-values for instance
	!you can integrate a function f of the same size over all those values.

	subroutine trapez(f, a, b, integral, n)
		implicit none
		integer		:: n, j, a, b, h
		double precision	:: f(n), integral
		h = (b - a + 1)/n !+1 because indexing starts at 1
		integral = 0.5*h*(f(a) + f(b))
		do j = 1, n-1
			integral = integral + h*f(a + j*h)
		end do
		return
	end subroutine trapez




