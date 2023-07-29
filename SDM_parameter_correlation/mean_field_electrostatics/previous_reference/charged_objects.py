class protein:
    def __init__(self, name, fasta, initial_radius, extra_charge=0.0, n_seq=1):
        self.name           = name
        self.fasta          = fasta
        self.initial_radius = initial_radius
        self.radius         = self.initial_radius
        self.extra_charge   = extra_charge
        self.n_seq          = n_seq
        return

    def get_seq(self):
        file = open(self.fasta, mode='r')
        contents = file.read().splitlines()
        seq = ""

        for seq_num in range(self.n_seq):
            for index in range(len(contents)):
                if index == 0:
                    pass
                else:
                    seq = seq + contents[index]
        self.seq = seq
        return

    def get_mass(self):
        mass = SeqUtils.molecular_weight(self.seq, seq_type='protein')
        self.mass = mass
        return

    def get_aa_count(self):
        # Creates a dictionary of amino acids and their occurrence number
        x = ProteinAnalysis(self.seq)
        aa_count = x.count_amino_acids()
        aa_count['N_term'] = 1
        aa_count['C_term'] = 1
        self.aa_count = aa_count
        return

    def get_charged_aa(self):
        # Creates a list of charged
        charged_aa = []
        for aa in self.aa_count:
            if aa in pK:
                charged_aa.append(aa)
        self.charged_aa = charged_aa
        return

    def get_surf_dens(self, aa):
        # Compute ligand density [mol m-2] from n [-] ligands on sphere of radius a [m]
        n = self.aa_count[aa]
        return n/(4*numpy.pi*(self.radius**2)*Na)

    def log_pK_surf_dens(self):
        self.log = {}
        for aa in self.charged_aa:
            self.log[aa] = pK[aa].copy()
            surf_dens    = self.get_surf_dens(aa)
            self.log[aa].append(surf_dens)
        return

    def m_log_pKeff_surf_dens(self, solution):
        # Log of charged aa:  pKa, z, surf_dens [mol m-2], pKa_eff, n ligands
        self.log_pK_surf_dens()
        self.log_pK_eff = m_pK_eff(self.log, solution)
        for aa in self.log_pK_eff:
            self.log_pK_eff[aa].append(self.aa_count[aa])
        return

    def m_get_z_protein(self, solution):
        self.m_log_pKeff_surf_dens(solution)
        ch = solution.ch
        z_acid = 0
        z_base = 0

        for aa in self.log_pK_eff:
            n = self.log_pK_eff[aa][4]
            pKa_eff = self.log_pK_eff[aa][3]
            Ka_eff  = 10**(-1.0*pKa_eff)
            if self.log_pK_eff[aa][1] == 0:  # acid
                z_acid += -1.0*n/(1.0+ch/Ka_eff)
            elif self.log_pK_eff[aa][1] == 1:  # base
                z_base += n/(1.0+Ka_eff/ch)
        self.m_z_acid = z_acid
        self.m_z_base = z_base
        self.m_z_protein = z_acid + z_base
        return

    def m_get_charge_dens(self, solution):
        # Charge density in [C m-2]

        # # Constant charge density (original model)
        # self.m_get_z_protein(solution)
        # self.m_charge_dens = self.m_z_protein*F/(4*numpy.pi*Na*self.radius**2)

        # # Charge regulation
        # self.m_log_pKeff_surf_dens(solution)
        # get_psi0(solution, self)
        # self.m_charge_dens = LHS_1([self.psi0], solution)

        # Charge density input for a given node
        self.m_charge_dens = self.m_charge_dens

        # print(self.m_charge_dens)
        return

    def m_get_charge_dens_dl(self, solution):
        self.m_get_charge_dens(solution)
        self.m_charge_dens_dl = solution.m_z*e*self.m_charge_dens*self.radius/(kT*solution.eps*eps0)
        return

    def m_get_theta_p0(self, solution):
        solution.get_kappa()
        self.m_get_charge_dens_dl(solution)
        fit = scipy.optimize.minimize(m_theta_p0_residual, x0=[0], args=(solution, self), bounds=[(-1000, 1000)], tol=1e-12)
        self.theta_0 = fit.x[0]
        # print(self.theta_0)
        return

    def m_get_capital_theta(self, solution, dpr):
        self.m_get_theta_p0(solution)
        self.capital_theta = m_capital_theta(dpr, solution, self)
        return

    def m_fp(self, solution):
        solution.get_kappa()
        kap_r = solution.kappa*self.radius
        num   = 1 - 1/kap_r + (1 + 1/kap_r)*numpy.exp(-2.0*kap_r)
        denom = 1 + 1/kap_r - (1 + 1/kap_r)*numpy.exp(-2.0*kap_r)     # 2016 Guelat:  1/kap_r is (+) {agrees with Hsu and Liu 2009};        2012 Guelat:  1/kap_r is (-)
        self.fp = num/denom
        return

class resin:
    def __init__(self, name, ligand, surf_dens):
        self.name      = name
        self.ligand    = ligand
        self.surf_dens = surf_dens # [mol m-2]
        return

    def log_pK_surf_dens(self):
        self.log = {}
        self.log[self.ligand] = pK[self.ligand].copy()
        self.log[self.ligand].append(self.surf_dens)
        return

    def m_log_pKeff(self, solution):
        self.log_pK_surf_dens()
        self.m_log_pK_eff = m_pK_eff(self.log, solution)
        return

    def m_get_charge_dens(self, solution):
        self.m_log_pKeff(solution)
        ch = solution.ch
        self.m_charge_dens = 0.0

        # Constant charge density (original model)
        for ligand in self.m_log_pK_eff:
            pKa_eff = self.m_log_pK_eff[ligand][3]
            Ka_eff = 10**(-1.0*pKa_eff)

            if self.m_log_pK_eff[ligand][1] == 0:  # acid
                sign = -1.0
                denom = (1.0+ch/Ka_eff)
            elif self.m_log_pK_eff[ligand][1] == 1:  # base
                sign = 1.0
                denom = (1.0+Ka_eff/ch)
            self.m_charge_dens += sign*F*self.surf_dens/denom       # Would need to change for multiple ligands

        # # Charge regulation
        # get_psi0(solution, self)
        # self.m_charge_dens = LHS_1([self.psi0], solution)

        return

    def m_get_charge_dens_dl(self, solution):
        self.m_get_charge_dens(solution)
        solution.get_kappa()
        self.m_charge_dens_dl = solution.m_z*e*self.m_charge_dens/(kT*eps0*solution.eps*solution.kappa)
        return

    def m_get_theta_r0(self, solution):
        solution.get_kappa()
        self.m_get_charge_dens_dl(solution)
        fit = scipy.optimize.minimize(m_theta_r0_residual, x0=[0], args=(solution, self), bounds=[(-1000, 1000)], tol=1e-12)
        self.theta_0 = fit.x[0]
        return

    def m_get_capital_theta(self, solution, dpr):
        self.m_get_theta_r0(solution)
        self.capital_theta = m_capital_theta(dpr, solution, self)
        return

class solution:
    def __init__(self, pH, ion_str, eps=80.1, m_z=1):
        self.pH      = pH
        self.ion_str = ion_str  # [M]
        self.eps     = eps      # [-]
        self.m_z     = m_z      # (Symmetrical) electorlyte valence
        return

    def get_ch(self):
        # Convert pH to bulk solution hydronium ion concentration [M]
        self.ch = 10.0**(-1.0*self.pH)
        return

    def get_kappa(self):
        # Convert ionic strength [M] into the inverse Debye length [m-1]
        self.kappa = numpy.sqrt(2*e**2*self.ion_str*Na*1.0e3/(kT*self.eps*eps0))
        return

    def get_cap_dif(self):
        # Compute the diffuse layer capacitance [C m-2 V-1]
        self.cap_dif = self.eps*eps0*self.kappa
        return
