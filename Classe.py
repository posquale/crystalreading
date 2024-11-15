class System:
    def __init__(self, name, trigger_up_x, trigger_up_y, trigger_down_x, trigger_down_y, trigger_updown_distance, triggerup_scint_distance,
                 rifrazione_ambiente, verticale, filtroA, filtroB):
        self.name = name
        
        self.tupx = trigger_up_x
        self.tupy = trigger_up_y
        self.tdownx = trigger_down_x
        self.tdowny = trigger_down_y
        self.dt = trigger_updown_distance
        
        self.tup_scint = triggerup_scint_distance
        self.rifr = rifrazione_ambiente

        self.v = verticale
        self.fa = filtroA
        self.fb = filtroB


class sipm:
    def __init__(self, sipmA, eff_geom_A, rifrazioneA, lambdamin_A, lambdamax_A, sipmB, eff_geom_B, rifrazioneB, lambdamin_B, lambdamax_B):
        self.na = sipmA
        self.nb = sipmB

        self.mina = lambdamin_A
        self.maxa = lambdamax_A
        self.minb = lambdamin_B
        self.maxb = lambdamax_B
        
        self.ra = rifrazioneA
        self.rb = rifrazioneB
        self.ega = eff_geom_A
        self.egb = eff_geom_B


class filtri:
    def __init__(self, filtroA, filtroB):
        self.na = filtroA
        self.nb = filtroB
        

class Scintillatore:
    def __init__(self, name, dim_z, dim_y, dim_x, density, light_yield, rifrazione, radlen, ene_loss, emissione, trasmittanza, lambdarange):
        self.name = name
        
        self.z = dim_z
        self.y = dim_y
        self.x = dim_x
        
        self.density = density
        self.ly = light_yield
        self.rifrazione = rifrazione
        self.radlen = radlen
        self.dedx = ene_loss
        
        self.e = emissione
        self.t = trasmittanza
        self.range = lambdarange


#Sistemi
sysVnof=System('demoverticaleNOFILTRI', 0.049, 0.049, 0.049, 0.049, 0.12, 0.025, 1.0003, True, False, False)
sysOnof=System('demoorizzontaleNOFILTRI', 0.049, 0.049, 0.049, 0.049, 0.12, 0.085, 1.0003, False, False, False)
sysV=System('demoverticale', 0.049, 0.049, 0.049, 0.049, 0.12, 0.025, 1.0003, True, True, True)
sysO=System('demoorizzontale', 0.049, 0.049, 0.049, 0.049, 0.12, 0.085, 1.0003, False, True, True)

sys4=System('quarta misura tesi', 0.049, 0.049, 0.049, 0.049, 0.12, 0.085, 1.0003, False, False, True)


#Scintillatori, dim x, y e z sono le dimensioni del cristallo in m, densità in Kg/m^3, light yield in #/MeV, radlen in m, ene_loss in MeV / (kg/ m^2)
dim_z = 0.05
dim_y = 0.012
dim_x = 0.012

BGO=Scintillatore('BGO', dim_z, dim_y, dim_x, 7130, 8200, 2.15, 0.011, 0.1272, 'dati/scintillatori/emissione/BGO.csv', 'dati/scintillatori/trasmittanza/BGO.asc', [200, 900])
PWO=Scintillatore('PWO', dim_z, dim_y, dim_x, 8280, 190, 2.16, 0.009, 0.1225, 'dati/scintillatori/emissione/PWO.csv', 'dati/scintillatori/trasmittanza/PWO.asc', [200, 900])
