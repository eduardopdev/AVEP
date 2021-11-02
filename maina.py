# -*- coding: utf-8 -*-
import iba
import ctc
import wtr

#Diâmetro interno da tubulação em m.
di = 0.254

#Diâmetro externo da tubulação em m.
de = 0.300

#Temperatura do fluido em K.
Ti = 6 + 273.15

#Temperatura do ar ambiente em K.
Ta = 29 + 273.15

#Condutividade térmica do material da tubulaçãoem SI.
lmd_tube = 60

#Velocidade do vento em m/s. Caso convecção natural, U = 0.
U = 3

#Altura da tubulação em m, caso vertical.
H = 0

#comprimento da tubulação em m.
z = 100

#Emissividade da superfície. Pode ser importada de ctc.py.
eps = ctc.emissividade('Alumínio')

#Vazão mássica de água em kg/s.
m = 138

#Coeficiente de convecção do escoamento interno à tubulação em SI.
h_fld = wtr.hc(m, di, Ti)

#Calor específico em J/(kg.K) da água. Pode ser importado de wtr.py.
c = wtr.cp(Ti)

#Variação máxima admissível de temperatura da água em °C. Ignorado se == 0.
Dt_max = 0

#Resitência térmica da seção de revestimento protetivo. m.K/W.
R_rev = 1e-4

#Umidade relativa.
RH = 0.78
TF = 6.5 + 273.15
(N, CEE, eta, COP, n, i, tm) = (400, 0.12, 0.80, 4, 15, 0.12, 0.01)

if True:
    result = iba.iso_tubes(di, de, Ti, Ta, h_fld, lmd_tube, U, H, z, eps, m, c, Dt_max, R_rev, RH, N, CEE, eta, COP, n, i, tm, TF)
    result.to_excel("resultado_a.xlsx", sheet_name='MERCK', index = False, startcol = 0, freeze_panes = (2,2))