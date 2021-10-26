# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import ctc
import rt
import wtr
#import ASTMC534
import espuma_elastomerica
import ASTMC552
import ASTMC591
#import ASTMC1728a
import aerogel_cryogel
import cst, custos

# =============================================================================
# Custos.
# =============================================================================

#Custo de energia perdida trazido para valor atual em $/(ano.m^2).
def CE(q, N, CEE, eta, COP):#, n, i, delta):

    CE = (q*N*CEE)/(1000*eta*COP)

    #j = ((1 + i)/(1 + delta)) - 1

    #f = (((1 + j)**n) - 1)/(j*((1 + j)**n))

    #CEVA = f*CE

    return (CE)

#Custo de manutenção do isolamento em $/(ano.m^2).
def CM(CI, tm):#, n, i):

    CM = tm*CI

    #f = (((1+i)**n)-1)/(i*((1+i)**n))

    #CMVA = f*CM

    return (CM)

# =============================================================================
# Função principal.
# =============================================================================

def iso_tubes(di, de, Ti, Ta, h_fld, lmd_tube, U, H, z, eps, m, c, Dt_max, R_rev, RH, N, CEE, eta, COP, n, i, delta, tm):

    Di = de #Diamentro interno do isolante = diâmetro externo da tubulação.

    #Lista de espessuras consideradas.
    #As listas LE2 e LE3 consideram convertem as espessuras de mm para metro.
    LE2 = [(6.35*x)/1000 for x in range(1, 17)]
    #LE3 = [(12.7*x)/1000 for x in range(1, 17)]

    #Lista de espessuras consideradas arredondadas.
    LE_Disp2 =  [(int(x*1000)) for x in LE2]
    #LE_Disp3 =  [(int(x*1000)) for x in LE3]

    #Lista de espessuras consideradas em polegadas.
    LE_Disp_Imp2 = [0.25*x for x in range (1, 17)]
    #LE_Disp_Imp3 = [0.50*x for x in range (1, 17)]

    #Lista de funções de condutividade térmica.

    #LLMD2 = [ASTMC1728a.lamed,
    #         ASTMC534.lamed]
    LLMD2 = [aerogel_cryogel.aerogel_condutividade_term,
            espuma_elastomerica.espuma_condutividade_term]
    #LLMD3 = [ASTMC552.lamed,
    #         ASTMC591.lamed]

    #Lista de nomes de isolantes.
    #Isolantes flexívies.
    Lnm2 = ['Aerogel',
            'Espuma Elastomérica']
    #Isolantes rígidos.
    #Lnm3 = ['Vidro Celular',
    #        'Poliisocianurato']

    #Número de espessuras consideradas.
    ne2 = len(LE2)
    #ne3 = len(LE3)

    #Número de isolantes considerados.
    ni2 = len(LLMD2)
    #ni3 = len(LLMD3)

    #Lista de espessuras para o DataFrame.
    #CONFIRMAR O VALOR DA ESPESSURA
    #DO ISOLANTE A SER SUBSTITUIDO
    #LE_Disp_DF = [0] + LE_Disp2*ni2 #+ LE_Disp3*ni3
    #LE_Disp_Imp_DF = [0] + LE_Disp_Imp2*ni2 #+ LE_Disp_Imp3*ni3
    LE_Disp_DF = LE_Disp2*ni2 #+ LE_Disp3*ni3
    LE_Disp_Imp_DF = LE_Disp_Imp2*ni2 #+ LE_Disp_Imp3*ni3

    #Lista de nomes para o DataFrame.
    LNM = ['Isolamento Atual']

    for d in Lnm2:
        LNM += [d]*ne2
    #for d in Lnm3:
    #    LNM += [d]*ne3

    #Lista de soluções para a temperatura na face externa.
    LTe = []

    #Lista de soluções para a condutividade térmica.
    Lslmd = [np.NaN] #OBS : Não será exibida/adicionada no dataframe

    #Lista de soluções para o taxa de transferência de calor.
    Lq = []

    #Lista de diâmetros externos.
    LDe = [de]

    #Lista de resistências térmicas.
    LR = [] #Lembrar de adicionar como coluna no DataFrame

    #Lista de variações de temperatura do fluido.
    LVT = []

    #Lista de temperaturas na saída.
    LTS = []

    if ((U != 0)):
        for flmd in LLMD2:
            for E in LE2:
                De = Di + 2*E
                #err = errf(di, de, Ti, Ta, h_fld, lmd_tube, U, H, eps, E, flmd, R_rev)
                Te = wtr.TDP(Ta, RH)
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+\
                rt.rt_cond_cili(di,de,lmd_tube)+\
                rt.rt_cond_cili(Di,De,flmd)+\
                R_rev+((((rt.rt_conv_cili(De,ctc.hc_m_ch(U,De,Te,Ta))))+((rt.rt_crad_cili(De,ctc.hr(eps,Te,Ta)))))))]
                qdp = (Ta - Ti) / LR[-1]
                #qdp = (ctc.qc_m_ch(U, De, Te, Ta) + ctc.qr(eps, Te, Ta))*(np.pi*(De))

                TDi = Ti - qdp*(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube))
                #Te = TDi + qdp * (rt.rt_cond_cili(Di, De, flmd) + R_rev)
                Lslmd = Lslmd + [flmd]
                Lq = Lq + [qdp]
                LTe = LTe + [Te]
                LDe = LDe + [De]
                '''A divisão por z e os expoentes (-1) foram eliminados por enquanto
                   LR = LR + [ (rt.rt_conv_cili(di,h_fld)+ \
                            rt.rt_cond_cili(di,de,lmd_tube)+ \
                            rt.rt_cond_cili(Di,De,flmd(Tiso))+ \
                            R_rev+ \
                            ((((rt.rt_conv_cili(De,ctc.hc_m_ch(U,De,Te,Ta)))**(-1))+((rt.rt_crad_cili(De,ctc.hr(eps,Te,Ta)))**(-1)))**(-1)))/z ]'''

    if ((U == 0)):
        for flmd in LLMD2:
            for E in LE2:
                De = Di + 2*E
                LR = LR + [(rt.rt_conv_cili(di,h_fld)+\
                           rt.rt_cond_cili(di,de,lmd_tube)+\
                           rt.rt_cond_cili(Di,De,flmd)+\
                           R_rev+\
                           ((((rt.rt_conv_cili(De,ctc.hc_n_ch(U,De,Te,Ta))))+((rt.rt_crad_cili(De,ctc.hr(eps,Te,Ta)))))))]
                qdp = (Ta - Ti) / LR[-1]
                #qdp = (ctc.qc_n_ch(U, De, Te, Ta) + ctc.qr(eps, Te, Ta))*(np.pi*(De))
                TDi = Ti - qdp*(rt.rt_conv_cili(di,h_fld)+rt.rt_cond_cili(di,de,lmd_tube))
                Te = TDi + qdp * (rt.rt_cond_cili(Di, De, flmd) + R_rev)
                Lslmd = Lslmd + [flmd]
                Lq = Lq + [qdp]
                LTe = LTe + [Te]
                LDe = LDe + [De]
                '''A divisão por z e os expoentes (-1) foram eliminados por enquanto
                    LR = LR + [(rt.rt_conv_cili(di,h_fld)+\
                            rt.rt_cond_cili(di,de,lmd_tube)+\
                            rt.rt_cond_cili(Di,De,flmd(Tiso))+\
                            R_rev+\
                            ((((rt.rt_conv_cili(De,ctc.hc_n_ch(U,De,Te,Ta)))**(-1))+((rt.rt_crad_cili(De,ctc.hr(eps,Te,Ta)))**(-1)))**(-1)))/z]'''

    #Diâmetros externos em milímetros.
    LDe_Disp = [1000*D for D in LDe]

    #Lista de soluções para a temperatura na face externa em °C.
    Lte = [(Te - 273.15) for Te in LTe]

    #Lista de q através de LMTD e R.
    LVT = list(map(lambda x: Ti - (Ta - (Ta - Ti)*np.exp(-1/(m*c*x))),LR))

    LTSK = list(map(lambda x: Ti - x, LVT))
    LLMTD = [((x - Ta) - (Ti - Ta))/(np.log((x - Ta)/(Ti - Ta))) for x in LTSK]
    ZTR = zip(LLMTD, LR)
    #Lq está em W/m, LR deve ser multiplicado por z.
    Lq = [x[0]/(x[1]) for x in ZTR]
    #Lq = [x[0]/(x[1]*z) for x in ZTR] Verificar se o z deve ser removido.

    #Análise econômica.
    #Custo de energia perdida atualizada.
    LCE = [CE(abs(q), N, CEE, eta, COP) for q in Lq]
    #Custo de investimento.

    #O investimento é o custo do material + custo mão de mao_de_obra
    #A lista LFCI2 vai pegar os dados do módulo (biblioteca) custos
    # custo_do_material * área_da_tubulação + custo da mão de obra *comprimento da tubulação
    LFCI2 = []
    LFCI2 += [custos.aerogel_custos[esp] for esp in LE_Disp_Imp2] #Alterar
    LFCI2 += LFCI2 + [custos.espuma_custos[esp] for esp in LE_Disp_Imp2]

    LCI = [np.nan]
    for ci in LFCI2:
        for E in LE2:
            LCI = LCI + [ci(Di, E)]
    '''for ci in LFCI3:
        for E in LE3:
            LCI = LCI + [ci(Di, E)]'''
    #Custo de manutenção.
    LCM = [CM(CI, tm) for CI in LCI]
    #Custo total.
    LCT = [x + y + z for (x,y,z) in zip(LCE, LCI, LCM)]

    #Exlusão de casos reprovados.

    if Dt_max > 0:
        LVT = list(map(lambda x: Ti - (Ta - (Ta - Ti)*np.exp(-1/(m*c*x))),LR))
        for i in range(len(LVT) - 1, 0, -1):
            if abs(LVT[i]) > Dt_max:
                del LVT[i]
                del Lq[i]
                del LNM[i]
                del LE_Disp_DF[i]
                del LE_Disp_Imp_DF[i]
                del LDe_Disp[i]
                del Lte[i]
                del Lslmd[i]
                del LR[i]
                del LCEVA[i]
                del LCI[i]
                del LCMVA[i]
                del LCT[i]

    LFL = [wtr.fl(x, Ta, RH) for x in LTe]
    for i in range(len(LFL) - 1,  0, -1):
        if (LFL[i]) == "Sim":
            del LVT[i]
            del Lq[i]
            del LNM[i]
            del LE_Disp_DF[i]
            del LE_Disp_Imp_DF[i]
            del LDe_Disp[i]
            del Lte[i]
            del Lslmd[i]
            del LR[i]
            del LCEVA[i]
            del LCI[i]
            del LCMVA[i]
            del LCT[i]

    #Organização dos dados em um DataFrame.
    Disp = pd.DataFrame({'Material' : LNM,
                         'Espessura \n [mm]' : LE_Disp_DF,
                         'Espessura \n [pol]' : LE_Disp_Imp_DF,
                         'Diâmetro \n Externo [mm]' : LDe_Disp,
                         'Temperatura na \n Face Externa [°C]' : Lte,
                         'Condutividade Térmica \n do Isolante [W/(m.k)]': Lslmd})

    LVT = list(map(lambda x: Ti - (Ta - (Ta - Ti)*np.exp(-1/(m*c*x))),LR))
    LVT_Disp = [abs(T) for T in LVT]
    Disp['Variação de Temperatura \n do Fluido [°C]'] = LVT_Disp
    LTS = list(map(lambda x: Ti - 273.15 - x, LVT))
    Disp['Temperatura do Fluido \n na Saída [°C]'] = LTS

    Disp['Fluxo de Calor \n [W/m]'] = Lq

    Disp['Custo de Energia \n [$/m]'] = LCEVA
    Disp['Custo de Investimento \n [$/m]'] = LCI
    Disp['Custo de Manutenção \n [$/m]'] = LCMVA
    Disp['Custo Total \n [$/m]'] = LCT

    #Disp = Disp.sort_values(by=['Custo Total \n [$/m]'])

    Disp = Disp.iloc[[0]].append(Disp.iloc[1:].sort_values(by=['Custo Total \n [$/m]']))

    if (len(Lq) >= 24):
        LA = ['Ti', 'Ponto de Orvalho', 'di', 'de', 'Ta', 'h_fld', 'lmd_tube', 'U', 'H', 'z', 'eps', 'm', 'c', 'Dt_max', 'R_rev', 'RH', 'N', 'CEE', 'eta', 'COP', 'n', 'i', 'delta', 'tm']
        LB = [Ti, wtr.TDP(Ta, RH), di, de, Ta, h_fld, lmd_tube, U, H, z, eps, m, c, Dt_max, R_rev, RH, N, CEE, eta, COP, n, i, delta, tm]
        while len(LA) < len(Lq):
            LA.append(np.nan)
            LB.append(np.nan)
        Disp.insert(loc = 0, column = 'Input', value = LB)
        Disp.insert(loc = 0, column = 'Varáivel', value = LA)

    return Disp
