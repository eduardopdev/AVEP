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
import vpl_revisado

# =============================================================================
# Custos.
# =============================================================================
def get_greater_y(xlist, ylist, material):
    if material=="Aerogel":
        for i in range(16):
            if i == 0:
                x = xlist[0]
                max = ylist[0]
            else:
                if ylist[i] > max:
                    x = i
    if material=="Espuma":
        for i in range(16,32):
            if i == 16:
                x = xlist[0]
                max = ylist[0]
            else:
                if ylist[i] > max:
                    x = i
    return x

def fesp_q(pltObj, LE, LQ):
    LQ_Aerogel = []
    LQ_Espuma = []
    for i in range(16):
        LQ_Aerogel.append(LQ[i])
    for i in range(16,32):
        LQ_Espuma.append(LQ[i])
    pltObj.plot(LE, LQ_Aerogel, "b.-", label = "Aerogel")
    pltObj.plot(LE, LQ_Espuma, "r.-", label = "Espuma")
    pltObj.title('Espessura (pol) x Taxa de Calor')
    plt.xlabel('polegadas')
    plt.ylabel('W')
    plt.legend() #exibir a legenda
    plt.show()

def fTe_q(pltObj, Te_aerogel, Te_espuma, q_aerogel, q_espuma):
    pltObj.plot(Te_aerogel, q_aerogel, "b", label = "Aerogel")
    pltObj.plot(Te_espuma, q_espuma, "r", label = "Espuma")
    pltObj.title('Temperatura x Calor')
    plt.xlabel('K')
    plt.ylabel('calor')
    plt.legend()
    plt.show()

def fTe_VPL(pltObj, Te_aerogel, Te_espuma, VPL_aerogel, VPL_espuma):
    pltObj.plot(Te_aerogel, VPL_aerogel, "b", label = "Aerogel")
    pltObj.plot(Te_espuma, VPL_espuma, "r", label = "Espuma")
    pltObj.title('Temperatura x VPL')
    plt.xlabel('K')
    plt.ylabel('VPL')
    plt.legend()
    plt.show()

def fTe_q_atual(pltObj, Te, q):
    pltObj.plot(Te, q, "b", label = "Atual")
    pltObj.title('Temperatura x Calor')
    plt.xlabel('K')
    plt.ylabel('VPL')
    plt.legend()
    plt.show()

def fTe_VPL_atual(pltObj, Te, VPL):
    pltObj.plot(Te, VPL, "b", label = "Atual")
    pltObj.title('Temperatura x VPL')
    plt.xlabel('K')
    plt.ylabel('VPL')
    plt.legend()
    plt.show()

def fesp_VPL(pltObj, LE, VPL):
    VPL_Aerogel = []
    VPL_Espuma = []
    for i in range(16):
        VPL_Aerogel.append(VPL[i])
    for i in range(16,32):
        VPL_Espuma.append(VPL[i])
    pltObj.plot(LE, VPL_Aerogel, "b.-", label = "Aerogel")
    pltObj.plot(LE, VPL_Espuma, "r.-", label = "Espuma")
    pltObj.title('Espessura x VPL')
    plt.xlabel('polegadas')
    plt.ylabel('VPL')
    plt.legend()
    plt.show()

def fesp_payback(pltObj, LE, payback):
    payback_Aerogel = []
    payback_Espuma = []
    for i in range(16):
        payback_Aerogel.append(payback[i])
    for i in range(16,32):
        payback_Espuma.append(payback[i])
    pltObj.plot(LE, payback_Aerogel, "b.-", label = "Aerogel")
    pltObj.plot(LE, payback_Espuma, "r.-", label = "Espuma")
    pltObj.title('Espessura x payback')
    plt.xlabel('polegadas')
    plt.ylabel('anos')
    plt.legend()
    plt.show()

def fn_vpl(pltObj, vpl_lista_aerogel, vpl_lista_espuma, n):
    n_lista = list(range(n+1))
    pltObj.plot(n_lista, vpl_lista_aerogel,"b.-", label="Aerogel")
    pltObj.plot(n_lista, vpl_lista_espuma,"r.-", label="Espuma")
    pltObj.title('Tempo (anos) x VPL')
    pltObj.xlabel('tempo')
    pltObj.ylabel('VPL')
    plt.legend()
    plt.show()

def fi_vpl(pltObj, vpl_aerogel, vpl_espuma, juros):
    pltObj.plot(juros, vpl_aerogel,"b", label="Aerogel")
    pltObj.plot(juros, vpl_espuma,"r", label="Espuma")
    pltObj.title('Taxa de Juros x VPL')
    pltObj.xlabel('Juros')
    pltObj.ylabel('VPL')
    plt.legend()
    plt.show()

def fn_vpl_atual(pltObj, vpl_lista_atual, n):
    n_lista = list(range(n+1))
    pltObj.plot(n_lista, vpl_lista_atual,"b", label="Material Atual")
    pltObj.title('Tempo (anos) x VPL')
    pltObj.xlabel('tempo')
    pltObj.ylabel('VPL')
    plt.legend()
    plt.show()

def fi_vpl_atual(pltObj, vpl_lista_atual, juros):
    pltObj.plot(juros, vpl_lista_atual,"b", label="Material Atual")
    pltObj.title('Taxa de Juros x VPL')
    pltObj.xlabel('Juros')
    pltObj.ylabel('VPL')
    plt.legend()
    plt.show()

#Custo de energia perdida trazido para valor atual em $/(ano.m^2).
def CE(q, N, CEE, eta, COP):
    CE = (q*N*CEE)/(1000*eta*COP)
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

def iso_tubes(di, de, Ti, Ta, h_fld, lmd_tube, U, H, z, eps, m, c, Dt_max, R_rev, RH, N, CEE, eta, COP, n, i, tm, TF):

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
    LNM = []
    for d in Lnm2:
        LNM += [d]*ne2
    #for d in Lnm3:
    #    LNM += [d]*ne3

    #Lista de soluções para a temperatura na face externa.
    LTe = []

    #Lista de soluções para a condutividade térmica.
    #Lslmd = [np.NaN] #OBS : Não será exibida/adicionada no dataframe
    Lslmd = []
    #Lista de soluções para o taxa de transferência de calor.
    Lq = []

    #Lista de diâmetros externos.
    #LDe = [de]
    LDe= []
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
                LR = LR + [ (rt.rt_conv_cili(di,h_fld,z)+ \
                    rt.rt_cond_cili(di,de,lmd_tube,z)+ \
                    rt.rt_cond_cili(Di,De,flmd,z)+ \
                    R_rev+ \
                    ((((rt.rt_conv_cili(De,ctc.hc_f_c(U,De,Te,Ta),z))**(-1))+((rt.rt_crad_cili(De,ctc.hr(eps,Te,Ta),z))**(-1)))**(-1)))]
                qdp = (Ta - Ti) / LR[-1]
                #qdp = (ctc.qc_m_ch(U, De, Te, Ta) + ctc.qr(eps, Te, Ta))*(np.pi*(De))

                TDi = Ti - qdp*(rt.rt_conv_cili(di,h_fld,z)+rt.rt_cond_cili(di,de,lmd_tube,z))
                #Te = TDi + qdp * (rt.rt_cond_cili(Di, De, flmd) + R_rev)
                Lslmd = Lslmd + [flmd]
                Lq = Lq + [qdp]
                LTe = LTe + [Te]
                LDe = LDe + [De]
                '''A divisão por z e os expoentes (-1) foram eliminados por enquanto
                   LR = LR + [ (rt.rt_conv_cili(di,h_fld,z)+ \
                            rt.rt_cond_cili(di,de,lmd_tube,z)+ \
                            rt.rt_cond_cili(Di,De,flmd,z)+ \
                            R_rev+ \
                            ((((rt.rt_conv_cili(De,ctc.hc_f_c(U,De,Te,Ta),z))**(-1))+((rt.rt_crad_cili(De,ctc.hr(eps,Te,Ta),z))**(-1)))**(-1))) ]'''
    
    #Diâmetros externos em milímetros.
    LDe_Disp = [1000*D for D in LDe]

    #Lista de soluções para a temperatura na face externa em °C.
    Lte = [(Te - 273.15) for Te in LTe]

    #Lista de q através de LMTD e R.
    LVT = list(map(lambda x: Ti - (Ta - (Ta - Ti)*np.exp(-1/(m*c*x))),LR))

    LTSK = list(map(lambda x: Ti - x, LVT))
    LLMTD = [((x - Ta) - (Ti - Ta))/(np.log((x - Ta)/(Ti - Ta))) for x in LTSK]
    ZTR = zip(LLMTD, LR)
    Lq = [x[0]/(x[1]) for x in ZTR]

    #Análise econômica das combinações Isolamento-Espessora
    #Custo de energia perdida atualizada.
    LCE = [CE(abs(q), N, CEE, eta, COP) for q in Lq]
    #Custo de investimento.

    #O investimento é o custo do material + custo mão de mao_de_obra
    #A lista LFCI2 vai pegar os dados do módulo (biblioteca) custos
    # custo_do_material * área lateral da tubulação + custo da mão de obra * comprimento da tubulação
    LFCI2 = []
    LFCI2 += [custos.aerogel_custos[esp] for esp in LE_Disp_Imp2]
    LFCI2 += [custos.espuma_custos[esp] for esp in LE_Disp_Imp2]
    #LCI = [np.nan]
    
    LCI = []
    for ci in LFCI2:
        teste = ci * 2 * np.pi * z * de/2 + custos.mao_de_obra_custos * z
        LCI += [teste]
    #Custo de manutenção.
    LCM = [CM(CI, tm) for CI in LCI]
    #Custo total.
    LCT = [x + y + z for (x,y,z) in zip(LCE, LCI, LCM)]
    
    #Isolamento Atual
    cp_media = (wtr.cp(TF) + wtr.cp(Ti))/2
    q_Atual = m * cp_media * (TF - Ti)
    CE_Atual = CE(abs(q_Atual), N, CEE, eta, COP)
    CM_Atual = custos.custo_manutencao
    CT_Atual = CM_Atual + CE_Atual
    
    #Para cada espessura dos materiais trabalhados, essa
    #lista armazena o quanto é economizado em relação
    #ao custo total do isolamento atual
    Dif_CTAtual_CTNovos = [CT_Atual - custo_LCT for custo_LCT in LCT] #Saving

    #Viabilidade de Projeto
    Valor_Pres_Liq = [] #VPL
    for index in range(0,32):
        Valor_Pres_Liq.append(-LCI[index] + vpl_revisado.VPL(Dif_CTAtual_CTNovos[index], n, i))

    #VPLAE
    Valor_Pres_LiqAE = [vpl * ((i*((1+i)**n))/(((1+n)**n)-1)) for vpl in Valor_Pres_Liq]

    payback = [] #Payback
    for index in range(0,32):
        payback.append(vpl_revisado.payback(-LCI[index], Dif_CTAtual_CTNovos[index], n, i))
    
    viabilidade_bool = [] #Viável?
    for index in range(0,32):
        if Valor_Pres_Liq[index] > 0 and Valor_Pres_LiqAE[index] > 0:
            viabilidade_bool.append("Sim")
        else:
            viabilidade_bool.append("Não")

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

   
    #Organização dos dados em um DataFrame.
    Disp = pd.DataFrame({'Material' : LNM,
                         'Espessura \n [mm]' : LE_Disp_DF,
                         'Espessura \n [pol]' : LE_Disp_Imp_DF,
                         'Diâmetro \n Externo [mm]' : LDe_Disp,
                         'Temperatura na \n Face Externa [°C]' : Lte,
                         'Resistência Térmica': LR})

    LVT = list(map(lambda x: Ti - (Ta - (Ta - Ti)*np.exp(-1/(m*c*x))),LR))
    LVT_Disp = [abs(T) for T in LVT]
    Disp['Variação de Temperatura \n do Fluido [°C]'] = LVT_Disp
    LTS = list(map(lambda x: Ti - 273.15 - x, LVT))
    Disp['Temperatura do Fluido \n na Saída [°C]'] = LTS
    Disp['Fluxo de Calor \n [W/m]'] = Lq
    Disp['Custo de Energia \n [$/m]'] = LCE
    Disp['Custo de Investimento \n [$/m]'] = LCI
    Disp['Custo de Manutenção \n [$/m]'] = LCM
    Disp['Custo Total \n [$/m]'] = LCT
    Disp['Saving'] = Dif_CTAtual_CTNovos
    Disp['VPL'] = Valor_Pres_Liq
    Disp['VPLAE'] = Valor_Pres_LiqAE
    Disp['Payback'] = payback
    Disp['Viável?'] = viabilidade_bool 
    #Disp = Disp.sort_values(by=['Custo Total \n [$/m]'])

    Disp = Disp.iloc[[0]].append(Disp.iloc[1:].sort_values(by=['VPL']))

    #Isolamento atual (DataFrame)
    Isolamento_Atual = {'Material' : 'Isolamento Atual',
                         'Espessura \n [mm]' : 38.1,
                         'Espessura \n [pol]' : 1.5, 
                         'Diâmetro \n Externo [mm]' : 304.8, 
                         'Temperatura na \n Face Externa [°C]' : np.NaN,
                         'Resistência Térmica': (Ta-Ti) / q_Atual}
    Isolamento_Atual['Variação de Temperatura \n do Fluido [°C]'] = TF-Ti
    Isolamento_Atual['Temperatura do Fluido \n na Saída [°C]'] = TF
    Isolamento_Atual['Fluxo de Calor \n [W/m]'] = q_Atual
    Isolamento_Atual['Custo de Energia \n [$/m]'] = CE_Atual
    Isolamento_Atual['Custo de Investimento \n [$/m]'] = 0
    Isolamento_Atual['Custo de Manutenção \n [$/m]'] = CM_Atual
    Isolamento_Atual['Custo Total \n [$/m]'] = CT_Atual
    Isolamento_Atual['Saving'] = 0
    Isolamento_Atual['VPL'] = 0
    Isolamento_Atual['VPLAE'] = 0
    Isolamento_Atual['Payback'] = 0 
    Isolamento_Atual['Viável?'] = np.NaN
    Isolamento_Atual_Lista = [Isolamento_Atual]
    Disp = pd.concat([pd.DataFrame(Isolamento_Atual_Lista), Disp], ignore_index=True)

    LA = ['Ti', 'Ponto de Orvalho', 'di', 'de', 'Ta', 'h_fld', 'lmd_tube', 'U', 'H', 'z', 'eps', 'm', 'c', 'Dt_max', 'R_rev', 'RH', 'N', 'CEE', 'eta', 'COP', 'n', 'i', 'tm', 'TF']
    LB = [Ti, wtr.TDP(Ta, RH), di, de, Ta, h_fld, lmd_tube, U, H, z, eps, m, c, Dt_max, R_rev, RH, N, CEE, eta, COP, n, i, tm, TF]
    while len(LA) < 33:
        LA.append(np.nan)
        LB.append(np.nan)
    Disp.insert(loc = 0, column = 'Input', value = LB)
    Disp.insert(loc = 0, column = 'Varáivel', value = LA)
#########################################################################################################3
    fesp_VPL(plt, LE_Disp_Imp2, Valor_Pres_Liq)
    print(Valor_Pres_Liq)
    """fesp_q(plt, LE_Disp_Imp2, Lq)
    fesp_payback(plt, LE_Disp_Imp2, payback)

    melhor_espessura_aerogel = get_greater_y(LE_Disp_Imp2, Valor_Pres_Liq, "Aerogel")
    melhor_espessura_espuma = get_greater_y(LE_Disp_Imp2, Valor_Pres_Liq, "Espuma")
    fct_aerogel = Dif_CTAtual_CTNovos[melhor_espessura_aerogel]
    custo_inicial_aerogel = -LCI[melhor_espessura_aerogel]
    fct_espuma = Dif_CTAtual_CTNovos[melhor_espessura_espuma]
    custo_inicial_espuma = -LCI[melhor_espessura_espuma]
    vpl_lista_aerogel = vpl_revisado.VPL_lista(custo_inicial_aerogel, fct_aerogel, n, i)
    vpl_lista_espuma = vpl_revisado.VPL_lista(custo_inicial_espuma, fct_espuma, n, i)
    vpl_lista_atual = vpl_revisado.VPL_lista(0, 0, n, i)
    fn_vpl(plt, vpl_lista_aerogel, vpl_lista_espuma, n)
    fn_vpl_atual(plt, vpl_lista_atual, n)
    
    juros = []
    vpl_aerogel = []
    vpl_espuma = []
    vpl_atual = []
    for x in np.linspace(4,15,1000):
        new_i = (x)/100
        juros.append(new_i)
        vpl_aerogel.append(vpl_revisado.VPL(fct_aerogel, n, new_i))
        vpl_espuma.append(vpl_revisado.VPL(fct_espuma, n, new_i))
        vpl_atual.append(vpl_revisado.VPL(0, n, new_i))
    fi_vpl(plt, vpl_aerogel, vpl_espuma, juros)
    fi_vpl_atual(plt, vpl_atual, juros)
    
    vpl_aerogel = []
    vpl_espuma = []
    vpl_atual = []
    temperaturas_aerogel = [] #variando de Te até Ta
    temperaturas_espuma = [] #variando de Te até Ta
    temperaturas_atual = []
    temperatura_aerogel = Lte[melhor_espessura_aerogel]
    temperatura_espuma = Lte[melhor_espessura_espuma]
    for t in np.linspace(temperatura_aerogel*100, Ta*100 + 1, 1000):
        temperaturas_aerogel.append(t/100)
    q_aerogel = []
    R_aerogel = 0
    for t in temperaturas_aerogel:
        R_aerogel = (rt.rt_conv_cili(di,h_fld)+\
                rt.rt_cond_cili(di,de,lmd_tube)+\
                rt.rt_cond_cili(Di,De,flmd)+\
                R_rev+((((rt.rt_conv_cili(De,ctc.hc_m_ch(U,De,t,Ta))))+((rt.rt_crad_cili(De,ctc.hr(eps,t,Ta)))))))
        q_aerogel.append((Ta - Ti) / R_aerogel)
        lce_aerogel = CE(abs(q_aerogel[-1]), N, CEE, eta, COP)
        lci_aerogel = LFCI2[melhor_espessura_aerogel] * 2 * np.pi * z * de/2 + custos.mao_de_obra_custos * z
        lcm_aerogel = CM(lci_aerogel, tm)
        lct_aerogel = lcm_aerogel + lci_aerogel + lce_aerogel
        fct_aerogel = CT_Atual - lct_aerogel
        vpl_aerogel.append(vpl_revisado.VPL(fct_aerogel, n, i))

    for t in np.linspace(temperatura_espuma*100, Ta*100 + 1, 1000):
        temperaturas_espuma.append(t/100)
    q_espuma = []
    R_espuma = 0
    for t in temperaturas_espuma:
        R_espuma = (rt.rt_conv_cili(di,h_fld)+\
                rt.rt_cond_cili(di,de,lmd_tube)+\
                rt.rt_cond_cili(Di,De,flmd)+\
                R_rev+((((rt.rt_conv_cili(De,ctc.hc_m_ch(U,De,t,Ta))))+((rt.rt_crad_cili(De,ctc.hr(eps,t,Ta)))))))
        q_espuma.append((Ta - Ti) / R_espuma)
        lce_espuma = CE(abs(q_espuma[-1]), N, CEE, eta, COP)
        lci_espuma = LFCI2[melhor_espessura_espuma] * 2 * np.pi * z * de/2 + custos.mao_de_obra_custos * z
        lcm_espuma = CM(lci_espuma, tm)
        lct_espuma = lcm_espuma + lci_espuma + lce_espuma
        fct_espuma = CT_Atual - lct_espuma
        vpl_espuma.append(vpl_revisado.VPL(fct_espuma, n, i))
    fTe_q(plt, temperaturas_aerogel, temperaturas_espuma, q_aerogel, q_espuma)
    fTe_VPL(plt, temperaturas_aerogel, temperaturas_espuma, vpl_aerogel, vpl_espuma)

    q_Atual = []
    for t in np.linspace(TF*100, Ta*100 + 1, 1000):
        temperaturas_atual.append(t/100)
    for t in temperaturas_atual:
        cp_media = (wtr.cp(t) + wtr.cp(Ti))/2
        q_Atual.append(m * cp_media * (Ta - Ti))
        CE_Atual = CE(abs(q_Atual[-1]), N, CEE, eta, COP)
        CM_Atual = custos.custo_manutencao
        CT_Atual = CM_Atual + CE_Atual
        vpl_atual.append(vpl_revisado.VPL(CT_Atual, n, i))
    fTe_q_atual(plt, temperaturas_atual, q_Atual)
    fTe_VPL_atual(plt, temperaturas_atual, vpl_atual)"""
    
    return Disp
