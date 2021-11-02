
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
def get_greater_saving(espessuras, saving, material):
    """
    Retorna o índice da espessura associada
    ao maior valor de saving para o materia
    passado como parâmetro
    """
    if material=="Aerogel":
        index_espessura = 0
        max_saving = saving[0]
        for index in range(1,16):
            if saving[index] > max_saving:
                    max_saving = saving[index]
                    index_espessura = index

    elif material=="Espuma":
        index_espessura = 16
        max_saving = saving[16]
        for index in range(16,32):
            if saving[index] > max_saving:
                    max_saving = saving[index]
                    index_espessura = index
    return index_espessura

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

def fTe_q(pltObj, Ts, q_aerogel, q_espuma):
    pltObj.plot(Ts, q_aerogel, "b.-", label = "Aerogel")
    pltObj.plot(Ts, q_espuma, "r.-", label = "Espuma")
    pltObj.title('Temperatura (°C) x Calor')
    plt.xlabel('K')
    plt.ylabel('calor')
    plt.legend()
    plt.show()

def fTe_VPL(pltObj, Ts, VPL_aerogel, VPL_espuma):
    pltObj.plot(Ts, VPL_aerogel, "b.-", label = "Aerogel")
    pltObj.plot(Ts, VPL_espuma, "r.-", label = "Espuma")
    pltObj.title('Temperatura (°C) x VPL')
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
    """
    Plota um gráfico que o valor de vpl a cada ano
    para aerogel e espuma levando em consideração
    a espessura de cada material que gerou saving
    em cada material.
    """
    n_lista = list(range(n+1))
    pltObj.plot(n_lista, vpl_lista_aerogel,"b.-", label="Aerogel")
    pltObj.plot(n_lista, vpl_lista_espuma,"r.-", label="Espuma")
    pltObj.title('Tempo (anos) x VPL')
    pltObj.xlabel('tempo')
    pltObj.ylabel('VPL')
    pltObj.ylim(-15000, 30000)
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
    
    LCI = []
    for ci in LFCI2:
        custo_de_investimento = ci * 2 * np.pi * z * de/2 + custos.mao_de_obra_custos * z
        #exit()
        LCI += [custo_de_investimento]
    #Custo de manutenção.
    LCM = [CM(CI, tm) for CI in LCI]
    #Custo total.
    LCT = [x + y for (x,y) in zip(LCE, LCM)]
    
    #Isolamento Atual
    cp_media = (wtr.cp(TF) + wtr.cp(Ti))/2
    q_Atual = m * cp_media * (TF - Ti)
    CE_Atual = CE(abs(q_Atual), N, CEE, eta, COP)
    CM_Atual = custos.custo_manutencao
    CT_Atual = CM_Atual + CE_Atual
    
    #Para cada espessura dos materiais trabalhados, essa
    #lista armazena o quanto é economizado em relação
    #ao custo total do isolamento atual
    Saving = [CT_Atual - custo_LCT for custo_LCT in LCT] #Saving

    #Viabilidade de Projeto
    Valor_Pres_Liq = [] #VPL
    for index in range(0,32):
        Valor_Pres_Liq.append(-LCI[index] + vpl_revisado.VPL(Saving[index], n, i))

    #VPLAE
    Valor_Pres_LiqAE = [vpl * ((i*((1+i)**n))/(((1+n)**n)-1)) for vpl in Valor_Pres_Liq]

    payback = [] #Payback
    for index in range(0,32):
        payback.append(vpl_revisado.payback(-LCI[index], Saving[index], n, i))
    
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
    Disp['Saving'] = Saving
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
    fesp_q(plt, LE_Disp_Imp2, Lq)
    fesp_payback(plt, LE_Disp_Imp2, payback)

    maior_saving_aerogel = get_greater_saving(LE_Disp_Imp2, Saving, "Aerogel")
    maior_saving_espuma = get_greater_saving(LE_Disp_Imp2, Saving, "Espuma")
    fct_aerogel = Saving[maior_saving_aerogel]
    custo_inicial_aerogel = -LCI[maior_saving_aerogel]
    fct_espuma = Saving[maior_saving_espuma]
    custo_inicial_espuma = -LCI[maior_saving_espuma]
    vpl_lista_aerogel = vpl_revisado.VPL_lista(custo_inicial_aerogel, fct_aerogel, n, i)
    vpl_lista_espuma = vpl_revisado.VPL_lista(custo_inicial_espuma, fct_espuma, n, i)
    fn_vpl(plt, vpl_lista_aerogel, vpl_lista_espuma, n)
    

    fct_aerogel = Saving[5] 
    custo_inicial_aerogel = -LCI[5]
    fct_espuma = Saving[21]
    custo_inicial_espuma = -LCI[21]
    vpl_lista_aerogel = vpl_revisado.VPL_lista(custo_inicial_aerogel, fct_aerogel, n, i)
    vpl_lista_espuma = vpl_revisado.VPL_lista(custo_inicial_espuma, fct_espuma, n, i)
    fn_vpl(plt, vpl_lista_aerogel, vpl_lista_espuma, n)

    fct_aerogel = Saving[maior_saving_aerogel]
    fct_espuma = Saving[maior_saving_espuma]
    juros = []
    vpl_aerogel = []
    vpl_espuma = []
    for x in np.linspace(4,15,1000):
        new_i = (x)/100
        juros.append(new_i)
        vpl_aerogel.append(-LCI[maior_saving_aerogel] + vpl_revisado.VPL(fct_aerogel, n, new_i))
        vpl_espuma.append(-LCI[maior_saving_espuma] + vpl_revisado.VPL(fct_espuma, n, new_i))
    fi_vpl(plt, vpl_aerogel, vpl_espuma, juros)

    fct_aerogel = Saving[5]
    fct_espuma = Saving[21]
    juros = []
    vpl_aerogel = []
    vpl_espuma = []
    for x in np.linspace(4,15,1000):
        new_i = (x)/100
        juros.append(new_i)
        vpl_aerogel.append(-LCI[5] + vpl_revisado.VPL(fct_aerogel, n, new_i))
        vpl_espuma.append(-LCI[21] + vpl_revisado.VPL(fct_espuma, n, new_i))
    fi_vpl(plt, vpl_aerogel, vpl_espuma, juros)

    # Gráficos de Temperatura por calor e VPL considerando a espessura do melhor saving
    vpl_aerogel = []
    vpl_espuma = []
    temperaturas = []
    temperatura_inicial = int(Lte[maior_saving_aerogel]*100)
    Ta_ajustado = Ta - 273.15
    for t in np.linspace(temperatura_inicial, Ta_ajustado + 1, 16):
        temperaturas.append(t/100)
    q_aerogel = []
    R_aerogel = []
    De_aerogel = Di + 2*LE_Disp_Imp2[maior_saving_aerogel]
    for t in temperaturas:
        R_aerogel = R_aerogel + [(rt.rt_conv_cili(di,h_fld,z)+ \
                    rt.rt_cond_cili(di,de,lmd_tube,z)+ \
                    rt.rt_cond_cili(Di,De_aerogel,flmd,z)+ \
                    R_rev+ \
                    ((((rt.rt_conv_cili(De_aerogel,ctc.hc_f_c(U,De_aerogel,t,Ta),z))**(-1))+\
                            ((rt.rt_crad_cili(De_aerogel,ctc.hr(eps,t,Ta),z))**(-1)))**(-1)))]

    LVT_aerogel = list(map(lambda x: Ti - (Ta - (Ta - Ti)*np.exp(-1/(m*c*x))),R_aerogel))
    LTSK_aerogel = list(map(lambda x: Ti - x, LVT_aerogel))
    LLMTD_aerogel = [((x - Ta) - (Ti - Ta))/(np.log((x - Ta)/(Ti - Ta))) for x in LTSK_aerogel]
    ZTR_aerogel = zip(LLMTD_aerogel, R_aerogel)
    q_aerogel = [x[0]/(x[1]) for x in ZTR_aerogel]
    LCE_aerogel = [CE(abs(q), N, CEE, eta, COP) for q in q_aerogel]
    LFCI2_aerogel = [custos.aerogel_custos[esp] for esp in LE_Disp_Imp2]
    LCI_aerogel = []
    for ci in LFCI2_aerogel:
        custo_de_investimento = ci * 2 * np.pi * z * de/2 + custos.mao_de_obra_custos * z
        LCI_aerogel += [custo_de_investimento]
    LCM_aerogel = [CM(CI, tm) for CI in LCI_aerogel]
    LCT_aerogel = [x + y for (x,y) in zip(LCE_aerogel, LCM_aerogel)]
    Saving_aerogel = [CT_Atual - custo_LCT for custo_LCT in LCT_aerogel]
    vpl_aerogel = []
    for index in range(0,16):
        vpl_aerogel.append(-LCI_aerogel[index] + vpl_revisado.VPL(Saving_aerogel[index], n, i))        

    q_espuma = []
    R_espuma = []
    De_espuma = Di + 2*LE_Disp_Imp2[maior_saving_espuma%16]
    for t in temperaturas:
        R_espuma = R_espuma + [(rt.rt_conv_cili(di,h_fld,z)+ \
                    rt.rt_cond_cili(di,de,lmd_tube,z)+ \
                    rt.rt_cond_cili(Di,De_espuma,flmd,z)+ \
                    R_rev+ \
                    ((((rt.rt_conv_cili(De_espuma,ctc.hc_f_c(U,De_espuma,t,Ta),z))**(-1))+\
                            ((rt.rt_crad_cili(De_espuma,ctc.hr(eps,t,Ta),z))**(-1)))**(-1)))]
    LVT_espuma = list(map(lambda x: Ti - (Ta - (Ta - Ti)*np.exp(-1/(m*c*x))),R_espuma))
    LTSK_espuma = list(map(lambda x: Ti - x, LVT_espuma))
    LLMTD_espuma = [((x - Ta) - (Ti - Ta))/(np.log((x - Ta)/(Ti - Ta))) for x in LTSK_espuma]
    ZTR_espuma = zip(LLMTD_espuma, R_espuma)
    q_espuma = [x[0]/(x[1]) for x in ZTR_espuma]

    LCE_espuma = [CE(abs(q), N, CEE, eta, COP) for q in q_espuma]
    LFCI2_espuma = [custos.espuma_custos[esp] for esp in LE_Disp_Imp2]
    LCI_espuma = []
    for ci in LFCI2_espuma:
        custo_de_investimento = ci * 2 * np.pi * z * de/2 + custos.mao_de_obra_custos * z
        LCI_espuma += [custo_de_investimento]
    LCM_espuma = [CM(CI, tm) for CI in LCI_espuma]
    LCT_espuma = [x + y for (x,y) in zip(LCE_espuma, LCM_espuma)]
    Saving_espuma = [CT_Atual - custo_LCT for custo_LCT in LCT_espuma]
    vpl_espuma = []
    for index in range(0,16):
        vpl_espuma.append(-LCI_espuma[index] + vpl_revisado.VPL(Saving_espuma[index], n, i))
    
    fTe_q(plt, temperaturas , q_aerogel, q_espuma)
    fTe_VPL(plt, temperaturas , vpl_aerogel, vpl_espuma)

    # Gráficos de Temperatura por calor e VPL considerando a espessura do isolamento atual
    vpl_aerogel = []
    vpl_espuma = []
    temperaturas = []

    for t in np.linspace(temperatura_inicial, Ta_ajustado + 1, 16):
        temperaturas.append(t/100)
    q_aerogel = []
    R_aerogel = []
    De_aerogel = Di + 2*LE_Disp_Imp2[5]
    for t in temperaturas:
        R_aerogel = R_aerogel + [(rt.rt_conv_cili(di,h_fld,z)+ \
                    rt.rt_cond_cili(di,de,lmd_tube,z)+ \
                    rt.rt_cond_cili(Di,De_aerogel,flmd,z)+ \
                    R_rev+ \
                    ((((rt.rt_conv_cili(De_aerogel,ctc.hc_f_c(U,De_aerogel,t,Ta),z))**(-1))+\
                            ((rt.rt_crad_cili(De_aerogel,ctc.hr(eps,t,Ta),z))**(-1)))**(-1)))]
        LVT_aerogel = list(map(lambda x: Ti - (Ta - (Ta - Ti)*np.exp(-1/(m*c*x))),R_aerogel))
        LTSK_aerogel = list(map(lambda x: Ti - x, LVT))
        LLMTD_aerogel = [((x - Ta) - (Ti - Ta))/(np.log((x - Ta)/(Ti - Ta))) for x in LTSK_aerogel]
        ZTR_aerogel = zip(LLMTD_aerogel, R_aerogel)
        q_aerogel = [x[0]/(x[1]) for x in ZTR_aerogel]
        LCE_aerogel = [CE(abs(q), N, CEE, eta, COP) for q in q_aerogel]
        LFCI2_aerogel = [custos.aerogel_custos[esp] for esp in LE_Disp_Imp2]
        LCI_aerogel = []
        for ci in LFCI2_aerogel:
            custo_de_investimento = ci * 2 * np.pi * z * de/2 + custos.mao_de_obra_custos * z
            LCI_aerogel += [custo_de_investimento]
        LCM_aerogel = [CM(CI, tm) for CI in LCI_aerogel]
        LCT_aerogel = [x + y for (x,y) in zip(LCE_aerogel, LCM_aerogel)]
        Saving_aerogel = [CT_Atual - custo_LCT for custo_LCT in LCT_aerogel]
        vpl_aerogel = []
        for index in range(0,16):
            vpl_aerogel.append(-LCI_aerogel[index] + vpl_revisado.VPL(Saving[index], n, i))        

    q_espuma = []
    R_espuma = []
    De_espuma = Di + 2*LE_Disp_Imp2[5]
    for t in temperaturas:
        R_espuma = R_espuma + [(rt.rt_conv_cili(di,h_fld,z)+ \
                    rt.rt_cond_cili(di,de,lmd_tube,z)+ \
                    rt.rt_cond_cili(Di,De_espuma,flmd,z)+ \
                    R_rev+ \
                    ((((rt.rt_conv_cili(De_espuma,ctc.hc_f_c(U,De_espuma,t,Ta),z))**(-1))+\
                            ((rt.rt_crad_cili(De_espuma,ctc.hr(eps,t,Ta),z))**(-1)))**(-1)))]
        LVT_espuma = list(map(lambda x: Ti - (Ta - (Ta - Ti)*np.exp(-1/(m*c*x))),R_espuma))
        LTSK_espuma = list(map(lambda x: Ti - x, LVT))
        LLMTD_espuma = [((x - Ta) - (Ti - Ta))/(np.log((x - Ta)/(Ti - Ta))) for x in LTSK_espuma]
        ZTR_espuma = zip(LLMTD_espuma, R_espuma)
        q_espuma = [x[0]/(x[1]) for x in ZTR_espuma]
        LCE_espuma = [CE(abs(q), N, CEE, eta, COP) for q in q_espuma]
        LFCI2_espuma = [custos.espuma_custos[esp] for esp in LE_Disp_Imp2]
        LCI_espuma = []
        for ci in LFCI2_espuma:
            custo_de_investimento = ci * 2 * np.pi * z * de/2 + custos.mao_de_obra_custos * z
            LCI_espuma += [custo_de_investimento]
        LCM_espuma = [CM(CI, tm) for CI in LCI_espuma]
        LCT_espuma = [x + y for (x,y) in zip(LCE_espuma, LCM_espuma)]
        Saving_espuma = [CT_Atual - custo_LCT for custo_LCT in LCT_espuma]
        vpl_espuma = []
        for index in range(0,16):
            vpl_espuma.append(-LCI_espuma[index] + vpl_revisado.VPL(Saving[index], n, i))
    fTe_q(plt, temperaturas , q_aerogel, q_espuma)
    fTe_VPL(plt, temperaturas , vpl_aerogel, vpl_espuma)
    
    return Disp
