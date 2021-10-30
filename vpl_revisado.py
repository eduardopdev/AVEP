# -*- coding: utf-8 -*-

import numpy as np

# VPL
# somatorio de: t =1 até t=n, sendo n = vida util do isolamento em anos
def VPL(FCt, n, i):
    vpl = 0
    for t in range(1, n+1):
        #if t == 1:
        #   vpl = vpl + x (x < 0)
        vpl = vpl + FCt/(1+i)**t
    return vpl

def VPL_lista(custo_inicial_aerogel, FCt, n, i):
    vpl = custo_inicial_aerogel
    vpl_lista = [vpl]
    for t in range(1, n+1):
        #if t == 1:
        #   vpl = vpl + x (x < 0)
        vpl = vpl + FCt/(1+i)**t
        vpl_lista.append(vpl)
    return vpl_lista

#FCt = LCT - LCR  #LCT é o custo total do isolamento bom e LCR custo total com isolamento ruim/atual

# VPLAE (Valor presente líquido anualizado equivalente) - dividir em partes
def VPLAE(FCt, n, i):
    vpl = VPL(FCt, n, i)
    vplae = vpl*((i*((1+i)**n))/(((1+n)**n)-1))
    return vplae

#PAYBACK
def payback(lci, FCt, n, i): # retorna o t necessesário para
    vpl = lci #que o investimento passe a gerar
    for t in range(1, n+1):          #retorno do capital investido.
        vpl = vpl + FCt/(1+i)**t
        if vpl >= 0:
            return t
        else:
            continue
#Qual valor de t para VPL=0

# CRIAR LISTA (VIB) para cada material e espessura
#CRIAR 4 LISTAS - Fct, VPL, VPLAE, PAYBACK
