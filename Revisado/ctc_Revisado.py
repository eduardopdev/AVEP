# -*- coding: utf-8 -*-

import numpy as np

import ar

# =============================================================================
# Coeficientes de transferência de calor e fluxo convectivo. Referências [016]
# - N-550 - Projeto de Isolamento Térmico a Alta Temperatura, Anexo E; e [017]
# - Incropera - Fundamentals of Heat & Mass Transfer - 8th Edition.
# =============================================================================



# =============================================================================
# Funções auxiliares.
# =============================================================================

# =============================================================================
# Emissividade de superfícies. Input em string, output 1.
# =============================================================================

def emissividade(nome):

    """
    Preencha com um nome entre aspas.
    Opções:
    'Alumínio';
    'Tinta Preta';
    'Tinta de Alumínio';
    'Chapa de Aço';
    'Tinta Branca';
    'Massa Asfáltica'.
    """

################# REVISAR VALORES DE EMISSIVIDADES ##############
    emissividades = {'Alumínio': 0.2,
                     'Tinta Preta': 0.97,
                     'Tinta de Alumínio': 0.5,
                     'Chapa de Aço' : 0.955,
                     'Tinta Branca' : 0.88,
                     'Massa Asfáltica' : 0.93}

    epsilon = emissividades[nome]

    return(epsilon)

# =============================================================================
# Número de Reynolds de um escoamento. Input em m/s, m e m^2/s, output em 1.
# =============================================================================

def Re(U, D, Te, Ta):

    Tf = (Te + Ta)/2

    nu = ar.nu(Tf)

    Re = U*D/nu

    return(Re)


# =============================================================================
# Radiação.
# =============================================================================

# =============================================================================
# Coeficiente de transferência de calor por radiação. Input em 1 e K, output em
# W/(m^2.K).
# =============================================================================

def hr(epsilon, Te, Ta):

    sigma = 5.670374e-08

    T = ((Te)**2+(Ta)**2)*(Te + Ta)

    hr = sigma*epsilon*T

    return(hr)

# Fluxo de Calor por Radiação - Verificar Manutenção dessa função

def qr(epsilon, Te, Ta):
    
    h = hr(epsilon, Te, Ta)
    
    return(h*(Te-Ta))


 # Convecção forçada.
# =============================================================================

# =============================================================================
# Coeficiente de transferência de calor por convecção forçada em um cilindro.
# Input em m/s, m e K, output em W/(m^2.K).
# =============================================================================

def hc_f_c(U, D, Te, Ta):
    
    Tf = (Te + Ta)/2
    
    PrL = ar.Pr(Tf)
    
    ReL = Re(U, D, Te, Ta)
    
    a = (ReL/282000)**(5.0/8.0)
    
    b = (1 + a)**(4.0/5.0)
    
    c = 0.62*(ReL**0.5)*(PrL**(1.0/3.0))*b
    
    d = (0.4/PrL)**(2.0/3.0)
    
    e = (1 + d)**(0.25)
    
    f = c/e
    
    NuL = 0.3 + f
    
    hc = (NuL*(ar.lamed(Tf)))/D
    
    return (hc)

# Fluxo de calor por convecção externa forçada

def qc_f_c(U, D, Te, Ta):
    
    hc = hc_f_c(U, D, Te, Ta)
    
    return (hc*(Te-Ta))
