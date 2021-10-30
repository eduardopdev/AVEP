# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Funções de analogia de circuitos elétricos. Referências [016] - N-550 -
# Projeto de Isolamento Térmico a Alta Temperatura, Anexo E; e [017] -
# Incropera - Fundamentals of Heat & Mass Transfer - 8th Edition.
# =============================================================================



# =============================================================================
# Resistências térmicas em cilindros.
# =============================================================================

# =============================================================================
# Resistências térmicas de condução  em um cilindro.
# Input em m e W/(m.K), output em K/W.
# =============================================================================

def rt_cond_cili(Di, De, lmd, z):
    
    a = np.log(De/Di)
    
    b = 2.0*np.pi*z*lmd
    
    return (a/b)

# =============================================================================
# Resistência térmica de convecção em um cilindro. Input em
# m e W/(m^2.K), output em K/W.
# =============================================================================

def rt_conv_cili(D, h, z):
    
    return (1.0/(h*np.pi*D*z))

# =============================================================================
# Resistência térmica de radiação em um cilindro. Input em
# m e W/(m^2.K), output em K/W.
# =============================================================================

def rt_crad_cili(D, h, z):
    
    return (1.0/(h*np.pi*D*z))