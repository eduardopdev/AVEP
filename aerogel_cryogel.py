# -*- coding: utf-8 -*-

import numpy as np

# =============================================================================
# Propriedades de tubos flexíveis de aerogel.
# Referência [018] - NIA Insulation Materials Spec Chart - ASTM C1728.
# =============================================================================

def lin(t, t2, t1, l2, l1):

    a = l1

    b = l2 - l1

    x = (t - t1)/(t2 - t1)

    return(a + b*x)

# =============================================================================
# Termoisolantes Tipo I, Grau 1B. Temperatura máxima de utilização de 125 °C.
# Intervalo de temperatura média para o cálculo da condutividade térmica de
# -73 °C a 93 °C.
# =============================================================================

# =============================================================================
# Condutividade térmica. Input em K, output em W/(m.K).
# =============================================================================
aerogel_condutividade_term = 0.017 # W/m.K 
# =============================================================================
# Intervalo de temperatura média para o cálculo da condutividade térmica.
# Output em K.
# =============================================================================

def tmdI():

    tmd = np.array([-73, 93])

    tmd = tmd + 273.15

    return(tmd)

# =============================================================================
# Temperatura máxima de utilização. Output em K.
# =============================================================================

def tmxI():

    tmx = 125 + 273.15

    return(tmx)
