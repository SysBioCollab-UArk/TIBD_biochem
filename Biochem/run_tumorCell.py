import numpy as np
from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as mp
import matplotlib.pyplot as plt
from pysb.util import alias_model_components
from tumorCell import model

print(model.monomers)

# 1D11 binds TGFb
Monomer('x1D11', ['tgfb'])
alias_model_components()
print(model.monomers)
Initial(x1D11(tgfb=None), Parameter('x1D11_init', 0))
# Observable('x1D11_TGFB', TGFB(r=1) % x1D11(tgfb=1))
Parameter('kf_1d11_tgfb_bind', 1000)
Parameter('kr_1d11_tgfb_bind', 10)
alias_model_components()
Rule('x1D11_binds_TGFB', TGFB(r=None) + x1D11(tgfb=None) | TGFB(r=1) % x1D11(tgfb=1), kf_1d11_tgfb_bind,
     kr_1d11_tgfb_bind)
Observable('mPTHrP_Obs', mPTHrP())

print(model.observables)

tspan = np.linspace(0, 72, 721)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
output = sim.run()

plt.plot(tspan, output.observables['mPTHrP_Obs'], lw=3, label='mPTHrP')
plt.xlabel('time (h)')
plt.ylabel('PTHrP mRNA')
plt.legend(loc=0)
plt.tight_layout()

plt.show()
