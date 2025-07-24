import matplotlib.pyplot as plt
import numpy as np


import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt( "declination vs ss.txt", skiprows = 1)
flux_norm = 1e-8 
declination =[i[0]  for i in data]
flux =[i[1]*flux_norm for i in data]
flux_2 = [i[2]*flux_norm*2 for i in data]

plt.figure(figsize=(8, 4.5), dpi = 180)
plt.scatter(declination, flux, label = r"$ \Phi_1 $")
plt.scatter(declination, flux_2, label = r"$ \Phi_2$")
plt.xlabel(r"$ \alpha $")
plt.ylabel(r"$ \Phi $ (sensitivity)")
plt.legend()
plt.show()
#neutrino telescope vuol dire che è indipendente quasi dal flusso che io injecto.Dipende quindi SOLO DAI TAGLI non tanto dal background

 


