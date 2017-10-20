##### PLOTEAR CATALOGO #####

from __future__ import division
import scipy as sp
import catalog_reader
import re
import matplotlib.pyplot as plt

path1 = "../../../../Naxo/Universidad/Magister/Semestre I/Taller de Investigacion III/SED-fitting/example/Photometry-run/"
path2 = "../../../../Naxo/Universidad/Magister/Semestre I/Taller de Investigacion III/SED-fitting/example/Filters/"

translator = catalog_reader.read_translate_file(path1 + 'goodss.translate')
#catalog = catalog_reader.read_catalog(path1 + 'goodss.cat')
catalog = catalog_reader.read_catalog(path1 + 'goodss_dani.cat')
n_gal = len(catalog['id'])
ind_filters = []
pattern_F = re.compile(r'F{1}[0-9]+$')
pattern_E = re.compile(r'E{1}[0-9]+$')
pattern_TOT = re.compile(r'T{1}O{1}T{1}[0-9]+$')

for ind, val in enumerate(list(catalog.columns)):
   try:
       if pattern_F.match(translator[val]):
           ind_filters.append(int(translator[val][1:]))
   except: continue
   
ind_filters = sp.array(ind_filters)
n_filt = len(ind_filters)

filters, lamb = catalog_reader.read_filters(path2 + 'FILTER.RES.latest', ind_filters)

fnu = []
efnu = []
for ind, val in enumerate(list(catalog.columns)):
    try:
        if pattern_F.match(translator[val]):
          fnu.append(catalog[val])
        elif pattern_E.match(translator[val]):
          efnu.append(catalog[val])
        print "No Falla"
    except: print "Falla"
fnu = sp.array(fnu)
efnu = sp.array(efnu)

cat = 0

# Obtener unidades flujo en catalogo en CGS
'''
m_AB_zp_25_flux = -2.5*sp.log10(fnu) + 25
m_AB_zp_25_eflux = -2.5*sp.log10(efnu) + 25
flux_cgs = 10 ** ((m_AB_zp_25_flux + 48.6) / (-2.5))  # erg s**-1 cm**-2 Hz**-1
eflux_cgs = 10 ** ((m_AB_zp_25_eflux + 48.6) / (-2.5))  # erg s**-1 cm**-2 Hz**-1
fnu = flux_cgs
efnu = eflux_cgs
'''
# Formula para unidades x Jy -> 10**-x <=> Z0 = 2.5 * (23 + x) - 48.6
# Entonces x = 6.44 para unidades con Z0 = 25 <=> indades son 10**(-x)Jy = 10**(-23) erg s**-1 cm**-2 Hz**-1 (unidades de fnu)
fnu = fnu * 10**(-6.44) * 10**(-23)
efnu = efnu  * 10**(-6.44) * 10**(-23)  # Al final parece q es lo mismo q hacer lo de arriba

# Convert fnu into flambda

mat_lamb = sp.ones([n_gal, n_filt])  # Matrix with lambda^2 values
mat_lamb[:] = (lamb * 1e-8)**2  # Lambda in cm
mat_lamb = sp.transpose(mat_lamb)
eflux = 3.0e10 * efnu / mat_lamb  # c en cgs es 3e10
flux = 3.0e10 * fnu / mat_lamb
rep = sp.where(fnu == -99)
n_rep = len(rep[0])
if n_rep >= 1:
  eflux[rep] = -99
  flux[rep] = -99
fnu = 0
efnu = 0


data = {"flux":flux, "lambda":lamb, "eflux":eflux}

num_gal = 1

### SOLO PARA CHI2 ###
f_obs = flux[:, num_gal-1]  / 3.828e33
ef_obs = eflux[:, num_gal-1]  / 3.828e33
#print len(f_obs)
######################

# Eliminar ratios error_flujo/flujo muy altos

ind_ratio_sobre_3 = sp.where((data["eflux"][:,num_gal-1]/data["flux"][:,num_gal-1] > 3) |
                             (sp.isnan(data["eflux"][:,num_gal-1]/data["flux"][:,num_gal-1])))

lamb = sp.delete(data["lambda"], ind_ratio_sobre_3)
flux = sp.delete(data["flux"][:,num_gal-1], ind_ratio_sobre_3)
eflux = sp.delete(data["eflux"][:,num_gal-1], ind_ratio_sobre_3)


# Cgs to L_sol/AA

# L_sol = 3.828e33 erg s**-1
# AA = 1e-8 cm

flux = flux / 3.828e33
eflux = eflux / 3.828e33

# Para plotear con lineas y no puntos

ind_sort = sp.argsort(lamb)
flux = flux[ind_sort]
eflux = eflux[ind_sort]
lamb = lamb[ind_sort]

fig1 = plt.figure(1)
fig1.clf()
plt.plot(lamb, flux, 'r*')
plt.errorbar(lamb, flux, yerr=eflux)
plt.xlabel(r"$\lambda \ [\AA]$")
plt.ylabel(r"$f_{\lambda} \ [L_{\odot}/\AA]$")
plt.title("Catalogue's galaxy SED")
plt.xlim(0,10000)
plt.tight_layout()


# HACER BARRAS DE ERROR DE LOS PUNTOS Y ELIMINAR AQUELLOS CON RATIO EFLUJO/FLUJO > 3, ADEMAS
# USAR LA FORMULA DE LA FOTO PARA OBTENER UNIDADES DEL FLUJO

# REVISAR UNIDADES DE FLUJO PQ ME DAN DEMASIADO PEQUENAS EN L_sol/AA

# Para arreglar las unidades tengo q minimizar la distancia en chi cuadrado para poder
# obtener el factor de escala que hace q se acerque el modelo a lo observado



##### PLOTEAR MODELOS #####

import scipy as sp
import matplotlib.pyplot as plt
from ised import read_bc_ised
import catalog_reader
import dust
import madau
import re
from scipy.interpolate import interp1d, interp2d
from scipy.integrate import simps

models_path = "../../../../Naxo/Universidad/Magister/Semestre I/Taller de Investigacion III/SED-fitting/example/Libraries/ised_exp.pr/"
model_name = "bc03_pr_ch_z004_ltau10.9.ised"
path1 = "../../../../Naxo/Universidad/Magister/Semestre I/Taller de Investigacion III/SED-fitting/example/Photometry-run/"
path2 = "../../../../Naxo/Universidad/Magister/Semestre I/Taller de Investigacion III/SED-fitting/example/Filters/"

n_gal = 1
z_peak = 3.815

# Read model from ised file
path = models_path + model_name 
m = read_bc_ised(path, post2012=False)
#print sp.shape(m['wavelengths']), sp.shape(m['seds'])

wv_model = m['wavelengths']
f_lambda_model = m['seds']

# Reducir seds por polvo
f_lambda_model_dust = dust.add_dust(f_lambda_model, wv_model, [1.0, 1.3, 1.5], 'calzetti')

# Reducir seds por madau
f_lambda_model_madau = madau.add_madau(f_lambda_model, wv_model, z_peak)

# Reducir seds por polvo y madau
f_lambda_model_madau_dust = madau.add_madau(f_lambda_model_dust, wv_model, z_peak)


# Agregar influencia del redshift a espectro (correr espectro a la derecha)
def lamb_redshift(z, lamb_emit):
    lamb_obs = (1+z) * lamb_emit  # Lamb emit es el del modelo y lamb obs el del catalogo en si
                                  # ie, el del modelo corrido al rojo
    return lamb_obs
#wv_model = wv_model + 3400  # Correrlo a mano ya que el primer filtro esta como a eso
                            # Tambien podria implementar aqui q sume el lamb del primero filtro
wv_model = lamb_redshift(z_peak, wv_model)  # z_peak=3.815 del cat de la dani

# Read filters from catalogue
translator = catalog_reader.read_translate_file(path1 + 'goodss.translate')
#catalog = catalog_reader.read_catalog(path1 + 'goodss.cat')
catalog = catalog_reader.read_catalog(path1 + 'goodss_dani.cat')
num_gal = len(catalog['id'])
ind_filters = []
pattern_F = re.compile(r'F{1}[0-9]+$')
pattern_E = re.compile(r'E{1}[0-9]+$')
pattern_TOT = re.compile(r'T{1}O{1}T{1}[0-9]+$')

for ind, val in enumerate(list(catalog.columns)):
   try:
       if pattern_F.match(translator[val]):
           ind_filters.append(int(translator[val][1:]))
   except: continue
   
ind_filters = sp.array(ind_filters)
n_filt = len(ind_filters)
filters, lamb = catalog_reader.read_filters(path2 + 'FILTER.RES.latest', ind_filters)

# Number of filter
num = filters[:,0]
# Wavelength
wv = filters[:,1]
# Transmission
tr = filters[:,2]

# Interpolar longitudes de onda continuas (SED) en discretas (catalogo)
tr_new = []
wv_new = []
f_lambda_new = []
for i in range(n_filt):
    ind = sp.where(num==i)
    f_ip_tr = interp1d(wv[ind], tr[ind])
    f_lambda_ip = f_lambda_model[n_gal][sp.where((wv_model>=sp.sort(wv[ind])[0]) & (wv_model<=sp.sort(wv[ind])[-1]))]
    wv_ip = wv_model[sp.where((wv_model>=sp.sort(wv[ind])[0]) & (wv_model<=sp.sort(wv[ind])[-1]))]
    tr_ip = f_ip_tr(wv_ip)
    tr_new.append(tr_ip)
    wv_new.append(wv_ip)
    f_lambda_new.append(f_lambda_ip)
tr_new = sp.array(tr_new)
wv_new = sp.array(wv_new)
f_lambda_new = sp.array(f_lambda_new)

def filtros(f_lambda_modelo, tr_lambda_modelo, lambda_modelo):
    def integral(lamb):
        h = 6.62e-27 / 3.828e33  # En AA
        c = 3e10 * 1e8  # En AA
        integral = f_lambda_modelo * tr_lambda_modelo * lamb / (h * c)
        return integral
    I1 = simps(integral(lambda_modelo), lambda_modelo)
    I2 = sp.sum(tr_lambda_modelo * lambda_modelo)
    return I1 / I2  # Retorna flujo del filtro en unidades de flujo de fotones

filtros_int = []
for ind, val in enumerate(tr_new):
    filtros_var = filtros(f_lambda_new[ind], tr_new[ind], wv_new[ind])
    # Supongo que debo pasar denuevo flujo de fotones a flujo de energia
    h = 6.62e-27 / 3.828e33  # En AA
    c = 3e10 * 1e8  # En AA
    filtros_var = filtros_var * (h * c / lamb[ind])
    filtros_int.append(filtros_var)
filtros_int = sp.array(filtros_int)

### SOLO PARA CHI2 ###    
f_model = filtros_int
#print len(filtros_int)
######################
    
# Plots (podria plotear la masa estelar mstar o sfr o mas cosas)
fig2 = plt.figure(2)
fig2.clf()
plt.plot(m['wavelengths'], m['seds'][n_gal], 'b-', label="Without redshift")
plt.plot(wv_model,f_lambda_model[n_gal], 'r-', label="With redshift")
#plt.plot(wv_new[0], (f_lambda_new * tr_new)[0],'g-', label="Flux per filters")
#for i in range(len(wv_new)):
#    plt.plot(wv_new[i], (f_lambda_new * tr_new)[i],'g-')  # Revisar si es asi o tb tengo q convertir a fotones
plt.xlabel(r"$\lambda \ [\AA]$")
plt.ylabel(r"$f_{\lambda} \ [L_{\odot}/\AA]$")
plt.title("SED model " + model_name)
plt.legend(loc='upper right')
plt.xlim(0,10000)
plt.tight_layout()

# Plots flujos integrados en filtros
fig3 = plt.figure(3)
fig3.clf()
ind_sort = sp.argsort(lamb)
plt.plot(lamb[ind_sort], filtros_int[ind_sort], 'g*')
plt.plot(lamb[ind_sort], filtros_int[ind_sort], 'g-')
plt.xlabel(r"$\lambda \ [\AA]$")
plt.ylabel(r"$f_{\lambda} \ [L_{\odot}/\AA]$")
plt.title("SED model " + model_name + "\n integrated flux in catalogue filters")
plt.xlim(0,10000)
plt.tight_layout()


# Hacer lo mismo con las unidades que en flambda_lambda_plotter.py


##### PLOT DATOS Y MODELO #####

# Eliminar nans y ratios mayores a 3 para grafico con plot datos y modelo
ind_nan = sp.where((ef_obs/f_obs > 3) | (sp.isnan(f_obs)) | (sp.isnan(ef_obs/f_obs)))
f_obs = sp.delete(f_obs, ind_nan)
f_model = sp.delete(f_model, ind_nan)
ef_obs = sp.delete(ef_obs, ind_nan)
lamb = sp.delete(lamb, ind_nan)
filtros_int = sp.delete(filtros_int, ind_nan)

def A_chi2(obs, err_obs, model):
    def chi2(A):
        comps_chi2 = []
        for ind, val in enumerate(obs):
            comp_i_chi2 = ((obs[ind] - A * model[ind])**2) / err_obs[ind]**2
            comps_chi2.append(comp_i_chi2)
        chi2 = sp.sum(sp.array(comps_chi2))
        return chi2
    p0 = obs[0] / model[0]
    from scipy.optimize import fmin
    return fmin(chi2, p0)[0]
A = A_chi2(f_obs, ef_obs, f_model)
print A

# Plots (podria plotear la masa estelar mstar o sfr o mas cosas)
print sp.shape(lamb), sp.shape(f_obs)

ind_sort = sp.argsort(lamb)

fig4 = plt.figure(4)
fig4.clf()
plt.plot(lamb[ind_sort], f_obs[ind_sort], 'r-', label="SED Catalogue")
plt.plot(lamb[ind_sort], f_obs[ind_sort], 'r*')
plt.plot(lamb[ind_sort], A * filtros_int[ind_sort], 'g*')
plt.plot(lamb[ind_sort], A * filtros_int[ind_sort], 'g-', label="SED Model")
plt.legend(loc='upper right')
plt.xlabel(r"$\lambda \ [\AA]$")
plt.ylabel(r"$f_{\lambda} \ [L_{\odot}/\AA]$")
plt.title(r"SED fitting (only scaling)")
plt.xlim(0,10000)
plt.tight_layout()

# Plot con polvo y madau y sin redshift
fig5 = plt.figure(5)
fig5.clf()
plt.plot(m['wavelengths'], m['seds'][n_gal], 'b-', label="Model")
plt.plot(m['wavelengths'],f_lambda_model_dust[0, n_gal], 'r-', label="Dust")
plt.plot(m['wavelengths'],f_lambda_model_madau[0, n_gal], 'g-', label="Madau")
plt.plot(m['wavelengths'],f_lambda_model_madau_dust[0, n_gal], 'c-', label="Both")
plt.legend(loc='upper right')
plt.xlabel(r"$\lambda \ [\AA]$")
plt.ylabel(r"$f_{\lambda} \ [L_{\odot}/\AA]$")
plt.title("SED model " + model_name)
plt.xlim(0,10000)
plt.tight_layout()
plt.show('all')
plt.close('all')

# ASUMIR Q TRANSMISION ES EN FOTONES NO EN ENERGIA ASI Q TENGO Q CONVERTIR FLAMBDA EN FLAMBDA DE FOTONES (DIVIDIR POR H/CLMABDA)
# TENGO QUE CORRER EL ESPECTRO SEGUN EL REDSHIFT FOTOMETRICO (LAMBD OBS - LAMBD EMITIDO BLA BLA)