# Título: Análisis de Series de Tiempo: Aprobación Presidencial y Economía en Uruguay y Venezuela (2000-2015)

Descripción: Este proyecto analiza la relación dinámica entre la aprobación presidencial neta y los indicadores económicos (crecimiento del PIB y desempleo) en Uruguay y Venezuela, utilizando datos trimestrales del paquete politicalds.

El objetivo es determinar el impacto de la economía en la popularidad presidencial a través de diversos modelos de series de tiempo, evaluando la estacionariedad de las variables y corrigiendo la autocorrelación en los errores.

## Metodología:

Preprocesamiento: Transformación de fechas y creación de variables con rezagos (lags) y primeras diferencias.

Diagnóstico: Pruebas de Dickey-Fuller Aumentada (ADF) para detección de raíces unitarias.

Modelado Econométrico: Estimación y comparación de modelos (Estático, Ajuste Parcial, FDL, ECM, ARDL, Curva de Phillips).

Validación: Prueba de Durbin-Watson para detectar autocorrelación serial.

Visualización y Estética: Para la representación gráfica se utilizó ggplot2 junto con la librería de paletas de colores artísticas MoMA Colors (scale_color_moma), utilizando esquemas inspirados en artistas como Kippenberger, Dalí, Palermo, Avedon, Koons, entre otros, para diferenciar claramente las series temporales de cada país.

## Herramientas: R, Tidyverse, Car, Tseries, MoMAColors.

## -------------------------------------------------------------------------
## Proyecto: Análisis de Series de Tiempo: Aprobación Presidencial y Economía en Uruguay y Venezuela (2000-2015)
## Autor: José César Romero Galván
## Descripción: Análisis comparado de series de tiempo para Uruguay y Venezuela.
## -------------------------------------------------------------------------

# 1. CONFIGURACIÓN DEL ENTORNO ---------------------------------------------

# Definición del directorio de trabajo (Ajustalo a tu propia ruta)
# setwd("/Ruta/A/Tu/Directorio") 

# Configuración regional y notación
Sys.setlocale("LC_ALL", "es_ES.UTF-8")
options(scipen = 999) # Desactivar notación científica

# Carga e instalación de paquetería ----------------------------------------
if (!require("pacman")) install.packages("pacman")

# Verificación e instalación de MoMAColors desde GitHub
# (Esta librería contiene las paletas artísticas y no se encuentra en CRAN)
if (!require("MoMAColors")) {
  if (!require("devtools")) install.packages("devtools")
  devtools::install_github("BlakeRMills/MoMAColors")
}

pacman::p_load(car, politicalds, tidyverse, zoo, tseries, tsibble, MoMAColors)

# 2. PREPARACIÓN DE DATOS --------------------------------------------------

# Carga y transformación inicial
datos <- approval %>% 
  mutate(fecha = paste0(year, "-", quarter)) %>% 
  mutate(fecha = as.yearqtr(fecha)) # Conversión a formato año-trimestre

# Creación de subsets por país
Uruguay <- datos %>% filter(country == "Uruguay")
Venezuela <- datos %>% filter(country == "Venezuela")

# 3. ANÁLISIS EXPLORATORIO (EDA) -------------------------------------------

# Función auxiliar para graficar series de tiempo
plot_series <- function(data, y_var, title, color_palette) {
  data %>% 
    ggplot(aes(x = fecha, y = .data[[y_var]], color = country, group = country)) +
    geom_line() +
    geom_point() +
    scale_color_moma_d(color_palette) + 
    theme_minimal() +
    labs(
      title = title,
      subtitle = "2000 - 2015",
      x = "Fecha",
      y = "Valor",
      caption = "Fuente: Politicalds"
    ) +
    theme(plot.caption = element_text(hjust = 0))
}

# Visualización: Aprobación Neta (Niveles)
# Análisis: Presencia de tendencias y posible raíz unitaria en ambos casos.
plot_series(Uruguay, "net_approval", "Aprobación Presidencial Neta - Uruguay", "Kippenberger")
plot_series(Venezuela, "net_approval", "Aprobación Presidencial Neta - Venezuela", "Dali")

# Visualización: Crecimiento del PIB
# Análisis: Posible cambio estructural 2003-2005. Parece estacionaria visualmente.
plot_series(Uruguay, "gdp_growth", "Crecimiento del PIB - Uruguay", "Palermo")
plot_series(Venezuela, "gdp_growth", "Crecimiento del PIB - Venezuela", "Avedon")

# Visualización: Desempleo
# Análisis: Tendencia a la baja clara, sugiere no estacionariedad.
plot_series(Uruguay, "unemployment", "Desempleo - Uruguay", "Lupi")
plot_series(Venezuela, "unemployment", "Desempleo - Venezuela", "Koons")

# 4. PRUEBAS DE ESTACIONARIEDAD (ADF) --------------------------------------

# Variables en niveles
# H0: La serie tiene raíz unitaria (no es estacionaria)
cat("--- Pruebas ADF (Niveles) ---\n")
print(adf.test(Uruguay$net_approval, k = 1)) # p > 0.05, No estacionaria
print(adf.test(Uruguay$gdp_growth, k = 1))   # p > 0.05, No estacionaria
print(adf.test(Uruguay$unemployment, k = 1)) # p > 0.05, No estacionaria

# Creación de Lags y Diferencias
create_lags_diffs <- function(df) {
  df %>% 
    mutate(
      net_approval_lag = lag(net_approval),
      net_approval_diff = difference(net_approval),
      gdp_growth_lag = lag(gdp_growth),
      gdp_growth_diff = difference(gdp_growth),
      unemployment_lag = lag(unemployment),
      unemployment_diff = difference(unemployment)
    )
}

Uruguay <- create_lags_diffs(Uruguay)
Venezuela <- create_lags_diffs(Venezuela)

# Variables en primeras diferencias
# Resultado: La diferenciación corrige la raíz unitaria en la mayoría de los casos.
cat("--- Pruebas ADF (Diferencias) ---\n")
print(adf.test(na.omit(Uruguay$net_approval_diff), k = 1)) # Estacionaria
print(adf.test(na.omit(Venezuela$net_approval_diff), k = 1)) # Marginalmente estacionaria (p=0.06)

# 5. MODELADO DE SERIES DE TIEMPO ------------------------------------------
# Variable dependiente: net_approval
# Variables independientes: gdp_growth, unemployment

# --- 5.1 Modelo Estático ---
# Yt = B0 + B1Xt + et
mod_static_uy <- lm(net_approval ~ gdp_growth + unemployment, data = Uruguay)
mod_static_vz <- lm(net_approval ~ gdp_growth + unemployment, data = Venezuela)

# Diagnóstico Durbin-Watson (Autocorrelación)
# Resultado: Alta autocorrelación en ambos modelos.
durbinWatsonTest(mod_static_uy)
durbinWatsonTest(mod_static_vz)

# --- 5.2 Modelo de Ajuste Parcial ---
# Inclusión del rezago de la variable dependiente.
mod_partial_uy <- lm(net_approval ~ net_approval_lag + gdp_growth + unemployment, data = Uruguay)
mod_partial_vz <- lm(net_approval ~ net_approval_lag + gdp_growth + unemployment, data = Venezuela)

durbinWatsonTest(mod_partial_uy) # Uruguay: Corrige autocorrelación (p > 0.05)
durbinWatsonTest(mod_partial_vz) # Venezuela: Persiste autocorrelación

# --- 5.3 Modelo FDL (Finite Distributed Lag) ---
# Inclusión de rezagos de las independientes.
mod_fdl_uy <- lm(net_approval ~ gdp_growth + gdp_growth_lag + unemployment + unemployment_lag, data = Uruguay)
mod_fdl_vz <- lm(net_approval ~ gdp_growth + gdp_growth_lag + unemployment + unemployment_lag, data = Venezuela)

# --- 5.4 Modelo de Diferencias ---
# Trabajo con variables diferenciadas para asegurar estacionariedad.
mod_diff_uy <- lm(net_approval_diff ~ gdp_growth_diff + unemployment_diff, data = Uruguay)
mod_diff_vz <- lm(net_approval_diff ~ gdp_growth_diff + unemployment_diff, data = Venezuela)

# --- 5.5 Modelo "Dead Start" ---
# Dinámica solo en las independientes (De Boef & Keele).
mod_dead_uy <- lm(net_approval ~ gdp_growth + gdp_growth_lag + unemployment + unemployment_lag, data = Uruguay)
mod_dead_vz <- lm(net_approval ~ gdp_growth + gdp_growth_lag + unemployment + unemployment_lag, data = Venezuela)

# --- 5.6 Modelo de Corrección de Errores (ECM) ---
# Captura relaciones de corto y largo plazo.
mod_ecm_uy <- lm(net_approval_diff ~ net_approval_lag + gdp_growth_diff + gdp_growth_lag + unemployment_diff + unemployment_lag, data = Uruguay)
mod_ecm_vz <- lm(net_approval_diff ~ net_approval_lag + gdp_growth_diff + gdp_growth_lag + unemployment_diff + unemployment_lag, data = Venezuela)

durbinWatsonTest(mod_ecm_uy) # Buen ajuste en Uruguay
durbinWatsonTest(mod_ecm_vz) # Problemas persistentes en Venezuela

# --- 5.7 Modelo ARDL (Autoregressive Distributed Lag) ---
mod_ardl_uy <- lm(net_approval ~ net_approval_lag + gdp_growth + gdp_growth_lag + unemployment + unemployment_lag, data = Uruguay)
mod_ardl_vz <- lm(net_approval ~ net_approval_lag + gdp_growth + gdp_growth_lag + unemployment + unemployment_lag, data = Venezuela)

# --- 5.8 Modelo Curva de Phillips ---
# Relación entre inflación (aprobación diff) y desempleo/actividad económica.
mod_phillips_uy <- lm(net_approval_diff ~ gdp_growth + gdp_growth_lag + unemployment_diff, data = Uruguay)
mod_phillips_vz <- lm(net_approval_diff ~ gdp_growth + gdp_growth_lag + unemployment_diff, data = Venezuela)

durbinWatsonTest(mod_phillips_uy) # Sin autocorrelación (D-W ~ 1.74)
durbinWatsonTest(mod_phillips_vz) # Con autocorrelación

# 6. CONCLUSIONES ----------------------------------------------------------

# Con base en los diagnósticos Durbin-Watson y la teoría (De Boef & Keele):
# - Para Uruguay: El modelo ECM y ARDL muestran buen comportamiento, eliminando la autocorrelación.
# - Para Venezuela: La persistencia de autocorrelación sugiere que las dinámicas políticas 
#   responden a factores estructurales no capturados solo por variables económicas.

# Guardar resultados
# saveRDS(mod_ecm_uy, "modelos/ecm_uruguay.rds")
