import pandas as pd
import streamlit as st
import evapo as ev

import pvlib
from pvlib.modelchain import ModelChain
from pvlib.location import Location
from pvlib.pvsystem import PVSystem
from pvlib.temperature import TEMPERATURE_MODEL_PARAMETERS

tmy = pd.read_csv(r'FPV_results.csv', index_col=0)
tmy.index = pd.date_range(start='2022-01-01 00:00', end='2022-12-31 23:00', freq="H", tz="America/Santiago")
module = pvlib.pvsystem.retrieve_sam('SandiaMod')["Canadian_Solar_CS5P_220M___2009_"]
inverter = pvlib.pvsystem.retrieve_sam('CECInverter')["iPower__SHO_1_1__120V_"]
temperature_parameters = {'u_c': 29.0, 'u_v': 0} #TEMPERATURE_MODEL_PARAMETERS['pvsyst']['freestanding']

altitude = 405

#c = st.slider('FPV Capacity (kWp): ', 0,1000, 300, 50)
lat = st.number_input("latitude", min_value=-50.0, max_value=-20.0, value=-33.8991)    
long = st.number_input("longitude", min_value=-80.0, max_value=-68.0, value=-70.7320)    
#creating a sample data consisting different points 

df_pos = pd.DataFrame([[lat, long]],columns=['latitude', 'longitude'])

#plotting a map with the above defined points

st.map(df_pos)



if st.button("Download TMY Data"):
    tmy_pvg_r  = pvlib.iotools.get_pvgis_tmy(lat, long, outputformat='json', usehorizon=True, userhorizon=None, startyear=None, endyear=None, url='https://re.jrc.ec.europa.eu/api/v5_2/', map_variables=None, timeout=30)[0]
    # move 3 first rows to back to convert to Chilean time 
    tmy_pvg = tmy_pvg_r.iloc[3:].append(tmy_pvg_r.iloc[:3])
    tmy_pvg.index = pd.date_range(start = "2022-01-01 00:00", end="2022-12-31 23:00", freq="h", tz="America/Santiago")
    # Rename for pvlib
    cols_to_use = ["T2m", "G(h)", "Gb(n)", "Gd(h)", "IR(h)", "WS10m", "RH", "SP"]
    pvlib_column_names = ["temp_air", "ghi", "dni", "dhi", "lwr_u", "wind_speed", "rh", "sp" ]
    tmy_pvg = tmy_pvg[cols_to_use]
    tmy_pvg.columns = pvlib_column_names
    tmy = tmy_pvg
    st.write('successfully updated tmy data')
else:
    st.write('Press to update tmy data')

tilt = st.slider('PV Surface tilt (deg): ', 0, 90, 10, 5)   
azimuth = st.slider('Azimut (deg): ', 0, 360, 0, 5)

system = PVSystem(surface_tilt = tilt, surface_azimuth=azimuth, module_parameters=module, inverter_parameters=inverter, temperature_model_parameters=temperature_parameters, modules_per_string=4, strings_per_inverter=1)
location = Location(lat, long, tz = 'America/Santiago', altitude = altitude)

m = ModelChain(system,location)
m.run_model(tmy)
chart_data = m.results.ac
st.line_chart(chart_data)
st.write("The total FPV Generation is " + str(round(chart_data.sum()/1000,2))+ " kWh in one year")

alb = st.slider('albedo %: :', 0,100, 8, 1)

def calculate_pet(surface_pressure_KPa,     # surface pressure KPa
                  temperature2m_C,          # Daily mean temperature at 2 m
                  rel_hum_p,                # Relative humidity
                  windspeed10m_m_s,         # Windspeed at 10 m
                  sw_down_radiation_W_m2,   # Daily downward shortwave/solar radiation W/m2
                  lw_up_radiation_W_m2,     # Daily downward longwave radiation
                  albedo,                   # Albedo of water
                  soil_hf,                  # factor used to get the soil heat flux
                  pet_time):                # 'daily' or 'hourly'  ETo value
    """
    This is the function that calculate the PET based on the PM method.
    """
    # Constants.
    lmbda = 2.45  # Latent heat of vaporization [MJ kg -1] (simplification in the FAO PenMon (latent heat of about 20°C)
    cp = 1.013e-3 # Specific heat at constant pressure [MJ kg-1 °C-1]
    eps = 0.622   # Ratio molecular weight of water vapour/dry air
    b = 243.04
    a = 17.625
    a_1 = 0.34
    b_1 = -0.14
    a_c = 1.35
    b_c = -0.35
    sbc = 4.903*10**-9

    # Wind Speed
    windspeed2m_m_s = windspeed10m_m_s*(4.87/(np.log(67.8*10-5.42)))

    # Dewpoint Temperature based on relative humidity
    dewpoint2m_C = (b * (np.log(rel_hum_p/100) + a * temperature2m_C/(b+temperature2m_C))) / (a - (np.log(rel_hum_p/100) + a * temperature2m_C/(b+temperature2m_C)))

    #  Soil heat flux density [MJ m-2 day-1] - set to 0 following eq 42 in FAO
    G = soil_hf    
   
    # Atmospheric pressure [kPa] eq 7 in FAO.
    P_kPa = surface_pressure_KPa/1000 #101.3*((293.0-0.0065*height_m) / 293.0)**5.26

    # Psychrometric constant (gamma symbol in FAO) eq 8 in FAO.
    psychometric_kPa_c = cp*P_kPa / (eps*lmbda)

    # Saturation vapour pressure, eq 11 in FAO.
    svp_kPa = 0.6108*np.exp((17.27*temperature2m_C) / (temperature2m_C+237.3))

    # Delta (slope of saturation vapour pressure curve) eq 13 in FAO.
    delta_kPa_C = 4098.0*svp_kPa / (temperature2m_C+237.3)**2

    # Actual vapour pressure, eq 14 in FAO.
    avp_kPa = 0.6108*np.exp((17.27*dewpoint2m_C) / (dewpoint2m_C+237.3))

    # Saturation vapour pressure deficit.
    svpdeficit_kPa = svp_kPa - avp_kPa

    # Net radiation
    net_sw_radiation_MJ_m2 = (1 - albedo) * sw_down_radiation_W_m2 *0.0036
    #out_lw_radiation_MJ_m2 = lw_up_radiation_W_m2 *0.0036
    out_lw_radiation_MJ_m2 =   sbc * ( ((temperature2m_C+273.3)**4)) * (a_1 + b_1 * avp_kPa**0.5) * (a_c * 1 + b_c) / 24
    net_radiation_MJ_m2 = net_sw_radiation_MJ_m2 - out_lw_radiation_MJ_m2


    if pet_time == 'daily':
        # Calculate ET0, equation 6 in FAO
        numerator = 0.408*delta_kPa_C*(net_radiation_MJ_m2 - G) + \
            psychometric_kPa_c*(900/(temperature2m_C+273))*windspeed2m_m_s*svpdeficit_kPa
        denominator = delta_kPa_C + psychometric_kPa_c*(1 + 0.34*windspeed2m_m_s)
    
        ET0_mm_day = numerator / denominator
        return ET0_mm_day
    
    elif pet_time == 'hourly':
        # Calculate ET0, equation 53 in FAO
        numerator = 0.408*delta_kPa_C*(net_radiation_MJ_m2 - G) + \
            psychometric_kPa_c*(37/(temperature2m_C+273))*windspeed2m_m_s*svpdeficit_kPa
        denominator = delta_kPa_C + psychometric_kPa_c*(1 + 0.34*windspeed2m_m_s)
    
        ET0_mm_hr = numerator / denominator

        return ET0_mm_hr
    
    else:
        raise ValueError("time only takes 'daily' or 'hourly'")

chart_data2 = calculate_pet(tmy["sp"], tmy["temp_air"], tmy["rh"], tmy["wind_speed"], tmy["ghi"], tmy["lwr_u"], alb/100, 0, "hourly")
#chart_data2.index = pd.date_range(start='2022-01-01 00:00', end='2022-12-31 23:00', freq="H", tz="America/Santiago")

st.line_chart(chart_data2)
st.write("The total evaporation is " + str(round(chart_data2.sum(),2))+ " mm/m2 in one year")