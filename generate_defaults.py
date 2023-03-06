import pandas as pd, urllib.request, os, yaml

df = pd.read_csv("defaults-initial.csv",
                 index_col=[0,1],
                 na_filter=False)

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)


years = config["tech_years"]


#get technology data
td = {}
for year in years:

    fn = f"costs_{year}.csv"
    url = f"https://raw.githubusercontent.com/PyPSA/technology-data/{config['tech_data_commit']}/outputs/{fn}"

    if not os.path.isfile(fn):
        print("downloading",fn)
        urllib.request.urlretrieve(url,fn)

    td[year] = pd.read_csv(fn,
                           index_col=[0,1])


#get traces efficiencies

fn = "efficiencies.csv"
url = f"https://raw.githubusercontent.com/euronion/trace/{config['trace_commit']}/data/{fn}"
if not os.path.isfile(fn):
    print("downloading",fn)
    urllib.request.urlretrieve(url,fn)

eff = pd.read_csv(fn,
                  index_col=[0,1,2])


#get traces shipping data

fn = "shipping.csv"
url = f"https://raw.githubusercontent.com/euronion/trace/{config['trace_commit']}/data/{fn}"
if not os.path.isfile(fn):
    print("downloading",fn)
    urllib.request.urlretrieve(url,fn)

shipping = pd.read_csv(fn,
                       index_col=[0,1])



for name,td_name,full_name in [("methanolisation","methanolisation","Methanol synthesis"),
                               ("methanation","methanation","Methanation"),
                               ("methane_liquefaction","CH4 liquefaction","Methane liquefaction"),
                               ("ft","Fischer-Tropsch","Fischer-Tropsch synthesis"),
                               ("hydrogen_liquefaction","H2 liquefaction","Hydrogen liquefaction"),
                               ("wind","onwind","Onshore wind"),
                               ("solar","solar-utility","Utility-scale solar PV"),
                               ("hydrogen_electrolyser","electrolysis","Hydrogen electrolyser"),
                               ("hydrogen_energy","hydrogen storage tank type 1","Hydrogen overground tank storage"),
                               ("methane_storage","methane storage tank incl. compressor","Methane overground tank storage"),
                               ("hydrogen_turbine","CCGT","Hydrogen combined cycle turbine"),
                               ("hydrogen_submarine_pipeline","H2 (g) submarine pipeline","Hydrogen submarine pipeline"),
                               ("hvdc_submarine_cable","HVDC submarine","HVDC submarine cable"),
                               ("hvdc_inverter_pair","HVDC inverter pair","HVDC AC-DC inverter pair"),
                               ("battery_energy","battery storage","Utility-scale battery energy"),
                               ("battery_power","battery inverter","Utility-scale battery converter power"),
                               ("dac","direct air capture","Direct air capture"),
                               ("heat_pump","industrial heat pump medium temperature","Industrial heat pump up to 125 C"),
                               ("liquid_carbonaceous_storage","General liquid hydrocarbon storage (product)","Liquid carbonaceous fuel storage tank"),
                               ("air_separation_unit","air separation unit","Air separation unit for nitrogen"),
                               ("haber_bosch","Haber-Bosch","Haber-Bosch ammonia synthesis"),
                               ("ammonia_storage","NH3 (l) storage tank incl. liquefaction","Liquid ammonia storage"),
                               ("co2_storage","CO2 storage tank","CO2 storage tank"),
                               ("hydrogen_compressor","hydrogen storage compressor","Hydrogen storage compressor")]:
    print(name,full_name)
    df.loc[(name + "_discount",""),:] = ["f",5,"percent",full_name + " discount rate",""]
    for year in years:
        value = td[year].loc[(td_name,"investment"),"value"]
        unit = td[year].loc[(td_name,"investment"),"unit"]

        df.loc[(name + "_cost",str(year)),:] = ["f",
                                                value,
                                                unit,
                                                full_name + " capital cost (overnight)",
                                                td[year].loc[(td_name,"investment"),"source"]]
        df.loc[(name + "_fom",str(year)),:] = ["f",
                                               td[year].loc[(td_name,"FOM"),"value"] if (td_name,"FOM") in td[year].index else 0,
                                               "percent of overnight cost per year",
                                               full_name + " fixed operation and maintenance costs",
                                               td[year].loc[(td_name,"FOM"),"source"] if (td_name,"FOM") in td[year].index else "default"]
        df.loc[(name + "_lifetime",str(year)),:] = ["f",
                                                    td[year].loc[(td_name,"lifetime"),"value"],
                                                    td[year].loc[(td_name,"lifetime"),"unit"],
                                                    full_name + " lifetime",
                                                    td[year].loc[(td_name,"lifetime"),"source"]]

#fix mistake in TD
for year in years:
    df.loc[("methanation_cost",str(year)),"unit"] = "EUR/kW_CH4"

for name,td_name,full_name in [("battery_power_efficiency_charging","battery inverter","Battery power charging efficiency"),
                               ("battery_power_efficiency_discharging","battery inverter","Battery power discharging efficiency"),
                               ("heat_pump_efficiency","industrial heat pump medium temperature","Industrial heat pump COP"),
                               ("hydrogen_electrolyser_efficiency","electrolysis","Hydrogen electrolyser efficiency"),
                               ("hydrogen_turbine_efficiency","CCGT","Hydrogen combined cycle turbine efficiency")]:

    for year in years:
        value = 100*td[year].loc[(td_name,"efficiency"),"value"]
        unit = "percent"
        if "battery" in name:
            value = 100*((value/100.)**0.5)
        elif "hydrogen" in name:
            unit ='"percent, LHV"'

        df.loc[(name,str(year)),:] = ["f",
                                      value,
                                      unit,
                                      full_name,
                                      td[year].loc[(td_name,"efficiency"),"source"]]


#hydrogen_submarine_pipeline_losses,,f,2.1,percent/1000km,Hydrogen submarine pipeline losses,

df.loc[("hydrogen_submarine_pipeline_losses",""),:] = ["f",
                                                       100*(1-eff.loc[("H2 (g) submarine pipeline","all","hydrogen (g) compressed submarine"),"efficiency"][0]),
                                                       "percent/1000km",
                                                       "Hydrogen submarine pipeline losses",
                                                       eff.loc[("H2 (g) submarine pipeline","all","hydrogen (g) compressed submarine"),"source"][0]]

df.loc[("hvdc_submarine_cable_losses",""),:] = ["f",
                                                100*(1-eff.loc[("HVDC submarine","all","hvdc submarine"),"efficiency"][0]),
                                                "percent/1000km",
                                                "HVDC submarine cable losses",
                                                eff.loc[("HVDC submarine","all","hvdc submarine"),"source"][0]]



df.loc[("methanolisation_efficiency",""),:] = ["f",
                                               eff.loc[("methanolisation","all","hydrogen (g)"),"to_amount"][0]/eff.loc[("methanolisation","all","hydrogen (g)"),"from_amount"][0],
                                               "MWh-MeOH-LHV/MWh-H2-LHV",
                                               "Methanol synthesis efficiency wrt hydrogen",
                                               eff.loc[("methanolisation","all","hydrogen (g)"),"source"][0]]

df.loc[("methanolisation_co2",""),:] = ["f",
                                        eff.loc[("methanolisation","all","CO2 (g)"),"from_amount"][0]/eff.loc[("methanolisation","all","CO2 (g)"),"to_amount"][0],
                                        "tCO2/MWh-MeOH-LHV",
                                        "Methanol synthesis carbon dioxide input",
                                        eff.loc[("methanolisation","all","CO2 (g)"),"source"][0]]

df.loc[("methanolisation_electricity",""),:] = ["f",
                                                eff.loc[("methanolisation","all","electricity"),"from_amount"][0]/eff.loc[("methanolisation","all","electricity"),"to_amount"][0],
                                                "MWhel/MWh-MeOH-LHV",
                                                "Methanol synthesis electricity input",
                                                eff.loc[("methanolisation","all","electricity"),"source"][0]]


df.loc[("ft_efficiency",""),:] = ["f",
                                  eff.loc[("Fischer-Tropsch","all","hydrogen (g)"),"to_amount"][0]/eff.loc[("Fischer-Tropsch","all","hydrogen (g)"),"from_amount"][0],
                                  "MWh-FT-LHV/MWh-H2-LHV",
                                  "Fischer-Tropsch synthesis efficiency wrt hydrogen",
                                  eff.loc[("Fischer-Tropsch","all","hydrogen (g)"),"source"][0]]

df.loc[("ft_co2",""),:] = ["f",
                           eff.loc[("Fischer-Tropsch","all","CO2 (g)"),"from_amount"][0]/eff.loc[("Fischer-Tropsch","all","CO2 (g)"),"to_amount"][0],
                           "tCO2/MWh-FT-LHV",
                           "Fischer-Tropsch synthesis carbon dioxide input",
                           eff.loc[("Fischer-Tropsch","all","CO2 (g)"),"source"][0]]

df.loc[("ft_electricity",""),:] = ["f",
                                   eff.loc[("Fischer-Tropsch","all","electricity"),"from_amount"][0]/eff.loc[("Fischer-Tropsch","all","electricity"),"to_amount"][0],
                                   "MWhel/MWh-FT-LHV",
                                   "Fischer-Tropsch synthesis electricity input",
                                   eff.loc[("Fischer-Tropsch","all","electricity"),"source"][0]]


df.loc[("methanation_efficiency",""),:] = ["f",
                                               eff.loc[("methanation","all","hydrogen (g)"),"to_amount"][0]/eff.loc[("methanation","all","hydrogen (g)"),"from_amount"][0],
                                               "MWh-CH4-LHV/MWh-H2-LHV",
                                               "Methanation efficiency wrt hydrogen",
                                               eff.loc[("methanation","all","hydrogen (g)"),"source"][0]]

df.loc[("methanation_co2",""),:] = ["f",
                                        eff.loc[("methanation","all","CO2 (g)"),"from_amount"][0]/eff.loc[("methanation","all","CO2 (g)"),"to_amount"][0],
                                        "tCO2/MWh-CH4-LHV",
                                        "Methanation carbon dioxide input",
                                        eff.loc[("methanation","all","CO2 (g)"),"source"][0]]

df.loc[("methane_liquefaction_electricity",""),:] = ["f",
                                                     eff.loc[("CH4 liquefaction","all","electricity"),"from_amount"][0]/eff.loc[("CH4 liquefaction","all","electricity"),"to_amount"][0]/config['mwh_per_t']['methane'],
                                                     "MWh-el/MWh-CH4-LHV",
                                                     "Methane liquefaction electricity input",
                                                     eff.loc[("CH4 liquefaction","all","electricity"),"source"][0]]


df.loc[("hydrogen_liquefaction_efficiency",""),:] = ["f",
                                                     eff.loc[("H2 liquefaction","all","hydrogen (g)"),"to_amount"][0]/eff.loc[("H2 liquefaction","all","hydrogen (g)"),"from_amount"][0],
                                                     "MWh-H2-liquid/MWh-H2-gas",
                                                     "Hydrogen liquefaction efficiency",
                                                     eff.loc[("H2 liquefaction","all","hydrogen (g)"),"source"][0]]


df.loc[("hydrogen_liquefaction_electricity",""),:] = ["f",
                                                     eff.loc[("H2 liquefaction","all","electricity"),"from_amount"][0]/eff.loc[("H2 liquefaction","all","electricity"),"to_amount"][0],
                                                     "MWh-el/MWh-H2-LHV",
                                                     "Hydrogen liquefaction electricity input",
                                                     eff.loc[("H2 liquefaction","all","electricity"),"source"][0]]


df.loc[("dac_electricity",""),:] = ["f",
                                    td[year].loc[("direct air capture","electricity-input"),"value"],
                                    td[year].loc[("direct air capture","electricity-input"),"unit"],
                                    "Direct air capture electricity consumption",
                                    td[year].loc[("direct air capture","electricity-input"),"source"]]

df.loc[("dac_heat",""),:] = ["f",
                             td[year].loc[("direct air capture","heat-input"),"value"],
                             td[year].loc[("direct air capture","heat-input"),"unit"],
                             "Direct air capture heat consumption",
                             td[year].loc[("direct air capture","heat-input"),"source"]]

df.loc[("air_separation_unit_efficiency",""),:] = ["f",
                                                   eff.loc[("air separation unit","all","electricity"),"to_amount"][0]/eff.loc[("air separation unit","all","electricity"),"from_amount"][0],
                                                   "tN2/MWhel",
                                                   "Air separation unit output",
                                                   eff.loc[("air separation unit","all","electricity"),"source"][0]]
df.loc[("haber_bosch_efficiency",""),:] = ["f",
                                           eff.loc[("Haber-Bosch","all","hydrogen (g)"),"to_amount"][0]/eff.loc[("Haber-Bosch","all","hydrogen (g)"),"from_amount"][0],
                                           "MWh-NH3-LHV/MWh-H2-LHV",
                                           "Haber-Bosch ammonia synthesis hydrogen input",
                                           eff.loc[("Haber-Bosch","all","hydrogen (g)"),"source"][0]]

df.loc[("haber_bosch_electricity",""),:] = ["f",
                                            eff.loc[("Haber-Bosch","all","electricity"),"from_amount"][0]/eff.loc[("Haber-Bosch","all","electricity"),"to_amount"][0],
                                            "MWhel/MWh-NH3-LHV",
                                            "Haber-Bosch ammonia synthesis electricity input",
                                            eff.loc[("Haber-Bosch","all","electricity"),"source"][0]]


df.loc[("haber_bosch_nitrogen",""),:] = ["f",
                                            eff.loc[("Haber-Bosch","all","nitrogen (g)"),"from_amount"][0]/eff.loc[("Haber-Bosch","all","nitrogen (g)"),"to_amount"][0],
                                            "tN2/MWh-NH3-LHV",
                                            "Haber-Bosch ammonia synthesis nitrogen input",
                                            eff.loc[("Haber-Bosch","all","nitrogen (g)"),"source"][0]]


df.loc[("hydrogen_compressor_electricity",""),:] = ["f",
                                                    eff.loc[("H2 storage compressor","all","electricity"),"from_amount"][0]/eff.loc[("H2 storage compressor","all","electricity"),"to_amount"][0],
                                                    "MWhel/MWh-H2-LHV",
                                                    "Hydrogen storage compressor electricity input",
                                                    eff.loc[("H2 storage compressor","all","electricity"),"source"][0]]



for name,td_name in [("methanol","MeOH"),("ammonia","NH3 (l)"),("methane","CH4 (l)"),("lh2","H2 (l)"),("ft","FT fuel")]:#,("lohc","LOHC")]


    df.loc[(name + "_ship_discount",""),:] = ["f",5,"percent",name + " shipping discount rate",""]


    fuel_df = shipping.loc[td_name + " transport ship"]


    for attr,td_attr,full in [("loading_loss","(un-) loading losses","unloading/loading loss"),
                              ("loading_time","(un-) loading time","unloading/loading time"),
                              ("average_speed","average speed","average speed"),
                              ("capacity_mwh","capacity","energy capacity"),
                              ("energy_demand","energy demand","energy demand")]:

        df.loc[(name + "_ship_" + attr,""),:] = ["f",
                                                 fuel_df.at[td_attr,"value"],
                                                 fuel_df.at[td_attr,"unit"],
                                                 name + " shipping " + full,
                                                 fuel_df.at[td_attr,"comment"]]


    for attr,td_attr,full in [("fom","FOM","fixed operation and maintenance costs"),
                              ("cost","investment","capital cost (overnight)"),
                              ("lifetime","lifetime","lifetime"),
                              ("capacity_t","capacity","mass capacity")]:

        for year in years:
            df.loc[(name + "_ship_" + attr,str(year)),:] = ["f",
                                                            td[year].loc[(td_name + " transport ship",td_attr),"value"],
                                                            td[year].loc[(td_name + " transport ship",td_attr),"unit"],
                                                            name + " shipping " + full,
                                                            td[year].loc[(td_name + " transport ship",td_attr),"source"]]

print(df)

df.to_csv("defaults.csv")
