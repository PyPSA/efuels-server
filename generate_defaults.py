import pandas as pd, urllib.request, os, yaml

df = pd.read_csv("defaults-initial.csv",
                 index_col=[0,1],
                 na_filter=False)

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)


years = config["tech_years"]


#get technology data

for year in years:
    fn = f"costs_{year}.csv"
    url = f"https://raw.githubusercontent.com/PyPSA/technology-data/{config['tech_data_commit']}/outputs/{fn}"

    if not os.path.isfile(fn):
        print("downloading",fn)
        urllib.request.urlretrieve(url,fn)

    td = {}
    for year in years:
        fn = f"costs_{year}.csv"
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
                               ("wind","onwind","Onshore wind"),
                               ("solar","solar-utility","Utility-scale solar PV"),
                               ("hydrogen_electrolyser","electrolysis","Hydrogen electrolyser"),
                               ("hydrogen_energy","hydrogen storage underground","Hydrogen underground storage"),
                               ("hydrogen_submarine_pipeline","H2 (g) submarine pipeline","Hydrogen submarine pipeline"),
                               ("battery_energy","battery storage","Utility-scale battery energy"),
                               ("battery_power","battery inverter","Utility-scale battery converter power")]:
    print(name,full_name)
    df.loc[(name + "_discount",""),:] = ["f",5,"percent",full_name + " discount rate",""]
    for year in years:
        value = td[year].loc[(td_name,"investment"),"value"]
        unit = td[year].loc[(td_name,"investment"),"unit"]
        if "EUR/MW" in unit and name != "hydrogen_submarine_pipeline":
            unit = unit.replace("EUR/MW","EUR/kW")
            value /= 1e3

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


for name,td_name,full_name in [("battery_power_efficiency_charging","battery inverter","Battery power charging efficiency"),
                               ("battery_power_efficiency_discharging","battery inverter","Battery power discharging efficiency"),
                               ("hydrogen_electrolyser_efficiency","electrolysis","Hydrogen electrolyser efficiency")]:


    for year in years:
        value = 100*td[year].loc[(td_name,"efficiency"),"value"]
        unit = "percent"
        if "battery" in name:
            value = 100*((value/100.)**0.5)
        elif "electrolyser" in name:
            unit ='"percent, LHV"'

        df.loc[(name,str(year)),:] = ["f",
                                      value,
                                      unit,
                                      full_name,
                                      td[year].loc[(td_name,"efficiency"),"source"]]


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


for name,td_name in [("methanol","MeOH")]:#,("lch4","CH4 (l)"),("ft","FT fuel"),("lh2","H2 (l)"),("lohc","LOHC"),("nh3","NH3 (l)")]


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
