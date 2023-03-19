## Copyright 2022-3 Tom Brown

## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU Affero General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero General Public License for more details.

## License and more information at:
## https://github.com/PyPSA/efuels-server



import pypsa

from pypsa.linopt import get_var, linexpr, define_constraints

from pypsa.geo import haversine

import pandas as pd
from pyomo.environ import Constraint
from rq import get_current_job

import json, os, hashlib, yaml

#required to stop wierd failures
import netCDF4

from atlite.gis import spdiag, compute_indicatormatrix

import xarray as xr

import scipy as sp

import numpy as np

from shapely.geometry import box, Point, Polygon, MultiPolygon

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

defaults = pd.read_csv("defaults.csv",index_col=[0,1],na_filter=False)


threshold = config["results_threshold"]


#based on mean deviation against renewables.ninja capacity factors for European countries for 2011-2013
solar_correction_factor = 0.926328


override_component_attrs = pypsa.descriptors.Dict(
        {k: v.copy() for k, v in pypsa.components.component_attrs.items()}
    )
override_component_attrs["Link"].loc["bus2"] = [
        "string",
        np.nan,
        np.nan,
        "2nd bus",
        "Input (optional)",
    ]
override_component_attrs["Link"].loc["bus3"] = [
        "string",
        np.nan,
        np.nan,
        "3rd bus",
        "Input (optional)",
    ]
override_component_attrs["Link"].loc["efficiency2"] = [
        "static or series",
        "per unit",
        1.0,
        "2nd bus efficiency",
        "Input (optional)",
    ]
override_component_attrs["Link"].loc["efficiency3"] = [
        "static or series",
        "per unit",
        1.0,
        "3rd bus efficiency",
        "Input (optional)",
    ]
override_component_attrs["Link"].loc["p2"] = [
        "series",
        "MW",
        0.0,
        "2nd bus output",
        "Output",
    ]
override_component_attrs["Link"].loc["p3"] = [
        "series",
        "MW",
        0.0,
        "3rd bus output",
        "Output",
    ]


def annuity(lifetime,rate):
    if rate == 0.:
        return 1/lifetime
    else:
        return rate/(1. - 1. / (1. + rate)**lifetime)



assumptions_df = pd.DataFrame(columns=["FOM","fixed","discount rate","lifetime","investment"],
                              dtype=float)



def error(message, jobid):
    with open('results-solve/results-{}.json'.format(jobid), 'w') as fp:
        json.dump({"jobid" : jobid,
                   "status" : "Error",
                   "error" : message
                   },fp)
    print("Error: {}".format(message))
    return {"error" : message}


def find_interval(interval_start,interval_length,value):
    return int((value-interval_start)//interval_length)


def get_octant_bounds(quadrant, hemisphere):

    x0 = -180 + quadrant*90.
    x1 = x0 + 90.

    y0 = -90. + hemisphere*90.
    y1 = y0 + 90.

    return x0,x1,y0,y1

def generate_octant_grid_cells(quadrant, hemisphere, mesh=0.5):

    x0,x1,y0,y1 = get_octant_bounds(quadrant, hemisphere)

    x = np.arange(x0,
                  x1 + mesh,
                  mesh)

    y = np.arange(y0,
                  y1 + mesh,
                  mesh)

    #grid_coordinates and grid_cells copied from atlite/cutout.py
    xs, ys = np.meshgrid(x,y)
    grid_coordinates = np.asarray((np.ravel(xs), np.ravel(ys))).T

    span = mesh / 2
    return [box(*c) for c in np.hstack((grid_coordinates - span, grid_coordinates + span))]


def get_octant(lon,lat):

    # 0 for lon -180--90, 1 for lon -90-0, etc.
    quadrant = find_interval(-180.,90,lon)

    #0 for lat -90 - 0, 1 for lat 0 - 90
    hemisphere = find_interval(-90,90,lat)

    print(f"octant is in quadrant {quadrant} and hemisphere {hemisphere}")

    rel_x = lon - quadrant*90 + 180.

    rel_y = lat - hemisphere*90 + 90.

    span = 0.5

    n_per_octant = int(90/span +1)

    i = find_interval(0-span/2,span,rel_x)
    j = find_interval(0-span/2,span,rel_y)

    position = j*n_per_octant+i

    print("position",position)

    #paranoid check
    if True:
        grid_cells = generate_octant_grid_cells(quadrant, hemisphere, mesh=span)
        assert grid_cells[position].contains(Point(lon,lat))

    return quadrant, hemisphere, position


def get_weather(lat, lon, year, cf_exponent):

    if lon < -180 or lon > 180 or lat > 90 or lat < -90:
        return "Point's coordinates not within lon*lat range of (-180,180)*(-90,90)", None

    quadrant, hemisphere, position = get_octant(lon,lat)

    pu = {}

    for tech in ["solar", "onwind"]:
        filename = os.path.join(config['octant_folder'],
                                f"octant-{year}-{quadrant}-{hemisphere}-{tech}.nc")
        o = xr.open_dataarray(filename)
        pu[tech] = o.loc[{"dim_0":position}].to_pandas()

    pu = pd.DataFrame(pu)

    if pu is not None:
        pu["solar"] = solar_correction_factor*pu["solar"]

    return pu, None


def export_time_series(n):

    bus_carriers = n.buses.carrier.unique()

    all_carrier_dict = {}

    for i in bus_carriers:
        bus_map = (n.buses.carrier == i)
        bus_map.at[""] = False

        carrier_df = pd.DataFrame(index=n.snapshots,
                                  dtype=float)

        for c in n.iterate_components(n.one_port_components):

            items = c.df.index[c.df.bus.map(bus_map).fillna(False)]

            if len(items) == 0:
                continue

            s = c.pnl.p[items].multiply(c.df.loc[items,'sign'],axis=1).groupby(c.df.loc[items,'carrier'],axis=1).sum()
            carrier_df = pd.concat([carrier_df,s],axis=1)

        for c in n.iterate_components(n.branch_components):

            for end in [col[3:] for col in c.df.columns if col[:3] == "bus"]:

                items = c.df.index[c.df["bus" + str(end)].map(bus_map,na_action=False)]

                if len(items) == 0:
                    continue

                s = (-1)*c.pnl["p"+end][items].groupby(c.df.loc[items,'carrier'],axis=1).sum()
                carrier_df = pd.concat([carrier_df,s],axis=1)

        all_carrier_dict[i] = carrier_df

    all_carrier_df = pd.concat(all_carrier_dict, axis=1)

    return all_carrier_df



def run_optimisation(assumptions, pu):
    """Needs cleaned-up assumptions and pu.
    return results_overview, results_series, error_msg"""


    Nyears = 1

    techs = [tech[:-5] for tech in assumptions if tech[-5:] == "_cost" and tech[-14:] != "_marginal_cost" and tech != "co2_cost"]

    print("calculating costs for",techs)

    for item in techs:
        assumptions_df.at[item,"discount rate"] = assumptions[item + "_discount"]/100.
        assumptions_df.at[item,"investment"] = assumptions[item + "_cost"]*1e3 if "EUR/kW" in defaults.loc[item + "_cost"]["unit"][0] else assumptions[item + "_cost"]
        assumptions_df.at[item,"FOM"] = assumptions[item + "_fom"]
        assumptions_df.at[item,"lifetime"] = assumptions[item + "_lifetime"]

    assumptions_df["fixed"] = [(annuity(v["lifetime"],v["discount rate"])+v["FOM"]/100.)*v["investment"]*Nyears for i,v in assumptions_df.iterrows()]


    distance = haversine([assumptions["destination_lng"],assumptions["destination_lat"]],
                         [assumptions["source_lng"],assumptions["source_lat"]])[0][0]

    print(f'distance between points is {distance} km')


    if assumptions["efuel"] in ["electricity_hvdc", "hydrogen_submarine_pipeline"]:
        distance_transported = assumptions["pipeline_distance_factor"]*distance
    elif assumptions["efuel"] in ["methanol","ammonia","methane","lh2","ft"]:
        distance_transported = assumptions["shipping_distance_factor"]*distance
    else:
        distance_transported = 0.

    print(f'distance transported between points is {distance_transported} km')

    print('Starting task with assumptions {}'.format(assumptions_df))

    network = pypsa.Network(override_component_attrs=override_component_attrs)

    snapshots = pd.date_range("{}-01-01".format(assumptions["year"]),"{}-12-31 23:00".format(assumptions["year"]),
                              freq=str(assumptions["frequency"])+"H")

    network.set_snapshots(snapshots)

    network.snapshot_weightings = pd.Series(float(assumptions["frequency"]),index=network.snapshots)

    network.add("Bus","electricity",
                carrier="electricity")

    if assumptions["solar"]:
        network.add("Generator","solar",
                    bus="electricity",
                    carrier="solar",
                    p_max_pu = pu["solar"],
                    p_nom_extendable = True,
                    p_nom_min = assumptions["solar_min"],
                    p_nom_max = assumptions["solar_max"],
                    marginal_cost = 0.1, #Small cost to prefer curtailment to destroying energy in storage, solar curtails before wind
                    capital_cost = assumptions_df.at['solar','fixed'])

    if assumptions["wind"]:
        network.add("Generator","wind",
                    bus="electricity",
                    carrier="wind",
                    p_max_pu = pu["onwind"],
                    p_nom_extendable = True,
                    p_nom_min = assumptions["wind_min"],
                    p_nom_max = assumptions["wind_max"],
                    marginal_cost = 0.2, #Small cost to prefer curtailment to destroying energy in storage, solar curtails before wind
                    capital_cost = assumptions_df.at['wind','fixed'])

    for i in range(1,3):
        name = "dispatchable" + str(i)
        if assumptions[name]:
            network.add("Carrier",name,
                        co2_emissions=assumptions[name+"_emissions"])
            network.add("Generator",name,
                        bus="electricity",
                        carrier=name,
                        p_nom_extendable=True,
                        marginal_cost=assumptions[name+"_marginal_cost"],
                        capital_cost=assumptions_df.at[name,'fixed'])

    if assumptions["battery"]:

        network.add("Bus","battery")

        network.add("Store","battery_energy",
                    bus = "battery",
                    carrier="battery storage",
                    e_nom_extendable = True,
                    e_cyclic=True,
                    capital_cost=assumptions_df.at['battery_energy','fixed'])

        network.add("Link","battery_power",
                    bus0 = "electricity",
                    bus1 = "battery",
                    carrier="battery inverter",
                    efficiency = assumptions["battery_power_efficiency_charging"]/100.,
                    p_nom_extendable = True,
                    capital_cost=assumptions_df.at['battery_power','fixed'])

        network.add("Link","battery_discharge",
                    bus0 = "battery",
                    bus1 = "electricity",
                    carrier="battery discharger",
                    p_nom_extendable = True,
                    efficiency = assumptions["battery_power_efficiency_discharging"]/100.)

        def extra_functionality(network,snapshots):

            link_p_nom = get_var(network, "Link", "p_nom")

            lhs = linexpr((1,link_p_nom["battery_power"]),
                          (-network.links.loc["battery_discharge", "efficiency"],
                           link_p_nom["battery_discharge"]))
            define_constraints(network, lhs, "=", 0, 'Link', 'charger_ratio')
    else:
        def extra_functionality(network,snapshots):
            pass

    network.add("Bus",
                "hydrogen",
                carrier="hydrogen")

    network.add("Bus",
                "water",
                carrier="water")

    network.add("Link",
                "hydrogen_electrolyser",
                carrier="hydrogen electrolyser",
                bus0="electricity",
                bus1="hydrogen",
                bus2="water",
                p_nom_extendable=True,
                efficiency=assumptions["hydrogen_electrolyser_efficiency"]/100.,
                efficiency2=-assumptions["hydrogen_electrolyser_water"]*assumptions["hydrogen_electrolyser_efficiency"]/100.,
                capital_cost=assumptions_df.at["hydrogen_electrolyser","fixed"])

    network.add("Link",
                "desalination",
                carrier="seawater desalination",
                bus0="electricity",
                bus1="water",
                p_nom_extendable=True,
                efficiency=1/assumptions["desalination_electricity"],
                capital_cost=assumptions_df.at["desalination","fixed"]/assumptions["desalination_electricity"])

    if assumptions["hydrogen"]:

        network.add("Bus",
                    "compressed hydrogen",
                    carrier="compressed hydrogen")

        network.add("Link",
                    "hydrogen_compressor",
                    carrier="hydrogen storing compressor",
                    bus0="hydrogen",
                    bus1="compressed hydrogen",
                    bus2="electricity",
                    p_nom_extendable=True,
                    efficiency=1,
                    efficiency2=-assumptions["hydrogen_compressor_electricity"],
                    capital_cost=assumptions_df.at["hydrogen_compressor","fixed"])

        network.add("Link",
                    "hydrogen_decompressor",
                    carrier="hydrogen storing decompressor",
                    bus0="compressed hydrogen",
                    bus1="hydrogen",
                    p_nom_extendable=True)

        network.add("Store",
                    "hydrogen_energy",
                    bus="compressed hydrogen",
                    carrier="hydrogen storage",
                    e_nom_extendable=True,
                    e_cyclic=True,
                    capital_cost=assumptions_df.at["hydrogen_energy","fixed"])

    #add CO2 source
    if assumptions["efuel"] in ["methanol","methane","ft"]:

        network.add("Bus",
                    "co2",
                    carrier="co2")

        if assumptions["dac"]:

            network.add("Bus",
                        "heat",
                        carrier="heat")

            network.add("Link",
                        "heat pump",
                        bus0="electricity",
                        bus1="heat",
                        carrier="heat pump",
                        p_nom_extendable=True,
                        capital_cost=assumptions_df.at["heat_pump","fixed"]*assumptions["heat_pump_efficiency"]/100.,
                        efficiency=assumptions["heat_pump_efficiency"]/100.)

            network.add("Link",
                        "dac",
                        bus0="electricity",
                        bus1="co2",
                        bus2="heat",
                        carrier="dac",
                        p_nom_extendable=True,
                        capital_cost=assumptions_df.at["dac","fixed"]/assumptions["dac_electricity"],
                        efficiency=1/assumptions["dac_electricity"],
                        efficiency2=-assumptions["dac_heat"]/assumptions["dac_electricity"])

            network.add("Store",
                        "co2",
                        bus="co2",
                        carrier="co2 storage",
                        e_nom_extendable=True,
                        e_cyclic=True,
                        capital_cost=assumptions_df.at["co2_storage","fixed"])
        else:
            network.add("Generator",
                        "co2",
                        bus="co2",
                        carrier="co2",
                        marginal_cost=assumptions["co2_cost"],
                        p_nom=8760*assumptions["efuels_load"])


    if assumptions["efuel"] == "hydrogen_submarine_pipeline":

        network.add("Bus",
                    "destination",
                    carrier="efuel")


        network.add("Load","hydrogen_load",
                    bus="destination",
                    carrier="efuel",
                    p_set=assumptions["efuels_load"])

        network.add("Link",
                    "hydrogen_submarine_pipeline",
                    bus0="hydrogen",
                    bus1="destination",
                    carrier="hydrogen submarine pipeline",
                    p_nom_extendable=True,
                    efficiency=(1-assumptions["hydrogen_submarine_pipeline_losses"]/100.)**(distance_transported/1000.),
                    capital_cost=assumptions_df.at["hydrogen_submarine_pipeline","fixed"]*distance_transported)

        print("pipeline efficiency is", (1-assumptions["hydrogen_submarine_pipeline_losses"]/100.)**(distance_transported/1000.))

        print("pieline cost per MW is",assumptions_df.at["hydrogen_submarine_pipeline","fixed"]*distance_transported)


    elif assumptions["efuel"] == "lh2":

        network.add("Bus",
                    "destination",
                    carrier="lh2")

        network.add("Bus",
                    "lh2",
                    carrier="lh2")

        network.add("Load","lh2_load",
                    bus="destination",
                    carrier="efuel",
                    p_set=assumptions["efuels_load"])

        network.add("Link",
                    "hydrogen liquefaction",
                    bus0="hydrogen",
                    bus1="lh2",
                    bus2="electricity",
                    carrier="hydrogen liquefaction",
                    p_nom_extendable=True,
                    efficiency=assumptions["hydrogen_liquefaction_efficiency"],
                    efficiency2=-assumptions["hydrogen_liquefaction_electricity"]*assumptions["hydrogen_liquefaction_efficiency"],
                    capital_cost=assumptions_df.at["hydrogen_liquefaction","fixed"])

    elif assumptions["efuel"] == "methanol":

        network.add("Bus",
                    "destination",
                    carrier="methanol")

        network.add("Bus",
                    "methanol",
                    carrier="methanol")

        network.add("Store",
                    "methanol",
                    bus="methanol",
                    carrier="methanol storage",
                    e_nom_extendable=True,
                    e_cyclic=True,
                    capital_cost=assumptions_df.at["liquid_carbonaceous_storage","fixed"]/config["mwh_per_m3"]["methanol"])

        network.add("Load","methanol_load",
                    bus="destination",
                    carrier="efuel",
                    p_set=assumptions["efuels_load"])

        network.add("Link",
                    "methanol synthesis",
                    bus0="hydrogen",
                    bus1="methanol",
                    bus2="electricity",
                    bus3="co2",
                    carrier="methanol synthesis",
                    p_nom_extendable=True,
                    p_min_pu=assumptions["methanolisation_min_part_load"]/100,
                    efficiency=assumptions["methanolisation_efficiency"],
                    efficiency2=-assumptions["methanolisation_electricity"]*assumptions["methanolisation_efficiency"],
                    efficiency3=-assumptions["methanolisation_co2"]*assumptions["methanolisation_efficiency"],
                    capital_cost=assumptions_df.at["methanolisation","fixed"]*assumptions["methanolisation_efficiency"]) #NB: cost is EUR/kW_MeOH

    elif assumptions["efuel"] == "ft":

        network.add("Bus",
                    "destination",
                    carrier="ft")

        network.add("Bus",
                    "ft",
                    carrier="ft")

        network.add("Store",
                    "ft",
                    bus="ft",
                    carrier="ft storage",
                    e_nom_extendable=True,
                    e_cyclic=True,
                    capital_cost=assumptions_df.at["liquid_carbonaceous_storage","fixed"]/config["mwh_per_m3"]["ft"])

        network.add("Load","ft_load",
                    bus="destination",
                    carrier="efuel",
                    p_set=assumptions["efuels_load"])

        network.add("Link",
                    "ft synthesis",
                    bus0="hydrogen",
                    bus1="ft",
                    bus2="electricity",
                    bus3="co2",
                    carrier="ft synthesis",
                    p_nom_extendable=True,
                    p_min_pu=assumptions["ft_min_part_load"]/100,
                    efficiency=assumptions["ft_efficiency"],
                    efficiency2=-assumptions["ft_electricity"]*assumptions["ft_efficiency"],
                    efficiency3=-assumptions["ft_co2"]*assumptions["ft_efficiency"],
                    capital_cost=assumptions_df.at["ft","fixed"]*assumptions["ft_efficiency"]) #NB: cost is EUR/kW_MeOH

    elif assumptions["efuel"] == "methane":

        network.add("Bus",
                    "destination",
                    carrier="methane")

        #this is liquefied methane
        network.add("Bus",
                    "methane",
                    carrier="methane")

        network.add("Bus",
                    "gaseous methane",
                    carrier="gaseous methane")

        network.add("Store",
                    "gaseous methane",
                    bus="gaseous methane",
                    carrier="methane storage",
                    e_nom_extendable=True,
                    e_cyclic=True,
                    capital_cost=assumptions_df.at["methane_storage","fixed"]/config["mwh_per_m3"]["methane"]/config["methane storage pressure"])

        network.add("Load","methane_load",
                    bus="destination",
                    carrier="efuel",
                    p_set=assumptions["efuels_load"])

        network.add("Link",
                    "methanation",
                    bus0="hydrogen",
                    bus1="gaseous methane",
                    bus2="co2",
                    carrier="methanation",
                    p_nom_extendable=True,
                    p_min_pu=assumptions["methanation_min_part_load"]/100,
                    efficiency=assumptions["methanation_efficiency"],
                    efficiency2=-assumptions["methanation_co2"]*assumptions["methanation_efficiency"],
                    capital_cost=assumptions_df.at["methanation","fixed"]*assumptions["methanation_efficiency"]) #NB: cost is EUR/MW_CH4

        network.add("Link",
                    "methane liquefaction",
                    bus0="gaseous methane",
                    bus1="methane",
                    bus2="electricity",
                    carrier="methane liquefaction",
                    p_nom_extendable=True,
                    efficiency=1,
                    efficiency2=-assumptions["methane_liquefaction_electricity"],
                    capital_cost=assumptions_df.at["methane_liquefaction","fixed"])

    elif assumptions["efuel"] == "ammonia":

        network.add("Bus",
                    "destination",
                    carrier="ammonia")

        network.add("Bus",
                    "ammonia",
                    carrier="ammonia")

        network.add("Store",
                    "ammonia",
                    bus="ammonia",
                    carrier="ammonia storage",
                    e_nom_extendable=True,
                    e_cyclic=True,
                    capital_cost=assumptions_df.at["ammonia_storage","fixed"])

        network.add("Load","ammonia_load",
                    bus="destination",
                    carrier="efuel",
                    p_set=assumptions["efuels_load"])

        network.add("Bus",
                    "nitrogen",
                    carrier="nitrogen")

        network.add("Link",
                    "air separation unit",
                    bus0="electricity",
                    bus1="nitrogen",
                    carrier="air separation unit",
                    p_nom_extendable=True,
                    efficiency=assumptions["air_separation_unit_efficiency"], # in tN2/MWhel
                    capital_cost=assumptions_df.at["air_separation_unit","fixed"]*assumptions["air_separation_unit_efficiency"]) #NB: cost is EUR/(tN2/h)

        network.add("Link",
                    "Haber Bosch",
                    bus0="hydrogen",
                    bus1="ammonia",
                    bus2="electricity",
                    bus3="nitrogen",
                    carrier="Haber-Bosch synthesis",
                    p_nom_extendable=True,
                    p_min_pu=assumptions["haber_bosch_min_part_load"]/100,
                    efficiency=assumptions["haber_bosch_efficiency"],
                    efficiency2=-assumptions["haber_bosch_electricity"]*assumptions["haber_bosch_efficiency"],
                    efficiency3=-assumptions["haber_bosch_nitrogen"]*assumptions["haber_bosch_efficiency"],
                    capital_cost=assumptions_df.at["haber_bosch","fixed"]*assumptions["haber_bosch_efficiency"]) #NB: cost is EUR/kW_NH3

    elif assumptions["efuel"] == "electricity_hvdc":

        network.add("Bus",
                    "destination",
                    carrier="efuel")

        network.add("Load","electricity_load",
                    bus="destination",
                    carrier="efuel",
                    p_set=assumptions["efuels_load"])

        network.add("Link",
                    "hvdc",
                    bus0="electricity",
                    bus1="destination",
                    carrier="hvdc",
                    p_nom_extendable=True,
                    efficiency=(1-assumptions["hvdc_submarine_cable_losses"]/100.)**(distance_transported/1000.),
                    capital_cost=assumptions_df.at["hvdc_submarine_cable","fixed"]*distance_transported + assumptions_df.at["hvdc_inverter_pair","fixed"])

        network.add("Link",
                    "hydrogen_turbine",
                    bus0="hydrogen",
                    bus1="electricity",
                    carrier="hydrogen turbine",
                    p_nom_extendable=True,
                    efficiency=assumptions["hydrogen_turbine_efficiency"]/100.,
                    capital_cost=assumptions_df.at["hydrogen_turbine","fixed"]*assumptions["hydrogen_turbine_efficiency"]/100.)  #NB: fixed cost is per MWel

    else:
        return None, None, "Efuel {} was not recognised".format(assumptions["efuel"])


    if assumptions["efuel"] in ["methanol","ammonia","methane","lh2","ft"]:
        efuel = assumptions["efuel"]
        loading_loss = 2*assumptions[efuel+"_ship_loading_loss"]/100
        transport_loss = assumptions[efuel+"_ship_energy_demand"]/assumptions[efuel+"_ship_capacity_mwh"]*(2*distance_transported)
        total_loss = loading_loss+transport_loss
        efficiency = (1-total_loss)
        print("efficiency",round(efficiency,3))

        travel_time = 2*assumptions[efuel+"_ship_loading_time"] + 2*distance_transported/assumptions[efuel+"_ship_average_speed"]
        trips_per_year = 8760/travel_time

        print("trips per year",round(trips_per_year,3))

        cost_per_t_capacity = assumptions_df.at[efuel+'_ship','fixed']/assumptions[efuel+"_ship_capacity_t"]

        cost_per_MWh = cost_per_t_capacity/config["mwh_per_t"][efuel]

        cost_per_MW = 8760*cost_per_MWh/trips_per_year
        print("cost per MW", cost_per_MW)

        network.add("Link",
                    efuel + " shipping",
                    bus0=efuel,
                    bus1="destination",
                    carrier=efuel + " shipping",
                    p_nom_extendable=True,
                    efficiency=efficiency,
                    capital_cost=cost_per_MW)


    network.consistency_check()

    solver_name = "cbc"
    solver_options = {}
    #solver_name = "gurobi"
    #solver_options = {"method": 2, # barrier
    #                  "crossover": 0}
                      #"BarConvTol": 1.e-5,
                      #"AggFill": 0,
                      #"PreDual": 0,
                      #"GURO_PAR_BARDENSETHRESH": 200}

    formulation = "kirchhoff"
    status, termination_condition = network.lopf(pyomo=False,
                                                 solver_name=solver_name,
                                                 solver_options=solver_options,
                                                 formulation=formulation,
                                                 extra_functionality=extra_functionality)

    print(status,termination_condition)

    if termination_condition in ["infeasible","infeasible or unbounded"]:
        return None, None, "Problem was infeasible"
    elif termination_condition in ["numeric"]:
        return None, None, "Numerical trouble encountered, problem could be infeasible"
    elif status == "ok" and termination_condition == "optimal":
        pass
    elif status == "warning" and termination_condition == "suboptimal":
        pass
    else:
        return None, None, "Job failed to optimise correctly"


    results_overview = pd.Series(dtype=float)
    results_overview["distance"] = distance
    results_overview["distance_transported"] = distance_transported
    results_overview["objective"] = network.objective/8760
    results_overview["average_price"] = network.buses_t.marginal_price.mean()["electricity"]
    results_overview["average_efuel_price"] = network.buses_t.marginal_price.mean()["destination"]
    if "electricity" not in assumptions["efuel"]:
        efuel = assumptions["efuel"]
        fuel = efuel if "_" not in efuel else efuel[:efuel.find("_")]
        results_overview["average_efuel_price_per_t"] = results_overview["average_efuel_price"]*config["mwh_per_t"][fuel]
    if assumptions["efuel"] in ["methanol","ft"]:
        results_overview["average_efuel_price_per_l"] = results_overview["average_efuel_price"]*config["mwh_per_m3"][fuel]/1e3
    results_overview["average_hydrogen_price"] = network.buses_t.marginal_price.mean()["hydrogen"]


    results_series = export_time_series(network)

    absmax = results_series.abs().max()

    to_drop = absmax.index[absmax < threshold*assumptions["efuels_load"]]
    results_series.drop(to_drop,
                        axis=1,
                        inplace=True)

    #results_series["electricity_price"] = network.buses_t.marginal_price["electricity"]


    results_overview["average_cost"] = sum([results_overview[s] for s in results_overview.index if s[-5:] == "_cost"])/assumptions["efuels_load"]


    stats = network.statistics(aggregate_time="sum").groupby(level=1).sum()

    stats["Total Expenditure"] = stats[["Capital Expenditure","Operational Expenditure"]].sum(axis=1)

    #exclude components contributing less than 0.1 EUR/MWh
    selection = stats.index[stats["Total Expenditure"]/assumptions["efuels_load"] > 100*threshold]
    stats = stats.loc[selection]

    for name,full_name in [("capex","Capital Expenditure"),("opex","Operational Expenditure"),("totex","Total Expenditure"),("capacity","Optimal Capacity")]:
        results_overview = pd.concat((results_overview,
                                      stats[full_name].rename(lambda x: x+ f" {name}")))

    results_overview = pd.concat((results_overview,
                                  (stats["Curtailment"]/(stats["Supply"]+stats["Curtailment"])).rename(lambda x: x+ " curtailment")))

    results_overview = pd.concat((results_overview,
                                  (stats["Total Expenditure"]/(stats["Supply"])).rename(lambda x: x+ " LCOE")))

    stats_mean = network.statistics(aggregate_time="mean").groupby(level=1).sum().loc[selection]
    results_overview = pd.concat((results_overview,
                                  stats_mean["Capacity Factor"].rename(lambda x: x+ " cf used")))
    results_overview = pd.concat((results_overview,
                                  ((stats_mean["Supply"]+stats_mean["Curtailment"])/stats_mean["Optimal Capacity"]).rename(lambda x: x+ " cf available")))

    #this is not a real capacity
    results_overview.drop("co2 capacity",
                          inplace=True,
                          errors='ignore')

    fn = 'networks/{}.nc'.format(assumptions['results_hex'])
    network.export_to_netcdf(fn)

    return results_overview, results_series, None


def solve(assumptions):

    job = get_current_job()
    jobid = job.get_id()

    job.meta['status'] = "Reading in data"
    job.save_meta()

    if assumptions["version"] != config["current_version"]:
        return error(f'Outdated version {assumptions["version"]} can no longer be calculated; please use version {config["current_version"]} instead', jobid)

    pu, error_msg = get_weather(assumptions["source_lat"], assumptions["source_lng"], assumptions["year"], assumptions['cf_exponent'])
    if error_msg is not None:
        return error(error_msg, jobid)
    pu = pu.round(3)

    #for test data stored monthly, make hourly again
    snapshots = pd.date_range("{}-01-01".format(assumptions["year"]),"{}-12-31 23:00".format(assumptions["year"]),freq="H")
    pu = pu.reindex(snapshots,method="nearest")

    carrier_series_csv = 'data/results-carrier-series-{}.csv'.format(assumptions['results_hex'])
    overview_csv = 'data/results-overview-{}.csv'.format(assumptions['results_hex'])

    print("Calculating results from scratch, saving as:", overview_csv, carrier_series_csv)
    job.meta['status'] = "Solving optimisation problem"
    job.save_meta()
    results_overview, results_series, error_msg = run_optimisation(assumptions, pu)
    if error_msg is not None:
        return error(error_msg, jobid)
    results_series = results_series.round(1)

    results_series.to_csv(carrier_series_csv)
    results_overview.to_csv(overview_csv,header=False)

    results_series.to_csv(carrier_series_csv)

    with open('data/results-assumptions-{}.json'.format(assumptions['results_hex']), 'w') as fp:
        json.dump(assumptions,fp)

    job.meta['status'] = "Processing and sending results"
    job.save_meta()

    with open('results-solve/results-{}.json'.format(jobid), 'w') as fp:
        json.dump({"jobid" : jobid,
                   "status" : "Finished",
                   "job_type" : assumptions["job_type"],
                   "average_cost" : results_overview["average_cost"],
                   "results_hex" : assumptions['results_hex']
                   },fp)

    #with open('results-{}.json'.format(job.id), 'w') as fp:
    #    json.dump(results,fp)

    return {"job_type" : "solve", "results_hex" : assumptions['results_hex']}
