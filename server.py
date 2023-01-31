## Copyright 2022 Tom Brown

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


from flask import Flask, request, jsonify, render_template, Markup

from redis import Redis

import rq
from rq.job import Job
from rq import Queue

import json, os, hashlib

import pandas as pd

import datetime


app = Flask(__name__)
app.jinja_env.filters['json'] = lambda v: Markup(json.dumps(v))

conn = Redis.from_url('redis://')

queue = Queue('efuels', connection=conn)


#na_filter leaves "" as "" rather than doing nan which confuses jinja2
defaults = pd.read_csv("defaults.csv",index_col=[0,1],na_filter=False)

for (n,t) in [("f",float),("i",int)]:
    defaults.loc[defaults["type"] == n, "value"] = defaults.loc[defaults["type"] == n,"value"].astype(t)

#work around fact bool("False") returns True
defaults.loc[defaults.type == "b","value"] = (defaults.loc[defaults.type == "b","value"] == "True")

defaults_t = {year: defaults.swaplevel().loc[year] for year in ["2020","2030","2050"]}
defaults = defaults.swaplevel().loc[""]

booleans = ["wind","solar","battery","hydrogen","dispatchable1","dispatchable2"]

floats = ["cf_exponent","efuels_load","wind_cost","solar_cost","battery_energy_cost",
          "source_lat","source_lng","destination_lat","destination_lng",
          "wind_fom","wind_lifetime","wind_discount",
          "solar_fom","solar_lifetime","solar_discount",
          "battery_energy_fom","battery_energy_lifetime","battery_energy_discount",
          "battery_power_fom","battery_power_lifetime","battery_power_discount",
          "battery_power_efficiency_charging","battery_power_efficiency_discharging",
          "hydrogen_energy_fom","hydrogen_energy_lifetime","hydrogen_energy_discount",
          "hydrogen_electrolyser_fom","hydrogen_electrolyser_lifetime","hydrogen_electrolyser_discount",
          "hydrogen_turbine_fom","hydrogen_turbine_lifetime","hydrogen_turbine_discount",
          "battery_power_cost","hydrogen_electrolyser_cost",
          "hydrogen_energy_cost",
          "hydrogen_electrolyser_efficiency",
          "hydrogen_turbine_cost",
          "hydrogen_turbine_efficiency",
	  "hydrogen_submarine_pipeline_cost",
	  "hydrogen_submarine_pipeline_losses",
	  "hydrogen_submarine_pipeline_fom",
	  "hydrogen_submarine_pipeline_lifetime",
	  "hydrogen_submarine_pipeline_discount",
          "dispatchable1_cost",
          "dispatchable1_marginal_cost",
          "dispatchable1_emissions",
          "dispatchable1_discount",
          "dispatchable1_fom","dispatchable1_lifetime",
          "dispatchable2_cost",
          "dispatchable2_marginal_cost",
          "dispatchable2_emissions",
          "dispatchable2_discount",
          "dispatchable2_fom","dispatchable2_lifetime",
          "wind_min",
          "solar_min",
          "wind_max",
          "solar_max"]


ints = ["year","frequency","version"]

strings = ["efuel"]


colors = {"wind":"#3B6182",
          "solar" :"#FFFF00",
          "battery" : "#999999",
          "battery_charge" : "#999999",
          "battery_discharge" : "#999999",
          "battery_power" : "#999999",
          "battery_energy" : "#666666",
          "hydrogen_turbine" : "red",
          "hydrogen_electrolyser" : "cyan",
          "hydrogen_energy" : "magenta",
          "dispatchable1" : "orange",
          "dispatchable2" : "lime",
}


years_available_start = 2011
years_available_end = 2012

float_upper_limit = 1e7



def sanitise_assumptions(assumptions):
    """
    Fix types of assumptions and check they are in correct
    range.

    Parameters
    ----------
    assumptions : dict
        Assumptions (location, technical and economic parameters)

    Returns
    -------
    error_message : None or string
        If there was an error, details of the error
    assumptions : dict
        If there was no error, the clean type-safe assumptions
    """
    for key in strings+ints+booleans+floats:
        if key not in assumptions:
            return f"{key} missing from assumptions", None

    for key in booleans:
        try:
            assumptions[key] = bool(assumptions[key])
        except:
            return "{} {} could not be converted to boolean".format(key,assumptions[key]), None

    for key in floats:
        try:
            assumptions[key] = float(assumptions[key])
        except:
            return "{} {} could not be converted to float".format(key,assumptions[key]), None

        if "lat" not in key and "lng" not in key:
            if assumptions[key] < 0 or assumptions[key] > float_upper_limit:
                return "{} {} was not in the valid range [0,{}]".format(key,assumptions[key],float_upper_limit), None

    for key in ints:
        try:
            assumptions[key] = int(assumptions[key])
        except:
            return "{} {} could not be converted to an integer".format(key,assumptions[key]), None

    for key in strings:
        assumptions[key] = str(assumptions[key])

    if assumptions["frequency"] < 1 or assumptions["frequency"] > 8760:
        return "Frequency {} is not in the valid range [1,8760]".format(assumptions["frequency"]), None

    if assumptions["year"] < years_available_start or assumptions["year"] > years_available_end:
        return "Year {} not in valid range".format(assumptions["year"]), None

    if assumptions["efuels_load"] == 0:
        return "No load", None

    if assumptions["efuel"] not in ["hydrogen_submarine_pipeline"]:
        return f"E-fuel {efuel} is not recognised", None

    return None, assumptions


def compute_results_hash(assumptions):
    results_string = ""
    for item in strings+ints+booleans+floats:
        results_string += "&{}={}".format(item,assumptions[item])
    print(results_string)
    return hashlib.md5(results_string.encode()).hexdigest()


def find_results(results_hash):

    assumptions_json = f'data/results-assumptions-{results_hash}.json'
    series_csv = f'data/results-series-{results_hash}.csv'
    overview_csv = f'data/results-overview-{results_hash}.csv'

    if not os.path.isfile(assumptions_json):
        return "Assumptions file is missing", {}
    if not os.path.isfile(series_csv):
        return "Series results file is missing", {}
    if not os.path.isfile(overview_csv):
        return "Overview results file is missing", {}

    print("Using preexisting results files:", assumptions_json, series_csv, overview_csv)
    with(open(assumptions_json, 'r')) as f:
        assumptions = json.load(f)
    results_overview = pd.read_csv(overview_csv,
                                   index_col=0,
                                   header=None,
                                   squeeze=True)
    results_series = pd.read_csv(series_csv,
                                 index_col=0,
                                 parse_dates=True)

    #fill in old results before dispatchable
    for i in range(1,3):
        g = "dispatchable" + str(i)
        if not assumptions[g]:
            results_overview[g+"_capacity"] = 0.
            results_overview[g+"_cost"] = 0.
            results_overview[g+"_marginal_cost"] = 0.
            results_overview[g+"_used"] = 0.
            results_overview[g+"_cf_used"] = 0.
            results_overview[g+"_rmv"] = 0.
            results_series[g] = 0.

    results = dict(results_overview)

    results["assumptions"] = assumptions

    results["snapshots"] = [str(s) for s in results_series.index]

    columns = {"positive" : ["wind","solar","battery_discharge","hydrogen_turbine","dispatchable1","dispatchable2"],
               "negative" : ["battery_charge","hydrogen_electrolyser"]}


    for sign, cols in columns.items():
        results[sign] = {}
        results[sign]["columns"] = cols
        results[sign]["data"] = results_series[cols].values.tolist()
        results[sign]["color"] = [colors[c] for c in cols]

    balance = results_series[columns["positive"]].sum(axis=1) - results_series[columns["negative"]].sum(axis=1)

    print(balance.describe())

    return None, results


#defaults to only listen to GET and HEAD
@app.route('/')
def root():
    return render_template('index.html',
                           defaults=defaults.T.to_dict(),
                           defaults_t={year: defaults_t[year].T.to_dict() for year in defaults_t})


@app.route('/jobs', methods=['GET','POST'])
def jobs_api():
    if request.method == "POST":
        if request.headers.get('Content-Type','missing') != 'application/json':
            return jsonify({"status" : "Error", "error" : "No JSON assumptions sent."})

        print(request.json)

        error_message, assumptions = sanitise_assumptions(request.json)

        if error_message is not None:
            return jsonify({"status" : "Error", "error" : error_message})

        assumptions["results_hex"] = compute_results_hash(assumptions)
        error_message, results = find_results(assumptions["results_hex"])
        if error_message is None:
            assumptions["timestamp"] = str(datetime.datetime.now())
            assumptions["jobid"] = hashlib.md5(assumptions["timestamp"].encode()).hexdigest()
            assumptions["queue_length"] = 0
            with open('assumptions/assumptions-{}.json'.format(assumptions["jobid"]), 'w') as fp:
                json.dump(assumptions, fp)
            mini_results = {"jobid" : assumptions["jobid"], "status" : "Finished",
                            "error" : None, "average_cost" : results["average_cost"]}
            with open('results/results-{}.json'.format(assumptions["jobid"]), 'w') as fp:
                json.dump(mini_results, fp)
            return jsonify(results)

        job = queue.enqueue("solve.solve", args=(assumptions,), job_timeout=300)
        result = {"jobid" : job.get_id()}
        assumptions.update({"jobid" : result["jobid"],
                            "timestamp" : str(datetime.datetime.now()),
                            "queue_length" : len(queue.jobs)})
        with open('assumptions/assumptions-{}.json'.format(assumptions["jobid"]), 'w') as fp:
            json.dump(assumptions, fp)
        print("jobid {} request:".format(result["jobid"]))
        print(assumptions)
        return jsonify(result)
    elif request.method == "GET":
        return "jobs in queue: {}".format(len(queue.jobs))

@app.route('/jobs/<jobid>')
def jobid_api(jobid):
    try:
        job = Job.fetch(jobid, connection=conn)
    except:
        return jsonify({"status" : "Error", "error" : "Failed to find job!"})

    if job.is_failed:
        return jsonify({"status" : "Error", "error" : "Job failed."})

    try:
        status = job.meta['status']
    except:
        status = "Waiting for job to run (jobs in queue: {})".format(len(queue.jobs))

    result = {"status" : status}

    if job.is_finished:
        if "error" in job.result:
            result["status"] = "Error"
            result["error"] = job.result["error"]
        else:
            result["status"] = "Finished"

            if job.result["job_type"] == "weather":
                error_message, results = find_weather(job.result["weather_hex"])
                if error_message is not None:
                    result["status"] = "Error"
                    result["error"] = error_message
            elif job.result["job_type"] == "solve":
                error_message, results = find_results(job.result["results_hex"])
                if error_message is not None:
                    result["status"] = "Error"
                    result["error"] = error_message

        if result["status"] == "Finished":
            results.update(result)
            result = results

        mini_results = {"jobid" : jobid,
                        "status" : result["status"],
                        "error" : result.get("error",None),
                        "average_cost" : result.get("average_cost",None)}

        print("jobid {} results:".format(jobid))
        print(mini_results)
        with open('results/results-{}.json'.format(jobid), 'w') as fp:
            json.dump(mini_results, fp)

    return jsonify(result)



if __name__ == '__main__':
    app.run(port='5002')
