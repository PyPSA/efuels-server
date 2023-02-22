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


from flask import Flask, request, jsonify, render_template, Markup

from redis import Redis

import rq
from rq.job import Job
from rq import Queue

import json, os, hashlib, yaml

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

booleans = defaults.index[defaults.type == "b"].to_list()

floats = defaults.index[defaults.type == "f"].union(defaults_t["2020"].index[defaults_t["2020"]["type"] == "f"]).to_list()

ints = defaults.index[defaults.type == "i"].to_list()

strings = defaults.index[defaults.type == "s"].to_list()

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)


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
            if assumptions[key] < 0 or assumptions[key] > config["float_upper_limit"]:
                return "{} {} was not in the valid range [0,{}]".format(key,assumptions[key],config["float_upper_limit"]), None

    for key in ints:
        try:
            assumptions[key] = int(assumptions[key])
        except:
            return "{} {} could not be converted to an integer".format(key,assumptions[key]), None

    for key in strings:
        assumptions[key] = str(assumptions[key])

    if assumptions["frequency"] < 1 or assumptions["frequency"] > 8760:
        return "Frequency {} is not in the valid range [1,8760]".format(assumptions["frequency"]), None

    if assumptions["year"] < config["years_available_start"] or assumptions["year"] > config["years_available_end"]:
        return "Year {} not in valid range".format(assumptions["year"]), None

    if assumptions["efuels_load"] == 0:
        return "No load", None

    if assumptions["efuel"] not in ["hydrogen_submarine_pipeline","methanol"]:
        return "E-fuel {} is not recognised".format(assumptions["efuel"]), None

    return None, assumptions


def compute_results_hash(assumptions):
    results_string = ""
    for item in strings+ints+booleans+floats:
        results_string += "&{}={}".format(item,assumptions[item])
    print(results_string)
    return hashlib.md5(results_string.encode()).hexdigest()


def find_results(results_hash):

    assumptions_json = f'data/results-assumptions-{results_hash}.json'
    overview_csv = f'data/results-overview-{results_hash}.csv'
    carrier_series_csv = f'data/results-carrier-series-{results_hash}.csv'

    if not os.path.isfile(assumptions_json):
        return "Assumptions file is missing", {}
    if not os.path.isfile(overview_csv):
        return "Overview results file is missing", {}
    if not os.path.isfile(carrier_series_csv):
        return "Carrier series results file is missing", {}

    print("Using preexisting results files:", assumptions_json, overview_csv, carrier_series_csv)
    with(open(assumptions_json, 'r')) as f:
        assumptions = json.load(f)
    results_overview = pd.read_csv(overview_csv,
                                   index_col=0,
                                   header=None,
                                   squeeze=True)
    carrier_series = pd.read_csv(carrier_series_csv,
                                 index_col=0,
                                 header=[0,1],
                                 parse_dates=True).round(1)

    #determine nice ordering of components
    current_order = results_overview.index[results_overview.index.str[-6:] == " totex"].str[:-6]
    preferred_order = pd.Index(config["preferred order"])
    new_order = preferred_order.intersection(current_order).append(current_order.difference(preferred_order))

    print("old:",current_order)
    print("new:",new_order)

    results = dict(results_overview)

    results["assumptions"] = assumptions

    results["order"] = list(new_order)

    results["snapshots"] = [str(s) for s in carrier_series.index]

    results["carrier_series"] = {}

    for carrier in config["balances_to_display"]:

        print("processing series for energy carrier", carrier)

        #group technologies
        df =  carrier_series[carrier]

        #sort into positive and negative
        separated = {}
        separated["positive"] = pd.DataFrame(index=df.index,
                                             dtype=float)
        separated["negative"] = pd.DataFrame(index=df.index,
                                             dtype=float)

        for col in df.columns:
            if df[col].min() > -1:
                separated["positive"][col] = df[col]
                separated["positive"][col][separated["positive"][col] < 0] = 0

            elif df[col].max() < 1:
                separated["negative"][col] = df[col]
                separated["negative"][col][separated["negative"][col] > 0] = 0

            else:
                separated["positive"][col] = df[col]
                separated["positive"][col][separated["positive"][col] < 0] = 0
                separated["negative"][col] = df[col]
                separated["negative"][col][separated["negative"][col] > 0] = 0

        separated["negative"] *= -1

        results["carrier_series"][carrier] = {}
        results["carrier_series"][carrier]["label"] = "power"
        results["carrier_series"][carrier]["units"] = "MW"

        for sign in ["positive","negative"]:
            results["carrier_series"][carrier][sign] = {}
            results["carrier_series"][carrier][sign]["columns"] = separated[sign].columns.tolist()
            results["carrier_series"][carrier][sign]["data"] = (separated[sign].values).tolist()
            print(sign,separated[sign].columns)
            results["carrier_series"][carrier][sign]["color"] = [config["colors"][i] for i in separated[sign].columns]

    return None, results


#defaults to only listen to GET and HEAD
@app.route('/')
def root():

    print("requests:",request.args)
    if "results" in request.args:
        results_hash = request.args.get("results",type=str)
        error_message, results = find_results(results_hash)
        if error_message is not None:
            print("error:",error_message)
            results = {}
    else:
        results = {}

    return render_template('index.html',
                           defaults=defaults.T.to_dict(),
                           defaults_t={year: defaults_t[year].T.to_dict() for year in defaults_t},
                           config=config,
                           results=results)


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
