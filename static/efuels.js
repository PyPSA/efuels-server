// Copyright 2022-3 Tom Brown

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation; either version 3 of the
// License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.

// License and more information at:
// https://github.com/PyPSA/efuels-server


var licenceText = '(Licence: <a href="https://creativecommons.org/licenses/by/4.0/">CC BY 4.0</a>, Attribution: <a href="https://model.energy">model.energy</a> & <a href="https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5">ECMWF ERA5</a> via <a href="https://cds.climate.copernicus.eu/">Copernicus Climate Change Service</a> (see <a href="https://apps.ecmwf.int/datasets/licences/copernicus/">weather data licence</a>))';

var parseDate = d3.timeParse("%Y-%m-%d %H:%M:00");

var formatDate = d3.timeFormat("%b %d %H:%M");


// Centered on Frankfurt
var map = L.map('mapid').setView([20,0], 2);

L.tileLayer('https://api.mapbox.com/styles/v1/{id}/tiles/{z}/{x}/{y}?access_token={accessToken}', {
    attribution: 'Map data &copy; <a href="https://www.openstreetmap.org/">OpenStreetMap</a> contributors, <a href="https://creativecommons.org/licenses/by-sa/2.0/">CC-BY-SA</a>, Imagery © <a href="https://www.mapbox.com/">Mapbox</a>',
    tileSize: 512,
    maxZoom: 18,
    zoomOffset: -1,
    id: 'mapbox/streets-v11',
    accessToken: 'pk.eyJ1IjoibndvcmJtb3QiLCJhIjoiY2prbWxibTUyMjZsMDNwcGp2bHR3OWZsaSJ9.MgSprgR6BEbBLXl5rPvXvQ'
}).addTo(map);


let balances = config["balances_to_display"];

let colors = config["colors"];


let assumptions = {};

for (let key in defaults){
    assumptions[key] = defaults[key]["value"];
}

let tech_assumptions = {};

for (let i = 0; i < config["tech_years"].length; i++){
    let year = String(config["tech_years"][i]);
    tech_assumptions[year] = {};
    for(let key in defaults_t[year]){
	tech_assumptions[year][key] = defaults_t[year][key]["value"];
    }
}


d3.select("#tech_scenario").on("change", function(){
    let scenario = this.value;
    console.log("tech scenario change to",scenario);
    for (let i = 0; i < Object.keys(tech_assumptions[scenario]).length; i++){
	let key = Object.keys(tech_assumptions[scenario])[i];
	let value = tech_assumptions[scenario][key];
	assumptions[key] = value;
	document.getElementsByName(key)[0].value = value;
    };
});




// overwrite assumptions if given already
if("assumptions" in results){
    assumptions = results["assumptions"];
};



let locations = ["destination", "source"];

for (let i=0; i<locations.length; i++){
    let location = locations[i];
    let lat = assumptions[location +"_lat"];
    let lng = assumptions[location +"_lng"];
    let marker = L.marker([lat,lng], {draggable:'true'});
    map.addLayer(marker);
    marker._icon.style.filter = (location == "source") ? "hue-rotate(120deg)" : "";
    marker.bindPopup('product ' + location);
    document.getElementsByName(location+"_location")[0].value = printLocation(lat,lng);
    marker.on('dragend', function(){
	lat = parseFloat(marker.getLatLng().lat.toFixed(1));
	lng = parseFloat(marker.getLatLng().lng.toFixed(1));
	assumptions[location +"_lat"] = lat;
	assumptions[location +"_lng"] = lng;
	document.getElementsByName(location+"_location")[0].value = printLocation(lat,lng);
	console.log(location, 'changed location to: lat:', lat, ", lon:", lng);
    });
};



function printLocation(lat,lng){
    return "latitude: " + lat + ", longitude: " + lng;
};




for (let i = 0; i < Object.keys(assumptions).length; i++){
    let key = Object.keys(assumptions)[i];
    let value = assumptions[key];
    console.log(key,value);
    if(["job_type","source_lat","source_lng","destination_lat","destination_lng","version","jobid","timestamp","queue_length","weather_hex","results_hex","cf_exponent","overland_fraction"].includes(key)){
    }
    else if(typeof value === "boolean"){
	document.getElementsByName(key)[0].checked = value;
	d3.selectAll("input[name='" + key + "']").on("change", function(){
	    assumptions[key] = this.checked;
	    console.log(key,"changed to",assumptions[key]);
	});
    }
    else if(typeof value === "string"){
	document.getElementsByName(key)[0].value = value;
	d3.selectAll("select[name='" + key + "']").on("change", function(){
	    assumptions[key] = this.value;
	    console.log(key,"changed to",assumptions[key]);
	});
    }
    else{
	document.getElementsByName(key)[0].value = value;
	d3.selectAll("input[name='" + key + "']").on("change", function(){
	    assumptions[key] = this.value;
	    console.log(key,"changed to",assumptions[key]);
	});
    }
}


var solveButton = d3.select("#solve-button");

var solveButtonText = {"before" : "Solve",
		       "after" : "Solving"};

var jobid = "";

var timer;
var timeout;
var timerStart;
var timerExpected = config["expected_time"];


// time between status polling in milliseconds
var poll_interval = 1000*config["poll_interval"];

// Shouldn't be divisible by poll_interval
var poll_timeout = 1000*config["poll_timeout"] + poll_interval/2;



function solve() {
    if (solveButton.text() == solveButtonText["before"]) {
	clear_results();
	var send_job = new XMLHttpRequest();
	send_job.open('POST', './jobs', true);
	send_job.setRequestHeader("Content-Type", "application/json");
	send_job.onload = function () {
	    var data = JSON.parse(this.response);
	    if("jobid" in data){
		jobid = data["jobid"];
		console.log("Jobid:",jobid);
		timer = setInterval(poll_result, poll_interval);
		timerStart = new Date().getTime();
		document.getElementById("countdown").innerHTML="Ready in around " + timerExpected + " seconds";
		console.log("timer",timer,"polling every",poll_interval,"milliseconds");
		timeout = setTimeout(poll_kill, poll_timeout);
		solveButton.text(solveButtonText["after"]);
		solveButton.attr("disabled","");
		document.getElementById("status").innerHTML="Sending job to solver";
	    } else if("status" in data && data["status"] == "Error") {
		console.log("results:", data);
		document.getElementById("status").innerHTML = data["status"] + ": " + data["error"];
	    } else {
		console.log("results:", data);
		document.getElementById("status").innerHTML = "Found previously calculated version";
		results = data;
		display_results();
	    };
	};
	assumptions["job_type"] = "solve";
	send_job.send(JSON.stringify(assumptions));
    };
};


solveButton.on("click", solve);


if(assumptions["job_type"] === "solve"){
    display_results();
};


function poll_result() {

    var poll = new XMLHttpRequest();

    poll.open('GET', './jobs/' + jobid, true);

    poll.onload = function () {
	results = JSON.parse(this.response);
	status = results["status"];
	document.getElementById("status").innerHTML=status;
	console.log("status is",status);

	document.getElementById("countdown").innerHTML = "Ready in around " + Math.round(timerExpected - (new Date().getTime() - timerStart)/1000.) + " seconds";

	if(status == "Error"){
	    clearInterval(timer);
	    clearTimeout(timeout);
	    document.getElementById("countdown").innerHTML = "Ready in around " + timerExpected + " seconds";
	    console.log("results:",results);
	    document.getElementById("status").innerHTML=status + ": " + results["error"];
	    solveButton.text(solveButtonText["before"]);
	    $('#solve-button').removeAttr("disabled");
	};
	if(status == "Finished"){
	    clearInterval(timer);
	    clearTimeout(timeout);
	    document.getElementById("countdown").innerHTML = "Solved in " + Math.round((new Date().getTime() - timerStart)/1000.) + " seconds";
	    console.log("results:",results);
	    solveButton.text(solveButtonText["before"]);
	    $('#solve-button').removeAttr("disabled");
	    display_results();
	};
    };
    poll.send();
};


function poll_kill() {
    clearInterval(timer);
    solveButton.text(solveButtonText["before"]);
    $('#solve-button').removeAttr("disabled");
    document.getElementById("status").innerHTML="Error: Timed out";
};




function clear_results(){
    document.getElementById("results_assumptions").innerHTML="";
    document.getElementById("average_cost").innerHTML="";
    document.getElementById("average_cost_per_t").innerHTML="";
    document.getElementById("average_cost_per_l").innerHTML="";
    document.getElementById("distance").innerHTML="";
    document.getElementById("distance_transported").innerHTML="";

    var table = document.getElementById('results_table');
    var rowCount = table.rows.length;
    for (let i=rowCount-1; i>0; i--){
	table.deleteRow(i);
    };

    d3.select("#average_cost_graph").selectAll("g").remove();
    d3.select("#power_capacity_bar").selectAll("g").remove();
    d3.select("#energy_capacity_bar").selectAll("g").remove();


    for (var k=0; k < balances.length; k++){
	let balance = balances[k];
	d3.select("#" + balance + "_power_graph").selectAll("g").remove();
	d3.select("#" + balance + "_power_graph_legend").selectAll("g").remove();
    };

    document.getElementById("results-overview-download").innerHTML = '';
    document.getElementById("results-series-download").innerHTML = '';
    document.getElementById("results-netcdf-download").innerHTML = '';
    document.getElementById("results-link").innerHTML = '';

};




function fill_results_table(){
    var table = document.getElementById('results_table');


    for(let i=0; i<results["order"].length; i++){
	let asset = results["order"][i];
	// escape for co2
	if(!(asset + " capacity" in results)){
	    continue;
	};
	let row = table.insertRow(table.rows.length);
	let cap = results[asset + " capacity"].toFixed(1) + " MW";
	if (asset.includes("storage")) cap += "h";
	if (asset.includes("co2 storage")) cap = cap.slice(0,-3) + "tCO2";
	let cfUsed = (100*results[asset + " cf used"]).toFixed(1);
	let cfAvailable = "";
	if (["wind","solar"].includes(asset)) cfAvailable = (100*results[asset + " cf available"]).toFixed(1);
	let curtailment = "";
	if (["wind","solar"].includes(asset)) curtailment = (100*results[asset + " curtailment"]).toFixed(1);
	let lcoe = "";
	if (["wind","solar"].includes(asset)) lcoe = results[asset + " LCOE"].toFixed(1);

	let content = [asset,
		       cap,
		       cfUsed,
		       cfAvailable,
		       curtailment,
		       lcoe];

	for(let j=0; j<content.length;j++){
	    let cell = row.insertCell(j);
	    cell.innerHTML = content[j];
	    if(j == 0){
		cell.className="tab_asset";
	    } else {
		cell.className="tab_data";
	    };
	};
    };
};


function display_results(){

    $('#collapseResults').addClass("show");

    document.getElementById("results_assumptions").innerHTML=" for " + config["efuels"][assumptions["efuel"]];
    document.getElementById("average_cost").innerHTML="Average cost [EUR/MWh]: " + results["average_efuel_price"].toFixed(1);
    document.getElementById("average_cost_per_t").innerHTML= "average_efuel_price_per_t" in results ? "Average cost [EUR/tonne]: " + results["average_efuel_price_per_t"].toFixed(1) : "";
    document.getElementById("average_cost_per_l").innerHTML= "average_efuel_price_per_l" in results ? "Average cost [EUR/litre]: " + results["average_efuel_price_per_l"].toFixed(2) : "";
    document.getElementById("distance").innerHTML="Distance as crow flies [km]: " + results["distance"].toFixed(0);
    document.getElementById("distance_transported").innerHTML="Distance transported [km]: " + results["distance_transported"].toFixed(0);


    for(var j=0; j < results.snapshots.length; j++) {
	results.snapshots[j] = parseDate(results.snapshots[j]);
    };

    fill_results_table();
    draw_power_capacity_bar();
    draw_energy_capacity_bar();
    draw_cost_stack();

    for (var k=0; k < balances.length; k++){
	let balance = balances[k];
	if(balance in results["carrier_series"]){
	    console.log("Drawing power time series for", balance);
	    let series = results["carrier_series"][balance];
	    draw_series(series, results.snapshots, balance);
	} else {
	    console.log("Balance data not available for", balance);
	};
    };


    document.getElementById("results-overview-download").innerHTML = '<a href="data/results-overview-' + results.assumptions.results_hex + '.csv">Download Comma-Separated-Variable (CSV) file of results overview</a> ' + licenceText;
    document.getElementById("results-series-download").innerHTML = '<a href="data/results-carrier-series-' + results.assumptions.results_hex + '.csv">Download Comma-Separated-Variable (CSV) file of results time series</a> ' + licenceText;
    document.getElementById("results-netcdf-download").innerHTML = '<a href="networks/' + results.assumptions.results_hex + '.nc">Download PyPSA NetCDF file</a> ' + licenceText;
    document.getElementById("results-link").innerHTML = '<a href="./?results=' + results.assumptions.results_hex + '#solve">Link to these results</a>';
};



// taken from pypsa-server and renamed results -> series to avoid name clash
function draw_series(series, snapshots, balance){

    let selection = [...Array(snapshots.length).keys()];

    // Inspired by https://bl.ocks.org/mbostock/3885211

    var svgGraph = d3.select("#" + balance + "_power_graph"),
	margin = {top: 20, right: 20, bottom: 110, left: 50},
	marginContext = {top: 430, right: 20, bottom: 30, left: 50},
	width = svgGraph.attr("width") - margin.left - margin.right,
	height = svgGraph.attr("height") - margin.top - margin.bottom,
	heightContext = svgGraph.attr("height") - marginContext.top - marginContext.bottom;

    // remove existing
    svgGraph.selectAll("g").remove();

    var x = d3.scaleTime().range([0, width]).domain(d3.extent(snapshots));
    var y = d3.scaleLinear().range([height, 0]);
    var xContext = d3.scaleTime().range([0, width]).domain(d3.extent(snapshots));
    var yContext = d3.scaleLinear().range([heightContext, 0]);

    var xAxis = d3.axisBottom(x),
	xAxisContext = d3.axisBottom(xContext),
	yAxis = d3.axisLeft(y);


    var brush = d3.brushX()
        .extent([[0, 0], [width, heightContext]])
        .on("start brush end", brushed);


    var zoom = d3.zoom()
        .scaleExtent([1, Infinity])
        .translateExtent([[0, 0], [width, height]])
        .extent([[0, 0], [width, height]])
        .on("zoom", zoomed);


    var data = [];

    // Custom version of d3.stack

    var previous = new Array(selection.length).fill(0);

    for (var j = 0; j < series["positive"].columns.length; j++){
	var item = [];
	for (var k = 0; k < selection.length; k++){
	    item.push([previous[k], previous[k] + series["positive"]["data"][selection[k]][j]]);
	    previous[k] = previous[k] + series["positive"]["data"][selection[k]][j];
	    }
	data.push(item);
    }
    var previous = new Array(selection.length).fill(0);

    for (var j = 0; j < series["negative"].columns.length; j++){
	var item = [];
	for (var k = 0; k < selection.length; k++){
	    item.push([-previous[k] - series["negative"]["data"][selection[k]][j],-previous[k]]);
	    previous[k] = previous[k] + series["negative"]["data"][selection[k]][j];
	    }
	data.push(item);
    }

    var ymin = 0, ymax = 0;
    for (var k = 0; k < selection.length; k++){
	if(data[series["positive"].columns.length-1][k][1] > ymax){ ymax = data[series["positive"].columns.length-1][k][1];};
	if(data[series["positive"].columns.length+series["negative"].columns.length-1][k][0] < ymin){ ymin = data[series["positive"].columns.length+series["negative"].columns.length-1][k][0];};
    };

    y.domain([ymin,ymax]);
    yContext.domain([ymin,ymax]);

    var area = d3.area()
        .curve(d3.curveMonotoneX)
        .x(function(d,i) { return x(snapshots[i]); })
        .y0(function(d) { return y(d[0]); })
        .y1(function(d) { return y(d[1]); });

    var areaContext = d3.area()
        .curve(d3.curveMonotoneX)
        .x(function(d,i) { return xContext(snapshots[i]); })
        .y0(function(d) { return yContext(d[0]); })
        .y1(function(d) { return yContext(d[1]); });


    svgGraph.append("defs").append("clipPath")
        .attr("id", "clip")
	.append("rect")
        .attr("width", width)
        .attr("height", height);

    var focus = svgGraph.append("g")
        .attr("class", "focus")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var context = svgGraph.append("g")
        .attr("class", "context")
        .attr("transform", "translate(" + marginContext.left + "," + marginContext.top + ")");

    var layer = focus.selectAll(".layer")
        .data(data)
        .enter().append("g")
        .attr("class", "layer");

    layer.append("path")
        .attr("class", "area")
        .style("fill", function(d, i) {if(i < series["positive"].color.length){ return series["positive"].color[i];} else{return series["negative"].color[i-series["positive"].color.length];} })
        .attr("d", area);


    focus.append("g")
        .attr("class", "axis axis--x")
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis);

    focus.append("g")
        .attr("class", "axis axis--y")
        .call(yAxis);


    var label = svgGraph.append("g").attr("class", "y-label");

    // text label for the y axis
    label.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 0)
        .attr("x",0 - (height / 2))
        .attr("dy", "1em")
        .style("text-anchor", "middle")
        .text(series["label"] + " [" + series["units"] + "]");


    var layerContext = context.selectAll(".layerContext")
        .data(data)
        .enter().append("g")
        .attr("class", "layerContext");

    layerContext.append("path")
        .attr("class", "area")
        .style("fill", function(d, i) {if(i < series["positive"].color.length){ return series["positive"].color[i];} else{return series["negative"].color[i-series["positive"].color.length];} })
        .attr("d", areaContext);

    context.append("g")
        .attr("class", "axis axis--x")
        .attr("transform", "translate(0," + heightContext + ")")
        .call(xAxisContext);

    var gBrush = context.append("g")
        .attr("class", "brush")
        .call(brush);

    // brush handle follows
    // https://bl.ocks.org/Fil/2d43867ba1f36a05459c7113c7f6f98a

    // following handle looks nicer
    // https://bl.ocks.org/robyngit/89327a78e22d138cff19c6de7288c1cf

    var brushResizePath = function(d) {
	var e = +(d.type == "e"),
	    x = e ? 1 : -1,
	    y = heightContext / 2;
	return "M" + (.5 * x) + "," + y + "A6,6 0 0 " + e + " " + (6.5 * x) + "," + (y + 6) + "V" + (2 * y - 6) + "A6,6 0 0 " + e + " " + (.5 * x) + "," + (2 * y) + "Z" + "M" + (2.5 * x) + "," + (y + 8) + "V" + (2 * y - 8) + "M" + (4.5 * x) + "," + (y + 8) + "V" + (2 * y - 8);
    }

    var handle = gBrush.selectAll(".handle--custom")
	.data([{type: "w"}, {type: "e"}])
	.enter().append("path")
        .attr("class", "handle--custom")
        .attr("stroke", "#000")
        .attr("cursor", "ew-resize")
        .attr("d", brushResizePath);

    gBrush.call(brush.move, x.range()); //this sets initial position of brush


    svgGraph.append("rect")
        .attr("class", "zoom")
        .attr("width", width)
        .attr("height", height)
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
        .call(zoom);

    function brushed() {
	if (d3.event.sourceEvent && d3.event.sourceEvent.type === "zoom") return; // ignore brush-by-zoom
	var s = d3.event.selection || xContext.range();
	x.domain(s.map(xContext.invert, xContext));
	layer.attr("d", area);
	focus.select(".axis--x").call(xAxis);
	svgGraph.select(".zoom").call(zoom.transform, d3.zoomIdentity
				      .scale(width / (s[1] - s[0]))
				      .translate(-s[0], 0));
	handle.attr("transform", function(d, i) { return "translate(" + [ s[i], - heightContext / 4] + ")"; });
    };

    function zoomed() {
	if (d3.event.sourceEvent && d3.event.sourceEvent.type === "brush") return; // ignore zoom-by-brush
	var t = d3.event.transform;
	x.domain(t.rescaleX(xContext).domain());
	layer.select(".area").attr("d", area);
	focus.select(".axis--x").call(xAxis);
	var newRange = x.range().map(t.invertX, t);
	context.select(".brush").call(brush.move, newRange);
	handle.attr("transform", function(d, i) { return "translate(" + [ newRange[i], - heightContext / 4] + ")"; });
    };


    //Legend

    //slice to make copy
    let labels = series["positive"]["columns"].slice().reverse().concat(series["negative"]["columns"]);
    let color = series["positive"]["color"].slice().reverse().concat(series["negative"]["color"]);

    //remove duplicate labels
    console.log(balance, labels);
    let toRemove = [];
    for(var i=0; i < labels.length; i++){
	let label = labels[i];
	if(i != labels.indexOf(label)) {
	    console.log(i, label, labels.indexOf(label));
	    toRemove.push(i);
	};
    };

    console.log("need to remove",toRemove);

    for(var i=0; i < toRemove.length; i++){
	// have to subtract i because array shifts to left each time
	labels.splice(toRemove[i]-i, 1);
	color.splice(toRemove[i]-i, 1);
    };

    let legendSVG = d3.select("#" + balance + "_power_graph_legend");

    let legend = legendSVG.selectAll("g")
	.data(labels)
	.enter()
	.append("g")
	.attr("transform", function (d, i) {  return "translate(0," + (5 + i * 20) + ")" });

    legend.append("rect")
	.attr("x",0)
	.attr("y",0)
	.attr("width", 10)
	.attr("height", 10)
	.style("fill", function (d, i) { return color[i] });

    legend.append("text")
	.attr("x",20)
	.attr("y",10)
	.text(function (d, i) { return d});

};



function draw_cost_stack(){

    let data = [];
    let color = [];
    let labels = [];

    for(let i=0; i<results["order"].length; i++){
	let asset = results["order"][i];
	data.push(results[asset + " totex"]/results["assumptions"]["efuels_load"]/8760.);
	color.push(colors[asset]);
	labels.push(asset);
    };

    draw_stack(data, labels, color, "Breakdown of average cost [EUR/MWh]", "#average_cost_graph", " EUR/MWh");
};


function draw_power_capacity_bar(){

    let data = [];
    let color = [];
    let labels = [];
    let units = [];

    data.push(results["assumptions"]["efuels_load"]);
    color.push(colors["efuels_load"]);
    labels.push("product demand");
    units.push("MW");

    for(let i=0; i<results["order"].length; i++){
	let asset = results["order"][i];
	if(!asset.includes("storage") && (asset + " capacity" in results)){
	    data.push(results[asset + " capacity"]);
	    color.push(colors[asset]);
	    labels.push(asset.replace("hydrogen","H2"));
	    units.push(asset.includes("co2") ? "tCO2/h" : "MW");
	};
    };

    draw_bar(data, labels, color, units, "Power capacity [MW]", "#power_capacity_bar");
};


function draw_energy_capacity_bar(){

    let data = [];
    let color = [];
    let labels = [];
    let units = [];

    data.push(results["assumptions"]["efuels_load"]*24/1000.);
    color.push(colors["efuels_load"]);
    labels.push("24h product demand");
    units.push("GWh");

    for(let i=0; i<results["order"].length; i++){
	let asset = results["order"][i];
	if(asset.includes("storage")){
	    data.push(results[asset + " capacity"]/1e3);
	    color.push(colors[asset]);
	    labels.push(asset.replace("hydrogen","H2"));
	    units.push(asset.includes("co2") ? "ktCO2" : "GWh");
	};
    };

    draw_bar(data, labels, color, units, "Storage capacity [GWh or ktCO2]", "#energy_capacity_bar");
};


function draw_stack(data, labels, color, ylabel, svgName, suffix){

    // Inspired by https://bl.ocks.org/mbostock/3885211 and
    // https://bl.ocks.org/mbostock/1134768

    let totals = [0.];

    for(let i=0; i < data.length; i++){
	totals.push(totals[i] + data[i]);
    };


    let svgGraph = d3.select(svgName),
	margin = {top: 20, right: 20, bottom: 30, left: 50},
	width = svgGraph.attr("width") - margin.left - margin.right,
	height = svgGraph.attr("height") - margin.top - margin.bottom;

    // remove existing
    svgGraph.selectAll("g").remove();
    let x = d3.scaleLinear().range([0, width]);
    let y = d3.scaleLinear().range([height, 0]);

    x.domain([0,1]);
    y.domain([0,totals[totals.length-1]]).nice();

    var g = svgGraph.append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    let tip = d3.tip()
	.attr('class', 'd3-tip')
	.offset([-8, 0])
	.html(function(d,i) {
	    return labels[i] + ": " + Math.abs(data[i]).toFixed(1) + suffix;
	});
    svgGraph.call(tip);

    var layer = g.selectAll("rect")
        .data(data)
        .enter().append("rect")
	.attr("x", x(0.1))
        .attr("y", function(d,i) { return y(totals[i+1]);})
        // following abs avoids rect with negative height e.g. -1e10
	.attr("height", function(d,i) { return Math.abs((y(totals[i]) - y(totals[i+1])).toFixed(2)); })
    	.attr("width", x(0.8))
        .style("fill", function(d, i) { return color[i];})
        .on('mouseover', tip.show)
        .on('mouseout', tip.hide);

    //g.append("g")
    //    .attr("class", "axis axis--x")
    //    .attr("transform", "translate(0," + height + ")")
    //    .call(d3.axisBottom(x));

    g.append("g")
        .attr("class", "axis axis--y")
        .call(d3.axisLeft(y));

    var label = svgGraph.append("g").attr("class", "y-label");

    // text label for the y axis
    label.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 0)
        .attr("x",0 - (height / 2))
        .attr("dy", "1em")
        .style("text-anchor", "middle")
        .text(ylabel);


    var label = svgGraph.append("g").attr("class", "column-total")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");


    // text label for the y axis
    label.append("text")
        .attr("y", y(totals[totals.length-1])-25)
        .attr("x",x(0.5))
        .attr("dy", "1em")
        .style("text-anchor", "middle")
        .text(totals[totals.length-1].toFixed(1));


};





function draw_bar(data, labels, color, units, ylabel, svgName){

    let svgGraph = d3.select(svgName),
	margin = {top: 20, right: 20, bottom: 30, left: 50},
	width = svgGraph.attr("width") - margin.left - margin.right,
	height = svgGraph.attr("height") - margin.top - margin.bottom;

    // remove existing
    svgGraph.selectAll("g").remove();
    let x = d3.scaleBand().range([0, width]).padding(0.1);
    let y = d3.scaleLinear().range([height,0]);

    x.domain(labels);
    y.domain([0,d3.max(data)]).nice();

    var g = svgGraph.append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    let tip = d3.tip()
	.attr('class', 'd3-tip')
	.offset([-8, 0])
	.html(function(d,i) {
	    return labels[i] + ": " + Math.abs(data[i]).toFixed(1) + " " + units[i];
	});
    svgGraph.call(tip);

    var layer = g.selectAll("rect")
        .data(data)
        .enter().append("rect")
	.attr("x", function(d,i) {return x(labels[i]);})
        .attr("y", function(d,i) { return Math.abs(y(d).toFixed(2));})
        // following abs avoids rect with negative height e.g. -1e10
	.attr("height", function(d,i) { return Math.abs((height -y(d)).toFixed(2)); })
    	.attr("width", x.bandwidth())
        .style("fill", function(d, i) { return color[i];})
        .on('mouseover', tip.show)
        .on('mouseout', tip.hide);

    //g.append("g")
    //    .attr("class", "axis axis--x")
    //    .attr("transform", "translate(0," + height + ")")
    //    .call(d3.axisBottom(x));

    g.append("g")
        .attr("class", "axis axis--y")
        .call(d3.axisLeft(y));

    var label = svgGraph.append("g").attr("class", "y-label");

    // text label for the y axis
    label.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 0)
        .attr("x",0 - (height / 2))
        .attr("dy", "1em")
        .style("text-anchor", "middle")
        .text(ylabel);

    // add the x Axis
    // https://stackoverflow.com/questions/20947488/d3-grouped-bar-chart-how-to-rotate-the-text-of-x-axis-ticks
    svgGraph.append("g")
        .attr("transform", "translate(" + margin.left + "," + (height + margin.top) + ")")
        .call(d3.axisBottom(x))
	.selectAll("text")
        .style("text-anchor", "middle")
        .attr("dx", "0em")
        .attr("dy", "1em")
        .attr("transform", "rotate(-10)");

};
