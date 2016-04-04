var EventEmitter = require('events').EventEmitter;
var d3 = require('d3');
var Globals = require('./globals.js')
var ANIMATION_DURATION = Globals.ANIMATION_DURATION;

var MuPlot = {};

MuPlot.mu = [];

MuPlot.padding = 30;
MuPlot.left_padding = 50;

MuPlot.create = function(el, props, state){
    console.log("CREATING MU PLOT")
    var svg = d3.select(el).append('svg')
      .attr('class', 'd3')
      .attr('width', props.width)
      .attr('height', props.height);
    
    svg.append('g')
      .attr('class', 'd3-mu-axes')

    svg.append('g')
      .attr('class', 'd3-mu-points')

    this.mu = props.mu
    
    this._update_points();
    var dispatcher = new EventEmitter();
    this.update(el, this.mu, state, dispatcher);
    return dispatcher;

};

MuPlot._update_points = function(){
    this.points = [];
    var mus= this.mu
    this.points = mus
        .map(function(d){return (
            {   x: d.numdate_given,
                y: d.xValue,
                name:d.name
            })
        }).filter(function(d){return d.numdate_given != 0.0 && d.xValue!=0});
};

MuPlot.update = function(el, mu, state, dispatcher){
    console.log("UPDATING MU");
    
    this.mu = mu;
    this._update_points();

    var selected_tip = state.selected_tip
        this.points.map(function(d){
      if (selected_tip && d.name == selected_tip){
          d.selected = true;
      }else{
          d.selected = false;
      }
    });
    var scales = this._scales(el);
    this._draw_axis(el, scales)
    this._draw_points(el, scales, dispatcher)
    
    // selected node
};

MuPlot._scales = function(el){

  
  var width = el.offsetWidth;
  var height = el.offsetHeight;
  var xs = this.points.map(function(d){return d.x});
  var ys = this.points.map(function(d){return d.y});
  var x = d3.scale.linear()
    .domain([d3.min(xs) , d3.max(xs)])
    .range([this.left_padding, width-this.padding]);

  var y = d3.scale.linear()
      .domain([d3.max(ys), d3.min(ys)])
      .range([this.padding,height-this.padding])
  return {x: x, y: y};

};

MuPlot._tipRadius = function(d){
    return d.selected ? 8.0 : 6.0;
};

MuPlot._draw_axis = function(el, scales){
    
    var width = el.offsetWidth;
    var height = el.offsetHeight;

    var xAxis = d3.svg.axis()
        .scale(scales.x)
        .orient("bottom")
        .ticks(5)
    
    var yAxis = d3.svg.axis()
        .scale(scales.y)
        .orient("left")
        .ticks(10)

    // function for the y grid lines
    function make_x_axis() {
    return d3.svg.axis()
      .scale(scales.x)
      .orient("bottom")
      .ticks(5)
    }

    function make_y_axis() {
    return d3.svg.axis()
      .scale(scales.y)
      .orient("left")
      .ticks(10)
    }

    var svg = d3.select(el).append("svg")
        .attr("class", "d3_mu_axis") 
        .attr("width", width)
        .attr("height", height)

    svg.append("g")
        .attr("class", "d3_mu_axis")
        .attr("transform", "translate(0," + (height - this.padding) + ")")
        .call(xAxis)

    svg.append("g")
        .attr("class", "d3_mu_grid")
        .attr("transform", "translate(0," + (this.padding) + ")")
        .call(make_x_axis()
            .tickSize(height-2*this.padding, 0, 0)
            .tickFormat("")
            )
    
    svg.append("g")
        .attr("transform", "translate(" + (this.left_padding) + ",0)")
        .call(yAxis);
    
    svg.append("g")
        .attr("class", "d3_mu_grid")
        .attr("transform", "translate(" + (width - this.padding) + ",0)")
        .call(make_y_axis()
            .tickSize(width-this.padding-this.left_padding, 0, 0)
            .tickFormat("")
            )
};

MuPlot._draw_points = function(el, scales, dispatcher){

    console.log("DRAW POINTS")
    var g = d3.select(el).selectAll('.d3-mu-points');

    var tip = g.selectAll('.d3-mu-point')
        .data(this.points);

    tip.enter()
      .append("circle")
      .attr("class", "d3-mu-point")
      .attr("id", function(d){return (d.name).replace(/\//g, "")})
    
    tip
      .attr("cx", function(d){return scales.x(d.x)})
      .attr("cy", function(d){return scales.y(d.y)})
      .attr("r", this._tipRadius)
      .style("fill", "#4D92BF")
      .on('mouseover', function(d) {
          dispatcher.emit('point:point_mouseover', d);
      })
      .on('mouseout', function(d) {
          dispatcher.emit('point:point_mouseout', d);
      })
};

export default MuPlot;