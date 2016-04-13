var EventEmitter = require('events').EventEmitter;
var d3 = require('d3');
var Globals = require('./globals.js')
var ANIMATION_DURATION = Globals.ANIMATION_DURATION;


var getLeafNodes = function(leafNodes, obj){
        if(obj.children){
            obj.children.forEach(function(child){getLeafNodes(leafNodes,child)});
        } else{
            leafNodes.push(obj);
        }
};

var MuPlot = {};

MuPlot.mu = [];

MuPlot.old_state = {};

MuPlot.padding_bottom = 80;

MuPlot.padding_text = 20;

MuPlot.padding_top = 10;

MuPlot.padding_right = 10;

MuPlot.padding_left = 120;

MuPlot.regression = {};

MuPlot.create = function(el, props, state){
    console.log("CREATING MU PLOT")
    var width = el.offsetWidth
    var height = el.offsetHeight

    var svg = d3.select(el).append('svg')
      .attr('class', 'd3_mu')
      .attr('width',  width)
      .attr('height', height);
    
    svg.append('g')
      .attr('class', 'd3_mu_axis')

    svg.append('g')
      .attr('class', 'd3_mu_points')

    var dispatcher = new EventEmitter();
    //this.update(el, props.root, state, dispatcher);
    return dispatcher;

};

MuPlot._set_points_from_root = function(dispatcher){

  console.log("ROOT updates LH...");
  if (!this.tree) {return ;}
  var tip_lhs = [];
  getLeafNodes(tip_lhs, this.tree);
  
 
  this.points = tip_lhs
      .map(function(d){
        return ({
          name: d.name,
          x: d.numdate,
          y: d.xvalue
        });
      }).filter(function(d){
        return (
          (d.numdate != 0.0) && (d.xvalue != 0)
          );
      });

  this._update_lin_regression(dispatcher);

};

MuPlot._update_lin_regression = function(dispatcher){
    
    console.log("Updateting linear regression for MU plot...");
    
    var n = this.points.length;
    
    if (n == 0) {
      this.regression = {};
      return;
    }
    
    var sum_x =  d3.sum(this.points.map(function(d){return d.x}));
    var sum_y =  d3.sum(this.points.map(function(d){return d.y}));
    var sum_xy = d3.sum(this.points.map(function(d){return d.x * d.y}));
    var sum_xx = d3.sum(this.points.map(function(d){return d.x * d.x}));
    var sum_yy = d3.sum(this.points.map(function(d){return d.y * d.y}));

    var slope = (n * sum_xy - sum_x * sum_y) / (n*sum_xx - sum_x * sum_x);
    var intercept =  (sum_y - slope * sum_x)/n;
    var r2 = Math.pow((n*sum_xy - sum_x*sum_y)/Math.sqrt((n*sum_xx-sum_x*sum_x)*(n*sum_yy-sum_y*sum_y)),2);

    this.regression = {
      'slope' :  slope,
      'intercept' :  intercept,
      'r2' :  r2
    };
    console.log("Molecular clock plot linear regression changed, emitting signal.")
    dispatcher.emit('mol_clock:regression_changed', this.regression);

};

MuPlot._draw_regression = function(el, scales) {

  if (!this.regression) return;
  if (!this.points ) return;
  console.log("MuPlot updating the molecular clock linear regression...")
  var max_x = d3.max(this.points.map(function(d){return d.x}));
  var min_x = d3.min(this.points.map(function(d){return d.x}));

  var max_y = this.regression.slope * max_x + this.regression.intercept
  var min_y = this.regression.slope * min_x + this.regression.intercept

  var svg = d3.select(el).select('.d3_mu').append('line')
    .attr('class', 'd3_mu_regression')
    .attr('x1', function(d){return scales.x(min_x)})
    .attr('y1', function(d){return scales.y(min_y)})
    .attr('x2', function(d){return scales.x(max_x)})
    .attr('y2', function(d){return scales.y(max_y)})
    .style("stroke", "#4D92BF")
    .style("stroke-width", '2px')

  console.log(svg);

};

MuPlot.update = function(el, root, state, dispatcher){
    
    console.log("UPDATING MU");

    if (this.tree != root){
      // update all points 
      console.log("MuPlot detected Tree changes, recreating the plot...")
      this.tree = root;
      this._set_points_from_root(dispatcher);
      this._update_lin_regression(dispatcher);
      var scales = this._scales(el);
      this._draw_axis(el, scales)
      this._draw_points(el, scales, dispatcher)
      this._draw_regression(el, scales)

    }

    if(this.old_state.selected_tip != state.selected_tip){
        console.log("MU: tip selection changed.")
        var selected_tip = state.selected_tip
            this.points.map(function(d){
          if (selected_tip && d.name == selected_tip){
              d.selected = true;
          }else{
              d.selected = false;
          }
        });
        this._refresh_selected_tip(el);
    }
        
    this.old_state = state;
    // selected node

};

MuPlot._refresh_selected_tip = function(el){
    
    var g = d3.select(el).selectAll('.d3_mu_points');
    var tip = g.selectAll('.d3_mu_point')
    tip.attr("r", this._tipRadius)

};

MuPlot._scales = function(el){

  
  var width = el.offsetWidth;
  var height = el.offsetHeight;
  var xs = this.points.map(function(d){return d.x});
  var ys = this.points.map(function(d){return d.y});
  var x = d3.scale.linear()
    .domain([d3.min(xs) , d3.max(xs)])
    .range([this.padding_left, width - this.padding_right]);

  var y = d3.scale.linear()
      .domain([d3.max(ys), d3.min(ys)])
      .range([this.padding_top,height-this.padding_bottom])
  return {x: x, y: y};

};

MuPlot._tipRadius = function(d){
    return d.selected ? 10.0 : 6.0;
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

    var svg = d3.select(el).select('.d3_mu_axis')
        

    svg.append("g")
        .attr("class", "d3_mu_x_axis")
        .attr("transform", "translate(0," + (height -  this.padding_bottom) + ")")
        .call(xAxis)

    svg.append("text")      // text label for the x axis
        .attr("x", width / 2 )
        .attr("y", height - this.padding_text )
        .style("text-anchor", "middle")
        .text("Sampling date");


   svg.append("g")
        .attr("class", "d3_mu_x_grid")
        .attr("transform", "translate(0," +  ( + this.padding_text) + ")")
        .call(make_x_axis()
            .tickSize(height-this.padding_bottom, 0, 0)
            .tickFormat("")
            )
    
    svg.append("g")
        .attr("class", "d3_mu_y_axis")
        .attr("transform", "translate(" + (this.padding_left) + ",0)")
        .call(yAxis);
    
    svg.append("g")
        .attr("class", "d3_mu_y_grid")
        .attr("transform", "translate(" + ( width) + ",0)")
        .call(make_y_axis()
            .tickSize(width-this.padding_left, 0, 0)
            .tickFormat("")
            )
    svg.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", this.padding_text)
        .attr("x",  - (height - this.padding_bottom)/2)
        .attr("dy", "1em")
        .style("text-anchor", "middle")
        .text("Distance to root");
};

MuPlot._draw_points = function(el, scales, dispatcher){

    console.log("MU DRAW POINTS...")
    var g = d3.select(el).selectAll('.d3_mu_points');

    
    var tip = g.selectAll('.d3_mu_point')
        .data(this.points);

    tip.enter()
      .append("circle")
      .attr("class", "d3_mu_point")
      .attr("id", function(d){
        return "NAME" //(d.name).replace(/\//g, "")
      })
    
    tip
      .attr("cx", function(d){return scales.x(d.x)})
      .attr("cy", function(d){return scales.y(d.y)})
      .attr("r", this._tipRadius)
      .style("fill", "#4D92BF")
      .style('stroke',"blue")
      .style('stroke-width',"2")

      .on('mouseover', function(d) {
          dispatcher.emit('point:point_mouseover', d);
      })
      .on('mouseout', function(d) {
          dispatcher.emit('point:point_mouseout', d);
      })
};

export default MuPlot;