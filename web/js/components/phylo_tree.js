var EventEmitter = require('events').EventEmitter;
var d3 = require('d3');
var d3tip = require('d3-tip');

var Globals = require('./globals.js')
var colors = Globals.colors;
var ANIMATION_DURATION = Globals.ANIMATION_DURATION;

var minimumAttribute = function (node, attr, min) {

  if (typeof node.children != "undefined") {
    for (var i=0, c=node.children.length; i<c; i++) {
      min = minimumAttribute(node.children[i], attr, min);
    }
  }
  else {
    if (attr(node) < min) {
      min = attr(node);
    }
  }
  return min;

};

var maximumAttribute = function (node, attr, max) {
  if (typeof node.children != "undefined") {
    for (var i=0, c=node.children.length; i<c; i++) {
      max = maximumAttribute(node.children[i], attr, max);
    }
  }
  else {
    if (attr(node) > max) {
      max = attr(node);
    }
  }
  return max;
};

var virusTooltip = d3tip()
  .direction('w')
  .attr('class', 'd3-tip')
  .offset([0, -12])
  .html(function(d){
    var string = "";
    string += "<h4>" + d.strain +"</h4>"
    "<div class=\"smallspacer\"></div>";
    string += "<div class=\"smallnote\">";

    string += "Numdate: " + d.numdate + " <br/>";
    string += "Time since Tmrca: " + d.tvalue + "<br/>";
    string += "Country: " + d.country + "<br/>" ;
    string += "Region: " + d.region + "<br/>" ;
    string += "</div>";

    return (string)});

var linkTooltip = d3tip()
  .direction(function(d){if (d.target.x > 500) {return 's';} else{return 's';};})
  .attr('class', 'd3-tip')
  .offset(function(d){if (d.source.x > 500) {return [0, -3]} else{return [0,3]};})
  .html(function(d) {
    //console.log(d)

    var string = ""
    if (typeof d.frequency != "undefined") {
      string += "Frequency: " + (100 * d.frequency).toFixed(1) + "%"
    }

    string += "<div class=\"smallspacer\"></div>";
    string += "<div class=\"smallnote\">";

    if (false) {//((typeof d.aa_muts !="undefined")&&(mutType=='aa')){
      var ncount = 0;
      for (tmp_gene in d.aa_muts) {ncount+=d.aa_muts[tmp_gene].length;}
      if (ncount) {string += "<b>Mutations:</b><ul>";}
      for (tmp_gene in d.aa_muts){
        if (d.aa_muts[tmp_gene].length){
          string+="<li>"+tmp_gene+":</b> "+d.aa_muts[tmp_gene].replace(/,/g, ', ') + "</li>";
        }
      }
    }

    else if ((typeof d.target.muts !="undefined")&&(d.target.muts.length)){
      var tmp_muts = d.target.muts.split(',');
      var nmuts = tmp_muts.length;
      tmp_muts = tmp_muts.slice(0,Math.min(10, nmuts))
      string += "<li>"+tmp_muts.join(', ');
      if (nmuts>10) {string+=' + '+ (nmuts-10) + ' more';}
      string += "</li>";
    }
    string += "</ul>";
    string += "click to zoom into clade"
    string += "</div>";
    return string;

  });


var PhyloTree = {};

PhyloTree.state = {};

PhyloTree.padding_bottom = 80;
PhyloTree.padding_text = 20;
PhyloTree.padding_top = 10;
PhyloTree.padding_right = 25;
PhyloTree.padding_left = 20;

PhyloTree.create = function(el, props, state){

  if (!props.root || typeof props.root == 'undefined') return null;
  //console.log("CREATING tree");
  this.width  = el.offsetWidth;
  this.height = el.offsetHeight;
  this.svg = d3.select(el).append('svg')
      .attr('class', 'd3-phylotree')
      .attr('width', this.width)
      .attr('height', this.height)


  this.svg.append('g')
      .attr('class', 'd3-tree_axis')
  this.svg.append('g')
      .attr('class', 'd3-links');
  this.svg.append('g')
      .attr('class', 'd3-tips');

  this.svg.call(virusTooltip);
  this.svg.call(linkTooltip);

  var dispatcher = new EventEmitter();
  this._create_data(props);
  this.state = state;

  this._update_node_colors(this.state);
  var scales = this._scales(el, state);
  this._drawTips(el,scales, dispatcher);
  this._drawLinks(el,scales,dispatcher);
  return dispatcher;

};

PhyloTree._create_data = function(props){

    ////console.log(props.root)
    this.root = props.root;
    console.log("Creating tree data...")
    console.log(this.root)
    this.layout = d3.layout.tree();
    this.vis_nodes = this.layout.nodes(this.root);
    this.vis_links = this.layout.links(this.vis_nodes);
    this.vis_links.push({"source":this.root, "target":this.root});
    this.vis_tips = []
    this.gatherTips(this.root, this.vis_tips);
    this.vis_tips.map(function(d){d.selected=false;});
    console.log(this.vis_tips)
    this._update_vis(this.root);
};

PhyloTree._update_vis = function(vis_root){

  this.vis_root = vis_root;
  //this.vis_nodes = this.layout.nodes(this.vis_root);
  //this.vis_links = this.layout.links(this.vis_nodes);
  //this.vis_links.push({"source":this.vis_root, "target":this.vis_root});
  //this.vis_tips = this.gatherTips(this.vis_root, []);
  //this.vis_tips.map(function(d){d.selected=false;});
};

PhyloTree.save_svg = function(el){

  var html = d3.select(el).select(".d3-phylotree")
        .attr("title", "test2")
        .attr("version", 1.1)
        .attr("xmlns", "http://www.w3.org/2000/svg")
        .node().parentNode.innerHTML;

  return html;
};

PhyloTree.resetLayout = function(el, dispatcher){
  this._zoom(this.root, el, dispatcher)
};

PhyloTree.gatherTips = function (node, tips) {

    if (typeof node.children != "undefined") {
        for (var i=0, c=node.children.length; i<c; i++) {
            this.gatherTips(node.children[i], tips);
        }
    }
    else {
        tips.push(node);
    }
    return tips;
};

PhyloTree.update = function(el, state, dispatcher) {
  
  if (this.width != el.offsetWidth || this.height != el.offsetHeight){
      this._hide_axis(el, state, dispatcher);
      this._update_scales(el, state, dispatcher);
      this.width = el.offsetWidth;
      this.height = el.offsetHeight;
      //this._update_axis(el, state, scales, dispatcher)
  }

  if (this.state.cscale != state.cscale || this.state.cvalue != state.cvalue){
    this._update_node_colors(state);
  }

  if (this.state.selected_tip != state.selected_tip){
      this._update_selected(el, state);
  }

  if (this.state.xUnit != state.xUnit){

    this._x_pos = function(d){
      return d[state.xUnit];
    };
    this._hide_axis(el, state, dispatcher)
    this._update_scales(el, state, dispatcher);

  }

  // update the tree state
  this.state = state;
};


PhyloTree._update_selected = function(el,state){

    this.vis_tips.map(function(d){
          if (state.selected_tip && d.strain == state.selected_tip){
              d.selected = true;
          }else{
              d.selected = false;
          }
    });

    var g = d3.select(el).selectAll('.d3-tips');
    var tip = g.selectAll('.d3-tip')
    tip.attr("r", this._tipRadius)
};

PhyloTree._update_scales = function(el, state, dispatcher){
    var scales = this._scales(el, state);
    //this._drawLinks(el, scales, dispatcher);
    this._moveLinks(el, scales, dispatcher);
    //this._drawTips(el, scales, dispatcher);
    this._moveTips(el,scales, dispatcher);
    this._update_axis (el, state, scales, dispatcher);
};

PhyloTree._scales = function(el, state) {

  var width = el.offsetWidth;
  var height = el.offsetHeight;
  this.svg.attr('width', width)
          .attr('height', height);
  var xpos = this._x_pos;
  var ypos = this._y_pos;

  var xMin = minimumAttribute(this.vis_root, xpos, xpos(this.vis_root));
  var xMax = maximumAttribute(this.vis_root, xpos, xpos(this.vis_root));

  console.log(xMin, xMax)
  var yMin = minimumAttribute(this.vis_root, ypos, ypos(this.vis_root));
  var yMax = maximumAttribute(this.vis_root, ypos, ypos(this.vis_root));

  var x = d3.scale.linear()
      .domain([xMin, xMax])
      .range([this.padding_left, width-this.padding_right]);

  var y = d3.scale.linear()
      .domain([yMin, yMax])
      .range([this.padding_top,height-this.padding_bottom])

  return {x: x, y: y};
};

/*
Get x position for the tree node
*/
PhyloTree._x_pos = function (d){
  return d["xvalue"];
};

/*
Get y position for the tree node
*/
PhyloTree._y_pos = function(d){
  return d["yvalue"];
};

PhyloTree._update_node_colors = function(state){

    this.vis_tips.forEach(function (d) {
        if (typeof(d) =='undefined') {return;}
        var cval = state.cvalue(d);
        d.color = state.cscale.get_color(cval);
    });

    var  _drawVirusToolTip = this._drawVirusToolTip
    var _hideVirusToolTip = this._hideVirusToolTip

    var g = d3.selectAll('.d3-tips');
    var tip = g.selectAll('.d3-tip')
    tip
      .style("fill", this._tipFillColor)
      .style("stroke", this._tipStrokeColor)
};

PhyloTree.get_color_scale = function(){
};

PhyloTree._tipFillColor = function(d) {

    var c =  d3.rgb(d.color).brighter([0.65]);
    return c;
};

PhyloTree._tipStrokeColor = function(d)  {return d.color;}

PhyloTree._branchStrokeColor = function (d) {return "#BBBBBB";}

PhyloTree._tipVisibility = function (d) { return d.current?"visible":"hidden";}

PhyloTree._branchStrokeWidth = function(d) {return 3;}

PhyloTree._tipRadius = function(d) {return d.selected ? 16.0 : 10.0;}

PhyloTree._drawTips = function(el, scales, dispatcher) {

    if (typeof this.vis_tips ==' undefined' || this.vis_tips.length < 2) return;
    var  _drawVirusToolTip = this._drawVirusToolTip
    var _hideVirusToolTip = this._hideVirusToolTip
    var xpos = this._x_pos;
    var ypos = this._y_pos;

    var g = d3.select(el).selectAll('.d3-tips');
    var tip = g.selectAll('.d3-tip')
        .data(this.vis_tips);

    console.log("VIS TIPS")
    console.log(this.vis_tips)
    console.log(this.vis_tips.map(function(d){return d.strain}))
    console.log(d3.sum(this.vis_tips.map(function(d){if (typeof d.strain == 'undefined') return 1; else return 0;}) ) )
    tip.enter()
      .append("circle")
      .attr("class", "d3-tip")
      .attr("id", function(d) { return (d.strain).replace(/\//g, ""); })
      .attr("cx", function(d) {

            return scales.x(xpos(d))
            //if (d.old_x){
            //    return d.old_x;
            //}
            //return d.x;
        })
    tip
      .attr("r", this._tipRadius)
      .attr("cy", function(d) {
        return scales.y(ypos(d))
      })
      .style("fill", this._tipFillColor)
      .style("stroke", this._tipStrokeColor)
      .on('mouseover', function(d) {
          virusTooltip.show(d);
          dispatcher.emit('point:tip_mouseover', d);
      })
      .on('mouseout', function(d) {
          //_hideVirusToolTip();
          virusTooltip.hide();
          dispatcher.emit('point:tip_mouseout', d);
      });
};

PhyloTree._moveTips = function(el, scales, dispatcher) {

    var xpos = this._x_pos;
    var ypos = this._y_pos;

    var g = d3.select(el).selectAll('.d3-tips');
    var tip = g.selectAll('.d3-tip')
    tip.transition()
      .duration(ANIMATION_DURATION)
      .attr('cx', function(d) {
        return scales.x(xpos(d))
      })
      .attr('cy', function(d){
        return scales.y(ypos(d))
      });
};

PhyloTree._drawLinks = function(el, scales, dispatcher){

    var bpoints = this._branchPoints;

    var g = d3.select(el).selectAll('.d3-links');
    var link = g.selectAll('.d3-link')
        .data(this.vis_links);



    var zoom = this._zoom;

    link.enter()
        .append("polyline")
        .attr("class", "d3-link")
        .attr("points", function(d){return bpoints(d, scales);})
        .style("stroke-width", 3)
        .style("stroke", "gray")
        .style("stroke-linejoin", "round")
        .style("cursor", "pointer")
        .on("mouseover", function(d){
            linkTooltip.show(d)
        })
        .on("mouseout", function(d){
            linkTooltip.hide()
        })
        .on("click", function(d){
          console.log("MOUSECLICK")
          zoom(d.target, el, dispatcher);
        });

    // /link.transition()
    // /        .attr("points", this._branchPoints)
    // /        .duration(ANIMATION_DURATION)
};

PhyloTree._moveLinks = function(el, scales, dispatcher){

    var bpoints = this._branchPoints;

    var g = d3.select(el).selectAll('.d3-links');
    var link = g.selectAll('.d3-link')

    link.transition()
            .attr("points",function(d){return bpoints(d, scales);})
            .duration(ANIMATION_DURATION)
};

PhyloTree._zoom = function (d, el, dispatcher){

  PhyloTree._hide_axis(el, PhyloTree.state, dispatcher);
  PhyloTree._update_vis(d);
  PhyloTree._update_scales(el, PhyloTree.state, dispatcher);

};

PhyloTree._branchPoints = function(d, scales) {


    var xpos = PhyloTree._x_pos;
    var ypos = PhyloTree._y_pos;

    var tmp =   scales.x(xpos(d.source)).toString() + "," + scales.y(ypos(d.target)).toString() + " "
              + scales.x(xpos(d.target)).toString() + "," + scales.y(ypos(d.target)).toString();

    if (typeof d.target.children != "undefined"){

        var child_ys = d.target.children.map(function (x){
          return scales.y(ypos(x));
        });

        tmp+= " "+ scales.x(xpos(d.target)).toString()+","+d3.min(child_ys).toString() + " "
                 + scales.x(xpos(d.target)).toString()+","+d3.max(child_ys).toString();

    }
    return tmp;
},

PhyloTree._branchPoints_old = function(d) {

    var tmp =   d.source.old_x.toString() + "," + d.target.old_y.toString() + " "
              + d.target.old_x.toString() + "," + d.target.old_y.toString();
    if (typeof d.target.children != "undefined"){
        var child_ys = d.target.children.map(function (x){return x.old_y;});
        tmp+= " "+ d.target.old_x.toString()+","+d3.min(child_ys).toString() + " "
                 + d.target.old_x.toString()+","+d3.max(child_ys).toString();
    }
    return tmp;
},

PhyloTree._update_axis = function(el, state, scales, dispatcher){

    if (state.xUnit == "numdate"){
        this._draw_axis(el, state, scales, dispatcher);
    }else{
        this._hide_axis(el, state, dispatcher);
    }
};

PhyloTree._draw_axis = function(el, state, scales, dispatcher){
    console.log("DRAW axis called")
    var width = el.offsetWidth;
    var height = el.offsetHeight;

    var xAxis = d3.svg.axis()
        .scale(scales.x)
        .orient("bottom")
        .ticks(10)

    // function for the x grid lines
    function make_x_axis() {
    return d3.svg.axis()
      .scale(scales.x)
      .orient("bottom")
      .ticks(10)
    }

    var svg = d3.select(el).select('.d3-tree_axis')

    svg.append("g")
        .attr("class", "d3_tree_x_axis")
        .attr("transform", "translate(0," + (height -  this.padding_bottom) + ")")
        .call(xAxis)

    svg.append("text")      // text label for the x axis
        .attr("x", width / 2 )
        .attr("y", height - this.padding_text )
        .style("text-anchor", "middle")
        .text("Date");


   svg.append("g")
        .attr("class", "d3_tree_x_grid")
        .attr("transform", "translate(0," +  ( + 0) + ")")
        .call(make_x_axis()
            .tickSize(height-this.padding_bottom, 0, 0)
            .tickFormat("")
            )

    // var bnd = scales.x.domain()
    // var yticks_xs = d3.range(bnd[0], bnd[1], 1);
    // var y_bnds = scales.y.range()
    // var glines = yticks_xs.map(function(d){return ({x1:d, y1:y_bnds[0], x2:d, y2:y_bnds[1]});})

    // var gridLine = g.selectAll('.d3-tree-axis-line')
    //     .data(glines);

    // gridLine.enter()
    //     .append('line')
    //     .attr("class","d3-tree-axis-line")
    //     .attr("x1", function(d){return scales.x(d.x1);})
    //     .attr("x2", function(d){return scales.x(d.x2);})
    //     .attr("y1", function(d){return d.y1;})
    //     .attr("y2", function(d){return d.y2;})
    //     .style("stroke-width", 1)
    //     .style("stroke", "#DDDDDD");
};

PhyloTree._hide_axis = function (el, state, dispatcher){

    
    var g = d3.select(el).select('.d3-tree_axis').selectAll("*");
    g.remove();
};

PhyloTree.destroy = function(el) {
};

export default PhyloTree;
