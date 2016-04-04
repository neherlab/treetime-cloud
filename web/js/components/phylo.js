module.exports = {
treeSetUp: function (start, end){
        calcFullTipCounts(rootNode);
        rootNode.branch_length=0.001;
        rootNode.tbranch_length=0.001;
        calcBranchLength(rootNode);
        nDisplayTips = displayRoot.fullTipCount;

        xValues = nodes.map(function(d) {return +d.xvalue;});
        tValues = nodes.map(function(d) {return +d.tvalue;});
        yValues = nodes.map(function(d) {return +d.yvalue;});

        if (timetree){
            currentXValues = tValues;
        }else{
            currentXValues = xValues;
        }

        xScale = d3.scale.linear()
            .domain([d3.min(currentXValues), d3.max(currentXValues)]);
        yScale = d3.scale.linear()
            .domain([d3.min(yValues), d3.max(yValues)]);
        drawGrid();

        canvas.selectAll(".link")
            .data(links)
            .enter().append("polyline")
            .attr("class", "link");

        canvas.selectAll(".tip")
            .data(tips)
            .enter()
            .append("circle")
            .attr("class", "tip")
            .attr("id", function(d) { return (d.strain).replace(/\//g, ""); });

    }

}
