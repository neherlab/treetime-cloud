var webpack = require('webpack');  

module.exports = {  
    entry: {
        welcome: "./js/components/welcome.js",
        progress: "./js/components/progress.js", 
        results: "./js/components/results.js", 
        terms: "./js/components/terms.js", 

    },

    output: {
        path: __dirname + '/static/js',
        filename: "[name].js"
    },
    
    module: {
        loaders: [
            { test: /\.js?$/, loaders: ['react-hot', 'babel'], exclude: /node_modules/ },
            { test: /\.js$/, exclude: /node_modules/, loader: 'babel-loader'},
            { test: /\.css$/, loader: "style!css" }
        ]
    },
    plugins: [
      new webpack.NoErrorsPlugin()
    ],

    node: {
        net: 'empty',
        tls: 'empty'
    }
};
