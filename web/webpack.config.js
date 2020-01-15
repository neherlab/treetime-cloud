var webpack = require('webpack');

module.exports = {
    entry: {
        welcome_1: "./js/components/welcome_1.js",
        welcome_anc: "./js/components/welcome_ancestral_reconstruction.js",
        welcome_treetime: "./js/components/welcome_treetime.js",
        progress_ancestral: "./js/components/progress_ancestral.js",
        progress_treetime: "./js/components/progress_treetime.js",
        results_treetime: "./js/components/results_treetime.js",
        about: "./js/components/about.js",

    },

    output: {
        path: __dirname + '/static/js',
        filename: "[name].js"
    },

    module: {
        loaders: [
            { test: /\.js?$/, loaders: ['react-hot', 'babel'], exclude: /node_modules/ },
            { test: /\.js$/,  exclude: /node_modules/, loader: 'babel-loader'},
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
